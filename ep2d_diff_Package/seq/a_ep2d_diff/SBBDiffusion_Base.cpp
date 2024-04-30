//-----------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2021. All Rights Reserved.
//-----------------------------------------------------------------------------


#include "MrProtSrv/Domain/CoreNative/MeasAqcDefs.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/KSpace/MrKSpace.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrTXSpec.h"
#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"

#include "MrVista/Ice/IceApplicationFunctors/IceDiffusion/DiffusionDef.h" // DIFFUSIONINTSTORAGE_...
#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Base.h"
#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"
#include "MrImagingFW/libSeqUTIF/libsequt.h"                          // SeqUT
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include <cmath>

#ifdef WIP
#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_BipolarPlus.h"
#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_StejskalPlus.h"
#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_STEAM.h"
#endif //   #ifdef WIP

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif
using namespace SEQ_NAMESPACE;

// ===========================================================================
/*!
\class  SBBDiffusion_Base

\brief  This virtual base class provides members common to all diffusion modes.

SBBDiffusion_Base is the virtual base class of this sequence
building block familiy which offers the possibility to easily implement
diffusion weighting into a sequence. The member functions of this
base class are not intended for direct use, but must be overloaded
by derived classes.

\author Michael.Zwanger@med.siemens.de
*/
// ============================================================================

// ***************************************************************************
///	The constructor initializes some members variables.
/**   In detail, the following variables are initialized:
- lDebugLevel
- allow maximum possible amplitudes for the spoiler gradients
- \ref m_dMaxAmpl and \ref m_dMinRiseTime are initialized depending
on the installed GPA and GX.
*/
SBBDiffusion_Base::SBBDiffusion_Base(SBBList* pSBBList)
    : SeqBuildBlock(pSBBList)
{
    setIdent("SBBDiffusion_Base");

    // ----------------------------------------------------
    // Configure SBB-specific dynamic adjustment properties
    // ----------------------------------------------------
    SeqBuildBlock::seteSliceAdjOptimizationMode(SLICEADJ::HOLD_OPTIMIZATION); // the diffusion SBB shall only hold the results from the pulses
    SeqBuildBlock::setsSliceAdjParametersRequestedBySBB(SLICEADJ::ADJNONE);   // all SliceAdj functionality is in the refocusing RF SBBs

    // Initialize diffusion directions for adjustment scans
    if(!m_AdjDidi.prepInternal(0, '\0', false))
    {
        SEQ_TRACE_ERROR.print("FATAL ERROR: AdjDidi preparation failed!");
        // Constructor failed: throw an exception
        throw "SBBDiffusion_Base::SBBDiffusion_Base: ERROR: AdjDidi.prepInternal() failed, this should never happen";
    }

    // Set maximum possible(!) gradient amplitudes (used during check())
    m_DSp1.setMaxMagnitude(SysProperties::getGradMaxAmplAbsolute());
    m_DSr1.setMaxMagnitude(SysProperties::getGradMaxAmplAbsolute());
    m_DSs1.setMaxMagnitude(SysProperties::getGradMaxAmplAbsolute());

    m_DSs1.setAxis(SEQ::AXIS_SLICE);
    m_DSp1.setAxis(SEQ::AXIS_PHASE);
    m_DSr1.setAxis(SEQ::AXIS_READOUT);


    // --------------------------------------------------------------------------
    // Select max. gradient amplitude
    // Take care: appropriate GPA models should exist.
    // WARNING: These values are only preliminary, as usually the method
    //          CalcMaximumAmplitude is called.
    //
    // During the sequence check method the gradient amplitude is checked for a
    // maximum value which is GradMaxAmplAbsolute * (1.0 - (fECCMargin / 100.0))
    // where fECCMargin = 5.
    // For the diffusion SBBs this sets 0.95 * GradMaxAmplAbsolute as absolute
    // upper limit for all gradients. The SBBs have to take care that this value
    // never is exceeded (including all libBalance dependencies and all correction
    // mechanisms like Maxwell corrections).
    // The bipolar SBB uses a Maxwell correction and has to take care of including 
    // a margin for the additional moments.
    // --------------------------------------------------------------------------
    m_dMaxAmpl = SysProperties::getGradMaxAmplAbsolute() * 0.95;

    const auto dMasterAmplitude = m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/max_grad_ampl", 0.0);
    if(dMasterAmplitude > 0.0 && dMasterAmplitude < 1.0)
    {
        m_dMaxAmpl = dMasterAmplitude;
    }

    // --------------------------------------------------------------------------
    // Select min. gradient risetime.
    //
    // If we are in the GSWD-mode we do not know, if the Diffusion-SBB or the
    // EPI-ReadOut-SBB causes the lookahead-stimulation level to exceed the allowed
    // limits. This is a problem we can't solve at the moment.
    // There are now the following scenarios:
    // (1) We could use the GSWD-risetime in the protocol both for the EPI readout
    // and the diffusion gradients. Actually we do NOT do this because, if the Diffusion
    // SBB causes the stimulation, this will cause a needless loss of performance
    // for the EPI-ReadOut which will increase the echospacing and therefore has
    // a major impact on image quality. This is much worse than an increase in TE
    // due to longer diffusion gradients.
    // (2) We restrict the minimum allowed risetime for the Diffusion-SBB
    // to a safe CONSTANT value, *if* the GSWD is searching an optimum risetime for the
    // sequence with which a measurement without stimulation is possible. This may
    // cause a needless loss of performance for the EPIReadout, if the SBBDiffusion
    // causes the stimulation. This is what we did until VA21B.
    // (3) For all cases (i.e. b values and gradient systems), we use a safe
    // slew rate for SBBDiffusion. This is what we do;-)
    //
    // --------------------------------------------------------------------------


    // This is an empiric value which is safe with repect to GSWD.
    // You must be REALLY sure that it is safe - else the sequence will abort.
    double dRiseTime = 16.0;

    // For Prisma, a longer risetime is required to avoid running into GSWD limitations.
    // Note that this value might not be completely safe: if all three gradients get
    // ramped from +Gmax to -Gmax simultaneously, a rise time of ~22us is required
    // in order to stay below stimulation limits. If the GSWD does not suggest a 
    // reasonable protocol modification in these rare cases, manually increasing TE
    // will help.
    if(SysProperties::isPrisma())
    {
        dRiseTime = 20.0;
    }

    // The Allegra has a head coil not stimulating as much as a body coil.
    // So we can use an empiric(!) value for faster performance
    if(SysProperties::isAllegra())
    {
        dRiseTime = 10.0;
    }

    // Ensure staying within hardware limits
    dRiseTime = std::max(dRiseTime, static_cast<double>(SysProperties::getGradMinRiseTime(SEQ::GRAD_NORMAL)));

    m_dMinRiseTime = dRiseTime;

    // This is an empiric value which is safe with respect to GSWD.
    // You must be REALLY sure that it is safe - else some protocols will abort.
    m_lStimoDelayus = 300;

    // m_dMaxAmpl contains the absolute maximum amplitude allowed for diffusion gradients and
    // must not be changed.
    // m_dAmpl denotes the maximum amplitude allowed for the current timing. If a GPA model
    // is available, this value will be adapted accordingly.
    m_dAmpl = m_dMaxAmpl;

#ifdef QUIETDWI
    // --------------------------------------------------------------------------
    // Select max. gradient amplitude for offset Diffusion gradients in qDWI
    // WARNING: These values are only preliminary
    // --------------------------------------------------------------------------
    m_DGoffp.setMaxMagnitude(20);
    m_DGoffr.setMaxMagnitude(20);
    m_DGoffs.setMaxMagnitude(20);

    // --------------------------------------------------------------------------
    // Select min. gradient risetime for offset Diffusion gradients in qDWI
    // WARNING: These values are only preliminary
    m_DGoffs.setMinRiseTime(1000.0 / 25.0);
    m_DGoffp.setMinRiseTime(1000.0 / 25.0);
    m_DGoffr.setMinRiseTime(1000.0 / 25.0);
#endif // QUIETDWI

    // Reset preparation status
    resetPrepared();
}


// ===========================================================================
/// Basic preparations and identification of required TE
// ===========================================================================
bool SBBDiffusion_Base::prepSBB(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo)
// ===========================================================================
{
    // Reset preparation status
    resetPrepared();

    // Set a default error exit code in case we should bail out unexpectedly
    // (This must be done before the first return statement)
    setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);

    // -------------------------------------------------------------------------------
    // Prohibit use of free diffusion mode if
    // a) the directory containing external diffusion vector sets does not exist AND
    // b) no free diffusion vector set is stored in the protocol
    // -------------------------------------------------------------------------------
    const bool bFreeModeAvailable = m_Didi.isExtDiffDirPresent() || (rMrProt.diffusion().getsFreeDiffusionData().getlDiffDirections() != 0);

    if((rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_FREE) && !bFreeModeAvailable)
    {
        // Free mode is not possible
        return false;
    }

    // ------------------------------------
    // Extract relevant protocol parameters
    // ------------------------------------
    if(!setParametersFromUI(rMrProt, rSeqLim, rSeqExpo))
    {
        SEQ_TRACE_ERROR.print("ERROR: SBBDiffusion_Base::setParameterFromUI returned false");
        return false;
    }

    // -------------------------
    // Basic preparations
    // -------------------------
    if(!prepParameters(rMrProt, rSeqLim, rSeqExpo))
    {
        if(!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("ERROR: SBBDiffusion_Base::prepParameters returned false");
        }
        return false;
    }

    // ---------------------------------------------------------
    // Prepare diffusion ordering (only used if thermal balancing is active, otherwise only for GSWD checking)
    // ---------------------------------------------------------
    if (!prepDiffusionOrder(rMrProt, rSeqLim, rSeqExpo))
    {
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ALWAYS.print("ERROR: Diffusion ordering preparation returned false");
        }
        return false;
    }
    

    //--------------------------------------------------
    // Get maximum b value (which determines the timing)
    //--------------------------------------------------
    double dMaxRequestedBValue = (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE) ? m_dQSpaceMaxBValue : *std::max_element(m_vdBValues.begin(), m_vdBValues.end());

    // Ensure a value > 0 (avoid division by zero)
    dMaxRequestedBValue = std::max<double>(1., dMaxRequestedBValue);

    // -------------------------------
    // Additional initial preparations
    // -------------------------------

    // Call ::prepInit of derived class. 

    // Anything required by ::prepTiming has to be prepared here. Usually, this 
    // includes spoiler gradients (PE-axis only) and the calculation of the reference 
    // spoil moment m_dRefSpoilMoment.
    if(!prepInit(rMrProt, rSeqLim, rSeqExpo))
    {
        if(!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("ERROR: ::prepInit returned false");
        }
        return false;
    }

    // ----------------
    // Calculate Timing
    // ----------------

    // Depending on m_bMinimizeTE, two strategies apply:
    //
    // false: Check whether the required maximum b-value can be realized with
    //        the current timing (i.e. TE). 
    //        If yes: prepare everything and return true. 
    //        If no:  increase TE up to the necessary value, prepare everything and return true
    // 
    // true:  Identify the minimum TE that is required in order realize
    //        the required protocol (i.e. maximum b-value). Prepare everything
    //        and return true.
    //
    // The sequence kernel receives the information on changes in TE
    // via the parameters m_lPreEchoTimeContrib and m_lPostEchoTimeContrib.

    // Initialize TE search limits and smallest search increment
    long       lMinTE;
    const long lMaxTE = m_vlTEMax_Limit[m_iTEArrayIndex];
    
    // The TE increment in the main cpp is 100. The UILink increases this increment by a factor of 10 for all TEs > 10ms (which is always the case for diffusion).
    // Hence m_lTEInc_Limit = 100, but the UI limit is 1000. By simulation of many protocols it turns out that half the increment (=500) is a good stop criterion
    // for the following binary search.
    const long lDeltaTE  = m_lTEInc_Limit * 5;  

    bool bStartBinarySearch = false;

    // with RESOLVE sequence perform a binary search in all contexts to avoid possible UI inconsistencies

    if ((!m_bMinimizeTE || !(rSeqLim.isContextPrepForBinarySearch() || rSeqLim.isContextPrepForScanTimeCalculation())) && !m_bResolve)
    {
        // If we are not supposed to minimize TE OR if we assume that
        // the minimum TE is already set, the current protocol value
        // should be used as the minimum.
        lMinTE = m_vlTE[m_iTEArrayIndex];

        // Try current TE (i.e. calculate maximum b-value that can be realized).
        // Note: m_dMaxPossibleBValue is calculated within prepTiming()
        // Note: pMrProt, pSeqLim and pSeqExpo are only required for some libRT interfaces
        if(!prepTiming(rMrProt, rSeqLim, rSeqExpo, lMinTE))
        {
            // prepTiming will ONLY return false if a severe problem occured (e.g. no solution
            // for the eddy current nulled gradient scheme or event preparations that failed).
            // Otherwise, it will prepare everything for the maximum possible diffusion
            // gradient amplitudes and calculate the related b-value.

            // prepTiming failed and we are not running in ContextPrepForBinarySearch
            // => something went wrong
            if(!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("ERROR: prepTiming returned false ");
                return false;
            }
            
            // prepTiming failed, but we are running in ContextPrepForBinarySearch 
            // => start search for minimum TE that allows for the given b-value
            bStartBinarySearch = true;
        }
        else
        {
            // prepTiming succeeded => maximum b-value for the current TE has been calculated

            if(m_dMaxPossibleBValue < dMaxRequestedBValue)
            {
                // Current TE does not allow for requested b-value
                // => start search for minimum TE that allows for the given b-value
                bStartBinarySearch = true;
            }
        }
    }
    else
    {
        // We are supposed to find the minimum TE that allows for the
        // given b-value => set lower search limit to zero.

        // Don't try a calcTiming() with this TE - it will fail. 
        // We have to start the optimization process anyway.
        lMinTE = 0;

        // Indicate that current timing is not possible.
        m_dMaxPossibleBValue = 0.;

        // => start search for minimum TE that allows for the given b-value
        bStartBinarySearch = true;
    }

    m_lActualTE = lMinTE;

    // Start the binary search if
    // - the minimum TE has to be identified or
    // - if the current TE does not allow for the required b-value
    if(bStartBinarySearch)
    {
        /*  Now we are in trouble: With the current timing, we cannot reach the
        required b value.
        The expected behavior of the SBB is now to calculate the TE which
        is necessary to prepare this b value successfully. This new TE will
        be returned to the calling sequence indirectly by setting the variable
        DurationPerRequest.
        As there is no easy algorithm to calculate the necessary TE from the
        requested b value, we do a binary search over TE.
        As TE is no "normal" variable, but always accessed as a part of the
        MrProt class, we use a copy of MrProt, called "KastriertesProtocol".
        It is now possible to write our TE value to this copy.
        */

        // Start value for binary search (mind that TE must be on the double gradient raster time)
        long lCurrentTE = static_cast<long>(
            fSDSDoubleRoundUp(0.0, m_MaxValueForRounding, static_cast<double>(lMaxTE), 2.0 * GRAD_RASTER_TIME));

        // Start increment x 2 for binary search
        long lInc       = lCurrentTE;

        // Calculate the b value for the TE starting value
        if(!prepTiming(rMrProt, rSeqLim, rSeqExpo, lCurrentTE))
        {
            // Even with lMaxTE, no valid timing can be calculated.
            if(!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("ERROR: Timing recalculation failed unexpectedly.");
            }
            setNLSStatus(MRI_SBB_SBB_NEGATIV_TEFILL);
            return false;
        }

        // variable to remember the last valid TE
        long lValidLastTE = lCurrentTE;

        // Start binary search loop now by adjusting the TE:
        while(lInc > lDeltaTE)
        {
            lInc /= 2;

            if(m_dMaxPossibleBValue >= dMaxRequestedBValue)
            {
                lCurrentTE -= lInc;
            }
            else
            {
                lCurrentTE += lInc;
            }

            if(lCurrentTE < lMinTE)
            {
                lCurrentTE = lMinTE;
            }

            // Mind that TE must be on the double gradient raster time!
            lCurrentTE = static_cast<int>(
                fSDSDoubleRoundUp(0.0, m_MaxValueForRounding, static_cast<double>(lCurrentTE), 2.0 * GRAD_RASTER_TIME));

            if(!prepTiming(rMrProt, rSeqLim, rSeqExpo, lCurrentTE))
            {
                // Indicate that current timing is not possible.
                m_dMaxPossibleBValue = 0.;
            }

            // remember last TE that was successful
            if(m_dMaxPossibleBValue >= dMaxRequestedBValue)
            {
                lValidLastTE = lCurrentTE;
            }
        }

        // one last check: if we are not on an integer ms, it might happen that the integer ms just below the current
        // value did not get tested, although it would work.
        // test it manually.
        const long lLastValidTERoundDown = static_cast<long>(
            fSDSDoubleRoundDown(0.0, m_MaxValueForRounding, static_cast<double>(lValidLastTE), 1000.0));

        if (lLastValidTERoundDown < lValidLastTE)
        {
            if (prepTiming(rMrProt, rSeqLim, rSeqExpo, lLastValidTERoundDown))
            {
                if (m_dMaxPossibleBValue >= dMaxRequestedBValue)
                {
                    lValidLastTE = lLastValidTERoundDown;
                }
            }
        }

        //set last TE which was possible
        lCurrentTE = lValidLastTE;


        // Although we are already on the gradient raster time, we might still be off the
        // increment raster of SeqLim. In principle, fSeqPrep would take care of this - however, 
        // increasing TE might in seldom cases increase (!) the GPA load and thus yield an
        // inconsistent protocol.
        // Therefore we increase TE stepwise and in accordance with MrUILink here until we have
        // both, a valid TE value and the required b-value.
        lCurrentTE = lCalcTEOnInc(lCurrentTE);

        if(!prepTiming(rMrProt, rSeqLim, rSeqExpo, lCurrentTE))
        {
            if(!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("ERROR: Timing recalculation failed unexpectedly.");
            }
            setNLSStatus(MRI_SBB_SBB_NEGATIV_TEFILL);
            return false;
        }

        while(m_dMaxPossibleBValue < dMaxRequestedBValue)
        {
            lCurrentTE += lDeltaTE;
            lCurrentTE = lCalcTEOnInc(lCurrentTE);

            if(!prepTiming(rMrProt, rSeqLim, rSeqExpo, lCurrentTE))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: Timing recalculation failed unexpectedly.");
                }
                setNLSStatus(MRI_SBB_SBB_NEGATIV_TEFILL);
                return false;
            }
        }

        m_lActualTE = lCurrentTE;
    }


    // -------------------------------
    // Ice exports
    // -------------------------------

    // IceProgramParam(ICE_PROGRAM_DIFF_THRESHOLD): Threshold for noise
    rSeqExpo.setICEProgramParam(ICE_PROGRAM_DIFF_THRESHOLD, m_lNoiseThreshold);

    // IceProgramParam(ICE_PROGRAM_DIFF_NO_OF_DIRECTIONS): No. of directions for the first b value
    // Note: here we assume that b values are in ascending order (as provided by UI)
    if (!(m_eDiffusionMode == SEQ::DIFFMODE_QSPACE) && (m_lDirections > 1) && (m_vdBValues[0] <= ONE_FOR_THREE_THRESHOLD))
    {
        rSeqExpo.setICEProgramParam(ICE_PROGRAM_DIFF_NO_OF_DIRECTIONS, 1);
    }
    else
    {
        rSeqExpo.setICEProgramParam(ICE_PROGRAM_DIFF_NO_OF_DIRECTIONS, m_lDirections);
    }

    // -------------------------------
    // Additional final preparations
    // -------------------------------

    // Call ::prepFinal  of derived class.
    if(!prepFinal(dMaxRequestedBValue, rSeqLim.isContextPrepForBinarySearch()))
    {
        if(!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("ERROR: ::prepInit returned false");
        }
        return false;
    }

    // Validate exports
    setPrepared();

    setNLSStatus(MRI_SBB_SBB_NORMAL);

    return true;
}


// ===========================================================================
/*!
\author PLM AW Neuro

\brief  Copy relevant UI (protocol) parameters to member variables
*/
// ===========================================================================
bool SBBDiffusion_Base::setParametersFromUI(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo & /* rSeqExpo */)
{
    // Reset preparation status
    resetPrepared();

    const MeasNucleus myMeasNucleus = rMrProt.txSpec().nucleusInfoArray()[0].gettNucleus().c_str();

    m_dGamma                   = myMeasNucleus.getLarmorConst() * 2. * M_PI * 1.e6;     // [1 / T s]
    m_lRepetitions             = rMrProt.repetitions();
    m_eRFPulseType             = rMrProt.txSpec().rfPulseType();
    m_eDiffusionMode           = rMrProt.diffusion().getulMode();
    m_eDynDistMode             = rMrProt.getsDynDistCorrFilter().getucMode();
    m_lDiffusionDirectionsMDDW = rMrProt.diffusion().getlDiffDirections();
    m_lBValueInc_Limit         = rSeqLim.getBValue().getInc();
    m_lTEInc_Limit             = rSeqLim.getTE()[0].getInc();
    m_lSlices                  = rMrProt.sliceSeries().getlSize();
    m_dSliceThickness          = rMrProt.sliceSeries().aFront().thickness();
    m_dReadoutMoment           = 1.0E6 * rMrProt.kSpace().getlBaseResolution()  / (myMeasNucleus.getLarmorConst() * rMrProt.sliceSeries().front().readoutFOV());

    m_strFreeUserComment.clear();
    m_vFreeDiffDir.resize(0);
    m_eFreeCoordinateSystem    = MrProtocolData::DIFFDIR_CS_XYZ;    // Dummy initialization

    if(m_eDiffusionMode == SEQ::DIFFMODE_FREE)
    {
        if(m_lDiffusionDirectionsMDDW != rMrProt.diffusion().getsFreeDiffusionData().getlDiffDirections())
        {
            // Should never happen
            SEQ_TRACE_ERROR.print("ERROR: inconsistent number of diffusion directions");
            return false;
        }
        m_eFreeCoordinateSystem = rMrProt.diffusion().getsFreeDiffusionData().getulCoordinateSystem();
        // Extract and store user comment (strip first line which contains the filename)
        const std::string strTemp     = rMrProt.diffusion().getsFreeDiffusionData().getsComment();
        m_strFreeUserComment    = strTemp.substr(strTemp.find_first_of('\n') + 1);
        m_vFreeDiffDir.resize(m_lDiffusionDirectionsMDDW);
        for(auto lI = 0; lI < m_lDiffusionDirectionsMDDW; ++lI)
        {
            m_vFreeDiffDir[lI].dx = rMrProt.diffusion().getsFreeDiffusionData().getasDiffDirVector()[lI].getdSag();
            m_vFreeDiffDir[lI].dy = rMrProt.diffusion().getsFreeDiffusionData().getasDiffDirVector()[lI].getdCor();
            m_vFreeDiffDir[lI].dz = rMrProt.diffusion().getsFreeDiffusionData().getasDiffDirVector()[lI].getdTra();
        }
    }

    m_vdBValues.resize(0);
    m_vlLocalAverages.resize(0);

    if(m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
    {
        // Store q-space related protocol parameters
        // (b-value related parameters are meaningless: do not use them in this case)


        m_eQSpaceCoverage  = rMrProt.diffusion().getulQSpaceCoverage();
        m_eQSpaceSampling  = rMrProt.diffusion().getulQSpaceSampling();
        m_dQSpaceMaxBValue = rMrProt.diffusion().getlQSpaceMaxBValue();
        m_lQSpaceSteps     = rMrProt.diffusion().getlQSpaceSteps();

    }
    else
    {
        // Store b-value related protocol parameters
        // (q-space related parameters are meaningless: do not use them in this case)
        m_vdBValues.reserve(rMrProt.diffusion().getlDiffWeightings());
        m_vlLocalAverages.reserve(rMrProt.diffusion().getlDiffWeightings());

        // Get the max value in average array:
        m_lMaxValueInAveArray = *std::max_element(rMrProt.getsDiffusion().getalAverages().begin(), rMrProt.getsDiffusion().getalAverages().end());

        for(auto lI = 0; lI < rMrProt.diffusion().getlDiffWeightings(); ++lI)
        {
            m_vdBValues.push_back(static_cast<double>(rMrProt.diffusion().getalBValue()[lI]));
            m_vlLocalAverages.push_back(std::max(static_cast<long>(1), static_cast<long>  (rMrProt.getsDiffusion().getalAverages()[lI])));
        }
    }

    m_vlTEMax_Limit.clear();
    m_vlTE.clear();
    for(auto lI = 0; lI < rSeqLim.getTE().getNoElements(); ++lI)
    {
        m_vlTEMax_Limit.push_back(rSeqLim.getTE()[m_iTEArrayIndex].getMax());
        m_vlTE.push_back(rMrProt.te()[m_iTEArrayIndex]);
    }

    // 20140320 DP: MR_00443009: we need to know if this protocol uses ZOOMit
    if(rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
    {
        m_bZoomedExcitationIsUsed = true;
    }
    else
    {
        m_bZoomedExcitationIsUsed = false;
    }

    return true;
}

// ===========================================================================
///	This method checks input and WIP parameters.

/**   This method should be called at the beginning of each prep function
of a derived class. It does the following general jobs:
- Set the Gradient ramp time
- Set m_NoOfWeightings and m_lDirections (using the Didi class)
- Asserts that the input parameters m_lSpinPrepTimeus, m_lADCusTillEcho
and rMrProt.te()[m_iTEArrayIndex] are on the gradient raster
- Display Diffusion Parameter Card
- Initializes the WIP parameters (if there are any).
*/
// ===========================================================================
bool SBBDiffusion_Base::prepParameters(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
// ===========================================================================
{
    SEQ::DiffusionMode eMode;                 // Diffusion Mode


    // ----------------------
    // Set gradient ramp time
    // ----------------------
    // This cannot yet be done in the constructor as the sequence might
    // modify m_dMinRiseTime and m_dMaxAmpl
    if(m_dMinRiseTime <= 0.001)
    {
        // Programming error: trace always
        SEQ_TRACE_ERROR.print("ERROR: No or negative rise time specified.");
        setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }

    m_lRampTime = fSDSRoundUpGRT(m_dMinRiseTime * m_dMaxAmpl);

    // Slice selection gradient amplitude is <= GradMaxAmplFast, thus
    // the ramps can be shorter.
    m_lRampTimeRF = fSDSRoundUpGRT(m_dMinRiseTime * SysProperties::getGradMaxAmpl(SEQ::GRAD_FAST));

#ifdef DEBUG
    bool bWorkAsDebug = true;
#else
    bool bWorkAsDebug = SeqUT.isUnitTestActive();
#endif
    if(bWorkAsDebug)
    {
        /* A individual rotation matrix must be valid for at least 300 us.
        So be sure that a triangular gradient is not shorter.
        For performance reasons, this check is only enabled for debug versions. */
        if(2 * m_lRampTime < 300)
        {
            if(!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("ERROR: m_lRampTime=%ld is too short", m_lRampTime);
            }
            setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
            return false;
        }
    }

    // --------------------------------
    // Prepare directions
    // --------------------------------

    eMode            = m_eDiffusionMode;
    m_NoOfWeightings = (eMode == SEQ::DIFFMODE_QSPACE) ? 1 : static_cast<int>(m_vdBValues.size());

    // PRS = Phase-Read-Slice coordinates = current slice rotation matrix
    // XYZ = x-y-z Magnet hardware coordinates = unity rotation matrix
    switch(eMode)
    {
        case SEQ::DIFFMODE_PHASE:
            m_lDirections = 1;
            if(!m_Didi.prepInternal(m_lDirections, 'P', rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepInternal returned false.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_READ:
            m_lDirections = 1;
            if(!m_Didi.prepInternal(m_lDirections, 'R', rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepInternal returned false.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_SLICE:
            m_lDirections = 1;
            if(!m_Didi.prepInternal(m_lDirections, 'S', rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepInternal returned false.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_DIAGONAL:
            m_lDirections = 1;
            if(!m_Didi.prepInternal(m_lDirections, 'D', rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepInternal returned false.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_ORTHOGONAL:
            m_lDirections = 3;
            if(!m_Didi.prepInternal(m_lDirections, 'O', rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepInternal returned false.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_THREE_SCAN_TRACE:   // runs with unity matrix
            m_lDirections = 3;
            if(!m_Didi.prepInternal(m_lDirections, 'T', rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepInternal returned false.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_FOUR_SCAN_TRACE:   // runs with unity matrix
            m_lDirections = 4;
            if(!m_Didi.prepInternal(m_lDirections, 'T', rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepInternal returned false.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_TENSOR:   // runs with unity matrix
            m_lDirections = m_lDiffusionDirectionsMDDW;
            if(!m_Didi.prepInternal(m_lDirections, '\0', rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepInternal returned false.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_QSPACE:   // runs with unity matrix
        {
            // Relevant for partial q-space coverage only: set vector that defines the
            // hemisphere that should get scanner preferably. Default: +z hemisphere.
            VectorStruct sHemisphere;
            sHemisphere.dx = m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/QSpaceHemisphereX", 0.0);
            sHemisphere.dy = m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/QSpaceHemisphereY", 0.0);
            sHemisphere.dz = m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/QSpaceHemisphereZ", 1.0);
            m_Didi.setQSpaceHemisphere(sHemisphere);

            // Prepare q-space diffusion vectors (number of directions depends on actual settings)
            m_lDirections = m_Didi.prepQSpace(m_eQSpaceCoverage, m_eQSpaceSampling, m_lQSpaceSteps, rSeqLim.isContextPrepForBinarySearch());
            if(m_lDirections == 0)
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepQSpace returned false");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
        }
            break;
        case SEQ::DIFFMODE_FREE:   // runs with unity matrix
            m_lDirections = m_lDiffusionDirectionsMDDW;
            if(!m_Didi.prepExplicit(m_lDirections, m_eFreeCoordinateSystem, m_strFreeUserComment, m_vFreeDiffDir, rSeqLim.isContextPrepForBinarySearch()))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERROR.print("ERROR: m_Didi.prepProtocol with %ld directions returned false.", m_lDirections);
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
            break;
        case SEQ::DIFFMODE_ONE_SCAN_TRACE:
            // the unloved child: just set a reasonable number of directions and survive
            m_lDirections = 1;
            break;
        default:
            SEQ_TRACE_ERROR.print("ERROR: Invalid diffusion mode %d", eMode);
            return false;

    }

    // ---------------------------------------------------------------------------
    // Assure that input timing parameters are on gradient raster
    // ---------------------------------------------------------------------------

    // It is required that m_lSpinPrepTimeus is on the gradient raster
    // This condition cannot be fixed because the total duration of an
    // event block must be on the gradient raster.
    if(m_lSpinPrepTimeus % GRAD_RASTER_TIME)
    {
        if(!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("ERROR: m_lSpinPrepTimeus=%ld not on gradient raster", m_lSpinPrepTimeus);
        }
        setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }

    // It is required that m_lADCusTillEcho is on the gradient raster
    if(m_lADCusTillEcho % GRAD_RASTER_TIME)
    {
        if(!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("ERROR: m_lADCusTillEcho=%ld not on gradient raster", m_lADCusTillEcho);
        }
        setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }

    // TE must have a 2*GRAD_RASTER_TIME raster
    // Check this (as the sequence programmer might try to ignore this):
    if(m_vlTE[m_iTEArrayIndex] % (2*GRAD_RASTER_TIME))
    {
        // this is no problem as TE has usually(!) an increment of > 100
        if(!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("ERROR: TE[%d]=%ld not on double gradient raster", m_iTEArrayIndex, m_vlTE[m_iTEArrayIndex]);
        }
        setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }


    // ---------------------------------------------------------------------------
    // Display Diffusion Parameter Card
    // ---------------------------------------------------------------------------

    rSeqExpo.setApplicationCard(SEQ::APPLICATION_CARD_DIFF);
    rSeqExpo.setApplicationCard(SEQ::APPLICATION_CARD_FMRI);
    rSeqExpo.setApplicationCardName(SEQ::APPLICATION_CARD_NAME_DIFF);


#ifdef WIP
    //. ---------------------------------------------------------------------------
    //. Initialize parameters on the Sequence/Special card
    //. ---------------------------------------------------------------------------
#ifdef WIN32
    // In this part (only performed on the host) we initialize the values
    //  for the parameters on the Sequence/Special card.
    // The concept of ContextPrepForBinarySearch and ContextPrepForMrProtUpdate
    //  is used here:
    //  When creating the default protocol, fSEQPrep() is called in
    //   ContextPrepForBinarySearch. The parameters are not yet initialized,
    //   hence a MRI_SEQ_SEQU_ERROR is returned. This causes a 2nd call of fSEQPrep(),
    //   this time in ContextPrepForMrProtUpdate,
    //   allowing for a modification of the protocol.
    // Now the parameters can be initialized.

    // if called within ContextPrepForMrProtUpdate, we can (and must) modify the protocol.
    if(rSeqLim.isContextPrepForMrProtUpdate())
    {
        PRINT0(DEBUG_INTERNAL, "WIP parameters initialized.");

        // set the WIP parameters to their default value:
        if(rMrProt.wipMemBlock().getadFree()[WIP_Anything] == 0.)
        {
            rMrProt.wipMemBlock().getadFree()[WIP_Anything] = 1.;
        }
    }


    // If any of the WIP parameters have not yet been initialized,
    //  return with error to induce a ContextPrepForMrProtUpdate
    if(rMrProt.wipMemBlock().getadFree()[WIP_Anything] == 0)
    {
#ifdef DEBUG
        if(!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("WIP parameter must be initialized.");
        }
#endif
        setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }

#endif   // of #ifdef WIN32

#endif   // of WIP


    setNLSStatus(MRI_SBB_SBB_NORMAL);
    return true;
}


// ===========================================
// This method performs the diffusion ordering
// if thermal balancing is active
// ===========================================
bool SBBDiffusion_Base::prepDiffusionOrder(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo)
{
    m_sDiffusionOrderInfo.setDiffusionMode(m_eDiffusionMode);
    m_sDiffusionOrderInfo.setDirections(m_lDirections);
    m_sDiffusionOrderInfo.setBValues(m_vdBValues);
    m_sDiffusionOrderInfo.setLocalAverages(m_vlLocalAverages);
    m_sDiffusionOrderInfo.setThermalBalancingReorder(m_bThermalBalancing);

    if (!m_sDiffusionOrderInfo.prepAll(rSeqLim.isContextPrepForBinarySearch(), &m_Didi))
    {
        SEQ_TRACE_ALWAYS.print("ERROR: m_sDiffusionOrderInfo.prepAll() failed");
        return false;
    }

    return m_sDiffusionOrderInfo.isFullyPrepared();
}


// ============================================================
// Diffusion loop counter for highest b-value, for kernel check
// ============================================================
long SBBDiffusion_Base::getDiffLoopCounterForHighestBValue()
{
    return m_sDiffusionOrderInfo.getDiffLoopCounterForHighestBValue();
}


// ============================================================
// Diffusion loop counter for scan with largest Read-axis component
// ============================================================
long SBBDiffusion_Base::getDiffLoopCounterForHighestReadComponent(sSLICE_POS* pSlice)
{
    return m_sDiffusionOrderInfo.getDiffLoopCounterForHighestReadComponent(pSlice, getDidiPointer());
}


// ============================================================
// Set thermal balancing flag
// ============================================================
void SBBDiffusion_Base::setThermalBalancing(bool bThermalBalancing)
{
    m_bThermalBalancing = bThermalBalancing;
}


// ============================================================
// Get thermal balancing flag
// ============================================================
bool SBBDiffusion_Base::getThermalBalancing() const
{
    return m_bThermalBalancing;
}

// ============================================================
// set parameters needed before the preparation of compensation gradients
// e.g. min rise time, max amplitude.
// ============================================================
void SBBDiffusion_Base::prePrepareCompGrad(MrProt& rMrProt, SeqLim& rSeqLim, long lRampTimeOutsideSBB)
{
    if (getUseGPABalance())
    {
        m_CompGrad.setbUseGPABalance(true);
    }
    else
    {
        m_CompGrad.setbUseGPABalance(false);
    }

    // MK: in cases where max. magnitude and min. rise time are close to the system limits (e.g. boost mode at
    // Free.Max): reduce magnitudes / increase rise times to avoid exceeding the system limit with oblique FoVs
    // solver" would have to intervene to often otherwise

    double         dGradSpecFactor = 1.2;
    SEQ::Gradients eGradMode       = rMrProt.getsGRADSPEC().getucMode();
    double         dMinRiseTime    = std::max(
        dGradSpecFactor * SysProperties::getGradMinRiseTimeAbsolute(), m_CompGrad.getMinRiseTime(eGradMode, 0));

    double dMaxMagnitude
        = std::min(SysProperties::getGradMaxAmplAbsolute() / dGradSpecFactor, m_CompGrad.getMaxMagnitude(eGradMode, 0));
    m_CompGrad.setMaxMagnitude(eGradMode, dMaxMagnitude, -1);
    m_CompGrad.setMinRiseTime(eGradMode, dMinRiseTime, -1);
    m_CompGrad.setCompGradParaLimit(50000, 1000000.0); // Max total time: 50ms, Max moment: 1000000.0

    m_CompGrad.setdMaxAmplitude(dMaxMagnitude); // always use Max Amplitude to calculated the duration
    m_CompGrad.setdMinRiseTime(dMinRiseTime);
    m_CompGrad.setCompensationPara(m_bCompensationDecay, m_dCompensationFraction, m_dEddycurrentTau);

    m_CompGrad.setsSliceAdjParametersRequestedBySequence(SLICEADJ::ADJNONE);
    m_CompGrad.setDesiredRampUpTime(lRampTimeOutsideSBB);
    m_CompGrad.setGSWDGradientPerformance(rMrProt, rSeqLim);
    m_CompGrad.setSATSpoiler(false);
    m_CompGrad.resetMoments();
    m_CompGrad.setMomentsFor10mm(/*RO:*/ 0.0, /*PE:*/ 0.0, /*SS:*/ 0.0);

    // assign these values for SBB
    m_CompGrad.setCalcMode(SeqBuildBlockSpoilGrad::eSpoilMode::eFixedDurationAndMoment);
}

// ===========================================================================
///   This method prepares the spoiler gradients (DSp1, DSr1, DSs1).
/**   It must be called explicitely in the prep() function .
From  the gradient properties, the shortest
gradient to produce the specifed moment will be calculated and prepared.
\pre    m_lRampTime		Duration of the gradient ramps
\pre    m_dMaxAmpl		Maximum gradient amplitude -
\return m_lSpoilerTotalTime

\remark For small moments, a triangular gradient will be generated with
the specified ramp times 2*m_lRampTime. In general, the amplitude will
be below the m_dMaxAmp.

\remark The moment actually applied is usually greater than the moment
specified here due to rounding the gradient duration to the gradient raster time.
In case of spoilers, this behaviour can be accepted.


*/
// ===========================================================================
bool SBBDiffusion_Base::prepSpoilGrad(
    double dMoment   /*!<  The gradient moment to be achieved, specified in mT/m*ms.
                     A negative sign is allowed and will be interpreted as a
                     negative gradient amplitude. */
                     )
                     // ===========================================================================
{
    // EPIRO with regridding causes an FID artefact.
    // It can be avoided by applying a moment which has a moment of the readout ramp, at least.


    // In the past, spoiler amplitude and ramp time were identical to those of the diffusion
    // gradients. However, with the introduction of variable diffusion gradient amplitudes
    // this yielded inconsistencies. Thus, the amplitude now is set to the maximum gradient
    // amplitude of the actual GPA, divided by sqrt(3) in order to take into account oblique
    // slice orientations. The ramp time is calculated using the defined minimum
    // rise time and the acutal gradient amplitude.

    double dSpoilerAmplitude = SysProperties::getGradMaxAmplAbsolute() / sqrt(3.);
    long   lFlattopDuration  = 0;
    long   lSpoilerRampTime  = fSDSRoundUpGRT(m_dMinRiseTime * dSpoilerAmplitude);

    const int    iSign      = (dMoment > 0 ? 1 : -1);
    const double dAbsMoment = (dMoment > 0 ? dMoment : -dMoment) * 1000; // ms -> us

    const auto   dSpoilerRampTime = static_cast<double>(lSpoilerRampTime);
    const double dRampMoment      = 2. * 0.5 *  dSpoilerAmplitude * dSpoilerRampTime;

    if(fGSLAlmEqual(dMoment, 0.0))
    {
        dSpoilerAmplitude = 0.0;
        lSpoilerRampTime  = 0;
    }
    else if(dAbsMoment > dRampMoment)
    {
        // we need some flat top to achieve the requested moment
        lFlattopDuration  = static_cast<long>(fSDSDoubleRoundUp(
            0.0, m_MaxValueForRounding, (dAbsMoment - dRampMoment) / dSpoilerAmplitude,
            GRAD_RASTER_TIME));
    }
    else
    {
        // we can scale down the amplitude of the triangle gradient pulse
        dSpoilerAmplitude = dAbsMoment / dSpoilerRampTime;
    }

    m_lSpoilerTotalTime = 2*lSpoilerRampTime + lFlattopDuration;
    dSpoilerAmplitude *= iSign;

    m_DSp1.setRampTimes(lSpoilerRampTime);
    m_DSp1.setDuration(lSpoilerRampTime + lFlattopDuration);
    m_DSp1.setStartTime(0);
    m_DSp1.setAmplitude(dSpoilerAmplitude);

    if(! m_DSp1.prep())
    {
        SEQ_TRACE_ERROR.print("ERROR: m_DSp1.prep() failed.");
        setNLSStatus(m_DSp1.getNLSStatus());
        return false;
    }

    if(! m_DSp1.check())
    {
        SEQ_TRACE_ERROR.print("ERROR: m_DSp1.check() failed.");
        printGradient(&m_DSp1);
        setNLSStatus(m_DSp1.getNLSStatus());
        return false;
    }

    m_DSr1.setRampTimes(lSpoilerRampTime);
    m_DSr1.setDuration(lSpoilerRampTime + lFlattopDuration);
    m_DSr1.setStartTime(0);
    m_DSr1.setAmplitude(dSpoilerAmplitude);

    if(! m_DSr1.prep())
    {
        SEQ_TRACE_ERROR.print("ERROR: m_DSr1.prep() failed.");
        setNLSStatus(m_DSr1.getNLSStatus());
        return false;
    }

    if(! m_DSr1.check())
    {
        SEQ_TRACE_ERROR.print("ERROR: m_DSr1.check() failed.");
        setNLSStatus(m_DSr1.getNLSStatus());
        return false;
    }

    m_DSs1.setRampTimes(lSpoilerRampTime);
    m_DSs1.setDuration(lSpoilerRampTime + lFlattopDuration);
    m_DSs1.setStartTime(0);
    m_DSs1.setAmplitude(dSpoilerAmplitude);

    if(! m_DSs1.prep())
    {
        SEQ_TRACE_ERROR.print("ERROR: m_DSs1.prep() failed.");
        setNLSStatus(m_DSs1.getNLSStatus());
        return false;
    }

    if(! m_DSs1.check())
    {
        SEQ_TRACE_ERROR.print("ERROR: m_DSs1.check() failed.");
        setNLSStatus(m_DSs1.getNLSStatus());
        return false;
    }

    return true;
}




// ===========================================================================
///   This method prepares the specified gradient to produce the specified moment.
/**   From  the gradient properties, the shortest
gradient to produce the specifed moment will be calculated and prepared.
As the spoiler gradients are applied in the phase-read-slice coordinate
system, the maximum amplitude is m_dMaxAmpl/sqrt(3).
\pre    m_lRampTime		Duration of the gradient ramps
\pre    m_dMaxAmpl		Maximum gradient amplitude -
\return Grad.RampUpTime, Grad.RampDownTime, Grad.Duration, Grad.Amplitude
*/
// ===========================================================================
bool SBBDiffusion_Base::prepGradMoment(
    sGRAD_PULSE_TRAP *sGrad, /*!<  A pointer to the gradient object to be prepared.*/
    double dMoment           /*!<  The gradient moment to be achieved, specified in mT/m*ms.
                             A negative sign is allowed and will be interpreted as a
                             negative gradient amplitude. */
                             )
                             // ===========================================================================
{
    double       dSpoilerAmplitude = m_dMaxAmpl / sqrt(3.0);
    long         lFlattopDuration  = 0;
    long         lSpoilerRampTime  = m_lRampTime;
    const double dSlewRate         = m_dMaxAmpl / static_cast<double>(m_lRampTime);

    const int    iSign      = (dMoment < 0 ? -1 : 1);
    const double dAbsMoment = (dMoment > 0 ? dMoment : -dMoment) * 1000; // ms -> us

    const double dRampMoment       = 2. * 0.5 *  dSpoilerAmplitude * static_cast<double>(lSpoilerRampTime);

    // Assert that we do not have a null pointer
    assert(sGrad);

    if(fGSLAlmEqual(dMoment, 0.0))
    {
        dSpoilerAmplitude = 0.0;
        lSpoilerRampTime  = 0;
    }
    else if(dAbsMoment > dRampMoment)
    {
        // we need some flat top to achieve the requested moment
        lFlattopDuration  = static_cast<long>(fSDSDoubleRoundUp(
            0.0, m_MaxValueForRounding, (dAbsMoment - dRampMoment) / dSpoilerAmplitude,
            GRAD_RASTER_TIME));
    }
    else
    {
        // we need a triangle gradient pulse and can scale the ramp times: M= RampTime * Amplitude
        // However we must stick to the slew rate specification: Slewrate=Ampl / RampTime
        // --> Moment = RampTime^2 * SlewRate
        lSpoilerRampTime = static_cast<long>(fSDSDoubleRoundUp(0.0, m_MaxValueForRounding, sqrt(dAbsMoment / dSlewRate),
            GRAD_RASTER_TIME));
        // Amplitude will be calculated in next section ...
    }

    // Due to the rounding of the ramp times and duration to the gradient raster time,
    // the generated moment will be a bit too high.
    // Therefore we have to adjust the gradient ampltitude .
    if(lSpoilerRampTime)   // avoid division by zero!
    {
        const auto dFlattopDuration  = static_cast<double>(lFlattopDuration);
        const auto dSpoilerRampTime  = static_cast<double>(lSpoilerRampTime);

        dSpoilerAmplitude = dAbsMoment / (2. * 0.5 *  dSpoilerRampTime + dFlattopDuration);
    }


    // Set amplitude polarity
    dSpoilerAmplitude *= iSign;


#ifdef DEBUG
    if(! fGSLAlmEqual(1000.0 * dMoment, dSpoilerAmplitude * (lSpoilerRampTime+lFlattopDuration)))
    {
        SEQ_TRACE_WARN.print("WARNING: Moment mismatch: requested=%f, calculated: %f",
                   1000.0 * dMoment, dSpoilerAmplitude * (lSpoilerRampTime+lFlattopDuration));
    }
#endif


    sGrad->setRampTimes(lSpoilerRampTime);
    sGrad->setDuration(lSpoilerRampTime + lFlattopDuration);
    sGrad->setAmplitude(dSpoilerAmplitude);

    if(! sGrad->prep())
    {
        SEQ_TRACE_ERROR.print("ERROR: Could not prepare gradient.");
        printGradient(*sGrad);
        setNLSStatus(sGrad->getNLSStatus());
        return false;
    }

    if(! sGrad->check())
    {
        SEQ_TRACE_ERROR.print("ERROR: Check of gradient with %f mT/m failed:", sGrad->getAmplitude());
        printGradient(*sGrad);
        setNLSStatus(sGrad->getNLSStatus());
        return false;
    }

    return true;
}



// ===========================================================================
///	Get total number of diffusion weighted scans (excluding prescans and the   
/// global average from MrProt->Averages(), but including local averages 
/// for individual BValue).
///
// ===========================================================================
long SBBDiffusion_Base::getTotalScans(bool bIncludeLocalAverages)
// ===========================================================================
{
    long lReallyTotalScans = 0;

    if(m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
    {
        // In q-space mode, only the direction index is relevant
        // => ignore local averages and b-values
        lReallyTotalScans = m_lDirections;
    }
    else
    {
        for(unsigned int i = 0; i < m_vdBValues.size(); i++)
        {
            const long lLocalAverages = bIncludeLocalAverages ?  m_vlLocalAverages[i] : 1;

            if(m_vdBValues[i] <= ONE_FOR_THREE_THRESHOLD)
            {
                lReallyTotalScans += lLocalAverages;
            }
            else
            {
                lReallyTotalScans += m_lDirections * lLocalAverages;
            }
        }
    }

    return lReallyTotalScans;
}


// ===========================================================================
///	Returns  the value of the average index, b-value index, direction index 
/// for the current diffusion scan.
// ===========================================================================
void SBBDiffusion_Base::getActualCounter(long  lDiffLoopCounter, long &lActualAverageCounter, long &lActualBValueCounter, long &lActualDirectionCounter)
{
    /*  What for hell is happening in here?
    Example:
    Orthogonal mode, b values = [0, 300, 1000]
    Local averages for individual b value = [1, 2, 3] (not including global average in MrProt->average() )
    Global average in MrProt->average() = 2;
    m_NoOfWeigthings has the value 3
    getTotalScans() returns 16 (b=0 is scanned only once in "dir" dimension)

    The measurement sequence is  0p - 0r - 0s - 300p - 300r - 300s - 1000p - 1000r - 1000s
    We really measure
    Global AVE 0:
    Local AVE 0:   0               1      2      3       4       5       6
    Local AVE 1:                   7      8      9       10      11      12
    Local AVE 2:                                         13      14      15
    Global AVE 1:
    Local AVE 0:   0               1      2      3       4       5       6
    Local AVE 1:                   7      8      9       10      11      12
    Local AVE 2:                                         13      14      15

    !!!lDiffLoopCounter is recounted from 0 for each global average

    */

    // In q-space mode, only the direction index is relevant
    // => ignore local averages and b-values
    if(m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
    {
        lActualAverageCounter   = 0;
        lActualBValueCounter    = 0;
        lActualDirectionCounter = lDiffLoopCounter;

        return;
    }

    const long lS = lDiffLoopCounter + 1; // Diffusion scan counter starts with 0 as first valid scan
    long       lIndexAve;
    long       lIndexBValue = 0;
    long       lIndexDir    = 0;
    long       lCounter     = 0;
    const long lNoOfBValues = static_cast<long>(m_vdBValues.size());

    // Loop all dimensions to find out the indexes of average, bValue 
    // and direction which are mapping to the lDiffLoopCounter
    for(lIndexAve = 0; lIndexAve < m_lMaxValueInAveArray; lIndexAve++)
    {
        for(lIndexBValue = 0; lIndexBValue < lNoOfBValues; lIndexBValue++)
        {
            if(m_vlLocalAverages[lIndexBValue] >= (lIndexAve + 1))
            {
                long lDir = m_lDirections;
                if(m_vdBValues[lIndexBValue] <= ONE_FOR_THREE_THRESHOLD)
                {
                    lDir = 1;
                }

                for(lIndexDir = 0; lIndexDir < lDir; lIndexDir++)
                {
                    lCounter ++;
                    if(lS == lCounter) break;
                }
            }
            if(lS == lCounter) break;
        }
        if(lS == lCounter) break;
    }

    lActualAverageCounter   = lIndexAve;
    lActualBValueCounter    = lIndexBValue;
    lActualDirectionCounter = lIndexDir;
}


// ===========================================================================
///	Returns the value of the average index for the current diffusion scan.
// ===========================================================================
long SBBDiffusion_Base::getActualAverageCounter(long lDiffLoopCounter)
// ===========================================================================
{
    long lActualAverageCounter   = 0;
    long lActualBValueCounter    = 0;
    long lActualDirectionCounter = 0;

    getActualCounter(lDiffLoopCounter, lActualAverageCounter, lActualBValueCounter, lActualDirectionCounter);

    return lActualAverageCounter;
}

// ===========================================================================
///	Returns the value of the b-value index for the current diffusion scan.
// ===========================================================================
long SBBDiffusion_Base::getActualBValueCounter(long lDiffLoopCounter)
// ===========================================================================
{
    long lActualAverageCounter   = 0;
    long lActualBValueCounter    = 0;
    long lActualDirectionCounter = 0;

    getActualCounter(lDiffLoopCounter, lActualAverageCounter, lActualBValueCounter, lActualDirectionCounter);

    return lActualBValueCounter;
}


// ===========================================================================
///	Returns the value of the direction index for the current diffusion scan.
// ===========================================================================
long SBBDiffusion_Base::getActualDirectionCounter(long lDiffLoopCounter)
// ===========================================================================
{
    long lActualAverageCounter   = 0;
    long lActualBValueCounter    = 0;
    long lActualDirectionCounter = 0;

    getActualCounter(lDiffLoopCounter, lActualAverageCounter, lActualBValueCounter, lActualDirectionCounter);

    return lActualDirectionCounter;
}


// ===========================================================================
///	This method sets the maximum amplitude used by the diffusion gradients of this SBB.

/**   The argument will be written into the member variable m_dAmpl.

If the maximum gradient amplitude of the system is exceeded, a warning
will be generated and the value will be set to
SysProperties::getGradMaxAmplAbsolute().
*/
// ===========================================================================
void SBBDiffusion_Base::setMaxAmplitude(double dMaxAmplitude)
// ===========================================================================
{
    if(dMaxAmplitude > SysProperties::getGradMaxAmplAbsolute())
    {
        SEQ_TRACE_WARN.print("SBBDiffusion_Base::setMaxAmplitude WARNING: "
                   "Requested GradAmpl=%0.1f clipped to MaxAmplAbs=%0.1f",
                   dMaxAmplitude, SysProperties::getGradMaxAmplAbsolute());
        dMaxAmplitude = SysProperties::getGradMaxAmplAbsolute();
    }

    if(dMaxAmplitude <= 0.)
    {
        SEQ_TRACE_ERROR.print("SBBDiffusion_Base::setMaxAmplitude ERROR: "
                   "No or negative MaxAmplitude specified: %0.1f. Using 5.",
                   dMaxAmplitude);
        dMaxAmplitude = 5.0;
    }

    if(dMaxAmplitude != m_dAmpl)
    {
        m_dAmpl = dMaxAmplitude;
        // Reset preparation status
        resetPrepared();
    }
}


// ===========================================================================
///	This method sets the minimum gradient rise time used by this SBB
/**
The argument will be written into the member variable m_dMinRiseTime.

If the maximum gradient slew rate of the system is exceeded, a warning
will be generated and the value will be set to
SysProperties::getGradClipRiseTime().
*/
// ===========================================================================
void SBBDiffusion_Base::setMinRiseTime(double dMinRiseTime)
// ===========================================================================
{
    /*  SEQ_TRACE_ALWAYS.print(
    "SBBDiffusion_Base::setMinRiseTime: Requested MinRiseTime %0.0f (line %d)",
    dMinRiseTime, __LINE__);
    */
    if(dMinRiseTime < SysProperties::getGradClipRiseTime())
    {
        SEQ_TRACE_ERROR.print("SBBDiffusion_Base::setMinRiseTime ERROR: Requested MinRiseTime %0.0f exceeds GradClipRiseTime=%f.",
                   dMinRiseTime, SysProperties::getGradClipRiseTime());
        dMinRiseTime = SysProperties::getGradClipRiseTime();
    }

    if(dMinRiseTime != m_dMinRiseTime)
    {
        m_dMinRiseTime = dMinRiseTime;
        resetPrepared();
    }
}

// ===========================================================================
///   Set the noise threshold used by the ICE program for calculating the ADC maps.
/**
The value will be writen to the member variable m_lNoiseThreshold.
A negative input value will be changed to a default of 20.
*/
// ===========================================================================
void SBBDiffusion_Base::setNoiseThreshold(long lNoiseThreshold)
// ===========================================================================
{
    if(lNoiseThreshold < 0)
    {
        SEQ_TRACE_ERROR.print("SBBDiffusion_Base::setNoiseThreshold ERROR: "
                   "negative threshold specified. Using 20.");

        lNoiseThreshold = 20;
    }

    m_lNoiseThreshold = lNoiseThreshold;
}


// ===========================================================================
///   Set the TE array index
/**
Some more sophisticated sequences may have several TE times.
This variable is used to select the TE array element which is relevant
for timing calculation. By default, the timing calculation of this SBB
is based on the first TE time of the array, i.e. MrProt.TE()[0].

The value will be writen to the member variable m_iTEArrayIndex.
If the index is outside the rage [0..9], zero will be used.
*/
// ===========================================================================
void SBBDiffusion_Base::setTEArrayIndex(int iIndex)
// ===========================================================================
{
    if((iIndex < 0) || (iIndex > 9))
    {
        SEQ_TRACE_ERROR.print("SBBDiffusion_Base::setTEArrayIndex ERROR: "
                   "Invalid array index specified. Using 0.");
        iIndex = 0;
    }

    if(iIndex != m_iTEArrayIndex)
    {
        m_iTEArrayIndex = iIndex;
        resetPrepared();
    }
}


// ===========================================================================
/*
\author   Michael.Zwanger@med.siemens.de

\brief This function transforms a XYZ vector to a PRS vector.

If m_Didi specifies a vector set in XYZ coordinates, this method transforms
this vector into a PRS vector by multiplying the XYZ vector with the inverse
(=transposed) rotation matrix of the current slice.

The function works also if the current vector set is already specified in
PRS coordinates.

*/
// ===========================================================================
void SBBDiffusion_Base::DidiXYZ2PRS(
    DiffusionDirections* pDidi, sSLICE_POS* pSLC,           long lDirectionCounter,     double *dDidiP,             double *dDidiR,             double *dDidiS
                                ) const
{
    assert(dDidiP);
    assert(dDidiR);
    assert(dDidiS);

    if (pDidi->getCoordinateSystem() == MrProtocolData::DIFFDIR_CS_PRS)
    {
        // If current vector set is already specified in PRS, we are already done:
        *dDidiP = pDidi->getX(lDirectionCounter);
        *dDidiR = pDidi->getY(lDirectionCounter);
        *dDidiS = pDidi->getZ(lDirectionCounter);
    }
    else
    {
        // Step 2a: For (eCoordinateFrame == XYZ), the sequence uses the unity rotation matrix.
        // Therefore transform vector from XYZ-Hardware coordinate system into Phase-Read-Slice system
        // This is done by applying the inverse roation matrix 'm_sROT_MATRIX'.
        // As the rotation matrix is orthogonal, the inverse matrix is equal to the transponed one.
        const double dX = pDidi->getX(lDirectionCounter);
        const double dY = pDidi->getY(lDirectionCounter);
        const double dZ = pDidi->getZ(lDirectionCounter);
        // Matrix multiplication with the inverse = transposed matrix
        // regarding the indices, refer to line 278 ff. in file
        // \n4\comp\Measurement\Sequence\libGSL\fGSLCalcRotMat.cpp@@\main\17,
        const sROT_MATRIX& rot_matrix = pSLC->getROT_MATRIX();

        *dDidiP =
            dX * rot_matrix.dMat[0][0] +
            dY * rot_matrix.dMat[1][0] +
            dZ * rot_matrix.dMat[2][0];
        *dDidiR =
            dX * rot_matrix.dMat[0][1] +
            dY * rot_matrix.dMat[1][1] +
            dZ * rot_matrix.dMat[2][1];
        *dDidiS =
            dX * rot_matrix.dMat[0][2] +
            dY * rot_matrix.dMat[1][2] +
            dZ * rot_matrix.dMat[2][2];
    }
}


// ===========================================================================
/*!
\author   PLM Neuro

\brief This function transforms a PRS vector to a XYZ vector.

If Didi specifies a vector set in PRS coordinates, this method transforms
this vector into a XYZ vector by multiplying the PRS vector with the
rotation matrix of the current slice.

The function works also if the current vector set is already specifed in
XYZ coorinates.

*/
// ===========================================================================
void SBBDiffusion_Base::DidiPRS2XYZ(
DiffusionDirections* pDidi, sSLICE_POS* pSLC,           long lDirectionCounter,     double *dDidiX,             double *dDidiY,             double *dDidiZ
                            ) const
{
    assert(dDidiX);
    assert(dDidiY);
    assert(dDidiZ);

    if(pDidi->getCoordinateSystem() == MrProtocolData::DIFFDIR_CS_XYZ)
    {
        // If current vector set is already specified in XYZ, we are already done:
        *dDidiX = pDidi->getX(lDirectionCounter);
        *dDidiY = pDidi->getY(lDirectionCounter);
        *dDidiZ = pDidi->getZ(lDirectionCounter);
    }
    else
    {
        // Transform vector from Phase-Read-Slice system into XYZ-Hardware coordinate system
        // This is done by applying the roation matrix 'm_sROT_MATRIX'.
        const double dX = pDidi->getX(lDirectionCounter);
        const double dY = pDidi->getY(lDirectionCounter);
        const double dZ = pDidi->getZ(lDirectionCounter);

        // Matrix multiplication with the rotation matrix
        const sROT_MATRIX& rot_matrix = pSLC->getROT_MATRIX();

        *dDidiX =
            dX * rot_matrix.dMat[0][0] +
            dY * rot_matrix.dMat[0][1] +
            dZ * rot_matrix.dMat[0][2];
        *dDidiY =
            dX * rot_matrix.dMat[1][0] +
            dY * rot_matrix.dMat[1][1] +
            dZ * rot_matrix.dMat[1][2];
        *dDidiZ =
            dX * rot_matrix.dMat[2][0] +
            dY * rot_matrix.dMat[2][1] +
            dZ * rot_matrix.dMat[2][2];
    }

    return;
}

// ===========================================================================
/*!
\author   PLM Neuro

\brief This function sets the actual loop counters for diffusion scans

The function is called at the beginning of the sequence runKernel and
sets the loop counters for IceProgramDiffusion2D and IceProgramDti.
Repetitions and diffusion loop counters are stored in member variables.
*/
// ===========================================================================

bool SBBDiffusion_Base::setLoopCounters
(
int       iAdjustmentScan,    /**< Imp:     Adjustment scan index                                                   */
long      lRepLoopCounter,    /**< Imp:     SeqLoop repetition counter--- global average number                     */
long      lDiffLoopCounter,   /**< Imp:     SeqLoop diffusion counter --- including local average of each b-value  */
sREADOUT* pADC                /**< Imp/Exp: Mdh                                                                     */
)
{
    m_lBValueCounter       = 0;                   // b-value index
    m_lDirectionCounter    = 0;                   // diffusion direction index
    m_lRepLoopCounter      = lRepLoopCounter;     // SeqLoop repetition index 
    m_lDiffLoopCounter     = lDiffLoopCounter;    // SeqLoop diffusion index
    long m_lAverageCounter = 0;

    // Set diffusion encoding properties
    if(iAdjustmentScan == 0)
    {
        // Imaging scans
        if (m_bThermalBalancing)
        {
            m_lAverageCounter = m_sDiffusionOrderInfo.getAverageIndex(m_lDiffLoopCounter);
            m_lBValueCounter = m_sDiffusionOrderInfo.getBValueIndex(m_lDiffLoopCounter);
            m_lDirectionCounter = m_sDiffusionOrderInfo.getDirectionIndex(m_lDiffLoopCounter);
        }
        else
        {
            getActualCounter(m_lDiffLoopCounter, m_lAverageCounter, m_lBValueCounter, m_lDirectionCounter);
        }

        // SEQ_TRACE_ALWAYS.print("DiffLoopCounter( %d),AveCounter(%d), BValueCounter ( %d), DirCounter( %d) ",m_lDiffLoopCounter, m_lAverageCounter, m_lBValueCounter, m_lDirectionCounter);
    }
    else
    {
        // Avoid division by zero
        if(m_AdjDidi.getNumberOfDirections() == 0)
        {
            SEQ_TRACE_ERROR.print("EROR: Division by zero");
            return false;
        }

        // Adjustment scans (for dynamic distortion correction)
        m_lBValueCounter = 0;

        switch(m_eDynDistMode)
        {
            case SEQ::DYN_DISTCORR_ADJ:
                // In adjust mode, m_iAdjScan runs from 1 to lAdjAverages * m_AdjDidi.getNumberOfDirections()
                m_lDirectionCounter = (m_iAdjScan - 1) % m_AdjDidi.getNumberOfDirections();
                break;
            case SEQ::DYN_DISTCORR_DIRECT:
                // In direct mode, there is only a single diffusion direction: use z-gradient (see didi.cpp)
                m_lDirectionCounter = 6;
                break;
            default:
            case SEQ::DYN_DISTCORR_NONE:
                SEQ_TRACE_ERROR.print("ERROR: direction index for adjustment scans undefined");
                return false;
        }
    }

    // Set MDH loop counters
    // All diffusion related Mdh loop counter entries for DWI and DTI Ice programs are set, 
    // adjustment scans (if there are any) for the dyanimc distortion correction
    // are appropriately considered.
    //
    // Image numbering: 
    // - Only Rep counter is used
    // - Rep counter runs over all adjustment- and imaging-scans
    // - Adjustment scans: Rep = 0 ... getNoOfAdjScans - 1 (different directions, b-values and averages)
    // - Imaging scans   : Rep = getNoOfAdjScans ...       (different directions, b-values and averages)

    if(m_iAdjScan == 0)
    {
        // Imaging scan

        // Patch repetitions counter
        pADC->getMDH().setCrep(static_cast<unsigned short>(m_lRepLoopCounter * getTotalScans() + m_lDiffLoopCounter + getNoOfAdjScans()));

        // Tag as imaging scan
        pADC->getMDH().deleteFromEvalInfoMask(MDH_TAGFLAG1);
    }
    else
    {
        // Adjustment scan

        // Patch repetitions counter
        pADC->getMDH().setCrep(static_cast<unsigned short>(m_iAdjScan - 1));

        // Tag as adjustment scan
        pADC->getMDH().addToEvalInfoMask(MDH_TAGFLAG1);
    }

    return true;
}



// ===========================================================================
/*!
\author   Thorsten.Feiweier@siemens.com

\brief This function exports the b-matrix for DICOM header (DTI Ice program).

The DICOM standard defines as "<b>Diffusion Gradient Orientation</b>" in tag
(0018,9089): "The direction cosines of the diffusion gradient vector with
respect to the patient."

This function exports the complete b-matrix to the DTI Ice program, which
takes care of the corresponding Dicom entries. The polarity of the diffusion
encoding direction is encoded in the sign of the diagonal matrix elements.
In addition, the nominal b-value and the direction index are passed to Ice.

Due to the limited capacity of this array and due to the fact it does not
provide a 'float' data type, the b-value information (which is a double)
is scaled to integer the following way:
\f[   MDHEntry = (unsigned short) (b-value + 16384.5)   \f]
Note that the ICE program must properly unpack this encoded information.

A similar approach is used for the actual diffusion gradient amplitude:
\f[   MDHEntry = (unsigned short) (amplitude * 10. + 16384.5)   \f]

Technically spoken, the following steps are performed:
-# Get slice position vector (in patient coordinates)
-# Get b-matrix information  (in xyz coordinates)
-# Transform b-matrix inforamtion into Patient Coordinate System
-# Rescale values to unsigned short type and export them to MDH
*/
// ===========================================================================


NLS_STATUS SBBDiffusion_Base::CalculateDicomHeaderInformation(
    const Slice &slice,           /**< Input: slice orientation vectors                */
    sSLICE_POS* pSLC,             /**< Input: getSliceIndex, m_sROT_MATRIX             */
    DiffusionDirections* pDidi,   /**< Input: Diffusion direction vector information   */
    double dAmpl,                 /**< Input: Diffusion gradient amplitude factor      */
    long lPolarity,               /**< Input: Direction sign                           */
    long lBValueCounter,          /**< Input: Current bvalue index                     */
    long lDirectionCounter,       /**< Input: Current direction vector index           */
    sREADOUT* pADC                /**< Output: Mdh.setIceProgramPara                   */
    )
{

    double dBxx = m_dBxx;
    double dByy = m_dByy;
    double dBzz = m_dBzz;
    double dBxy = m_dBxy;
    double dBxz = m_dBxz;
    double dByz = m_dByz;

    long   lNominalBValue = 0;
    bool   bDumpBMatrix   = m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_BMatrix", false);

    sROT_MATRIX rotMatrix;

    // rotMatrix transforms from GCS/DCS to PCS
    NLS_STATUS lStatus = CalculateGradientsToPCS(pSLC, rotMatrix);
    if(NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        return setNLSStatus(lStatus);
    }

    sROT_MATRIX covB;
    sROT_MATRIX matTmp;

    sROT_MATRIX covBRot;

    if(bDumpBMatrix)
    {
        SEQ_TRACE_ALWAYS.print("B-Matrix0 xx:%f xy:%f xz:%f yy:%f yz:%f zz:%f",
                   dBxx, dBxy, dBxz, dByy, dByz, dBzz);
    }


    covB.dMat[0][0] = dBxx;
    covB.dMat[0][1] = dBxy;
    covB.dMat[0][2] = dBxz;

    covB.dMat[1][0] = dBxy;
    covB.dMat[1][1] = dByy;
    covB.dMat[1][2] = dByz;

    covB.dMat[2][0] = dBxz;
    covB.dMat[2][1] = dByz;
    covB.dMat[2][2] = dBzz;

    // compute rotated Bmatrix covBRot = rotMatrix * covB * rotMatrix'
    // rotMatrix gives the rotation from PRS into PCS coordinate system
    // and thus covBRot is finally stored in PCS

    MatMultTrans(covB, rotMatrix, matTmp);       // matTmp = covB * rotMatrix'
    MatMult(rotMatrix, matTmp, covBRot);    //covBRot = rotMatrix * matTmp

    dBxx = covBRot.dMat[0][0];
    dBxy = covBRot.dMat[0][1];
    dBxz = covBRot.dMat[0][2];
    dByy = covBRot.dMat[1][1];
    dByz = covBRot.dMat[1][2];
    dBzz = covBRot.dMat[2][2];

    if(pDidi != nullptr)
    {
        // Encode diffusion direction polarity in b-matrix
        // ===============================================
        // Set sign of the b-matrix diagonal elements depending on the 'polarity' of the
        // diffusion encoding direction. Hint: for a well-defined b-matrix, all diagonal
        // elements are positive - this allows to use the sign for our purpose.
        //
        // Conventions:
        // - diffusion direction coordinates are provided in PCS coordinates
        // - polarity of Sag-direction is encoded as the sign of b-matrix[0][0]
        //               Cor-direction            as the sign of b-matrix[1][1]
        //               Tra-direction            as the sign of b-matrix[2][2]
        //
        // Before any matrix operations take place, Ice stores the signs of the
        // diagonal elements e1 = sign(b[0][0]), e2 = sign(b[1][1]), e3 = sign(b[2][2])
        // and applies a positive sign to every diagonal element of the matrix.
        // Given a diffusion encoding direction d = (d1, d2, d3) used in the sequence,
        // Ice then calculates a diffusion direction d' = (d1', d2', d3') from the
        // b-matrix. Either d' or -d' is similar (though not exactly identical!) to d.
        // The  correct sign of d' can be assigned by calculating the vector product
        // d' * (e1 * abs(d1'), e2 * abs(d2'), e3 * abs(d3')). If it is positive, d'
        // points in the same direction as d. If it is negative, d' has to be inverted.

        // 1st step: get diffusion direction in PRS (GCS) coordinates
        double dDidiP = 0.;
        double dDidiR = 0.;
        double dDidiS = 0.;

        DidiXYZ2PRS(pDidi, pSLC, lDirectionCounter, &dDidiP, &dDidiR, &dDidiS);

        if(lPolarity < 0)
        {
            // Invert diffusion direction
            dDidiP *= -1.;
            dDidiR *= -1.;
            dDidiS *= -1.;
        }

        // 2nd step: convert to PCS coordinates
        double dDidiSag = 0.;
        double dDidiCor = 0.;
        double dDidiTra = 0.;

        transformGCSToPCS     // Library: libGSL
            (
            dDidiP,                 // Import: Phase encoding component
            dDidiR,                 // Import: Readout component
            dDidiS,                 // Import: Slice selection component
            dDidiSag,               // Export: Sagittal component
            dDidiCor,               // Export: Coronal component
            dDidiTra,               // Export: Transverse component
            slice.normal().sag(),   // Import: Sagittal component of slice normal vector
            slice.normal().cor(),   // Import: Coronal component of slice normal vector
            slice.normal().tra(),   // Import: Transverse component of slice normal vector
            slice.rotationAngle()   // Import: Slice rotation angle ("swap Fre/Pha")
            );

        // 3rd step: encode direction

        // there might be very small negative numbers due to rounding errors (CHARM 443674)
        if(dBxx < 0.0 && dBxx > -1E-10)
            dBxx = 0.0;
        if(dByy < 0.0 && dByy > -1E-10)
            dByy = 0.0;
        if(dBzz < 0.0 && dBzz > -1E-10)
            dBzz = 0.0;

        if((dBxx < 0.) || (dByy < 0.) || (dBzz < 0.))
        {
            // If any of the diagonal elements is negative, the encoding will not
            // work. This should never happen, but who knows ...
            SEQ_TRACE_ERROR.print("ERROR: negative b-matrix diagonal element");
        }

        if(dDidiSag < 0.)
        {
            dBxx *= -1.;
        }
        if(dDidiCor < 0.)
        {
            dByy *= -1.;
        }
        if(dDidiTra < 0.)
        {
            dBzz *= -1.;
        }
    }   // End     if ( pDidi != NULL )

    if(bDumpBMatrix)
    {
        SEQ_TRACE_ALWAYS.print("B-Matrix1 xx:%f xy:%f xz:%f yy:%f yz:%f zz:%f",
                   dBxx, dBxy, dBxz, dByy, dByz, dBzz);

        SEQ_TRACE_ALWAYS.print("Rot \n                   %f %f %f\n                   %f %f %f\n                   %f %f %f\n",
                   rotMatrix.dMat[0][0], rotMatrix.dMat[0][1], rotMatrix.dMat[0][2],
                   rotMatrix.dMat[1][0], rotMatrix.dMat[1][1], rotMatrix.dMat[1][2],
                   rotMatrix.dMat[2][0], rotMatrix.dMat[2][1], rotMatrix.dMat[2][2]
                   );
    }

    // Store b-matrix in MDH
    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_BXX, static_cast<unsigned short>(dBxx + 16384.5));
    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_BYY, static_cast<unsigned short>(dByy + 16384.5));
    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_BZZ, static_cast<unsigned short>(dBzz + 16384.5));
    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_BXY, static_cast<unsigned short>(dBxy + 16384.5));
    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_BXZ, static_cast<unsigned short>(dBxz + 16384.5));
    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_BYZ, static_cast<unsigned short>(dByz + 16384.5));

    // Store nominal b-value and direction index in MDH
    if((m_eDiffusionMode == SEQ::DIFFMODE_FREE)
       || (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
       || (m_iAdjScan       != 0))
    {
        // Free mode, q-space mode and adjustment scans can include more than one b-value
        // => calculate nominal b-value from b-matrix
        long lBValueStep = m_lBValueInc_Limit;

        lNominalBValue = static_cast<long>((std::lround((m_dBxx + m_dByy + m_dBzz) / static_cast<double>(lBValueStep))) * lBValueStep);

        // SEQ_TRACE_ALWAYS.print("lBValueStep = %d, lNominalBValue = %d", lBValueStep, lNominalBValue) ;
    }
    else
    {
        // the nominal b value for image text is taken from the UI b value list
        lNominalBValue = std::lround(m_vdBValues[lBValueCounter]);
    }

    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_NOMINALB, static_cast<unsigned short>(static_cast<double>(lNominalBValue)+ 16384.5));   // nominal b value for image text
    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_DIRINDEX, static_cast<unsigned short>(lDirectionCounter));   // direction index
    pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_BVALUEINDEX, static_cast<unsigned short>(lBValueCounter));   // b-value index

    if(pDidi != nullptr)
    {
        // Calculate diffusion gradient direction in xyz coordinates
        // =========================================================
        double dDidiX = 0.;
        double dDidiY = 0.;
        double dDidiZ = 0.;

        DidiPRS2XYZ(pDidi, pSLC, lDirectionCounter, &dDidiX, &dDidiY, &dDidiZ);

        if(lPolarity < 0)
        {
            // Invert diffusion direction
            dDidiX *= -1.;
            dDidiY *= -1.;
            dDidiZ *= -1.;
        }

        // Store diffusion gradient in MDH
        // - Maximum amplitude:   +/- 1600.0 mT/m
        // - Amplitude precision:        0.1 mT/m
        pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_DIFFGRADX, static_cast<unsigned short>(dAmpl * dDidiX * 10. + 16384.5)); // 10 * x-gradient [mT/m]
        pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_DIFFGRADY, static_cast<unsigned short>(dAmpl * dDidiY * 10. + 16384.5)); // 10 * x-gradient [mT/m]
        pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_DIFFGRADZ, static_cast<unsigned short>(dAmpl * dDidiZ * 10. + 16384.5)); // 10 * x-gradient [mT/m]
    }
    else
    {
        // No explicit submission of diffusion gradient if pDidi is not provided
        pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_DIFFGRADX, 0);
        pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_DIFFGRADY, 0);
        pADC->getMDH().setIceProgramPara(DIFFUSIONINTSTORAGE_DIFFGRADZ, 0);
    }

    if(m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_ICE", false))
    {
        uint16_t iI;

        SEQ_TRACE_ALWAYS.print("MDH loop counter: REP=%03d ", pADC->getMDH().getCrep());

        for(iI = 0; iI < DIFFUSIONINTSTORAGE_SIZE; iI++)
        {
            switch(iI)
            {
                case DIFFUSIONINTSTORAGE_BXX:
                case DIFFUSIONINTSTORAGE_BYY:
                case DIFFUSIONINTSTORAGE_BZZ:
                case DIFFUSIONINTSTORAGE_BXY:
                case DIFFUSIONINTSTORAGE_BXZ:
                case DIFFUSIONINTSTORAGE_BYZ:
                case DIFFUSIONINTSTORAGE_NOMINALB:
                {
                    SEQ_TRACE_ALWAYS.print("IceProgramPara[%d]=%d   (i.e. b=%d)", iI, pADC->getMDH().getIceProgramPara(iI), pADC->getMDH().getIceProgramPara(iI) - 16384);
                    break;
                }
                case DIFFUSIONINTSTORAGE_DIRINDEX:
                case DIFFUSIONINTSTORAGE_BVALUEINDEX:
                {
                    SEQ_TRACE_ALWAYS.print("IceProgramPara[%d]=%d", iI, pADC->getMDH().getIceProgramPara(iI));
                    break;
                }
                case DIFFUSIONINTSTORAGE_DIFFGRADX:
                case DIFFUSIONINTSTORAGE_DIFFGRADY:
                case DIFFUSIONINTSTORAGE_DIFFGRADZ:
                {
                    SEQ_TRACE_ALWAYS.print("IceProgramPara[%d]=%d   (i.e. g=%d)", iI, pADC->getMDH().getIceProgramPara(iI), (pADC->getMDH().getIceProgramPara(iI) - 16384) / 10);
                    break;
                }
                default:
                    // Do nothing
                    break;
            }
        }
    }

    setNLSStatus(MRI_SBB_SBB_NORMAL);
    return true;
}


// ===========================================================================
/*!
\author   Thorsten.Feiweier@siemens.com

\brief This ensures that TE is set to a value compatible with MrUILink

Note: we assume something about TE-limit-handling in MrUILink!
*/
long SBBDiffusion_Base::lCalcTEOnInc
(
long lTENotOnInc                    /**< Imp: TE              */
)
{
    long   lInc = m_lTEInc_Limit;

    if(lTENotOnInc >   10000) lInc *= 10;
    if(lTENotOnInc > 1000000) lInc *= 10;

    if((lTENotOnInc % lInc) == 0)
    {
        return lTENotOnInc;
    }

    const auto dInc = static_cast<double>(lInc);

    return std::lround(dInc * ceil(static_cast<double>(lTENotOnInc) / dInc));
}


// ===========================================================================
/*!
\author   Stefan.Huwer@med.siemens.de

\brief Calculate the transformation from GradientCoordinatesystem depending on the
coordinate system used for gradients (PRS=GCS or XYZ=DCS)


\return NLS_STATUS NLS_SUCCESS

*/
// ===========================================================================
NLS_STATUS SBBDiffusion_Base::CalculateGradientsToPCS(
    sSLICE_POS* pSLC,
    sROT_MATRIX& rotMatrix)
{
    // get the rotation matrix from DCS to PCS ...
    // -------------------------------------------

    LROTMAT          sPLAMatrix;
    const NLS_STATUS lStatus = fGSLCalcPLAMatrix(
        m_iPatDirection,
        m_iPatPosition,
        &(sPLAMatrix));

    if(NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        return setNLSStatus(lStatus);
    }

    // B_matrix is internally stored in PRS coordinate system thus 
    // a transformation  into the PCS coordinate system is necessary
    // to get the dicom attribute B_matrix in the correct  coordinate system
    //
    //      Rotation from PRS to XYZ via pSLC_sROT_MATRIX
    //      Rotation from XYZ to PPCS via sPLAMatrix
    //
    // rotMatrix = sPLAMatrix * pSLC_sROT_MATRIX
    const sROT_MATRIX& rot_matrix = pSLC->getROT_MATRIX();

    for(int r=0;r<3;r++)
    {
        for(int c=0;c<3;c++)
        {
            rotMatrix.dMat[r][c] =
                static_cast<double>(sPLAMatrix.lMat[0][r]) * rot_matrix.dMat[0][c] +
                static_cast<double>(sPLAMatrix.lMat[1][r]) * rot_matrix.dMat[1][c] +
                static_cast<double>(sPLAMatrix.lMat[2][r]) * rot_matrix.dMat[2][c];
        }
    }

    return NLS_SUCCESS;
}

// ===========================================================================
/*!
\author   Stefan.Huwer@med.siemens.de

\brief multiplies the matrices given: matC = matA * matB

*/
// ===========================================================================
void SBBDiffusion_Base::MatMult(sROT_MATRIX const& matA, sROT_MATRIX const& matB, sROT_MATRIX& matC) const
{
    for(int r=0;r<3;r++)
    {
        for(int c=0;c<3;c++)
        {
            matC.dMat[r][c] =
                matA.dMat[r][0] * matB.dMat[0][c] +
                matA.dMat[r][1] * matB.dMat[1][c] +
                matA.dMat[r][2] * matB.dMat[2][c];
        }
    }
}

// ===========================================================================
/*!
\author   Stefan.Huwer@med.siemens.de

\brief multiplies the matrices given: matC = matA * matB'

*/
// ===========================================================================
void SBBDiffusion_Base::MatMultTrans(sROT_MATRIX const& matA, sROT_MATRIX const& matB, sROT_MATRIX& matC) const
{
    for(int r=0;r<3;r++)
    {
        for(int c=0;c<3;c++)
        {
            matC.dMat[r][c] =
                matA.dMat[r][0] * matB.dMat[c][0] + //transposed matB
                matA.dMat[r][1] * matB.dMat[c][1] +
                matA.dMat[r][2] * matB.dMat[c][2];
        }
    }
}


// ===========================================================================
/*!
\brief Calculation of maximum amplitude for diffusion encoding gradients

The following strategy applies:

1. Check if the GPA allows to run the gradient event series
<acquisition> - <diffusion> once with the absolute maximum amplitude.
If this is the case: continue with step 3.
Note: Only one gradient axis is considered in this step. Diffusion
encoding is assumed to exclusively take place on this axis.

2. Perform a binary search to identify the maximum allowable diffusion
gradient amplitude for running the event series <acquisition> - <diffusion>
once.
Note: Only one gradient axis is considered in this step. Diffusion
encoding is assumed to exclusively take place on this axis.

3. Calculate minimum pause duration required to run the event series
<acquisition> - <diffusion> - <pause> with the calculated diffusion
gradient amplitude forever.
Note: Only one gradient axis is considered in this step. Diffusion
encoding is assumed to exclusively take place on this axis.

4. Calculate minimum pause duration required to run the event series
<acquisition> - <diffusion> - <pause> with the calculated diffusion
gradient amplitude forever. Here, the focus is on long (thermal)
time constants.
If the pause calculated here is longer than the one calculated in
step 3, we are finished.
Note: All gradient axes are considered in this step. All diffusion
directions are considered (distribution of thermal energy
among the different axes).

5. Perform a binary search to identify the minimum allowable pause
duration for running the event series <acquisition> - <diffusion> -
<pause> forever.
Note: Only one gradient axis is considered in this step. All diffusion
directions are considered (reduced average load on the axis
considered)

The calculated maximum diffusion gradient amplitude and the minimum
pause duration are returned.

CalcMaximumAmplitude and CalcTRIncrement are the wrappers used
internally to calculate only the maximum amplitude or the required
TR increment, respectively.
*/
// ===========================================================================
bool SBBDiffusion_Base::CalcMaximumAmplitude(MrProt &rMrProt, double &dMaxAmplitude)
{
    long         lTRIncDummy = 0;   // Dummy value
    const double dUpperLimit = -1.; // Indicate maximum amplitude calculation, no TR increment calculation

    const bool status = CalcMaxAmpAndTRInc(dMaxAmplitude, lTRIncDummy, dUpperLimit);

    // 20140320 DP: MR_00443009: in case of ZOOMit, the maximum grad amplitude will be reduced due to duty cycle problems. The reduction factor is empirical so far.
#if defined ZOOM_2DRF
    if(m_bZoomedExcitationIsUsed)
    {
        dMaxAmplitude = m_dDiffGradMaxAmpReduction * dMaxAmplitude;
    }
#endif
    
    SMSProperties::calcMaximumAmplitude(rMrProt, dMaxAmplitude);

    return status;
}

bool SBBDiffusion_Base::CalcTRIncrement(double dAmplitude, long &lTRIncrement)
{
    if(dAmplitude < 0.)
    {
        SEQ_TRACE_ERROR.print("Error: negative amplitude %fmT/m. Setting TR increment to 0.",  dAmplitude);
        lTRIncrement = 0;

        return false;
    }

    double dMaxAmpDummy = 0.; // Dummy value, calculate TR increment only

    return CalcMaxAmpAndTRInc(dMaxAmpDummy, lTRIncrement, dAmplitude);
}

bool SBBDiffusion_Base::CalcMaxAmpAndTRInc(double &dMaxAmplitude, long &lTRIncrement, double dUpperLimit)
{

    GPABalance::GPABalanceValue sBalanceIn;
    GPABalance::GPABalanceValue sBalanceOut;
    GPABalance::GPABalanceValue sBalanceMin;
    GPABalance::GPABalanceValue sBalance;

    double     dAmplitude    = 0.;
    double     dNewAmplitude = 0.;
    double     dUpperAmp     = 0.;
    double     dLowerAmp     = 0.;
    long       lResult       = 0;
    long       lDuration     = 0;
    bool       bReturn       = true;
    bool       bCalcMaxAmp   = false;
    bool       bCalcTRInc    = false;


    // Scaling of calculated diffusion gradient amplitude
    // (maximum = Gmax)
    auto     dGradScale = m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/diff_grad_scale", 1.0);


    const double dAmpSearchLimit   = 0.5;    // Search maximum amplitude down to this precision [mT/m]
    // For Gmax = 45mT/m, this yields a maximum of 7 iterations.

    const long   lTRIncSearchLimit = 2000;   // Search minimum TR increment down to this precision [us]


    // Initialize return values
    dMaxAmplitude = 0.;
    lTRIncrement  = 0;

    // Set upper limit for diffusion gradient amplitude (and start value for the binary search)
    if(dUpperLimit <= 0.)
    {
        // Set upper limit depending on GPA restrictions
        dUpperAmp = std::min(m_sBalanceDiff.dGetMaximumAmplitude(), m_dMaxAmpl);

        // Calculate maximum amplitude only
        bCalcMaxAmp = true;
        bCalcTRInc  = false;
    }
    else
    {
        // Set maximum amplitude provided as input parameter
        // If it is possible to play out the gradient events with this amplitude,
        // it will be used - even if higher gradient amplitudes would in principle
        // be possible. This feature can be used to incidentally apply a reduced
        // gradient amplitude and trade this for a shorter TR increment.
        dUpperAmp = dUpperLimit;

        // Calculate TR increment only
        // (However, it is checked whether it is possible to apply the given amplitude)
        bCalcMaxAmp = false;
        bCalcTRInc  = true;
    }

    dAmplitude = dUpperAmp;

    // ------------------------------------------------------------------------
    // Check whether we are supposed to use the GPA load model. If not, the stored default
    // value for the diffusion gradient amplitudes will be used. lGetStatus will return
    // an error status (and dump a message) e.g. if the current GPA is not supported.

    if(!m_bUseGPABalance || (m_sBalanceDiff.lGetStatus() != 0))
    {
        setTRIncrement(0);

        if(m_bUseGPABalance && (m_sBalanceDiff.lGetStatus() != 0))
        {
            // We are supposed to use the GPA model, but something failed 
            // => dump status of GPA model
            m_sBalanceDiff.lGetStatus(true);

            // Since something failed, return a meaningless amplitude
            dMaxAmplitude = 0.;

            return false;
        }
        else
        {
            // GPA model not used. Return the original maximum amplitude
            // (see SBBDiffusion base class).
            dMaxAmplitude = m_dMaxAmpl;

            return true;
        }
    }

    // ------------------------------------------------------------------------
    // 1.
    // ------------------------------------------------------------------------
    // Prepare diffusion gradient events for balance calculations (m_sBalanceDiff).
    // Note: The balance of the readout module (m_sBalanceAcq) has been provided
    //       by the sequence
    if(!prepGPALoadDiff(dAmplitude))
    {
        SEQ_TRACE_ERROR.print("ERROR: preparation of diffusion GPA load failed");

        // Since something failed, return a meaningless amplitude
        dMaxAmplitude = 0.;

        return false;
    }

    // Prepare compensation gradient events for balance calculations
    // Note: The compensation gradient should be already prepared
    if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
    {
        double dECCGAmplitude = m_CompGrad.getdMaxAmplitude() * dAmplitude / m_dMaxAmpl;
        if (!m_CompGrad.prepGPALoad(dECCGAmplitude))
        {
            SEQ_TRACE_ERROR.print("ERROR: preparation of ECCG GPA load failed");

            // Since something failed, return a meaningless amplitude
            dMaxAmplitude = 0.;

            return false;
        }
    }

    // No history
    GPABalance::ResetBalanceValue(sBalanceIn);
    GPABalance::ResetBalanceValue(sBalance);
    GPABalance::ResetBalanceValue(sBalanceMin);

    // Calculate total time of gradient events
    if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
    {
        lDuration
            = m_sBalanceDiff.lGetDuration() + m_sBalanceAcq.lGetDuration() + m_CompGrad.getGPALoad().lGetDuration();
    }
    else
    {
        lDuration = m_sBalanceDiff.lGetDuration() + m_sBalanceAcq.lGetDuration();
    }

    // ------------------------------------------------------------------------
    // Before we start to search: check whether the maximum possible gradient
    // amplitude can be used. If so, there is no need to search. If not, we
    // have an upper limit for the binary search.

    // It is important to perform the check of the all modules in the
    // right (worst case) order. The repeatability check (see below) 
    // calculates the required TR increment based on the balance of 
    // the concatenated modules. Typically, the maximum balance can
    // be observed immediately after the diffusion module. Thus, the
    // concatenated balance has to be calculated for the order
    // [Readout - CompGrad - Diff].

    // Check readout module
    lResult = m_sBalanceAcq.lCalc(GPABALANCE_X_AXIS, sBalanceIn, sBalanceOut);

    // Check diffusion module
    if(lResult == 0)
    {
        sBalanceIn = sBalanceOut;

        if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
        {
            // Check compensation module
            lResult = m_CompGrad.getGPALoad().lCalc(GPABALANCE_X_AXIS, sBalanceIn, sBalanceOut);
            if (lResult == 0)
            {
                sBalanceIn = sBalanceOut;
                lResult    = m_sBalanceDiff.lCalc(GPABALANCE_X_AXIS, sBalanceIn, sBalanceOut);
            }
        }
        else
        {
            lResult = m_sBalanceDiff.lCalc(GPABALANCE_X_AXIS, sBalanceIn, sBalanceOut);
        } 
    }

    // Are we supposed to calculate the maximum amplitude? 
    if(!bCalcMaxAmp && (lResult != 0))
    {
        // No - but the given amplitude appears too high => error
        SEQ_TRACE_ERROR.print("diffusion gradient amplitude %fmT/m is too high!", dAmplitude);

        return false;
    }

    // ------------------------------------------------------------------------
    // 2.
    // ------------------------------------------------------------------------
    // Search for maximum amplitude (precision given above)
    //
    // Strategy: binary search
    // Note:     only diffusion gradient amplitudes are scaled
    //           (spoiler amplitudes are constant)

    if(lResult == 0)
    {
        // It is possible to apply the diffusion module and a succeeding readout module at least once.
        // This means that with a sufficiently long pause, it is possible to apply this forever.

        // Store balance of this series of gradient events
        sBalanceMin = sBalanceOut;
    }
    else
    {
        // Maximum amplitude is not possible - start binary search
        while((dUpperAmp - dLowerAmp) > dAmpSearchLimit)
        {
            dNewAmplitude = (dUpperAmp + dLowerAmp) / 2.;

            // Scale amplitude and remember error status
            if(!scaleGPALoadDiff(dNewAmplitude / dAmplitude))
            {
                SEQ_TRACE_ERROR.print("scaling error");

                // Since something failed, return a meaningless amplitude
                dMaxAmplitude = 0.;

                return false;
            }

            if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
            {
                // Scale amplitude and remember error status
                if (!m_CompGrad.scaleGPALoad(dNewAmplitude / dAmplitude))
                {
                    SEQ_TRACE_ERROR.print("scaling error");

                    // Since something failed, return a meaningless amplitude
                    dMaxAmplitude = 0.;

                    return false;
                }
            }

            // No history
            GPABalance::ResetBalanceValue(sBalanceIn);

            // Check readout module
            lResult = m_sBalanceAcq.lCalc(GPABALANCE_X_AXIS, sBalanceIn, sBalanceOut);

            // Check diffusion module
            if(lResult == 0)
            {
                sBalanceIn = sBalanceOut;
                if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
                {
                    // Check compensation module
                    lResult = m_CompGrad.getGPALoad().lCalc(GPABALANCE_X_AXIS, sBalanceIn, sBalanceOut);
                    if (lResult == 0)
                    {
                        sBalanceIn = sBalanceOut;
                        lResult    = m_sBalanceDiff.lCalc(GPABALANCE_X_AXIS, sBalanceIn, sBalanceOut);
                    }
                }
                else
                {
                    lResult = m_sBalanceDiff.lCalc(GPABALANCE_X_AXIS, sBalanceIn, sBalanceOut);
                }  
            }

            // Update search boundaries depending on result
            if(lResult != 0)
            {
                // Further reduction necessary
                dUpperAmp = dNewAmplitude;
            }
            else
            {
                // Further increase possible
                dLowerAmp = dNewAmplitude;

                // Remember balance change of highest amplitude that works
                sBalanceMin = sBalanceOut;
            }

            dAmplitude = dNewAmplitude;
        }

        // Remember the lowest amplitude that works
        dAmplitude = dLowerAmp;
    }

    // ------------------------------------------------------------------------
    // 3.
    // ------------------------------------------------------------------------
    // With this amplitude, it is possible to apply at least one repetition ...
    if(dUpperLimit <= 0.)
    {
        // Option: overload the gradient system by providing a dGradScale > 1
        dMaxAmplitude = dAmplitude * dGradScale;    // Export
    }
    else
    {
        // No overload if a gradient amplitude has been explicitly provided
        dMaxAmplitude = dAmplitude;                 // Export
    }

    // Are we supposed to calculate the required TR increment? 
    if(!bCalcTRInc)
    {
        // No - we have finished

        return true;
    }

    // ... however, in order to be able to run 'infinitely', an additional cooling pause is required.
    // This dynamic cooling pause has to be applied by the sequence once for each effective TR.
    lTRIncrement = m_sBalanceDiff.lCalcPause(sBalanceMin, lDuration);      // Export

    // ------------------------------------------------------------------------
    // 4.
    // ------------------------------------------------------------------------
    // Before we try to reduce the pause duration by taking into account
    // the different diffusion directions, an absolute minimum pause
    // duration based on thermal long-term effects (mean square gradient
    // amplitude limitation) is calculated.

    GPABalance               DiffMemBalance;
    GPABalance               CompGradMemBalance;
    GPABalance::GrmsValue    sGrms;

    int          iAmplAvgMax         = 0;
    long         lI                  = 0;
    long         lDirections         = m_Didi.getNumberOfDirections();
    long         lTotalScans         = 0;
    long         lMinPause           = 0;
    double       dAmplAvgX           = 0.;
    double       dAmplAvgY           = 0.;
    double       dAmplAvgZ           = 0.;

    lTotalScans = m_bThermalBalancing ? getTotalScans(true) : lDirections;

    // For long term thermal effects, all axes have to be considered
    if(!prepGPALoadDiff(dAmplitude, dAmplitude, dAmplitude))
    {
        SEQ_TRACE_ERROR.print("preparation of diffusion GPA load failed");

        // Since something failed, return a meaningless amplitude
        dMaxAmplitude = 0.;

        return false;
    }

    // Store prepared balance
    DiffMemBalance = m_sBalanceDiff;

    if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
    {
        double dECCGAmplitude = m_CompGrad.getdMaxAmplitude() * dAmplitude / m_dMaxAmpl;
        if (!m_CompGrad.prepGPALoad(dECCGAmplitude, dECCGAmplitude, dECCGAmplitude))
        {
            SEQ_TRACE_ERROR.print("ERROR: preparation of ECCG GPA load failed");

            // Since something failed, return a meaningless amplitude
            dMaxAmplitude = 0.;

            return false;
        }
        // Store prepared balance
        CompGradMemBalance = m_CompGrad.getGPALoad();
    }
    // No history
    GPABalance::ResetRMSValue(sGrms);

    // ------------------------------------------------------------------------
    // Identify the axis with the highest average amplitude 
    // (For isotropic diffusion directions, all axes should be similar)
    for(lI = 0; lI < lDirections; ++lI)
    {
        dAmplAvgX += fabs(m_Didi.getX(lI));
        dAmplAvgY += fabs(m_Didi.getY(lI));
        dAmplAvgZ += fabs(m_Didi.getZ(lI));
    }

    if((dAmplAvgX >= dAmplAvgY) && (dAmplAvgX >= dAmplAvgZ))
    {
        iAmplAvgMax = 0;    // X
    }
    else if((dAmplAvgY >= dAmplAvgX) && (dAmplAvgY >= dAmplAvgZ))
    {
        iAmplAvgMax = 1;    // Y
    }
    else
    {
        iAmplAvgMax = 2;    // Z
    }

    // maximum b-value
    double dMaxBValue = (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE) ? m_dQSpaceMaxBValue : *std::max_element(m_vdBValues.begin(), m_vdBValues.end());

    // Ensure a value > 0 (avoid division by zero)
    dMaxBValue = std::max<double>(1., dMaxBValue);

    for(lI = 0; lI < lTotalScans; ++lI)
    {
        double dAmpScaleX = m_Didi.getX(lI);
        double dAmpScaleY = m_Didi.getY(lI);
        double dAmpScaleZ = m_Didi.getZ(lI);

        if (m_bThermalBalancing)
        {
            double dScaleForBValue = (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE) ? 1.0 : sqrt(m_vdBValues[m_sDiffusionOrderInfo.getBValueIndex(lI)] / dMaxBValue);

            dAmpScaleX = m_Didi.getX(m_sDiffusionOrderInfo.getDirectionIndex(lI)) * dScaleForBValue;
            dAmpScaleY = m_Didi.getY(m_sDiffusionOrderInfo.getDirectionIndex(lI)) * dScaleForBValue;
            dAmpScaleZ = m_Didi.getZ(m_sDiffusionOrderInfo.getDirectionIndex(lI)) * dScaleForBValue;
        }

        // Restore reference diffusion module (which gets modified by scaleGPALoadDiff below)
        m_sBalanceDiff = DiffMemBalance;

        if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
        {
            // Restore reference compensation module (which gets modified by scaleGPALoadDiff below)
            m_CompGrad.setGPALoad(CompGradMemBalance);
        }

        // Make sure that highest average amplitude gets applied to x-axis
        // (worst case assumption - relevant only for anisotropic vector sets)
        switch(iAmplAvgMax)
        {
            case 0:
            default:
                break;
            case 1:
            {
                // Exchange X with Y amplitude
                double dTempValue = dAmpScaleX;
                dAmpScaleX = dAmpScaleY;
                dAmpScaleY = dTempValue;
                break;
            }
            case 2:
            {
                // Exchange X with Z amplitude
                double dTempValue = dAmpScaleX;
                dAmpScaleX = dAmpScaleZ;
                dAmpScaleZ = dTempValue;
                break;
            }
        }

        // Scale according to actual diffusion direction
        if(!scaleGPALoadDiff(dAmpScaleX, dAmpScaleY, dAmpScaleZ))
        {
            SEQ_TRACE_ERROR.print("scaling error");

            // Since something failed, return a meaningless amplitude
            dMaxAmplitude = 0.;

            return false;
        }

        if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
        {
            // Scale according to actual diffusion direction
            if (!m_CompGrad.scaleGPALoad(dAmpScaleX, dAmpScaleY, dAmpScaleZ))
            {
                SEQ_TRACE_ERROR.print("scaling error");

                // Since something failed, return a meaningless amplitude
                dMaxAmplitude = 0.;

                return false;
            }
        }  

        // Readout
        m_sBalanceAcq.bCalcRMS(sGrms, sGrms);

        // Compensation gradients
        if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
        {
            m_CompGrad.getGPALoad().bCalcRMS(sGrms, sGrms);
        }

        // Diffusion
        m_sBalanceDiff.bCalcRMS(sGrms, sGrms);
    }

    // Required pause per diffusion direction (=> per effective TR)
    //
    // Check sum of all gradient axes
    lMinPause = DiffMemBalance.lCalcPauseRMS(sGrms) / lTotalScans;

    if(lMinPause < 0)
    {
        // No additional pause required
        lMinPause = 0;
    }

    if(lMinPause > lTRIncrement)
    {
        // Thermal load determines minimum pause duration
        lTRIncrement = lMinPause;

        return true;
    }

    if(lTRIncrement == 0)
    {
        // No additional pause required
        return true;
    }

    // ------------------------------------------------------------------------
    // 5.
    // ------------------------------------------------------------------------
    // So far, we assumed that all diffusion gradient activity takes place
    // on the readout axis for all the time - usually, the diffusion direction
    // will change and the load will be somewhat distributed among the three
    // gradient axes. Let's check whether this helps in order to decrease the
    // pause duration.

    // First we need to calculate the balance for all diffusion directions.
    // While in principle it would also be possible to explicitly consider
    // the actual b-values, this would yield an extremely (!) complicated UI
    // behaviour. Thus, only the dependency on the direction is investigated.

    // Instead of considering the actual gradient amplitudes of each direction,
    // we assign each amplitude to one of a limited number N of predefined
    // amplitudes. Only for those N amplitudes, the actual balance evolution
    // is actually calculated. This helps to speed up the calculation for
    // large diffusion direction sets.

    // Number of predefined amplitudes (sign is included with thermal balancing, so the number increases)
    int    iAmpGrid = m_bThermalBalancing ? 41 : 21;

    // we want to check all 3 axes with Thermal Balancing
    int	 iAxesToCheck = m_bThermalBalancing ? 3 : 1;

    GPABalance::GPABalanceValue sBalanceTemp;
    std::vector<GPABalance>		DiffBalance(iAmpGrid);
    std::vector<GPABalance>     CompGradBalance(iAmpGrid);
    GPABalance   PauseBalance;

    long         lPauseID            = -1;
    long         lSlices             = m_lSlices;
    int          iI                  = 0;
    double       dAmpScale           = 0.;

    std::vector<std::vector<int>  >		viAmpAssign(lTotalScans, std::vector<int>(iAxesToCheck, 0));
    std::vector<std::vector<double> >	vdAmpTable(lTotalScans, std::vector<double>(iAxesToCheck, 0.));
    std::vector<double>					vdAmpGrid(iAmpGrid);

    // Prepare diffusion gradient events (single axis only)
    if(!prepGPALoadDiff(dAmplitude))
    {
        SEQ_TRACE_ERROR.print("preparation of diffusion GPA load failed");

        // Since something failed, return a meaningless amplitude
        dMaxAmplitude = 0.;

        return false;
    }

    // Prepare compensation gradient events (single axis only)
    if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
    {
        double dECCGAmplitude = m_CompGrad.getdMaxAmplitude() * dAmplitude / m_dMaxAmpl;
        if (!m_CompGrad.prepGPALoad(dECCGAmplitude))
        {
            SEQ_TRACE_ERROR.print("ERROR: preparation of ECCG GPA load failed");

            // Since something failed, return a meaningless amplitude
            dMaxAmplitude = 0.;

            return false;
        }
        CompGradMemBalance = m_CompGrad.getGPALoad();
    }

    // Store prepared balance
    DiffMemBalance = m_sBalanceDiff;

    // ------------------------------------------------------------------------
    // Prepare balances for a set of predefined amplitudes
    for(iI = 0; iI < iAmpGrid; ++iI)
    {
        if (m_bThermalBalancing)
        {
            vdAmpGrid[iI] = -1. + 2. * static_cast<double>(iI) / static_cast<double>(iAmpGrid - 1);
        }
        else
        {
            vdAmpGrid[iI] = static_cast<double>(iI) / static_cast<double>(iAmpGrid - 1);
        }

        // Copy reference balance object
        m_sBalanceDiff  = DiffMemBalance;

        if(!scaleGPALoadDiff(vdAmpGrid[iI]))
        {
            SEQ_TRACE_ERROR.print("scaling error");

            // Since something failed, return a meaningless amplitude
            dMaxAmplitude = 0.;

            return false;
        }

        DiffBalance[iI] = m_sBalanceDiff;

        if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
        {
            // Copy reference balance object
            m_CompGrad.setGPALoad(CompGradMemBalance);

            if (!m_CompGrad.scaleGPALoad(vdAmpGrid[iI]))
            {
                SEQ_TRACE_ERROR.print("scaling error");

                // Since something failed, return a meaningless amplitude
                dMaxAmplitude = 0.;

                return false;
            }

            CompGradBalance[iI] = m_CompGrad.getGPALoad();
        } 
    }

    // ------------------------------------------------------------------------
    // Compile amplitude table that covers all scans.
    for (int iAxis = 0; iAxis < iAxesToCheck; iAxis++)
    {
        for (lI = 0; lI < lTotalScans; ++lI)
        {
            // Get actual diffusion scaling for the actual axis, or fo the worst-case axis
            int iSwitchAxis = m_bThermalBalancing ? iAxis : iAmplAvgMax;

            // depending on thermal balancing, use the appropriate direction index
            long lDirectionIndex = m_bThermalBalancing ? m_sDiffusionOrderInfo.getDirectionIndex(lI) : lI;

            switch (iSwitchAxis)
            {
            case 0:
                dAmpScale = m_Didi.getX(lDirectionIndex);
                break;
            case 1:
                dAmpScale = m_Didi.getY(lDirectionIndex);
                break;
            case 2:
                dAmpScale = m_Didi.getZ(lDirectionIndex);
                break;
            default:
                SEQ_TRACE_ERROR.print("unknown gradient axis.");
                // Since something failed, return a meaningless amplitude
                dMaxAmplitude = 0.;

                return true;
            }

            // set the relative amplitude according to the actual b-values for thermal balancing
            if (m_bThermalBalancing && (m_eDiffusionMode != SEQ::DIFFMODE_QSPACE))
                dAmpScale *= sqrt(m_vdBValues[m_sDiffusionOrderInfo.getBValueIndex(lI)] / dMaxBValue);

            // without thermal balancing, take the absolute value
            if (!m_bThermalBalancing)
                vdAmpTable[lI][iAxis] = std::fabs(dAmpScale);
            else
                vdAmpTable[lI][iAxis] = dAmpScale;
        }
    }
 
    // ------------------------------------------------------------------------
    // Since the order of amplitudes might be different for each axis (and we
    // only consider one axis without thermal balancing), we assume a worst case ordering of the
    // amplitude series. This means: sort the amplitudes in increasing order.
    if (!m_bThermalBalancing)
    {
        std::sort(vdAmpTable.begin(), vdAmpTable.end(), [](auto const& elem1, auto const& elem2) {
            return elem1[0] < elem2[0];
        });
    }

    // ------------------------------------------------------------------------
    // For each scan, assign an index of the amplitude grid depending on the 
    // actual amplitude.
    for (int iAxis = 0; iAxis < iAxesToCheck; iAxis++)
    {
        for (lI = 0; lI < lTotalScans; ++lI)
        {
            dAmpScale = vdAmpTable[lI][iAxis];

            int iIndex;

            if (dAmpScale >= 0.)
            {
                if (m_bThermalBalancing)
                    // Round up to nearest bin within boundaries   (considers worst case)
                    iIndex = std::min(iAmpGrid / 2 + static_cast<int>(ceil(dAmpScale * static_cast<double>(iAmpGrid / 2))), iAmpGrid - 1);
                else
                    // we only have positive values without thermal balancing
                    iIndex = static_cast<int>(dAmpScale * static_cast<double>(iAmpGrid - 1) + (1. - 1. / (static_cast<double>(iAmpGrid) - 1)));
            }
            else
            {
                if (m_bThermalBalancing)
                    // Round down to nearest bin within boundaries (considers worst case)
                    iIndex = std::max(iAmpGrid / 2 - static_cast<int>(ceil(-dAmpScale * static_cast<double>(iAmpGrid / 2))), 0);
                else
                {
                    //We should not have a negative value here, something failed. return a meaningless amplitude
                    SEQ_TRACE_ERROR.print("unexpected negative amplitude scale %7f in scan %li.", dAmpScale, lI);
                    dMaxAmplitude = 0.;
                    return true;
                }
            }

            // check and make sure the final value is meaningful
            if (iIndex >= iAmpGrid)
            {
                iIndex = iAmpGrid;
            }
            if (iIndex < 0)
            {
                iIndex = 0;
            }

            viAmpAssign[lI][iAxis] = iIndex;
        }
    }


    // ------------------------------------------------------------------------
    // Perform binary search: identify minimum pause duration
    long lTRIncrementMax = lMinPause;

    for (int iAxis = 0; iAxis < iAxesToCheck; iAxis++)
    {
        long lTRIncrementTemp = lTRIncrement; // start from the same value for all axes
        long lTRIncUpper  = lTRIncrementTemp;
        long lTRIncLower  = lMinPause;
        long lTRIncNew    = 0;

        // Insert dummy pause for first loop iteration
        lPauseID = PauseBalance.lAddEmptyEvent(0);

        while((lTRIncUpper - lTRIncLower) > lTRIncSearchLimit)
        {
            // Prepare new pause
            lTRIncNew = static_cast<long>(lTRIncUpper + lTRIncLower) / 2;
            PauseBalance.bRemoveEvent(lPauseID);
            lPauseID  = PauseBalance.lAddEmptyEvent(lTRIncNew);
            if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
            {
                lDuration = m_sBalanceAcq.lGetDuration() + m_CompGrad.getGPALoad().lGetDuration()
                            + m_sBalanceDiff.lGetDuration() + lTRIncNew;
            }
            else
            {
                lDuration = m_sBalanceAcq.lGetDuration() + m_sBalanceDiff.lGetDuration() + lTRIncNew;
            }

            // Reset values
            GPABalance::ResetBalanceValue(sBalanceIn);
            lResult = 0;
            lI      = 0;

            // Check the whole chain of gradient events!
            while((lI < lTotalScans) && (lResult == 0))
            {
                // Calculate balance change for ONE slice
                GPABalance::ResetBalanceValue(sBalanceTemp);

                // Readout
                if(!m_sBalanceAcq.bCalcBalanceChange(GPABALANCE_X_AXIS, sBalanceTemp, sBalanceTemp))
                {
                    lResult = -1;
                }
                else
                {
                    if (m_CompGrad.isPrepared() && m_CompGrad.getbUseGPABalance())
                    {
                        if (!CompGradBalance[viAmpAssign[lI][iAxis]].bCalcBalanceChange(
                                GPABALANCE_X_AXIS, sBalanceTemp, sBalanceTemp))
                        {
                            lResult = -1;
                        }
                        else
                        {
                            // Pause
                            if (!PauseBalance.bCalcBalanceChange(GPABALANCE_X_AXIS, sBalanceTemp, sBalanceTemp))
                            {
                                lResult = -1;
                            }
                            else
                            {
                                // Diffusion
                                if (!DiffBalance[viAmpAssign[lI][iAxis]].bCalcBalanceChange(
                                        GPABALANCE_X_AXIS, sBalanceTemp, sBalanceTemp))
                                {
                                    lResult = -1;
                                }
                            }
                        }
                    }
                    else
                    {
                        // Pause
                        if (!PauseBalance.bCalcBalanceChange(GPABALANCE_X_AXIS, sBalanceTemp, sBalanceTemp))
                        {
                            lResult = -1;
                        }
                        else
                        {
                            // Diffusion
                            if (!DiffBalance[viAmpAssign[lI][iAxis]].bCalcBalanceChange(
                                    GPABALANCE_X_AXIS, sBalanceTemp, sBalanceTemp))
                            {
                                lResult = -1;
                            }
                        }
                    }
                }

                // Extrapolate to actual number of slices
                if(lResult == 0)
                {
                    if(!m_sBalanceDiff.bCalcBalanceChange(sBalanceIn, sBalanceOut, sBalanceTemp, lDuration, lSlices))
                    {
                        lResult = -1;
                    }
                    else
                    {
                        sBalanceIn = sBalanceOut;
                    }
                }

                ++lI;
            }

            // Check for infinite repeatability
            if(lResult == 0)
            {
                lDuration *= lSlices * lTotalScans;

                if(m_sBalanceDiff.lRepeatabilityCheck(sBalanceOut, lDuration) < lInfRepeatability)
                {
                    lResult = -1;
                }
            }

            // Update search boundaries depending on result
            if(lResult != 0)
            {
                // Further increase necessary
                lTRIncLower = lTRIncNew;
            }
            else
            {
                // Further reduction possible
                lTRIncUpper = lTRIncNew;
            }

            // Remember the longest pause that works
            lTRIncrementTemp = lTRIncUpper;
        }

        if (lTRIncrementTemp > lTRIncrementMax)
            lTRIncrementMax = lTRIncrementTemp;
    }

    // storing the maximum increment in the output variable
    lTRIncrement = lTRIncrementMax;

    return bReturn;
}



// ===========================================================================
/*!
\author   Michael.Zwanger@med.siemens.de

\brief This function returns the sign of the argument

A trivial helper method which returns the sign of the argument.

\return +1 or -1

*/
// ===========================================================================


int SBBDiffusion_Base::sign(double x)
{
    if(x < 0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

bool SBBDiffusion_Base::prep_SBBMultibandRF(MrProt &rMrProt, SeqLim  &rSeqLim, SeqExpo &rSeqExpo, SBBMultibandRF & SBBMultiband, IRF_PULSE * pBaseRF)
{
    // pointer for MrProtFacade for easier protocol queries
    MrProtFacade protFacade(rMrProt);

    KernelCalculationLimits DummyCalcLimits;
    DummyCalcLimits.resetAllLimits();

    if(!SBBMultiband.setPointerToCalculationLimits(&DummyCalcLimits))
        return false;

    if(!SBBMultiband.setMultibandParam(SMSProperties::getMultiBandFactor(rMrProt), SMSProperties::getSliceSeparation(rMrProt)))
        return false;

    const sRFPulseProperties myRFPulseProperties = m_RFPulseLibrary.getPulsePropertiesRefocusing(rMrProt);
    if(!SBBMultiband.setVERSEParam(myRFPulseProperties.bIsVERSE, myRFPulseProperties.fFactorVERSE, myRFPulseProperties.fRelPlateauLengthVERSE))
        return false;

    if(!SBBMultiband.setBandwidthOptimizations(true))
        return false;

    SBBMultiband.setSMSPulsePreDampening(SMSProperties::isSuitableForPulsePredampening());

    if(!SBBMultiband.setRFPulse(rMrProt, rSeqLim, rSeqExpo, pBaseRF))
        return false;

    // set gradient performance for * modi
    SBBMultiband.setGSWDGradientPerformance(rMrProt, rSeqLim);

    // the flat top time of the slice selection gradient should be a multiple of 2 * RGT
    SBBMultiband.setSliceSelectionGradientDurationOnDoubleRasterTime(true);

    SBBMultiband.setIgnoreRFLeadTime(true);

    // set dynamic adjustment data
    SBBMultiband.seteSliceAdjOptimizationMode(SLICEADJ::HOLD_OPTIMIZATION);
    SBBMultiband.setsSliceAdjParametersRequestedBySequence(getsSliceAdjParametersRequestedBySequence());

    if(protFacade.isSliceAdj())
    {
        // SliceAdj requires that all SBBs are in single band mode
        setRunMode(SINGLE_BAND);
        SBBMultiband.setRunMode(SINGLE_BAND);

        // the refocusing pulse SBBs are in HOLD_OPTIMIZATION mode such that we need to 
        // pass the cuboids from this SBB to the refocusing SBBs
        std::vector<SLICEADJ::sCuboidGeometry> vsSliceAdjCuboids;

        if(!getSliceAdjCuboids(vsSliceAdjCuboids))
            return false;

        if(!SBBMultiband.prep(rMrProt, rSeqLim, rSeqExpo, vsSliceAdjCuboids))
            return false;
    }
    else
    {
        if(!SBBMultiband.prep(rMrProt, rSeqLim, rSeqExpo))
            return false;
    }

    if(!SBBMultiband.checkGradients(rMrProt, rSeqLim))
    {
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ALWAYS.print("ERROR: m_SBBMultibandRFRefoc.checkGradients failed.");
        }
        return false;
    }

    //ensure RF pulse is on double the gradient raster time
    SBBMultiband.adjustDurationToGradRasterTime();


    //check for pulse clipping 
    if(!rSeqLim.isContextPrepForBinarySearch())
    {
        if(protFacade.isSliceAcceleration())
        {
            SBBMultiband.setRunMode(MULTI_BAND);
            if(SBBMultiband.isRFClipped())
            {
                SEQ_TRACE_WARN.print("INFO: m_SBBMultibandRFRefoc is clipping");

                // write clipped flip angle to protocol
                if(SBBMultiband.getRFPulsePointer()->hasIdent())
                {
                    // get a protocol entry for this RF pulse ...
                    MrProtocolData::MrRFPulseData::Pointer  sProtRfPulse;

                    if(rMrProt.txSpec().rfPulse().getbyName(SBBMultiband.getRFPulsePointer()->getIdent(), sProtRfPulse))
                    {
                        // ... and write new flip angle into 'userInfo' field
                        sProtRfPulse->setdUserInfo(SBBMultiband.getRFPulsePointer()->getActualFlipAngle());
                    }
                }
            }
            SBBMultiband.setRunMode(m_eSliceAccelRFRunMode);
        }
    }

    return true;
}

bool SEQ_NAMESPACE::SBBDiffusion_Base::setADCforDiffusionMDHentries(sREADOUT* pADC)
{
    
    if(!pADC)
    {
        SEQ_TRACE_ERROR.print("ERROR: received NULL pointer");
        setNLSStatus(MRI_SBB_SBB_ERROR);
        return false;
    }

    m_pADC = pADC;
    return true;
}

long SBBDiffusion_Base::getIVIMIncrement() const
{
    return IVIM_B_VALUE_INCREMENT;
}

long SBBDiffusion_Base::getMaxBValueSmallIVIMIncrement() const
{
    return IVIM_B_VALUE_MAX_SMALL_INC;
}

