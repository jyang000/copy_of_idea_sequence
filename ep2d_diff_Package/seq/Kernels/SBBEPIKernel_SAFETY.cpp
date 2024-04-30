//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2021  All Rights Reserved.
//    -----------------------------------------------------------------------------

#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrRXSpec.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrSysSpec.h"
#include "MrProtSrv/Domain/CoreNative/MrPat.h"

#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"

#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"     // for fGSLGetPEDeltaMoment etc.
#include "MrImagingFW/libSeqUTIF/libsequt.h"                // for mSEQTest
#include "MrMeasSrv/SeqIF/libRT/sFREQ_PHASE.h" // for sFREQ_PHASE
#include "MrMeasSrv/SeqIF/libRT/RTController.h" // for getAbsTimeOfEventBlockMSec
#include "MrImagingFW/libSeqUtilFW/KernelCalculationLimits.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrImagingFW/libSeqSysProp/SysProperties.h"
#include "MrImaging/libSeqPTX/C2DGradientsBlippedEPI.h"
#include "MrImaging/libSeqPTX/SBB2DPtx.h"     
#include "MrVista/Ice/IceCommonFunctorsAndAlgos/IceSFC/SFCDefinitions.h"           // SeqToIceSliceAdjustData
#include "MrImaging/libSeqUtil/SliceAccelerationUtils.h"

//---------------------------------------------------------------------------
// Debug flags
//---------------------------------------------------------------------------
//#define DEBUG_calcEffEchoSpacingAndBWPerPixelPE


//---------------------------------------------------------------------------
// Requirement keys
//---------------------------------------------------------------------------
//
// EGA Requirement Key: As shown on the following lines:
//
//   Abbrev.   Translation                                        Relevant for
//   -------   -----------                                        ------------
//   EGA-All   One, some or all of the following keys:            Any EGA requirement
//   EGA-01    {:IMPLEMENT:000_EGA_BildOri_SW_SequenzROVz::}      GR/GP   polarity
//   EGA-02    {:IMPLEMENT:000_EGA_BildPos_SW_SequenzSSelVz::}    GS      polarity
//   EGA-03    {:IMPLEMENT:000_EGA_BildMass_SW_SequenzROPC::}     GR/GP   amplitude
//   EGA-04    {:IMPLEMENT:000_EGA_BildPos_SW_SequenzSSel::}      GS      amplitude
//   EGA-05    {:IMPLEMENT:000_EGA_BildPos_SW_NCOFrequenzSSel::}  SRF     frequency
//   EGA-06    {:IMPLEMENT:000_EGA_BildPos_SW_NCOFrequenzRO::}    Readout frequency
//   EGA-07    {:IMPLEMENT:000_EGA_BildOri_SW_OrientierungTest::} Image orientation
//
//
#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif

using namespace SEQ_NAMESPACE;
using namespace std;


SeqBuildBlockEPIKernel::SeqBuildBlockEPIKernel(SBBList* pSBBList)
    : SeqBuildBlockEPIReadOut(pSBBList)
{
    setUseEPILikePhaseCorrection();
    SeqBuildBlockEPIReadOut::setPerformPhaseCorrection(true);//Call explicitly to avoid PCLint Warning 1506

    m_aGradMoments[0].setIdent("MO section 0");
    m_aGradMoments[1].setIdent("MO section 1");
    m_aGradMoments[2].setIdent("MO section 2");
    m_aGradMoments[3].setIdent("MO section 3");

    m_pSBBExcite = nullptr;

    // Preset: no history
    GPABalance::ResetBalanceValue(m_sBalanceIn);
    GPABalance::ResetRMSValue(m_sBalanceRmsIn);

    for (int i = 0; i < 64; i++)
    {
        m_tIdent[i] = '\0';
    }

    // ----------------------------------------------------
    // Configure SBB-specific dynamic adjustment properties
    // ----------------------------------------------------
    // Note: SeqBuildBlockEPIKernel applies a number of SBB's, each with its own RTEB
    //       - Excitation
    //       - Phase correction
    //       - PlugIn (e.g. refocusing)
    //       - Readout
    //       An update of dynamic adjustments shall take place only before the excitation
    //       module - control parameters have to stay constant afterwards. This has to be
    //       considered when preparing or running the SBB's internally.
    seteSliceAdjOptimizationMode(SLICEADJ::LOCAL_OPTIMIZATION);  // Selective excitation     => local optimization
    setsSliceAdjParametersRequestedBySBB(SLICEADJ::ADJALL);  // Applies excitation pulse => consider all parameters


    //-------------------------------------------------------------------------------------
    // change default gradient performance of excitation SBB
    //-------------------------------------------------------------------------------------
    const double dRiseTimeMin = SysProperties::getGradMinRiseTimeAbsolute();
    double adRiseTimes[3] = { dRiseTimeMin, dRiseTimeMin, dRiseTimeMin };

    const double dGradMaxMagn = SysProperties::getGradMaxAmpl(SEQ::GRAD_FAST);
    double adMagnitudes[3] = { dGradMaxMagn, dGradMaxMagn, dGradMaxMagn };

    m_SBBExcitation.setMinRiseTimes(adRiseTimes, SBBBinomialPulses_GRAD_PERF_BINOMIAL);
    m_SBBExcitation.setMinRiseTimes(adRiseTimes, SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
    m_SBBExcitation.setMaxMagnitudes(adMagnitudes, SBBBinomialPulses_GRAD_PERF_BINOMIAL);
    m_SBBExcitation.setMaxMagnitudes(adMagnitudes, SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);

    m_SBBExcitation.setIdent("SBBExcitation");

    //-------------------------------------------------------------------------------------
    // protocol independent configuration of binomial pulse excitation
    //-------------------------------------------------------------------------------------       
    m_SBBExcitationBinomial.setMinRiseTimes(adRiseTimes, SBBBinomialPulses_GRAD_PERF_BINOMIAL);
    m_SBBExcitationBinomial.setMinRiseTimes(adRiseTimes, SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
    m_SBBExcitationBinomial.setMaxMagnitudes(adMagnitudes, SBBBinomialPulses_GRAD_PERF_BINOMIAL);
    m_SBBExcitationBinomial.setMaxMagnitudes(adMagnitudes, SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);

    m_SBBExcitationBinomial.setBandwidthTimeProduct(5.2);
    m_SBBExcitationBinomial.setUsePossibleBandwidthTimeProduct(true, 1.6);     // Allow BWT reduction down to the given limit
    m_SBBExcitationBinomial.setTBWAdaption(true);                              // Consider internal filtering of Sinc-pulses => slice gradient correction
}

bool SeqBuildBlockEPIKernel::calcSingleSliceAdjSBBCuboid(
    MrProt                  &rMrProt,       //< IMP: points to the protocol structure.
    SeqLim                  &rSeqLim,       //< points to the sequence limits structure.
    SeqExpo                 &rSeqExpo,      //< points to the sequence exports structure.
    const sSLICE_POS*        pSLC,          //< IMP: points to the slice position structure of the current slice
    SLICEADJ::sCuboidGeometry &sSliceAdjCuboid  //< EXP: cuboid geometry (single)
)
{
    return SeqBuildBlock::calcSliceAdjSliceCuboid(rMrProt, rSeqLim, rSeqExpo, pSLC, sSliceAdjCuboid);
}

bool SeqBuildBlockEPIKernel::calcSliceAdjSBBCuboids(
    MrProt                               & rMrProt,         //< IMP: points to the protocol structure.
    SeqLim                               & rSeqLim,         //< points to the sequence limits structure.
    SeqExpo                              & rSeqExpo,        //< points to the sequence exports structure.
    std::vector<SLICEADJ::sCuboidGeometry> &vsSliceAdjCuboids   //< EXP: cuboid geometries (multiple elements)
)
{
    vsSliceAdjCuboids.clear();
    return SeqBuildBlock::calcSliceAdjSliceSeriesCuboids(rMrProt, rSeqLim, rSeqExpo, vsSliceAdjCuboids);
}


//===============================================================================
//
// Function:    setUseSyncBits()
//
// Description: Tells the kernel whether sync-bits should be sent before
//              the excitation SBB or not. If so, an osc bit and an external
//              trigger bit are prepared. Member m_lMaxSyncBitDuration is set.
//                If the function executed with success, true is returned.
//                The execution of the osc-bit can be disabled/enabled during
//              run-time of the sequence using the method setDoNotSendOscBit.
//                The execution of the external trigger bit can be disabled/enabled
//              during run-time of the sequence using the method
//              setDoNotSendExtTrigger.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::setUseSyncBits(bool bValue, long lOscChannel, long lOscDuration, long lOscStartTime, long lExtTrigDuration, long lExtTrigStartTime)
{
    m_bUseSyncBits = bValue;

    if (m_bUseSyncBits)
    {
        if (!m_OscBit.prep(lOscChannel, lOscDuration))
        {
            setNLSStatus(m_OscBit.getNLSStatus(), "SeqBuildBlockEPIKernel::setUseSyncBits", "m_OscBit.prep failed");
            return false;
        }
        m_OscBit.setStartTime(lOscStartTime);

        if (!m_ExtTrig.prep(0, lExtTrigDuration))
        {
            setNLSStatus(m_ExtTrig.getNLSStatus(), "SeqBuildBlockEPIKernel::setUseSyncBits", "m_ExtTrig.prep failed");
            return false;
        }
        m_ExtTrig.setStartTime(lExtTrigStartTime);

        m_lMaxSyncBitDuration = std::max(m_OscBit.getStartTime() + m_OscBit.getDuration(),
            m_ExtTrig.getStartTime() + m_ExtTrig.getDuration()
        );
    }
    else
    {
        m_lMaxSyncBitDuration = 0;
    }

    resetPrepared();
    return true;
}


//===============================================================================
//
// Function:    setUseEchoShifting()
//
// Description: Activates or deactivates echo-shifting. If echo shifting
//              should be activated
//                It is important for the kernel to know the number of
//              (lines res. partitions), counters per segment and the
//              (line res. partition) counter within the center
//                segment with which the k-space center is measured.
//                For deactivation only false needs to be passed to this
//              function.
//
//===============================================================================

void SeqBuildBlockEPIKernel::setUseEchoShifting(bool bValue, long lCountersPerSegment, long lCounterInSegmentWithEcho)
{
    SeqBuildBlockEPIReadOut::setUseEchoShifting(bValue, lCountersPerSegment);

    if (bValue) m_lCounterInSegmentWithEcho = lCounterInSegmentWithEcho;
    else        m_lCounterInSegmentWithEcho = 0;

    resetPrepared();
}


//===============================================================================
//
// Function:    setAdditionalPhase()
//
// Description: Can be used to specify an additional phase to the NCO
//              during the measurement.
//                Note that this will overwrite any additional phase
//              specified for the excitation SBB within the sequence before
//              calling this method!
//
//===============================================================================

void SeqBuildBlockEPIKernel::setAdditionalPhase(double dValue)
{
    //---------------------------------------------------------------------------
    // set additional phase for EPI read out
    //---------------------------------------------------------------------------
    m_dAdditionalPhase = dValue;

    //---------------------------------------------------------------------------
    // set additional phase for excitation SBB
    //---------------------------------------------------------------------------
    if (m_pSBBExcite) m_pSBBExcite->setAdditionalPhase(m_dAdditionalPhase);
}


//===============================================================================
//
// Function:    updateGradientPerformance()
//
// Description: Internal function to update the performance data of all
//              gradients involved within the kernel.
//
//===============================================================================

void SeqBuildBlockEPIKernel::updateGradientPerformance(SEQ::Gradients eGradMode)
{
    //---------------------------------------------------------------------------
    // update performance for EPI read out
    //---------------------------------------------------------------------------
    SeqBuildBlockEPIReadOut::updateGradientPerformance(eGradMode);

    //---------------------------------------------------------------------------
    // update performance for gradient moment sections
    //---------------------------------------------------------------------------
    for (auto& gradMoment : m_aGradMoments)
    {
        gradMoment.setMinRiseTime(getMinRiseTime(eGradMode, SBBEPIKernel_GRAD_PERF_GRAD_MOMENTS));
        gradMoment.setMaxMagnitude(getMaxMagnitude(eGradMode, SBBEPIKernel_GRAD_PERF_GRAD_MOMENTS));
    }
}


//===============================================================================
//
// Function:    calculateTimingOfGradMoments()
//
// Description: Internal function to calculate the timing of the gradients
//              realizing gradient moments between the SBBs used by the kernel.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::calculateTimingOfGradMoments(MrProt &rMrProt, SeqLim &rSeqLim)
{
    double adMaxGradMoment[4];

    double dMr = 0., dMp = 0., dMs = 0.;

    //---------------------------------------------------------------------------
    // m_pSBBExcite must have been set
    //---------------------------------------------------------------------------
    if (m_pSBBExcite == nullptr)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "missing pointer to excitation SBB");
        return false;
    }

    //---------------------------------------------------------------------------
    // This is IMPORTANT!
    // Otherwise getPEPrePhasingMoment and get3DPrePhasingMoment may return zero
    // when called next time, if EPIReadOut is in phase correction mode.
    //---------------------------------------------------------------------------
    setNextExecutionIsImaging(false); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // calculate basic frequency and phase encoding moments
    //---------------------------------------------------------------------------
    double dMax3DPrephasingMoment = 0.0;
    double dMaxROPrephasingMoment = 0.0;
    double dMaxPEPrephasingMoment = 0.0;

    const char* nucleus = rMrProt.getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str();

    // 3D imaging
    if (rMrProt.kSpace().partitions() > 1)
    {
        const double dThick      = rMrProt.sliceSeries().front().getdThickness();
        const double d3dOSFactor = rMrProt.kSpace().sliceOversampling();

        dMax3DPrephasingMoment = fabs(fGSLGet3DDeltaMoment(dThick, d3dOSFactor, nucleus) * (double)rMrProt.kSpace().partitions() / 2.0 * m_pCalcLimits->getFactorForPixelSize3D(rMrProt));
    }
    // Multiband imaging
    else if (m_bIsSliceAcceleration && m_lFOVShiftFactor > 1)
    {
        dMax3DPrephasingMoment = fabs( SliceAccelerationUtils::calcSliceAccelerationDeltaMoment(rMrProt) * static_cast<double>( rMrProt.getsSliceAcceleration().getlFOVShiftFactor() - 1 ) / 2.0 );
    }
    // 2D imaging
    else
    {
        dMax3DPrephasingMoment = 0.0;
    }

    const double dFOVPH         = rMrProt.sliceSeries().front().getdPhaseFOV();
    const double dPhaseOSFactor = rMrProt.phaseOversampling();

    dMaxROPrephasingMoment = fabs(getROPrePhasingMoment(rMrProt, rSeqLim) * m_pCalcLimits->getFactorForPixelSizeRO(rMrProt));
    dMaxPEPrephasingMoment = fabs(fGSLGetPEDeltaMoment(dFOVPH, dPhaseOSFactor, nucleus) * (double)rMrProt.kSpace().phaseEncodingLines() / 2.0 * m_pCalcLimits->getFactorForPixelSizePE(rMrProt));

    // PEPrephaser is enlarged by the maximum possible m_lPlaceCurrent
    long lPlaceCurrent = 0;
    // choose starting point as MIN of vector: this does not ensure that overflow of phaseEncodingLines() will occur!
    lPlaceCurrent = m_lPlaceCurrentMin;   // (0) default 

    //---------------------------------------------------------------------------
    // calculate PLACE settings: PEPrephaser is enlarged by the maximum possible m_lPlaceCurrent
    //---------------------------------------------------------------------------
    if ((lPlaceCurrent < 0) || (lPlaceCurrent > rMrProt.kSpace().phaseEncodingLines()))
    {
        dMaxPEPrephasingMoment = fabs(fGSLGetPEDeltaMoment(dFOVPH, dPhaseOSFactor, nucleus)
            * (static_cast<double>(rMrProt.kSpace().phaseEncodingLines()) / 2.0 - static_cast<double>(lPlaceCurrent))
            * m_pCalcLimits->getFactorForPixelSizePE(rMrProt));
    }

    //---------------------------------------------------------------------------
    // calculate basic moments for REphasing the excitation SBB
    //---------------------------------------------------------------------------
    double dMaxGSRephasingMoment = 0.0;
    double dMaxGRRephasingMoment = 0.0;
    double dMaxGPRephasingMoment = 0.0;

    if (m_pSBBExcite->getGSData().getbAxisActive())
    {
        dMaxGSRephasingMoment = fabs(m_pSBBExcite->getGSData().getdRequiredRefocusingMoment());
    }

    if (m_pSBBExcite->getGPData().getbAxisActive())
    {
        dMaxGPRephasingMoment = fabs(m_pSBBExcite->getGPData().getdRequiredRefocusingMoment());
    }

    if (m_pSBBExcite->getGRData().getbAxisActive())
    {
        dMaxGRRephasingMoment = fabs(m_pSBBExcite->getGRData().getdRequiredRefocusingMoment());
    }

    // Apply kernel calculation limits only to GS refocusing moment,
    // because on other axes there is no useful kernel calculation limit
    // available:
    //
    dMaxGSRephasingMoment *= m_pCalcLimits->getFactorForSliceThickness(rMrProt);


    //---------------------------------------------------------------------------
    // calculate basic moments for PREphasing the excitation SBB
    //---------------------------------------------------------------------------
    double dMaxGSPrephasingMoment = 0.0;
    double dMaxGRPrephasingMoment = 0.0;
    double dMaxGPPrephasingMoment = 0.0;

    if (m_pSBBExcite->getGSData().getbAxisActive())
    {
        dMaxGSPrephasingMoment = fabs(m_pSBBExcite->getGSData().getdMomentBeforeFirstRFCenter());
    }

    if (m_pSBBExcite->getGPData().getbAxisActive())
    {
        dMaxGPPrephasingMoment = fabs(m_pSBBExcite->getGPData().getdMomentBeforeFirstRFCenter());
    }

    if (m_pSBBExcite->getGRData().getbAxisActive())
    {
        dMaxGRPrephasingMoment = fabs(m_pSBBExcite->getGRData().getdMomentBeforeFirstRFCenter());
    }

    // Apply kernel calculation limits only to GS Prefocusing moment,
    // because on other axes there is no useful kernel calculation limit
    // available:
    //
    dMaxGSPrephasingMoment *= m_pCalcLimits->getFactorForSliceThickness(rMrProt);

    //---------------------------------------------------------------------------
    // calculate moment for section 0
    //---------------------------------------------------------------------------
    dMr = dMp = dMs = 0.0;

    if (m_bInternalPhaseCorrection || !m_bPrephaseROAfterRTEBPlugIn)
    {
        dMr = dMaxROPrephasingMoment;
    }

    if (!m_bInternalPhaseCorrection && !m_bPrephaseBlipsAfterRTEBPlugIn)
    {
        dMp = dMaxPEPrephasingMoment;
        dMs = dMax3DPrephasingMoment;
    }

    dMs += dMaxGSRephasingMoment;
    dMp += dMaxGPRephasingMoment;
    dMr += dMaxGRRephasingMoment;

    if (m_bRotationProofGradientMoments)
    {
        adMaxGradMoment[0] = sqrt(dMr*dMr + dMp * dMp + dMs * dMs);
    }
    else
    {
        adMaxGradMoment[0] = std::max(std::max(dMr, dMp), dMs);
    }

    //---------------------------------------------------------------------------
    // calculate moment for section 1
    //---------------------------------------------------------------------------
    dMr = dMp = dMs = 0.0;

    if (m_bPrephaseROAfterRTEBPlugIn && m_bInternalPhaseCorrection)
    {
        dMr = dMaxROPrephasingMoment;
    }

    if (!m_bPrephaseBlipsAfterRTEBPlugIn && m_bInternalPhaseCorrection)
    {
        dMp = dMaxPEPrephasingMoment;
        dMs = dMax3DPrephasingMoment;
    }

    if (m_bRotationProofGradientMoments)
    {
        adMaxGradMoment[1] = sqrt(dMr*dMr + dMp * dMp + dMs * dMs);
    }
    else
    {
        adMaxGradMoment[1] = std::max(std::max(dMr, dMp), dMs);
    }

    //---------------------------------------------------------------------------
    // calculate moment for section 2
    //---------------------------------------------------------------------------
    dMr = dMp = dMs = 0.0;

    if (m_bPrephaseROAfterRTEBPlugIn)
    {
        dMr = dMaxROPrephasingMoment;
    }

    if (m_bPrephaseBlipsAfterRTEBPlugIn)
    {
        dMp = dMaxPEPrephasingMoment;
        dMs = dMax3DPrephasingMoment;
    }

    if (m_bRotationProofGradientMoments)
    {
        adMaxGradMoment[2] = sqrt(dMr*dMr + dMp * dMp + dMs * dMs);
    }
    else
    {
        adMaxGradMoment[2] = std::max(std::max(dMr, dMp), dMs);
    }

    //---------------------------------------------------------------------------
    // calculate moment for section 3
    //---------------------------------------------------------------------------
    dMr = dMp = dMs = 0.0;

    if (m_bRewindRO)
    {
        dMr = dMaxROPrephasingMoment;
    }
    else if ( m_bSpoilRO )
    {
        dMr = fabs( m_dSpoilROFactor ) * dMaxROPrephasingMoment;
    }

    if (m_bRewindBlips)
    {
        dMp = dMaxPEPrephasingMoment;
        dMs = dMax3DPrephasingMoment;
    }

    if (m_bPrefaceExcitationSBB)
    {
        dMs += dMaxGSPrephasingMoment;
        dMp += dMaxGPPrephasingMoment;
        dMr += dMaxGRPrephasingMoment;
    }

    if (m_bRotationProofGradientMoments)
    {
        adMaxGradMoment[3] = sqrt(dMr*dMr + dMp * dMp + dMs * dMs);
    }
    else
    {
        adMaxGradMoment[3] = std::max(std::max(dMr, dMp), dMs);
    }

    //---------------------------------------------------------------------------
    // calculate timing for all sections
    //---------------------------------------------------------------------------
    for (long lI = 0; lI < 4; lI++)
    {
        if (adMaxGradMoment[lI] > 0.001) // adMaxGradMoment[lI] is always >= 0.0
        {
            if (!m_aGradMoments[lI].prepSymmetricTOTShortestTime(adMaxGradMoment[lI]))
            {
                // unexpected error
                // => trace also if pSeqLim->isContextPrepForBinarySearch()
                setNLSStatus(m_aGradMoments[lI].getNLSStatus(), __FUNCTION__, "m_GradMoments[lI].prepSymmetricTOTShortestTime failed");
                return false;
            }
        }
        else
        {
            m_aGradMoments[lI].set(0, 0, 0, 0.0);
            if (!m_aGradMoments[lI].prep())
            {
                setNLSStatus(m_aGradMoments[lI].getNLSStatus(), __FUNCTION__, "m_GradMoments[lI].prep failed");
                return false;
            }
        }
    }

    //---------------------------------------------------------------------------
    // finished
    //---------------------------------------------------------------------------
    return true;
}


//===============================================================================
//
// Function:    prep()
//
// Description: Prepares the kernel.
//                At least an excitation SBB configuration function must
//              have been registered.
//                All other configuration steps have to be performed.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::prepSBB(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
{

#ifdef ZOOM_EXTENDED
    //---------------------------------------------------------------------------
    // Check whether to use debug settings
    //---------------------------------------------------------------------------
    m_bUseDebugSettings = SysProperties::ReadSeqSettingGeneral("EPI_ZOOMIT/USE_ZOOMIT_DEBUG_SETTINGS", false, true);
    //SEQ_TRACE_ALWAYS.print("DebugXML: m_bUseDebugSettings = %d", m_bUseDebugSettings);
#endif // ZOOM_EXTENDED

    //---------------------------------------------------------------------------
    // if we continue now, we are no longer prepared
    //---------------------------------------------------------------------------
    resetPrepared();

    //---------------------------------------------------------------------------
    // check whether simultaneous multislice imaging is activated,
    // will be needed in getRFInfoPerMeasurement, calculateTimingOfGradMoments etc.
    //---------------------------------------------------------------------------
    MrProtFacade protFacade(rMrProt);
    m_bIsSliceAcceleration = protFacade.isSliceAcceleration();

    // ---------------------------------------------------------------------------
    // set GSWD gradient performance for excitation SBB
    // ---------------------------------------------------------------------------
    m_SBBExcitationBinomial.setGSWDGradientPerformance(rMrProt, rSeqLim);
    m_SBBExcitation.setGSWDGradientPerformance(rMrProt, rSeqLim);



    // dummy calculation limits for m_SBBExcitation
    KernelCalculationLimits DummyCalcLimitsEx;
    DummyCalcLimitsEx.resetAllLimits();
    m_SBBExcitationBinomial.setPointerToCalculationLimits(&DummyCalcLimitsEx);
    m_SBBExcitation.setPointerToCalculationLimits(&DummyCalcLimitsEx);


    // ---------------------------------------------------------------------------
    // set and prepare B1 control loop in the excitation SBBs,
    // this activates or deactivates whether the SBB marks an excitation pulse
    // for measurement of the reflected power of the BC 
    // ---------------------------------------------------------------------------
    m_SBBExcitationBinomial.setB1ControlLoopActive(rMrProt.getsTXSPEC().getB1CorrectionParameters().getbActive());
    m_SBBExcitation.setB1ControlLoopActive(rMrProt.getsTXSPEC().getB1CorrectionParameters().getbActive());

    //---------------------------------------------------------------------------
    // prepare excitation SBB and set pointer to SBB
    //---------------------------------------------------------------------------
    if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
    {
#ifdef ZOOM_2DRF
        if (!prepZOOMitExcitation(rMrProt, rSeqLim))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
                setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "prepZOOMitExciation() failed");
            return false;
        }
#endif //ZOOM_2DRF
    }
    else
    {
        // Enlarge the max amplitude of Gradient Slice of ep2d_diff to allow thinner slice thickness (3mm) for low-field
        if (SysProperties::isLowField())
        {
            increaseSliceSelectionGradientMaxAmplitude(rMrProt);
        }

        if (!initExcitation(rMrProt, rSeqLim, rSeqExpo))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
                setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "initExcitation() failed");
            return false;
        }

        // for Slice Adjust the runmode always needs to be set to SINGLE_BAND, otherwise the 
        // cuboid calculation will fail
        if (protFacade.isSliceAdj())
        {
            setRunMode(SINGLE_BAND);

            SBBMultibandRF* pSBBMultibandRF = dynamic_cast<SBBMultibandRF*>(m_pSBBExcite);
            if (pSBBMultibandRF != nullptr)
                pSBBMultibandRF->setRunMode(m_eRunMode);
        }
    }

    //---------------------------------------------------------------------------
    // reset local fill times
    //---------------------------------------------------------------------------
    m_lLocalTEFill = 0; // fill time in EPIReadout::run before echo train
    m_lLocalTRFill = 0; // fill time in EPIReadout::run after  echo train
    m_lTRFill = 0;

    // ---------------------------------------------------------------------------
    // update gradient performance
    // ---------------------------------------------------------------------------
    if (didGradPerfChange(rMrProt.gradSpec().mode())) updateGradientPerformance(rMrProt.gradSpec().mode());

    //---------------------------------------------------------------------------
    // configure and prepare excitation SBB
    //---------------------------------------------------------------------------
    m_pSBBExcite->setUseOwnEventBlock(false);

    if (!m_pSBBExcite->setPointerToCalculationLimits(m_pCalcLimits))
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(m_pSBBExcite->getNLSStatus(), __FUNCTION__, "m_pSBBExcite->setPointerToCalculationLimits failed.");
        return false;
    }

    // updateGradientPerformance() in GSWD mode for excitation SBB 
    m_pSBBExcite->setGSWDGradientPerformance(rMrProt, rSeqLim);

    // m_pSBBExcite is part of GREKernel => use identical dynamic adjustment parameter settings
    m_pSBBExcite->seteSliceAdjOptimizationMode(SLICEADJ::HOLD_OPTIMIZATION);
    m_pSBBExcite->setsSliceAdjParametersRequestedBySequence(getsSliceAdjParametersRelevantForUpdate());

    // Get the cuboids of the kernel, because all included SBB are in HOLD_OPTIMIZATION mode and thus
    // do not calculate their cuboids on their own. They require the cuboids as parameter during prep.
    std::vector<SLICEADJ::sCuboidGeometry> vsSliceAdjCuboidsFromKernel;

    if (protFacade.isSliceAdj())
    {
        if (!calcSliceAdjSBBCuboids(rMrProt, rSeqLim, rSeqExpo, vsSliceAdjCuboidsFromKernel))
        {
            // Preparation of excitation module preparation might fail in binary search
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                setNLSStatus(m_pSBBExcite->getNLSStatus(), __FUNCTION__, "calcSliceAdjSBBCuboids failed.");
            }
            return false;
        }

        m_sRFInfoStorage.setCuboids(vsSliceAdjCuboidsFromKernel);

        if (!m_pSBBExcite->prep(rMrProt, rSeqLim, rSeqExpo, vsSliceAdjCuboidsFromKernel))
        {
            // Preparation of excitation module preparation might fail in binary search
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                setNLSStatus(m_pSBBExcite->getNLSStatus(), __FUNCTION__, "m_pSBBExcite->prep failed.");
            }
            return false;
        }
    }
    else
    {
        if (!m_pSBBExcite->prep(rMrProt, rSeqLim, rSeqExpo))
        {
            // Preparation of excitation module preparation might fail in binary search
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                setNLSStatus(m_pSBBExcite->getNLSStatus(), __FUNCTION__, "m_pSBBExcite->prep failed.");
            }
            return false;
        }
    }

#ifdef ZOOM_2DRF
#ifdef ZOOM_EXTENDED
    if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
    {
        //---------------------------------------------------------------------------
        // After successfull preparation of the excitation trajectory,
        // the optimal rotation angle can be calculated
        //---------------------------------------------------------------------------
        if (!updateZOOMitRotation(rMrProt, rSeqLim, rSeqExpo))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "updateZOOMitRotation() failed");
            }
            return false;
        }
    }
#endif // ZOOM_EXTENDED
#endif // ZOOM_2DRF

    //---------------------------------------------------------------------------
    // we do not support TGSE-like phase correction
    //
    // NOTE this:
    // Although we check it here the user of the kernel can still trick us by
    // calling the according set-method of SeqBuildBlockEPIReadOut which does not
    // reset the prepare-status after this prep. We can get problems during run.
    // That's the user's own problem ;-).
    //---------------------------------------------------------------------------
    if (m_bTGSELikePhaseCorr)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_bTGSELikePhaseCorr currently is not supported");
        return false;
    }

    //---------------------------------------------------------------------------
    // The whole moment calculation stuff in RO direction depends on the fact
    // that previously 3 echos were acquired for phase correction.
    // Now the echo-train-length for phase correction is variable, but it must
    // still be an odd number.
    // More then 3 phase correction echos is only supported for internal phase
    // correction.
    //
    // NOTE this:
    // Although we check it here the user of the Kernel can still trick us by
    // calling the according set-method of SeqBuildBlockEPIReadOut which does not
    // reset the prepare-status after this prep. We can get problems during run.
    // That's the user's own problem ;-).
    //---------------------------------------------------------------------------
    if (getEchoTrainLengthPhaseCorrScan() % 2 == 0)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "getEchoTrainLengthPhaseCorrScan() must be an odd number currently");
        return false;
    }

    if (!m_bInternalPhaseCorrection && getEchoTrainLengthPhaseCorrScan() != 3)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_bInternalPhaseCorrection false and getEchoTrainLengthPhaseCorrScan()!=3 currently not allowed");
        return false;
    }


    //---------------------------------------------------------------------------
    // calculate timing of EPI read out
    //---------------------------------------------------------------------------
    if (!SeqBuildBlockEPIReadOut::calculateTiming(rMrProt, rSeqLim))
    {
        // timing calculation might fail in binary search
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "SeqBuildBlockEPIReadOut::calculateTiming failed!");
        return false;
    }

    //---------------------------------------------------------------------------
    // prepare EPI read out
    //---------------------------------------------------------------------------
    // Important: We have to call ::prepSBB() here and not ::prep()!
    // (SBBEPIKernel takes care of preparations related to dynamic adjustments.)
    if (!SeqBuildBlockEPIReadOut::prepSBB(rMrProt, rSeqLim, rSeqExpo))
    {
        // unexpected error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        SEQ_TRACE_ERROR.print("SeqBuildBlockEPIReadOut::prep failed!");
            return false;
    }

    //---------------------------------------------------------------------------
    // calculate timing of gradient moment sections
    //---------------------------------------------------------------------------
    if (!calculateTimingOfGradMoments(rMrProt,rSeqLim))
    {
        // unexpected error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        SEQ_TRACE_ERROR.print("calculateTimingOfGradMoments failed!");
        return false;
    }

    //---------------------------------------------------------------------------
    // Prepares the additional Gradients,e.g. flow compensated gradients
    //---------------------------------------------------------------------------
    if (!prepAdditionalGradients(rMrProt, rSeqLim, rSeqExpo))
    {
        // unexpected error
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "prepAdditionalGradients failed!");

        return false;
    }

    //---------------------------------------------------------------------------
    // calculate TE contribution before and after RTEBPlugIn
    //---------------------------------------------------------------------------
    m_lTEContributionBeforeRTEBPlugIn = 0;
    m_lTEContributionAfterRTEBPlugIn  = 0;
    long lEchoShiftingDelay              = 0;
    long lEchoShiftingFillEnd            = 0;

    if (!getEchoShiftingData (m_lCounterInSegmentWithEcho, &lEchoShiftingDelay, &lEchoShiftingFillEnd))
    {
        // unexpected error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        SEQ_TRACE_ERROR.print("getEchoShiftingData failed!");
        return false;
    }

    //---------------------------------------------------------------------------
    // Check consistency of numbers of contrasts before/after PlugIn
    //---------------------------------------------------------------------------
    if ( m_lNumberOfContrastsBeforeRTEBPlugIn >= m_lNumberOfContrasts )
    {
        setNLSStatus( MRI_SEQ_SEQU_ERROR, __PRETTY_FUNCTION__, "m_lNumberOfContrastsBeforeRTEBPlugIn must be smaller than m_lNumberOfContrasts!" );
        return false;
    }

    if ( !m_bPlugInAvailable && ( m_lNumberOfContrastsBeforeRTEBPlugIn != 0 ) )
    {
        setNLSStatus( MRI_SEQ_SEQU_ERROR, __PRETTY_FUNCTION__, "Without an RTEB-PlugIn, m_lNumberOfContrastsBeforeRTEBPlugIn has to be zero!" );
        return false;
    }

    if (  m_bPlugInAvailable && ( m_lTEContrastIndex != m_lNumberOfContrasts - 1 ) && ( m_lTEContrastIndex != m_lNumberOfContrastsBeforeRTEBPlugIn ) )
    {
        setNLSStatus( MRI_SEQ_SEQU_ERROR, __PRETTY_FUNCTION__, "With an RTEB-PlugIn, the first or the last contrast index after the RTEB-PlugIn has to realize TE!" );
        return false;
    }

    if ( !m_bPlugInAvailable && ( m_lTEContrastIndex != m_lNumberOfContrasts - 1 ) && ( m_lTEContrastIndex != 0) )
    {
        setNLSStatus( MRI_SEQ_SEQU_ERROR, __PRETTY_FUNCTION__, "Without an RTEB-PlugIn, the first or the last contrast index has to realize TE!" );
        return false;
    }

    //---------------------------------------------------------------------------
    // Assemble TE contributions before PlugIn
    //---------------------------------------------------------------------------
    m_lTEContributionBeforeRTEBPlugIn = 
          m_pSBBExcite->getTEContribution() + m_pSBBExcite->getRequiredHoldTime()           // Excitation
        +                                        m_aGradMoments[0].getTotalTime()           // Excitation refocusing
        + ( m_bInternalPhaseCorrection ? getDurationPhaseCorrScanPerRequest() : 0 )         // Internal phase correction
        +                                        m_aGradMoments[1].getTotalTime()           // Phase correction refocusing
        + m_lNumberOfContrastsBeforeRTEBPlugIn * m_aGradMoments[2].getTotalTime()           // Prephasers
        + m_lNumberOfContrastsBeforeRTEBPlugIn * getDurationEPIReadOutPerRequest()          // Echo trains
        + m_lNumberOfContrastsBeforeRTEBPlugIn * m_aGradMoments[3].getTotalTime();          // Rewinders

    //---------------------------------------------------------------------------
    // Assemble TE contributions after PlugIn
    //---------------------------------------------------------------------------
    m_lTEContributionAfterRTEBPlugIn = 
          ( m_lTEContrastIndex - m_lNumberOfContrastsBeforeRTEBPlugIn ) * m_aGradMoments[2].getTotalTime()      // Prephasers
        + ( m_lTEContrastIndex - m_lNumberOfContrastsBeforeRTEBPlugIn ) * getDurationEPIReadOutPerRequest()     // Echo trains
        + ( m_lTEContrastIndex - m_lNumberOfContrastsBeforeRTEBPlugIn ) * m_aGradMoments[3].getTotalTime()      // Rewinders
        +                                                                  m_aGradMoments[2].getTotalTime()      // Prephaser final echo train
        + m_lCenterSegment * m_lEchoSpacing + m_lEchoSpacing / 2 + lEchoShiftingDelay;                           // Half of final echo train

    m_lRTEBPlugInToCompGradTime
        = (m_lTEContrastIndex - m_lNumberOfContrastsBeforeRTEBPlugIn) * m_aGradMoments[2].getTotalTime()   // Prephasers
          + (m_lTEContrastIndex - m_lNumberOfContrastsBeforeRTEBPlugIn) * m_aGradMoments[3].getTotalTime() // Rewinders
          + m_aGradMoments[2].getTotalTime()                        // Prephaser final echo train
          + m_lTEContrastIndex * getDurationEPIReadOutPerRequest(); // Echo trains

    //---------------------------------------------------------------------------
    // adjust TE contributions due to the additional Gradients applied
    // e.g. flow compensated gradients
    //---------------------------------------------------------------------------
    if (!adjustTEContributions(rMrProt))
    {
        // unexpected error
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "AdjustTEContributionsAfterRTEBPlugIn failed!");

        return false;
    }

    m_lTEContributionBeforeRTEBPlugIn = fSDSRoundDownGRT(m_lTEContributionBeforeRTEBPlugIn);
    m_lTEContributionAfterRTEBPlugIn  = fSDSRoundDownGRT(m_lTEContributionAfterRTEBPlugIn );

    //---------------------------------------------------------------------------
    // Check if the readout train stays within gradient limitations
    // Prerequisite: EPI readout is prepared, timing of gradient moment sections
    //               is calculated.
    //---------------------------------------------------------------------------

    // Reset required TR increment
    m_lTRIncrement = 0;

    if (getUseGPABalance())
    {
        // Some notes on the GPA balance implementation:
        //
        // The basic idea behind the GPA balance implementation of the EPI kernel
        // is to allow for a flexible compromise between readout gradient amplitude
        // (=> pixel bandwidth, FOV, echo spacing) and required cooling pauses.
        // Elevated amplitudes up to GradMaxAmplAbsolute can be used - it is
        // assumed that:
        //
        // 1. Only the logical readout and phase encoding gradients
        //    might constructively overlap on a single physical gradient axis (i.e.
        //    no additional gradients along the slice direction exist).
        // 2. The ramp-up of the phase encoding gradients does not start not before
        //    the ramp-down of the readout gradient.
        // 3. The slewrate used for the phase encoding gradients is not higher
        //    than used for the readout gradients.
        //
        // An intermediate mode (using only up to GradMaxAmplAbsolute / sqrt(2))
        // is supported: here only the condition 1. has to be fulfilled.
        //
        // In a first step, it is checked that the readout gradient events can
        // be applied at all. Phase encoding blips are neglected due to their
        // (presumably) low gradient load.
        // In a second step, this EPI readout balance contribution is passed
        // to the plug-in (if there is one). The plug-in can use this information
        // to set up its gradient events appropriately (e.g. diffusion module)
        // and exports the required cooling pause per request (taking into
        // account the gradient load of the EPI readout). If no plug-in
        // is registered, the required cooling pause is calculated considering
        // the EPI readout contribution including the corrsponding fill times
        // only.
        // The required cooling pause is exported (::getTRIncrement) and the
        // sequence has to take care to apply it appropriately.
        if (!checkBalance())
        {
            // balance check might fail in binary search
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "SeqBuildBlockEPIKernel::checkBalance failed!");
            return false;
        }

        // Indicate that TR increment calculation is required
        // (after plug-in preparation)
        m_lTRIncrement = -1;
    }

    //---------------------------------------------------------------------------
    // prepare RTEB-PlugIn
    //---------------------------------------------------------------------------
    m_lRTEBPlugInTEContribution = 0;
    m_lRTEBPlugInStorageTime = 0;
    m_bRTEBPlugInInvertsMagnetization = false;
    m_lRTEBPlugInDurationPerRequest = 0;

    if (m_bPlugInAvailable)
    {
        if (!prepPlugIn(rMrProt, rSeqLim, rSeqExpo))
        {
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "prepPlugIn failed!");
            return false;
        }
    }

    if (m_lRTEBPlugInTEContribution%GRAD_RASTER_TIME)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_lRTEBPlugInTEContribution%GRAD_RASTER_TIME must be zero!");
        return false;
    }

    if (m_lRTEBPlugInStorageTime%GRAD_RASTER_TIME)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_lRTEBPlugInStorageTime%GRAD_RASTER_TIME must be zero!");
        return false;
    }

    if (m_lRTEBPlugInDurationPerRequest%GRAD_RASTER_TIME)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_lRTEBPlugInDurationPerRequest%GRAD_RASTER_TIME must be zero!");
        return false;
    }

    //---------------------------------------------------------------------------
    // wanted TE should be on gradient raster time if devided by two
    //---------------------------------------------------------------------------
    if (m_lWantedTE % (2 * GRAD_RASTER_TIME))
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_lWantedTE%(2*GRAD_RASTER_TIME) must be zero!");
        return false;
    }

    //---------------------------------------------------------------------------
    // calculate TE-Fill-Times and required TE
    //---------------------------------------------------------------------------
    long lNeedMoreTEBefore = 0;
    long lNeedMoreTEAfter = 0;

    if (m_lRTEBPlugInDurationPerRequest != m_lRTEBPlugInTEContribution)
    {
        // Try satisfying spin-echo condition
        m_lRTEBPlugInTEFillBefore = m_lWantedTE / 2
            - m_lTEContributionBeforeRTEBPlugIn
            - ((m_lRTEBPlugInDurationPerRequest - m_lRTEBPlugInStorageTime) - m_lRTEBPlugInTEContribution);

        m_lRTEBPlugInTEFillAfter = m_lWantedTE / 2
            - m_lRTEBPlugInTEContribution
            - m_lTEContributionAfterRTEBPlugIn;

        if (m_lRTEBPlugInTEFillBefore < 0)
        {
            lNeedMoreTEBefore = -m_lRTEBPlugInTEFillBefore;
            m_lRTEBPlugInTEFillBefore = 0;
        }

        if (m_lRTEBPlugInTEFillAfter < 0)
        {
            lNeedMoreTEAfter = -m_lRTEBPlugInTEFillAfter;
            m_lRTEBPlugInTEFillAfter = 0;
        }

        if (lNeedMoreTEAfter > lNeedMoreTEBefore)
        {
            m_lRTEBPlugInTEFillBefore += lNeedMoreTEAfter - lNeedMoreTEBefore;
        }
        else
        {
            m_lRTEBPlugInTEFillAfter += lNeedMoreTEBefore - lNeedMoreTEAfter;
        }

        m_lNeededTE = m_lWantedTE + 2 * max(lNeedMoreTEBefore, lNeedMoreTEAfter);
    }
    else
    {
        m_lRTEBPlugInTEFillBefore = 0;

        m_lRTEBPlugInTEFillAfter = m_lWantedTE
            - m_lTEContributionBeforeRTEBPlugIn
            - (m_lRTEBPlugInDurationPerRequest - m_lRTEBPlugInStorageTime)
            - m_lTEContributionAfterRTEBPlugIn;

        if (m_lRTEBPlugInTEFillAfter < 0)
        {
            lNeedMoreTEAfter = -m_lRTEBPlugInTEFillAfter;
            m_lRTEBPlugInTEFillAfter = 0;
        }

        m_lNeededTE = m_lWantedTE + lNeedMoreTEAfter;
    }

    //---------------------------------------------------------------------------
    // Calculate required TR increment
    //---------------------------------------------------------------------------

    // Required only if
    // 1. GPA balance model should be used and
    // 2. the plug-in did not take care of the TR increment calculation
    if (m_lTRIncrement < 0)
    {
        // Cooling pause calculation, based on EPI readout and TE fill time only
        // (the latter is considered as an implicit pause contribution)
        long lTEFillTime = m_lRTEBPlugInTEFillBefore + m_lRTEBPlugInTEFillAfter;

        if (lTEFillTime < 0)
        {
            lTEFillTime = 0;
        }

        if (!calcTRIncrement(lTEFillTime))
        {
            // TR increment calculation might fail in binary search
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "calcTRIncrement failed!");
            return false;
        }
    }

    //---------------------------------------------------------------------------
    // calculate start time of excitation SBB within even block
    //---------------------------------------------------------------------------
    m_lSBBExciteStartTime = m_pSBBExcite->getRequiredLeadTime();

    if (m_bUseSyncBits)
    {
        m_lSBBExciteStartTime += m_lMaxSyncBitDuration;
    }

    //---------------------------------------------------------------------------
    // calculate duration per request
    //---------------------------------------------------------------------------
    setSBBDurationPerRequest(m_lSBBExciteStartTime
        + m_pSBBExcite->getDurationPerRequest()
        + m_pSBBExcite->getRequiredHoldTime()
        + m_aGradMoments[0].getTotalTime()
        + (m_bInternalPhaseCorrection ? getDurationPhaseCorrScanPerRequest() : 0)
        + m_aGradMoments[1].getTotalTime()
        + m_lRTEBPlugInTEFillBefore
        + m_lRTEBPlugInDurationPerRequest
        + m_lRTEBPlugInTEFillAfter
        + m_lNumberOfContrasts * m_aGradMoments[2].getTotalTime()
        + m_lNumberOfContrasts * getDurationEPIReadOutPerRequest()
        + m_lNumberOfContrasts * m_aGradMoments[3].getTotalTime() );

    //---------------------------------------------------------------------------
    // calculate remaining duration after echo was measured
    //---------------------------------------------------------------------------
    m_lDurationAfterEcho = getSBBDurationPerRequest()
        - m_lNeededTE
        - fSDSRoundDownGRT(m_pSBBExcite->getDurationPerRequest() - m_pSBBExcite->getTEContribution())
        - m_lSBBExciteStartTime;

    //---------------------------------------------------------------------------
    // calculate start time of navigator (can be used in EPI phase correction)
    //---------------------------------------------------------------------------
    // "Duration per request" is misleading here: what is required is the time between
    // center and end of the excitation module => use duration without dynamic adjustment
    // update.
    m_lNavigatorStartTime = m_pSBBExcite->getDurationPerRequest(SLICEADJ::ADJNONE) / 2 
        + m_pSBBExcite->getRequiredHoldTime() 
        + m_aGradMoments[0].getTotalTime();

    //---------------------------------------------------------------------------
    // calculate actual TE's
    //---------------------------------------------------------------------------
    m_vlActualTE.assign( m_lNumberOfContrasts, 0 );

    const long lTEOffsetBeforeRTEBPlugIn = 
          m_pSBBExcite->getTEContribution() + m_pSBBExcite->getRequiredHoldTime()           // Excitation
        + m_aGradMoments[0].getTotalTime()                                                  // Excitation refocusing
        + ( m_bInternalPhaseCorrection ? getDurationPhaseCorrScanPerRequest() : 0 )         // Internal phase correction
        + m_aGradMoments[1].getTotalTime();                                                 // Phase correction refocusing

    for ( long lI = 0; lI < m_lNumberOfContrasts; ++lI )
    {
        if ( lI < m_lNumberOfContrastsBeforeRTEBPlugIn )
        {
            m_vlActualTE[lI] = lTEOffsetBeforeRTEBPlugIn 
                + lI * m_aGradMoments[2].getTotalTime()                                             // Prephasers
                + lI * getDurationEPIReadOutPerRequest()                                            // Echo trains
                + lI * m_aGradMoments[3].getTotalTime()                                             // Rewinders
                +      m_aGradMoments[2].getTotalTime()                                             // Prephaser final echo train
                + m_lCenterSegment * m_lEchoSpacing + m_lEchoSpacing / 2 + lEchoShiftingDelay;      // Half of final echo train
        }
        else
        {
            m_vlActualTE[lI] = m_lNeededTE
                - ( m_lTEContrastIndex - lI ) * m_aGradMoments[2].getTotalTime()     // Prephasers
                - ( m_lTEContrastIndex - lI ) * getDurationEPIReadOutPerRequest()    // Echo trains
                - ( m_lTEContrastIndex - lI ) * m_aGradMoments[3].getTotalTime();    // Rewinders
        }
    }

    //---------------------------------------------------------------------------
    // decide who has to apply the TE/TR fill times
    //---------------------------------------------------------------------------
    if (m_bPlugInAvailable)
    {
        m_bEPIReadOutAppliesTEFill = false;
        m_bGradMoment2AppliesTEFill = false;
    }
    else
    {
        if (m_aGradMoments[2].getTotalTime() != 0)
        {
            m_bEPIReadOutAppliesTEFill = false;
            m_bGradMoment2AppliesTEFill = true;
        }
        else
        {
            m_bEPIReadOutAppliesTEFill = true;
            m_bGradMoment2AppliesTEFill = false;
        }
    }

    if (m_aGradMoments[3].getTotalTime() != 0)
    {
        m_bEPIReadOutAppliesTRFill = false;
        m_bGradMoment3AppliesTRFill = true;
    }
    else
    {
        m_bEPIReadOutAppliesTRFill = true;
        m_bGradMoment3AppliesTRFill = false;
    }

    //---------------------------------------------------------------------------
    // EPIReadOut has to apply fill times if phase correction is active
    //---------------------------------------------------------------------------
    m_bApplyFillTimesIfPhaseCorrection = true;

    //---------------------------------------------------------------------------
    // adjust additional timing config due to the applied additional Gradients
    // e.g.flow compensated gradients
    //---------------------------------------------------------------------------
    if (!adjustAdditionalTiming(rMrProt))
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "AdjustSBBDurationPerRequest failed!");
        return false;
    }


    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    setPrepared();
    return true;
}


//===============================================================================
//
// Function:    checkGradients()
//
// Description: Checks all gradients involved for gradient specification
//              violations in the LOGICAL coordinate system.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::checkGradients(MrProt &rMrProt, SeqLim &rSeqLim)
{
    //---------------------------------------------------------------------------
    // m_pSBBExcite must have been set
    //---------------------------------------------------------------------------
    if (!m_pSBBExcite)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "missing pointer to excitation SBB");
        return false;
    }

    //---------------------------------------------------------------------------
    // check gradients of excitation SBB
    //---------------------------------------------------------------------------
    if (!m_pSBBExcite->checkGradients(rMrProt, rSeqLim))
    {
        if (rSeqLim.isContextPrepForBinarySearch()) setNLSStatus(m_pSBBExcite->getNLSStatus());
        else                                         setNLSStatus(m_pSBBExcite->getNLSStatus(), __FUNCTION__, "m_pSBBExcite->checkGradients failed!");
        return false;
    }

    //---------------------------------------------------------------------------
    // check gradients of EPI read out
    //---------------------------------------------------------------------------
    if (!SeqBuildBlockEPIReadOut::checkGradients(rMrProt, rSeqLim))
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "SeqBuildBlockEPIReadOut::checkGradients failed!");
        return false;
    }

    //---------------------------------------------------------------------------
    // check gradient moment sections
    //---------------------------------------------------------------------------
    for (long lI = 0; lI < 4; lI++)
    {
        if (!m_aGradMoments[lI].check())
        {
            if (rSeqLim.isContextPrepForBinarySearch()) setNLSStatus(m_aGradMoments[lI].getNLSStatus());
            else                                         setNLSStatus(m_aGradMoments[lI].getNLSStatus(), __PRETTY_FUNCTION__, "m_aGradMoments[lI].check() failed!");
            return false;
        }
    }


    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    return true;
}


//===============================================================================
//
// Function:    getDurationPerRequest()
//
// Description: Returns duration of the SBB of one call of the run function.
//
//===============================================================================

long SeqBuildBlockEPIKernel::getSBBDurationPerRequest()
{
    // Note: Do not call the corresponding SBBEPIReadout method here!
    return SeqBuildBlock::getSBBDurationPerRequest();
}


//===============================================================================
//
// Function:    getRFInfoPerRequest()
//
// Description: Returns energy of the SBB of one call of the run function.
//              Method needs to be overloaded for specific kernel flavors.
//===============================================================================

MrProtocolData::SeqExpoRFInfo SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getRFInfoPerRequest()
{
    if (m_pSBBExcite && isPrepared())
        return m_pSBBExcite->getRFInfoPerRequest();
    else
        return MrProtocolData::SeqExpoRFInfo();
}



//===============================================================================
//
// Function:    increaseTE()
//
// Description: Increases the current TE to a certain value.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::increaseTE(long lNewTE)
{

    //---------------------------------------------------------------------------
    // are we prepared ?
    //---------------------------------------------------------------------------
    if (!isPrepared())
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "SBB is not prepared.");
        return false;
    }

    //---------------------------------------------------------------------------
    // check new TE
    //---------------------------------------------------------------------------
    if (lNewTE < m_lNeededTE)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "can't decrease TE.");
        return false;
    }

    //---------------------------------------------------------------------------
    // check size of actual TE array
    //---------------------------------------------------------------------------
    if ( m_vlActualTE.size() < static_cast<size_t>( m_lNumberOfContrasts ) )
    {
        setNLSStatus( MRI_SEQ_SEQU_ERROR, __PRETTY_FUNCTION__, "m_vlActualTE array size too small." );
        return false;
    }

    //---------------------------------------------------------------------------
    // recalculate fill times and needed TE
    //---------------------------------------------------------------------------
    if (m_lRTEBPlugInDurationPerRequest != m_lRTEBPlugInTEContribution)
    {
        const long lAddFillBeforePlugInTime = fSDSRoundDownGRT(static_cast<long>((lNewTE-m_lNeededTE)/2));
        const long lAddFillAfterPlugInTime  = lAddFillBeforePlugInTime;
        // Update member variables correspondingly
        m_lRTEBPlugInTEFillBefore   += lAddFillBeforePlugInTime;
        m_lRTEBPlugInTEFillAfter    += lAddFillAfterPlugInTime;
        m_lNeededTE                 += lAddFillBeforePlugInTime + lAddFillAfterPlugInTime;
        addSBBDurationPerRequest( lAddFillBeforePlugInTime + lAddFillAfterPlugInTime );

        for ( long lI = m_lNumberOfContrastsBeforeRTEBPlugIn; lI < m_lNumberOfContrasts; ++lI )
        {
            m_vlActualTE[lI] += lAddFillBeforePlugInTime + lAddFillAfterPlugInTime;
        }
    }
    else
    {
        const long lAddFillTime = fSDSRoundDownGRT(lNewTE - m_lNeededTE);

        // Update member variables correspondingly
        m_lRTEBPlugInTEFillAfter += lAddFillTime;
        m_lNeededTE += lAddFillTime;
        addSBBDurationPerRequest(lAddFillTime);

        for ( long lI = m_lNumberOfContrastsBeforeRTEBPlugIn; lI < m_lNumberOfContrasts; ++lI )
        {
            m_vlActualTE[lI] += lAddFillTime;
        }
    }

    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    return true;
}


//===============================================================================
//
// Function:    setTRFill()
//
// Description: Sets the TR-fill that the kernel should apply.
//
//===============================================================================

void SeqBuildBlockEPIKernel::setTRFill(long lTRFill)
{
    if (isPrepared())
    {
        m_lTRFill = fSDSRoundDownGRT(lTRFill);
    }
    else
    {
        m_lTRFill = 0;
    }
}


//    Internal function to run one gradient moment section.
bool SeqBuildBlockEPIKernel::runGradMoments(long lSectionIndex, MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo , sSLICE_POS* pSLC)

{

    double dMoment[3];

    long   lEventBlockEnd = 0;
    long   lFillTimeBeforeGradients = 0;

    //---------------------------------------------------------------------------
    // m_pSBBExcite must have been set
    //---------------------------------------------------------------------------
    if (m_pSBBExcite == nullptr)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "missing pointer to excitation SBB");
        return false;
    }

    //---------------------------------------------------------------------------
    // This is IMPORTANT!
    // Otherwise getPEPrePhasingMoment and get3DPrePhasingMoment may return zero
    // when called next time, if EPIReadOut is in phase correction mode.
    //---------------------------------------------------------------------------
    setNextExecutionIsImaging(false); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // check index
    //---------------------------------------------------------------------------
    if (lSectionIndex < 0 || lSectionIndex>3)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "lSectionIndex out of range");
        return false;
    }

    //---------------------------------------------------------------------------
    // check, if there is something to do
    //---------------------------------------------------------------------------
    if (m_aGradMoments[lSectionIndex].getTotalTime() == 0)
    {
        return true;
    }

    //---------------------------------------------------------------------------
    // determine sign of first RO-gradient of echo train
    //---------------------------------------------------------------------------
    bool bLocalStartImagingReadOutWithNegativeGradient = m_bStartImagingReadOutWithNegativeGradient;

    if (m_bExecuteKernelAsPhaseCorrectionScan && !m_bExternalEPIPhaseCorrection && (m_lCenterSegment - 1) % 2)
    {
        bLocalStartImagingReadOutWithNegativeGradient = !bLocalStartImagingReadOutWithNegativeGradient;
    }

    //---------------------------------------------------------------------------
    // calculate moments ...
    //---------------------------------------------------------------------------
    dMoment[m_lGP] = dMoment[m_lGR] = dMoment[m_lGS] = 0.0;


    //---------------------------------------------------------------------------
    // ... for section 0
    //---------------------------------------------------------------------------
    if (lSectionIndex == 0)
    {
        if (m_bInternalPhaseCorrection)
        {
            dMoment[m_lGR] = -fabs(getROPrePhasingMoment(rMrProt, rSeqLim)); // we do not want the result to depend on
            // SeqBuilBlockEPIReadOut::m_bStartWithNegativeROGrad
            if (!m_bPrephaseROAfterRTEBPlugIn)
            {
                dMoment[m_lGR] *= -1.0; // we know: internal phase correction scan has 3 echoes

                if (m_bRTEBPlugInInvertsMagnetization)
                {
                    dMoment[m_lGR] *= -1.0;
                }

                if (bLocalStartImagingReadOutWithNegativeGradient)
                {
                    dMoment[m_lGR] *= -1.0;
                }
            }
        }
        else if (m_bEarlyFIDPhaseCorrection && m_bExecuteKernelAsPhaseCorrectionScan)
        {
            dMoment[m_lGR] = -fabs(getROPrePhasingMoment(rMrProt, rSeqLim));
        }
        else
        {
            if (!m_bPrephaseROAfterRTEBPlugIn)
            {
                dMoment[m_lGR] = -fabs(getROPrePhasingMoment(rMrProt, rSeqLim)); // we do not want the result to depend on
                // SeqBuilBlockEPIReadOut::m_bStartWithNegativeROGrad
                if (m_bRTEBPlugInInvertsMagnetization)
                {
                    dMoment[m_lGR] *= -1.0;
            }

                if (bLocalStartImagingReadOutWithNegativeGradient)
                {
                    dMoment[m_lGR] *= -1.0;
                }
        }
    }

        if (SysProperties::isPhaseEncodingEnabled()
            && (!m_bInternalPhaseCorrection && !m_bPrephaseBlipsAfterRTEBPlugIn)
            && !m_bExecuteKernelAsPhaseCorrectionScan
            )
        {
            if (m_bRTEBPlugInInvertsMagnetization)
            {
                dMoment[m_lGP] -= getPEPrePhasingMoment(rMrProt, rSeqLim);
                dMoment[m_lGS] -= get3DPrePhasingMoment(rMrProt, rSeqLim);
            }
            else
            {
                dMoment[m_lGP] += getPEPrePhasingMoment(rMrProt, rSeqLim);
                dMoment[m_lGS] += get3DPrePhasingMoment(rMrProt, rSeqLim);
            }
        }


        if ( m_pSBBExcite->getGSData().getbAxisActive()) dMoment[m_lGS] += m_pSBBExcite->getGSData().getdRequiredRefocusingMoment();
        if ( m_pSBBExcite->getGPData().getbAxisActive()) dMoment[m_lGP] += m_pSBBExcite->getGPData().getdRequiredRefocusingMoment();
        if ( m_pSBBExcite->getGRData().getbAxisActive()) dMoment[m_lGR] += m_pSBBExcite->getGRData().getdRequiredRefocusingMoment();
}


    //---------------------------------------------------------------------------
    // ... for section 1
    //---------------------------------------------------------------------------
    if (lSectionIndex == 1)
    {
        if (m_bInternalPhaseCorrection)
        {
            if (m_bPrephaseROAfterRTEBPlugIn)
            {
                dMoment[m_lGR] = -fabs(getROPrePhasingMoment(rMrProt, rSeqLim)); // we do not want the result to depend on
                // SeqBuilBlockEPIReadOut::m_bStartWithNegativeROGrad
            }

            if (SysProperties::isPhaseEncodingEnabled()
                && !m_bPrephaseBlipsAfterRTEBPlugIn
                && !m_bExecuteKernelAsPhaseCorrectionScan
                )
            {
                if (m_bRTEBPlugInInvertsMagnetization)
                {
                    dMoment[m_lGP] -= getPEPrePhasingMoment(rMrProt, rSeqLim);
                    dMoment[m_lGS] -= get3DPrePhasingMoment(rMrProt, rSeqLim);
                }
                else
                {
                    dMoment[m_lGP] += getPEPrePhasingMoment(rMrProt, rSeqLim);
                    dMoment[m_lGS] += get3DPrePhasingMoment(rMrProt, rSeqLim);
                }
            }
        }
    }


    //---------------------------------------------------------------------------
    // ... for section 2
    //---------------------------------------------------------------------------
    if (lSectionIndex == 2)
    {
        if (m_bPrephaseROAfterRTEBPlugIn)
        {
            dMoment[m_lGR] = -fabs(getROPrePhasingMoment(rMrProt, rSeqLim)); // we do not want the result to depend on
            // SeqBuilBlockEPIReadOut::m_bStartWithNegativeROGrad
            if (bLocalStartImagingReadOutWithNegativeGradient)
            {
                dMoment[m_lGR] *= -1.0;
            }
        }


        if (SysProperties::isPhaseEncodingEnabled()
            && m_bPrephaseBlipsAfterRTEBPlugIn
            && !m_bExecuteKernelAsPhaseCorrectionScan
            )
        {
            dMoment[m_lGP] = getPEPrePhasingMoment(rMrProt, rSeqLim);
            dMoment[m_lGS] = get3DPrePhasingMoment(rMrProt, rSeqLim);



        }
    }


    //---------------------------------------------------------------------------
    // ... for section 3
    //---------------------------------------------------------------------------
    if (lSectionIndex == 3)
    {
        if ( m_bRewindRO || m_bSpoilRO )
        {
            dMoment[m_lGR] = -fabs(getROPrePhasingMoment(rMrProt, rSeqLim)); // we do not want the result to depend on
            // SeqBuilBlockEPIReadOut::m_bStartWithNegativeROGrad

            if (!m_bInternalPhaseCorrection && m_bExecuteKernelAsPhaseCorrectionScan && !m_bExternalEPIPhaseCorrection)
            {
                if (!m_bEarlyFIDPhaseCorrection)
                {
                    if ((m_lCenterSegment - 1) % 2)
                    {
                        dMoment[m_lGR] *= -1.0;
                    }

                    if (m_bStartImagingReadOutWithNegativeGradient)
                    {
                        dMoment[m_lGR] *= -1.0;
                    }
                }
            }
            else
            {
                if ((m_lEchoTrainLength % 2) == 0)
                {
                    dMoment[m_lGR] *= -1.0;
                }

                if (m_bStartImagingReadOutWithNegativeGradient)
                {
                    dMoment[m_lGR] *= -1.0;
                }
            }

            // Spoiling rather than rewinding
            if ( m_bSpoilRO )
            {
                dMoment[m_lGR] *= -m_dSpoilROFactor;
            }
        }

        if (m_bRewindBlips && !m_bExecuteKernelAsPhaseCorrectionScan)
        {
            dMoment[m_lGP] = getPERePhasingMoment(rMrProt, rSeqLim);
            dMoment[m_lGS] = get3DRePhasingMoment(rMrProt, rSeqLim);

        }

        if (m_bPrefaceExcitationSBB)
        {
            if (m_pSBBExcite->getGSData().getbAxisActive()) dMoment[m_lGS] -= m_pSBBExcite->getGSData().getdMomentBeforeFirstRFCenter();
            if (m_pSBBExcite->getGPData().getbAxisActive()) dMoment[m_lGP] -= m_pSBBExcite->getGPData().getdMomentBeforeFirstRFCenter();
            if (m_pSBBExcite->getGRData().getbAxisActive()) dMoment[m_lGR] -= m_pSBBExcite->getGRData().getdMomentBeforeFirstRFCenter();
        }
    }

    //---------------------------------------------------------------------------
    // prepare and check gradients to be switched
    //---------------------------------------------------------------------------
    for (long lI = 0; lI < 3; lI++)
    {
        sprintf(m_tIdent, "%s [%ld]", m_aGradMoments[lSectionIndex].getIdent(), lI);

        m_sGrad[lI].setIdent(m_tIdent);
        m_sGrad[lI].setMinRiseTime(m_aGradMoments[lSectionIndex].getMinRiseTime());
        m_sGrad[lI].setMaxMagnitude(m_aGradMoments[lSectionIndex].getMaxMagnitude());
        m_sGrad[lI].setRampUpTime(m_aGradMoments[lSectionIndex].getRampUpTime());
        m_sGrad[lI].setDuration(m_aGradMoments[lSectionIndex].getDuration());
        m_sGrad[lI].setRampDownTime(m_aGradMoments[lSectionIndex].getRampDownTime());

        if (!m_sGrad[lI].prepMomentumTOT(dMoment[lI]))
        {
            // unexpected error
            // => trace also if pSeqLim->isContextPrepForBinarySearch()
            setNLSStatus(m_sGrad[lI].getNLSStatus(), __FUNCTION__, "m_sGrad[lI].prepMomentumTOT(dMoment[lI]) failed!");
            return false;
        }


        if (!m_sGrad[lI].check())
        {
            // caused by programming error
            // => trace also if pSeqLim->isContextPrepForBinarySearch()
            setNLSStatus(m_sGrad[lI].getNLSStatus(), __FUNCTION__, "unexpected, but true: m_sGrad[lI].check() failed!");
            return false;
        }
    }

    //---------------------------------------------------------------------------
    // do we have to apply any fill times ?
    //---------------------------------------------------------------------------
    if (lSectionIndex == 2)
    {
        if (m_bGradMoment2AppliesTEFill)
        {
            lFillTimeBeforeGradients = m_lRTEBPlugInTEFillAfter;
        }
    }

    lEventBlockEnd = lFillTimeBeforeGradients + m_aGradMoments[lSectionIndex].getTotalTime();

    if (lSectionIndex == 3)
    {
        if (m_bGradMoment3AppliesTRFill)
        {
            lEventBlockEnd += m_lTRFill;
        }
    }


    //---------------------------------------------------------------------------
    // update additional gradient moment applied
    // e.g. flow compensated gradients
    //---------------------------------------------------------------------------
    if (!updateAdditionalGradMoments(lSectionIndex, rMrProt, rSeqLim, lFillTimeBeforeGradients, lEventBlockEnd))
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "updatedAdditonalGradMoments failed");
        return false;
    }


    //---------------------------------------------------------------------------
    // run gradients prepared above
    //---------------------------------------------------------------------------

    if (!runFinalGradMoments(lSectionIndex, rMrProt, rSeqLim, rSeqExpo, pSLC, lFillTimeBeforeGradients, lEventBlockEnd))
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "runFinalGradMoments failed");
        return false;
    }


    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    return true;
}


//    Internal function to run the sync-bits.

bool SeqBuildBlockEPIKernel::runSyncBits(MrProt&, SeqLim&, SeqExpo&, sSLICE_POS*)
{
    if (!m_bDoNotSendOscBit)
    {
        m_OscBit.run();
    }

    if (!m_bDoNotSendExtTrigger)
    {
        m_ExtTrig.run();
    }

    return true;
}


//===============================================================================
//
// Function:    runExcitationAndSyncBits()
//
// Description: Internal function to run the sync-bit and the excitation SBB.
//              Calls runSync Bits.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::runExcitationAndSyncBits(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
{

    //---------------------------------------------------------------------------
    // m_pSBBExcite must have been set
    //---------------------------------------------------------------------------
    if (!m_pSBBExcite)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "missing pointer to excitation SBB");
        return false;
    }

    //---------------------------------------------------------------------------
    // when to close the event block
    //---------------------------------------------------------------------------
    long lEventBlockEnd = m_lSBBExciteStartTime
        + m_pSBBExcite->getDurationPerRequest()
        + m_pSBBExcite->getRequiredHoldTime();

    //---------------------------------------------------------------------------
    // put start time in event block into excitation SBB
    //---------------------------------------------------------------------------
    m_pSBBExcite->setStartTimeInEventBlock(m_lSBBExciteStartTime);

    //---------------------------------------------------------------------------
    // open event block
    //---------------------------------------------------------------------------
    fRTEBInit(pSLC->getROT_MATRIX()); /*! EGA-07 !*/


    //---------------------------------------------------------------------------
    // Provide slice-specific adjustment information to the Ice world:
    // Certain functors (e.g. static field correction) might require this
    // information for proper operation.
    //---------------------------------------------------------------------------
    if ( getsSliceAdjParametersRelevantForUpdate().isAdjFre() ||
         getsSliceAdjParametersRelevantForUpdate().isAdjShim() )
    {
        if ( m_bSendSliceAdjustData )
        {
            // Obtain cuboid for current slice object
            SLICEADJ::sCuboidGeometry sCurrentCuboid;
            if ( !calcSliceAdjSliceCuboid( rMrProt, rSeqLim, rSeqExpo, pSLC, sCurrentCuboid ) )
            {
                setNLSStatus( MRI_SEQ_SEQU_ERROR, __FUNCTION__, "Cuboid geometry for current slice not available." );
                return false;
            }
            // Get parameter set for current cuboid
            SLICEADJ::sControlParameters sCurrentControlParameters;
            if ( !SLICEADJ::ControlParameterBox::getParameters( rMrProt, sCurrentCuboid, getsSliceAdjParametersRelevantForUpdate(), sCurrentControlParameters ) )
            {
                setNLSStatus( MRI_SEQ_SEQU_ERROR, __FUNCTION__, "Adjustment data for current cuboid not available." );
                return false;
            }

            SFC::SeqToIceSliceAdjustData sSliceAdjustData;

            sSliceAdjustData.usRepIndex   = static_cast<unsigned short>( m_ADC.getMDH().getCrep() );
            sSliceAdjustData.usSlcIndex   = static_cast<unsigned short>( m_ADC.getMDH().getCslc() );
            sSliceAdjustData.lFrequency   = +sCurrentControlParameters.lFrequency;          // Consider applied frequency offset (rather than compensation field => positive sign)
            sSliceAdjustData.dGradOffsetX = -sCurrentControlParameters.dGradOffsetX_mTm;    // Consider applied x-axis compensation gradient
            sSliceAdjustData.dGradOffsetY = -sCurrentControlParameters.dGradOffsetY_mTm;    // Consider applied y-axis compensation gradient
            sSliceAdjustData.dGradOffsetZ = -sCurrentControlParameters.dGradOffsetZ_mTm;    // Consider applied z-axis compensation gradient


            if ( !m_sSliceAdjustSEQData.setID( SFC::SEQ2ICE_SLICEADJUSTDATA_ID.c_str() ) )
            {
                setNLSStatus( MRI_SEQ_SEQU_ERROR, __FUNCTION__, "Setting the ID of SEQData object failed." );
                return false;
            }

            if ( !m_sSliceAdjustSEQData.setData( &sSliceAdjustData, sizeof( sSliceAdjustData ) ) )
            {
                setNLSStatus( MRI_SEQ_SEQU_ERROR, __FUNCTION__, "Setting the content of SEQData object failed." );
                return false;
            }

            NLSStatus lStatus = m_sSliceAdjustSyncData.setData( m_sSliceAdjustSEQData );
            if ( setNLSStatus( lStatus, __FUNCTION__, "Setting the content of SyncData object failed." ) )
            {
                return false;
            }

            lStatus = m_sSliceAdjustSyncData.prep();
            if ( setNLSStatus( lStatus, __FUNCTION__, "Preparation of SyncData object failed." ) )
            {
                return false;
            }

            // Send SyncData object to Ice
            fRTEI( 0, &m_sSliceAdjustSyncData );
        }
    }

    //---------------------------------------------------------------------------
    // sync-bits
    //---------------------------------------------------------------------------
    if (m_bUseSyncBits)
    {
        if (!runSyncBits(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
        {
            SEQ_TRACE_ERROR.print("runSyncBits failed!");
            return false;
        }
    }

    //---------------------------------------------------------------------------
    // excitation SBB
    //---------------------------------------------------------------------------

    // In case of multiband scans set runmode accordingly
    SBBMultibandRF* pSBBMultibandRF = dynamic_cast<SBBMultibandRF*>(m_pSBBExcite);
    if (pSBBMultibandRF != nullptr)
        pSBBMultibandRF->setRunMode(m_eRunMode);

    if (!m_pSBBExcite->run(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
    {
        setNLSStatus(m_pSBBExcite->getNLSStatus(), __FUNCTION__, "m_pSBBExcite->run failed!");
        return false;
    }

    //---------------------------------------------------------------------------
    // run additional gradients in excitation event block
    // e.g. flow compensated slice refocusing gradient  
    //---------------------------------------------------------------------------
    if (!runAdditionalGradExcitation(lEventBlockEnd))
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "runAdditionalGradExcitation failed");
        return false;
    }

    //---------------------------------------------------------------------------
    // close event block
    //---------------------------------------------------------------------------
    fRTEI(lEventBlockEnd, 0, 0, 0, 0, 0, 0, 0); /*! EGA-All !*/
    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 0, 0, pSLC->getSliceIndex(), 0, 0); /*! EGA-All !*/

    const NLS_STATUS lStatus = fRTEBFinish();

    if (setNLSStatus(lStatus, __FUNCTION__, "fRTEBFinish for SBBExcitation failed!"))
    {
        // lStatus is error code, ergo:
        //
        return false;
    }

    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    return true;
}


//===============================================================================
//
// Function:    runInternalPhaseCorrection()
//
// Description: Internal function to run the internal phase correction echos.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::runInternalPhaseCorrection(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
{
    const long lStoreCountersPerSegmentForEchoShifting = m_lCountersPerSegmentForEchoShifting;
    bool bSuccess = false;

    //---------------------------------------------------------------------------
    // determine sign of first RO-gradient
    //---------------------------------------------------------------------------
    bool bFirstRONegative = false;

    if (m_bPrephaseROAfterRTEBPlugIn)
    {
        bFirstRONegative = false;
    }
    else
    {
        if (m_bRTEBPlugInInvertsMagnetization)
        {
            bFirstRONegative = m_bStartImagingReadOutWithNegativeGradient;
        }
        else
        {
            bFirstRONegative = !m_bStartImagingReadOutWithNegativeGradient;
        }
    }

    setNextExecutionIsPhaseCorrection(bFirstRONegative); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // take care that EPIReadOut does not apply echo-shifting fill times
    //---------------------------------------------------------------------------
    m_lCountersPerSegmentForEchoShifting = 0;

    //---------------------------------------------------------------------------
    // run internal phase correction scans
    //---------------------------------------------------------------------------
    // Important: We have to call ::runSBB() here and not ::run()!
    // (SBBEPIKernel takes care of preparations related to dynamic adjustments.)
    bSuccess = SeqBuildBlockEPIReadOut::runSBB(rMrProt, rSeqLim, rSeqExpo, pSLC); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // reset m_lCountersPerSegmentForEchoShifting
    //---------------------------------------------------------------------------
    m_lCountersPerSegmentForEchoShifting = lStoreCountersPerSegmentForEchoShifting;

    //---------------------------------------------------------------------------
    // check success of internal phase correction scans
    //---------------------------------------------------------------------------
    if (!bSuccess)
    {
        SEQ_TRACE_ERROR.print("SeqBuildBlockEPIReadOut::run failed.");
        return false;
    }

    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    return true;
}


//===============================================================================
//
// Function:    runEarlyFIDPhaseCorrection()
//
// Description: Internal function to run the early FID phase correction echos.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::runEarlyFIDPhaseCorrection(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
{
    const long lStoreCountersPerSegmentForEchoShifting = m_lCountersPerSegmentForEchoShifting;
    bool bSuccess = false;

    //---------------------------------------------------------------------------
    // m_pSBBExcite must have been set
    //---------------------------------------------------------------------------
    if (m_pSBBExcite == nullptr)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "missing pointer to excitation SBB");
        return false;
    }

    //---------------------------------------------------------------------------
    // determine sign of first RO-gradient
    //---------------------------------------------------------------------------
    setNextExecutionIsPhaseCorrection(false); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // calculate fill time after phase correction scans
    //---------------------------------------------------------------------------
    m_lLocalTRFill = getSBBDurationPerRequest()
        - m_lSBBExciteStartTime
        - m_pSBBExcite->getDurationPerRequest()
        - m_pSBBExcite->getRequiredHoldTime()
        - m_aGradMoments[0].getTotalTime()
        - getDurationPhaseCorrScanPerRequest()
        - m_aGradMoments[3].getTotalTime();

    
    //---------------------------------------------------------------------------
    // adjust local TR fill due to the additional Gradients applied 
    // e.g. flow compensated gradients
    //---------------------------------------------------------------------------
    if (!adjustLocalTRFill(rMrProt))
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "AdjustLocalTRFill failed");
        return false;
    }


    if (m_lLocalTRFill < 0)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "can not insert negative fill time after early FID phase correction scan!");
        return false;
    }

    //---------------------------------------------------------------------------
    // decide, if EPIReadOut has to apply TR-Fill
    //---------------------------------------------------------------------------
    if (m_bEPIReadOutAppliesTRFill)
    {
        m_lLocalTRFill += m_lTRFill;
    }

    //---------------------------------------------------------------------------
    // check correct configuration of EPIReadOut
    //---------------------------------------------------------------------------
    if (!m_bApplyFillTimesIfPhaseCorrection)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_bApplyFillTimesIfPhaseCorrection must be true!");
        return false;
    }

    //---------------------------------------------------------------------------
    // take care that EPIReadOut does not apply echo-shifting
    //---------------------------------------------------------------------------
    m_lCountersPerSegmentForEchoShifting = 0;

    //---------------------------------------------------------------------------
    // run early fid phase correction scans
    //---------------------------------------------------------------------------
    // Important: We have to call ::runSBB() here and not ::run()!
    // (SBBEPIKernel takes care of preparations related to dynamic adjustments.)
    bSuccess = SeqBuildBlockEPIReadOut::runSBB(rMrProt, rSeqLim, rSeqExpo, pSLC); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // reset m_lCountersPerSegmentForEchoShifting
    //---------------------------------------------------------------------------
    m_lCountersPerSegmentForEchoShifting = lStoreCountersPerSegmentForEchoShifting;

    //---------------------------------------------------------------------------
    // check success of internal phase correction scans
    //---------------------------------------------------------------------------
    if (!bSuccess)
    {
        SEQ_TRACE_ERROR.print("SeqBuildBlockEPIReadOut::run failed.");
        return false;
    }

    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    return true;
}


//===============================================================================
//
// Function:    runEPIReadOutPhaseCorrection()
//
// Description: Internal function to run the phase correction EPI read-out and
//              the appropriate fill times, if the kernel is in the phase
//              correction mode.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::runEPIReadOutPhaseCorrection(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
{
    bool bSuccess = false;

    //---------------------------------------------------------------------------
    // check line in segment for echo-shifting
    //---------------------------------------------------------------------------
    if (m_lCounterInSegmentForEchoShifting != 0)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "expect m_lCounterInSegmentForEchoShifting==0, this is not the case.");
        return false;
    }

    //---------------------------------------------------------------------------
    // take care that EPIReadOut applies correct echo-shifting delay
    //---------------------------------------------------------------------------
    m_lCounterInSegmentForEchoShifting = m_lCounterInSegmentWithEcho;

    //---------------------------------------------------------------------------
    // determine fill times
    //---------------------------------------------------------------------------
    if (m_bEPIReadOutAppliesTEFill) m_lLocalTEFill = m_lRTEBPlugInTEFillAfter;
    else                            m_lLocalTEFill = 0;

    if (m_bEPIReadOutAppliesTRFill) m_lLocalTRFill = m_lTRFill;
    else                            m_lLocalTRFill = 0;

    m_lLocalTEFill += (m_lCenterSegment - 1)*getEchoSpacing();
    m_lLocalTRFill += (m_lEchoTrainLength - m_lCenterSegment - 2)*getEchoSpacing();

    if (m_lLocalTEFill < 0)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "can not insert negative fill time before EPIReadOut!");
        return false;
    }

    if (m_lLocalTRFill < 0)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "can not insert negative fill time after EPIReadOut!");
        return false;
    }

    //---------------------------------------------------------------------------
    // determine sign of first RO-gradient
    //---------------------------------------------------------------------------
    bool bFirstRONegative = m_bStartImagingReadOutWithNegativeGradient;

    if ((m_lCenterSegment - 1) % 2)
    {
        bFirstRONegative = !bFirstRONegative;
    }

    setNextExecutionIsPhaseCorrection(bFirstRONegative); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // check correct configuration of EPIReadOut
    //---------------------------------------------------------------------------
    if (!m_bApplyFillTimesIfPhaseCorrection)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_bApplyFillTimesIfPhaseCorrection must be true!");
        return false;
    }


    //---------------------------------------------------------------------------
    // run internal phase correction scans
    //---------------------------------------------------------------------------
    // Important: We have to call ::runSBB() here and not ::run()!
    // (SBBEPIKernel takes care of preparations related to dynamic adjustments.)
    bSuccess = SeqBuildBlockEPIReadOut::runSBB(rMrProt, rSeqLim, rSeqExpo, pSLC); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // reset line in segment for echo-shifting
    //---------------------------------------------------------------------------
    m_lCounterInSegmentForEchoShifting = 0;

    //---------------------------------------------------------------------------
    // check success of internal phase correction scans
    //---------------------------------------------------------------------------
    if (!bSuccess)
    {
        SEQ_TRACE_ERROR.print("SeqBuildBlockEPIReadOut::run failed.");
            return false;
    }

    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    return true;
}


//===============================================================================
//
// Function:    runEPIReadOut()
//
// Description: Internal function to run the imaging EPI read-out.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::runEPIReadOut(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
{
    //---------------------------------------------------------------------------
    // determine sign of first RO-gradient
    //---------------------------------------------------------------------------
    setNextExecutionIsImaging(m_bStartImagingReadOutWithNegativeGradient); /*! EGA-07 !*/

    //---------------------------------------------------------------------------
    // decide, if EPIReadOut has to apply TE/TR-Fill
    //---------------------------------------------------------------------------
    if (m_bEPIReadOutAppliesTEFill) m_lLocalTEFill = m_lRTEBPlugInTEFillAfter;
    else                            m_lLocalTEFill = 0;

    if (m_bEPIReadOutAppliesTRFill) m_lLocalTRFill = m_lTRFill;
    else                            m_lLocalTRFill = 0;

    //---------------------------------------------------------------------------
    // run imaging EPI read out
    //---------------------------------------------------------------------------
    // Important: We have to call ::runSBB() here and not ::run()!
    // (SBBEPIKernel takes care of preparations related to dynamic adjustments.)
    if (!SeqBuildBlockEPIReadOut::runSBB(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
    {
        SEQ_TRACE_ERROR.print("SeqBuildBlockEPIReadOut::run failed.");
        return false;
    }

    //---------------------------------------------------------------------------
    // if this point is reached successfully
    //---------------------------------------------------------------------------
    return true;
}


//    Internal function to run the RTEB-plug-in.

bool SeqBuildBlockEPIKernel::runRTEBPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
{
    // nothing to do here: all variants which use a plugin have to overload this method and play out events there
    // this implementation does not do anything 
    SEQ_TRACE_ALWAYS.print("SeqBuildBlockEPIKernel::runRTEBPlugIn is not implemented");
    return false;
}


//===============================================================================
//
// Function:    run()
//
// Description: Executes the real-time part of the SBB.
//
//===============================================================================

bool SeqBuildBlockEPIKernel::runSBB(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
{
    //---------------------------------------------------------------------------
    // execute run function only, if SBB is prepared
    //---------------------------------------------------------------------------
    if (!isPrepared())
    {
        SEQ_TRACE_ERROR.print("WARNING: SBB not prepared => nothing to do");
        return true;
    }

    //---------------------------------------------------------------------------
    // we do not like other people to access m_lLocalTEFill and m_lLocalTRFill
    // of SeqBuildBlockSBBEPIReadOut
    //---------------------------------------------------------------------------
    if (m_lLocalTEFill)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_lLocalTEFill was set outside SeqBuildBlockEPIKernel. We're not gonna take it.");
        return false;
    }

    if (m_lLocalTRFill)
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "m_lLocalTRFill was set outside SeqBuildBlockEPIKernel. We're not gonna take it.");
        return false;
    }

    //---------------------------------------------------------------------------
    // send clock check event for unit test (SBBBinomialPulses does not do it)
    //---------------------------------------------------------------------------
    // Note: Has to get applied before / in the RTEB that plays out the excitation pulse.
    //       (Take care of transit RTEB applied by dynamic adjustments!)
    //if (SeqUT.isUnitTestActive())
    //{
    //    mSEQTest(rMrProt,rSeqLim,rSeqExpo,RTEB_ClockCheck,27,-1,pSLC->getSliceIndex(),0,0); /*! EGA-All !*/
    //}

    //---------------------------------------------------------------------------
    // Disable unsupported test cases for current kernel execution
    //---------------------------------------------------------------------------
    if (SeqUT.isUnitTestActive())
    {
        if ( m_bIsSliceAcceleration && ( m_eRunMode == MULTI_BAND ) )
        {
            SeqUT.DisableTestCase( lFrequNotEquOmegaRFErr, RTEB_ORIGIN_fSEQRunKernel,
                                   "Disabled for slice accelerated sequence. "
                                   "Image orientation and image geometry are checked in a functional scanner test." );
        }

        if ( rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED )
        {
            SeqUT.DisableTestCase( lAmplSignRFErr, RTEB_ORIGIN_fSEQRunKernel,
                                   "ZOOM_2DRF slice selection has a 2D trajectory" );
            SeqUT.DisableTestCase( lRFAmplValErr, RTEB_ORIGIN_fSEQRunKernel,
                                   "ZOOM_2DRF slice-select gradient amplitude is alternating with trajectory"
                                   "Image orientation and image geometry are checked in a functional scanner test" );
        }
    }

    //---------------------------------------------------------------------------
    // excitation SBB and Osc-Bit
    //---------------------------------------------------------------------------
    if (!runExcitationAndSyncBits(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
    {
        SEQ_TRACE_ERROR.print("runExcitationAndSyncBits failed.");
        return false;
    }

    //---------------------------------------------------------------------------
    // Re-enable test cases
    //---------------------------------------------------------------------------
    if (SeqUT.isUnitTestActive())
    {
        if ( rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED )
        {
            SeqUT.EnableTestCase( lAmplSignRFErr, RTEB_ORIGIN_fSEQRunKernel );
            SeqUT.EnableTestCase( lRFAmplValErr,  RTEB_ORIGIN_fSEQRunKernel );
        }
    }

    //---------------------------------------------------------------------------
    // gradient moment section 0
    //---------------------------------------------------------------------------
    if (!runGradMoments(0, rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
    {
        SEQ_TRACE_ERROR.print("runGradMoments for section 0 failed.");
        return false;
    }

    //---------------------------------------------------------------------------
    // run internal phase correction scans
    //---------------------------------------------------------------------------
    if (m_bInternalPhaseCorrection)
    {
        if (!runInternalPhaseCorrection(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
        {
            SEQ_TRACE_ERROR.print("runInternalPhaseCorrection failed.");
            return false;
        }
    }

    //---------------------------------------------------------------------------
    // run early fid phase correction scans, gradient moment section 3 and exit
    //---------------------------------------------------------------------------
    if (m_bEarlyFIDPhaseCorrection && m_bExecuteKernelAsPhaseCorrectionScan)
    {
        if (!runEarlyFIDPhaseCorrection(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
        {
            SEQ_TRACE_ERROR.print("runEarlyFIDPhaseCorrection failed.");
            return false;
        }

        if (!runGradMoments(3, rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
        {
            SEQ_TRACE_ERROR.print("runGradMoments for section 3 failed.");
            return false;
        }

        m_lLocalTEFill = 0;
        m_lLocalTRFill = 0;

        return true;
    }

    //---------------------------------------------------------------------------
    // gradient moment section 1
    //---------------------------------------------------------------------------
    if (!runGradMoments(1, rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
    {
        SEQ_TRACE_ERROR.print("runGradMoments for section 1 failed.");
        return false;
    }

    //---------------------------------------------------------------------------
    // contrasts loop
    //---------------------------------------------------------------------------
    // Store fill time and status of LastScanInConcat and LastScanInMeas
    const bool bIsLastScanInMeas   = m_ADC.getMDH().isLastScanInMeas();
    const bool bIsLastScanInConcat = m_ADC.getMDH().isLastScanInConcat();
    const long lTRFill             = m_lTRFill;
    const long lTEFill             = m_lRTEBPlugInTEFillAfter;
    // Remember status of readout events
    bool bIsReadoutEnabled   = fRTIsReadoutEnabled();

    for ( long lContrastIndex = 0; lContrastIndex < m_lNumberOfContrasts; ++lContrastIndex )
    {
        //---------------------------------------------------------------------------
        // Apply LastScan-flags and TR fill time for the last contrast only
        //---------------------------------------------------------------------------
        const bool bIsLastContrast = (lContrastIndex == (m_lNumberOfContrasts - 1));

        m_ADC.getMDH().setLastScanInMeas(bIsLastContrast ? bIsLastScanInMeas : false);
        m_ADC.getMDH().setLastScanInConcat(bIsLastContrast ? bIsLastScanInConcat : false);
        m_lTRFill = bIsLastContrast ? lTRFill : 0;

        //---------------------------------------------------------------------------
        // Apply TE fill time for the first contrast only
        //---------------------------------------------------------------------------
        const bool bIsFirstContrast = (lContrastIndex == 0);

        // If no PlugIn is specified, TE fill time will get applied for the first contrast only
        m_lRTEBPlugInTEFillAfter = (m_bPlugInAvailable || bIsFirstContrast) ? lTEFill : 0;

        if (SeqUT.isUnitTestActive())
        {
            if (m_bPlugInAvailable)
            {
                // Spin-echo condition should be valid for the last contrast (without fill time redistribution) only
                const bool bExpectSERefocused
                    = (lContrastIndex == m_lTEContrastIndex);

                if (bExpectSERefocused)
                {
                    SeqUT.EnableTestCase(lSERefocErr, RTEB_ORIGIN_fSEQRunKernel);
                }
                else
                {
                    SeqUT.DisableTestCase(
                        lSERefocErr, RTEB_ORIGIN_fSEQRunKernel, "Spin-echo condition deliberately violated.");
                }
            }
        }

        //---------------------------------------------------------------------------
        // acquire external phase correction data for the first contrast only
        //---------------------------------------------------------------------------
        if (m_bExecuteKernelAsPhaseCorrectionScan && m_bExternalEPIPhaseCorrection && bIsReadoutEnabled)
        {
            // Without plug-in => use first contrast
            // With    plug-in => use first contrast after plug-in
            if (lContrastIndex == m_lNumberOfContrastsBeforeRTEBPlugIn)
            {
                fRTSetReadoutEnable(1);
            }
            else
            {
                fRTSetReadoutEnable(0);
            }
        }

        //---------------------------------------------------------------------------
        // execute RTEBPlugIn after having acquired all preceding contrasts
        //---------------------------------------------------------------------------
        if (lContrastIndex == m_lNumberOfContrastsBeforeRTEBPlugIn)
        {
            if (m_bPlugInAvailable)
            {
                if (!runRTEBPlugIn(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
                {
            		SEQ_TRACE_ERROR.print("runRTEBPlugIn failed.");
                    return false;
                }
            }
        }

        //---------------------------------------------------------------------------
        // gradient moment section 2
        //---------------------------------------------------------------------------
        if (!runGradMoments(2, rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
        {
            SEQ_TRACE_ERROR.print("runGradMoments for section 2 failed.");
            return false;
        }

        //---------------------------------------------------------------------------
        // EPI read out train
        //---------------------------------------------------------------------------
        if (m_bExecuteKernelAsPhaseCorrectionScan && !m_bExternalEPIPhaseCorrection)
        {
            if (!runEPIReadOutPhaseCorrection(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
            {
                SEQ_TRACE_ERROR.print("runEPIReadOutPhaseCorrection failed.");
                return false;
            }
        }
        else
        {
            if (!runEPIReadOut(rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
            {
                SEQ_TRACE_ERROR.print("runEPIReadOut failed.");
                return false;
            }
        }

        //---------------------------------------------------------------------------
        // gradient moment section 3
        //---------------------------------------------------------------------------
        if (!runGradMoments(3, rMrProt, rSeqLim, rSeqExpo, pSLC)) /*! EGA-07 !*/
        {
            SEQ_TRACE_ERROR.print("runGradMoments for section 3 failed.");
            return false;
        }
        //---------------------------------------------------------------------------
        // advance reordering table to next contrast index
        //---------------------------------------------------------------------------
        // Note: This uses implicit knowledge about the reordering table layout.
        m_pRI->increaseReorderIndexOffset(m_pRI->getEchoTrainLength());
    }

    // Restore original status of readout events
    fRTSetReadoutEnable( bIsReadoutEnabled );

    // Restore original fill times
    m_lTRFill                = lTRFill;
    m_lRTEBPlugInTEFillAfter = lTEFill;

    //---------------------------------------------------------------------------
    // Re-enable test cases
    //---------------------------------------------------------------------------
    if (SeqUT.isUnitTestActive())
    {
        if ( m_bIsSliceAcceleration && ( m_eRunMode == MULTI_BAND ) )
        {
            SeqUT.EnableTestCase( lFrequNotEquOmegaRFErr, RTEB_ORIGIN_fSEQRunKernel );
            SeqUT.EnableTestCase( lFqNCOPhaseRoErr,       RTEB_ORIGIN_fSEQRunKernel );
        }
    }

    // --------------------------------------------------------------------------
    // reset local TE-fill
    // --------------------------------------------------------------------------
    m_lLocalTEFill = 0;
    m_lLocalTRFill = 0;

    // --------------------------------------------------------------------------
    // if this point is reached successfully
    // --------------------------------------------------------------------------
    return true;
}


//===============================================================================
//
// Function:    calcEffEchoSpacingAndBWPerPixelPE()
//
// Calculates: (1) An effective echo-spacing, which takes the number of
//                 interleaves and the PAT acceleration factor into account
//
//             (2) The bandwidth per pixel in the phase-encoding direction
//                 FOR THE RECONSTRUCTED IMAGE
//
//===============================================================================

bool SeqBuildBlockEPIKernel::calcEffEchoSpacingAndBWPerPixelPE(MrProt &rMrProt, long& lEffectiveEchoSpacing, double& dBandwidthPerPixelPE)
{

    // --------------------------------------------------------------------------
    // read parameters from protocol
    // --------------------------------------------------------------------------
    const bool   b2DInterpolation   = rMrProt.kSpace().Interpolation2D();
    const long   lPatAccFactor      = rMrProt.PAT().getlAccelFactPE();
    const double dFovRo             = rMrProt.sliceSeries().aFront().readoutFOV();
    const double dFovPc             = rMrProt.sliceSeries().aFront().phaseFOV();
    const double dBaseResolution    = rMrProt.kSpace().getlBaseResolution();
    const double dPhaseOversampling = rMrProt.phaseOversampling();


    SEQ_TRACE_DEBUG.print("b2DInterpolation = %d", b2DInterpolation);
    SEQ_TRACE_DEBUG.print("lPatAccFactor = %ld", lPatAccFactor);
    SEQ_TRACE_DEBUG.print("dFovRo = %f", dFovRo);
    SEQ_TRACE_DEBUG.print("dFovPc = %f", dFovPc);
    SEQ_TRACE_DEBUG.print("dBaseResolution = %f", dBaseResolution);
    SEQ_TRACE_DEBUG.print("dPhaseOversampling = %f", dPhaseOversampling);


    // --------------------------------------------------------------------------
    // check that number of interleaves has been set
    // --------------------------------------------------------------------------
    if (m_lNumInterleaves < 1)
    {
        SEQ_TRACE_ERROR.print("m_lNumInterleaves has not been set.");
        return false;
    }

    // --------------------------------------------------------------------------
    // make sure that a non-zero PAT acceleration factor is read from MrProt
    // --------------------------------------------------------------------------
    if (lPatAccFactor < 1)
    {
        SEQ_TRACE_ERROR.print("lPatAccFactor is less than 1.");
        return false;
    }

    // --------------------------------------------------------------------------
    // determine effective echo-spacing in us
    // --------------------------------------------------------------------------
    const double dEffectiveEchoSpacingUs = (double)getEchoSpacing() / (double)(m_lNumInterleaves * lPatAccFactor);

    lEffectiveEchoSpacing = long(dEffectiveEchoSpacingUs + 0.5);

    // --------------------------------------------------------------------------
    // determine bandwidth per pixel in phase-encoding direction
    // --------------------------------------------------------------------------

    // interpolation factor with value 1 (interpolation off) or 2 (interpolation on)

    const double dInterpolFactor = b2DInterpolation ? 2.0 : 1.0;

    // phase-oversampling factor with range 1 (no OS) to 2 (100% OS)

    const double dPhaseOsFactor = 1.0 + dPhaseOversampling;

    // number of pixels in image in phase-encoding direction

    const double dNumPixelsPE = dBaseResolution * dFovPc / dFovRo;

    // bandwidth per pixel (Hz) in phase-encoding direction

    const double dTemp = dNumPixelsPE * dPhaseOsFactor * dInterpolFactor * dEffectiveEchoSpacingUs;

    if (dTemp < DBL_EPSILON)
    {
        SEQ_TRACE_ERROR.print("Denominator too small in phase-encode bandwidth per pixel calculation.");
        SEQ_TRACE_DEBUG.print("dNumPixelsPE = %f", dNumPixelsPE);
        SEQ_TRACE_DEBUG.print("dPhaseOsFactor = %f", dPhaseOsFactor);
        SEQ_TRACE_DEBUG.print("dInterpolFactor = %f", dInterpolFactor);
        SEQ_TRACE_DEBUG.print("dEffectiveEchoSpacingUs = %f", dEffectiveEchoSpacingUs);

        return false;
    }

    dBandwidthPerPixelPE = 1.0e6 / dTemp;

    // --------------------------------------------------------------------------
    // finish
    // --------------------------------------------------------------------------
    return true;
}

//--------------------------------------------------------------------
//	Overloaded base class method
//  Set gradient performance:
//  - If no GPA balance model is used: call base class method
//  - If GPA model is activated: call base class method and change
//    maximum amplitude afterwards
void SeqBuildBlockEPIKernel::setDefaultGradientPerformance()
{

    // Call base class method
    SeqBuildBlockEPIReadOut::setDefaultGradientPerformance();

    if (getUseGPABalance())
    {
        // If the balance model is active, we allow gradient amplitudes up to the
        // maximum absolute value divided by sqrt(2). Basic assumption: 
        // gradient activity only on two axes (RO and PE).
        double dMaxMagnitude = SysProperties::getGradMaxAmplAbsolute() / sqrt(2.);

        if (m_dGradMaxAmpl > 0.)
        {
            // If m_dGradMaxAmpl has been set, we allow gradient amplitudes up to this
            // value (taking into account absolute limits). Basic assumptions:
            //  a) PE gradients do not overlap with RO flat top
            //  b) Slewrate of PE gradients is not higher than that of RO gradients
            // If these conditions apply, even slice rotation / tilting will never
            // yield a single axis gradient amplitude higher than this.
            //
            // Note: It is not recommended to use values above 90% of GradMaxAmplAbsolute - 
            //       the combination with high slewrates might lead to a GPA overshoot.
            dMaxMagnitude = std::min<double>(m_dGradMaxAmpl, std::min<double>(m_sSBBBalance.dGetMaximumAmplitude(), SysProperties::getGradMaxAmplAbsolute()));
    }

        double adMagnitudes[3] = { dMaxMagnitude, dMaxMagnitude, dMaxMagnitude };

        // Ensure that we don't make things worse ...
        if (getMaxMagnitude(SEQ::GRAD_FAST) < dMaxMagnitude)
        {
            setMaxMagnitudes(adMagnitudes, SBBEPIReadOut_GRAD_PERF_RO);
            setMaxMagnitudes(adMagnitudes, SBBEPIReadOut_GRAD_PERF_BLIPS);
        }

        // additional 'don't make it worse' setting if system has an ULTRAFAST gradient mode
        if (SysProperties::hasGradientMode(SEQ::GRAD_ULTRAFAST))
        {
            if (getMaxMagnitude(SEQ::GRAD_ULTRAFAST) < dMaxMagnitude)
            {
                setMaxMagnitude(SEQ::GRAD_ULTRAFAST, dMaxMagnitude, SBBEPIReadOut_GRAD_PERF_RO);
                setMaxMagnitude(SEQ::GRAD_ULTRAFAST, dMaxMagnitude, SBBEPIReadOut_GRAD_PERF_BLIPS);
            }
        }

}

    return;
}

//--------------------------------------------------------------------
//	Check whether balance allows to apply the readout train
//  at least once
//  Prerequisite: EPI readout has to be prepared
bool SeqBuildBlockEPIKernel::checkBalance()
{
    GPABalance::GPABalanceValue sBalanceOut;

    // Clear content
    m_sSBBBalance.Reset();

    if ((m_sSBBBalance.lGetStatus() != 0) || (!getUseGPABalance()))
    {
        // This should not happen => trace always
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "ERROR: Balance model initialisation failed!");
        return false;
    }

    // Frequency content of readout module [Hz]
    // Note: providing a 'typical worst case' (e.g. 1000Hz) instead of the 
    //       actual value might yield a more consistent UI behaviour ...
    const long lFreqAcq = 1000000 / (2 * m_lEchoSpacing);
    m_sSBBBalance.bSetFrequency(lFreqAcq);

    // Set up a 'worst case' series of kernel gradient events:
    // [ Readout - Excitation - Phase Correction ]
    // Note: plugins are handled separately 
    // Note: only the worst case axis (readout) is evaluated

    const double dGROAmp      = m_GRO.getAmplitude();
    const long   lGRORampUp   = m_GRO.getRampUpTime();
    const long   lGRORampDown = m_GRO.getRampDownTime();
    const long   lGROFlatTop  = m_GRO.getFlatTopTime();

    int      iI = 0;
    long     lTimeUS = 0;
    double   dSign = 1.;

    const long lEchoTrainDuration
        = m_aGradMoments[2].getTotalTime() + getDurationEPIReadOutPerRequest() + m_aGradMoments[3].getTotalTime();


    // Note: It is assumed that all echo trains get applied sequentially without
    //       any intermittent delays. For elaborate PlugIns (e.g. diffusion
    //       preparation), this assumption does not apply.
    for ( long lJ = 0; lJ < m_lNumberOfContrasts; ++lJ )
    {
        // Add readout gradients (plus one for initial dephasing)
        for ( long lI = 0; lI < m_lEchoTrainLength + 1; ++lI )
        {
            m_sSBBBalance.lAddGradient(lTimeUS, dSign * dGROAmp, lGRORampUp, lGRORampDown, lGROFlatTop, GPABALANCE_X_AXIS);

            // Alternate the sign of the readout gradient
            dSign *= -1.;

            // Time increment
            lTimeUS += m_lEchoSpacing;
        }
        lTimeUS += std::max( 0L, lEchoTrainDuration - ( m_lEchoTrainLength + 1 ) * m_lEchoSpacing );
    }

    // Add excitation (consider as pause)
    lTimeUS += m_pSBBExcite->getDurationPerRequest();

    if (m_bInternalPhaseCorrection)
    {

        // Add three phase correction scans (plus one for initial dephasing)
        for (iI = 0; iI < 4; ++iI)
        {
            m_sSBBBalance.lAddGradient(lTimeUS, dSign * dGROAmp, lGRORampUp, lGRORampDown, lGROFlatTop, GPABALANCE_X_AXIS);

            // Alternate the sign of the readout gradient
            dSign *= -1.;

            // Time increment
            lTimeUS += m_lEchoSpacing;
        }
    }
    else
    {
        m_sSBBBalance.lAddEmptyEvent(lTimeUS);
    }

    // Note: m_sBalanceIn contains no history if not explicitely set
    const long lResult = m_sSBBBalance.lCalc(GPABALANCE_X_AXIS, m_sBalanceIn, sBalanceOut);

    if (lResult != 0)
    {
        // This might happen within binary search => no trace
        // SEQ_TRACE_ERROR.print("Balance calculation failed at %li", lResult );
        return false;
    }

    // EPI readout can be applied at least once
    return true;
}


//--------------------------------------------------------------------
//	Calculate TR increment necessary to repeat the readout train
//  forever
//  Prerequisite: EPI readout has to be prepared
bool SeqBuildBlockEPIKernel::calcTRIncrement(long lTEFillTime)
{

    if ((m_sSBBBalance.lGetStatus() != 0) || (!getUseGPABalance()))
    {
        // This should not happen => trace always
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "ERROR: Balance model initialization failed!");
        return false;
    }

    GPABalance::GPABalanceValue sBalanceOut;

    // ------------------------------------------------------------------------
    // 1. Calculate TR increment required by the GPA
    // ------------------------------------------------------------------------

    // Calculate balance contribution of EPI readout
    // Note: m_sBalanceIn contains no history if not explicitly set
    const long lResult = m_sSBBBalance.lCalc(GPABALANCE_X_AXIS, m_sBalanceIn, sBalanceOut);

    if (lResult != 0)
    {
        // This should not happen, since ::checkBalance() has been called before
        // => trace always
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "ERROR: Balance calculation failed unexpectedly!");
        return false;
    }

    // Duration of the readout events
    const long lDuration = m_sSBBBalance.lGetDuration();

    // Calculate required TR increment
    long lTRIncrement = m_sSBBBalance.lCalcPause(sBalanceOut, lDuration);

    if (lTRIncrement < 0)
    {
        // This should not happen => trace always
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "ERROR: Pause calculation failed!");
        return false;
    }

    m_lTRIncrement = lTRIncrement;


    // ------------------------------------------------------------------------
    // 2. Calculate TR increment required by longterm thermal effects
    // ------------------------------------------------------------------------

    GPABalance::GrmsValue sBalanceRmsOut;

    // Calculate longterm RMS contribution of EPI readout
    if (!m_sSBBBalance.bCalcRMS(m_sBalanceRmsIn, sBalanceRmsOut))
    {
        // This should not happen => trace always
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "ERROR: Balance rms calculation failed!");
        return false;
    }

    // Calculate required TR increment
    lTRIncrement = m_sSBBBalance.lCalcPauseRMS(sBalanceRmsOut);

    if (lTRIncrement < 0)
    {
        // This should not happen => trace always
        setNLSStatus(MRI_SEQ_SEQU_ERROR, __FUNCTION__, "ERROR: Pause calculation failed!");
        return false;
    }

    if (lTRIncrement > m_lTRIncrement)
    {
        // Use worst case
        m_lTRIncrement = lTRIncrement;
    }

    // Consider TE fill time as pause
    m_lTRIncrement -= lTEFillTime;

    if (m_lTRIncrement < 0)
    {
        m_lTRIncrement = 0;
    }

    return true;
}

bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::updateSliceAdjLocalRFInfo(MrProt &rMrProt, /*< IMP: points to the protocol structure. */ SeqLim &rSeqLim, /*< IMP: points to the sequence limits structure. */ SeqExpo &rSeqExpo /*< IMP: points to the sequence exports structure */)
{
    // The kernel does not store any energy at all but always asks the incorporated SBBs for their energy
    // and passes this information ==> we can overload the method and do nothing here
    return true;
}

bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::prepPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
{
    // this implementation does not do anything 
    SEQ_TRACE_ALWAYS.print("SeqBuildBlockEPIKernel::prepPlugIn is not implemented");
    return false;
}

bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::initExcitation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
{
    // this implementation does not do anything but needs to be overloaded
    SEQ_TRACE_ALWAYS.print("SeqBuildBlockEPIKernel::initExcitation is not implemented");
    return false;
}



#ifdef ZOOM_2DRF
bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::prepZOOMitExcitation(MrProt &rMrProt, const SeqLim &rSeqLim)
{

    // initialize ZOOM_2DRF parameters
    if (!setZOOMitParameterFromUI(rMrProt, rSeqLim))
    {
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("setZOOMitParameterFromUI() failed.");
        }
        return false;
    }

    // select excitation SBB based on PTXTrajectoryType
    if (m_ExcType == MrProtocolData::PTXTrajectoryType_EPI_1D)
    {
        if (!configureSBB2DExc(rMrProt, rSeqLim))
        {
            SEQ_TRACE_ERROR.print("onfigureSBB2DExc() failed.");
            m_pSBBExcite = nullptr;
        }
        else
        {
            // set and return pointer to class
            m_pSBBExcite = &m_2DExcite;
        }
    }
    else if (m_ExcType == MrProtocolData::PTXTrajectoryType_ExternalIni)
    {
        // do not use OptPTXVolume
        m_vdOptPTXVolumeThickness.clear();
        m_vdOptPTXVolumeShift.clear();

        if (m_2DExciteIniMultiCh.getIniRFPulseChannels() > 1)
        {
            m_2DExciteIniMultiCh.setPhysDesignFile("PPDconfig_ep2d_ini"); // input file to PPD
            m_2DExciteIniMultiCh.setUseOwnEventBlock(false);          // ep2d requires open RTEB: done in Kernel

            // return pointer to class
            m_pSBBExcite = &m_2DExciteIniMultiCh;

        }
        else if (m_2DExciteIniOneCh.getIniRFPulseChannels() == 1)
        {
            m_2DExciteIniOneCh.setUseOwnEventBlock(false);            // ep2d requires open RTEB: done in Kernel

            // return pointer to class
            m_pSBBExcite = &m_2DExciteIniOneCh;

        }
        else
        {
            SEQ_TRACE_ERROR.print("wrong external INI pulse file");
            m_pSBBExcite = nullptr;
            return false;
        }
    }
    else
    {
        SEQ_TRACE_ERROR.print("wrong ExcType == <%d>", m_ExcType);
        m_pSBBExcite = nullptr;
        return false;
    }
    return true;
}

bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setZOOMitParameterFromUI(MrProt &rMrProt, const SeqLim &rSeqLim)
{

#ifdef ZOOM_EXTENDED
    // Enable transmit partial Fourier
    double dPulsePartialFourierFactor = 6. / 8.;
    m_Param.dPulsePartialFourierFactor = dPulsePartialFourierFactor;

    // Enable rotation with default value
    m_Param.dRotationAngleDeg = m_dRotationAngleDegDefault;

    // update parameters based on debug xml
    if (m_bUseDebugSettings)
    {
        // rotation angle
        double dRotationAngleDeg = SysProperties::ReadSeqSettingGeneral("EPI_ZOOMIT/RotationAngleDeg", m_dRotationAngleDegDefault, true);
        SEQ_TRACE_ALWAYS.print("DebugXML: dRotationAngleDeg = %f", dRotationAngleDeg);
        m_Param.dRotationAngleDeg = dRotationAngleDeg;

        // pTx partial Fourier factor
        double dPulsePartialFourierFactorDefault = 6. / 8.;
        dPulsePartialFourierFactor = SysProperties::ReadSeqSettingGeneral("EPI_ZOOMIT/PulsePartialFourierFactor", dPulsePartialFourierFactorDefault, true);
        SEQ_TRACE_ALWAYS.print("DebugXML: dPulsePartialFourierFactor = %f", dPulsePartialFourierFactor);
        m_Param.dPulsePartialFourierFactor = dPulsePartialFourierFactor;
    }
#endif // ZOOM_EXTENDED

    // determine some parameters from the *first* PTXRFPulse (index 0)
    const size_t uiCurVol = 0;
    long lPTXRFPulseSize = rMrProt.getsTXSPEC().getaPTXRFPulse().size();

    if (lPTXRFPulseSize > 0)
    {
        m_ExcType = rMrProt.getsTXSPEC().getaPTXRFPulse()[uiCurVol].getlTrajectoryType();
        m_B0CorrectionType = rMrProt.getsTXSPEC().getaPTXRFPulse()[uiCurVol].getlB0CorrectionType();
        m_Param.dPulseAcceleration = rMrProt.getsTXSPEC().getaPTXRFPulse()[uiCurVol].getdPulseAcceleration();
    }
    else  // set to default
    {
        m_ExcType = MrProtocolData::PTXTrajectoryType_EPI_1D;
        m_B0CorrectionType = MrProtocolData::PTXB0CorrectionType_Off;
        m_Param.dPulseAcceleration = 1.0;
    }

    // take phase oversampling into account
    const double dPhaseOversampling = rMrProt.phaseOversampling();
    const double dPhaseOsFactor     = 1.0 + dPhaseOversampling;

    double dReducedFoE        = rMrProt.sliceSeries().aFront().phaseFOV() * dPhaseOsFactor * m_Param.dRelExcitationSize;
    double dMinFoEOuterBorder = rMrProt.sliceSeries().aFront().phaseFOV() * dPhaseOsFactor;

    double dFoEOuterBorder    = -1; // [mm]  

#ifdef ZOOM_EXTENDED
    // for rotated ZOOMit, pTx volume no longer available/applicable
    // therefore, set default value
    dFoEOuterBorder = RF2DArbGenerator::DefaultSideLobeDistance;
#else
    // retrieve FoEOuterBorder from the PhaseFOV of the *first* PTXVolume with the property 'Optimization'
    long lPTXVolumeSize = rMrProt.getsPTXData().getasPTXVolume().size();

    if (lPTXVolumeSize > 0)
    {
        for (long iVol = 0; ((iVol < lPTXVolumeSize) && (dFoEOuterBorder < 0)); iVol++)
        {
            if (rMrProt.getsPTXData().getasPTXVolume()[iVol].getlVolProperty() == MrProtocolData::PTxVolProp_Optimization)
            {
                dFoEOuterBorder = rMrProt.getsPTXData().getasPTXVolume()[iVol].getsSliceData().getdPhaseFOV();  // [mm]
            }
        }
    }

    if (dFoEOuterBorder < 0)
    {
        dFoEOuterBorder = RF2DArbGenerator::DefaultSideLobeDistance;
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ALWAYS.print("No Optimization PTXVolume found. Using default value for FoEOuterBorder = %.2f.", dFoEOuterBorder);
        }
    }
#endif

    if ((dFoEOuterBorder < dMinFoEOuterBorder))
    {
#ifdef ZOOM_EXTENDED
        // this check is not applicable for rotated ZOOMit since dFoEOuterBorder is only a dummy value
        // therefore, patch it to dMinFoEOuterBorder
        dFoEOuterBorder = dMinFoEOuterBorder;
#else
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ALWAYS.print("PTXVolume phaseFOV (dFoEOuterBorder) = %.2f is too small (< %.2f).", dFoEOuterBorder, dMinFoEOuterBorder);
        }
        return false;
#endif
    }

    if (m_Param.dXReadBWT < 0.1)   // check for DivByZero
    {
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ALWAYS.print("dXReadBWT (%.2f) is too small .", m_Param.dXReadBWT);
        }
        return false;
    }

    // set length of slice selection gradient   
    //     empirical slice thickness correction
    // depending on define flags COMPILE_EP2D_FID or COMPILE_EP2D_DIFF
    double dSliceThkFactor = 1.0;
#if defined COMPILE_EP2D_FID  
    dSliceThkFactor = 1.1;
#elif defined COMPILE_EP2D_DIFF
    dSliceThkFactor = 1.0;
#endif
    m_Param.dXReadResolution = rMrProt.sliceSeries().aFront().thickness() * 2 / m_Param.dXReadBWT / dSliceThkFactor;

    // calculate field-of-excitation FoE (slidelobe distance)
    // in general: FOE = width_passband + width_stopband = 1/area_blip_phase
    // in our implementation:
    //    "Reduced FoE" is the width_passband
    //    "FoE outer border" taken from GSP/UI parameter, is the distance from the left to the right side of suppression (includes the passband)
    //    "XPhaseResolution": the edges of the pulse is included to guarantee better FoE definition on both sides
    //    "SidelobeDistance" is the trajectory design parameter, the effective FoE, determines the area_blip_phase
    //
    //    note: the factor 0.5 results from the extension of the FoEOuterBorder to *both* sides of the passband:
    //          FoEOuterBorder = edge + width_stopband + width_passband + width_stopband + edge:
    m_Param.lSidelobeDistance = static_cast<long>(0.5*(dFoEOuterBorder + dReducedFoE + 2 * m_Param.dXPhaseResolution) + 0.5);

#ifdef ZOOM_EXTENDED
    // Overwrite default sidelobe distance by automatically optimized sidelobe distance
    bool bOptimizeSidelobeDistance = true;

    // update parameters based on debug xml
    if (m_bUseDebugSettings)
    {
        // optimize sidelobe distance
        bOptimizeSidelobeDistance = SysProperties::ReadSeqSettingGeneral("EPI_ZOOMIT/OptimizeSideLobeDistance", true, true);
        SEQ_TRACE_ALWAYS.print("DebugXML: bOptimizeSidelobeDistance = %d", bOptimizeSidelobeDistance);
    }

    if (bOptimizeSidelobeDistance)
    {
        // Rotating the field of excitation allows automatically choosing the minimal possible sidelobe distance while still avoiding infolding artifacts

        // However, The final rotation angle depends on the realized excitation trajectory, which in turn depends on the sidelobe distance
        // To solve this circular dependency, several prepares would have to be required, which is not the case on the imager
        // Therfore, a worst case scenario for the sidelobe distance is assumed and used to calculate the excitation trajectory
        // Later on, the rotation angle is calculated based on this information

        // For this worst case scenario, the sidelobe distance is defined as the sidelobedistance which would be required for m_dRotationAngleDegLowerLimit and m_dRotationAngleDegUpperLimit

        const double sinAngleLowerLimit = sin(m_dRotationAngleDegLowerLimit * M_PI / 180.);
        const double cosAngleLowerLimit = cos(m_dRotationAngleDegLowerLimit * M_PI / 180.);

        const double sinAngleUpperLimit = sin(m_dRotationAngleDegUpperLimit * M_PI / 180.);
        const double cosAngleUpperLimit = cos(m_dRotationAngleDegUpperLimit * M_PI / 180.);

        const double dProfileDimFast_mm = rMrProt.sliceSeries().aFront().phaseFOV() * dPhaseOsFactor; // [mm]
        const double dProfileDimSlow_mm = rMrProt.sliceSeries().aFront().thickness(); // [mm]
        const double dTransitionRegion_mm = m_Param.dXPhaseResolution; // [mm]

        // calculate sidelobeDistance1
        // required sidelobe distance to avoid infolding (overlap of side excitations with imaging plane)
        const double sidelobeDistance1_RotAngleLowerLimit = dProfileDimSlow_mm / sinAngleLowerLimit;
        const double sidelobeDistance1_RotAngleUpperLimit = dProfileDimSlow_mm / sinAngleUpperLimit;

        // calculate sidelobeDistance2
        // required sidelobe distance to avoid interference of side excitations with slice stack
        const double sidelobeDistance2_RotAngleLowerLimit = (dProfileDimFast_mm + dTransitionRegion_mm) / cosAngleLowerLimit;
        const double sidelobeDistance2_RotAngleUpperLimit = (dProfileDimFast_mm + dTransitionRegion_mm) / cosAngleUpperLimit;

        // for both rotation angles, get the maximum sidelobe distance
        const double sidelobeDistanceRotAngleLowerLimit = std::max(sidelobeDistance1_RotAngleLowerLimit, sidelobeDistance2_RotAngleLowerLimit);
        const double sidelobeDistanceRotAngleUpperLimit = std::max(sidelobeDistance1_RotAngleUpperLimit, sidelobeDistance2_RotAngleUpperLimit);

        // optimized sidelobe distance is maximum of sidelobeDistanceRotAngleLowerLimit and sidelobeDistanceRotAngleUpperLimit
        m_Param.lSidelobeDistance = static_cast<long>(std::max(sidelobeDistanceRotAngleLowerLimit, sidelobeDistanceRotAngleUpperLimit));
    }
#endif // ZOOM_EXTENDED

    // To avoid infolding artifacts of the TX ghost (occur at half sidelobe distance) we double the effective FoE
    //    NOTE: this can be looked at as a PSEUDO flyback trajectory mode. 
    //          Disadvantage is the double length of TX trajectory, i.e. double pulse duration, double off-resonance artifacts
    m_Param.lSidelobeDistance *= 2;


#ifdef ZOOM_EXTENDED
    // Update dFoEOuterBorder based on automatically calculated sidelobe distance
    dFoEOuterBorder = 2 * static_cast<double>(m_Param.lSidelobeDistance) - 2 * m_Param.dXPhaseResolution - dReducedFoE;

    // repeat check from above
    if ((dFoEOuterBorder < dMinFoEOuterBorder))
    {
        SEQ_TRACE_ALWAYS.print("PTXVolume phaseFOV (dFoEOuterBorder) = %.2f is too small (< %.2f).", dFoEOuterBorder, dMinFoEOuterBorder);
        dFoEOuterBorder = dMinFoEOuterBorder;
    }
#endif // ZOOM_EXTENDED

    // store rel. size parameter
    m_Param.dFoERelativeSize = dReducedFoE / dFoEOuterBorder;

    // update (dependent) PTXPulse fields of pulseID 0 - for display only
    if (rMrProt.getsTXSPEC().getaPTXRFPulse().size() > 0)
    {
        rMrProt.getsTXSPEC().getaPTXRFPulse()[0].setdFlipAngleDegree(static_cast<float>(rMrProt.flipAngle()));
    }

    // calculate OptPTXVolume parameters
    m_vdOptPTXVolumeThickness.clear();
    m_vdOptPTXVolumeShift.clear();

    // Target definition of OptVolume thickness and shift:  cover the OptPTXVolume up to the outer border
    double dOptPTXVolumeThickness = (dFoEOuterBorder - dReducedFoE) * 0.5;
    double dOptPTXVolumeShift = (dFoEOuterBorder + dReducedFoE) * 0.25;

    // DISABLED: shift in addition by the edge width (dXPhaseResolution) not to touch the inner volume excitation
    ////dOptPTXVolumeShift += dXPhaseResolution;

    bool bEnableOptPTXVolume = true;

#ifdef ZOOM_EXTENDED
    // With rotated ZOOMit, these hidden saturation pulses are not necessary anymore
    bEnableOptPTXVolume = false;
#endif // ZOOM_EXTENDED


    // OptPTXVolume is activated for TX Acceleration (all sequences)
    // For EP2D_DIFF: OptPTXVolume is linked to the UI parameter MR_TAG_ADJ_VOL_COUPLE_TO
    if (bEnableOptPTXVolume)
    {
        if (isOptPTXVolumeCondition(rMrProt))
        {
            // possibly disable OptVolume through registry flag: PARAM_EPI_Disable_OptPTXVolume
            long lRegistryKeyValue = getMaskFromRegistry("PARAM_EPI_Disable_OptPTXVolume", false);
            if (lRegistryKeyValue == 0)
            {
                m_vdOptPTXVolumeThickness.push_back(dOptPTXVolumeThickness);
                m_vdOptPTXVolumeShift.push_back(+1 * dOptPTXVolumeShift);   // positive direction

                m_vdOptPTXVolumeThickness.push_back(dOptPTXVolumeThickness);
                m_vdOptPTXVolumeShift.push_back(-1 * dOptPTXVolumeShift);   // reverse direction
            }
        }
    }

    return true;
}

bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::configureSBB2DExc(MrProt &rMrProt, const SeqLim &)
{

    // settings for SBBExcitation
    m_2DExcite.setFlipAngle(rMrProt.flipAngle());
    m_2DExcite.setUseOwnEventBlock(false);                // ep2d requires open RTEB: done in Kernel
    m_2DExcite.setStartTimeInEventBlock(0);
    m_2DExcite.setOptPTXVolume(m_vdOptPTXVolumeThickness, m_vdOptPTXVolumeShift);

    // struct to configure RF2DArbGenerator
    m_Param.dFlipAngle = rMrProt.flipAngle();

    // transfer parameter struct to RF2DArbGenerator
    if (!m_2DExcite.setRFParams_Epi(m_Param))
    {
        SEQ_TRACE_ALWAYS.print("m_2DExcite.setRFParams_Epi() failed.");
        return false;
    }

    // set type of B0 correction for PTX pulse calculation
    B0CorrMode eB0Corr;
    switch (m_B0CorrectionType)
    {
    case MrProtocolData::PTXB0CorrectionType_Standard:
        eB0Corr = B0CORR_STD;
        break;
    case MrProtocolData::PTXB0CorrectionType_Advanced:
        eB0Corr = B0CORR_ADV;
        break;
    case MrProtocolData::PTXB0CorrectionType_Off:
    default:
        eB0Corr = B0CORR_OFF;
    }

    bool bErrorWithQCCheck = true;

    m_2DExcite.setPhysParams_Epi(
        m_Param.dXPhaseResolution,                          // dFoEResolution
        m_Param.dFoERelativeSize,                           // dFoERelativeSize: relative to FoE outer border
        static_cast<double>(m_Param.lSidelobeDistance),   // dFoE
        m_Param.dFlipAngle,                                 // dFlipAngle
        eB0Corr,                                            // B0 correction
        bErrorWithQCCheck);                      // sequence error, if QC control is not passed

    return true;
}

#ifdef ZOOM_EXTENDED
bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::updateZOOMitRotation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
{
    const double dEps = std::numeric_limits<double>::epsilon();

    // Overwrite default rotation distance by automatically optimized sidelobe distance
    bool bOptimizeRotationAngle = true;

    // update parameters based on debug xml
    if (m_bUseDebugSettings)
    {
        // optimize sidelobe distance
        bOptimizeRotationAngle = SysProperties::ReadSeqSettingGeneral("EPI_ZOOMIT/OptimizeRotationAngle", true, true);
        SEQ_TRACE_ALWAYS.print("DebugXML: bOptimizeRotationAngle = %d", bOptimizeRotationAngle);
    }

    if (bOptimizeRotationAngle)
    {
        // In order to avoid infolding artifacts from side excitations, the EPI excitation trajectory is rotated.
        // However, in the offresonant case, the rotation introduces a shift of offresonant signal components along the slice direction and therefore 
        // reduces the overlap of these shifted components with the refocussing pulse.
        // As a consequence, a certain B0 sensitivity is introduced. 
        // While this B0 sensitivity improves fat suppression and can be beneficial (to some extent), it can also lead to unwanted signal decrease in offresonant regions.
        // The shift of offresonances depends on the parameters of the 2D excitation pulse, and is for example aggravated for thin slices.

        // The exact shift can be calculated following Eq. [9] in "Alley_1997_Angiographic Imaging with 2D RF Pulses":
        // Shift = DeltaFrequency * Nblip * DurationLine * FOVphase / BWTP

        // In order to achieve a constant shift (and therefore a constant B0 sensitivity) across different protocol settings such as slice thickness, etc.,
        // a logic is implemented to automatically calculate the rotation angle.

        // In particular, the desired spatial shift of signal with an offresonance of dFrequencyShiftHz to be dAllowedRelativeShift times the slice thickness is defined
        // and the corresponding rotation angle to achieve such shift is calculated.

        double dOptimizedRotationAngleRad = m_dRotationAngleDegDefault / 180. * M_PI;
        const double dSliceThickness = rMrProt.sliceSeries().aFront().thickness(); // [mm]

        // To generate the same effect in offresonance sensitivity, the desired frequency shift should be field strength dependent
        const double dFrequencyShiftHz = 35. * SysProperties::getNominalB0(); // [Hz], empirical value

        const double dAllowedRelativeShift = 4.; // Empirically set value to limit offresonant shift [units of slice thickness]
        const double dAllowedAbsoluteShift = dAllowedRelativeShift * dSliceThickness; // [mm]

        if (RF2DArbGenerator* pRF2DArbGenerator = m_2DExcite.getRFArbGenerator())
        {
            if (C2DGradientsBlippedEPI* pC2DGradientsBlippedEPI = dynamic_cast<C2DGradientsBlippedEPI*>(pRF2DArbGenerator->getTrajectory()))
            {
                sGRAD_PULSE sLobe_fast = pC2DGradientsBlippedEPI->getXLine();
                const int iTrajectoryDuration_us = pC2DGradientsBlippedEPI->getTrajectoryDuration();
                const int iDuration_lobe_fast_us = sLobe_fast.getTotalTime();
                const int iNumberOfKLines = iDuration_lobe_fast_us == 0 ? 0 : iTrajectoryDuration_us / iDuration_lobe_fast_us;
                const int iNBlip = iNumberOfKLines - 1;
                const double dBWTP_slow = m_Param.dXReadBWT;

                // Calculate offresonant shift along line direction, according to Eq. [9] in "Alley_1997_Angiographic Imaging with 2D RF Pulses"
                // Shift = DeltaFrequency * Nblip * DurationLine * FOVphase / BWTP
                // In the rotated case: Shift = DeltaFrequency * Nblip * DurationLine * FOVphase / BWTP * sin(rotationAngle)
                const double dPhaseOversampling = rMrProt.phaseOversampling();
                const double dPhaseOsFactor     = 1.0 + dPhaseOversampling;

                const double dPhaseFOV    = rMrProt.sliceSeries().aFront().phaseFOV() * dPhaseOsFactor;
                const double dNominator = dAllowedAbsoluteShift * dBWTP_slow;
                const double dDenominator = iNBlip * dFrequencyShiftHz * (1e-6 * iDuration_lobe_fast_us) * dPhaseFOV;

                double dSinOptimizedRotationAngle = 0;
                if (dDenominator > dEps)
                {
                    dSinOptimizedRotationAngle = dNominator / dDenominator;
                }
                // check for valid value range of asin
                if (dSinOptimizedRotationAngle >= -1 && dSinOptimizedRotationAngle <= 1)
                {
                    dOptimizedRotationAngleRad = asin(dSinOptimizedRotationAngle);
                }

                // rotation angle should not exceed or fall below certain limits
                const double dRotationUpperLimit = m_dRotationAngleDegUpperLimit * M_PI / 180.;
                const double dRotationLowerLimit = m_dRotationAngleDegLowerLimit * M_PI / 180.;

                dOptimizedRotationAngleRad = std::min(dOptimizedRotationAngleRad, dRotationUpperLimit);
                dOptimizedRotationAngleRad = std::max(dOptimizedRotationAngleRad, dRotationLowerLimit);

                const double dOptimizedRotationAngleDeg = dOptimizedRotationAngleRad * 180. / M_PI;

                // update parameter struct with the optimized rotation angle
                m_Param.dRotationAngleDeg = dOptimizedRotationAngleDeg;

                // transfer updated parameter struct to RF2DArbGenerator
                if (!m_2DExcite.setRFParams_Epi(m_Param))
                {
                    if (!rSeqLim.isContextPrepForBinarySearch())
                    {
                        SEQ_TRACE_ALWAYS.print("m_2DExcite.setRFParams_Epi() failed.");
                    }
                    return false;
                }

                // reprepare everything that updated rotation angle is played out
                if (!m_pSBBExcite->prep(rMrProt, rSeqLim, rSeqExpo))
                {
                    // Preparation of excitation module preparation might fail in binary search
                    if (!rSeqLim.isContextPrepForBinarySearch())
                    {
                        setNLSStatus(m_pSBBExcite->getNLSStatus(), __FUNCTION__, "m_pSBBExcite->prep failed.");
                    }
                    return false;
                }
            }
        }
    }

    return true;
}
#endif

bool SeqBuildBlockEPIKernel::isOptPTXVolumeCondition(MrProt& /*rMrProt*/)
{
    return ((m_ExcType == MrProtocolData::PTXTrajectoryType_EPI_1D) && (m_Param.dPulseAcceleration > 1.0));
}

#endif //ZOOM_2DRF


//===============================================================================
//
// Function:    prepAdditionalGradients()
//
// Description: Prepares additional gradients, e.g. flow compensated gradients.
//               
// Return: true, real implementation is in SBBEPIFCKernel
//===============================================================================

bool SeqBuildBlockEPIKernel::prepAdditionalGradients(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
{
    return true;
}

//===============================================================================
//
// Function:    adjustTEContributions()
//
// Description: Adjust TE contributions because of additional gradients
// e.g. flow compensated gradients.
//               
// Return: true, real implementation is in SBBEPIFCKernel
//===============================================================================

bool SeqBuildBlockEPIKernel::adjustTEContributions(MrProt &rMrProt)
{
    return true;
}

//===============================================================================
//
// Function:    adjustSBBDurationPerRequest()
//
// Description: Adjust additional timing config because of additional gradients
// e.g. flow compensated gradients.
//               
// Return: true, real implementation is in SBBEPIFCKernel
//===============================================================================
bool SeqBuildBlockEPIKernel::adjustAdditionalTiming(MrProt &rMrProt)
{
    return true;
}


//===============================================================================
//
// Function:    adjustLocalTRFill()
//
// Description: Adjust local TR Fill because of additional gradients
// e.g. flow compensated gradients.
//               
// Return: true, real implementation is in SBBEPIFCKernel
//===============================================================================
bool SeqBuildBlockEPIKernel::adjustLocalTRFill(MrProt &rMrProt)
{
    return true;
}

//===============================================================================
//
// Function:    runAdditionalGradExcitation()
//
// Description: run additional gradients in Excitation event block
// e.g. flow compensated slice refocusing gradient  
//               
// Return: true, real implementation is in SBBEPIFCKernel
//===============================================================================
bool SeqBuildBlockEPIKernel::runAdditionalGradExcitation(long &lEventBlockEnd)
{
    return true;
}
//===============================================================================
//
// Function:    updateAdditionalGradMoments()
//
// Description: updated gradient moments for flow compensated gradients.
//               
// Return: true, real implementation is in SBBEPIFCKernel
//===============================================================================
bool SeqBuildBlockEPIKernel::updateAdditionalGradMoments(long lSectionIndex, MrProt &rMrProt, SeqLim &rSeqLim, long &lFillTimeBeforeGradients, long &lEventBlockEnd)
{
    return true;
}

//===============================================================================
//
// Function:    runFinalGradMoments()
//
// Description: run gradients.
//               
// Return: 
//===============================================================================
bool SeqBuildBlockEPIKernel::runFinalGradMoments(long lSectionIndex, MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC, long lFillTimeBeforeGradients, long lEventBlockEnd)
{
    //---------------------------------------------------------------------------
    // switch them
    //---------------------------------------------------------------------------
    fRTEBInit(pSLC->getROT_MATRIX()); /*! EGA-07 !*/

    fRTEI(lFillTimeBeforeGradients, 0, 0, 0, &m_sGrad[m_lGP], &m_sGrad[m_lGR], &m_sGrad[m_lGS], 0);          /*! EGA-All !*/
    fRTEI(lEventBlockEnd, 0, 0, 0, 0, 0, 0, 0);          /*! EGA-All !*/


    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 01, 0, pSLC->getSliceIndex(), 0, 0); /*! EGA-All !*/

    const NLS_STATUS lStatus = fRTEBFinish();

    if (setNLSStatus(lStatus, __FUNCTION__, "fRTEBFinish for dMoment[lI] failed!"))
    {
        // lStatus is error code, ergo:
        //
        return false;
    }

    return true;
}

void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::increaseSliceSelectionGradientMaxAmplitude(MrProt& rMrProt)
{
    // functionality only in derived Diffusion kernel
}
