//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 2010  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\SBBDiffusion_Trace.cpp
//
//      Author: PLM AW NEUR
//
//        Lang: C++
//
//     Descrip: Implementation of the class Diffusion_Trace.
//
//     Classes: Diffusion_Trace
//
//    -----------------------------------------------------------------------------
//
// EGA Requirement Key: As shown on the following lines:
//
//   Abbrev.   Translation                                        Relevant for
//   -------   -----------                                        ------------
//   EGA-05    {:IMPLEMENT:000_EGA_BildPos_SW_NCOFrequenzSSel::}  SRF     frequency
//
//    -----------------------------------------------------------------------------


/**
***************************************************************************

\changed     1-Sep-2002; M.Zwanger; 4a21a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- ancient monster file split up in several files
- class 'Diffusion_Stejskal' added
- comments changed to doxygen format (http://www.doxygen.org)
- some .h includes removed
- removed member m_Images
- SBBDiffusion_Base::prepSpoilGrad() - Spoiler is now specified in uTs
- m_lSpoilerDuration renamed in m_lSpoilerTotalTime
- m_Norm removed
- Spoilers moved next to the RF pulses in Diffusion_Trace::run()

***************************************************************************

\changed    31-Okt-2003; M.Zwanger; 4a21a; CHARM: n.a.
\description
- 'm_bDiffusionGradientsEnabled' support added
- include paths adapted for VA25A archive

***************************************************************************
*/

#include "MrProtSrv/Domain/MrProtData/MrProt/Application/Application.h"
#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Trace.h"
#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"
#include "MrImagingFW/libSeqUTIF/libsequt.h"            // mSEQTest
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"
#include "DiffusionRFPulseProperties.h"

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif
using namespace SEQ_NAMESPACE;



// ***************************************************************************
// class Diffusion_Trace
// ***************************************************************************



// ===========================================================================
///  The constructor initializes the starting time of the diffusion gradients 
///  with 0 and the maximum possible gradient amplitudes.
SBBDiffusion_Trace::SBBDiffusion_Trace(SBBList* pSBBList):
SBBDiffusion_Base(pSBBList),
m_lPreFill(0),
m_lPostFill(0),
m_bRunOnlyOnce(false),
m_bDumpBMatrix(false),
m_dPreparedSlcThk(0.0),
m_DRF("RTEIdentDRF"),
m_DFPset("RTEIdentDFPset"),
m_DFPneg("RTEIdentDFPneg"),
m_DGSS("RTEIdentDGSS")
// ===========================================================================
{
    setIdent("Diffusion_Trace");

    seteSliceAdjOptimizationMode(SLICEADJ::HOLD_OPTIMIZATION);  // the diffusion SBB shall only hold the results from the pulses
    setsSliceAdjParametersRequestedBySBB(SLICEADJ::ADJALL); //  SliceAdj functionality is in this SBB

    int  iI        = 0;
    char istring[] = "0000";
    char tRTEIdent[32];

    // Set maximum possible(!) gradient amplitudes (used during check())
    double dAmpl = SysProperties::getGradMaxAmplAbsolute();

    // Is there an elegant way to initialize member arrays???
    int iPolarityPhase[16] ={-1, +1, +1, -1, -1, -1, -1, -1, +1, +1, -1, -1, -1, -1, +1, +1};
    int iPolarityRead[16] ={-1, -1, +1, +1, -1, +1, +1, -1, -1, -1, +1, +1, -1, -1, -1, -1};
    int iPolaritySlice[16] ={-1, -1, -1, -1, +1, +1, +1, +1, +1, +1, +1, +1, -1, +1, +1, -1};

    for(iI = 0; iI < 16; iI++)
    {
        m_PhaseGradSign[iI] = iPolarityPhase[iI];
        m_ReadGradSign[iI] = iPolarityRead[iI];
        m_SliceGradSign[iI] = iPolaritySlice[iI];

        m_DGP[iI].setStartTime(0);

        m_DGP[iI].setMaxMagnitude(dAmpl);
        m_DGR[iI].setMaxMagnitude(dAmpl);
        m_DGS[iI].setMaxMagnitude(dAmpl);

        // Set identifier
        sprintf(istring, "%02d", iI);
        strcpy(tRTEIdent, "RTEIdentDGP");
        strcat(tRTEIdent, istring);
        m_DGP[iI].setIdent(tRTEIdent);

        strcpy(tRTEIdent, "RTEIdentDGR");
        strcat(tRTEIdent, istring);
        m_DGR[iI].setIdent(tRTEIdent);

        strcpy(tRTEIdent, "RTEIdentDGS");
        strcat(tRTEIdent, istring);
        m_DGS[iI].setIdent(tRTEIdent);
    }
}


// ===========================================================================
///   This destructor does nothing.
SBBDiffusion_Trace::~SBBDiffusion_Trace()
// ===========================================================================
{
}


bool SBBDiffusion_Trace::calcSliceAdjSBBRFInfo(
    MrProt                                     &rMrProt,        //< IMP: points to the protocol structure.
    SeqLim                                     &rSeqLim,        //< IMP: points to the sequence limits structure.
    SeqExpo                                    &rSeqExpo,       //< IMP: points to the sequence exports structure
    const SLICEADJ::sCuboidGeometry              &sSliceAdjCuboid,  //< IMP: cuboid geometry (single)
    std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo        //< EXP: RF info
    )
{
    // Initialize export
    vsRFInfo.clear();

    if(! SLICEADJ::calcRFInfoForSliceAdjCuboid(rMrProt, rSeqLim, rSeqExpo, &m_DRF, sSliceAdjCuboid, getsSliceAdjParametersRelevantForUpdate(), vsRFInfo))
    {
        MRTRACE("SLICEADJ::calcRFInfoForSliceAdjCuboid() failed");
        setNLSStatus(MRI_SBB_SBB_ERROR);
        return false;
    }

    return true;
}



// ===========================================================================
/// Implementation of the pure virtual base class method
// ===========================================================================
bool SBBDiffusion_Trace::prepInit(MrProt & /* &rMrProt */, SeqLim & /* &rSeqLim */, SeqExpo & /* &rSeqExpo */)
// ===========================================================================
{
    // -------------------
    // Set Gradient limits
    // -------------------
    // Set max gradient Amplitude
    // "Trace" has a more critical duty cycle behaviour than the other modes
    setMaxAmplitude(SysProperties::getGradMaxAmpl(SEQ::GRAD_FAST));

    // -------------------------
    // Prepare spoiler gradients
    // -------------------------
    // In this single shot trace mode, spoilers are ONLY applied for b=0.
    // But the sequence timing is governed by the gradients for b>0.
    // There is a guarantee that there will be enough time to play out these spoilers.
    // Therefore the preparation of these gradients will not take place
    // in prepTiming(), but here.
    if(! prepSpoilGrad(-0.70 * m_dAmpl * static_cast<double> (m_lRampTime)* 7. / 1000.))
    {
        SEQ_TRACE_ERROR.print("ERROR: Preparation of spoiler gradients failed.");
        return false;
    }

    return true;
}

// ===========================================================================
/// Implementation of the pure virtual base class method
// ===========================================================================
bool SBBDiffusion_Trace::prepFinal(double dMaxRequestedBValue, bool bIsContextForBinarySeach)
// ===========================================================================
{
    // ---------------------------
    // Prepare diffusion gradients
    // ---------------------------
    // Scale Amplitude exactly to requested maximum b value
    double dAmpl = m_dAmpl * sqrt(dMaxRequestedBValue / m_dMaxPossibleBValue);

    for(int iI = 0; iI < 16; iI++)
    {
        m_DGP[iI].setAmplitude(dAmpl * m_PhaseGradSign[iI]);
        m_DGP[iI].setDuration(m_DGP[0].getDuration());
        m_DGP[iI].setAxis(SEQ::AXIS_PHASE);
        m_DGP[iI].setRampTimes(m_lRampTime);

        if(! m_DGP[iI].prep())
        {
            if(!bIsContextForBinarySeach)
            {
                SEQ_TRACE_ERROR.print("ERROR: m_DGP[%d].prep() returned false ", iI);
            }
            setNLSStatus(m_DGP[iI].getNLSStatus());
            return false;
        }

        if(! m_DGP[iI].check())
        {
            if(!bIsContextForBinarySeach)
            {
                SEQ_TRACE_ERROR.print("ERROR: m_DGP[%d].check() returned false ", iI);
            }
            setNLSStatus(m_DGP[iI].getNLSStatus());
            return false;
        }

        m_DGR[iI].setAmplitude(dAmpl * m_ReadGradSign[iI]);
        m_DGR[iI].setDuration(m_DGP[0].getDuration());
        m_DGR[iI].setAxis(SEQ::AXIS_READOUT);
        m_DGR[iI].setRampTimes(m_lRampTime);
        m_DGR[iI].setStartTime(m_DGP[iI].getStartTime());

        if(! m_DGR[iI].prep())
        {
            if(!bIsContextForBinarySeach)
            {
                SEQ_TRACE_ERROR.print("ERROR: m_DGR[%d].prep() returned false ", iI);
            }
            setNLSStatus(m_DGR[iI].getNLSStatus());
            return false;
        }

        if(! m_DGR[iI].check())
        {
            if(!bIsContextForBinarySeach)
            {
                SEQ_TRACE_ERROR.print("ERROR: m_DGR[%d].check() returned false ", iI);
            }
            setNLSStatus(m_DGR[iI].getNLSStatus());
            return false;
        }

        m_DGS[iI].setAmplitude(dAmpl * m_SliceGradSign[iI]);
        m_DGS[iI].setDuration(m_DGP[0].getDuration());
        m_DGS[iI].setAxis(SEQ::AXIS_SLICE);
        m_DGS[iI].setRampTimes(m_lRampTime);
        m_DGS[iI].setStartTime(m_DGP[iI].getStartTime());

        if(! m_DGS[iI].prep())
        {
            if(!bIsContextForBinarySeach)
            {
                SEQ_TRACE_ERROR.print("ERROR: m_DGS[%d].prep() returned false ", iI);
            }
            setNLSStatus(m_DGS[iI].getNLSStatus());
            return false;
        }

        if(! m_DGS[iI].check())
        {
            if(!bIsContextForBinarySeach)
            {
                SEQ_TRACE_ERROR.print("ERROR: m_DGS[%d].check() returned false ", iI);
            }
            setNLSStatus(m_DGS[iI].getNLSStatus());
            return false;
        }

    }


    //----------------------------
    // Calculate Export Parameters
    //----------------------------
    // Store global RF info
    m_RFInfoPerRequest = m_DRF.getRFInfo();

    // The times must be calculated from the events itself
    // (MrProt->TE must not be used, as SBBDurationPerRequest contains the time *requested* 
    // to run the SBB for all b values. (Sorry, don't blame me, I didn't design this)
    m_lPreFill  =  m_lActualTE/2 - m_lSpinPrepTimeus                   - m_DGSS.getTotalTime()/2 - 10 * m_DGP[0].getTotalTime();
    m_lPostFill =  m_lActualTE/2 - m_lADCusTillEcho  - m_lStimoDelayus - m_DGSS.getTotalTime()/2 -  6 * m_DGP[0].getTotalTime();

    // in case of te()[m_iTEArrayIndex] on single gradient raster or even worse not on raster, time before plus
    // time after te()[m_iTEArrayIndex]/2 available for diffsbb will differ from the primarily offered time
    // due to rounding error if we do a fSDSRoundUpGRT for both  m_lPreEchoTimeContrib and
    // m_lPostEchoTimeContrib. To ensure that both times are ON the raster and to avoid the
    // above mentioned differences we round one time UP and one time DOWN!
    // Additionally m_lADCusTillEcho may also be not on gradient raster.
    // Note: a_ep2d_diff.h does neither use m_lPreEchoTimeContrib nor m_lPostEchoTimeContrib for
    // internal calculations.

    m_lPreEchoTimeContrib =  m_DGSS.getTotalTime()/2 + 10 * m_DGP[0].getTotalTime() + m_lPreFill;
    m_lPostEchoTimeContrib = m_DGSS.getTotalTime()/2 +  6 * m_DGP[0].getTotalTime() + m_lPostFill + m_lStimoDelayus;

    setSBBDurationPerRequest(m_lActualTE - m_lADCusTillEcho - m_lSpinPrepTimeus);

    m_bIsMagnetizationInverted = true;

    return true;
}

// ===========================================================================
bool SBBDiffusion_Trace::runSBB(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
// ===========================================================================
{
    long       lStartTime   = 0;
    bool       bSpoiler     = false;        // Triggers application of spoilers instead of diffusion encoding gradients
    int        iI           = 0;

    // Since by definition there is no direction linked to the trace diffusion 
    // encoding, the unity rotation matrix can be used.
    sROT_MATRIX theUnityRotMatrix;

    // Initialize error return code in case of unexpected bail-outs
    setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);

    if(! isPrepared())
    {
        SEQ_TRACE_ERROR.print("ERROR: Module is not prepared");

        return false;
    }

    if(m_pADC == NULL)
    {
        SEQ_TRACE_ERROR.print("ERROR: NULL pointer for m_pADC");
        return false;
    }

    // Check for valid loop counters
    if((m_iAdjScan != 0) || (m_lBValueCounter < 0)  || (m_lDirectionCounter < 0) || (m_lDiffLoopCounter < 0))
    {
        SEQ_TRACE_ERROR.print("ERROR: invalid loop counters (probably setLoopCounters has not been called)");
        return false;
    }

    if(! m_bRunOnlyOnce)
    {
        m_bRunOnlyOnce = true;
        m_bDumpBMatrix = m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_BMatrix", false);

        if(m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_diffusion_vector", false))
        {
            // Dump the rotation matrix
            SEQ_TRACE_INFO.print("               = ( %4.1f  %4.1f  %4.1f  ) ",
                       pSLC->getROT_MATRIX().dMat[0][0], pSLC->getROT_MATRIX().dMat[1][0], pSLC->getROT_MATRIX().dMat[2][0]);
            SEQ_TRACE_INFO.print("RotationMatrix = ( %4.1f  %4.1f  %4.1f  ) ",
                       pSLC->getROT_MATRIX().dMat[0][1], pSLC->getROT_MATRIX().dMat[1][1], pSLC->getROT_MATRIX().dMat[2][1]);
            SEQ_TRACE_INFO.print("               = ( %4.1f  %4.1f  %4.1f  ) ",
                       pSLC->getROT_MATRIX().dMat[0][2], pSLC->getROT_MATRIX().dMat[1][2], pSLC->getROT_MATRIX().dMat[2][2]);
        }
    }


    // -------------------------------------
    // Prepare diffusion gradient amplitudes 
    // -------------------------------------

    // scale as  b \propto g^2
    double dAmpl = m_dAmpl * sqrt(m_vdBValues[m_lBValueCounter] / m_dMaxPossibleBValue);

    for(iI = 0; iI < 16; ++iI)
    {
        m_DGP[iI].prepAmplitude(dAmpl * m_PhaseGradSign[iI]);
        m_DGR[iI].prepAmplitude(dAmpl * m_ReadGradSign[iI]);
        m_DGS[iI].prepAmplitude(dAmpl * m_SliceGradSign[iI]);
    }

    // --------------------------------------------------
    // Check whether spoiler gradients have to be applied
    // --------------------------------------------------

    // For small b-values, the implicit spoiling of the diffusion gradients will
    // not be sufficient to spoil any undesired coherence pathways. In that
    // case, spoiler gradients will be applied instead of diffusion gradients.
    bSpoiler = (rMrProt.diffusion().getalBValue()[m_lBValueCounter] < SPOILER_THRESHOLD);

    if(bSpoiler)
    {
        dAmpl = 0.;

        for(iI = 0; iI < 16; ++iI)
        {
            m_DGP[iI].prepAmplitude(0.);
            m_DGR[iI].prepAmplitude(0.);
            m_DGS[iI].prepAmplitude(0.);
        }
    }


    // ------------------------------------------
    // Prepare NCO objects for refocussing pulses
    // ------------------------------------------
    m_DFPset.prepSet(*pSLC, m_DRF);   /*! EGA-05 !*/
    m_DFPneg.prepNeg(*pSLC, m_DRF);   /*! EGA-05 !*/


    // ---------------------------------------------------------------------------
    // Calculate b matrix and prepare header information
    // ---------------------------------------------------------------------------

    calcBMatrix(bSpoiler);

    if(m_bDumpBMatrix)
    {
        SEQ_TRACE_INFO.print("               |   ( %5.1f )|      ( %7.1f   %7.1f  %7.1f )",
                   dAmpl, m_dBxx, m_dBxy, m_dBxz);
        SEQ_TRACE_INFO.print("g[%2ld]=%6.2f * | = ( %5.1f )|    b=( %7.1f   %7.1f  %7.1f )",
                   m_lDirectionCounter, dAmpl, dAmpl, m_dBxy, m_dByy, m_dByz);
        SEQ_TRACE_INFO.print("-              |D  ( %5.1f )|P   = ( %7.1f   %7.1f  %7.1f )",
                   dAmpl, m_dBxz, m_dByz, m_dBzz);
        SEQ_TRACE_INFO.print("tr b[%ld]=%6.1f", m_lDirectionCounter, m_dBxx + m_dByy + m_dBzz);
    }

    NLS_STATUS lStatus = CalculateDicomHeaderInformation
        (
        rMrProt.sliceSeries()[static_cast<int32_t>(pSLC->getSliceIndex())], // Input: sliceSeries()
        pSLC,            // Input: getSliceIndex, getROT_MATRIX()
        NULL,            // Input: Diffusion direction vector information (not available)
        dAmpl,           // Input: Diffusion amplitude information
        1,               // Input: Direction sign
        m_lBValueCounter,     // Input: Current bvalue index 
        m_lDirectionCounter,  // Input: Current direction vector index 
        m_pADC             // Output: Mdh.setIceProgramPara 
        );

    if(NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        return setNLSStatus(lStatus);
    }

    // --------------------------
    // Now play out the events...
    // --------------------------

    // Event table #1: TE fill time in the first half
    // ----------------------------------------------
    if(m_lPreFill)
    {
        fRTEBInit(theUnityRotMatrix);        /* GP->GX, GR->GY, GS->GZ */
        fRTEI(m_lPreFill, 0, 0, 0, 0, 0, 0, 0);

        if(SeqUT.isUnitTestActive())
        {
            mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRun, 'R', 0, pSLC->getSliceIndex(), 0, 0);
        }

        if(setNLSStatus(fRTEBFinish()))
        {
            SEQ_TRACE_ERROR.print("ERROR: Execution of event block 1 failed.");
            return false;
        }
    }

    // Event table #2: diffusion gradients
    // ----------------------------------------------
    if(bSpoiler)
    {
        // -------------
        // Apply spoiler instead of trace train
        // -------------
        long lEventBlockDuration = 10 * m_DGP[0].getTotalTime();

        fRTEBInit(pSLC->getROT_MATRIX());
        m_DSr1.setStartTime(lEventBlockDuration - m_DSr1.getTotalTime());
        runGradient(&m_DSr1);
        fRTEI(lEventBlockDuration, 0, 0, 0, 0, 0, 0, 0);

        if(SeqUT.isUnitTestActive())
        {
            mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'Q', 0, pSLC->getSliceIndex(), 0, 0);
        }

        if(setNLSStatus(fRTEBFinish()))
        {
            SEQ_TRACE_ERROR.print("ERROR: Execution of event block 2 (spoiler) failed.");
            return false;
        }
    }
    else
    {
        // -----------------------------
        // Apply trace train, first part
        // -----------------------------
        fRTEBInit(theUnityRotMatrix);        /* GP->GX, GR->GY, GS->GZ */

        if(m_bDiffusionGradientsEnabled)
        {
            for(iI = 0; iI < 10; iI++)
            {
                lStartTime = iI * m_DGP[0].getTotalTime();

                fRTEI(lStartTime, 0, 0, 0, &m_DGP[iI], 0, 0, 0);
                fRTEI(lStartTime, 0, 0, 0, 0, &m_DGR[iI], 0, 0);
                fRTEI(lStartTime, 0, 0, 0, 0, 0, &m_DGS[iI], 0);
            }
        }

        fRTEI(10 * m_DGP[0].getTotalTime(), 0, 0, 0, 0, 0, 0, 0);

        if(SeqUT.isUnitTestActive())
        {
            mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRun, 'R', 0, pSLC->getSliceIndex(), 0, 0);
        }

        if(setNLSStatus(fRTEBFinish()))
        {
            SEQ_TRACE_ERROR.print("ERROR: Execution of event block 2 failed.");
            return false;
        }
    }

    // Event table #3: inversion pulse
    // -------------------------------------
    fRTEBInit(pSLC->getROT_MATRIX());
    /*   Start Time              |    NCO   |  SRF         |  ADC  |    Gradient Events    | Sync   */
    /*     (usec)                |   Event  | Event        | Event | phase | read  | slice | Event  */
    fRTEI(m_DFPset.getStartTime(), &m_DFPset, 0, 0, 0, 0, 0, 0);
    fRTEI(m_DRF.getStartTime(), 0, &m_DRF, 0, 0, 0, 0, 0);
    fRTEI(m_DFPneg.getStartTime(), &m_DFPneg, 0, 0, 0, 0, 0, 0);

   runGradient(&m_DGSS);

   fRTEI(m_DGSS.getStartTime()+m_DGSS.getTotalTime(), 0, 0, 0, 0, 0, 0, 0);

    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'S', 0, pSLC->getSliceIndex(), 0, 0);

    if(setNLSStatus(fRTEBFinish()))
    {
        SEQ_TRACE_ERROR.print("ERROR: Execution of event block 3 failed.");
        return false;
    }


    // Event table #4: diffusion gradients / stimulation delay
    // -------------------------------------------------------
    if(bSpoiler)
    {
        // -------------
        // Apply spoiler instead of trace train
        // -------------
        fRTEBInit(pSLC->getROT_MATRIX());
        m_DSr1.setStartTime(0);
        runGradient(&m_DSr1);
        fRTEI(6 * m_DGP[0].getTotalTime() + m_lStimoDelayus, 0, 0, 0, 0, 0, 0, 0);

        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'U', 0, pSLC->getSliceIndex(), 0, 0);

        if(setNLSStatus(fRTEBFinish()))
        {
            SEQ_TRACE_ERROR.print("ERROR: Execution of event block 4 (spoiler) failed.");
            return false;
        }
    }
    else
    {
        // ------------------------------
        // Apply trace train, second part
        // ------------------------------
        fRTEBInit(theUnityRotMatrix);        /* GP->GX, GR->GY, GS->GZ */

        if(m_bDiffusionGradientsEnabled)
        {
            for(iI = 10; iI < 16; iI++)
            {
                lStartTime = (iI - 10) * m_DGP[0].getTotalTime();

                fRTEI(lStartTime, 0, 0, 0, &m_DGP[iI], 0, 0, 0);
                fRTEI(lStartTime, 0, 0, 0, 0, &m_DGR[iI], 0, 0);
                fRTEI(lStartTime, 0, 0, 0, 0, 0, &m_DGS[iI], 0);
            }
        }
        fRTEI(6 * m_DGP[0].getTotalTime() + m_lStimoDelayus, 0, 0, 0, 0, 0, 0, 0);

        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRun, 'T', 0, pSLC->getSliceIndex(), 0, 0);

        if(setNLSStatus(fRTEBFinish()))
        {
            SEQ_TRACE_ERROR.print("ERROR: Execution of event block 4 failed.");
            return false;
        }
    }

    // Event table #5: TE fill time after the second half
    // --------------------------------------------------
    if(m_lPostFill)
    {
        fRTEBInit(theUnityRotMatrix);        /* GP->GX, GR->GY, GS->GZ */
        fRTEI(m_lPostFill, 0, 0, 0, 0, 0, 0, 0);

        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRun, 'R', 0, pSLC->getSliceIndex(), 0, 0);

        if(setNLSStatus(fRTEBFinish()))
        {
            SEQ_TRACE_ERROR.print("ERROR: Execution of event block 5 failed.");
            return false;
        }
    }

    setNLSStatus(MRI_SBB_SBB_NORMAL);
    return true;
}



// ===========================================================================
bool SBBDiffusion_Trace::prepTiming(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, long lActualTE)
// ===========================================================================
{
    MrProtFacade::Pointer protFacade(new MrProtFacade(rMrProt));

    // Flag: Print a message before every return statement
    bool bDebugReturn = ((!(rSeqLim.isContextPrepForBinarySearch() || rSeqLim.isContextPrepForScanTimeCalculation())) && !m_bResolve);

    long lt1          = 0;
    long lt2          = 0;
    long lDuration    = 0;
    int  iI           = 0;

    //  ----------------------------------------------------
    /** Set up timing and calculate maximum possible b value */
    //  ----------------------------------------------------

    sRFPulseProperties myRFPulseProperties = m_RFPulseLibrary.getPulsePropertiesRefocusing(rMrProt);

    // Mind that m_lRFGradDuration must be on the double gradient raster time!
    long lRFGradDuration = long(fSDSDoubleRoundUp(
        0.0, m_MaxValueForRounding, static_cast<double>(myRFPulseProperties.lDuration_us), 2. * GRAD_RASTER_TIME));

    // These durations have to be set first: they are required
    // in order to calculate lt1 / lt2
    m_DGSS.setRampTimes(m_lRampTimeRF);

    m_DRF.setDuration(myRFPulseProperties.lDuration_us);
    m_DGSS.setDuration(lRFGradDuration + m_lRampTimeRF);   // Total time must be on double gradient raster!

    // Time available for one diffusion encoding gradient before / after refocussing
    lt1 = ((lActualTE/2 - m_lSpinPrepTimeus                   - m_DGSS.getTotalTime()/2) / 10) - m_lRampTime;
    lt2 = ((lActualTE/2 - m_lADCusTillEcho  - m_lStimoDelayus - m_DGSS.getTotalTime()/2) /  6) - m_lRampTime;

    lDuration = fSDSRoundDownGRT(std::min(lt1, lt2));

    // Do some cross checks
    // --------------------
    if(m_lRampTime > lDuration)
    {
        setNLSStatus(MRI_SBB_SBB_NEGATIV_TEFILL);

        if(bDebugReturn)
        {
            SEQ_TRACE_ERROR.print("ERROR: lDuration = %ldus is too short or negative (ramp time: %ldus)", lDuration, m_lRampTime);
        }

        return false;
    }


    // -------------------------------------------------
    // 2nd event block (diffusion gradients)
    // -------------------------------------------------
    for(iI = 0; iI < 16; ++iI)
    {
        m_DGP[iI].setAxis(SEQ::AXIS_PHASE);
        m_DGP[iI].setAmplitude(m_dAmpl * m_PhaseGradSign[iI]);
        m_DGP[iI].setStartTime(0);
        m_DGP[iI].setDuration(lDuration);  // Remember: gradient pulse duration = RampUp + FlatTop
        m_DGP[iI].setRampTimes(m_lRampTime);
    }


    // -------------------------------------------------
    // 3rd event block (refocussing pulse)
    // -------------------------------------------------

    // For details on the definition of the RF pulses, refer to the comments in a_ep2d.cpp.
    long   lDRFStartTime           = m_lRampTimeRF + (lRFGradDuration - myRFPulseProperties.lDuration_us) / 2;

    // Duration of m_DRF has already been set above
    m_DRF.setFamilyName("SE2560A180.SE180_12A2_2");
    m_DRF.setTypeRefocussing();
    m_DRF.setFlipAngle(180.0);
    m_DRF.setInitialPhase(0);
    m_DRF.setThickness(m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_DRF.setStartTime(lDRFStartTime);  // Symmetric insertion on m_DGSS flat top

    // Apply gradient reversal to refocussing pulse (just for ep2d_diff but not for RESOLVE)
    if(protFacade->isGradientReversalDiffusion() && !m_bResolve)
        m_DRF.setRequiredGSPolarity(-1.0);
    else
        m_DRF.setRequiredGSPolarity(+1.0);

    // If dRFPulseThickness is deliberately set unequal 1.0, SeqUT needs to know this.
    SeqUT.setRFThicknessInfo(&m_DRF, m_dRFPulseThicknessFactor * m_dSliceThickness);

    if((m_dPreparedSlcThk != m_dRFPulseThicknessFactor * m_dSliceThickness) || (!rSeqLim.isContextPrepForBinarySearch()))
    {
        if(! m_DRF.prepExternal(rMrProt, rSeqExpo))
        {
            if(bDebugReturn)
            {
                SEQ_TRACE_ERROR.print("ERROR: DRF.prepExternal failed.");
            }
            setNLSStatus(m_DRF.getNLSStatus());
            return false;
        }
        m_dPreparedSlcThk = m_dRFPulseThicknessFactor * m_dSliceThickness;
    }

    m_DFPset.setStartTime(m_DRF.getStartTime());
    m_DFPneg.setStartTime(m_DRF.getStartTime() + m_DRF.getDuration());

    // Duration and ramp times of m_DGSS have already been set above
    m_DGSS.setAmplitude(m_DRF.getGSAmplitude());
    m_DGSS.setAxis(SEQ::AXIS_SLICE);
    m_DGSS.setStartTime(0);

    if(! m_DGSS.prep())
    {
        if(bDebugReturn)
        {
            SEQ_TRACE_ERROR.print("ERROR: m_DGSS.prepAmplitude failed.");
        }
        setNLSStatus(m_DGSS.getNLSStatus());
        return false;
    }

    if(! m_DGSS.check())
    {
        if(bDebugReturn)
        {
            SEQ_TRACE_ERROR.print("ERROR: m_DGSS.check failed.");
        }
        setNLSStatus(m_DGSS.getNLSStatus());
        return false;
    }

    // Calculate fill times (required for b-value calculation)
    m_lPreFill  =  lActualTE/2 - m_lSpinPrepTimeus                   - m_DGSS.getTotalTime()/2 - 10 * m_DGP[0].getTotalTime();
    m_lPostFill =  lActualTE/2 - m_lADCusTillEcho  - m_lStimoDelayus - m_DGSS.getTotalTime()/2 -  6 * m_DGP[0].getTotalTime();

    // ------------------------------------------------------------
    /** Calculate the b value that will be obtained with this timing
        and store it in SBBDiffusion_Base::m_dMaxPossibleBValue.
        Remember that all gradients have been prepared with maximum amplitude.
        Only spoilers on the P-axis are considered in this b value estimation.
        */
    // ------------------------------------------------------------

    m_dMaxPossibleBValue = calcBValue();

    return true;

}



// ===========================================================================
double SBBDiffusion_Trace::calcBValue(void)
// ===========================================================================

{
    // This function is called within protocol preparation, thus no actual
    // diffusion direction is defined. At least the diffusion encoding pulses
    // in phase encoding direction (see prepTiming) are
    // prepared with maximum amplitude. For the calculation of the maximum
    // b-value, it is sufficient to look at the corresponding diagonal element
    // (Bxx) and multiply this b-value with the square (since the
    // b-value scales with the square of the gradient amplitude) of the 
    // maximum norm of the diffusion direction.
    //
    // e.g. orthogonal diffusion directions: norm = 1
    // e.g. room diagonals (3-scan-trace)  : norm = sqrt(3)

    double dBValue = 0.;

    // Set gamma of actual nucleus
    m_theBMatrix.setGamma(m_dGamma);

    // Prepare for b-value calcuation (set up diffusion encoding events):
    //  b-value only, no spoilers
    PrepBMatrix(true, false);

    // Calculate b-value
    if(!m_theBMatrix.bCalcBValue(dBValue))
    {
        SEQ_TRACE_WARN.print("WARNING: Non-null moment");
    }

    // Factor 3 due to trace
    return (dBValue * 3.);
}


// ===========================================================================
/// Calculate the b matrix which will be obtained by the current timing.
/** Each component is calculated individually by summing up the moments
    generated by the gradients called on the axis.

    \pre For the calculation, it is supposed that all diffusion gradients have
    been prepared completely.

    \return The result will be returned in the member variables
    SBBDiffusion_Base::m_dBxx, m_dByy, m_dBzz, m_dBxy, m_dBxz, and m_dByz.

    In the DEBUG version, there is a check that the gradient moments are
    zero at the end of the kernel. If this check fails, a warning message
    will be generated.
    */
void SBBDiffusion_Trace::calcBMatrix(bool bApplySpoiler)
// ===========================================================================
{
    // Set gamma of actual nucleus
    m_theBMatrix.setGamma(m_dGamma);

    // Prepare for b-matrix calcuation (set up diffusion encoding events):
    //  complete b-matrix, spoilers on demand
    PrepBMatrix(false, bApplySpoiler);

    // Calculate b-matrix
    if(!m_theBMatrix.bCalcBMatrix(m_dBxx, m_dBxy, m_dBxz, m_dByy, m_dByz, m_dBzz))
    {
        SEQ_TRACE_WARN.print("WARNING: Non-null moment");
    }

    // Calculate b-value
    m_theBMatrix.bCalcBValue(m_dBValue);
}

void SBBDiffusion_Trace::calcBMatrix(void)
// ===========================================================================
{
    // Prepare for b-matrix calcuation (set up diffusion encoding events):
    //  complete b-matrix, no spoilers
    PrepBMatrix(false, false);

    // Calculate b-matrix with default gamma
    if(!m_theBMatrix.bCalcBMatrix(m_dBxx, m_dBxy, m_dBxz, m_dByy, m_dByz, m_dBzz))
    {
        SEQ_TRACE_WARN.print("WARNING: Non-null moment");
    }


    // Calculate b-value
    m_theBMatrix.bCalcBValue(m_dBValue);
}

bool SBBDiffusion_Trace::PrepBMatrix(bool bBValueOnly, bool bApplySpoiler)
// ===========================================================================
{
    // Initialize b-matrix calculation
    m_theBMatrix.Reset();

    // Fill in basic diffusion encoding events
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DGP[0]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DGP[1]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad3_Start), m_DGP[2]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad4_Start), m_DGP[3]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad5_Start), m_DGP[4]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad6_Start), m_DGP[5]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad7_Start), m_DGP[6]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad8_Start), m_DGP[7]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad9_Start), m_DGP[8]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad10_Start), m_DGP[9]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad11_Start), m_DGP[10]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad12_Start), m_DGP[11]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad13_Start), m_DGP[12]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad14_Start), m_DGP[13]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad15_Start), m_DGP[14]);       // PE-axis diffusion
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad16_Start), m_DGP[15]);       // PE-axis diffusion

    m_theBMatrix.bAddEvent(getEventTime(RefocRF1_Center), m_DRF);       // refocussing RF

    if(!bBValueOnly)
    {
        // Full b-matrix calculation
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DGR[0]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DGR[1]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad3_Start), m_DGR[2]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad4_Start), m_DGR[3]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad5_Start), m_DGR[4]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad6_Start), m_DGR[5]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad7_Start), m_DGR[6]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad8_Start), m_DGR[7]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad9_Start), m_DGR[8]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad10_Start), m_DGR[9]);   // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad11_Start), m_DGR[10]);  // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad12_Start), m_DGR[11]);  // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad13_Start), m_DGR[12]);  // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad14_Start), m_DGR[13]);  // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad15_Start), m_DGR[14]);  // RO-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad16_Start), m_DGR[15]);  // RO-axis diffusion

        m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DGS[0]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DGS[1]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad3_Start), m_DGS[2]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad4_Start), m_DGS[3]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad5_Start), m_DGS[4]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad6_Start), m_DGS[5]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad7_Start), m_DGS[6]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad8_Start), m_DGS[7]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad9_Start), m_DGS[8]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad10_Start), m_DGS[9]);   // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad11_Start), m_DGS[10]);  // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad12_Start), m_DGS[11]);  // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad13_Start), m_DGS[12]);  // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad14_Start), m_DGS[13]);  // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad15_Start), m_DGS[14]);  // SL-axis diffusion
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad16_Start), m_DGS[15]);  // SL-axis diffusion

        m_theBMatrix.bAddEvent(getEventTime(SliceGrad1_Start), m_DGSS);     // SL-axis selection

        if(bApplySpoiler)
        {
            m_theBMatrix.bAddEvent(getEventTime(SpoilGrad1_Start), m_DSr1);   // RO-axis spoiler 1
            m_theBMatrix.bAddEvent(getEventTime(SpoilGrad2_Start), m_DSr1);   // RO-axis spoiler 2
        }
    }

    if(m_theBMatrix.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}

bool SBBDiffusion_Trace::prepGPALoadDiff(double /* dAmplitudeX */)
{
    SEQ_TRACE_ERROR.print("ERROR: Not implemented.");

    return false;
}

bool SBBDiffusion_Trace::prepGPALoadDiff(double /* dAmplitudeX */, double /* dAmplitudeY */, double /* dAmplitudeZ */)
{
    SEQ_TRACE_ERROR.print("ERROR: Not implemented.");

    return false;
}

bool SBBDiffusion_Trace::scaleGPALoadDiff(double /* dScaleX */)
{
    SEQ_TRACE_ERROR.print("ERROR: Not implemented.");

    return false;
}

bool SBBDiffusion_Trace::scaleGPALoadDiff(double /* dScaleX */, double /* dScaleY */, double /* dScaleZ */)
{
    SEQ_TRACE_ERROR.print("ERROR: Not implemented.");

    return false;
}


long SBBDiffusion_Trace::getEventTime(EnumEventTime eEvent)
{
    // Start of 1st diffusion encoding gradient
    long lTime0 = m_lPreFill;
    // Start time of the 2nd diffusion encoding gradient
    long lTime1 = lTime0 + m_DGP[0].getTotalTime();
    // Start time of the 3rd diffusion encoding gradient
    long lTime2 = lTime1 + m_DGP[0].getTotalTime();
    // Start time of the 4th diffusion encoding gradient
    long lTime3 = lTime2 + m_DGP[0].getTotalTime();
    // Start time of the 5th diffusion encoding gradient
    long lTime4 = lTime3 + m_DGP[0].getTotalTime();
    // Start time of the 6th diffusion encoding gradient
    long lTime5 = lTime4 + m_DGP[0].getTotalTime();
    // Start time of the 7th diffusion encoding gradient
    long lTime6 = lTime5 + m_DGP[0].getTotalTime();
    // Start time of the 8th diffusion encoding gradient
    long lTime7 = lTime6 + m_DGP[0].getTotalTime();
    // Start time of the 9th diffusion encoding gradient
    long lTime8 = lTime7 + m_DGP[0].getTotalTime();
    // Start time of the 10th diffusion encoding gradient
    long lTime9 = lTime8 + m_DGP[0].getTotalTime();

    // Start time of slice selection gradient
    long lTime10 = lTime9 + m_DGP[0].getTotalTime() + m_DGSS.getStartTime();
    // Refocussing time
    long lTime11 = lTime10 + m_DGSS.getTotalTime()/2;

    // Start time of the 11th diffusion encoding gradient
    long lTime12 = lTime10 + m_DGSS.getTotalTime();
    // Start time of the 12th diffusion encoding gradient
    long lTime13 = lTime12 + m_DGP[0].getTotalTime();
    // Start time of the 13th diffusion encoding gradient
    long lTime14 = lTime13 + m_DGP[0].getTotalTime();
    // Start time of the 14th diffusion encoding gradient
    long lTime15 = lTime14 + m_DGP[0].getTotalTime();
    // Start time of the 15th diffusion encoding gradient
    long lTime16 = lTime15 + m_DGP[0].getTotalTime();
    // Start time of the 16th diffusion encoding gradient
    long lTime17 = lTime16 + m_DGP[0].getTotalTime();

    // Start time of first spoiler (applied if b-value is below threshold)
    long lTime18 = 10 * m_DGP[0].getTotalTime() - m_DSr1.getTotalTime();
    // Start time of second spoiler (applied if b-value is below threshold)
    long lTime19 = lTime10 + m_DGSS.getTotalTime();


    switch(eEvent)
    {
        case DiffGrad1_Start:
            return lTime0;
        case DiffGrad2_Start:
            return lTime1;
        case DiffGrad3_Start:
            return lTime2;
        case DiffGrad4_Start:
            return lTime3;
        case DiffGrad5_Start:
            return lTime4;
        case DiffGrad6_Start:
            return lTime5;
        case DiffGrad7_Start:
            return lTime6;
        case DiffGrad8_Start:
            return lTime7;
        case DiffGrad9_Start:
            return lTime8;
        case DiffGrad10_Start:
            return lTime9;

        case DiffGrad11_Start:
            return lTime12;
        case DiffGrad12_Start:
            return lTime13;
        case DiffGrad13_Start:
            return lTime14;
        case DiffGrad14_Start:
            return lTime15;
        case DiffGrad15_Start:
            return lTime16;
        case DiffGrad16_Start:
            return lTime17;

        case SpoilGrad1_Start:
            return lTime18;
        case SpoilGrad2_Start:
            return lTime19;

        case SliceGrad1_Start:
            return lTime10;
        case RefocRF1_Center:
            return lTime11;
        default:
            SEQ_TRACE_ERROR.print("unknown event enumerator: %i", eEvent);
            return -1;
    }
}

long SEQ_NAMESPACE::SBBDiffusion_Trace::getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch)
{
    return m_lBValueInc_Limit;
}

