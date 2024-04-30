//  [ Compilation unit *******************************************************
//
// Project: NUMARIS/4
//
// File: \src\MrImaging\seq\a_ep2d.cpp
//
//  Author: Clinical
//          Thomas Kluge;    Siemens AG Med MRIA/Seq;  (09131) 84-8049 (01/1999-02/2002)
//          Michael Zwanger; Siemens AG Med MREA-fMRI; (09131) 84-2672 (only ep2d_diff)
//          Josef Pfeuffer;  Siemens AG Med MR PLM AW Neuro; (only ep2d_asl)
//          PLM AW Neuro
//
// Lang        : cpp
//
// Descrip: Source Code for Standard NUMARIS/4 Sequences of type EPI.
//

/** \requirement N4_elh_DIFF_FLAIR,
N4_elh_DTI_measurement_Parameter_Range_Matrix_PAT


\file   a_ep2d.cpp
\brief  File containing source code for the sequences
- ep2d_diff
- ep2d_perf
- ep2d_se
- ep2d_fid
- ep2d_bold
- ep2d_pace
- ep2d_asl
- tgse_asl

This file contains the implementation of the class Ep2d.

***************************************************************************

\changed    19-Mar-2007; T.Feiweier; 4c11a; CHARM: n.a.
\description
Adaptive IR slice thickness for SE and DIFF variants

If an inversion pulse is used with diffusion epi, the doubled inversion thickness
can generate crosstalk issues (e.g. STIR: imperfect fat suppression) if SeqLoop
dictates a nesting of inversions (e.g. invert slice 1 - invert slice 3 - scan
slice 1 - scan slice 3). Thus, an inversion thickness identical to the slice
thickness seems reasonable.
However, an inversion pulse is also used for fluid attenuation (FLAIR). Due to
inflow effects, an increased inversion thickness significantly enhances the quality
of fluid attenuation.
As a compromise, a doubled inversion thickness will be only used if the inversion
time exceeds 500ms, which is a reasonable indication of a FLAIR protocol.


\changed    19-Mar-2007; T.Feiweier; 4c11a; CHARM: n.a.
\description
Monopolar water excitation pulse

Although the bipolar WE-pulse allows for better slice profiles, the monopolar
variant provides more robust water excitation quality. At the same time, the
gradient polarity for the SE-variant has to be inverted in order to avoid
33rd arm artefacts (gradient amplitudes of excitation and refocussing must not
be the identical).

***************************************************************************

*/

// ------------------------------------------------------------------------------
// General includes
// ------------------------------------------------------------------------------

// MrProt
#include "MrProtSrv/Domain/CoreNative/MeasAqcDefs.h" // MOSAIC image size limit
#include "MrProtSrv/Domain/CoreNative/MrFastImaging.h"
#include "MrProtSrv/Domain/CoreNative/MrNavigator.h"
#include "MrProtSrv/Domain/CoreNative/MrPat.h"
#include "MrProtSrv/Domain/CoreNative/MrPreparationPulses.h"
#include "MrProtSrv/Domain/CoreNative/Slice.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Application/Application.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Filter/MrFilter.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrRXSpec.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrSysSpec.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrTXSpec.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrCoilInfo.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrSliceGroup.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Physiology/MrPhysiology.h"
#include "MrProtSrv/Domain/CoreNative/MrChannelMatrix.h"
// MrProt

// MrProt KSpace Wrapper
#include "MrImaging/seq/SystemProperties.h" // * Siemens system properties *
#include "MrMeasSrv/SeqIF/csequence.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/CoilSelect/MrRxCoilSelect.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpoRFBlockInfo.h"
#ifdef WIN32
#include "MrMeasSrv/SeqIF/Sequence/ISequence.h"
#endif

#include "MrImaging/libSBB/SBBOptfs.h"
#include "MrImagingFW/libSBBFW/libSBBFW.h" // for fSBBMeasRepetDelayGetDurationForZeroPauseUs()
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrMeasSrv/MeasUtils/MrTimeStamp.h" // getTimeus
// EPIPhaseCorrPEFunctor feedback format
#include "MrVista/Ice/IceScanFunctors/EPIDorkFeedbackData.h"
// EPIPhaseCorrPEFunctor parametrization defins
#include "MrVista/Ice/IceScanFunctors/EPIPhaseCorrPEFunctorDefs.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/JustVol/MrProtJustVol.h"           // calcAdjSlice()
#include "MrAdjustSrv/AdjAccessIF/AdjAccessIF.h"                                // isLocalShimRequested()

#if (defined SUPPORT_iPAT_a_ep2d) || (defined SUPPORT_iPAT_TGSE)
//  #pragma message ("NOTE: This single shot EPI sequence supports iPAT.")
#include "MrImaging/seq/common/iPAT/iPAT.h"
#endif

#include "MrImagingFW/libSeqUTIF/libsequt.h" // mSEQTest, SeqUT

#ifdef ASL
#include "MrProtSrv/Domain/CoreNative/MrAsl.h"
#endif

#ifdef SUPPORT_PACE
#include "MrMeasSrv/SeqIF/libMES/SEQData.h"
//  Needed for PACE::fInit
#include "MrImaging/SequenceLibraries/libPace/PACE.h"
#endif

// SMS support
#include "MrImaging/libSeqUtil/SMSProperties.h"
#include "MrImaging/libSeqUtil/SliceAccelerationUtils.h"
// SMS support

// trace support
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

//#ifdef PACE3D
//#include "MrVista\include\Ice\IceImagePostProcFunctors\MotionCorrFeedbackData.h"
//#endif

// sequence messages
#include "MrMeasSrv/SeqIF/Sequence/sequmsg.h"

#if defined COMPILE_EP2D_SE
#include "MrImaging/seq/greMRE/ParameterCheck.h"
#endif

// ------------------------------------------------------------------------------
// Application includes
// ------------------------------------------------------------------------------
#include "MrImaging/seq/a_ep2d.h"
#include "MrImaging/seq/common/IterativeDenoisingUIParameter/IterativeDenoisingUIParameter.h"

#ifdef EP2D_MS
#include "MrImaging/seq/a_ep2d_se_ms/KSpaceHelpers.h"
#include "MrImaging/seq/a_ep2d_se_ms/msEPIDefs.h"
#endif

//-------------------------------------------------------------------------------------
// macro for error-handling
//-------------------------------------------------------------------------------------
#define mSBBErrGotoFinish(A, B)                             \
    {                                                       \
        lStatus = (A).getNLSStatus();                       \
        if (!rSeqLim.isContextPrepForBinarySearch())        \
        {                                                   \
            SEQ_TRACE_ERROR.print("%s: 0x%lx", B, lStatus); \
        }                                                   \
        goto FINISHED;                                      \
    }


#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif


//  --------------------------------------------------------------------------
//
//  Name        :  SEQIF_DEFINE
//
//  Description :
///  \brief        Create instance of the sequence
//
//  Return      :  SeqIF *
//
//  --------------------------------------------------------------------------
#ifdef SEQUENCE_CLASS_EP2D
#ifndef COMPILE_EP2D_DIFF
SEQIF_DEFINE(SEQ_NAMESPACE::Ep2d)
#endif
#endif

using namespace SEQ_NAMESPACE;

#ifdef WIN32

//  ----------------------------------------------------------------------
//
//  Name        :  getUI
//
//  Description :
/// \brief         Returns the pointer to the UI class
///
//
//  Return      :  EpUI*
//
//  ----------------------------------------------------------------------
EpUI* SEQ_NAMESPACE::getUI(MrUILinkBase* const pThis)
{
#if defined COMPILE_EP2D_SE || defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_FID
    return (static_cast<Ep2d*>(pThis->sequence().getSeq())->getUI());
#else
    return nullptr;
#endif
}

//  ----------------------------------------------------------------------
//
//  Name        :  getUI
//
//  Description :
/// \brief         Returns the pointer to the UI class
///
//
//  Return      :  EpUI*
//
//  ----------------------------------------------------------------------
EpUI* SEQ_NAMESPACE::getUI(MrMeasSrv::ISequence* const pSeq)
{
#if defined COMPILE_EP2D_SE || defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_FID
    return (static_cast<Ep2d*>(pSeq->getSeq())->getUI());
#else
    return nullptr;
#endif
}

//  ----------------------------------------------------------------------
//
//  Name        :  getSeq
//
//  Description :
/// \brief         Returns the pointer to the sequence Ep2d
//
//  Return      :  Ep2d*
//
//  ----------------------------------------------------------------------
Ep2d* SEQ_NAMESPACE::getSeq(MrUILinkBase* const pThis)
{
    return (static_cast<Ep2d*>(pThis->sequence().getSeq()));
}

#endif // #ifdef WIN32

//  ----------------------------------------------------------------------
//
//  Name        :  useB1ControlLoop
//
//  Description :
/// \brief         Decides (based on protocol parameters) whether
//                 the B1 control loop is used or not. The method does not change the
//                 values in MrProt. This task has to be done in e.g. the UILink handlers.
//                 Values for the two referenced variables fCorrectionFactorMax and
//                 fPeakReserveFactor are determined from the protocol. If the B1 control loop
//                 is not used the defauls values for the variables (1.0 and 0.0) are provided.
//
//  Return      :  true  - if B1 control loop is used
//                 false - else
//
//  ----------------------------------------------------------------------
bool SEQ_NAMESPACE::useB1ControlLoop(MrProt& rMrProt, float& fCorrectionFactorMax, float& fPeakReserveFactor)
{
    MrProtFacade protFacade(rMrProt);

    // set default values
    bool bUseB1ControlLoop = false;
    fCorrectionFactorMax   = 1.0f;
    fPeakReserveFactor     = 0.0f;

    // the B1 control loop shall only be activated in the ep2d_diff
    // conditions for B1 control loop active:
    // 1) the protocol is a DTI protocol ==> diffusion modes free, MDDW or q-Space

    // if conditions are fulfilled ==> 10% for B1 control loop
    if (protFacade.iSDTI())
    {
        bUseB1ControlLoop    = true;
        fCorrectionFactorMax = 1.10f;
        fPeakReserveFactor   = 0.10f;
    }

    // activate B1 control loop also for ep2d_bold scans (charm #445264)
#if (defined BOLD && !defined PACE3D)
    bUseB1ControlLoop    = true;
    fCorrectionFactorMax = 1.10f;
    fPeakReserveFactor   = 0.10f;
#endif

    // return value
    return bUseB1ControlLoop;
}

/*###########################################################################
SeqLoopEP2D members
--> Are now located here:
"MrServers/MrImaging/seq/common/SeqLoopEPI/SeqLoopEPI.h"
###########################################################################*/

Ep2d::Ep2d()
{
    // ---------------------------------------------------------------------------
    // initialize osc-bit control-flags
    // ---------------------------------------------------------------------------
    initializeOscBitFlags();

#ifdef EP2D_MS
    m_isForcedFastGRERefscanWithoutPAT
        = SysProperties::ReadSeqSettingGeneral("FORCED_REFSCAN_WITHOUT_PAT", false, true);
#endif
}

Ep2d::~Ep2d()
{
    //  ----------------------------------------------------------------------
    ///  Delete existing UI instance
    //  ----------------------------------------------------------------------
    if (m_pUI)
    {
        delete m_pUI;
        m_pUI = nullptr;
    }
}

/*[ Function ****************************************************************\
*
* Name        : calculateTRTIFillTimes
*
* Description : Calculates TR- and TI-fill times for SeqLoop and required TR
*               and TI values.
*
* Return      : true, if success
*
\****************************************************************************/
bool Ep2d::calculateTRTIFillTimes(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long* plNeededTI, long* plNeededTR)
{
#undef DEBUG_ORIGIN
#define DEBUG_ORIGIN 0x00000200

    // ---------------------------------------------------------------------------
    // calculate the basic scan time and time needed for sats.  Needed for spir calculation
    // it is necessary to add the time spent for the spir pulse but first this must be calculated
    // ---------------------------------------------------------------------------
    long lScanTimeSatsEtc = m_mySeqLoop.getlScanTimeAllSats();   // SBB scan time
    long lScanTimeBasic   = m_EPIKernel.getDurationPerRequest(); // Kernel scan time excluding mandatory fill time
    long lScanTimeStore   = lScanTimeSatsEtc;                    // Helper

#ifdef COMPILE_EP2D_DIFF
    if (m_EPIKernel.getbCompensationEnable())
    {
        if (m_EPIKernel.getPointerCompGrad()->isPrepared())
        {
            lScanTimeBasic += m_EPIKernel.getPointerCompGrad()->getDurationPerRequest();
        }
    }
#endif

    // The sequence asks for a pointer to the SPAIR SBB from SeqLoop.  Using this pointer it then uses the SBB to calculate the required inversion time.
    // This is set in the SBB from within the SBB.
    if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_Spair)
    {
        SeqBuildBlockOptfs*     pSBBOptfs     = m_mySeqLoop.getpOptfs();
        SeqBuildBlockOptfsPrep* pSBBOptfsPrep = m_mySeqLoop.getpOptfsPrep();

        // Add SPAIR time (depending on protocol parameters, e.g. TR) to lScanTimeSatsEtc
        if (pSBBOptfs == nullptr || !pSBBOptfs->calcSPIRTime(rMrProt, rSeqLim, rSeqExpo, lScanTimeBasic, lScanTimeSatsEtc, 0, SeqBuildBlockOptfs::SPIR_CALC_TYPE_EPI, pSBBOptfsPrep))
        {
            SEQ_TRACE_ERROR.print("The calc SPIR time has failed, prob due to an error in the tickle pulse");
            return false;
        }
    }

    // ---------------------------------------------------------------------------
    // consider mandatory pause
    // ---------------------------------------------------------------------------
    long lPause = m_lCoolPauseTotal;

    // Reset implicit cooling contribution
    m_lCoolPauseImplicit = 0;

    if (m_EPIKernel.getUseGPABalance())
    {
        // since gradients during the IR block (spoiling, slice selection) are small, the duration is considered for the cooling time
        considerIRBlockForImplicitCoolingPause(rSeqLim, rMrProt, rSeqExpo, lScanTimeSatsEtc, lScanTimeBasic);

        if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_FatSaturation)
        {
            // Consider duration of fat saturation module as pause
            m_lCoolPauseImplicit += fSDSRoundUpGRT(m_mySeqLoop.getlFSDuration());
        }
        else if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_Spair)
        {
            // Consider duration of SPAIR module as pause
            m_lCoolPauseImplicit += fSDSRoundUpGRT(m_mySeqLoop.getlOptFSDuration());
        }

        if (rMrProt.preparationPulses().getucMTC())
        {
            // Consider duration of MTC module as pause
            m_lCoolPauseImplicit += fSDSRoundUpGRT(m_mySeqLoop.getlMTCDuration());
        }

        for (long lI = 0; lI < rMrProt.satList().size(); ++lI)
        {
            // Consider duration of RSat module as pause
            m_lCoolPauseImplicit += fSDSRoundUpGRT(m_mySeqLoop.getlRSatDuration(lI));
        }

        // 20140320 DP: MR_00443009: due to gradient duty cycle limitations, we will not handle the excitation using ZOOMit as cooling pause any more.

        lPause -= m_lCoolPauseImplicit;

        // Restrict pause duration to meaningful values
        if (lPause < 0)
        {
            lPause               = 0;
            m_lCoolPauseImplicit = m_lCoolPauseTotal;
        }
    }
    
    bool bSuccess = false;
    long lTIMinAdd1        = 0;
    long lTIMinAdd2        = 0;
    long lNegativeFillTime = 0;

    if ( !isIIRSchemeStandard( rMrProt, rSeqLim ) )
    {
        // ---------------------------------------------------------------------------
        // calculate the TR/TI fill times
        // ---------------------------------------------------------------------------
        if ( rMrProt.preparationPulses().getucInversion() != SEQ::INVERSION_OFF )
        {
            lTIMinAdd1 = lScanTimeSatsEtc;
            lTIMinAdd2 = lScanTimeSatsEtc;
        }
        else
        {
            lTIMinAdd1 = lTIMinAdd2 = 0;
        }

        // Add fill time to kernel
        lScanTimeBasic += lPause;

        bSuccess = m_mySeqLoop.TrTiFillTimes
        (
            rMrProt,
            rSeqLim,
            rSeqExpo,
            lScanTimeBasic + lScanTimeSatsEtc, // Minimum TR (assuming no IR)
            1,                                 // Number of TR fills
            lTIMinAdd1,                        // Additional minimum TI, when not interleaved
            lTIMinAdd2,                        // Additional minimum TI, when     interleaved
            lScanTimeBasic,                    // Scan time
            lScanTimeSatsEtc,                  // SBB time
            &lNegativeFillTime
        );
        

        // The following section is required by TSE solve handlers - not needed for EPI?
        /*
        m_lSliceAcqDuration = lScanTimeBasic;
        if ( ( rMrProt.getsPrepPulses().getucInversion() != SEQ::INVERSION_OFF ) &&
            ( *plNeededTI > lScanTimeSatsEtc ) )
        {
            m_lSliceAcqDuration += *plNeededTI;
        }
        else
        {
            m_lSliceAcqDuration += lScanTimeSatsEtc;
        }
        */
    }
    else
    {

        // Update SBB scan time including SPAIR
        m_mySeqLoop.setlSBBScanTime(lScanTimeSatsEtc);

    #ifdef ASL
        // add fat sat time for strong mode
        if (rMrProt.preparationPulses().getucFatSatMode() == SEQ::FAT_SAT_STRONG)
        {
            m_mySeqLoop.setlKernelScanTime(lScanTimeBasic + lPause + m_CSatFat.getDurationPerRequest() + m_SpoilGrad.getDurationPerRequest());
        }
        else
        {
            m_mySeqLoop.setlKernelScanTime(lScanTimeBasic + lPause);
        }
    #else

        
        long lScanTimeOptPTX = 0;
    #if defined ZOOM_2DRF
        if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
        {
            lScanTimeOptPTX = m_OptPTXVolume.getDurationPerRequest();
        }
    #endif

        m_mySeqLoop.setlKernelScanTime(lScanTimeBasic + lPause + lScanTimeOptPTX);
    #endif

    #ifdef SUPPORT_FAST_IR
        m_mySeqLoop.setCoolPauseWithinKernelTime_us(static_cast<int>(lPause));
    #endif

        bSuccess = m_mySeqLoop.calcFillTimes(rMrProt, rSeqLim, rSeqExpo);
    }

    if (!bSuccess)
    {
        // What follows here is basically a copy of fUICEvaluateSeqLoopTrTiFillTimes. However, a special
        // handling is required in case that SPAIR is enabled: if the sequence demands to increase TR by
        // a certain amount, the SPAIR time needs to be updated which in turn might require an additional
        // TR increase. This hen-egg-problem needs to be explicitly resolved.

        if (m_mySeqLoop.getlTRneeded() != 0 || m_mySeqLoop.getlTIneeded() != 0)
        {
            // We could proceed, if we change TR and/or TI !

            if (m_mySeqLoop.getlTRneeded() != 0)
            {
                // NOTE: In contrast to TI mySeqLoop does NOT take care of the TR
                //       increment. So we have to do it here.

                double dInc = static_cast<double>(rSeqLim.getTR()[0].getInc());

                // Avoid division by zero
                if (fabs(dInc) < 1.e-10)
                {
                    // Unexpected error
                    // => trace also if rSeqLim.isContextPrepForBinarySearch()
                    SEQ_TRACE_ERROR.print("ERROR in calculateTRTIFillTimes(): dInc == 0.");
                    return false;
                }

                // ("..........NOTE!!!!! we assume something about TR-limit-handling in MrUILink !")
                if (m_mySeqLoop.getlTRneeded() > 10000)
                {
                    dInc *= 10.0;
                }
                if (m_mySeqLoop.getlTRneeded() > 1000000)
                {
                    dInc *= 10.0;
                }

                *plNeededTR = static_cast<long>(0.5 + dInc * ceil(static_cast<double>(m_mySeqLoop.getlTRneeded()) / dInc));

                // SPAIR:
                // in this case we have to consider NeededTR as the actual TR (because it will be applied by the solve handler)
                // and thus perform the calculation again!
                if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_Spair)
                {
                    // Make a copy of MrProt and set the needed TR
                    MrProt sTmp(*rMrProt.clone());
                    sTmp.tr()[0] = static_cast<int32_t>(*plNeededTR);

                    // Revert original SBB scan time excluding SPAIR
                    lScanTimeSatsEtc = lScanTimeStore;

                    SeqBuildBlockOptfs*     pSBBOptfs     = m_mySeqLoop.getpOptfs();
                    SeqBuildBlockOptfsPrep* pSBBOptfsPrep = m_mySeqLoop.getpOptfsPrep();
                    // Add SPAIR time (depending on protocol parameters, e.g. TR) to lScanTimeSatsEtc
                    if (pSBBOptfs == nullptr || !pSBBOptfs->calcSPIRTime(sTmp, rSeqLim, rSeqExpo, lScanTimeBasic, lScanTimeSatsEtc, 0, SeqBuildBlockOptfs::SPIR_CALC_TYPE_EPI, pSBBOptfsPrep))
                    {
                        SEQ_TRACE_ERROR.print("The calc SPIR time has failed, prob due to an error in the tickle pulse");
                        return false;
                    }

                    lPause               = m_lCoolPauseTotal;
                    m_lCoolPauseImplicit = 0;
                    if (m_EPIKernel.getUseGPABalance())
                    {
                        // Caution: We might be in trouble here ... ?
                        //
                        // We have just recalculated the SPAIR time based on the needed TR and are now
                        // going to correspondingly recalculate the needed TR. However, the required
                        // cool pause is also affected by this recalculation process...
                        //
                        // A full solution might require an iterative procedure here that seeks for
                        // a complete set of valid and consistent parameters:
                        // a) needed TR
                        // b) SPAIR time
                        // c) pause
                        //
                        // For the moment, try this direct approach and check whether any related
                        // internal errors show up.

                        // ---------------------------------------------------------------------------
                        // consider mandatory pause
                        // ---------------------------------------------------------------------------
                        lPause = m_lCoolPauseTotal;

                        // Consider duration of SPAIR module as pause
                        // Note: Since neither standard FatSat nor IR can be combined with SPAIR, they
                        //       don't need to be considered here.
                        m_lCoolPauseImplicit = fSDSRoundUpGRT(m_mySeqLoop.getlOptFSDuration());
                                                
                        if (rMrProt.preparationPulses().getucMTC())
                        {
                            // Consider duration of MTC module as pause
                            m_lCoolPauseImplicit += fSDSRoundUpGRT(m_mySeqLoop.getlMTCDuration());
                        }

                        for (long lI = 0; lI < rMrProt.satList().size(); ++lI)
                        {
                            // Consider duration of RSat module as pause
                            m_lCoolPauseImplicit += fSDSRoundUpGRT(m_mySeqLoop.getlRSatDuration(lI));
                        }

                        lPause -= m_lCoolPauseImplicit;

                        // Restrict pause duration to meaningful values
                        if (lPause < 0)
                        {
                            lPause               = 0;
                            m_lCoolPauseImplicit = m_lCoolPauseTotal;
                        }
                    }

                    if (!isIIRSchemeStandard(rMrProt, rSeqLim))
                    {
                        // Add fill time to kernel
                        lScanTimeBasic += lPause;

                        bSuccess = m_mySeqLoop.TrTiFillTimes(
                            rMrProt,
                            rSeqLim,
                            rSeqExpo,
                            lScanTimeBasic + lScanTimeSatsEtc, // Minimum TR (assuming no IR)
                            1,                                 // Number of TR fills
                            lTIMinAdd1,                        // Additional minimum TI, when not interleaved
                            lTIMinAdd2,                        // Additional minimum TI, when     interleaved
                            lScanTimeBasic,                    // Scan time
                            lScanTimeSatsEtc,                  // SBB time
                            &lNegativeFillTime);
                    }
                    else
                    {
#if defined ZOOM_2DRF
                        m_mySeqLoop.setlKernelScanTime(lScanTimeBasic + lPause + m_OptPTXVolume.getDurationPerRequest());
#else
                        // Update kernel scan time including mandatory pause
                        m_mySeqLoop.setlKernelScanTime(lScanTimeBasic + lPause);
#endif
                    
                        // Update SBB scan time including SPAIR
                        m_mySeqLoop.setlSBBScanTime(lScanTimeSatsEtc);

                        bSuccess = m_mySeqLoop.calcFillTimes(sTmp, rSeqLim, rSeqExpo);

                    }

                    dInc = static_cast<double>(rSeqLim.getTR()[0].getInc());

                    // Avoid division by zero
                    if (fabs(dInc) < 1.e-10)
                    {
                        // Unexpected error
                        // => trace also if rSeqLim.isContextPrepForBinarySearch()
                        SEQ_TRACE_ERROR.print("ERROR in calculateTRTIFillTimes(): dInc == 0.");
                        return false;
                    }

                    // ("..........NOTE!!!!! we assume something about TR-limit-handling in MrUILink !")
                    if (m_mySeqLoop.getlTRneeded() > 10000)
                    {
                        dInc *= 10.0;
                    }
                    if (m_mySeqLoop.getlTRneeded() > 1000000)
                    {
                        dInc *= 10.0;
                    }

                    long lNeededTR_new = static_cast<long>(0.5 + dInc * ceil(static_cast<double>(m_mySeqLoop.getlTRneeded()) / dInc));

                    *plNeededTR = std::max(lNeededTR_new, *plNeededTR);
                }
            }
            else
            {
                *plNeededTR = rMrProt.tr()[0];
            }

            if (m_mySeqLoop.getlTIneeded() != 0)
            {
                *plNeededTI = m_mySeqLoop.getlTIneeded();
            }
            else
            {
                *plNeededTI = rMrProt.ti()[0];
            }

            bSuccess = true;
        }
        else
        {
            // no chance
            //
            bSuccess = false;
        }
    }
    else
    {
        *plNeededTR = rMrProt.tr()[0];
        *plNeededTI = rMrProt.ti()[0];
    }

    return bSuccess;
}


//   --------------------------------------------------------------------------
//
//   Name        :  Ep2d::initialize
//
//   Description :
///  \brief        Initialization of the sequence
///
///                On the host, the object m_pUI will actually contain sensible
///                  data after Ep2d::initialize. On the measurement system, it
///                  is basically an empty object behind it.
///
//   Return      :  NLS status
//
//   --------------------------------------------------------------------------
NLSStatus Ep2d::initialize(SeqLim& rSeqLim)
{
    if (!isSystemCompatibleWithFlavor())
    {
        return MRI_SEQ_SEQU_SEQU_SYSTEM_INCOMPATIBLE;
    }

    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    //   ----------------------------------------------------------------------
    //
    //   Definition of sequence hard limits
    //
    //   ----------------------------------------------------------------------
#undef DEBUG_ORIGIN
#define DEBUG_ORIGIN 0x00000100

    //  ----------------------------------------------------------------------
    //  Instantiate of UI class
    //  ----------------------------------------------------------------------
    if ((NLS_SEV & (lStatus = createUI(rSeqLim))) == NLS_SEV)
    {
        SEQ_TRACE_ERROR.print("Instantiation of UI class failed: Gre::createUI(SeqLim&)");
        return (lStatus);
    }

    //  ----------------------------------------------------------------------
    //  set solve handlers according to thermal balancing feature toggle (only in the derived Ep2d_diff)
    //  ----------------------------------------------------------------------
#if defined WIN32
    setUIThermalBalancing();
#endif

    //-------------------------------------------------------------------------------------
    // Init the pace feedback class
    //-------------------------------------------------------------------------------------
#ifdef PACE3D
    // JR inserted PACE3d parameter (starting feedback in ICE)
    // rSeqLim.setBold3dPace(1, 0);
    rSeqLim.setBold3dPace(SEQ::ON, SEQ::OFF);

    // init 3D PACE class
    m_PaceFeedback.Init();

#ifndef ASL
    // handle time required for wakeup-eventblock within TR
    m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(TWAKEUP);
#else
    m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(TWAKEUP + m_ASL_SBB.getDurationPerRequest());
#endif
#else
    rSeqLim.setBold3dPace(SEQ::OFF);
#ifndef ASL
    // default no plugin within TR
    m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(0);
#else
    m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(m_ASL_SBB.getDurationPerRequest());
#endif
#endif

    //-------------------------------------------------------------------------------------
    // sequence should not be executed on Magnetom OPEN
    //-------------------------------------------------------------------------------------
    rSeqLim.setNotSupportedSystemTypes("007");

    //-------------------------------------------------------------------------------------
    // sequence hint text
    //-------------------------------------------------------------------------------------
    {
        const std::string ptVariant = getSequenceVariantText();
        char t[512];
        sprintf(t, "%s (compile date: %s)", ptVariant.c_str(), (char*)__DATE__);
        rSeqLim.setSequenceHintText(t);
    }

    //-------------------------------------------------------------------------------------
    // set sequence type
    //-------------------------------------------------------------------------------------
    rSeqLim.setSequenceType(SEQ::SEQUENCE_TYPE_EPI);
    rSeqLim.getSequenceType().setDisplayMode(SEQ::DM_OFF);

    //-------------------------------------------------------------------------------------
    // general settings
    //-------------------------------------------------------------------------------------
    rSeqLim.setMyOrigFilename((char*)__FILE__);

    //-------------------------------------------------------------------------------------
    // Set gradient mode options with default "fast"
    // The "normal" mode option has been introduced to allow a reduced slew-rate
    // to be set for under-voltage situations (CHARM 333764) - see fSEQPrep.
    //-------------------------------------------------------------------------------------
    rSeqLim.setGradients(SEQ::GRAD_FAST, SEQ::GRAD_FAST_GSWD_RISETIME, SEQ::GRAD_NORMAL, SEQ::GRAD_NORMAL_GSWD_RISETIME, SEQ::GRAD_ULTRAFAST, SEQ::GRAD_ULTRAFAST_GSWD_RISETIME);

    rSeqLim.setAdjustmentMode(AdjustmentMode_Standard);
    rSeqLim.setAdjSliceBySliceFirstOrderShim(SEQ::OFF);
    rSeqLim.setAdjSliceBySliceFrequency(SEQ::OFF);
    rSeqLim.setAdjSliceBySliceTxRef(SEQ::OFF);
    rSeqLim.setAdjSliceBySlicePtx(SEQ::OFF);


    // Allow distortion correction 2D, 3D, and OFF (OFF is needed for MR Neurology in combination with DWI)
    rSeqLim.setDistortionCorrMode(SEQ::DISTCORR_DIS2D, SEQ::DISTCORR_DIS3D, SEQ::DISTCORR_NDIS);

    //-------------------------------------------------------------------------------------
    // EPI kernel initialization: configure protocol-independent settings
    //-------------------------------------------------------------------------------------
    if(!initializeEPIKernel())
        return m_EPIKernel.getNLSStatus();


    //-------------------------------------------------------------------------------------
    // protocol independent configuration of SeqLoop
    //-------------------------------------------------------------------------------------
    m_mySeqLoop.setPerformTRFill(false);
#ifndef ASL
    m_mySeqLoop.setPerformSATs(true);
#else
    // we need to do the SATs ourselves
    m_mySeqLoop.setPerformSATs(false);
#endif
    m_mySeqLoop.setEffectiveTRForR_CSat(false, false);
    m_mySeqLoop.setbSpoilGradAfterKernel(false);
    m_mySeqLoop.setdDistFacMinIfConcYes(-1.0);
    m_mySeqLoop.setbHandleTRTIConflict(true);
    m_mySeqLoop.initRegistryEntries();
    m_mySeqLoop.setTrigHaltSingleShot(false);
    m_mySeqLoop.setPerformOscBit(false);

#if (defined COMPILE_EP2D_SE) || (defined COMPILE_EP2D_DIFF)
    m_mySeqLoop.setScaleIRThickness(2.0);
#endif

    //-------------------------------------------------------------------------------------
    // standard epi hard limits
    //-------------------------------------------------------------------------------------
    if (!m_pUI->fEPIStdInit(rSeqLim, (SeqBuildBlockEPIReadOut*)&m_EPIKernel, (ReorderInfo*)&m_REOInfo))
    {
        SEQ_TRACE_ERROR.print("fEPIStdInit failed");
        return MRI_SEQ_SEQU_ERROR;
    }
    rSeqLim.setAdjFreProtRelated(SEQ::ON);

    //-------------------------------------------------------------------------------------
    // standard single-shot epi hard limits
    //-------------------------------------------------------------------------------------
    rSeqLim.setMultiSliceMode(SEQ::MSM_INTERLEAVED);
    rSeqLim.setDelayTimeInTR(0, 30000000, 1000, 0);
    rSeqLim.setSliceSeriesMode(SEQ::INTERLEAVED, SEQ::ASCENDING, SEQ::DESCENDING);
    rSeqLim.setMultipleSeriesMode(SEQ::MULTIPLE_SERIES_OFF, SEQ::MULTIPLE_SERIES_EACH_MEASUREMENT);
    if (SysProperties::isUHFSystem())
        rSeqLim.setSliceThickness(0.100, 10.000, 0.100, 5.000);
    else
        rSeqLim.setSliceThickness(0.800, 10.000, 0.100, 5.000);


    // base resolution with flexible adjustment of increment
    rSeqLim.setBaseResolution(64, 512, SEQ::INC_NORMAL, 128);

#ifdef SUPPORT_iPAT_a_ep2d
    {
        // set default SeqLim for PAT: PATMode, RefScanMode, AccelFactorPE, RefLinesPE
        fPATSetDefaultSeqLim(rSeqLim);
        // EPI uses the ExtraRefScanMode. With VD15A, it is also possible to use a Flash reference scan (SEQ::PAT_REF_SCAN_EXTRA)
        rSeqLim.setRefScanMode(SEQ::PAT_REF_SCAN_EXTRA_EPI, SEQ::PAT_REF_SCAN_EXTRA);
    }
#endif

#ifdef SUPPORT_iPAT_TGSE
    {
        // set default SeqLim for MrProtocolData::MrPatData: PATMode, RefScanMode, AccelFactorPE, RefLinesPE
        // second argument 'true': AccelFactor3D and RefLines3D will be initialized as well
        fPATSetDefaultSeqLim(rSeqLim, true);

        // set to EXTRA mode: ref.lines are scanned separately for each slice by SBBPATRefScan (within SeqLoop)
        rSeqLim.setRefScanMode(SEQ::PAT_REF_SCAN_EXTRA);
    }
#endif

    //-------------------------------------------------------------------------------------
    // variant specific limits
    //-------------------------------------------------------------------------------------
    setVariantSpecificHardLimits(rSeqLim);


    //-------------------------------------------------------------------------------------
    // single shot epi sequences have problems with rephasing first PE-moment
    //-------------------------------------------------------------------------------------
#ifdef WIN32
    if (SeqUT.isUnitTestActive())
    {
        SeqUT.DisableTestCase(lGpFirstMomentNotRephasedErr, RTEB_ORIGIN_fSEQRunKernel, "single shot sequence, k-space center line not in first measured segment");
    }
#endif

#ifdef SUPPORT_PACE
    lStatus = PACE::fInit(&rSeqLim, SEQ::RESP_COMP_OFF, SEQ::RESP_COMP_BREATH_HOLD, SEQ::RESP_COMP_TRIGGER);
    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        SEQ_TRACE_ERROR.print("PACE::fInit() failed.");
        return lStatus;
    }
#endif

    // default TrueformC behaviour for all systems (is been overriden by pTX systems)
    // GUI logic activates the parameter on the System/Adjustments card for 3T systems only
    rSeqLim.setB1ShimMode(SEQ::TX_B1_SHIM_TRUEFORM, SEQ::TX_B1_SHIM_TRUEFORM_CP);

    // default case for excitation
    rSeqLim.setExcitationPulse(SEQ::EXCITATION_PULSE_STANDARD);

#if defined ZOOM_2DRF
    // ZOOMit is available on pTX and non-pTX systems
    if (SysProperties::isPTxSystem())
    {
        rSeqLim.setExcitationPulse(SEQ::EXCITATION_PULSE_STANDARD, SEQ::EXCITATION_ZOOMED);
        rSeqLim.setB1ShimMode(SEQ::TX_B1_SHIM_TRUEFORM, SEQ::TX_B1_SHIM_TRUEFORM_CP, SEQ::TX_B1_SHIM_PAT_SPEC, SEQ::TX_B1_SHIM_VOL_SEL);
        rSeqLim.setpTxVolProperty(MrProtocolData::PTxVolProp_B1Shim, MrProtocolData::PTxVolProp_Optimization);
        rSeqLim.setPTXPulses(1, 1, 1, 1); // index, min, max, inc, def
#ifdef ZOOM_EXTENDED
                                          // disable coupling of adjust volume to ptx volume
        rSeqLim.setAdjVolCoupling(MrProtocolData::AdjVolCoupling_ImagingVolume);

        // disable pulse acceleration for rotated ZOOMit
        rSeqLim.setPTXPulseAcceleration(0, 1.0, 1.0, 0.1, 1.0); // index, min, max, inc, def
#else
        rSeqLim.setAdjVolCoupling(MrProtocolData::AdjVolCoupling_ImagingVolume, MrProtocolData::AdjVolCoupling_PTxVolume);
        rSeqLim.setPTXPulseAcceleration(0, 1.0, 2.0, 0.1, 1.0); // index, min, max, inc, def
#endif
    }
#ifdef ZOOM_EXTENDED
    else
    {
        rSeqLim.setExcitationPulse(SEQ::EXCITATION_PULSE_STANDARD, SEQ::EXCITATION_ZOOMED);
        rSeqLim.setpTxVolProperty(MrProtocolData::PTxVolProp_Optimization);
        rSeqLim.setPTXPulses(1, 1, 1, 1);                       // index, min, max, inc, def
        rSeqLim.setPTXPulseAcceleration(0, 1.0, 1.0, 0.1, 1.0); // index, min, max, inc, def
    }
#endif // ZOOM_EXTENDED
#endif // ZOOM_2DRF

#ifndef EP2D_SE_MRE
    IterativeDenoisingUIParameter::initialize(rSeqLim);
#endif // !EP2D_SE_MRE
   
    if (SysProperties::is3TOrHigherSystem())
        rSeqLim.setStaticFieldCorrection(SEQ::OFF, SEQ::ON);

    //  ----------------------------------------------------------------------
    //  Declaration of pointer to UI parameter classes
    //  ----------------------------------------------------------------------
    lStatus = m_pUI->registerUI(rSeqLim);

    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        SEQ_TRACE_ERROR.print("Initialization of UI failed : 0x%lx", lStatus);
        return (lStatus);
    }

    return (lStatus);
}

NLSStatus SEQ_NAMESPACE::prepareError(const SeqLim& rSeqLim, std::string const& sErrorText)
{
    SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), sErrorText.c_str());

    return MRI_SEQ_SEQU_ERROR;
}

bool Ep2d::isNumberOfCoilsSufficient(MrProt& rMrProt, SeqLim& rSeqLim)
{
    // If coil selection is available, and the number of selected coil elements is smaller than the total
    // acceleration factor, inform the user with a popup. (do not do it during a test)
    if (!SysProperties::isSystemOperatingModeHostOnly() && rMrProt.rxSpec().isCoilSelectAvailable(0)
        && rSeqLim.isContextPrepForMeasurement() && SMSProperties::isSMS(rMrProt))
    {
        const MrRxCoilSelect CoilSelection = rMrProt.coilInfo().Meas().getaRxCoilSelectData()[0];

        if (rMrProt.getsSliceAcceleration().getlMultiBandFactor() * rMrProt.PAT().getlAccelFactPE()
            > CoilSelection.getNumOfUsedADCChan())
        {
            return false;
        }
    }

    return true;
}

void Ep2d::setPhaseCorrScansAndAdjustInitialDummyScans(MrProt& rMrProt)
{
    m_lPhaseCorrPrepScans = 0;
    if (rMrProt.preparationPulses().getlPhaseCorrectionMode() == MrProtocolData::PHASECORR_EXTERNAL)
    {
        m_lPhaseCorrPrepScans = 1;
        // If possible, reduce the number of initial preparing scans accordingly
        long numberOfMinimalPrepScans = 2;
        if (rMrProt.getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST)
            numberOfMinimalPrepScans = 0;

        if (m_lInitialDummyScans > numberOfMinimalPrepScans)
        {
            m_lInitialDummyScans--;
        }
    }
}

#ifdef EPI_SUPPORT_FREQ_FEEDBACK
long Ep2d::getCurrentVolumeForB0Correction()
{
    if (!fRTIsReadoutEnabled())
        return -1;
    
    // Increase volume counter for every first slice
    // Note: Only consider slices that are actually acquired
    // Note: On the Ice side (EPIPhaseCorrPEFunctor), the volumes are counted
    //       independently: sequence and Ice counting have to match!
    if (m_mySeqLoop.getlInnerSliceCounter() == 0)
    {
        m_lVolumeCounter++;
    }
    // Set current volume (counting starts at 0)
    long currentVolume = m_lVolumeCounter - 1;

    if (currentVolume < 0)
    {
        SEQ_TRACE_WARN.print("Volume counter %ld < 0 - this should never happen", currentVolume);
    }
    
    return currentVolume;
}
#endif

bool Ep2d::isSystemCompatibleWithFlavor() const
{
#if (defined COMPILE_EP2D_SE || defined ASL || defined BOLD || defined EP2D_MS)
    if (SysProperties::isFreeDot())
    {
        return false;
    }
#endif

    return true;
}

void Ep2d::setForcedFastGRERefscanWithoutPAT(bool isForcedFastGRERefScanWithoutPAT)
{
    m_isForcedFastGRERefscanWithoutPAT = isForcedFastGRERefScanWithoutPAT;
}

bool Ep2d::isForcedFastGRERefscanWithoutPAT() const
{
    return m_isForcedFastGRERefscanWithoutPAT;
}

NLSStatus SEQ_NAMESPACE::Ep2d::SeverePrepareErrorReturn(NLSStatus lStatus) const
{
    if (NLS_SEVERITY(lStatus) == NLS_SUCCESS)
        lStatus = MRI_SEQ_SEQU_ERROR; // make sure that an ERROR code is returned

    m_pUI->fEPIStdResetSolveHandlerControlTETITR(); // even if TE,TI or TR is changed, we can't help

    return (lStatus);
}

bool Ep2d::isFastGreRefScan(const MrProt& rMrProt)
{
    return rMrProt.getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST;
}

bool Ep2d::isGreRefScanType(const MrProt& rMrProt)
{
    return rMrProt.getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST || rMrProt.getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA;
}

uint16_t Ep2d::getFastGreRefScanBaseRes()
{
    return 24; // MZ: Temporary until we have a better possibility for switching
}

uint16_t Ep2d::getFastGreRefScanBandwidth()
{
    return 780; // MZ: Temporary until we have a better possibility for switching
}

long Ep2d::getMaxPATFactorForDLRecon() const
{
    return 3;
}

long Ep2d::getMinPATFactorForDLRecon() const
{
    return 2;
}


NLSStatus Ep2d::prePrepare(const MrProt& rMrProt, const SeqLim& rSeqLim, SeqExpo& rSeqExpo)
{
    MrProtFacade protFacade(rMrProt);

    // ==============================================================================
    // forbidden combinations for SMS
    // ==============================================================================

#if defined COMPILE_EP2D_DIFF || (defined BOLD && !defined PACE3D) || defined EP2D_MS
    if (protFacade.isSliceAcceleration())
    {
        const long lMultibandFactor = SMSProperties::getMultiBandFactor(rMrProt);
        // ---------------------------------------------------------------------------
        // Number of slices must be multiple of slice acceleration factor
        // ---------------------------------------------------------------------------
        if ((rMrProt.getsSliceArray().getlSize() + lMultibandFactor) % lMultibandFactor != 0)
            return prepareError(rSeqLim, "number of slices is NOT a multiple of slice acceleration factor");

        // ---------------------------------------------------------------------------
        // Number of concatenations must be 1
        // ---------------------------------------------------------------------------
        if (rMrProt.concatenations() != 1)
            return prepareError(rSeqLim, "number of concatenations must be 1 for slice acceleration");

        // ---------------------------------------------------------------------------
        // Number of slice groups must be 1
        // ---------------------------------------------------------------------------
        if (rMrProt.sliceGroupList().size() != 1)
            return prepareError(rSeqLim, "number of slice groups must be 1 for slice acceleration");

        // ---------------------------------------------------------------------------
        // Product of in-plane and slice acceleration must be MAX_MULTIBAND_FACTOR
        // ---------------------------------------------------------------------------
        if (lMultibandFactor * rMrProt.getsPat().getlAccelFactPE() > SMSProperties::MAX_MULTIBAND_FACTOR)
            return prepareError(rSeqLim, "the product of in-plane and slice acceleration must be <= MAX_MULTIBAND_FACTOR");

        // ---------------------------------------------------------------------------
        // GRE ref scan only allowed if PE acceleration factor > 1
        // ---------------------------------------------------------------------------
        if (rMrProt.getsPat().getlAccelFactPE() < 2 && rMrProt.getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA)
            return prepareError(rSeqLim, "GRE ref scan only allowed if iPAT > 1");

        // ---------------------------------------------------------------------------
        // Prohibit use of water excitation
        // ---------------------------------------------------------------------------
        if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_WaterExcitation)
            return prepareError(rSeqLim, "water excitation is not supported for slice acceleration");

        // ---------------------------------------------------------------------------
        // Prohibit certain diffusion modes
        // ---------------------------------------------------------------------------
        if (rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_ONE_SCAN_TRACE)
            return prepareError(rSeqLim, "diffusion mode DIFFMODE_ONE_SCAN_TRACE is not supported for slice acceleration");

#ifndef EP2D_MS
        // ---------------------------------------------------------------------------
        // Prohibit use of IR
        // ---------------------------------------------------------------------------
        if (rMrProt.getsPrepPulses().getucInversion() != SEQ::INVERSION_OFF)
            return prepareError(rSeqLim, "IR is not supported for slice acceleration");

        // ---------------------------------------------------------------------------
        // Prohibit use of min TE strategy
        // ---------------------------------------------------------------------------
        if (rMrProt.TOM() == SEQ::TOM_MINIMIZE_TE)
            return prepareError(rSeqLim, "min TE strategy is not supported for slice acceleration");

        // ---------------------------------------------------------------------------
        // With selective IR enabled, only fast FLASH reference scans are supported
        // ---------------------------------------------------------------------------
        // Note: In principle, SeqLoopEPI does support other reference scans for the
        //       combination SMS + IR (this may not apply to SeqLoopIIR, though).
        //       However, in this case IR modules get applied for the EPI reference
        //       scans as well, which requires corresponding adaptions to scan time
        //       and energy calculations.
        if (rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE
            && rMrProt.getsPat().getucRefScanMode() != SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST)
        {
            return prepareError(rSeqLim, "IR-preparation requires FLASH reference scans (only) for slice acceleration");
        }

                // ---------------------------------------------------------------------------
        // Number of slices must be an even multiple of the slice acceleration factor
        // in case of interleaved, slice-selective preparations
        // ---------------------------------------------------------------------------
        if (rMrProt.sliceSeries().mode() == SEQ::INTERLEAVED
            && (rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE
                || rMrProt.getsPrepPulses().getucSatRecovery() == SEQ::SATREC_SLICE_SELECTIVE))
        {
            if (rMrProt.getsSliceArray().getlSize() % (2 * lMultibandFactor) != 0)
                return prepareError(rSeqLim, "number of slices is NOT an even multiple of slice acceleration factor");
        }

#endif

        // ---------------------------------------------------------------------------
        // Prohibit ZOOMit
        // ---------------------------------------------------------------------------
        if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
            return prepareError(rSeqLim, "ZOOMit is not supported for slice acceleration");

        // ---------------------------------------------------------------------------
        // Prohibit PACE
        // ---------------------------------------------------------------------------
        if (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF)
            return prepareError(rSeqLim, "PACE functionality is not supported for slice acceleration");

        if (rMrProt.getucReconstructionPrio())
            return prepareError(rSeqLim, "Prio Recon is not supported for slice acceleration");


    }
#endif // COMPILE_EP2D_DIFF || (BOLD && !PACE3D) || defined EP2D_MS
    // ==============================================================================

    // ---------------------------------------------------------------------------
    // Only allow Abdomen/Thorax optimization with SPAIR
    // ---------------------------------------------------------------------------
#ifdef COMPILE_EP2D_DIFF
    if (!SysProperties::isLowField())
    {
        if ((rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_ABDOMEN)
            && !protFacade.isSPAIRFatSat())
            return prepareError(rSeqLim, "only SPAIR allowed with abdomen optimization");
    }
#else
    if ((rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_ABDOMEN)
        && !protFacade.isSPAIRFatSat())
        return prepareError(rSeqLim, "only SPAIR allowed with abdomen optimization");
#endif // COMPILE_EP2D_DIFF

    if ((rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_THORAX) && !protFacade.isSPAIRFatSat())
        return prepareError(rSeqLim, "only SPAIR allowed with thorax optimization");

    if ((rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_BREAST) && !protFacade.isSPAIRFatSat())
        return prepareError(rSeqLim, "only SPAIR allowed with breast optimization");

    // ---------------------------------------------------------------------------
    // Only allow Brain optimization with FatSat
    // ---------------------------------------------------------------------------
    if ((rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_BRAIN)
        && (rMrProt.preparationPulses().getlFatWaterContrast() != FatWaterContrast_FatSaturation))
        return prepareError(rSeqLim, "only FatSat allowed with brain optimization");

    // ---------------------------------------------------------------------------
    // POCS requires adaptive coil combination
    // ---------------------------------------------------------------------------
    if ((rMrProt.getsKSpace().getucPOCS() != SEQ::POCS_OFF)
        && (rMrProt.getucCoilCombineMode() != SEQ::COILCOMBINE_ADAPTIVE_COMBINE))
    {
        return prepareError(rSeqLim, "POCS requires adaptive coil combination");
    }

    // ---------------------------------------------------------------------------
    // check that dwell time is not less than hardware limit of 1000ns
    // ---------------------------------------------------------------------------
    for (long lI = 0; lI < rMrProt.contrasts(); lI++)
    {
        if (rMrProt.rxSpec().realDwellTime()[lI] < 1000)
            return prepareError(rSeqLim, "hardware sampling-rate limit exceeded");
    }

    // ---------------------------------------------------------------------------
    // Prohibit simultaneous use of IR and SPAIR
    // ---------------------------------------------------------------------------
    if (protFacade.isSPAIRFatSat() && !(rMrProt.preparationPulses().getucInversion() == SEQ::INVERSION_OFF))
        return prepareError(rSeqLim, "SPAIR and inversion preparation are incompatible");

    // ---------------------------------------------------------------------------
    // Prohibit IVIM and bipolar
    // ---------------------------------------------------------------------------
    if (protFacade.isIVIM() && (rMrProt.diffusion().getdsScheme() != SEQ::DIFFSCHEME_MONOPOLAR))
        return prepareError(rSeqLim, "IVIM is only allowed with monopolar diffusion mode");

    // UHF systems allow all B0 shim modes with ADJSHIM_TUNEUP as default. Since SliceAdj is not enabled by default, this would always fail.
    // ADJSHIM_TUNEUP should still be available for test purposes and the order of B0 shim modes should be identical for all sequences.
    if (!SysProperties::isUHFSystem())
    {
        // ---------------------------------------------------------------------------
        // tune up shim is only allowed if SliceAdj is enabled
        // ---------------------------------------------------------------------------
        if (rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_TUNEUP && !protFacade.isSliceAdj())
            return prepareError(rSeqLim, "Tuneup B0 shim is only allowed if SliceAdj is enabled");
    }

    // ---------------------------------------------------------------------------
    // whole body shim is only allowed if SliceAdj is enabled
    // ---------------------------------------------------------------------------
    if (rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_WHOLE_BODY && !protFacade.isSliceAdj())
        return prepareError(rSeqLim, "Whole body B0 shim is only allowed if SliceAdj is enabled");

    // ---------------------------------------------------------------------------
    // prohibit whole body B0 shim and local shim
    // ---------------------------------------------------------------------------
    if (rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_WHOLE_BODY && rMrProt.getsAdjData().getuiLocalShim() != LocalShim_OFF)
        return prepareError(rSeqLim, "Whole body B0 shim and local shim are incompatible");

    // ---------------------------------------------------------------------------
    // slice selective IR is only allowed with interleaved slice order
    // ---------------------------------------------------------------------------
    if (rMrProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE && rMrProt.sliceSeries().mode() != SEQ::INTERLEAVED)
        return prepareError(rSeqLim, "slice selective IR is only allowed with interleaved slice order");

    // ==============================================================================
    // forbidden combinations for SliceAdj
    // ==============================================================================
    if (protFacade.isSliceAdj())
    {
        // ---------------------------------------------------------------------------
        // Prohibit SliceAdj and slice acceleration
        // ---------------------------------------------------------------------------
        if (protFacade.isSliceAcceleration())
            return prepareError(rSeqLim, "SliceAdj and slice acceleration preparation are incompatible");

        // ---------------------------------------------------------------------------
        // Prohibit SliceAdj and ZOOMit
        // ---------------------------------------------------------------------------
        if (protFacade.isZOOMit())
            return prepareError(rSeqLim, "SliceAdj and ZOOMit are incompatible");

        // ---------------------------------------------------------------------------
        // Only allow tune up, standard or whole body B0 shim
        // ---------------------------------------------------------------------------
        if (!SysProperties::isUHFSystem())
        {
            if (!(rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_TUNEUP
                || rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_STANDARD
                || rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_WHOLE_BODY
                || rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_ABSOLUTE))
            {
                return prepareError(rSeqLim, "SliceAdj and B0 shim mode are incompatible");
            }
        }
    }

    // ---------------------------------------------------------------------------
    // CHARM 305893:
    //
    // TR should ALWAYS be the physically correct TR for single shot EPI. Also and
    // especially over multiple measurements. This leads to the following restrictions:
    //
    // - multi-slice mode must be interleaved
    //
    // - We do not allow repetition delay times, because this would disturb the
    //   steady state and change the true TR. This restriction can be overcome
    //   with a new parameter which allows to force SeqLoop to use a certain TR-fill
    //   at the end of the concatenation. So the measurement delay is included in
    //   TR.
    //
    // - For the same reason the measPause must be zero.
    //
    // - It is not allowed to use multiple concatenations and multiple measurements
    //   at the same time. Again this would lead to physically incorrect TRs for the
    //   first acquisition of a certain slice within the measurement for all measurements
    //   except the first, i.e. for all repetitions.
    // ---------------------------------------------------------------------------
    if (rMrProt.kSpace().getucMultiSliceMode() != SEQ::MSM_INTERLEAVED)
        return prepareError(rSeqLim, "rMrProt.kSpace().getucMultiSliceMode() == SEQ::MSM_INTERLEAVED required!");

    for (long lRepetition = 0; lRepetition < K_NO_REP_TIMES; lRepetition++)
    {
        if (rMrProt.repetitionDelayTime()[lRepetition])
            return prepareError(rSeqLim, "All rMrProt.repetitionDelayTime()[x] must be zero!");
    }

    if (rMrProt.measPause())
        return prepareError(rSeqLim, "rMrProt.measPause() must be zero!");

        // -------------------------------------------------------------------------------
        // Prohibit use of PTX when this is not defined by the sequence
        // -------------------------------------------------------------------------------
#ifndef ZOOM_2DRF
    if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
        return prepareError(rSeqLim, "PTX is not supported by this sequence");
#endif

#ifdef ASL
    // -------------------------------------------------------------------------------
    // In PCASL mode, labeling plain thickness is hard-coded and should not be modified
    // -------------------------------------------------------------------------------
    if ((rMrProt.getsGroupArray().getsPSat().getdThickness() != 10.0) && (rMrProt.getsAsl().getulMode() == SEQ::ASL_PSEUDOCASL))
        return prepareError(rSeqLim, "In PCASL mode, Sat-band thickness should not be modified.");
#endif


#ifdef EP2D_MS

    #ifdef COMPILE_EP2D_SE
    // ---------------------------------------------------------------------------
    // Prohibit simultaneous use of flow attenuation and multiple contrasts
    // ---------------------------------------------------------------------------
    if ((rMrProt.getsPrepPulses().getlFlowAttenuation() != MrProtocolData::FLOW_ATTENUATION_OFF)
        && (rMrProt.getlContrasts() > 1))
    {
        return prepareError(rSeqLim, "Flow attenuation is supported for a single contrast only");
    }

    // ---------------------------------------------------------------------------
    // Prohibit simultaneous use of IR-preparation and TE minimization
    // ---------------------------------------------------------------------------
    if ((rMrProt.getsPrepPulses().getucInversion() != SEQ::INVERSION_OFF) && (rMrProt.getlTOM() != SEQ::TOM_OFF))
    {
        return prepareError(rSeqLim, "IR preparation is not compatible with TE minimization");
    }
    #endif

        // ---------------------------------------------------------------------------
    // Prohibit simultaneous use of EPI reference scans and multiple segments
    // ---------------------------------------------------------------------------
    if ((rMrProt.getsPat().getucPATMode() != SEQ::PAT_MODE_NONE)
        && (rMrProt.getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI)
        && (rMrProt.getsFastImaging().getlSegments() > 1))
    {
        return prepareError(rSeqLim, "Multi-shot EPI requires FLASH reference scan");
    }

    // ---------------------------------------------------------------------------
    // Prohibit simultaneous use of EPI reference scans and multiple contrasts
    // ---------------------------------------------------------------------------
    if ((rMrProt.getsPat().getucPATMode() != SEQ::PAT_MODE_NONE)
        && (rMrProt.getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI) && (rMrProt.getlContrasts() > 1))
    {
        return prepareError(rSeqLim, "Multi-contrast EPI requires FLASH reference scan");
    }

    // ---------------------------------------------------------------------------
    // Non-standard reconstruction requires FLASH reference scans
    // ---------------------------------------------------------------------------
    if (isDLReconMode(rMrProt)
        && rMrProt.getsPat().getucPATMode() != SEQ::PAT_MODE_NONE
        && rMrProt.getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI)
    {
        return prepareError(rSeqLim, "Non-standard reconstructions require FLASH reference scans");
    }

    // ---------------------------------------------------------------------------
    // Non-standard reconstruction handles coil compression internally
    // ---------------------------------------------------------------------------
    if (isDLReconMode(rMrProt) && rMrProt.getsChannelMatrix().getucChannelMixingMode() != SEQ::ChannelMixingMode::CMM_OFF)
    {
        return prepareError(rSeqLim, "Non-standard reconstruction requires matrix optimization off");
    }

    // ---------------------------------------------------------------------------
    // Non-standard reconstruction requires adaptive coil combination
    // ---------------------------------------------------------------------------
    // This ensures the availability of CoilSensServerFunctor in the pipeline
    if (isDLReconMode(rMrProt)
        && (rMrProt.getucCoilCombineMode() != SEQ::COILCOMBINE_ADAPTIVE_COMBINE))
    {
        return prepareError(rSeqLim, "Non-standard reconstructions require adaptive coil combination");
    }

    // ---------------------------------------------------------------------------
    // Non-standard reconstruction requires some sort of PAT
    // ---------------------------------------------------------------------------
    if (isDLReconMode(rMrProt) && (rMrProt.getsPat().getucPATMode() == SEQ::PAT_MODE_NONE))
    {
        return prepareError(rSeqLim, "Non-standard reconstructions requires active PAT");
    }

    // ---------------------------------------------------------------------------
    // Averages not supported
    // ---------------------------------------------------------------------------
    if (isDLReconMode(rMrProt) && (rMrProt.getlAverages()>1))
    {
        return prepareError(rSeqLim, "Non-standard reconstructions does not support multiple averages.");
    }

    // ---------------------------------------------------------------------------
    // Non-standard reconstruction requires PAT factors 2, 3, or 4
    // ---------------------------------------------------------------------------
    const long minPATFactorForDL = getMinPATFactorForDLRecon();
    const long maxPATFactorForDL = getMaxPATFactorForDLRecon();
    const long currentPATFactor  = static_cast<long>(rMrProt.getsPat().getlAccelFactPE());

    if (isDLReconMode(rMrProt) && (currentPATFactor < minPATFactorForDL || currentPATFactor > maxPATFactorForDL))
    {
        SEQ_TRACE_ERROR_COND(
            !rSeqLim.isContextPrepForBinarySearch(),
            "Non-standard reconstructions requires PAT factors between %ld and %ld, current PAT factor it's %ld",
            minPATFactorForDL,
            maxPATFactorForDL,
            currentPATFactor);

        return MRI_SEQ_SEQU_ERROR;
    }

    // ---------------------------------------------------------------------------
    // check for k-space sampling limitations (adapted from DefaultKSpace.cpp)
    // ---------------------------------------------------------------------------
    // Note: Ideally, every parameter change affecting k-space coverage would
    //       go along with corresponding consistency checks and required
    //       adaptions (as realized within DefaultKSpace.cpp). However, the
    //       combination of partial Fourier, PAT and segmentation currently
    //       leads to inconsistencies of UI parameters.
    //       As a workaround, we mimic the _checkPhOSAndPEFTLen function
    //       here in order to explicitly trigger corresponding solve handlers.
    if (!KSpace::KSpaceHelpers::checkPhOSAndPEFTLen(rMrProt, rSeqLim))
    {
        SEQ_TRACE_ERR << "checkPhOSAndPEFTLen failed!";
        return MRI_SEQ_SEQU_ERROR;
    }

#endif

    if (!IterativeDenoisingUIParameter::ProtocolOkForIterativeDenoising(rMrProt, rSeqLim))
    {
        return MRI_SEQ_SEQU_ERROR;
    }

    return MRI_SEQ_SEQU_NORMAL;
}

//  --------------------------------------------------------------------------
//
//  Name        :  Ep2d::prepare
//
//  Description :
/// \brief         Preparation of the sequence during binary search and prior
///                 to sequence execution
//
//  Return      :  NLS status
//
//  --------------------------------------------------------------------------
NLSStatus Ep2d::prepare(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL; // * Return status *

    //   ----------------------------------------------------------------------
    //
    //   Preparation of sequence
    //
    //   ----------------------------------------------------------------------
#undef DEBUG_ORIGIN
#define DEBUG_ORIGIN 0x00000200

    MrProtocolData::SeqExpoRFInfo RFInfo;
    MrRxCoilSelect                CoilSelection;
    double                        dBandwidthPerPixelPE         = 0.0;
    long                          lEffectiveEchoSpacing        = 0;
    long                          lDefaultNoOfRXChannels       = 256;
    long                          lNeededTE                    = -1;
    long                          lNeededTI                    = -1;
    long                          lNeededTR                    = -1;
    long                          lRequiredPrepScans           = 0; // Total number of prepscans (including PAT reference, adjust, etc.)
    long                          lRequiredPrepScansWithoutADC = 0; // number of prepscans withot ADC event
    long                          lIceUserCtrlMask             = 0;
    bool                          bSuccess                     = 0;

    // pointer for MrProtFacade for easier protocol queries
    MrProtFacade protFacade(rMrProt);

#ifdef ASL
    // check if TE, TR and TI are in the same precision as the UI handlers
    // if not, let prepare fail so that protocol conversion is enforced
    if (static_cast<int>(rMrProt.te()[0] / 100.0 + 0.5) * 100 != rMrProt.te()[0])
    {
        SEQ_TRACE_ERROR.print("TE in protocol has incorrect precision, protocol has to be converted");
        return MRI_SEQ_SEQU_ERROR;
    }
    if (static_cast<int>(rMrProt.tr()[0] / 1000.0 + 0.5) * 1000 != rMrProt.tr()[0])
    {
        SEQ_TRACE_ERROR.print("TR in protocol has incorrect precision, protocol has to be converted");
        return MRI_SEQ_SEQU_ERROR;
    }
    if (static_cast<int>(rMrProt.ti()[0] / 1000.0 + 0.5) * 1000 != rMrProt.ti()[0])
    {
        SEQ_TRACE_ERROR.print("TI in protocol has incorrect precision, protocol has to be converted");
        return MRI_SEQ_SEQU_ERROR;
    }
#endif // ASL

    // check member pointer
    if (m_pUI == nullptr)
    {
        SEQ_TRACE_ERROR.print("ERROR: m_pUI == nullptr");
        return MRI_SEQ_SEQU_ERROR;
    }

    // load SliceAdj values from ini file (if configured) and set all m_sSliceAdjParametersRequestedBySequence params here
    if (protFacade.isSliceAdj())
    {
        loadSliceAdjData(rSeqLim, rMrProt);

        m_sSliceAdjParametersRequestedBySequence.setAdjFre(rMrProt.getsAdjData().getuiAdjSliceBySliceFrequency());
        m_sSliceAdjParametersRequestedBySequence.setAdjShim(rMrProt.getsAdjData().getuiAdjSliceBySliceFirstOrderShim());
        m_sSliceAdjParametersRequestedBySequence.setAdjTra(rMrProt.getsAdjData().getuiAdjSliceBySliceTxRef());
        if (SysProperties::isPTxSystem())
            m_sSliceAdjParametersRequestedBySequence.setAdjTxScale(rMrProt.getsAdjData().getuiAdjSliceBySlicePtx());
        else
            m_sSliceAdjParametersRequestedBySequence.setAdjTxScale(false);
    }
    else
    {
        m_sSliceAdjParametersRequestedBySequence.setAdjFre(false);
        m_sSliceAdjParametersRequestedBySequence.setAdjShim(false);
        m_sSliceAdjParametersRequestedBySequence.setAdjTra(false);
        m_sSliceAdjParametersRequestedBySequence.setAdjTxScale(false);
    }

#ifdef EP2D_MS 
    if (isDLReconMode(rMrProt))
    {
        // DL reconstruction requires noise adjust even for no PAT
        // => enable noise adjust (UI handlers present, this here is just for failsafety)
        rMrProt.setucEnableNoiseAdjust(true);
    }
#endif

#ifdef ZOOM_2DRF
    // determine some parameters from the *first* PTXRFPulse (index 0)
    const size_t                      uiCurVol        = 0;
    long                              lPTXRFPulseSize = rMrProt.getsTXSPEC().getaPTXRFPulse().size();
    MrProtocolData::PTXTrajectoryType eExcType        = MrProtocolData::PTXTrajectoryType_EPI_1D;

    if (lPTXRFPulseSize > 0)
    {
        eExcType = rMrProt.getsTXSPEC().getaPTXRFPulse()[uiCurVol].getlTrajectoryType();
    }

#ifdef ZOOM_EXTENDED
    // to enable support of pre-VA20 protocols where lPTXRFPulseSize can be equal to 0,
    // a protocol conversion is forced to set valid default values for getaPTXRFPulse()[0]
    if (lPTXRFPulseSize == 0)
    {
        return MRI_SEQ_SEQU_ERROR;
    }
#endif

#endif

    // Variables for simultaneous multislice imaging
    long                          lNMeasPrepScanTR               = 0;
    long                          lNMeasTotal                    = 0;
    long                          lNMeasImagingScanTR            = 0;
    long                          lNMeasImagingScanTR_SecondMeas = 0;
    MrProtocolData::SeqExpoRFInfo RFInfoSeqLoopPerSlice;
    long                          lTR_SliceAccPrepMeas = 0;
    long                          lCoolPauseExplicit   = 0;
    long                          lScanTimeSatsEtc     = 0;
    long                          lScanTimeBasic       = 0;
    long                          lReducedSlices       = 0;
    long                          lTotalSlices         = 0;

#ifdef EP2D_SE_MRE
    // Check validity of MRE parameters
    if (m_dMEGFrequency <= 0.0)
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "m_dMEGFrequency <= 0.0");

        return MRI_SEQ_SEQU_ERROR;
    }
    const long   lTR_us                             = rMrProt.tr()[0];
    const long   lCyclesPerTR                       = static_cast<long>((1.0E-6 * m_dMEGFrequency) * static_cast<double>(lTR_us) + 0.5);
    const long   lNominalDriverBurstDuration_cycles = m_bIsDefaultMode ? 3 : lCyclesPerTR;
    const double dDriverBurstDuration_us            = 1.0E6 * static_cast<double>(lNominalDriverBurstDuration_cycles) / m_dMEGFrequency;
    const long   lDriverBurstDuration_us            = static_cast<long>(dDriverBurstDuration_us + 0.5);

    const int32_t nSlc      = rMrProt.sliceSeries().size();
    const long    lEffNoSlc = protFacade.isSliceAcceleration() ? nSlc / SMSProperties::getMultiBandFactor(rMrProt) : nSlc;
    if (lEffNoSlc <= 0)
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "lEffNoSlc <= 0");
        return MRI_SEQ_SEQU_ERROR;
    }
    const long lSliceExcitationSpacing_us = lTR_us / lEffNoSlc;
    if (lSliceExcitationSpacing_us % lDriverBurstDuration_us != 0)
    {
       SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "[defaultMode] lSliceExcitationSpacing_us %% lDriverBurstDuration_us != 0");
        return MRI_SEQ_SEQU_ERROR;
    }

    // For inline ICE, the matrix size has to be 256, but we want flexible acquisition matrix sizes
    // -> choose float interpolation factor accordingly
    // - Regular 2D interpolation must be off (see limits)
    const bool bInterpolation2D = rMrProt.kSpace().Interpolation2D();
    if (bInterpolation2D)
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "[Inline ICE] Regular 2D interpolation is on");
        return MRI_SEQ_SEQU_ERROR;
    }
    const long  lBaseMatrix           = rMrProt.kSpace().baseResolution();
    const float flInterpolationFactor = 256.0f / static_cast<float>(lBaseMatrix);
    if ((flInterpolationFactor < 1.0f) || (flInterpolationFactor > 4.0f))
    {
SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "[Inline ICE] necessary interpolation factor (%f) out of range", flInterpolationFactor);
        return MRI_SEQ_SEQU_ERROR;
    }
    rMrProt.kSpace().setfl2DInterpolation(flInterpolationFactor);

    // Set the other SBB's members for MRE
    m_pSBBRefocSE = m_EPIKernel.getpSBBRefoc();
    if (m_pSBBRefocSE == nullptr)
    {
        SEQ_TRACE_ERROR.print("ERROR: m_pSBBRefocSE == nullptr");
        return MRI_SEQ_SEQU_ERROR;
    }
    // MEG frequency
    m_pSBBRefocSE->m_dMEGFrequencyHz = m_dMEGFrequency;
    // Configure running trigger
    m_pSBBRefocSE->m_lExtTriggerDuration  = EXT_TRIGGER_DURATION_US;
    m_pSBBRefocSE->m_lFirstExtTriggerTime = 0;
    m_pSBBRefocSE->m_lExtTriggerSpacing   = lDriverBurstDuration_us;

    m_EPIKernel.setbDefaultMode(m_bIsDefaultMode);
    // - initialize kernel phase offset counter here in prepare
    m_EPIKernel.setlMEGPhaseOffsetCounter(0);
    // MEG frequency
    m_EPIKernel.m_dMEGFrequencyHz = m_dMEGFrequency;
    // Configure running trigger
    m_EPIKernel.m_lExtTriggerDuration  = EXT_TRIGGER_DURATION_US;
    m_EPIKernel.m_lFirstExtTriggerTime = 0;
    m_EPIKernel.m_lExtTriggerSpacing   = lDriverBurstDuration_us;

    m_mySeqLoop.setlExtTriggerSpacing(lDriverBurstDuration_us);

    // Free loop controls both phase offset and MEG polarity
    const long lfreeLoopLength_seqLoop = NMAXMEGS * MEG_DEFAULT_TRIGGER_STEPS;
    m_mySeqLoop.setFreeLoopLength(lfreeLoopLength_seqLoop);
    m_EPIKernel.setlNumMEGPhaseOffsets(MEG_DEFAULT_TRIGGER_STEPS);

    // for visual distinction between the two regions in the graphical limits
    const long gapBetweenFract_And_NonFract_us = 1000;

    const long lFractionalEncPerc = static_cast<long>(MRE_ALLOWED_FRACTIONAL_FACTOR * 100.0);

    // for MRE application, single slice group is more appropriate.So we forbidden the user to add the slice group number
    if (rMrProt.sliceGroupList().size() > 1)
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "MRE sequences only for one slice group");
        return MRI_SEQ_SEQU_ERROR;
    }

    // The MRE processing will not know what to do with the uncombined images, and only(combined) magnitude and phase difference images are allowed.
    // so we have to make sure this is not selected for both SoS and ACC.
    if (rMrProt.uncombImages() == true)
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "'save uncombined images' is not allowed in MRE sequence");
        return MRI_SEQ_SEQU_ERROR;
    }

    // check timing, The total time of RSATs and CSATs should be less than trigger spacing 50ms, short enough to not need external trigger event inside.
    long lPrepPulseTime = 0; // us
    if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_FatSaturation)
    {
        // Consider duration of fat saturation module
        lPrepPulseTime += fSDSRoundUpGRT(m_mySeqLoop.getlFSDuration());
    }

    for (long lI = 0; lI < rMrProt.satList().size(); ++lI)
    {
        // Consider duration of RSat module
        lPrepPulseTime += fSDSRoundUpGRT(m_mySeqLoop.getlRSatDuration(lI));
    }

    if (lPrepPulseTime > lDriverBurstDuration_us)
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "The total time of RSATs and CSATs is %ld us, longer than 50000us is not allowed in MRE", lPrepPulseTime);
        return MRI_SEQ_SEQU_ERROR;
    }

    // It is used to transport the info about “I am an MRE sequence” to the AddIn
    rMrProt.applicationDetails(SEQ::APPLICATION_MRE);

#endif

    // ---------------------------------------------------------------------------
    // initialize UI class
    // ---------------------------------------------------------------------------
    lStatus = m_pUI->initializeUI(rMrProt, rSeqLim);

    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        // SEQ_TRACE_ERROR.print("%s : Initialization of UI failed : 0x%lx"  , lStatus);
        return (lStatus);
    }

    // --------------------------------------------------------------------------------------------
    // Pass all slice acceleration parameters to underline utility classes together
    // and patch all Unit Test special cases
    // --------------------------------------------------------------------------------------------
    if (protFacade.isSliceAcceleration())
    {
        if (!m_mySeqLoop.setMultibandFactor(SMSProperties::getMultiBandFactor(rMrProt)))
            return MRI_SEQ_SEQU_ERROR;

        if (isIIRSchemeStandard(rMrProt, rSeqLim))
        {
            // Account for trigger halt duration in case of triggered SMS acquisition
            m_mySeqLoop.setPutTriggerDelayIntoTR(SMSProperties::isSMS(rMrProt));
        }

        // ---------------------------------------------------------------------------
        // configure SBBEPIKernel
        // ---------------------------------------------------------------------------
        // Setup slice acceleration properties of SBBEPIKernel
        if (!m_EPIKernel.setMultiBandFactor(SMSProperties::getMultiBandFactor(rMrProt)))
            return MRI_SEQ_SEQU_ERROR;

        if (!m_EPIKernel.setFOVShiftFactor(SMSProperties::getFOVShiftFactor(rMrProt)))
            return MRI_SEQ_SEQU_ERROR;
    }

#ifdef BOLD
    // enable real-time processing mode by default for some systems
    m_bIsRealtimeProcessingEnabled = m_debugSettings.getDefaultSetting<bool>("EP2D_BOLD/EnableRealtimeProcessing", SysProperties::isUHFSystem() ? true : false);
#else
    m_bIsRealtimeProcessingEnabled = m_debugSettings.getDefaultSetting<bool>("EP2D_BOLD/EnableRealtimeProcessing", false);
#endif
    SEQ_TRACE_DEBUG_COND(!rSeqLim.isContextPrepForBinarySearch(), "m_bIsRealtimeProcessingEnabled: %d", m_bIsRealtimeProcessingEnabled);

    // disable B0 drift correction if frequalizer is active for some systems
    m_bB0Correction = m_debugSettings.getDefaultSetting<bool>(
        "EPI_GENERAL/UseB0Correction",
        SysProperties::isUHFSystem() && SysProperties::bIsTempToFreqB0CorrectionEnabled() ? false : true);

    // -------------------------------------------------------------------------------
    // Configure B0 correction
    // -------------------------------------------------------------------------------
    // B0 correction requires internal phase correction scans
    m_bB0Correction = m_bB0Correction && (rMrProt.preparationPulses().getlPhaseCorrectionMode() != MrProtocolData::PHASECORR_EXTERNAL);

    // -------------------------------------------------------------------------------
    // Prohibit use of B0 correction with incompatible loop structure
    // -------------------------------------------------------------------------------
    // B0 correction requires innermost acquisition of complete volumes.
    // This is not the case for
    // - respiration compensation (concatenations loop outside diffusion loop)
    // - multiple concatenations (if long TR triggering mode is inactive)
    m_bB0Correction = m_bB0Correction && isLoopStructureCompatibleWithB0Correction(rMrProt);
    SEQ_TRACE_DEBUG_COND(!rSeqLim.isContextPrepForBinarySearch(), "m_bB0Correction: %d", m_bB0Correction);
    m_bSequentialVolumeAcquisition = isLoopStructureCompatibleWithB0Correction(rMrProt);
#ifdef WIN32
#if defined BOLD && !defined PACE3D
    // Set B1 control loop parameters for BOLD
    float fCorrectionFactorMax = 1.0f;
    float fPeakReserveFactor   = 0.0f;

    // get B1 control loop structure
    MrTXSpecData&               rTXSpec                 = rMrProt.getsTXSPEC();
    B1CorrectionParametersType& rB1CorrectionParameters = rTXSpec.getB1CorrectionParameters();

    // Has the B1 control loop been deactivated in the INI file?
    bool bB1ControlLoopActive = m_debugSettings.getDefaultSetting<bool>("EP2D_BOLD/is_B1_control_active", true);

    // If not, check for internal activation condition
    if (bB1ControlLoopActive)
        bB1ControlLoopActive = useB1ControlLoop(rMrProt, fCorrectionFactorMax, fPeakReserveFactor);

    // retrieve information if and how the control loop is used
    rB1CorrectionParameters.setbActive(bB1ControlLoopActive);
    if (bB1ControlLoopActive)
    {
        rB1CorrectionParameters.setflCorrectionFactorMax(fCorrectionFactorMax);
        rB1CorrectionParameters.setflPeakReserveFactor(fPeakReserveFactor);
        rB1CorrectionParameters.setbValid(true); // this flag has to be set true, otherwise the values above will be set to the values from the SeqLims later on
    }

#endif
#endif

    // ---------------------------------------------------------------------------
    // reset solve handler control
    // ---------------------------------------------------------------------------
    m_pUI->fEPIStdResetSolveHandlerControlTETITR();

    // ---------------------------------------------------------------------------
    // set default value for maximum numnber of receiver channels
    // ---------------------------------------------------------------------------
    rSeqExpo.setMaxReceiverChannels(static_cast<int32_t>(lDefaultNoOfRXChannels));

    // ---------------------------------------------------------------------------
    // Calculate the rotation matrices and offsets for slices
    // ---------------------------------------------------------------------------
    lStatus = fSUPrepSlicePosArray(rMrProt, rSeqLim, m_asSLC);
    
    // error can only be caused by programming error
    // => trace also if rSeqLim.isContextPrepForBinarySearch()
    CheckStatusPB(lStatus, "fSUPrepSlicePosArray");

    // set 3D phase offsets for SMS cases
    if (protFacade.isSliceAcceleration())
    {
        prepSliceAccelPhaseOffcenter3D(rMrProt, m_asSLC);
    }

    // ---------------------------------------------------------------------------
    // Specify type of reference scan for PAT
    //
    // For PAT factors greater than 2 a segmented reference scan is used
    // (this behavior does not hold for the Flash reference scan). The
    // primary reason for this is that the same timing is currently used for the
    // PAT reference scans as for the imaging scans. Consequently with a high PAT
    // factor there is a restriction on the number of reference lines that can be
    // acquired with a single-shot reference scan.
    // ---------------------------------------------------------------------------
#ifdef SUPPORT_iPAT_a_ep2d
    m_bSegmentedRefLines = isSegmentedPATRefLinesCondition(rMrProt);
#endif


#ifdef SUPPORT_iPAT_TGSE
    {
        if (m_pUI->calculateKspace(rMrProt) == MRI_SEQ_SEQU_ERROR)
        {
            lStatus = MRI_SEQ_SEQU_ERROR;
            return SeverePrepareErrorReturn(lStatus);
        }

        // * ---------------------------------------------------------------------- *
        // * iPAT: switch on/off SBBPATRefScan for Ref.Scan Mode 'extra'            *
        // * ---------------------------------------------------------------------- *
        if (rMrProt.PAT().getucPATMode() != SEQ::PAT_MODE_NONE)
        {
            // switch on/off SBBPATRefScan for RefScanMode 'extra'
            if (isGreRefScanType(rMrProt))
            {
                // advice SeqLoop to run SBBPATRefScan once (i.e. only for first measurement)
                m_mySeqLoop.setePATRefScanLoopMode(SeqLoop::PATRefScanLoopMode_ONCE);
                m_mySeqLoop.setRepetitionValueForMdh(0);
            }
            else
            {
                // inplace mode => no SBBPATRefScan from SeqLoop required
                m_mySeqLoop.setePATRefScanLoopMode(SeqLoop::PATRefScanLoopMode_NEVER);
            }
        }
        else
        {
            m_mySeqLoop.setePATRefScanLoopMode(SeqLoop::PATRefScanLoopMode_NEVER);
        }
    }
#endif // SUPPORT_iPAT_TGSE

    // ---------------------------------------------------------------------------
    // calculate EPI reordering
    // note: in case of changes, please adapt them in EpCommonUINS::adaptRefLinesPE accordingly
    // ---------------------------------------------------------------------------
    
    // configure reorder info structure
    configureReorderInfo(rMrProt);

    // prepare reordering
    if (!m_REOInfo.reorderEPI(rMrProt, rSeqLim))
        mSBBErrGotoFinish(m_REOInfo, "m_REOInfo.reorderEPI");

        // ---------------------------------------------------------------------------
        // segmented EPI (phase direction)
        // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Set gradient performance for EPI readout
    //
    // UltraFast gradient mode: set slew rate to "UltraFast" value specified
    //                          in MeasPerm section.
    //
    // Fast gradient mode:      set slew rate to 95% of "fast" value specified
    //                          in MeasPerm section (CHARM 324468).
    //
    // Normal gradient mode:    set slew rate to 80% of "fast" value specified
    //                          in MeasPerm section (CHARM 333764);
    //                          this option is provided to reduce Nyquist ghosts
    //                          under low supply voltage conditions.
    // ---------------------------------------------------------------------------
    switch (rMrProt.gradSpec().gradModeBeforeGSWD())
    {
    case SEQ::GRAD_ULTRAFAST:

        // update maximum gradient amplitude used in GPA balance checks
        // (only for sequences that support GPA balance model)
#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE || defined EP2D_MS 

        m_EPIKernel.setUseGPABalance(true, SysProperties::getGradMaxAmpl(SEQ::GRAD_ULTRAFAST));

#endif

        // set scaling factor for rise time
        // note that this call also sets the default gradient performance (SBBEPIReadout)
        m_EPIKernel.setMinRiseTimeScaleFactor(1.0 + DBL_EPSILON);

        break;

    case SEQ::GRAD_FAST:

        // update maximum gradient amplitude used in GPA balance checks
        // (only for sequences that support GPA balance model)
#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE || defined EP2D_MS 

        m_EPIKernel.setUseGPABalance(true, SysProperties::getGradMaxAmpl(SEQ::GRAD_FAST));

#endif

        // set scaling factor for rise time
        // note that this call also sets the default gradient performance (SBBEPIReadout)
        m_EPIKernel.setMinRiseTimeScaleFactor(1.05);

        break;

    case SEQ::GRAD_NORMAL:

        // update maximum gradient amplitude used in GPA balance checks
        // (only for sequences that support GPA balance model)
#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE || defined EP2D_MS

        m_EPIKernel.setUseGPABalance(true, SysProperties::getGradMaxAmpl(SEQ::GRAD_FAST));

#endif

        // set scaling factor for rise time
        // note that this call also sets the default gradient performance (SBBEPIReadout)
        m_EPIKernel.setMinRiseTimeScaleFactor(1.25);
        break;

    default:

        SEQ_TRACE_ERROR.print(
            "Gradient mode %d not supported by EPI sequence", rMrProt.gradSpec().gradModeBeforeGSWD());
        return MRI_SEQ_SEQU_ERROR;
    }

    // ---------------------------------------------------------------------------
    // configure SBBEPIKernel
    // ---------------------------------------------------------------------------
    lStatus = configureEPIKernel(rSeqLim, rMrProt);

    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
        return SeverePrepareErrorReturn(lStatus);


#ifdef SUPPORT_iPAT_TGSE
    {
        // forbid phasePartialFourierFactor in non-iPAT mode (recon/reordering issue)
        if ((rMrProt.PAT().getucPATMode() == SEQ::PAT_MODE_NONE) && (rMrProt.kSpace().phasePartialFourierFactor() != SEQ::PF_OFF))
        {
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "For non-iPAT phasePartialFourierFactor is not possible! ");
            return MRI_SEQ_SEQU_ERROR;
        }

        // forbid slicePartialFourierFactor in non-iPAT mode (recon/reordering issue)
        if ((rMrProt.PAT().getucPATMode() == SEQ::PAT_MODE_NONE) && (rMrProt.kSpace().slicePartialFourierFactor() != SEQ::PF_OFF))
        {
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "For non-iPAT slicePartialFourierFactor is not possible! ");
            return MRI_SEQ_SEQU_ERROR;
        }

        // forbid rectangular FOV in non-iPAT mode (backwards compatibility)
        if ((rMrProt.PAT().getucPATMode() == SEQ::PAT_MODE_NONE) && (rMrProt.getsSliceArray().getasSlice()[0].getdPhaseFOV() != rMrProt.getsSliceArray().getasSlice()[0].getdReadoutFOV()))
        {
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "For non-iPAT rectangular FOV is not possible! ");
            return MRI_SEQ_SEQU_ERROR;
        }

        //---------------------------------------------------------------------------
        // call fPATPrepPost, which checks some MrProtocolData::MrPatData related restrictions
        //---------------------------------------------------------------------------
        lStatus = fPATPrepPost(rMrProt, rSeqLim, rSeqExpo, &m_REOInfo);

        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            CheckStatusPB(lStatus, "fPATPrepPost");
        }
        else
        {
            CheckStatusB(lStatus);
        }
    }
#endif

#ifdef EP2D_SE_MRE
    // Fractional: mimic the solution of greMRE
    // - see MrServers\MrImaging\seq\a_ep_seg_therm\a_ep_seg_therm_UI.cpp
    //   and MrServers\MrImaging\seq\a_ep_seg_therm\a_ep_seg_therm.cpp
    // 1. Prepare the kernel with and without fractional MEG
    // 2. Log minTE for both cases
    // 3. Calculate the upper limit for fractional -> "gap" in UI TE limit range
    m_pSBBRefocSE->m_lMEGFractionalEncPerc = 100;
    m_EPIKernel.setWantedTE(0);
    // Don't mSBBErrGotoFinish(m_EPIKernel, "m_EPIKernel.prep");
    // -> We need values for m_lMinTE_NonFractional_us, m_lMinTE_Fractional_us
    if (!m_EPIKernel.prep(rMrProt, rSeqLim, rSeqExpo) && !rSeqLim.isContextPrepForBinarySearch())
    {
        SEQ_TRACE_ERROR.print("m_EPIKernel.prep(...) failed for non-fractional MEG");
    }
    m_lMinTE_NonFractional_us = m_EPIKernel.getNeededTE();

    m_pSBBRefocSE->m_lMEGFractionalEncPerc = lFractionalEncPerc;
    m_EPIKernel.setWantedTE(0);
    if (!m_EPIKernel.prep(rMrProt, rSeqLim, rSeqExpo) && !rSeqLim.isContextPrepForBinarySearch())
    {
        SEQ_TRACE_ERROR.print("m_EPIKernel.prep(...) failed for fractional MEG of %ld%%", lFractionalEncPerc);
    }
    m_lMinTE_Fractional_us = m_EPIKernel.getNeededTE();

    m_EPIKernel.setWantedTE(rMrProt.te()[0]);

    m_lFractionalTE_LimitsRange_us = std::max(0l, m_lMinTE_NonFractional_us - gapBetweenFract_And_NonFract_us - m_lMinTE_Fractional_us);

    if (rMrProt.te()[0] > (m_lMinTE_Fractional_us + m_lFractionalTE_LimitsRange_us))
    {
        m_pSBBRefocSE->m_lMEGFractionalEncPerc = 100;
        if (!m_EPIKernel.prep(rMrProt, rSeqLim, rSeqExpo))
            mSBBErrGotoFinish(m_EPIKernel, "m_EPIKernel.prep");
    }
#endif

    // ---------------------------------------------------------------------------
    // pass dynamic adjustment data to all SBBs
    // ---------------------------------------------------------------------------
    m_EPIKernel.setsSliceAdjParametersRequestedBySequence(m_sSliceAdjParametersRequestedBySequence);

    // ---------------------------------------------------------------------------
    // prepare SBBEPIKernel
    // ---------------------------------------------------------------------------
    if (!m_EPIKernel.prep(rMrProt, rSeqLim, rSeqExpo))
        mSBBErrGotoFinish(m_EPIKernel, "m_EPIKernel.prep");

    // --------------------------------------------------------------------------
    // Export gradient timing to UI
    // --------------------------------------------------------------------------
#ifdef WIN32
    if (!exportDiffusionTimingToUI(rSeqLim, rMrProt))
        return MRI_SEQ_SEQU_ERROR;
#endif

    // ---------------------------------------------------------------------------
    // check gradients of the EPI-kernel for grad spec violations
    // ---------------------------------------------------------------------------
    if (!m_EPIKernel.checkGradients(rMrProt, rSeqLim))
    {
        lStatus = m_EPIKernel.getNLSStatus();
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "m_EPIKernel.checkGradients failed: 0x%lx", m_EPIKernel.getNLSStatus());
        return SeverePrepareErrorReturn(lStatus);
    }

#ifdef ZOOM_2DRF
    if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
    {
        MrGradSpec sMrGradSpec(rMrProt.getsGRADSPEC());
        if (sMrGradSpec.isGSWDMode() && !rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ALWAYS.print(
                "ZOOM_2DRF: rt=%.2f/%.2f, g=%.2f/%.2f (isGSWDMode rt=%.2f)\n",
                m_EPIKernel.getMinRiseTime(rMrProt.gradSpec().mode(), SBBEPIKernel_EPIRO_GRAD_PERF_BLIPS),
                m_EPIKernel.getMinRiseTime(rMrProt.gradSpec().mode(), SBBEPIKernel_EPIRO_GRAD_PERF_RO),
                m_EPIKernel.getMaxMagnitude(rMrProt.gradSpec().mode(), SBBEPIKernel_EPIRO_GRAD_PERF_BLIPS),
                m_EPIKernel.getMaxMagnitude(rMrProt.gradSpec().mode(), SBBEPIKernel_EPIRO_GRAD_PERF_RO),
                rMrProt.getsGRADSPEC().getflGSWDMinRiseTime());
        }

#ifdef WIN32
        if (SeqUT.isUnitTestActive())
        {
            setZoomItSeqUTExpectations(rMrProt);

        }
#endif
    }
#endif

    // ---------------------------------------------------------------------------
    // check TE and set it to minimum if needed
    // ---------------------------------------------------------------------------
    // Note: If TOM is enabled or if it is not possible to apply the requested
    // b-value with the current protocol TE value, the 'needed TE' will be updated here.
    lStatus = CheckAndAdaptTE(rSeqLim, rMrProt, rSeqExpo, lNeededTE);

    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
        return SeverePrepareErrorReturn(lStatus);


    // ---------------------------------------------------------------------------
    // handle setting of thickness of IR-pulse
    // ---------------------------------------------------------------------------
    // NOTE: mySeqLoop.setScaleIRThickness(2.0)) was already passed to mySeqLoop in
    // fSEQInit.
#if (defined COMPILE_EP2D_DIFF) || (defined COMPILE_EP2D_SE) || (defined EP2D_MS)

    m_mySeqLoop.setScaleIRThickness(2.0);
    m_mySeqLoop.setInterleavedIRAllowed(rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL ? false : true);

    // Reduced inversion thickness for short inversion times (e.g. STIR) - avoid crosstalk
    // for nested inversion schemes.
    if (rMrProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE && rMrProt.ti()[0] < 500000)
    {
        setReducedIRThicknessForShortTI();
    }

#ifdef EP2D_MS
#if defined COMPILE_EP2D_SE

    {
        const double DEFAULT_IR_PULSE_THICKNESS_FACTOR_PERC = 77.0;
        const double DEFAULT_IR_PULSE_THICKNESS_PERC_THICK = 200.0;
        const double DEFAULT_IR_PULSE_THICKNESS_PERC_THIN = 125.0;

        double dRelIRPulseThicknessFactorPercent = DEFAULT_IR_PULSE_THICKNESS_PERC_THICK;

        if (rMrProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
        {

            // The following internal logic (adapted from a_tse.cpp) appears reasonable for all IR protocols:
            if (rMrProt.ti()[0] < 500000) // [us]
            {
                dRelIRPulseThicknessFactorPercent = DEFAULT_IR_PULSE_THICKNESS_PERC_THIN;
            }
            else
            {
                // ------------------------------------------------------------------------------------------
                // For FLAIR contrast we want to suppress inflow. Therefore we choose the IR pulse thickness
                // as thick as possible
                // ------------------------------------------------------------------------------------------
                //  search slice group with smallest distance factor
                MrProtSliceGroup sGroup = rMrProt.sliceGroupList()[0];
                double dDistFact_min = sGroup.distFactor();

                for (long lCntr = 1; lCntr < rMrProt.sliceGroupList().size(); ++lCntr)
                {
                    sGroup = rMrProt.sliceGroupList()[static_cast<int32_t>(lCntr)];
                    if (sGroup.distFactor() < dDistFact_min)
                    {
                        dDistFact_min = sGroup.distFactor();
                    }
                }
                const double dIRThickFact_max = SMSProperties::getNReducedSlices(rMrProt) == 1
                                                    ? 10.
                                                    : (1 + dDistFact_min) * rMrProt.concatenations();

                dRelIRPulseThicknessFactorPercent = dIRThickFact_max * DEFAULT_IR_PULSE_THICKNESS_FACTOR_PERC;
            }


            m_mySeqLoop.setScaleIRThickness(0.01 * dRelIRPulseThicknessFactorPercent);
        }
    }
#elif defined COMPILE_EP2D_FID
    if (rMrProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
    {
        m_mySeqLoop.setScaleIRThickness(1.0);
    }
#endif
#endif



    //  Align sign of the slice selection gradients of IR pulse and excitation pulse
    if (SeqBuildBlockIRsel* pIRsel = const_cast<SeqBuildBlockIRsel*>(&m_mySeqLoop.getSBBIRsel()))
    {
        bool bReverseGSAmplitude = false;

#ifdef ZOOM_2DRF
        if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
        {
            SBB2DPtx* pSBBExcite = dynamic_cast<SBB2DPtx*>(m_EPIKernel.getExcitationPointer());
            if (pSBBExcite != nullptr)
            {
                bReverseGSAmplitude = false;
            }
            else
            {
                SEQ_TRACE_ERROR.print("ERROR - nullptr pointer found");
                lStatus = MRI_SEQ_SEQU_ERROR;
                return SeverePrepareErrorReturn(lStatus);
            };
        }
        else

#endif
        {
            sRF_PULSE* pRFExcit;

            pRFExcit = dynamic_cast<SeqBuildBlockExcitationRFPulse*>(m_EPIKernel.getExcitationPointer())->getExcitationRFPulse();

            if (pRFExcit != nullptr)
            {
                bReverseGSAmplitude = (pRFExcit->getRequiredGSPolarity() < 0);
            }
            else
            {
                SEQ_TRACE_ERROR.print("ERROR - nullptr pointer found");
                lStatus = MRI_SEQ_SEQU_ERROR;
                return SeverePrepareErrorReturn(lStatus);
            }
        }
        pIRsel->setReverseGSAmplitude(bReverseGSAmplitude);
    }

#ifdef EP2D_MS
#ifdef SUPPORT_IIR
    // Configure continuous saturation
    m_mySeqLoop.ActivateContinuousFatSatPulsing(false);
    m_mySeqLoop.ActivateContinousMTSatPulsing(false);
    m_mySeqLoop.ActivateContinousRSatPulsing(false);

    if (rMrProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
    {
        if (rMrProt.preparationPulses().getlFatWaterContrast() == MrProtocolData::FatWaterContrast_FatSaturation)
        {
            m_mySeqLoop.ActivateContinuousFatSatPulsing(true);
        }

        if (rMrProt.preparationPulses().getucMTC())
        {
            m_mySeqLoop.ActivateContinousMTSatPulsing(true);
        }

        if (rMrProt.getsRSatArray().getlSize() > 0)
        {
            m_mySeqLoop.ActivateContinousRSatPulsing(true);
        }
    }
#endif
    // This avoids, that m_mySeqLoop switches the interleaving scheme across concats
    m_mySeqLoop.setInterleavedIR_SameOuterLoopsForAllConcats(true);   
#endif

    // ---------------------------------------------------------------------------
    // handle SPAIR pulse settings
    // ---------------------------------------------------------------------------
    if (protFacade.isSPAIRFatSat())
    {
        if (SysProperties::isLowField()) // Optimized SPAIR for low-field
        {
            if (!m_mySeqLoop.seteSPAIRPulse(SeqBuildBlockOptfs::LOW_FIELD_OPTIMIZED))
            {
                SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "ERROR: seteSPAIRPulse to LOW_FIELD_OPTIMIZED failed.");
                lStatus = MRI_SEQ_SEQU_ERROR;
                return SeverePrepareErrorReturn(lStatus);
            }
        }
        else
        {
            if (rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_ABDOMEN)
            {
                if (!m_mySeqLoop.seteSPAIRPulse(SeqBuildBlockOptfs::ABDOMEN_OPTIMIZED))
                {
                    SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "ERROR: seteSPAIRPulse to ABDOMEN_OPTIMIZED failed.");
                    lStatus = MRI_SEQ_SEQU_ERROR;
                    return SeverePrepareErrorReturn(lStatus);
                }
            }
            else if (rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_THORAX)
            {
                if (!m_mySeqLoop.seteSPAIRPulse(SeqBuildBlockOptfs::THORAX_OPTIMIZED))
                {
                    SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "ERROR: seteSPAIRPulse to THORAX_OPTIMIZED failed.");
                    lStatus = MRI_SEQ_SEQU_ERROR;
                    return SeverePrepareErrorReturn(lStatus);
                }
            }
            else if (rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_BREAST)
            {
                if (!m_mySeqLoop.seteSPAIRPulse(SeqBuildBlockOptfs::BREAST_OPTIMIZED))
                {
                    SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "ERROR: seteSPAIRPulse to BREAST_OPTIMIZED failed.");
                    lStatus = MRI_SEQ_SEQU_ERROR;
                    return SeverePrepareErrorReturn(lStatus);
                }
            }
            else // reset settings
            {
                if (!m_mySeqLoop.seteSPAIRPulse(SeqBuildBlockOptfs::STANDARD))
                {
                    SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "ERROR: seteSPAIRPulse to STANDARD failed.");
                    lStatus = MRI_SEQ_SEQU_ERROR;
                    return SeverePrepareErrorReturn(lStatus);
                }
            }
        }
    }

    // ---------------------------------------------------------------------------
    // configure fat suppression
    // ---------------------------------------------------------------------------
    // Switch between different FatSat modes in case FatSat is active in protocol
    if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_FatSaturation)
    {
        // Defect 658379 - The default fatsat should be improved for low-field
        if (SysProperties::isLowField())
        {
            m_mySeqLoop.seteFatSatPulse(SeqBuildBlockCSat::Sinc);
            m_mySeqLoop.setFatSatOffcenterFrequency(0); // offset of the sinc pulse  for 0.55T will be set in the
                                                        // SBBCSat, since it is consistant for all sequence.
        }
        else
        {
            if (rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_BRAIN)
            {
                m_mySeqLoop.seteFatSatPulse(SeqBuildBlockCSat::Sinc);
                m_mySeqLoop.setFatSatOffcenterFrequency(-50);
            }
            else
            {
                m_mySeqLoop.seteFatSatPulse(SeqBuildBlockCSat::Gauss);
                m_mySeqLoop.setFatSatOffcenterFrequency(0);
            }
        }

    }

    // ---------------------------------------------------------------------------
    // configure magnetization preparation
    // ---------------------------------------------------------------------------
    if (rMrProt.getsPrepPulses().getucMTC())
    {
        long   lMTReplications    = 0;
        switch (rMrProt.getsPrepPulses().getlMTCMode())
        {
            case MrProtocolData::MTC_MODE_STRONG:
                lMTReplications = 5; // => 6 MT-prep RF-pulses with alternating frequency offset
                break;
            case MrProtocolData::MTC_MODE_MEDIUM:
                lMTReplications = 3; // => 4 MT-prep RF-pulses with alternating frequency offset (to reduce SAR)
                break;
            case MrProtocolData::MTC_MODE_WEAK:
            case MrProtocolData::MTC_MODE_OFF: // Default if not visible
            default:
                break;
        }

        m_mySeqLoop.getpMSat()->setAdditionalReplications(lMTReplications);
    }

#endif // COMPILE_EP2D_DIFF || COMPILE_EP2D_SE || EP2D_MS

    // ---------------------------------------------------------------------------
    // set loop parameters different from standard settings
    // ---------------------------------------------------------------------------
    setVariantSpecificLoopSettings(rSeqLim, rMrProt);

           
    if (m_REOInfo.isPATActive() && (m_REOInfo.getPATAccelerationFactorPE() > 1))
    {
        m_mySeqLoop.setLinesToMeasure(m_REOInfo.getPATLinesPerSegment());
    }
    else
    {
        m_mySeqLoop.setLinesToMeasure(m_REOInfo.getLinesPerSegment());
    }


    //----------------------------------------------------------------------------
    // determine number of prep-scans required
    //----------------------------------------------------------------------------
    //
    // compare SeqLoop.cpp:
    //
    // - PreparingScans = (m_lPreparingTime + (TR*Phases) - 1) / (TR*Phases)
    // - total prep scans per concatenations = PreparingScans * lSlicesInConc * Phases
    //
    //----------------------------------------------------------------------------


    // basic initial dummy scans (without PAT, SMS, ...)
    setInitialDummyScansBasic(rMrProt);

    setPhaseCorrScansAndAdjustInitialDummyScans(rMrProt);



#ifdef ASL

    m_bEnableFirstPrepScanAsM0Scan = true;

#ifdef SUPPORT_iPAT_a_ep2d
    // with iPAT a separate M0 scan BEFORE the iPAT scan is not possible
    if (m_REOInfo.isPATActive() && m_REOInfo.getPATAccelerationFactorPE() > 1)
    {
        m_bEnableFirstPrepScanAsM0Scan = false;
    }
#endif // PAT

    // 	if (protFacade.isSliceAcceleration())
    // 	{
    // 		// with SliceAccel. a separate M0 scan BEFORE the reference scan is not possible
    // 		m_bEnableFirstPrepScanAsM0Scan = false;
    // 	}

    if (rMrProt.getlRepetitions() > 0)
    {
        // minimum PrepScanTime [ms] is taken from ASL UI
        m_lInitialDummyScans = (long)(ceil(m_ASL_SBB.getPrepScanTime() * 1e3 / rMrProt.tr()[0]));
    }
    else
    {
        // no preparing for 'real' single shot
        m_lInitialDummyScans = 0;

        // no separate M0 scan for single shot
        m_bEnableFirstPrepScanAsM0Scan = false;
    }

    if (m_bEnableFirstPrepScanAsM0Scan)
    {
        // ensure at least a single M0 scan: >= 1;
        m_lInitialDummyScans = std::max(1L, m_lInitialDummyScans);

        // decrease repetitions by one (first repetition is M0 reference scan measured as PrepScan)
        m_mySeqLoop.setRepetitionsToMeasure(rMrProt.getlRepetitions() - 1);
    }
    else // No M0 scan
    {
        m_mySeqLoop.setRepetitionsToMeasure(rMrProt.getlRepetitions());
    }

    m_mySeqLoop.setRepetitionValueForMdh(0); // After calling setRepetitionsToMeasure(), the repetition Mdh should be set to 0 for gre ref scan

#endif // ASL

#ifdef SUPPORT_iPAT_a_ep2d

    m_lPATRefScans = 0;   
    if (m_REOInfo.isPATActive() && m_REOInfo.getPATAccelerationFactorPE() > 1 && !isGreRefScanType(rMrProt))
    {
        // Make sure that a minimum number of prep scans is applied before the
        // PAT reference scans (evolution of steady state)
        if (m_lMinPrepScansNoPATRefScans > m_lInitialDummyScans + m_lPhaseCorrPrepScans)
        {
            m_lInitialDummyScans = m_lMinPrepScansNoPATRefScans - m_lPhaseCorrPrepScans;
        }

        if (m_bSegmentedRefLines)
        {
            m_lPATRefScans = m_REOInfo.getPATAccelerationFactorPE();
        }
        else
        {
            m_lPATRefScans = 1;
        }
    }
    
#endif

    // Overview of prep scans for slice accelerated case
    //
    // m_lInitialDummyScans---m_lPhaseCorrPrepScans---m_lSliceAccelRefScans---m_lPATPrepScans------m_lSliceAccelDummyScans---m_lSliceAccelPhaseCorrScans---m_lAdjPrepScans------Imaging----------
    //
    // ----------n---------------------1/0---------------------1-------------------1/0-----------------------m-----------------------1/0-------------------------1/0-----------------------------
    //
    // --------dummy-------------------DAQ---------------------DAQ-----------------DAQ----------------------dummy--------------------DAQ------------------------DAQ---------------DAQ------------
    //
    // --------all slices------------all slices---------------all slices----------all slices-----------reduced slices------------reduced slices--------------reduced slices----reduced slices----


    if (!isNumberOfCoilsSufficient(rMrProt, rSeqLim))
    {
        return MRI_SEQ_SEQU_SMS_PAT_NOT_ENOUGH_COILS;
    }

    // For slice acceleration two dummy scan types are added:
    // Type 1 acquired in singleband mode before the ACS scans
    // Type 2 acquired in multiband mode before the imaging scans
    
    m_lSliceAccelPhaseCorrScans = 0;
    if (protFacade.isSliceAcceleration())
    {
        setInitialDummyScansSMS(rMrProt);

        // If only FLASH reference scans get acquired: use a slice-accelerated scan for
        // obtaining (collapsed) phase-correction information
        if ((rMrProt.preparationPulses().getlPhaseCorrectionMode() == MrProtocolData::PHASECORR_EXTERNAL)
            && (isFastGreRefScan(rMrProt)))
        {
            m_lSliceAccelDummyScans     = 0;
            m_lPhaseCorrPrepScans       = 0;
            m_lSliceAccelPhaseCorrScans = 1;
        }

#if (defined BOLD && !defined PACE3D) || defined COMPILE_EP2D_DIFF /*|| (defined ASL && !defined PACE3D)*/

        // Check for custom dummy scan parameters in Ini file
        long lMBDummyDebug = static_cast<long>(m_debugSettings.getDefaultSetting<int32_t>("EPI_GENERAL/number_of_multiband_dummy_scans", -1));
        long lSBDummyDebug = static_cast<long>(m_debugSettings.getDefaultSetting<int32_t>("EPI_GENERAL/number_of_singleband_dummy_scans", -1));

        if (lMBDummyDebug != -1)
        {
            SEQ_TRACE_ALWAYS.print("Config File: Overwriting m_lSliceAccelDummyScans (%ld) with value: %ld.", m_lSliceAccelDummyScans, lMBDummyDebug);
            m_lSliceAccelDummyScans = lMBDummyDebug;
        }
        if (lSBDummyDebug != -1)
        {
            SEQ_TRACE_ALWAYS.print("Config File: Overwriting m_lInitialDummyScans (%ld) with value: %ld.", m_lInitialDummyScans, lSBDummyDebug);
            m_lInitialDummyScans = lSBDummyDebug;
        }
#endif
    }

#if (defined BOLD && !defined PACE3D) || defined COMPILE_EP2D_DIFF /*|| (defined ASL && !defined PACE3D)*/
    if (m_debugSettings.getDefaultSetting<bool>("EPI_GENERAL/dump_dummy_scan_timing", false) && !rSeqLim.isContextPrepForBinarySearch())
    {
        m_mySeqLoop.setRunMode(MULTI_BAND);
        long ltCoolPauseExplicit = m_lCoolPauseTotal - m_lCoolPauseImplicit;
        long ltScanTimeSatsEtc   = m_mySeqLoop.getlScanTimeAllSats();   // SBB scan time
        long ltScanTimeBasic     = m_EPIKernel.getDurationPerRequest(); // Kernel scan time excluding mandatory fill time
        m_mySeqLoop.setRunMode(SINGLE_BAND);
        long ltTR_SliceAccPrepMeas
            = (ltScanTimeBasic + m_mySeqLoop.getlTRFillInConcat(0) + ltScanTimeSatsEtc + ltCoolPauseExplicit) * (rMrProt.sliceSeries().getlSize()) + m_mySeqLoop.getlTRFillEndInConcat(0);

        SEQ_TRACE_ALWAYS.print("TR for preparation: %ld", ltTR_SliceAccPrepMeas);
        SEQ_TRACE_ALWAYS.print("m_lInitialDummyScans: %ld", m_lInitialDummyScans);
        SEQ_TRACE_ALWAYS.print("Total: %ld", ltTR_SliceAccPrepMeas * m_lInitialDummyScans);
        SEQ_TRACE_ALWAYS.print("TR for imaging: %d", rMrProt.tr()[0]);
        SEQ_TRACE_ALWAYS.print("m_lSliceAccelDummyScans: %ld", m_lSliceAccelDummyScans);
        SEQ_TRACE_ALWAYS.print("Total: %ld", rMrProt.tr()[0] * m_lSliceAccelDummyScans);
    }
#endif

    // ---------------------------------------------------------------------------
    // setting slice-accelerated prep scans
    // ---------------------------------------------------------------------------
    if (protFacade.isSliceAcceleration() && !isFastGreRefScan(rMrProt))
    {
        m_lSliceAccelRefScans = 1;
    }
    else
    {
        m_lSliceAccelRefScans = 0; // Fast GRE ref scan provides data for slice GRAPPA
    }

    // ---------------------------------------------------------------------------
    // calculate total number of required prep scans, including dummy scans
    // ---------------------------------------------------------------------------
    lRequiredPrepScans = calcRequiredPrepScans(rMrProt);

    // ---------------------------------------------------------------------------
    // Store the number of prep scans without ADC separately
    // ---------------------------------------------------------------------------

    lRequiredPrepScansWithoutADC = m_lInitialDummyScans;

    if (protFacade.isSliceAcceleration())
    {
        lRequiredPrepScansWithoutADC += m_lSliceAccelDummyScans;
    }

    // ---------------------------------------------------------------------------
    // force SeqLoop to use a specific time for prep scans
    // ---------------------------------------------------------------------------
    m_mySeqLoop.setlPreparingTime(lRequiredPrepScans * rMrProt.tr()[0] * rMrProt.physiology().phases());
    m_mySeqLoop.setPreparingWithoutADCTime(lRequiredPrepScansWithoutADC * rMrProt.tr()[0] * rMrProt.physiology().phases());

    // ---------------------------------------------------------------------------
    // Advise SeqLoop to perform prep-scans only within the first measurement.
    // ---------------------------------------------------------------------------
    if (m_bPrepScansOnlyInFirstMeasurement)
    {
        m_mySeqLoop.setePerformPreparingScans(OnlyFirstRepetition);
    }
    else
    {
        m_mySeqLoop.setePerformPreparingScans(Always);
    }

    // ---------------------------------------------------------------------------
    // Set long TR triggering mode
    // ---------------------------------------------------------------------------

    setLongTRTrigMode(rMrProt);


    // -------------------------------------------------------------------------
    // Apply restrictions for multiple concatenations
    // -------------------------------------------------------------------------
    if (rMrProt.concatenations() > 1)
    {
        // Return with error if multiple concatenations are not allowed
        if (!isMultiConcatsAllowed(rMrProt))
        {
            lStatus = MRI_SEQ_SEQU_ERROR;

            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "Multiple concatenations are not possible ");

            return SeverePrepareErrorReturn(lStatus);
        }

        // check dynamic field correction with long TR trigger mode
        if (!checkDFCWithLongTRTrigMode(rSeqLim, rMrProt))
        {
            lStatus = MRI_SEQ_SEQU_ERROR;
            return SeverePrepareErrorReturn(lStatus);
        }
    }

    // navigator triggering

    if ((rMrProt.NavigatorParam().getlRespComp() == SEQ::RESP_COMP_TRIGGER) || (rMrProt.NavigatorParam().getlRespComp() == SEQ::RESP_COMP_TRIGGER_AND_FOLLOW))
    {
        m_mySeqLoop.setOptfsPrepflag(true);
    }
    else
    {
        m_mySeqLoop.setOptfsPrepflag(false);
    }

    m_mySeqLoop.setpCppSequence(this);

#if (defined SUPPORT_iPAT_a_ep2d) || (defined SUPPORT_iPAT_TGSE)
    {
        bool bAcquisitionRequiresFLASHRefScan = isGreRefScanType(rMrProt)
                                            && (rMrProt.PAT().getucPATMode() != SEQ::PAT_MODE_NONE)
                                            && (rMrProt.PAT().getlAccelFactPE() > 1);
        bool bUserRequestsFLASHRefScan = isFastGreRefScan(rMrProt);

        long lPATFlashRefScanBaseRes = m_lPATFlashRefScanBaseRes;
        long lPATFlashRefScanBandwidth = m_lPATFlashRefScanBandwidth;

        if (bUserRequestsFLASHRefScan)
        {
            lPATFlashRefScanBaseRes   = std::min(getFastGreRefScanBaseRes(), static_cast<uint16_t>(rMrProt.getsKSpace().getlBaseResolution()));
            lPATFlashRefScanBandwidth = getFastGreRefScanBandwidth();
        }

        // Configure separate Flash reference scan
        if (bAcquisitionRequiresFLASHRefScan || bUserRequestsFLASHRefScan)
        {
            
            m_mySeqLoop.setePATRefScanLoopMode(SeqLoop::PATRefScanLoopMode_ONCE);
            m_mySeqLoop.setlSBBPATRefScanBaseRes(lPATFlashRefScanBaseRes);
            m_mySeqLoop.setlSBBPATRefScanBandwidth(lPATFlashRefScanBandwidth);
            m_mySeqLoop.setSkipOnlinePhaseCorrFlag();
            m_mySeqLoop.setSkipRegriddingFlag();
        }
        else
        {
            m_mySeqLoop.setePATRefScanLoopMode(SeqLoop::PATRefScanLoopMode_NEVER);
        }
    }

#endif

    // Set required dynamic adjustments
    m_mySeqLoop.setsSliceAdjParametersRequestedBySequence(m_sSliceAdjParametersRequestedBySequence);
    // In some cases, TI counting does not end at the beginning of the corresponding
    // excitation SBB but after a certain offset only. This information is required for
    // the fill time calculations.
    //
    // 1. The "TimeToCriticalEvent" mechanism gets explicitly employed. If the
    //    corresponding property is set, the SBB provides information about the
    //    actual start time. Note that this time will include the transit-RTEB
    //    duration if dynamic adjustments are enabled.
    // 2. With dynamic adjustments, the transit-RTEB duration has to be considered.
    if (m_EPIKernel.getbUseTimeToCriticalEvent_us())
    {
        // Note that this contains contributions from dynamic adjustments already
        m_mySeqLoop.setlExcOffsetForTICalc_us(m_EPIKernel.getlTimeToCriticalEvent_us());
    }
    else
    {
        m_mySeqLoop.setlExcOffsetForTICalc_us(m_EPIKernel.getSliceAdjUpdateDuration());
    }

    // ---------------------------------------------------------------------------
    // prepare standard loop
    // ---------------------------------------------------------------------------

#ifdef SUPPORT_iPAT_a_ep2d
    m_mySeqLoop.setNumberOfPATPrepScans(m_lPATRefScans);
    m_mySeqLoop.setNumberOfInitialDummyScans(m_lInitialDummyScans);
    m_mySeqLoop.setNumberOfPhaseCorrectionScans(m_lPhaseCorrPrepScans);
    if (protFacade.isSliceAcceleration())
    {
        m_mySeqLoop.setNumberOfSliceAccelDummyScans(m_lSliceAccelDummyScans);
        m_mySeqLoop.setNumberOfSliceAccelPrepScans(m_lSliceAccelRefScans);
        m_mySeqLoop.setNumberOfSliceAccelPhaseCorrectionScans( m_lSliceAccelPhaseCorrScans );
    }
#endif

#ifdef COMPILE_EP2D_DIFF
    if (m_EPIKernel.getbCompensationEnable() && m_mySeqLoop.getpCSatFat() != nullptr)
    {
        if (m_EPIKernel.getPointerCompGrad()->isPrepared())
        {
            if (rMrProt.satList().size() < 1)
            {
                m_mySeqLoop.getpCSatFat()->setSpoilerFront(10, 10, 10, 0);
            }

            m_mySeqLoop.getpCSatFat()->setSpoilerBack(10, 10, 10, 0);
        }
        else
        {
            m_mySeqLoop.getpCSatFat()->resetSpoiler(true);
        }
    }
    else if (m_mySeqLoop.getpCSatFat() != nullptr)
    {
        m_mySeqLoop.getpCSatFat()->resetSpoiler(true);
    }
#endif

    setSPAIRSpoilingType();

    if (!m_mySeqLoop.prep(rMrProt, rSeqLim, rSeqExpo))
        mSBBErrGotoFinish(m_mySeqLoop, "m_mySeqLoop.prep");

    // Do a separate preparation of the Pat ref scan
    if ( ( m_mySeqLoop.getePATRefScanLoopMode() != SeqLoop::PATRefScanLoopMode_NEVER ) && SMSProperties::isSMS( rMrProt ) )
    {
        
        if (!isFastGreRefScan(rMrProt))
        {
            // ... prepare the reference scan as if SMS was disabled
            // Note: This will just remove the corresponding flag from the reference scan mdh.
            long lPreviousSMSFactor = SMSProperties::getMultiBandFactor(rMrProt);
            rMrProt.getsSliceAcceleration().setlMultiBandFactor(1);
            m_mySeqLoop.getPATRefScan()->prepSBB(rMrProt, rSeqLim, rSeqExpo);
            // Reset SMS factor
            rMrProt.getsSliceAcceleration().setlMultiBandFactor(lPreviousSMSFactor);
        }
    }

#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
    if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_FatSaturation)
    {
        MrProtocolData::SeqExpoRFInfo rfInfoInSBBs;

        if (rMrProt.getlRepetitions() > 0)
        {
            // FatSat / Spoil gradient preparation - FATSAT_ALL_SLICES
            m_CSatFat.setRequestsPerMeasurement(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas));
        }
        else
        {
            // calc request for FirstMeas without PrepScans
            m_CSatFat.setRequestsPerMeasurement(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, FirstMeas));
        }

        m_CSatFat.setCSatMode(SBBCSatCode_Fat);
        m_CSatFat.setIdent("AddCSat");
        m_CSatFat.setGSWDGradientPerformance(rMrProt, rSeqLim);

        m_SpoilGrad.setRequestsPerMeasurement(m_CSatFat.getRequestsPerMeasurement());
        m_SpoilGrad.setGSWDGradientPerformance(rMrProt, rSeqLim);

        if (!m_SBB.prepSBBAll(rMrProt, rSeqLim, rSeqExpo, &rfInfoInSBBs, m_sSliceAdjParametersRequestedBySequence))
        {
            return (m_SBB.getpSBBLastPrep()->getNLSStatus());
        }
    }

#endif // #if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE

#if defined ZOOM_2DRF
    if ((rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED) && (eExcType == MrProtocolData::PTXTrajectoryType_EPI_1D))
    {
        std::vector<double> vdThickness;
        std::vector<double> vdShift;

        // access to SBBExcitation only AFTER m_EPIKernel.prep !
        SBB2DPtx* pSBBExcite = dynamic_cast<SBB2DPtx*>(m_EPIKernel.getExcitationPointer());
        if (pSBBExcite != nullptr)
        {
            pSBBExcite->getOptPTXVolume(vdThickness, vdShift);
        }

        MrProtocolData::SeqExpoRFInfo rfInfoInSBBs;

        if (rMrProt.getlRepetitions() > 0)
        {
            m_OptPTXVolume.setRequestsPerMeasurement(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas));
        }
        else
        {
            // calc request for FirstMeas without PrepScans
            m_OptPTXVolume.setRequestsPerMeasurement(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, FirstMeas));
        }

        m_OptPTXVolume.setIdent("OptPTXVolume");
        m_OptPTXVolume.setGSWDGradientPerformance(rMrProt, rSeqLim);

        m_OptPTXVolume.init(vdThickness, vdShift);

        if (!m_SBBPTX.prepSBBAll(rMrProt, rSeqLim, rSeqExpo, &rfInfoInSBBs))
        {
            return (m_SBBPTX.getpSBBLastPrep()->getNLSStatus());
        }
    }
#endif

#ifndef ASL
#ifdef PACE3D
    m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(TWAKEUP);
#elif defined EPI_SUPPORT_FREQ_FEEDBACK
    if (m_bB0Correction)
    {
        // no wait-for-wakeup when in real-time mode
        m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(m_sFreqFeedback.getWakeUpDuration(m_bIsRealtimeProcessingEnabled ? false : true));
    }
    else
    {
        m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(0);
    }
#else
    m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(0);
#endif // #ifdef PACE3D
#else

    if (rMrProt.preparationPulses().getucFatSatMode() == SEQ::FAT_SAT_STRONG)
    {
        MrProtocolData::SeqExpoRFInfo rfInfoInSBBs;

        if (rMrProt.getlRepetitions() > 0)
        {
            // FatSat / Spoil gradient preparation - FATSAT_ALL_SLICES
            m_CSatFat.setRequestsPerMeasurement(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas));
        }
        else
        {
            // calc request for FirstMeas without PrepScans
            m_CSatFat.setRequestsPerMeasurement(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, FirstMeas) / (lRequiredPrepScans + 1));
        }

        m_CSatFat.setCSatMode(SBBCSatCode_Fat);
        m_CSatFat.setIdent("PaslFS2");
        m_CSatFat.setGSWDGradientPerformance(rMrProt, rSeqLim);

        m_SpoilGrad.setRequestsPerMeasurement(m_CSatFat.getRequestsPerMeasurement());
        m_SpoilGrad.setGSWDGradientPerformance(rMrProt, rSeqLim);

        if (!m_SBB.prepSBBAll(rMrProt, rSeqLim, rSeqExpo, &rfInfoInSBBs, m_sSliceAdjParametersRequestedBySequence))
        {
            return (m_SBB.getpSBBLastPrep()->getNLSStatus());
        }
    }

    // set Requests for FatSat in ASL
    m_ASL_SBB.setRequestsPerMeasurement(1);

    // set time between start of SBB and excitation of the first image slice
    m_ASL_SBB.setTimeToImageSliceExc_us(m_EPIKernel.getExcitationPointer()->getlRFCenterTime()); // ASL does not support water excitation ==> no support for m_SBBExcitationBinomial

    m_ASL_SBB.setASLMode(rMrProt.Asl().getulMode());

    // Provide requested dynamic adjustments
    m_ASL_SBB.setsSliceAdjParametersRequestedBySequence(m_sSliceAdjParametersRequestedBySequence);

    if (!m_ASL_SBB.prep(rMrProt, rSeqLim, rSeqExpo))
    {
        return (m_ASL_SBB.getNLSStatus());
    }

    if (rMrProt.preparationPulses().getucFatSatMode() == SEQ::FAT_SAT_STRONG)
    {
#ifdef PACE3D
        // FATSAT_ALL_SLICES
        // tell seqloop that first FatSat is not to be taken into account
        m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(TWAKEUP + m_ASL_SBB.getDurationPerRequest() - m_CSatFat.getDurationPerRequest() - m_SpoilGrad.getDurationPerRequest());
#endif // PACE3D
    }
    else
    {
#ifdef PACE3D
        m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(TWAKEUP + m_ASL_SBB.getDurationPerRequest());
#else
        m_mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(m_ASL_SBB.getDurationPerRequest());
#endif // PACE3D
    }

#endif // ASL

    //----------------------------------------------------------------------------
    // Check, if we got the correct number of preparing scans from SeqLoop.
    // This is important for PAT-measurements, because reference lines are acquired
    // during the last prep-scans.
    //----------------------------------------------------------------------------
#ifdef SUPPORT_iPAT_a_ep2d
    {
        if (m_mySeqLoop.getlPreparingScans() != lRequiredPrepScans)
        {
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "NEED %ld PrepScans from SeqLoop, BUT got %ld", lRequiredPrepScans, m_mySeqLoop.getlPreparingScans());
            lStatus = MRI_SEQ_SEQU_ERROR;
            return SeverePrepareErrorReturn(lStatus);
        }
    }
#endif

    // ---------------------------------------------------------------------------
    // calculate mandatory fill times
    // ---------------------------------------------------------------------------
    // For some Magnetoms, an elevated gradient performance is used during the
    // EPI readout that requires the introduction of an additional pause afterwards
    // in order to stay within hardware limitations. This is realized by introducing
    // a mandatory fill time after each kernel. This time is calculated based on the
    // duration of the echo train.
    //
    // Exception #1: for the EP2D_SE variant no additional fill time is neccessary
    // (no gradient activity between excitation and refocussing)
    //
    // Exception #2: for the EP2D_DIFF variant, an additional dynamic component is
    // required which is calculated based on the gradient load of the diffusion
    // module. Note: In this case, the same pause duration has to be used in
    // fPrepDIFFPlugInForEPIKernel!

    // m_lCoolPauseTotal    is required (and possibly adapted) by calculateTRTIFillTimes
    // m_lCoolPauseImplicit is calculated by calculateTRTIFillTimes

    m_lCoolPauseTotal    = 0;
    m_lCoolPauseImplicit = 0;

    if (m_EPIKernel.getUseGPABalance())
    {
        // Receive required cool pause per execution of EPI kernel (including plugin)
        // Note: EPIKernel has to be prepared
        m_lCoolPauseTotal = fSDSRoundDownGRT(m_EPIKernel.getTRIncrement());
    }
    else
    {
        // All other variants: calculate cooling pause as a fraction of the readout train
        // Note: EPIKernel and REOInfo have to be prepared
        //
        // Strategy:    (t_readout * G_readout^2) / (t_readout + t_pause) = G_nominal^2
        //           =>  t_pause = t_readout * (G_readout^2 / G_nominal^2 - 1)
        //
        // - Quadratic increase of losses with gradient amplitude is assumed
        // - Ramps are not considered yet

        double dAmplReadout = m_EPIKernel.getGRO().getAmplitude();
        double dAmplNominal = SysProperties::getGradMaxAmplNominal();
        double dTRIncFactor = (dAmplReadout * dAmplReadout) / (dAmplNominal * dAmplNominal) - 1.;
        long   lEchoSpacing = m_EPIKernel.getEchoSpacing();
        long   lEchoNumber  = m_REOInfo.getEchoTrainLength();

        if (dTRIncFactor < 0.)
        {
            dTRIncFactor = 0.;
        }

        m_lCoolPauseTotal = fSDSRoundDownGRT(dTRIncFactor * static_cast<double>(lEchoNumber * lEchoSpacing));
    }

    // ---------------------------------------------------------------------------
    // calculate the TR/TI fill times
    // ---------------------------------------------------------------------------
    if (!calculateTRTIFillTimes(rMrProt, rSeqLim, rSeqExpo, &lNeededTI, &lNeededTR))
    {
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            CheckStatusPB(MRI_SEQ_SEQU_ERROR, "calculateTRTIFillTimes failed");
        }
        else
        {
            CheckStatusB(MRI_SEQ_SEQU_ERROR);
        }
    }

#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
    // Decide whether an additional CSat is played out at the beginning of the
    // cool pause.
    m_bApplyExtraCSat = isExtraCSatApplied(rMrProt);


#endif // #if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE

    // ---------------------------------------------------------------------------
    // Prepare frequency stabilization
    // ---------------------------------------------------------------------------
#ifdef EPI_SUPPORT_FREQ_FEEDBACK
    if (m_bB0Correction)
    {
        // Initialize volume counter:
        // In contrast to the repetitions counter, all acquired scans including
        // iPAT reference and adjustment scans are considered.
        m_lVolumeCounter = 0;

        if (rSeqLim.isContextNormal())
        {
            long lTotalNoOfVolumes = calcTotalNumberOfVolumesForFreqFeedback(rMrProt);

            // Feedback preparation - synchronization and incorporation at the end of each volume acquisition
            if (!m_sFreqFeedback.Prep(rMrProt, rSeqLim, rSeqExpo, lTotalNoOfVolumes, false))
            {
                // Error can only be caused by programming error
                SEQ_TRACE_ERROR.print("m_sFreqFeedback.Prep() failed.");
                lStatus = MRI_SEQ_SEQU_ERROR;
                return SeverePrepareErrorReturn(lStatus);
            }
        }
    }
#endif

    // ---------------------------------------------------------------------------
    // prep pace stuff:
    // copy slice position data
    // ---------------------------------------------------------------------------
#ifdef PACE3D
    if (!m_PaceFeedback.Prep(rMrProt, m_asSLC))
    {
        // error can only be caused by programming error
        // => trace also if rSeqLim.isContextPrepForBinarySearch()
        SEQ_TRACE_ERROR.print("PaceFeedback.Prep failed.");
        lStatus = MRI_SEQ_SEQU_ERROR;
        return SeverePrepareErrorReturn(lStatus);
    }
#endif

    // ---------------------------------------------------------------------------
    // calculate energy
    // ---------------------------------------------------------------------------
    if (protFacade.isSliceAcceleration())
    {
        // ---------------------------------------------------------------------------
        // Slice-accelerated case
        // ---------------------------------------------------------------------------

        lTotalSlices   = SMSProperties::getNSlices(rMrProt);
        lReducedSlices = SMSProperties::getNReducedSlices(rMrProt);

        // Returns preparation pulses energy
        RFInfo = m_mySeqLoop.getRFInfo(rMrProt);

        lNMeasPrepScanTR = 0;
#ifdef SUPPORT_iPAT_a_ep2d
        lNMeasPrepScanTR += m_lPATRefScans;
#endif
        lNMeasPrepScanTR += m_lSliceAccelRefScans;
        lNMeasPrepScanTR += m_lInitialDummyScans;
        lNMeasPrepScanTR += m_lPhaseCorrPrepScans;

        // getKernelRequestsPerMeasurement multiplies with the reduced number of slices as it is passed to by the
        // prep method of the SeqLoopMultiband. The division accounts for that
        lNMeasTotal         = m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, FirstMeas) / lReducedSlices;
        lNMeasImagingScanTR = lNMeasTotal - lNMeasPrepScanTR;

        // Repetitions do not include reference and dummy scans anymore
        // same logic for the division here
        lNMeasImagingScanTR_SecondMeas = (m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas) / lReducedSlices) * rMrProt.getlRepetitions();

        // RFInfoPerSlice is average energy per unit time throughout the whole acquisition
        // Todo: Redistribute into different RFInfoBlocks for preparation and measurement scans with differing TRs
        RFInfoSeqLoopPerSlice = RFInfo / static_cast<double>((lNMeasTotal + lNMeasImagingScanTR_SecondMeas) * lReducedSlices);

        RFInfo.clear();

        // Single-band reference scans
        RFInfo += RFInfoSeqLoopPerSlice * (double)(lNMeasPrepScanTR * lTotalSlices);
        m_EPIKernel.setRunMode(SINGLE_BAND);
        RFInfo += m_EPIKernel.getRFInfoPerRequest() * (double)(lNMeasPrepScanTR * lTotalSlices);

        // Multi-band mode
        m_EPIKernel.setRunMode(MULTI_BAND);

        // Slice-accelerated prep scans and first imaging scan
        RFInfo += RFInfoSeqLoopPerSlice * (double)(lNMeasImagingScanTR * lReducedSlices);
        RFInfo += m_EPIKernel.getRFInfoPerRequest() * (double)(lNMeasImagingScanTR * lReducedSlices);

        // Repeated slice-accelerated imaging scans
        RFInfo += m_EPIKernel.getRFInfoPerRequest() * (double)(lNMeasImagingScanTR_SecondMeas * lReducedSlices);
        RFInfo += RFInfoSeqLoopPerSlice * (double)(lNMeasImagingScanTR_SecondMeas * lReducedSlices);

#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
        if (m_bApplyExtraCSat)
        {
            // Additional fat saturation pulse applied explicitly after each kernel execution
            RFInfo += m_CSatFat.getRFInfoPerRequest() * (double)(lNMeasPrepScanTR * lTotalSlices + lNMeasImagingScanTR * lReducedSlices + lNMeasImagingScanTR_SecondMeas * lReducedSlices);
        }
#endif // defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
    }
    else
    {
        // ---------------------------------------------------------------------------
        // Non-accelerated case
        // ---------------------------------------------------------------------------

        RFInfo = m_mySeqLoop.getRFInfo(rMrProt);

        RFInfo += m_EPIKernel.getRFInfoPerRequest()
                  * (double)(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, FirstMeas) + m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas) * rMrProt.getlRepetitions());

#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
        if (m_bApplyExtraCSat)
        {
            // Additional fat saturation pulse applied explicitly after each kernel execution
            RFInfo += static_cast<double>(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, FirstMeas)) * m_CSatFat.getRFInfoPerRequest()
                      + static_cast<double>(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas)) * m_CSatFat.getRFInfoPerRequest() * (double)rMrProt.getlRepetitions();
        }
#endif // defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
    }

    //-------------------------------------------------------------------------------------------------
    // Inform Unit Test exceptions for slice acceleration EPI sequences ref scan behavior
    //-------------------------------------------------------------------------------------------------
#ifdef WIN32
    if (protFacade.isSliceAcceleration())
    {
        if (SeqUT.isUnitTestActive())
        {
            //-------------------------------------------------------------------------------------------------
            // Slice accelerated Ref scans have different TR than that specified in the sequence
            //-------------------------------------------------------------------------------------------------
            long lExpectedNotOK = rMrProt.getsSliceArray().getlSize() * (m_lPATRefScans + m_lInitialDummyScans + m_lPhaseCorrPrepScans) + isFastGreRefScan(rMrProt) ? 0 : SMSProperties::getNReducedSlices(rMrProt);
            /* MZ: Todo: Currently disabled down below. Reactivate
             * if (lTRClockErr>0)
                SeqUT.SetExpectedNotOk(lTRClockErr, RTEB_ClockCheck, lExpectedNotOK, "Slice accelerated Ref scans have different TR than that specified in the sequence");
                */
            //-------------------------------------------------------------------------------------------------
            // Inform UT about the slice grappa Ref scan mismatch.
            // Root cause is how the UT calculate the line range of slice grappa reference scans: it is calculated based on the the line number of k-space center and number of in-plane grappa
            // reference scans. Please check fillAllMCEs in \src\MrImagingFW\ut\libSeqUTII\DataCollector\MomentChangeEventCollector.cpp The working around below cannot cover all cases. The problem
            // should be fixed from UT side in the future.
            //-------------------------------------------------------------------------------------------------
            long lExpectedNotOKPerSlice = 0;

            if (rMrProt.PAT().getlAccelFactPE() < 2)
            {
                long lSMSRefScans = m_REOInfo.getEchoTrainLength();
                long lCenterLinNo = m_REOInfo.getKSCenterLin();

                long lUTMinExtraRefLin = std::max(lCenterLinNo - lSMSRefScans / 2, 0L);
                long lUTMaxExtraRefLin = lCenterLinNo - lSMSRefScans / 2 + lSMSRefScans;

                lExpectedNotOKPerSlice = lSMSRefScans - rMrProt.PAT().getlRefLinesPE();
                lExpectedNotOKPerSlice -= (lSMSRefScans - lUTMaxExtraRefLin);
                lExpectedNotOKPerSlice -= lUTMinExtraRefLin;
            }
            else
            {
                if (rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI)
                {
                    lExpectedNotOKPerSlice = std::max((m_REOInfo.getEchoTrainLength() - (rMrProt.PAT().getlRefLinesPE() + 1)) / rMrProt.PAT().getlAccelFactPE(), 0L);
                }
                else if (isFastGreRefScan(rMrProt))
                {
                    lExpectedNotOKPerSlice = 0;
                }
                else
                {
                    long lUTMinExtraRefLin = m_REOInfo.getKSCenterLin() - rMrProt.PAT().getlRefLinesPE() / 2;
                    long lUTMaxExtraRefLin = lUTMinExtraRefLin + rMrProt.PAT().getlRefLinesPE();

                    lExpectedNotOKPerSlice = std::max((lUTMinExtraRefLin - 0), 0L) + std::abs(lUTMaxExtraRefLin - rMrProt.PAT().getlRefLinesPE());
                     // SeqUT counting is not prepared for this case. Separate ref lines for PAT ref scan are labelled 1..NPE-1 but SeqUT checks center region
                     // SMS EPI Ref scan lines are labeled same as imaging scan lines but SeqUT detects some central lines as duplicates of PAT ref lines
                     // and peripheral lines as invalid pat lines (every line in beginning of k-space but only every Nth right to center but stopping way before end)
                     // SeqUT shall check different ref scan types separately in the future
                }
            }

            lExpectedNotOK = lExpectedNotOKPerSlice * rMrProt.getsSliceArray().getlSize();
            if (lExpectedNotOK > 0)
            {
                SeqUT.SetExpectedNotOk(lPPAMissingExtraReferenceSample, RTEB_ORIGIN_fSEQRunFinish, lExpectedNotOK, "Slice accelerated Ref scans have different number of lines than GRAPPA Ref scans");
            }
        }
    }
#endif // End of WIN32

    // Add energy from iPAT reference scan  (SBBPATRefScan)
    if ( m_mySeqLoop.getePATRefScanLoopMode() == SeqLoop::PATRefScanLoopMode_ONCE )
    {
        RFInfo += m_mySeqLoop.getTotalRFInfoPATRefScan();
    }

#if defined ZOOM_2DRF
    if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
    {
        // Additional OptPTXVolume
        RFInfo += static_cast<double>(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, FirstMeas)) * m_OptPTXVolume.getRFInfoPerRequest()
                  + static_cast<double>(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas)) * m_OptPTXVolume.getRFInfoPerRequest() * (double)rMrProt.getlRepetitions();
    }
#endif

#ifdef ASL

    RFInfo += (static_cast<double>(m_ASL_SBB.getRequestsPerMeasurement()) * m_ASL_SBB.getRFInfoPerRequest()) * static_cast<double>(rMrProt.getlRepetitions() + 1 + lRequiredPrepScans);

    // Add energy of Strong FatSat
    if (rMrProt.preparationPulses().getucFatSatMode() == SEQ::FAT_SAT_STRONG)
    {
        if ((m_CSatFat.getRequestsPerMeasurement() > 1))
        {
            RFInfo += static_cast<double>(m_CSatFat.getRequestsPerMeasurement() - 1) * m_CSatFat.getRFInfoPerRequest() * static_cast<double>(rMrProt.getlRepetitions() + 1 + lRequiredPrepScans);
        }
    }

    if (m_bEnableFirstPrepScanAsM0Scan)
    {
        // Add energy of reference scan which is not accounted for
        RFInfo += m_ASL_SBB.getReferenceScanRFInfo();

        // one volume is explicitly taken out of SeqLoop (pMrProt->getlRepetitions() - 1) and acquired as PrepScan (M0 volume)
        m_EPIKernel.setRequestsPerMeasurement(m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas));
        RFInfo -= m_EPIKernel.getRFInfoPerRequest() * m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas);

        // subtract one volume in SeqLoop and one REFERENCE (RF switched off) from PrepScans
        RFInfo -= static_cast<double>(2 * m_ASL_SBB.getRequestsPerMeasurement()) * m_ASL_SBB.getRFInfoPerRequest();

        // subtract energy of Strong FatSat (taken out of SeqLoop)
        if (rMrProt.preparationPulses().getucFatSatMode() == SEQ::FAT_SAT_STRONG)
        {
            if (m_CSatFat.getRequestsPerMeasurement() > 1)
            {
                RFInfo -= static_cast<double>(m_CSatFat.getRequestsPerMeasurement() - 1) * m_CSatFat.getRFInfoPerRequest();
            }
        }
    }

#endif // ASL

    //---------------------------------------------------------------------------
    // activate online ice-process
    //---------------------------------------------------------------------------
    m_EPIKernel.getReadOutAddress()->getMDH().addToEvalInfoMask(MDH_ONLINE);

    // ---------------------------------------------------------------------------
    // set the gain of the receiver
    // ---------------------------------------------------------------------------
    lStatus = m_pSSL->setRxGain(K_RX_GAIN_CODE_HIGH, rMrProt, rSeqLim);

    // error can only be caused by programming error
    // => trace also if rSeqLim.isContextPrepForBinarySearch()
    CheckStatusPB(lStatus, "SSLProfileStatic::setRxGain");

    // ---------------------------------------------------------------------------
    // calculate effective echo-spacing
    // and bandwidth per pixel in phase-encode direction
    // ---------------------------------------------------------------------------
    if (!m_EPIKernel.calcEffEchoSpacingAndBWPerPixelPE(rMrProt, lEffectiveEchoSpacing, dBandwidthPerPixelPE))
    {
        SEQ_TRACE_ERROR.print("ERROR: calcEffEchoSpacingAndBWPerPixelPE() failed");
        lStatus = MRI_SEQ_SEQU_ERROR;
        return SeverePrepareErrorReturn(lStatus);
    }


    // ---------------------------------------------------------------------------
    // prepare common exports
    // ---------------------------------------------------------------------------
    int32_t iB0CorrectionMask;

    iB0CorrectionMask = 0;

    if (m_bB0Correction)
    {
        iB0CorrectionMask |= EPIPC_FREQ_CORRECTION;
        iB0CorrectionMask |= EPIPC_LOGGING;
#ifdef ASL
        iB0CorrectionMask |= EPIPC_ASL_MODE;
#endif //   #ifdef ASL
#ifdef EPI_SUPPORT_FREQ_FEEDBACK
        iB0CorrectionMask |= EPIPC_SEND_FEEDBACK;
        iB0CorrectionMask |= EPIPC_REALTIME;
#endif
    }

    // Slice acceleration functors use EPIPhaseCorrPE as their hook
    // => Enable at least logging to ensure availability of functor in the pipeline
    if (SMSProperties::isSMS(rMrProt) && (iB0CorrectionMask == 0))
    {
        iB0CorrectionMask |= EPIPC_LOGGING;
    }

    rSeqExpo.clearScanningSequence();
    rSeqExpo.addtoScanningSequence("EP");
    rSeqExpo.setPhaseCorScans(2);
    rSeqExpo.setTotalMeasureTimeUsec(m_mySeqLoop.getTotalMeasTimeUsec(rMrProt, rSeqLim));
    rSeqExpo.setRFInfo(RFInfo);
    rSeqExpo.setMeasuredPELines(static_cast<int32_t>(m_REOInfo.getLinesToMeasure()));
    rSeqExpo.setMeasured3dPartitions(static_cast<int32_t>(m_REOInfo.getPartitionsToMeasure()));
    rSeqExpo.setEchoSpacing(static_cast<int32_t>(m_EPIKernel.getEchoSpacing()));
    rSeqExpo.setEffectiveEpiEchoSpacing(static_cast<int32_t>(lEffectiveEchoSpacing));
    rSeqExpo.setBandwidthPerPixelPhaseEncode(static_cast<float>(dBandwidthPerPixelPE));

    {
        const double dReadFOV = static_cast<double>(rMrProt.sliceSeries().aFront().readoutFOV());
        const double dReadMatrix
            = static_cast<double>(rMrProt.kSpace().getlBaseResolution() * (rMrProt.kSpace().Interpolation2D() ? 2 : 1));
        const double dPixelSize = (dReadMatrix > 0.) ? (dReadFOV / dReadMatrix) : 0.;
        // Calculate pixel bandwidths [Hz/mm]
        // Note: dBandwidthPerPixelPE refers to the pixel size in the reconstructed image (quadratic, interpolated pixels)
        const double dPhaseBW = (dPixelSize > 0.) ? (dBandwidthPerPixelPE / dPixelSize) : 0.;

        rSeqExpo.setBandwidth_Hz_per_mm_PhaseEncode(static_cast<float>(dPhaseBW));
    }


    rSeqExpo.setPCAlgorithm(SEQ::PC_ALGORITHM_NONE);
    rSeqExpo.setRelevantReadoutsForMeasTime(static_cast<int32_t>(m_mySeqLoop.getNumberOfRelevantADCs()));
    if (lEffectiveEchoSpacing <= 0)
    {
        // Disable Maxwell correction if no valid effective echo spacing is provided (e.g. TGSE).
        rSeqExpo.setMaxwellCorrection(false);
    }
    else
    {
        rSeqExpo.setMaxwellCorrection(true);
        rSeqExpo.setMaxwellIntegralROGradient(static_cast<float>(m_EPIKernel.getGRO().getMaxwellIntegral()));
    }

    // Set B0 correction mode and add EPIPhaseCorrPEFunctor to Ice pipeline
    rSeqExpo.setB0Correction(iB0CorrectionMask);
    rSeqExpo.AddAdditionalIceProgramFileName ("%SiemensIceProgs%\\IceDecoratorEPIPhaseCorrPE" );

    // specify that measurement duration for image text is read from protocol
    lIceUserCtrlMask = rSeqExpo.getICEProgramParam(ICE_PROGRAM_PARA_USER_CTRL_MASK);
    lIceUserCtrlMask |= ICE_PROGRAM_MSK_USE_TOTALSCANTIME_FROM_UI;
    rSeqExpo.setICEProgramParam(ICE_PROGRAM_PARA_USER_CTRL_MASK, lIceUserCtrlMask);

    // ---------------------------------------------------------------------------
    // prepare exports concerning measurement time per measurement
    // ---------------------------------------------------------------------------
    if (rMrProt.getlRepetitions())
    {
        rSeqExpo.setMeasureTimeUsec(m_mySeqLoop.getMeasurementTimeUsec(rMrProt, rSeqLim, SecondMeas));
        rSeqExpo.setPreparingTimeInFirstMeasUSec(m_mySeqLoop.getPreparingTimeInFirstMeasUSec());
    }
    else
    {
        rSeqExpo.setMeasureTimeUsec(m_mySeqLoop.getMeasurementTimeUsec(rMrProt, rSeqLim, FirstMeas));
        rSeqExpo.setPreparingTimeInFirstMeasUSec(0);
    }

    // ---------------------------------------------------------------------------
    // set acquisition contrasts
    // ---------------------------------------------------------------------------
    setDICOMAcquisitionContrast(rSeqExpo);


    // In case of slice accelerated scans, recalculate the required times
    if (protFacade.isSliceAcceleration())
    {
        // the 200 fill time for calibration scans is required to have a 200 us gap between ADC and RF if SPAIR is selected
        // m_EPIKernel.setTRFill(200);
        m_mySeqLoop.setRunMode(MULTI_BAND);

        lCoolPauseExplicit = m_lCoolPauseTotal - m_lCoolPauseImplicit;
        lScanTimeSatsEtc   = m_mySeqLoop.getlScanTimeAllSats();   // SBB scan time
        lScanTimeBasic     = m_EPIKernel.getDurationPerRequest(); // Kernel scan time excluding mandatory fill time

        // The sequence asks for a pointer to the SPAIR SBB from SeqLoop.  Using this pointer it then uses the SBB to calculate the required inversion time.
        // This is set in the SBB from within the SBB.
        if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_Spair)
        {
            SeqBuildBlockOptfs*     pSBBOptfs     = m_mySeqLoop.getpOptfs();
            SeqBuildBlockOptfsPrep* pSBBOptfsPrep = m_mySeqLoop.getpOptfsPrep();

            // Add SPAIR time (depending on protocol parameters, e.g. TR) to lScanTimeSatsEtc
            if (pSBBOptfs == nullptr || !pSBBOptfs->calcSPIRTime(rMrProt, rSeqLim, rSeqExpo, lScanTimeBasic, lScanTimeSatsEtc, 0, SeqBuildBlockOptfs::SPIR_CALC_TYPE_EPI, pSBBOptfsPrep))
            {
                SEQ_TRACE_ERROR.print("The calc SPIR time has failed, prob due to an error in the tickle pulse");
                lStatus = MRI_SEQ_SEQU_ERROR;
                return SeverePrepareErrorReturn(lStatus);
            }
        }

        // lScanTimeBasic      - kernel time
        // lScanTimeSatsEtc    - time for fatsat/spatial sats
        // lCoolPauseExplicit  - Coolpause applied after EPI kernel
        m_mySeqLoop.setRunMode(SINGLE_BAND);
        lTR_SliceAccPrepMeas = (lScanTimeBasic + m_mySeqLoop.getlTRFillInConcat(0) + lScanTimeSatsEtc + lCoolPauseExplicit) * (rMrProt.sliceSeries().getlSize()) + m_mySeqLoop.getlTRFillEndInConcat(0);
        const long lTokTokTokTime    = m_mySeqLoop.getTokTokTokTime();
        const long lNoiseMeasTime    = m_mySeqLoop.getNoiseMeasTime();
        const long lFlashRefScanTime = m_mySeqLoop.getlDurationPATRefScan();
        const long lPrepScanMeasTime = lNMeasPrepScanTR * lTR_SliceAccPrepMeas;
        m_mySeqLoop.setRunMode(MULTI_BAND);
        rSeqExpo.setTotalMeasureTimeUsec(
            static_cast<double>(lNMeasImagingScanTR + lNMeasImagingScanTR_SecondMeas)
                * static_cast<double>(rMrProt.tr()[0])
            + static_cast<double>(lFlashRefScanTime + lPrepScanMeasTime + lNoiseMeasTime + lTokTokTokTime));

        if (rMrProt.getlRepetitions())
        {
            long lAdjPrepScans = getDiffusionAdjPrepScans();

            rSeqExpo.setPreparingTimeInFirstMeasUSec(static_cast<double>(
                lPrepScanMeasTime + lFlashRefScanTime + lNoiseMeasTime
                + (m_lSliceAccelDummyScans + m_lSliceAccelPhaseCorrScans + lAdjPrepScans) * rMrProt.tr()[0]));
            rSeqExpo.setMeasureTimeUsec(static_cast<double>(
                (lNMeasImagingScanTR - (m_lSliceAccelDummyScans + m_lSliceAccelPhaseCorrScans + lAdjPrepScans))
                * rMrProt.tr()[0]));
        }
        else
        {
            rSeqExpo.setPreparingTimeInFirstMeasUSec(0);
            rSeqExpo.setMeasureTimeUsec(rSeqExpo.getTotalMeasureTimeUsec());
        }
    }

    //  If GRE-PatRefScan or PACE triggering is used the RF exposure of the sequence is non-uniform.
    //  Hence it is necessary to work with 'RFBlockInfos' to avoid that SAR look ahead stops the sequence at run-time.
    if (rSeqLim.isContextNormal())
    {
        rSeqExpo.resetAllRFBlockInfos();

        MrProtocolData::SeqExpoRFInfo sRemainingEnergy_Ws(RFInfo);
        double                        dRemainingMeasureTime_us = rSeqExpo.getTotalMeasureTimeUsec();

        // *************************************
        // 1st block zero energy:
        // TokTokTok and / or noise measurement
        // *************************************
        if (rMrProt.intro() || (m_mySeqLoop.getPerformNoiseMeas() && (m_mySeqLoop.getNoiseMeasTime() > 0)))
        {
            MrProtocolData::SeqExpoRFInfo sEnergy1_Ws;
            int                           iTime1_us = 0;

            if (rMrProt.intro())
            {
                iTime1_us += static_cast<int>(m_mySeqLoop.getTokTokTokTime());
            }
            if (m_mySeqLoop.getPerformNoiseMeas())
            {
                iTime1_us += static_cast<int>(m_mySeqLoop.getNoiseMeasTime());
            }

            if (!rSeqExpo.addRFBlockInfo(
                    1E-6 * iTime1_us,                                             // duration in sec.
                    SeqExpoRFBlockInfo::VALUETYPE_ACTUAL,                         // exact energy
                    rMrProt.txSpec().nucleusInfoArray()[0].gettNucleus().c_str(), // nucleus type
                    sEnergy1_Ws,                                                  // energy
                    SeqExpoRFBlockInfo::VALUETYPE_ACTUAL                          // exact energy
                    ))
            {
                SEQ_TRACE_ERROR.print("ERROR: cannot add energy block");
                return MRI_SEQ_SEQU_ERROR;
            }
            sRemainingEnergy_Ws -= sEnergy1_Ws;
            dRemainingMeasureTime_us -= (double)iTime1_us;
        }

        // *************************************
        // 2nd block:
        // iPAT pre scan and SMS pre scan
        // *************************************
        if (((rMrProt.getsPat().getucPATMode() != SEQ::PAT_MODE_NONE) && (m_mySeqLoop.getePATRefScanLoopMode() != SeqLoop::PATRefScanLoopMode_NEVER)) || protFacade.isSliceAcceleration() || isFastGreRefScan(rMrProt))
        {
            MrProtocolData::SeqExpoRFInfo sEnergy2_Ws;
            int                           iTime2_us = 0;

            if (protFacade.isSliceAcceleration())
            {
                sEnergy2_Ws = RFInfoSeqLoopPerSlice * static_cast<double>(lNMeasPrepScanTR * lTotalSlices);

#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
                if (m_bApplyExtraCSat)
                    sEnergy2_Ws += calcEnergyOfExtraCSat(lNMeasPrepScanTR * lTotalSlices);
#endif

                m_EPIKernel.setRunMode(SINGLE_BAND);
                sEnergy2_Ws += m_EPIKernel.getRFInfoPerRequest() * static_cast<double>(lNMeasPrepScanTR * lTotalSlices);
                iTime2_us = static_cast<int>(lNMeasPrepScanTR * lTR_SliceAccPrepMeas);
            }
            // Add GRE ref scan energy and timing if required
            if (m_mySeqLoop.getePATRefScanLoopMode() != SeqLoop::PATRefScanLoopMode_NEVER )
            {
                sEnergy2_Ws += m_mySeqLoop.getTotalRFInfoPATRefScan();
                iTime2_us += static_cast<int>(m_mySeqLoop.getlDurationPATRefScan());
            }
            
            if (!rSeqExpo.addRFBlockInfo(
                    1e-6 * iTime2_us, SeqExpoRFBlockInfo::VALUETYPE_ACTUAL, rMrProt.getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str(), sEnergy2_Ws, SeqExpoRFBlockInfo::VALUETYPE_ACTUAL))
            {
                SEQ_TRACE_ERROR.print("ERROR: cannot add energy block");
                return MRI_SEQ_SEQU_ERROR;
            }
            sRemainingEnergy_Ws -= sEnergy2_Ws;
            dRemainingMeasureTime_us -= (double)iTime2_us;
        }

#ifdef SUPPORT_PACE
        // *************************************
        // 3rd block:
        // navigator learning phase
        // *************************************
        if ((rMrProt.NavigatorParam().getlRespComp() & (SEQ::RESP_COMP_TRIGGER | SEQ::RESP_COMP_TRIGGER_AND_FOLLOW)) != 0)
        {
#ifdef SUPPORT_PACE2
            MrProtocolData::SeqExpoRFInfo sEnergy3_Ws;
            m_mySeqLoop.getPACEStatPhaseRFInfo(sEnergy3_Ws, rMrProt, rSeqLim, rSeqExpo);
#else
            MrProtocolData::SeqExpoRFInfo sEnergy3_Ws(m_mySeqLoop.getPACEStatPhaseRFInfo(rMrProt));
#endif
            const int iTime3_us = static_cast<int>(m_mySeqLoop.getPACEStatPhaseDuration_us(rMrProt));

            if (!rSeqExpo.addRFBlockInfo(
                    1e-6 * iTime3_us, SeqExpoRFBlockInfo::VALUETYPE_ACTUAL, rMrProt.getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str(), sEnergy3_Ws, SeqExpoRFBlockInfo::VALUETYPE_ACTUAL))
            {
                SEQ_TRACE_ERROR.print("ERROR: cannot add energy block");
                return MRI_SEQ_SEQU_ERROR;
            }
            sRemainingEnergy_Ws -= sEnergy3_Ws;
            dRemainingMeasureTime_us -= (double)iTime3_us;
        }
#endif //  SUPPORT_PACE

        // *************************************
        // 4th block:
        // remaining time and energy
        // *************************************
        if (!rSeqExpo.addRFBlockInfo(
                1e-6 * dRemainingMeasureTime_us,
                SeqExpoRFBlockInfo::VALUETYPE_ACTUAL,
                rMrProt.getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str(),
                sRemainingEnergy_Ws,
                SeqExpoRFBlockInfo::VALUETYPE_ACTUAL))
        {
            SEQ_TRACE_ERROR.print("ERROR: cannot add energy block");
            return MRI_SEQ_SEQU_ERROR;
        }
    } //  Context Normal

    // ----------------------------------------------------------------------
    // Set echo train length (-> used to set corresponding DICOM attributes)
    // ----------------------------------------------------------------------
    rSeqExpo.setEchoTrainLength(static_cast<int32_t>(m_REOInfo.getEchoTrainLength()));
    rSeqExpo.setGradientEchoTrainLength(static_cast<int16_t>(m_REOInfo.getEchoTrainLength()));
#if (defined COMPILE_EP2D_SE) || (defined COMPILE_EP2D_DIFF)
    rSeqExpo.setRFEchoTrainLength(static_cast<int16_t>(1));
#else
    rSeqExpo.setRFEchoTrainLength(static_cast<int16_t>(0));
#endif

    // ---------------------------------------------------------------------------
    // determine default ice-program
    // ---------------------------------------------------------------------------
    bSuccess = m_mySeqLoop.setIceProgram(rMrProt, rSeqLim, rSeqExpo, SEQ::RS_LIN_IN_PAR);
    if (!bSuccess)
        mSBBErrGotoFinish(m_mySeqLoop, "m_mySeqLoop.setIceProgram failed");

    // ---------------------------------------------------------------------------
    // activate auto-crosscorrelation / across segments and primary mode
    // phase-correction in ICE program
    // ---------------------------------------------------------------------------
    rSeqExpo.setOnlinePhaseCorrectionAlgo(ICE_ONLINEPC_AUTOCROSSCORR_ACROSSSEGMENTS | ICE_ONLINEPC_PRIMARYMODE);
    // Pre-VD13A phase correction for EPI variants other than ep2d_diff and ep2d_se
    // rSeqExpo.setOnlinePhaseCorrectionAlgo(ICE_ONLINEPC_AUTOCORR|ICE_ONLINEPC_CROSSCORR_ACROSSSEGMENTS|ICE_ONLINEPC_PRIMARYMODE);

    // ---------------------------------------------------------------------------
    // set adaptive coil combine algorithm (CHARM 437817)
    // ---------------------------------------------------------------------------
    rSeqExpo.setAdaptiveCoilCombineAlgo(ACC_ALGO_PSNSENS);

    

    // ---------------------------------------------------------------------------
    // variant specific exports, ice program selection, etc.
    // ---------------------------------------------------------------------------
    lStatus = setVariantSpecificExports(rSeqLim, rMrProt, rSeqExpo);
    
    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
        return lStatus;


    IterativeDenoisingUIParameter::activateNoZeroFillingAndIcePATIfNeeded(rMrProt, rSeqExpo);

    // ---------------------------------------------------------------------------
    // specify submatrix for PC algorithm
    // ---------------------------------------------------------------------------
    if (rSeqExpo.getPCAlgorithm() != SEQ::PC_ALGORITHM_NONE)
    {
        rSeqExpo.setNoOfPhaseCorrLines(16 * rMrProt.kSpace().getlBaseResolution() / 256);
        rSeqExpo.setLinSlopeLength(16 * rMrProt.kSpace().getlBaseResolution() / 256);
        rSeqExpo.setNoOfPhaseCorrColumns(128 * rMrProt.kSpace().getlBaseResolution() / 256);
        rSeqExpo.setColSlopeLength(128 * rMrProt.kSpace().getlBaseResolution() / 256);
    }

    SEQ_TRACE_DEBUG_COND(!rSeqLim.isContextPrepForBinarySearch(), "rSeqExpo.getMaxReceiverChannels() = %d", rSeqExpo.getMaxReceiverChannels());

#ifdef SUPPORT_iPAT_a_ep2d
    //---------------------------------------------------------------------------
    // call fPATPrepPost, which checks some PAT related restrictions
    //---------------------------------------------------------------------------
    lStatus = fPATPrepPost(rMrProt, rSeqLim, rSeqExpo, &m_REOInfo);

#ifdef BOLD
    // Forbid channel reduction for GRAPPA
    if (rMrProt.PAT().getucPATMode() == SEQ::PAT_MODE_GRAPPA)
    {
        // Read ICE control flag
        long lIceControl = rSeqExpo.getICEProgramParam(ICE_PROGRAM_PARA_CTRL_MASK);

        // Set the "no channel reduction" bit
        lIceControl = lIceControl | ICE_PROGRAM_MSK_iPAT_NO_CHANNEL_REDUCTION;

        // Write back
        rSeqExpo.setICEProgramParam(ICE_PROGRAM_PARA_CTRL_MASK, lIceControl);
    }
#endif // endif BOLD

    if (!rSeqLim.isContextPrepForBinarySearch())
    {
        CheckStatusPB(lStatus, "fPATPrepPost");
    }
    else
    {
        CheckStatusB(lStatus);
    }
#endif

    // ---------------------------------------------------------------------------
    // Check the required TE, TI and TR values against the current ones, adapt
    // those values, if rSeqLim.isContextPrepForMrProtUpdate() is true.
    // ---------------------------------------------------------------------------
    if (!m_pUI->fEPIStdUICheckNeededTETITR(rMrProt, rSeqLim, (int32_t)lNeededTE, (int32_t)lNeededTI, (int32_t)lNeededTR, getTEContrastIndex( rMrProt.getlContrasts() )))
    {
        lStatus = MRI_SEQ_SEQU_ERROR; // no specific error-message available,
        // not nice but SEQU__NEGATIVE_TE_FILL would also be misleading
    }

    // check acq window when resp trigger is used
    if (!checkAcqWindowForRespTriggering(rMrProt))
    {
        lStatus = MRI_SEQ_SEQU_ERROR;
    }

    // ----------------------------------------------------------------------------------
    // Protocol related final actions for static field correction SFC
    // ----------------------------------------------------------------------------------
    if (!rSeqLim.isContextPrepForBinarySearch()
        && (rMrProt.getucStaticFieldCorrection()))
    {
        // A consistency check integrated into SFC requires information about
        // the actual adjustment volume. Since the corresponding protocol
        // entries will not always contain valid data, we add this manually.
        // Note: In a product implementation, one might consider implementing
        //       an Ice interface which provides access to the actual
        //       adjustment volume (which is used by AdjFre).
        if (!rMrProt.getsAdjData().getuiAdjVolumeValid())
        {
            MrProtJustVol sMrProtAdjVol;
            Slice         sAdjSlice;
            AdjAccessIF   sAdjAccessIF;

            sMrProtAdjVol.calcAdjSlice(rMrProt.getProtData(), sAdjSlice, sAdjAccessIF.isLocalShimRequested(rMrProt));

            // Note: Do not touch the 'valid' flag, since this might affect adjustments!
            rMrProt.getsAdjData().getsAdjVolume().setdPhaseFOV(sAdjSlice.getdPhaseFOV());
            rMrProt.getsAdjData().getsAdjVolume().setdReadoutFOV(sAdjSlice.getdReadoutFOV());
            rMrProt.getsAdjData().getsAdjVolume().setdThickness(sAdjSlice.getdThickness());
            rMrProt.getsAdjData().getsAdjVolume().setdInPlaneRot(sAdjSlice.getdInPlaneRot());
            rMrProt.getsAdjData().getsAdjVolume().getsNormal().setdSag(sAdjSlice.getsNormal().getdSag());
            rMrProt.getsAdjData().getsAdjVolume().getsNormal().setdCor(sAdjSlice.getsNormal().getdCor());
            rMrProt.getsAdjData().getsAdjVolume().getsNormal().setdTra(sAdjSlice.getsNormal().getdTra());
            rMrProt.getsAdjData().getsAdjVolume().getsPosition().setdSag(sAdjSlice.getsPosition().getdSag());
            rMrProt.getsAdjData().getsAdjVolume().getsPosition().setdCor(sAdjSlice.getsPosition().getdCor());
            rMrProt.getsAdjData().getsAdjVolume().getsPosition().setdTra(sAdjSlice.getsPosition().getdTra());
        }
    }

    dumpSliceAdjData(rSeqLim, rMrProt);

    if (!IterativeDenoisingUIParameter::SeqExpoOkForIterativeDenoising(rMrProt, rSeqExpo, rSeqLim))
    {
        return MRI_SEQ_SEQU_ERROR;
    }

    // ---------------------------------------------------------------------------
    // finished
    // ---------------------------------------------------------------------------

    return (lStatus);

    // ---------------------------------------------------------------------------
    // in the error case:
    // ---------------------------------------------------------------------------
FINISHED:
    if (NLS_SEVERITY(lStatus) == NLS_SUCCESS)
        lStatus = MRI_SEQ_SEQU_ERROR; // make sure that an ERROR code is returned

    m_pUI->fEPIStdResetSolveHandlerControlTETITR(); // even if TE,TI or TR is changed, we can't help

    return (lStatus);
}

//  --------------------------------------------------------------------------
//
//  Name        :  Ep2d::check
//
//  Description :
/// \brief         Check of the sequence for gradient stimulation
///
///                This method is called by the framework prior to a
///                 measurement on the host to ensure, that
///                 - no gradient overflow occurs
///                 - the stimulation will not exceed the threshold
///
//  Return      :  NLS status
//
//  --------------------------------------------------------------------------
NLSStatus Ep2d::check(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, SEQCheckMode*)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    //  ----------------------------------------------------------------------
    //
    //  Sequence check:
    //     - max. gradient amplitude
    //     - GSWD look ahead
    //
    //  ----------------------------------------------------------------------
#undef DEBUG_ORIGIN
#define DEBUG_ORIGIN 0x00000400

    // pointer for MrProtFacade for easier protocol queries
    MrProtFacade protFacade(rMrProt);

#if (defined ZOOM_2DRF)
    // determine some parameters from the *first* PTXRFPulse (index 0)
    const size_t                      uiCurVol        = 0;
    const long                        lPTXRFPulseSize = rMrProt.getsTXSPEC().getaPTXRFPulse().size();
    MrProtocolData::PTXTrajectoryType eExcType        = MrProtocolData::PTXTrajectoryType_EPI_1D;
    if (lPTXRFPulseSize > 0)
    {
        eExcType = rMrProt.getsTXSPEC().getaPTXRFPulse()[uiCurVol].getlTrajectoryType();
    }

    if ((rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED) && (eExcType == MrProtocolData::PTXTrajectoryType_EPI_1D))
    {
        // trace 2DRF pulse details to log file
        std::string strTmp;
        SBB2DPtx*   pSBBExcite = dynamic_cast<SBB2DPtx*>(m_EPIKernel.getExcitationPointer());
        if (pSBBExcite != nullptr)
        {
            strTmp = pSBBExcite->getToolTipInfo();
            SEQ_TRACE_ALWAYS.print("ZOOM_2DRF Info:\n%s", strTmp.c_str());
        }
    }
#endif

    // prepare bookeeping for the dynamic field correction
    if (protFacade.isBookkeepingConditionForDFC())
    {
        prepareDFCBookkeeping();
    }

    //---------------------------------------------------------------------------
    // Set the looping parameters
    //---------------------------------------------------------------------------
    m_mySeqLoop.setlinesToCheck(1);
    m_mySeqLoop.setlineNoToCheck(0, 0);

    m_mySeqLoop.setpartitionsToCheck(1);
    m_mySeqLoop.setparNoToCheck(0, 0);

    //---------------------------------------------------------------------------
    // Execute the check loops
    //---------------------------------------------------------------------------
    // the check should always run in single band mode as the gradients are the
    // same for single band and multi band
    // MZ: Todo: Really? Slice blips are not considered
    if (protFacade.isSliceAcceleration() && isIIRSchemeStandard(rMrProt, rSeqLim))
    {
        m_mySeqLoop.setRunMode(SINGLE_BAND);
        m_EPIKernel.setRunMode(SINGLE_BAND);
    }

    const bool bSuccess = m_mySeqLoop.check(rMrProt, rSeqLim, rSeqExpo, m_asSLC, m_EPIKernel.getReadOutAddress());
    if (!bSuccess)
        mSBBErrGotoFinish(m_mySeqLoop, " m_mySeqLoop.check");

    //---------------------------------------------------------------------------
    // ready.
    //---------------------------------------------------------------------------
FINISHED:
    return (lStatus);
}

//  --------------------------------------------------------------------------
//
//  Name        :  Ep2d::run
//
//  Description :
///     \brief     Execution of the sequence
//
//  Return      :  NLS status
//
//  --------------------------------------------------------------------------
NLSStatus Ep2d::run(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

#ifdef EP2D_SE_MRE
    if (SeqUT.isUnitTestActive())
    {
        if (rMrProt.intro())
        {
            if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_Spair)
            {
                SeqUT.SetExpectedNotOk(lExtTriggerSpacingErr, RTEB_ORIGIN_fSBBOptfs, 1, "First external trigger too late due to the TokTokTokTime.");
            }
            else
            {
                SeqUT.SetExpectedNotOk(lExtTriggerSpacingErr, RTEB_ORIGIN_fSEQRunKernel, 1, "First external trigger too late due to the TokTokTokTime.");
            }
        }
    }
#endif

    // enable real-time processing mode for some systems
    if (m_bIsRealtimeProcessingEnabled)
    {
        getRTController().setRealtimeProcessing(true);
        getRTController().setRealtimeProcessingMargin(m_debugSettings.getDefaultSetting<double>("EP2D_BOLD/RealtimeProcessingMargin_s", 0.010)); // [s]
        SEQ_TRACE_ALWAYS.print("getRTController().isRealtimeProcessingEnabled()/Margin: %d/%f s", getRTController().isRealtimeProcessingEnabled(), getRTController().getRealtimeProcessingMargin());
    }

    //  ----------------------------------------------------------------------
    //
    //  Execution of sequence timing
    //
    //  ----------------------------------------------------------------------
#undef DEBUG_ORIGIN
#define DEBUG_ORIGIN 0x00002000

    // pointer for MrProtFacade for easier protocol queries
    MrProtFacade protFacade(rMrProt);

    // prepare bookkeeping for the dynamic field correction
    if (protFacade.isBookkeepingConditionForDFC())
    {
        prepareDFCBookkeeping();
    }

    // ---------------------------------------------------------------------------
    // provide unit test with special information:
    // ---------------------------------------------------------------------------
#ifdef WIN32

    if (SeqUT.isUnitTestActive() && SMSProperties::isSMS(rMrProt))
    {
        // SeqUT uses the EPI factor from MrProt to calculate
        // expected SMS-related moments along the slice-encoding axis - but
        // EPI sequences do not use or set this protocol value. Thus, we
        // write it to the protocol here in order to assist the SeqUT.
        // Note: In principle, one could update the EPI factor in MrProt
        // when changing corresponding UI parameters. However, a value
        // other than 1 currently leads to unexpected UI side effects,
        // probably due to libUILink internal handling.
        rMrProt.getsFastImaging().setlEPIFactor(static_cast<int32_t>(m_REOInfo.getEchoTrainLength()));
    }

    if (SeqUT.isUnitTestActive())
    {
        configureDiffusionSpecificSeqUTSettings(rSeqLim, rMrProt);
    }
#endif

    // ---------------------------------------------------------------------------
    // Disable unit test TR calculation when using multiple
    // concatenations with long TR triggering mode
    // ---------------------------------------------------------------------------
#ifdef WIN32

    if (SeqUT.isUnitTestActive())
    {
#ifndef SUPPORT_IIR
        if (m_mySeqLoop.isLongTRTrigMode() && (rMrProt.concatenations() > 1))
        {
            SeqUT.DisableTestCase(lTRClockErr, RTEB_ClockCheck, "The unit test TR calculations are not valid for multiple concatenations in long TR triggering mode");
        }
#endif

#ifdef EP2D_SE_MRE
        SeqUT.setSizeOfDimSet(NMAXMEGS);
        SeqUT.setSizeOfDimFree(MEG_DEFAULT_TRIGGER_STEPS);
        SeqUT.setRequiredExtTrigSpacingForMRE_us(static_cast<long>(MRE_DRIVER_DEFFAULT_NUMBEROFBURSTS * 1.0E6 / MRE_DRIVER_DEFAULT_FREQUENCY));
#endif

#ifdef SUPPORT_PACE
        if (rMrProt.NavigatorParam().getlRespComp() == SEQ::RESP_COMP_TRIGGER || rMrProt.NavigatorParam().getlRespComp() == SEQ::RESP_COMP_TRIGGER_AND_FOLLOW)
        {
            SeqUT.DisableTestCase(lTRClockErr, RTEB_ClockCheck, "PACE respiratory triggering: The TR is determined by the breathing cycle and not by the protocol TR");
        }
#endif
    }

#endif

    //---------------------------------------------------------------------------
    // initialize m_alPrepScanCounter
    //---------------------------------------------------------------------------
#ifdef SUPPORT_iPAT_a_ep2d
    {
        std::fill(begin(m_alPrepScanCounter),end(m_alPrepScanCounter),0);
    }
#endif

    //---------------------------------------------------------------------------
    // initialize alASLPrepScanCounter, alASLScanCounter
    //---------------------------------------------------------------------------
    // TODO: Replace this with stl::<vector>
#ifdef ASL
    {
        std::fill(begin(m_alASLPrepScanCounter),end(m_alASLPrepScanCounter),0);
        std::fill(begin(m_alASLlScanCounter),end(m_alASLlScanCounter),0);
    }
#endif

    // ---------------------------------------------------------------------------
    // initialize osc-bit control-flags
    // ---------------------------------------------------------------------------
    initializeOscBitFlags();

    // Initialize runmode: starting with single band mode is important for SeqLoop,
    // since SBBPatRefScan relies on this setting. For the succeeding EPI acquisitions,
    // SeqLoopMultiBand::runOuterSliceLoop takes care of applying the desired mode
    // (which gets applied to EPIKernel later on).
    if (SMSProperties::isSMS(rMrProt))
    {
        m_mySeqLoop.setRunMode(SINGLE_BAND);
        m_EPIKernel.setRunMode(SINGLE_BAND);
    }


    //---------------------------------------------------------------------------
    // Initialization of the unit test function
    //---------------------------------------------------------------------------
    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunStart, 0, 0, 0, 0, 0);

#ifdef EP2D_SE_MRE
    // Initialize the running trigger: "-1" signals that the trigger has never been sent
    m_mySeqLoop.setlExtTrigStartTime(-1);
#endif

    //---------------------------------------------------------------------------
    // Execute the measurement loops
    //---------------------------------------------------------------------------
    const bool bSuccess = m_mySeqLoop.run_new(rMrProt, rSeqLim, rSeqExpo, m_asSLC, m_EPIKernel.getReadOutAddress());
    if (!bSuccess)
        mSBBErrGotoFinish(m_mySeqLoop, " m_mySeqLoop.run");

    //---------------------------------------------------------------------------
    // ready.
    //---------------------------------------------------------------------------
FINISHED:
    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunFinish, 0, 0, 0, 0, 0);

    return (lStatus);
}

//   --------------------------------------------------------------------------
//
//   Name        :  Ep2d::runKernel
//
//   Description :
///                 Executes the basic timing of the real-time sequence.
//
//   Return      :  NLS status
//
//   --------------------------------------------------------------------------
NLS_STATUS Ep2d::runKernel(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lKernelMode, long lSlice, long /* lPartition */, long lCurrentShot /*long lLine*/) // MZ: Todo: Align signature
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    //  ----------------------------------------------------------------------
    //
    //  Excecution of sequence kernel
    //
    //  ----------------------------------------------------------------------
#undef DEBUG_ORIGIN
#define DEBUG_ORIGIN 0x00004000

    // pointer for MrProtFacade for easier protocol queries
    MrProtFacade protFacade(rMrProt);

    // Send data to SFC functor
    const bool isSendSliceAdjustData = ( lKernelMode == KERNEL_IMAGE ) && protFacade.isSliceAdj();
    m_EPIKernel.setSendSliceAdjustData(isSendSliceAdjustData);

    //---------------------------------------------------------------------------
    // to send phase-FT flags correctly:
    //---------------------------------------------------------------------------
    m_REOInfo.setIsLastAcquisition(m_EPIKernel.getReadOutAddress()->getMDH().getCacq() == rMrProt.averages() - 1);

    //---------------------------------------------------------------------------
    // send clock check event for unit test (SBBBinomialPulses does not do it)
    //---------------------------------------------------------------------------
    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ClockCheck, 27, -1, m_asSLC[lSlice].getSliceIndex(), 0, 0);

    //---------------------------------------------------------------------------
    // remaining kernel configurations
    //---------------------------------------------------------------------------
#ifdef ASL
    // disable seqloop TR fill to shift slices together, add lTRfillASL after last slice
    //
    // The problem is that running EPI without fat sat and zero fill may cause the ADC and RF pulse to be too close
    //
    long lTRfillASL = 0;

    // This is a conservative setting, does not account for time after last ADC and GS ramp up time

    MrRXSpec MrSpecWrapper(rMrProt.getsRXSPEC());
    long     lMinDurationBetweenADCAndRF = fSDSRoundUpGRT(static_cast<long>(SysProperties::getMinDurationBetweenReadoutAndRFPulse(MrSpecWrapper.realDwellTime()[0] / 1000.0)));
    lTRfillASL                           = m_mySeqLoop.getlTRFill() - lMinDurationBetweenADCAndRF;
    m_EPIKernel.setTRFill(lMinDurationBetweenADCAndRF);

#else
    m_EPIKernel.setTRFill(m_mySeqLoop.getlTRFill());
#endif

    // set external trigger and Osc Bit
    setTriggerAndOscBit(lKernelMode);

#ifdef EP2D_SE_MRE
    // Free counter for both phase offsets and MEG instances
    // - ToDo: consistent settings should really happen inside the (MRE-) kernel
    // -> distribute to ICE SET and IDA dimensions
    m_EPIKernel.getReadOutAddress()->getMDH().setCset(static_cast<unsigned short>(m_mySeqLoop.getlFreeLoopCounter() % NMAXMEGS));
    m_EPIKernel.getReadOutAddress()->getMDH().setCida(static_cast<unsigned short>(m_mySeqLoop.getlFreeLoopCounter() / NMAXMEGS));

#endif

    //----------------------------------------------------------------------------
    // adapt Crep counter and Mdh flags when using ASL reference scan
    // NOTE: needs to be done BEFORE gPaceFeedback.SyncAndIncorporateFeedback !
    //----------------------------------------------------------------------------
    // if this is the first slice, execute m_ASL_SBB (which includes FatSat+Spoiler)
    //----------------------------------------------------------------------------
#ifdef ASL
    if (m_bEnableFirstPrepScanAsM0Scan)
    {
        if (lKernelMode == KERNEL_PREPARE)
        {
            if (m_alASLPrepScanCounter[lSlice] == 0)
            {
                //-------------------------------------------------------------------------------------
                // enable readout for very first M0 scan (this is even before PAT: lPrepScansNoPATRefScans)
                //-------------------------------------------------------------------------------------
                fRTSetReadoutEnable(1);
                m_mySeqLoop.setRepetitionValueForMdh(0);
                m_EPIKernel.getReadOutAddress()->getMDH().setCrep((unsigned short)(0));

                if (false)
                {
                    SEQ_TRACE_ALWAYS.print("%li: => m_alASLPrepScanCounter[lSlice]", m_alASLPrepScanCounter[lSlice]);
                    SEQ_TRACE_ALWAYS.print("%li: => lSlice", lSlice);
                    SEQ_TRACE_ALWAYS.print("%li: => m_mySeqLoop.getNumberOfRelevantADCs()", m_mySeqLoop.getNumberOfRelevantADCs());
                    SEQ_TRACE_ALWAYS.print("%i: => rMrProt.sliceSeries().getlSize()-1", rMrProt.sliceSeries().getlSize() - 1);
                }

                // Correction for physio triggering - M0 has to be counted properly
                // CHARM MR_00397864
                m_EPIKernel.getReadOutAddress()->setRelevantForMeasTime(false); // Is this necessarry?

                if (m_mySeqLoop.getNumberOfRelevantADCs() > 1 && m_mySeqLoop.getNumberOfRelevantADCs() >= rMrProt.sliceSeries().getlSize())
                {
                    if (lSlice == (rMrProt.sliceSeries().getlSize() - 1))
                    {
                        m_EPIKernel.setADCRelevant(false, true);
                    }
                    else
                    {
                        m_EPIKernel.setADCRelevant(true, false);
                    }
                }
                else if (m_mySeqLoop.getNumberOfRelevantADCs() > 1 && m_mySeqLoop.getNumberOfRelevantADCs() < rMrProt.sliceSeries().getlSize())
                {
                    if (lSlice + 1 >= rMrProt.sliceSeries().getlSize() - m_mySeqLoop.getNumberOfRelevantADCs() + 2)
                    {
                        if (lSlice == (rMrProt.sliceSeries().getlSize() - 1))
                        {
                            m_EPIKernel.setADCRelevant(true, true);
                        }
                        else
                        {
                            m_EPIKernel.setADCRelevant(true, false);
                        }
                    }
                }

                // Last slice - start counting measurement time for the physio trigger
                if (lSlice >= (rMrProt.sliceSeries().getlSize() - 1))
                {
                    m_EPIKernel.getReadOutAddress()->getMDH().setLastScanInMeas(true);
                    m_EPIKernel.getReadOutAddress()->getMDH().setLastScanInConcat(true);

                    if (m_mySeqLoop.getNumberOfRelevantADCs() > 0)
                    {
                        m_EPIKernel.getReadOutAddress()->setRelevantForMeasTime(true);
                    }
                }
            }

            m_alASLPrepScanCounter[lSlice]++;
        }
        else if (lKernelMode == KERNEL_IMAGE)
        {
            // increase counter first (start at crep=1): first crep is M0 reference scan
            m_alASLlScanCounter[lSlice]++;

            m_mySeqLoop.setRepetitionValueForMdh(m_alASLlScanCounter[lSlice]);
            m_EPIKernel.getReadOutAddress()->getMDH().setCrep((unsigned short)(m_alASLlScanCounter[lSlice]));
        }
    }
    else if (lKernelMode == KERNEL_IMAGE)
    {
        m_mySeqLoop.setRepetitionValueForMdh(m_alASLlScanCounter[lSlice]);
        m_EPIKernel.getReadOutAddress()->getMDH().setCrep((unsigned short)(m_alASLlScanCounter[lSlice]));
        m_alASLlScanCounter[lSlice]++; // We want label state to be first
    }
#endif // ASL

    //-------------------------------------------------------------------------------------
    // synchronize PACE feedbacks and sequence here; i.e. for every new measurement starting
    //-------------------------------------------------------------------------------------
#ifdef PACE3D
    if (!m_PaceFeedback.SyncAndIncorporateFeedback(rMrProt, rSeqLim, rSeqExpo, m_asSLC, &m_EPIKernel, lKernelMode))
    {
        SEQ_TRACE_ERROR.print("PaceFeedback.SyncAndIncorporateFeedback failed.");
        return MRI_SEQ_SEQU_ERROR;
    }
#endif

    //-------------------------------------------------------------------------------------
    // execute m_OptPTXVolume module
    //-------------------------------------------------------------------------------------

#if defined ZOOM_2DRF
    if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
    {
        lStatus = m_OptPTXVolume.run(rMrProt, rSeqLim, rSeqExpo, &m_asSLC[lSlice]);
        if (!lStatus)
        {
            SEQ_TRACE_ERROR.print("m_OptPTXVolume.run failed: 0x%lx", lStatus);
            return lStatus;
        }
    }
#endif

    //-------------------------------------------------------------------------------------
    // execute ASL module before first slice and handle alternate label/control state
    //-------------------------------------------------------------------------------------
#ifdef ASL
    if (lSlice == 0)
    {
        // set ASL to label/control scan -> first is M0 ASL_REFERENCE, 2nd is "label", 3rd is "control" etc etc
        //    in PostProc: ASL_LABEL scan is subtracted from ASL_CONTROL (3-2, 5-4, 7-6 ...)

        if (lKernelMode == KERNEL_IMAGE)
        {
            if (((m_EPIKernel.getReadOutAddress()->getMDH().getCrep() % 2) != 0) || rMrProt.getlRepetitions() == 0)
            {
                m_ASL_SBB.setLabelState(LIBASL::SYMBOL_ASL::ASL_LABEL);
                // SEQ_TRACE_ALWAYS.print("%s setLabelState(ASL_LABEL) crep=%d, cslc=%d (LabelingScheme=%d)", m_EPIKernel.getReadOutAddress()->getMDH().getCrep(), lSlice, iLabelingScheme);
            }
            else
            {
                m_ASL_SBB.setLabelState(LIBASL::SYMBOL_ASL::ASL_CONTROL);
                // SEQ_TRACE_ALWAYS.print("%s setLabelState(ASL_CONTROL) crep=%d, cslc=%d (LabelingScheme=%d)", m_EPIKernel.getReadOutAddress()->getMDH().getCrep(), lSlice, iLabelingScheme);
            }
        }

        else // PrepScans
        {
            // only a single reference scan
            // TODO: Should this be m_alASLPrepScanCounter[lSlice] < 1
            if (m_bEnableFirstPrepScanAsM0Scan && m_alASLPrepScanCounter[lSlice] <= 1)
            {
                m_ASL_SBB.setLabelState(LIBASL::SYMBOL_ASL::ASL_REFERENCE);
            }
            else // this is a prepscan
            {
                m_ASL_SBB.setLabelState(LIBASL::SYMBOL_ASL::ASL_CONTROL);
            }
        }

        if (!m_ASL_SBB.run(rMrProt, rSeqLim, rSeqExpo, &m_asSLC[lSlice]))
        {
            SEQ_TRACE_ERROR.print("m_ASL_SBB.run failed.");
            return (m_ASL_SBB.getNLSStatus());
        }

#if (defined DEBUG) && (defined ASL_DEBUG_COUT)
        if (m_EPIKernel.getReadOutAddress()->getMDH().getCrep() == 1)
            std::cout << m_ASL_SBB;
#endif
    }
    else
    { // FATSAT_ALL_SLICES
        // run fatsat for each slice EXCEPT the first slice, which is included in the ASL SBB
        if (rMrProt.preparationPulses().getucFatSatMode() == SEQ::FAT_SAT_STRONG)
        {
            if (!m_CSatFat.run(rMrProt, rSeqLim, rSeqExpo, &m_asSLC[lSlice]))
            {
                SEQ_TRACE_ERROR.print("m_CSatFat.run failed.");
                return MRI_SEQ_SEQU_ERROR;
            }
            if (!m_SpoilGrad.run(rMrProt, rSeqLim, rSeqExpo, &m_asSLC[lSlice]))
            {
                SEQ_TRACE_ERROR.print("SpoilGrad.run failed.");
                return MRI_SEQ_SEQU_ERROR;
            }
        }
    }
#endif // ASL

    //-------------------------------------------------------------------------------------
    // Settings for diffusion module.
    //-------------------------------------------------------------------------------------
    setDiffusionAdjustmentScans(rMrProt, lSlice, lKernelMode);

    
    //--------------------------------------------------------------------------
    // Set diffusion and repetitions loop counter for diffusion module
    //--------------------------------------------------------------------------
    if (!setDiffusionLoopCounters(lSlice, lKernelMode))
        return MRI_SEQ_SEQU_ERROR;

    
    // Each volume is handled internally as a separate repetition => set
    // LastScanInMeas flag correspondingly. Exception: PatRefScans.
    if (m_bSequentialVolumeAcquisition)
    {
        setLastScanInMeasFlagForB0Correction(rMrProt, lSlice, lCurrentShot);
    }

    //-------------------------------------------------------------------------------------
    // Settings for slice accelerated reference scans.
    // For ep2d_diff the b-value is set to zero for PAT ref scans
    //-------------------------------------------------------------------------------------
#ifdef SUPPORT_iPAT_a_ep2d
    if (SMSProperties::isSMS(rMrProt) && (lKernelMode == KERNEL_PREPARE))
    {
        long lFirstSliceAccRefScan = m_lInitialDummyScans + m_lPhaseCorrPrepScans;
        long lLastSliceAccRefScan  = lFirstSliceAccRefScan + m_lSliceAccelRefScans - 1;

        // Settings for Slice Acc ACS scans
        if ((m_alPrepScanCounter[lSlice] >= lFirstSliceAccRefScan)
            && (m_alPrepScanCounter[lSlice] <= lLastSliceAccRefScan))
        {
            // Enable readouts in KERNEL_PREPARE mode
            fRTSetReadoutEnable(1);
            // Calculate shot-index corresponding to the current preparing scan
            if (m_REOInfo.isMultiShot())
            {
                lCurrentShot = m_alPrepScanCounter[lSlice] - (m_lInitialDummyScans + m_lPhaseCorrPrepScans);
            }
            disableDiffusionForPrepScan();
        }
    }
#endif

    // set last scan in meas when all slices have been acquired for one specific repetition counter
    // this is required if dynamic field correction is used with multiple concats
    if (protFacade.isBookkeepingConditionForDFC())
    {
        if (NLS_SEVERITY(setLastScanInMeasFlagForDFCBookkeeping(rMrProt, lKernelMode)) != NLS_SUCCESS)
            return MRI_SEQ_SEQU_ERROR;
    }


    //---------------------------------------------------------------------------
    // Disable diffusion encoding for external phase correction scan
    //---------------------------------------------------------------------------
    handleExternalPhaseCorrectionRun(rMrProt, lKernelMode, lSlice);


    //-------------------------------------------------------------------------------------
    // Settings for PAT reference scans.
    // For ep2d_diff the b-value is set to zero for PAT ref scans
    //-------------------------------------------------------------------------------------
#ifdef SUPPORT_iPAT_a_ep2d

    long lFirstPATRefScan = m_lInitialDummyScans + m_lPhaseCorrPrepScans;
    long lLastPATRefScan  = lFirstPATRefScan + m_lPATRefScans - 1;

    // In case of slice acceleration shift PAT reference scans
    if (protFacade.isSliceAcceleration())
    {
        lFirstPATRefScan += m_lSliceAccelRefScans;
        lLastPATRefScan += m_lSliceAccelRefScans; // lFirstPATRefScan + m_lPATPrepScans - 1;
    }

    if (m_REOInfo.isPATActive() && m_REOInfo.getPATAccelerationFactorPE() > 1 && (rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI))
    {
        if ((lKernelMode == KERNEL_PREPARE) && (m_alPrepScanCounter[lSlice] >= lFirstPATRefScan) && (m_alPrepScanCounter[lSlice] <= lLastPATRefScan))
        {
            fRTSetReadoutEnable(1);

            const long lPATRefScanCounterInSegment = m_alPrepScanCounter[lSlice] - lFirstPATRefScan;

            m_REOInfo.setPATReorderIndexOffsetForRefScans(lPATRefScanCounterInSegment);
            m_EPIKernel.setCounterInSegmentForEchoShifting(lPATRefScanCounterInSegment);
            m_EPIKernel.setBlindImagingADCs(m_REOInfo.getPATBlindADCsBeforeRefScans(), m_REOInfo.getPATBlindADCsAfterRefScans());

#ifdef ASL
            // The reference scan must have repetition 0
            m_EPIKernel.getReadOutAddress()->getMDH().setCrep(0);
#endif

            // disable diffusion gradients for PAT reference scans
            disableDiffusionForPrepScan();

        }
        else
        {
            // Why we don't use
            // m_EPIKernel.setCounterInSegmentForEchoShifting(m_REOInfo.getPATRefCounterInSegmentWithKSCenter());
            // in the following line is explained where m_EPIKernel.setUseEchoShifting is used in fSEQPrep.
            //
            m_REOInfo.setPATReorderIndexOffsetForImagingScans();
            m_EPIKernel.setCounterInSegmentForEchoShifting(0);
            m_EPIKernel.setBlindImagingADCs(0, 0);
        }
    }
    else
    {
        long lCounterInSegmentForEchoShifting = 0;

        // Multi-shot acquisition: consider echo shifting
        if (m_REOInfo.isMultiShot())
        {
            if (m_REOInfo.isPATActive())
            {
                lCounterInSegmentForEchoShifting = lCurrentShot % m_REOInfo.getPATLinesPerSegment();
            }
            else
            {
                lCounterInSegmentForEchoShifting = lCurrentShot % m_REOInfo.getLinesPerSegment();
            }
        }

        // Shift reordering table offset to first contrast (ECO) of current shot
        m_REOInfo.setReorderIndexOffset(0, 0, lCurrentShot, 0);
        // Configure echo-shifting of EPI kernel for current shot
        m_EPIKernel.setCounterInSegmentForEchoShifting(lCounterInSegmentForEchoShifting);
        // All ADC's of EPI kernel are used
        m_EPIKernel.setBlindImagingADCs(0, 0);
    }

#endif

#ifdef EP2D_MS
    //---------------------------------------------------------------------------
    // Write shot index to IceProgramPara section of MDH
    //---------------------------------------------------------------------------
    m_EPIKernel.getReadOutAddress()->getMDH().setIceProgramPara(
        ICEPROGRAMPARA_SHOT_INDEX, static_cast<uint16_t>(lCurrentShot));
    m_EPIKernel.getReadOutAddress()->getMDH().setCset(0); // Must not be overwritten in case of EP2D_SE_MRE 
#endif

    if (protFacade.isSliceAcceleration())
    {
        const SliceAccelRFRunMode eRunMode = m_mySeqLoop.getRunMode();
        m_EPIKernel.setRunMode(eRunMode);

        if (eRunMode == MULTI_BAND)
        {
            m_REOInfo.setRunModeMultiBand();
        }
        else
        {
            m_REOInfo.setRunModeSingleBand();
        }

        if (eRunMode == MULTI_BAND || SMSProperties::getMultiBandFactor(rMrProt) == 1)
            m_EPIKernel.setFlagPCforRTFeedback(true);
        else
            m_EPIKernel.setFlagPCforRTFeedback(false);
    }

        //---------------------------------------------------------------------------
    // Disable unsupported test cases for current kernel execution
    //---------------------------------------------------------------------------
#ifdef WIN32
    if (SeqUT.isUnitTestActive())
    {
#ifndef SUPPORT_IIR
        if (m_mySeqLoop.isLongTRTrigMode() && (rMrProt.concatenations() > 1))
        {
            SeqUT.DisableTestCase(
                lTRClockErr,
                RTEB_ClockCheck,
                "The unit test TR calculations are not valid for multiple concatenations in long TR triggering mode");
        }
#endif

#ifdef SUPPORT_PACE
        if (rMrProt.NavigatorParam().getlRespComp() == SEQ::RESP_COMP_TRIGGER
            || rMrProt.NavigatorParam().getlRespComp() == SEQ::RESP_COMP_TRIGGER_AND_FOLLOW)
        {
            SeqUT.DisableTestCase(
                lTRClockErr,
                RTEB_ClockCheck,
                "PACE respiratory triggering: The TR is determined by the breathing cycle and not by the protocol TR");
        }
#endif

        if (protFacade.isSliceAcceleration())
        {
            if (m_mySeqLoop.getRunMode() != MULTI_BAND)
            {
                // Single-band reference scans use a longer TR
                SeqUT.DisableTestCase(
                    lTRClockErr, RTEB_ClockCheck, "TR of single-band preparation scans does not match protocol value");
            }

            if (lKernelMode == KERNEL_PREPARE)
            {
                // Start TR clock with the preparation scan preceding the multi-band imaging scans
                mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ClockInitTR, 27, -1, m_asSLC[lSlice].getSliceIndex(), 0, 0);
            }
        }
    }
#endif


    //---------------------------------------------------------------------------
    // execute EPI Kernel
    //---------------------------------------------------------------------------

#ifdef EP2D_SE_MRE
    // Set private members for MRE in the Kernel and RefocRTEB
    const long lLcFree = m_mySeqLoop.getlFreeLoopCounter();
    const long lLcMEG  = lLcFree % (NMAXMEGS * MEG_DEFAULT_TRIGGER_STEPS);
    const long lLcPol  = lLcMEG % NMAXMEGS;

    m_EPIKernel.setlMEGPhaseOffsetCounter(lLcMEG / NMAXMEGS);
    m_pSBBRefocSE->setMEGInstanceCounter(lLcPol);
#endif

#ifdef COMPILE_EP2D_DIFF
    if (m_EPIKernel.getbCompensationEnable() && m_EPIKernel.getPointerCompGrad()->isPrepared())
    {
        m_EPIKernel.setEPIReadOutAppliesTRFill(false);
        m_EPIKernel.getPointerCompGrad()->setAdditionalFillTime(m_mySeqLoop.getlTRFill());
    }
#endif

    if (!m_EPIKernel.run(rMrProt, rSeqLim, rSeqExpo, &m_asSLC[lSlice]))
    {
        SEQ_TRACE_ERROR.print("m_EPIKernel.run failed: 0x%lx", m_EPIKernel.getNLSStatus());
        return m_EPIKernel.getNLSStatus();
    }

#ifdef COMPILE_EP2D_DIFF
    if (m_EPIKernel.getbCompensationEnable() && m_EPIKernel.getPointerCompGrad()->isPrepared())
    {
        if (!updateCompGrad(rSeqLim, rMrProt, rSeqExpo))
        {
            SEQ_TRACE_ERROR.print("updateSpoilNull failed");
            return (m_EPIKernel.getPointerCompGrad()->getNLSStatus());
        }

        if (!m_EPIKernel.getPointerCompGrad()->run(rMrProt, rSeqLim, rSeqExpo, &m_asSLC[lSlice]))
        {
            SEQ_TRACE_ERROR.print("m_EPIKernel.getPointerCompGrad()->run(...) failed");
            return (m_EPIKernel.getPointerCompGrad()->getNLSStatus());
        }
    }
#endif
    //---------------------------------------------------------------------------
    // Re-enable test cases
    //---------------------------------------------------------------------------
    if (SeqUT.isUnitTestActive())
    {
        SeqUT.EnableTestCase(lTRClockErr, RTEB_ClockCheck);
    }

    // Duration of explicitly applied cool pause
    long lCoolPauseExplicit = m_lCoolPauseTotal - m_lCoolPauseImplicit;

#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
    //---------------------------------------------------------------------------
    // execute extra fat sat and spoiler within explicit cool pause
    //---------------------------------------------------------------------------
    if (m_bApplyExtraCSat)
    {
        lStatus = m_CSatFat.run(rMrProt, rSeqLim, rSeqExpo, &m_asSLC[lSlice]);
        if (!lStatus)
        {
            SEQ_TRACE_ERROR.print("m_CSatFat.run failed: 0x%lx", lStatus);
            return lStatus;
        }
        lStatus = m_SpoilGrad.run(rMrProt, rSeqLim, rSeqExpo, &m_asSLC[lSlice]);
        if (!lStatus)
        {
            SEQ_TRACE_ERROR.print("m_SpoilGrad.run failed: 0x%lx", lStatus);
            return lStatus;
        }

        // Reduce explicit cool pause duration (applied below) by event durations
        // Note: m_bApplyExtraCSat is set to 'true' only if lCoolPauseExplicit is long enough
        //       to include both events.
        lCoolPauseExplicit -= m_CSatFat.getDurationPerRequest() + m_SpoilGrad.getDurationPerRequest();
    }
#endif

    if (isCoolTimeExecutedWithIR(rMrProt))
    {
        //---------------------------------------------------------------------------
        // execute mandatory fill time
        //---------------------------------------------------------------------------
        lStatus = fSBBCoolTimeRun(rMrProt, lCoolPauseExplicit);
        if (!lStatus)
        {
            SEQ_TRACE_ERROR.print("fSBBFillTimeRun failed: 0x%lx", lStatus);
            return lStatus;
        }
    } // Either no IR or standard IR scheme

    //------------------------------------------------------------------------------------------
    // Synchronize frequency feedback and sequence here; i.e. at the end of each volume
    //------------------------------------------------------------------------------------------
    // Synchronization events are incorporated after the EPI kernel of the last slice
    // has been played out. By doing so, interference with SeqLoop events (e.g. inversion
    // pulses / TI fill time) is avoided.
#ifdef EPI_SUPPORT_FREQ_FEEDBACK
    if (m_bB0Correction)
    {


        
        const bool isKernelRunSuitableForB0Feedback = !protFacade.isSliceAcceleration() || (protFacade.isSliceAcceleration() && m_mySeqLoop.getRunMode() == MULTI_BAND);        
        if (isKernelRunSuitableForB0Feedback)  
        {
            const auto currentVolume = getCurrentVolumeForB0Correction();

            // Feedback synchronization and bookkeeping; no wait-for-wakeup when in real-time mode
            if (!m_sFreqFeedback.SyncFeedback(rMrProt, rSeqLim, rSeqExpo, currentVolume, m_mySeqLoop.getInnerSliceCounter(), getRTController().isRealtimeProcessingEnabled() ? false : true))
            {
                SEQ_TRACE_ERROR.print("m_sFreqFeedback.SyncFeedback failed.");
                return MRI_SEQ_SEQU_ERROR;
            }
        }
    }

#endif

    //---------------------------------------------------------------------------
    // increase prep-scan-counter for current slice
    //---------------------------------------------------------------------------
#ifdef SUPPORT_iPAT_a_ep2d
    {
        if (lKernelMode == KERNEL_PREPARE)
        {
            m_alPrepScanCounter[lSlice]++;

            if (m_alPrepScanCounter[lSlice] == m_mySeqLoop.getlPreparingScans())
            {
                m_alPrepScanCounter[lSlice] = 0;
            }
        }
    }
#endif

    //-------------------------------------------------------------------------------------------
    // disable read-outs again (may have been switched on for PAT reference or adjustment scans)
    //-------------------------------------------------------------------------------------------
    disableReadoutForPrepScans(lKernelMode);


#ifdef ASL
    // disable seqloop TR fill to shift slices together, add lTRfillPASL after last slice
    if (lSlice >= (rMrProt.sliceSeries().getlSize() - 1))
    {
        lStatus = fSBBFillTimeRun(lTRfillASL * rMrProt.sliceSeries().getlSize());
        if (lStatus != MRI_SEQ_SEQU_NORMAL)
            return (lStatus);
    }
#endif

    //---------------------------------------------------------------------------
    // finished
    //---------------------------------------------------------------------------
    return (lStatus);
}

#ifdef WIN32
// ------------------------------------------------------------------------------
// Name      : Ep2d::convProt
// ------------------------------------------------------------------------------
//
// Description : try to convert protocols from previous software versions
//
//
// Return      : MRI_SEQ_SEQU_NORMAL for success
//               MRI_SEQ_SEQU_ERROR  for error
//
// ------------------------------------------------------------------------------
NLSStatus Ep2d::convProt(const MrProt& rMrProtSrc, MrProt& rMrProtDst)
{
    const NLS_STATUS lRet = MRI_SEQ_SEQU_NORMAL;

    // charm MR_00348199: actively set coilCombineMode for protocol conversion
    //      automatic protocol conversion to COILCOMBINE_ADAPTIVE_COMBINE is here disabled
    // convert to new parameters if protocol version is earlier than first VB15
    if (rMrProtSrc.getConvFromVersion() < 21510000)
    {
        rMrProtDst.coilCombineMode(SEQ::COILCOMBINE_SUM_OF_SQUARES);
    }

    return lRet;
}
#endif // #ifdef WIN32

#if (defined EPI_SUPPORT_FREQ_FEEDBACK) || (defined SUPPORT_PACE)

#if defined SUPPORT_PACE
NLSStatus Ep2d::receive(SeqLim& rSeqLim, SeqExpo& rSeqExpo, const SEQData& rSEQData)
#elif defined EPI_SUPPORT_FREQ_FEEDBACK
NLSStatus Ep2d::receive(SeqLim& /* rSeqLim */, SeqExpo& /* rSeqExpo */, const SEQData& rSEQData)
#endif
{
    //---------------------------------------------------------------------------
    // Reserve semaphore
    //---------------------------------------------------------------------------
    m_sFreqFeedback.GetFeedbackSemaphore()->acquire(30);

    MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback).print("called with ID '%s.'", rSEQData.getID());

    //---------------------------------------------------------------------------
    // No action if feedback data has inappropriate ID
    //---------------------------------------------------------------------------
    if (strcmp(rSEQData.getID(), "MRIR:DORK") != 0)
    {
        // Release semaphore
        m_sFreqFeedback.GetFeedbackSemaphore()->release();

#ifdef SUPPORT_PACE
        // Maybe someone else wants this feedback?
        if (!m_mySeqLoop.receive(rSeqLim, rSeqExpo, rSEQData))
        {
            SEQ_TRACE_ERROR.print("Error calling m_mySeqLoop.receive()");
            return m_mySeqLoop.getNLSStatus();
        }
#else
        SEQ_TRACE_ALWAYS.print("Feedback data has wrong ID (%s)", rSEQData.getID());
#endif
        return MRI_SEQ_SEQU_NORMAL;
    }

    //--------------------------------------
    // Read feedback and copy data to buffer
    //--------------------------------------

    // Interface to data stored by ICE (see EPIPhaseCorrPEFunctor.cpp)
    EPIDorkFBData* pFreqFBData = (EPIDorkFBData*)rSEQData.getData();

    // Volume counter
    const long lCurrFBVolumeNo = pFreqFBData->lRepetition;

    // Frequency offset
    const double dCurrFBFreq = pFreqFBData->dDorkRelFreq;

    //----------------------
    // Send data to FB class
    //----------------------
    // Update frequency update only if required
    if (m_bB0Correction)
    {
        // Note: EPIPhaseCorrPEFunctor starts counting with volume = 1,
        //       here we start with 0
        m_sFreqFeedback.StoreCurrFBData(lCurrFBVolumeNo - 1, dCurrFBFreq);
    }

    //------------------
    // Release semaphore
    //------------------
    m_sFreqFeedback.GetFeedbackSemaphore()->release();

    // Dump feedback data - always out of semaphore area
    MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback).print("feedback received. lCurrFBVolumeNo=%li. m_bB0Correction=%i", lCurrFBVolumeNo, m_bB0Correction);

    //---------------------------------------------------------------------------
    // Finished
    //---------------------------------------------------------------------------
    return MRI_SEQ_SEQU_NORMAL;
}

NLSStatus Ep2d::cancel()
{
    // Release semaphore
    m_sFreqFeedback.GetFeedbackSemaphore()->release();

    if (!m_mySeqLoop.cancel())
    {
        return m_mySeqLoop.getNLSStatus();
    }
    return MRI_SEQ_SEQU_NORMAL;
}
#endif // #if (defined EPI_SUPPORT_FREQ_FEEDBACK) || (defined SUPPORT_PACE)

#if (defined COMPILE_EP2D_DIFF) || (defined SUPPORT_PACE)
NLSStatus Ep2d::resume(const SEQData& rSEQData)
{
    if (!m_mySeqLoop.resume(rSEQData))
    {
        SEQ_TRACE_ERROR.print("Ep2d::resume(%s) failed.", rSEQData.getID());
        return m_mySeqLoop.getNLSStatus();
    }
    SEQ_TRACE_ALWAYS.print("*** Ep2d::resume(%s) ***", rSEQData.getID());
    return MRI_SEQ_SEQU_NORMAL;
}
#endif // #if (defined COMPILE_EP2D_DIFF) || (defined SUPPORT_PACE)

#if defined ZOOM_2DRF
//  --------------------------------------------------------------------------
//
//  Name        :  Ep2d::calculatePTX
//
//  Description :
///     \brief     Calculates the sequence's pTX pulses.
/// The sequence limits must NOT be modified by the sequence.
///
/// This method is optional, it should be implemented by
/// sequences which use pTX RF pulses.
//
//  Return      :  NLS status
//
//  --------------------------------------------------------------------------
NLSStatus Ep2d::calculatePTX(MrProt& rMrProt, const SeqLim& rSeqLim)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    // timing start
    const int64_t uiStartStamp = MrTimeStamp::getTimeus();
    const clock_t uiStartClock = clock();

    // determine some parameters from the *first* PTXRFPulse (index 0)
    const size_t                      uiCurVol        = 0;
    const long                        lPTXRFPulseSize = rMrProt.getsTXSPEC().getaPTXRFPulse().size();
    MrProtocolData::PTXTrajectoryType eExcType        = MrProtocolData::PTXTrajectoryType_EPI_1D;
    if (lPTXRFPulseSize > 0)
    {
        eExcType = rMrProt.getsTXSPEC().getaPTXRFPulse()[uiCurVol].getlTrajectoryType();
    }

    if ((rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED) && ((eExcType == MrProtocolData::PTXTrajectoryType_EPI_1D) || (eExcType == MrProtocolData::PTXTrajectoryType_ExternalIni)))
    {
        // prepZOOMitExcitation has to be called to initialize the SBB and receive a valid
        // pointer from the kernel
        // the kernel prepare is called after Ep2d::calculatePTX()
        m_EPIKernel.prepZOOMitExcitation(rMrProt, rSeqLim);
        // prepare is called after Ep2d::calculatePTX()
        SBB2DPtx* pSBBExcite = dynamic_cast<SBB2DPtx*>(m_EPIKernel.getExcitationPointer());

        if (pSBBExcite == nullptr)
        {
            return MRI_SEQ_SEQU_NORMAL; // for classes other than SBB2DPtx(): skip calculatePTX()
        }

        m_myCalcLimits.resetAllLimits();
        if (!pSBBExcite->setPointerToCalculationLimits(&m_myCalcLimits))
        {
            // caused by programming error
            // => trace also if rSeqLim.isContextPrepForBinarySearch()
            SEQ_TRACE_ERROR.print("pSBBExcite->setPointerToCalculationLimits(&m_myCalcLimits) failed.");
            lStatus = MRI_SEQ_SEQU_ERROR;
            return (lStatus);
        }

        const bool bOK = pSBBExcite->calculatePTX(rMrProt, rSeqLim);
        if (!bOK)
        {
            lStatus = pSBBExcite->getNLSStatus();
            if (MrSucceeded(lStatus))
            {
                lStatus = MRI_SEQ_SEQU_ERROR;
            }

            SEQ_TRACE_ERROR.print("calculatePTX failed: " MEAS_FMT_HEXINT64, lStatus);

            return (lStatus);
        }
    }

    // timing stop
    const int64_t uiStopStamp = MrTimeStamp::getTimeus();
    const clock_t uiStopClock = clock();

    const double dCalcTime  = static_cast<double>(uiStopStamp - uiStartStamp) * 1e-6;                                // [s]
    const double dClockTime = static_cast<double>(uiStopClock - uiStartClock) / static_cast<double>(CLOCKS_PER_SEC); // [s]

    SEQ_TRACE_ALWAYS.print("ZOOM_2DRF Info: used time (clock time) = %.3f (%.3f) s.", dCalcTime, dClockTime);

    return (lStatus);
}
#endif // ZOOM_2DRF

NLS_STATUS Ep2d::createUI(SeqLim&)
{
    //  ----------------------------------------------------------------------
    //  Delete existing instance if necessary
    //  ----------------------------------------------------------------------
    if (m_pUI)
    {
        delete m_pUI;
        m_pUI = nullptr;
    }

    //  ----------------------------------------------------------------------
    //  Variant specific instantiation of the UI class
    //  ----------------------------------------------------------------------

    try
    {
        m_pUI = new EpUI();
    }

    catch (...)
    {
        SEQ_TRACE_ERROR.print("Cannot instantiate UI class !");
        return (MRI_SEQ_SEQU_ERROR);
    }

    return (MRI_SEQ_SEQU_NORMAL);
} // end: Ep2d::createUI

EpUI* Ep2d::getUI() const
{
    return (m_pUI);
}

void Ep2d::loadSliceAdjData(SeqLim& rSeqLim, MrProt& rMrProt)
{
    // functionality only in derived Ep2d_diff sequence
}

void Ep2d::dumpSliceAdjData(SeqLim& rSeqLim, MrProt& rMrProt)
{
    // functionality only in derived Ep2d_diff sequence
}

#ifdef WIN32
bool Ep2d::exportDiffusionTimingToUI(SeqLim& rSeqLim, MrProt& rMrProt)
{
    // functionality only in derived Ep2d_diff sequence
    return true;
}
#endif

void Ep2d::setDICOMAcquisitionContrast(SeqExpo& rSeqExpo) const
{
#if (defined PERF || defined ASL)
    rSeqExpo.setDICOMAcquisitionContrast("PERFUSION"); // Acquisition Contrast (0008,9209)
#endif

#ifdef BOLD
    rSeqExpo.setDICOMAcquisitionContrast("UNKNOWN"); // Acquisition Contrast (0008,9209)
    rSeqExpo.setDICOMImageFlavor("FMRI");            // Value 3 of Image Type (0008,0008) and Frame Type (0008,9007)
#endif
}

void Ep2d::setVariantSpecificHardLimits(SeqLim& rSeqLim)
{

#ifdef EP2D_MS
    // Parallel imaging
    rSeqLim.setPATMode(SEQ::PAT_MODE_NONE, SEQ::PAT_MODE_GRAPPA, SEQ::PAT_MODE_SLICE_ACCELERATION);
    rSeqLim.setSliceAccMultiBandFactor(1, 4, 1, 1);
    rSeqLim.setSliceAccFOVShiftFactor(1, 4, 1, 1);
    rSeqLim.setRefScanMode(SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST, SEQ::PAT_REF_SCAN_EXTRA_EPI, SEQ::PAT_REF_SCAN_EXTRA);

    // Slice Adjust
    rSeqLim.setAdjustmentMode(AdjustmentMode_Standard, AdjustmentMode_SliceBySlice, AdjustmentMode_FastView);
    rSeqLim.setAdjSliceBySliceFirstOrderShim(SEQ::OFF, SEQ::ON);
    rSeqLim.setAdjSliceBySliceFrequency(SEQ::OFF, SEQ::ON);
    rSeqLim.setAdjSliceBySliceTxRef(SEQ::OFF, SEQ::ON);
    if (SysProperties::isPTxSystem())
        rSeqLim.setAdjSliceBySlicePtx(SEQ::OFF, SEQ::ON);
    else
        rSeqLim.setAdjSliceBySlicePtx(SEQ::OFF);

    // Segmentation
    rSeqLim.setSegments(1, 32, 1, 2);
    rSeqLim.setBaseResolution(64, 512, SEQ::INC_16, 128);
    rSeqLim.setPELines(32, 512, SEQ::INC_SEGMENTED_IPAT, 128);
    // Just show actual EPI-factor
    rSeqLim.getEPIFactor().setDisplayMode(SEQ::DM_SHOW);
    // Multi-shot might benefit from lower bandwidth
    rSeqLim.setBandWidthPerPixel(0, 250, 10000, 2, 750);
    
    // Multiple contrasts
    rSeqLim.setContrasts(1, 16, 1, 1);
    for (long lI = rSeqLim.getContrasts().getMin(); lI <= rSeqLim.getContrasts().getMax();
         lI += rSeqLim.getContrasts().getInc())
    {
        rSeqLim.setTE(lI - 1, 1000, 400000, 100, 1000);
    }

    rSeqLim.setRFSpoiling(SEQ::OFF);
    rSeqLim.getRFSpoiling().setDisplayMode(SEQ::DM_OFF);

#else
    // Must only be set for non-EP2D_MS flavored sequences
    rSeqLim.setEPIFactor(1, 512, SEQ::INC_SINGLESHOT, 128);
#endif

    //set up SMS-related limits
#if defined BOLD && !defined PACE3D
    rSeqLim.setPATMode(SEQ::PAT_MODE_NONE, SEQ::PAT_MODE_GRAPPA, SEQ::PAT_MODE_SENSE, SEQ::PAT_MODE_SLICE_ACCELERATION);
    rSeqLim.setSliceAccMultiBandFactor(1, 8, 1, 1);
    rSeqLim.setSliceAccFOVShiftFactor(1, 4, 1, 1);

    rSeqLim.setRefScanMode(SEQ::PAT_REF_SCAN_EXTRA_EPI, SEQ::PAT_REF_SCAN_EXTRA);
#endif

    //-------------------------------------------------------------------------------------
#ifdef COMPILE_EP2D_FID
    //-------------------------------------------------------------------------------------
    rSeqLim.setTD(0, 0, 10000000, 100, 0);
    rSeqLim.setPhasePartialFourierFactor(SEQ::PF_OFF, SEQ::PF_7_8, SEQ::PF_6_8);
    rSeqLim.setInversion(SEQ::INVERSION_OFF);
    if (SysProperties::isUHFSystem())
        rSeqLim.setRFPulseType(SEQ::RF_FAST, SEQ::RF_NORMAL, SEQ::RF_LOW_SAR);
    else
        rSeqLim.setRFPulseType(SEQ::RF_NORMAL);
    rSeqLim.setRFSpoiling(SEQ::OFF);
    rSeqLim.getRFSpoiling().setDisplayMode(SEQ::DM_OFF);

    // default number of repetitions must be consistent with EVA protocol (CHARM 339267)
    rSeqLim.setRepetitions(0, 4095, 1, 14);

#ifdef BOLD
    if (SysProperties::isUHFSystem())
        rSeqLim.setRFPulseType(SEQ::RF_FAST, SEQ::RF_NORMAL, SEQ::RF_LOW_SAR);
    else
        rSeqLim.setRFPulseType(SEQ::RF_NORMAL, SEQ::RF_LOW_SAR);
    rSeqLim.set2DInterpolation(SEQ::OFF);
    rSeqLim.setRepetitions(0, 4095, 1, 19);

    // mosaic image format does not support multi-slice-multi-angle
    rSeqLim.disableMSMA();

    // long-term averaging is not supported by IceProgramOnline2D (CHARM 308391)
    rSeqLim.setAverages(1, 1, 1, 1);

    // set flag to exclude SET and REP dimensions from data amount calculation in Sequence.cpp this is appropriate
    // for ep2d_bold and ep2d_pace sequences, which perform inline image calculation without raw data storage
    rSeqLim.setInteractiveRealtime(SEQ::ON);

    rSeqLim.setMultipleSeriesMode(SEQ::MULTIPLE_SERIES_OFF);
    rSeqLim.getMultipleSeriesMode().setDisplayMode(SEQ::DM_OFF);

    // Allow asymmetric sat regions
    rSeqLim.setRSatShapeMode(SEQ::SAT_SHAPE_STANDARD, SEQ::SAT_SHAPE_ASYMMETRIC);

    // file containing the default postprocessing protocol (EVAProtocol)
#ifdef WIN32
    rSeqLim.setDefaultEVAProt(_T("%SiemensEvaDefProt%\\BOLD\\t-test_10B10A_moco.evp"));
#endif

    // sequence ep2d_pace can only be run with fatsat on (CHARM 309540)
#ifdef PACE3D

    rSeqLim.setFatWaterContrast(FatWaterContrast_FatSaturation);

#endif // PACE3D

#endif // BOLD

#ifdef ASL
    rSeqLim.setBaseResolution(64, 256, SEQ::INC_NORMAL, 64); // enabled for ASL
    rSeqLim.setPELines(32, 256, 1, 64);                      // enabled for ASL
    rSeqLim.set2DInterpolation(SEQ::OFF);
    rSeqLim.setSlices(2, 64, 1, 2);         // minimum of 2
    rSeqLim.setRepetitions(2, 4095, 2, 20); // should be only odd !!
    rSeqLim.setSliceSeriesMode(SEQ::ASCENDING, SEQ::DESCENDING, SEQ::INTERLEAVED);
    rSeqLim.setTR(0, 10000, 30000000, 1000, 5000000); // pushed up for asl
    rSeqLim.setIntro(SEQ::ON, SEQ::OFF);              // for pace
    rSeqLim.setDelayTimeInTR(0, 0, 1000, 0);          // disabled

    // mosaic image format does not support multi-slice-multi-angle
    rSeqLim.disableMSMA();

    // long-term averaging is not supported by IceProgramOnline2D (CHARM 308391)
    rSeqLim.setAverages(1, 1, 1, 1);

    // sequence ep2d_pace can only be run with fatsat on (CHARM 309540)
    rSeqLim.setFatWaterContrast(FatWaterContrast_FatSaturation);

    // allow two modes of fat saturation, weak - only one fat sat before all slices, strong - fat sat every slice
    rSeqLim.setFatSatMode(SEQ::FAT_SAT_WEAK, SEQ::FAT_SAT_STRONG); // (CHARM 366331)

    // set Rsat for display of labelling region
    rSeqLim.setRSats(0, 2, 1, 0);
    rSeqLim.setRSatThickness(1.000, 200.000, 1.000, 100.000); // enabled for ASL
    rSeqLim.setPSatMode(SEQ::PSAT_SINGLE_REG, SEQ::PSAT_DOUBLE_REG);
    rSeqLim.setPSatThickness(1.000, 200.000, 1.000, 100.000); // enabled for ASL
    rSeqLim.setPSatGapToSlice(0.000, 100.000, 0.100, 25.000); // enabled for ASL
    rSeqLim.setMTC(SEQ::OFF);

    // KAH Incorporating VD11 changes to ASL parameter card - commented options are for 3D
    rSeqLim.setAslBolusDuration(0, 2000000, 25000, 700000);
    rSeqLim.setAslInversionTime(0, 0, 4000000, 100, 1800000); // Array!!

    // PCASL settings
    rSeqLim.setAslLabelingDuration(100000, 4000000, 10000, 1800000);
    rSeqLim.setAslDelayArraySize(1, MAX_ASL_INFLOW_PHASES, 1, 1);
    for (int32_t i = 0; i < MAX_ASL_INFLOW_PHASES; i++)
        rSeqLim.setAslPostLabelingDelay(i, 100000, 5000000, 10000, 1800000); // array

    rSeqLim.setAslMode(SEQ::ASL_PICOREQ2TIPS, SEQ::ASL_PSEUDOCASL);
    rSeqLim.setAslFlowLimit(0.1, CRUSHGRAD_MAX_VELOCITY_ENC, 0.1, CRUSHGRAD_MAX_VELOCITY_ENC);
    rSeqLim.setAslInversionArraySize(1, 1, 1, 1);

    rSeqLim.setInversion(SEQ::INVERSION_OFF);

    rSeqLim.getMTC().setDisplayMode(SEQ::DM_OFF);

    // file containing the default postprocessing protocol (EVAProtocol)
#ifdef WIN32
    rSeqLim.setDefaultEVAProt(_T("%SiemensEvaDefProt%\\ASL\\ASL_moco.evp"));
#endif
#endif // ASL

#ifdef PERF
    // file containing the default postprocessing protocol (EVAProtocol)
#ifdef WIN32
    rSeqLim.setDefaultEVAProt(_T("%SiemensEvaDefProt%\\Perfusion\\LOCALAIF_DEFAULT.evp"));
#endif

    // perfusion and motion correction calculations in ICE do not support multi-slice-multi-angle
    rSeqLim.disableMSMA();

#endif // PERF

    //-------------------------------------------------------------------------------------
#endif // COMPILE_EP2D_FID

    const float fExtendedMaxFoV = SysProperties::getExtendedFoVMax();
    const float fMaxFoV         = SysProperties::getFoVMax();
    rSeqLim.setPhaseFOV(25.0, fExtendedMaxFoV, 1.0, fMaxFoV);
    rSeqLim.setReadoutFOV(25.0, fExtendedMaxFoV, 1.0, fMaxFoV);

    //-------------------------------------------------------------------------------------
#ifdef COMPILE_EP2D_SE
    //-------------------------------------------------------------------------------------
    rSeqLim.setTD(0, 0, 10000000, 100, 0);
    rSeqLim.setInversion(SEQ::INVERSION_OFF, SEQ::SLICE_SELECTIVE);
    rSeqLim.setTI(0, 10000, 5000000, 100, 2500000);
    rSeqLim.setRFPulseType(SEQ::RF_NORMAL, SEQ::RF_LOW_SAR);
    rSeqLim.setFlipAngle(90.000, 90.000, 5.000, 90.000);
    rSeqLim.getFlipAngle().setDisplayMode(SEQ::DM_OFF);

    // remove option for partial Fourier 4/8 due to poor image quality (CHARM 310630)
    rSeqLim.setPhasePartialFourierFactor(SEQ::PF_OFF, SEQ::PF_5_8, SEQ::PF_6_8, SEQ::PF_7_8);
#ifdef EP2D_SE_MRE
    rSeqLim.setFatWaterContrast(
        FatWaterContrast_Spair,
        FatWaterContrast_FatSaturation,
        FatWaterContrast_WaterExcitation,
        FatWaterContrast_Standard);
    rSeqLim.setFatSatMode(SEQ::FAT_SAT_STRONG, SEQ::FAT_SAT_WEAK);
    rSeqLim.setCoilCombineMode(SEQ::COILCOMBINE_ADAPTIVE_COMBINE); // COILCOMBINE_SUM_OF_SQUARES is not compatible with
                                                                   // current complex-valued solution for defect 324693
#else
    rSeqLim.setFatWaterContrast(
        FatWaterContrast_FatSaturation,
        FatWaterContrast_WaterExcitation,
        FatWaterContrast_Spair,
        FatWaterContrast_Standard);
    // Allow two modes of fat saturation
    //  weak:   standard
    //  strong: additional gradient reversal (excitation vs. refocusing slice selection gradient)
    rSeqLim.setFatSatMode(SEQ::FAT_SAT_WEAK, SEQ::FAT_SAT_STRONG);
#endif

#ifdef EP2D_SE_MRE
    rSeqLim.setTE(0, 1000, 200000, 100, 200000);
    // Beyond 1000ms, the TR limits increments progress x10... the "500" below
    //   translate into 50ms increments for values > 1000ms. This must be an overloaded getLimits
    //   handler.
    rSeqLim.setTR(0, 500000, 10000000, 500, 1500000);
    rSeqLim.setConcatenations(1, 1, 1, 1);

    // Restrict number of RSATs so that 50ms trigger spacing is kept
    rSeqLim.setRSats(0, 2, 1, 0);

    // Hide incompatible settings
    rSeqLim.getMTC().setDisplayMode(SEQ::DM_OFF);
    rSeqLim.getInversion().setDisplayMode(SEQ::DM_OFF);

    // Resolution limits for MRE
    // - original: rSeqLim.setBaseResolution           (64, 512, SEQ::INC_NORMAL, 128                  );
    // - default 100 led to some strange effects in the default protocol -> stick with the regular 128
    rSeqLim.setBaseResolution(64, 256, SEQ::INC_NORMAL, 128);

    // Use (hidden) flexible resolution to make MRE images compatible with inline inversion
    // - this also currently requires a restriction to square FoVs
    // - Interpolation checkbox is currently hidden
    rSeqLim.get2DInterpolation().setDisplayMode(SEQ::DM_OFF);
    rSeqLim.setSquareFOVOnly(true);

    // No triggering yet
    rSeqLim.getPhysioModes().unset(SEQ::SIGNAL_ALL, SEQ::METHOD_TRIGGERING);
    rSeqLim.getPhysioModes().setDisplayMode(SEQ::DM_OFF);

    // Necessary for reconstruction of PC Angio in the right manner for MRE
    rSeqLim.setReconstructionMode(SEQ::RECONMODE_MAGNITUDE);
    rSeqLim.getReconstructionMode().setDisplayMode(SEQ::DM_OFF);
    rSeqLim.setPhaseImages(SEQ::YES, SEQ::NO);
    rSeqLim.setMagnitudeImages(SEQ::OFF, SEQ::ON);
    rSeqLim.setRephasedImage(SEQ::ON, SEQ::OFF);
    rSeqLim.setMagnitudeSum(SEQ::OFF, SEQ::ON);

    rSeqLim.setAverages(1, 1, 1, 1);
    rSeqLim.setRepetitions(0, 0, 0, 0);
    rSeqLim.getRepetitions().setDisplayMode(SEQ::DM_OFF);
    rSeqLim.getDelayTimeInTR().setDisplayMode(SEQ::DM_OFF);
    rSeqLim.getMultipleSeriesMode().setDisplayMode(SEQ::DM_OFF);

    // It is used to transport the info about “I am an MRE sequence” to the AddIn
    rSeqLim.setApplicationDetails(SEQ::APPLICATION_MRE);
#endif

    // file containing the default postprocessing protocol (EVAProtocol)
#ifdef WIN32
    rSeqLim.setDefaultEVAProt(_T("%SiemensEvaDefProt%\\Inline\\Inline.evp"));
#endif
    //-------------------------------------------------------------------------------------
#endif // COMPILE_EP2D_SE

    // Limits for multi-shot configurations
#ifdef EP2D_MS

    // Inversion
    rSeqLim.setInversion(SEQ::INVERSION_OFF, SEQ::SLICE_SELECTIVE);
    rSeqLim.setFreezeSuppressedTissue(SEQ::OFF, SEQ::ON);
    rSeqLim.setTI(0, 10000, 5000000, 100, 2500000);
    rSeqLim.setIRScheme(SEQ::IR_SCHEME_AUTO);

    rSeqLim.setPhasePartialFourierFactor(SEQ::PF_OFF, SEQ::PF_7_8, SEQ::PF_6_8, SEQ::PF_5_8);
    rSeqLim.setPhaseCorrectionMode(MrProtocolData::PHASECORR_INTERNAL, MrProtocolData::PHASECORR_EXTERNAL);
    if (SysProperties::getNominalB0() > 2.5)
    {
        rSeqLim.setFatSupOpt(MrProtocolData::FATSUPOPT_DEFAULT, MrProtocolData::FATSUPOPT_BRAIN);
    }

    rSeqLim.setCoilCombineMode(SEQ::COILCOMBINE_SUM_OF_SQUARES, SEQ::COILCOMBINE_ADAPTIVE_COMBINE);

    // Allow automatic optimization of TE
    rSeqLim.setTOM(SEQ::TOM_OFF, SEQ::TOM_MINIMIZE_TE);

    // Allow image filter
    rSeqLim.setFilterType(
        SEQ::FILTER_NONE,
        SEQ::FILTER_RAW,
        SEQ::LARGE_FOV,
        SEQ::ELLIPTICAL,
        SEQ::HAMMING,
        SEQ::NORMALIZE,
        SEQ::PRESCAN_NORMALIZE,
        SEQ::FILTER_IMAGE);

    // Allow disabling of enforced frequency adjustment
    rSeqLim.setAdjFreProtRelated(SEQ::ON, SEQ::OFF);

    // Support POCS
    rSeqLim.setPOCS(SEQ::POCS_OFF, SEQ::POCS_READ_PHASE);

    
    // SE FLAIR multishot specifics
#ifdef COMPILE_EP2D_SE
    rSeqLim.setMTC(SEQ::OFF, SEQ::ON);
    rSeqLim.setMTCMode(
        MrProtocolData::MTC_MODE_OFF,
        MrProtocolData::MTC_MODE_WEAK,
        MrProtocolData::MTC_MODE_MEDIUM,
        MrProtocolData::MTC_MODE_STRONG);
    rSeqLim.setFlowAttenuation(
        MrProtocolData::FLOW_ATTENUATION_OFF,
        MrProtocolData::FLOW_ATTENUATION_WEAK,
        MrProtocolData::FLOW_ATTENUATION_STRONG);

    rSeqLim.setFlipAngle(130.000, 180.000, 1.000, 180.000);

#endif

#ifdef COMPILE_EP2D_FID
    rSeqLim.setRepetitions(0, 4095, 1, 0);
#endif
#endif

}

NLSStatus Ep2d::setVariantSpecificExports(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo)
{
    MrProtFacade protFacade(rMrProt);

#ifdef COMPILE_EP2D_FID
    fSUSetSequenceString(rMrProt, rSeqLim, rSeqExpo, "epfid");

#ifdef PERF

#ifdef WIN32
    rSeqExpo.setApplicationCard(SEQ::APPLICATION_CARD_PERF);
    rSeqExpo.setApplicationCardName(SEQ::APPLICATION_CARD_NAME_PERF);
#endif // WIN32

    //
    //  Due to the new concept of EVA-protocols the perfusion-data was removed from the MrProt.
    //  Now we can not choose between IceProgramOnline2D and IceProgramStandard any more.
    //  So we are forced to use the IceProgramStandard in any case and have to live with the
    //  following restriction for the moment (or forever?, CHARM 305697):
    //
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  SEQUENCE ep2d_fid CAN ONLY BE EXECUTED FOR IMAGING WITHOUT PERF-POSTPROCESSING WITH
    //  OFFLINE RECONSTRUCTION.
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    //  This is due to the fact that IceProgramStandard does not support online-reconstruction,
    //  but is required to do perfusion postprocessing.
    //
    rSeqExpo.setOnlineFFT(SEQ::ONLINE_FFT_NONE);
    rSeqExpo.setICEProgramFilename("%SiemensIceProgs%\\IceProgramStandard");
    // Insert additional SMS ICE program
    // 	if (protFacade.isSliceAcceleration())
    // 	{
    // 		rSeqExpo.AddAdditionalIceProgramFileName("%SiemensIceProgs%\\IceProgramSMSAcc");
    // 	}

    rSeqExpo.setICEProgramParam(ICE_PROGRAM_PARA_SHOW_OFFLINE, SEQ::SO_SHOW_YES);
    rSeqExpo.setICEProgramParam(ICE_PROGRAM_PERF_THRESHOLD, 50);
    rSeqExpo.setICEProgramParam(ICE_PROGRAM_PERF_BASE_START, 5);
    rSeqExpo.setICEProgramParam(ICE_PROGRAM_PERF_BASE_END, 9);
    rSeqExpo.setICEProgramParam(ICE_PROGRAM_PERF_TTP_SCALING, 100000);
    rSeqExpo.setICEProgramParam(ICE_PROGRAM_PERF_TTP_START, 1000);
#endif // PERF

#if defined BOLD

#ifdef WIN32
    rSeqExpo.setApplicationCard(SEQ::APPLICATION_CARD_FMRI);
    rSeqExpo.setApplicationCardName(SEQ::APPLICATION_CARD_NAME_BOLD);
#endif // WIN32

    //
    //  Due to the new concept of EVA-protocols the fMRI-data was removed from the MrProt.
    //  Now we can not choose between IceProgramOnline2D and IceProgramStandard any more.
    //  So we are forced to use the IceProgramOnline2D in any case and have to live with the
    //  following restriction for the moment (or forever?, see CHARM 305697):
    //
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  SEQUENCE ep2d_bold CAN ONLY BE EXECUTED FOR IMAGING WITHOUT fMRI-POSTPROCESSING WHEN
    //  ONLINE RECONSTRUCTION IS POSSIBLE.
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    //  This is due to the fact that IceProgramOnline2D does not support OFFLINE-reconstruction,
    //  but is required to do fMRI postprocessing.
    //
    rSeqExpo.setOnlineFFT(SEQ::ONLINE_FFT_PHASE);
    // JR 4 lines commented out
    // #ifdef PACE3D
    rSeqExpo.setICEProgramFilename("%SiemensIceProgs%\\IceProgramStandard");

    // Insert additional SMS ICE program
    if (protFacade.isSliceAcceleration())
    {
        rSeqExpo.AddAdditionalIceProgramFileName("%SiemensIceProgs%\\IceProgramSMSAcc");
    }

#ifdef WIN32
    if (rSeqExpo.getOnlineFFT() != SEQ::ONLINE_FFT_PHASE)
    {
        // m_mySeqLoop.setIceProgram says offline!
        // we can't help, if IceProgramOnline2D supports only online,
        // ergo:
        SEQ_TRACE_WARN_COND(!rSeqLim.isContextPrepForBinarySearch(), "m_mySeqLoop.setIceProgram says offline");

        if (m_mySeqLoop.getCTRL_SeqLoop_SwitchToOffline())
        {
            SEQ_TRACE_WARN.print(
                "m_mySeqLoop.setIceProgram says offline, because it is forced to do so.\nCheck registry-key: "
                "SOFTWARE\\Siemens\\Numaris4\\Config\\Modality\\Sequence\\CTRL_SeqLoop_SwitchToOffline");
        }

        return SeverePrepareErrorReturn(MRI_SEQ_SEQU_ERROR);
    }
#endif // WIN32

#endif // BOLD

#if defined ASL

#ifdef WIN32
    rSeqExpo.setApplicationCard(SEQ::APPLICATION_CARD_PERF);
    rSeqExpo.setApplicationCardName(SEQ::APPLICATION_CARD_NAME_PERF);
#endif // WIN32

    rSeqExpo.setOnlineFFT(SEQ::ONLINE_FFT_PHASE);
    rSeqExpo.setICEProgramFilename("%SiemensIceProgs%\\IceProgramStandard");

#ifdef WIN32
    if (rSeqExpo.getOnlineFFT() != SEQ::ONLINE_FFT_PHASE)
    {
        SEQ_TRACE_WARN_COND(!rSeqLim.isContextPrepForBinarySearch(), "m_mySeqLoop.setIceProgram says offline");

        if (m_mySeqLoop.getCTRL_SeqLoop_SwitchToOffline())
        {
            SEQ_TRACE_WARN.print(
                "m_mySeqLoop.setIceProgram says offline, because it is forced to do so.\nCheck registry-key: "
                "SOFTWARE\\Siemens\\Numaris4\\Config\\Modality\\Sequence\\CTRL_SeqLoop_SwitchToOffline");
        }

        return SeverePrepareErrorReturn(MRI_SEQ_SEQU_ERROR);
    }
#endif // WIN32

#endif // ASL

#endif // COMPILE_EP2D_FID

#ifdef COMPILE_EP2D_SE
    if (rMrProt.preparationPulses().getucInversion() == SEQ::INVERSION_OFF)
        fSUSetSequenceString(rMrProt, rSeqLim, rSeqExpo, "epse");
    else
        fSUSetSequenceString(rMrProt, rSeqLim, rSeqExpo, "epir");

#ifdef EP2D_SE_MRE

    fSUSetSequenceString(rMrProt, rSeqLim, rSeqExpo, "epseMRE");

    rSeqExpo.setNSet(NMAXMEGS);
    rSeqExpo.setNIda(MEG_DEFAULT_TRIGGER_STEPS);

    rSeqExpo.setICEProgramFilename("%SiemensIceProgs%\\IceProgramElastography2D");

    // PCAngio settings
    rSeqExpo.setICEProgramParam(
        ICE_PROGRAM_PARA_CTRL_MASK,
        (rSeqExpo.getICEProgramParam(ICE_PROGRAM_PARA_CTRL_MASK) | ICE_PROGRAM_PARA_IS_ANGIO));

    // Zero out the Maxwell Coefficients for now
    for (int iI = 0; iI < NMAXMEGS; iI++)
    {
        rSeqExpo.setMaxwellCoefficients(0 + (iI)*4, 0.0);
        rSeqExpo.setMaxwellCoefficients(1 + (iI)*4, 0.0);
        rSeqExpo.setMaxwellCoefficients(2 + (iI)*4, 0.0);
        rSeqExpo.setMaxwellCoefficients(3 + (iI)*4, 0.0);
    }

    if (!rSeqLim.isContextPrepForBinarySearch())
    {
        if (!IsTRAdequateForMRE(rMrProt.getalTR()[0], MREProperties::epiMRE))
        {
            SEQ_TRACE_ERROR.print("IsTRAdequateForMRE(...) failed");
            return MRI_SEQ_SEQU_ERROR;
        }
    }
#endif

    if (m_REOInfo.getPhasePartialFourierFactor() < 0.70)
    {
        rSeqExpo.setPCAlgorithm(SEQ::PC_ALGORITHM_MARGOSIAN);
    }

#ifdef WIN32
#ifdef EP2D_SE_MRE
    rSeqExpo.setApplicationCard(SEQ::APPLICATION_CARD_NONE);
#else
    rSeqExpo.setApplicationCard(SEQ::APPLICATION_CARD_INLINE);
#endif
#endif
#endif

#ifdef EP2D_MS
    if (isDLReconMode(rMrProt))
    {
        rSeqExpo.AddAdditionalIceProgramFileName("%SiemensIceProgs%\\IceProgramMsEpiAdd");
        rSeqExpo.AddAdditionalIceProgramFileName("%SiemensIceProgs%\\IceChannelCompression");
        const auto targetChannels = 12;
        rSeqExpo.setGeometricChannelCompression_TargetChannels(targetChannels);
    }
    if (protFacade.isSliceAcceleration())
    {
        rSeqExpo.AddAdditionalIceProgramFileName("%SiemensIceProgs%\\IceProgramSMSAcc");
    }
    if ((m_REOInfo.getPhasePartialFourierFactor() < 0.90) && (rMrProt.getsKSpace().getucPOCS() == SEQ::POCS_READ_PHASE))
    {
        rSeqExpo.setPCAlgorithm(SEQ::PC_ALGORITHM_POCS_PE);
    }
#endif

    return MRI_SEQ_SEQU_NORMAL;
}

void Ep2d::setInitialDummyScansBasic(MrProt& rMrProt)
{
    if ((rMrProt.getlRepetitions() > 0) || (rMrProt.getlAverages() > 1))
    {
        // we want to prepare at least for 3 sec; i.e.
        m_lInitialDummyScans = (long)(3000000.0 / (rMrProt.tr()[0] * rMrProt.physiology().phases()) + 1.0);
    }

    else
    {
        // no preparing for 'real' single shot
        m_lInitialDummyScans = 0;
    }

#ifdef EP2D_SE_MRE
    // fixed number of prep scans for MRE
    m_lInitialDummyScans = 1;
#endif
}

void Ep2d::setInitialDummyScansSMS(MrProt& rMrProt)
{
    // BOLD and ms-EPI need longer to reach steady state as TR is typically shorter in multi-band case
    // thus the value of m_lInitialPrepScans calculated above is utilized.
    // The number of initial dummy scans however is set to zero as the TR of the preparation scan
    // typically is much longer than the TR of the multi-band accelerated scan on which m_lInitialDummyScans is based
    // and would thus unnecessarily increase the measurement time
    m_lSliceAccelDummyScans = m_lInitialDummyScans;
    m_lInitialDummyScans    = 0;
}

long Ep2d::calcRequiredPrepScans(MrProt& rMrProt)
{
    long lRequiredPrepScans = 0;

    // Add up total number of prep scans
    lRequiredPrepScans = m_lInitialDummyScans + m_lPhaseCorrPrepScans;

#ifdef SUPPORT_iPAT_a_ep2d
    lRequiredPrepScans += m_lPATRefScans;
#endif

    // ---------------------------------------------------------------------------
    // For slice acceleration one prep scan for SliceAcc ACS is acquired in singleband mode
    // ---------------------------------------------------------------------------

    lRequiredPrepScans += m_lSliceAccelRefScans;
    lRequiredPrepScans += m_lSliceAccelDummyScans;
    lRequiredPrepScans += m_lSliceAccelPhaseCorrScans;

    return lRequiredPrepScans;
}

NLSStatus Ep2d::CheckAndAdaptTE(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo, long& lNeededTE)
{
    // ---------------------------------------------------------------------------
    // check TE
    // ---------------------------------------------------------------------------
    // Note: If TOM is enabled or if it is not possible to apply the requested
    // b-value with the current protocol TE value, the 'needed TE' will be updated here.
    if (m_EPIKernel.getNeededTE() != rMrProt.te()[getTEContrastIndex( rMrProt.getlContrasts() )])
    {
        const long lNeededTENotOnInc = m_EPIKernel.getNeededTE();
        double     dInc              = static_cast<double>(rSeqLim.getTE()[getTEContrastIndex( rMrProt.getlContrasts())].getInc());

        // ///////////////////////////////////////////////////////////////////////////////////////////////
        // #pragma message ("..........NOTE!!!!! we assume something about TE-limit-handling in MrUILink !")
        // ///////////////////////////////////////////////////////////////////////////////////////////////
        if (lNeededTENotOnInc > 10000)
            dInc *= 10.0;
        if (lNeededTENotOnInc > 1000000)
            dInc *= 10.0;

        lNeededTE = static_cast<long>(0.5 + dInc * ceil(static_cast<double>(lNeededTENotOnInc) / dInc));

        if (!m_EPIKernel.increaseTE(lNeededTE))
        {
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "%m_EPIKernel.increaseTE() failed with status: 0x%lx", m_EPIKernel.getNLSStatus());
            return SeverePrepareErrorReturn(m_EPIKernel.getNLSStatus());
        }

    }
    else
    {
        lNeededTE = rMrProt.te()[getTEContrastIndex( rMrProt.getlContrasts())];
    }

    
  // Check that all TE's (including those calculated automatically) stay within the boundaries
    // Note: We are introducing UI logic here. As an alternative, one might consider
    //       adding suitable getLimits-handlers for all parameters which affect TE.
    for (int32_t iI = 0; iI < rMrProt.getlContrasts(); ++iI)
    {
        if ((m_EPIKernel.getActualTE(iI) < rSeqLim.getTE()[iI].getMin())
            || (m_EPIKernel.getActualTE(iI) > rSeqLim.getTE()[iI].getMax()))
        {
            return prepareError(rSeqLim, "TE out of bounds");
        }
    }

    return MRI_SEQ_SEQU_NORMAL;
}

void Ep2d::setDiffusionAdjustmentScans(MrProt& rMrProt, long lSlice, long lKernelMode)
{
    // functionality only in derived Ep2d_diff sequence
}

bool Ep2d::setDiffusionLoopCounters(long lSlice, long lKernelMode)
{
    // functionality only in derived Ep2d_diff sequence
    return true;
}

std::string Ep2d::getSequenceVariantText() const
{
//-------------------------------------------------------------------------------------
// sequence variant name
//-------------------------------------------------------------------------------------
#if defined COMPILE_EP2D_FID
#ifdef PERF
    return {"Single-shot EPI 2D FID sequence with perfusion post-processing"};
#endif
#ifdef BOLD
    return {"Single-shot EPI 2D FID sequence with fMRI post-processing"};
#endif
#ifdef ASL
    return {"Single-shot EPI 2D FID sequence with pulsed ASL preparation"};
#endif
#ifdef EP2D_MS
    return {"Multi-shot EPI 2D FID sequence with phase-encoding segmentation"};
#endif
#elif defined COMPILE_EP2D_SE
#ifdef EP2D_MS
    return {"Multi-shot EPI 2D SE sequence with phase-encoding segmentation"};
#else
    return {"Single-shot EPI 2D SE sequence"};
#endif
#endif
    return {""};
}

NLSStatus Ep2d::configureEPIKernel(SeqLim& rSeqLim, MrProt& rMrProt)
{
#if defined EP2D_MS && defined COMPILE_EP2D_SE

    const bool isGradientRewindingAndPrephasingRequired = rMrProt.getlContrasts() > 1;

    // For multiple contrasts, rewinders are required
    m_EPIKernel.setRewindBlips(isGradientRewindingAndPrephasingRequired);
    m_EPIKernel.setRewindRO(isGradientRewindingAndPrephasingRequired);

    // Apply pre-phasing immediately before EPI readout
    // => Mandatory for multi-contrast
    // => Recommended for diffusion-weighted (single- or multi-contrast)
    // => Optional otherwise
    m_EPIKernel.setPrephaseROAfterRTEBPlugIn(isGradientRewindingAndPrephasingRequired);
    m_EPIKernel.setPrephaseBlipsAfterRTEBPlugIn(isGradientRewindingAndPrephasingRequired);
    
    // set flow attenuation in EPI kernel
    switch (rMrProt.getsPrepPulses().getlFlowAttenuation())
    {
        case MrProtocolData::FLOW_ATTENUATION_STRONG:
            m_EPIKernel.setFlowAttenuationDirection(SEQ::AXIS_SLICE);
            m_EPIKernel.setFlowAttenuationStrength(1300);
            break;
        case MrProtocolData::FLOW_ATTENUATION_WEAK:
            m_EPIKernel.setFlowAttenuationDirection(SEQ::AXIS_SLICE);
            m_EPIKernel.setFlowAttenuationStrength(500);
            break;
        case MrProtocolData::FLOW_ATTENUATION_OFF:
            m_EPIKernel.deactivateFlowAttenuation();
        default:
            break;
    }
#endif

      // Set number of contrasts before / after PlugIn
    const auto numberOfContrasts = rMrProt.getlContrasts();
    m_EPIKernel.setNumberOfContrasts(numberOfContrasts);
    m_EPIKernel.setNumberOfContrastsBeforeRTEBPlugIn(getNumberOfContrastsBeforeRTEBPlugIn(numberOfContrasts));
    m_EPIKernel.setTEContrastIndex(getTEContrastIndex(numberOfContrasts));

    if (rMrProt.fastImaging().getucFreeEchoSpacing())
    {
        m_EPIKernel.setUseFixedEchoSpacing(rMrProt.fastImaging().getlEchoSpacing());
    }
    else
    {
        m_EPIKernel.setUseShortestEchoSpacing();
    }

    if (!m_EPIKernel.setPointerToReorderInfo(&m_REOInfo))
    {
        SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "m_EPIKernel.setPointerToReorderInfo() failed with status: 0x%lx", m_EPIKernel.getNLSStatus());
        return SeverePrepareErrorReturn(m_EPIKernel.getNLSStatus());
    }

    m_EPIKernel.setGSWDGradientPerformance(rMrProt, rSeqLim);
    m_EPIKernel.setIgnoreForbiddenEchoSpacingRange(false); // important here in fSEQPrep, do not move to fSEQInit
    m_EPIKernel.setEchoTrainLength(m_REOInfo.getEchoTrainLength());
    m_EPIKernel.setCenterSegment(m_REOInfo.getKSpaceCenterSegment());

    // External phase correction scan
    if (rMrProt.preparationPulses().getlPhaseCorrectionMode() == MrProtocolData::PHASECORR_EXTERNAL)
    {

        // External phase correction requires at least one echo before and after k-space center
        if ((m_REOInfo.getSegmentsBeforeKSpaceCenter() < 1)
            || (m_REOInfo.getMeasuredSegments() - m_REOInfo.getSegmentsBeforeKSpaceCenter() < 2))
        {
            SEQ_TRACE_ERROR_COND(!rSeqLim.isContextPrepForBinarySearch(), "Insufficient number of echoes for external phase correction");
            return MRI_SEQ_SEQU_ERROR;
        }

        m_EPIKernel.setInternalPhaseCorrection(false);
        m_EPIKernel.setExternalEPIPhaseCorrection(true);
    }
    else
    {
        m_EPIKernel.setInternalPhaseCorrection(true);
        m_EPIKernel.setExternalEPIPhaseCorrection(false);
    }

#ifdef EPI_SUPPORT_FREQ_FEEDBACK
    // If B0 correction is active: enable realtime feedback of phase correction scans
    m_EPIKernel.setFlagPCforRTFeedback(m_bB0Correction);
#endif // #ifdef EPI_SUPPORT_FREQ_FEEDBACK

// unfortunately this "if not diffusion" needs to be here since the "m_EPIKernel.setWantedTE(rMrProt.te()[0]);"
// below should not apply to diffusion. it can be removed once all non-TGSE variants are refactored to derived classes.
#ifndef COMPILE_EP2D_DIFF
    m_EPIKernel.setWantedTE(rMrProt.te()[0]);
#endif // not diffusion

#ifdef EP2D_MS
    if (rMrProt.TOM() != SEQ::TOM_MINIMIZE_TE)
    {
        // Prepare kernel with the determining TE
        m_EPIKernel.setWantedTE(rMrProt.te()[getTEContrastIndex(rMrProt.getlContrasts())]);
    }
    else
    {
        // Prepare kernel with shortest possible TE
        m_EPIKernel.setWantedTE(0);
    }
#endif

#ifdef SUPPORT_iPAT_a_ep2d
    {
        if (m_REOInfo.isPATActive() && m_REOInfo.getPATAccelerationFactorPE() > 1  && !m_REOInfo.isPATGRERefScans() )
        {
            // To achieve that the segmented reference scans have exactly the same effective
            // TE as the imaging scans we would write the following line:
            //
            // m_EPIKernel.setUseEchoShifting(true, m_REOInfo.getPATAccelerationFactorPE(),
            // m_REOInfo.getPATRefCounterInSegmentWithKSCenter());
            //
            // This leads to a non-convex reference line parameter space, because the value
            // m_REOInfo.getPATRefCounterInSegmentWithKSCenter() can change depending on the
            // number of reference lines like this:
            //
            // PATRefCounterInSegmentWithKSCenter = (PATRefLinesPE/2)%PATAccelerationFactorPE
            //
            // Therefore we currently use always the minimum possible TE for the imaging scans.
            // The maximum error for the effective TE of the extra reference lines is EchoSpacing/2.
            //
            m_EPIKernel.setUseEchoShifting(true, m_REOInfo.getPATAccelerationFactorPE(), 0);
            m_EPIKernel.setExpectedMaxLinesForPEBlip(m_REOInfo.getMaxLinIncrementBetweenEchoes());
            m_EPIKernel.setExpectedMaxPartitionsFor3DBlip(m_REOInfo.getMaxParIncrementBetweenEchoes());
        }
        else
        {
            m_EPIKernel.setUseEchoShifting(false);
            if (m_REOInfo.getEchoTrainLength() > 1)
            {
                m_EPIKernel.setExpectedMaxLinesForPEBlip(m_REOInfo.getMaxLinIncrementBetweenEchoes());
                m_EPIKernel.setExpectedMaxPartitionsFor3DBlip(m_REOInfo.getMaxParIncrementBetweenEchoes());
            }
            else
            {
                m_EPIKernel.setExpectedMaxLinesForPEBlip(0);
                m_EPIKernel.setExpectedMaxPartitionsFor3DBlip(0);
            }

            // Multi-shot acquisition with multiple echoes per train: EPI kernel requires knowledge about echo shift
            // between segments
            if (m_REOInfo.isMultiShot() && (m_REOInfo.getEchoTrainLength() > 1))
            {
                if (m_REOInfo.isPATActive())
                {
                    m_EPIKernel.setUseEchoShifting(
                        true, m_REOInfo.getPATLinesPerSegment(), m_REOInfo.getCounterInSegmentWithKSpaceCenter());
                }
                else
                {
                    m_EPIKernel.setUseEchoShifting(
                        true, m_REOInfo.getLinesPerSegment(), m_REOInfo.getCounterInSegmentWithKSpaceCenter());
                }
            }
        }
    }
#endif

    // NOTE: the following part was not directly after the parts above; but after the many protocol checks for TGSE, currently being after the configureEPIKernel() function call in ep2d.cpp
    // it may or may not change behavior (tests does not indicate this to be a case, but maybe in some edge cases...)

    // set the number of interleaves to 1 for single-shot sequence
    // this information is required by the m_EPIKernel ONLY for the function calcEffEchoSpacingAndBWPerPixelPE()
    // the parameter has NO effect on the EPI readout gradient waveform
    m_EPIKernel.setNumInterleaves(1);

#ifdef EP2D_MS
    m_EPIKernel.setNumInterleaves(rMrProt.getsFastImaging().getlSegments());
#endif

    return MRI_SEQ_SEQU_NORMAL;
}

#ifdef WIN32
void Ep2d::setUIThermalBalancing()
{
    // functionality only in derived Ep2d_diff
}
#endif

void Ep2d::setVariantSpecificLoopSettings(SeqLim& rSeqLim, MrProt& rMrProt)
{
    m_mySeqLoop.setdDistFacMinIfConcNo(rSeqLim.getSliceDistanceFactor().getMin());
}

bool Ep2d::initializeEPIKernel()
{
#if defined COMPILE_EP2D_SE || defined EP2D_MS 
    // Enable balance model
    // (if current gradient system supports GPA balance models)
    m_EPIKernel.setUseGPABalance(true, SysProperties::getGradMaxAmpl(SEQ::GRAD_FAST));
#endif

    // set polarity of first RO pulse in echo-train to positive
    m_EPIKernel.setStartImagingReadOutWithNegativeGradient(false);

    // type of phase correction is set later
    m_EPIKernel.setInternalPhaseCorrection(false);

    // activate RO ramp-sampling
    m_EPIKernel.setUseRegriddingForRO(true);

    // set slew rate to 95% of "fast" value specified in MeasPerm section (CHARM 32446)
    // note that this can be changed to 80% by switching to "normal" gradient mode (see fSEQPrep)
    // Also note that this call also sets the default gradient performance (SBBEPIReadout)
    m_EPIKernel.setMinRiseTimeScaleFactor(1.05);

    //-------------------------------------------------------------------------------------
    // we do not use calculation limits
    // => full available gradient performance is used all the time
    //-------------------------------------------------------------------------------------
    m_myCalcLimits.resetAllLimits();

    if (!m_EPIKernel.setPointerToCalculationLimits(&m_myCalcLimits))
    {
        SEQ_TRACE_ERROR.print("m_EPIKernel.setPointerToCalculationLimits failed: 0x%lx", m_EPIKernel.getNLSStatus());
        return false;
    }

    //-------------------------------------------------------------------------------------
    // preparation of osc-bit
    //-------------------------------------------------------------------------------------
#ifdef EP2D_SE_MRE
    // Extend the duration for the driver to catch the TTL trigger.
    if (!m_EPIKernel.setUseSyncBits(true, 0, 10, 0, EXT_TRIGGER_DURATION_US, 0))
#else
    if (!m_EPIKernel.setUseSyncBits(true)) // using default arguments from SeqBuildBlockEPIKernel::setUseSyncBits
#endif
    {
        SEQ_TRACE_ERROR.print("m_EPIKernel.setUseOscBit failed: 0x%lx", m_EPIKernel.getNLSStatus());
        return false;
    }

#ifdef ASL
    m_EPIKernel.setPrephaseBlipsAfterRTEBPlugIn(true);
#endif

    return true;
}

void Ep2d::configureReorderInfo(MrProt& rMrProt)
{
    if (rMrProt.getsFastImaging().getlSegments() > 1)
    {
        m_REOInfo.setModeMultiShot();
    }
    else
    {
        m_REOInfo.setModeSingleShot();
    }


#ifdef SUPPORT_iPAT_a_ep2d

    if (isGreRefScanType(rMrProt))
    {
        // Switch ReorderInfoEPI to external reference scan mode
        m_REOInfo.setPATMultiShotRefScans(false);
        m_REOInfo.setPATGRERefScans(true);
    }
    else
    {
        // Internal reference scans
        m_REOInfo.setPATMultiShotRefScans(m_bSegmentedRefLines);
        m_REOInfo.setPATGRERefScans(false);
    }

#endif
}

void Ep2d::setPartialFourierToReorderInfo(MrProt& rMrProt, ReorderInfo* pReorderInfo) const
{
    // functionality only in derived Ep2d_diff
}

void Ep2d::initializeOscBitFlags()
{
    std::fill(m_abOscBitSentForMeas.begin(), m_abOscBitSentForMeas.end(), false);
}

void Ep2d::prepareDFCBookkeeping()
{
    // functionality only in derived Ep2d_diff
}

void Ep2d::handleExternalPhaseCorrectionRun(MrProt& rMrProt, long lKernelMode, long lSlice)
{
    if (rMrProt.getsPrepPulses().getlPhaseCorrectionMode() == MrProtocolData::PHASECORR_EXTERNAL)
    {
        long lFirstPhaseCorrScan = m_lInitialDummyScans;
        long lLastPhaseCorrScan  = lFirstPhaseCorrScan + m_lPhaseCorrPrepScans - 1;

        // Consider slice-accelerated (collapsed) phase correction scans
        if (SMSProperties::isSMS(rMrProt))
        {
            if (isFastGreRefScan(rMrProt))
            {
                lFirstPhaseCorrScan = m_lInitialDummyScans + m_lPhaseCorrPrepScans + m_lSliceAccelRefScans
                                      + m_lPATRefScans + m_lSliceAccelDummyScans;
                lLastPhaseCorrScan = lFirstPhaseCorrScan + m_lSliceAccelPhaseCorrScans - 1;
            }
        }

        if ((lKernelMode == KERNEL_PREPARE) && (m_alPrepScanCounter[lSlice] >= lFirstPhaseCorrScan)
            && (m_alPrepScanCounter[lSlice] <= lLastPhaseCorrScan))
        {
            // Re-enable readouts for phase correction scans
            fRTSetReadoutEnable(1);
            // Disable diffusion gradients for phase correction scans
            disableDiffusionForPrepScan();
            // Configure kernel for phase correction
            m_EPIKernel.setExecuteKernelAsPhaseCorrectionScan(true);
            m_EPIKernel.setExecuteExternalEPIPhaseCorrectionScan(true);

            // Repetitions loop counter of phase correction scans is always zero
            // (setLoopCounters might have set a different value)
            m_EPIKernel.getReadOutAddress()->getMDH().setCrep(0);
            // No LastScanInMeas flag for phase correction scans
            // (might have been set above)
            m_EPIKernel.getReadOutAddress()->getMDH().deleteFromEvalInfoMask(MDH_LASTSCANINMEAS);
        }
        else
        {
            m_EPIKernel.setExecuteKernelAsPhaseCorrectionScan(false);
            m_EPIKernel.setExecuteExternalEPIPhaseCorrectionScan(false);
        }
    }
    else
    {
        m_EPIKernel.setExecuteKernelAsPhaseCorrectionScan(false);
        m_EPIKernel.setExecuteExternalEPIPhaseCorrectionScan(false);
    }
}

void Ep2d::setReducedIRThicknessForShortTI()
{
    m_mySeqLoop.setScaleIRThickness(1.25);
}

MrProtocolData::SeqExpoRFInfo Ep2d::calcEnergyOfExtraCSat(long lNumberOfPulses)
{
    // functionality only in derived Ep2d_diff. here return zero block (generated by default constructor of the RFInfo)
    return {};
}

void SEQ_NAMESPACE::Ep2d::setLastScanInMeasFlagForB0Correction(MrProt& rMrProt, long lSlice, long lShot)
{
    // functionality only in derived Ep2d_diff
}

void Ep2d::disableDiffusionForPrepScan()
{
    // functionality only in derived Ep2d_diff
}

NLS_STATUS Ep2d::setLastScanInMeasFlagForDFCBookkeeping(MrProt& rMrProt, long lKernelMode)
{
    // functionality only in derived Ep2d_diff
    return MRI_SEQ_SEQU_NORMAL;
}

void Ep2d::configureDiffusionSpecificSeqUTSettings(SeqLim& rSeqLim, MrProt& rMrProt)
{
    // functionality only in derived Ep2d_diff
}

long Ep2d::getDiffusionAdjPrepScans() const
{
    // gets the number of diffusion adjustment preparation scans in the derived Ep2d_diff, returns 0 in other flavors.
    return 0;
}

void Ep2d::setTriggerAndOscBit(long lKernelMode)
{
    {
        // CHARM 304696: send osc bit only once per measurement
        //
        if (m_EPIKernel.getReadOutAddress()->getMDH().getCrep() > m_lMaxOscBitSentFlags - 1)
        {
            SEQ_TRACE_ALWAYS.print("m_EPIKernel.getReadOutAddress()->getMDH().getCrep() > m_lMaxOscBitSentFlags-1");
            SEQ_TRACE_ALWAYS.print("=> sending osc-bit uncontrolled");

            m_EPIKernel.setDoNotSendOscBit(false);
            m_EPIKernel.setDoNotSendExtTrigger(false);
        }
        else
        {
            if (m_abOscBitSentForMeas[m_EPIKernel.getReadOutAddress()->getMDH().getCrep()]
                || lKernelMode == KERNEL_PREPARE)
            {
                m_EPIKernel.setDoNotSendOscBit(true);
                m_EPIKernel.setDoNotSendExtTrigger(true);
            }
            else
            {
                m_EPIKernel.setDoNotSendOscBit(false);
                m_EPIKernel.setDoNotSendExtTrigger(false);
                m_abOscBitSentForMeas[m_EPIKernel.getReadOutAddress()->getMDH().getCrep()] = true;
            }
        }
    }

#ifdef EP2D_SE_MRE
    // Re-enable the External trigger for all TRs - Override dependency on the REP counter
    m_EPIKernel.setDoNotSendExtTrigger(false);
#endif
}

void Ep2d::setSPAIRSpoilingType()
{
    // functionality only in derived Ep2d_diff sequence
}

void Ep2d::considerIRBlockForImplicitCoolingPause(
    SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo, long lScanTimeSatsEtc, long lScanTimeBasic)
{
    // SBB's with low gradient activity that are played out with each kernel
    // might be considered as contributions to the pause. SPAIR duration
    // is known only after calcSPIRTime has been executed.
    if ((rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
        && (m_mySeqLoop.getInterleavedIRAllowed() == false))
    {
        //  We need the TI-Fill time, to determine the implicit pause.
        //  However the fill time has not been calculated, yet.
        //  Since the TI-Fill time is independent of the kernel duration,
        //  we first call SeqLoop::calcFillTimes with zero pause ...
        m_mySeqLoop.setlSBBScanTime(lScanTimeSatsEtc);
        m_mySeqLoop.setlKernelScanTime(lScanTimeBasic);
        if (m_mySeqLoop.calcFillTimes(rMrProt, rSeqLim, rSeqExpo))
        {
            const SeqConcat* pSeqConcat = m_mySeqLoop.getsSeqConcat(0);
            if (pSeqConcat != 0)
            {
                //  Note: in the case of sequential IR the TI-Fill time is independent of the actual concatenation.
                m_lCoolPauseImplicit += pSeqConcat->m_alTIFill[0];
                //  The IR pulse is comparatively long (at the time of implementation: 20640 us)
                //  and the slice selection gradient is usually small. Hence the RF time is considered
                //  The spoiler is moderate (at the time of implemenation 8 mT/m) and therefore the spoiler time is
                //  not considered.
                m_lCoolPauseImplicit += const_cast<SeqBuildBlockIRsel*>(&m_mySeqLoop.getSBBIRsel())->getIRRFDuration();
            }
        }
    }
}

bool Ep2d::isSegmentedPATRefLinesCondition(MrProt& rMrProt) const
{
    if ((rMrProt.PAT().getucPATMode() != SEQ::PAT_MODE_NONE) && (rMrProt.PAT().getlAccelFactPE() > 2)
        && (rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI))
    {
        return true;
    }
    return false;
}

void Ep2d::disableReadoutForPrepScans(long lKernelMode, long lSlice/*=0*/)
{
#ifdef SUPPORT_iPAT_a_ep2d
    {
        if (lKernelMode == KERNEL_PREPARE)
        {
            fRTSetReadoutEnable(0);
        }
    }
#endif

#ifdef ASL

    if (m_bEnableFirstPrepScanAsM0Scan)
    {
        // disable readout for further PrepScans after M0 reference scan
        if (lKernelMode == KERNEL_PREPARE)
        {
            if (m_alASLPrepScanCounter[lSlice] == 1)
            {
                fRTSetReadoutEnable(0);
            }
        }
    }

#endif // ifdef ASL
}

bool Ep2d::isLoopStructureCompatibleWithB0Correction(MrProt& rMrProt) const
{
    if (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF)
    {
        return false;
    }

    // no long TR triggering mode for ep2d_se, ep2d_fid and ep2d_asl sequences
    // (see m_mySeqLoop configuration below)
#if ((defined COMPILE_EP2D_SE) || (defined COMPILE_EP2D_FID) || (defined ASL))
    if (rMrProt.concatenations() > 1)
    {
        return false;
    }
#endif

    return true;
}

void Ep2d::setLongTRTrigMode(MrProt& rMrProt)
{   
    // switch on by default
    m_mySeqLoop.setLongTRTrigMode(true);

    // no long TR triggering mode if navigator triggering is switched on
#ifdef SUPPORT_PACE

    if (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF)
    {
        m_mySeqLoop.setLongTRTrigMode(false);
    }
#endif

    // no long TR triggering mode for ep2d_se, ep2d_fid and ep2d_asl sequences
#if ((defined COMPILE_EP2D_SE) || (defined COMPILE_EP2D_FID) || (defined ASL))

    m_mySeqLoop.setLongTRTrigMode(false);

#endif
}

bool Ep2d::isCoolTimeExecutedWithIR(MrProt& rMrProt) const
{
    // actual check only in derived Ep2d_diff; in non-diffusion flavors it's always true
    return true;
}

bool Ep2d::isMultiConcatsAllowed(MrProt& rMrProt) const
{
    bool bMultiConcatsAllowed = false;

    // Multiple concatenations are allowed if long TR triggering mode is enabled
    // and standard triggering is active

    if (m_mySeqLoop.isLongTRTrigMode())
    {
        SEQ::PhysioSignal FirstSignal;
        SEQ::PhysioMethod FirstMethod;
        SEQ::PhysioSignal SecondSignal;
        SEQ::PhysioMethod SecondMethod;

        rMrProt.physiology().getPhysioMode(FirstSignal, FirstMethod, SecondSignal, SecondMethod);

        if (FirstMethod == SEQ::METHOD_TRIGGERING)
        {
            bMultiConcatsAllowed = true;
        }
    }

    // Multiple concatenations are allowed if navigator triggering is active

#ifdef SUPPORT_PACE
    if (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF)
    {
        bMultiConcatsAllowed = true;
    }
#endif

    // Multiple concatenations are allowed with ep2d_se and ep2d_fid
    // sequences if there is a single repetition

#if ((defined COMPILE_EP2D_SE) || (defined COMPILE_EP2D_FID))

    if (rMrProt.getlRepetitions() == 0)
        bMultiConcatsAllowed = true;

#endif

    return bMultiConcatsAllowed;
}

bool Ep2d::isExtraCSatApplied(MrProt& rMrProt)
{
    // Decide whether an additional CSat is played out at the beginning of the
    // cool pause.
    //    Condition #1: explicit cool pause >= 2 * CSat duration
    //    Condition #2: ratio explicit cool pause vs. effective TR >= 0.30 (empirical value)
    // Note: the duration of the additional CSat will be considered as a
    //       contribution to the cool pause.
    // Note: No additional CSat will be applied if respiration triggering is enabled.
#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
    if ((rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_FatSaturation)
        && (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_TRIGGER)
        && (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_TRIGGER_AND_FOLLOW)
#ifdef COMPILE_EP2D_DIFF
        && !m_EPIKernel.getbCompensationEnable()
#endif // COMPILE_EP2D_DIFF
    )
    {
        const double dCoolPauseExplicit = static_cast<double>(m_lCoolPauseTotal - m_lCoolPauseImplicit);
        const double dCSatDuration
            = static_cast<double>(m_CSatFat.getDurationPerRequest() + m_SpoilGrad.getDurationPerRequest());
        const double dTReff = static_cast<double>(m_mySeqLoop.getlKernelScanTime());

        if (dTReff > 0.)
        {
            if ((dCoolPauseExplicit >= 2. * dCSatDuration) && (dCoolPauseExplicit / dTReff >= 0.30))
            {
                return true;
            }
        }
    }
#endif

    return false;
}

void Ep2d::setZoomItSeqUTExpectations(MrProt& rMrProt)
{
    SeqUT.SetExpectedOk(
        lAmplSignRFErr, RTEB_ORIGIN_fSEQRunKernel, 0, "ZOOM_2DRF slice selection has a 2D trajectory");
    SeqUT.SetExpectedOk(
        lRFAmplValErr,
        RTEB_ORIGIN_fSEQRunKernel,
        0,
        "ZOOM_2DRF slice-select gradient amplitude is alternating with trajectory. Image orientation and image "
        "geometry are checked in a functional scanner test");

    SeqUT.SetExpectedOk(
        lMoreGrInXRFErr,
        RTEB_ORIGIN_fSEQRunKernel,
        0,
        "ZOOM_2DRF has more than one active gradient during the RF_PULSE event ()");
}

bool Ep2d::checkDFCWithLongTRTrigMode(SeqLim& rSeqLim, MrProt& rMrProt) const
{
    // functionality only in derived Ep2d_diff
    return true;
}

bool Ep2d::checkAcqWindowForRespTriggering(MrProt& rMrProt) const
{
    // functionality only in derived Ep2d_diff
    return true;
}

long Ep2d::calcTotalNumberOfVolumesForFreqFeedback(MrProt& rMrProt)
{
    // TODO: check this; in the original implementation, lTotalNoOfVolumes do not exist outside of BOLD and DIFF, but is needed if FreqFeedback is active!
    long lTotalNoOfVolumes = 0;

#ifdef BOLD
    // Total number of expected volumes, including preparation and adjustment scans
    lTotalNoOfVolumes = rMrProt.measurements() + m_lPATRefScans;
#endif

    return lTotalNoOfVolumes;
}

void Ep2d::prepSliceAccelPhaseOffcenter3D(MrProt& rMrProt, sSLICE_POS asSLC[])
{
    MrProtFacade protFacade(rMrProt);

    if (!protFacade.isSliceAcceleration())
        return;

    const long lNumberOfSlices = rMrProt.getsSliceArray().getlSize();

    for (long lI = 0; lI < lNumberOfSlices; ++lI)
    {
        const double dPhaseOffcenter3D
            = SliceAccelerationUtils::calcSliceAccelerationPhaseOffcenter3D(rMrProt, asSLC[lI]);

        asSLC[lI].setPhaseOffCenter3D(dPhaseOffcenter3D);
    }
}

long Ep2d::getTEContrastIndex(long lNumberOfContrasts) const
{
#if defined EP2D_MS && defined COMPILE_EP2D_SE
    // User specifies value of the longest TE
    return lNumberOfContrasts - 1;
#else
    // User specifies value of the shortest TE
    return 0;
#endif
}

long Ep2d::getNumberOfContrastsBeforeRTEBPlugIn(long lNumberOfContrasts) const
{
#if defined EP2D_MS && defined COMPILE_EP2D_SE
    // Number of contrasts after PlugIn: equal or larger by one than the number before PlugIn
    return lNumberOfContrasts - (lNumberOfContrasts + 1) / 2;
#else
    // No contrasts before the PlugIn
    return 0;
#endif
}

bool Ep2d::isIIRSchemeStandard(MrProt& rMrProt, SeqLim& rSeqLim)
{
#if defined SUPPORT_IIR
    return m_mySeqLoop.IsIIRSchemeStd( rMrProt, rSeqLim );
#else
    return true;
#endif
}

#ifdef COMPILE_EP2D_DIFF
bool Ep2d::updateCompGrad(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo)
{
    // scaled factors are with RPS ordering
    std::vector<double> vdScaleFactorinRun = m_EPIKernel.getvdScaleFactorinRun();
    // reset compensation gradient amplitude on all axes according to the applied diffusion gradients

    m_EPIKernel.getPointerCompGrad()->scaleAmplitude(
        /*RO:*/ -1.0 * vdScaleFactorinRun[1],
        /*PE:*/ -1.0 * vdScaleFactorinRun[0],
        /*SS:*/ -1.0 * vdScaleFactorinRun[2]);

    return true;
}
#endif
