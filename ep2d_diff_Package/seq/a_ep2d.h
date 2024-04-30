//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d.h
//     Version: \main\52
//      Author: Clinical
//        Date: 2015-03-06 08:04:09 +01:00
//
//        Lang: C++
//
//     Descrip: Declarations for a_ep2d.cpp
//
//
/// \brief  File containing declarations for the sequences
///         - ep2d_diff
///         - ep2d_perf
///         - ep2d_se
///         - ep2d_fid
///         - ep2d_bold
///         - ep2d_pace
///
/// This file contains the declaration of the class Ep2d.
///
//  ***************************************************************************/
///

#pragma once

#ifndef a_ep2d_h
#define a_ep2d_h

#define ICEPROGRAMPARA_SHOT_INDEX      ( MDH_NUMBEROFICEPROGRAMPARA - 1 )

//========================================
// defines for ep2d_bold and ep2d_pace
//========================================
#ifdef BOLD

//----------------------------------------------------------------------------------------------------
// Default setting for limit to the number of receiver channels for different matrix sizes.
// The 4 elements correspond to 64, 128, 192 and 256 matrix sizes respectively.
// These values can be overridden using the sequence registry key PARAM_EPI_BOLD_ReceiverChannelLimits.
//----------------------------------------------------------------------------------------------------
#define EPI_RECEIVER_CHANNEL_LIMITS_DEFAULT \
    {                                       \
        16, 8, 1, 1                         \
    }

//----------------------------------------------------------------------------------------------------
// This define causes the 'Save uncombined' switch to be removed.
// Reason: 'Save uncombined' option has two potential pitfalls:
// 1. It is incompatible with mosaic images.
// 2. It may cause trouble in the database when many images are measured with many channels enabled.
//----------------------------------------------------------------------------------------------------
#define EPI_DISABLE_SAVE_UNCOMBINED

//========================================
// defines for ep2d_bold only
//========================================
#ifndef PACE3D
//----------------------------------------------------------------------------------------------------
// This define activates the frequency feedback.
// EPI navigator data is used to estimate the scanner frequency drift
// (relative ot the first acquired volume). A realtime feedback mechanism
// is used to provide the sequence with the evaluated frequency offset:
// this information is used to update the scanner frequency (involves a
// temporal filtering).
//----------------------------------------------------------------------------------------------------
#define EPI_SUPPORT_FREQ_FEEDBACK
#endif

#endif


//========================================
// defines for ep2d_fid
//========================================
#ifdef PERF

//----------------------------------------------------------------------------------------------------
// This define causes the 'Save uncombined' switch to be removed.
// Reason: incompatible with perfusion post-processing (CHARM 388035)
//----------------------------------------------------------------------------------------------------
#define EPI_DISABLE_SAVE_UNCOMBINED

#endif

//  --------------------------------------------------------------------------
//  General Includes
//  --------------------------------------------------------------------------
#include "MrGlobalDefinitions/MrResult.h"
#include "MrImagingFW/libSBBFW/StdSeqIF.h"
#ifdef SUPPORT_FAST_IR
#include "MrImaging/seq/common/IR/SeqLoopFastIR.h"
#elif defined SUPPORT_IIR
#include "MrImaging/libSBB/SeqLoopIIR.h"
#include "MrImaging/seq/a_ep2d_se_ms/SeqLoopIIR_msEPI.h"
#else
#include "MrImaging/seq/common/SeqLoopLongTRTrig/SeqLoopLongTRTrig.h"
#endif

#include "MrMeasSrv/SeqFW/libSSL/SSLProfile.h"


//  --------------------------------------------------------------------------
//  Include centralized SeqLoopEP2D structure
//  --------------------------------------------------------------------------
#ifdef EP2D_SE_MRE
#include "MrImaging/seq/a_ep2d_se_mre/SeqLoopMRE.h"
#else
#include "MrImaging/seq/common/IR/SeqLoopFastIR.h"
#include "MrImaging/seq/common/SeqLoopEPI/SeqLoopEP2D.h"
#include "MrImaging/seq/common/SeqLoopEPI/SeqLoopMultiband.h"
#endif

// SMS support
#include "MrImaging/libSBB/SBBMultibandRF.h"
// SMS support

//  --------------------------------------------------------------------------
//  Application includes
//  --------------------------------------------------------------------------
#ifdef PACE3D
#include "MrImaging/seq/a_ep2d_pace/PaceFeedback.h"
#endif

#ifdef EPI_SUPPORT_FREQ_FEEDBACK
#include "MrImaging/seq/a_ep2d_diff/FreqFeedback.h"
#endif // #ifdef EPI_SUPPORT_FREQ_FEEDBACK

#if defined COMPILE_EP2D_DIFF
// Variant specific includes
#include "MrImaging/seq/a_ep2d_diff/DiffusionSBBContainer.h"
#include "MrImaging/seq/a_ep2d_diff/SBBEPIKernelDiffusion.h"

// Variant specific UI handlers
#include "MrImaging/seq/a_ep2d_diff/a_ep2d_diff_UI.h"

namespace SEQ_NAMESPACE
{
class Ep2d_diff_UI;        // Forward declaration
typedef Ep2d_diff_UI EpUI; // Variant specific UI hook
} // namespace SEQ_NAMESPACE
#elif defined COMPILE_EP2D_SE

#if defined EP2D_MS && !defined EP2D_MS_WIP // multishot SE EPI
#include "MrImaging/seq/a_ep2d_se_ms/SBBEPIKernelSE_ms.h"
#include "MrImaging/seq/a_ep2d_se_ms/a_ep2d_se_ms_UI.h"

namespace SEQ_NAMESPACE
{
class Ep2d_se_ms_UI;        // Forward declaration
typedef Ep2d_se_ms_UI EpUI; // Variant specific UI hook
} // namespace SEQ_NAMESPACE

#else
#include "MrImaging/seq/a_ep2d_se/SBBEPIKernelSE.h"
#include "MrImaging/seq/a_ep2d_se/a_ep2d_se_UI.h"

namespace SEQ_NAMESPACE
{
class Ep2d_se_UI;        // Forward declaration
typedef Ep2d_se_UI EpUI; // Variant specific UI hook
} // namespace SEQ_NAMESPACE
#endif //if defined EP2D_MS

#elif defined ASL

#include "MrImaging/SequenceLibraries/libASL/ASLCommon.h"
#include "MrImaging/SequenceLibraries/libASL/SBBAsl.h"

#ifdef COMPILE_EP2D_FID
#include "MrImaging/seq/a_ep2d_asl/SBBEPIKernelASL2D.h"
#endif

#include "MrImaging/seq/a_ep2d_asl/a_ep2d_asl_UI.h"

namespace SEQ_NAMESPACE
{
class ep2d_ASL_UI;        // Forward declaration
typedef ep2d_ASL_UI EpUI; // Variant specific UI hook
} // namespace SEQ_NAMESPACE

#elif defined COMPILE_EP2D_FID

#if defined PERF || defined EP2D_MS_WIP
#include "MrImaging/seq/a_ep2d_fid/SBBEPIKernelFID.h"
#include "MrImaging/seq/a_ep2d_fid/a_ep2d_fid_UI.h"

namespace SEQ_NAMESPACE
{
class Ep2d_fid_UI;        // Forward declaration
typedef Ep2d_fid_UI EpUI; // Variant specific UI hook
} // namespace SEQ_NAMESPACE

#elif defined EP2D_MS && !defined EP2D_MS_WIP // multishot EPI FID
#include "MrImaging/seq/a_ep2d_fid/SBBEPIKernelFID.h"
#include "MrImaging/seq/a_ep2d_fid_ms/a_ep2d_fid_ms_UI.h"

namespace SEQ_NAMESPACE
{
    class Ep2d_fid_ms_UI;        // Forward declaration
    typedef Ep2d_fid_ms_UI EpUI; // Variant specific UI hook
} // namespace SEQ_NAMESPACE

#else // BOLD

#include "MrImaging/seq/a_ep2d_bold/a_ep2d_bold_UI.h"
#include "MrImaging/seq/a_ep2d_fid/SBBEPIKernelFID.h"

namespace SEQ_NAMESPACE
{
class Ep2d_bold_UI;        // Forward declaration
typedef Ep2d_bold_UI EpUI; // Variant specific UI hook
} // namespace SEQ_NAMESPACE
#endif

#else // unflavored epi
#include "MrImaging/seq/a_ep_CommonUI.h"

namespace SEQ_NAMESPACE
{
class EpCommonUI;        // Forward declaration
typedef EpCommonUI EpUI; // Variant specific UI hook
} // namespace SEQ_NAMESPACE
#endif

#ifdef ZOOM_2DRF
// Variant specific includes: 2D RF Excitation
#include "MrImaging/libSeqPTX/SBBOptPTXVolume.h"
#endif

#include "MrImaging/libSBB/SBBBinomialPulses.h"
#include "MrImaging/libSeqUtil/ReorderInfoEPI.h"
#ifdef EP2D_SE_MRE
#include "MrImaging/seq/a_ep2d_se_mre/SBBEPIKernelSEMRE.h"
#else
#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#endif
#include "MrImaging/seq/SystemProperties.h"
#include "MrMeasSrv/MeasUtils/MeasMath.h" // minimum/maximum

// MRProtFacade
#include "MrImaging/seq/a_ep2d_diff/SequenceDebugSettings.h"
#include "MrImaging/seq/common/MrProtFacade/IMrProtFacade.h"
#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

namespace SEQ_NAMESPACE
{
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
//                 is not used the default values for the variables (1.0 and 0.0) are provided.
//
//  Return      :  true  - if B1 control loop is used
//                 false - else
//
//  ----------------------------------------------------------------------
bool __IMP_EXP useB1ControlLoop(MrProt& rMrProt, float& fCorrectionFactorMax, float& fPeakReserveFactor);

/*###########################################################################

C l a s s:   S e q l o o p E P 2 D
     --> Is now located here:
    "MrServers/MrImaging/seq/common/SeqLoopEPI/SeqLoopEP2D.h"
###########################################################################*/

//  --------------------------------------------------------------------------
//
/// \brief <b> Class definition of Ep2d. This class is used by the sequences
///         - ep2d_diff
///         - ep2d_perf
///         - ep2d_se
///         - ep2d_fid
///         - ep2d_bold
///         - ep2d_pace    </b>
///
/// This file contains the declaration of the class Ep2d.
///
//  --------------------------------------------------------------------------
class __IMP_EXP Ep2d : public StdSeqIF
{
  public:
    //  ------------------------------------------------------------------
    //
    //  Name        :  Ep2d::Ep2d
    //
    //  Description :
    /// \brief         Initialization of class members
    //
    //  Return      :  void
    //
    //  ------------------------------------------------------------------
    Ep2d();

    //  ------------------------------------------------------------------
    //
    //  Name        :  Ep2d::~Ep2d
    //
    //  Description :
    /// \brief         Destructor. Deletes existing Ep2dUI instances.
    //
    //  Return      :  void
    //
    //  ------------------------------------------------------------------
    virtual ~Ep2d();

    //  ------------------------------------------------------------------
    ///  Copy constructor not implemented
    //  ------------------------------------------------------------------
    Ep2d(const Ep2d& right) = delete;

    //  ------------------------------------------------------------------
    ///  Assignment operator not implemented
    //  ------------------------------------------------------------------
    Ep2d& operator=(const Ep2d& right) = delete;

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
    virtual NLSStatus initialize(SeqLim& rSeqLim);

    //  --------------------------------------------------------------------------
    //
    //  Name        :  Ep2d::prePrepare
    //
    //  Description :
    /// \brief <b>     Check on protocol, if there are forbidden combinations of  parameters.
    ///                This speeds up the binary search, because the method is called very early in
    ///                the prepare process and therefore can save lots of lengthy calculations.
    ///                The method must not prepare any objects like pulses or SBBs.</b>
    //
    //  Return      :  NLS status
    //
    //  --------------------------------------------------------------------------
    virtual NLSStatus prePrepare(const MrProt& rMrProt, const SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    //  --------------------------------------------------------------------------
    //
    //  Name        :  Ep2d::prepare
    //
    //  Description :
    /// \brief <b>     Preparation of the sequence during binary search and prior
    ///                 to sequence execution  </b>
    //
    //  Return      :  NLS status
    //
    //  --------------------------------------------------------------------------
    virtual NLSStatus prepare(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    //  --------------------------------------------------------------------------
    //
    //  Name        :  Ep2d::check
    //
    //  Description :
    /// \brief  <b>    Check of the sequence for gradient stimulation </b>
    ///
    ///                This method is called by the framework prior to a
    ///                 measurement on the host to ensure, that
    ///                 - no gradient overflow occurs
    ///                 - the stimulation will not exceed the threshold
    ///
    //  Return      :  NLS status
    //
    //  --------------------------------------------------------------------------
    virtual NLSStatus check(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, SEQCheckMode* pSEQCheckMode);

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
    virtual NLSStatus run(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    //   --------------------------------------------------------------------------
    //
    //   Name        :  Ep2d::runKernel
    //
    //   Description :
    ///  \brief <b>     Executes the basic timing of the real-time sequence.   </b>
    ///
    ///                 The method runKernel plays out a sequence "Kernel",
    ///                  consisting of one or more lines in k-Space.
    ///                 It is called by SeqLoop.
    ///
    //   Return      :  NLS status
    //
    //   --------------------------------------------------------------------------
    virtual NLS_STATUS runKernel(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lKernelMode, long lSlice, long lPartition, long lLine);

#ifdef WIN32
    // ------------------------------------------------------------------------------
    //   Name        : Ep2d::convProt
    // ------------------------------------------------------------------------------
    //
    //   Description :
    ///  \brief <b>    Try to convert protocols from previous software versions. </b>
    //
    //
    //   Return      : MRI_SEQ_SEQU_NORMAL for success
    //                 MRI_SEQ_SEQU_ERROR  for error
    //
    // ------------------------------------------------------------------------------
    virtual NLSStatus convProt(const MrProt& rMrProtSrc, MrProt& rMrProtDst);
#endif

#if (defined EPI_SUPPORT_FREQ_FEEDBACK) || (defined SUPPORT_PACE) || (defined PACE3D)
    // ------------------------------------------------------------------------------
    // Function    : Ep2d::receive
    // ------------------------------------------------------------------------------
    //
    // Description : Receives and processes feedback data.
    //
    // Return      : NLSStatus.
    //
    // ------------------------------------------------------------------------------
    virtual NLSStatus receive(SeqLim&, SeqExpo&, const SEQData&) override;

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EP2D::cancel                                        *
    // *                                                                    *
    // * Description :  Requests to release all semaphores ...              *
    // *                                                                    *
    // * Return      :  NLS status                                          *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    virtual NLSStatus cancel();
#endif

#if (defined COMPILE_EP2D_DIFF) || (defined SUPPORT_PACE)
    // ------------------------------------------------------------------------------
    // Function    : Ep2d::resume
    // ------------------------------------------------------------------------------
    //
    // Description : Triggers to continue with data acquisition after breath-hold command.
    //
    // Return      : NLSStatus.
    //
    // ------------------------------------------------------------------------------
    virtual NLSStatus resume(const SEQData& rSEQData) override;
#endif

#if defined ZOOM_2DRF
    //  --------------------------------------------------------------------------
    //
    //  Name        :  Ep2d::calculatePTX
    //
    //  Description :
    ///     \brief     Calculates the sequence's pTX pulses.
    ///                The sequence limits must NOT be modified my the sequence.
    ///
    ///                This method is optional, it should be implemented by
    ///                sequences which use pTX RF pulses.
    //
    //  Return      :  NLS status
    //
    //  --------------------------------------------------------------------------
    virtual NLSStatus calculatePTX(MrProt& rMrProt, const SeqLim& rSeqLim);
#endif // ZOOM_2DRF

    //  --------------------------------------------------------------
    //
    //  Name        :  getUI
    //
    //  Description :
    /// \brief <b>     Returns the pointer to the Ep2d UI class  </b>
    ///
    ///                This method is only sensible on the host.
    ///                On the measurement system, it will return an nearly empty object.
    ///
    //  Return      :  EpUI* (variant specific)
    //
    //  --------------------------------------------------------------
    EpUI* getUI() const;

    //  ------------------------------------------------------------------
    //  Declare additional P U B L I C member functions, e.g. functions
    //  that will be accessed in UI handlers (outside the class).
    //  ------------------------------------------------------------------
    bool calculateTRTIFillTimes(MrProt&, SeqLim&, SeqExpo&, long*, long*);

    // set PF factors to reorder info (functionality only in derived Ep2d_diff)
    virtual void setPartialFourierToReorderInfo(MrProt& rMrProt, ReorderInfo* pReorderInfo) const;

#ifdef WIN32
    // Exporting diffusion timing (small and large delta) to UI and checking for IVIM condition
    virtual bool exportDiffusionTimingToUI(SeqLim& rSeqLim, MrProt& rMrProt);

    // set thermal balancing flag in the UI (only in derived Ep2d_diff)
    virtual void setUIThermalBalancing();
#endif

    virtual long getTEContrastIndex(long lNumberOfContrasts) const;

    virtual uint16_t getFastGreRefScanBaseRes();

    virtual long getMaxPATFactorForDLRecon() const;

    virtual long getMinPATFactorForDLRecon() const;

    void setForcedFastGRERefscanWithoutPAT(bool isForcedFastGRERefScanWithoutPAT);

    bool isForcedFastGRERefscanWithoutPAT() const;

protected:

    virtual bool isSystemCompatibleWithFlavor() const;

    // since gradients during the IR block (spoiling, slice selection) are small, the duration is considered for the
    // cooling time. different in diffusion and others.
    virtual void considerIRBlockForImplicitCoolingPause(
        SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo, long lScanTimeSatsEtc, long lScanTimeBasic);

    virtual void loadSliceAdjData(SeqLim& rSeqLim, MrProt& rMrProt);
    virtual void dumpSliceAdjData(SeqLim& rSeqLim, MrProt& rMrProt);

    virtual std::string getSequenceVariantText() const;

    // sequence variant-specific hard limit settings in initialize
    virtual void setVariantSpecificHardLimits(SeqLim& rSeqLim);

    // sequence variant-specific exports to SeqExpo, including some parameter combination checks
    virtual NLSStatus setVariantSpecificExports(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo);

    // severe error encountered; return error
    NLSStatus SeverePrepareErrorReturn(NLSStatus lStatus) const;

    // setting the appropriate, variant-specific DICOM contrast text to SeqExpo
    virtual void setDICOMAcquisitionContrast(SeqExpo& rSeqExpo) const;

    // setting initial dummy scans at the beginning, without considering PAT, SMS, ...
    virtual void setInitialDummyScansBasic(MrProt& rMrProt);

    // set external phase correction and modify initial dummy scans accordingly
    virtual void setPhaseCorrScansAndAdjustInitialDummyScans(MrProt& rMrProt);

    // modify initial dummy scans considering SMS
    virtual void setInitialDummyScansSMS(MrProt& rMrProt);

    // calculate number of required prep scans, including dummy scans
    virtual long calcRequiredPrepScans(MrProt& rMrProt);

    // check TE and increase it if necessary / required by protocol parameters (e.g. TOM)
    virtual NLSStatus CheckAndAdaptTE(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo, long& lNeededTE);

    // diffusion adjustment scan settings in RunKernel (only in derived Ep2d_diff sequence)
    virtual void setDiffusionAdjustmentScans(MrProt& rMrProt, long lSlice, long lKernelMode);

    // setting diffusion loop counters in RunKernel (only in derived Ep2d_diff sequence)
    virtual bool setDiffusionLoopCounters(long lSlice, long lKernelMode);

    // configure EPI kernel depending protocol parameters and sequence variant
    virtual NLSStatus configureEPIKernel(SeqLim& rSeqLim, MrProt& rMrProt);

    // setting variant-specific settings in the loop structure
    virtual void setVariantSpecificLoopSettings(SeqLim& rSeqLim, MrProt& rMrProt);

    // configure protocol-independent settings in EPI kernel in initialize()
    virtual bool initializeEPIKernel();

    // configure settings of the ReorderInfo structure
    virtual void configureReorderInfo(MrProt& rMrProt);

    // initialize osc-bit control-flags
    void initializeOscBitFlags();

    // setting external trigger in EPI kernel, and Osc bit
    virtual void setTriggerAndOscBit(long lKernelMode);

    // setting the spoilers before and after SPAIR pulse (only in derived Ep2d_diff)
    virtual void setSPAIRSpoilingType();

    // prepare bookkeeping for Dynamic Field Correction (only in derived Ep2d_diff)
    virtual void prepareDFCBookkeeping();

    // disable diffusion gradients in the runkernel for external phase correction (only in derived Ep2d_diff)
    virtual void handleExternalPhaseCorrectionRun(MrProt& rMrProt, long lKernelMode, long lSlice);

    // Reduced inversion thickness for short inversion times (e.g. STIR) - avoid crosstalk
    // for nested inversion schemes. different for diffusion and others.
    virtual void setReducedIRThicknessForShortTI();

    // calculate the energy of an extra chemical saturation pulse to the RF block (only in derived Ep2d_diff)
    virtual MrProtocolData::SeqExpoRFInfo calcEnergyOfExtraCSat(long lNumberOfPulses);

    // If B0 correction is enabled, complete volumes are acquired innermost.
    // Each volume is handled internally as a separate repetition => set
    // LastScanInMeas flag correspondingly. Exception: PatRefScans.
    // (only in derived Ep2d_diff)
    virtual void setLastScanInMeasFlagForB0Correction(MrProt& rMrProt, long lSlice, long lShot);

    // disabling diffusion gradients for PAT and SMS reference scans (only in derived Ep2d_diff)
    virtual void disableDiffusionForPrepScan();

    // set last scan in meas when all slices have been acquired for one specific repetition counter
    // this is required if dynamic field correction is used with multiple concats
    // (only in derived Ep2d_diff)
    virtual NLS_STATUS setLastScanInMeasFlagForDFCBookkeeping(MrProt& rMrProt, long lKernelMode);

    // set some diffusion-specific SeqUT settings in Run (only in derived Ep2d_diff)
    virtual void configureDiffusionSpecificSeqUTSettings(SeqLim& rSeqLim, MrProt& rMrProt);

    // disable readouts for preparation scans in runKernel
    virtual void disableReadoutForPrepScans(long lKernelMode, long lSlice=0);

    // check acq window when resp trigger is used (only in derived Ep2d_diff)
    virtual bool checkAcqWindowForRespTriggering(MrProt& rMrProt) const;

    // calculate otal number of expected volumes for frequency feedback, including preparation and adjustment scans
    virtual long calcTotalNumberOfVolumesForFreqFeedback(MrProt& rMrProt);

    // get diffusion adjustment prep scan number in Ep2d_diff, returns 0 in other flavors.
    virtual long getDiffusionAdjPrepScans() const;

    // check the combination "dynamic field correction with multiple concats and respiratory compensation"
    // (only in derived Ep2d_diff)
    virtual bool checkDFCWithLongTRTrigMode(SeqLim& rSeqLim, MrProt& rMrProt) const;

    // set Zoomit-specific expectations for SeqUT
    virtual void setZoomItSeqUTExpectations(MrProt& rMrProt);

    // check restrictions on multiple concatenations
    virtual bool isMultiConcatsAllowed(MrProt& rMrProt) const;

    // check whether the loop structure allows the complete volumes necessary for B0 correction.
    virtual bool isLoopStructureCompatibleWithB0Correction(MrProt& rMrProt) const;

    // check if segmented PAT reference lines should be used.
    virtual bool isSegmentedPATRefLinesCondition(MrProt& rMrProt) const;

    // set long TR triggering mode depending on protocol parameters
    virtual void setLongTRTrigMode(MrProt& rMrProt);

    // Decide whether an additional CSat is played out at the beginning of the cool pause.
    virtual bool isExtraCSatApplied(MrProt& rMrProt);

    // check whether cooling time is executed with the given IR scheme in the protocol. (always true in non-diffusion)
    virtual bool isCoolTimeExecutedWithIR(MrProt& rMrProt) const;

    virtual bool isFastGreRefScan(const MrProt& rMrProt);
    virtual bool isGreRefScanType(const MrProt& rMrProt);

    virtual uint16_t getFastGreRefScanBandwidth();

    virtual bool isNumberOfCoilsSufficient(MrProt& rMrProt, SeqLim& rSeqLim);
    virtual long getNumberOfContrastsBeforeRTEBPlugIn(long lNumberOfContrasts) const;

    virtual bool isIIRSchemeStandard(MrProt& rMrProt, SeqLim& rSeqLim);

#ifdef EPI_SUPPORT_FREQ_FEEDBACK
    virtual long getCurrentVolumeForB0Correction();
#endif 

  public:
    //-------------------------------------------------------------------------------------
    // standard loop structure realized within class
    //-------------------------------------------------------------------------------------
#ifdef EP2D_SE_MRE
    SeqLoopMRE m_mySeqLoop;
#elif defined SUPPORT_FAST_IR
    SeqLoopEP2D<SeqLoopMultiBand<SeqLoopFastIR>> m_mySeqLoop;
    #elif defined SUPPORT_IIR
    SeqLoopEP2D<SeqLoopMultiBand<SeqLoopIIR_msEPI>> m_mySeqLoop;
#else
    SeqLoopEP2D<SeqLoopMultiBand<SeqLoopLongTRTrig>> m_mySeqLoop;
#endif

    //-------------------------------------------------------------------------------------
    // Slice position information (rotation matrices and shifts)
    //-------------------------------------------------------------------------------------
    sSLICE_POS m_asSLC[K_NO_SLI_MAX];

    //-------------------------------------------------------------------------------------
    // kernel calculation limits
    //-------------------------------------------------------------------------------------
    KernelCalculationLimits m_myCalcLimits;

#ifdef COMPILE_EP2D_DIFF
    SBBEPIKernelDiffusion m_EPIKernel{nullptr}; // the EPI kernel class
#elif defined COMPILE_EP2D_SE
#ifdef EP2D_SE_MRE
    SBBEPIKernelSEMRE m_EPIKernel{nullptr}; // the EPI kernel class
#elif defined EP2D_MS && !defined EP2D_MS_WIP
    SBBEPIKernelSE_ms m_EPIKernel{nullptr}; // the EPI kernel class
#else
    SBBEPIKernelSE m_EPIKernel{nullptr}; // the EPI kernel class
#endif

#elif defined COMPILE_EP2D_FID

#ifdef ASL
    SBBEPIKernelASL2D m_EPIKernel{nullptr}; // the EPI kernel class
#else
    SBBEPIKernelFID m_EPIKernel{nullptr}; // the EPI kernel class
#endif

#endif

    //-------------------------------------------------------------------------------------
    // reordering information data
    //-------------------------------------------------------------------------------------
    ReorderInfoEPI m_REOInfo;

    //-------------------------------------------------------------------------------------
    // control execution of osc-bit
    //-------------------------------------------------------------------------------------
    enum {m_lMaxOscBitSentFlags = 4096};
    std::array<bool, m_lMaxOscBitSentFlags> m_abOscBitSentForMeas;

    //  ------------------------------------------------------------------
    //  Declare additional member variables of the sequence, here.
    //
    //  Initialization should be done in the constructor of the class.
    //  Maybe, that these variables have to be deleted in the destructor.
    //  ------------------------------------------------------------------
#ifdef PACE3D
    PaceFeedback m_PaceFeedback;
#endif
    bool m_bPrepScansOnlyInFirstMeasurement{true};
    /// Number of initial dummy scans
    /** No data is acquired during these scans - they just serve to
    drive the magnetization into a steady state. Initial dummy
    scans are played out with the preparation scans.
    */
    long m_lInitialDummyScans{0};

    long m_lPhaseCorrPrepScans{0};

    /// This flag indicates whether B0 correction (frequency feedback and/or image shift) is applied
    bool m_bB0Correction{false};
    bool m_bSequentialVolumeAcquisition{true};
    bool m_bIsRealtimeProcessingEnabled{false}; // sequence runs in real-time mode if enabled

#ifdef EPI_SUPPORT_FREQ_FEEDBACK
    // Frequency feedback is supported by:
    // - ep2d_diff
    // - ep2d_bold

    /// This FreqFeedback object is used for synchronization and incorporation of frequency feedback
    FreqFeedback m_sFreqFeedback;

    /// Counter that indicates the currently acquired volume
    /** Similar to the repetitions counter, but all acquired volumes including
    averages, iPAT reference scans and adjustment scans are considered.
    */
    long m_lVolumeCounter{0};
#endif // #ifdef EPI_SUPPORT_FREQ_FEEDBACK


#if defined COMPILE_EP2D_DIFF || defined COMPILE_EP2D_SE
    /// Flag to indicate application of additional CSat
    /** If elevated gradient amplitudes are used for during diffusion encoding, a mandatory
    cool pause is inserted. If fat saturation is enabled and this cool pause is
    sufficiently long, an additional CSat is applied at the beginning of the pause
    to improve fat suppression.
    */
    bool m_bApplyExtraCSat{false};

    /// SBBs used for the additional CSat
    SBBList                m_SBB{};
    SeqBuildBlockCSat      m_CSatFat{&m_SBB};
    SeqBuildBlockSpoilGrad m_SpoilGrad{&m_SBB};
#endif

#if defined ZOOM_2DRF
    SBBList         m_SBBPTX{};
    SBBOptPTXVolume m_OptPTXVolume{&m_SBBPTX};
#endif

#if (defined SUPPORT_iPAT_a_ep2d) || (defined SUPPORT_iPAT_TGSE)
    long m_lMinPrepScansNoPATRefScans{2};
    /// Number of PAT reference scans
    /** PAT reference scans are played out with the preparation scans.
     */
    long m_lPATRefScans{0};
    long m_alPrepScanCounter[K_NO_SLI_MAX];
    bool m_bSegmentedRefLines{false};
    long m_lPATFlashRefScanBaseRes{64};
    long m_lPATFlashRefScanBandwidth{260};

#endif

#ifdef ASL
    // Fat sat all slices for Strong Fat Sat Mode
    // We need to execute fat sat ourselves, instead of SeqLoop executing them because SeqLoop
    // will execute them before the ASL label pulses, and we want the fatsat after the label pulses
    SBBList                m_SBB{};
    SeqBuildBlockCSat      m_CSatFat{&m_SBB};
    SeqBuildBlockSpoilGrad m_SpoilGrad{&m_SBB};

    SBBList                        m_SBBListASL{};
    LIBASL::SeqBuildBlockAsl       m_ASL_SBB{&m_SBBListASL};
    LIBASL::SeqBuildBlockCrushGrad m_PaslCrushGrad{nullptr}; // PCASL: libASL
    bool                           m_bEnableFirstPrepScanAsM0Scan{true};
    // TODO: change these to stl::<vector>
    long m_alASLPrepScanCounter[K_NO_SLI_MAX];
    long m_alASLlScanCounter[K_NO_SLI_MAX];
#endif

    /// Mandatory cool pause (per kernel) [us]
    long m_lCoolPauseTotal{0};
    /// Implicit cool pause (due to SBB's etc.) [us]
    long m_lCoolPauseImplicit{0};

    /// Mask specifying requested dynamic adjustments
    SLICEADJ::sAdjParametersMask m_sSliceAdjParametersRequestedBySequence{}; // Default: Nothing specified (note that this is different from explicitly specifying ADJNONE)

    /// access to SSL Rx Gain methods
    SSLProfileStatic::Pointer m_pSSL{SSLProfileStatic::create()};


    //  ------------------------------------------------------------------
    //  Declare additional member functions of the sequence
    //  ------------------------------------------------------------------
    bool fEPIStdInit(SeqLim& rSeqLim, SeqBuildBlockEPIReadOut* _pEPIRO, ReorderInfo* _pREOInfo);

    //  --------------------------------------------------------------
    /// \brief <b> UI class for Ep2d
    ///
    ///         This class is basically empty on the measurement system
    //  --------------------------------------------------------------
    EpUI* m_pUI{nullptr};

    //  ------------------------------------------------------------------
    //
    //  Name        :  Ep2d::createUI
    //
    //  Description :
    /// \brief <b>     Instantiation of UI classes   </b>
    //
    //  Return      :  NLS status
    //
    //  ------------------------------------------------------------------
    virtual NLS_STATUS createUI(SeqLim& rSeqLim);

#ifdef COMPILE_EP2D_DIFF
    virtual bool updateCompGrad(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo);
#endif

  protected:
    long                                         m_lSliceAccelRefScans{0}; // Number of slice acceleration prepare scans
    long                                         m_lSliceAccelDummyScans{0}; // Number of slice acceleration dummy scans
    long                                            m_lSliceAccelPhaseCorrScans{0};
    SequenceDebugSettings::SequenceDebugSettings m_debugSettings = SequenceDebugSettings::SequenceDebugSettings("USE_EPI_DEBUG_SETTINGS");

    bool m_isForcedFastGRERefscanWithoutPAT{false};

#ifdef EP2D_SE_MRE
    bool        m_bIsDefaultMode{true}; // 60Hz w/running trigger
    double      m_dMEGFrequency{MEG_DEFAULT_HZ}; ///< MRE wobbler frequency
    long        m_lMEGFractionalEnc{100};
    SBBRefocSE* m_pSBBRefocSE{nullptr};

  public:
    long   m_lMinTE_Fractional_us{0};
    long   m_lMinTE_NonFractional_us{0};
    long   m_lFractionalTE_LimitsRange_us{0};
    double getMinTE_Fractional_ms();
    double getMinTE_NonFractional_ms();
    double getFractionalTE_LimitsRange_ms();
#endif



  private:

    void prepSliceAccelPhaseOffcenter3D(MrProt& rMrProt, sSLICE_POS asSLC[]);

};

#ifdef EP2D_SE_MRE
inline double Ep2d::getMinTE_Fractional_ms()
{
    return m_lMinTE_Fractional_us / 1000.;
}

inline double Ep2d::getMinTE_NonFractional_ms()
{
    return m_lMinTE_NonFractional_us / 1000.;
}

inline double Ep2d::getFractionalTE_LimitsRange_ms()
{
    return m_lFractionalTE_LimitsRange_us / 1000.;
}
#endif

#ifdef WIN32

//  ----------------------------------------------------------------------
//
//  Name        :  getUI
//
//  Description :
/// \brief         Returns the pointer to the UI class
///
//
//  Return      :  EpCommonUI*
//
//  ----------------------------------------------------------------------
EpUI* getUI(MrUILinkBase* const pThis);

//  ----------------------------------------------------------------------
//
//  Name        :  getUI
//
//  Description :
/// \brief         Returns the pointer to the UI class
///
//
//  Return      :  EpCommonUI*
//
//  ----------------------------------------------------------------------
EpUI* getUI(MrMeasSrv::ISequence* const pSeq);

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
Ep2d* getSeq(MrUILinkBase* const pThis);

#endif // #ifdef WIN32

// logging error if not in binary search and return error status
NLSStatus prepareError(const SeqLim& rSeqLim, std::string const& sErrorText);

} // namespace SEQ_NAMESPACE

#endif
