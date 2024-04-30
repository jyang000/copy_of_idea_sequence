//----------------------------------------------------------------------------------
// <copyright file="a_ep2d_diff.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

//----------------------------------------------------------------------------------------------------
// This define causes the 'Save uncombined' switch to be removed.
// Reason: incompatible with IceProgramDiffusion2D (CHARM 327088)
//----------------------------------------------------------------------------------------------------
//#define EPI_DISABLE_SAVE_UNCOMBINED           now defined in makefile

//----------------------------------------------------------------------------------------------------
// This define activates the frequency feedback.
// EPI navigator data is used to estimate the scanner frequency drift
// (relative ot the first acquired volume). A realtime feedback mechanism
// is used to provide the sequence with the evaluated frequency offset:
// this information is used to update the scanner frequency (involves a
// temporal filtering).
//----------------------------------------------------------------------------------------------------
//#define EPI_SUPPORT_FREQ_FEEDBACK             now defined in makefile

#include "MrImaging/seq/a_ep2d.h"

#ifdef BUILD_WIPParameterTool
#include "MrImagingFW/WIPParameterTool/WIPParameterTool.h"
#endif

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"


namespace SEQ_NAMESPACE
{
class __IMP_EXP Ep2d_diff : public Ep2d
{
  public:
    // constructor
    Ep2d_diff();

    // destructor
    virtual ~Ep2d_diff() = default;

    // no copy and move constructor and assignment operator
    Ep2d_diff(const Ep2d_diff& right) = delete;
    Ep2d_diff& operator=(const Ep2d_diff& right) = delete;
    Ep2d_diff(Ep2d_diff&& right)                 = delete;
    Ep2d_diff& operator=(Ep2d_diff&& right) = delete;

    // init of the sequence
    virtual NLSStatus initialize(SeqLim& rSeqLim) override;

    // Check on protocol, if there are forbidden combinations of  parameters.
    // This speeds up the binary search, because the method is called very early in
    // the prepare process and therefore can save lots of lengthy calculations.
    // The method must not prepare any objects like pulses or SBBs.</b>
    virtual NLSStatus prePrepare(const MrProt& rMrProt, const SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    // preparation of the sequence during binary search and prior to sequence execution
    virtual NLSStatus prepare(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    // checking the sequence for gradient stimulation
    virtual NLSStatus check(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, SEQCheckMode* pSEQCheckMode) override;

    // execution of the sequence
    virtual NLSStatus run(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

#ifdef WIN32
    // conversion rules from protocols of previous versions
    NLSStatus convProt(const MrProt& rMrProtSrc, MrProt& rMrProtDst) override;

    // Exporting diffusion timing (small and large delta) to UI and checking for IVIM condition
    bool exportDiffusionTimingToUI(SeqLim& rSeqLim, MrProt& rMrProt) override;

    // set thermal balancing flag in the UI (only in derived Ep2d_diff)
    void setUIThermalBalancing() override;
#endif

    // set PF factors to reorder info
    void setPartialFourierToReorderInfo(MrProt& rMrProt, ReorderInfo* pReorderInfo) const override;
  
#ifdef BUILD_WIPParameterTool
    WPT_NAMESPACE::WIPParameterTool m_WIPParamTool{*this};
#endif

  protected:

    bool isSystemCompatibleWithFlavor() const override;

    // since gradients during the IR block (spoiling, slice selection) are small, the duration is considered for the
    // cooling time. different in diffusion and others.
    void considerIRBlockForImplicitCoolingPause(
        SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo, long lScanTimeSatsEtc, long lScanTimeBasic) override;

    void loadSliceAdjData(SeqLim& rSeqLim, MrProt& rMrProt) override;

    void dumpSliceAdjData(SeqLim& rSeqLim, MrProt& rMrProt) override;

    std::string getSequenceVariantText() const override;

    // Export DICOM contrast text "Diffusion" to SeqExpo
    void setDICOMAcquisitionContrast(SeqExpo& rSeqExpo) const override;

    // sequence variant-specific (here diffusion) hard limit settings in initialize
    void setVariantSpecificHardLimits(SeqLim& rSeqLim) override;

    // sequence variant-specific (here diffusion) exports to SeqExpo, including some parameter combination checks
    NLSStatus setVariantSpecificExports(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo) override;

    // setting initial dummy scans at the beginning, without considering PAT, SMS, ...
    void setInitialDummyScansBasic(MrProt& rMrProt) override;

    // modify initial dummy scans considering SMS
    void setInitialDummyScansSMS(MrProt& rMrProt) override;

    // calculate number of required prep scans, including dummy scans
    long calcRequiredPrepScans(MrProt& rMrProt) override;

    // check TE and increase it if necessary / required by protocol parameters (e.g. TOM)
    NLSStatus CheckAndAdaptTE(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo, long& lNeededTE) override;

    // diffusion adjustment scan settings in RunKernel
    void setDiffusionAdjustmentScans(MrProt& rMrProt, long lSlice, long lKernelMode) override;

    // setting diffusion loop counters in RunKernel
    bool setDiffusionLoopCounters(long lSlice, long lKernelMode) override;

    // configure EPI kernel depending protocol parameters and sequence variant
    NLSStatus configureEPIKernel(SeqLim& rSeqLim, MrProt& rMrProt) override;

    // setting variant-specific (here diffusion) settings in the loop structure
    void setVariantSpecificLoopSettings(SeqLim& rSeqLim, MrProt& rMrProt) override;

    // configure protocol-independent settings in EPI kernel in initialize()
    bool initializeEPIKernel() override;

    // configure settings of the ReorderInfo structure
    void configureReorderInfo(MrProt& rMrProt) override;

    // setting external trigger in EPI kernel, and Osc bit
    void setTriggerAndOscBit(long lKernelMode) override;

    // setting the spoilers before and after SPAIR pulse
    void setSPAIRSpoilingType() override;

    // prepare bookkeeping for Dynamic Field Correction
    void prepareDFCBookkeeping() override;

    // Reduced inversion thickness for short inversion times (e.g. STIR) - avoid crosstalk
    // for nested inversion schemes.
    void setReducedIRThicknessForShortTI() override;

    // calculate the energy of an extra chemical saturation pulse to the RF block
    MrProtocolData::SeqExpoRFInfo calcEnergyOfExtraCSat(long lNumberOfPulses) override;

    // If B0 correction is enabled, complete volumes are acquired innermost.
    // Each volume is handled internally as a separate repetition => set
    // LastScanInMeas flag correspondingly. Exception: PatRefScans.
    void setLastScanInMeasFlagForB0Correction(MrProt& rMrProt, long lSlice, long lShot) override;

    // disabling diffusion gradients for PAT and SMS reference scans
    void disableDiffusionForPrepScan() override;

    // set last scan in meas when all slices have been acquired for one specific repetition counter
    // this is required if dynamic field correction is used with multiple concats
    NLS_STATUS setLastScanInMeasFlagForDFCBookkeeping(MrProt& rMrProt, long lKernelMode) override;

    // set some diffusion-specific SeqUT settings in Run (only in derived Ep2d_diff)
    void configureDiffusionSpecificSeqUTSettings(SeqLim& rSeqLim, MrProt& rMrProt) override;

    // disable readouts for preparation scans
    void disableReadoutForPrepScans(long lKernelMode, long lSlice=0) override;

    // check acq window when resp trigger is used
    bool checkAcqWindowForRespTriggering(MrProt& rMrProt) const override;

    // calculate otal number of expected volumes for frequency feedback, including preparation and adjustment scans
    long calcTotalNumberOfVolumesForFreqFeedback(MrProt& rMrProt) override;

    // get diffusion adjustment prep scan number in Ep2d_diff, returns 0 in other flavors.
    long getDiffusionAdjPrepScans() const override;

    // check the combination "dynamic field correction with multiple concats and respiratory compensation"
    // (only in derived Ep2d_diff)
    bool checkDFCWithLongTRTrigMode(SeqLim& rSeqLim, MrProt& rMrProt) const override;

    // set Zoomit-specific expectations for SeqUT
    void setZoomItSeqUTExpectations(MrProt& rMrProt) override;

    // check restrictions on multiple concatenations
    bool isMultiConcatsAllowed(MrProt& rMrProt) const override;

    // check whether the loop structure allows the complete volumes necessary for B0 correction
    bool isLoopStructureCompatibleWithB0Correction(MrProt& rMrProt) const override;

    // check if segmented PAT reference lines should be used.
    bool isSegmentedPATRefLinesCondition(MrProt& rMrProt) const override;

    // set long TR triggering mode depending on protocol parameters
    void setLongTRTrigMode(MrProt& rMrProt) override;

    // Decide whether an additional CSat is played out at the beginning of the cool pause.
    bool isExtraCSatApplied(MrProt& rMrProt) override;

    // check whether cooling time is executed with the given IR scheme in the protocol. (always true in non-diffusion)
    bool isCoolTimeExecutedWithIR(MrProt& rMrProt) const override;

    // set up compensation related parameters
    NLSStatus setCompensationPara(MrProt& rMrProt, SeqLim& rSeqLim);

  private:

    // thermal balancing flag
    bool m_bThermalBalancing{false};

    // maxwell correction flag (correcting the EPI readout)
    bool m_bMaxwellCorrection{false};

    /// this parameter handles the the setting of the external trigger during the sequence
    /** If the parameter is set, external trigger pulses are active. This can be handled
    controlled with an ini file.
    Default value is 'false', meaning that the external trigger pulses are not played out
    during the standard diffusion measurements.*/
    bool m_bSetExternalTrigger{false};

    /// Number of adjustment scans for dynamic distortion correction
    /** Adjustment scans are played out with the preparation scans. Value is provided by
    the diffusion module and includes averages.
    */
    long m_lAdjPrepScans{0};

    //  Length of the Array is equal to number of repetitions + m_lAdjPrepScans
    //  Each bit corresponds to one slice
    //  Used to set Last Scan in meas flag if all slices of a particular repetition have been acquired.
    std::vector<std::pair<uint64_t, uint64_t>> m_asAcqInRep;

    // index used to distinguish calls for check() for different diffusion directions
    int m_iCheckIndex{0};
};

} // namespace SEQ_NAMESPACE
