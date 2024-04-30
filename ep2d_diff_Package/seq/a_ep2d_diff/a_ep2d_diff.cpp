//----------------------------------------------------------------------------------
// <copyright file="a_ep2d_diff.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#include "a_ep2d_diff.h"
#include "MrProtSrv/Domain/CoreNative/MrNavigator.h"
#include "MrImaging/SequenceLibraries/libPace/PACE.h"
#include "MrGlobalDefinitions/MrResult.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Physiology/MrPhysiology.h"
#include "MrNFramework/MrFeatureToggle/Toggle.h"

#ifdef BUILD_WIPParameterTool
// Eddy current compensation
#include "MrImaging/seq/common/Common_TSE/EddyCurrentComp_UI.h"
#include "MrProtSrv/Domain/CoreNative/MrWipMemBlock.h"
#endif

#ifdef SEQUENCE_CLASS_EP2D
SEQIF_DEFINE(SEQ_NAMESPACE::Ep2d_diff)
#endif

using namespace SEQ_NAMESPACE;


Ep2d_diff::Ep2d_diff() : Ep2d()
{
    // ---------------------------------------------------------------------------
    // check diffusion thermal balancing feature toggle
    // ---------------------------------------------------------------------------
    try
    {
        const MrFeatureToggle::Toggle toggle("DiffThermalBalancing");
        m_bThermalBalancing = toggle.isOn();
    }
    // default: feature is off
    catch (const std::exception&)
    {
        m_bThermalBalancing = false;
    }

    // ---------------------------------------------------------------------------
    // set maxwell correction from sequence settings xml
    // ---------------------------------------------------------------------------
    m_bMaxwellCorrection = m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/enable_maxwell_correction_readout", true);

    SEQ_TRACE_DEBUG.print("Maxwell correction enabled: %d", m_bMaxwellCorrection);
}


NLSStatus Ep2d_diff::initialize(SeqLim& rSeqLim)
{
    // ----------------------------------------------------------------------
    // for diffusion, set the flag of including prep scans to
    // measurement-time-relevant scans, to decrease the delay in scan time
    // display in respiratory triggered diffusion measurements
    // ----------------------------------------------------------------------
    m_mySeqLoop.setPrepScansRelevantForMeasTime(true);

    // call base class initialization
    NLSStatus lStatus = Ep2d::initialize(rSeqLim);

    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
        return lStatus;

    if (SysProperties::isLowField())
    {
        rSeqLim.setFatSupOpt(MrProtocolData::FATSUPOPT_DEFAULT, MrProtocolData::FATSUPOPT_ABDOMEN);
    }

#ifdef BUILD_WIPParameterTool
    const bool bSeqSpecialCard
        = SysProperties::ReadSeqSettingGeneral("EDDYCURRCOMP_WIP_PARAMS_SPECIALCARD_VISIBLE", false, true);
    if (bSeqSpecialCard)
    {
        rSeqLim.setEddyCurrentComp(SEQ::OFF, SEQ::ON);
    }
#ifdef WIN32
    if (!initUIEddyCurrentComp(rSeqLim, m_WIPParamTool))
    {
        SEQ_TRACE_ERROR.print("Error at %s(%d).", __FILE__, __LINE__);
        return MRI_SEQ_SEQU_ERROR;
    }
#endif // WIN32
#endif // BUILD_WIPParameterTool

    return lStatus;
}

NLSStatus Ep2d_diff::prePrepare(const MrProt &rMrProt, const SeqLim &rSeqLim, SeqExpo &rSeqExpo)
{
    // call base class prePrepare
    NLSStatus lStatus = Ep2d::prePrepare(rMrProt, rSeqLim, rSeqExpo);

    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
        return lStatus;

    MrProtFacade protFacade(rMrProt);

    FatWaterContrast eFatSuppression  = rMrProt.preparationPulses().getlFatWaterContrast();
    FatSupOptRegion  eFatSupOptRegion = rMrProt.preparationPulses().getlFatSupOpt();

    bool bCombinationsAllowed    = true;
    bool bCompensationNotAllowed = false;

    bool bContextNormal = rSeqLim.isContextNormal();

    // only provide fat sat option in low field
    // use abdomen mode to activate the eddy current compensation
    if (SysProperties::isLowField())
    {
        bCombinationsAllowed
            = ((eFatSuppression == MrProtocolData::FatWaterContrast_FatSaturation
                && eFatSupOptRegion == MrProtocolData::FATSUPOPT_DEFAULT)
               || (eFatSuppression == MrProtocolData::FatWaterContrast_FatSaturation
                   && eFatSupOptRegion == MrProtocolData::FATSUPOPT_ABDOMEN)
               || (eFatSuppression == MrProtocolData::FatWaterContrast_Spair
                   && eFatSupOptRegion == MrProtocolData::FATSUPOPT_DEFAULT)
               || (eFatSuppression == MrProtocolData::FatWaterContrast_Spair
                   && eFatSupOptRegion == MrProtocolData::FATSUPOPT_ABDOMEN)
               || (eFatSuppression == MrProtocolData::FatWaterContrast_WaterExcitation
                   && eFatSupOptRegion == MrProtocolData::FATSUPOPT_DEFAULT)
               || (eFatSuppression == MrProtocolData::FatWaterContrast_WaterExcitation
                   && eFatSupOptRegion == MrProtocolData::FATSUPOPT_ABDOMEN)
               || (eFatSuppression == MrProtocolData::FatWaterContrast_Standard
                   && eFatSupOptRegion == MrProtocolData::FATSUPOPT_DEFAULT));

        bCompensationNotAllowed
            = ((rMrProt.diffusion().getdsScheme() != SEQ::DIFFSCHEME_MONOPOLAR)
               && ((eFatSuppression == MrProtocolData::FatWaterContrast_FatSaturation
                    && eFatSupOptRegion == MrProtocolData::FATSUPOPT_ABDOMEN)
                   || (eFatSuppression == MrProtocolData::FatWaterContrast_Spair
                       && eFatSupOptRegion == MrProtocolData::FATSUPOPT_ABDOMEN)
                   || (eFatSuppression == MrProtocolData::FatWaterContrast_FastWaterExcitation
                       && eFatSupOptRegion == MrProtocolData::FATSUPOPT_ABDOMEN)));

        if (!bCombinationsAllowed)
        {
            if (bContextNormal)
            {
                SEQ_TRACE_ERROR.print("ERROR: FatSat Mode currently not supported.");
            }
            return MRI_SEQ_SEQU_ERROR;
        }

        // ---------------------------------------------------------------------------
        // EC compensation is only for monopolar diffusion scheme
        // ---------------------------------------------------------------------------
        if (bCompensationNotAllowed)
        {
            if (bContextNormal)
            {
                SEQ_TRACE_ERROR.print("ERROR: EC compensation is only for monopolar diffusion scheme");
            }
            return MRI_SEQ_SEQU_ERROR;
        }
    }
    else
    {
#ifdef BUILD_WIPParameterTool
        // ---------------------------------------------------------------------------
        // EC compensation is only for monopolar diffusion scheme
        // ---------------------------------------------------------------------------
        if (rMrProt.getucEddyCurrentComp() && (rMrProt.diffusion().getdsScheme() != SEQ::DIFFSCHEME_MONOPOLAR))
            return prepareError(rSeqLim, "EC compensation is only for monopolar diffusion scheme");
#endif
    }

    // -------------------------------------------------------------------------------
    // Prohibit use of dynamic distortion correction in diffusion mode one_scan_trace
    // -------------------------------------------------------------------------------
    if (rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_ONE_SCAN_TRACE
        && rMrProt.getsDynDistCorrFilter().getucMode() != SEQ::DYN_DISTCORR_NONE)
        return prepareError(
            rSeqLim, "Dynamic distortion correction and ONE_SCAN_TRACE diffusion mode are incompatible");

    // -------------------------------------------------------------------------------
    // Do not allow distortion correction for DTI
    // -------------------------------------------------------------------------------
    if (protFacade.iSDTI())
    {
        if (rMrProt.getsDistortionCorrFilter().getucMode() != SEQ::DISTCORR_NDIS)
            return prepareError(rSeqLim, "DTI can only be measured without distortion correction");
    }

    // ---------------------------------------------------------------------------
    // Trigger solve handler to show dvs import / export error message
    // ---------------------------------------------------------------------------
#ifdef WIN32
    if (m_pUI->m_bImportExportError)
        return prepareError(rSeqLim, "error in diffusion direction import");
#endif

    if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
    {
        // -------------------------------------------------------------------------------
        // Prohibit use of PTX acceleration without bipolar diffusion or gradient reversal
        // -------------------------------------------------------------------------------
        if ((rMrProt.diffusion().getdsScheme() != SEQ::DIFFSCHEME_BIPOLAR || !protFacade.isGradientReversalDiffusion())
            && rMrProt.getsTXSPEC().getaPTXRFPulse()[0].getdPulseAcceleration() > 1.0)
            return prepareError(
                rSeqLim, "PTX acceleration is available only with bipolar diffusion or gradient reversal");

        // ---------------------------------------------------------------------------
        // Prohibit use of IR
        // ---------------------------------------------------------------------------
        if (rMrProt.getsPrepPulses().getucInversion() != SEQ::INVERSION_OFF)
            return prepareError(rSeqLim, "IR is not supported for ZOOMit");

        // -------------------------------------------------------------------------------
        // Prohibit use of PTX acceleration when rotated trajectory is used
        // -------------------------------------------------------------------------------
        if (rMrProt.getsTXSPEC().getaPTXRFPulse()[0].getdPulseAcceleration() > 1.0)
            return prepareError(rSeqLim, "PTX acceleration is not available when excitation trajectory is rotated");

        // -------------------------------------------------------------------------------
        // Prohibit use of FatWaterContrast_WaterExcitation with ZOOMit option
        // -------------------------------------------------------------------------------
        if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_WaterExcitation)
            return prepareError(
                rSeqLim, "FatWaterContrast_WaterExcitation and EXCITATION_ZOOMED mode are incompatible");

        // -------------------------------------------------------------------------------
        // Prohibit use of GRE Ref Scan with ZOOMit option
        // -------------------------------------------------------------------------------
        if (isGreRefScanType(rMrProt))
            return prepareError(rSeqLim, "GRE PAT ref scan and EXCITATION_ZOOMED mode are incompatible");

        // -------------------------------------------------------------------------------
        // Prohibit selection of slice acceleration with ZOOMit option
        // -------------------------------------------------------------------------------
        if (rMrProt.PAT().getucPATMode() == SEQ::PAT_MODE_SLICE_ACCELERATION)
            return prepareError(rSeqLim, "PAT_MODE_SLICE_ACCELERATION and EXCITATION_ZOOMED mode are incompatible");
    }

    return MRI_SEQ_SEQU_NORMAL;
}

NLSStatus Ep2d_diff::prepare(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo)
{
    NLSStatus lStatus = MRI_SEQ_SEQU_NORMAL;

    lStatus = setCompensationPara(rMrProt, rSeqLim);
    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
        return lStatus;

    // call base class prepare
    return Ep2d::prepare(rMrProt, rSeqLim, rSeqExpo);
}

NLSStatus Ep2d_diff::check(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, SEQCheckMode* pSEQCheckMode)
{
    // call base class check. for diffusion, 2 different directions are checked.
    NLSStatus lStatus = MRI_SEQ_SEQU_NORMAL;

    for (m_iCheckIndex = 0; m_iCheckIndex < 2; ++m_iCheckIndex)
    {
         lStatus = Ep2d::check(rMrProt, rSeqLim, rSeqExpo, pSEQCheckMode);

         if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
             return lStatus;
    }

    return lStatus;
}

NLSStatus Ep2d_diff::run(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo)
{
    // call base class run
    return Ep2d::run(rMrProt, rSeqLim, rSeqExpo);
}

#ifdef WIN32
NLSStatus Ep2d_diff::convProt(const MrProt& rMrProtSrc, MrProt& rMrProtDst)
{
    // Charm MR_00391990:
    // Set fat saturation mode depending on RF pulse type if source
    // protocol version is earlier than VD11A.
    // RF pulse type 'low SAR' activates fat saturation mode 'strong',
    // otherwise this mode is set to 'weak'.
    if (rMrProtSrc.getConvFromVersion() < 41110000)
    {
        if (rMrProtSrc.getsTXSPEC().getucRFPulseType() == SEQ::RF_LOW_SAR)
        {
            rMrProtDst.getsPrepPulses().setucFatSatMode(SEQ::FAT_SAT_STRONG);
        }
        else
        {
            rMrProtDst.getsPrepPulses().setucFatSatMode(SEQ::FAT_SAT_WEAK);
        }
    }
    if (rMrProtSrc.getConvFromVersion() <= 41310000)
    {
        //  Early VD13A version or older
        const int32_t i32bValueSize = rMrProtSrc.getsDiffusion().getlDiffWeightings();
        int32_t       i32Pos        = 0;
        for (; i32Pos != i32bValueSize; ++i32Pos)
        {
            if (rMrProtSrc.getsDiffusion().getalAverages()[i32Pos] < 1)
            {
                rMrProtDst.getsDiffusion().getalAverages()[i32Pos] = 1;
            }
        }

        if ((rMrProtSrc.getsNavigatorPara().getlRespComp() != SEQ::RESP_COMP_OFF) && (rMrProtSrc.getlAverages() > 1))
        {
            //  PACE can not use repetition loop for multiple averages since PACE usually requires multiple
            //  concatenations. To reduce code complexity the average loop is not supported any more and the b-value
            //  specific #average mechanism is used.
            const int32_t i32NAvg = rMrProtSrc.getlAverages();
            for (i32Pos = 0; i32Pos != i32bValueSize; ++i32Pos)
            {
                rMrProtDst.getsDiffusion().getalAverages()[i32Pos] *= i32NAvg;
            }
            rMrProtDst.setlAverages(1);
        }
    }

    // CHARM 472775 in VE11C and previous versions the diffusion averages got multiplied by the number of repetitions.
    // In VA10A the diffusion averages are directly written to the protocol and are not multiplied by the number of
    // repetitions.
    if (rMrProtSrc.getConvFromVersion() <= 51130001) // VE11C
    {
        const int32_t i32bValueSize = rMrProtSrc.getsDiffusion().getlDiffWeightings();
        int32_t       i32Pos        = 0;
        for (; i32Pos != i32bValueSize; ++i32Pos)
        {
            rMrProtDst.getsDiffusion().getalAverages()[i32Pos]
                = rMrProtDst.getsDiffusion().getalAverages()[i32Pos] * (rMrProtSrc.getlRepetitions() + 1);
        }
    }

    const NLS_STATUS u64Status = PACE::fConvProt(*rMrProtSrc.getProtData(), *rMrProtDst.getProtData());
    if (NLS_SEVERITY(u64Status) != NLS_SUCCESS)
    {
        SEQ_TRACE_ERROR.print("Error at %s(%d).", __FILE__, __LINE__);
        return u64Status;
    }

    return MRI_SEQ_SEQU_NORMAL;
}
#endif

void Ep2d_diff::loadSliceAdjData(SeqLim& rSeqLim, MrProt& rMrProt)
{
#ifdef WIN32
    if (rSeqLim.isContextNormal())
    {
        if (m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/load_sliceAdjust_data", false))
        {
            string strFileName = getenv("CustomerSeq");
            strFileName.append("\\SliceAdjDump.txt");
            SLICEADJ::LoadSliceAdjData(rMrProt, strFileName.c_str());
        }
    }
#endif // WIN32
}

void Ep2d_diff::dumpSliceAdjData(SeqLim& rSeqLim, MrProt& rMrProt)
{
#ifdef WIN32
    if (rSeqLim.isContextPrepForMeasurement())
    {
        if (m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_final_sliceAdjust_data", false))
        {
            string strFileName = getenv("CustomerSeq");
            strFileName.append("\\SliceAdjDump.txt");
            SLICEADJ::DumpSliceAdjData(rMrProt, strFileName.c_str(), true);
        }
    }

    if (rSeqLim.isContextPrepForCuboidCalculation())
    {
        if (m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_default_sliceAdjust_data", false))
        {
            string strFileName = getenv("CustomerSeq");
            strFileName.append("\\SliceAdjDumpDefault.txt");
            SLICEADJ::DumpSliceAdjData(rMrProt, strFileName.c_str(), true);
        }
    }
#endif // WIN32
}

#ifdef WIN32
bool Ep2d_diff::exportDiffusionTimingToUI(SeqLim& rSeqLim, MrProt& rMrProt)
{
    MrProtFacade protFacade(rMrProt);

    // ---------------------------------------------------------------------------
    // check IVIM condition
    // ---------------------------------------------------------------------------
    if (protFacade.isIVIM())
    {
        // check if smallest b value (>0) in protocol is larger than the smallest possible b value
        // calculated by the diffusion SBB
        const long lSmallestPossibleBValue
            = m_EPIKernel.getSmallestIVIMbValuePossible(rSeqLim.isContextPrepForBinarySearch());

        for (long lBValue = 0; lBValue < rMrProt.getsDiffusion().getlDiffWeightings(); lBValue++)
        {
            if ((rMrProt.getsDiffusion().getalBValue()[lBValue] > 0)
                && (rMrProt.getsDiffusion().getalBValue()[lBValue] < 50))
            {
                if (lSmallestPossibleBValue > rMrProt.getsDiffusion().getalBValue()[lBValue])
                {
                    if (!rSeqLim.isContextPrepForBinarySearch())
                    {
                        SEQ_TRACE_ERROR.print(
                            "ERROR: IVIM b value set in protocol (%d) but cannot be realized. Smallest possible value: %ld",
                            rMrProt.getsDiffusion().getalBValue()[lBValue],
                            lSmallestPossibleBValue);
                    }
                    return false;
                }
            }
        }
    }

    // get pointer to diffusion gradient object and make it available to UI functions
    m_pUI->fUILinkRegisterDidi(m_EPIKernel.getDidiPointer());

    // Get diffusion gradient duration (little delta) and make it available to UI functions
    long lIndex                    = m_pUI->ToolTipParamDiffGradDuration;
    m_pUI->m_dToolTipParam[lIndex] = m_EPIKernel.getDiffGradDuration_ms();

    // Get diffusion gradient spacing (big delta) and make it available to UI functions
    lIndex                         = m_pUI->ToolTipParamDiffGradSpacing;
    m_pUI->m_dToolTipParam[lIndex] = m_EPIKernel.getDiffGradSpacing_ms();

    // Assemble diffusion vector set tooltip and make it available to UI functions
    lIndex = m_pUI->ToolTipStringDVSInfo;
    if (m_EPIKernel.getDidiPointer()->isDirectionInternal(m_EPIKernel.getDidiPointer()->getNumberOfDirections()))
    {
        // Tooltip for Siemens internal vector set (does not include filename)
        m_pUI->m_sToolTipString[lIndex] = m_EPIKernel.getDidiPointer()->getComment();
    }
    else
    {
        // Tooltip for user defined vector sets (including filename)
        m_pUI->m_sToolTipString[lIndex]
            = m_EPIKernel.getDidiPointer()->getVectorFileName() + "\n" + m_EPIKernel.getDidiPointer()->getComment();
    }

    return true;
}
#endif

void Ep2d_diff::setDICOMAcquisitionContrast(SeqExpo& rSeqExpo) const
{
    rSeqExpo.setDICOMAcquisitionContrast("DIFFUSION"); // Acquisition Contrast (0008,9209)
}

void Ep2d_diff::setVariantSpecificHardLimits(SeqLim& rSeqLim)
{
    rSeqLim.setEPIFactor(1, 512, SEQ::INC_SINGLESHOT, 128);
    
    //set up SMS-related limits
    rSeqLim.setPATMode(SEQ::PAT_MODE_NONE, SEQ::PAT_MODE_GRAPPA, SEQ::PAT_MODE_SENSE, SEQ::PAT_MODE_SLICE_ACCELERATION);
    rSeqLim.setSliceAccMultiBandFactor(1, 4, 1, 1);
    rSeqLim.setSliceAccFOVShiftFactor(1, 4, 1, 1);

    rSeqLim.setRefScanMode(SEQ::PAT_REF_SCAN_EXTRA_EPI, SEQ::PAT_REF_SCAN_EXTRA, SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST);

    stringstream ssSequenceHintText;
    ssSequenceHintText << "Sequence: ep2d_diff: ";
#ifdef WIP
    ssSequenceHintText << "WIP sequence. Not for clinical use. ";
#endif
    ssSequenceHintText << "Application: EPI-based diffusion weighted imaging and DTI. Basics: 2D single shot EPI; "
                          "asymmetric sampling in phase direction; Diffusion encoding with four bipolar "
                          "gradients and double spin echo; pixelwise phase correction using three measured projection "
                          "lines; read-out module optimized for gradient performance; Build: "
                       << __DATE__ << "  " << __TIME__;
    rSeqLim.setSequenceHintText(ssSequenceHintText.str().c_str());

    rSeqLim.setPhasePartialFourierFactor(SEQ::PF_6_8, SEQ::PF_7_8, SEQ::PF_OFF, SEQ::PF_5_8);
    rSeqLim.setInversion(SEQ::INVERSION_OFF, SEQ::SLICE_SELECTIVE); /* FLAIR ;-) */
    rSeqLim.setIRScheme(SEQ::IR_SCHEME_AUTO, SEQ::IR_SCHEME_SEQUENTIAL);
    rSeqLim.setTI(0, 0, 9000000, 100, 2500000);
    rSeqLim.setRFPulseType(SEQ::RF_NORMAL, SEQ::RF_LOW_SAR);
    rSeqLim.setFlipAngle(90.0, 90.0, 1.0, 90.0);
    rSeqLim.setNoiseLevel(0, 4095, 1, 20);
    rSeqLim.setAdjFreProtRelated(SEQ::ON, SEQ::OFF);
    if (m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/turn_AdjFreProtRelated_off", false))
    {
        rSeqLim.setAdjFreProtRelated(SEQ::OFF);
    }

    rSeqLim.setSliceDistanceFactor(-0.500, 8.000, 0.010, 0.500);

    // There is a service version of the EPI diffusion sequence which requires that there is no
    // AjdFre which is enforced by the prot related flag. The following define is used by the
    // service sequence to compile a corresponding target.
#ifdef COMPILE_AS_SERVICE_SEQUENCE
    rSeqLim.setAdjFreProtRelated(SEQ::OFF);
#endif

    rSeqLim.setFatWaterContrast(
        FatWaterContrast_FatSaturation,
        FatWaterContrast_WaterExcitation,
        FatWaterContrast_Spair,
        FatWaterContrast_Standard);
    // Allow two modes of fat saturation
    //  weak:   standard
    //  strong: additional gradient reversal (excitation vs. refocusing slice selection gradient)
    rSeqLim.setFatSatMode(SEQ::FAT_SAT_WEAK, SEQ::FAT_SAT_STRONG);

    // 20150105 DP: allow "abdomen" mode for SPAIR for 3T systems
    if (SysProperties::getNominalB0() > 2.5)
    {
        rSeqLim.setFatSupOpt(
            MrProtocolData::FATSUPOPT_DEFAULT,
            MrProtocolData::FATSUPOPT_BRAIN,
            MrProtocolData::FATSUPOPT_ABDOMEN,
            MrProtocolData::FATSUPOPT_THORAX,
            MrProtocolData::FATSUPOPT_BREAST);
    }

    // Default combination mode: adaptive combine
    rSeqLim.setCoilCombineMode(SEQ::COILCOMBINE_ADAPTIVE_COMBINE, SEQ::COILCOMBINE_SUM_OF_SQUARES);

    // Allow Noise Mask
    rSeqLim.setFilterPrescanNormalizeNoiseMask(SEQ::OFF, SEQ::ON);

    // Allow automatic optimization of TE
    rSeqLim.setTOM(SEQ::TOM_OFF, SEQ::TOM_MINIMIZE_TE);

    // Allow dynamic distortion correction
    rSeqLim.setFilterType(
        SEQ::FILTER_NONE,
        SEQ::FILTER_RAW,
        SEQ::LARGE_FOV,
        /* SEQ::NORMALIZE, */ SEQ::ELLIPTICAL,
        /* SEQ::FILTER_IMAGE, */ SEQ::PRESCAN_NORMALIZE,
        SEQ::DYNAMIC_DIST_CORR /* , SEQ::FILTER_BIFIC */);
    rSeqLim.setDynamicDistortionCorrMode(SEQ::DYN_DISTCORR_NONE, SEQ::DYN_DISTCORR_ADJ, SEQ::DYN_DISTCORR_DIRECT);

    // modify maximum TE and TR after introduction of modified PAT selection solve-handler (CHARM 315966)
    // previous values did not allow solve-handler to work correctly with large number of slices
    rSeqLim.setTE(0, 1000, 400000, 100, 200000);
    rSeqLim.setTR(0, 10000, 30000000, 1000, 500000);

    // Do not display the following UI parameters:
    rSeqLim.getFlipAngle().setDisplayMode(SEQ::DM_OFF);
    // rSeqLim.setMultipleSeriesMode(SEQ::MULTIPLE_SERIES_EACH_MEASUREMENT);
    // rSeqLim.getMultipleSeriesMode().setDisplayMode(SEQ::DM_OFF);
    // check whether debug settings are used for the sequence and if so

    // mark sequence as "user sequence"
    if (m_debugSettings.areDebugSettingsUsed())
    {
        SEQ_TRACE_WARN << "This sequences uses debug settings and its behavior might deviate from the product "
                          "sequence! It is not released for clinical use!";
        rSeqLim.setSequenceOwner("USER");
    }

    // Initialize SBBDiffusion

    m_EPIKernel.getDiffusionSeqLims(rSeqLim);

    //  Set by the sequence and not by SBBDiffusion (like the other diffusion parameters)
    //  since the feature is implemented by the sequence and not by SBBDiffusion.
    rSeqLim.setDiffAverages(1, 64, 1, 1); // Local averages  (individual for each b-value)
    rSeqLim.getAverages().setDisplayMode(SEQ::DM_OFF);

    rSeqLim.setAdjustmentMode(AdjustmentMode_Standard, AdjustmentMode_SliceBySlice, AdjustmentMode_FastView);
    rSeqLim.setAdjSliceBySliceFirstOrderShim(SEQ::OFF, SEQ::ON);
    rSeqLim.setAdjSliceBySliceFrequency(SEQ::OFF, SEQ::ON);
    rSeqLim.setAdjSliceBySliceTxRef(SEQ::OFF, SEQ::ON);
    if (SysProperties::isPTxSystem())
        rSeqLim.setAdjSliceBySlicePtx(SEQ::OFF, SEQ::ON);
    else
        rSeqLim.setAdjSliceBySlicePtx(SEQ::OFF);

    if (!SysProperties::isUHFSystem())
    {
        rSeqLim.getAdjShim().setBitMask(rSeqLim.getAdjShim().getBitMask() | SEQ::ADJSHIM_WHOLE_BODY);
    }

    // file containing the default post processing protocol (EVAProtocol)
#ifdef WIN32
    rSeqLim.setDefaultEVAProt(_T("%SiemensEvaDefProt%\\DTI\\DTI.evp"));
#endif

    rSeqLim.enableCaptureCycle();
    rSeqLim.setConcatenationsSelectModeResp(SEQ::CONCAT_MANUAL, SEQ::CONCAT_AUTOMATIC);
    rSeqLim.setAcquisitionWindowSelectModeResp(SEQ::ACQUISITION_WINDOW_MS, SEQ::ACQUISITION_WINDOW_PERCENT);
    rSeqLim.setAcqusitionWindowPercentResp(1, 100, 1, 25);

    rSeqLim.setPhaseCorrectionMode(MrProtocolData::PHASECORR_INTERNAL, MrProtocolData::PHASECORR_EXTERNAL);

    // call base class implementation for common part
    Ep2d::setVariantSpecificHardLimits(rSeqLim);
}

NLSStatus Ep2d_diff::setVariantSpecificExports(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo)
{
    MrProtFacade protFacade(rMrProt);

    rSeqExpo.setAdditionalScaleFactor(1);
    rSeqExpo.setAdditionalRawDataSizeFactor(1.); // Default: No additional raw data size scale factor

    // Enable dedicated adaptive combine mode for diffusion imaging (also affects calculation of
    // raw data object size, see Sequence::checkRawDataSizeOfMeasurement).
    // Note: the mode is activated explicitely within IceDiffusion / IceDti configurator
    rSeqExpo.setICEProgramParam(
        ICE_PROGRAM_PARA_CTRL_MASK,
        rSeqExpo.getICEProgramParam(ICE_PROGRAM_PARA_CTRL_MASK) | ICE_PROGRAM_MSK_AC_USEFIRSTMAPONLY);

    // Export number of adjustment prep scans
    rSeqExpo.setAdjScans(static_cast<int32_t>(m_lAdjPrepScans));

    // Provide Ice with total number of repetitions
    const long lTotalNumberOfRepetitions = m_EPIKernel.getTotalScans() * rMrProt.measurements() + m_lAdjPrepScans;
    rSeqExpo.setICEProgramParam(ICE_PROGRAM_DIFF_REPETITIONS, lTotalNumberOfRepetitions);

    if (protFacade.iSDTI())
    {
        // Prohibit simultaneous use of respiration compensation (PACE) and tensor evaluation
        if (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF)
        {
            if (rSeqLim.isContextNormal())
            {
                SEQ_TRACE_ALWAYS.print("Error at %s(%d)", __FILE__, __LINE__);
            }
            return MRI_SEQ_SEQU_ERROR;
        }

        // Prohibit simultaneous use of multiple slice groups (MSMA) and tensor evaluation (Mosaic!)
        if (rMrProt.sliceGroupList().size() > 1)
        {
            if (rSeqLim.isContextNormal())
            {
                SEQ_TRACE_ALWAYS.print("Error at %s(%d)", __FILE__, __LINE__);
            }
            return MRI_SEQ_SEQU_ERROR;
        }

        rSeqExpo.setICEProgramFilename("%SiemensIceProgs%\\IceProgramDti2D");
        rSeqExpo.setOnlineFFT(SEQ::ONLINE_FFT_PHASE);
        rSeqExpo.setICEProgramParam(ICE_PROGRAM_PARA_SHOW_OFFLINE, SEQ::SO_SHOW_NO);

        // Additional settings for raw data object size calculation (see Sequence::checkRawDataSizeOfMeasurement)
        // ======================================================================================================
        // Averages, diffusion directions and b-values are handled as repetitions in the Dti Ice world.
        // The raw object memory allocated for each repetition gets deallocated once it is no
        // longer needed. By default, the raw data size estimation considers a data amount
        // corresponding to the number of measurements (protocol value). Here, this gets scaled down
        // to the number of volumes that actually have to be kept in memory.
        //
        // - If adjustment scans are acquired, we assume that all adjustment scans and an empirical maximum
        //   number of additional imaging scans have to be kept in memory simultaneously.
        // - If no adjustment scans are acquired, we assume that 4 imaging scans have to be kept in memory
        //   simultaneously.
        double dVolumesInMemory = 4.;

        if (m_lAdjPrepScans > 0)
        {
            dVolumesInMemory = static_cast<double>(m_lAdjPrepScans + 20);
        }

        rSeqExpo.setAdditionalRawDataSizeFactor(
            std::min(1., dVolumesInMemory / static_cast<double>(lTotalNumberOfRepetitions)));
    }
    else
    {
        if ((rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF)
            && (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_BREATH_HOLD))
        {
            rSeqExpo.setICEProgramFilename("%SiemensIceProgs%\\IceProgramDiffusion2D+PACE");
        }
        else
        {
            rSeqExpo.setICEProgramFilename("%SiemensIceProgs%\\IceProgramDiffusion2D");
        }

        // Additional settings for raw data object size calculation (see Sequence::checkRawDataSizeOfMeasurement)
        // ======================================================================================================
        // Averages, diffusion directions and b-values are handled as repetitions in the Diffusion Ice world:
        // - If PACE is disabled, the raw object memory allocated for each repetition gets deallocated once it is no
        //   longer needed. Thus, it is not neccessary to consider averages for the raw data object size
        //   calculation.
        // - If PACE is enabled, the raw object includes the repetitions dimension and averages have to be
        // considered. Additional adjustment scans have to be considered in any case.

        long lNumberOfScansForRawObj;

        if (rMrProt.NavigatorParam().getlRespComp() == SEQ::RESP_COMP_OFF)
        {
            // Get total number of scans excluding averages
            lNumberOfScansForRawObj = m_EPIKernel.getTotalScans(false) + m_lAdjPrepScans;
        }
        else
        {
            // Get total number of scans including averages
            lNumberOfScansForRawObj = lTotalNumberOfRepetitions;
        }

        rSeqExpo.setAdditionalRawDataSizeFactor(
            rSeqExpo.getAdditionalRawDataSizeFactor() * static_cast<double>(lNumberOfScansForRawObj)
            / static_cast<double>(rMrProt.measurements()));
    }

    // Insert additional SMS ICE program
    if (protFacade.isSliceAcceleration())
    {
        rSeqExpo.AddAdditionalIceProgramFileName("%SiemensIceProgs%\\IceProgramSMSAcc");
    }

    rSeqExpo.setPCAlgorithm(SEQ::PC_ALGORITHM_NONE);

    rSeqExpo.setApplicationCard(SEQ::APPLICATION_CARD_DIFF);
    rSeqExpo.setApplicationCardName(SEQ::APPLICATION_CARD_NAME_DIFF);

    if (rMrProt.getsTXSPEC().getucExcitMode() != SEQ::EXCITATION_ZOOMED)
    {
        if (rMrProt.preparationPulses().getucInversion() == SEQ::INVERSION_OFF)
        {
            fSUSetSequenceString(rMrProt, rSeqLim, rSeqExpo, "epse");
        }
        else
        {
            fSUSetSequenceString(rMrProt, rSeqLim, rSeqExpo, "epir");
        }
    }
    else
    {
        // for ZOOMit mark 2nd char with "z" (diffusion functor overwrites the sequence string from 3rd character)
        if (rMrProt.preparationPulses().getucInversion() == SEQ::INVERSION_OFF)
        {
            fSUSetSequenceString(rMrProt, rSeqLim, rSeqExpo, "ezse");
        }
        else
        {
            fSUSetSequenceString(rMrProt, rSeqLim, rSeqExpo, "ezir");
        }
    }

    // define last scan in meas as "trigger" for image looper
    if (protFacade.isBookkeepingConditionForDFC())
    {
        rSeqExpo.setOnlineFFT(SEQ::ONLINE_FFT_LASTSCANINMEAS);
    }

    return MRI_SEQ_SEQU_NORMAL;
}

void Ep2d_diff::setInitialDummyScansBasic(MrProt& rMrProt)
{
    if ((m_EPIKernel.getTotalScans(true) * (rMrProt.getlRepetitions() + 1)) > 1)
    {
        double prepScanTime_us = 3000000.0;
        if (isFastGreRefScan(rMrProt))
        {
            // Prepare at least for 2800ms => a protocol with TR 3000ms will need a single preparing scan only
            prepScanTime_us = 2800000.0;
        }
        // we want to prepare at least for 3 sec; i.e.
        m_lInitialDummyScans = static_cast<long>(prepScanTime_us / (rMrProt.tr()[0] * rMrProt.physiology().phases()) + 1.0);
    }
    else
    {
        // min. one prep scan for diffusion
        m_lInitialDummyScans = 1;
    }
}

 void Ep2d_diff::setInitialDummyScansSMS(MrProt& rMrProt)
{
    // An initial dummy scan is necessary for SMS acquisitions with an preceding GRE PAT ref scan but not for the EPI
    // and the fast GRE (skipping additional EPI) PAT ref scan
    if (rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI || isFastGreRefScan(rMrProt))
    {
        m_lInitialDummyScans = 0;
    }
    else
    {
        m_lInitialDummyScans = 1;
    }
    // Number of slice acceleration dummy scans must be at least 1 due to internal SeqLoop tracking
    m_lSliceAccelDummyScans = 1;
}

long Ep2d_diff::calcRequiredPrepScans(MrProt& rMrProt)
{
    // call base implementation
    long lRequiredPrepScans = Ep2d::calcRequiredPrepScans(rMrProt);

    // add diffusion-specific adjustment scans for DFC
    lRequiredPrepScans += m_lAdjPrepScans;

    return lRequiredPrepScans;
}

NLSStatus Ep2d_diff::CheckAndAdaptTE(SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo, long& lNeededTE)
{
    // call base implementation
    const NLSStatus lStatus = Ep2d::CheckAndAdaptTE(rSeqLim, rMrProt, rSeqExpo, lNeededTE);
    
    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
        return SeverePrepareErrorReturn(lStatus);

    // do diffusion-specific settings if TE was changed in the base class
    if (m_EPIKernel.getNeededTE() != rMrProt.te()[0])
    {
        // ---------------------------------------------------------------------------
        // Now we are in trouble: TE has been increased without preparing the kernel
        // again (EPIKernel.increaseTE() only updates some fill times). However, with
        // a different TE, the diffusion gradients will change and thus also the
        // dynamic TR fill time! If we would just continue without taking care of this,
        // the needed TR which will be calculated below gets inconsistent!
        //
        // That's how we handle this situation: create a copy of the original
        // protocol (we are not allowed to change the original one), set TE to the
        // new (increased) value and prepare the kernel (and thus the diffusion
        // module) with TOM disabled. The latter is required to ensure that the
        // diffusion module preparation behaves the same way it does in the
        // final prep (context normal).

        // ---------------------------------------------------------------------------
        // prepare SBBEPIKernel
        // ---------------------------------------------------------------------------
        MrProt NewProt(*rMrProt.clone());

        NewProt.te()[0] = static_cast<int>(lNeededTE);
        NewProt.TOM(SEQ::TOM_OFF);

        if (!m_EPIKernel.prep(NewProt, rSeqLim, rSeqExpo))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("%s: 0x%lx", "m_EPIKernel.prep", m_EPIKernel.getNLSStatus());
            }
            return SeverePrepareErrorReturn(m_EPIKernel.getNLSStatus());
        }
#ifdef WIN32
        m_pUI->fUILinkRegisterDidi(m_EPIKernel.getDidiPointer());
#endif
    }

    return MRI_SEQ_SEQU_NORMAL;
}

void Ep2d_diff::setDiffusionAdjustmentScans(MrProt& rMrProt, long lSlice, long lKernelMode)
{
    //---------------------------------------------------------------------------------------
    // Identify prep-scans, iPAT reference scans, internal adjustment scans and imaging scans
    //---------------------------------------------------------------------------------------
    int iAdjustmentScan = 0; // Default: tag as imaging scan

    m_EPIKernel.setDiffusionGradientsEnabled(true); // Default: apply diffusion gradients

    if (m_lAdjPrepScans > 0)
    {
        // Series of prep scans if dynamic distortion correction is enabled:
        // initial preps       : ignore
        // intermediate preps  : iPAT reference
        // final preps         : dynamic distortion correction adjustment

        long lFirstAdjScan = m_lInitialDummyScans + m_lPhaseCorrPrepScans;
        long lLastAdjScan  = lFirstAdjScan + m_lAdjPrepScans - 1;

        // If PAT is enabled, adjustment scans are acquired after
        // PAT reference scans.
        lFirstAdjScan += m_lPATRefScans;
        lLastAdjScan += m_lPATRefScans;

        MrProtFacade protFacade(rMrProt);

        if (protFacade.isSliceAcceleration())
        {
            lFirstAdjScan += (m_lSliceAccelRefScans + m_lSliceAccelDummyScans + m_lSliceAccelPhaseCorrScans);
            lLastAdjScan += (m_lSliceAccelRefScans + m_lSliceAccelDummyScans + m_lSliceAccelPhaseCorrScans);
        }

        if ((lKernelMode == KERNEL_PREPARE) && (m_alPrepScanCounter[lSlice] >= lFirstAdjScan)
            && (m_alPrepScanCounter[lSlice] <= lLastAdjScan))
        {
            // Enable readout
            fRTSetReadoutEnable(1);
            // Send data to SFC functor
            m_EPIKernel.setSendSliceAdjustData(true);
            // Set adjustment scan number (number 0 is reserved for imaging scans)
            iAdjustmentScan = static_cast<int>(m_alPrepScanCounter[lSlice] - lFirstAdjScan + 1);
        }
    }

    m_EPIKernel.setAdjustmentScan(iAdjustmentScan);
}

bool Ep2d_diff::setDiffusionLoopCounters(long lSlice, long lKernelMode)
{
    long lDiffLoopCounter;
    const long lRepLoopCounter  = m_mySeqLoop.getRepetitionLoopCounter();

    if ((lKernelMode & KERNEL_CHECK) == KERNEL_CHECK)
    {
        if (m_iCheckIndex == 0)
        {
            // Check the diffusion direction with the largest read component
            lDiffLoopCounter = m_EPIKernel.getDiffLoopCounterForHighestReadComponent(&m_asSLC[lSlice]);
        }
        else
        {
            // Check the last occurrence of the highest b-value scan
            lDiffLoopCounter = m_EPIKernel.getDiffLoopCounterForHighestBValue();
        }
    }
    else
    {
        lDiffLoopCounter = m_mySeqLoop.getDiffusionLoopCounter();
    }

    // Provide diffusion module with information about current scan:
    //   - Adjustment scan index (0 if it's a regular imaging scan)
    //   - Repetition loop counter (average)
    //   - Diffusion loop counter (runs over b-values and directions)
    //   - Pointer to ADC (used to modify Mdh entries)
    // Note: iPAT reference scans are not handled appropriately -
    //       repetitions loop counter will get patched below.
    // Note: Correct values are required before the EPI kernel is called -
    //       setting the values within the diffusion plugin is too late!
    if (!m_EPIKernel.setLoopCounters(m_EPIKernel.getAdjustmentScan(), lRepLoopCounter, lDiffLoopCounter))
    {
        SEQ_TRACE_ERROR.print("ERROR: Setting of Mdh entries failed");
        return false;
    }

    SEQ_TRACE_DEBUG.print(
        "Scan SeqLoop #%ld, Crep=%ld    (Kernelmode=%ld)", lDiffLoopCounter, lRepLoopCounter, lKernelMode);

    return true;
}

std::string Ep2d_diff::getSequenceVariantText() const
{
    return {"N/X EPI 2D DIFFUSION Sequence"};
}

NLSStatus Ep2d_diff::configureEPIKernel(SeqLim& rSeqLim, MrProt& rMrProt)
{
    // set thermal balancing flag
    m_EPIKernel.setThermalBalancing(m_bThermalBalancing);

    if (rMrProt.TOM() != SEQ::TOM_MINIMIZE_TE)
    {
        m_EPIKernel.setWantedTE(rMrProt.te()[0]);
    }
    else
    {
        m_EPIKernel.setWantedTE(0);
    }

    MrProtFacade protFacade(rMrProt);

    // Enable Maxwell correction if SMS is not selected
    if (!protFacade.isSliceAcceleration())
        m_EPIKernel.setMaxwellCorrSB_Slice(m_bMaxwellCorrection);
    else
        m_EPIKernel.setMaxwellCorrSB_Slice(false);


    // call base class implementation for common settings among sequence variants
   return Ep2d::configureEPIKernel(rSeqLim, rMrProt);
}

#ifdef WIN32
void Ep2d_diff::setUIThermalBalancing()
{
    // set thermal balancing flag in UI for solve handler registration
    m_pUI->setThermalBalancing(m_bThermalBalancing);
}
#endif

void Ep2d_diff::setVariantSpecificLoopSettings(SeqLim& rSeqLim, MrProt& rMrProt)
{
    //  Ep2d_diff allows multiple concatenations and Distance factor < 0
    m_mySeqLoop.setdDistFacMinIfConcNo(0);

    // Set diffusion loop length.
    // Note: Diff has to be prepared.
    m_mySeqLoop.setDiffusionLoopLength(m_EPIKernel.getTotalScans());

    if (rMrProt.getsPhysioImaging().getlSignal1() == SEQ::SIGNAL_RESPIRATION)
    {
        // causes a shift of the TRFill time to TRFillEnd within setFillTimes
        // used within TSE for Resp. Triggering in order to reach the shortest possible acq. window
        m_mySeqLoop.setShiftTRFillToTRFillEnd(true);
    }
    else
    {
        m_mySeqLoop.setShiftTRFillToTRFillEnd(false);
    }

    m_mySeqLoop.setPhaseCorScans(false);

    // Depending on the filter mode, additional PrepScans are acquired for the
    // dynamic field correction for each slice (see SBBDiffusion)
    m_lAdjPrepScans = m_EPIKernel.getNoOfAdjScans();
}

bool Ep2d_diff::initializeEPIKernel()
{
    // Enable balance model
    // (if current gradient system supports GPA balance models)
    m_EPIKernel.setUseGPABalance(true, SysProperties::getGradMaxAmpl(SEQ::GRAD_FAST));

    // set prephaser after diffusion block
    if (m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/prephase_after_SBBDiffusion", true))
    {
        m_EPIKernel.setPrephaseROAfterRTEBPlugIn(true);
        m_EPIKernel.setPrephaseBlipsAfterRTEBPlugIn(true);
    }
    else
    {
        SEQ_TRACE_ALWAYS.print("m_EPIKernel.setPrephaseAfterRTEBPlugIn (false)");
    }
    
    // call base implementation
    return Ep2d::initializeEPIKernel();
}

void Ep2d_diff::configureReorderInfo(MrProt& rMrProt)
{
    // set PF info to the reorderinfo structure
    setPartialFourierToReorderInfo(rMrProt, &m_REOInfo);

    // call base class (common) implementation
    Ep2d::configureReorderInfo(rMrProt);
}

void Ep2d_diff::setPartialFourierToReorderInfo(MrProt& rMrProt, ReorderInfo* pReorderInfo) const
{
    switch (rMrProt.kSpace().phasePartialFourierFactor())
    {
        case SEQ::PF_AUTO:
        case SEQ::PF_OFF:
            pReorderInfo->usePrivatePartialFourierFactors(1.00);
            break;
        case SEQ::PF_7_8:
            pReorderInfo->usePrivatePartialFourierFactors(0.87);
            break;
        case SEQ::PF_6_8:
            pReorderInfo->usePrivatePartialFourierFactors(0.75);
            break;
        case SEQ::PF_5_8:
            pReorderInfo->usePrivatePartialFourierFactors(0.66);
            break;
        case SEQ::PF_HALF:
            pReorderInfo->usePrivatePartialFourierFactors(0.50);
            break;
    }
}

void Ep2d_diff::setTriggerAndOscBit(long lKernelMode)
{
    // switch off blips for phase corr scan (only if not PAT)
    // if (pMrProt->MrProtocolData::MrPatData().getucPATMode()==SEQ::PAT_MODE_NONE)
    // m_EPIKernel.setExecuteKernelAsPhaseCorrectionScan (lKernelMode==KERNEL_PHASECOR);

    // if the sending of external trigger signals is activated using an ini file, the extern trigger is sent
    if (m_bSetExternalTrigger)
        m_EPIKernel.setDoNotSendExtTrigger(false);

    // otherwise: not (standard case, as this interferes with MR Elastography measurements in the liver case)
    else
        m_EPIKernel.setDoNotSendExtTrigger(true);
}

void Ep2d_diff::setSPAIRSpoilingType()
{
    // ------- ini option needs to be removed after feedback is positive --------------
    // set SPAIR SBB to new mode
    SeqBuildBlockOptfs* pSPAIR_SBB = nullptr;

    pSPAIR_SBB = m_mySeqLoop.getpOptfs();

    if (pSPAIR_SBB)
    {
        pSPAIR_SBB->seteSpoilingType(SeqBuildBlockOptfs::SPOILING_BEFORE_AND_AFTER_RF);
    }
    // ------- ini option needs to be removed after feedback is positive --------------
}

void Ep2d_diff::prepareDFCBookkeeping()
{
    const size_t uNSize = (m_mySeqLoop.getRepetitionsToMeasure() + 1) * m_EPIKernel.getTotalScans() + m_lAdjPrepScans;
    m_asAcqInRep.clear();
    m_asAcqInRep.resize(uNSize, std::pair<uint64_t, uint64_t>(0, 0));
}

void Ep2d_diff::setReducedIRThicknessForShortTI()
{
    m_mySeqLoop.setScaleIRThickness(m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/IR_slice_thickness", 1.25));
}

bool Ep2d_diff::isSystemCompatibleWithFlavor() const
{
    // all systems are compatible with diffusion
    return true;
}

void SEQ_NAMESPACE::Ep2d_diff::setLastScanInMeasFlagForB0Correction(MrProt& rMrProt, long lSlice, long lShot)
{
    MrProtFacade protFacade(rMrProt);

    if (protFacade.isSliceAcceleration())
    {
        if ((m_mySeqLoop.getInnerSliceCounter() == SMSProperties::getNReducedSlices(rMrProt) - 1)
            && ( lShot                       == rMrProt.getsFastImaging().getlSegments()  - 1 ) )
        {
            m_EPIKernel.getReadOutAddress()->getMDH().setLastScanInMeas(true);
        }
    }
    else
    {
        if ((lSlice == rMrProt.sliceSeries().getlSize() - 1)
            && ( lShot                       == rMrProt.getsFastImaging().getlSegments()  - 1 ) )
        {
            m_EPIKernel.getReadOutAddress()->getMDH().setLastScanInMeas(true);
        }
    }
}

void Ep2d_diff::disableDiffusionForPrepScan()
{
    m_EPIKernel.setDiffusionGradientsEnabled(false);
    // Repetitions loop counter of iPAT and SMS reference scans is always zero
    // (setLoopCounters might have set a different value)
    m_EPIKernel.getReadOutAddress()->getMDH().setCrep(0);
    // No LastScanInMeas flag for iPAT and SMS reference scans
    // (might have been set above)
    m_EPIKernel.getReadOutAddress()->getMDH().deleteFromEvalInfoMask(MDH_LASTSCANINMEAS);
}

NLS_STATUS Ep2d_diff::setLastScanInMeasFlagForDFCBookkeeping(MrProt& rMrProt, long lKernelMode)
{
    if ((lKernelMode & KERNEL_CHECK) == 0)
    {
        if (fRTIsReadoutEnabled())
        {
            MdhProxy&            rMDH    = m_EPIKernel.getReadOutAddress()->getMDH();
            const unsigned short ushCslc = rMDH.getCslc(), ushCrep = rMDH.getCrep();
            if (ushCrep >= m_asAcqInRep.size() || (ushCslc >= 128))
            {
                SEQ_TRACE_ERROR.print("Error at %s(%d).", __FILE__, __LINE__);
                return MRI_SEQ_SEQU_ERROR;
            }
            //  In the case of PAT reference scans (which do not get a separate repetition index) readout should be
            //  still disabled In the case of AdjPrepScans (which get a separte repetition index in
            //  SBBDiffusion_Base::setLoopCounters) the readout is already enabled
            if (ushCslc < 64)
            {
                m_asAcqInRep[ushCrep].first |= (static_cast<uint64_t>(1) << ushCslc);
            }
            else
            {
                m_asAcqInRep[ushCrep].second |= (static_cast<uint64_t>(1) << (ushCslc - 64));
            }

            const uint64_t u64NSlc           = static_cast<uint64_t>(rMrProt.sliceSeries().getlSize()),
                u64MaskFirst     = u64NSlc < 64 ? (static_cast<uint64_t>(1) << u64NSlc) - 1 : 0xffffffffffffffff,
                           u64MaskSecond     = u64NSlc < 64 ? 0 : (static_cast<uint64_t>(1) << (u64NSlc - 64)) - 1,
                           u64AcqInRepFirst  = m_asAcqInRep[ushCrep].first,
                           u64AcqInRepSecond = m_asAcqInRep[ushCrep].second;
            rMDH.setLastScanInMeas(
                ((u64AcqInRepFirst & u64MaskFirst) == u64MaskFirst)
                && ((u64AcqInRepSecond & u64MaskSecond) == u64MaskSecond));
            if (rMDH.isLastScanInMeas())
            {
                //  Temporary
                SEQ_TRACE_ALWAYS.print("Last Scan In Meas: REP(%d)", static_cast<int>(ushCrep));
            }
        }
    }

    return MRI_SEQ_SEQU_NORMAL;
}

MrProtocolData::SeqExpoRFInfo Ep2d_diff::calcEnergyOfExtraCSat(long lNumberOfPulses)
{
    return m_CSatFat.getRFInfoPerRequest() * static_cast<double>(lNumberOfPulses);
}

void Ep2d_diff::configureDiffusionSpecificSeqUTSettings(SeqLim& rSeqLim, MrProt& rMrProt)
{
    // The sequences uses the Free loop to step through b-values and diffusion
    // directions. Each completely acquired volume receives a separate
    // Mdh repetition index. Thus, a) the actual number of repetitions does not
    // match the number given in the protocol and b) concatenations get
    // distributed among multiple repetitions.

    // Provide UT with actual number or repetition indices (including DFC adjustment scans)
    if (m_bSequentialVolumeAcquisition)
    {
        SeqUT.setSizeOfDimRep(m_EPIKernel.getTotalScans() * rMrProt.measurements() + m_lAdjPrepScans - 1);
        long lExpectNotOK = m_EPIKernel.getTotalScans() + m_lAdjPrepScans - 1;
        if (rMrProt.getProtData()->getlRepetitions() > 0)
        {
            lExpectNotOK *= rMrProt.measurements();
        }
        // Disable check of LastScanInConcat flags
        SeqUT.SetExpectedNotOk(
            lNoOfLastScanInConcatErr,
            RTEB_ORIGIN_fSEQRunFinish,
            lExpectNotOK,
            "SeqLoop distributes concatenations among multiple repetitions, thus no LastScanInConcat-flags are "
            "set");
    }

    // WE + Bipolar diffusion scheme will enable diffusion gradient reversal, so do not have to check the GS ratio
    // between excitation and refocusing
    MrProtFacade protFacade(rMrProt);

    if (SysProperties::isLowField())
    {
        long const lExpectNotOK
            = m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, FirstMeas)
              + m_mySeqLoop.getKernelRequestsPerMeasurement(rSeqLim, SecondMeas) * rMrProt.getlRepetitions();
        if ((rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_WaterExcitation)
            && (rMrProt.diffusion().getdsScheme() == SEQ::DIFFSCHEME_BIPOLAR))
        {
            SeqUT.SetExpectedNotOk(
                lRatioExciteRefocErr,
                RTEB_ORIGIN_fSEQRunKernel,
                lExpectNotOK,
                "WE + Bipolar diffusion scheme will enable diffusion gradient reversal, so do not have to check "
                "the GS ratio between excitation and refocusing");
        }
        else if (
            (!protFacade.isGradientReversalDiffusion()) && (rMrProt.getsTXSPEC().getucRFPulseType() != SEQ::RF_LOW_SAR)
            && (rMrProt.sliceSeries()[0].getdThickness() < 3.1)
            && (rMrProt.diffusion().getdsScheme() == SEQ::DIFFSCHEME_MONOPOLAR))
        {
            SeqUT.SetExpectedNotOk(
                lRatioExciteRefocErr,
                RTEB_ORIGIN_fSEQRunKernel,
                lExpectNotOK,
                "With monopolar and thin slices without gradient reversal, the GS ratio is intentionally smaller "
                "at low field");
        }
    }
}

long Ep2d_diff::getDiffusionAdjPrepScans() const
{
    return m_lAdjPrepScans;
}

void Ep2d_diff::considerIRBlockForImplicitCoolingPause(
    SeqLim& rSeqLim, MrProt& rMrProt, SeqExpo& rSeqExpo, long lScanTimeSatsEtc, long lScanTimeBasic)
{
    if ((rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
        && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        //  The IR pulse is comparatively long (at the time of implementation: 20640 us)
        //  and the slice selection gradient is usually small. Hence the RF time is considered
        //  The spoiler is moderate (at the time of implemenation 8 mT/m) and therefore the spoiler time is not
        //  considered.
        m_lCoolPauseImplicit += const_cast<SeqBuildBlockIRsel*>(&m_mySeqLoop.getSBBIRsel())->getIRRFDuration();
    }
}

void Ep2d_diff::disableReadoutForPrepScans(long lKernelMode, long lSlice/*=0*/)
{
    if (m_lAdjPrepScans > 0)
    {
        if (lKernelMode == KERNEL_PREPARE)
        {
            fRTSetReadoutEnable(0);
        }
    }

    // call base class implementation in addition
    Ep2d::disableReadoutForPrepScans(lKernelMode);
}

bool Ep2d_diff::isCoolTimeExecutedWithIR(MrProt& rMrProt) const
{
    return ((rMrProt.getsPrepPulses().getucInversion() != SEQ::SLICE_SELECTIVE)
           || (rMrProt.getsPrepPulses().getucIRScheme() != SEQ::IR_SCHEME_SEQUENTIAL));
}

bool Ep2d_diff::isExtraCSatApplied(MrProt& rMrProt)
{
    // check the diffusion-specific forbidden condition
    if ((rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
        && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        return false;
    }

    // otherwise base implementation decides.
    return Ep2d::isExtraCSatApplied(rMrProt);
}

void Ep2d_diff::setLongTRTrigMode(MrProt& rMrProt)
{
    // call base implementation
    Ep2d::setLongTRTrigMode(rMrProt);

    // also execute diffusion-specific logic
    if ((rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
        && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        m_mySeqLoop.setLongTRTrigMode(false);
    }
}

bool Ep2d_diff::isSegmentedPATRefLinesCondition(MrProt& rMrProt) const
{
    // check the diffusion-specific forbidden condition
    if ((rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
        && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        //  Intended for free-breathing WBDWI: Segmentation cannot work!
        return false;
    }

    // otherwise base implementation decides
    return Ep2d::isSegmentedPATRefLinesCondition(rMrProt);
}

bool Ep2d_diff::isMultiConcatsAllowed(MrProt& rMrProt) const
{
    // call base implementation for common EPI multicontrast checks
    bool bMultiConcatsAllowed = Ep2d::isMultiConcatsAllowed(rMrProt);

    if ((rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
        && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        bMultiConcatsAllowed = true;
    }

    return bMultiConcatsAllowed;
}

bool Ep2d_diff::isLoopStructureCompatibleWithB0Correction(MrProt& rMrProt) const
{
    // call base implementation, if that already forbids B0 correction we are done
    if (!Ep2d::isLoopStructureCompatibleWithB0Correction(rMrProt))
        return false;

    if ((rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
        && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        //  Results in TR, TI error
        return false;
    }

    return true;
}

void Ep2d_diff::setZoomItSeqUTExpectations(MrProt& rMrProt)
{
    const long lExpectNotOK
        = rMrProt.getsSliceArray().getlSize() * rMrProt.measurements() + rMrProt.getsSliceArray().getlSize();

    SeqUT.SetExpectedNotOk(
        lAmplSignRFErr, RTEB_ORIGIN_fSEQRunKernel, lExpectNotOK, "ZOOM_2DRF slice selection has a 2D trajectory");
    SeqUT.SetExpectedNotOk(
        lRFAmplValErr,
        RTEB_ORIGIN_fSEQRunKernel,
        lExpectNotOK,
        "ZOOM_2DRF slice-select gradient amplitude is alternating with trajectory. Image orientation and image "
        "geometry are checked in a functional scanner test.");

        SeqUT.SetExpectedOk(
        lMoreGrInXRFErr,
        RTEB_ORIGIN_fSEQRunKernel,
        0,
        "ZOOM_2DRF has more than one active gradient during the RF_PULSE event ()");
}

bool Ep2d_diff::checkDFCWithLongTRTrigMode(SeqLim& rSeqLim, MrProt& rMrProt) const
{
    MrProtFacade protFacade(rMrProt);

    // the combination "dynamic field correction with multiple concats and respiratory compensation"
    // is enabled, otherwise keep the old behavior
    if (!protFacade.isBookkeepingConditionForDFC()
        && (rMrProt.getsDynDistCorrFilter().getucMode() != SEQ::DYN_DISTCORR_NONE)
        && (m_mySeqLoop.isLongTRTrigMode() == false))
    {
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            TEXT_TR(rSeqLim, "Dynamic distortion correction requires complete volumes")
        }
        return false;
    }

    return true;
}

bool Ep2d_diff::checkAcqWindowForRespTriggering(MrProt& rMrProt) const
{
    if (rMrProt.getsPhysioImaging().getlSignal1() == SEQ::SIGNAL_RESPIRATION)
    {
        const int32_t lAcqScanWindow_ms = rMrProt.physiology().scanWindow(SEQ::SIGNAL_RESPIRATION);
        const long    lMinTR            = static_cast<long>(ceil(m_mySeqLoop.getMaxMinTR() * 1e-03));

        if (lAcqScanWindow_ms < lMinTR)
            return false;
    }
    return true;
}

long Ep2d_diff::calcTotalNumberOfVolumesForFreqFeedback(MrProt& rMrProt)
{
    // Total number of expected volumes, including preparation and adjustment scans
    long lTotalNoOfVolumes = m_EPIKernel.getTotalScans() * rMrProt.measurements() + m_lAdjPrepScans + m_lPATRefScans;

    MrProtFacade protFacade(rMrProt);
    if (protFacade.isSliceAcceleration())
        lTotalNoOfVolumes = m_EPIKernel.getTotalScans() * rMrProt.measurements() + m_lAdjPrepScans;

    return lTotalNoOfVolumes;
}

NLSStatus Ep2d_diff::setCompensationPara(MrProt& rMrProt, SeqLim& rSeqLim)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL; // * Return status *

    // Eddy current compensation
    if (SysProperties::isLowField())
    {
        if ((rMrProt.preparationPulses().getlFatWaterContrast() == MrProtocolData::FatWaterContrast_FatSaturation
             || rMrProt.preparationPulses().getlFatWaterContrast() == MrProtocolData::FatWaterContrast_Spair
             || rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_WaterExcitation)
            && rMrProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_ABDOMEN)
        {
            double dEddyCurrentTau = 1000000. * SysProperties::getEddyCurrentCompGradTimeConstant();

            if (dEddyCurrentTau < 0.01)
            {
                SEQ_TRACE_ERROR.print("Eddy current decay time for EC compensation has not been set! ");
                return (MRI_SBB_SBB_ERROR);
            }
            m_EPIKernel.setbCompensationEnable(true);
            m_EPIKernel.setCompensationPara(true, 0.8, dEddyCurrentTau);
        }
        else
        {
            m_EPIKernel.setbCompensationEnable(false);
        }
    }
#ifdef BUILD_WIPParameterTool
    // Eddy current compensation
    const bool bSeqSpecialCard
        = SysProperties::ReadSeqSettingGeneral("EDDYCURRCOMP_WIP_PARAMS_SPECIALCARD_VISIBLE", false, true);
    if (bSeqSpecialCard)
    {
        // prepare sequence-special card
        if (!m_WIPParamTool.prepare(rMrProt, rSeqLim))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("ERROR: Failed to prepare WIPParamTool sequence!");
            }
            return MRI_SEQ_SEQU_ERROR;
        }

        if (rMrProt.getucEddyCurrentComp())
        {
            m_EPIKernel.setbCompensationEnable(true);
            bool   bCompensationDecay    = rMrProt.getsWipMemBlock().getalFree()[POS_bCompSpoilDecay];
            double dCompensationFraction = rMrProt.getsWipMemBlock().getadFree()[POS_dCompSpoilPercentRead] / 100.;
            double dEddycurrentTau       = rMrProt.getsWipMemBlock().getadFree()[POS_dECTau] * 1000.0;
            m_EPIKernel.setCompensationPara(bCompensationDecay, dCompensationFraction, dEddycurrentTau);
        }
        else
        {
            m_EPIKernel.setbCompensationEnable(false);
        }
    }
#endif // BUILD_WIPParameterTool
    return lStatus;
}