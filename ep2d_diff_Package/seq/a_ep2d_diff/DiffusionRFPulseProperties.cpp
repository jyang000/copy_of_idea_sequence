/*! 
***************************************************************************
\file   DiffusionPulseProperties.cpp 

\brief  Class for storing all pulse properties

\author Uvo Hoelscher

\b Language: C++

\b Copyright: &copy; Siemens AG (http://www.siemensmedical.com).
All rights reserved.   

***************************************************************************

*/
 
#include "MrImaging/seq/a_ep2d_diff/DiffusionRFPulseProperties.h"

// MRProtFacade
#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"
// MRProtFacade

// MrProt
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProtSliceSeries.h"
// MrProt

// Reference values for gradient and RF calculations
#include "MrMeasSrv/SeqIF/csequence.h"
#include "MrImaging/SequenceLibraries/libPace/MODULE/MODULE_Routines.h"

using namespace SEQ_NAMESPACE;

// Some comments on the RF pulse design of the diffusion sequence:
//
// For (double) spin echoes, the slice profile is determined by a complex interplay
// of the excitation and refocussing pulses. With the current implementation, the
// slice profile is predominantly defined by the refocussing pulses. Note that
// profile changes with the number of refocussing pulses (i.e. a twice refocused
// spin echo ends up with thinner slices than a single refocused spin echo if 
// identical RF pulse shapes are used): thus each diffusion encoding module has
// to ensure a proper selection of RF pulse properties.
//
// It should be possible to realize a slice thickness of 2mm (for 1.5T) and 1mm (for
// 3T), respectively. Considering the gradient performance of current Magnetoms and
// taking into account oblique slice orientations, slice selection gradients must not
// exceed GradMaxAmplFast for the given thickness. We aim for a slice selection gradient
// not higher than 24mT/m for a slice thickness of 1mm.
//
// In order to avoid unambiguity ('3rd arm') artifacts, the slice selection gradients of 
// the refocussing pulses have to deviate from those of the excitation pulses by at least
// 20% (this is also checked by the SeqUT for any spin-echo sequence). Using even larger
// deviations helps to reduce signal from undesired, chemically shifted spin species
// (e.g. fat signal). Note however, that this signal reduction occurs for any off-resonance
// spins (e.g. static or dynamic field inhomogeneities).
//
// With low slice selection gradients, slice bending (due to static or dynamic field 
// inhomogeneities) increases. We aim for a slice selection gradient of at least 3mT/m 
// for a slice thickness of 5mm.
//
// RF pulses should not be clipped because of RFPA limitations. Since diffusion 
// imaging is also used in the body, peak amplitudes should not exceed
// 20uT (1.5T) and 16uT (3T), respectively.
//
// Since the duration of the RF pulses contributes to the minimum possible echo time
// TE, they should be kept as short as reasonable.

// We decided to use a symmetric external Sinc-like RF pulse with one side lobe 
// on either side. It provides a reasonable compromise between slice profile
// and excitation duration. For improved slice profile quality (which is e.g.
// required for gradient reversal fat saturation) we can also use a sinc pulse
// with a higher bandwidth time product and two side lobes on either side.
// For field strengths higher than 4T, the pulse duration is extended in order 
// to reduce SAR.


// All pulse durations are defined in this class. They have previously (VE11A and
// earlier) been defined in multiple places:
// * Ep2dDiffCallbackConfigExcitation::fConfigureExcitationForEPIKernel (excitation)
// * SBBDiffusion_Base::prepParameters (refocussing)
// * Diffusion_Stejskal::prepTiming (refocussing)
// * Diffusion_Bipolar::prepTiming (refocussing)


// Here is an overview of the different pulse durations. Durations will double
// if we chose the 'use better slice profile' option

//                \            <2.0T     |     <4.0T     |   >=4.0T    
//                 \ Exc    #1      #2   |  #1      #2   |  #1      #2
//             Refoc\      1792    3584  | 2048    4608  | 2560    4864
//          ------------------------------------------------------------        Ratio
//   <2.0T  #1  2560 |     0.700   1.400 |               |                      Exc : Refoc
//          #2  4864 |     0.368   0.737 |               |
//          ------------------------------------------------------------
//   <4.0T  #1  3328 |                   | 0.615   1.385 |
//          #2  5632 |                   | 0.333   0.818 |
//          ------------------------------------------------------------
//  >=4.0T  #1  6144 |                   |               | 0.417  0.792
//          #2 10240 |                   |               | 0.250  0.475
//
//
// For the selected durations, the following conditions are always fulfilled:
// (G_refoc / G_exc) is >= (max.empirical slice thickness factor) * 1.2  OR 
//                      <= (min.empirical slice thickness factor) / 1.2, respectively.
//
//      Maximum empirical slice thickness factor (double spin-echo schemes): 1.15
//      Minimum empirical slice thickness factor (single spin-echo schemes): 1.00
//   => Ratio of Excitation and refocusing duration larger than (1.2 * 1.15) = 1.380 
//                                              or smaller than (1.00 / 1.2) = 0.833 required


// optimization for eMeRge: 
// boundary condition: eMeRge Peak B1: 23.5uT; eMeRge gradient: 24/22/22

// 1. Reduce pulse duration of RF Normal to use higher B1(result of sequence simulaiton(uT): Singleband Excit: 17.5; Single band Refoc: 20.8; Multiband Excite(factor 2): 19.3 ; Multiband Refoc (factor 2): 19.6);
// 2. eMeRge low SAR mode = 1.5T normal mode
// 3. In case of slice selection gradient reversal(SSGR), do not extend the pulse duration; Reason:due to the water-fat chemical shift (in Hz) is low, SSGR is not as important as at high field.


// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

DiffusionRFPulseProperties::DiffusionRFPulseProperties()
{
    // set field strength
    m_eFieldStrength = undefined;

    if (SysProperties::isLowField())
        m_eFieldStrength = system_LowField;
    else if (SysProperties::is1_5TSystem())
        m_eFieldStrength = system_1_5T;
    else if (SysProperties::is3TSystem())
        m_eFieldStrength = system_3T;
    else if (SysProperties::is7TSystem())
        m_eFieldStrength = system_7T;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================


long DiffusionRFPulseProperties::getPulseDurationExcitationMultiBand_us(MrProt &rMrProt)
{
    // ==============================================================
    // preparation of all information
    // ==============================================================
    MrProtFacade     protFacade(rMrProt);
    SEQ::RFPulseType eRFPulseType = rMrProt.getsTXSPEC().getucRFPulseType();

    uint32_t iDuration = m_debugSettings.getDefaultSetting<int32_t>("EP2D_DIFF/pulse_duration_ex_multi", 0);
    if(iDuration > 0)
        return iDuration;

    if ((eRFPulseType == SEQ::RF_OPTIMIZED) && rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
    {
        const long lBaseDuration_us         = 100;
        const long lIRRFDuration_us         = 12800; // duration of IR pulse when SMS is used, defined in a_resolve.cpp
        long       lExcitationRFDuration_us = 0;

        if (getMatchedPulseDurationExcitation(rMrProt, lExcitationRFDuration_us, lBaseDuration_us, lIRRFDuration_us))
        {
            return lExcitationRFDuration_us;
        }
    }

    // RF_NORMAL
    if (eRFPulseType == SEQ::RF_NORMAL)
    {
        switch (m_eFieldStrength)
        {
        case system_LowField:
            return 2000;
        case system_1_5T:
            return 4000;
        case system_3T:
            return 2000;
        case system_7T:
            return 6000;
        case undefined:
            return 0;
        }
    }

    // RF_LOW_SAR
    // For RF_OPTIMIZED, use setting from RF_LOW_SAR
    else if (eRFPulseType == SEQ::RF_LOW_SAR || eRFPulseType == SEQ::RF_OPTIMIZED)
    {
        switch (m_eFieldStrength)
        {
        case system_LowField:
            return 4000;
        case system_1_5T:
            return 6000;
        case system_3T:
            return 4000;
        case system_7T:
            return 8000;
        case undefined:
            return 0;
        }
    }

    return 0;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

long DiffusionRFPulseProperties::getPulseDurationExcitationSingleBand_us(MrProt &rMrProt, double dMaxGradientStrength)
{

    uint32_t iDuration = m_debugSettings.getDefaultSetting<int32_t>("EP2D_DIFF/pulse_duration_ex_single", 0);
    if(iDuration > 0)
        return iDuration;

    // ==============================================================
    // preparation of all information
    // ==============================================================

    // Deal with different gradient strength per slice thickness. The following conditions need to 
    // be fulfilled for the pulses:
    // 1) SE2560A90.SE90_12A2_2:    10 mm * 1.2 mT/m * 5120 us = const
    // 2) SE3D5120A90.SE90_24A2_1:  10 mm * 2.4 mT/m * 5120 us = const
    const double    dConstForSliceThickness     = 10.0 * 1.2 * 5120.0;    

    const long      lBaseDuration_us            = 256;
    const double    dSliceThickness             = rMrProt.sliceSeries().aFront().thickness();

    long            lMultiplicationFactor       = 0;

    SEQ::RFPulseType eRFPulseType = rMrProt.getsTXSPEC().getucRFPulseType();

    if ((eRFPulseType == SEQ::RF_OPTIMIZED) && rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE)
    {
        long       lExcitationRFDuration_us = 0;
        const long lIRRFDuration_us = 10240; // default duration of non-SMS IR pulse, defined in SBBIRsel
        if (getMatchedPulseDurationExcitation(rMrProt, lExcitationRFDuration_us, lBaseDuration_us, lIRRFDuration_us))
        {
            return lExcitationRFDuration_us;
        }
    }

    MrProtFacade protFacade(rMrProt);

    // ==============================================================
    // determine the multiplication factor
    // ==============================================================

    // gradient reversal is not active
    if (!protFacade.isGradientReversalDiffusion())
    {
        if (m_eFieldStrength == system_LowField)
            lMultiplicationFactor = 5;
        else if (m_eFieldStrength == system_1_5T)
            lMultiplicationFactor = 7;
        else if (m_eFieldStrength == system_3T)
            lMultiplicationFactor = 8;
        else if (m_eFieldStrength == system_7T)
            lMultiplicationFactor = 15;

        // check for thin slice profile ==> longer pulses
        if (dSliceThickness * dMaxGradientStrength * (double)(lBaseDuration_us * lMultiplicationFactor) < dConstForSliceThickness)
        {
            if (m_eFieldStrength == system_LowField)
                lMultiplicationFactor = 8;
            else if (m_eFieldStrength == system_1_5T)
                lMultiplicationFactor = 14;
            else if (m_eFieldStrength == system_3T)
                lMultiplicationFactor = 18;
            else if (m_eFieldStrength == system_7T)
                lMultiplicationFactor = 19;
        }
    }
    // gradient reversal is active
    else
    {
        if (m_eFieldStrength == system_LowField)
            lMultiplicationFactor = 5;
        else if (m_eFieldStrength == system_1_5T)
            lMultiplicationFactor = 14;
        else if (m_eFieldStrength == system_3T)
            lMultiplicationFactor = 18;
        else if (m_eFieldStrength == system_7T)
            lMultiplicationFactor = 19;
    }

    // check for better slice profile
    if(protFacade.isBetterSliceProfileDiffusion())
        lMultiplicationFactor *= 2;

    return lBaseDuration_us * lMultiplicationFactor;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

long DiffusionRFPulseProperties::getPulseDurationExcitation_us(MrProt &rMrProt, double dMaxGradientStrength)
{
    MrProtFacade protFacade(rMrProt);

    // ==============================================================
    // there is completely different logic for single-band and 
    // multi-band (= slice acceleration) pulses
    // ==============================================================
    if(protFacade.isSliceAcceleration())
        return getPulseDurationExcitationMultiBand_us(rMrProt);
    else
        return getPulseDurationExcitationSingleBand_us(rMrProt, dMaxGradientStrength);
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

string DiffusionRFPulseProperties::getPulseFamilyTypeExcitation(MrProt &rMrProt)
{
    string sFamilyType;
    MrProtFacade protFacade(rMrProt);

    // default
    sFamilyType = "SE2560A90.SE90_12A2_2";

    // better slice profile
    if(protFacade.isBetterSliceProfileDiffusion())
        sFamilyType = "SE3D5120A90.SE90_24A2_1";  

    // multi-band pulses for slice acceleration
    if(protFacade.isSliceAcceleration())
        sFamilyType = "sms_multiband.90_excitation";  

    return sFamilyType;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

long DiffusionRFPulseProperties::getPulseDurationRefocusingSingleBand_us(MrProt &rMrProt)
{
    uint32_t iDuration = m_debugSettings.getDefaultSetting<int32_t>("EP2D_DIFF/pulse_duration_refoc_single", 0);
    if(iDuration > 0)
        return iDuration;

    // ==============================================================
    // preparation of all information
    // ==============================================================
    const long          lBaseDuration_us        = 256;
    long                lMultiplicationFactor   = 0;
    SEQ::RFPulseType    eRFPulseType            = rMrProt.getsTXSPEC().getucRFPulseType();

    MrProtFacade protFacade(rMrProt);

    // ==============================================================
    // determine the multiplication factor
    // ==============================================================
    if (m_eFieldStrength == system_LowField)
    {
        if (!(eRFPulseType == SEQ::RF_LOW_SAR || eRFPulseType == SEQ::RF_OPTIMIZED))
            lMultiplicationFactor = 9;
        else
            lMultiplicationFactor = 10;
    }
    else if (m_eFieldStrength == system_1_5T)
    {
        if (!(eRFPulseType == SEQ::RF_LOW_SAR || eRFPulseType == SEQ::RF_OPTIMIZED))
            lMultiplicationFactor = 10;
        else
            lMultiplicationFactor = 19;
    }
    else if (m_eFieldStrength == system_3T)
    {
        if (!(eRFPulseType == SEQ::RF_LOW_SAR || eRFPulseType == SEQ::RF_OPTIMIZED))
            lMultiplicationFactor = 13;
        else
            lMultiplicationFactor = 22;
    }
    else if (m_eFieldStrength == system_7T)
    {
        if (!(eRFPulseType == SEQ::RF_LOW_SAR || eRFPulseType == SEQ::RF_OPTIMIZED))
            lMultiplicationFactor = 24;
        else
            lMultiplicationFactor = 40;
    }

    // check for better slice profile
    if(protFacade.isBetterSliceProfileDiffusion())
            lMultiplicationFactor *= 2;

    return lBaseDuration_us * lMultiplicationFactor;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

long DiffusionRFPulseProperties::getPulseDurationRefocusingMultiBand_us(MrProt &rMrProt)
{
    // ==============================================================
    // preparation of all information
    // ==============================================================
    SEQ::RFPulseType eRFPulseType = rMrProt.getsTXSPEC().getucRFPulseType();

    uint32_t iDuration = m_debugSettings.getDefaultSetting<int32_t>("EP2D_DIFF/pulse_duration_refoc_multi", 0);
    if(iDuration > 0)
        return iDuration;

    // RF_NORMAL
    if (eRFPulseType == SEQ::RF_NORMAL)
    {
        switch (m_eFieldStrength)
        {
        case system_LowField:
            return 6000;
        case system_1_5T:
            return 8000;
        case system_3T:
            return 6000;
        case system_7T:
            return 12000;
        case undefined:
            return 0;
        }
    }

    // RF_LOW_SAR
    else if (eRFPulseType == SEQ::RF_LOW_SAR || eRFPulseType == SEQ::RF_OPTIMIZED)
    {
        switch (m_eFieldStrength)
        {
        case system_LowField:
            return 8000;
        case system_1_5T:
            return 12000;
        case system_3T:
            return 8000;
        case system_7T:
            return 16000;
        case undefined:
            return 0;
        }
    }

    return 0;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

long DiffusionRFPulseProperties::getPulseDurationRefocusing_us(MrProt &rMrProt)
{
    MrProtFacade protFacade(rMrProt);

    // ==============================================================
    // there is completely different logic for single-band and 
    // multi-band (= slice acceleration) pulses
    // ==============================================================
    if(protFacade.isSliceAcceleration())
        return getPulseDurationRefocusingMultiBand_us(rMrProt);
    else
        return getPulseDurationRefocusingSingleBand_us(rMrProt);
}


bool DiffusionRFPulseProperties::getMatchedPulseDurationExcitation(
    MrProt& rMrProt, long& lExcitationRFDuration_us, long lBaseDuration_us, long lIRRFDuration_us)
{
    // to excite identical fat slice with STIR and excitation pulse,
    // the bandwidth and therefore the amplitude of the corresponding slice selection gradients should be matched

    const double dSliceThickness = rMrProt.sliceSeries().aFront().thickness();
    bool         bIsSTIR         = false;
    if (rMrProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE && rMrProt.ti()[0] < 500000)
    {
        bIsSTIR = true;
    }

    double dScaleIRThickness = 2.;

    if (bIsSTIR)
    {
        dScaleIRThickness = 1.25;
    }

    const double dIRThickness_mm = dSliceThickness * dScaleIRThickness;

    //---------------------------------------------------------------------------
    // get the filename of the external pulse file
    //---------------------------------------------------------------------------

    // IR10240H180.IR180_36B1_1 is used in IRPrep as default setting, currently not reset in any product sequences
    const double dIRGSRefGrad_mT_m = 3.6; // RefGrad of IR10240H180.IR180_36B1_1, taken from PulseTool
    double       dExGSRefGrad_mT_m = 1.2; // values taken from PulseTool, entry RefGrad

    std::string sExcFamilyName = getPulseFamilyTypeExcitation(rMrProt);

    if (strcmp(sExcFamilyName.c_str(), "SE2560A90.SE90_12A2_2") == 0)
    {
        dExGSRefGrad_mT_m = 1.2;
    }
    else if (strcmp(sExcFamilyName.c_str(), "SE3D5120A90.SE90_24A2_1") == 0)
    {
        dExGSRefGrad_mT_m = 2.4;
    }
    else if (strcmp(sExcFamilyName.c_str(), "sms_multiband.90_excitation") == 0)
    {
        dExGSRefGrad_mT_m = 1.638198;
    }
    else
    {
        SEQ_TRACE_ERROR << "Error: Unknown family name of excitation RF";
        return false;
    }

    double dIRGSAmpl_mT_m = dIRGSRefGrad_mT_m * (double(REFSLICETHICKmm) / dIRThickness_mm)
                                * (double(EXTSRFREFLEN) / lIRRFDuration_us);

    const double dExcictationDuration_us
            = double(EXTSRFREFLEN) * (dExGSRefGrad_mT_m / dIRGSAmpl_mT_m) * (double(REFSLICETHICKmm) / dSliceThickness);
        lExcitationRFDuration_us = MODULE::IMULT(dExcictationDuration_us, int(lBaseDuration_us));

    return true;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

std::string DiffusionRFPulseProperties::getPulseFamilyTypeRefocusing(MrProt &rMrProt)
{
    string sFamilyType;
    MrProtFacade protFacade(rMrProt);

    // default
    sFamilyType = "SE2560A180.SE180_12A2_2";

    if (m_eFieldStrength == system_7T)
    {
        sFamilyType = "GAUSS5120.B375";
    }
    else
    {
        // check for better slice profile
        if (protFacade.isBetterSliceProfileDiffusion())
            sFamilyType = "SE3D5120A180.SE180_24A2_1";
    }

    if (protFacade.isSliceAcceleration())
        sFamilyType = "sms_multiband.180_refocussing";

    return sFamilyType;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

sRFPulseProperties DiffusionRFPulseProperties::getPulsePropertiesExcitation(MrProt &rMrProt, double dMaxGradientStrength)
{
    sRFPulseProperties sReturnProperties;

    sReturnProperties.lDuration_us  = getPulseDurationExcitation_us(rMrProt, dMaxGradientStrength);
    sReturnProperties.sFamilyName   = getPulseFamilyTypeExcitation(rMrProt);

    getVERSESettingsExcitation(rMrProt, sReturnProperties);

    return sReturnProperties;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

sRFPulseProperties DiffusionRFPulseProperties::getPulsePropertiesRefocusing(MrProt &rMrProt)
{
    sRFPulseProperties sReturnProperties;

    sReturnProperties.lDuration_us  = getPulseDurationRefocusing_us(rMrProt);
    sReturnProperties.sFamilyName   = getPulseFamilyTypeRefocusing(rMrProt);

    getVERSESettingsRefocusing(rMrProt, sReturnProperties);

    return sReturnProperties;
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

void SEQ_NAMESPACE::DiffusionRFPulseProperties::getVERSESettingsExcitation(MrProt &rMrProt, sRFPulseProperties &sProperties)
{
    MrProtFacade protFacade(rMrProt);

    if(protFacade.isSliceAcceleration())
        getVERSESettingsExcitationMultiBand(rMrProt, sProperties);
    else
        getVERSESettingsExcitationSingleBand(rMrProt, sProperties);
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

void SEQ_NAMESPACE::DiffusionRFPulseProperties::getVERSESettingsExcitationMultiBand(MrProt &rMrProt, sRFPulseProperties &sProperties)
{
    SEQ::RFPulseType eRFPulseType = rMrProt.getsTXSPEC().getucRFPulseType();

    // RF_NORMAL
    if (eRFPulseType == SEQ::RF_NORMAL)
    {
        switch (m_eFieldStrength)
        {
        case system_LowField:
            sProperties.fFactorVERSE = 1.9f;
            sProperties.fRelPlateauLengthVERSE = 0.4f;
            break;
        case system_1_5T:
            sProperties.fFactorVERSE            = 1.9f;
            sProperties.fRelPlateauLengthVERSE  = 0.4f;
            break;
        case system_3T:
            sProperties.fFactorVERSE            = 1.6f;
            sProperties.fRelPlateauLengthVERSE  = 0.5f;
            break;
        case system_7T:
            sProperties.fFactorVERSE            = 2.5f;
            sProperties.fRelPlateauLengthVERSE  = 0.4f;
            break;
        case undefined:
            break;
        }
    }

    // RF_LOW_SAR
    else if (eRFPulseType == SEQ::RF_LOW_SAR || eRFPulseType == SEQ::RF_OPTIMIZED)
    {
        switch (m_eFieldStrength)
        {
        case system_LowField:
			sProperties.fFactorVERSE = 1.9f;
			sProperties.fRelPlateauLengthVERSE = 0.4f;
            break;
        case system_1_5T:
            sProperties.fFactorVERSE            = 1.6f;
            sProperties.fRelPlateauLengthVERSE  = 0.6f;
            break;
        case system_3T:
            sProperties.fFactorVERSE            = 1.6f;
            sProperties.fRelPlateauLengthVERSE  = 0.6f;
            break;
        case system_7T:
            sProperties.fFactorVERSE            = 2.5f;
            sProperties.fRelPlateauLengthVERSE  = 0.4f;
            break;
        case undefined:
            break;
        }
    }

    float fVerseFactor = (float)m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/verse_factor_ex_multi", 0.0);
    float fRelPlateau = (float)m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/verse_plateau_ex_multi", 0.0);

    if( fVerseFactor > 0.0)
    {
        sProperties.fFactorVERSE            = fVerseFactor;
        sProperties.fRelPlateauLengthVERSE  = fRelPlateau; 
    }

    // set VERSE bool
    sProperties.bIsVERSE = (sProperties.fFactorVERSE > 1.0f);
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

void SEQ_NAMESPACE::DiffusionRFPulseProperties::getVERSESettingsExcitationSingleBand(MrProt &rMrProt, sRFPulseProperties &sProperties)
{
    // do not change default

    float fVerseFactor = (float)m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/verse_factor_ex_single", 0.0);
    float fRelPlateau = (float)m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/verse_plateau_ex_single", 0.0);

    if( fVerseFactor > 0.0)
    {
        sProperties.fFactorVERSE            = fVerseFactor;
        sProperties.fRelPlateauLengthVERSE  = fRelPlateau; 
    }

    // set VERSE bool
    sProperties.bIsVERSE = (sProperties.fFactorVERSE > 1.0f);
}

void SEQ_NAMESPACE::DiffusionRFPulseProperties::getVERSESettingsRefocusing(MrProt &rMrProt, sRFPulseProperties &sProperties)
{
    MrProtFacade protFacade(rMrProt);

    if(protFacade.isSliceAcceleration())
        getVERSESettingsRefocusingMultiBand(rMrProt, sProperties);
    else
        getVERSESettingsRefocusingSingleBand(rMrProt, sProperties);

}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

void SEQ_NAMESPACE::DiffusionRFPulseProperties::getVERSESettingsRefocusingMultiBand(MrProt &rMrProt, sRFPulseProperties &sProperties)
{
    SEQ::RFPulseType eRFPulseType = rMrProt.getsTXSPEC().getucRFPulseType();

    // RF_NORMAL
    if (eRFPulseType == SEQ::RF_NORMAL)
    {
        switch (m_eFieldStrength)
        {
        case system_LowField:
            sProperties.fFactorVERSE = 1.9f;
            sProperties.fRelPlateauLengthVERSE = 0.4f;
            break;
        case system_1_5T:
            sProperties.fFactorVERSE            = 1.9f;
            sProperties.fRelPlateauLengthVERSE  = 0.4f;
            break;
        case system_3T:
            sProperties.fFactorVERSE            = 1.6f;
            sProperties.fRelPlateauLengthVERSE  = 0.5f;
            break;
        case system_7T:
            sProperties.fFactorVERSE            = 1.6f;
            sProperties.fRelPlateauLengthVERSE  = 0.5f;
            break;
        case undefined:
            break;
        }
    }

    // RF_LOW_SAR
    else if (eRFPulseType == SEQ::RF_LOW_SAR || eRFPulseType == SEQ::RF_OPTIMIZED)
    {
        switch (m_eFieldStrength)
        {
        case system_LowField:
			sProperties.fFactorVERSE = 1.9f;
			sProperties.fRelPlateauLengthVERSE = 0.4f;
            break;
        case system_1_5T:
            sProperties.fFactorVERSE            = 1.6f;
            sProperties.fRelPlateauLengthVERSE  = 0.6f;
            break;
        case system_3T:
            sProperties.fFactorVERSE            = 1.6f;
            sProperties.fRelPlateauLengthVERSE  = 0.6f;
            break;
        case system_7T:
            sProperties.fFactorVERSE            = 1.6f;
            sProperties.fRelPlateauLengthVERSE  = 0.6f;
            break;
        case undefined:
            break;
        }
    }

    float fVerseFactor = (float)m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/verse_factor_refoc_multi", 0.0);
    float fRelPlateau = (float)m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/verse_plateau_refoc_multi", 0.0);

    if( fVerseFactor > 0.0)
    {
        sProperties.fFactorVERSE            = fVerseFactor;
        sProperties.fRelPlateauLengthVERSE  = fRelPlateau; 
    }

    // set VERSE bool
    sProperties.bIsVERSE = (sProperties.fFactorVERSE > 1.0f);
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================

void SEQ_NAMESPACE::DiffusionRFPulseProperties::getVERSESettingsRefocusingSingleBand(MrProt &rMrProt, sRFPulseProperties &sProperties)
{
    // do not change default

    float fVerseFactor = (float)m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/verse_factor_refoc_single", 0.0);
    float fRelPlateau = (float)m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/verse_plateau_refoc_single", 0.0);

    if( fVerseFactor > 0.0)
    {
        sProperties.fFactorVERSE            = fVerseFactor;
        sProperties.fRelPlateauLengthVERSE  = fRelPlateau; 
    }

    // set VERSE bool
    sProperties.bIsVERSE = (sProperties.fFactorVERSE > 1.0f);
}

// =======================================================================================================
// *******************************************************************************************************
// =======================================================================================================