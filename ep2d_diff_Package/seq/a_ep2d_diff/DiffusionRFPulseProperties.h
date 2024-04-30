/*! 
***************************************************************************
\file   DiffusionPulseProperties.h 

\brief  Class for storing all pulse properties

\author Uvo Hoelscher

\b Language: C++

\b Copyright: &copy; Siemens AG (http://www.siemensmedical.com).
All rights reserved.   

***************************************************************************

*/
 
#pragma once

#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrTXSpec.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrRXSpec.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Application/Application.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrSysSpec.h"
#include "MrProtSrv/Domain/CoreNative/MrPreparationPulses.h"

#include "MrImagingFW/libSeqSysProp/SysProperties.h"

#include "MrImaging/seq/a_ep2d_diff/SequenceDebugSettings.h"

#include <string>
using namespace std;

// dll export
#if defined BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"
// dll export




namespace SEQ_NAMESPACE
{
    struct __IMP_EXP sRFPulseProperties 
    {
        long    lDuration_us;      
        string  sFamilyName;
        bool    bIsVERSE;
        float   fFactorVERSE;
        float   fRelPlateauLengthVERSE;

        sRFPulseProperties()
            : lDuration_us(0)
            , sFamilyName()
            , bIsVERSE(false)
            , fFactorVERSE(1.0f)
            , fRelPlateauLengthVERSE(1.0f)

        {};

        sRFPulseProperties(
            long lDurationInput_us, 
            string sFamilyNameInput,
            bool bIsVERSEInput,
            float fFactorVERSEInput,
            float fRelPlateauLengthVERSEInput)
        {
            lDuration_us            = lDurationInput_us;
            sFamilyName             = std::move(sFamilyNameInput);
            bIsVERSE                = bIsVERSEInput;
            fFactorVERSE            = fFactorVERSEInput;
            fRelPlateauLengthVERSE  = fRelPlateauLengthVERSEInput;
        }
    };

    enum eSystemFieldStrength
    {
        system_LowField,
        system_1_5T,
        system_3T,
        system_7T,
        undefined
    };

    class __IMP_EXP  DiffusionRFPulseProperties
    {
    public:
        // general
        DiffusionRFPulseProperties();

        ~DiffusionRFPulseProperties(){}

        // methods
        sRFPulseProperties getPulsePropertiesExcitation(MrProt &rMrProt, double dMaxGradientStrength);
        sRFPulseProperties getPulsePropertiesRefocusing(MrProt &rMrProt);

    protected:
        long getPulseDurationExcitation_us(MrProt &rMrProt, double dMaxGradientStrength);
        long getPulseDurationExcitationMultiBand_us(MrProt &rMrProt);
        long getPulseDurationExcitationSingleBand_us(MrProt &rMrProt, double dMaxGradientStrength);

        long getPulseDurationRefocusing_us(MrProt &rMrProt);
        long getPulseDurationRefocusingSingleBand_us(MrProt &rMrProt);
        long getPulseDurationRefocusingMultiBand_us(MrProt &rMrProt);

        bool getMatchedPulseDurationExcitation(MrProt& rMrProt, long& lExcitationRFDuration_us, long lBaseDuration_us, long lIRRFDuration_us);

        void getVERSESettingsExcitation(MrProt &rMrProt, sRFPulseProperties &sProperties);
        void getVERSESettingsExcitationMultiBand(MrProt &rMrProt, sRFPulseProperties &sProperties);
        void getVERSESettingsExcitationSingleBand(MrProt &rMrProt, sRFPulseProperties &sProperties);

        void getVERSESettingsRefocusing(MrProt &rMrProt, sRFPulseProperties &sProperties);
        void getVERSESettingsRefocusingMultiBand(MrProt &rMrProt, sRFPulseProperties &sProperties);
        void getVERSESettingsRefocusingSingleBand(MrProt &rMrProt, sRFPulseProperties &sProperties);

        string getPulseFamilyTypeExcitation(MrProt &rMrProt);
        string getPulseFamilyTypeRefocusing(MrProt &rMrProt);

        eSystemFieldStrength m_eFieldStrength;
        SequenceDebugSettings::SequenceDebugSettings m_debugSettings = SequenceDebugSettings::SequenceDebugSettings("USE_EPI_DEBUG_SETTINGS");
    };
}
