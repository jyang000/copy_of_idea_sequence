//----------------------------------------------------------------------------------
// <copyright file="MREProperties.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2016-2021. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

// Default Values
#define MRE_DRIVER_DEFAULT_FREQUENCY         60.0  // Hz
#define MRE_DRIVER_DEFFAULT_NUMBEROFBURSTS   3
#define MRE_MAX_NUMBEROFMEGs                 2
#define MRE_DRIVER_DEFAULT_AMPLITUDE         50.0  // % of max amplitude

#define MRE_SETTINGS_DEFAULT_MEGLENGTH       16.66667 // ms
#define MRE_SETTINGS_DEFAULT_RAPIDFACTOR     0.5  //

#define MRE_MIN_FRACTIONAL_FACTOR            0.50 //
#define MRE_ALLOWED_FRACTIONAL_FACTOR        0.65 // Currently only one fractional factor is allowed

#define MEG_AMPLITUDE_REDUCTION_FACTOR       0.95

#define MEG_DEFAULT_TRIGGER_STEPS            4
#define MEG_RISETIME_HARD_DEFAULT            10
#define MEG_RISETIME_HARD_LOW               25
#define MEG_AMPLITUDE_HARD_HIGH              40

#define ITRLIMIT                             30  // defines the TR maximum value for rapid MRE

#include <iostream>
#include "MrCommon/UTrace/Macros.h"
//#include "MrGlobalDefinitions/MrBasicTypes.h"

/// contains the parameters that are used for setting the MR
/// Elastography hardware (active driver) and the corresponding 
/// set and get methods
/// \todo set up driver communication
class MREDriverSettings
{
public:
    double getFrequency();
    void setFrequency(double Freq);

    double getAmplitude();
    void setAmplitude(double Ampl);

    int32_t getNumberBursts();
    void setNumberBursts(int32_t NBurst);


    MREDriverSettings();
    MREDriverSettings(double dFrequency, int32_t iNumberOfBursts, double dAmplitude);
    ~MREDriverSettings();

    friend UTrace::ITraceStream& operator<<(UTrace::ITraceStream&, MREDriverSettings&);

private:
    double m_dFrequency;            /// driver frequency [Hz], for liver typically 60.1 Hz
    int32_t m_iNumberBursts;        /// number of bursts which means: how many sinus waves are played after one trigger. Default value is '3'
    double m_dAmplitude;            /// driver amplitude [%]. Default value is '50'
};

/// class for setting up the MRE configuration
/// e.g. defining rapid MRE, fractional MRE, ....
class MREProperties
{
public:
    enum eMEGAxis { MEG_Slice, MEG_Read, MEG_Phase, MEG_None };
    enum eSeqType { greMRE, epiMRE };

    bool getAllowRapid();
    void setAllowRapid(bool RapidAllowed);

    bool getAllowFractional();
    void setAllowFractional(bool FractionalAllowed);

    double getRapidFactor();
    void setRapidFactor(double RapFact);

    double getFractionalFactor();
    void setFractionalFactor(double FractFact);
    double getMinFractionalFactor();
    // do not allow setting of the minimal fraction MRE factor from outside
    //void setFractionalFactor(double FractFact);

    double getAmpReductionFactor();
    void setAmpReductionFactor(double ReductFact);

    double getMEGWavelength();
    void setMEGWavelength(double Wavelength);

    MREDriverSettings getDriverSettings();
    void setMREDriverSettings(MREDriverSettings);

    eMEGAxis getMEGAxis();
    void setMEGAxis(eMEGAxis MEGAxisOrientation);


    MREProperties();
    MREProperties(bool bIsFractionalAllowed, bool bIsRapidAllowed);
    virtual ~MREProperties();

    friend UTrace::ITraceStream& operator<<(UTrace::ITraceStream&, MREProperties&);


protected:
    bool m_bAllowRapid;            /// allow for rapid MRE
    bool m_bAllowFractional;       /// allow for fractional MRE

    int32_t m_MaxNumbersOfDifferentMEGs;

    double m_dFractionalFactor;    /// current fractional MRE factor
    double m_dMinFractionalFactor; /// the minimally allowed MEG reduction 
    double m_dRapidFactor;         /// the rapid factor; default value is 0.5
    double m_dAmpReductionFactor;  /// the MEGs amplitude reduction factor

    double m_dWavePeriod; 
    eMEGAxis m_eMEGAxis;

    MREDriverSettings m_DriverSettings;
};
