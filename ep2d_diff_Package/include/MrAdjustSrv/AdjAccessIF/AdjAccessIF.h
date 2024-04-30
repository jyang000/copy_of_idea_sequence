//------------------------------------------------------------------------------
//  Copyright (C) Siemens AG 2009  All Rights Reserved.  Confidential
//------------------------------------------------------------------------------
//
// Project: NUMARIS/4
//    File: 
// Version: 
//  Author: cc_adjust
//    Date: 2014-05-21 16:28:05 +02:00
//
//    Lang: C++
//
// Descrip:
//
// Classes:
//
//------------------------------------------------------------------------------

#ifndef ADJACCESSIF_H_INCLUDED
#define ADJACCESSIF_H_INCLUDED

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------

#include <complex>
#include "MrMeasSrv/MeasCompilerDefs.h"
#include "MrMeasSrv/MeasNuclei/Base/MeasNucleus.h"
#include "MrAdjustSrv/AdjDefines.h"
#include "MrAdjustSrv/AdjData/AdjDynAdjust/AdjDynAdjustVolumeData.h"

//------------------------------------------------------------------------------
// declarations
//------------------------------------------------------------------------------

typedef std::complex<double> Complex;
class MrProt;
class AdjContext;
class AdjResult;
class AdjVector;
class AdjDicoResult;
class AdjDicoResultData;
class AdjDynAdjustContext;
class AdjDynAdjustResult;
class AdjMDSResult;
class AdjMDSResultData;
class AdjTraResultData;

//------------------------------------------------------------------------------
// import/export control
//------------------------------------------------------------------------------

#ifdef BUILD_AdjAccessIF
  #define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

//------------------------------------------------------------------------------
// class AdjAccessIF
//------------------------------------------------------------------------------

//Comfortable interface to the adjustments
class __IMP_EXP AdjAccessIF
{
// friends
friend class AdjUI4ExpertLogic;
friend class AdjUIConfirmDlgLogic;

public:
  //Default constructor
  AdjAccessIF();

  //Destructor
  virtual ~AdjAccessIF();

private:
  //Copy constructor
  AdjAccessIF(const AdjAccessIF& rSource);

  //Semantic assignment
  AdjAccessIF& operator=(const AdjAccessIF& rSource);

private:
  //Get context and result for the given protocol.
  //Returns false, if no such data is available.
  static bool get(AdjContext& rContext, AdjResult& rResult, const MrProt& rProt);

  //Get vector for the given protocol.
  //Returns false, if no such data is available.
  static bool get(AdjVector& rVector, const MrProt& rProt);

  //Help function -  read volume data of AdjDynAdjust Context and Result into AdjDynAdjustVolumeData
  static bool readDynAdjustVolumesFromContextAndResult(AdjDynAdjustVolumeData* asDynAdjustVolumeData[ADJ_DYNADJUST_MAX_VOLUMES], size_t& riDynAdjustVolumeDataSize, const AdjDynAdjustContext& rContext, const AdjDynAdjustResult& rResult);

public:
  //Local shim current result structure (see getLocalShimCurrents)
  struct LocalShimCurrent
  {
    LocalShimCurrent() : dCurrent(0.0), bUseShimChannel(false) {}
    double dCurrent;
    bool bUseShimChannel;
  };

  //Get the variable capacity voltages for the given protocol [mV].
  //Returns false, if no such voltages are available.
  static bool getVariCapVoltages(double adVariCapVoltage[ADJ_TUNELC_MAX_CHANNELS], const MrProt& rProt);

  //Returns the last used frequency for the given nucleus [Hz].
  static double getFrequency(const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Get the frequency for the given protocol and nucleus [Hz].
  //Returns false, if no such frequency is available.
  static bool getFrequency(double& rdFrequency, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);
  static bool getFrequency(long&   rlFrequency, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Get the transmitter reference amplitude for the given protocol and nucleus [V].
  //Returns false, if no such amplitude is available.
  static bool getTraRefAmpl(AdjTraResultData& rTraResult, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Get the transmitter reference amplitude for the given protocol and nucleus [V].
  //Returns false, if no such amplitude is available.
  static bool getTraRefAmpl(double& rdTraRefAmpl, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Get the transmitter reference amplitude and the adjustment excitation mode for the given protocol and nucleus [V].
  //Returns false, if no such amplitude is available.
  static bool getTraRefAmpl(double& rdTraRefAmpl, int& riExcitationMode, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);
  
  //Get the transmitter reference amplitude from the full FoV for the given protocol and nucleus [V].
  //Returns false, if no such amplitude is available.
  static bool getTraRefAmplFullFoV(double& rdTraRefAmplFullFoV, bool& rbIsMeasured, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Get the transmitter reference amplitude and the adjustment excitation mode from the full FoV for the given protocol and nucleus [V].
  //Returns false, if no such amplitude is available.
  static bool getTraRefAmplFullFoV(double& rdTraRefAmplFullFoV, int& riExcitationMode, bool& rbIsMeasured, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Get the Dico test results for the given protocol and nucleus.
  //Returns false, if no such results are available.
  static bool getDicoTest(AdjDicoResult& rDicoResult, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Get the Dico test results for the given protocol and nucleus.
  //Returns false, if no such results are available.
  static bool getDicoTestData(AdjDicoResultData& rDicoResult, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Calculate gain from broadband average power values for given Tx channel and pulse
  //The gain describes the ratio of sqrt((forward Tx power magnitude - reflected Tx power magnitude) * 50) / requested Tx voltage magnitude
  //Returns false if no values for this signature were measured or calculation is not possible
  static bool calcGain(double& rdGain, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H, int iTxChannel = 0, int iPulse = ADJ_DICO_TRAREF_AMPLITUDE, bool bOnlyForward = false);

  //Calculate reflection from narrowband average voltage values for given Tx channel and pulse
  //The reflection describes the ration of reflected Tx voltage magnitude / forward Tx voltage magnitude
  //Returns false if no values for this signature were measured or calculation is not possible
  static bool calcReflection(double& rdReflection, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H, int iTxChannel = 0, int iPulse = ADJ_DICO_TRAREF_AMPLITUDE);

  //Calculate reflection from narrowband average voltage values at RFPA level for given Tx channel and pulse
  //The reflection describes the ration of reflected Tx voltage magnitude / forward Tx voltage magnitude
  //Returns false if no values for this signature were measured or calculation is not possible
  static bool calcRFPAReflection(double& rdReflection, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H, int iTxChannel = 0, int iPulse = ADJ_DICO_TRAREF_AMPLITUDE);

  //Get the UID of the fieldmap for the given protocol [].
  //Returns false, if no such UID is available.
  static bool getFieldMapUID(long& rlFieldMapUID, const MrProt& rProt, bool bB0MapAcq = false);

  //Get the UID of the static field correction fieldmap for the given protocol [].
  //Returns false, if no such UID is available.
  static bool getSFCB0MapUID(long& rlB0MapUID, const MrProt& rProt);
 
  //Get the gradient offsets for the given protocol [DAC units].
  //Returns false, if no such offsets are available.
  static bool getGradientOffsets(double adGradientOffset[3], const MrProt& rProt, int iGPANumber = 0);
  static bool getGradientOffsets(long   alGradientOffset[3], const MrProt& rProt, int iGPANumber = 0);

  //Get the shim currents for the given protocol [mA].
  //Returns false, if no such currents are available.
  static bool getShimCurrents(double adShimCurrent[ADJ_SHIM_MAX_CHANNELS], const MrProt& rProt, int iGPANumber = 0);
  static bool getShimCurrents(float  afShimCurrent[ADJ_SHIM_MAX_CHANNELS], const MrProt& rProt, int iGPANumber = 0);

  //Get the local shim currents for the given protocol:
  // asCurrent[i].dCurrent       : local shim current for channel i [mA]
  // asCurrent[i].bUseShimChannel: use local shim channel i
  //Returns false, if no such currents are available.
  static bool getLocalShimCurrents(LocalShimCurrent asCurrent[ADJ_LOCAL_SHIM_MAX_CHANNELS], const MrProt& rProt);

  //Get the UID of the RF map for the given protocol [].
  //Returns false, if no such UID is available.
  static bool getRFMapUID(long& rlRFMapUID, const MrProt& rProt);

  //Get the Tx scale factors for the given protocol []. 
  //Returns false, if no such factors are available.
  static bool getTxScaleFactors(Complex axTxScaleFactor[ADJ_B1SHIM_MAX_CHANNELS], const MrProt& rProt); 

  //Get the correction factor (for TraRefAmpl) for the given protocol and nucleus [].
  //Returns false, if no such factor is available.
  static bool getCorrFac(double& rdCorrectionFactor, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Get the UID of the channel mixing data for the given protocol [].
  //Returns false, if no such UID is available.
  static bool getChannelMixingMatrixUID(long& rlChannelMixingMatrixUID, const MrProt& rProt);

  //Get the UID of the sensitivity map for the given protocol [].
  //Returns false, if no such UID is available.
  static bool getSensMapUID(long& rlSensMapUID, const MrProt& rProt);

  //Get the dynamic adjustment results for the given protocol.
  //Returns false, if no such results are available.
  static bool getDynamicAdjustments(AdjDynAdjustVolumeData* asDynAdjustVolumeData[ADJ_DYNADJUST_MAX_VOLUMES], size_t& riDynAdjustVolumeDataSize, const MrProt& rProt);

  //Get the MDS prescan results for the given protocol.
  //Returns false, if no such results are available.
  static bool getMDSPrescan(AdjMDSResult& rMDSResult, const MrProt& rProt);

  //Get the MDS prescan results for the given protocol.
  //Returns false, if no such results are available.
  static bool getMDSPrescanData(AdjMDSResultData& rMDSResult, const MrProt& rProt);

public:
  //Returns true if a frequency adjustment measurement is required before each imaging measurement.
  static bool isAdjFrePerformAlways();

  //Returns true if RFMap is required 
  static bool isAdjRFMapRequested(const MrProt& rProt);

  //Returns true if ChannelMixing is required 
  static bool isAdjChannelMixingRequested(const MrProt& rProt);

  //Returns true if CoilSens is required 
  static bool isAdjCoilSensRequested(const MrProt& rProt);

  //Returns true if fieldmap is required 
  static bool isAdjFieldMapRequested(const MrProt& rProt, bool bB0MapAcq = false);

  //Returns true if static field correction B0Map is requested by the protocol
  static bool isAdjSFCB0MapRequested(const MrProt& rProt);

  //Returns true, if local shim is required. The following conditions must be fullfilled
  // 1. local shim is enabled in protocol
  // 2. b0 shim is standard or standard neck
  // 3. local shim elements are selected
  // 4a. protocol is opened in DOT configurator or
  // 4b. position is inside local shim working range
  static bool isLocalShimRequested(const MrProt& rProt, bool bDump = false);

  //Returns true, if a coil position adjustment for the given protocol would be performed
  static bool isCoilPosAdjustmentRequired(const MrProt& rProt);

public:
  //Invalidate all adjustment results inside the adjustment server.
  static bool invalidateAll();

  //Invalidate all Dico test results inside the adjustment server.
  static bool invalidateAllDicoTests();

public:
  //Validate variable capacity voltages inside the adjustment server for the given protocol [mV].
  static bool validateVariCapVoltages(const double adVariCapVoltage[ADJ_TUNELC_MAX_CHANNELS], const MrProt& rProt);

  //Validate system frequency inside the adjustment server for the given protocol and nucleus [Hz].
  static bool validateFrequency(double dFrequency, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);

  //Validate transmitter reference amplitude inside the adjustment server for the given protocol and nucleus [V].
  static bool validateTraRefAmpl(double dTraRefAmpl, const MrProt& rProt, const MeasNucleus& rNucleus = NUCLEUS_1H);
  static bool MEAS_DEPRECATED validateTraRefAmpl(double dTraRefAmpl, const MrProt& rProt, const MeasNucleus& rNucleus, bool bDefExcitationMode);

  //Validate UID of fieldmap inside the adjustment server for the given protocol [].
  //Note: Take care that a fieldmap with the corresponding UID really exists!
  static bool validateFieldMapUID(long lFieldMapUID, const MrProt& rProt);

  //Validate gradient offsets and shim currents inside the adjustment server for the given protocol [DAC units] / [mA] / [mA].
  static bool validateShim(const double adGradientOffset[3], const double adShimCurrent[ADJ_SHIM_MAX_CHANNELS], const LocalShimCurrent asCurrent[ADJ_LOCAL_SHIM_MAX_CHANNELS], const MrProt& rProt);

  // Same as upper, but local shim currents are given in a vector of values with unused channels leaved out
  static bool validateShim(const double adGradientOffset[3], const double adShimCurrent[ADJ_SHIM_MAX_CHANNELS], const double adLocalShimCurrents[ADJ_LOCAL_SHIM_MAX_CHANNELS], int iNumLocalShimCurrents, const MrProt& rProt);

  //Validate UID of RF map inside the adjustment server for the given protocol [].
  //Note: Take care that a RF map with the corresponding UID really exists!
  static bool validateRFMapUID(long lRFMapUID, const MrProt& rProt);

  //Validate Tx scale factors inside the adjustment server for the given protocol [].
  static bool validateTxScaleFactors(const Complex axTxScaleFactor[ADJ_B1SHIM_MAX_CHANNELS], unsigned int iNumScaleFactors, const MrProt& rProt); 

  //Validate correction factor (for TraRefAmpl) inside the adjustment server for the given protocol [].
  static bool validateCorFactor(double dCorrectionFactor, const MrProt& rProt);

  //Validate channel mixing data inside the adjustment server for the given protocol [].
  //Note: Take care that data with the corresponding UID really exists!
  static bool validateChannelMixingMatrixUID(long lChannelMixingMatrixUID, const MrProt& rProt);

  //Validate UID of sensitivity map inside the adjustment server for the given protocol [].
  //Note: Take care that a sensitivity map with the corresponding UID really exists!
  static bool validateSensMapUID(long lSensMapUID, const MrProt& rProt);

  // Save diagnostic data
  // NOTE: call this method only in error case for further investigations 
  //       since the method is time consuming !
  static bool saveDiagnosticData(const MrProt& rProt);

  // Callback for pat not pos / pat not reg dialog user confirmation
  // bPatientRegisterDataIsStillValid: Patient registration data is still related to the current patient
  // bPatientHasMoved: Patient has moved on the table (NOT: table has moved)
  static bool handlePatRegTableHomeEvent(bool bPatientRegisterDataIsStillValid, bool bPatientHasMovedOnTable);
};

//------------------------------------------------------------------------------

#endif
