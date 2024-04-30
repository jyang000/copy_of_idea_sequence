//-----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential
//-----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//    File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProtocol\UILink\MrUILinkPhysio.h
// Version: \main\24
//  Author: Comp_ProBe
//    Date: 2013-04-11 12:52:40 +02:00
//
//    Lang: C++
//
// Descrip: 
//
// Classes: -
//
//-----------------------------------------------------------------------------
#ifndef _MRUILINKPHYSIO_H
#define _MRUILINKPHYSIO_H

#include <vector>

#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkSelection.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkString.h"
#include "MrProtSrv/Domain/MrProtocol/UILink/macros.h"

#ifdef BUILD_MrUILink
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

namespace MrProtocolData
{
    class MrProtData;
}                      
                                    
int32_t _maximizeAcqWnd(MrUILinkBase* const pThis);

int32_t _minimizeAcqWnd(MrUILinkBase* const pThis, int32_t lOldAcqWnd);

bool _solveAcqWnd(MrUILinkBase* const pThis);

void getScanWindowAndDelay(const MrUILinkBase* const pThis, int32_t *scanWindow_ms, int32_t *delay_ms);

//  ---------------------------------------------------------------------------
//  physio delay [ms]
//
//  Cardiac-Triggering         TriggerDelay				
//  Cardiac-Gating             GateOn
//  Cardiac-Retrogating        0
//  Respiration                0	
//
bool fUILinkPhysioDelayIsAvailable(LINK_DOUBLE_TYPE* const pThis, int32_t SignalIndex);
double fUILinkPhysioDelayGetValue(LINK_DOUBLE_TYPE* const pThis, int32_t SignalIndex);

//  ---------------------------------------------------------------------------
//  Physio Acq-Part delay == TRfill [ms]
//
//  All signal/methods         (< 0: left; == 0: not available; > 0 right)
//
bool fUILinkPhysioPartDelayIsAvailable(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
double fUILinkPhysioPartDelayGetValue(LINK_DOUBLE_TYPE* const pThis, int32_t pos);

//  ---------------------------------------------------------------------------
//  Physio Acq-Part == TR [ms]
//
//  Cardiac-Triggering         TR
//  Cardiac-Gating             TR
//	Cardiac-Retrogating        TR
//  Respiratory-Triggering     TR
//  Respiratory-Gating         TR
//
bool fUILinkPhysioPartIsAvailable(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
double fUILinkPhysioPartGetValue(LINK_DOUBLE_TYPE* const pThis, int32_t pos);

//  -----------------------------------------------------------------
//  Repetition of physio acq parts
//                                      
//  Cardiac-Triggering         Phases * AqcWndFactor
//  Cardiac-Gating             ScanWindow / AcqPart
//  Cardiac-Retrogating        Phases * AqcWndFactor
//  Respiratory-Triggering     Phases * AqcWndFactor
//  Respiratory-Gating         0
//                             (specified by GATE-OFF-Signal 
//                              from PMU, defined with 
//                              threshold and resp. phase!!)
//
bool fUILinkPhysioRepetitionIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);
int32_t fUILinkPhysioRepetitionGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Physio Event count 
//  (defines which n-th PMU event is to be used to start measurement)
//
//  Cardiac-Triggering         Trigger Pulse
//  All other signals          1
//
bool fUILinkPhysioEventCountIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);
int32_t fUILinkPhysioEventCountGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

//  ---------------------------------------------------------------------------
//  physio signal 
//
//  SEQ::SIGNAL_NONE
//  SEQ::SIGNAL_EKG
//  SEQ::SIGNAL_PULSE
//  SEQ::SIGNAL_BEATSENSOR_CARDIAC
//  SEQ::SIGNAL_EXT
//  SEQ::SIGNAL_RESPIRATION
//
bool fUILinkPhysioSignalIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);
int32_t fUILinkPhysioSignalGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

//  ---------------------------------------------------------------------------
//  physio method 
//
//  SEQ::METHOD_NONE
//  SEQ::METHOD_TRIGGERING
//  SEQ::METHOD_GATING
//  SEQ::METHOD_RETROGATING
//  SEQ::METHOD_SOPE			ROPE or COPE
//
bool fUILinkPhysioMethodIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);
int32_t fUILinkPhysioMethodGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

//  /////////////////////////////////////////////////////////////////
//  Physio Card
//  /////////////////////////////////////////////////////////////////
bool fUILinkPhysioCardGetValue(LINK_BOOL_TYPE* const pThis, int32_t pos);

// Dummy heartbeats
unsigned fUILinkPhysioDummyHeartBeatsGetToolTipId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t pos);
unsigned fUILinkPhysioDummyHeartBeatsGetLabelId (LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos);
bool fUILinkPhysioDummyHeartBeatsIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);
int32_t fUILinkPhysioDummyHeartBeatsGetValue(LINK_LONG_TYPE* const pThis, int32_t pos);
int32_t fUILinkPhysioDummyHeartBeatsSetValue(LINK_LONG_TYPE* const pThis, int32_t value, int32_t pos);
bool fUILinkPhysioDummyHeartBeatsGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Signal and Mode
//
unsigned fUILinkSignalModeGetLabelId(LINK_SELECTION_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkSignalModeIsAvailable(LINK_SELECTION_TYPE* const pThis,int32_t SignalIndex);

unsigned fUILinkSignalModeGetValue(LINK_SELECTION_TYPE* const pThis,int32_t SignalIndex);

bool fUILinkSignalModeGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t SignalIndex);

unsigned fUILinkSignalModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned newVal, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Acquisition Window
//
unsigned fUILinkAcquisitionWindowGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkAcquisitionWindowGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkAcquisitionWindowGetToolTipId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkAcquisitionWindowIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkAcquisitionWindowGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkAcquisitionWindowSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t SignalIndex);

bool fUILinkAcquisitionWindowGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);
int32_t fUILinkAcquisitionWindowInternalSetValue(LINK_LONG_TYPE* const pThis, int32_t x, int32_t SignalIndex);


//  -----------------------------------------------------------------
//  Trigger Pulses
//
unsigned fUILinkTriggerPulsesGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkTriggerPulsesGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkTriggerPulsesIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkTriggerPulsesGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkTriggerPulsesSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t SignalIndex);

bool fUILinkTriggerPulsesGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Gate On
//
unsigned fUILinkGateOnGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkGateOnGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkGateOnIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkGateOnGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkGateOnSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t SignalIndex);

bool fUILinkGateOnGetLimits( LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Gate Ratio
//
unsigned fUILinkGateRatioGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkGateRatioGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkGateRatioIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkGateRatioGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkGateRatioSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t SignalIndex);

bool fUILinkGateRatioGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Trigger Delay
//
unsigned fUILinkTriggerDelayGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkTriggerDelayGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkTriggerDelayIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkTriggerDelayGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkTriggerDelaySetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t SignalIndex);

bool fUILinkTriggerDelayGetLimits( LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Respiratory Threshold
//
unsigned fUILinkThresholdGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkThresholdGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkThresholdIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkThresholdGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkThresholdSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t SignalIndex);

bool fUILinkThresholdGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Gate Off
//
unsigned fUILinkGateOffGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkGateOffGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkGateOffIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkGateOffGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkGateOffSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t SignalIndex);

bool fUILinkGateOffGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Physiological Calculated Phases (= retrogated Images);
//
unsigned fUILinkCalculatedPhasesGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkCalculatedPhasesGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkCalculatedPhasesIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkCalculatedPhasesGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkCalculatedPhasesSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t SignalIndex);

bool fUILinkCalculatedPhasesGetLimits( LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Respiration Phase
//
unsigned fUILinkRespirationPhaseGetLabelId(LINK_SELECTION_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkRespirationPhaseIsAvailable(LINK_SELECTION_TYPE* const pThis,int32_t SignalIndex);

unsigned fUILinkRespirationPhaseGetValue(LINK_SELECTION_TYPE* const pThis,int32_t SignalIndex);

bool fUILinkRespirationPhaseGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t SignalIndex);

unsigned fUILinkRespirationPhaseSetValue(LINK_SELECTION_TYPE* const pThis, unsigned newVal, int32_t SignalIndex);

bool fUILinkRespirationPhaseIDIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkRespirationPhaseIDGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Physiological Measured Phases
//
unsigned fUILinkPhysioPhasesGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t pos);

bool fUILinkPhysioPhasesIsAvailable(LINK_LONG_TYPE* const pThis, int32_t pos);

int32_t fUILinkPhysioPhasesGetValue(LINK_LONG_TYPE* const pThis, int32_t pos);

int32_t fUILinkPhysioPhasesSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t pos);

bool fUILinkPhysioPhasesGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t pos);

unsigned fUILinkPhysioPhasesSolve(LINK_LONG_TYPE* const pThis, char* arg_list[], const void* pVoid, const MrProtocolData::MrProtData*, int32_t pos);

bool fUILinkPhysioPhasesTry(LINK_LONG_TYPE* const pThis, void* pClientMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lPos);

//  -----------------------------------------------------------------
//  Arrythmia Detection
//
unsigned fUILinkArrythmiaDetectionGetLabelId(LINK_SELECTION_TYPE* const pThis, char* arg_list[], int32_t pos);

bool fUILinkArrythmiaDetectionIsAvailable(LINK_SELECTION_TYPE* const pThis,int32_t pos);

unsigned fUILinkArrythmiaDetectionGetValue(LINK_SELECTION_TYPE* const pThis, int32_t pos);

bool fUILinkArrythmiaDetectionGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);

unsigned fUILinkArrythmiaDetectionSetValue(MrUILinkSelection<unsigned>* const pThis, unsigned newVal, int32_t pos);

//  -----------------------------------------------------------------
//  Trigger Window
//
unsigned fUILinkTriggerWindowGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t pos);

bool fUILinkTriggerWindowIsAvailable(LINK_LONG_TYPE* const pThis, int32_t pos);

int32_t fUILinkTriggerWindowGetValue(LINK_LONG_TYPE* const pThis, int32_t pos);

int32_t fUILinkTriggerWindowSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t pos);

bool fUILinkTriggerWindowGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t pos);

//  -----------------------------------------------------------------
//  Average Resp./Cardiac Cycle to be displayed as text
//
bool fUILinkAvgCycleStringIsAvailable(MrUILinkString* const pThis, int32_t SignalIndex);

unsigned fUILinkAvgCycleStringGetLabelId(MrUILinkString* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkAvgCycleStringGetUnitId(MrUILinkString* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkAvgCycleStringGetValue(MrUILinkString* const pThis, char* arg_list[], int32_t SignalIndex);



unsigned fUILinkShortAvgCycleStringGetValue(MrUILinkString* const pThis, char* arg_list[], int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Average Cycle in ms
//
unsigned fUILinkAvgCycleMsGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkAvgCycleMsGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkAvgCycleMsIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkAvgCycleMsGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);


int32_t fUILinkShortAvgCycleMsGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  first and second average cycle
//
int32_t fUILinkTypAvgCycleGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

//  -----------------------------------------------------------------
//  Standard deviation of Average Cycle in ms
//
unsigned fUILinkStdDevMsGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

unsigned fUILinkStdDevMsGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t SignalIndex);

bool fUILinkStdDevMsIsAvailable(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);

int32_t fUILinkStdDevMsGetValue(LINK_LONG_TYPE* const pThis, int32_t SignalIndex);


//  -----------------------------------------------------------------
//  Scan acquisition window in ms
//
unsigned fUILinkScanWindowMsGetLabelId(LINK_DOUBLE_TYPE* const pThis, char* arg_list[], int32_t pos);

unsigned fUILinkScanWindowMsGetUnitId(LINK_DOUBLE_TYPE* const pThis, char* arg_list[], int32_t pos);

bool fUILinkScanWindowMsIsAvailable(LINK_DOUBLE_TYPE* const pThis, int32_t pos);

double fUILinkScanWindowMsGetValue(LINK_DOUBLE_TYPE* const pThis, int32_t pos);

//  -----------------------------------------------------------------
//  Capture RR
// 
void fUILinkCaptureRRCTor(MrUILinkBase* const pThis);

void fUILinkCaptureRRDTor(MrUILinkBase* const pThis);


bool fUILinkCaptureRRIsAvailable(LINK_BOOL_TYPE* const ,int32_t );


unsigned fUILinkCaptureRRLabelId(LINK_BOOL_TYPE* const pThis, char **, int32_t);

bool fUILinkCaptureRRGetValue(LINK_BOOL_TYPE* const pThis, int32_t);

unsigned
fUILinkCaptureRRGetUnitId(LINK_BOOL_TYPE* const pThis, char**, int32_t);

int
fUILinkCaptureRRFormat(LINK_BOOL_TYPE* const pThis, bool nID, char* arg_list[], int32_t /*pos*/);

bool fUILinkCaptureRRSetValue(LINK_BOOL_TYPE* const pThis, bool value, int32_t);
__IMP_EXP bool fUILinkCaptureRRSetValue(LINK_BOOL_TYPE* const pThis, bool value, int32_t, int32_t minimumAcqWindow, short percentageOfAcqWindow);

//  -----------------------------------------------------------------
//  Sampling rate of PMU signals
//
bool fUILinkPMUSamplingRateIsAvailable(LINK_LONG_TYPE* const pThis, int32_t lSignalIndex);

int32_t fUILinkPMUSamplingRateGetValue(LINK_LONG_TYPE* const pThis, int32_t lSignalIndex);

//  -----------------------------------------------------------------
//  Cine Mode
//
DECLARE_STD_UILINK_HANDLERS_SELECTION(CineMode)

// ------------------------------------------------------------------
//  NATIVE functions
//

//  ---------------------------------------------------------------------------
//  Physio "physio_native_mode"
//
bool fUILinkPhysioNativeModeIsAvailable (LINK_SELECTION_TYPE* const pThis, int32_t lPos);
unsigned fUILinkPhysioNativeModeGetLabelId (LINK_SELECTION_TYPE* const pThis, char* arg_list[], int32_t lPos);
unsigned fUILinkPhysioNativeModeGetValue (LINK_SELECTION_TYPE* const pThis, int32_t lPos);
unsigned fUILinkPhysioNativeModeSetValue (LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos);
bool fUILinkPhysioNativeModeGetOptions (LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t lPos);

//  ---------------------------------------------------------------------------
//  Physio "physio_native_flow_sensitivity"
//
bool fFlowSensitivityIsAvailable (LINK_SELECTION_TYPE* const pThis, int32_t lPos);
unsigned fFlowSensitivityGetLabelId (LINK_SELECTION_TYPE* const pThis, char* arg_list[], int32_t lPos);
unsigned fFlowSensitivityGetValue (LINK_SELECTION_TYPE* const pThis, int32_t lPos);
unsigned fFlowSensitivitySetValue (LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos);
bool fFlowSensitivityGetOptions (LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t lPos);
unsigned fFlowSensitivityGetToolTipId(LINK_SELECTION_TYPE* const pThis, char* arg_list[], int32_t pos);


//  ---------------------------------------------------------------------------
//  Physio "physio_native_time1"
//
bool fUILinkPhysioNativeTime1IsAvailable(LINK_LONG_TYPE* const pThis, int32_t lPos);
unsigned fUILinkPhysioNativeTime1GetLabelId (LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos);
unsigned fUILinkPhysioNativeTime1GetUnitId (LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos);
int32_t fUILinkPhysioNativeTime1GetValue  (LINK_LONG_TYPE* const pThis, int32_t lPos);
int32_t fUILinkPhysioNativeTime1SetValue (LINK_LONG_TYPE* const pThis, int32_t lNewVal_ms, int32_t lPos);
bool fUILinkPhysioNativeTime1GetLimits (LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t lPos);

//  ---------------------------------------------------------------------------
//  Physio "physio_native_time2"
//
bool fUILinkPhysioNativeTime2IsAvailable(LINK_LONG_TYPE* const pThis, int32_t lPos);
unsigned fUILinkPhysioNativeTime2GetLabelId (LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos);
unsigned fUILinkPhysioNativeTime2GetUnitId (LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos);
int32_t fUILinkPhysioNativeTime2GetValue  (LINK_LONG_TYPE* const pThis, int32_t lPos);
int32_t fUILinkPhysioNativeTime2SetValue (LINK_LONG_TYPE* const pThis, int32_t lNewVal_ms, int32_t lPos);
bool fUILinkPhysioNativeTime2GetLimits (LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t lPos);

//  ---------------------------------------------------------------------------
//  Physio "physio_native_measurements"
//
bool fUILinkPhysioNativeMeasurementsIsAvailable (LINK_LONG_TYPE* const pThis, int32_t lPos);
unsigned fUILinkPhysioNativeMeasurementsGetLabelId (LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos);
int32_t fUILinkPhysioNativeMeasurementsGetValue  (LINK_LONG_TYPE* const pThis, int32_t lPos);
int32_t fUILinkPhysioNativeMeasurementsSetValue (LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos);
bool fUILinkPhysioNativeMeasurementsGetLimits (LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t lPos);
bool fUILinkPhysioNativeMeasurementsTry(LINK_LONG_TYPE* const pThis, void* pClientMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lPos);
unsigned fUILinkPhysioNativeMeasurementsSolve(LINK_LONG_TYPE* const pThis, char* arg_list[], const void* pClientMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lPos);


//  ---------------------------------------------------------------------------
//  Physio "adaptive triggering"
//
bool fUILinkAdaptiveTriggerIsAvailable(LINK_BOOL_TYPE* const ,int32_t );
unsigned fUILinkAdaptiveTriggerGetLabelId(LINK_BOOL_TYPE* const pThis, char **, int32_t);
bool fUILinkAdaptiveTriggerGetValue(LINK_BOOL_TYPE* const pThis, int32_t);
bool fUILinkAdaptiveTriggerGetOptions(LINK_BOOL_TYPE* const pThis, std::vector<unsigned>&,uint32_t&, int32_t);
bool fUILinkAdaptiveTriggerSetValue(LINK_BOOL_TYPE* const pThis, bool newVal, int32_t);

//  ---------------------------------------------------------------------------
//  Physio "trigger_lock_time"
//
unsigned fUILinkTriggerLockTimeGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos);
unsigned fUILinkTriggerLockTimeGetUnitId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos);
bool fUILinkTriggerLockTimeIsAvailable(LINK_LONG_TYPE* const pThis, int32_t lPos);
int32_t fUILinkTriggerLockTimeGetValue(LINK_LONG_TYPE* const pThis, int32_t lPos);
int32_t fUILinkTriggerLockTimeSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t lPos);
bool fUILinkTriggerLockTimeGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t lPos);

__IMP_EXP bool fUILinkBeatSensorTrainingStepIsAvailable(LINK_BOOL_TYPE* const, int32_t);
__IMP_EXP bool fUILinkBeatSensorTrainingStepGetValue(LINK_BOOL_TYPE* const pThis, int32_t);
__IMP_EXP bool fUILinkBeatSensorTrainingStepSetValue(LINK_BOOL_TYPE* const pThis, bool newVal, int32_t);
__IMP_EXP unsigned fUILinkBeatSensorTrainingStepGetLabelId(LINK_BOOL_TYPE* const pThis, char **, int32_t);
__IMP_EXP bool fUILinkBeatSensorTrainingStepGetOptions(LINK_BOOL_TYPE* const pThis, std::vector<unsigned>&, uint32_t&, int32_t);
__IMP_EXP unsigned fUILinkBeatSensorTrainingStepGetToolTip(LINK_BOOL_TYPE* const pThis, char* arg_list[], int32_t iElementIndex);

//  ---------------------------------------------------------------------------
//  MR_TAG_ACQUISITION_WINDOW_SELECT_MODE_RESP
//
bool fUILinkAcquisitionWindowSelectModeRespIsAvailable(LINK_SELECTION_TYPE* const pThis, int32_t lPos);
unsigned fUILinkAcquisitionWindowSelectModeRespGetValue(LINK_SELECTION_TYPE* const pThis, int32_t lPos);
unsigned fUILinkAcquisitionWindowSelectModeRespSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos);
bool fUILinkAcquisitionWindowSelectModeRespGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t lPos);

//  ---------------------------------------------------------------------------
//  MR_TAG_ACQUISITION_POSITION
//
unsigned fUILinkAcquisitionPositionGetLabelId(LINK_SELECTION_TYPE* const pThis, char* arg_list[], int32_t pos);
bool fUILinkAcquisitionPositionIsAvailable(LINK_SELECTION_TYPE* const pThis, int32_t pos);
unsigned fUILinkAcquisitionPositionGetValue(LINK_SELECTION_TYPE* const pThis, int32_t pos);
unsigned fUILinkAcquisitionPositionSetValue(LINK_SELECTION_TYPE* const pThis, unsigned newVal, int32_t pos);
bool fUILinkAcquisitionPositionGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);

#endif // _MRUILINKPHYSIO_H

//  -----------------------------------------------------------------
//  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential
//  -----------------------------------------------------------------