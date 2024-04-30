//  -----------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//  -----------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\common\libPace\SyncRC.h
//     Version: \main\14
//      Author: stemal8q
//        Date: 2014-03-17 15:46:04 +01:00
//
//        Lang: C++
//
//     Descrip: Class SyncRC can be used to synchronize the sequence with
//              the patient's (R)espiratory (C)ycle.
//
//     Classes:
//
//  -----------------------------------------------------------------

#ifndef _SyncRC_h
#define _SyncRC_h

//  Definition of NLS_STATUS
#include "MrGlobalDefinitions/MrResult.h"
//  Definition of SEQ::PhysioSignal
#include "MrProtSrv/Domain/CoreNative/PhysioDefines.h"


//  Energy structure
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpoRFBlockInfo.h"
//  SEQData object
#include "MrMeasSrv/SeqIF/libMES/SEQData.h"

class MrProt;
class SeqLim;
class SEQCheckMode;
class SEQData;
class Slice;

template <class TYPE>
class VectorPat;

namespace MrProtocolData
{
  class SeqExpo;
}

//  -----------------------------------------------------------------
//  Import - Export - Control
//  -----------------------------------------------------------------

#ifdef __IMP_EXP
#undef __IMP_EXP
#endif 

#ifdef LOCAL_PACE_BUILD
#define __IMP_EXP
#pragma message( "Local PACE build" )
#else  //   LOCAL_PACE_BUILD not defined
#ifdef BUILD_libPACE
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control
#endif  //  LOCAL_PACE_BUILD


namespace PACE
{
    struct SyncRCData;
    struct GRESCOUT_PARA;
    struct B0SCOUT_PARA;
    struct NAVSCOUT_PARA;

    enum 
    {
        //  Single respiratory triggering (default)
        DBLTRIGGER_MODE_OFF                = 0x2
    };
    class __IMP_EXP SyncRC
    {
    public:

        struct __IMP_EXP PARA
        {
            //  Either SEQ::RESP_COMP_TRIGGER or SEQ::RESP_COMP_TRIGGER_AND_FOLLOW
            int32_t              m_iTriggerMode; 
            
            //  Repetition time of the navigator during navigator phase
            int32_t              m_i32ScoutTR_us;
            
            //  Minimum recovery time between imaging and navigator phase
            int32_t              m_i32TRecoveryBeforeScout_us;
            
            //  The number of respiratory intervals to be skipped between two
            //  adjacent trigger events.
            int32_t              m_i32TriggerPulse;
            
            //  The minimum temporal distance between two profile images send
            //  to the host.
            int32_t              m_i32MinSendIntrv_us;
            
            //  The duration of the statistic/learning phase in units of the respiratory cycle.
            int32_t              m_i32NStatRespCycles;

            //  Maximum number of statistic scans.
            //  If the actual number of statistic scans exceeds m_i32NStatScansMax
            //  the learning phase stop's, even if the number of detected respiratory
            //  intervals is less than m_i32NStatRespCycles.
            int32_t              m_i32NStatScansMax;

            //  Number of navigator preparing scans performed after each reset to ensure steady state.
            int32_t              m_i32NPrepScans;

            //  After the scout scan has been send to the DSP's the sequence stops feeding 
            //  the DSP's with meaningful data as the feedback is pending.
            //  It feeds them with a series of synchronization events, each followed
            //  by a pause time. Whenever the GSU sends the echo of the synchronization event
            //  it is checked whether the feedback is still pending. If so another WAKEUP/Pause
            //  interval is inserted, otherwise the sequence breaks out of the loop. 
            //  The return value of the static member function determines the time of a single
            //  pause interval in [us]. The duration of the first pause interval is returned by
            //  sleepIntervAfterScout_us. To retrieve the duration off all following intervals
            //  pollIntervAfterScout_us must be used. If one of the pause 
            //  intervals is choosen too short a FIFO-underflow error may occur. A long pause time
            //  worses the time resolution. 
            int32_t              m_i32SleepIntervAfterScout_us;
            int32_t              m_i32PollIntervAfterScout_us;

            
            //  Double Trigger Mode: Either DBLTRIGGER_MODE_OFF, DBLTRIGGER_MODE_SIMULTANEOUSLY or DBLTRIGGER_MODE_CARDIAC_AFTER_RESP
            int32_t              m_i32DblTriggerMode;

            //  If double triggering is used, the cardiac trigger dealy.
            //  This value is only used to enforce a minimum trigger delay time at
            //  prepare time. The actual value is passed at run time.
            int32_t              m_i32CardiacTriggerDelay_us;

            //  Refrence signal from PMU send to ICE (Default SEQ:: SIGNAL_NONE).
            int32_t              m_i32RefSignal;

            //  Maximum simulated translation [mm] in simulation mode
            double               m_dMaxSimTra_mm;

            // In the simulation mode the navigator image is acquired at a position that is shifted
            //  compared with the position stored in the protocol.
            //  See PACE::SimMode for potential values.
            unsigned char        m_ucSimMode;

            //  Maximum simulated translation [mm] in simulation mode
            double               m_dSimCycleDeviation_Percent;

            //  Channel ID used for communication with ICE (default zero)
            unsigned char        m_ucChaID;

            //  Boolean flag, indicating whether images should be send to the InlineDisplay
            bool                 m_bSendIma;

            //  Boolean flag, indicating whether a noise adjust scan is required,
            //  for improved channel combination
            bool                 m_bNoiseAdjScan;

            //  Boolean flag, indicating whether a phase correction scan is required,
            //  for fast channel combination in complex space (prior to ROFT)
            bool                 m_bPhsCorrScan;
            
            //  Boolean value which indicates whether the trigger threshold is set
            //   manually or automatically. In the manual mode the trigger threshold is
            //  determined by the return value of triggerThreshold, in the
            //  automatic mode the optimal trigger threshold is determined
            //  retrospectively at the end of the learning phase of the
            //  gating algorithm.
            bool                 m_bManualTriggerThreshold;

            //  Boolean value, that indicates whether the sequence runs in scout mode
            //  (no real image are acquired) or not.
            //  If prepared in scout mode, the subsequent call to statPhase is expected
            //  but calls to the waitOn... methods are unexpected.
            bool                 m_bScoutMode;

            //  True announces that trace will be called during the imaging portion of the sequence
            //  However, trace can be called at run-time whether the flag has been set or not.
            //  The actual effect of the flag is that it suppresses the reset of the trigger
            //  after each trigge event.
            bool                 m_bContinuousMode;
            
            //  Usually true, if prepared in context NORMAL
            bool                 m_bContextNormal;
            
            //  Default constructor
            PARA();

            //  Default initialization from protocol
            NLS_STATUS init(const MrProt* pProt, SeqLim* pSeqLim);
        };

        struct TIMING_PARA
        {
            //  Parameters needed to configure trigger

            //  Constant part of the anatomical image acquisition
            //  period, per respiratory trigger halt
            int32_t                    m_i32AcqPeriod_us;
        
            //  Variable part of the anatomical image acquisition period.
            //  Usually zero. Different from zero only if multiple RR-intervals
            //  are acquired per respiratory trigger halt.
            //  Example with m_i32AcqPeriodRR=1:

            //  R-1        C-1                       C-2
            //   |          |                         |
            //   v          v                         v      
            
            //              <--- Acq-Period-1 --->    <--- Acq-Period-2 --->
            //   <- ECG-TD ->
            //              <--------- Avg-RR -------->            
            int32_t                    m_i32AcqPeriodRR;


        
            //  Time delay between the trigger event and the time
            //  point at which anatomical image acquisition starts.
            //  (Needed, if the trigger threshold is determined
            //  automatically. May be zero.).
            //  May contain the time of a non-selective preparation pulse, for example.
            int32_t                   m_i32RespTriggerDelay_us;

            //  If double triggering is used. Time interval between ECG/Pulse trigger
            //  and start of acquisition period
            int32_t                   m_i32CardTriggerDelay_us;


            //  Total number of trigger halts during acquisition phase.
            int32_t                   m_i32NHalts;

            //  Number of concatenations (needed only, if timing is changed at run time)
            int32_t                   m_i32NConc;

            //  Default constructor.
            TIMING_PARA()
                : m_i32AcqPeriod_us(1000000)
                , m_i32AcqPeriodRR(0)
                , m_i32RespTriggerDelay_us(0)
                , m_i32CardTriggerDelay_us(0)
                , m_i32NHalts(1)
                , m_i32NConc(1)
            {}

        };

        class IModifyTiming
        {
        public:
            //  Operator() is called, whenever the static data (# of detected respiratory intervals, median respiratory interval) changes during the statistic phase.
            //  Return values are the boolean variable rbUpdateTiming and the struct TIMING_PARA.
            //  If rbUpdateTiming recievs true the timing parameters written to rNewTiming are send to the trigger algorithm.
            virtual NLS_STATUS operator()(bool& rbUpdateTiming, TIMING_PARA& rNewTiming, int32_t i32NRespIntrv, double dMediaRespIntrv_us) = 0;
        };

        class IContinuousSat
        {
        public:
            //  Called before each scout. Returns the duration of the event bloc(s) to ri32Duarion_us.
            //  The client must ensure that this duration + getTraceDuration_us(0,0) is less
            //  than or equal to ScoutTR

            //  Duration
            //  <---->
            //  +----+----+-------+
            //  | FS |    | Scout |
            //  +----+----+-------+
            //  <---- Scout TR --->
            virtual NLS_STATUS run(int32_t& ri32Duarion_us) = 0;

            //  Returns the energy per Request
            virtual MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest() = 0;
        };


        //  Preparation of the calss.
        //  Either of the four functions must be invoked once before configure
        //  The fourth overloaded function first reads the parameter setting from the protocol
        //  and then calls either the first, the second or the third variant.
        virtual NLS_STATUS prep(const PARA* pPara, const GRESCOUT_PARA* pScoutPara);
        virtual NLS_STATUS prep(const PARA* pPara, const  B0SCOUT_PARA* pScoutPara);
#ifdef SUPPORT_NAV_SCOUT
        virtual NLS_STATUS prep(const PARA* pPara, const  NAVSCOUT_PARA* pScoutPara, MrProt* pProt, SeqLim* pSeqLim, MrProtocolData::SeqExpo* pSeqExpo);
#endif
        virtual NLS_STATUS prep(MrProt* pProt, SeqLim* pSeqLim, MrProtocolData::SeqExpo* pSeqExpo);
        virtual NLS_STATUS prep(const PARA* pPara, MrProt* pProt, SeqLim* pSeqLim, MrProtocolData::SeqExpo* pSeqExpo);

        //  Returns duration of statistic phase when excuted in simulation mode.
        //  The actual duration of the static phase depends on the breathing cycle
        //  of the patient and is unknown.
        //  If the class has not been prepared the return value is undefined.
        virtual int32_t getStatPhaseDuration_us() const;

        //  Returns energy exposure during the statistic phase when excuted in simulation mode.
        //  The actual enegy exposure during the static phase depends on the breathing cycle
        //  of the patient and is unknown.
        //  If the class has not been prepared the return value is undefined.
        virtual MrProtocolData::SeqExpoRFInfo getStatPhaseRFInfo(IContinuousSat* pIContinuousSat = 0) const;

        //  Returns duration of one call to waitOnTrigger when excuted in simulation mode.
        //  The actual duration of the call depends on the breathing cycle
        //  of the patient and is unknown.
        //  From the second trigger halt onwards there is an additional pause
        //  between the imaging scan and the navigator.
        //  Therefore the argument 'bFirstHalt'
        //  If the class has not been prepared the return value is undefined.
        virtual int32_t getWaitOnTriggerDuration_us(int32_t i32DelayPriorToTigger_us, SEQ::PhysioSignal iCardiacSignal, int32_t i32CardiacTriggerDelay) const;

        //  Returns energy exposure of one call to waitOnTrigger when excuted in simulation mode.
        //  The actual enegy exposure during the call depends on the breathing cycle
        //  of the patient and is unknown.
        //  If the class has not been prepared the return value is undefined.
        virtual MrProtocolData::SeqExpoRFInfo getWaitOnTriggerRFInfo(int32_t i32DelayPriorToTigger_us, SEQ::PhysioSignal iCardiacSignal, int32_t i32CardiacTriggerDelay) const;
        virtual MrProtocolData::SeqExpoRFInfo getWaitOnTriggerRFInfo(int32_t i32DelayPriorToTigger_us, IContinuousSat* pIContinuousSat);

        //  The method trace plays out a single navigator to track the respiratory cycle and to
        //  sustain steady state of the navigator sequence.
        //  The intention is that this method is called during the imaging portion of the
        //  sequence (e.g. between different slices in an interleaved multi-slice sequence).
        //  Ideally the time between subsequent navigators is the same during the imaging
        //  and the navigator phase of the sequence.
        //  The input parameter pPrevSlice is optional (may be 0).
        //  The object describes the position and orientation of the previous imaging slice.
        //  This information is used to consider saturation stripes in the naviagtor image
        //  caused by the excitation which may disturb the result otherwise.
        virtual NLS_STATUS trace(int32_t i32DelayPriorToNav_us, int32_t i32DelayPastNav, const Slice* pPrevSlice);
        virtual int32_t getTraceDuration_us(int32_t i32DelayPriorToNav_us, int32_t i32DelayPastNav) const;
        virtual MrProtocolData::SeqExpoRFInfo getTraceRFInfo() const;


        //  Configuration of trigger algorithm.
        //  Must be invoked past preparation
        virtual NLS_STATUS configure(const TIMING_PARA* pTimingPara, MrProtocolData::SeqExpo* pSeqExpo);


        //  Optional modification of trigger configuration at run_time (usually called past statPhase and
        //  before trigger phase.
        //  i32Duration_us is the duration of the event block, used to update the trigger algorithm,
        //  which must be at least the value returned by the method getMinModifyDuration_us
        virtual NLS_STATUS modify(const TIMING_PARA* pTimingPara, int32_t i32Duration_us);
        virtual int32_t getMinModifyDuration_us() const;


        //  Performs learning phase of the trigger algorithm.
        //  Must be invoked once at run time prior to the first call to waitOnTrigger.
        //  If pi32NRespIntrv is different from 0 the integer varaible pi32NRespIntrv points
        //  to receives the number of respiratory intervals detected during the statistic phase.
        //  If pi32NRespIntrv is different from 0 the integer varaible pi32NRespIntrv points
        //  to receives the number of respiratory intervals detected during the statistic phase.

        //  If pdMediaRespIntrv_us is different from 0 the double variable pdMediaRespIntrv_us
        //  points to receive the median respiratory interval (breathing cycle) found during the
        //  learning phase. If *pi32NRespIntrv is less than 1 the median respiratory interval
        //  cannot be calculated and *pdMediaRespIntrv_us is undefined.
        virtual NLS_STATUS statPhase(int32_t* pi32NRespIntrv, double* pdMediaRespIntrv_us, IContinuousSat* pIContinuousSat = 0);
        virtual NLS_STATUS statPhase(IModifyTiming* pIModifyTiming);
        virtual NLS_STATUS statPhase(IModifyTiming* pIModifyTiming,PACE::SyncRC::IContinuousSat* pIContinuousSat);


        //  Intended to be used during sequence check.
        //  Executes the navigator i32NScoutCalls times.
        virtual NLS_STATUS check(int32_t i32NScoutCalls);


        //  Executes the navigator sequence until a respiartory trigger is generated.
        //  i32DelayPriorToTigger_us specifies an optional time interval during which
        //  no trigger will be genarted. The Scout will be executed i32DelayPriorToTigger_us/ScouTR
        //  times during this interval. Thereby the expression i32DelayPriorToTigger_us/ScouTR
        //  is truncated to the next integer.
        //  iCardiacSignal != SEQ::SIGNAL_NONE specifies double triggering
        //  If iCardiacSignal == SEQ::SIGNAL_NONE the i32CardiacTriggerDelay argument is ignored.
        virtual NLS_STATUS waitOnTrigger(int32_t i32DelayPriorToTigger_us, SEQ::PhysioSignal iCardiacSignal, int32_t i32CardiacTriggerDelay);

        //  Similar to previous function, except that pIContinuousSat->operator()() is executed once
        //  per scout. Fails, if pIContinuousSat->duration_us()+ScoutDuration is less than scout TR.
        //  Doesn't support double triggering.
        virtual NLS_STATUS waitOnTrigger(int32_t i32DelayPriorToTigger_us, IContinuousSat* pIContinuousSat);

        virtual NLS_STATUS waitOnRespiratoryTrigger(long i32DelayPriorToTigger_us, long lContinuousSatTR, IContinuousSat& pIContinuousSat);
        virtual NLS_STATUS runContinuousSatAndFillTime(long lContinuousSatTR, IContinuousSat& pIContinuousSat);
        virtual bool isTriggerFound(SEQ::PhysioSignal iPhysioSignal);
        virtual MrProtocolData::SeqExpoRFInfo getWaitOnRespiratoryTriggerRFInfo(long i32DelayPriorToTigger_us, long lContinuousSatTR, IContinuousSat& rIContinuousSat);

        //  The member function is usually invoked immediately after waitOnResume returns.
        //  It excutes a single event block of duration i32Duration_us, which triggers
        //  the countdown of the current breath-hold.
        //  Parameters are
        //  i32Duration_us           ... Duration of the event block in microseconds (at least GRAD_RASTER_TIME but usually significantly longer (e.g. 10 ms))
        //  i32CBH                   ... Breath-hold counter of the actual breath-hold (1 based)
        //  i32NBH                   ... Total number of breath-holds for all measurements.
        //  i32BreathHoldDuration_us ... For non triggered measurements: Duration of the current breath-hold in us.
        //                               Must be zero for triggered measurements.
        //  i32NRelevantRxPerBH      ... For triggered measurements: Number of 'relevant ADCs' per breath-hold.
        //                               Must be zero for non-triggered measurements.
        //NLS_STATUS startCountdown(int32_t i32Duration_us, int32_t i32CBH, int32_t i32NBH, int32_t i32BreathHoldDuration_us, int32_t i32NRelevantRxPerBH);

        //  The member function must be invoked by the receive method of the sequence.
        virtual NLS_STATUS receive(const SEQData& rSEQData);

        //  The member function must be invoked by the cancel method of the sequence.
        virtual NLS_STATUS cancel();

        //  If slice tracking is used. Returns the latest shift vector.
        //  If the information is not available the null vector is returned
        virtual void shift_mm(VectorPat<double>& rShift);

        //  Default constructor
        SyncRC();

        //  Destructor
        virtual ~SyncRC();

        //  Returns the prep input parameter m_i32ScoutTR_us.
        //  The return value is undefined unless the class has been prepared.
        //  The method may be useful for clients which work with the default
        //  preparation (prep(MrProt*,SeqLim*,MrProtocolData::SeqExpo*)).
        int32_t getScoutTR_us() const;
    protected:
        SyncRCData* m_pData;

        SEQData m_sModifyData;
        //if log in special is activated writes timestamps to textfile %MRIEXCHANGE%
        virtual NLS_STATUS  dumpTriggerTimeStamp();
        virtual NLS_STATUS waitOnRespTrigger();
        virtual NLS_STATUS waitOnDblTrigger_sd(SEQ::PhysioSignal iCardiacSignal, int32_t i32CardiacTriggerDelay);
        virtual NLS_STATUS waitOnDblTrigger_ld(SEQ::PhysioSignal iCardiacSignal, int32_t i32CardiacTriggerDelay);
    private:
        /// Actually measured time passed since last trigger
        unsigned long m_lActualTimeStamp = 0;

        /// Previously measured time passed since last 
        /// trigger                                    
        unsigned long m_lLastTimeStamp = 0;

    };

    //  To allow forward declarations
    struct SYNCRC_PARA        : public SyncRC::PARA{};
    struct SYNCRC_TIMING_PARA : public SyncRC::TIMING_PARA{};

} //  namespace PACE

#ifdef __IMP_EXP_DATA
#undef __IMP_EXP_DATA
#endif

#ifdef __IMP_EXP
#undef __IMP_EXP
#endif

#endif // _SyncRC_h

