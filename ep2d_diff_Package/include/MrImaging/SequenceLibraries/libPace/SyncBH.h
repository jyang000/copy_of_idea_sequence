//  -----------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//  -----------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4\pkg\MrServers\MrImaging\seq\common\libPace\SyncBH.h
//     Version:
//      Author: stemal8q
//        Date: n.a.
//
//        Lang: C++
//
//     Descrip: Class SyncBH can be used to synchronize the sequence with
//              breath-hold commands.
//
//     Classes:
//
//  -----------------------------------------------------------------

#ifndef _SyncBH_h
#define _SyncBH_h


//  Definition of NLS_STATUS
#include "MrGlobalDefinitions/MrResult.h"



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

class SEQData;

namespace PACE
{
    struct SyncBHData;
    class __IMP_EXP SyncBH
    {
    public:
        //  Feeds an event block of duration max(i32MinTimeToResume_us,i32PollIntervall_us).
        //  Checks wether resume method was called as soon as this event block 
        //  is seen on the DSP's. If it was called waitOnResume returns SEQU_NORMAL.
        //  Otherwise it feeds another event block of duration i32PollIntervall_us ...
        //  Return value different from SEQU_NORMAL indicates an internal error.
        NLS_STATUS waitOnResume(int32_t i32MinTimeToResume_us, int32_t i32PollIntervall_us);

        //  Returns the minimum DSP-time of one call to waitOnResume.
        int32_t getMinWaitOnResumeDuration_us(int32_t i32MinTimeToResume_us, int32_t i32PollIntervall_us) const;

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
        NLS_STATUS startCountdown(int32_t i32Duration_us, int32_t i32CBH, int32_t i32NBH, int32_t i32BreathHoldDuration_us, int32_t i32NRelevantRxPerBH);

        //  The member function must be invoked by the resume method of the sequence.
        NLS_STATUS resume(const SEQData& rSEQData);

        //  The member function must be invoked by the resume method of the sequence.
        NLS_STATUS cancel();

        //  Default constructor
        SyncBH();

        //  Destructor
        ~SyncBH();
    protected:
        SyncBHData* m_pData;
    };


} //  namespace PACE

#ifdef __IMP_EXP_DATA
#undef __IMP_EXP_DATA
#endif

#ifdef __IMP_EXP
#undef __IMP_EXP
#endif

#endif // _SyncBH_h

