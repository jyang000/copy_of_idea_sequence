//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4\pkg\MrServers\MrImaging\libSBB\SeqLoopPARA.h
//	 Version: \main\23
//	  Author: stemal8q 
//	    Date: 2011-10-24 05:15:57 +02:00
//
//	    Lang: C++
//
//	 Descrip: Definition of structures that encapsulate members used by SeqLoop
//
//	 Classes:
//
//	-----------------------------------------------------------------------------

#pragma once

#ifndef SeqLoopPARA_h
#define SeqLoopPARA_h


#include "MrGlobalDefinitions/MrBasicTypes.h"
#include "MrMeasSrv/SeqIF/libRT/SEQSemaphore.h"

//  Definition of RT_MDSUPDATE_ADJ_NONE
#include "MrMeasSrv/SeqIF/libRT/libRT.h"

#include <limits>
#include <vector>

//  Definition of class SyncBH
#include "MrImaging/SequenceLibraries/libPace/SyncBH.h"

//  Definition of class SyncRC
#include "MrImaging/SequenceLibraries/libPace/SyncRC.h"


struct PACE_PARA
{
    //  Additional concatenation information
    //  Will be moved to SeqConcat
    struct CONC
    {
        //  Measurement time in first/remaining measurements
        long  m_lTA1_us;
        long  m_lTA2_us;
        // sum of all triggerlock times in the concat, needed to more accurately estimate SAR
        long  m_lTotalTriggerLock_us;
        // sum of all triggerdelays plus physio halts (20ms) in the concat, needed to pass UT, and to more accurately estimate SAR
        // note that the physio halt of 20 ms is included in the trigger delay, unless the trigger delay is zero;
        // i.e. for triggered sequences m_lTotalTriggerDelayAndPhysioHalt_us = 20 us for trigger delay = 0,
        // and  m_lTotalTriggerDelayAndPhysioHalt_us = trigger delay otherwise
        long  m_lTotalTriggerDelayAndPhysioHalt_us;
        //  (N)umber of (R)elevant Readouts in first/remaining measurements
        long  m_lNRRx;

        //  Boolean value indicating whether a breath-hold command is given
        //  before this particular breath-hold.
        bool  m_bHalt;

        //  Default Constructor
        CONC()
            : m_lTA1_us(0)
            , m_lTA2_us(0)
            , m_lTotalTriggerLock_us(0)
            , m_lTotalTriggerDelayAndPhysioHalt_us(0)
            , m_lNRRx(0)
            , m_bHalt(true)
        {}

        //  Restores state after construction
        void clear()
        {
            this->m_lTA1_us = 0;
            this->m_lTA2_us = 0;
            this->m_lTotalTriggerLock_us = 0;
            this->m_lTotalTriggerDelayAndPhysioHalt_us = 0;
            this->m_lNRRx   = 0;
            this->m_bHalt   = true;
        }

    };

    struct TRTIFILLTIMES_ARG
    {
        //  Stores arguments passed to SeqLoop member function TrTiFillTimes
        long m_lScanTime;
        long m_lMultiplier;
        long m_lTIMinAdd1;
        long m_lTIMinAdd2;
        long m_lDummyScanTime;
        long m_lDummySBBTime;

        TRTIFILLTIMES_ARG()
            : m_lScanTime(0)
            , m_lMultiplier(0)
            , m_lTIMinAdd1(0)
            , m_lTIMinAdd2(0)
            , m_lDummyScanTime(0)
            , m_lDummySBBTime(0)
        {}

        //  Return value true indicates that TrTiFillTimes was never called
        bool virgin() const
        {
            return this->m_lScanTime      == 0
                && this->m_lMultiplier    == 0
                && this->m_lTIMinAdd1     == 0
                && this->m_lTIMinAdd2     == 0
                && this->m_lDummyScanTime == 0
                && this->m_lDummySBBTime  == 0
                ;
        }
    };

    TRTIFILLTIMES_ARG m_sTRTIFillTimesArg;


    //  Size of the array is equal to the number of concatenations
    std::vector<CONC> m_asConc;
 
    //  The object is used to synchronize the sequence with breath-hold commands
    PACE::SyncBH  m_sSyncBH;

    //  Passed to last two arguments of PACE::fAdaptRSatPos
    //  See MrImaging/seq/common/ibPace/SyncBH.h for explanation.
    int32_t       m_i32PSatNmbr1;
    int32_t       m_i32PSatNmbr2;


    //  Parameters for respiratory Triggering:

    //  Number of respiratory cycles detected during the learning phase of the algorithm
    int32_t m_i32NRespIntrv;
    //  Median respiratory cycle found during the learning phase
    double  m_dRespCycle_us;

    //  The PACE trigger halt
    PACE::SyncRC  m_sSyncRC;


    //  Default constructor
    PACE_PARA()
        : m_sSyncBH()
        , m_i32PSatNmbr1(std::numeric_limits<int32_t>::max())
        , m_i32PSatNmbr2(std::numeric_limits<int32_t>::max())
        , m_i32NRespIntrv(0)
        , m_dRespCycle_us(0)
        , m_sSyncRC()
    {}
};

#endif  //  SeqLoopPARA_h
