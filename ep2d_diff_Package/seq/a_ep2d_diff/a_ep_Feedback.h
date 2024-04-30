//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2010  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\a_ep_Feedback.h
//	 Version:
//	  Author: PLM AW Neuro
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: Generic class for volume-by-volume feedback synchronization.
//            Actual functionality should be implemented in derived classes.
//
//	 Classes: template <class ep_Feedback_Data> class ep_Feedback
//
//	-----------------------------------------------------------------------------
//  Example #1: good case, Ice keeps track with sequence, 3 slices per volume,
//              synchronization at last slice of volume
//
//  - Acquisition Vol#0, Slc#0
//  - Acquisition Vol#0, Slc#1
//  - Acquisition Vol#0, Slc#2
//  - Ice calc & send FB (Vol#0)
//  - SeqReceive StoreCurrFBData (Vol#0)
//      m_lCurrFBVolumeNo     = 0
//      m_bNewFeedbackOccurred = true
//  - SeqRun SyncFeedback (Vol#0)
//      m_lLastVolumeCalledWithSync        = 0
//      m_bNewFeedbackOccurred             = false
//      m_lLastVolumeWithFeedback          = 0
//      m_vbFeedbackPerformedOnVolumeNo[1] = true
//  => Use feedback from Vol#0 for next volume (#1)
//
//  - Acquisition Vol#1, Slc#0
//  - Acquisition Vol#1, Slc#1
//  - Acquisition Vol#1, Slc#2
//  - Ice calc & send FB (Vol#1)
//  - SeqReceive StoreCurrFBData (Vol#1)
//      m_lCurrFBVolumeNo     = 1
//      m_bNewFeedbackOccurred = true
//  - SeqRun SyncFeedback (Vol#1)
//      m_lLastVolumeCalledWithSync        = 1
//      m_bNewFeedbackOccurred             = false
//      m_lLastVolumeWithFeedback          = 1
//      m_vbFeedbackPerformedOnVolumeNo[2] = true
//  => Use feedback from Vol#1 for next volume (#2)
//
//  - Acquisition Vol#2, Slc#0
//  - Acquisition Vol#2, Slc#1
//  - Acquisition Vol#2, Slc#2
//  - Ice calc & send FB (Vol#2)
//  - SeqReceive StoreCurrFBData (Vol#2)
//      m_lCurrFBVolumeNo     = 2
//      m_bNewFeedbackOccurred = true
//  - SeqRun SyncFeedback (Vol#2)
//      m_lLastVolumeCalledWithSync        = 2
//      m_bNewFeedbackOccurred             = false
//      m_lLastVolumeWithFeedback          = 2
//      m_vbFeedbackPerformedOnVolumeNo[3] = true
//  => Use feedback from Vol#2 for next volume (#3)
//
//	-----------------------------------------------------------------------------
//  Example #2: bad case, Ice lags behind sequence at some point, 3 slices per volume,
//              synchronization at last slice of volume
//
//  - Acquisition Vol#0, Slc#0
//  - Acquisition Vol#0, Slc#1
//  - Acquisition Vol#0, Slc#2
//  - Ice calc & send FB (Vol#0)
//  - SeqReceive StoreCurrFBData (Vol#0)
//      m_lCurrFBVolumeNo     = 0
//      m_bNewFeedbackOccurred = true
//  - SeqRun SyncFeedback (Vol#0)
//      m_lLastVolumeCalledWithSync        = 0
//      m_bNewFeedbackOccurred             = false
//      m_lLastVolumeWithFeedback          = 0
//      m_vbFeedbackPerformedOnVolumeNo[1] = true
//  => Use feedback from Vol#0 for next volume (#1)
//
//  - Acquisition Vol#1, Slc#0
//  - Acquisition Vol#1, Slc#1
//  - Acquisition Vol#1, Slc#2
//  - >>> Ice calc & send FB (Vol#1) is delayed, will be called below <<<
//  - SeqRun SyncFeedback (Vol#1)
//      m_lLastVolumeCalledWithSync        = 1
//  => NO update of feedback for next volume (#2)
//
//  - Acquisition Vol#2, Slc#0
//  - Acquisition Vol#2, Slc#1
//  - Acquisition Vol#2, Slc#2
//  - Ice calc & send FB (Vol#1 !!!)
//  - SeqReceive StoreCurrFBData (Vol#1 !!!)
//      m_lCurrFBVolumeNo     = 1
//      m_bNewFeedbackOccurred = true
//  - SeqRun SyncFeedback (Vol#2)
//      m_lLastVolumeCalledWithSync        = 2
//      m_bNewFeedbackOccurred             = false
//      m_lLastVolumeWithFeedback          = 2
//      m_vbFeedbackPerformedOnVolumeNo[3] = true
//  => Use feedback from Vol#1 for next volume (#3) !!!
//
//  - Acquisition Vol#3, Slc#0
//  - Acquisition Vol#3, Slc#1
//  - Acquisition Vol#3, Slc#2
//  - Ice calc & send FB (Vol#2 !!!)
//  - SeqReceive StoreCurrFBData (Vol#2 !!!)
//      Since m_vbFeedbackPerformedOnVolumeNo[2] = false, m_bNewFeedbackOccurred is not set to true
//  - SeqRun SyncFeedback (Vol#3)
//      m_lLastVolumeCalledWithSync        = 3
//  => Feedback from Vol#2 is completely omitted!!!
//  => Next feedback that will be considered originates from Vol#3 (m_vbFeedbackPerformedOnVolumeNo[3] = true)!!!
//	-----------------------------------------------------------------------------

#pragma once

#include "MrMeasSrv/SeqIF/libRT/sSYNC.h"
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"
#include "MrMeasSrv/SeqIF/libRT/SEQSemaphore.h"           // for synchronisation between run and receive functions on MCIR
#include "MrMeasSrv/SeqIF/libRT/libRT.h"                  // fRTEB...
#include "MrImagingFW/libSeqUTIF/libsequt.h"                        // mSEQTest
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProtSliceSeries.h"           // SliceSeries
#include "MrImaging/libSeqUtil/SMSProperties.h"
#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"

#include <vector>

// ---------------------------------------------------------------------------
// Definitions
// ---------------------------------------------------------------------------
static const unsigned long SeqTraceMarker_Feedback = 0x00001000; ///< feedback specific traces
const long                 minEventDuration_us         = 10;
const long                 startTimeOfFeedbackEvent_us = 300;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

// ===========================================================================
/*!
\class ep_Feedback

\brief This template class implements all interfaces required to
       synchronize a volume-by-volume feedback between Ice and sequence.
       The format of the actual feedback data is user definable using
       the template mechanism.

       As a prerequisite, an Ice functor is required that provides
       the desired feedback values for each acquired volume.
       The standard feedback mechanism (fSEQReceive) is used to provide
       the sequence with the current feedback values and the index
       of the corresponding volume. Index counting starts with one
       for the first volume, each acquired volume (including
       e.g. reference or adjustment scans) is considered.

       Before using the central synchronization methods of the
       ep_Feedback class, the instance has to be provided with the
       expected number of volumes using the ::Prep method. Afterwards,
       ::SyncFeedback has to be called once per acquisition with
       the current slice and volume index. Once per volume (for the
       provided slice index), a real time synchronization takes place. This
       mechanism allows to incorporate feedback exactly once
       per volume and in the correct temporal order. Note that the
       volume index passed to ::SyncFeedback has to match the counting
       used within the corresponding Ice functor.

       The sequence passes the data acquired within the real time
       feedback (fSEQReceive) to this class using the method
       ::StoreCurrFBData. A semaphore is provided by ::GetFeedbackSemaphore
       to synchronize the access to relevant member variables.

       Actual feedback functionality should be implemented using a
       derived class. A plugin mechanism is available to support
       dedicated additional preparations and feedback functionality
       (::PrepPlugin and ::SyncFeedbackPlugin).

       Actual feedback values are stored in corresponding member
       variables within ::StoreCurrFBData and can be accessed without
       semaphore protection using ::GetCurrFBData.

*/
// ===========================================================================
template <class ep_Feedback_Data> class ep_Feedback
{
public:
    /// The constructor initializes the member variables.
    ep_Feedback();

    /// This destructor does nothing.
    virtual ~ep_Feedback() = default;

    // no copy and move constructor and assignment operator
    ep_Feedback(const ep_Feedback& right) = delete;
    ep_Feedback& operator=(const ep_Feedback& right) = delete;
    ep_Feedback(ep_Feedback&& right) = delete;
    ep_Feedback& operator=(ep_Feedback&& right) = delete;

    /// Prepare instance
    /** \b Input:
        \n rMrProt, rSeqLim, rSeqExpo, lTotalVolumeNo

        \b Output:
        \n n.a.

        \b Return value:
        \n true = success, false = error

        Prepare everything for feedback, initialize bookkeeping. The total
        number of expected volumes includes every volume for which a valid
        feedback is sent by Ice. It used to set up an array of flags that
        indicate whether feedback has been applied to a certain volume.
        By setting bSyncWithFirstSlice to false, synchronization takes place
        when the last slice of a volume is acquired (feedback will be used
        for the next volume). By default, synchronization takes place when
        the first slice of a volume is acquired (feedback will be used for
        the current volume).

        At the end ::PrepPlugin is called which can be used to implement
        dedicated preparations within derived classes.
    */
    virtual bool Prep(
        MrProt   &rMrProt,                    /**< Imp: The protocol (required for SeqUT)                      */
        SeqLim   &rSeqLim,                    /**< Imp: The sequence limits (required for SeqUT)               */
        SeqExpo  &rSeqExpo,                   /**< Imp: The sequence exports (required for SeqUT)              */
        long      lTotalVolumeNo,             /**< Imp: Total number of expected volumes                       */
        bool      bSyncWithFirstSlice = true  /**< Imp: Synchronization with first or last slice in volume     */
        );

    /// Feedback synchronization and bookkeeping
    /** \b Input:
        \n rMrProt, rSeqLim, rSeqExpo, lVolumeToBeAcquired, lSliceToBeAcquired, bRunSyncEvent

        \b Output:
        \n n.a.

        \b Return value:
        \n true = success, false = error

        Synchronize feedback with real time by playing out an appropriate
        event block. The duration of this block (::getWakeUp) has to be
        considered by the timing calculation of the calling sequence.

        Typically, this method is called within the sequence ::run method
        once before or after each kernel execution. Afterwards, the currently
        valid feedback data can be accessed using ::GetCurrFBData and used to
        update the succeeding real time events.

        If a new and usable feedback is available ::SyncFeedbackPlugin is
        called which can be used to implement dedicated functionality
        within derived classes.

        If lCurrVolume is negative, only synchronization RT events are played out
        (no update of feedback bookkeeping). This is useful if SyncFeedback is
        called for kernel executions that actually don't generate feedback (e.g.
        preparation scans or iPAT reference scans).

        By setting bRunSyncEvent to false, multiple feedback instances can
        be concatenated without the need to incorporate multiple sync
        wait times. In that case, the first call of ::SyncFeedback should
        employ the real time synchronization.
    */
    virtual bool SyncFeedback(
        MrProt   &rMrProt,                  /**< Imp: The protocol (required for SeqUT)                      */
        SeqLim   &rSeqLim,                  /**< Imp: The sequence limits (required for SeqUT)               */
        SeqExpo  &rSeqExpo,                 /**< Imp: The sequence exports (required for SeqUT)              */
        long     lCurrVolume,               /**< Imp: Current volume - consistent counting with Ice required */
        long     lCurrSlice,                /**< Imp: Current slice - apply feedback only once per volume    */
        bool     bRunSyncEvent = true       /**< Imp: Flag to disable synchronization RT events              */
        );

    /// Set total duration of feedback real time events (wake up)
    /** \b Input:
        \n lWakeupTime

        \b Output:
        \n n.a.

        \b Return value:
        \n m_lWakeupTime

        For real time synchronization, a real time event block
        (using the wakeup mechanism) is applied. The duration
        of this block can be set here (default is set within
        constructor).
    */
    virtual long setWakeUpDuration (
        long lWakeupTime                    /**< Imp: Duration of synchronization event block [us]            */
        );

    /// Get duration of feedback real time events
    /** \b Input:
        \n bRunSyncEvent

        \b Output:
        \n n.a.

        \b Return value:
        \n Event block duration [us]

        For real time synchronization, a real time event block
        (using the wakeup mechanism) is applied. The duration
        of this block has to be considered by sequence timing
        calculation.

        If the RT synchronization is disabled in ::SyncFeedback,
        this has to be considered here by setting bRunSyncEvent
        to false.
    */
    virtual long getWakeUpDuration (
        bool bRunSyncEvent = true           /**< Imp: Flag to indicate disabled synchronization RT events    */
        );

    /// Provide pointer to semaphore
    /** \b Input:
        \n n.a.

        \b Output:
        \n n.a.

        \b Return value:
        \n Semaphore

        Some member variables might get accessed in parallel
        by the calling sequence (within fSEQReceive and
        fSEQRun, respectively). This semaphore can be used
        to serialize the access.
        */
    virtual SEQSemaphore *GetFeedbackSemaphore();

    /// Store current feedback data if they are acceptable
    /** \b Input:
        \n lCurrFBVolumeNo, sCurrFBData

        \b Output:
        \n n.a.

        \b Return value:
        \n n.a.

        Store current feedback information. This is usually
        called within fSEQReceive to pass the values received
        from the Ice world.
    */
    virtual void StoreCurrFBData(
        long             lCurrFBVolumeNo,         /**< Imp: Current feedback volume index (from Ice) */
        ep_Feedback_Data sCurrFBData              /**< Imp: Current feedback              (from Ice) */
        );

    /// Access current feedback data
    /** \b Input:
        \n n.a.

        \b Output:
        \n sCurrFBData

        \b Return value:
        \n n.a.

        Get most recent valid feedback data
        */
    virtual void GetCurrFBData(
        ep_Feedback_Data &sCurrFBData             /**< Exp: Most recent valid feedback */
        );

protected:

    /// Virtual prepare plugin
    /** \b Input:
        \n lTotalVolumeNo

        \b Output:
        \n n.a.

        \b Return value:
        \n true = success, false = error

        Called at the end of ::Prep - can be used to implement all dedicated
        feedback preparations. Should be overloaded if required.
    */
    virtual bool PrepPlugIn(
        MrProt   & /* rMrProt  */,                  /**< Imp: The protocol                     */
        SeqLim   & /* rSeqLim  */,                  /**< Imp: The sequence limits              */
        SeqExpo  & /* rSeqExpo */,                  /**< Imp: The sequence exports             */
        long       /* lTotalVolumeNo */             /**< Imp: Total number of expected volumes */
        ) {return true;}

    /// Virtual feedback synchronization plugin
    /** \b Input:
        \n rMrProt, rSeqLim, rSeqExpo, lVolumeToBeAcquired, lSliceToBeAcquired

        \b Output:
        \n n.a.

        \b Return value:
        \n true = success, false = error

        Called within ::SyncFeedback if a new and usable feedback is available.
        Can be used to implement all dedicated feedback activities. Should be
        overloaded if required.
    */
    virtual bool SyncFeedbackPlugIn(
        MrProt   & /* rMrProt  */,                  /**< Imp: The protocol (required for SeqUT)                      */
        SeqLim   & /* rSeqLim  */,                  /**< Imp: The sequence limits (required for SeqUT)               */
        SeqExpo  & /* rSeqExpo */,                  /**< Imp: The sequence exports (required for SeqUT)              */
        long       /* lCurrVolume */,               /**< Imp: Current volume - consistent counting with Ice required */
        long       /* lCurrSlice  */                /**< Imp: Current slice - frequency update only once per volume  */
        ) {return true;}

    void runEmptyEventBlock(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) const;

    /// Invalid counter index for initialization purposes.
    const long m_clInvalidVolume{-19222};

    /// Wakeup time for realtime stuff [us]
    long m_lWakeupTime{20000};

    /// Flag to indicate that ::Prep has been called with success
    bool m_bPrepPerformed{false};

    /// Array of flags for each volume indicating whether a feedback has been applied
    std::vector<bool> m_vbFeedbackPerformedOnVolumeNo;

    /// Indicator of a new feedback - semaphore recommended for access
    bool m_bNewFeedbackOccurred{false};

    /// Last volume for which a feedback has been applied
    long m_lLastVolumeWithFeedback{m_clInvalidVolume};

    /// Last volume for which SyncFeedback has been called
    long m_lLastVolumeCalledWithSync{m_clInvalidVolume};

    /// Total number of expected volumes
    long m_lTotalVolumeNo{m_clInvalidVolume};

    /// Synchronization with first (true) or last (false) slice of each volume
    bool m_bSyncWithFirstSlice{true};

    /// Semaphore to synchronise run and receive
    SEQSemaphore m_sFeedbackSemaphore;

    /// Data from last feedback: volume index - semaphore recommended for access
    long   m_lCurrFBVolumeNo{m_clInvalidVolume};
    /// Same as before, but safe for access without semaphore
    long   m_lLocalCurrFBVolumeNo{m_clInvalidVolume};

    /// Data from last feedback: actual feedback data - semaphore recommended for access
    ep_Feedback_Data m_sCurrFBData;
    /// Same as before, but safe for access without semaphore
    ep_Feedback_Data m_sLocalCurrFBData;

    /// Trace feedback
    bool   m_bDumpFB{true};

    /// Synchronization event
    sSYNC_WAKEUP m_sWakeUp;
};


// ===========================================================================
// The constructor initializes the member variables
// ===========================================================================
template<class ep_Feedback_Data>
ep_Feedback<ep_Feedback_Data>::ep_Feedback() : m_sCurrFBData(), m_sLocalCurrFBData()
{
}

// ===========================================================================
// Prepare instance
// ===========================================================================
template<class ep_Feedback_Data>
bool ep_Feedback<ep_Feedback_Data>::Prep(
    MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lTotalVolumeNo, bool bSyncWithFirstSlice)
{
    MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback) << __FUNCTION__ << " called";

    // Init curr FB data
    m_lCurrFBVolumeNo      = m_clInvalidVolume;
    m_lLocalCurrFBVolumeNo = m_clInvalidVolume;

    // flags ...
    m_bPrepPerformed            = false;
    m_lLastVolumeWithFeedback   = m_clInvalidVolume;
    m_lLastVolumeCalledWithSync = m_clInvalidVolume;
    m_bNewFeedbackOccurred      = false;

    m_lTotalVolumeNo      = std::max<long>(1L, lTotalVolumeNo);
    m_bSyncWithFirstSlice = bSyncWithFirstSlice;

    // Init array for feedback security control
    m_vbFeedbackPerformedOnVolumeNo.clear();
    m_vbFeedbackPerformedOnVolumeNo.resize(m_lTotalVolumeNo, false);

    m_vbFeedbackPerformedOnVolumeNo[0] = true; // first/reference volume - gets no feedback
    // Note: If the corresponding Ice functor does not send a feedback event for
    //       the first (reference) volume, it is required that
    //       m_vbFeedbackPerformedOnVolumeNo[1] is also set to true. This preferably
    //       takes place within ::PrepPlugIn of the actual implementation.

    // Call plugin: dedicated preparations of derived class
    const bool bSuccess = PrepPlugIn(rMrProt, rSeqLim, rSeqExpo, lTotalVolumeNo);

    // Remember that we were here
    m_bPrepPerformed = true;

    return bSuccess;
}

// ===========================================================================
// Feedback synchronization and bookkeeping
// ===========================================================================
template<class ep_Feedback_Data>
bool ep_Feedback<ep_Feedback_Data>::SyncFeedback(
    MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lCurrVolume, long lCurrSlice, bool bRunSyncEvent)
{
    bool bSuccess = true;

    MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback) << __FUNCTION__ << " called";

    // Error handling
    if (m_bPrepPerformed == false)
    {
        SEQ_TRACE_INFO.print("Error: called before Prep.");
        return false;
    }

    //-------------------------------------
    // Check if sync event block must be run
    //-------------------------------------

    // Synchronization with first slice of each volume?
    long lSyncSliceNo = 0;
    if (!m_bSyncWithFirstSlice)
    {
        // pointer for MrProtFacade for easier protocol queries
        MrProtFacade protFacade(rMrProt);

        if (protFacade.isSliceAcceleration())
        {
            lSyncSliceNo = SMSProperties::getNReducedSlices(rMrProt) - 1;
        }
        else
        {
            // No - synchronization with last slice of each volume
            lSyncSliceNo = rMrProt.sliceSeries().getlSize() - 1;
        }
    }

    // Execute following event block only for the desired synchronization slice
    if (lCurrSlice != lSyncSliceNo)
    {
        MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback)
            .print("%s called. Volume=%li , Slice=%li", __FUNCTION__, lCurrVolume, lCurrSlice);
        return true;
    }

    if (bRunSyncEvent)
    {
        m_sWakeUp->setDuration(minEventDuration_us);
        m_sWakeUp->setIdent("epFBsWakeup");

        const long endTimeOfEventBlock_us = m_lWakeupTime - minEventDuration_us;

        //----------------------------------
        // Synchronize sequence with reality
        //----------------------------------
        // Event block for synchronization of MaRS calculation and current time
        fRTEBInit(sROT_MATRIXUnity);
        fRTEI(startTimeOfFeedbackEvent_us, 0, 0, 0, 0, 0, 0, &m_sWakeUp);
        fRTEI(endTimeOfEventBlock_us, 0, 0, 0, 0, 0, 0, 0);

        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRun, 'P', 0, 0, 0, 0);
        fRTEBFinish();

        // Now wait until we are really here (in realtime)
        fRTWaitForWakeup();
    }

    // We only want to incorporate feedback once per measurement
    // but we need the upper event block for every acquisition to maintain TR.
    // => Check if we have already been called for the current volume
    // => Check if for the current volume data will be acquired at all
    if ((lCurrVolume == m_lLastVolumeCalledWithSync) || (lCurrVolume < 0))
    {
        MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback)
            .print("%s called. Volume=%li , Slice=%li", __FUNCTION__, lCurrVolume, lCurrSlice);

        runEmptyEventBlock(rMrProt, rSeqLim, rSeqExpo);

        return true;
    }

    // Update
    m_lLastVolumeCalledWithSync = lCurrVolume;

    //--------------------------------------------
    // Process feedback data: administrative stuff
    //--------------------------------------------

    // Reserve semaphore: we want to be sure that no further feedback will interrupt us now
    m_sFeedbackSemaphore.acquire(30);

    MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback).print("CurrFBVolumeNo=%li", m_lCurrFBVolumeNo);
    MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback)
        .print(
            "NewFeedbackOccurred=%i, (VolumeToBeAcquired=%li > LastVolumeWithFeedback=%li)?",
            m_bNewFeedbackOccurred,
            lCurrVolume,
            m_lLastVolumeWithFeedback);
    if ((m_lCurrFBVolumeNo >= 0) && (m_lCurrFBVolumeNo < m_lTotalVolumeNo))
    {
        MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback)
            .print(
                "FeedbackPerformedOnVolumeNo[m_lCurrFBVolumeNo]=%i",
                static_cast<int>(m_vbFeedbackPerformedOnVolumeNo[m_lCurrFBVolumeNo]));
    }

    // Is there an interesting (new + usable) feedback ?
    if ((m_bNewFeedbackOccurred == true)             // Is there a NEW feedback available ?
        && (lCurrVolume > m_lLastVolumeWithFeedback) // Security against messages overtaking each other
        && (m_vbFeedbackPerformedOnVolumeNo[m_lCurrFBVolumeNo]
            == true) // Accept feedback only if basing image data was itself affected by a feedback
    )
    {
        // So there is a feedback which we want to use. Update status data.
        m_bNewFeedbackOccurred    = false;
        m_lLastVolumeWithFeedback = lCurrVolume;

        if (m_bSyncWithFirstSlice)
        {
            // Feedback gets applied to current volume
            if (((unsigned long)lCurrVolume) < m_vbFeedbackPerformedOnVolumeNo.size())
                m_vbFeedbackPerformedOnVolumeNo[lCurrVolume] = true;
        }
        else
        {
            // Feedback gets applied to next volume
            if(((unsigned long)lCurrVolume + 1) < m_vbFeedbackPerformedOnVolumeNo.size())
                m_vbFeedbackPerformedOnVolumeNo[lCurrVolume + 1] = true;
        }

        // Copy feedback data to local variables
        m_lLocalCurrFBVolumeNo = m_lCurrFBVolumeNo;
        m_sLocalCurrFBData     = m_sCurrFBData;

        // Release semaphore: data is now safe in local variables so we can give up the sema
        m_sFeedbackSemaphore.release();

        // Dump feedback info if required
        if (m_bDumpFB == true)
        {
            if (m_bSyncWithFirstSlice)
            {
                // Feedback gets applied to current volume
                SEQ_TRACE_INFO.print(
                    "Received a valid feedback [%li]. Synchronization with volume no.: %li",
                    m_lLocalCurrFBVolumeNo,
                    lCurrVolume);
            }
            else
            {
                // Feedback gets applied to next volume
                SEQ_TRACE_INFO.print(
                    "Received a valid feedback [%li]. Synchronization with volume no.: %li",
                    m_lLocalCurrFBVolumeNo,
                    lCurrVolume + 1);
            }
        }

        // Call plugin: dedicated synchronization functionality of derived class
        // There, the minimum event duration also has to be considered
        bSuccess = SyncFeedbackPlugIn(rMrProt, rSeqLim, rSeqExpo, lCurrVolume, lCurrSlice);
    }
    else
    {
        runEmptyEventBlock(rMrProt, rSeqLim, rSeqExpo);

        // Release semaphore in case of not using feedback
        m_sFeedbackSemaphore.release();

        MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback)
            .print("Feedback NOT used. Volume = %li, Slice = %li", lCurrVolume, lCurrSlice);
    }

    //--------------------------------
    // Finished if we reach this point
    //--------------------------------
    return bSuccess;
}

// ===========================================================================
// Stores current FB data variables if they are acceptable
// ===========================================================================
template<class ep_Feedback_Data>
void ep_Feedback<ep_Feedback_Data>::StoreCurrFBData(long lCurrFBVolumeNo, ep_Feedback_Data sCurrFBData)
{
    MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback) << __FUNCTION__ << " called";

    // Accept feedback only if preparation has been successful
    if (!m_bPrepPerformed)
    {
        SEQ_TRACE_INFO.print("Feedback data not set (m_bPrepPerformed = %i)", m_bPrepPerformed);
        return;
    }

    // Accept only meaningful volume indices
    if (lCurrFBVolumeNo < 0)
    {
        SEQ_TRACE_INFO.print("Feedback data not set (lCurrFBVolumeNo %li < 0)", lCurrFBVolumeNo);
        return;
    }

    // Handling of unexpected additional volumes
    if (lCurrFBVolumeNo >= m_lTotalVolumeNo)
    {
        // We did not expect this volume - extend storage range correspondingly
        SEQ_TRACE_INFO.print(
            "Unexpected volume index %li >= %li. Cannot store feedback data.", lCurrFBVolumeNo, m_lTotalVolumeNo);

        return;
    }

    // Accept feedback only if basing imagedata was itself affected by a feedback !
    // Note: no outputs here as this function probably runs with acquired semaphore
    if (m_vbFeedbackPerformedOnVolumeNo[lCurrFBVolumeNo] == true)
    {
        // Feedback o.k. -> copy received data to class members
        m_lCurrFBVolumeNo      = lCurrFBVolumeNo;
        m_sCurrFBData          = sCurrFBData;
        m_bNewFeedbackOccurred = true;
    }
}

// ===========================================================================
// Access current feedback data
// ===========================================================================
template<class ep_Feedback_Data>
void ep_Feedback<ep_Feedback_Data>::GetCurrFBData(ep_Feedback_Data& sCurrFBData)
{
    // Access most recent feedback data
    sCurrFBData = m_sCurrFBData;

    return;
}

// ===========================================================================
// Provide pointer to semaphore
// ===========================================================================
template<class ep_Feedback_Data>
SEQSemaphore* ep_Feedback<ep_Feedback_Data>::GetFeedbackSemaphore()
{
    return (&m_sFeedbackSemaphore);
}

// ===========================================================================
// Get duration of feedback real time events
// ===========================================================================
template<class ep_Feedback_Data>
long ep_Feedback<ep_Feedback_Data>::getWakeUpDuration(bool bRunSyncEvent)
{
    if (bRunSyncEvent)
    {
        return m_lWakeupTime;
    }

    return minEventDuration_us;
}

// ===========================================================================
/// Set total duration of feedback real time events (wake up)
// ===========================================================================
template<class ep_Feedback_Data>
long ep_Feedback<ep_Feedback_Data>::setWakeUpDuration(long lWakeupTime)
{
    // Headroom of minEventDuration_us is needed for both
    //  - FreqPhase update object in FreqFeedback_SAFETY, and
    //  - to play out the feedback sync event in SyncFeedback
    const long minimumTimeNeeded = startTimeOfFeedbackEvent_us + 2 * minEventDuration_us;

    if (lWakeupTime >= minimumTimeNeeded)
    {
        m_lWakeupTime = lWakeupTime;
    }
    else
    {
        SEQ_TRACE_WARN << "Specified wakeup time of " << lWakeupTime
                       << " us is too short. Setting to minimum duration of " << minimumTimeNeeded << " us.";
        m_lWakeupTime = minimumTimeNeeded;
    }

    return m_lWakeupTime;
}

template<class ep_Feedback_Data>
void ep_Feedback<ep_Feedback_Data>::runEmptyEventBlock(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) const
{
    // dummy RTEB to account for possible freq phase adaptations in FreqFeedback
    fRTEBInit(sROT_MATRIXUnity);
    fRTEI(minEventDuration_us, 0, 0, 0, 0, 0, 0, 0);
    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRun, 'P', 0, 0, 0, 0);
    fRTEBFinish();
}
