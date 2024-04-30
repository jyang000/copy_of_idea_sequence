//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2010  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\FreqFeedback.h
//	 Version:
//	  Author: PLM AW Neuro
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: Class for volume-by-volume synchronization of frequency feedback.
//
//	 Classes: FreqFeedbackBase, FreqFeedback
//
//	-----------------------------------------------------------------------------

/// Double include protection:
#ifndef FreqFeedback_h
#define FreqFeedback_h 1

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------
#include "MrImaging/seq/a_ep2d_diff/a_ep_Feedback.h"
                          
// ===========================================================================
/*!
\class FreqFeedback

\brief Uses synchronization mechanism of ep_Feedback template
       class to implement a volume-by-volume frequency feedback.

*/
// ===========================================================================

/// Frequency feedback base class
//  Feedback data consists of a single double value (frequency offset [Hz])
typedef class ep_Feedback<double> FreqFeedbackBase;

class FreqFeedback : public FreqFeedbackBase
{
public:
    // ------------------------------------
    // Public methods
    // ------------------------------------

    // ------------------------------------
    // Constructors and deconstructors
    // ------------------------------------

    /// The constructor initializes the member variables.
    FreqFeedback();

    /// The destructor does nothing.
    virtual ~FreqFeedback(){};

    /// Get duration of feedback real time events
    /** \b Input:
        \n n.a.

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

    // ---------------------------------------
    // Protected methods and member variables
    // ---------------------------------------
protected:

    /// Prepare plugin
    /** \b Input:
        \n lTotalVolumeNo

        \b Output:
        \n n.a.

        \b Return value:
        \n true = success, false = error

        Called at the end of ::Prep.
        Set initial feedback frequency to zero.
    */
    virtual bool PrepPlugIn(
        MrProt   &rMrProt,                  /**< Imp: The protocol                     */
        SeqLim   &rSeqLim,                  /**< Imp: The sequence limits              */
        SeqExpo  &rSeqExpo,                 /**< Imp: The sequence exports             */
        long      lTotalVolumeNo            /**< Imp: Total number of expected volumes */     
        );

    /// Feedback synchronization plugin
    /** \b Input:
        \n rMrProt, rSeqLim, rSeqExpo, lVolumeToBeAcquired, lSliceToBeAcquired

        \b Output:
        \n n.a.

        \b Return value:
        \n true = success, false = error

        Called within ::SyncFeedback if a new and usable feedback is available. 
        Apply new frequency offset.
    */
    virtual bool SyncFeedbackPlugIn( 
        MrProt   &rMrProt,                  /**< Imp: The protocol (required for SeqUT)                      */
        SeqLim   &rSeqLim,                  /**< Imp: The sequence limits (required for SeqUT)               */
        SeqExpo  &rSeqExpo,                 /**< Imp: The sequence exports (required for SeqUT)              */
        long      lCurrVolume,              /**< Imp: Current volume - consistent counting with Ice required */
        long      lCurrSlice                /**< Imp: Current slice - frequency update only once per volume  */
        );

    double m_dAbsFreqOffset;
    double m_dPrevAbsFreqOffset;
    double m_dPrevAbsFreqOffsetFiltered;

};

#endif  // #ifndef FreqFeedback_h

