//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2010  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\FreqFeedback.cpp
//	 Version: \main\3
//	  Author: PLM AW Neuro
//	    Date: 2011-02-08 12:40:53 +01:00
//
//	    Lang: C++
//
//	 Descrip: Class for volume-by-volume synchronization of frequency feedback.
//
//	 Classes: FreqFeedbackBase, FreqFeedback
//
// EGA Requirement Key: As shown on the following lines:
//
//   Abbrev.   Translation                                        Relevant for
//   -------   -----------                                        ------------
//   EGA-All   One, some or all of the following keys:            Any EGA requirement
//   EGA-01    {:IMPLEMENT:000_EGA_BildOri_SW_SequenzROVz::}      GR/GP   polarity
//   EGA-02    {:IMPLEMENT:000_EGA_BildPos_SW_SequenzSSelVz::}    GS      polarity
//   EGA-03    {:IMPLEMENT:000_EGA_BildMass_SW_SequenzROPC::}     GR/GP   amplitude
//   EGA-04    {:IMPLEMENT:000_EGA_BildPos_SW_SequenzSSel::}      GS      amplitude
//   EGA-05    {:IMPLEMENT:000_EGA_BildPos_SW_NCOFrequenzSSel::}  SRF     frequency
//   EGA-06    {:IMPLEMENT:000_EGA_BildPos_SW_NCOFrequenzRO::}    Readout frequency
//   EGA-07    {:IMPLEMENT:000_EGA_BildOri_SW_OrientierungTest::} Image orientation
//
//	-----------------------------------------------------------------------------

#include "MrImaging/seq/a_ep2d_diff/FreqFeedback.h"
#include "MrMeasSrv/SeqIF/libRT/sFREQ_PHASE.h"        // sFREQ_PHASE
#include "MrImagingFW/libSBBFW/SliceAdjSupport.h"       // ControlParameterBox
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

// ===========================================================================
// The constructor initializes the member variables
// ===========================================================================
FreqFeedback::FreqFeedback() 
:  FreqFeedbackBase()
, m_dAbsFreqOffset(0.)
, m_dPrevAbsFreqOffset(0.)
, m_dPrevAbsFreqOffsetFiltered(0.)
{
}

// ===========================================================================
// Prepare plugin
// ===========================================================================
bool FreqFeedback::PrepPlugIn( MrProt & /* rMrProt */, SeqLim & /* rSeqLim */, SeqExpo & /* rSeqExpo */, long /* lTotalVolumeNo */ )
{
    // Initialize feedback data (frequency offset)
    m_sCurrFBData      = 0.;
    m_sLocalCurrFBData = 0.;

    m_dAbsFreqOffset             = 0.;
    m_dPrevAbsFreqOffset         = 0.;
    m_dPrevAbsFreqOffsetFiltered = 0.;

    return true;
}

// ===========================================================================
// Feedback synchronization plugin
// ===========================================================================
bool FreqFeedback::SyncFeedbackPlugIn( MrProt & /* rMrProt */, SeqLim & /* rSeqLim */, SeqExpo & /* rSeqExpo */, long /* lCurrVolume */, long /* lCurrSlice */ )
{
    sFREQ_PHASE sFreqPhase;
    double      dRelFreqOffset = 0.;

    // Get relative frequency offset from most recent feedback
    GetCurrFBData ( dRelFreqOffset );

    // Keep track of absolute frequency offset
    // Note: This cannot be done in Ice, since only the sequence knows
    //       which frequency offset has been actually applied. Since
    //       feedback takes place asynchronously, it's not guaranteed
    //       that each feedback is considered.
    m_dAbsFreqOffset += dRelFreqOffset;

    // Simple infinite impulse response (IIR) filter
    static const double a1 = 0.726542528005;
    static const double b0 = 0.136728735997;
    static const double b1 = 0.136728735997;

    double dAbsFreqOffsetFiltered =   
          a1 * m_dPrevAbsFreqOffsetFiltered
        + b0 * m_dAbsFreqOffset 
        + b1 * m_dPrevAbsFreqOffset;

    m_dPrevAbsFreqOffset         = m_dAbsFreqOffset;
    m_dPrevAbsFreqOffsetFiltered = dAbsFreqOffsetFiltered;

    MARKER_SEQ_TRACE_DEBUG(SeqTraceMarker_Feedback)
        .print(
            "Relative offset from feedback %f Hz\n"
            "Accumulated absolute offset   %f Hz\n"
            "Applied filtered offset       %f Hz",
            dRelFreqOffset,
            m_dAbsFreqOffset,
            dAbsFreqOffsetFiltered);

    // Update system frequency (relative to original frequency adjustment)
    // Note: The frequency offset gets activated only within the next
    //       'real' FreqPhase event. A minimum event block duration is required

    sFreqPhase.setFrequency( dAbsFreqOffsetFiltered, FP_FREQ_MODE_OFFS_SEQ );       /*! EGA-05; EGA-06 !*/

    fRTEBInit( sROT_MATRIXUnity );
    fRTEI   ( 0,                    sFreqPhase    );                                                  /*! EGA-05; EGA-06 !*/
    fRTEI   ( minEventDuration_us,  0, 0, 0, 0, 0, 0, 0);
    fRTEBFinish();

    // Consider frequency offset in dynamic adjustments
    SLICEADJ::ControlParameterBox::setFrequencyCorrection( static_cast<long>( dAbsFreqOffsetFiltered ) );

    return true;
}

// ===========================================================================
// Get duration of feedback real time events
// ===========================================================================
inline long FreqFeedback::getWakeUpDuration( bool bRunSyncEvent )
{
    // No additional time required
    return FreqFeedbackBase::getWakeUpDuration( bRunSyncEvent );
}
