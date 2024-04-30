//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 2013  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4\pkg\MrServers\MrImaging\seq\common\SeqLoopEPI\SeqLoopEP2D.h
//
//     Lang: C++
//     Authors: Virtual Neuro Team:
//					Dingxin (Guilong) Wang
//					Himanshu Bhat
//					Thomas Beck
//                  Uvo Hoelscher
//                  Mario Zeller
//
//     Descrip: This file is a merge of the SeqLoopEP2D and the SeqLoopMultiBand functionality
//
//              For multi-band imaging:
//              - runOuterLoop handles switching between single-band (preparation) and multi-band scans
//              - Slices indices are re-calculated to reflect the grouping of multi-band slices
//
//     Classes: SeqLoopEP2D
//
//    -----------------------------------------------------------------------------

#pragma once

//------------------------------------------------------------
// Debug
//------------------------------------------------------------
#ifdef DEBUG_ORIGIN
#undef DEBUG_ORIGIN
#endif
#define DEBUG_ORIGIN DEBUG_SEQLOOP

#include "MrImaging/libSBB/SEQLoop.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

//-------------------------------------------------------------------------------------
// CHARM 305893 : TR for single shot EPI sequences should always be physically correct!
//
// This requires new functionality of SeqLoop which is realized in a class derived from
// SeqLoop. The interface of SeqLoop was adapded so that we are able to realize the
// following things within this module:
//
// - TR check of libSeqUT is performed also between two measurements to guarantee that
//   the modifications made are also tested by the sequence UT.
// - The time for SUBFINI / SUBSTRT between two repetitions is taken into account for
//   TR and measuerement time calculations.
// - The time for a halt in the case of single-shot triggering is now correctly handled
//   within the TR-period concerning TR and measurement time calculations.
// - SeqLoop can be asked not to distribute the slices evenly over a TR-period but use
//   a certain TR-fill at the end of each concatenation. This time period can be
//   specified with a new protocol parameter. Choosing the  maximum value of this time
//   period for a given TR results in the measurement of all slices without any TR-fill
//   times between them.
//-------------------------------------------------------------------------------------
//
// CHARM 354879 : Adaptive IR slice thickness for SE and DIFF variants
//
// If an inversion pulse is used with diffusion epi, the doubled inversion thickness
// can generate crosstalk issues (e.g. STIR: imperfect fat suppression) if SeqLoop
// dictates a nesting of inversions (e.g. invert slice 1 - invert slice 3 - scan
// slice 1 - scan slice 3). Thus, an inversion thickness identical to the slice
// thickness seems reasonable.
// However, an inversion pulse is also used for fluid attenuation (FLAIR). Due to
// inflow effects, an increased inversion thickness significantly enhances the quality
// of fluid attenuation.
// As a compromise, a doubled inversion thickness will be only used if the inversion
// time exceeds 500ms, which is a reasonable indication of a FLAIR protocol.
//-------------------------------------------------------------------------------------
//
// SeqLoopLongTRTrig:
// SeqLoopEP2D was previously derived from SeqLoopLongTRTrig, until the creation of 
// SeqLoopMultiBand. This SeqLoopLongTRTrig class supports the possibility
// of activating a 'Long TR Triggering' mode, in which the slices are organized into
// different slice groups, which are acquired during different RR intervals. The number
// of slice groups (and hence the number of RR intervals per volume) is set using the
// number of concatenations parameter. Long TR Triggering mode makes it possible to
// use a long TR and still acquire all slices during diastole. This can be useful
// in diffusion imaging of the brain to avoid artifacts associated with CSF pulsation.
// Note, that navigator triggering and Long TR Triggering mode
// cannot be used at the same time. Long TR triggering mode is not used, if navigator
// triggering or Multi-breath-hold mode is selected.
// Long TR triggering mode is activated for some sequence types in fSEQPrep().
//-------------------------------------------------------------------------------------

// after the inclusion of multiband acquisition, SeqLoopEP2D was usually derived from
// SeqLoopMultiBand. That class is refactored to be a template to enable instance-specific
// inheritance hierarchy, thus SeqLoopEP2D was also reworked in that way.

namespace SEQ_NAMESPACE
{
template <class Seqloop_BASE_TYPE>
class SeqLoopEP2D : public Seqloop_BASE_TYPE
{
  public:
    SeqLoopEP2D() = default;

    virtual bool runOuterSliceLoop(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC);

    virtual bool   calcFillTimesOnly(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lWantedTR = -1, long lTimeForOSC = -1);
    virtual void   calcMeasurementTimeUsec(MrProt& rMrProt, SeqLim& rSeqLim);
    virtual double getTotalMeasTimeUsec(MrProt& rMrProt, SeqLim& rSeqLim);
    virtual void   doClockInitTRBetweenRepetitions(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lSliceIndex, long lEcho);
    virtual bool   insertTRFillEnd(long lFillTime);
    virtual void   setPutTriggerDelayIntoTR(bool bValue);
    virtual void   setAdditionalTimeForExtraEBinTRUsec(long lAdditionalTimeForExtraEBinTRUsec);

    // Sets multi-band/single runmode
    bool setRunMode(SliceAccelRFRunMode eRunMode) override;

    bool TrTiFillTimes(
        MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lScanTime, long lMultiplier, long lTIMinAdd1, long lTIMinAdd2, long lDummyScanTime, long lDummySBBTime, long* plNegativeFillTime) override;
    bool calcFillTimes(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;
    // --------------------------------------------------------------
    // \brief <b> get FatSat duration
    //
    // - Returns duration of FatSat module [us]
    //   (excluding spoiler gradients)
    // - Required for implicit cooling time calculations
    //  --------------------------------------------------------------
    virtual long getlFSDuration();

    // --------------------------------------------------------------
    // \brief <b> get OptFS duration
    //
    // - Returns duration of OptFatSat module [us]
    //   (excluding spoiler gradients)
    // - Required for implicit cooling time calculations
    //  --------------------------------------------------------------
    virtual long getlOptFSDuration();

    // --------------------------------------------------------------
    // \brief <b> get MTC duration
    //
    // - Returns duration of MTC module [us]
    //   (excluding spoiler gradients)
    // - Required for implicit cooling time calculations
    //  --------------------------------------------------------------
    virtual long getlMTCDuration();

    // --------------------------------------------------------------
    // \brief <b> get RSat duration
    //
    // - Returns duration of RSat module [us]
    //   (excluding spoiler gradients)
    // - Required for implicit cooling time calculations
    //  --------------------------------------------------------------
    virtual long getlRSatDuration(long lIndex);

    inline bool setNumberOfPATPrepScans(long lScans)
    {
        m_lPATPrepScans = lScans;
        return true;
    }

    inline long getNumberOfPATPrepScans() const
    {
        return m_lPATPrepScans;
    }

    inline bool setNumberOfSliceAccelPrepScans(long lScans)
    {
        m_lSliceAccelPrepScans = lScans;
        return true;
    }

    inline long getNumberOfSliceAccelPrepScans() const
    {
        return m_lSliceAccelPrepScans;
    }

    inline bool setNumberOfInitialDummyScans(long lScans)
    {
        // Currently at least 0 dummy scans are needed before actual measurement starts
        if (lScans < 0)
        {
            SEQ_TRACE_ERROR.print("Number of dummy scans (%ld) needs to be at least 0", lScans);
            return false;
        }

        m_lInitialPrepScans = lScans;
        return true;
    }

    inline long getNumberOfInitialDummyScans() const
    {
        return m_lInitialPrepScans;
    }

    inline bool setNumberOfSliceAccelDummyScans(long lScans)
    {
        if (lScans < 0)
        {
            SEQ_TRACE_ERROR.print("Number of dummy scans (%ld) needs to be at least 0", lScans);
            return false;
        }

        m_lSliceAccelDummyScans = lScans;
        return true;
    }

    inline long getNumberOfSliceAccelDummyScans() const
    {
        return m_lSliceAccelDummyScans;
    }

    inline bool setNumberOfPhaseCorrectionScans(long lScans)
    {
        m_lPhaseCorrPrepScans = lScans;
        return true;
    }

    inline long getNumberOfPhaseCorrectionScans() const
    {
        return m_lPhaseCorrPrepScans;
    }
    
    inline bool setNumberOfSliceAccelPhaseCorrectionScans(long lScans)
    {
        m_lSliceAccelPhaseCorrectionScans = lScans;
        return true;
    }

    void setFatSatOffcenterFrequency(long lFrequency);
    // --------------------------------------------------------------
    // \brief <b> get repetitions loop counter
    //
    // - Returns current repetition
    // - Required for setting Mdh entries in diffusion module
    //  --------------------------------------------------------------
    virtual long getRepetitionLoopCounter() const;

#if defined SUPPORT_iPAT_a_ep2d || defined SUPPORT_iPAT_TGSE

    void setlSBBPATRefScanBaseRes(long lBaseRes);
    void setlSBBPATRefScanBandwidth(long lBandwidth);
    void setSkipOnlinePhaseCorrFlag();
    void setSkipRegriddingFlag();

#endif

  protected:
    bool m_bPutTriggerDelayIntoTR{false};
    long m_lAdditionalTimeForExtraEBinTRUsec{0};

    // Number of PAT reference scans
    long m_lPATPrepScans{0};
    // Number of slice acceleration reference scans
    long m_lSliceAccelPrepScans{0};
    // Number of single-band dummy scans
    long m_lInitialPrepScans{1};
    // Number of multi-band dummy scans
    long m_lSliceAccelDummyScans{0};
    // Number of external phase correction scans (applied after initial dummy scans)
    long m_lPhaseCorrPrepScans{0};
    
    long m_lSliceAccelPhaseCorrectionScans{0};

    // Total number of outer loops for single-band mode
    long m_lTotalNumberOfOuterLoopsSingleBand{0};
    // Total number of outer loops for multi-band mode
    long m_lTotalNumberOfOuterLoopsMultiBand{0};
    // Total number of slices for single-band mode
    long m_lSlicesToMeasureSingleBand{0};
    // Total number of slices for single-band mode
    long m_lSlicesToMeasureMultiBand{0};
};

} // namespace SEQ_NAMESPACE




namespace SEQ_NAMESPACE
{
    template <class SeqLoop_BASE_TYPE>
    bool SeqLoopEP2D<SeqLoop_BASE_TYPE>::runOuterSliceLoop(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC)
    {
        //-----------------------------------------------------------------------
        // Modify scan settings in correspondence to run mode
        //-----------------------------------------------------------------------
        if (SMSProperties::isSMS(rMrProt))
        {
            // Calculates number of single-band prep scans
            long lSingleBandPrepScans = m_lPATPrepScans + m_lSliceAccelPrepScans + m_lInitialPrepScans + m_lPhaseCorrPrepScans;

            SeqLoop_BASE_TYPE::m_lSliceOffset = 0;
            SeqLoop_BASE_TYPE::m_lSliceOffsetConc = 0;

            if (SeqLoop_BASE_TYPE::m_lPrepareLoopCounter < lSingleBandPrepScans)
            {
                // settings for single-band reference scan
            setRunMode(SINGLE_BAND);

            }
            else
            {
                // settings for multi-band scan
            setRunMode(MULTI_BAND);

            }

            for (long lCounter = 0; lCounter < SeqLoop_BASE_TYPE::m_lConcatenationCounter; lCounter++)
            {
                SeqLoop_BASE_TYPE::m_lSliceOffset += SeqLoop_BASE_TYPE::m_SeqConcat[lCounter].m_Slices;
                SeqLoop_BASE_TYPE::m_lSliceOffsetConc += SeqLoop_BASE_TYPE::m_SeqConcat[lCounter].m_Slices;
            }
        }
        else
        {
            // settings for single-band scan
            SeqLoop_BASE_TYPE::m_eRunMode = SINGLE_BAND;
        }

        //-----------------------------------------------------------------------
        // call base class
        //-----------------------------------------------------------------------
        if (!SeqLoop_BASE_TYPE::runOuterSliceLoop(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC))
        {
            SEQ_TRACE_ERROR.print("Error encountered in runOuterSliceLoop(...).");
            SeqLoop_BASE_TYPE::setNLSStatus(MRI_SEQ_SEQU_ERROR);
            return false;
        }

        return true;
    }

    template <class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::setPutTriggerDelayIntoTR(bool bValue)
    {
        m_bPutTriggerDelayIntoTR = bValue;
    }

    template <class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::setAdditionalTimeForExtraEBinTRUsec(long lAdditionalTimeForExtraEBinTRUsec)
    {
        m_lAdditionalTimeForExtraEBinTRUsec = lAdditionalTimeForExtraEBinTRUsec;
    }

    template <class SeqLoop_BASE_TYPE>
    bool SeqLoopEP2D<SeqLoop_BASE_TYPE>::calcFillTimesOnly(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lWantedTR_orig, long lTimeForOSC_orig)
    {
        //---------------------------------------------------------------------------
        // we expect that we are called from SeqLoop with the default parameters:
        //---------------------------------------------------------------------------
        if (lWantedTR_orig != -1 || lTimeForOSC_orig != -1)
        {
            // caused by programming error
            // => trace also if rSeqLim.isContextPrepForBinarySearch()
            TEXT_TR(rSeqLim, "unexpected ERROR: lWantedTR_orig!=-1 || lTimeForOSC_orig!=-1")
                SeqLoop_BASE_TYPE::setNLSStatus(MRI_SEQ_SEQU_ERROR);
            return false;
        }

#ifndef TGSE
        //---------------------------------------------------------------------------
        // handling of m_lAdditionalTimeForExtraEBinTRUsec and m_TrigHaltDuration
        // may (or will) not work when multiple phases are selected
        //---------------------------------------------------------------------------
        if (SeqLoop_BASE_TYPE::m_PhasesToMeasure > 1)
        {
            // caused by programming error
            // => trace also if rSeqLim.isContextPrepForBinarySearch()
            TEXT_TR(rSeqLim, "m_PhasesToMeasure>1 is currently not supported")
                SeqLoop_BASE_TYPE::setNLSStatus(MRI_SEQ_SEQU_ERROR);
            return false;
        }
#endif

        //---------------------------------------------------------------------------
        // We want to force SeqLoop::calcFillTimesOnly to use a minimum for TRFillEnd
        // so we reduce the wanted TR for calculation of fill times.
        //
        // For multiple measurements we want at least fSBBMeasRepetDelayGetDurationForZeroPauseUs() to be able
        // to compensate the increase of TR caused by the SUBFINI/SUBSTRT-event-block
        // inserted by SeqLoop. Compensation is done later by simply reducing the last
        // TRFillEnd by fSBBMeasRepetDelayGetDurationForZeroPauseUs(). Therefore TRFilleEnd must have at least this value.
        //
        // Additional times must be put into TR:
        // - pMrProt->delayTimeInTR()
        // - m_lAdditionalTimeForExtraEBinTRUsec
        // - m_TrigHaltDuration
        //---------------------------------------------------------------------------
        long lWantedTR = rMrProt.tr()[0];
        long lReduceActualTRForCalculation = 0;

        lReduceActualTRForCalculation = maximum((long)fSBBMeasRepetDelayGetDurationForZeroPauseUs(), (long)rMrProt.delayTimeInTR());
        lReduceActualTRForCalculation += m_lAdditionalTimeForExtraEBinTRUsec;

        if (m_bPutTriggerDelayIntoTR && !SeqLoop_BASE_TYPE::m_TrigHaltSingleShot)
        {
            lReduceActualTRForCalculation += SeqLoop_BASE_TYPE::m_TrigHaltDuration;
        }

        lWantedTR -= lReduceActualTRForCalculation;

        if (lWantedTR < 0)
        {
            lWantedTR = 0;
        }

        //---------------------------------------------------------------------------
        // By increasing the time for the OSC-bit we force SeqLoop to put the ECG-fill
        // into TR, if m_TrigHaltSingleShot is true.
        //
        // Compare graphs of SeqLoop inner loop structure: fSBBECGFillTimeRun is called
        // directly before SBBOscBitRun.
        //---------------------------------------------------------------------------
        long lTimeForOSC = SeqLoop_BASE_TYPE::getScanTimeOscBit();

        if (m_bPutTriggerDelayIntoTR && SeqLoop_BASE_TYPE::m_TrigHaltSingleShot)
        {
            lTimeForOSC += SeqLoop_BASE_TYPE::getTrigHaltDuration();
        }

        //---------------------------------------------------------------------------
        // Do standard fill time calculation with modified wanted TR and time for OSC
        //---------------------------------------------------------------------------
        bool bRet = false;

        // In case multi-band acquisition is active calculate fill times for single-band reference scans first
        if (SMSProperties::isSMS(rMrProt))
        {
            // Modify number of slices for ::calcFillTimesOnly in SeqLoop
            SeqLoop_BASE_TYPE::m_SlicesToMeasure = SeqLoop_BASE_TYPE::m_SlicesToMeasure * SeqLoop_BASE_TYPE::m_lMultibandFactor;

#ifdef SUPPORT_FAST_IR
            if ((rMrProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
            {
                lWantedTR = rMrProt.tr()[0];
                lReduceActualTRForCalculation = 0;
                lTimeForOSC = 0;
            }
#endif
            bRet = SeqLoop_BASE_TYPE::calcFillTimesOnly(rMrProt, rSeqLim, rSeqExpo, lWantedTR, lTimeForOSC);

            if (!bRet)
            {
                SEQ_TRACE_ERROR.print("ERROR: calcFillTimesOnly failed.");
                return false;
            }

        m_lTotalNumberOfOuterLoopsSingleBand = SeqLoop_BASE_TYPE::m_lTotalNumberOfOuterLoops;
        m_lSlicesToMeasureSingleBand         = SeqLoop_BASE_TYPE::m_SlicesToMeasure;
            SeqLoop_BASE_TYPE::copySeqConcat(SeqLoop_BASE_TYPE::m_SeqConcatSingleBand, SeqLoop_BASE_TYPE::m_SeqConcat);

            for (long lCounter = 0; lCounter < SeqLoop_BASE_TYPE::m_lConcatenations; lCounter++)
            {
                SeqLoop_BASE_TYPE::m_SeqConcatSingleBand[lCounter].m_TRFill = 200;                                                                 // needed for SPAIR to avoid error: ADCToRFTimeTooShort
            SeqLoop_BASE_TYPE::m_SeqConcatSingleBand[lCounter].m_TRFillEnd
                = SeqLoop_BASE_TYPE::m_bShiftTRFillToTRFillEnd ? m_lAdditionalTimeForExtraEBinTRUsec : 0; // Normally 0, special case for triggered SMS protocols;
        }

            // Restore previous number of slices
            SeqLoop_BASE_TYPE::m_SlicesToMeasure
                = SeqLoop_BASE_TYPE::m_SlicesToMeasure / SeqLoop_BASE_TYPE::m_lMultibandFactor;
        }

        // The following section applies for the unaccelerated and the multi-band case
        // In the multi-band case the reduced slice number, in the standard case the full slice number is used here
#ifdef SUPPORT_FAST_IR
        if ((rMrProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
        {
            lWantedTR = rMrProt.tr()[0];
            lReduceActualTRForCalculation = 0;
            lTimeForOSC = 0;
        }
#endif
        bRet = SeqLoop_BASE_TYPE::calcFillTimesOnly(rMrProt, rSeqLim, rSeqExpo, lWantedTR, lTimeForOSC);

        //---------------------------------------------------------------------------
        // correct members due to manipulations of lWantedTR
        //---------------------------------------------------------------------------
        for (long lI = 0; lI < SeqLoop_BASE_TYPE::m_lConcatenations; ++lI)
        {
            SeqLoop_BASE_TYPE::m_SeqConcat[lI].m_TRFillEnd += lReduceActualTRForCalculation;
        }

        SeqLoop_BASE_TYPE::m_lActualTR += lReduceActualTRForCalculation;

        if (SeqLoop_BASE_TYPE::m_bHandleTRTIConflict && SeqLoop_BASE_TYPE::m_lTRneeded != 0)
        {
            SeqLoop_BASE_TYPE::m_lTRneeded += lReduceActualTRForCalculation;
        }

        //---------------------------------------------------------------------------
        // Save m_SeqConcat for multi-band scan
        //---------------------------------------------------------------------------
        if (SMSProperties::isSMS(rMrProt))
        {
            m_lTotalNumberOfOuterLoopsMultiBand = SeqLoop_BASE_TYPE::m_lTotalNumberOfOuterLoops;
            m_lSlicesToMeasureMultiBand         = SeqLoop_BASE_TYPE::m_SlicesToMeasure;
            SeqLoop_BASE_TYPE::copySeqConcat(SeqLoop_BASE_TYPE::m_SeqConcatMultiBand, SeqLoop_BASE_TYPE::m_SeqConcat);
        }

        //---------------------------------------------------------------------------
        // finished
        //---------------------------------------------------------------------------
        return bRet;
}

template<class SeqLoop_BASE_TYPE>
bool SeqLoopEP2D<SeqLoop_BASE_TYPE>::setRunMode(SliceAccelRFRunMode eRunMode)
{
#ifdef SUPPORT_IIR
    // Interleaved schemes do not support switching of the run mode
    if (SeqLoop_BASE_TYPE::IIRScheme() != SeqLoop_BASE_TYPE::PROT_MASK_IIR_SCHEME_STD)
    {
        return true;
    }
#endif
    if (!SeqLoop_BASE_TYPE::setRunMode(eRunMode))
        return false;

    switch (eRunMode)
    {
        case SINGLE_BAND:
            SeqLoop_BASE_TYPE::m_lTotalNumberOfOuterLoops = m_lTotalNumberOfOuterLoopsSingleBand;
            SeqLoop_BASE_TYPE::m_SlicesToMeasure                      = m_lSlicesToMeasureSingleBand;
                break;
        case MULTI_BAND:
            SeqLoop_BASE_TYPE::m_lTotalNumberOfOuterLoops = m_lTotalNumberOfOuterLoopsMultiBand;
            SeqLoop_BASE_TYPE::m_SlicesToMeasure                      = m_lSlicesToMeasureMultiBand;
                break;
        default:
            // Throw error if run mode is unknown
            SEQ_TRACE_ERROR.print("ERROR: RunMode not supported.");
            return false;
    }

    SeqLoop_BASE_TYPE::SBBIRsel.setRunMode(eRunMode);

    return true;
}

template<class SeqLoop_BASE_TYPE>
bool SeqLoopEP2D<SeqLoop_BASE_TYPE>::TrTiFillTimes(
    MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lScanTime, long lMultiplier, long lTIMinAdd1, long lTIMinAdd2, long lDummyScanTime, long lDummySBBTime, long* plNegativeFillTime)
{
    // Call corresponding base class implementation
    bool bResult = SeqLoop_BASE_TYPE::TrTiFillTimes(rMrProt, rSeqLim, rSeqExpo, lScanTime, lMultiplier, lTIMinAdd1, lTIMinAdd2, lDummyScanTime, lDummySBBTime, plNegativeFillTime);

    // With SMS, prohibit SpecialInterleaveMode_INTERLEAVED_IR
    // Note: There have been good reasons for introducing this special interleaving mode. However, so far this
    //       does not consider multi-band specialities, which might e.g. lead to successive inversion of
    //       (virtually) adjacent slices. In the future, one might envision more sophisticated interleaving
    //       patterns which consider multi-band excitation.
    if (SMSProperties::isSMS(rMrProt))
    {
        if (SeqLoop_BASE_TYPE::m_eSpecialSliceInterleaveMode == SpecialInterleaveMode_INTERLEAVED_IR)
        {
            // With SMS, prohibit SpecialInterleaveMode_INTERLEAVED_IR
            // Note: There have been good reasons for introducing this special interleaving mode. However, so far this
            //       does not consider multi-band specialities, which might e.g. lead to successive inversion of
            //       (virtually) adjacent slices. In the future, one might envision more sophisticated interleaving
            //       patterns which consider multi-band excitation.
            SeqLoop_BASE_TYPE::m_eSpecialSliceInterleaveMode = SpecialInterleaveMode_OFF;
        }

#ifdef SUPPORT_IIR
        if (!SeqLoop_BASE_TYPE::IsIIRSchemeStd(rMrProt, rSeqLim))
        {
            // Store concatenation-related information for use in multi-band mode
            m_lTotalNumberOfOuterLoopsMultiBand = SeqLoop_BASE_TYPE::m_lTotalNumberOfOuterLoops;
            m_lSlicesToMeasureMultiBand         = SeqLoop_BASE_TYPE::m_SlicesToMeasure;
            SeqLoop_BASE_TYPE::copySeqConcat(SeqLoop_BASE_TYPE::m_SeqConcatMultiBand, SeqLoop_BASE_TYPE::m_SeqConcat);

            // SeqLoopIIR does not support single-band acquisitions
            m_lTotalNumberOfOuterLoopsSingleBand = 0;
            m_lSlicesToMeasureSingleBand = SeqLoop_BASE_TYPE::m_SlicesToMeasure * SeqLoop_BASE_TYPE::m_lMultibandFactor;
            SeqLoop_BASE_TYPE::copySeqConcat(SeqLoop_BASE_TYPE::m_SeqConcatSingleBand, SeqLoop_BASE_TYPE::m_SeqConcat);
        }
#endif
    }

    return bResult;
}

template<class SeqLoop_BASE_TYPE>
bool SeqLoopEP2D<SeqLoop_BASE_TYPE>::calcFillTimes(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo)
{
    // Call corresponding base class implementation
    bool bResult = SeqLoop_BASE_TYPE::calcFillTimes(rMrProt, rSeqLim, rSeqExpo);

    // With SMS, prohibit SpecialInterleaveMode_INTERLEAVED_IR
    // Note: There have been good reasons for introducing this special interleaving mode. However, so far this
    //       does not consider multi-band specialities, which might e.g. lead to successive inversion of
    //       (virtually) adjacent slices. In the future, one might envision more sophisticated interleaving
    //       patterns which consider multi-band excitation.
    if (SMSProperties::isSMS(rMrProt)
        && (SeqLoop_BASE_TYPE::m_eSpecialSliceInterleaveMode == SpecialInterleaveMode_INTERLEAVED_IR))
    {
        // Switch back to standard interleaving pattern
        SeqLoop_BASE_TYPE::m_eSpecialSliceInterleaveMode = SpecialInterleaveMode_OFF;
    }

    return bResult;
}

    template <class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::calcMeasurementTimeUsec(MrProt& rMrProt, SeqLim& rSeqLim)
    {
        SeqLoop_BASE_TYPE::calcMeasurementTimeUsec(rMrProt, rSeqLim);

        // we took care that ECG-fill-time was put into TR during calcFillTimesOnly,
        // standard SeqLoop-timing calculation does not recognize this, so we have to
        // do a correction here
        if (m_bPutTriggerDelayIntoTR)
        {
            SeqLoop_BASE_TYPE::m_dMeasureTimeInFirstMeasUsec -= SeqLoop_BASE_TYPE::m_dTrigHaltTimeInFirstMeasUsec;
            SeqLoop_BASE_TYPE::m_dTrigHaltTimeInFirstMeasUsec = 0;

            if (SeqLoop_BASE_TYPE::m_RepetitionsToMeasure)
            {
                SeqLoop_BASE_TYPE::m_dMeasureTimeInSecondMeasUsec -= SeqLoop_BASE_TYPE::m_dTrigHaltTimeInSecondMeasUsec;
                SeqLoop_BASE_TYPE::m_dTrigHaltTimeInSecondMeasUsec = 0;
            }
        }
    }

    template <class SeqLoop_BASE_TYPE>
    double SeqLoopEP2D<SeqLoop_BASE_TYPE>::getTotalMeasTimeUsec(MrProt& rMrProt, SeqLim& rSeqLim)
    {
        SeqLoop_BASE_TYPE::getTotalMeasTimeUsec(rMrProt, rSeqLim);

        // we take care that TR is always exact, so we consider the time for
        // SUBFINI / SUBSTRT between two repetitions during TR-calculation and during
        // execution of the sequence timing in calcFillTimesOnly
#ifdef SUPPORT_IIR
        if (SeqLoop_BASE_TYPE::IsIIRSchemeStd(rMrProt, rSeqLim))
        {
            SeqLoop_BASE_TYPE::m_dTotalMeasureTimeUsec
                -= fSBBMeasRepetDelayGetDurationForZeroPauseUs()
                   * static_cast<double>(SeqLoop_BASE_TYPE::m_RepetitionsToMeasure);
        }
#else

        SeqLoop_BASE_TYPE::m_dTotalMeasureTimeUsec -= fSBBMeasRepetDelayGetDurationForZeroPauseUs() * static_cast<double>(SeqLoop_BASE_TYPE::m_RepetitionsToMeasure);
#endif

#if defined ASL && !defined TGSE

        bool bEnableFirstPrepScanAsM0Scan = true;

        if (rMrProt.getlRepetitions() <= 0)
        {
            bEnableFirstPrepScanAsM0Scan = false;
        }

#ifdef SUPPORT_iPAT_a_ep2d

        // with iPAT a separate M0 scan BEFORE the iPAT scan is not possible
        if ((rMrProt.PAT().getucPATMode() != SEQ::PAT_MODE_NONE) && (rMrProt.PAT().getlAccelFactPE() > 1))
        {
            bEnableFirstPrepScanAsM0Scan = false;
        }

#endif

        if (bEnableFirstPrepScanAsM0Scan)
        {
            SeqLoop_BASE_TYPE::m_dTotalMeasureTimeUsec += fSBBMeasRepetDelayGetDurationForZeroPauseUs();
        }
#endif

        return SeqLoop_BASE_TYPE::m_dTotalMeasureTimeUsec;
    }

    template <class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::doClockInitTRBetweenRepetitions(MrProt&, SeqLim&, SeqExpo&, long, long)
    {
        // We never do a clock initialization for TR for the unit test between repetitions,
        // because TR must be physically correct for single shot EPI even over repetitions!
        return;
    }

    template <class SeqLoop_BASE_TYPE>
    bool SeqLoopEP2D<SeqLoop_BASE_TYPE>::insertTRFillEnd(long lFillTime)
    {
        // SeqLoop will insert an event block for SUBFINI/SUBSTRT sync events which will disturb
        // our TR for multiple measurements. So we have to shorten the last TR-fill-end of the
        // measurement to compensate this effect.
        //
        if (SeqLoop_BASE_TYPE::m_bIsLastScanInMeas && SeqLoop_BASE_TYPE::m_lRepetitionCounter < SeqLoop_BASE_TYPE::m_RepetitionsToMeasure && SeqLoop_BASE_TYPE::m_RepetitionsToMeasure)
        {
            lFillTime -= static_cast<long>(fSBBMeasRepetDelayGetDurationForZeroPauseUs());
        }

        // adapt TRFillEnd for extra eventblocks in TR blocks
        if (SeqLoop_BASE_TYPE::m_lRepetitionCounter <= SeqLoop_BASE_TYPE::m_RepetitionsToMeasure)
        {
            lFillTime -= m_lAdditionalTimeForExtraEBinTRUsec;
        }

        // adapt TRFillEnd due to Trig-Halt
        if (m_bPutTriggerDelayIntoTR && !SeqLoop_BASE_TYPE::m_TrigHaltSingleShot)
        {
            lFillTime -= SeqLoop_BASE_TYPE::m_TrigHaltDuration;
        }

        // check fill time:
        if (lFillTime < 0)
        {
            // caused by programming error
            // => trace also if rSeqLim.isContextPrepForBinarySearch()
            TEXT_TR(rSeqLim, "FATAL: lFillTime<0")
                SeqLoop_BASE_TYPE::setNLSStatus(MRI_SEQ_SEQU_ERROR);
            return false;
        }

        // setNLSStatus returns false for success
        return !SeqLoop_BASE_TYPE::setNLSStatus(fSBBFillTimeRun(lFillTime));
    }


    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::getlFSDuration
    // ------------------------------------------------------------------------------
    //
    // Description : return duration of FatSat module
    //
    // Return      : FatSsat duration [us]
    //
    // ------------------------------------------------------------------------------
    template <class SeqLoop_BASE_TYPE>
    long SeqLoopEP2D<SeqLoop_BASE_TYPE>::getlFSDuration(void)
    {
        return SeqLoop_BASE_TYPE::CSatFat.getDurationPerRequest();
    }

    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::getlOptFSDuration
    // ------------------------------------------------------------------------------
    //
    // Description : return duration of OptFatSat module
    //
    // Return      : OptFatSat duration [us]
    //
    // ------------------------------------------------------------------------------
    template <class SeqLoop_BASE_TYPE>
    long SeqLoopEP2D<SeqLoop_BASE_TYPE>::getlOptFSDuration()
    {
        return SeqLoop_BASE_TYPE::m_SBBOptfs.getRFDuration() + SeqLoop_BASE_TYPE::m_SBBOptfs.getDwellTime();
    }

    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::getlMTCDuration
    // ------------------------------------------------------------------------------
    //
    // Description : return duration of MTC module
    //
    // Return      : MTC duration [us]
    //
    // ------------------------------------------------------------------------------
    template<class SeqLoop_BASE_TYPE>
    long SeqLoopEP2D<SeqLoop_BASE_TYPE>::getlMTCDuration()
    {
        return SeqLoop_BASE_TYPE::MSat.getDurationPerRequest();
    }

    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::getlRSatDuration
    // ------------------------------------------------------------------------------
    //
    // Description : return duration of RSat module of given index
    //
    // Return      : RSat duration [us]
    //
    // ------------------------------------------------------------------------------
    template <class SeqLoop_BASE_TYPE>
    long SeqLoopEP2D<SeqLoop_BASE_TYPE>::getlRSatDuration(long lIndex)
    {
        return SeqLoop_BASE_TYPE::RSat[lIndex]->getDurationPerRequest();
    }

    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::getRepetitionLoopCounter
    // ------------------------------------------------------------------------------
    //
    // Description : return current repetition loop counter
    //
    // Return      : repetition loop counter
    //
    // ------------------------------------------------------------------------------
    template <class SeqLoop_BASE_TYPE>
    long SeqLoopEP2D<SeqLoop_BASE_TYPE>::getRepetitionLoopCounter(void) const
    {
        return SeqLoop_BASE_TYPE::m_lRepetitionCounter;
    }


    template<class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::setFatSatOffcenterFrequency(long lFrequency)
    {
        SeqLoop_BASE_TYPE::CSatFat.setUserDefinedOffsetHz(lFrequency);
    }

#if defined SUPPORT_iPAT_a_ep2d || defined SUPPORT_iPAT_TGSE

    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::setlSBBPATRefScanBaseRes
    // ------------------------------------------------------------------------------
    //
    // Description : set base resolution of external iPAT reference scans
    //
    // Return      : n.a.
    //
    // ------------------------------------------------------------------------------
    template <class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::setlSBBPATRefScanBaseRes(long lBaseRes)
    {
        SeqLoop_BASE_TYPE::SBBPATRefScan.getGRERefScan()->setlBaseRes(lBaseRes);
    }

    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::setlSBBPATRefScanBandwidth
    // ------------------------------------------------------------------------------
    //
    // Description : set bandwidth of external iPAT reference scans
    //
    // Return      : n.a.
    //
    // ------------------------------------------------------------------------------
    template <class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::setlSBBPATRefScanBandwidth(long lBandwidth)
    {
        SeqLoop_BASE_TYPE::SBBPATRefScan.getGRERefScan()->setlBandwidth(lBandwidth);
    }

    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::setSkipOnlinePhaseCorrFlag
    // ------------------------------------------------------------------------------
    //
    // Description : add  MDH_SKIP_ONLINE_PHASCOR flag to MDH
    //
    // Return      : n.a.
    //
    // ------------------------------------------------------------------------------
    template <class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::setSkipOnlinePhaseCorrFlag()
    {
        SeqLoop_BASE_TYPE::SBBPATRefScan.getGRERefScan()->getADC()->getMDH().addToEvalInfoMask(MDH_SKIP_ONLINE_PHASCOR);
    }

    // ------------------------------------------------------------------------------
    // Name      : SeqLoopEP2D<SeqLoop_BASE_TYPE>::setSkipRegriddingFlag
    // ------------------------------------------------------------------------------
    //
    // Description : add  MDH_SKIP_REGRIDDING flag to MDH
    //
    // Return      : n.a.
    //
    // ------------------------------------------------------------------------------
    template <class SeqLoop_BASE_TYPE>
    void SeqLoopEP2D<SeqLoop_BASE_TYPE>::setSkipRegriddingFlag()
    {
        SeqLoop_BASE_TYPE::SBBPATRefScan.getGRERefScan()->getADC()->getMDH().addToEvalInfoMask(MDH_SKIP_REGRIDDING);
    }

#endif

} // namespace SEQ_NAMESPACE
