//----------------------------------------------------------------------------------
// <copyright file="SeqLoopLongTRTrig.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 2008-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015. All Rights Reserved. Confidential.
// </copyright>
// <description>
//   Modified SeqLoop with modified triggering mode for long TR.
//
//   Slices are grouped so that each group is executed during a
//     different RR-interval, allowing all slices to be acquired during diastole.
//
//   This is realized by repositioning the concatenations loop.
// </description>
//----------------------------------------------------------------------------------

#include "MrImaging/seq/common/SeqLoopLongTRTrig/SeqLoopLongTRTrig.h"

#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrProtSrv/Domain/CoreNative/MrNavigator.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Physiology/MrPhysiology.h"

//------------------------------------------------------------
// Debug
//------------------------------------------------------------
#define DEBUG_ORIGIN DEBUG_SEQLOOP

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif
using namespace SEQ_NAMESPACE;

//-----------------------------------------------------------------------
// function:    runConcatenationLoop()
//
// description: Overloaded SeqLoop function.
//              Standard concatenations loop fixed to a single iteration.
//              Concatenations now performed outside 'outer slices' loop.
//-----------------------------------------------------------------------
bool SeqLoopLongTRTrig::runConcatenationLoop(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC)
{
    NLS_STATUS lStatus = MRI_SBB_SBB_NORMAL;

    // * -------------------------------------------------------------------------- *
    // * Run base class version if long TR triggering mode is switched off          *
    // * or if navigator triggering is active                                       *
    // * -------------------------------------------------------------------------- *
    if (!m_bLongTRTrigMode || (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF))
    {
        return SeqLoop_BASE_TYPE::runConcatenationLoop(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC);
    }

    // * -------------------------------------------------------------------------- *
    // initialization
    // * -------------------------------------------------------------------------- *
    m_ADCCounter            = 0;
    m_lConcatenationCounter = m_lSliceOffsetConc = 0;

    // * -------------------------------------------------------------------------- *
    // * Old concatenations loop                                                    *
    // * Now always run with a single iteration                                     *
    // * A new concatenations loop is exectuted in runOuterSliceLoop()              *
    // * -------------------------------------------------------------------------- *
    {
        for (m_lConcatenationCounter = m_lSliceOffsetConc = 0; m_lConcatenationCounter < m_lConcatenations; m_lConcatenationCounter++) // charm453626
        {
            // * ---------------------------------------------------------------------- *
            // * run SBBPATRefScan for current concatenation                            *
            // *                                                                        *
            // * Not used by current EPI sequence. If this is required in the future    *
            // * a modified version of runPATRefScanConcatenated() will be required for *
            // * compatability with the new concatenations loop.                        *
            // * ---------------------------------------------------------------------- *
            if ((m_ePATRefScanLoopMode == PATRefScanLoopMode_EACH_REPETION) || ((m_ePATRefScanLoopMode == PATRefScanLoopMode_ONCE) && (m_lRepetitionCounter == 0)))
            {
#ifdef SHOW_LOOP_STRUCTURE
                if (IS_TR_LAND(rSeqLim))
                {
                    ShowEvent((char*)"SBBPATRefScan (concatenated)");
                }
#endif
                if (!runPATRefScanConcatenated(rMrProt, rSeqLim, rSeqExpo, pSlcPos))
                {
                    if (!rSeqLim.isContextPrepForBinarySearch())
                    {
                        SEQ_TRACE_ALWAYS.print("Error encountered in running runPATRefScanConcatenated()");
                    }
                    return (false); // NLS Status was set in runPATRefScanConcatenated()
                }
            }
            m_lSliceOffsetConc += m_SeqConcat[m_lConcatenationCounter].m_Slices; // charm453626
        }

        // * ---------------------------------------------------------------------- *
        // * Initialize counter for quick sats and triggering                       *
        // * ---------------------------------------------------------------------- *
        m_RawLinesCount                = 0;
        m_CountFirstScansForTriggering = 0;

        // reset loop counters
        m_FreeLoopCounter          = 0;
        m_lOuterAcquisitionCounter = 0;
        m_lLineCounter             = 0;
        m_lPartitionCounter        = 0;
        m_lInnerAcquisitionCounter = 0;

        // * ---------------------------------------------------------------------- *
        // * Reset mdh index values and control flags for actual concatenation      *
        // * ---------------------------------------------------------------------- *
        psADC->getMDH().setCacq(0);
        psADC->getMDH().setClin(0);
        psADC->getMDH().setCpar(0);
        psADC->getMDH().setFirstScanInSlice(m_bIsFirstScanInSlice = false);
        psADC->getMDH().setLastScanInSlice(m_bIsLastScanInSlice = false);

        if (IS_TR_INT(rSeqLim))
        {
            SEQ_TRACE_ALWAYS.print(
                "Concatenation %ld/%ld: %ld slices; %ld outer loops; %ld (max.) inner loops\n",
                m_lConcatenationCounter + 1,
                m_lConcatenations,
                m_SeqConcat[m_lConcatenationCounter].m_Slices,
                m_SeqConcat[m_lConcatenationCounter].m_LoopsOuter,
                m_SeqConcat[m_lConcatenationCounter].m_alLoopsInner[0]);
            SEQ_TRACE_ALWAYS.print("Preparing scans for this concatenation %ld\n", m_SeqConcat[m_lConcatenationCounter].m_PreparingScans);
        }

        // * ---------------------------------------------------------------------- *
        // * Execute real time events at the beginning of the concatenation loop    *
        // *                                                                        *
        // * Currently an unused dummy function                                     *                                                    *
        // * ---------------------------------------------------------------------- *
        if (!runConcatenationKernelPre(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ALWAYS.print("Error encountered in runConcatenationKernelPre(...).");
            }
            setNLSStatus(MRI_SBB_SBB_ERROR);
            return (false);
        }

        // * ---------------------------------------------------------------------- *
        // * Execute local version of preparingscans                                             *
        // * ---------------------------------------------------------------------- *
        if (!runPreparingScans(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ALWAYS.print("Error encountered in runPreparingScans(...).");
            }
            setNLSStatus(MRI_SBB_SBB_ERROR);
            return (false);
        }

        // * ---------------------------------------------------------------------- *
        // * Execute the phase correction and/or phase stabilization scans          *
        // *                                                                        *
        // * Not used by current EPI sequence. If this is required in the future    *
        // * a modified version of runPATRefScanConcatenated() will be required for *
        // * compatability with the new concatenations loop.                        *
        // * ---------------------------------------------------------------------- *
        if (!runPhaseCorScans(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ALWAYS.print("Error encountered in runPhaseCorScans(...).");
            }
            setNLSStatus(MRI_SBB_SBB_ERROR);
            return (false);
        }

        // * ---------------------------------------------------------------------- *
        // * Call of the free loop (Used e.g. for diffusion or Ciss)                *
        // * ---------------------------------------------------------------------- *
        if (!runFreeLoop(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ALWAYS.print("Error encountered in runFreeLoop(...).");
            }
            setNLSStatus(MRI_SBB_SBB_ERROR);
            return (false);
        }

        // * ---------------------------------------------------------------------- *
        // * Execute real time events at the end of the concatenation loop          *
        // *                                                                        *
        // * Currently an unused dummy function                                     *
        // * ---------------------------------------------------------------------- *
        if (!runConcatenationKernelPost(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC))
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ALWAYS.print("Error encountered in runConcatenationKernelPost(...).");
            }
            setNLSStatus(MRI_SBB_SBB_ERROR);
            return (false);
        }

    } // end of old concatenations loop

    setNLSStatus(lStatus);

    return (true);
}

//--------------------------------------------------------------------------------
// function:    runOuterSliceLoop()
//
// description: Overloaded SeqLoop function.
//              Modified to add new concatenations loop outside outer slices loop.
//--------------------------------------------------------------------------------
bool SeqLoopLongTRTrig::runOuterSliceLoop(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC)
{
    NLS_STATUS lStatus = MRI_SBB_SBB_NORMAL;

    // * -------------------------------------------------------------------------- *
    // * Run base class version if long TR triggering mode is switched off          *
    // * or if navigator triggering is active                                       *
    // * -------------------------------------------------------------------------- *
    if (!m_bLongTRTrigMode || (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF))
    {
        return SeqLoop_BASE_TYPE::runOuterSliceLoop(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC);
    }

    // * -------------------------------------------------------------------------- *
    // * New concatenations loop
    // * -------------------------------------------------------------------------- *
    m_ADCCounter = 0;

    for (m_lConcatenationCounter = m_lSliceOffsetConc = 0; m_lConcatenationCounter < m_lConcatenations; m_lConcatenationCounter++)
    {
#ifdef SHOW_LOOP_STRUCTURE
        if (IS_TR_LAND(rSeqLim))
        {
            ShowLoopStart((char*)"concatenations", m_lConcatenationCounter, m_lConcatenations);
        }
#endif

        if (IS_TR_INT(rSeqLim))
        {
            SEQ_TRACE_ALWAYS.print(
                "Concatenation %ld/%ld: %ld slices; %ld outer loops; %ld (max.) inner loops\n",
                m_lConcatenationCounter + 1,
                m_lConcatenations,
                m_SeqConcat[m_lConcatenationCounter].m_Slices,
                m_SeqConcat[m_lConcatenationCounter].m_LoopsOuter,
                m_SeqConcat[m_lConcatenationCounter].m_alLoopsInner[0]);
            SEQ_TRACE_ALWAYS.print("Preparing scans for this concatenation %ld\n", m_SeqConcat[m_lConcatenationCounter].m_PreparingScans);
        }

        m_lSliceOffset = m_lSliceOffsetConc;

        // * -------------------------------------------------------------------------- *
        // * Loop over the outer loops of the concatenation                             *
        // * (used for long TI periods)                                                 *
        // * -------------------------------------------------------------------------- *
        if ((IS_TR_INP(rSeqLim)) && (m_lOuterSliceCounter == 0) && (m_lLineCounter == 0) && (m_lPartitionCounter == 0))
        {
            SEQ_TRACE_ALWAYS.print("m_PhasesToMeasure       = %d\n", (int)m_PhasesToMeasure);
            SEQ_TRACE_ALWAYS.print("m_LoopsOuter            = %d\n", (int)m_SeqConcat[m_lConcatenationCounter].m_LoopsOuter);
        }

        for (m_lOuterSliceCounter = 0; m_lOuterSliceCounter < m_SeqConcat[m_lConcatenationCounter].m_LoopsOuter; m_lOuterSliceCounter++)
        {
#ifdef SHOW_LOOP_STRUCTURE
            if (IS_TR_LAND(rSeqLim))
            {
                ShowLoopStart((char*)"outer slices", m_lOuterSliceCounter, m_SeqConcat[m_lConcatenationCounter].m_LoopsOuter);
            }
#endif

            // get index that determines whether the current OuterLoop contains reduced number of inner slices
            // (this happens, if the number of slices within the concat is not a multiple of TotalNumberOfOuterLoops).
            // All OuterLoops < m_LoopOuterSwitch slice = m_alLoopsInner[0],
            // the others contain a reduced slicenumber = m_alLoopsInner[1].
            if (m_eSpecialSliceInterleaveMode == SpecialInterleaveMode_INTERLEAVED_IR)
            {
                // this is a special slice sorting used for interleaved IR to avoid crosstalk problems
                // (more details are given within mapLoopCounterToSliceIndex())
                // -> instead of taking OuterSliceCounter directly, it has to be mapped to an 'InterleavedOuterCounter'
                //    as it is used by mapLoopCounterToSliceCouter() during run()
                if (getInterleavedOuterCounter(m_lTotalNumberOfOuterLoops, m_lOuterSliceCounter) < m_SeqConcat[m_lConcatenationCounter].m_LoopOuterSwitch)
                {
                    m_lConcIndex = 0;
                }
                else
                {
                    m_lConcIndex = 1;
                }
            }
            else
            {
                if (m_lOuterSliceCounter < m_SeqConcat[m_lConcatenationCounter].m_LoopOuterSwitch)
                {
                    m_lConcIndex = 0;
                }
                else
                {
                    m_lConcIndex = 1;
                }
            }

            // * ---------------------------------------------------------------------- *
            // * Explaination of m_lInnerSliceNumber                                    *
            // *                                                                        *
            // * Interleaved IR   : Slices inside the actual outer loop                 *
            // * No Interleaved IR: Slices inside the actual concatenation (outer       *
            // *                    loop is only passed once)                           *
            // * ---------------------------------------------------------------------- *
            m_lInnerSliceNumber = m_SeqConcat[m_lConcatenationCounter].m_alLoopsInner[m_lConcIndex];

            // * ---------------------------------------------------------------------- *
            // * Calculate fill times                                                   *
            // * ---------------------------------------------------------------------- *
            setFillTimes(rMrProt, rSeqLim);

            // * ---------------------------------------------------------------------- *
            // * Call of the Kernel of the OuterSliceLoop                               *
            // * ---------------------------------------------------------------------- *
            if (!runOuterSliceLoopKernel(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC))
            {
                if (!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ALWAYS.print("Error encountered in runOuterSliceLoopKernel(...).");
                }
                setNLSStatus(MRI_SBB_SBB_ERROR);
                return (false);
            }

#ifdef SHOW_LOOP_STRUCTURE
            if (IS_TR_LAND(rSeqLim))
            {
                ShowLoopEnd((char*)"loop outer", m_lOuterSliceCounter);
            }
#endif

            // * ---------------------------------------------------------------------- *
            // * Increase slice offset by the number of slices that have been executed  *
            // * within the previous inner slice loop.                                  *
            // * ---------------------------------------------------------------------- *
            m_lSliceOffset += m_lInnerSliceNumber;

            //  Retrogating:
            //  When the sequence kernel reports a trigger, break out of the outer slice loop
            if (m_bIsRetroGating)
            {
                if (m_RetroGating.wasPMUTriggerFound())
                {
                    if (m_RetroGating.isMPCUTraceDeeply())
                    {
                        SEQ_TRACE_ALWAYS.print("BREAK OUT OF Outer Slice Loop\n");
                    }
                    break;
                }
            }

        } // End of loop m_lOuterSliceCounter

        // Retrogating:
        // Clear the flag that indicates that a PMU trigger was detected in preparation for the next trigger
        if (m_bIsRetroGating)
        {
            m_RetroGating.setPMUTriggerFound(false);
        }

        if (m_lTRFillEnd && (m_PhasesToMeasure == 1))
        {
            if (!insertTRFillEnd(m_lTRFillEnd))
            {
                if (!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ALWAYS.print("Error encountered in insertTRFillEnd.");
                }

                if ((m_NLSStatus & NLS_SEV) != NLS_SEV)
                {
                    setNLSStatus(MRI_SBB_SBB_ERROR);
                }

                return (false);
            }
        }

        // * ---------------------------------------------------------------------- *
        // * Add delay time before the next concatenation: TD value from protocol   *
        // * ---------------------------------------------------------------------- *
        if ((m_TDFill) && (m_lConcatenations > 1) && (m_lConcatenationCounter < m_lConcatenations - 1))
        {
            lStatus = fSBBFillTimeRun(m_TDFill);
            if ((lStatus & NLS_SEV) != NLS_SUCCESS)
            {
                if (!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ALWAYS.print("Error encountered in fSBBFillTimeRun(...).");
                }
                setNLSStatus(MRI_SBB_SBB_ERROR);
                return (false);
            }
        }

        m_lSliceOffsetConc += m_SeqConcat[m_lConcatenationCounter].m_Slices;

#ifdef SHOW_LOOP_STRUCTURE
        if (IS_TR_LAND(rSeqLim))
        {
            ShowLoopEnd((char*)"concatenation", m_lConcatenationCounter);
        }
#endif

    } // end of new concatenations loop

    setNLSStatus(lStatus);

    return (true);
}

//-----------------------------------------------------------------------------------------------------------------
// function:    runPreparingScans()
//
// description: Overloaded SeqLoop function.
//              Modified to use m_lPreparingScans instead of m_SeqConcat[m_lConcatenationCounter].m_PreparingScans.
//-----------------------------------------------------------------------------------------------------------------
bool SeqLoopLongTRTrig::runPreparingScans(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC)
{
    NLS_STATUS lStatus = MRI_SBB_SBB_NORMAL;

    // * -------------------------------------------------------------------------- *
    // * Run base class version if long TR triggering mode is switched off          *
    // * or if navigator triggering is active                                       *
    // * -------------------------------------------------------------------------- *
    if (!m_bLongTRTrigMode || (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF))
    {
        return SeqLoop_BASE_TYPE::runPreparingScans(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC);
    }

    // * -------------------------------------------------------------------------- *
    // * Helpful MPCU trace for retrogating                                         *
    // * -------------------------------------------------------------------------- *
    if (m_bIsRetroGating)
    {
        SEQ_TRACE_ALWAYS.print("# Start of preparing scan; Slice %ld/%ld #\n", m_lSliceOffsetConc + 1, (long int)rMrProt.sliceSeries().getlSize());
    }

    // * -------------------------------------------------------------------------- *
    // * Execute preparing scans if required                                        *
    // * -------------------------------------------------------------------------- *
    if ((m_ePerformPreparingScans == Always) || ((m_ePerformPreparingScans == OnlyFirstRepetition) && (m_lRepetitionCounter == 0)))
    {
        // * ---------------------------------------------------------------------- *
        // * Execute the preparing scans                                            *
        // * ---------------------------------------------------------------------- *

        if (m_lPreparingScans > 0)
        {
            // * ------------------------------------------------------------------ *
            // * Disable ADC events the preparing scans                             *
            // * ------------------------------------------------------------------ *
            fRTSetReadoutEnable(0); // disable ADC events
            m_FreeLoopCounter = 0;  // set free loop counter to a well-defined value

            for (m_lPrepareLoopCounter = 0; m_lPrepareLoopCounter < m_lPreparingScans; m_lPrepareLoopCounter++)
            {
#ifdef SHOW_LOOP_STRUCTURE
                if (IS_TR_LAND(rSeqLim))
                {
                    SEQ_TRACE_ALWAYS.print("##############################################\n");
                    SEQ_TRACE_ALWAYS.print("#    Start of preparing scan: %ld            #\n", m_lPrepareLoopCounter + 1);
                    SEQ_TRACE_ALWAYS.print("##############################################\n");
                    ShowLoopStart((char*)"preparation scans", m_lPrepareLoopCounter, m_SeqConcat[m_lConcatenationCounter].m_PreparingScans);
                }
#endif

                long ADCCounterOld = m_ADCCounter;

                m_lKernelMode = KERNEL_PREPARE;
                m_bExecuteIR  = true;
                if (!runOuterSliceLoop(rMrProt, rSeqLim, rSeqExpo, pSlcPos, psADC))
                {
                    if (!rSeqLim.isContextPrepForBinarySearch())
                    {
                        SEQ_TRACE_ALWAYS.print("Error encountered in runOuterSliceLoop(...).");
                    }
                    setNLSStatus(MRI_SBB_SBB_ERROR);
                    return (false);
                }

                m_ADCCounter = ADCCounterOld;

#ifdef SHOW_LOOP_STRUCTURE
                if (IS_TR_LAND(rSeqLim))
                {
                    ShowLoopEnd((char*)"preparation scans", m_lPrepareLoopCounter);
                    SEQ_TRACE_ALWAYS.print("##############################################\n");
                    SEQ_TRACE_ALWAYS.print("#    End of preparing scan                   #\n");
                    SEQ_TRACE_ALWAYS.print("##############################################\n");
                }
#endif

            } // End of loop m_lPrepareLoopCounter

            // * ------------------------------------------------------------------- *
            // * Re-enable ADC events at the end of the preparing scans              *
            // * ------------------------------------------------------------------- *
            fRTSetReadoutEnable(1);

        } // End of the preparing scans
    }

    // Helpful MPCU trace for retrogating
    if (m_bIsRetroGating)
    {
        SEQ_TRACE_ALWAYS.print("#   End of preparing scan; Slice %ld/%ld #\n", m_lSliceOffsetConc + 1, (long int)rMrProt.sliceSeries().getlSize());
    }

    setNLSStatus(lStatus);

    return (true);
}
