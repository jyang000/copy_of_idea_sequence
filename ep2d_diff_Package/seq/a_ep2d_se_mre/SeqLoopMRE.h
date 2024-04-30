//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2016  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \src\MrImaging\seq\a_ep2d_se_mre\SeqLoopMRE.h
//
//     Lang: C++
//     Authors: Bolster Jr, Bradley
//              Kannengiesser, Stephan
//              Fang, Dong
//
//     Descrip: This file contains SeqLoopMRE, a derived SeqLoop variant for the MRE sequence.
//              It derives from SeqLoopEP2D.
//
//     Classes: SeqLoopMRE
//
//    -----------------------------------------------------------------------------


#ifndef LocalSeqLoopEP2D_h
#define LocalSeqLoopEP2D_h

#include "MrImaging/seq/a_ep2d_se_mre/SeqBuildBlockOptfsWithExtTrig.h"
#include "MrImaging/seq/common/SeqLoopLongTRTrig/SeqLoopLongTRTrig.h"
#include "MrImaging/seq/common/SeqLoopEPI/SeqLoopMultiband.h"
#include "MrImaging/seq/common/SeqLoopEPI/SeqLoopEP2D.h"

// For MRE: derive a SeqLoop with derived OptFS class as member
//   (for running external trigger implementation)
class SeqLoopMRE : public SeqLoopEP2D<SeqLoopMultiBand<SeqLoopLongTRTrig>>
{
public:
    SeqLoopMRE();

    ///    Preparation of the SeqLoop class and following SBB
    virtual bool prep (MrProt &rMrProt,     ///< Pointer to the protocol structure
        SeqLim &rSeqLim,     ///< Pointer to the sequence limits structure
        SeqExpo &rSeqExpo    ///< Pointer to the sequence export structure
        );

    ///    Executes the kernel calls loop.
    virtual bool runKernelCallsLoop (
        MrProt &rMrProt,     ///< Pointer to the protocol class
        SeqLim &rSeqLim,     ///< Pointer to the sequence limits class
        SeqExpo &rSeqExpo,     ///< Pointer to the sequence export class
        sSLICE_POS* pSlcPos,     ///< Pointer to the array holding the slice position information (rotation
        ///< matrices and shifts)
        sREADOUT* psADC    ///< Pointer to the sequence ADC
        );
    //  Provides local access to the Repetitions counter
    long getlRepetitionCounter() const;
    void setRepetitionsToMeasure (long);

    //  Provides local access to the Running trigger parameters
    long getlExtTriggerSpacing();
    void setlExtTriggerSpacing (long);
    long getlExtTrigStartTime();
    void setlExtTrigStartTime (long);

protected:

    /// For running trigger implementation
    SeqBuildBlockOptfsWithExtTrig m_localSBBOptfs;

    long m_lExtTriggerSpacing;
    long m_lExtTrigStartTime;
};

//------------------------------------------------------------
// inline functions
//------------------------------------------------------------
inline long SeqLoopMRE::getlRepetitionCounter() const
{
    return m_lRepetitionCounter;
}
inline void SeqLoopMRE::setRepetitionsToMeasure (long repetitions)
{
  //Don't invoke this flag as the incorrect number of repetitions result
  m_RtMSetBySequence     = true;
  m_RepetitionsToMeasure = maximum(0L, repetitions);
}
inline long SeqLoopMRE::getlExtTriggerSpacing()
{
    return m_lExtTriggerSpacing;
}
inline void SeqLoopMRE::setlExtTriggerSpacing (long T)
{
    m_lExtTriggerSpacing = T;
}
inline long SeqLoopMRE::getlExtTrigStartTime()
{
    return m_lExtTrigStartTime;
}
inline void SeqLoopMRE::setlExtTrigStartTime (long tstart)
{
    m_lExtTrigStartTime = tstart;
}

#endif
