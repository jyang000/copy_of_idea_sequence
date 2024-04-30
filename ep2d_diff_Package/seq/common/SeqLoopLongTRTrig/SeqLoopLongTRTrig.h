//----------------------------------------------------------------------------------
// <copyright file="SeqLoopLongTRTrig.h" company="Siemens Healthcare GmbH">
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

#pragma once

#ifndef SeqLoopLongTRTrig_h
#define SeqLoopLongTRTrig_h

//------------------------------------------------------------
// Specify SeqLoop version to be used as base class
//------------------------------------------------------------

// standard SeqLoop
#include "MrImaging/libSBB/SEQLoop.h"
typedef SeqLoop SeqLoop_BASE_TYPE;

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif
namespace SEQ_NAMESPACE
{
class SeqLoopLongTRTrig : public SeqLoop_BASE_TYPE
{
  public:
    // constructor
    SeqLoopLongTRTrig() = default;

    // destructor
    virtual ~SeqLoopLongTRTrig() = default;

    // switch diffusion triggering mode on/off
    void setLongTRTrigMode(bool bSwitch);

    // get status of diffusion triggering mode
    bool isLongTRTrigMode() const;

  protected:
    // overloaded SeqLoop function: runConcatenationLoop
    bool runConcatenationLoop(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC) override;

    // overloaded SeqLoop function: runOuterSliceLoop
    bool runOuterSliceLoop(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC) override;

    // overloaded SeqLoop function: runPreparingScans
    bool runPreparingScans(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC) override;

    // flag to specify whether long TR triggering mode is active
    // default value set to true by constructor
    bool m_bLongTRTrigMode{true};
};

} // namespace SEQ_NAMESPACE

using namespace SEQ_NAMESPACE;

inline void SeqLoopLongTRTrig::setLongTRTrigMode(bool bSwitch)
{
    m_bLongTRTrigMode = bSwitch;
}

inline bool SeqLoopLongTRTrig::isLongTRTrigMode() const
{
    return m_bLongTRTrigMode;
}

#endif
