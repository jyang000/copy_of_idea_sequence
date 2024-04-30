//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2016  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \src\MrImaging\seq\a_ep2d_se_mre\SeqBuildBlockOptfsWithExtTrig.h
//     Lang: C++
//     Authors: Bolster Jr, Bradley
//              Kannengiesser, Stephan
//              Fang, Dong
//
//     Descrip: MR::Measurement::Sequence::libSBB
//
//     Classes:
//
//    -----------------------------------------------------------------------------

#pragma once

#include "MrImaging/libSBB/SBBOptfs.h"
#include "MrMeasSrv/SeqIF/libRT/sSYNC.h"
#include "MrImaging/seq/greMRE/a_gre_mre_def.h"

// Derived class from SeqBuildBlockOptfs with continuous external trigger (for MRE)
class SeqBuildBlockOptfsWithExtTrig : public SeqBuildBlockOptfs
{

public:
    //    Constructor
    SeqBuildBlockOptfsWithExtTrig (SBBList* pSBBList);

    //    Destructor
    virtual ~SeqBuildBlockOptfsWithExtTrig();

    ///    Prepares the SBB-functions for Run. Returns TRUE, if successful.
    virtual bool prepSBB(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

    ///    run inside the Event block.
    virtual bool runInnerEvents(long lValue);

    virtual void setlTimeAfterSBB(long T);
    virtual void setlExtTriggerSpacing (long T);
    virtual long getlExtTrigStartTime();
    virtual void setlExtTrigStartTime (long tstart);

protected:

    long m_lExtTriggerSpacing;
    long m_lExtTrigStartTime;
    long m_lTimeAfterSBB;


    // External trigger
    sSYNC_EXTTRIGGER m_ExtTrigOptfs;

private:
    SeqBuildBlockOptfsWithExtTrig(const SeqBuildBlockOptfsWithExtTrig &right);
    const SeqBuildBlockOptfsWithExtTrig & operator=(const SeqBuildBlockOptfsWithExtTrig &right);
};
inline void SeqBuildBlockOptfsWithExtTrig::setlTimeAfterSBB(long T)
{
    m_lTimeAfterSBB = T;
}
inline void SeqBuildBlockOptfsWithExtTrig::setlExtTriggerSpacing (long T)
{
    m_lExtTriggerSpacing = T;
}
inline long SeqBuildBlockOptfsWithExtTrig::getlExtTrigStartTime()
{
    return m_lExtTrigStartTime;
}
inline void SeqBuildBlockOptfsWithExtTrig::setlExtTrigStartTime (long tstart)
{
    m_lExtTrigStartTime = tstart;
}
