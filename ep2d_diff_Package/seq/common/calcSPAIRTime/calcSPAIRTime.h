//----------------------------------------------------------------------------------
// <copyright file="calcSPAIRTime.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 2006-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#ifndef calcSPAIRTime_h
#define calcSPAIRTime_h

#include "MrImaging/libSBB/SBBOptfs.h"
#include "MrImaging/libSBB/SBBOptfsPrep.h"

/// This class calculates the SPAIR inversion time for the EPI sequence
class calcSPAIRTime
{
  public:
    calcSPAIRTime();
    
    virtual ~calcSPAIRTime() = default;

    calcSPAIRTime(const calcSPAIRTime& right) = delete;
    const calcSPAIRTime& operator=(const calcSPAIRTime& right) = delete;
    calcSPAIRTime(calcSPAIRTime&& right)                       = delete;
    const calcSPAIRTime& operator=(calcSPAIRTime&& right) = delete;

    static bool calcSPAIRTimeEPI(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lScanTimeBasic, long& lScanTimeSatsEtc, SeqBuildBlockOptfsPrep* pSBBOptfsPrep, SeqBuildBlockOptfs* pSBBOptfs);

  protected:
 
  private:

};


#endif
