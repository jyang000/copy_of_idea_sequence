//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/X
//        File: \src\MrImaging\seq\a_tgse_asl\ReorderInfoGrase3D.h
//      Author: PLM AW NERUO
//        Date: 2018-08-07 12:23:43 +02:00
//
//        Lang: C++
//
//
//
///  \file   ReorderInfoGrase3D.h
///  \brief  File containing declaraion of the Reorder class
///        
///
///  This file contains the implementation of the class ReorderInfoGrase3D.
///
//    -----------------------------------------------------------------------------

#pragma once

// ReorderInfoEPI
#include "MrImaging/libSeqUtil/ReorderInfoEPI.h"
#include "MrImaging/seq/SeqDebug.h"
#define DEBUG_ORIGIN 0x00800000

class ReorderInfoGrase3D : public ReorderInfoEPI
{
public:
    ReorderInfoGrase3D (long lMaxNoOfReorderIndices = -1, long lCalcBufferSize = -1);

    virtual ~ReorderInfoGrase3D(){};

    virtual void setGsSegments(long seg);
    virtual long getGsSegments();
    virtual long getKCenterParIndex();
    bool isPATRefScan (long DeltaReorderIndex = 0) override;
    bool isLastScanInSlice (long DeltaReorderIndex=0) override;

    long getLinesToMeasure()      {return m_lLinesToMeasure;     };
    long getPartitionsToMeasure() {return m_lPartitionsToMeasure;};
    long getMeasuredSegments_Lines() {return m_lMeasuredSegments_Line;};
    long getMeasuredSegments_Partitions() {return m_lMeasuredSegments_Partition;};

  private:

    bool prepareCalculation_Standard(MrProt& rMrProt, SeqLim& rSeqLim) override;
    bool reorderEPI_Standard(MrProt& rMrProt, SeqLim& rSeqLim) override;

#ifdef SUPPORT_iPAT_TGSE
    bool   prepareCalculation_PATRefScanEPI(MrProt& rMrProt, SeqLim& rSeqLim) override;
    bool   reorderEPI_PATRefScanEPI(MrProt& rMrProt, SeqLim& rSeqLim) override;
    void   setPATReorderIndexOffsetForRefScans(long lPATRefScanCounterInSegment) override;
    double getPartialFourierFactor(SEQ::PartialFourierFactor pffPFF);

    long   m_lMinPATLineNo{0};
    long   m_lMinPATParNo{0};
    long   m_lLastScanInSliceIndex{0};
    bool   m_bIsPAT{false};
#endif

    int32_t m_lLinesToMeasure{-1};
    int32_t m_lPartitionsToMeasure{-1};
    long    m_lGsSegments{-1};
    long    m_lPartitionsPerSegment{-1};
    long    m_lMeasuredSegments_Line{-1};
    long    m_lMeasuredSegments_Partition{-1};
};


#ifdef SUPPORT_iPAT_TGSE
inline double ReorderInfoGrase3D::getPartialFourierFactor(SEQ::PartialFourierFactor pffPFF) 
{
    switch (pffPFF) 
    {
    case SEQ::PF_HALF: 
      return 0.5;
      break;
    case SEQ::PF_5_8: 
      return 5./8.;
      break;
    case SEQ::PF_6_8: 
      return 0.75;
      break;
    case SEQ::PF_7_8: 
      return 7./8.;
      break;
    case SEQ::PF_OFF: 
    default:
      return 1.;
      break;
    }
};
#endif //SUPPORT_iPAT_TGSE
