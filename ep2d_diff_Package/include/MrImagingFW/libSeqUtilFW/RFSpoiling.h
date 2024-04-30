//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2015  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \src\MrImagingFW\libSeqUtilFW\RFSpoiling.h
//     Version:
//      Author: KUEHN
//        Date: n.a.
//
//        Lang: C++
//
//     Descrip: MR::Measurement::Sequence::libSeqUtilFW
//
//     Classes:
//
//    -----------------------------------------------------------------------------

#ifndef RFSpoiling_h
#define RFSpoiling_h 1


// MrProtocolData::MrProtData
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
// SeqLim
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"

//-----------------------------------------------------------------------------
// import/export control
//-----------------------------------------------------------------------------
#ifdef BUILD_libSeqUtilFW
  #define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

//    Wrapper class to comfortably use the old C-function fSUVerifyRFSpoiling
//    without declaring lots of static variables in the sequence.
//
//    Do it like this:
//    - create static object within your sequence:
//      static RFSpoiling myRFSpoil;
//    - in fSEQRunKernel do the following:
//      if (!myRFSpoil.verify(pMrProt,pSeqLim,lSlice))
//      {
//          // error handling !!!
//          // ...
//      }
//      if (!myRFSpoil.setPhase())
//      {
//          // error handling !!!
//          // ...
//      }
//    - then you can access the actual spoil phase using the method
//      myRFSpoil.getPhase()


class __IMP_EXP RFSpoiling
{
  public:
      RFSpoiling();

      virtual ~RFSpoiling();

      bool verify (MrProt &rMrProt, SeqLim &rSeqLim, long lSlice);

      virtual bool setPhase ();

      double getPhase ();

      NLS_STATUS getNLSStatus ();


  protected:

      virtual bool doSUVerifyRFSpoil (MrProt &rMrProt, SeqLim &rSeqLim, long lSlice);

    // Data Members for Class Attributes

      bool m_bSpoilingActive;

      double m_dRFSpoilPhase;

      double m_dRFSpoilPhasePrevSlice;

      double m_dRFSpoilIncrement;

      double m_dRFSpoilIncrementPrevSlice;

      double m_dRFPrevSlicePosSag;

      double m_dRFPrevSlicePosCor;

      double m_dRFPrevSlicePosTra;

      double m_dRFPrevSliceNormalSag;

      double m_dRFPrevSliceNormalCor;

      double m_dRFPrevSliceNormalTra;

      double m_dRFSPOIL_INCREMENTdeg;

      NLS_STATUS m_NLSStatus;


  private:
      RFSpoiling(const RFSpoiling &right);

      RFSpoiling & operator=(const RFSpoiling &right);
};


#endif
