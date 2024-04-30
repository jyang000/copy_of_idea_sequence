//  -----------------------------------------------------------------------------
//    Copyright (C) Siemens Healthcare GmbH 2015  All Rights Reserved.
//  -----------------------------------------------------------------------------
//
//   Project: NUMARIS/4
//      File: src\MrImagingFW\libSeqUtilFW\IReorderInfo.h
//   Version: \main\4
//    Author: Heitland
//      Date: 2014-12-12 10:28:04 +01:00
//
//      Lang: C++
//
//   Descrip: MR::Measurement::Sequence::libSeqUtilFW
//
//   Classes:
//
//  -----------------------------------------------------------------------------

/*!
\file ReorderInfo.h
\brief Base class for reordering schemes
*/


#ifndef IReorderInfo_h
#define IReorderInfo_h 1

#include <algorithm>

// MrProt KSpace Wrapper
#include "MrProtSrv/Domain/MrProtData/MrProt/KSpace/MrKSpace.h"

// SeqLim
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"


//-----------------------------------------------------------------------------
// import/export control
//-----------------------------------------------------------------------------
#ifdef BUILD_libSeqUtilFW
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

//  ======================================================================
//
//  Class : IReorderInfo
//  ======================================================================

class __IMP_EXP IReorderInfo
{
public:
    IReorderInfo (void) {}
    virtual ~IReorderInfo (void) {}

    virtual long getTotalNoOfReorderIndices () = 0;
    virtual long getLinNo (long DeltaReorderIndex = 0) = 0;
    virtual long getParNo (long DeltaReorderIndex = 0) = 0;
    virtual bool isPhaseFT (long DeltaReorderIndex = 0) = 0;
    virtual bool isPartitionFT (long DeltaReorderIndex = 0) = 0;
    virtual bool isPATRefScan (long DeltaReorderIndex = 0) = 0;
    virtual bool isPATRefAndImaScan (long DeltaReorderIndex = 0) = 0;
    virtual void setIsLastAcquisition (bool bValue) = 0;
    virtual long getKSCenterLin () = 0;
    virtual bool isPATActive () = 0;
    virtual long getPATAccelerationFactorPE () = 0;
    virtual long getEchoTrainLength () = 0;

private:
    IReorderInfo(const IReorderInfo &rRright);
    IReorderInfo& operator=(const IReorderInfo &rRright);
};


#endif
