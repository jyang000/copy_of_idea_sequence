// -----------------------------------------------------------------------------
//   Copyright (C) Siemens Healthcare GmbH 2015  All Rights Reserved.
// -----------------------------------------------------------------------------
//
//  Project: NUMARIS/4
//     File: \n4_servers1\pkg\MrServers\MrMeasSrv\SeqIF\libRT\sSYNCDATA.h
//  Version: \main\12
//   Author: CC_MEAS SCHOSTZF
//     Date: 2013-12-06 11:12:20 +01:00
//
//     Lang: C++
//
//  Model:   $subsystem
//
//  Descrip: <add a module description here>
//
// -----------------------------------------------------------------------------
//

#ifdef _MSC_VER
#pragma once
#endif

#ifndef sSYNCDATA_h
#define sSYNCDATA_h 1

#include "MrMeasSrv/SeqIF/libRT/ISYNCData.h"
#include "MrMeasSrv/SeqIF/libRT/RTEvent.h"
#include "MrMeasSrv/SeqIF/libRT/libRTDefines.h"

// --------------------------------------------------------------------------------------------
/// \ingroup libRT_oldAPI
/// \todo write documentation
/// \brief Class sSYNCDATA (synchroneous data transfer to MRIR)
class sSYNCDATA : public RTEvent_PIMPL<sSYNCDATA, ISYNCDATA>
{
  LIBRT_PIMPL_BOILERPLATE( sSYNCDATA )

public:


  MrResult prep()
  {
    return get()->prep();
  }

  MrResult setData(const SEQData& rSeqData)
  {
    return get()->setData(rSeqData) ;
  }

  const SEQData* getData() const
  {
    return get()->getData();
  }

  long getDataLength() const
  {
    return get()->getDataLength();
  }

  double getDuration() const
  {
    return get()->getDuration() ;
  }

  long getRoundedDuration(long lRoundingIncrement = 1) const
  {
    return get()->getRoundedDuration(lRoundingIncrement);
  }

  RTController& getRTController() const
  {
    return get()->getRTController();
  }

  /// compares rRhs and returns true if equal
  bool compare(const ISYNCDATA& rRhs) const
  {
    return get()->compare(rRhs);
  }

  /// copy content of rRhs to this
  ISYNCDATA& copy(const ISYNCDATA& rRhs)
  {
    return get()->copy(rRhs);
  }
};

#endif
