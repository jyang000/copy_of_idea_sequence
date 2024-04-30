// -----------------------------------------------------------------------------
//   Copyright (C) Siemens Healthcare GmbH 2015  All Rights Reserved.
// -----------------------------------------------------------------------------
//
//  Project: NUMARIS/4
//     File: \n4_servers1\pkg\MrServers\MrMeasSrv\SeqIF\libRT\include\ISYNCData.h
//  Version: \main\6
//   Author: CC_MEAS SCHOSTZF
//     Date: 2013-12-03 14:50:48 +01:00
//
//     Lang: C++
//
//  Model:   $subsystem
//
//  Descrip: <add a module description here>
//
// -----------------------------------------------------------------------------
//

#ifndef ISYNCDATA_h
#define ISYNCDATA_h 1

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include <math.h>
#include "MrGlobalDefinitions/MrResult.h"
#include "MrMeasSrv/SeqIF/libRT/IRTEvent.h"                
#include "MrMeasSrv/SeqIF/libRT/libRTDefines.h"

class SeqTime;
class SEQData;
class RTEventProcessor;

/// \ingroup libRT
/// \todo write documentation
/// \brief SyncData object
class PARC_STD_INTERFACE ISYNCDATA: public IRTEvent
{
  DECLARE_PARC_INTERFACE(ISYNCDATA)

public:

  virtual MrResult prep() = 0;
  virtual MrResult setData(const SEQData& rSeqData) = 0;
  virtual const SEQData* getData() const = 0;
  virtual long getDataLength() const = 0;

  /// Returns the total duration [us] of the readout event.
  virtual double getDuration() const = 0;

  /// \brief Returns the total duration [us] of the readout event.
  ///
  /// Always rounds up (when rounding is necessary).
  virtual long getRoundedDuration(long lRoundingIncrement = 1) const = 0;

  /// compares this with rRhs and returns true if equal
  virtual bool compare(const ISYNCDATA& rRhs) const = 0;

  /// copy content of rRhs to this
  virtual ISYNCDATA& copy(const ISYNCDATA& rRhs) = 0;

  /// assignment operator
  ISYNCDATA& operator=(const ISYNCDATA& rRhs) { return copy(rRhs); };
};


//
//------------------------------------------------------------------------------
//
inline bool operator==(const ISYNCDATA& rLhs, const ISYNCDATA& rRhs)
{
  return rLhs.compare(rRhs);
}

//
//------------------------------------------------------------------------------
//
inline bool operator!=(const ISYNCDATA& rLhs, const ISYNCDATA& rRhs)
{
  return ! (rLhs == rRhs);
}

#endif
