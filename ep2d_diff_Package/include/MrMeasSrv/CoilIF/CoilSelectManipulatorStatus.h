// -----------------------------------------------------------------------------
//   Copyright (C) Siemens Healthcare GmbH 2015  All Rights Reserved.
// -----------------------------------------------------------------------------
//     File: \n4_servers1\pkg\MrMeasSrv\CoilIF\CoilSelectManipulatorStatus.h
//  Descrip: MR::MrServers::MrMeasSrv::SpuSer::CoilInterface
// -----------------------------------------------------------------------------

#pragma once


#ifndef CoilSelectManipulatorStatus_h
#define CoilSelectManipulatorStatus_h 1

#include "MrVista/Parc/Reflection/IObject.h"
#include "MrMeasSrv/MeasUtils/MeasRealTimeStatic.h"
#include <string>

class MeasDumpEnv;
class CoilSelectManipulatorStatus;

namespace MeasSrv
{

/// Class CoilSelectManipulator
class ICoilSelectManipulatorStatusStatic : public Parc::Interface
{
public:

  DECLARE_PARC_INTERFACE( ICoilSelectManipulatorStatusStatic );

protected:

  static ICoilSelectManipulatorStatusStatic* create();

  virtual bool isGood(const CoilSelectManipulatorStatus*) const = 0;

  virtual void dump(const CoilSelectManipulatorStatus*) const = 0;

  virtual void dump(const CoilSelectManipulatorStatus*, MeasDumpEnv&) const = 0;


  //Returns true, if the Status
  //ICoilSelectManipulator::COILSel_OK
  //Ignore following Bits:
  //- ICoilSelectManipulator::COMP_ELEM_EXCHANGED
  //- ICoilSelectManipulator::UNNEEDED_ELEM_COMBINATION
  //- ICoilSelectManipulator::RX_CHAN_INVAL
  virtual bool isModifySilence(const CoilSelectManipulatorStatus* pObj) = 0;

  //Returns true, if the Status
  //ICoilSelectManipulator::COILSel_OK
  //Ignore following Bits:
  //- ICoilSelectManipulator::CONNECTED_ELEM_CHANGED
  //- ICoilSelectManipulator::UNNEEDED_ELEM_COMBINATION
  virtual bool isValidForMeasurement(const CoilSelectManipulatorStatus* pObj) = 0;


  /// Get the current status as string
  virtual std::string getStatusStr(const CoilSelectManipulatorStatus* pObj) = 0;

  friend class ::CoilSelectManipulatorStatus;
};


MEAS_REALTIME_STATIC_VAR_DECLARE(ICoilSelectManipulatorStatusStatic*,ICoilSelectManipulatorStatusStatic, nullptr);


//
//-----------------------------------------------------------------------------
//
inline ICoilSelectManipulatorStatusStatic* ICoilSelectManipulatorStatusStatic::create()
{
  // cast volatile away
  ICoilSelectManipulatorStatusStatic* pCached = const_cast<ICoilSelectManipulatorStatusStatic*> (MEAS_REALTIME_STATIC_LOAD_ACQUIRE(ICoilSelectManipulatorStatusStatic));

  if (pCached)
    return pCached;

  // object is a parc singleton. no need to keep the ref count
  Pointer ptr = Pointer::Create("CoilSelectManipulatorStatusStaticImpl@CoilManipulator");
  MEAS_REALTIME_STATIC_STORE_RELEASE(ICoilSelectManipulatorStatusStatic, ptr.get());

  return ptr.get();
}

}

class CoilSelectManipulatorStatus
{

public:

  /// Status bit constants
  enum
  {
    COILSEL_OK                    = 0,
    COILSEL_CHANGED               = 0x00000001, // Info
    COILSEL_EMPTY                 = 0x00000002, // Error
    ELEM_DESELECTED               = 0x00000004, // Info
    ELEM_NOT_CONNECTED            = 0x00000008, // Error
    NO_RX_CHAN_AVAIL              = 0x00000010, // Error
    RX_CHAN_INVAL                 = 0x00000020, // Error
    MUX_CHAN_INVAL                = 0x00000040, // Error
    NO_TX_ELEM                    = 0x00000080, // Error
    CONNECTED_ELEM_CHANGED        = 0x00000100, // Info
    COMP_ELEM_EXCHANGED           = 0x00000200, // Info
    UNNEEDED_ELEM_COMBINATION     = 0x00000400, // Error
    ELEM_IS_NOT_SELECTED          = 0x00000800, // Error
    NUCLEUS_NOT_SUPPORTED_BY_ELEM = 0x00001000, // Error
    NUCLEUS_DOES_NOT_MATCH        = 0x00002000, // Error
    PAT_ELEM_DESELECTED           = 0x00004000, // Info
    CONNECTED_ELEM_INVAL          = 0x00008000, // Error
    ELEM_COULD_NOT_BE_SELECTED    = 0x00010000, // Error
    ELEM_COULD_NOT_BE_DESELECTED  = 0x00020000, // Error
    BUILD_REQUIRED                = 0x00040000, // Error
    CSM_INTERNAL_ERROR            = 0x00080000, // Error
    COILENV_INVALID               = 0x00100000, // Error
    NOT_ENOUGH_RX_CHAN_AVAIL      = 0x00200000, // Info
    RX_CHAN_ALREADY_ASSIGNED      = 0x00400000, // Info
    MUX_CONSTRAINT_CONFLICT       = 0x00800000, // info
    TX_CHAN_ALREADY_ASSIGNED      = 0x01000000, // Info
    BCCOMBINEMATRIX_COULD_NOT_SET = 0x02000000, // Error
    TX_CHAN_INVAL                 = 0x04000000, // Error
    RX_CHAN_CHANGED               = 0x10000000, // Info
    TX_CHAN_CHANGED               = 0x20000000, // Info
    LC_ELEM_DESELECTED            = 0x40000000  // Info
  };


  CoilSelectManipulatorStatus(unsigned long ulStatus);

  unsigned long getStatus() const;

  //Returns true, if the Status
  //CoilSelectManipulator::COILSel_OK
  //Ignore following Bits:
  //- CoilSelectManipulator::CONNECTED_ELEM_CHANGED
  //- CoilSelectManipulator::UNNEEDED_ELEM_COMBINATION
  bool isValidForMeasurement() const
  {
    bool bIsValidForMeasurement = m_pImpl->isValidForMeasurement(this);
    return bIsValidForMeasurement;
  }

  //Return true, if the status flag PAT_ELEM_DESELECTED is
  //set.
  bool isPATElementDeselected() const;

  //Returns true, if the Status
  //CoilSelectManipulator::COILSel_OK
  //Ignore following Bits:
  //- CoilSelectManipulator::COMP_ELEM_EXCHANGED
  //- CoilSelectManipulator::UNNEEDED_ELEM_COMBINATION
  //- CoilSelectManipulator::RX_CHAN_INVAL
  bool isModifySilence() const
  {
    bool bIsModifySilence = m_pImpl->isModifySilence(this);
    return bIsModifySilence;
  }

  bool isNucleusNotSupportedByElement() const;

  bool isNucleusDoesNotMatch() const;

  void dump()
  {
    m_pImpl->dump(this);
  }

  void dump(MeasDumpEnv& rMeasDumpEnv)
  {
    m_pImpl->dump(this, rMeasDumpEnv);
  }

  const CoilSelectManipulatorStatus &  operator |=(const CoilSelectManipulatorStatus& rCoilSelectManipulatorStatus);

  const CoilSelectManipulatorStatus& operator |=(unsigned long ulStatus);

  const CoilSelectManipulatorStatus& operator =(unsigned long ulStatus);

  bool isElementIsNotSelected() const;

  bool isElementCouldNotBeSelected() const;

  bool isElementCouldNotBeDeselected() const;

  CoilSelectManipulatorStatus();

  operator bool()
  {
    bool bGood = m_pImpl->isGood(this);
    return bGood;
  }

  bool isCoilSelectChanged() const;

  bool isBuildRequired() const;

  bool isCoilSelectEmpty() const;

  bool isElementDeselected() const;

  bool isLCElementDeselected() const;

  bool isElementNotConnected() const;

  bool isNoRxChannelAvailable() const;

  bool isRxChannelInvalid() const;

  bool isMuxChannelInvalid() const;

  bool isNoTxElementAvailable() const;

  bool isConnectedElementChanged() const;

  bool isCompatibleElementExchanged() const;

  bool isUnneededElementCombination() const;

  bool isConnectedElementInvalid() const;

  bool isInternalError() const;

  bool isCoilEnvironmentInvalid() const;

  bool isNotEnoughRxChannelAvailable() const;

  bool isRxChannelAlreadyAssigned() const;

  bool isMuxConstraintConflictAvailable() const;

  bool isTxChannelAlreadyAssigned() const;

  bool isBCCombineMatrixCouldNotSet() const;

  bool isTxChannelInvalid() const;

  bool isRxChannelChanged() const;

  bool isTxChannelChanged() const;

  /// Get the current status as string
  std::string getStatusStr()
  {
    return m_pImpl->getStatusStr(this);
  }

protected:



private:
  MeasSrv::ICoilSelectManipulatorStatusStatic* m_pImpl;

  unsigned long m_ulStatus;

  //---------------------------------------------------------------------------
  // Friends
  //---------------------------------------------------------------------------
  friend class CoilSelectManipulator;

};

//
//-----------------------------------------------------------------------------
//
inline CoilSelectManipulatorStatus::CoilSelectManipulatorStatus() 
  :m_pImpl(MeasSrv::ICoilSelectManipulatorStatusStatic::create())
  ,m_ulStatus(COILSEL_OK)
{}

//
//-----------------------------------------------------------------------------
//
inline CoilSelectManipulatorStatus::CoilSelectManipulatorStatus(unsigned long ulStatus)
  :m_pImpl(MeasSrv::ICoilSelectManipulatorStatusStatic::create())
  ,m_ulStatus(ulStatus)
{}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isCoilSelectChanged() const
{
  return ((m_ulStatus & COILSEL_CHANGED) == COILSEL_CHANGED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isBuildRequired() const
{
  return ((m_ulStatus & BUILD_REQUIRED) == BUILD_REQUIRED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isCoilSelectEmpty() const
{
  return ((m_ulStatus & COILSEL_EMPTY) == COILSEL_EMPTY);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isElementDeselected() const
{
  return ((m_ulStatus & ELEM_DESELECTED) == ELEM_DESELECTED);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isLCElementDeselected() const
{
  return ((m_ulStatus & LC_ELEM_DESELECTED) == LC_ELEM_DESELECTED);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isPATElementDeselected() const
{
  return ((m_ulStatus & PAT_ELEM_DESELECTED) == PAT_ELEM_DESELECTED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isElementNotConnected() const
{
  return ((m_ulStatus & ELEM_NOT_CONNECTED) == ELEM_NOT_CONNECTED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isNoRxChannelAvailable() const
{
  return ((m_ulStatus & NO_RX_CHAN_AVAIL) == NO_RX_CHAN_AVAIL);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isRxChannelInvalid() const
{
  return ((m_ulStatus & RX_CHAN_INVAL) == RX_CHAN_INVAL);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isMuxChannelInvalid() const
{
  return ((m_ulStatus & MUX_CHAN_INVAL) == MUX_CHAN_INVAL);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isNoTxElementAvailable() const
{
  return ((m_ulStatus & NO_TX_ELEM) == NO_TX_ELEM);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isConnectedElementChanged() const
{
  return ((m_ulStatus & CONNECTED_ELEM_CHANGED) == CONNECTED_ELEM_CHANGED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isConnectedElementInvalid() const
{
  return ((m_ulStatus & CONNECTED_ELEM_INVAL) == CONNECTED_ELEM_INVAL);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isCompatibleElementExchanged() const
{
  return ((m_ulStatus & COMP_ELEM_EXCHANGED) == COMP_ELEM_EXCHANGED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isUnneededElementCombination() const
{
  return ((m_ulStatus & UNNEEDED_ELEM_COMBINATION) == UNNEEDED_ELEM_COMBINATION);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isElementIsNotSelected() const
{
  return ((m_ulStatus & ELEM_IS_NOT_SELECTED) == ELEM_IS_NOT_SELECTED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isElementCouldNotBeSelected() const
{
  return ((m_ulStatus & ELEM_COULD_NOT_BE_SELECTED) == ELEM_COULD_NOT_BE_SELECTED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isElementCouldNotBeDeselected() const
{
  return ((m_ulStatus & ELEM_COULD_NOT_BE_DESELECTED) == ELEM_COULD_NOT_BE_DESELECTED);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isNucleusNotSupportedByElement() const
{
  return ((m_ulStatus & NUCLEUS_NOT_SUPPORTED_BY_ELEM) == NUCLEUS_NOT_SUPPORTED_BY_ELEM);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isNucleusDoesNotMatch() const
{
  return ((m_ulStatus & NUCLEUS_DOES_NOT_MATCH) == NUCLEUS_DOES_NOT_MATCH);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isInternalError() const
{
  return ((m_ulStatus & CSM_INTERNAL_ERROR) == CSM_INTERNAL_ERROR);
}


//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isCoilEnvironmentInvalid() const
{
  return ((m_ulStatus & COILENV_INVALID) == COILENV_INVALID);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isNotEnoughRxChannelAvailable() const
{
  return ((m_ulStatus & NOT_ENOUGH_RX_CHAN_AVAIL) == NOT_ENOUGH_RX_CHAN_AVAIL);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isRxChannelAlreadyAssigned() const
{
  return ((m_ulStatus & RX_CHAN_ALREADY_ASSIGNED) == RX_CHAN_ALREADY_ASSIGNED);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isMuxConstraintConflictAvailable() const
{
  return ((m_ulStatus & MUX_CONSTRAINT_CONFLICT) == MUX_CONSTRAINT_CONFLICT);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isTxChannelAlreadyAssigned() const
{
  return ((m_ulStatus & TX_CHAN_ALREADY_ASSIGNED) == TX_CHAN_ALREADY_ASSIGNED);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isBCCombineMatrixCouldNotSet() const
{
  return ((m_ulStatus & BCCOMBINEMATRIX_COULD_NOT_SET) == BCCOMBINEMATRIX_COULD_NOT_SET);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isTxChannelInvalid() const
{
  return ((m_ulStatus & TX_CHAN_INVAL) == TX_CHAN_INVAL);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isRxChannelChanged() const
{
  return ((m_ulStatus & RX_CHAN_CHANGED) == RX_CHAN_CHANGED);
}

//
//-----------------------------------------------------------------------------
//
inline bool CoilSelectManipulatorStatus::isTxChannelChanged() const
{
  return ((m_ulStatus & TX_CHAN_CHANGED) == TX_CHAN_CHANGED);
}


//
//-----------------------------------------------------------------------------
//
inline const CoilSelectManipulatorStatus &  CoilSelectManipulatorStatus::operator |=(const CoilSelectManipulatorStatus& rCoilSelectManipulatorStatus)
{
  m_ulStatus |= rCoilSelectManipulatorStatus.m_ulStatus;

  return (*this);
}


//
//-----------------------------------------------------------------------------
//
inline const CoilSelectManipulatorStatus&CoilSelectManipulatorStatus:: operator|= (unsigned long ulStatus)
{
  m_ulStatus |= ulStatus;

  return (*this);
}


//
//-----------------------------------------------------------------------------
//
inline const CoilSelectManipulatorStatus&CoilSelectManipulatorStatus:: operator= (unsigned long ulStatus)
{
  m_ulStatus = ulStatus;

  return (*this);
}

//
//-----------------------------------------------------------------------------
//
inline unsigned long CoilSelectManipulatorStatus::getStatus() const
{
  return (m_ulStatus);
}



#endif
