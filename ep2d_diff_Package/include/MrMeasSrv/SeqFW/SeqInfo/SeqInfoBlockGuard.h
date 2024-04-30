//-----------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH -2018 
//  All Rights Reserved. Restricted
//-----------------------------------------------------------------------------

#pragma once

#ifndef SeqInfoBlockGuard_h
#define SeqInfoBlockGuard_h


#include "MrMeasSrv/SeqIF/Sequence/ISequence.h"                       // ISequence

//-----------------------------------------------------------------------------
// Prototypes
//-----------------------------------------------------------------------------

namespace MeasSrv
{

  class ISeqInfoBlock;


  /// Class SeqInfoBlockGuard
  /// Stack guard for puhsing SeqInfoBlockData on TLS
  /// ion case the CallStack doesnt hit the sequence implementation
  /// CTOR and DTOR must be called from the same thread
  class SeqInfoBlockGuard
  {

  public:
    SeqInfoBlockGuard() = delete;
    SeqInfoBlockGuard(const SeqInfoBlockGuard&) = delete;
    SeqInfoBlockGuard& operator=(const SeqInfoBlockGuard&) = delete;

    SeqInfoBlockGuard(const ::MrMeasSrv::ISequence::Pointer&);
    SeqInfoBlockGuard(::MrMeasSrv::ISequence *);
    SeqInfoBlockGuard(::MrMeasSrv::ISequence &);

    ~SeqInfoBlockGuard();
  private:
    ::MrMeasSrv::ISequence::Pointer m_pSequence;
    ISeqInfoBlock * m_pInfoBlock = nullptr;
  };

  inline SeqInfoBlockGuard::SeqInfoBlockGuard(const ::MrMeasSrv::ISequence::Pointer& ptr)
    :m_pSequence(ptr)
  {
    if (m_pSequence)
      m_pInfoBlock = m_pSequence->createSeqInfoBlock();
  }

  inline SeqInfoBlockGuard::SeqInfoBlockGuard(::MrMeasSrv::ISequence* ptr)
    : m_pSequence(ptr)
  {
    if (m_pSequence)
      m_pInfoBlock = m_pSequence->createSeqInfoBlock();
  }

  inline SeqInfoBlockGuard::SeqInfoBlockGuard(::MrMeasSrv::ISequence& ref)
    : m_pSequence(&ref)
  {
    if (m_pSequence)
      m_pInfoBlock = m_pSequence->createSeqInfoBlock();
  }

  inline SeqInfoBlockGuard::~SeqInfoBlockGuard()
  {
    if (m_pSequence)
      m_pSequence->releaseSeqInfoBlock(m_pInfoBlock);
  }
}

#define MR_SEQUENCE_EXPORT_TLS(_sequence_) ::MeasSrv::SeqInfoBlockGuard _seqInfoBlockGuard(_sequence_)

#endif // SeqInfoBlockGuard_h
