//----------------------------------------------------------------------------------
// <copyright file="SBBEPIKernelFID.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2019. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#ifndef SBBEPIKernelFID_h
#define SBBEPIKernelFID_h

#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include "MrImaging/seq/a_ep2d_diff/SequenceDebugSettings.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

namespace SEQ_NAMESPACE
{
class __IMP_EXP SBBEPIKernelFID : public SeqBuildBlockEPIKernel
{
  public:
    SBBEPIKernelFID(SBBList* pSBBList);

    virtual ~SBBEPIKernelFID() = default;

    // no copy and move constructor and assignment operator
    SBBEPIKernelFID(const SBBEPIKernelFID& right) = delete;
    SBBEPIKernelFID& operator=(const SBBEPIKernelFID& right) = delete;
    SBBEPIKernelFID(SBBEPIKernelFID&& right) = delete;
    SBBEPIKernelFID& operator=(SBBEPIKernelFID&& right) = delete;

    // init the excitation SBBs
    bool initExcitation(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    // energy export
    MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest() override;

  protected:
    // excitation pulse
    sRF_PULSE_SINC                               m_sExcitationRF;
    SequenceDebugSettings::SequenceDebugSettings m_debugSettings = SequenceDebugSettings::SequenceDebugSettings("EP2D_BOLD/USE_CUSTOM_PULSE_SETTINGS");
};

}; // namespace SEQ_NAMESPACE

#endif
