//----------------------------------------------------------------------------------
// <copyright file="SBBEPIKernelSE_ms.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2021. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#include "MrImaging/seq/a_ep2d_se/SBBEPIKernelSE.h"
#include "SBBRefocSE_ms.h"

namespace SEQ_NAMESPACE
{

class SBBEPIKernelSE_ms : public SBBEPIKernelSE
{
  public:
    SBBEPIKernelSE_ms(SBBList* pSBBList);

    virtual ~SBBEPIKernelSE_ms() = default;

    SBBRefocSE_ms& getSBBRefocSE_ms();

    virtual void setRunMode(SliceAccelRFRunMode eRunMode);

    virtual void activateFlowAttenuationAllDirWithStrength(long lFlowAttennuationStrength);

    virtual void setFlowAttenuationDirection(SEQ::GradientAxis eAxis);

    virtual void deactivateFlowAttenuation();

    /// Set flow attenuation strength in [mT/m ms^2]
    virtual void setFlowAttenuationStrength(long lFlowAttennuationStrength);

    bool initExcitation(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest() override;

  protected:

    // instantiating refocusing SBB
    void initSBBRefocSE_ms();

    // Initialize excitation: Sinc-type RF-pulse (same slice profile before and after refocusing)
    virtual bool initExcitationSinc(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    bool setMultibandParametersForExcitation(MrProt& rMrProt, SeqLim& rSeqLim);

    // Initialize excitation: External  RF-pulse (suitable for a single contrast after refocusing)
    virtual bool initExcitationOpt(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    // Check whether 'optimized' RF-pulse settings should be used
    virtual bool isExcitationOpt(const MrProt& rMrProt) const;

    bool configureRefocusingRFPulse(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    bool configureRefocusingRFPulseSinc(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    void applyGradientReversalForRefocusingRFPulse(MrProt& rMrProt, IRF_PULSE* pRFPulse);

    bool setMultibandParametersForRefocusing(MrProt& rMrProt, SeqLim& rSeqLim);

    bool configureRefocusingRFPulseOpt(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    bool prepSBBRefocus(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    void configureSpoilerGradientPolarity() override;

    // excitation pulse
    sRF_PULSE_SINC m_sExcitationRFsinc{"SincExciteRF"};
    sRF_PULSE_EXT  m_sExcitationRFopt{"OptExciteRF"};

    // refocusing pulse
    sRF_PULSE_SINC m_sRefocusingRFopt{"OptRefocRF"};

    // Factor used to scale the refocusing pulse thickness
    double m_dRefocThicknessFactor{1.2};
};

}//end of namespace SEQ_NAMESPACE
