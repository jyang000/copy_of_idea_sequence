//----------------------------------------------------------------------------------
// <copyright file="SBBRefocSE_ms.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2021. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#include "MrImaging/seq/a_ep2d_se/SBBRefocSE.h"
#include "MrImaging/libSeqUtil/SMSProperties.h"
#include "MrImaging/libSBB/SBBMultibandRF.h"
#include "MrImagingFW/libSeqUTIF/libsequt.h"

namespace SEQ_NAMESPACE
{

enum class EnumFlowAttenuationDirection
{
    FlowAttenuationDirection_None,
    FlowAttenuationDirection_Read,
    FlowAttenuationDirection_Phase,
    FlowAttenuationDirection_Slice,
    FlowAttenuationDirection_All,
};

class SBBRefocSE_ms : public SBBRefocSE
{
  public:
    SBBRefocSE_ms(SBBList* pSBBList = nullptr);

    virtual ~SBBRefocSE_ms() = default;

    bool prepSBB(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;
    bool runSBB(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC) override;

    IRF_PULSE* getRFPulsePointer() override;

    // Provide RF-pulse
    bool setRFPulseForRefocusing(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, IRF_PULSE* pRF);
    // Setup multi band RF pulse design parameters
    virtual bool setMultibandParam(long lMultibandFactor, double dMultibandDistance);
    // Setup VERSE RF pulse design parameters
    virtual bool setVERSEParam(bool bIsVerse, double dVerseFactor, double dRelativePlateauLength = 0.4);
    // Setup flag for bandwidth optimized base band of multi band RF pulse
    virtual bool setBandwidthOptimizations(bool bOptimize);

    // Set polarity of RO spoil moment
    void setIsROSpoilMomentNegative(bool bIsNegative);
    // Set polarity of PE spoil moment
    void setIsPESpoilMomentNegative(bool bIsNegative);

    /// Set flow attenuation strength [mT/m ms^2]
    void setFlowAttenuationStrength(long lFlowAttenuationStrength);

    // interface for setting flow attenuation direction
    void setFlowAttenuationDirection(SEQ::GradientAxis eAxis);

    void setFlowAttenuationDirectionNone();

    void setFlowAttenuationDirectionAll();

    /// See SeqBuildBlock for detailed information: Calculate RF info
    //  Note: Instead of adding a complete list of cuboids to this PlugIn, we just
    //        provide the current amplitude factor (used to scale the flip angles).
    //        Usually, this factor is available from the superior instance (i.e. SBBEPIKernel).
    virtual bool calcSliceAdjSBBRFInfo(
        MrProt&                                     rMrProt,         //< IMP: points to the protocol structure.
        SeqLim&                                     rSeqLim,         //< IMP: points to the sequence limits structure.
        SeqExpo&                                    rSeqExpo,        //< IMP: points to the sequence exports structure
        const SLICEADJ::sCuboidGeometry&            sSliceAdjCuboid, //< IMP: cuboid geometry (single)
        std::vector<MrProtocolData::SeqExpoRFInfo>& vsRFInfo         //< EXP: RF info
    );

    /// Get stored RF info for a certain geometry
    virtual bool getSliceAdjRFInfo(
        const SLICEADJ::sCuboidGeometry&            sSliceAdjCuboid, //< IMP: cuboid geometry
        std::vector<MrProtocolData::SeqExpoRFInfo>& vsRFInfo         //< EXP: RF info
    );

    virtual MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest();
    virtual MrProtocolData::SeqExpoRFInfo getRFInfoPerRequestMB();
    virtual long                          getDurationPerRequest();
    virtual long                          getTEContributionPerRequest();

    /// set the run mode: either single band or multi band
    virtual void setRunMode(SliceAccelRFRunMode eRunMode);


  protected:

    void setGradientMinRiseTimesAndMaxAmplitudes() override;
    bool prepareRFPulse(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;
    long getSliceSelectionGradientRampUpTime() override;
    void adaptSliceSelectionGradientTiming() override;
    void calcSliceSelectionGrad(MrProt& rMrProt) override;
    void adaptToDoubleGRT() override;
    bool prepareGradients(SeqLim& rSeqLim, double dGSSpoilMoment, double dGRSpoilMoment, double dGPSpoilMoment) override;
    bool prepareFlowAttenuationGradients();

    void setExports() override;
    void setEventStartTimes() override;
    bool checkGradients(SeqLim& rSeqLim) override;
    bool checkFlowAttenuationGradients(SeqLim& rSeqLim);
    bool rePrepareSpoilerGradient() override;
    void prepareNCOEvents(sSLICE_POS* pSLC) override;
    bool run_insertEvents(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC) override;
    bool runRTEvents(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC) override;
    long calcStartTimeOfFirstSpoilerGradient() override;
    long calcStartTimeOfSecondSpoilerGradient() override;
    bool runRefocusing(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC) override;
    bool runFlowAttenuationGradients();

    // prepares all properties of the slice accelerated refocusing SBB
    virtual bool prep_SBBMultibandRF(
        MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, SBBMultibandRF& SBBMultiband);

    //    Energy applied during one execution of the run-function [Ws].
    MrProtocolData::SeqExpoRFInfo m_RFInfoPerRequestMB;

    sGRAD_PULSE_TRAP m_GS_FlowAttenuation[2];
    sGRAD_PULSE_TRAP m_GR_FlowAttenuation[2];
    sGRAD_PULSE_TRAP m_GP_FlowAttenuation[2];

    // Invert polarity of PE spoil moment
    bool m_bPESpoilMomentNegative{false};

    /// Flow attenuation direction
    EnumFlowAttenuationDirection m_eFlowAttenuationDirection{EnumFlowAttenuationDirection::FlowAttenuationDirection_None};
    /// Flow attenuation strength [mT/m ms^2]
    long m_lFlowAttenuationStrength{0};
    /// Total time for flow attenuation gradients (as calculated during preparation)
    long m_lFlowAttenuationTotalTime{0};


    // Multiband run mode
    SliceAccelRFRunMode m_eSliceAccelRFRunMode{SINGLE_BAND};
    // SBBs for refocusing pulse which supports SMS
    SBBMultibandRF m_SBBMultibandRFRefoc{nullptr};
};

inline bool SBBRefocSE_ms::setRFPulseForRefocusing(
    MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, IRF_PULSE* pRF)
{
    resetPrepared();

    if (!m_SBBMultibandRFRefoc.setRFPulse(rMrProt, rSeqLim, rSeqExpo, pRF))
    {
        return false;
    }

#ifdef WIN32
    // If the RF-pulse thickness is scaled deliberately, SeqUT needs to know this.
    if (pRF->getThickness() != rMrProt.sliceSeries().aFront().thickness())
    {
        m_SBBMultibandRFRefoc.setRunMode(SINGLE_BAND);
        SeqUT.setRFThicknessInfo(m_SBBMultibandRFRefoc.getRFPulsePointer(), pRF->getThickness());
        m_SBBMultibandRFRefoc.setRunMode(MULTI_BAND);
        SeqUT.setRFThicknessInfo(m_SBBMultibandRFRefoc.getRFPulsePointer(), pRF->getThickness());
        m_SBBMultibandRFRefoc.setRunMode(m_eSliceAccelRFRunMode);
    }
#endif

    return true;
}

inline void SBBRefocSE_ms::setIsROSpoilMomentNegative(bool bIsNegative)
{
    m_bROSpoilMomentNegative = bIsNegative;
}

inline void SBBRefocSE_ms::setIsPESpoilMomentNegative(bool bIsNegative)
{
    m_bPESpoilMomentNegative = bIsNegative;
}

inline void SBBRefocSE_ms::setFlowAttenuationStrength( long lFlowAttenuationStrength )
{
    m_lFlowAttenuationStrength = lFlowAttenuationStrength;
}

inline void SBBRefocSE_ms::setRunMode(SliceAccelRFRunMode eRunMode)
{
    m_eSliceAccelRFRunMode = eRunMode;
}

inline long SBBRefocSE_ms::getDurationPerRequest()
{
    if (isPrepared())
    {
        return m_lSBBDurationPerRequest_us;
    }
    else
    {
        return 0;
    }
}

inline long SBBRefocSE_ms::getTEContributionPerRequest()
{
    if (isPrepared())
    {
        return m_lTEContribution;
    }
    else
    {
        return 0;
    }
}

inline bool SBBRefocSE_ms::setMultibandParam(long lMultibandFactor, double dMultibandDistance)
{
    resetPrepared();
    return m_SBBMultibandRFRefoc.setMultibandParam(lMultibandFactor, dMultibandDistance);
}

inline bool SBBRefocSE_ms::setVERSEParam(bool bIsVerse, double dVerseFactor, double dRelativePlateauLength)
{
    resetPrepared();
    return m_SBBMultibandRFRefoc.setVERSEParam(bIsVerse, dVerseFactor, dRelativePlateauLength);
}

inline bool SBBRefocSE_ms::setBandwidthOptimizations(bool bOptimize)
{
    resetPrepared();
    return m_SBBMultibandRFRefoc.setBandwidthOptimizations(bOptimize);
}

}//end of namespace SEQ_NAMESPACE
