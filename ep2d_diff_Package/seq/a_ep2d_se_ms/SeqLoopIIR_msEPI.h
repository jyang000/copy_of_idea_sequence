//----------------------------------------------------------------------------------
// <copyright file="SeqLoopIIR_msEPI.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2021. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

//  Definition of base class
#include "MrImaging/libSBB/SeqLoopIIR.h"

namespace SEQ_NAMESPACE
{
class SeqLoopIIR_msEPI : public SeqLoopIIR
{
  public:
    SeqLoopIIR_msEPI();

    virtual ~SeqLoopIIR_msEPI() = default;

    // Control continuous pulsing for MT saturation
    virtual void ActivateContinousMTSatPulsing(bool bContinousMTSatPulsing);
    virtual bool IsContinousMTSatPulsing() const;

    // Control continuous pulsing for regional saturation
    virtual void ActivateContinousRSatPulsing(bool bContinousRSatPulsing);
    virtual bool IsContinousRSatPulsing() const;

    // Get energy of all active continuous saturations
    MrProtocolData::SeqExpoRFInfo getRFInfoContinuousSat(MrProt& rMrProt, SeqLim& rSeqLim) override;

    //  Check whether the scheme PROT_MASK_IIR_SCHEME_STD will get used
    //  Note: This information is available even before prepartion!
    virtual bool IsIIRSchemeStd(MrProt& rMrProt, SeqLim& rSeqLim);

    virtual void setLongTRTrigMode(bool bSwitch);
    virtual bool isLongTRTrigMode() const;

  protected:

    void adaptNumberOfExecutionsInCheck() override;

    bool executeSeparatePATRefScans(
        MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, long lNExe_run) override;

    long getInnerSliceIndex(MrProt& rMrProt, long lCChronPos, long lEChronPos = 0, long lCSweep = 0) const override;

    bool isTimeSufficientForSatPulses_blockwise(long lContinuousSatDuration) override;

    bool isTimeSufficientForSatPulses_iiee(long lTFill_us, long lContinuousSatDuration) override;

    long calcTotalContinuousSatDuration() override;

    bool runAllSatPulses(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSlicePos) override;

    long setInitialPositionForOneConcat(MrProt& rMrProt, const SEQ::SeriesMode iExcitOrder) override;

    bool setChronologicalPositionForOneConcat(
        MrProt&                      rMrProt,
        long&                        lAPos,
        const SEQ::SeriesMode        iExcitOrder,
        std::vector<long>::iterator& pChronPos,
        std::vector<CONC>::iterator  pConc,
        const SliceSeries&           rSeries) override;

    bool setFinalPositionForOneConcat(MrProt& rMrProt, const SEQ::SeriesMode iExcitOrder, long& lAPos) override;

    void setPositionForMultipleConcats(
        MrProt& rMrProt, long& lCIR, long& lAPos, std::vector<long>::iterator& pChronPos, const SliceSeries& rSeries) override;

    double calcNumberOfBlocksForSat(MrProt& rMrProt) override;

    double calcNumberOfSatModules(MrProt& rMrProt, double dNumBlocks) override;

    MrProtocolData::SeqExpoRFInfo calcTotalEnergyFromAllSatPulses(double dNumberOfSatModules) override;

    MrProtocolData::SeqExpoRFInfo calcTotalEnergyFromMSat(double dNumberOfSatModules);

    MrProtocolData::SeqExpoRFInfo calcTotalEnergyFromRSat(double dNumberOfSatModules);

    void printEnergyRelevantInfos(
        MrProt& rMrProt, const double dNumBlocks, MrProtocolData::SeqExpoRFInfo& rfInfo_AddEnergy) override;

    // Control continuous pulsing for magnetization transfer saturation
    bool m_bContinousMTSatPulsing{false};

    // Control continuous pulsing for regional saturation
    bool m_bContinousRSatPulsing{false};

};

} // namespace SEQ_NAMESPACE