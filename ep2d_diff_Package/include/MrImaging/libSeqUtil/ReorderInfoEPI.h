//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2021  All Rights Reserved.
//    -----------------------------------------------------------------------------

#ifndef ReorderInfoEPI_h
#define ReorderInfoEPI_h 1

#include "MrImagingFW/libSeqUtilFW/ReorderInfo.h"

//-----------------------------------------------------------------------------
// import/export control
//-----------------------------------------------------------------------------
#ifdef BUILD_libSeqUtil
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

enum class EPIAcquisitionMode
{
    SingleShot, MultiShot
};

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4266)
#endif


class __IMP_EXP ReorderInfoEPI : public ReorderInfo
{

    friend struct ReorderInfoEPI_UT;

  public:

    ReorderInfoEPI(long lMaxNoOfReorderIndices = -1, long lCalcBufferSize = -1);


    virtual ~ReorderInfoEPI() = default;

    ReorderInfoEPI(const ReorderInfoEPI& right) = delete;
    ReorderInfoEPI& operator=(const ReorderInfoEPI& right) = delete;
    ReorderInfoEPI(ReorderInfoEPI&& right)                 = delete;
    ReorderInfoEPI& operator=(ReorderInfoEPI&& right) = delete;

    // Rationale why "MultiShot" is used rather than "Segmented":
    // The term "segments" is already introduced in the base class. It does NOT refer to
    // the distribution of k-space lines among multiple shots, but rather denotes regions
    // of k-space with similar magnetization evolution. For multi-echo acquisitions
    // (e.g. TSE), each echo is assigned to a segment: thus, the first "segment" comprises
    // data from all first echoes, and the last "segment" comprises data from all last
    // echoes.
    // For the sake of a consistent nomenclature (and since internally both, "shot index"
    // and "segment index" will be required), the term "shot" will be used for denoting the
    // index of a multi-shot acquisition.

    void setModeMultiShot();

    void setModeSingleShot();

    bool prepareCalculation(MrProt& rMrProt, SeqLim& rSeqLim) override;

    bool reorderEPI(MrProt& rMrProt, SeqLim& rSeqLim);

    bool isMultiShot() const;

    bool isSingleShot() const;

    // returns number of contrasts to measure
    long getContrastsToMeasure() const;

    long getKSpaceCenterSegment() const;

    long getCounterInSegmentWithKSpaceCenter() const;

    bool isPATRefAndImaScan(long = 0) override;

    bool isPATRefScan(long DeltaReorderIndex = 0) override;

    void setPATMultiShotRefScans(bool bValue);

    long getPATEchoTrainLengthRefScans() const;

    long getPATRefKSpaceCenterSegment() const;

    long getPATRefCounterInSegmentWithKSCenter() const;

    virtual long getPATBlindADCsBeforeRefScans() const;

    virtual long getPATBlindADCsAfterRefScans() const;

    virtual void setPATReorderIndexOffsetForRefScans(long lPATRefScanCounterInSegment);

    virtual void setPATReorderIndexOffsetForImagingScans();

    // set reordering offset to first echo of given partition / shot / contrast </summary>
    // PrepareCalculation() is required before calling this method.
    void setReorderIndexOffset(long lEchoNo, long lContrastNo, long lShotNo, long lPartitionNo);

    // compatibility for one-argument base implementation (non-virtual)
    void setReorderIndexOffset(long lReorderInfoOffset);

    // increase reordering offset by given partition / shot / contrast </summary>
    // ::PrepareCalculation() is required before calling this method.
    void increaseReorderIndexOffset(long lEchoInc, long lContrastInc, long lShotInc, long lPartitionInc);

    // compatibility for one-argument base implementation (non-virtual)
    void increaseReorderIndexOffset(long ReorderIndexOffsetIncrement);

    virtual void setPATGRERefScans(bool bValue);

    virtual bool isPATGRERefScans() const;

    // Number of lines per segment, divided by the PAT factor
    virtual long getPATLinesPerSegment() const;

    bool isLastScanInSlice(long lDeltaReorderIndex = 0) override;

    bool isFirstScanInSlice(long lDeltaReorderIndex = 0) override;

    // Reverse phase encoding direction for gradient-reversal distortion correction
    virtual void setReversePE(bool bValue);

    // retrieves the current state of reverse phase encoding
    virtual bool getIsReversePE() const;

    // Note: Slice-acceleration related methods should get moved to the base class.
    //       This would enable access to the corresponding information in methods
    //       which support access to the base class only (e.g. within SBBEPIReadout)

    /// returns true if the reordering is prepared for slice-acceleration
    virtual bool isSliceAccelerationActive() const;

    /// returns the slice-acceleration factor of the prepared reordering
    virtual long getSliceAccelerationFactor() const;

    /// returns the FOV-shift factor of the prepared reordering
    virtual long getFOVShiftFactor() const;

    /// set multi-band mode (similar to enabling the phase-correction mode in the base class)
    virtual void setRunModeMultiBand();

    /// set single-band mode (similar to disabling the phase-correction mode in the base class)
    virtual void setRunModeSingleBand();

    long getParNo(long lDeltaReorderIndex) override;

    long getParNoCenterZero(long lDeltaReorderIndex) override;

    long getKSCenterPar() override;

    bool isKSCenterPar(long lDeltaReorderIndex = 0) override;

    long getMinParNoCenterZero() override;

    long getMaxParNoCenterZero() override;

    // get maximum line increment between adjacent echoes (determines maximum PE-blip moment)
    virtual long getMaxLinIncrementBetweenEchoes() const;

    // get maximum partition increment between adjacent echoes (determines maximum 3D-blip moment)
    virtual long getMaxParIncrementBetweenEchoes() const;

    // printing ReorderInfo data to a string, which can be bassed to a trace.
    std::string printDataToString(const char* tText);


  protected:

    void resetBasic() override;

    virtual bool prepareCalculation_Standard(MrProt& rMrProt, SeqLim& rSeqLim);

    virtual bool prepareCalculation_PATRefScanEPI(MrProt& rMrProt, SeqLim& rSeqLim);

    virtual bool reorderEPI_Standard(MrProt&, SeqLim& rSeqLim);

    virtual bool reorderEPI_PATRefScanEPI(MrProt&, SeqLim& rSeqLim);

    /// Internal helper which calculates the reorder index for a given set of specifiers
    virtual long calcReorderIndex(long lEchoNo, long lContrastNo, long lShotNo, long lPartitionNo);

    void setParNo(long lDeltaReorderIndex, long lParNo);

    void calcMaxLinAndParIncrementBetweenEchoes();

    bool calcReorderingTable(long lShot = 0, long lPartition = 0, long lContrast = 0);

    void setSliceAccelerationParams(MrProt& rMrProt);

    EPIAcquisitionMode m_eAcquisitionMode{EPIAcquisitionMode::SingleShot};

    long m_lMeasuredLinesPerSegment{0};

    long m_lKSpaceCenterSegment{0};

    long m_lCounterInSegmentWithKSpaceCenter{0};

    bool m_bPATMultiShotRefScans{true};

    bool m_bPATGRERefScans{false};

    long m_lPATEchoTrainLengthRefScans{0};

    long m_lPATRefKSpaceCenterSegment{0};

    long m_lPATRefCounterInSegmentWithKSCenter{0};

    bool m_bIsReversePE{false};

    /// true if the reordering is prepared for slice-acceleration
    bool m_bSliceAccelerationActive{false};

    /// slice-acceleration factor considered within preparation
    long m_lSliceAccelerationFactor{0};

    /// FOV-shift factor considered within preparation
    long m_lFOVShiftFactor{0};

    /// maximum line increment between adjacent echoes
    long m_lMaxLinIncrementBetweenEchoes{0};

    /// maximum partition increment between adjacent echoes
    long m_lMaxParIncrementBetweenEchoes{0};

    /// true if the current run-mode is set to single-band
    bool m_bRunModeSingleBandActive{false};


  private:

    // checks for forbidden 3D parameters
    bool check3DParams(MrProt& rMrProt, SeqLim& rSeqLim);
};

#ifdef _MSC_VER
#pragma warning(pop)
#endif

inline void ReorderInfoEPI::setModeMultiShot()
{
    m_eAcquisitionMode = EPIAcquisitionMode::MultiShot;
}

inline void ReorderInfoEPI::setModeSingleShot()
{
    m_eAcquisitionMode = EPIAcquisitionMode::SingleShot;
}

inline bool ReorderInfoEPI::isMultiShot() const
{
    return m_eAcquisitionMode == EPIAcquisitionMode::MultiShot;
}

inline bool ReorderInfoEPI::isSingleShot() const
{
    return m_eAcquisitionMode == EPIAcquisitionMode::SingleShot;
}

inline long ReorderInfoEPI::getContrastsToMeasure() const
{
    return m_lMaxEchoNumber + 1;
}

inline long ReorderInfoEPI::getKSpaceCenterSegment() const
{
    return m_lKSpaceCenterSegment;
}

inline long ReorderInfoEPI::getCounterInSegmentWithKSpaceCenter() const
{
    return m_lCounterInSegmentWithKSpaceCenter;
}

inline bool ReorderInfoEPI::isPATRefScan(long DeltaReorderIndex)
{
    if (m_bPATActive && isSingleShot() && !isPATGRERefScans())
    {
        if (getValidDataIndex(DeltaReorderIndex) > m_lEchoTrainLength - 1)
            return true;
        else
            return false;
    }
    else
    {
        return false;
    }
}

inline void ReorderInfoEPI::setPATMultiShotRefScans(bool bValue)
{
    m_bPATMultiShotRefScans = bValue;
}

inline long ReorderInfoEPI::getPATEchoTrainLengthRefScans() const
{
    return m_lPATEchoTrainLengthRefScans;
}

inline long ReorderInfoEPI::getPATRefKSpaceCenterSegment() const
{
    return m_lPATRefKSpaceCenterSegment;
}

inline long ReorderInfoEPI::getPATRefCounterInSegmentWithKSCenter() const
{
    return m_lPATRefCounterInSegmentWithKSCenter;
}

inline void ReorderInfoEPI::setPATGRERefScans(bool bValue)
{
    m_bPATGRERefScans = bValue;
}

inline bool ReorderInfoEPI::isPATGRERefScans() const
{
    return m_bPATGRERefScans;
}

inline long ReorderInfoEPI::getPATLinesPerSegment() const
{
    if (m_bPATActive && (m_lPATAccelerationFactorPE > 0))
    {
        return m_lLinesPerSegment / m_lPATAccelerationFactorPE;
    }
    else
    {
        return m_lLinesPerSegment;
    }
}

inline bool ReorderInfoEPI::isFirstScanInSlice(long lDeltaReorderIndex)
{
    if (m_bPhaseCorrectionActive)
        return false;

    const long lIndex                 = getValidDataIndex(lDeltaReorderIndex);
    const long lFirstIndexInFirstShot = calcReorderIndex(0, 0, 0, 0);
    const long lLastIndexInFirstShot  = calcReorderIndex(getEchoTrainLength() - 1, getContrastsToMeasure() - 1, 0, 0);

    const bool bIsInFirstShot      = (lIndex >= lFirstIndexInFirstShot) && (lIndex <= lLastIndexInFirstShot);
    const bool bIsFirstEchoInTrain = (lIndex % getEchoTrainLength() == 0);

    return bIsInFirstShot && bIsFirstEchoInTrain;
};

inline bool ReorderInfoEPI::isLastScanInSlice(long lDeltaReorderIndex)
{
    if (m_bPhaseCorrectionActive)
        return false;

    const long lIndex                = getValidDataIndex(lDeltaReorderIndex);
    const long lNumberOfShots        = isPATActive() ? getPATLinesPerSegment() : getLinesPerSegment();
    const long lFirstIndexInLastShot = calcReorderIndex(0, 0, lNumberOfShots - 1, getPartitionsToMeasure() - 1);
    const long lLastIndexInLastShot  = calcReorderIndex(
        getEchoTrainLength() - 1, getContrastsToMeasure() - 1, lNumberOfShots - 1, getPartitionsToMeasure() - 1);

    const bool bIsInLastShot      = (lIndex >= lFirstIndexInLastShot) && (lIndex <= lLastIndexInLastShot);
    const bool bIsLastEchoInTrain = (lIndex % getEchoTrainLength() == getEchoTrainLength() - 1);

    return bIsInLastShot && bIsLastEchoInTrain;
}

inline void ReorderInfoEPI::setReversePE(bool bValue)
{
    m_bIsReversePE = bValue;
}

inline bool ReorderInfoEPI::getIsReversePE() const
{
    return m_bIsReversePE;
}

inline bool ReorderInfoEPI::isSliceAccelerationActive() const
{
    return m_bSliceAccelerationActive;
}

inline long ReorderInfoEPI::getSliceAccelerationFactor() const
{
    if (m_bSliceAccelerationActive)
    {
        return m_lSliceAccelerationFactor;
    }

    return 1;
}

inline long ReorderInfoEPI::getFOVShiftFactor() const
{
    if (m_bSliceAccelerationActive)
    {
        return m_lFOVShiftFactor;
    }

    return 1;
}

inline void ReorderInfoEPI::setParNo(long lDeltaReorderIndex, long lParNo)
{
    if (m_bSliceAccelerationActive)
    {
        // Limit partition index to the range [0, FOV-shift factor]
        m_aAllInfo[getValidDataIndex(lDeltaReorderIndex)].uiCPar
            = (unsigned short)std::max(0L, std::min(m_lFOVShiftFactor, lParNo));
    }
    else
    {
        // Call base class implementation
        ReorderInfo::setParNo(lDeltaReorderIndex, lParNo);
    }
}

inline void ReorderInfoEPI::setRunModeMultiBand()
{
    m_bRunModeSingleBandActive = false;
}

inline void ReorderInfoEPI::setRunModeSingleBand()
{
    m_bRunModeSingleBandActive = true;
}

inline long ReorderInfoEPI::getParNo(long lDeltaReorderIndex)
{
    // No slice partition encoding in single-band mode
    if (m_bSliceAccelerationActive && m_bRunModeSingleBandActive)
    {
        return 0;
    }

    // Call base class implementation
    return ReorderInfo::getParNo(lDeltaReorderIndex);
}

inline long ReorderInfoEPI::getParNoCenterZero(long lDeltaReorderIndex)
{
    // No slice partition encoding in single-band mode
    if (m_bSliceAccelerationActive && m_bRunModeSingleBandActive)
    {
        return 0;
    }

    // Call base class implementation
    return ReorderInfo::getParNoCenterZero(lDeltaReorderIndex);
}

inline long ReorderInfoEPI::getKSCenterPar()
{
    // No slice partition encoding in single-band mode
    if (m_bSliceAccelerationActive && m_bRunModeSingleBandActive)
    {
        return 0;
    }

    // Call base class implementation
    return ReorderInfo::getKSCenterPar();
}

inline bool ReorderInfoEPI::isKSCenterPar(long lDeltaReorderIndex)
{
    // No slice partition encoding in single-band mode
    if (m_bSliceAccelerationActive && m_bRunModeSingleBandActive)
    {
        return true;
    }

    // Call base class implementation
    return ReorderInfo::isKSCenterPar(lDeltaReorderIndex);
}

inline long ReorderInfoEPI::getMinParNoCenterZero()
{
    // No slice partition encoding in single-band mode
    if (m_bSliceAccelerationActive && m_bRunModeSingleBandActive)
    {
        return 0;
    }

    // Call base class implementation
    return ReorderInfo::getMinParNoCenterZero();
}

inline long ReorderInfoEPI::getMaxParNoCenterZero()
{
    // No slice partition encoding in single-band mode
    if (m_bSliceAccelerationActive && m_bRunModeSingleBandActive)
    {
        return 0;
    }

    // Call base class implementation
    return ReorderInfo::getMaxParNoCenterZero();
}

inline long ReorderInfoEPI::getMaxLinIncrementBetweenEchoes() const
{
    return m_lMaxLinIncrementBetweenEchoes;
}

inline long ReorderInfoEPI::getMaxParIncrementBetweenEchoes() const
{
    return m_lMaxParIncrementBetweenEchoes;
}

//} // end of namespace SEQ_NAMESPACE

#endif
