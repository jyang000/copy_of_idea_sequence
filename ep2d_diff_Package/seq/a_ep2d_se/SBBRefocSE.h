//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2016  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/X
//        File: \src\MrServers\MrImaging\seq\a_ep2d_se\SBBRefoc.h
//      Author: koellner
//              Thomas Kluge; Siemens AG Med MRIA/Seq; (09131) 84-8049
//              Uvo Hoelscher 
//
//        Lang: C++
//
//     Descrip: Declares SBB used in standard SE EPI sequences for slice selective
//            refocusing and FID-spoiling.
//
//    -----------------------------------------------------------------------------

#pragma once

#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpoRFBlockInfo.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sSLICE_POS.h"
#include "MrMeasSrv/SeqIF/libRT/sFREQ_PHASE.h"
#include "MrMeasSrv/SeqIF/libRT/sSYNC.h"
#include "MrImagingFW/libSeqUtilFW/KernelCalculationLimits.h"
#include "MrImagingFW/libSBBFW/SeqBuildBlock.h"

#ifdef EP2D_SE_MRE
#include "MrImaging/seq/greMRE/SBBMREMEG.h"
#endif

namespace SEQ_NAMESPACE
{



class SBBRefocSE: public SeqBuildBlock
{
public:
    SBBRefocSE  (SBBList* pSBBList = NULL);

    virtual ~SBBRefocSE() = default;

    void         setRFPulseForRefocusing(IRF_PULSE*);
    virtual bool prepSBB(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);
    virtual bool runSBB(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC);

    virtual IRF_PULSE* getRFPulsePointer();

    void setCalculationLimits(KernelCalculationLimits* pCalcLimits);

    double m_dROSpoilMomentFor10mmPixel{575.0};
    double m_dPESpoilMomentFactor{0.5};
    double m_dSSSpoilMomentFor10mmSlice{2000.0};
    double m_d3DSpoilMomentFor10mmPixel{575.0};
    bool   m_bROSpoilMomentNegative{false};
    bool   m_bRotationProofSpoilers{true};

    double m_dMinRiseTime{100.0};
    double m_dMaxMagnitude{1.0};

    long m_lTEContribution{0};
    long m_lStartTimeInEventBlock{0};
    long m_lFillEnd{0};

#ifdef EP2D_SE_MRE
    // * ------------------------------------------------------------------ *
    // * array of MEG SBB instances                                         *
    // * ------------------------------------------------------------------ *
    SBBMREMEG  m_aSBBMREMEG[NMAXMEGS];

    SBBMREMEG::eGradientAxis m_lMEGDirection{SBBMREMEG::Slice};
    long                     m_lMEGFractionalEncPerc{100};
    double                   m_dMEGFrequencyHz{0.0};

    long       getMEGTEContribution();

    //Get/Set functionality to move freeloop counter in and out of RefocRTEB
    long       getMEGInstanceCounter();
    void       setMEGInstanceCounter(long lctr);

    // For running external trigger implementation
    sSYNC_EXTTRIGGER m_ExtTrig{"ExTrig1"};
    long             m_lExtTriggerDuration {0};
    long             m_lFirstExtTriggerTime{-1};
    long             m_lExtTriggerSpacing{-1};

private:
    long       m_lMEGInstanceCounter;
#endif

protected:

    bool                                       checkCalcLimitsPointer();
    virtual void                               setGradientMinRiseTimesAndMaxAmplitudes();
    virtual bool                               prepareRFPulse(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);
    virtual void                               calcSliceSelectionGrad(MrProt& rMrProt);
    virtual void                               adaptToDoubleGRT();
    virtual bool                               checkAndAdaptForGradientReversal();
    virtual std::tuple<double, double, double> calcSpoilMoments(MrProt& rMrProt) const;
    virtual double                             calcMaxMomentForSpoilerTiming(
                                  MrProt& rMrProt, double dGPSpoilMoment, double dGRSpoilMoment, double dGSSpoilMoment);
    virtual void adaptSpoilMomentsForTGSE(double& dGPSpoilMoment, double& dGRSpoilMoment, double& dGSSpoilMoment);
    virtual bool setupSpoilerGradients(double dCalcMoment);
    virtual long getSliceSelectionGradientRampUpTime();
    virtual void adaptSliceSelectionGradientTiming();
#ifdef EP2D_SE_MRE
    bool setupMEGSBB(SeqLim& rSeqLim);
    bool prepareMEGSBB(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);
    bool checkMEGSBB(SeqLim& rSeqLim);
    bool prepareExternalTriggerForMRE();
    bool runSBBMREMEG(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC);
#endif
    virtual bool prepareGradients(SeqLim& rSeqLim, double dGSSpoilMoment, double dGRSpoilMoment, double dGPSpoilMoment);
    virtual bool prepareSpoilerGradients(SeqLim& rSeqLim, double dGPSpoilMoment, double dGRSpoilMoment, double dGSSpoilMoment);
    virtual void setExports();
    virtual void setEventStartTimes();
    virtual bool checkGradients(SeqLim& rSeqLim);
    virtual bool checkSpoilerGradients(SeqLim& rSeqLim);

    virtual bool run_insertEvents(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC);
    virtual bool checkRFPulse();
    virtual bool rePrepareSpoilerGradient();
    virtual void prepareNCOEvents(sSLICE_POS* pSLC);
    virtual bool runRTEvents(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC);
    virtual long calcStartTimeOfFirstSpoilerGradient();
    virtual long calcStartTimeOfSecondSpoilerGradient();
    void         runSpoilerGradients(long lFirstSpoilerStartTime, long lSecondSpoilerStartTime);
    virtual bool runRefocusing(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC);

    IRF_PULSE*  m_pRF{nullptr};
    sGRAD_PULSE m_GS{"SERefocRTEB::m_GS"};
    sGRAD_PULSE m_GSSpoil{"SERefocRTEB::m_GSSpoil"};
    sGRAD_PULSE m_GRSpoil{"SERefocRTEB::m_GRSpoil"};
    sGRAD_PULSE m_GPSpoil{"SERefocRTEB::m_GPSpoil"};
    long        m_lShiftRF{0};

    //Static variables in the old code. They're now declared as member variables
    //to avoid vilation between threads.
    sFREQ_PHASE m_sRFSet{"sRFSet"};
    sFREQ_PHASE m_sRFNeg{"sRFNeg"};

};

inline void SBBRefocSE::setRFPulseForRefocusing (IRF_PULSE* pRF)
{
    m_pRF = pRF;
    resetPrepared();
}

#ifdef EP2D_SE_MRE
//------------------------------------------------------------
// SERefocRTEB class: getMEGTEContribution()
//------------------------------------------------------------
inline long SBBRefocSE::getMEGTEContribution()
{
    return m_aSBBMREMEG[0].getDurationPerRequest();
}
//------------------------------------------------------------
// SERefocRTEB class: getMEGInstanceCounter()
//------------------------------------------------------------
inline long SBBRefocSE::getMEGInstanceCounter()
{
    return m_lMEGInstanceCounter;
}
//------------------------------------------------------------
// SERefocRTEB class: setMEGInstanceCounter()
//------------------------------------------------------------
inline void SBBRefocSE::setMEGInstanceCounter(long lctr)
{
    m_lMEGInstanceCounter = lctr;
}
#endif

}//end of namespace SEQ_NAMESPACE
