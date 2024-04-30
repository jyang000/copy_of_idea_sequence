
//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \src\MrImaging\seq\a_ep2d_se_mre\SBBEPIMREKernel.h
//     Lang: C++
//     Authors: Bolster Jr, Bradley
//              Kannengiesser, Stephan
//              Fang, Dong
//
//     Descrip: MR::MrServers::MrImaging::seq::a_ep2d_se_mre
//
//     Classes: SBBEPIKernelSEMRE
//              Derived from SBBEPIKernelSE
//
//    -----------------------------------------------------------------------------

#pragma once

#include "MrImaging/seq/a_ep2d_se/SBBEPIKernelSE.h"

#ifdef BUILD_SEQU
  #define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

namespace SEQ_NAMESPACE
{

class __IMP_EXP SBBEPIKernelSEMRE : public SBBEPIKernelSE
{
  public:

    //--------------------------------------------------------------------
    //  Constructor
      SBBEPIKernelSEMRE(SBBList* pSBBList);

    //--------------------------------------------------------------------
    //  Destructor
    virtual ~SBBEPIKernelSEMRE() = default;

    SBBEPIKernelSEMRE(const SBBEPIKernelSEMRE& right) = delete;
    SBBEPIKernelSEMRE& operator=(const SBBEPIKernelSEMRE& right) = delete;
    SBBEPIKernelSEMRE(SBBEPIKernelSEMRE&& right)                 = delete;
    SBBEPIKernelSEMRE& operator=(SBBEPIKernelSEMRE&& right) = delete;

    //--------------------------------------------------------------------
    //  Tells the kernel whether sync-bits should be sent
    //  before the excitation SBB or not. If so, an osc bit
    //  and an external trigger bit are prepared.
    //  Member m_lMaxSyncBitDuration is set.
    //  If the function executed with success, true is
    //  returned. The execution of the osc-bit can be dis-/en-abled
    //  during run-time of the sequence using the method setDoNotSendOscBit.
    //  The execution of the external trigger bit can be
    //  dis-/en-abled during run-time of the sequence using the method
    virtual bool setUseSyncBits (bool bValue, long lOscChannel = 0,     // channel for the osc bit
        long lOscDuration = 10,         // duration of osc bit
        long lOscStartTime = 0,         // start time of osc bit in event block
        long lExtTrigDuration = 10,     // duration of external trigger bit
        long lExtTrigStartTime = 0    // start time of external trigger bit in event block
        );


    //--------------------------------------------------------------------
    //  Prepares the kernel. At least an excitation SBB configuration
    //  function must have been registered. All other configuration
    //  steps have to be performed.
    virtual bool prepSBB (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

    //--------------------------------------------------------------------
    //  Internal/External function to run the sync-bit and the
    //  excitation SBB. Calls runSync Bits.
    virtual bool runExcitationAndSyncBits (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

    // run plug in
    virtual bool runRTEBPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

    /// run inside the Event block.
    virtual bool runInnerEvents(long lValue);

    // Methods for the prototype mre sequence
    // Set and get index of phase offsets from the kernel
    long getlMEGPhaseOffsetCounter();
    void setlMEGPhaseOffsetCounter(long lCtr);

    // Set and get number of MEG phase offsets from the kernel
    long getlNumMEGPhaseOffsets();
    void setlNumMEGPhaseOffsets(long lValue);

    void setbDefaultMode(bool);

    // Members for prototype mre sequence
    long m_lMEGPhaseOffset;

    // Added for communication between prep and run methods
    long m_lMEGPeriodBy4{0};

    double m_dMEGFrequencyHz{0.0};

    // For running external trigger implementation
    // negative values signal no external triggering
    sSYNC_EXTTRIGGER m_ExtTrigReadout{"ExTrig2"};

    double m_dReferenceTimeForExtTrigger_ms{0.0};
    long   m_lRTEBPlugInExtTrigStartTime;
    long   m_lExtTriggerDuration{0};
    long   m_lFirstExtTriggerTime{-1};
    long   m_lExtTriggerSpacing{-1};

    SBBRefocSE* getpSBBRefoc();

  protected:
    //--------------------------------------------------------------------
    // Data Members for Class Attributes

    // Index mechanical phase offsets
    long m_lNumMEGPhaseOffsets{0};
    long m_lMEGPhaseOffsetCounter{0};

  private:

    bool m_bIsDefaultMode{true};
};

// Members for prototype mre sequence
inline void SBBEPIKernelSEMRE::setbDefaultMode(bool flag)
{
    m_bIsDefaultMode = flag;
}

//-----------------------------------------------------------------
//  Set and get functions for passing repetition counter in and out
//  of kernel
//-----------------------------------------------------------------
inline long SBBEPIKernelSEMRE::getlMEGPhaseOffsetCounter()
{
    return m_lMEGPhaseOffsetCounter;
}
inline void SBBEPIKernelSEMRE::setlMEGPhaseOffsetCounter(long lctr)
{
    m_lMEGPhaseOffsetCounter = lctr;
}

//-----------------------------------------------------------------
//  Set and get functions for passing Number of MEG phase offsets
//  in and out of kernel
//-----------------------------------------------------------------
inline long SBBEPIKernelSEMRE::getlNumMEGPhaseOffsets()
{
    return m_lNumMEGPhaseOffsets;
}
inline void SBBEPIKernelSEMRE::setlNumMEGPhaseOffsets(long lvalue)
{
    m_lNumMEGPhaseOffsets = lvalue;
}

inline SBBRefocSE* SBBEPIKernelSEMRE::getpSBBRefoc()
{
    return m_pSBBRefocus.get();
}

}
