//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/X
//        File: \src\MrImaging\seq\a_tgse_asl\SBBGrase3DKernel.h
//      Author: PLM AW NERUO
//        Date: 2018-08-07 12:24:16 +02:00
//
//        Lang: C++
//
//
//
///  \file   SBBGrase3DKernel.h
///  \brief  File containing declaraion of the SBBGrase3DKernel class
///        
///
///  This file contains the implementation of the class SBBGrase3DKernel.
///
//    -----------------------------------------------------------------------------

#pragma once

// SBBEPIKernel
#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include "MrImaging/seq/SeqDebug.h"
#include "MrImagingFW/libSeqUTIF/libsequt.h"                // for mSEQTest
#include "MrImaging/seq/a_ep2d_se/SBBRefocSE.h"

#define DEBUG_ORIGIN 0x00800000

using namespace SEQ_NAMESPACE;

class SBBGrase3DKernel : public SeqBuildBlockEPIKernel
{
public:
    SBBGrase3DKernel (SBBList* pSBBList);
    
    //## Destructor (generated)
    virtual ~SBBGrase3DKernel();

    //SBBEPIKernel Section
    //overide--------------------------------------------------------------
    //    Prepares the kernel. At least an excitation SBB configuration
    //  function must have been registered. All other configuration
    //  steps have to be performed.
    bool prep (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

    //overide--------------------------------------------------------------
    //    Executes the real time part of the SBB.
    bool run (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC, long counter = 0);

    //overide--------------------------------------------------------------
    //    Internal function to run the phase correction EPI read-out and the
    //    appropriate fill times, if the kernel is in the phase correction mode.
    bool runEPIReadOutPhaseCorrection (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

    //overide-------------------------------------------------------------
    //    Internal function to run the sync-bit and the
    //  excitation SBB. Calls runSync Bits.
    bool runExcitationAndSyncBits (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

    // preparation of the SBBExcitation
    virtual bool initExcitation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

    // preparation of the plug in SBB
    virtual bool prepPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);
	
    //overload-----------------------------------------------------------
    //    Internal function to run the RTEB-plug-in.
    bool runRTEBPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC, long lCounter);

    //overide-------------------------------------------------------------
    //  Calculate: (1) An effective echo-spacing, which takes number of interleaves
    //                 and PAT acceleration factor into account
    //             (2) The bandwidth per pixel in the phase-encoding direction
    //                 FOR THE RECONSTRUCTED IMAGE
    //
    virtual bool calcEffEchoSpacingAndBWPerPixelPE(MrProt &/*rMrProt*/, long& /*lEffectiveEchoSpacing*/, double& /*dBandwidthPerPixelPE */)
    {
      return true;
    };

    //SBBReadout Section
    //overide-------------------------------------------------------------
    //  Returns the 2d phase encoding moment necessary to 
    //  RE-phase the magnetization
    //  after the measurement of the last ADC.
    double getPERePhasingMoment (MrProt &rMrProt, SeqLim& );

    //overide-------------------------------------------------------------
    //  Returns the 2d phase encoding moment necessary to 
    //  PRE-phase the magnetization before the measurement of the first ADC.
    double getPEPrePhasingMoment (MrProt &rMrProt, SeqLim& );

    //overide-------------------------------------------------------------
    //  Returns the number of RO-gradients applied for the 
    //  EPI-like phase correction scan.
    //  This number may be higher than the number of echos 
    //  acquired for phase correction.
    long getEchoTrainLengthPhaseCorrScan ();

    //overide-------------------------------------------------------------
    //  Internal function called by the run-function to set the Mdh-entries which
    //  change from ADC to ADC.
    bool run_setADCMdh (long lADCCounter, MdhProxy* pMdh, bool bIsLastADC);

    //overide-------------------------------------------------------------
    //  Returns the duration of the SBB in us when executed as imaging read out.
    long getDurationEPIReadOutPerRequest ();

    //New methods
    void set_CycleLength(long len)
    {
      m_CycleLength = len;
    }

#ifdef SUPPORT_CSL
    void setOverwriteLastScanInMeasAndConcat(bool bVal)   { m_bOverwriteIsLastScanInMeasAndConcat = bVal; };
    void setIsLastScanInMeas(bool bVal)                   { m_bIsLastScanInMeas = bVal; };
    void setIsLastScanInConcat(bool bVal)                 { m_bIsLastScanInConcat = bVal; };
#endif

    long getTurboFactor() const { return m_lTurboFactor; }
    void setTurboFactor(long lTurboFactor);

    // energy export
    MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest (); 

    // excitation pulse
    sRF_PULSE_SINC m_sExcitationRF;

    // refocusing pulse
    sRF_PULSE_EXT m_sRefocusingRF;

    // refocusing SBB
    SBBRefocSE m_SBBRefocus;
	
private:
    long m_CycleLength; // that's not nice in here!!
    long m_lCounter;
	  long m_lTurboFactor;
  #ifdef SUPPORT_CSL
    bool m_bOverwriteIsLastScanInMeasAndConcat;
    bool m_bIsLastScanInMeas;
    bool m_bIsLastScanInConcat;
  #endif
};

inline long SBBGrase3DKernel::getDurationEPIReadOutPerRequest ()
{
  //## begin SeqBuildBlockEPIReadOut::getDurationEPIReadOutPerRequest%36C4259F03D4.body preserve=yes
  // *****************************************************************************
  // **                                                                         **
  // **                                                                         **
  // ** FUNCTION:  getDurationEPIReadOutPerRequest                              **
  // **                                                                         **
  // **                                                                         **
  // *****************************************************************************
  long lRet = m_lEchoSpacing*m_lEchoTrainLength - m_lRampTimeOutsideSBB_us;

  if (m_lCountersPerSegmentForEchoShifting)
  {
      long lEchoShiftingDelay  =0;
      long lEchoShiftingFillEnd=0;

      if (getEchoShiftingData (0, &lEchoShiftingDelay, &lEchoShiftingFillEnd))
      {
          lRet += lEchoShiftingDelay+lEchoShiftingFillEnd;
      }
  }

  return lRet;
  //## end SeqBuildBlockEPIReadOut::getDurationEPIReadOutPerRequest%36C4259F03D4.body
}


inline double SBBGrase3DKernel::getPERePhasingMoment (MrProt &rMrProt, SeqLim& )
{
	double dFovPH = rMrProt.sliceSeries().front().getdPhaseFOV(); 
	double dPhaseOSFactor = rMrProt.phaseOversampling();
	const char* nucleus = rMrProt.getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str();
    return - fGSLGetPEDeltaMoment(dFovPH, dPhaseOSFactor, nucleus) * static_cast<double>(m_pRI->getLinNoCenterZero(m_lEchoTrainLength-1));
//  return fGSLGetPEDeltaMoment(pMrProt) * (m_pRI->getLinNoCenterZero(0));
//  return fGSLGetRODeltaMoment(pMrProt) * (pMrProt->fastImaging().EPIFactor()) - getROPrePhasingMoment(pMrProt,pSeqLim);
}

inline double SBBGrase3DKernel::getPEPrePhasingMoment (MrProt &rMrProt, SeqLim& )
{
	double dFovPH = rMrProt.sliceSeries().front().getdPhaseFOV();
	double dPhaseOSFactor = rMrProt.phaseOversampling();
	const char* nucleus = rMrProt.getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str();
	return fGSLGetPEDeltaMoment(dFovPH, dPhaseOSFactor, nucleus) * static_cast<double>(m_pRI->getLinNoCenterZero(0));
}

inline long SBBGrase3DKernel::getEchoTrainLengthPhaseCorrScan ()
{
  return 3+m_lAddROGradientsBeforeEPIPhaseCorrScans+m_lAddROGradientsAfterEPIPhaseCorrScans;
}

// we overwrite the behaviour of the original SBBEPIKernel to avoid an excitation pulse for each EPI echo train
inline bool SBBGrase3DKernel::runExcitationAndSyncBits (MrProt &/*rMrProt*/, SeqLim &/*rSeqLim*/, SeqExpo &/*rSeqExpo*/, sSLICE_POS* /*pSLC*/)
{
    return true;
}

