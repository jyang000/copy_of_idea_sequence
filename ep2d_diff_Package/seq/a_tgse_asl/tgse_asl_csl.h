//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2013  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/X
//	    File: \src\MrImaging\seq\a_tgse_asl\tgse_asl_csl.h
//	  Author: pfeujodj
//	    Date: 2018-08-02 08:20:30 +02:00
//
//	    Lang: C++
//
//	 Descrip: 3D ASL sequence with CompositeSeqLoop support
//
//	-----------------------------------------------------------------------------
#pragma once

#include "MrImagingFW/libSBBFW/StdSeqIF.h"
#include "MrGlobalDefinitions/MrResult.h"
#include "MrImaging/libSeqUtil/libSeqUtil.h"

#include "MrImaging/libSBB/libSBBmsg.h"                  // SBB_... error codes
#include "MrImaging/libSL/StdSL.h"
#include "MrImaging/libSL/StdSL_ID.h"
#include "MrImaging/libSL/StdMediator.h"
#include "MrImaging/libSL/StdLoopIdents.h"

#ifdef WIN32
#include "TCHAR.h"
#endif

#include "MrImaging/libSBB/SBBBinomialPulses.h"
//#include "MrImaging/seq/epi_Config.h"
#include "MrImagingFW/libCSL/Key.h"

#include "MrImaging/seq/a_tgse_asl/EpiKernelLeafAdapter.h"
#include "MrImaging/seq/a_tgse_asl/AslLeafAdapter.h"
#include "MrImaging/seq/a_tgse_asl/ReorderInfoGrase3D.h"
#include "MrImaging/seq/a_tgse_asl/SBBGrase3DKernel.h"
#include "MrImaging/seq/a_tgse_asl/a_tgse_asl_UI.h"
#include "MrImaging/seq/a_tgse_asl/ASLLoop.h"
#include "MrImaging/seq/a_tgse_asl/AslSL.h"
#include "MrImaging/seq/a_tgse_asl/ASLDefines.h"   // SEQ_ASL

// Forward declarations
class MrProt;
class SeqLim;


namespace SEQ_NAMESPACE
{
  class TGSEAslCsl
  {
  public:
    TGSEAslCsl(StdSeqIF *sequence);
    virtual ~TGSEAslCsl();

    virtual NLSStatus initialize(SeqLim &rSeqLim,Tgse_asl_UI* m_pUI);

    virtual NLSStatus prepare(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo,Tgse_asl_UI* m_pUI);

    virtual NLSStatus check(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo, SEQCheckMode *  pSEQCheckMode);

    virtual NLSStatus run(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo);

    virtual NLS_STATUS runKernel(MrProt &rMrProt,SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo, long lKernelMode, long lSlice, long lPartition, long lLine);

  protected:
    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  TGSEAslCsl::prepareSetLimits                            *
    // *                                                                        *
    // * Description :  Sets the parameters for the hard limits.                *
    // * Parameter   : - MrProt& rMrProt : Reference of the protocol strcture   *
    // *               - SeqLim& rSeqLim : Reference of the sequence limits     *
    // *                                   structure                            *
    // *               - SeqExpo& rSeqExpo : Reference of the sequence export   *
    // *                                   structure                            *
    // *               - bool bPerformM0Scan: bool specifying if M0Scan should  *
    // *                                    be performed                        *
    // *                                                                        *
    // * Return      :  NLSStatus                                               *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    virtual NLSStatus prepareSetLimits(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  TGSEAslCsl::prepareSeqExpo                              *
    // *                                                                        *
    // * Description :  Sets the parameters that get exported for the           *
    // *                measurement with the given protocol.                    *
    // * Parameter   : - MrProt& rMrProt : Reference of the protocol strcture   *
    // *               - SeqLim& rSeqLim : Reference of the sequence limits     *
    // *                                   structure                            *
    // *               - SeqExpo& rSeqExpo : Reference of the sequence export   *
    // *                                   structure                            *
    // *               - bool bPerformM0Scan: bool specifying if M0Scan should  *
    // *                                    be performed                        *
    // *                                                                        *
    // * Return      :  NLSStatus                                               *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    virtual NLSStatus prepareSeqExpo(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, bool bPerformM0Scan);

    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  TGSEAslCsl::prepareCompositeSeqLoop                     *
    // *                                                                        *
    // * Description :  Prepares the mediator keys according to the protocol    *
    // *                and calls prepareAslSL afterwards.                      *
    // * Parameter   : - MrProt& rMrProt : Reference of the protocol strcture   *
    // *               - SeqLim& rSeqLim : Reference of the sequence limits     *
    // *                                   structure                            *
    // *               - SeqExpo& rSeqExpo : Reference of the sequence export   *
    // *                                   structure                            *
    // *               - bool bPerformM0Scan: bool specifying if M0Scan should  *
    // *                                    be performed                        *
    // *                                                                        *
    // * Return      :  NLSStatus                                               *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    virtual NLSStatus prepareCompositeSeqLoop (MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, bool bPerformM0Scan);

    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  TGSEAslCsl::prepareAslSL                                *
    // *                                                                        *
    // * Description :  Parameterizes the AslSL and adds the readout and        *
    // *                preparation LeafAdapters                                *  
    // * Parameter   : - MrProt& rMrProt : Reference of the protocol strcture   *
    // *               - SeqLim& rSeqLim : Reference of the sequence limits     *
    // *                                   structure                            *
    // *               - SeqExpo& rSeqExpo : Reference of the sequence export   *
    // *                                   structure                            *
    // *               - bool bPerformM0Scan: bool specifying if M0Scan should  *
    // *                                    be performed                        *
    // *                                                                        *
    // * Return      :  NLSStatus                                               *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    virtual NLSStatus prepareAslSL(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, bool bPerformM0Scan);

    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  TGSEAslCsl::setIceProgram                               *
    // *                                                                        *
    // * Description :  Sets the correct IceProgram for given protocol.         *  
    // * Parameter   : - MrProt& rMrProt : Reference of the protocol strcture   *
    // *               - SeqLim& rSeqLim : Reference of the sequence limits     *
    // *                                   structure                            *
    // *               - SeqExpo& rSeqExpo : Reference of the sequence export   *
    // *                                   structure                            *
    // *                                                                        *
    // * Return      :  NLSStatus                                               *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    virtual NLSStatus setIceProgram(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  TGSEAslCsl::initializeReadoutKernel                     *
    // *                                                                        *
    // * Description :  Initializes the given EpiKernelLeafAdapter.             *  
    // * Parameter   : - EpiKernelLeafAdapter* epiKernelLeafAdapter: Pointer to *
    // *                              EpiKernelLeafAdapter to be initialized    *
    // *               - TgseCallbackConfigExcitation &excitationCallbackConfig:*
    // *                              Reference to excitisation callback config *
    // *               - SeqBuildBlockBinomialPulses& binomialPulses: Reference *
    // *                              to binomial pulses SBB.                   *
    // *               - TgseCallbackRTEBPlugIn &RTEBCallBackPlugin: Refence to *
    // *                              RTEBCallBackPlugin.                       *
    // *               - SBBRefocSE &refocRTEB: Reference to SBBRefocSE         *
    // *                                                                        *
    // * Return      :  bool                                                    *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    virtual NLSStatus initializeReadoutKernel(SL_NAME_SPACE::EpiKernelLeafAdapter* epiKernelLeafAdapter);

    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  TGSEAslCsl::performM0Scan                              *
    // *                                                                        *
    // * Description :  Returns if M0Scans are activated in the given protocol. *
    // *                                                                        *
    // * Return      :  bool                                                    *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    virtual bool performM0Scan(MrProt & rMrProt) const;

    StdSeqIF *m_pSequence;

    // --------------------------------------------------------------------------------------------
    ///  \brief    Standard mediator used for parameter exchange in CSL
    // --------------------------------------------------------------------------------------------
    SL_NAME_SPACE::StdMediator m_SLM;

    // --------------------------------------------------------------------------------------------
    ///  \brief    Composite SeqLoop, default implementation
    // --------------------------------------------------------------------------------------------
    SL_NAME_SPACE::AslSL m_AslSL;

    //-------------------------------------------------------------------------------------
    // Members for Epi-Readout
    //-------------------------------------------------------------------------------------
    SL_NAME_SPACE::EpiKernelLeafAdapter* m_pEpiKernelAdapter;
    SL_NAME_SPACE::EpiKernelLeafAdapter* m_pEpiKernelAdapter_M0;

    SL_NAME_SPACE::AslLeafAdapter* m_pAslLeafAdapter;
    SL_NAME_SPACE::AslLeafAdapter* m_pAslLeafAdapter_M0;

    // reordering information data
    ReorderInfoGrase3D m_REOInfo;
    ReorderInfoGrase3D m_REOInfo_M0;

    // kernel calculation limits
    KernelCalculationLimits m_myCalcLimits;

    //-------------------------------------------------------------------------------------
    //long     m_lLines;
    //long     m_lPartitions;

    int32_t  m_lKSpacePointsToMeasure;
    int32_t  m_lKSpacePointsToCheck;
    long     m_lSlicesToMeasure;
    double   m_dMinRiseTime;
    double   m_dGradMaxAmpl;
    SBBList m_dummySBBList;
    SBBList m_dummySBBList_M0;

    // Slice position information (rotation matrices and shifts)
    sSLICE_POS m_asSLC [K_NO_SLI_MAX];
    //epiConfig m_configInfo;

    MrProt m_m0MrProt;

  private:
    //  Copy constructor not implemented
    TGSEAslCsl (const TGSEAslCsl &right);

    // Assignment operator not implemented  
    TGSEAslCsl & operator=(const TGSEAslCsl &right);
  };

} // SEQ_NAMESPACE

