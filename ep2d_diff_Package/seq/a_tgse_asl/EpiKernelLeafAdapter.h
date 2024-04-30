//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2013  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/X
//	    File: \src\MrImaging\seq\a_tgse_asl\EpiKernelLeafAdapter.h
//	  Author: pfeujodj
//	    Date: 2018-08-06 13:59:07 +02:00
//
//	    Lang: C++
//
//	 Descrip: 3D ASL sequence with CompositeSeqLoop support
//
//	-----------------------------------------------------------------------------

#pragma once

#include "MrImagingFW/libCSL/LeafAdapter.h"
#include "MrImagingFW/libCSL/Counter.h"
#include "MrImagingFW/libCSL/Key.h"
#include "MrImaging/libSL/StdSL_ID.h"

#ifdef TGSE
  #include "MrImaging/seq/a_tgse_asl/ReorderInfoGrase3D.h"
  #include "MrImaging/seq/a_tgse_asl/SBBGrase3DKernel.h"
  #define EPIKERNEL SBBGrase3DKernel
#else
  #include "MrImaging/seq/Kernels/SBBEPIKernel.h"
  #define EPIKERNEL SeqBuildBlockEPIKernel
#endif
                                                  
#include "MrGlobalDefinitions/ImpExpCtrl.h"    // import/export control

class ReorderInfo;
class MrProt;
class SeqLim;

namespace MrProtocolData { class SeqExpo; }

namespace SL_NAME_SPACE  
{

// ----------------------------------------------------------------------------
/// \brief EpiKernelLeafAdapter
///
/// Leafadapter that encapsulates the SBBGrase3DKernel setting the required parameters read from the mediator.
/// 
/// In the context \c SL::PHASE_CORR the loop length is automatically set to one and providing the interface of the SBB.

///
/// \par Default Values
/// - none
///
/// \par Exported Values
/// The phase loop exports its loop counter using the mediator key 
/// \c m_LoopCounterKey. This key can be specified as parameter of the
/// constructor. Its default value is \c SLK::l_PHASE_LOOP_COUNTER.
///
/// \par Imported Values
/// The phase loop imports the loop length from the mediator using the key 
/// given by \c m_LoopLengthKey. This key can be specified as parameter of the
/// constructor. Its default value is \c SLK::l_PHASE_LOOP_LENGTH.
///
// ----------------------------------------------------------------------------
class EpiKernelLeafAdapter : public LeafAdapter<EPIKERNEL>  
{
  public:
    // * ------------------------------------------------------------------ *
    // * Keys used in the EpiKernelLeafAdapter                                     *
    // * ------------------------------------------------------------------ *
    static const Key<long> l_COOL_PAUSE_EXPLICIT;

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::EpiKernelLeafAdapter          *
    // *                                                                    *
    // * Description :  Constructor                                         *
    // *                                                                    *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    EpiKernelLeafAdapter (const std::string& sIdent, Mediator* pMediator,ReorderInfoGrase3D &reorderInfo);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::EpiKernelLeafAdapter          *
    // *                                                                    *
    // * Description :  Constructor                                         *
    // *                                                                    *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    EpiKernelLeafAdapter (const std::string& sScope, const std::string& sIdent, Mediator* pMediator,ReorderInfoGrase3D &reorderInfo);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::~EpiKernelLeafAdapter         *
    // *                                                                    *
    // * Description :  Destructor                                          *
    // *                                                                    *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    virtual ~EpiKernelLeafAdapter();

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::prep                          *
    // *                                                                    *
    // * Description :  Calculates the timing of the kernel and             *
    // *                prepares all real time events of the SBB            *
    // *                                                                    *
    // * Parameter   : - MrProt& rMrProt : Pointer to the protocol strcture *
    // *               - SeqLim& rSeqLim : Pointer to the sequence limits   *
    // *                                   structure                        *
    // *               - SeqExpo& rSeqExpo : Pointer to the sequence export *
    // *                                     structure                      *
    // *                                                                    *
    // * Return      :  bool                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    virtual bool prep (MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::run                           *
    // *                                                                    *
    // * Description :  Executes the Grase3D kernel                         *
    // *                                                                    *
    // * Parameter   : - MrProt& rMrProt : Pointer to the protocol strcture *
    // *               - SeqLim& rSeqLim : Pointer to the sequence limits   *
    // *                                   structure                        *
    // *               - SeqExpo& rSeqExpo : Pointer to the sequence export *
    // *                                     structure                      *
    // *               - sSLICE_POS* pSLC : Pointer to the rotation matrix  *
    // *                                    and slice position information  *
    // *                                                                    *
    // * Return      :  bool                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    virtual bool run (MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, sSLICE_POS* pSLC);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::initKernel                    *
    // *                                                                    *
    // * Description :  Initializes the Grase3D kernel SBB.                 *
    // *                                                                    *
    // * Parameter   :                                                      *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    void initKernel();

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::setPrepareParameters          *
    // *                                                                    *
    // * Description :  Sets the parameters required during prepation       *
    // *                Returns true on success, false otherwise            *
    // *                                                                    *
    // * Parameter   : -MrProt& rMrProt : Reference of the protocol strcture*
    // *               -SeqLim& rSeqLim : Reference of the sequence limits  *
    // *                                   structure                        *
    // *               -SeqExpo& rSeqExpo : Reference of the sequence export*
    // *                                     structure                      *
    // *               -ReorderInfoGrase3D &rREOInfo: Reference of          *
    // *                                               ReorderInfoGrase3D   *
    // *                                                                    *
    // * Return      :  bool                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    bool setPrepareParameters(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, ReorderInfoGrase3D &rREOInfo);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::initializeRun                 *
    // *                                                                    *
    // * Description :  Updates the parameter of the SBB that have to be    *
    // *                set before SBB.run is executed.                     *
    // *                                                                    *
    // * Parameter   : -MrProt& rMrProt : Reference of the protocol strcture*
    // *               -SeqLim& rSeqLim : Reference of the sequence limits  *
    // *                                   structure                        *
    // *               -SeqExpo& rSeqExpo : Reference of the sequence export*
    // *                                     structure                      *
    // *               -ReorderInfoGrase3D &rREOInfo: Reference of          *
    // *                                               ReorderInfoGrase3D   *
    // *                                                                    *
    // *              -long& lCurSegment: Reference to current Segment      *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    void initializeRun(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, long lSlice,
      ReorderInfoGrase3D &rREOInfo, long& lCurSegment);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::postRun                       *
    // *                                                                    *
    // * Description :  Executes code after SBB is executed, i.e.           *
    // *                fSBBCooltime, handling of freq_feedback, if required*
    // *                                                                    *
    // * Parameter   : -MrProt& rMrProt : Reference of the protocol strcture*
    // *               -SeqLim& rSeqLim : Reference of the sequence limits  *
    // *                                   structure                        *
    // *               -SeqExpo& rSeqExpo: Reference of the sequence export *
    // *                                     structure                      *
    // *               - long lSlice    : Current slice                     *
    // *               -ReorderInfoGrase3D &rREOInfo: Reference of          *
    // *                                               ReorderInfoGrase3D   *
    // *                                                                    *
    // *              -long& lCurSegment: Reference to current Segment      *
    // *              -long lCoolPauseExplicit: Cooling pause in ms.        *
    // *              -int64_t lTRfillTime: Additional fill time            *
    // * Return      :  NLS_STATUS                                          *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    NLS_STATUS postRun(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, long lSlice,
      ReorderInfoGrase3D &rREOInfo, long lCurSegment,
      long lCoolPauseExplicit,int64_t lTRFillTime);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::initializeReorderInfo         *
    // *                                                                    *
    // * Description :  Initializes the given ReorderInfoGrase3D for the    *
    // *                given protocol.                                     *
    // *                                                                    *
    // * Parameter   : -MrProt& rMrProt : Reference of the protocol strcture*
    // *               -SeqLim& rSeqLim : Reference of the sequence limits  *
    // *                                   structure                        *
    // *               -SeqExpo& rSeqExpo: Reference of the sequence export *
    // *                                     structure                      *
    // *               -ReorderInfoGrase3D &rREOInfo: Reference of          *
    // *                                               ReorderInfoGrase3D   *
    // *                                                                    *
    // * Return      :  NLS_STATUS                                          *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    NLS_STATUS initializeReorderInfo(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo,ReorderInfoGrase3D &rREOInfo);

    // * Following functions are just used to provide the interface of the  *
    // * SBB have a look at documentation of SBB for details. @{            *

    inline long getNeededTE()
    {
      return m_SBB.getNeededTE();
    }

    inline void set_CycleLength(long len)
    {
      m_SBB.set_CycleLength(len);
    }

    inline bool getUseGPABalance()
    {
      return m_SBB.getUseGPABalance();
    }

    inline sGRAD_PULSE_RO getGRO()
    {
      return m_SBB.getGRO();
    }

    inline long getEchoSpacing()
    {
      return m_SBB.getEchoSpacing();
    }

    inline bool setRequestsPerMeasurement(long Requests)
    {
      return m_SBB.setRequestsPerMeasurement(Requests);
    }

//    inline  MrProtocolData::SeqExpoRFInfo getRFInfoPerMeasurement (long lSlicesPerMeasurement)
//    {
//      return m_SBB.getRFInfoPerMeasurement(lSlicesPerMeasurement);
//    }

    inline long getTRIncrement()
    {
      return m_SBB.getTRIncrement();
    }

    inline sREADOUT* getReadOutAddress ()
    {
      return m_SBB.getReadOutAddress();
    }

    inline bool calcEffEchoSpacingAndBWPerPixelPE( MrProt &rMrProt, long& lEffectiveEchoSpacing, double& dBandwidthPerPixelPE )
    {
      return m_SBB.calcEffEchoSpacingAndBWPerPixelPE(rMrProt, lEffectiveEchoSpacing,dBandwidthPerPixelPE);;
    };

    inline void clearOscBits()
    {
      for (long lI=0; lI<m_lMaxOscBitSentFlags; lI++)
      {
        m_abOscBitSentForMeas[lI]=false;
      }
    }

    inline bool setPointerToCalculationLimits(KernelCalculationLimits* pPointer)
    {
      return m_SBB.setPointerToCalculationLimits(pPointer);
    }

    inline NLS_STATUS getKernelNLSStatus()
    {
      return m_SBB.getNLSStatus();
    }

    inline bool setUseSyncBits(
      bool bValue, 
      long lOscChannel = 0,     // channel for the osc bit
      long lOscDuration = 10,         // duration of osc bit
      long lOscStartTime = 0,         // start time of osc bit in event block
      long lExtTrigDuration = 10,     // duration of external trigger bit
      long lExtTrigStartTime = 0    // start time of external trigger bit in event block
      )
    {
      return m_SBB.setUseSyncBits(bValue,lOscChannel,lOscDuration,lOscStartTime,
        lExtTrigDuration, lExtTrigStartTime);
    }

//    inline bool setEPICallbackConfigExcitation(EpiCallbackConfigExcitation *pEpiCallbackConfigExcitation)
//    {
//      return m_SBB.setEPICallbackConfigExcitation(pEpiCallbackConfigExcitation);
//    }
//    inline void setEPICallbackRTEBPlugIn(EpiCallbackRTEBPlugIn *pEpiCallbackRTEBPlugIn)
//    {
//      m_SBB.setEPICallbackRTEBPlugIn(pEpiCallbackRTEBPlugIn);
//    }

    inline void setPrephaseBlipsAfterRTEBPlugIn(bool bValue)
    {
      m_SBB.setPrephaseBlipsAfterRTEBPlugIn(bValue);
    }

    inline void setPrephaseROAfterRTEBPlugIn(bool bValue)
    {
      m_SBB.setPrephaseROAfterRTEBPlugIn(bValue);
    }
    //@}

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::getSBB                        *
    // *                                                                    *
    // * Description :  Returns pointer to SBBGrase3DKernel                 *
    // *                                                                    *
    // * Parameter   :                                                      *
    // *                                                                    *
    // * Return      :  SBBGrase3DKernel*                                   *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    inline SBBGrase3DKernel* getSBB()
    {
      return &m_SBB;
    }

  protected:
    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  EpiKernelLeafAdapter::initialize                    *
    // *                                                                    *
    // * Description :  Initializes the Grase3D kernel                      *
    // *                                                                    *
    // * Parameter   :                                                      *
    // *                                                                    *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    void initialize();

    ReorderInfoGrase3D &m_rReorderInfo;

    //-------------------------------------------------------------------------------------
    // control execution of osc-bit
    //-------------------------------------------------------------------------------------
    enum{ m_lMaxOscBitSentFlags = 4096 };
    bool m_abOscBitSentForMeas[m_lMaxOscBitSentFlags];

  private:
    EpiKernelLeafAdapter (const EpiKernelLeafAdapter& right);   // Copy constructor not implemented 
    EpiKernelLeafAdapter& operator= (const EpiKernelLeafAdapter &right);  // Assignment operator not implemented
};

}   // namespace SL_NAME_SPACE
