//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2013  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/X
//	    File: \src\MrImaging\seq\a_tgse_asl\ASLLoop.h
//	  Author: pfeujodj
//	    Date: 2017-12-11 09:56:08 +01:00
//
//	    Lang: C++
//
//	 Descrip: 3D ASL sequence with CompositeSeqLoop support
//
//	-----------------------------------------------------------------------------

#pragma once

#include "MrImaging/seq/a_tgse_asl/AslSL_ID.h"
#include "MrImagingFW/libCSL/CompositeLoop.h"
#include "MrImaging/libSL/StdSL_ID.h"

// Forward declaration
class MrProt;
class SeqLim;

namespace MrProtocolData { class SeqExpo; }

namespace SL_NAME_SPACE  
{

// ----------------------------------------------------------------------------
/// \brief ASLLoop
///
/// The length of this loop is retrieved from the mediator using the key 
/// specified as a parameter of the constructor. If there is no value given the
/// default value \c SL::AslKeys::::l_ASL_LOOP_LENGTH is used. The current loop counter
/// can be accessed using the key specified as a parameter of the constructor. 
/// If there is no value given the default value \c SL::AslKeys::::l_ASL_LOOP_COUNTER
/// is used. \n
/// 
/// \par Default Values
/// - none
///
/// \par Exported Values
/// The asl loop exports its loop counter using the mediator key 
/// \c m_LoopCounterKey. This key can be specified as parameter of the
/// constructor. Its default value is \c SL::AslKeys::::l_ASL_LOOP_COUNTER.
///
/// \par Imported Values
/// The phase loop imports the loop length from the mediator using the key 
/// given by \c m_LoopLengthKey. This key can be specified as parameter of the
/// constructor. Its default value is \c SL::AslKeys::::l_ASL_LOOP_LENGTH.
///
// ----------------------------------------------------------------------------
class ASLLoop : public CompositeLoop
{
  public:
    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::ASLLoop
    //
    // Description :
    ///              \brief     Constructor
    ///
    ///                         \c SL::AslKeys::::l_ASL_LOOP_LENGTH and \c SL::AslKeys::::l_ASL_LOOP_COUNTER 
    ///                         will be used as default keys for the loop length and
    ///                         loop counter respectively
    //
    // Parameters  :
    ///              \param     sIdent              String identifier of the component
    ///              \param     pMediator           Pointer to Mediator class
    ///              \param     rLoopLengthKey      Key specifying the loop length
    ///              \param     rLoopCounterKey     Key specifying the loop counter
    //
    // Return      :
    ///              \return    void
    // 
    // ------------------------------------------------------------------------------
    ASLLoop(const std::string&              sIdent,
      Mediator*              pMediator,
      const Key<int32_t>&    rLoopLengthKey  = AslKeys::l_ASL_LOOP_LENGTH,
      const Key<int32_t>&    rLoopCounterKey = AslKeys::l_ASL_LOOP_COUNTER);

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::ASLLoop
    //
    // Description :
    ///              \brief     Constructor
    ///
    ///                         \c SL::AslKeys::::l_ASL_LOOP_LENGTH and \c SL::AslKeys::::l_ASL_LOOP_COUNTER 
    ///                         will be used as default keys for the loop length and
    ///                         loop counter respectively
    //
    // Parameters  :
    ///              \param     sScope              String identifier of the scope of this component
    ///              \param     sIdent              String identifier of the component
    ///              \param     pMediator           Pointer to Mediator class
    ///              \param     rLoopLengthKey      Key specifying the loop length
    ///              \param     rLoopCounterKey     Key specifying the loop counter
    //
    // Return      :
    ///              \return    void
    // 
    // ------------------------------------------------------------------------------
    ASLLoop(const std::string&              sScope,
      const std::string&                    sIdent,
      Mediator*              pMediator,
      const Key<int32_t>&    rLoopLengthKey  = AslKeys::l_ASL_LOOP_LENGTH,
      const Key<int32_t>&    rLoopCounterKey = AslKeys::l_ASL_LOOP_COUNTER);
    ~ASLLoop(void);

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::configure
    //
    // Description :
    ///              \brief     Configures the asl loop, i.e. registers
    ///                         the loop counter at the mediator
    //
    // Parameters  :
    ///              \param     rMrProt     Reference to the protocol class
    ///              \param     rSeqLim     Reference to the sequence limits class
    //
    // Return      :
    ///              \return    \arg \c     true  : successful execution \n
    ///                         \arg \c     false : error occurred
    //
    // ------------------------------------------------------------------------------
    virtual bool configure (const MrProt& rMrProt, const SeqLim& rSeqLim);

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::prep
    //
    // Description :
    ///              \brief     Prepares the asl loop
    ///
    ///                         Reads the loop length from the mediator using the 
    ///                         key specified in the constructor (\c rLoopLengthKey)
    //
    // Parameters  :
    ///              \param     rMrProt     Reference to the protocol class
    ///              \param     rSeqLim     Reference to the sequence limits class
    ///              \param     rSeqExpo    Reference to the sequence exports class
    //
    // Return      :
    ///              \return    \arg \c     true  : successful execution \n
    ///                         \arg \c     false : error occurred
    //
    // ------------------------------------------------------------------------------
    virtual bool prep (MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo);

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::run
    //
    // Description :
    ///              \brief     Executes the asl loop
    //
    // Parameters  :
    ///              \param     rMrProt     Reference to the protocol class
    ///              \param     rSeqLim     Reference to the sequence limits class
    ///              \param     rSeqExpo    Reference to the sequence exports class
    ///              \param     pSLC        Pointer to the rotation matrix and slice
    ///                                     position information
    //
    // Return      :
    ///              \return    \arg \c     true  : successful execution \n
    ///                         \arg \c     false : error occurred
    //
    // ------------------------------------------------------------------------------
    virtual bool run (MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, sSLICE_POS* pSLC);

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::aslRun
    //
    // Description :
    ///              \brief     Executes one asl state in the asl loop
    //
    // Parameters  :
    ///              \param     rMrProt     Reference to the protocol class
    ///              \param     rSeqLim     Reference to the sequence limits class
    ///              \param     rSeqExpo    Reference to the sequence exports class
    ///              \param     pSLC        Pointer to the rotation matrix and slice
    ///                                     position information
    ///              \param     aslState    Current aslState.
    //
    // Return      :
    ///              \return    \arg \c     true  : successful execution \n
    ///                         \arg \c     false : error occurred
    // 
    // ------------------------------------------------------------------------------
    virtual bool aslRun(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, sSLICE_POS* pSLC, long aslState );

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::updateASLParameters
    //
    // Description :
    ///              \brief     Updates mediator keys for given asl state.
    //
    // Parameters  :
    ///              \param     rMrProt     Reference to the protocol class
    ///              \param     rSeqLim     Reference to the sequence limits class
    ///              \param     rSeqExpo    Reference to the sequence exports class
    ///              \param     pSLC        Pointer to the rotation matrix and slice
    ///                                     position information
    ///              \param     aslState    Current aslState.
    //
    // Return      :
    ///              \return    \arg \c     true  : successful execution \n
    ///                         \arg \c     false : error occurred
    // 
    // ------------------------------------------------------------------------------
    bool updateASLParameters(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, sSLICE_POS* pSLC, long aslState );

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::updateASLState
    //
    // Description :
    ///              \brief     Updates the mediator key AslLeafAdapter::e_LABEL_STATE used in the asl specific leaf adapters. 
    ///                         The mode is determined by the scanContext, the mediator key SL::AslKeys::b_IS_M0_SCAN and the iteration
    ///                         of the ASLLoop.
    ///                         To support an aritrary number of labels a long is provided. The last iteration is used as control scan.
    ///                         Returns true if successfull.
    //
    // Parameters  :
    ///              \param     rMrProt     Reference to the protocol class
    ///              \param     aslState    Asl mode specified as long to support arbitrary number of labels. 
    //
    // Return      :
    ///              \return    Duration in us
    // 
    // ------------------------------------------------------------------------------
    bool updateASLState(MrProt& rMrProt, long aslState);

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::getDuration
    //
    // Description :
    ///              \brief     Returns the total duration of the preparing scans
    //
    // Parameters  :
    ///              \param     eRange  Specifies the range (ELEMENT or TREE)
    //
    // Return      :
    ///              \return    Duration in us
    // 
    // ------------------------------------------------------------------------------
    virtual int64_t getDuration (Range eRange);

    // ------------------------------------------------------------------------------
    //
    // Name        :            ASLLoop::getEnergy
    //
    // Description :
    ///              \brief     Returns the total energy of the preparing scan
    //
    // Parameters  :
    ///              \param     eRange      Specifies the range (ELEMENT or TREE)
    //
    // Return      :
    ///              \return    Energy in Ws
    // 
    // ------------------------------------------------------------------------------
    virtual MrProtocolData::SeqExpoRFInfo getRFInfo (Range eRange);

  private:
    ASLLoop (const ASLLoop &right);   // Copy constructor not implemented 
    ASLLoop & operator= (const ASLLoop &right);   // Assignment operator not implemented 
};

} // namespace SL_NAME_SPACE 


