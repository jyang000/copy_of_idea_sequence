//  -----------------------------------------------------------------------------
//    Copyright (C) Siemens AG 2013  All Rights Reserved.
//  -----------------------------------------------------------------------------
//
//   Project: NUMARIS/X
//      File: \src\MrImaging\seq\a_tgse_asl\AslLeafAdapter.h
//    Author: pfeujodj
//      Date: 2018-08-07 12:23:01 +02:00
//
//      Lang: C++
//
//   Descrip: 3D ASL sequence with CompositeSeqLoop support
//
//  -----------------------------------------------------------------------------

#pragma once

#include "MrImaging/SequenceLibraries/libASL/SBBAsl.h"
#include "MrImagingFW/libCSL/LeafAdapter.h"
#include "MrImagingFW/libCSL/Counter.h"
#include "MrImagingFW/libCSL/Key.h"
#include "MrImaging/libSL/StdSL_ID.h"

class ReorderInfo;
class MrProt;
class SeqLim;

namespace MrProtocolData { class SeqExpo; }

// Import/export control
#ifdef BUILD_libSL
  #define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

#ifndef LIB_NAMESPACE
  #error LIB_NAMESPACE not defined
#endif
using namespace LIB_NAMESPACE;
namespace SL_NAME_SPACE 
{


class AslLeafAdapter : public LeafAdapter<SeqBuildBlockAsl>  
{
  public:
    // * ------------------------------------------------------------------ *
    // * Keys used in the AslLeafAdapter                                   *
    // * ------------------------------------------------------------------ *
    static const Key<SYMBOL_ASL::AslLabelState> e_LABEL_STATE;
    static const Key<long>                  l_CUR_TI_PHASE;
    static const Key<long>                  l_TIME_TO_IMAGE_SLICE_EXC_US;

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  AslLeafAdapter::AslLeafAdapter                            *
    // *                                                                    *
    // * Description :  Constructor                                         *
    // *                                                                    *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    AslLeafAdapter (const std::string& sIdent, Mediator* pMediator, SBBList* pSBBList);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  AslLeafAdapter::AslLeafAdapter                            *
    // *                                                                    *
    // * Description :  Constructor                                         *
    // *                                                                    *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    AslLeafAdapter (const std::string& sScope, const std::string& sIdent, Mediator* pMediator, SBBList* pSBBList);

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  AslLeafAdapter::~AslLeafAdapter                           *
    // *                                                                    *
    // * Description :  Destructor                                          *
    // *                                                                    *
    // * Return      :  void                                                *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    virtual ~AslLeafAdapter ();

    // * ------------------------------------------------------------------ *
    // *                                                                    *
    // * Name        :  AslLeafAdapter::prep                                   *
    // *                                                                    *
    // * Description :  Calculates the timing of the TRUFI kernel and       *
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
    // * Name        :  AslLeafAdapter::run                                    *
    // *                                                                    *
    // * Description :  Executes the TRUFI kernel                           *
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
    // * Name        :  AslLeafAdapter::setdPhaseIncr                          *
    // *                                                                    *
    // * Description :  Increments the phase of the Medic kernel. May be    *
    // *                used to implement phase cycles, rf spoiling, ...    *
    // *                                                                    *
    // * Parameter   :  void                                                *
    // *                                                                    *
    // * Return      :  double                                              *
    // *                                                                    *
    // * ------------------------------------------------------------------ *
    virtual void setdPhaseIncr (double dPhaseIncr);

    // ------------------------------------------------------------------------------
    //
    // Name        :            AslLeafAdapter::getEnergy
    //
    // Description :
    ///              \brief     Returns the total energy of the Medic kernel
    //
    // Parameters  :
    ///              \param     Range eRange : Specifies the range (ELEMENT or TREE)
    //
    // Return      :
    ///              \return    Energy in Ws
    // 
    // ------------------------------------------------------------------------------
    virtual MrProtocolData::SeqExpoRFInfo getRFInfo (Range eRange);

    int64_t getDuration (Range eRange);

  private:
    // * ------------------------------------------------------------------ *
    // * Copy constructor not implemented                                   *
    // * ------------------------------------------------------------------ *
    AslLeafAdapter (const AslLeafAdapter& right);

    // * ------------------------------------------------------------------ *
    // * Assignment operator not implemented                                *
    // * ------------------------------------------------------------------ *
    AslLeafAdapter& operator= (const AslLeafAdapter &right);
};

} // namespace SL_NAME_SPACE

