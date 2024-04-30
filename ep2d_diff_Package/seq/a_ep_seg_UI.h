//----------------------------------------------------------------------------------
// <copyright file="a_ep_seg_UI.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 1998-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015-2021. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrImaging/seq/a_ep_CommonUI.h"
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"

#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/UILink/MrStdNameTags.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkSelection.h"
#endif

//  -------------------------------------------------------------------------- *
//  Forward declarations                                                       *
//  -------------------------------------------------------------------------- *
class MrProt;
class SeqLim;
class SeqExpo;

namespace MrMeasSrv
{
class ISequence;
}

namespace SEQ_NAMESPACE
{
//  --------------------------------------------------------------------------
//
//  Name        : Ep_segUI
//
//  Description :
/// \brief        This class basically is a storage for the pointers to the
///                original setValue / getValue / solve - handlers.
///
///               The sequence registers new UI handlers, which usually do
///                something, then call the original UI handler, and then
///                do something else. To keep the information of the original
///                UI handlers, the TemplUI class stores the pointers
///
///               It also provides the method registerUI to execute the
///                registration of all new handlers (and the storage of
///                 the original pointers)
///
//  --------------------------------------------------------------------------

class Ep_segUI : public EpCommonUI
{
  public:
    Ep_segUI()          = default;
    virtual ~Ep_segUI() = default;

    //  --------------------------------------------------------------------------
    //
    //  Name        : Ep2dUI::registerUI
    //
    //  Description :
    /// \brief        This function initializes the UI functions and
    ///                registers all given set / get / Solve - handlers
    ///
    ///               It can be executed on the measurement system, too, but is empty there.
    ///
    ///               On the host, it executes these steps
    ///               - Declaration of pointers to UI classes
    ///               - Registration of overloaded set value handlers
    ///
    ///               It returns an NLS status
    ///
    virtual NLS_STATUS registerUI(SeqLim& pSeqLim);

#ifdef WIN32

    //  --------------------------------------------------------------
    ///  \brief Helper class instances for UI handlers
    ///         - register new handler functions
    ///         - save pointer to original handler function
    ///         These classes exit only on the host.
    ///
    ///  The following line is an example which can be removed for
    ///  other sequences.
    //  --------------------------------------------------------------
    UI_ELEMENT_LONG      m_Averages;
    UI_ELEMENT_LONG      m_Measurements;
    UI_ELEMENT_SELECTION m_RespComp;

#ifdef EP_SEG_FID
    UISelectionElement<LINK_BOOL_TYPE> m_SWI;
    // For flow comp mode
    UISelectionElement<LINK_SELECTION_TYPE> m_FlowCompMode;
#endif // SUPPORT_SWI

#endif
};

namespace Ep_segUINS
{
#ifdef WIN32
//  ----------------------------------------------------------------------
//
//  Name        :  getSeq
//
//  Description :
/// \brief         Returns the pointer to the sequence Ep2d
//
//  Return      :  Ep2d*
//
//  ----------------------------------------------------------------------
SEQ_NAMESPACE::Ep_seg* getSeq(MrUILinkBase* const pThis);

#ifdef SUPPORT_PACE
// ===========================================================================
// Handlers needed for PACE
long     fAvgGetValue(LINK_LONG_TYPE* const, long);
long     fAvgSetValue(LINK_LONG_TYPE* const, long, long);
unsigned fRespCompSetValue(LINK_SELECTION_TYPE* const, unsigned, long);
#endif
#endif // #ifdef WIN32

#ifdef NEED_a_ep_seg_calculateTRTIFillTimes
bool pcalculateTRTIFillTimes(
    MrProt& pMrProt, SeqLim& pSeqLim, SeqExpo& pSeqExpo, MrMeasSrv::ISequence* pSeq, long* plNeededTI, long* plNeededTR);
#endif

} // namespace Ep_segUINS
}; // namespace SEQ_NAMESPACE
