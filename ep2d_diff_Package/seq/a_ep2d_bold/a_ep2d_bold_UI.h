//----------------------------------------------------------------------------------
// <copyright file="a_ep2d_bold_UI.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 2010-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#ifndef a_ep2d_bold_UI_h
#define a_ep2d_bold_UI_h

#include "MrImaging/seq/a_ep_CommonUI.h"
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"

#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#endif

namespace MrProtocolData
{
class MrProtData;
}
class SeqLim;
class SeqExpo;
class Sequence;

namespace SEQ_NAMESPACE
{
//  --------------------------------------------------------------------------
//
//  Name        : Ep2d_bold_UI
//
//  Description :
/// \brief        This class basically is a storage for the pointers to the
///               original setValue / getValue / solve - handlers.
///
///               The sequence registers new UI handlers, which usually do
///               something, then call the original UI handler, and then
///               do something else. To keep the information of the original
///               UI handlers, the Ep2d_bold_UI class stores the pointers
///
///               It also provides the method registerUI to execute the
///               registration of all new handlers (and the storage of
///               the original pointers)
///
//  --------------------------------------------------------------------------

class Ep2d_bold_UI : public EpCommonUI
{
  public:
    Ep2d_bold_UI()          = default;
    virtual ~Ep2d_bold_UI() = default;

    // no copy and move constructor and assignment operator
    Ep2d_bold_UI(const Ep2d_bold_UI& right) = delete;
    Ep2d_bold_UI& operator=(const Ep2d_bold_UI& right) = delete;
    Ep2d_bold_UI(Ep2d_bold_UI&& right)                 = delete;
    Ep2d_bold_UI& operator=(Ep2d_bold_UI&& right) = delete;

    //  --------------------------------------------------------------------------
    //
    //  Name        : Ep2d_bold_UI::registerUI
    //
    //  Description :
    /// \brief        This function initializes the UI functions and
    ///               registers all given set / get / Solve - handlers
    ///
    ///               It can be executed on the measurement system, too, but is empty there.
    ///
    ///               On the host, it executes these steps
    ///               - Declaration of pointers to UI classes
    ///               - Registration of overloaded set value handlers
    ///
    ///               It returns an NLS status
    ///
    NLS_STATUS registerUI(SeqLim& rSeqLim) override;

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
    UI_ELEMENT_BOOL m_FilterNorm;
    UI_ELEMENT_BOOL m_FilterNormPreScan;
    UI_ELEMENT_BOOL m_FilterNormBific;
    UI_ELEMENT_BOOL m_FilterDistCorr;
    // UI_ELEMENT_BOOL      m_FilterImage;

    UI_ELEMENT_SELECTION m_RFPulseType;
#endif // WIN32
};

} // namespace SEQ_NAMESPACE
#endif // #ifndef a_ep2d_bold_UI_h
