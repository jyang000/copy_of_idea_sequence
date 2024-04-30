//-----------------------------------------------------------------------------
// <copyright file="a_ep2d_fid_UI.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 2010-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
// <description>This file contains the implementation of the class Ep2d_fid_UI</description>
//-----------------------------------------------------------------------------

#pragma once

#ifndef a_ep2d_fid_UI_h
#define a_ep2d_fid_UI_h

//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrImaging/seq/a_ep_CommonUI.h"
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"

//  -------------------------------------------------------------------------- *
//  Forward declarations                                                       *
//  -------------------------------------------------------------------------- *
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
//  Name        : Ep2d_fid_UI
//
//  Description :
/// \brief        This class basically is a storage for the pointers to the
///                original setValue / getValue / solve - handlers.
///
///               The sequence registers new UI handlers, which usually do
///                something, then call the original UI handler, and then
///                do something else. To keep the information of the original
///                UI handlers, the Ep2d_fid_UI class stores the pointers
///
///               It also provides the method registerUI to execute the
///                registration of all new handlers (and the storage of
///                 the original pointers)
///
//  --------------------------------------------------------------------------

class Ep2d_fid_UI : public EpCommonUI
{
  public:
    //  --------------------------------------------------------------
    //
    //  Name        :  Ep2d_fid_UI::Ep2d_fid_UI
    //
    //  Description :
    /// \brief         Initialization of class members
    //
    //  Return      :
    //
    //  --------------------------------------------------------------
    Ep2d_fid_UI() = default;

    //  --------------------------------------------------------------
    //
    //  Name        :  Ep2d_fid_UI::~Ep2d_fid_UI
    //
    //  Description :
    /// \brief         Destructor
    //
    //  Return      :
    //
    //  --------------------------------------------------------------
    virtual ~Ep2d_fid_UI() = default;

    //  --------------------------------------------------------------------------
    //
    //  Name        : Ep2d_fid_UI::registerUI
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
    UI_ELEMENT_BOOL m_FilterImage;

#endif // WIN32
};

} // end of namespace SEQ_NAMESPACE
#endif // of #ifndef a_ep2d_diff_UI_h
