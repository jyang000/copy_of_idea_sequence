//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_se\a_ep2d_se_UI.h
//     Version: \main\3
//      Author: Clinical
//        Date: 2011-09-07 11:36:50 +02:00
//
//        Lang: C++
//
//
//
///  \file   a_ep2d_se_UI.h
///  \brief  File containing declaraion of the UI class
///         - ep2d_diff
///
///  This file contains the implementation of the class Ep2d_se_UI.
///
//    -----------------------------------------------------------------------------



#ifndef a_ep2d_se_UI_h
#define a_ep2d_se_UI_h 1



//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrImaging/seq/a_ep_CommonUI.h"
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"


#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkSelection.h"
#endif

//  -------------------------------------------------------------------------- *
//  Defines and typedefs                                                       *
//  -------------------------------------------------------------------------- *


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
//  Name        : Ep2d_se_UI
//
//  Description :
/// \brief        This class basically is a storage for the pointers to the
///                original setValue / getValue / solve - handlers.
///
///               The sequence registers new UI handlers, which usually do
///                something, then call the original UI handler, and then
///                do something else. To keep the information of the original
///                UI handlers, the Ep2d_se_UI class stores the pointers
///
///               It also provides the method registerUI to execute the
///                registration of all new handlers (and the storage of
///                 the original pointers)
///
//  --------------------------------------------------------------------------

class Ep2d_se_UI: public EpCommonUI
{

public:

    //  --------------------------------------------------------------
    //
    //  Name        :  Ep2d_se_UI::Ep2d_se_UI
    //
    //  Description :
    /// \brief         Initialization of class members
    //
    //  Return      :
    //
    //  --------------------------------------------------------------
    Ep2d_se_UI();


    //  --------------------------------------------------------------
    //
    //  Name        :  Ep2d_se_UI::~Ep2d_se_UI
    //
    //  Description :
    /// \brief         Destructor
    //
    //  Return      :
    //
    //  --------------------------------------------------------------
    virtual ~Ep2d_se_UI();

    //  --------------------------------------------------------------------------
    //
    //  Name        : Ep2d_se_UI::registerUI
    //
    //  Description :
    /// \brief        This function initializes the UI functions and
    ///                registers all given set / get / Solve - handlers
    ///
    ///               It can be executed on the measuement system, too, but is empty there.
    ///
    ///               On the host, it executes these steps
    ///               - Declaration of pointers to UI classes
    ///               - Registration of overloaded set value handlers
    ///
    ///               It returns an NLS status
    ///
    virtual NLS_STATUS registerUI (SeqLim &rSeqLim);


protected:

#ifdef WIN32

    //  --------------------------------------------------------------
    ///  \brief Helper class instances for UI handlers
    ///         - register new handler functions
    ///         - save pointer to original handler function
    ///         These classes exist only on the host.
    //  --------------------------------------------------------------
    UI_ELEMENT_SELECTION m_RFPulseType;
    UI_ELEMENT_SELECTION m_FatSatMode;

#ifdef EP2D_SE_MRE
    // Need access from overloaded handlers
public:
    UI_ELEMENT_DOUBLE m_TR;
    UI_ELEMENT_DOUBLE m_TE;
    UI_ELEMENT_LONG m_SGSize;
#endif

#endif // #ifdef WIN32

};

namespace Ep2d_se_UINS
{


} // end of namespace Ep2d_se_UINS
} // end of namespace SEQ_NAMESPACE
#endif  // of #ifndef a_ep2d_se_UI_h

