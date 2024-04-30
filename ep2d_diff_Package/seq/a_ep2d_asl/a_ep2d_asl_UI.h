//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_asl\a_ep2d_asl_UI.h
//     Version: \main\4
//      Author: PLM AW NERUO
//        Date: 2011-09-07 11:36:34 +02:00
//
//        Lang: C++
//
//
//
///  \file   a_ep2d_asl_UI.h
///  \brief  File containing declaraion of the UI class
///         - ep2d_asl
///
///  This file contains the implementation of the class ep2d_ASL_UI.
///
//    -----------------------------------------------------------------------------



#ifndef a_ep2d_asl_UI_h
#define a_ep2d_asl_UI_h 1



//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrImaging/seq/a_ep_CommonUI.h"
#include "MrImaging/seq/a_ep2d.h"
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
//  Name        : ep2d_ASL_UI
//
//  Description :
/// \brief        This class basically is a storage for the pointers to the
///                original setValue / getValue / solve - handlers.
///
///               The sequence registers new UI handlers, which usually do
///                something, then call the original UI handler, and then
///                do something else. To keep the information of the original
///                UI handlers, the ep2d_ASL_UI class stores the pointers
///
///               It also provides the method registerUI to execute the
///                registration of all new handlers (and the storage of
///                 the original pointers)
///
//  --------------------------------------------------------------------------

class ep2d_ASL_UI: public EpCommonUI
{

public:

    //  --------------------------------------------------------------
    //
    //  Name        :  ep2d_ASL_UI::ep2d_ASL_UI
    //
    //  Description :
    /// \brief         Initialization of class members
    //
    //  Return      :
    //
    //  --------------------------------------------------------------
    ep2d_ASL_UI();


    //  --------------------------------------------------------------
    //
    //  Name        :  ep2d_ASL_UI::~ep2d_ASL_UI
    //
    //  Description :
    /// \brief         Destructor
    //
    //  Return      :
    //
    //  --------------------------------------------------------------
    virtual ~ep2d_ASL_UI();

    //  --------------------------------------------------------------------------
    //
    //  Name        : ep2d_ASL_UI::registerUI
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

#ifdef WIN32

	UI_ELEMENT_SELECTION m_FatSatMode;
	UI_ELEMENT_SELECTION m_AslMode;
    UI_ELEMENT_BOOL      m_SaveUncombined;
	UI_ELEMENT_BOOL      m_FilterNorm;
	UI_ELEMENT_BOOL      m_FilterPrescan;
	UI_ELEMENT_BOOL      m_FilterDiscor;
	UI_ELEMENT_BOOL      m_FilterImage;
	UI_ELEMENT_BOOL      m_FilterNormBific;

#endif

protected:

};

namespace ep2d_ASL_UINS
{

#ifdef WIN32

// ------------------------------------------------------------------------------
// Function    : _bandwidthGetLimits
// ------------------------------------------------------------------------------
//
// Description : calls original getLimits-handler for bandwidth, but sets search
//               mode to VERIFY_SCAN_ALL
// Return      : whatever original getLimits-handler says
//
// ------------------------------------------------------------------------------
bool _bandwidthGetLimits (LINK_DOUBLE_TYPE * const pThis, std::vector< MrLimitDouble >& rLimitVector, uint32_t& rulVerify, int32_t lIndex);

// ------------------------------------------------------------------------------
// Function    : _UnavailableOption
// ------------------------------------------------------------------------------
//
// Overloaded is-available handler for generic removal of UI options
//
// ------------------------------------------------------------------------------
bool _UnavailableOption(LINK_BOOL_TYPE* const /*pThis*/, int32_t /*pos*/);

// ------------------------------------------------------------------------------
// Function    : AslMode_GetToolTip
// ------------------------------------------------------------------------------
//
// remove irrevalent tooltip
//
// ------------------------------------------------------------------------------
unsigned _AslMode_GetToolTip(LINK_SELECTION_TYPE* const /*pThis*/, char * /*arg_list*/[], int32_t /*pos*/);

bool fPostLabelingDelayArraySizeIsVisible   (LINK_LONG_TYPE* const /*pThis*/, int32_t /*pos*/);
bool fInversionTimeArraySizeIsVisible       (LINK_LONG_TYPE* const /*pThis*/, int32_t /*pos*/);

#endif   // of #ifdef WIN32

} // end of namespace ep2d_ASL_UINS

} // end of namespace SEQ_NAMESPACE

#endif  //of #ifndef a_ep2d_ASL_UI_h