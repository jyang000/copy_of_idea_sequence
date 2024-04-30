//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2019  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4\pkg\MrServers\MrImaging\seq\a_NetworkUI.h
//	 Version:
//	  Author: NEUR
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: Common UI for networks
//            - Provides a container for available network configurations
//              (considering PAT-factor- and segmentation-specific networks)
//            - Provides UI-handlers to select corresponding reconstruction
//              options (and to resolve parameter conflicts)
//
//	 Classes:
//
//	-----------------------------------------------------------------------------

#pragma once

//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrGlobalDefinitions/MrResult.h"                   // NLS_STATUS

#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"   // UI_ELEMENT_...
#endif


#if defined BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

#include <set>
#include <string>
#include <vector>

// Forward declarations
namespace MrProtocolData
{
    class SeqExpo;
}
class  SeqLim;
class  Sequence;
class  MrProt;

namespace SEQ_NAMESPACE
{
    //  --------------------------------------------------------------------------
    //
    //  Name        : NetworkUI
    //
    //  Description :
    /// \brief        This class basically is a storage for the pointers to the
    ///                original setValue / getValue / solve - handlers.
    ///         
    ///               All common handlers which are relevant in the context of 
    ///               network-enabled reconstructions are considered here.
    //  --------------------------------------------------------------------------

    class NetworkUI
    {
    public:

#ifdef WIN32

        //  ------------------------------------------------------------------
        /// Resolve dependencies between reconstruction mode and PAT parameters 
        //  (to be used by PAT mode setValue-handler)
        //  ------------------------------------------------------------------
        void setValuePATModeDependencies(
            MrUILinkBase*       pThis,
            unsigned            uNewPATMode, 
            int32_t             iOrigNumberOfRefLines
        );
#endif

        NetworkUI() = default;

        virtual ~NetworkUI() = default;



        //  --------------------------------------------------------------------------
        //  Name        : NetworkUI::registerUI
        //  Description : Register UI methods
        /// \brief        This function initializes the UI functions and
        ///               registers all given set / get / Solve - handlers
        //  --------------------------------------------------------------------------
        virtual NLS_STATUS registerUI (SeqLim &rSeqLim );

#ifdef WIN32
        //  --------------------------------------------------------------
        ///  \brief Helper class instances for UI handlers
        ///         - register new handler functions
        ///         - save pointer to original handler function
        ///         These classes exist only on the host.
        //  --------------------------------------------------------------
        UI_ELEMENT_SELECTION  m_ChannelCompression;
        UI_ELEMENT_SELECTION  m_CoilCombineMode;
        UI_ELEMENT_SELECTION  m_POCS;

        UI_ELEMENT_LONG       m_PATReferenceLines;

        UI_ELEMENT_SELECTION m_AdvancedReconMode;
        UI_ELEMENT_SELECTION m_DenoisingMethod;
#endif

    private:
    };

#ifdef WIN32

    namespace NetworkUINS
    {
        //  --------------------------------------------------------------
        //   Overloaded UI handlers
        //  --------------------------------------------------------------
        unsigned fCoilCombineModeSetValue         ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        bool     fCoilCombineModeTry              ( LINK_SELECTION_TYPE* const pThis, void* pClientMem, const MrProtocolData::MrProtData* pOrigProt,                      int32_t lIndex );
        unsigned fPOCSSetValue                    ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        bool     fPATReferenceLinesIsAvailable    ( LINK_LONG_TYPE* const pThis,                                                                                          int32_t lIndex );
        unsigned fAdvancedReconModeSetValue       ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fDenoisingMethodSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewValue, int32_t lIndex);

    } // end of namespace NetworkUINS
#endif // #ifdef WIN32

} // end of namespace SEQ_NAMESPACE
