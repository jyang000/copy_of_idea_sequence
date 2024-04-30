//  -----------------------------------------------------------------
//  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential
//  -----------------------------------------------------------------
//
//  Project: NUMARIS/4
//     File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProtocol\libUILink\UILinkMap.h
//  Version: \main\6
//   Author: KOELLNER
//     Date: 2011-01-25 08:59:42 +01:00
//
//     Lang: C++
//
//  Descrip: 
//
//  Classes: 
//
//  -----------------------------------------------------------------

//  -----------------------------------------------------------------
//  Used interfaces
//
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkBase.h"

//  -----------------------------------------------------------------
//  Additional Declarations
//
#ifndef __UILINKMAP_H
#define __UILINKMAP_H

//  -----------------------------------------------------------------
//  Import-Export-Control
//
#ifdef BUILD_libUILink
  #define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

class IProtocolConfigRepository;
class ICoilHardware;
class IPeripheryInterface;
class IAdjustmentInterface;
class ISeriesBlockManagerInterface;
namespace MrMeasSrv
{
	class ISequence;
}
//  -----------------------------------------------------------------
//  class MrUILinkMap
//
class __IMP_EXP MrUILinkMap : public MrUILinkBase
{
public:
    typedef MrUILinkMap MYTYPE;

    //  /////////////////////////////////////////////////////////////
    //  Definition of the interfaces to the various handler functions:
    //  /////////////////////////////////////////////////////////////
	
    //  See MrUILinkBase::PFctGetLabelId for a detailed description.
    typedef unsigned (*PFctGetLabelId)(
        MrUILinkMap* const pThis,
        char** arg_list,
        int32_t lPos
		);

    //  See MrUILinkBase::PFctGetToolTipId for a detailed description.
    typedef unsigned (*PFctGetToolTipId)(
        MrUILinkMap* const,
        char**,
        int32_t
		);

    //  See MrUILinkBase::PFctIsAvailable for a detailed description.
	typedef bool (*PFctIsAvailable)(
		    MrUILinkMap* const pThis,
            int32_t lPos
		);

    //  /////////////////////////////////////////////////////////////
    //  Member functions
    //  /////////////////////////////////////////////////////////////

private:
	MrUILinkMap(const MrUILinkMap&) = delete;
	MrUILinkMap& operator=(const MrUILinkMap&) = delete;

public:
	explicit MrUILinkMap(const UILinkSharedData&);

    //  Destructor
//    ~MrUILinkMap();

    //  The member function allows the sequence programmer to specify
    //  its own get-label-handler-function. The return value is a 
    //  pointer to the previous get-label-handler function registered
    //  by registerGetLabelIdHandler. If no previous function has been set,
    //  the return value is NULL.
    PFctGetLabelId registerGetLabelIdHandler(PFctGetLabelId pCallBackFct);
    PFctGetLabelId registerGetLabelIdHandler(MrUILinkBase::PFctGetLabelId pCallBackFct);

    //  The member function allows the sequence programmer to specify
    //  its own get-tool-tip-handler-function. The return value is a 
    //  pointer to the previous get-tool-tip-handler function registered
    //  by registerGetToolTipIdHandler. If no previous function has been set,
    //  the return value is NULL.
    PFctGetToolTipId registerGetToolTipIdHandler(PFctGetToolTipId pCallBackFct);
    PFctGetToolTipId registerGetToolTipIdHandler(MrUILinkBase::PFctGetToolTipId pCallBackFct);

    //  The member function allows the sequence programmer to specify
    //  its own is-available-handler-function. The return value is a 
    //  pointer to the previous is-available-handler function registered
    //  by registerIsAvailableHandler. If no previous function has been set,
    //  the return value is NULL.
    PFctIsAvailable registerIsAvailableHandler(PFctIsAvailable pCallBackFct);
    PFctIsAvailable registerIsAvailableHandler(MrUILinkBase::PFctIsAvailable pCallBackFct);
};

typedef MrUILinkMap LINK_MAP_TYPE;

template<>
class _search<MrUILinkMap> : public _search_base
  {
    public:
      _search(SeqLim* pSeqLim,const char* pszKey) 
        : _search_base(pSeqLim,pszKey)
      {}

	  _search(MrUILinkBase* pBase, const char* pszKey, int32_t pos = 0, MrUILinkBase::SEARCH_MODE nSearchMode = MrUILinkBase::SEARCH_MODE::SEARCH_DEFAULT)
        : _search_base(pBase,pszKey,pos, nSearchMode)
      {}

      operator MrUILinkMap*()
      {
         if (! dataValidForSearch() )
           return 0;

         MrUILinkMap* pParam = dynamic_cast<MrUILinkMap*>( execute_search() );

         if ( pParam == 0 )
           return 0;

         if ( m_Data.m_pSCData == 0)
           return 0;

         switch (m_Data.m_nSearchMode)
         {
			case MrUILinkBase::SEARCH_MODE::SEARCH_AVAILABLE:
             return pParam->isAvailable(m_Data.m_pos) ? pParam : 0;

			case MrUILinkBase::SEARCH_MODE::SEARCH_EDITABLE:
             return 0;
         }
         return pParam;
      }
  };

__IMP_EXP_TEMPLATE template class __IMP_EXP _searchElm< MrUILinkMap >;
__IMP_EXP_TEMPLATE template class __IMP_EXP _create< MrUILinkMap >;
__IMP_EXP_TEMPLATE template class __IMP_EXP _createArray< MrUILinkMap >;

#endif // __UILINKMAP_H

//  -----------------------------------------------------------------
//  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential
//  -----------------------------------------------------------------
