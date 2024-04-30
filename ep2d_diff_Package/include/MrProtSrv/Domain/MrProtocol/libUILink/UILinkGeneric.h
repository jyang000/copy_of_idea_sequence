//  -----------------------------------------------------------------
//  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential
//  -----------------------------------------------------------------
//
//  Project: NUMARIS/4
//     File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProtocol\libUILink\UILinkGeneric.h
//  Version: \main\8
//   Author: SCHAAL7H
//     Date: 2013-02-14 18:24:16 +01:00
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

#ifndef __UILINKGENERIC_H
#define __UILINKGENERIC_H

#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkBase.h"
#include "MrProtSrv/Common/MrGenericDC/ExternalInterface.h"
#include "MrProtSrv/Domain/MrProtocol/UILink/MrStdNameTags.h"

//  -----------------------------------------------------------------
//  Import-Export-Control
//
#ifdef BUILD_libUILink
  #define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

class __IMP_EXP MrUILinkGeneric : public MrUILinkBase
{
public:    
    typedef MrUILinkGeneric MYTYPE;
	
    //  /////////////////////////////////////////////////////////////
    //  Definition of the interfaces to the various handler functions:
    //  /////////////////////////////////////////////////////////////

  	//  See MrUILinkBase::PFctGetLabelId for a detailed description.
    typedef unsigned (*PFctGetLabelId)(
        MYTYPE* const pThis,
        char* arg_list[],
        INDEX_TYPE lPos
		);

	//  See MrUILinkBase::PFctGetToolTipId for a detailed description.
    typedef unsigned (*PFctGetToolTipId)(
        MYTYPE* const pThis,
        char* arg_list[],
        INDEX_TYPE lPos
		);

    //  The interface to a get-value-handler is defined by the 
    //  PFctGetValue procedure prototype. The return value of the
    //  get-value-handler determines the numeric value displayed
    //  in the user interface.
    //  Whenever the get-value-handler is invoked it is passed two
    //  argumants:
    //  pThis
    //    A pointer to the UILink-object related to the parameter.
    //  lPos
    //    Receives the array index, if the parameter itself is
    //    an array element. Otherwise the value is ambiguous.
	typedef MrPtr<MrGenericDC::IValueNode> (*PFctGetValue)(
        MYTYPE* const pThis,
        INDEX_TYPE lPos
		);

    //  The interface to a set-value-handler is defined by the 
    //  PFctSetValue procedure prototype. The set-value-handler
    //  is invoked by the UI-software, whenever the user types
    //  a new value and within binary-search.
    //  Three arguments are passed to the set-value-handler:
    //  pThis
    //    A pointer to the UILink-object related to the parameter
    //  newVal
    //    Receives the desired new value. The set-value-handler is
    //    responsible for the correct update of the actual protocol.
    //  lPos
    //    Receives the array index, if the parameter itself is
    //    an array element. Otherwise the value is ambiguous.
    //  
    //  The set-value-handler must return the actual new value.
    //
	typedef MrPtr<MrGenericDC::IValueNode> (*PFctSetValue)(
		MYTYPE* const pThis,
		const MrPtr<MrGenericDC::IValueNode>& newVal,
        INDEX_TYPE lPos
		);
	
	
	//  See MrUILinkBase::PFctSolve for a detailed description.
    typedef unsigned (*PFctSolve)(
        MYTYPE* const,
        char* arg_list[],
        const void*,
        const MrProtocolData::MrProtData*,
        INDEX_TYPE
        );

    //  See MrUILinkBase::PFctTry for a detailed description.
    typedef bool (*PFctTry)(
        MYTYPE* const,
        void*,
        const MrProtocolData::MrProtData*,
        INDEX_TYPE
        );

    
private:

    class UILinkGenericData
    {
    public:
		PFctGetValue        m_getValue;
        PFctSetValue        m_setValue;
    	   
        UILinkGenericData()
        {
            m_getValue = NULL;
            m_setValue = NULL;
        }
    private:
        UILinkGenericData(const UILinkGenericData& /*rhs*/);
    };

    
    mutable UILinkGenericData m_GenericData;
public:
    
	//  The return value is false, unless a get-size and a get-max-size
    //  handler has been set. If both size-handler and an is-available-
    //  handler have been set MrUILinkBase::isAvailable is called.
    virtual bool isAvailable(int32_t pos) const;

	bool isEditable(INDEX_TYPE pos) const;
	
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
    //  its own get-value-handler-function. The return value is a 
    //  pointer to the previous get-value-handler function registered
    //  by registerGetValueHandler. If no previous function has been set,
    //  the return value is NULL.
    PFctGetValue registerGetValueHandler(PFctGetValue _pCallBackFct);

    //  The member function allows the sequence programmer to specify
    //  its own set-value-handler-function. The return value is a 
    //  pointer to the previous set-value-handler function registered
    //  by registerSetValueHandler. If no previous function has been set,
    //  the return value is NULL.
    PFctSetValue registerSetValueHandler(PFctSetValue _pCallBackFct);
		
	//  Invokes get-value-handler, if isAvailable returns true.
    //  Otherwise the behavior is undefined.
    //unsigned value(char* arg_list[] newValue, int32_t pos) const;
    MrPtr<MrGenericDC::IValueNode> value(MrPtr<MrGenericDC::IValueNode> newValue, int32_t pos);
    MrPtr<MrGenericDC::IValueNode> value(int32_t pos);

	//  The member function allows the sequence programmer to specify
    //  its own is-available-handler-function. The return value is a 
    //  pointer to the previous is-available-handler function registered
    //  by registerIsAvailableHandler. If no previous function has been set,
    //  the return value is NULL.
    MrUILinkBase::PFctIsAvailable registerIsAvailableHandler(MrUILinkBase::PFctIsAvailable pCallBackFct);
    //  The member function allows the sequence programmer to specify
    //  its own try-handler. The return value is a 
    //  pointer to the previous try-handler-function registered
    //  by registerClearMemoryHandler. If no previous function has been set,
    //  the return value is NULL.
    PFctTry registerTryHandler(PFctTry pCallBackFct);
    PFctTry registerTryHandler(MrUILinkBase::PFctTry pCallBackFct);

    //  The member function allows the sequence programmer to specify
    //  its own solve-handler-function. The return value is a 
    //  pointer to the previous solve-handler function registered
    //  by registerSolveHandler. If no previous function has been set,
    //  the return value is NULL.
    PFctSolve registerSolveHandler(PFctSolve pCallBackFct);
    PFctSolve registerSolveHandler(MrUILinkBase::PFctSolve pCallBackFct);

private:
	MrUILinkGeneric(const MrUILinkGeneric&) = delete;
	MrUILinkGeneric& operator=(const MrUILinkGeneric&) = delete;

public:
	explicit MrUILinkGeneric(const UILinkSharedData&);

    //  Destructor (only for internal use)
    ~MrUILinkGeneric();

	virtual void unregister();
};

// typedefs
typedef MrUILinkGeneric LINK_GENERIC_TYPE;

__IMP_EXP_TEMPLATE template class __IMP_EXP _search< MrUILinkGeneric >;
__IMP_EXP_TEMPLATE template class __IMP_EXP _searchElm< MrUILinkGeneric >;
__IMP_EXP_TEMPLATE template class __IMP_EXP _create< MrUILinkGeneric >;


#endif // __UILINKGENERIC_H
//  -----------------------------------------------------------------
//  Copyright (C) Siemens AG 2008  All Rights Reserved.  Confidential
//  -----------------------------------------------------------------
