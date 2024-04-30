//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProt\Application\Application.h
//	 Version:
//	  Author: CC_PROTSRV KAIMRAJR SCHABES6 MAIECOD6 LISTSTR5
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: MR::Measurement::CSequence::Prot::Application
//
//	 Classes:
//
//	-----------------------------------------------------------------------------
#pragma once

#include "MrProtSrv/Domain/CoreNative/MrApplication.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrVector.h"

class Slice;

#ifdef BUILD_MrProt
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"



class __IMP_EXP Flow: public MrProtocolData::MrFlowDataDelegate
{
    typedef MrProtocolData::MrFlowDataDelegate BasicImplementation;

public:
    Flow();
    Flow(const Flow& rSource);
    Flow& operator=(const Flow& rSource);
    explicit Flow(MrProtocolData::MrFlowData& data);
    explicit Flow(MrProtocolData::MrFlowData::Pointer data);
    virtual ~Flow(void);

    operator MrProtocolData::MrFlowData*() { return m_pData.get(); }
    operator const MrProtocolData::MrFlowData*() const { return m_pData.get(); }

    //	Return velocity [cm/s].
    int16_t velocity (void ) const;
    //	Sets velocity to value [cm/s].
    int16_t velocity (int16_t value);

    //	Return enum FlowDir, which describes direction.
    SEQ::FlowDir dir (void ) const;
    //	Sets flow direction.
    SEQ::FlowDir dir (SEQ::FlowDir value);

    //	Return const reference to vector, which specifies the flow direction in the
    //	case dir() == SEQ::FLOW_DIR_FREE.
    const VectorPat<double>& freeDir (void ) const;
    //	Return reference to vector, which specifies the flow direction in the case
    //	dir() == SEQ::FLOW_DIR_FREE.
    VectorPat<double>& freeDir ();

    //	The member function returns an enum SEQ::FlowDirDisplay that describes the
    //	image text. On error the return value is SEQ::FLOW_DIR_INVALID;
    SEQ::FlowDirDisplay dirDisplay(const Slice& rSlice) const;

private:
    VectorPat<double>           m_FreeDir;
};


#undef __IMP_EXP
#undef __IMP_EXP_DATA