//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProt\MDS\MDS.h
//	 Version:
//	  Author: CC_PROTSRV PRINCAH5 SCHABES6  LISTSTR5 GIRALUJP
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: MR::Measurement::CSequence::Prot::FastImaging
//
//	 Classes:
//
//	-----------------------------------------------------------------------------
#pragma once

#ifdef _MSC_VER
#pragma once
#endif

#include <cassert>
#include <vector>

#include "MrGlobalDefinitions/MrBasicTypes.h"
#include "MrProtSrv/Domain/CoreNative/MdsDefines.h"
#include "MrProtSrv/Domain/CoreNative/MrMds.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrVector.h"

#undef __IMP_EXP

#ifdef BUILD_MrProt
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"


#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4263 4264)
#endif


class  __IMP_EXP CMds: public MrProtocolData::MrMdsDataDelegate
{
    typedef MrProtocolData::MrMdsDataDelegate BasicImplementation;

public:
    CMds(const CMds& rSource);
    CMds& operator=(const CMds& rSource);
	explicit CMds(MrProtocolData::MrMdsData& rhs);

    virtual ~CMds();

    operator MrProtocolData::MrMdsData*() { return m_pData.get(); }
    operator const MrProtocolData::MrMdsData*() const { return m_pData.get(); }

    //=========================================================================================
    uint32_t mdsModeMask(void) const;
    uint32_t mdsModeMask(uint32_t _ulMdsModeMask);
    //=========================================================================================
    uint32_t mdsVariableResolution(void) const;
    uint32_t mdsVariableResolution(uint32_t _ulMdsVariableResolution);
    //=========================================================================================
    uint32_t mdsReconMode(void) const;
    uint32_t mdsReconMode(uint32_t _ulMdsReconMode);
    //=============================================================================

    double mdsRangeExtension(void) const;
    double mdsRangeExtension(double _dMdsRangeExtension);

    //=========================================================================================
    //=========================================================================================
    const MrProtocolData::VectorPatDbl& mdsEndPosSBCS(void) const;
    MrProtocolData::VectorPatDbl& mdsEndPosSBCS(void);

    const VectorPat<double> mdsEndPosSBCS(const MrProtocolData::VectorPatDbl& _vpat);
    const VectorPat<double> mdsEndPosSBCS(const VectorPat<double>& _vpat);

    //=========================================================================================
    //=========================================================================================
    const MrProtocolData::MrPTabMotionIntervalData& ptabMotionInterval(uint32_t nIndex) const;

    MrProtocolData::MrPTabMotionIntervalData& ptabMotionInterval(uint32_t nIndex);

    int32_t getalFree(uint32_t nIndex) const;

    double getadFree(uint32_t nIndex) const;

    int32_t setalFree(int32_t nVal, uint32_t nIndex);

    double setalFree(double dVal, uint32_t nIndex);

    bool isMDSOn(void) const;

private:
};


#ifdef _MSC_VER
#pragma warning(pop)
#endif


#undef __OWNER