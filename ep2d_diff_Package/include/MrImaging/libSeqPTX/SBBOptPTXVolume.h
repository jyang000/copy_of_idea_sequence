//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 2013  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4\pkg\MrServers\MrImaging\libSeqPTX\SBBOptPTXVolume.h
//     Version: \main\1
//      Author: pfeujodj
//        Date: 2013-02-11 09:18:27 +01:00
//
//        Lang: C++
//
//     Descrip: optimized rf pulse for PTXVolume
//
//     Classes:
//
//    -----------------------------------------------------------------------------


#ifndef SBBOptPTXVolume_h
#define SBBOptPTXVolume_h 1


#include "MrImagingFW/libSBBFW/SeqBuildBlock.h"


//-----------------------------------------------------------------------------
// import/export control
//-----------------------------------------------------------------------------
#ifdef __IMP_EXP
#undef __IMP_EXP
#endif 

#ifdef LOCAL_PTX_BUILD
#define __IMP_EXP
#pragma message( "Local PTX build" )
#else  //   LOCAL_PTX_BUILD not defined
#ifdef BUILD_libSeqPTX
  #define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control
#endif  //  LOCAL_PTX_BUILD

class SeqBuildBlockOSat_Internal;

class __IMP_EXP SBBOptPTXVolume : public SeqBuildBlock
{
public:
    SBBOptPTXVolume (SBBList* pSBBList = nullptr);
    virtual ~SBBOptPTXVolume();

	SBBOptPTXVolume(const SBBOptPTXVolume &right) = delete;
	const SBBOptPTXVolume & operator=(const SBBOptPTXVolume &right) = delete;

	SBBOptPTXVolume(SBBOptPTXVolume&& right) = delete;
	const SBBOptPTXVolume & operator=(SBBOptPTXVolume&& right) = delete;

    ///    Prepares the SBB-functions for Run. Returns true, if successful.
    bool prep (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;

    ///    Executes the real time part of the SBB. Returns true, if successful.
    bool run (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC) override;

    // both vectors must have same size: zero size deactivates OptPTXVolume
    void init(std::vector<double> vdThickness, std::vector<double> vdShift);

private:

	SeqBuildBlockOSat_Internal * m_OSat;
    std::vector<double> m_vdOptPTXVolumeThickness{};
    std::vector<double> m_vdOptPTXVolumeShift{};
	size_t              m_uiOptPTXVolumeCount{ 0 };
};


#endif
