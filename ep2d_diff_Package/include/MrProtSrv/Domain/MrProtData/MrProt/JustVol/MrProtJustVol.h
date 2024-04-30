//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProt\JustVol\MrProtJustVol.h
//	 Version: \main\17
//	  Author: CC_PROTSRV PRINCAH5 SCHABES6 LISTSTR5 GIRALUJP
//	    Date: 2012-10-17 15:21:47 +02:00
//
//	    Lang: C++
//
//	 Descrip: MR::Measurement::CSequence::Prot::JustVol
//
//	 Classes:
//
//	-----------------------------------------------------------------------------

#pragma once

#include "MrProtSrv/Domain/CoreNative/MrProtData.h"

#ifdef BUILD_MrProtJustVol
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"


struct MrProtJustVolData;

class  MrProtJustVol
{
public:
    __IMP_EXP MrProtJustVol();
    __IMP_EXP explicit MrProtJustVol(bool bRegardNavigators);

    __IMP_EXP virtual ~MrProtJustVol();

    //	calculate adjustment slice
    __IMP_EXP void calcAdjSlice (const MrProtocolData::MrProtData::Pointer pProt, MrProtocolData::MrSliceData &adjSlice, bool localShimEnabled = false);
    __IMP_EXP void calcAdjSlice(
        const double dDefaultFovSag,
        const double dDefaultFovCor,
        const double dDefaultFovTra,
        const double dDefaultPosSag,
        const double dDefaultPosCor,
        const double dDefaultPosTra,
        const MrProtocolData::MrProtData::Pointer pProt,
        MrProtocolData::MrSliceData &adjSlice,
        bool localShimEnabled);

    __IMP_EXP static void calcAdjMDSCoilPosSlice(MrProtocolData::MrSliceData &adjSlice);

    __IMP_EXP static void calcAdjMDSScoutSlice(MrProtocolData::MrSliceData &adjSlice);

    __IMP_EXP static void calcAdjMDSFreSlice(MrProtocolData::MrSliceData &adjSlice);

    //	calculate adjustment slice by fitting a axis aligned bounding box around all
    //	slices and the spectro voi (if present)
    __IMP_EXP void calcAdjSliceAsBBox (const MrProtocolData::MrProtData::Pointer pProt, MrProtocolData::MrSliceData &adjSlice);

    //  calculate bounding box around imaging FoV. If Imaging FoV consists of multiple slice groups, the bounding box is axis aligned
    __IMP_EXP void calcImagingFoVAsBBox (const MrProtocolData::MrProtData::Pointer pProt, MrProtocolData::MrSliceData &eImagingFoVBB);

protected:
    void accountNavigator(bool bAccount);
    void adaptTraSize(bool bAdapt);

protected:
    //	calculate min/max corners
    void calcCornerMinMax (const MrProtocolData::MrSliceData &slice);

    //	calculate adjustment slice by getting the (valid !) adj. volume from the
    //	protocoll (and checking/correcting if phaseFOV <= readoutFOV)
    void calcAdjSliceFromProt (const MrProtocolData::MrProtData::Pointer pProt, MrProtocolData::MrSliceData &adjSlice);

    //	calculate adjustment slice by fitting it to the first slice group (that has
    //	to be in the protocoll)
    void calcAdjSliceFromOneGroup (const MrProtocolData::MrProtData::Pointer pProt, MrProtocolData::MrSliceData &adjSlice);

    //	calculate adjustment slice by copying the spectro voi that has to be valid
    //	(checks/corrects if phaseFOV <= readoutFOV)
    void calcAdjSliceFromSpectro (const MrProtocolData::MrProtData::Pointer pProt, MrProtocolData::MrSliceData &adjSlice);

    //	calculate adjustment slice by inserting the default volume (70% of GC, i.e.
    //	shim-fieldmap-volume)
    void calcAdjSliceAsDefaultVolume (double dFovSag, double dFovCor, double dFovTra, double dPosSag, double dPosCor, double dPosTra, MrProtocolData::MrSliceData &adjSlice) const;

    //	shrink adjustemt slice for neck shim (25% in AP direction)
    void shrinkAdjVolume4NeckShim (MrProtocolData::MrSliceData &adjSlice);

    //	shrink adjustemt slice for local shim (110 mm in RL direction)
    void shrinkAdjVolume4LocalShim (MrProtocolData::MrSliceData &adjSlice);

    //	shrink adjustemt slice for whole body shim (280 mm in RL and 350 mm in AP direction)
    void shrinkAdjVolume4WholeBodyShim (MrProtocolData::MrSliceData &adjSlice);

    void calcCornerAbsMin(void);

private:
    MrProtJustVolData*    m_pData;

private:
    MrProtJustVol(const MrProtJustVol &right);
    const MrProtJustVol& operator=(const MrProtJustVol &right);
};