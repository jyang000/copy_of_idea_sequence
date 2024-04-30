//-----------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2021. All Rights Reserved.
//-----------------------------------------------------------------------------

#pragma once

#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrMeasSrv/SeqIF/libRT/sSLICE_POS.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProtSliceSeries.h"
#include "MrProtSrv/Domain/CoreNative/MrSliceGroupData.h"
#include "MrProtSrv/Domain/CoreNative/MrPat.h"
#include "MrProtSrv/Domain/CoreNative/Slice.h"
#include "MrProtSrv/Domain/CoreNative/MrTXSpec.h"

#ifdef BUILD_libSeqUtil
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"



namespace SliceAccelerationUtils
{
struct __IMP_EXP SliceAccelerationParams
{
    double      dSliceThickness;
    double      dDistFactor;
    long        lFOVShiftFactor;
    long        lMultiBandFactor;
    long        lNumberOfSlices;
    const char* pcNucleus;

    void setSliceAccelerationParamsFromProtocol(MrProt& rMrProt)
    {
        dSliceThickness  = rMrProt.sliceSeries().front().getdThickness();
        dDistFactor      = rMrProt.getsGroupArray().getasGroup()[0].getdDistFact();
        lFOVShiftFactor  = rMrProt.getsSliceAcceleration().getlFOVShiftFactor();
        lMultiBandFactor = rMrProt.getsSliceAcceleration().getlMultiBandFactor();
        lNumberOfSlices  = rMrProt.getsSliceArray().getlSize();
        pcNucleus        = rMrProt.getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str();
    }
};

    // Note: The delta moment dM will introduce a phase shift of 2 pi between the
    //       outermost slices (similar to 3D encoding). For slice acceleration,
    //          D = ( Gap + Thickness) * ( N / MB ) * F
    //            = ( Gap + Thickness) * N * ( 1 + OS ),
    //       where N denotes the  total number of slices.
    //       With this (and the gyromagnetic ratio g), we get
    //         dM = (2 pi/g) / D
    //            = (2 pi/g) / ( ( Gap + Thickness) * ( N / MB ) * F )

    __IMP_EXP double calcSliceBandDistance(SliceAccelerationParams const& sliceAccelerationParams);

    __IMP_EXP double calcSliceAccelerationDeltaMoment(MrProt& rMrProt);

    // PhaseOffCenter3D is currently calculated according to
    // (see also the comments on fGSLGetSliceAccelerationDeltaMoment):
    //      360.0 * SliceShift / ( Thickness * ( 1 + OverSampling) )
    // For slice-acceleration, this translates to:
    //      360.0 * SliceShift / ( ( Thickness + Gap ) * NumberOfSlices / ( MultiBandFactor * FOVShiftFactor
    //      ) )
    //    = 360.0 * SliceShift / ( SliceBlockThickness * ( 1 + OverSampling ) )
    // where SliceBlockThickness = NumberOfSlices  * ( Thickness + Gap )
    // and   FOVShiftFactor      = MultiBandFactor * ( 1 + OverSampling).
    //
    // Note: Check whether adding a correction phase for even number of slices
    //       would be desirable here (compare sSLICE_POSImpl.cpp)
    __IMP_EXP double calcSliceAccelerationPhaseOffcenter3D(MrProt& rMrProt, sSLICE_POS const& sSLC);

    __IMP_EXP long performModuloWithFovShift(long lFOVShiftFactor, long line);

    __IMP_EXP long calcSlicePartitionNumber(long lFOVShiftFactor, long lLinNo);

} // namespace SliceAccelerationUtils
