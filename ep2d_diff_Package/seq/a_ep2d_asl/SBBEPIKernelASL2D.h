//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/X
//      Author: Uvo Hoelscher, SYS APPL
//
//        Lang: C++
//
//    -----------------------------------------------------------------------------

#ifndef SBBEPIKernelASL2D_h
#define SBBEPIKernelASL2D_h 1

#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include "MrImaging/SequenceLibraries/libASL/SBBPasl_CrushGrad.h" 

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control




namespace SEQ_NAMESPACE
{

    class __IMP_EXP SBBEPIKernelASL2D:public SeqBuildBlockEPIKernel
    {
    public:

        //--------------------------------------------------------------------
        //    Constructor
        SBBEPIKernelASL2D(SBBList* pSBBList);

        //--------------------------------------------------------------------
        //  Destructor
        virtual ~SBBEPIKernelASL2D() {};

        // init the excitation SBBs
        virtual bool initExcitation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        // preparation of the SBBDiffusion
        virtual bool prepPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        // run the SBBDiffusion
        virtual bool runRTEBPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

        // Calculate PlugIn RF info (required for dynamic adjustments)
        virtual NLS_STATUS calcSliceAdjPlugInRFInfo(MrProt &, SeqLim &, SeqExpo &, const SLICEADJ::sCuboidGeometry &, const SLICEADJ::sAdjParametersMask &, std::vector<MrProtocolData::SeqExpoRFInfo> &);

	protected:
		LIB_NAMESPACE::SeqBuildBlockCrushGrad      m_CrushGradASL;
        // excitation pulse
        sRF_PULSE_SINC              m_sExcitationRF;

    };

};

#endif
