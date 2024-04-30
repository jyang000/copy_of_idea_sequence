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


#ifndef SBBEPIKernelFID_h
#define SBBEPIKernelFID_h 1

#include "MrImaging/seq/Kernels/SBBEPIFCKernel.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control



namespace SEQ_NAMESPACE
{

    class __IMP_EXP SBBEPIsegKernelFID:public SeqBuildBlockEPIFCKernel
    {
    public:

        //--------------------------------------------------------------------
        //    Constructor
        SBBEPIsegKernelFID(SBBList* pSBBList);

        //--------------------------------------------------------------------
        //  Destructor
        virtual ~SBBEPIsegKernelFID() {};

        // init the excitation SBBs
        virtual bool initExcitation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);
     

    protected:
        // excitation pulse
        sRF_PULSE_SINC m_sExcitationRF;
    };

};

#endif
