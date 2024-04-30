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


#ifndef SBBEPIKernelSE_h
#define SBBEPIKernelSE_h 1

#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"


#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control
#include "../a_ep2d_se/SBBRefocSE.h"



namespace SEQ_NAMESPACE
{

    class __IMP_EXP SBBEPIsegKernelSE:public SeqBuildBlockEPIKernel
    {
    public:

        //--------------------------------------------------------------------
        //    Constructor
        SBBEPIsegKernelSE(SBBList* pSBBList);

        //--------------------------------------------------------------------
        //  Destructor
        virtual ~SBBEPIsegKernelSE() {};

        // init the excitation SBBs
        virtual bool initExcitation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        // preparation of the plug in SBB
        virtual bool prepPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        // run plug in
        virtual bool runRTEBPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

        // RF info
        virtual MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest();

    protected:
        // refocusing pulse
        sRF_PULSE_EXT m_sExcitationRF;

        // refocusing pulse
        sRF_PULSE_EXT m_sRefocusingRF;

        // refocusing SBB
        SBBRefocSE m_SBBRefocus;


    };

};

#endif
