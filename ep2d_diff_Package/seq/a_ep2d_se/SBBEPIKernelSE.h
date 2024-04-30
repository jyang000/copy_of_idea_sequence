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

#pragma once

#include "SBBRefocSE.h"

#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control
#include <memory>

#include <memory>

namespace SEQ_NAMESPACE
{

    class __IMP_EXP SBBEPIKernelSE:public SeqBuildBlockEPIKernel
    {
    public:

        //--------------------------------------------------------------------
        //    Constructor
         SBBEPIKernelSE(SBBList* pSBBList, bool bInitBaseSBBRefocSE = true);

        //--------------------------------------------------------------------
        //  Destructor
        virtual ~SBBEPIKernelSE() = default;
        ;

    // init the excitation SBBs
    virtual bool initExcitation(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    // preparation of the plug in SBB
    virtual bool prepPlugIn(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    // run plug in
    virtual bool runRTEBPlugIn(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC);

    // RF info
    virtual MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest();

    protected:

        // instantiating refocusing SBB
        void initSBBRefocSE();

        virtual bool configureRefocusingRFPulse(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

        virtual bool prepSBBRefocus(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

        virtual void configureSpoilerGradientPolarity();

        // excitation pulse
        sRF_PULSE_EXT m_sExcitationRF;

    // refocusing pulse
    sRF_PULSE_EXT m_sRefocusingRF;

        // refocusing SBB
        std::unique_ptr<SBBRefocSE> m_pSBBRefocus;


    };

}; // namespace SEQ_NAMESPACE