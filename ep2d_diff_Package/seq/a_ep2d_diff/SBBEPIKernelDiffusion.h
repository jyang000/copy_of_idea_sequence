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

#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include "DiffusionSBBContainer.h"
#include "MrImaging/libSeqUtil/SMSProperties.h"
#include "didi.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"

// Eddy current compensation gradients
#include "MrImaging/seq/a_ep2d_diff/SBBCompGrad.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

namespace SEQ_NAMESPACE
{

    class __IMP_EXP SBBEPIKernelDiffusion:public SeqBuildBlockEPIKernel
    {
    public:

        //--------------------------------------------------------------------
        //    Constructor
        SBBEPIKernelDiffusion(SBBList* pSBBList);

        //--------------------------------------------------------------------
        //  Destructor
        virtual ~SBBEPIKernelDiffusion() {};

        // preparation of the SBBDiffusion
        bool prepPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;

        // init the excitation SBBs
        bool initExcitation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;

        // run the SBBDiffusion
        bool runRTEBPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC) override;

        MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest() override;

        virtual long getTotalScans(bool bIncludeLocalAverages = true)
        {
            return m_Diffusion->getTotalScans(bIncludeLocalAverages);
        }

        virtual void setAdjustmentScan(int iAdjScan)
        {
            m_Diffusion->setAdjustmentScan(iAdjScan);
        }

        virtual long getAdjustmentScan()
        {
            return m_Diffusion->getAdjustmentScan();
        }

        virtual long getNoOfAdjScans(void)
        {
            return m_Diffusion->getNoOfAdjScans();
        }
        
        virtual bool getDiffusionSeqLims(SeqLim &rSeqLim)
        {
            return m_Diffusion.getDiffusionSeqLims(rSeqLim);
        }

        virtual void setDiffusionGradientsEnabled(bool bEnabled)
        {
            m_Diffusion->setDiffusionGradientsEnabled(bEnabled);
        }

        void setRunMode(SliceAccelRFRunMode eRunMode) override
        {
            m_eRunMode = eRunMode;
            m_Diffusion->setRunMode(m_eRunMode);
        }

        virtual DiffusionDirections * getDidiPointer()
        {
            return m_Diffusion->getDidiPointer();
        }

        // get diffusion gradient duration in ms for display in UI
        virtual double getDiffGradDuration_ms()
        {
            return m_Diffusion->getDiffGradDuration_ms();
        }

        // get diffusion gradient spacing in ms for display in UI
        virtual double getDiffGradSpacing_ms()
        {
            return m_Diffusion->getDiffGradSpacing_ms();
        }

        // set loop counters
        virtual bool setLoopCounters
            (
            int       iAdjustmentScan,    /**< Imp:     Adjustment scan index            */
            long      lRepLoopCounter,    /**< Imp:     SeqLoop repetition counter       */
            long      lDiffLoopCounter   /**< Imp:     SeqLoop diffusion counter        */
            )
        {
            return m_Diffusion->setLoopCounters(iAdjustmentScan, lRepLoopCounter, lDiffLoopCounter, &m_ADC);
        }

        // get smallest b value for which IVIM is possible
        long getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch = false);

        // get IVIM increment
        long getIVIMIncrement();

        // get max b value for small IVI increments
        virtual long getMaxBValueSmallIVIMIncrement();

        // get diffusion loop counter for highest b-value for kernel check
        virtual long getDiffLoopCounterForHighestBValue()
        {
            return m_Diffusion->getDiffLoopCounterForHighestBValue();
        }

        virtual long getDiffLoopCounterForHighestReadComponent(sSLICE_POS* pSlice)
        {
            return m_Diffusion->getDiffLoopCounterForHighestReadComponent(pSlice);
        }

        // set flag for thermal balancing
        virtual void setThermalBalancing(bool bThermalBalancing)
        {
            m_Diffusion->setThermalBalancing(bThermalBalancing);
        }

        // get flag for thermal balancing
        virtual bool getThermalBalancing()
        {
            return m_Diffusion->getThermalBalancing();
        }

        // set flag for TR filling, when compensation gradients used, TR fill will be done in compensation gradient SBB.
        virtual void setEPIReadOutAppliesTRFill(bool bTRFill) { m_bEPIReadOutAppliesTRFill = bTRFill; }

        // get scale factor for diffusion gradients in runtime, inculding the sign: dUsedDiffAmpl/m_dMaxAmpl;
        virtual std::vector<double> getvdScaleFactorinRun() { return m_Diffusion->getvdScaleFactorinRun(); }

        // set up compensation realted parameters
        void setbCompensationEnable(bool bComEnable) { m_bCompEnable = bComEnable; }

        bool getbCompensationEnable() { return m_bCompEnable; }

        void setCompensationPara(bool bCompensationDecay, double dCompensationFraction, double dEddycurrentTau)
        {
            m_bCompensationDecay    = bCompensationDecay;
            m_dCompensationFraction = dCompensationFraction;
            m_dEddycurrentTau       = dEddycurrentTau;
        };
        SeqBuildBlockCompGrad* getPointerCompGrad() { return m_Diffusion->getPointerCompGrad(); }

    protected:

#ifdef ZOOM_2DRF
      // check if the conditions of OptPTXVolume are present (different in diffusion and others)
      bool isOptPTXVolumeCondition(MrProt& rMrProt) override;
#endif

      // Enlarge the max amplitude of Gradient Slice of ep2d_diff to allow thinner slice thickness (3mm) for low-field
      void increaseSliceSelectionGradientMaxAmplitude(MrProt& rMrProt) override;

        DiffusionSBBContainer m_Diffusion;

        // excitation pulse
        sRF_PULSE_EXT m_sExcitationRF;

        // RF Pulse library which holds all all information about the used RF pulses
        DiffusionRFPulseProperties m_RFPulseLibrary;

        // compensation related parameters
        bool   m_bCompEnable{false};
        bool   m_bCompensationDecay{false};
        double m_dCompensationFraction{1.0};
        double m_dEddycurrentTau{0.0};

    };

};
