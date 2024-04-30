//-----------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2020. All Rights Reserved.
//-----------------------------------------------------------------------------

#include "SBBEPIKernelDiffusion.h"
#include "DiffusionSBBContainer.h"

#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"

#include "DiffusionRFPulseProperties.h"
#include "MrImaging/libSBB/SBBBinomialPulses.h"

#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

using namespace std;

namespace SEQ_NAMESPACE
{

    SBBEPIKernelDiffusion::SBBEPIKernelDiffusion(SBBList* pSBBList)
        :SeqBuildBlockEPIKernel(pSBBList)
    {
        m_bPlugInAvailable = true;
        m_sExcitationRF->setIdent("ExtExciteRF");
    }

    bool SBBEPIKernelDiffusion::prepPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
    {
        MrProtFacade protFacade(rMrProt);

        // --------------------------------------------------------------------------
        // configure SBBDiffusion
        // --------------------------------------------------------------------------
        if(!m_Diffusion.create(rMrProt))
        {
            SEQ_TRACE_ERROR << "ERROR m_Diffusion.create(rMrProt) failed";
            return false;
        }

        m_Diffusion->setSpinPrepTimeus(m_lTEContributionBeforeRTEBPlugIn);
        m_Diffusion->setADCusTillEcho(m_lTEContributionAfterRTEBPlugIn);
        m_Diffusion->setNoiseThreshold(rMrProt.diffusion().getlNoiseLevel());

        // set up compensation related parameters
        if (getbCompensationEnable())
        {
            m_Diffusion->setbCompensationEnable(true);
            m_Diffusion->setCompensationPara(m_bCompensationDecay, m_dCompensationFraction, m_dEddycurrentTau);
            m_Diffusion->setlPlugInToCompGradTime(m_lRTEBPlugInToCompGradTime);
        }
        else
        {
            m_Diffusion->setbCompensationEnable(false);
        }

        // Enable TE minimization
        if(rMrProt.TOM() != SEQ::TOM_OFF)
        {
            m_Diffusion->setMinimizeTE(true);
        }
        else
        {
            m_Diffusion->setMinimizeTE(false);
        }

        if(getUseGPABalance())
        {
            // Provide diffusion module with the readout gradient events
            m_Diffusion->setKernelGPALoad(m_sSBBBalance);

            // Enable GPA balance calculations
            m_Diffusion->setUseGPABalance(true);
        }
        else
        {
            // Disable GPA balance calculations
            m_Diffusion->setUseGPABalance(false);
        }

        // set parameters needed before the preparation of compensation gradients
        // e.g. min rise time, max amplitude.
        if (getbCompensationEnable())
        {
            m_Diffusion->prePrepareCompGrad(rMrProt, rSeqLim);
        }

        // --------------------------------------------------------------------------
        // prepare SBBDiffusion
        // --------------------------------------------------------------------------
        if(protFacade.isSliceAdj())
        {
            // set requested dynamic adjustments
            m_Diffusion->setsSliceAdjParametersRequestedBySequence(getsSliceAdjParametersRequestedBySequence());

            // prepare with given cuboids (because SBB is in HOLD_OPTIMIZATION mode)
            std::vector<SLICEADJ::sCuboidGeometry> vsSliceAdjCuboidsFromKernel;
            getSliceAdjCuboids(vsSliceAdjCuboidsFromKernel);
            if(! m_Diffusion->prep(rMrProt, rSeqLim, rSeqExpo, vsSliceAdjCuboidsFromKernel))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERR << "ERROR m_Diffusion->prep failed with " << m_Diffusion->getNLSStatus();
                }

                return false;
            }
        }
        else
        {
            if(! m_Diffusion->prep(rMrProt, rSeqLim, rSeqExpo))
            {
                if(!rSeqLim.isContextPrepForBinarySearch())
                {
                    SEQ_TRACE_ERR << "ERROR m_Diffusion->prep failed with " << m_Diffusion->getNLSStatus();
                }

                return false;
            }
        }

        m_lRTEBPlugInDurationPerRequest      = m_Diffusion->getDurationPerRequest();
        m_lRTEBPlugInTEContribution          = m_Diffusion->getTEContributionPerRequest();
        m_lRTEBPlugInStorageTime             = m_Diffusion->getStorageTimePerRequest();
        m_bRTEBPlugInInvertsMagnetization    = m_Diffusion->isMagnetizationInverted();
        m_lTRIncrement                       = m_Diffusion->getTRIncrement();

        return true;
    }

    bool SBBEPIKernelDiffusion::initExcitation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
    {
        // ===============================================
        // set pulse properties
        // ===============================================
        sRFPulseProperties myRFPulseProperties;

        if(rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_WaterExcitation)
            myRFPulseProperties = m_RFPulseLibrary.getPulsePropertiesExcitation(rMrProt, m_SBBExcitationBinomial.getMaxMagnitude(SEQ::GRAD_FAST));
        else
            myRFPulseProperties = m_RFPulseLibrary.getPulsePropertiesExcitation(rMrProt, m_SBBExcitation.getMaxMagnitude(SEQ::GRAD_FAST));

        m_sExcitationRF.setFamilyName(myRFPulseProperties.sFamilyName.c_str());
        m_sExcitationRF.setDuration(myRFPulseProperties.lDuration_us);
        m_sExcitationRF.setThickness(rMrProt.sliceSeries().aFront().thickness());
        m_sExcitationRF.setTypeExcitation();
        m_sExcitationRF.setFlipAngle(rMrProt.flipAngle());
        m_sExcitationRF.setInitialPhase(90.0);

        // Use positive slice selection gradient polarity for excitation pulses, 
        // negative polarity for refocusing pulses (if gradient reversal is active) will be set
        // in the SBB playing out the refocusing pulse.
        m_sExcitationRF.setRequiredGSPolarity(+1.0);

        // For water excitation we have to use a separate SBB which supports 
        // binomial pulses. All other excitation modes use the SBB which is able
        // to use multi-band RF pulses.
        if(rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_WaterExcitation)
        {
            // tell the excitation SBB about slice thickness and gradient polarity
            m_SBBExcitationBinomial.setThickness(rMrProt.sliceSeries().aFront().thickness());
            m_SBBExcitationBinomial.setRequiredGSPolarity(+1.0);

            // pass the pulse to the excitation SBB
            if(!m_SBBExcitationBinomial.setExcitationRFPulse(&m_sExcitationRF, rMrProt, rSeqExpo))
            {
                SEQ_TRACE << "m_SBBExcitationBinomial.setRFPulse " << m_SBBExcitationBinomial.getNLSStatus();
                return false;
            }

            m_pSBBExcite = &m_SBBExcitationBinomial;
        }
        else
        {
            if(!m_SBBExcitation.setVERSEParam(myRFPulseProperties.bIsVERSE, myRFPulseProperties.fFactorVERSE, myRFPulseProperties.fRelPlateauLengthVERSE))
            {
                SEQ_TRACE << "m_SBBExcitation.setVERSEParam " << m_SBBExcitation.getNLSStatus();
                return false;
            }

            if(!m_SBBExcitation.setMultibandParam(SMSProperties::getMultiBandFactor(rMrProt), SMSProperties::getSliceSeparation(rMrProt)))
            {
                SEQ_TRACE << "m_SBBExcitation.setMultibandParam " << m_SBBExcitation.getNLSStatus();
                return false;
            }

            m_SBBExcitation.setSMSPulsePreDampening(SMSProperties::isSuitableForPulsePredampening());

            if(!m_SBBExcitation.setBandwidthOptimizations(true))
            {
                SEQ_TRACE << "m_SBBExcitation.setBandwidthOptimizations " << m_SBBExcitation.getNLSStatus();
                return false;
            }



            if(!m_SBBExcitation.setRFPulse(rMrProt, rSeqLim, rSeqExpo, m_sExcitationRF))
            {
                SEQ_TRACE << "m_SBBExcitation.setRFPulse " << m_SBBExcitation.getNLSStatus();
                return false;
            }

            m_pSBBExcite = &m_SBBExcitation;
        }

        return true;
    }


    bool SBBEPIKernelDiffusion::runRTEBPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
    {
        // --------------------------------------------------------------------------
        // SBBDiffusion completely consumes the TE-fill times available, also it has
        // no methods to specify additional fill times.
        // So, if the RTEBPlugIn is told to insert fill times, we have a problem.
        // --------------------------------------------------------------------------
        if(m_lRTEBPlugInTEFillBefore!=0)
        {
            SEQ_TRACE_WARN <<  "WARNING: lRTEBPlugInTEFillBefore=" << m_lRTEBPlugInTEFillBefore << " (should be 0).";
        }
        if(m_lRTEBPlugInTEFillAfter!=0)
        {
            SEQ_TRACE_WARN << "WARNING: lRTEBPlugInTEFillAfter=" << m_lRTEBPlugInTEFillAfter << " (should be 0).";
        }

        if(! m_Diffusion->setADCforDiffusionMDHentries(&m_ADC))
        {
            SEQ_TRACE_ERROR << "setADCforDiffusionMDHentries failed with " << m_Diffusion->getNLSStatus();
            return false;
        }

        //---------------------------------------------------------------------------
        // store start time of NEXT RTEB to be inserted into the sequence timing
        //---------------------------------------------------------------------------
        const double dStartOfRTEBPlugIn = RTController::getInstance().getAbsTimeOfEventBlockMSec();

        if(! m_Diffusion->run(rMrProt, rSeqLim, rSeqExpo, pSLC))
        {
            // m_pSeq->m_lDiffLoopCounter and m_pSeq->m_lRepLoopCounter are global variables which have been set in fSeqRunKernel
            SEQ_TRACE_ERROR << "m_Diffusion failed with " << m_Diffusion->getNLSStatus();
            return false;
        }

        //---------------------------------------------------------------------------
        // Check duration of inserted event block(s):
        //
        // We spent lots of effort to do a correct calculation of TE-fill-times and
        // kernel duration. Therfore we expect the RTEB-plug-in to insert the correct
        // amount of time into the sequence timing.
        //---------------------------------------------------------------------------
        const long lRTEBsDuration    = static_cast<long>((RTController::getInstance().getAbsTimeOfEventBlockMSec() - dStartOfRTEBPlugIn) * 1000.0 + 0.5);

        const long lExpectedDuration =   m_lRTEBPlugInTEFillBefore + m_lRTEBPlugInDurationPerRequest + m_lRTEBPlugInTEFillAfter;

        if(lRTEBsDuration != lExpectedDuration)
        {
            SEQ_TRACE << "ERROR: runRTEBPlugIn did not insert the correct amount of time into sequence timing:\n"
                << "      expected duration" << lExpectedDuration << "\n" 
                << "      got duration" << lRTEBsDuration << "\n" 
                << "possible reasons for this error:\n" 
                << "- SeqBuildBlockEPIKernel::prep calculated wrong TE-fill-times:\n"
                << "  m_lRTEBPlugInTEFillBefore" << m_lRTEBPlugInTEFillBefore << "\n"
                << "  m_lRTEBPlugInTEFillAfter" << m_lRTEBPlugInTEFillAfter << "\n"
                << "- runRTEBPlugIn did not calculate correct duration:\n"
                << "  m_lRTEBPlugInDurationPerRequest" << m_lRTEBPlugInDurationPerRequest << "\n"
                << "- runRTEBPlugIn did not insert TE-fill times";

                setNLSStatus(MRI_SEQ_SEQU_ERROR);
            return false;
        }

        return true;
    }

    MrProtocolData::SeqExpoRFInfo SBBEPIKernelDiffusion::getRFInfoPerRequest()
    {
        if(m_pSBBExcite)
        {
            if(m_bIsSliceAcceleration && m_eRunMode == MULTI_BAND)
                return static_cast<SBBMultibandRF*>(m_pSBBExcite)->getRFInfoPerRequestMB() + m_Diffusion->getRFInfoPerRequestMB();
            else
                return m_pSBBExcite->getRFInfoPerRequest() + m_Diffusion->getRFInfoPerRequest();
        }
        else
            return MrProtocolData::SeqExpoRFInfo();
    }

    long SBBEPIKernelDiffusion::getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch)
    {
        return m_Diffusion->getSmallestIVIMbValuePossible(bIsContextPrepForBinarySearch);
    }
    
    long SBBEPIKernelDiffusion::getIVIMIncrement()
    {
        return m_Diffusion->getIVIMIncrement();
    }

    long SBBEPIKernelDiffusion::getMaxBValueSmallIVIMIncrement()
    {
        return m_Diffusion->getMaxBValueSmallIVIMIncrement();
    }

    void SBBEPIKernelDiffusion::increaseSliceSelectionGradientMaxAmplitude(MrProt& rMrProt)
    {
        MrProtFacade protFacade(rMrProt);

        double dGradMaxMagn = getMaxMagnitude(rMrProt.gradSpec().mode(), SBBBinomialPulses_GRAD_PERF_BINOMIAL);
        double dAdditionalScaleFactor = 1.0;
        if (protFacade.isGradientReversalDiffusion() && !protFacade.isSliceAcceleration())
        {
            dAdditionalScaleFactor = 1.1;
        }
        else
        {
            dAdditionalScaleFactor = 1.0;
        }

        double adMagnitudes[3]
            = {dAdditionalScaleFactor * dGradMaxMagn,
               dAdditionalScaleFactor * dGradMaxMagn,
               dAdditionalScaleFactor * dGradMaxMagn};
        m_SBBExcitation.setMaxMagnitudes(adMagnitudes, SBBBinomialPulses_GRAD_PERF_BINOMIAL);
        m_SBBExcitation.setMaxMagnitudes(adMagnitudes, SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
        m_SBBExcitationBinomial.setMaxMagnitudes(adMagnitudes, SBBBinomialPulses_GRAD_PERF_BINOMIAL);
        m_SBBExcitationBinomial.setMaxMagnitudes(adMagnitudes, SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
    }

#ifdef ZOOM_2DRF
    bool SBBEPIKernelDiffusion::isOptPTXVolumeCondition(MrProt& rMrProt)
    {
        return (
            SeqBuildBlockEPIKernel::isOptPTXVolumeCondition(rMrProt)
            || (rMrProt.getsAdjData().getlCoupleAdjVolTo() == MrProtocolData::AdjVolCoupling_PTxVolume));
    }
#endif

    } // namespace SEQ_NAMESPACE