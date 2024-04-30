//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 2013  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\libSBB\SBBMultibandRF.h
//     Version: \main\2
//     Authors: Virtual Neuro Team:
//					Himanshu Bhat
//					Dingxin (Guilong) Wang
//					Thomas Beck
//                  Mario Zeller
//                  Uvo Hoelscher
//        Date: 2015-01-28 17:03:29 +01:00
//
//        Lang: C++
//
//     Descrip: MR::MrServers::MrImaging::seq::common::SliceAcceleration
//              SBB for multi band RF pulses -> uses sRF_PULSE_MB to create multiband RF pulses
//				Supports VERSE
//				Holds and manages multiple RF pulses to enable RF CAIPIRHINA based FOV shifting
//
//     Classes: SBBMultibandRF
//
//    -----------------------------------------------------------------------------


#ifndef SBBMultibandRF_h
#define SBBMultibandRF_h 1

#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"
#include "MrImaging/libSBB/SBBExcitation.h"
#include "MrImaging/libSBB/libSBBmsg.h" 
#include "MrMeasSrv/SeqIF/csequence.h"        
#include "MrGlobalDefinitions/ImpExpCtrl.h" 
#include "MrImaging/libSeqUtil/sRF_PULSE_MB.h"
#include "MrImaging/libSeqUtil/SMSProperties.h"

#ifdef BUILD_libSBB
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h" 
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpoRFBlockInfo.h"

class SBBTestHelperMultibandRF;      // SBB iTest Helper Class to Dump Class Member

class __IMP_EXP_EXPORT_DECL SBBMultibandRF : public SeqBuildBlockExcitationRFPulse
{

public:

    SBBMultibandRF (SBBList* pSBBList = NULL);

    virtual ~SBBMultibandRF();

    // set basic RF pulse pointer, copy all RF samples to the internally used single- and multi-band pulses
    virtual bool setRFPulse (MrProt & rMrProt, SeqLim & rSeqLim, SeqExpo & rSeqExpo, IRF_PULSE* pRF);

    //Setup multi band RF pulse design parameters
    virtual bool setMultibandParam(long lMultibandFactor, double dMultibandDistance);

    //Setup VERSE RF pulse design parameters
    virtual bool setVERSEParam(bool bIsVerse, double dVerseFactor, double dRelativePlateauLength = 0.4);

    // prepare the MB pulse object, VERSE is turned off for binary search
    virtual bool prep_insertEvents (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo& );

    // run the MB pulse object
    virtual bool run_insertEvents (MrProt& rMrProt, SeqLim &rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC);

    //checks if the duration is on double gradient raster time and if not increases it
    virtual bool adjustDurationToGradRasterTime();

    //Check if the peak RF power is over amplifier limit
    virtual bool isRFClipped() const;

    // return VERSE value
    virtual bool isVerse() const;

    // return pointer to the used RF object
    virtual IRF_PULSE* getRFPulsePointer(void);  

    //Setup RF run mode
    virtual void setRunMode(SliceAccelRFRunMode eRunMode);

    virtual void setB1ControlLoopActiveForSMSSingleBandScans(bool bActive);

    // returns the RF info for MB run mode, SB RF info is handled by the base SBB
    virtual MrProtocolData::SeqExpoRFInfo getRFInfoPerRequestMB();

    // returns pointer to used slice selection gradient; ID is used for VERSEd gradients
    virtual sGRAD_PULSE_TRAP* getGradientPointer(int GradPulseID = 0);

    // flag for bandwidth optimized base band
    virtual bool setBandwidthOptimizations(bool bOptimize);

    /// Overloaded base class method: Calculate RF info for given cuboid geometry
    virtual bool calcSliceAdjSBBRFInfo(
        MrProt                                     &rMrProt,        //< IMP: points to the protocol structure.
        SeqLim                                     &rSeqLim,        //< IMP: points to the sequence limits structure.
        SeqExpo                                    &rSeqExpo,       //< IMP: points to the sequence exports structure
        const SLICEADJ::sCuboidGeometry              &sSliceAdjCuboid,  //< IMP: cuboid geometry (single)
        std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo        //< EXP: RF info
        );

	// sets the flag for SMS pulse pre-distortion in the SMS pulse member variable m_sMultiBandRF
	virtual void setSMSPulsePreDampening(bool bSMSPulsePreDampening);

protected: 

    virtual void modifySlicePosition(MrProt & rMrProt, sSLICE_POS* pSLC);

    // copy RF samples to single-band pulse and prepare the pulse for VERSE
    virtual bool copyParametersAndPrepareForVERSE (MrProt & rMrProt, SeqExpo & rSeqExpo, sRF_PULSE * pInputRF, sRF_PULSE_ARB & rOutputRF, double dVerseFactor, double dRelativePlateauLength = 0.4);

    //Slice acceleration factor equals to number of RF bands
    long m_lMultibandFactor;

    //Distance between RF bands in mm
    double m_dMultibandDistance;

    //VERSE related parameters
    bool m_bIsVerse;
    double m_dVerseFactor;
    double m_dVerseRelativePlateauLength;

    //RF run mode that decides which RF pulse to be executed
    SliceAccelRFRunMode m_eRunMode;

    // multi-band RF pulse
    sRF_PULSE_MB m_sMultiBandRF;

    // single-band RF pulse
    sRF_PULSE_ARB m_sSingleBandRF;

    // RF samples of the single-band pulse (1st half for non-versed version; 2nd half for versed version)
    sSample	m_sSingleBandRFSamples[IRF_MAX_ENVELOPE_ELEMENTS * 2];

    // flag for bandwidth optimized base band
    bool m_bOptimizeBandwidth;

    // for SliceAdj we need to store the MB energy in a separate variable, SB energy is handled by the base SBB mechanism
    // (SMS and SliceAdj are not designed to be used at the same time)
    MrProtocolData::SeqExpoRFInfo m_RFInfoPerRequestMB;

    bool m_bB1ControlLoopActiveForSMSSingleBandScans;

private: 
    SBBMultibandRF(const SBBMultibandRF &right);

    SBBMultibandRF & operator=(const SBBMultibandRF &right);

    friend class SBBTestHelperMultibandRF;      // SBB iTest Helper Class to Dump Class Member

};// Class SBBMultibandRF

#endif
