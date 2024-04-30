//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 2012  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\common\Excitation\a_ep2d_zoom_UINS.h
//     Version: \main\11
//      Author: pfeujodj
//
//        Lang: C++
//
//
//
///  \file   a_ep2d_zoom_UINS.h
///  \brief  File declaring the UI class for EP2D with 2D RF pulses 
///
//    -----------------------------------------------------------------------------


#ifndef a_ep2d_zoom_UINS_h
#define a_ep2d_zoom_UINS_h 1


#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"
#ifdef WIN32
    #include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkMap.h"
#endif

namespace Ep2d_zoom_UINS
{

#ifdef WIN32

    //-----------------------------------------------------------------------------------
    // register standard and WIP solve handlers for ZOOM_2DRF option
    //-----------------------------------------------------------------------------------
    void fEPIRegisterZoomHandlers(SeqLim &rSeqLim);

    //-----------------------------------------------------------------------------------
    // initializes UI functions and members for ZOOM_2DRF option
    //-----------------------------------------------------------------------------------
    NLS_STATUS initializeZoomUI(MrProt &rMrProt, SeqLim &rSeqLim);

    //	calculate PTXVolume by fitting it to the first slice group (that has to be in the protocol)
    void calcPTXVolFromSliceGroup (MrProt &rMrProt);

    //	calculate PTXMPR slice group by fitting it to the slice group and PhaseFOV of the PTXVolume
    void calcPTXMPRFromSliceGroup (MrProt &rMrProt);

    //	calculate outer-volume Sats (ovSat) from slice group
    bool calcOvSatFromSliceGroup (MrProt &rMrProt, SeqLim &rSeqLim);

    //  enables PTX pulses
    void setPTXCalculation (MrProt &rMrProt);

    void _refreshPTXVol(MrUILinkBase* const pThis);
  

    //-----------------------------------------------------------------------------------
    // implementation of overloaded solve handlers 
    //-----------------------------------------------------------------------------------
    //  MR_TAG_EXCIT_PULSE: Excitation Pulse
    //
    void _fUILinkInsertOptimizationVolume             (MrUILinkBase* const pThis);
    unsigned fUILinkExcitationPulseSetValueNew        (LINK_SELECTION_TYPE* const pThis, unsigned val, int32_t pos);

    //-----------------------------------------------------------------------------------
    //  MR_TAG_PTX_VOLUME_ARRAY: PTX Volumes
    //
    bool fUILinkPTXVolPropGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolArrayCanInsertNew(MrUILinkArray* const pThis, int32_t newPos, int32_t dummy);
    bool fUILinkPTXVolArrayCanEraseNew (MrUILinkArray* const pThis, int32_t doomed, int32_t pos);
    
	// make PTXVolume phase valid for binary search or make field not editable in case of rotated trajectory
	bool fUILinkPTXVolPDimGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
         
    // make PTXVolume fields NOT editable for type optimization
    bool fUILinkPTXVolRDimGetLimitsNew       (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolSDimGetLimitsNew       (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolRotGetLimitsNew        (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolPosSagGetLimitsNew     (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolPosSag_SBCSGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolPosCorGetLimitsNew     (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolPosCor_SBCSGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolPosTraGetLimitsNew     (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolPosTra_SBCSGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolOriAlphaGetLimitsNew   (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolOriBetaGetLimitsNew    (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos);

    bool fUILinkPTXVolPosGetOptionsNew       (LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolPos_SBCSGetOptionsNew  (LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolOriGetOptionsNew       (LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);
    bool fUILinkPTXVolOriHistoryGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);

#ifdef ZOOM_EXTENDED
	bool fUILinkPTXVolVisibilityGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);
#endif // ZOOM_EXTENDED

#ifdef ZOOM_EXTENDED
	// the following fields are not required for the rotated trajectory -> make them not available
	//bool fUILinkPTXVolPropIsAvailableNew(LINK_SELECTION_TYPE* const pThis, int32_t pos);
	bool fUILinkPTXVolPDimIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
	bool fUILinkPTXVolRotIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
	bool fUILinkPTXVolRDimIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
	bool fUILinkPTXVolSDimIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
	bool fUILinkPTXVolPosIsAvailableNew(LINK_SELECTION_TYPE* const pThis, int32_t pos);
	bool fUILinkPTXVolOriIsAvailableNew(LINK_SELECTION_TYPE* const pThis, int32_t pos);
	bool fUILinkPTXVolVisibilityIsAvailableNew(LINK_SELECTION_TYPE* const pThis, int32_t pos);
#endif // ZOOM_EXTENDED

    // update PTXVolume fields by calling _refreshPTXVol() / calcPTXVolFromSliceGroup()
    double fUILinkPTXVolRotGetValueNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
    double fUILinkPTXVolPDimGetValueNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
    double fUILinkPTXVolRDimGetValueNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
    double fUILinkPTXVolSDimGetValueNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);

    // UI handlers for GSP parameters
    GradDirPat         fUILinkPTXVolGSPRotGetValueNew(LINK_GRAD_DIR_TYPE* const pThis, int32_t lPos);
    VectorPat<double>  fUILinkPTXVolGSPPosGetValueNew(LINK_VECTOR_TYPE* const pThis, int32_t pos);

    // update PTXPulse fields / make them non-zero-sized with EXCITATION_ZOOMED
    int32_t fUILinkPTXPulseArraySizeNew(MrUILinkArray*    const pThis, int32_t pos);
    bool fUILinkPTXPulseIsAvailableNew(MrUILinkMap*      const pThis, int32_t pos);
    double fUILinkPTXPulseTxAccSetValueNew(LINK_DOUBLE_TYPE* const pThis, double newVal, int32_t pos);
    bool fUILinkPTXPulseTxAccGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t);
#ifdef ZOOM_EXTENDED
	bool fUILinkPTXPulseTxAccIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos);
#endif

    // make PTXPulse fields NOT editable
    bool fUILinkPTXPulseFlipAngleGetLimitsNew       (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify,int32_t pos);
    bool fUILinkPTXPulsePhaseFoEGetLimitsNew        (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify,int32_t pos);
    bool fUILinkPTXPulsePhaseMatrixSizeGetLimitsNew (LINK_LONG_TYPE*   const pThis, std::vector<MrLimitLong>&   rLimitVector, uint32_t& rulVerify,int32_t pos);

    unsigned fUILinkPhaseFOV_GetToolTip(LINK_DOUBLE_TYPE* const pThis, char * arg_list[], int32_t /*pos*/);

#endif   // of #ifdef WIN32

} // end of namespace Ep2d_diff_UINS

#endif  // of #ifndef a_ep2d_zoom_UINS_h

