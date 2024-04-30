//----------------------------------------------------------------------------------
// <copyright file="a_ep_CommonUI.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 1998-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

// --------------------------------------------------------------------------
// General Includes
// --------------------------------------------------------------------------
#include "MrMeasSrv/SeqIF/Sequence/ISequence.h"
#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/UILink/MrStdNameTags.h"
#endif

#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkArray.h"
#include <vector>
#endif

#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrRXSpec.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Application/Application.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrCoilInfo.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/CoilSelect/MrRxCoilSelect.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrSysSpec.h"
#include "MrProtSrv/Domain/CoreNative/MrPat.h"
#include "MrProtSrv/Domain/CoreNative/MrFastImaging.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Filter/MrFilter.h"
#include "MrProtSrv/Domain/CoreNative/MrPreparationPulses.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Physiology/MrPhysiology.h"
#include "MrImaging/seq/common/IterativeDenoisingUIParameter/IterativeDenoisingUIParameter.h"
#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"

#include "MrMeasSrv/MeasUtils/MrTimeStamp.h"
#include "MrImaging/libSeqUtil/SMSProperties.h"
//  --------------------------------------------------------------------------
//  Application includes
//  --------------------------------------------------------------------------
#include "MrImaging/seq/a_ep_CommonUI.h"

#ifdef ZOOM_2DRF
#include "MrImaging/seq/common/Excitation/a_ep2d_zoom_UINS.h"   // UI definitions, WipMemBlock prot
#endif

//  --------------------------------------------------------------------------
//  Sequence specific includes
//  --------------------------------------------------------------------------
#if (defined SEQUENCE_CLASS_EP2D || defined SEQUENCE_CLASS_BASED_ON_EP2D)
    #include "MrImaging/seq/a_ep2d.h"
#elif defined SEQUENCE_CLASS_EP_SEG
    #include "MrImaging/seq/a_ep_seg.h"
#elif defined SEQUENCE_CLASS_AslCsl   // SUPPORT_iPAT_TGSE
#include "MrImaging/seq/a_tgse_asl/asl_csl.h"
#else
    #error SEQUENCE CLASS not recognised
#endif

//  --------------------------------------------------------------------------
//  Specify namespace
//  --------------------------------------------------------------------------
using namespace SEQ_NAMESPACE;

#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/SetProtocolParameter.h"
using namespace UICtrl;
#include "MrImagingFW/libSeqSysProp/SysProperties.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrMeasSrv/CoilIF/CoilSelectManipulator.h"
#include <limits>

// * ---------------------------------------------------------------------- *
// *                                                                        *
// * Name        :  SolveSliceAccelerationSettings                          *
// *                                                                        *
// * Description : The function sets all UI parameters such that            *
// *               slice acceleration can be turned on                      *
// *                                                                        *
// * Return      : MRI_STD_CONFIRMATION_MSG for success, else 0             *
// *                                                                        *
// * ---------------------------------------------------------------------- *
template<class UILinkParameterType>
unsigned EpCommonUINS::solveSliceAccelerationSettings(UILinkParameterType* const pThis)
{
    MrTimeStamp theTimeStamp(TIMESTAMP_SEQ, "EpCommonUINS::solveSliceAccelerationSettings");


    // when switching on SMS we need to turn off several other settings

    // general info
    MrProt rMrProt(&pThis->prot());

#ifdef EP2D_MS
    long lMultibandFactor = SMSProperties::isSMS( rMrProt) ? SMSProperties::getMultiBandFactor( rMrProt ) : 1;
    long lNumberOfSlices  = rMrProt.getsSliceArray().getlSize();

    // For (possibly interleaved) slice-selective preparations, limit the number of slices to even multiples of the multi-band factor.
    // See corresponding comments in ::fUILinkNumberOfSlicesGetLimits().
    if ( ( lMultibandFactor > 1 ) && ( rMrProt.sliceSeries().mode() == SEQ::INTERLEAVED ) )
    {
        if ( ( rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE ) || ( rMrProt.getsPrepPulses().getucSatRecovery() == SEQ::SATREC_SLICE_SELECTIVE ) )
        {
            lMultibandFactor *= 2;

            // first try to increase slices 
            long lRemainderSlices   = lNumberOfSlices % lMultibandFactor;
            long lRequiredSlices    = ((long)ceil((double)lNumberOfSlices / (double)lMultibandFactor)) * lMultibandFactor;
            if(lRemainderSlices != 0)
                setProtocolParameterElm<LINK_LONG_TYPE>(pThis, MR_TAG_SLICE_GROUP_LIST, MR_TAG_SG_SIZE, lRequiredSlices, 0);

            if(pThis->sequence().prepareForBinarySearch(&pThis->prot()))
                return (MRI_STD_CONFIRMATION_MSG);
        }
    }
#endif
    
    // concatenations
    setProtocolParameter<LINK_LONG_TYPE>(pThis, MR_TAG_CONCATENATIONS, 1, 0);

    // no optimization strategy
    setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_TOM, MRI_STD_TOM_OFF, 0);

    // no PTX
    setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_EXCIT_PULSE, MRI_STD_EXCITATION_STANDARD, 0);

    // no IR
#ifndef EP2D_MS
    setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_INVERSION, MRI_STD_MAGN_PREP_NONE, 0);
#endif

    // PACE
    setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_RESP_COMP, MIR_STD_NAV_OFF, 0);

    // turn off water excitation
    if (rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_WaterExcitation || rMrProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_FastWaterExcitation)
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_FAT_WATER_CONTRAST, MRI_STD_FATWATERCONTRAST_FAT_SATURATION, 0);
        
    adaptRefLinesPE(pThis, false);

    // diffusion mode
    if(rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_ONE_SCAN_TRACE ||
       rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_DIAGONAL)
    setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_DIFF_MODE, MRI_STD_DIFFUSION_FOUR_SCAN_TRACE, 0);

    // slice adjust
    setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_ADJUSTMENT_MODE, MRI_STD_ADJUSTMENT_MODE_STANDARD, 0);

    // now check if protocol is consistent
    if(pThis->sequence().prepareForBinarySearch(&pThis->prot()))
        return MRI_STD_CONFIRMATION_MSG;

    // adjust TE
    long lTEContrastIndex = 0;
#if (defined SEQUENCE_CLASS_EP2D || defined SEQUENCE_CLASS_BASED_ON_EP2D)
    lTEContrastIndex = SEQ_NAMESPACE::getSeq( pThis )->getTEContrastIndex( pThis->prot().getlContrasts() );
#endif

    long lNeededTE = getUI(pThis)->EpCommonUI::m_lNeededTE;
    if(lNeededTE > 0)
    {
        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TE, NULL, lNeededTE / 1000.0, lTEContrastIndex );
        if(pThis->sequence().prepareForBinarySearch(&pThis->prot()))
            return MRI_STD_CONFIRMATION_MSG;
    }

    // adjust TR
    long lNeededTR = getUI(pThis)->EpCommonUI::m_lNeededTR;
    if(lNeededTR > 0)
    {
        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR, NULL, lNeededTR / 1000.0, 0);
        if(pThis->sequence().prepareForBinarySearch(&pThis->prot()))
            return MRI_STD_CONFIRMATION_MSG;
    }

    return 0;
}

#ifdef SUPPORT_IIR

// Copy/Paste from a_tse_ui.cpp
struct EpCommonUINS::TE_TR_TI_HANDLER
{
    int32_t  m_mem;
    int32_t  m_index;
    MrProt*   m_pMrProt;
    MrMeasSrv::ISequence* m_pSeq;

    TE_TR_TI_HANDLER( MrProt* _pTempMrProt, MrMeasSrv::ISequence* pSeq )
        : m_mem( 0 )
        , m_index( 0 )
        , m_pMrProt( _pTempMrProt )
        , m_pSeq( pSeq )
    {}

    void storeInMemory()
    {
        switch ( m_index )
        {
        default:  m_mem = m_pMrProt->te()[m_index];  break;
        case -1:  m_mem = m_pMrProt->tr()[0];        break;
        case -2:  m_mem = m_pMrProt->ti()[0];        break;
        }
    }

    void recallMemory()
    {
        switch ( m_index )
        {
        default:  m_pMrProt->te()[m_index] = m_mem;  break;
        case -1:  m_pMrProt->tr()[0] = m_mem;        break;
        case -2:  m_pMrProt->ti()[0] = m_mem;        break;
        }
    }

    void clearMemory() {}

    int32_t value( int32_t val )
    {
        switch ( m_index )
        {
        default:  return m_pMrProt->te()[m_index] = val;  break;
        case -1:  return m_pMrProt->tr()[0] = val;        break;
        case -2:  return m_pMrProt->ti()[0] = val;        break;
        }
    }

    bool accept()
    {
        try
        {
            return ( m_pSeq->prepareForBinarySearch( *m_pMrProt ) != false );
        }
        catch ( ... )
        {
            std::cerr
                << "ERROR: unknown exception in file '" << __FILE__ << ", line " << __LINE__ << std::endl;
            return false;
        }
    }
};

// Copy/Paste from a_tse_ui.cpp
unsigned EpCommonUINS::fUISearchTR( MrUILinkBase* const _this, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex )
{
    unsigned _ret   = 0;
    int32_t  _inv;
    int32_t  _val;
    long     lMaxTR = 0;

    SEQ::PhysioSignal TrigSignalHigh = SEQ::SIGNAL_NONE;   // * Selected trigger modes    *
    SEQ::PhysioSignal TrigSignalLow  = SEQ::SIGNAL_NONE;   // * and signals               *
    SEQ::PhysioMethod TrigMethodHigh = SEQ::METHOD_NONE;
    SEQ::PhysioMethod TrigMethodLow  = SEQ::METHOD_NONE;

    MrProt sNewProt ( &_this->prot() );    // MrProt wrapper
    MrProt rOrigProt( pOrigProt      );

    SeqLim               &rSeqLim = _this->sequence().getSEQ_BUFF()->getSeqLim();
    MrMeasSrv::ISequence* pSeq    = &_this->sequence();

    // get Physio Mode
    rOrigProt.physiology().getPhysioMode( TrigSignalHigh, TrigMethodHigh, TrigSignalLow, TrigMethodLow );

    TE_TR_TI_HANDLER _handler( &sNewProt, pSeq );

    _inv = sNewProt.tr()[0];  // this tr does not work

    // -----------------------------------------------
    // try with maximum tr
    // -----------------------------------------------
    if ( TrigMethodHigh != SEQ::METHOD_NONE )
    {
        lMaxTR = rOrigProt.physiology().scanWindow( TrigSignalHigh ) * 1000 - rOrigProt.physiology().triggerDelay( TrigSignalHigh );
    }
    else
    {
        lMaxTR = static_cast<long>( rSeqLim.getTR()[0].getMax() );
    }

    sNewProt.tr()[0] = lMaxTR;

    _ret = _this->tryProt( const_cast<void*>( pVoid ), &_this->prot(), lIndex );  // note: &_this_prot() is intentionally used instead of pOrigProt! Otherwise the _try_template() may reject the protocol. This problem should be fixed!

    if ( !_ret )
    {
        //bNeedOtherTETITR = true;
        _ret = fUICStandardSolveHandler( _this, arg_list, pOrigProt, &_getTENeeded, &_getTINeeded,/*&_getTINeeded*/NULL );
        if ( !_ret )
        {
            sNewProt.tr()[0] = static_cast<int32_t>( _inv ); //restore old protocol value
            return ( 0 ); // can't solve with any tr
        }
    }

    // * -------------------------------------------------------------------- *
    // * Solve the problem:                                                   *
    // * -------------------------------------------------------------------- *
    // if ( arg_list != 0 )
    {
        _val = sNewProt.tr()[0];

        _handler.m_index = -1;
        //  Note: We assume something about TR-limit-handling in MrUILink!
        int32_t _inc = rSeqLim.getTR()[0].getInc();
        if ( _inv > 10000   ) _inc *= 10;
        if ( _inv > 1000000 ) _inc *= 10;
        binary_search( _handler, _inv, _val, _inc, BINARY_SEARCH_TRY_INVALID );
        double dValTR_ms = _val / 1000.;
        //  Note: The maximum TR is in general different for triggered and non-triggered measurements
        fUICSetWithinLimits( _this, dValTR_ms, MR_TAG_TR, UIC_SET_WITHIN_LIMITS_ROUNDUP, 0 );
        sNewProt.tr()[0] = static_cast<int32_t>( 1000 * dValTR_ms + 0.5 );
    }
    return ( MRI_STD_CONFIRMATION_MSG );
}

#endif

// ------------------------------------------------------------------------------
// Function    : _getTENeeded
// ------------------------------------------------------------------------------
//
// Description : callback function needed for libUICtrl's standard solve handlers
//               to get the required TE
// Return      : true , if parameter conflict can be solved with another TE
//               false, if not
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::_getTENeeded(MrProtocolData::MrProtData* pMrProt, SeqLim* rSeqLim, SeqExpo*, MrMeasSrv::ISequence* pSeq, long* alNeededTE)
{
    MrProt rMrProt( pMrProt );

    // Identify contrast index which determines TE
    long lTEContrastIndex=0;
#if (defined SEQUENCE_CLASS_EP2D || defined SEQUENCE_CLASS_BASED_ON_EP2D)
     lTEContrastIndex = static_cast<Ep2d*>( pSeq->getSeq() )->getTEContrastIndex( rMrProt.getlContrasts() );
#endif  

    if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"fGetTENeeded: returned ";

    if(    getUI(pSeq)->EpCommonUI::m_bNeedOtherTETITR
        && getUI(pSeq)->EpCommonUI::m_lNeededTE>=rSeqLim->getTE()[lTEContrastIndex].getMin()
        && getUI(pSeq)->EpCommonUI::m_lNeededTE<=rSeqLim->getTE()[lTEContrastIndex].getMax()
        )
    {
        for ( long lI = 0; lI < rMrProt.getlContrasts(); ++lI )
        {
            if ( lI == lTEContrastIndex )
            {
                alNeededTE[lI] = getUI( pSeq )->EpCommonUI::m_lNeededTE;
            }
            else
            {
                // Just to satisfy the calling UI methods
                alNeededTE[lI] = rMrProt.getalTE()[lI];
            }
        }
        if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"true"<<std::endl;
        return true;
    }
    else
    {
        if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"false"<<std::endl;
        return false;
    }
}


// ------------------------------------------------------------------------------
// Function    : _getTRNeeded
// ------------------------------------------------------------------------------
//
// Description : callback function needed for libUICtrl's standard solve handlers
//               to get the required TR
// Return      : true , if parameter conflict can be solved with another TR
//               false, if not
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::_getTRNeeded(MrProtocolData::MrProtData*, SeqLim* rSeqLim, SeqExpo*, MrMeasSrv::ISequence* pSeq, long* alNeededTR)
{
    if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"fGetTRNeeded: returned ";

    if(    getUI(pSeq)->EpCommonUI::m_bNeedOtherTETITR
        && getUI(pSeq)->EpCommonUI::m_lNeededTR>=rSeqLim->getTR()[0].getMin()
        && getUI(pSeq)->EpCommonUI::m_lNeededTR<=rSeqLim->getTR()[0].getMax()
        )
    {
        alNeededTR[0] = getUI(pSeq)->EpCommonUI::m_lNeededTR;
        if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"true"<<std::endl;
        return true;
    }
    else
    {
        if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"false"<<std::endl;
        return false;
    }
}


// ------------------------------------------------------------------------------
// Function    : _getTINeeded
// ------------------------------------------------------------------------------
//
// Description : callback function needed for libUICtrl's standard solve handlers
//               to get the required TI
// Return      : true , if parameter conflict can be solved with another TI
//               false, if not
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::_getTINeeded(MrProtocolData::MrProtData*, SeqLim* rSeqLim, SeqExpo*, MrMeasSrv::ISequence* pSeq, long* alNeededTI)
{
    if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"_getTINeeded: returned ";

    if(    getUI(pSeq)->EpCommonUI::m_bNeedOtherTETITR
        && getUI(pSeq)->EpCommonUI::m_lNeededTI>=rSeqLim->getTI()[0].getMin()
        && getUI(pSeq)->EpCommonUI::m_lNeededTI<=rSeqLim->getTI()[0].getMax()
        )
    {
        alNeededTI[0] = getUI(pSeq)->EpCommonUI::m_lNeededTI;
        if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"true"<<std::endl;
        return true;
    }
    else
    {
        if(getUI(pSeq)->EpCommonUI::m_lDebug_SEQ_UILink&1) std::cout<<"false"<<std::endl;
        return false;
    }
}

// --------------------------------------------------------------------------
// Function    : pfGetTRTIFillTimes
// --------------------------------------------------------------------------
//
// Description : Serve as a delegation of the function calculateTRTIFillTimes
//               
// Return      : true , if successful
//               false, if not
// --------------------------------------------------------------------------
bool EpCommonUINS::pfGetTRTIFillTimes
    (
    MrProtocolData::MrProtData  *rMrProt,
    SeqLim                      *rSeqLim,
    MrProtocolData::SeqExpo     *rSeqExpo,
    MrMeasSrv::ISequence        *pSeq,
    long                        *plNeededTI,
    long                        *plNeededTR
    )
{ 
    // local copy of MrProt

    MrProt sMrProt(rMrProt);

    // sqeuence class EP2D
#if (defined SEQUENCE_CLASS_EP2D || defined SEQUENCE_CLASS_BASED_ON_EP2D)
    auto pSeqImpl = static_cast<Ep2d*>(pSeq->getSeq());
    return pSeqImpl->calculateTRTIFillTimes(sMrProt, *rSeqLim, *rSeqExpo, plNeededTI, plNeededTR);
#endif

#if defined SEQUENCE_CLASS_EP_SEG
    auto pSeqImpl = static_cast<Ep_seg*>(pSeq->getSeq());
    return pSeqImpl->calculateTRTIFillTimes(sMrProt, *rSeqLim, *rSeqExpo, plNeededTI, plNeededTR);
#endif
#ifdef SEQUENCE_CLASS_AslCsl   // SUPPORT_iPAT_TGSE
    return true;
#endif

}


// *****************************************************************************
//
// Name        : AdaptRefLinesPE
//
// Description : This function checks if reference lines are within the range
///              If this is not the case, it adapts the no. of ref lines
///
// *****************************************************************************
void EpCommonUINS::adaptRefLinesPE(MrUILinkBase* const pThis, bool bAlwaysSetToMax)
{
    MrProt rMrProt(pThis->prot());

    if (LINK_LONG_TYPE* pPATRefLinesPE = _search<LINK_LONG_TYPE>(pThis, MR_TAG_PAT_LINES_PE))
    {
        if ((pPATRefLinesPE != NULL) && (pPATRefLinesPE->isAvailable(0)) && (rMrProt.PAT().getlAccelFactPE() > 1))
        {
            MrProt rMrProtCopy = rMrProt.clone();

            rMrProtCopy.PAT().setlRefLinesPE(pThis->seqLimits().getRefLinesPE().getMin());

#ifdef SEQUENCE_CLASS_EP2D
            // set partial fourier factors to reorderinfo if this is an ep2d_diff sequence
            if (getUI(pThis)->isDiffusion())
            {
                Ep2d* pSeq = getSeq(pThis);

                if (pSeq)
                    pSeq->setPartialFourierToReorderInfo(rMrProt, getUI(pThis)->m_pREOInfo);
            }
#endif

            getUI(pThis)->m_pREOInfo->prepareCalculation(rMrProtCopy, pThis->seqLimits());

            int32_t lNoOfPATLinesInLowerKSpace = getUI(pThis)->m_pREOInfo->getNoOfPATLinesInLowerKSpace();
            int32_t lNoOfPATLinesInUpperKSpace = getUI(pThis)->m_pREOInfo->getNoOfPATLinesInUpperKSpace();

            int32_t lMaxNoOfPATLinesInOneSegment = 2 * std::min<int32_t>(lNoOfPATLinesInLowerKSpace, lNoOfPATLinesInUpperKSpace);

            int32_t lMaxNoOfPATLines = lMaxNoOfPATLinesInOneSegment;
            int32_t lAccelFactPE = rMrProt.PAT().getlAccelFactPE();

            if (lAccelFactPE > 2)
            {
                lMaxNoOfPATLines = lMaxNoOfPATLinesInOneSegment * lAccelFactPE;
            }

            if (rMrProt.PAT().getlRefLinesPE() > lMaxNoOfPATLines || (bAlwaysSetToMax && rMrProt.getlTOM() != SEQ::TOM_MINIMIZE_TE))
            {
                int32_t lRoundedMaxLimit = (pThis->seqLimits().getRefLinesPE().getMax() / rMrProt.PAT().getlAccelFactPE()) *rMrProt.PAT().getlAccelFactPE();
                pPATRefLinesPE->value(std::min<int32_t>(lMaxNoOfPATLines, lRoundedMaxLimit), 0);
            }
        }
    }
}


// ------------------------------------------------------------------------------
// Function    : fEPIStdSolveBoolParamConflict
// ------------------------------------------------------------------------------
//
// Description : invokes libUICtrl's standard solve handler for bool parameters
// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::fEPIStdSolveBoolParamConflict(MrUILinkSelection<bool> * const pThis, char **arg_list, const void * pMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex)
{
    return (fUICSolveBoolParamConflict(pThis, arg_list, pMem, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded));
}


// ------------------------------------------------------------------------------
// Function    : fEPIStdSolveLongParam
// ------------------------------------------------------------------------------
//
// Description : invokes libUICtrl's standard solve handler for int32_t parameters
// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::fEPIStdSolveLongParam(LINK_LONG_TYPE* const pThis, char** arg_list, const void* pMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex)
{
    return (fUICSolveLongParamConflict(pThis, arg_list, pMem, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded));
}


// ------------------------------------------------------------------------------
// Function    : fEPIStdSolveDoubleParam
// ------------------------------------------------------------------------------
//
// Description : invokes libUICtrl's standard solve handler for double parameters
// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::fEPIStdSolveDoubleParam(LINK_DOUBLE_TYPE* const pThis, char** arg_list, const void* pMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    MrProt rOrigProt(pOrigProt); MrProt rNewProt(pThis->prot());

    adaptRefLinesPE(pThis, false);
    if (pThis->tryProt((void*)pMem, pOrigProt, lIndex))  {
        return (MRI_STD_CONFIRMATION_MSG);
    }

    return (fUICSolveDoubleParamConflict(pThis, arg_list, pMem, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded));
}


// ------------------------------------------------------------------------------
// Function    : fEPIStdSolveSelection
// ------------------------------------------------------------------------------
//
// Description : invokes libUICtrl's standard solve handler for selection parameters
//               also checks if SPAIR and IR have both been selected and if this is the case then turns off one of them
// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::fEPIStdSolveSelection (LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    MrProt rOrigProt(pOrigProt); MrProt rNewProt(pThis->prot());     //

#ifdef EP2D_MS
    // Take care of slice-acceleration specific adaptions
    MrProtFacade protFacade(pThis->prot());

    if ( protFacade.isSliceAcceleration() )
    {
        unsigned uResult = solveSliceAccelerationSettings<LINK_SELECTION_TYPE>( pThis );
        if ( uResult != 0 )
            return uResult;
    }
#endif

    // 20150105 DP: in case of NOT SPAIR switch back to default FatSat Optimization Mode since other modes are not supported yet
    bool bFatSupOptNotAllowed = ((rOrigProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_ABDOMEN)
                                 || (rOrigProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_THORAX)
                                 || (rOrigProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_BREAST))
                                && (rNewProt.preparationPulses().getlFatWaterContrast() != FatWaterContrast_Spair);
#ifdef COMPILE_EP2D_DIFF
    // Abdomen mode is available in low field, to activate eddy current compensation
    if (SysProperties::isLowField())
    {
        bFatSupOptNotAllowed = ((rOrigProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_THORAX)
                                || (rOrigProt.preparationPulses().getlFatSupOpt() == MrProtocolData::FATSUPOPT_BREAST))
                               && (rNewProt.preparationPulses().getlFatWaterContrast() != FatWaterContrast_Spair);
    }
#endif

    if (bFatSupOptNotAllowed)
    {
        MrUILinkSelection<unsigned>* _FatSatOpt
            = _search<MrUILinkSelection<unsigned>>(pThis, MR_TAG_FAT_SUP_OPTIMIZATION);
        if (!_FatSatOpt || !_FatSatOpt->isEditable(0)
            || !pThis->seqLimits().getFatSupOpt().hasOption(MrProtocolData::FATSUPOPT_DEFAULT))
        {
            return 0;
        }
        if (_FatSatOpt->value(MRI_STD_STRING_DEFAULT, 0L) != MRI_STD_STRING_DEFAULT)
        {
            return 0;
        }
        //  calls prepareForBinarySearch by default
        if (pThis->tryProt((void*)pVoid, pOrigProt, lIndex))
        {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }


#ifdef EP2D_MS
    // In case of NOT FatSat switch back to default FatSat Optimization Mode since other modes are not supported yet
    if ( ( rOrigProt.preparationPulses().getlFatSupOpt()       == MrProtocolData::FATSUPOPT_BRAIN ) && 
         ( rNewProt.preparationPulses().getlFatWaterContrast() != FatWaterContrast_FatSaturation  ) )
    {
        MrUILinkSelection<unsigned> *_FatSatOpt	= _search<  MrUILinkSelection<unsigned>>(pThis, MR_TAG_FAT_SUP_OPTIMIZATION);
        if(!_FatSatOpt || !_FatSatOpt->isEditable(0) || !pThis->seqLimits().getFatSupOpt().hasOption(MrProtocolData::FATSUPOPT_DEFAULT)) { return 0; }
        if(_FatSatOpt->value(MRI_STD_STRING_DEFAULT, 0L) != MRI_STD_STRING_DEFAULT) { return 0; }
        //  calls prepareForBinarySearch by default
        if(pThis->tryProt((void*)pVoid, pOrigProt, lIndex))  {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }
#endif

    // IR is on and the user wished to select SPAIR.  The IR should be turned off
    if((rNewProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_Spair) && !(rOrigProt.preparationPulses().getucInversion()==SEQ::INVERSION_OFF))
    {

        LINK_SELECTION_TYPE* pInversionMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_INVERSION);

        // * Check whether the UI parameter "inversion" exists, is editable and has the option "off" *
        if(!pInversionMode || !pInversionMode->isEditable(0) || !pThis->seqLimits().getInversion().hasOption(SEQ::INVERSION_OFF))  { return (0); }

        // * Switch off inversion *
        if(pInversionMode->value(MRI_STD_MAGN_PREP_NONE, 0) != MRI_STD_MAGN_PREP_NONE)  { return (0); }

        //  calls prepareForBinarySearch by default
        if(pThis->tryProt((void*)pVoid, pOrigProt, lIndex))  {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }
    else if((rOrigProt.preparationPulses().getlFatWaterContrast() == FatWaterContrast_Spair) && !(rNewProt.preparationPulses().getucInversion()==SEQ::INVERSION_OFF))
    {
        // In this option SPAIR is on and the user wishes to select IR

        LINK_SELECTION_TYPE* pFatWaterContrast = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_FAT_WATER_CONTRAST);


        // * Check whether the UI parameter "fat water contrast" exists, is editable and has the option "standard" *
        if(!pFatWaterContrast || !pFatWaterContrast->isEditable(0) || !pThis->seqLimits().getFatWaterContrast().hasOption(FatWaterContrast_Standard))  { return (0); }

        // * Switch off fat suppression *
        if(pFatWaterContrast->value(MRI_STD_FATWATERCONTRAST_STANDARD, 0) != MRI_STD_FATWATERCONTRAST_STANDARD)  { return (0); }

        //  calls prepareForBinarySearch by default
        if(pThis->tryProt((void*)pVoid, pOrigProt, lIndex))  {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }

#ifdef EP2D_MS
    // The user has slice series mode different from INTERLEAVED but wants to select IR.  The slice series mode should be set to INTERLEAVED
    if ( ( rOrigProt.sliceSeries().mode()                != SEQ::INTERLEAVED       ) && 
         ( rOrigProt.sliceSeries().mode()                != SEQ::INTERLEAVED_IN_BH ) && 
         ( rNewProt.preparationPulses().getucInversion() != SEQ::INVERSION_OFF     ) )
    {
        setProtocolParameter<LINK_SELECTION_TYPE>( pThis, MR_TAG_SERIES_MODE, MRI_STD_SERIES_INTERLEAVED );
        if ( pThis->sequence().prepareForBinarySearch( &pThis->prot() ) )
        {
            return MRI_STD_CONFIRMATION_MSG;
        }
    }
#endif

    adaptRefLinesPE(pThis, false);
    if (pThis->tryProt((void*)pVoid, pOrigProt, lIndex))  {
        return (MRI_STD_CONFIRMATION_MSG);
    }
#ifdef SUPPORT_IIR
    // SeqLoopIIR does not properly support the needed TR/TI mechanism => adapt solution from a_tse_UI.cpp
    fUISearchTR( pThis, arg_list, pVoid, pOrigProt, lIndex );
    if ( pThis->tryProt( ( void* ) pVoid, pOrigProt, lIndex ) )
    {
        return ( MRI_STD_CONFIRMATION_MSG );
    }
#endif

    return (fUICSolveSelectionConflict(pThis, arg_list, pVoid, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded));
}



#ifdef SUPPORT_iPAT_a_ep2d
// ------------------------------------------------------------------------------
// Function    : fEPIPATModeSolveSelection
// ------------------------------------------------------------------------------
//
// Description : invokes libUICtrl's standard solve handler for selction parameters
// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// CHARM 315966: This PAT-specific selection solve-handler was originally not
//               correctly copied from VA21B to VA25A archive.
//               Call to fUICSolveSelectionConflict() mdified by replacing
//               pointers to functions giving "needed" TE, TI, TR with
//               NULL pointers. This forces fUICSolveSelectionConflict() to
//               solve conflicts by executing a binary search. This avoids the
//               problem that these "needed" functions are not correctly defined
//               when the diffusion SBB fails to prepare at short TE's. One
//               consequence of this was that PAT mode could not be set to none
//               once a minimum TE had been selected in GRAPPA or SENSE mode.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::fEPIPATModeSolveSelection (LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    MrProtFacade protFacade(pThis->prot());

    if(protFacade.isSliceAcceleration())
    {        
        unsigned uResult = solveSliceAccelerationSettings<LINK_SELECTION_TYPE>(pThis);
        if(uResult != 0)
            return uResult;
    }

    adaptRefLinesPE(pThis, false);
    if (pThis->tryProt((void*)pVoid, pOrigProt, lIndex))  {
        return (MRI_STD_CONFIRMATION_MSG);
    }

    // All the standard PAT related issues are handled by the standard UI software.
    // For EPI we have to add one case: If TE=min, we might not be able to decrease the PAT factor.
    // Therefore we try to solve tetitr here.
    if(fUICSolveSelectionConflict(pThis, arg_list, NULL, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded ))
    {
        if(arg_list)
        {
            fUICFormatTETRTIChanges(pThis, arg_list, NULL, pOrigProt, 0);
        }
    }

    if(pThis->sequence().prepareForBinarySearch(&pThis->prot()))
        return MRI_STD_CONFIRMATION_MSG;

    return 0;
}

// ---------------------------------------------------------------------------------
// Function    : fUILinkPATModeSetValue
// ---------------------------------------------------------------------------------
// Description : Invokes standard set handler and adapts the number of PAT reference
//               lines 
// ---------------------------------------------------------------------------------
unsigned EpCommonUINS::fUILinkPATModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    unsigned uResult = getUI(pThis)->EpCommonUI::m_PATMode.getOrigSetValueHandler()(pThis, uNewVal, lPos);

    // we only want the maximum reference lines if the refscan mode is not Fast GRE
    const bool bAlwaysSetToMax = pThis->prot().getsPat().getucRefScanMode() != SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST;
    adaptRefLinesPE(pThis, bAlwaysSetToMax);

    return uResult;
}

// ------------------------------------------------------------------------------
// Function    : fEPIRefLinesPEGetLimits
// ------------------------------------------------------------------------------
//
// Description : Enforces that the limits for PAT reference lines is a multiple
//               of the acceleration factor
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::fEPIRefLinesPEGetLimits(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t /* pos */)
{
    rulVerify = LINK_LONG_TYPE::VERIFY_BINARY_SEARCH;

    rLimitVector.clear();

    MrProtocolData::MrProtData* pProt   = &pThis->prot();
    SeqLim*                     pSeqLim = &pThis->seqLimits();

    const ParLim<int32_t>& _seqLimits = pSeqLim->getRefLinesPE();

    rulVerify = LINK_LONG_TYPE::VERIFY_SCAN_ALL;
    //  Reference lines must be an integer multiple of acceleration factor
    rLimitVector.resize(1);
    const int32_t lAccelPE = pProt->getsPat().getlAccelFactPE();

    return rLimitVector[0].setEqualSpaced(lAccelPE * ((_seqLimits.getMin()+lAccelPE-1L) / lAccelPE)
        , _seqLimits.getMax()
        , lAccelPE
        );
}

// ------------------------------------------------------------------------------
// Function    : fEPIAccPESetValue
// ------------------------------------------------------------------------------
//
// Description : Enforces that the limits for PAT reference lines is a multiple
//               of the acceleration factor
//
// ------------------------------------------------------------------------------
int32_t EpCommonUINS::fEPIAccPESetValue(LINK_LONG_TYPE* const pThis, int32_t val, int32_t pos)
{
    MrProtocolData::MrProtData*  pProt = &pThis->prot();
    const long lAccelPE_prev = pProt->getsPat().getlAccelFactPE();

    // set PAT ref scan to EPI if in-plane acceleration is deactivated
#ifndef EP2D_MS
    if (lAccelPE_prev > 1 && val == 1 && (pProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA || pProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST))
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_PAT_REF_SCAN_MODE, MRI_STD_PAT_REF_SCAN_EXTRA_EPI, 0);
#endif
    
    pProt->getsPat().setlAccelFactPE(val);
    
    // Make ref lines multiple of PE factor
    pProt->getsPat().setlRefLinesPE(val*(pProt->getsPat().getlRefLinesPE() / lAccelPE_prev));

    // set recommended FOV shift factor
    // In the current implementation the UI element for the FOV shift factor is hidden thus the factor is written to the protocol directly.
    // If the UI parameter is supposed to be used again (see fUILinkFOVShiftFactorIsAvailable) use the line which is currently commented out.
    // We need to specify if this is a diffusion sequence as we get different values for 1.5T there
    MrProt	rMrProt(&pThis->prot());

    SMSProperties::setFOVShiftFactor(rMrProt, SMSProperties::getRecommendedFOVShiftFactorCumulative(pThis->prot(), getUI(pThis)->isSpinEcho()));

    //setProtocolParameter<LINK_LONG_TYPE>(pThis, MR_TAG_SLICE_ACCEL_FOV_SHIFT_FACTOR, SMSProperties::getRecommendedFOVShiftFactorCumulative(pThis->prot()), 0, false, false);


    _addTimingDependencyPtr(pThis);

    if (rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI)
        adaptRefLinesPE(pThis, true);
    else
        adaptRefLinesPE(pThis, false);

    return pThis->value(pos);
}

// ------------------------------------------------------------------------------
// Function    : _addTimingDependencyPtr
// ------------------------------------------------------------------------------
//
// Description : adds timing dependency to setvalue handler, see MrUILink for 
//               original
//
// ------------------------------------------------------------------------------
void EpCommonUINS::_addTimingDependencyPtr(MrUILinkBase* const pThis)
{
    for(int nTiming = 0; nTiming < 3; nTiming++)
    {
        char szNameTag[64];
        switch(nTiming)
        {
        case 0:
                strcpy(szNameTag, MR_TAG_TR);
            break;
        case 1:
                strcpy(szNameTag, MR_TAG_TI);
            break;
        case 2:
                strcpy(szNameTag, MR_TAG_TE);
            break;
        }

        MrUILinkArray* pArray = _search<MrUILinkArray>(pThis, szNameTag);
        if(pArray)
        {
            LINK_DOUBLE_TYPE* pT = _searchElm<LINK_DOUBLE_TYPE>(pThis, szNameTag);

            int32_t lSize = pArray->size(0);
            for (int32_t lIndex = 0; lIndex < lSize; lIndex++)
            {
                pThis->addDependentParamPtr(pT, lIndex);
            }
        }
    }
}

#endif //#ifdef SUPPORT_iPAT_a_ep2d


// ------------------------------------------------------------------------------
// Function    : _solveBaseResolution
// ------------------------------------------------------------------------------
//
// Description : invokes libUICtrl's standard solve handler for base resolution,
//               which needs the original base resolution solve handler.

// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::_solveBaseResolution (LINK_LONG_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    //---------------------------------------------------------------------------
    // start with std. solve handler -> if it works, we like it
    //---------------------------------------------------------------------------
    const UI_ELEMENT_LONG&  rBaseResolution = getUI(pThis)->EpCommonUI::m_BaseResolution;

    //Before calling the standard solve handler, adapt the number of ref lines!
    adaptRefLinesPE(pThis, false);
    if (pThis->tryProt((void*)pVoid, pOrigProt, lIndex))  {
        return (MRI_STD_CONFIRMATION_MSG);
    }

    unsigned bResult =  fUICSolveBaseResolutionConflict(pThis, arg_list, pVoid, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded, rBaseResolution.getOrigSolveHandler());

    return bResult;
}


// ------------------------------------------------------------------------------
// Function    : _epi_BaseResolutionTry
// ------------------------------------------------------------------------------
//
// Description : resets getUI(pSeq)->m_bNeedOtherTETITR and calls original try-handler for base-
//               resolution
// Return      : whatever original try-handler says
//
// ------------------------------------------------------------------------------
static bool EpCommonUINS::_epi_BaseResolutionTry (LINK_LONG_TYPE* const pThis, void* pVoid, const MrProtocolData::MrProtData*  pOrig, int32_t lIndex)
{
    getUI(pThis)->EpCommonUI::m_bNeedOtherTETITR = false;

    const UI_ELEMENT_LONG&  rBaseResolution = getUI(pThis)->EpCommonUI::m_BaseResolution;

    if(rBaseResolution.getOrigTryHandler())
    {
        return((*(rBaseResolution.getOrigTryHandler()))(pThis, pVoid, pOrig, lIndex));
    }
    return false;
}


// ------------------------------------------------------------------------------
// Function    : _solveTETITR_TE
// ------------------------------------------------------------------------------
//
// Description : solves conflicts TE vs. TR/TI
// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::_solveTETITR_TE (LINK_DOUBLE_TYPE* const pThis, char** arg_list, const void*, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    if(pThis->prot().getalTE()[0]!=getUI(pThis)->EpCommonUI::m_lNeededTE)
    {
        // can't solve a conflict with myself
        //
        return 0;
    }

    unsigned _ret = fUICSolveDoubleParamConflict(pThis, arg_list, NULL, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded);

    // reformat TE-changes (we don't want to see the changed TE)
    //
    if(arg_list)
    {
        fUICFormatTETRTIChanges(pThis, arg_list, NULL, pOrigProt, 0, 0);
    }

    return _ret;
}


// ------------------------------------------------------------------------------
// Function    : _solveTETITR_PPF_EPIFactorConflict
// ------------------------------------------------------------------------------
//
// Description : Small EPI factors have restrictions concerning phase partial
//               fourier factors. So this solves handler tries to switch off
//               phase partial fourier, if the desired EPI factor could not be
//               selected.
//
//               An additional restriction for segmented EPI sequences in A15A
//               exists. EPI-factor 3 is not supported with phase partial fourier.
//
// Return      : MRI_STD_CONFIRMATION_MSG (if success) or 0
//
// ------------------------------------------------------------------------------
static unsigned EpCommonUINS::_solveTETITR_PPF_EPIFactorConflict (MrUILinkLimited<int32_t> * const pThis, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    bool      bChangedPPF    = false;
    bool      bChangedTETITR = false;
    unsigned  _ret           = 0;



    if((pThis->prot().getsKSpace().getucPhasePartialFourier() != getUI(pThis)->EpCommonUI::m_pREOInfo->getPhasePartialFourierFactor_x_8())
        || (pThis->prot().getsFastImaging().getlEPIFactor()==3 && pThis->prot().getsKSpace().getucPhasePartialFourier() !=SEQ::PF_OFF)
        )
    {
        pThis->prot().getsKSpace().setucPhasePartialFourier(SEQ::PF_OFF);

        _ret = pThis->sequence().prepareForBinarySearch(&pThis->prot());

        if(!_ret)
        {
            _ret = fUICSolveLongParamConflict(pThis, arg_list, pVoid, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded);

            if(_ret) { bChangedTETITR = true; }
        }
        else
        {
            (void)fUICSolveLongParamConflict(pThis, arg_list, pVoid, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded);
        }

        if(_ret) { bChangedPPF = true; }
    }
    else
    {
        _ret = fUICSolveLongParamConflict(pThis, arg_list, pVoid, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded);

        if(_ret)
        {
            bChangedTETITR = true;
        }
    }

    if(_ret && arg_list && bChangedPPF)
    {
        arg_list[20] = (char*)MRI_STD_PHASE_PARTIAL_FOURIER_LABEL;

        switch(pOrigProt->getsKSpace().getucPhasePartialFourier())
        {
        default:
            case   SEQ::PF_OFF:  arg_list[24]  = (char*)MRI_STD_PARTIAL_FOURIER_OFF; break;
            case   SEQ::PF_7_8:  arg_list[24]  = (char*)MRI_STD_PARTIAL_FOURIER_7_8; break;
            case   SEQ::PF_6_8:  arg_list[24]  = (char*)MRI_STD_PARTIAL_FOURIER_6_8; break;
            case   SEQ::PF_5_8:  arg_list[24]  = (char*)MRI_STD_PARTIAL_FOURIER_5_8; break;
        case   SEQ::PF_HALF:  arg_list[24]  = (char*)MRI_STD_PARTIAL_FOURIER_HALF; break;
        }

        arg_list[30]  = (char*)MRI_STD_PARTIAL_FOURIER_OFF;
        arg_list[36]  = (char*)MRI_STD_EMPTY;
        arg_list[39]  = (char*)(bChangedTETITR ? '\n' : '\0');

        fUICFormatTETRTIChanges(pThis, arg_list, NULL, pOrigProt, 1);
    }

    if(_ret) return MRI_STD_CONFIRMATION_MSG;
    else      return 0;
}


// ------------------------------------------------------------------------------
// Function    : _epi_EchoSpacingIsAvailable
// ------------------------------------------------------------------------------
//
// Description : determines, if parameter echo spacing can be seen on the UI
// Return      : true or false
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::_epi_EchoSpacingIsAvailable
    (
    LINK_DOUBLE_TYPE* const /*pThis*/,
 int32_t /*lIndex*/
    )
{
    return true;
}


// ------------------------------------------------------------------------------
// Function    : _epi_EchoSpacingGetValue
// ------------------------------------------------------------------------------
//
// Description : determines value of echo spacing to be shown on the UI
// Return      : echo sapcing [ms]
//
// ------------------------------------------------------------------------------
double EpCommonUINS::_epi_EchoSpacingGetValue (LINK_DOUBLE_TYPE* const pThis, int32_t /*lIndex*/)
{
    double dVal = 0.0;

    if(! pThis->prot().getsFastImaging().getucFreeEchoSpacing())
    {
        if(pThis->sequence().prepareForBinarySearch(&(pThis->prot())))
        {
            dVal = pThis->sequence().getSEQ_BUFF()->getSeqExpo().getEchoSpacing()/1000.0;
        }
    }
    else
    {
        dVal = pThis->prot().getsFastImaging().getlEchoSpacing()/1000.0;
    }

    char _buffer[32];
    if(dVal < 10.0) sprintf(_buffer, "%.2f", dVal);
    else             sprintf(_buffer, "%.1f", dVal);

    return atof(_buffer);
}


// ------------------------------------------------------------------------------
// Function    : _epi_EchoSpacingSetValue
// ------------------------------------------------------------------------------
//
// Description : Called when user modifies echo spacing in the UI. Puts appropriate
//               value into the protocol.
// Return      : new value of echo spacing
//
// ------------------------------------------------------------------------------
double EpCommonUINS::_epi_EchoSpacingSetValue (LINK_DOUBLE_TYPE* const pThis, double dNewVal_ms, int32_t lIndex)
{
    pThis->prot().getsFastImaging().setucFreeEchoSpacing(true);
    pThis->prot().getsFastImaging().setlEchoSpacing(static_cast<int32_t>(0.5+1000*dNewVal_ms));
    return pThis->value(lIndex);
}


// ------------------------------------------------------------------------------
// Function    : _epi_EchoSpacingGetLimits
// ------------------------------------------------------------------------------
//
// Description : Calculates limits for echo spacing. Needed, because the echo
//               spacing has certain forbidden ranges due to resonances of the gradient
//               coil.
// Return      : false, if echospacing may not be edited
//                      or if no valid limits can be found
//               true, if echo spacing has valid limits
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::_epi_EchoSpacingGetLimits (LINK_DOUBLE_TYPE* const pThis, std::vector< MrLimitDouble >& rLimitVector, uint32_t& rulVerify, int32_t /*lIndex*/)
{
    if(! pThis->prot().getsFastImaging().getucFreeEchoSpacing()) return false; //  read-only in this case

    if(! pThis->sequence().prepareForBinarySearch(&(pThis->prot())))
    {
        return false;
    }

    getUI(pThis)->EpCommonUI::m_pEPIRO->setIgnoreForbiddenEchoSpacingRange(true);
    getUI(pThis)->EpCommonUI::m_pEPIRO->setUseShortestEchoSpacing();

    MrProt sMrProtWrapper(&pThis->prot());

    if(! getUI(pThis)->EpCommonUI::m_pEPIRO->calculateTiming(sMrProtWrapper, pThis->seqLimits()))
    {
        return false;
    }

    const double dMinEchoSpacing_ms               = getUI(pThis)->EpCommonUI::m_pEPIRO->getEchoSpacing()/1000.;
    const double dMaxEchoSpacing_ms               = 3.0;
    const double dIncrEchoSpacing_ms              = 0.01;
    const unsigned int lForbiddenRangesSize       = getUI(pThis)->EpCommonUI::m_pEPIRO->getForbiddenEchoSpacingRangeSize();

    MrLimitDouble sInterval;


    // activate binary search
    //
    rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;

    // clear limits - we have to set this up on our own
    //






    rLimitVector.clear();

    // Our task: identify the allowable echo spacing intervals
    //  min                                                           max
    //   |                                                             |
    //          ]XXXXXXXXXXXXXXXXX[   ...      ]XXXXXXXXXXXXXXXXX[
    //           forbidden range 1              forbidden range N


    // Are there forbidden ranges at all?
    //
    if(lForbiddenRangesSize == 0)
    {
        //   |                              |
        //
        //  min                            max
        if(sInterval.setEqualSpaced(dMinEchoSpacing_ms, dMaxEchoSpacing_ms, dIncrEchoSpacing_ms))
        {
            rLimitVector.push_back(sInterval);
        }

        // SEQ_TRACE_INFO.print("No forbidden ranges - Interval from %7f to %7f ms", dMinEchoSpacing_ms, dMaxEchoSpacing_ms);

        return !rLimitVector.empty();
    }

    double dCurrentEchoSpacing_ms = 0.;         // Denotes the lower border of the currently considered interval
    double dNextEchoSpacing_ms    = 0.;         // Denotes the upper border of the currently considered interval
    bool bForbiddenRange          = false;      // Indicates wether the currently considered interval is forbidden or not
    bool bFinish                  = false;      // Indicates that the iteration ran up to the maximum echo spacing
    unsigned int iI               = 0;          // Counter

    // Search for the first forbidden range that affects allowable echo spacings
    // Note: forbidden ranges are stored in pEPIRO in increasing order
    //
    while((iI <= lForbiddenRangesSize) && !bFinish)
    {
        if(dMinEchoSpacing_ms >= getUI(pThis)->EpCommonUI::m_pEPIRO->getForbiddenEchoSpacingRangeMax(iI) / 1000.)
        {
            //                                         min          max
            //                                          |            |
            //  ...              ]XXXXXXXXXXXXXXXXX[   ...
            //                    forbidden range m

            // Don't care about current forbidden range - advance to next one
            ++iI;
        }
        else if(dMinEchoSpacing_ms > getUI(pThis)->EpCommonUI::m_pEPIRO->getForbiddenEchoSpacingRangeMin(iI) / 1000.)
        {
            //                           min        max
            //                            |     |---------|
            //  ...              ]XXXXXXXXXXXXXXXXX[   ...
            //                    forbidden range m
            dCurrentEchoSpacing_ms = dMinEchoSpacing_ms;
            dNextEchoSpacing_ms    = std::min(dMaxEchoSpacing_ms, getUI(pThis)->EpCommonUI::m_pEPIRO->getForbiddenEchoSpacingRangeMax(iI) / 1000.);
            bForbiddenRange = true;     // Inside a forbidden range
            bFinish         = true;     // First relevant range identified
        }
        else
        {
            //              min             max
            //               |  |-------------------------|
            //  ...              ]XXXXXXXXXXXXXXXXX[   ...
            //                    forbidden range m
            dCurrentEchoSpacing_ms = dMinEchoSpacing_ms;
            dNextEchoSpacing_ms    = std::min(dMaxEchoSpacing_ms, getUI(pThis)->EpCommonUI::m_pEPIRO->getForbiddenEchoSpacingRangeMin(iI) / 1000.);
            bForbiddenRange = false;    // Inside an allowed range
            bFinish         = true;     // First relevant range identified
        }
    }

    if(!bFinish)
    {
        //                                         min          max
        //                                          |            |
        //  ...              ]XXXXXXXXXXXXXXXXX[
        //                    forbidden range N

        // There is no forbidden range that affects allowable echo spacings
        if(sInterval.setEqualSpaced(dMinEchoSpacing_ms, dMaxEchoSpacing_ms, dIncrEchoSpacing_ms))
        {
            rLimitVector.push_back(sInterval);
        }

        return !rLimitVector.empty();
    }


    bFinish = false;

    // Iterate over all forbidden ranges
    while((iI <= lForbiddenRangesSize) && !bFinish)
    {
        // SEQ_TRACE_INFO.print("Iteration %ld: Current %7f Next %7f Forbidden %ld", iI, dCurrentEchoSpacing_ms, dNextEchoSpacing_ms, bForbiddenRange);

        if(bForbiddenRange)
        {
            // Currently, we are inside a forbidden range

            //                           cur  --> next
            //                            |        |
            //  ...              ]XXXXXXXXXXXXXXXXX[   ...
            //                    forbidden range m

            // Advance
            ++iI;
            bForbiddenRange        = false;    // Next interval is inside an allowed range
            dCurrentEchoSpacing_ms = dNextEchoSpacing_ms;

            if(iI < lForbiddenRangesSize)
            {
                dNextEchoSpacing_ms = std::min(dMaxEchoSpacing_ms, getUI(pThis)->EpCommonUI::m_pEPIRO->getForbiddenEchoSpacingRangeMin(iI) / 1000.);
            }
            else
            {
                dNextEchoSpacing_ms = dMaxEchoSpacing_ms;
            }
        }
        else
        {
            // Currently, we are inside an allowed range

            //         cur  --> next
            //          |        |
            //  ...              ]XXXXXXXXXXXXXXXXX[   ...
            //                    forbidden range m

            // Push the allowed interval
            // SEQ_TRACE_INFO.print("Allowed interval from %7f to %7f ms", dCurrentEchoSpacing_ms, dNextEchoSpacing_ms);

            if(sInterval.setEqualSpaced(dCurrentEchoSpacing_ms, dNextEchoSpacing_ms, dIncrEchoSpacing_ms))
            {
                rLimitVector.push_back(sInterval);
            }

            // Advance - no need to increment iI here
            // (first, the upper limit of the current forbidden range has to be considered)
            bForbiddenRange        = true; // Next interval is inside a forbidden range
            dCurrentEchoSpacing_ms = dNextEchoSpacing_ms;
            dNextEchoSpacing_ms    = std::min(dMaxEchoSpacing_ms, getUI(pThis)->EpCommonUI::m_pEPIRO->getForbiddenEchoSpacingRangeMax(iI) / 1000.);
        }

        if(dCurrentEchoSpacing_ms == dMaxEchoSpacing_ms)
        {
            bFinish = true;
        }
    }

    // SEQ_TRACE_INFO.print("Result: %ld allowed intervals", rLimitVector.size());

    return !rLimitVector.empty();
}


// ------------------------------------------------------------------------------
// Function    : _epi_EchoSpacing_GetToolTip
// ------------------------------------------------------------------------------
//
// Description : Tooltip for echo spapcing. Shows the user the fat/water-shift in
//               PE direction caused by phases aquired during echo train.
// Return      : 0, if no tooltip should/can be displayed,
//               MRI_STD_BW_TOOLTIP, else i.e. Fat-water-shift:\t%1!r!\t[Px]
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::_epi_EchoSpacing_GetToolTip(LINK_DOUBLE_TYPE* const pThis, char* arg_list[], int32_t)
{
    if(!pThis->sequence().prepareForBinarySearch(&(pThis->prot())))
    {
        // can't continue, should never happen
        return 0;
    }

    if(getUI(pThis)->EpCommonUI::m_pREOInfo->getEchoTrainLength()<2)
    {
        // no shift in PE, if no multi-echo measurement
        return 0;
    }

    const double bandWidthPerPixelPE = pThis->sequence().getSEQ_BUFF()->getSeqExpo().getBandwidthPerPixelPhaseEncode();

    if (bandWidthPerPixelPE < std::numeric_limits<double>::epsilon())
    {
        // bandwidth per pixel in PE direction <= 0, should never happen. avoid division by 0.
        return 0;
    }

    const MeasNucleus myNucleus(pThis->prot().getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str());
    const double      dGAMMA     = myNucleus.getLarmorConst();
    const double      dDeltaFreq = fabs(dGAMMA * SysProperties::getNominalB0() * CHEMICAL_SHIFT_FAT_PPM); // Hz

    const double dShift_Px = dDeltaFreq / bandWidthPerPixelPE;

    sprintf(getUI(pThis)->m_tLine_epi_EchoSpacing_GetToolTip, "%.2f", dShift_Px);

    arg_list[0] = (char*)MRI_STD_STRING; // %1!s!
    arg_list[1] = getUI(pThis)->EpCommonUI::m_tLine_epi_EchoSpacing_GetToolTip;

    return MRI_STD_BW_TOOLTIP;
}


// ------------------------------------------------------------------------------
// Function    : _epi_FreeEchoSpacing_GetLabelId
// ------------------------------------------------------------------------------
//
// Description : determines NLS supported text shown on the UI
// Return      : ressource-ID
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::_epi_FreeEchoSpacing_GetLabelId(LINK_BOOL_TYPE* const, char* /*arg_list*/[], int32_t /*lIndex*/)
{
    return MRI_STD_FREE_ECHO_SPACING;
}


// ------------------------------------------------------------------------------
// Function    : _epi_FreeEchoSpacing_GetOptions
// ------------------------------------------------------------------------------
//
// Description : determines options the user may select
// Return      : true
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::_epi_FreeEchoSpacing_GetOptions(LINK_BOOL_TYPE* const /*pThis*/, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t /*lIndex*/)
{
    rulVerify = LINK_BOOL_TYPE::VERIFY_ON;
    rOptionVector.resize(2);
    rOptionVector[0] = false;
    rOptionVector[1] = true;
    return true;
}


// ------------------------------------------------------------------------------
// Function    : _epi_FreeEchoSpacing_GetValue
// ------------------------------------------------------------------------------
//
// Description : determines current setting of free echo spacing
// Return      : true or false
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::_epi_FreeEchoSpacing_GetValue(LINK_BOOL_TYPE* const pThis, int32_t /*lIndex*/)
{
    return pThis->prot().getsFastImaging().getucFreeEchoSpacing() != 0;
}


// ------------------------------------------------------------------------------
// Function    : _epi_FreeEchoSpacing_SetValue
// ------------------------------------------------------------------------------
//
// Description : Activates or deactivates free echo spacing. If activated the
//               value for echo spacing parameter is initialized with the current
//               possible minimum.
//
// Return      : new setting of free echo spacing
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::_epi_FreeEchoSpacing_SetValue(LINK_BOOL_TYPE* const pThis, bool value, int32_t /*lIndex*/)
{
    if(value)
    {
        getUI(pThis)->EpCommonUI::m_pEPIRO->setIgnoreForbiddenEchoSpacingRange(false);
        getUI(pThis)->EpCommonUI::m_pEPIRO->setUseShortestEchoSpacing();

        MrProt sMrProtWrapper(&pThis->prot());

        if(! getUI(pThis)->EpCommonUI::m_pEPIRO->calculateTiming(sMrProtWrapper, pThis->seqLimits()))
        {
            pThis->prot().getsFastImaging().setlEchoSpacing(0);
        }
        else
        {
            pThis->prot().getsFastImaging().setlEchoSpacing(getUI(pThis)->EpCommonUI::m_pEPIRO->getEchoSpacing());
        }
    }
    else
    {
        pThis->prot().getsFastImaging().setlEchoSpacing(0);
    }

    pThis->prot().getsFastImaging().setucFreeEchoSpacing(value);

    return pThis->prot().getsFastImaging().getucFreeEchoSpacing() != 0;
}


// ------------------------------------------------------------------------------
// Function    : _epi_FreeEchoSpacing_Solve
// ------------------------------------------------------------------------------
//
// Description : free echo spacing solve handler
// Return      : 0, i.e. we can't help
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::_epi_FreeEchoSpacing_Solve(LINK_BOOL_TYPE* const /*pThis*/, char* /*arg_list*/[], const void* /*pVoid*/, const MrProtocolData::MrProtData *  /*rMrProt*/, int32_t /*lIndex*/)
{
    return 0;
}


// ------------------------------------------------------------------------------
// Function    : _TIGetLimits
// ------------------------------------------------------------------------------
//
// Description : Invokes libUICtrl's standard get-limit handler for TI depending
//               on help of the sequence itself supplied with call-back function
//               m_pfCalculateTRTIFill.
// Return      : true (if valid limits can be found) or false
//
// ------------------------------------------------------------------------------
static bool EpCommonUINS::_TIGetLimits (LINK_DOUBLE_TYPE * const pThis, std::vector< MrLimitDouble >& rLimitVector, uint32_t& rulVerify, int32_t lIndex)
{
    const UI_ELEMENT_DOUBLE&  rTI = getUI(pThis)->EpCommonUI::m_TI;

    return fUICGetTILimits(pThis, rLimitVector, rulVerify, lIndex, rTI.getOrigGetLimitsHandler(), getUI(pThis)->EpCommonUI::m_pfCalculateTRTIFill);
}


// ------------------------------------------------------------------------------
// Function    : _solveTETITR_TI
// ------------------------------------------------------------------------------
//
// Description : solves conflicts TI vs. TR/TE
// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// ------------------------------------------------------------------------------
static unsigned EpCommonUINS::_solveTETITR_TI (MrUILinkLimited<double>* const pThis, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex)
{
    if(getUI(pThis)->EpCommonUI::m_lNeededTI != pThis->prot().getalTI()[0])
    {
        // can't solve a conflict with myself
        //
        return 0;
    }

    unsigned _ret = fUICSolveDoubleParamConflict(pThis, arg_list, pVoid, pOrigProt, lIndex, &_getTENeeded, &_getTINeeded, &_getTRNeeded);

    // we don't want to see the changed TI in the list of adapted parameters
    //
    if(arg_list)
    {
        MrProt MyOrigProt(pOrigProt->clone());
        MyOrigProt.getalTI()[0]=pThis->prot().getalTI()[0];

        fUICFormatTETRTIChanges(pThis, arg_list, pVoid, MyOrigProt, 0);
    }

    return _ret;
}


// ------------------------------------------------------------------------------
// Function    : fEPIGradientModeGetOptions()
// ------------------------------------------------------------------------------
//
// Overloaded get-options handler for "Gradient mode".
//
// Removes the restriction in the standard function fUILinkGradientModeGetOptions()
// that the user can't switch between normal and fast modes if the corresponding
// MeasPerm specifications are the same (as for example with Z-Engine).
//
// By removing this restriction it is possible for the EPI sequence to set
// it's own specification for "normal", which is used to offer a slightly reduced
// slew rate for application in the under-voltage case. This reduces problems with
// high-spatial-frequency Nyquist ghosts caused by imperfect k-space trajectories
// during readout gradient ramps (CHARM 333764).
//
// The "normal" and "fast" behaviour is set in fSEQPrep() in a_ep2d.cpp.
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::fEPIGradientModeGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t /*pos*/)
{
    rulVerify = LINK_SELECTION_TYPE::VERIFY_ON;
    rOptionVector.clear();

    const ParLimOption<SEQ::Gradients>& _seqLimits = pThis->seqLimits().getGradients();

    bool bGradOptionDetected = false;

    if (_seqLimits.hasOption(SEQ::GRAD_WHISPER))
    {
        bGradOptionDetected = true;
        rOptionVector.push_back(MRI_STD_GRADIENT_MODE_WHISPER);
    }

    if (_seqLimits.hasOption(SEQ::GRAD_NORMAL))
    {
        bGradOptionDetected = true;
        rOptionVector.push_back(MRI_STD_GRADIENT_MODE_NORMAL);
    }

    if (_seqLimits.hasOption(SEQ::GRAD_FAST))
    {
        bGradOptionDetected = true;
        rOptionVector.push_back(MRI_STD_GRADIENT_MODE_FAST);
    }

    if (_seqLimits.hasOption(SEQ::GRAD_ULTRAFAST))
    {
        bGradOptionDetected = true;
        rOptionVector.push_back(MRI_STD_GRADIENT_MODE_ULTRAFAST);
    }

    if (!bGradOptionDetected)
    {
        SEQ_TRACE_INFO.print("No gradient options found");
    }

    //  Note: The user can never select one of the GSWD-modes explicitly.
    switch (pThis->prot().getsGRADSPEC().getucMode())
    {
    case SEQ::GRAD_ULTRAFAST_GSWD_RISETIME:
        rOptionVector.push_back(MRI_STD_STRING);
        SET_MODIFIER(rOptionVector.back(), MRI_STD_GRADIENT_MODE_ULTRAFAST);
        break;

    case SEQ::GRAD_FAST_GSWD_RISETIME:
        rOptionVector.push_back(MRI_STD_STRING);
        SET_MODIFIER(rOptionVector.back(), MRI_STD_GRADIENT_MODE_FAST);
        break;

    case SEQ::GRAD_NORMAL_GSWD_RISETIME:
        rOptionVector.push_back(MRI_STD_STRING);
        SET_MODIFIER(rOptionVector.back(), MRI_STD_GRADIENT_MODE_NORMAL);
        break;

    case SEQ::GRAD_WHISPER_GSWD_RISETIME:
        rOptionVector.push_back(MRI_STD_STRING);
        SET_MODIFIER(rOptionVector.back(), MRI_STD_GRADIENT_MODE_WHISPER);
        break;

    case SEQ::GRAD_FAST:                // FALL-THROUGH
    case SEQ::GRAD_NORMAL:              // FALL-THROUGH
    case SEQ::GRAD_WHISPER:             // FALL-THROUGH
    case SEQ::GRAD_ULTRAFAST:           // FALL-THROUGH
    case SEQ::GRAD_BOOST:               // FALL-THROUGH
    case SEQ::GRAD_BOOST_GSWD_RISETIME: // FALL-THROUGH
        break;                          // nothing to do

    default:
        SEQ_TRACE_INFO.print("Gradient mode '%04X' not recognised", pThis->prot().getsGRADSPEC().getucMode());
    }

    return !rOptionVector.empty();
}


// ------------------------------------------------------------------------------
// Function    : fEPIGradientModeGetValue()
// ------------------------------------------------------------------------------
//
// Overloaded get-value handler for "Gradient mode".
//
// Removes the restriction in the standard function fUILinkGradientModeGetValue()
// that the user can't switch between normal and fast modes if the corresponding
// MeasPerm specifications are the same (as for example with Z-Engine).
//
// By removing this restriction it is possible for the EPI sequence to set
// it's own specification for "normal", which is used to offer a slightly reduced
// slew rate for application in the under-voltage case. This reduces problems with
// high-spatial-frequency Nyquist ghosts caused by imperfect k-space trajectories
// during readout gradient ramps (CHARM 333764).
//
// The "normal" and "fast" behaviour is set in fSEQPrep() in a_ep2d.cpp.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::fEPIGradientModeGetValue(LINK_SELECTION_TYPE* const pThis, int32_t /*pos*/)
{
    unsigned nVal = 0;

    switch (pThis->prot().getsGRADSPEC().getucMode())
    {
    case SEQ::GRAD_ULTRAFAST:
        nVal = MRI_STD_GRADIENT_MODE_ULTRAFAST;
        break;

    case SEQ::GRAD_ULTRAFAST_GSWD_RISETIME:
        nVal = MRI_STD_STRING;
        SET_MODIFIER(nVal, MRI_STD_GRADIENT_MODE_ULTRAFAST);
        break;

    case SEQ::GRAD_FAST:
        nVal = MRI_STD_GRADIENT_MODE_FAST;
        break;

    case SEQ::GRAD_FAST_GSWD_RISETIME:
        nVal = MRI_STD_STRING;
        SET_MODIFIER(nVal, MRI_STD_GRADIENT_MODE_FAST);
        break;

    case SEQ::GRAD_NORMAL:
        nVal = MRI_STD_GRADIENT_MODE_NORMAL;
        break;

    case SEQ::GRAD_NORMAL_GSWD_RISETIME:
        nVal = MRI_STD_STRING;
        SET_MODIFIER(nVal, MRI_STD_GRADIENT_MODE_NORMAL);
        break;

    case SEQ::GRAD_WHISPER:
        nVal = MRI_STD_GRADIENT_MODE_WHISPER;
        break;

    case SEQ::GRAD_WHISPER_GSWD_RISETIME:
        nVal = MRI_STD_STRING;
        SET_MODIFIER(nVal, MRI_STD_GRADIENT_MODE_WHISPER);
        break;

    case SEQ::GRAD_BOOST:               // FALL-THROUGH
    case SEQ::GRAD_BOOST_GSWD_RISETIME: // FALL-THROUGH
        break;                          // nothing to do

    default:
        SEQ_TRACE_INFO.print("Gradient mode '%04X' not recognised", pThis->prot().getsGRADSPEC().getucMode());
        nVal = MRI_STD_INVALID;
    }

    return nVal;
}

// ------------------------------------------------------------------------------
// Function    : fEPIGradientModeSetValue
// ------------------------------------------------------------------------------
//
// Overloaded set-value handler for "Gradient mode".
//
// Removes the restriction in the standard function fUILinkGradientModeSetValue()
// that the user can't switch between normal and fast modes if the corresponding
// MeasPerm specifications are the same (as for example with Z-Engine).
//
// By removing this restriction it is possible for the EPI sequence to set
// it's own specification for "normal", which is used to offer a slightly reduced
// slew rate for application in the under-voltage case. This reduces problems with
// high-spatial-frequency Nyquist ghosts caused by imperfect k-space trajectories
// during readout gradient ramps (CHARM 333764).
//
// The "normal" and "fast" behaviour is set in fSEQPrep() in a_ep2d.cpp.
//
// ------------------------------------------------------------------------------
unsigned EpCommonUINS::fEPIGradientModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned newVal, int32_t pos)
{
    if(GET_MODIFIER(newVal))
    {
        return pThis->value(pos);
    }

    // const ParLimOption<SEQ::Gradients>& _seqLimits = pThis->seqLimits().getGradients();

    switch(newVal)
    {
    case MRI_STD_GRADIENT_MODE_ULTRAFAST:
        pThis->prot().getsGRADSPEC().setucMode(SEQ::GRAD_ULTRAFAST);
        break;

    case MRI_STD_GRADIENT_MODE_FAST:
        pThis->prot().getsGRADSPEC().setucMode(SEQ::GRAD_FAST);
        break;

    case MRI_STD_GRADIENT_MODE_NORMAL:
        pThis->prot().getsGRADSPEC().setucMode(SEQ::GRAD_NORMAL);
        break;

    case MRI_STD_GRADIENT_MODE_WHISPER:
        pThis->prot().getsGRADSPEC().setucMode(SEQ::GRAD_WHISPER);
        break;

    default:
        SEQ_TRACE_INFO.print("Gradient mode '%u' not recognised", newVal);
        return MRI_STD_INVALID;
    }

    return pThis->value(pos);
}


// ------------------------------------------------------------------------------
// Function    : fEPISaveUncombinedIsAvailable()
// ------------------------------------------------------------------------------
//
// Overloaded is-available handler for "Save uncombined".
//
// This overloaded version is activated by defining
// EPI_DISABLE_SAVE_UNCOMBINED in a_ep2d.h. It is used to remove the option
// "Save uncombined" in the following cases:
//
// 1. EPI BOLD sequences (ep2d_bold, ep2d_pace). In the case of many images
//    (as with epi), "Save uncombined" may overload the database. Furthermore,
//    "Save uncombined" in not compatible with mosaic mode.
//
// 2. EPI Diffusion sequence (ep2d_diff). "Save uncombined" not supported by
//    IcePrograDiffusion2D (CHARM 327088).
//
// ------------------------------------------------------------------------------
bool EpCommonUINS::fEPISaveUncombinedIsAvailable(LINK_BOOL_TYPE* const /* pThis */, int32_t)
{
    return false;
}


// ------------------------------------------------------------------------------
// Function    : fEPIFilterBoolIsAvailable()
// ------------------------------------------------------------------------------
//
// Overloaded is-available handler for "Save unfiltered".
//
// This overloaded version is used to remove the option
// "Save unfiltered" from various filters. It always returns false.
// ------------------------------------------------------------------------------
bool EpCommonUINS::fEPIFilterBoolIsAvailable(LINK_BOOL_TYPE* const /*pThis*/, int32_t /*pos*/)
{
    return false;
}

// ------------------------------------------------------------------------------
//  Function    : fEPITriggeringSolve()
// ------------------------------------------------------------------------------
//  Description : Solve handler for physio mode to allow triggering to be
//                switched off when number of concatenations > 1. This
//                is necessary with int32_t TR triggering mode, which uses
//                a repositioned concatenation loop outside of 'outer slices'
//                loop. Concatenations are only allowed when triggering is on,
//                leading to a conflict when switching triggering off.
//
//  Return      : 0              if no solution possible
//                MRI_STD_STRING on success. The text in arg_list is then used
//                                to format the confirmation message.
// ------------------------------------------------------------------------------

unsigned EpCommonUINS::fEPITriggeringSolve
    (
    LINK_SELECTION_TYPE*  const pThis,   // ...
    char**                arg_list,      // receives confirmation message
    const void*         /*  pToAddMemory */,
    const MrProtocolData::MrProtData*          pOrigProt,     // Original protocol
 int32_t                  lIndex          // Array index reserved
    )
{
    const SEQ::PhysioMethod FirstMethod = pThis->prot().getsPhysioImaging().getlMethod1();

    // check for known parameter conflict
    if((pThis->prot().getsSliceArray().getlConc() > 1) && (FirstMethod == SEQ::METHOD_NONE))
    {
        // set number of concats to 1

        pThis->prot().getsSliceArray().setlConc(1);

        // formulate confirmation message for user

        if(arg_list)
        {
            //arg_list[0] = (char*) "$&OK$&Undo$$Switching triggering off \nwill set concatenations to 1";

            // new triggering mode

            arg_list[0] = (char*)MRI_STD_FIRST_SIGNAL_MODE_LABEL;
            unsigned nID = pThis->value(lIndex);
            arg_list[10] = (char*)static_cast<int64_t>(RES_ID(nID));
            pThis->format(nID, arg_list+11, lIndex);
            arg_list[16] = (char*)static_cast<int64_t>(pThis->unitId(arg_list+17, lIndex));

            // old triggering mode

            MrProt _storage;
            void* pVoid     = pThis->storeInMemory(_storage, lIndex);
            pThis->recallMemory(NULL, pOrigProt, lIndex);
            nID = pThis->value(lIndex);
            arg_list[4]     = (char*)static_cast<int64_t>(RES_ID(nID));
            pThis->format(nID, arg_list+5, lIndex);

            // restore new state

            pThis->recallMemory(pVoid, _storage, lIndex);
            pThis->clearMemory(pVoid, lIndex);

            // concatenations changes

            int32_t lConcatsOld = pOrigProt->getsSliceArray().getlConc();
            int32_t lConcatsNew = pThis->prot().getsSliceArray().getlConc();

            arg_list[20]  = (char*)MRI_STD_CONC_LABEL;
            arg_list[21]  = (char*) '\0';
            arg_list[22]  = (char*) '\0';

            arg_list[24] = (char*)MRI_STD_INT;

            arg_list[25] = (char*)static_cast<int64_t>(lConcatsOld);
            //arg_list[25] = (char*) _this->prot().concatenations();
            //arg_list[25] = (char*) pOrigProt->concatenations();
            //arg_list[24]  = (char*) MRI_STD_ON;
            arg_list[30]  = (char*)MRI_STD_INT;

            arg_list[31] = (char*)static_cast<int64_t>(lConcatsNew);
            arg_list[36]  = (char*)MRI_STD_EMPTY;
            arg_list[39]  = (char*) '\n';
        }

        // if new protocol with one concatenation is not valid, TR may need to be increased

        if(! pThis->sequence().prepareForBinarySearch(&pThis->prot()))
        {
            if(getUI(pThis)->EpCommonUI::m_lNeededTR > pThis->prot().getalTR()[0])
            {
                // set new minimum TR value in new protocol

                double dTRMin_ms = static_cast<double>(getUI(pThis)->EpCommonUI::m_lNeededTR / 1000.0);     // * TR in ms *
                fUICSetWithinLimits(pThis, dTRMin_ms, MR_TAG_TR, UIC_SET_WITHIN_LIMITS_ROUNDUP | UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, 0);
            }

            else
            {
                // can't help further

                return 0;
            }
        }

        // finish

        return MRI_STD_CONFIRMATION_MSG;
    }

    // sorry - no help for other conflict types

    return(0);

}



#endif // #ifdef WIN32


EpCommonUI::EpCommonUI()
    : m_lNeededTR(0)
    , m_lNeededTI(0)
    , m_lNeededTE(0)
    , m_bNeedOtherTETITR(false)
    , m_lDebug_SEQ_UILink(0)
    , m_pEPIRO(NULL)
    , m_pREOInfo(NULL)
    , m_pfCalculateTRTIFill(NULL)
{
    int i;
    for(i=0;i<100;i++)
    {
        m_tLine_epi_EchoSpacing_GetToolTip[i]='\0';
    }
}

EpCommonUI::~EpCommonUI()
{
}


// ------------------------------------------------------------------------------
// Function    : fEPIStdResetSolveHandlerControlTETITR
// ------------------------------------------------------------------------------
//
// Description : resets solve handler control flag and needed TE,TR,TI.
// Return      : void
//
// ------------------------------------------------------------------------------
void EpCommonUI::fEPIStdResetSolveHandlerControlTETITR()
{
    m_bNeedOtherTETITR = false;
    m_lNeededTE        = -1;    // should now be != rMrProt.te()[0]
    m_lNeededTR        = -1;    // should now be != rMrProt.tr()[0]
    m_lNeededTI        = -1;    // should now be != rMrProt.ti()[0]
}


// ------------------------------------------------------------------------------
// Function    : fEPIStdUICheckNeededTETITR
// ------------------------------------------------------------------------------
//
// Description : - Checks, if TE, TR, TI required by the sequence are identical
//                 to the ones in the current protocol.
//               - If required: Puts the required values for TE, TR, TI into the
//                 current protocol so that it becomes consistent.
// Return      : false, if at least one needed value is not equal to the one in
//                      the current protocol and if the sequence is not allowed
//                      to update the protocol
//               true,  else
//
// ------------------------------------------------------------------------------
bool EpCommonUI::fEPIStdUICheckNeededTETITR(MrProt &rMrProt, SeqLim &rSeqLim, int32_t lSEQNeedsTE, int32_t lSEQNeedsTI, int32_t lSEQNeedsTR, int32_t lTEContrastIndex)
{
    // round TE, TR and TE to UI precision - otherwise this code can set values to the protocol which
    // cannot be realized by UI interaction ==> this causes inconsistencies in UI handling
    // exclude STIR with sequential scheme as this requires the sequence to set TR value with higher precision
    if((rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rMrProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        m_lNeededTE = lSEQNeedsTE;
        m_lNeededTR = lSEQNeedsTR;
        m_lNeededTI = lSEQNeedsTI;
    }
    else
    {
    m_lNeededTE = static_cast<int>(lSEQNeedsTE/ 100.0 + 0.5) * 100;
    m_lNeededTR = static_cast<int>(lSEQNeedsTR/ 1000.0 + 0.5) * 1000;
    m_lNeededTI = static_cast<int>(lSEQNeedsTI/ 1000.0 + 0.5) * 1000;
#ifdef EP2D_MS
    if ( rMrProt.getsPrepPulses().getucFreezeSupTissue() )
    {
        m_lNeededTI = lSEQNeedsTI;
    }
#endif
    }

    bool bOtherTRNeeded = false;
    bool bOtherTENeeded = false;
    bool bOtherTINeeded = false;

    if(m_lNeededTR != rMrProt.tr()[0])
    {
        if(!rSeqLim.isContextPrepForBinarySearch() || (m_lDebug_SEQ_UILink&2))
        {
            SEQ_TRACE_INFO.print("=============================================================");
            SEQ_TRACE_INFO.print("WantedTR(=%7d) !=  NeededTR(=%7d)", rMrProt.tr()[0], m_lNeededTR);
            SEQ_TRACE_INFO.print("=============================================================");
        }
        bOtherTRNeeded = true;
    }

    if(m_lNeededTI != rMrProt.ti()[0])
    {
        if(!rSeqLim.isContextPrepForBinarySearch() || (m_lDebug_SEQ_UILink&2))
        {
            SEQ_TRACE_INFO.print("=============================================================");
            SEQ_TRACE_INFO.print("WantedTI(=%7d) !=  NeededTI(=%7d)", rMrProt.ti()[0], m_lNeededTI);
            SEQ_TRACE_INFO.print("=============================================================");
        }
        bOtherTINeeded = true;
    }

    if(rMrProt.te()[lTEContrastIndex] != m_lNeededTE)
    {
        if(!rSeqLim.isContextPrepForBinarySearch() || (m_lDebug_SEQ_UILink&2))
        {
            SEQ_TRACE_INFO.print("=============================================================");
            SEQ_TRACE_INFO.print("WantedTE(=%7d) !=  NeededTE(=%7d)", rMrProt.te()[lTEContrastIndex], m_lNeededTE);
            SEQ_TRACE_INFO.print("=============================================================");
        }
        bOtherTENeeded=true;
    }

    if(bOtherTENeeded || bOtherTINeeded || bOtherTRNeeded)
    {
        if(rSeqLim.isContextPrepForMrProtUpdate())
        {
            if(bOtherTENeeded)
            {
                rMrProt.te()[lTEContrastIndex] = static_cast<int>(m_lNeededTE);
                bOtherTENeeded = false;
            }

            if(bOtherTINeeded)
            {
                rMrProt.ti()[0] = static_cast<int>(m_lNeededTI);
                bOtherTINeeded = false;
            }

            if(bOtherTRNeeded)
            {
                rMrProt.tr()[0] = static_cast<int>(m_lNeededTR);
                bOtherTRNeeded = false;
            }
        }
    }

    if(bOtherTENeeded || bOtherTINeeded || bOtherTRNeeded)
    {
        m_bNeedOtherTETITR = true;
    }

    return !m_bNeedOtherTETITR;
}

// ------------------------------------------------------------------------------
// Function    : EpCommonUI::fEPIStdInit
// ------------------------------------------------------------------------------
//
// Description : - sets static variables of the module
//               - specifies standard hard limits for epi sequences
// Return      : true (if success) or false
//
// ------------------------------------------------------------------------------
bool EpCommonUI::fEPIStdInit(SeqLim &rSeqLim, SeqBuildBlockEPIReadOut* _pEPIRO, ReorderInfo* _pREOInfo)
{
    //-------------------------------------------------------------------------------------
    // check the pointers we've got
    //-------------------------------------------------------------------------------------
    if(_pEPIRO==NULL || _pREOInfo==NULL)
    {
        return false;
    }
    m_pEPIRO   = _pEPIRO;
    m_pREOInfo = _pREOInfo;

    //-------------------------------------------------------------------------------------
    // general settings
    //-------------------------------------------------------------------------------------
    rSeqLim.setSequenceOwner(SEQ_OWNER_SIEMENS);
    rSeqLim.setAllowedFrequency(5000000, 500000000);
    rSeqLim.setRequiredGradAmpl(1.00);
    rSeqLim.setRequiredGradSlewRate(1.00);

    //-------------------------------------------------------------------------------------
    // image reconstruction
    //-------------------------------------------------------------------------------------
    rSeqLim.setReconstructionMode(SEQ::RECONMODE_MAGNITUDE);

    //-------------------------------------------------------------------------------------
    // image resolution                  (     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------
    rSeqLim.setBaseResolution(64, 256, SEQ::INC_64, 128);
    rSeqLim.setPhaseOversampling(0.000, 0.500, 0.010, 0.000);
    rSeqLim.setPELines(32, 256, 1, 128);
    rSeqLim.set2DInterpolation(SEQ::OFF, SEQ::ON);

    //-------------------------------------------------------------------------------------
    // bandwidth per pixel [Hz]     (ADC#,     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------
    rSeqLim.setBandWidthPerPixel(0, 750, 10000, 2, 750);

    //-------------------------------------------------------------------------------------
    // timing                            (     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------

    // modify maximum TE and TR after introduction of modified PAT selection solve-handler (CHARM 315966)
    // previous values did not allow solve-handler to work correctly with large number of slices
    rSeqLim.setTE(0, 1000, 400000, 100, 1000);
    rSeqLim.setTR(0, 10000, 30000000, 100, 10000);

    //-------------------------------------------------------------------------------------
    // slices and their attributes       (     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------
    rSeqLim.setSlices(1, K_NO_SLI_MAX, 1, 1);
    rSeqLim.setSliceThickness(0.100, 10.000, 0.100, 5.000);
    rSeqLim.setSliceDistanceFactor(0.000, 8.000, 0.010, 0.500);
    rSeqLim.setSliceSeriesMode(SEQ::INTERLEAVED);
    rSeqLim.enableSliceShift();
    rSeqLim.enableMSMA();
    rSeqLim.enableOffcenter();
    rSeqLim.setAllowedSliceOrientation(SEQ::DOUBLE_OBLIQUE);

    //-------------------------------------------------------------------------------------
    // RF                                (     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------
    rSeqLim.setExtSrfFilename           (                          "%MEASCONST%/extrf/extrf.dat");
    rSeqLim.setFlipAngle                (  10.000,     130.000,       1.000,      90.000);
    //rSeqLim.setB1Corr                   (                              SEQ::OFF, SEQ::ON);

    // ------------------------------------------------------------------------------------
    // loop control
    // ------------------------------------------------------------------------------------
    rSeqLim.setIntro(SEQ::OFF, SEQ::ON);
    rSeqLim.setAveragingMode(SEQ::OUTER_LOOP);
    rSeqLim.setEllipticalScanning(SEQ::OFF);

    //-------------------------------------------------------------------------------------
    // preparation pulses                (     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------
    rSeqLim.setRSats(0, 6, 1, 0);
    rSeqLim.setRSatThickness(3.000, 150.000, 1.000, 50.000);
    rSeqLim.setMTC(SEQ::OFF, SEQ::ON);
    rSeqLim.setPSatMode(SEQ::PSAT_NONE, SEQ::PSAT_SINGLE_REG, SEQ::PSAT_DOUBLE_REG);
    rSeqLim.setPSatThickness(3.000, 150.000, 1.000, 50.000);
    rSeqLim.setPSatGapToSlice(5.000, 50.000, 1.000, 10.000);
    rSeqLim.setFatWaterContrast(FatWaterContrast_FatSaturation, FatWaterContrast_WaterExcitation, FatWaterContrast_Standard);

    //-------------------------------------------------------------------------------------
    // Averages/Repetitions              (     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------
    rSeqLim.setAverages(1, 32, 1, 1);
    rSeqLim.setRepetitions(0, 511, 1, 0);
    rSeqLim.setRepetitionsDelayTime(0, 2100000000, 100000, 0);

    //-------------------------------------------------------------------------------------
    // Concatenations                    (     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------
    rSeqLim.setConcatenations(1, K_NO_SLI_MAX, 1, 1);

    //-------------------------------------------------------------------------------------
    // Physiologic measurements          (     min,         max,         inc,         def);
    //-------------------------------------------------------------------------------------
    rSeqLim.addPhysioMode(SEQ::SIGNAL_CARDIAC, SEQ::METHOD_TRIGGERING);
    rSeqLim.addPhysioMode(SEQ::SIGNAL_EXT_2, SEQ::METHOD_TRIGGERING);
    rSeqLim.addPhysioMode(SEQ::SIGNAL_RESPIRATION, SEQ::METHOD_TRIGGERING);
    rSeqLim.setPhases(1, 1, 1, 1);

    //-------------------------------------------------------------------------------------
    // shimming
    //-------------------------------------------------------------------------------------
    if (!SysProperties::isUHFSystem())
        rSeqLim.getAdjShim().setDef(SEQ::ADJSHIM_STANDARD);

    // ------------------------------------------------------------------------------------
    // database control
    // ------------------------------------------------------------------------------------
    rSeqLim.setMultipleSeriesMode(SEQ::MULTIPLE_SERIES_EACH_MEASUREMENT, SEQ::MULTIPLE_SERIES_OFF);

    // ------------------------------------------------------------------------------------
    // filters
    // ------------------------------------------------------------------------------------
    rSeqLim.setFilterType(SEQ::FILTER_NONE,
        SEQ::FILTER_RAW,
        SEQ::LARGE_FOV,
        SEQ::ELLIPTICAL,
        SEQ::HAMMING,
        SEQ::PRESCAN_NORMALIZE
        /* , SEQ::FILTER_BIFIC */
        );


    // ------------------------------------------------------------------------------------
    // CHARM 348199/356618: set coilCombineMode to default SUM_OF_SQUARES
    // ------------------------------------------------------------------------------------
#ifndef SUPPORT_iPAT_TGSE
    rSeqLim.setCoilCombineMode(SEQ::COILCOMBINE_SUM_OF_SQUARES, SEQ::COILCOMBINE_ADAPTIVE_COMBINE);
#else
    rSeqLim.setCoilCombineMode(SEQ::COILCOMBINE_ADAPTIVE_COMBINE);    // IcePAT only works with ACC
#endif

    //-------------------------------------------------------------------------------------
    // etc.
    //-------------------------------------------------------------------------------------
    rSeqLim.getEllipticalScanning().setDisplayMode(SEQ::DM_OFF);
    rSeqLim.setSarSolverStrategy(SARSolverStrategy_TR);
    rSeqLim.getReadoutFOV().setDisplayMode(SEQ::DM_EDIT);
    rSeqLim.getPhaseFOV().setDisplayMode(SEQ::DM_EDIT);

    return true;
}

// ------------------------------------------------------------------------------
// Function    : fEPIStdRegisterTIHandlers
// ------------------------------------------------------------------------------
//
// Description : registers TI handlers defined above
// Return      : true (if success) or false
//
// ------------------------------------------------------------------------------
#ifdef WIN32
bool EpCommonUI::fEPIStdRegisterTIHandlers(SeqLim &rSeqLim, fEPIStdInit_pfCalculateTRTIFillTimes _pfCalculateTRTIFill)
#else
bool EpCommonUI::fEPIStdRegisterTIHandlers(SeqLim& /* rSeqLim */, fEPIStdInit_pfCalculateTRTIFillTimes _pfCalculateTRTIFill)
#endif
{
    //-------------------------------------------------------------------------------------
    // check the pointers we've got
    //-------------------------------------------------------------------------------------
    if(_pfCalculateTRTIFill==NULL)
    {
        return false;
    }
    m_pfCalculateTRTIFill = _pfCalculateTRTIFill;

#ifdef WIN32
    //-------------------------------------------------------------------------------------
    // register solve-Handler for inversion
    //-------------------------------------------------------------------------------------
    m_Inversion.registerSolveHandler(rSeqLim, MR_TAG_INVERSION, EpCommonUINS::fEPIStdSolveSelection);

    //-------------------------------------------------------------------------------------
    // register get-Limit and solve-Handler for TI
    //-------------------------------------------------------------------------------------
    char tMrTagTI[64];
    sprintf(tMrTagTI, "%s.0\0", MR_TAG_TI);

    m_TI.registerGetLimitsHandler(rSeqLim, tMrTagTI, EpCommonUINS::_TIGetLimits);
    m_TI.registerSolveHandler(rSeqLim, tMrTagTI, EpCommonUINS::_solveTETITR_TI);
#endif
    return true;
}


// ------------------------------------------------------------------------------
// Function    : fEPIStdRegisterEPIFactorHandlers
// ------------------------------------------------------------------------------
//
// Description : registers EPI factor handlers defined above
// Return      : true (if success) or false
//
// ------------------------------------------------------------------------------
#ifdef WIN32
bool EpCommonUI::fEPIStdRegisterEPIFactorHandlers(SeqLim &rSeqLim)
#else
bool EpCommonUI::fEPIStdRegisterEPIFactorHandlers(SeqLim& /* rSeqLim */)
#endif
{
#ifdef WIN32
    //-------------------------------------------------------------------------------------
    // register solve-Handler for EPI factor
    //-------------------------------------------------------------------------------------
    if (MrUILinkLimited<int32_t> *_epiFactor = _search< MrUILinkLimited<int32_t> >(rSeqLim, MR_TAG_EPI_FACTOR))
    {
        _epiFactor->registerSolveHandler(EpCommonUINS::_solveTETITR_PPF_EPIFactorConflict);
    }
#endif
    return true;
}


//  --------------------------------------------------------------------------
//
//  Name        : registerUI
//
//  Description :
/// \brief        This method registers all given set / get / Solve - handlers
///
///               It can be executed on the measurement system, too, but is empty there.
///
///
///               It returns an NLS status
///
//  Return      : int32_t
//
//  --------------------------------------------------------------------------

#ifdef WIN32
NLS_STATUS EpCommonUI::registerUI(SeqLim &rSeqLim)
#else
NLS_STATUS EpCommonUI::registerUI(SeqLim& /* rSeqLim */)
#endif
{
#ifdef WIN32
    //-------------------------------------------------------------------------------------
    // read debug-masks
    //-------------------------------------------------------------------------------------
    m_lDebug_SEQ_UILink   = getMaskFromRegistry("DEBUG_SEQ_UILink");

    //-------------------------------------------------------------------------------------
    // get standard imaging TD/TI/TR/TE handling for the UI
    //-------------------------------------------------------------------------------------
    fStdImagingInitPost(rSeqLim);

    //-------------------------------------------------------------------------------------
    // create bandwidth name tag
    //-------------------------------------------------------------------------------------
    char tMrTagBW[64], tMrTagTE[64];
    sprintf(tMrTagBW, "%s.0\0", MR_TAG_BANDWIDTH);
    sprintf(tMrTagTE, "%s.0\0", MR_TAG_TE);

    //-------------------------------------------------------------------------------------
    // register (standard) solve handlers
    //-------------------------------------------------------------------------------------
    m_Bandwidth.registerSolveHandler(rSeqLim, tMrTagBW, EpCommonUINS::fEPIStdSolveDoubleParam);
    m_TE.registerSolveHandler(rSeqLim, tMrTagTE, EpCommonUINS::_solveTETITR_TE);
    m_ReadFOV.registerSolveHandler(rSeqLim, MR_TAG_READOUT_FOV, EpCommonUINS::fEPIStdSolveDoubleParam);
    m_PhaseFOV.registerSolveHandler(rSeqLim, MR_TAG_PHASE_FOV, EpCommonUINS::fEPIStdSolveDoubleParam);
    m_SliceThick.registerSolveHandler(rSeqLim, MR_TAG_SLICE_THICKNESS, EpCommonUINS::fEPIStdSolveDoubleParam);
    m_PhaseResolution.registerSolveHandler(rSeqLim, MR_TAG_PHASE_RESOLUTION, EpCommonUINS::fEPIStdSolveDoubleParam);
    m_PhaseOS.registerSolveHandler(rSeqLim, MR_TAG_PHASE_OVERSAMPLING, EpCommonUINS::fEPIStdSolveDoubleParam);
    m_PhasePartF.registerSolveHandler(rSeqLim, MR_TAG_PHASE_PARTIAL_FOURIER, EpCommonUINS::fEPIStdSolveSelection);
    m_Dimension.registerSolveHandler(rSeqLim, MR_TAG_DIMENSION, EpCommonUINS::fEPIStdSolveSelection);
    m_AdjustmentMode.registerSolveHandler(rSeqLim, MR_TAG_ADJUSTMENT_MODE, EpCommonUINS::fUILinkAdjustmentModeSolve);
    m_AdjustmentMode.registerSetValueHandler(rSeqLim, MR_TAG_ADJUSTMENT_MODE, EpCommonUINS::fUILinkAdjustmentModeSetValue);

    //-------------------------------------------------------------------------------------
    // register advanced solve handler for base resolution
    //-------------------------------------------------------------------------------------
    m_BaseResolution.registerSolveHandler(rSeqLim, MR_TAG_BASE_RESOLUTION, EpCommonUINS::_solveBaseResolution);
    m_BaseResolution.registerTryHandler(rSeqLim, MR_TAG_BASE_RESOLUTION, EpCommonUINS::_epi_BaseResolutionTry);

    //-------------------------------------------------------------------------------------
    // register handlers for echo spacing
    //-------------------------------------------------------------------------------------
    m_EchoSpacing.registerIsAvailableHandler(rSeqLim, MR_TAG_ECHO_SPACING, EpCommonUINS::_epi_EchoSpacingIsAvailable);
    m_EchoSpacing.registerGetValueHandler(rSeqLim, MR_TAG_ECHO_SPACING, EpCommonUINS::_epi_EchoSpacingGetValue);
    m_EchoSpacing.registerSetValueHandler(rSeqLim, MR_TAG_ECHO_SPACING, EpCommonUINS::_epi_EchoSpacingSetValue);
    m_EchoSpacing.registerGetLimitsHandler(rSeqLim, MR_TAG_ECHO_SPACING, EpCommonUINS::_epi_EchoSpacingGetLimits);
    m_EchoSpacing.registerToolTipHandler(rSeqLim, MR_TAG_ECHO_SPACING, EpCommonUINS::_epi_EchoSpacing_GetToolTip);
    m_EchoSpacing.registerSolveHandler(rSeqLim, MR_TAG_ECHO_SPACING, EpCommonUINS::fEPIStdSolveDoubleParam);

    if(LINK_BOOL_TYPE* pBool = _create< LINK_BOOL_TYPE >(rSeqLim, MR_TAG_MANUAL_ECHO_SPACING))
    {
        pBool->registerGetLabelIdHandler(EpCommonUINS::_epi_FreeEchoSpacing_GetLabelId);
        pBool->registerGetOptionsHandler(EpCommonUINS::_epi_FreeEchoSpacing_GetOptions);
        pBool->registerGetValueHandler(EpCommonUINS::_epi_FreeEchoSpacing_GetValue);
        pBool->registerSetValueHandler(EpCommonUINS::_epi_FreeEchoSpacing_SetValue);
        pBool->registerSolveHandler(EpCommonUINS::_epi_FreeEchoSpacing_Solve);
    }

    //-------------------------------------------------------------------------------------
    // register modified handlers  for "Gradient mode"
    //-------------------------------------------------------------------------------------
    m_GradMode.registerGetOptionsHandler(rSeqLim, MR_TAG_GRADIENT_MODE, EpCommonUINS::fEPIGradientModeGetOptions);
    m_GradMode.registerGetValueHandler(rSeqLim, MR_TAG_GRADIENT_MODE, EpCommonUINS::fEPIGradientModeGetValue);
    m_GradMode.registerSetValueHandler(rSeqLim, MR_TAG_GRADIENT_MODE, EpCommonUINS::fEPIGradientModeSetValue);
    m_GradMode.registerSolveHandler(rSeqLim, MR_TAG_GRADIENT_MODE, EpCommonUINS::fEPIStdSolveSelection);

    //-------------------------------------------------------------------------------------
    // Register solve handler for triggering
    //-------------------------------------------------------------------------------------

    if(LINK_SELECTION_TYPE* pSelect  = _search<LINK_SELECTION_TYPE >(rSeqLim, MR_TAG_FIRST_SIGNAL_MODE))
    {
        pSelect->registerSolveHandler(EpCommonUINS::fEPITriggeringSolve);
    }

    //-----------------------------------------------------------------------------------
    // do not display 'Pause after Measurement'
    //-----------------------------------------------------------------------------------
    if(MrUILinkArray* _pPause = _search< MrUILinkArray >(rSeqLim, MR_TAG_MEASUREMENT_DELAY_TIMES))
    {
        _pPause->unregister();
    }

    //-----------------------------------------------------------------------------------
    // register solve handler for measurements:
    //
    // This solve handler became necessary due to CHARM 305893 (TR should always be 
    // physically correct). The sequence has to take into account the time for SUBFINI
    // and SUBSTRT between two repetitions for TR and measuerement time calculations.
    // If this is done only for multiple measurements, then TR may have to be increased,
    // if the number of measurements is increased from 1 to any other value.
    //
    // Unfortuneatly the UI already provides a solve-handler for measurements dealing with
    // problems of the EVA-protocols. Incorporating an original solve-handler into an
    // overloaded one is no trivial task, because reformatting of the message box is
    // necessary.
    //
    // So the problem is solved within the method SeqLoopEP2D::calcFillTimesOnly  by taking
    // into account the 340us for SUBFINI and SUBSTRT also for one measurement only.
    // Another possibility (which would be more transparent for the user) for solving the
    // problem would be to set the minimum for delayTimeInTR to 500us.
    //
    // SO THE SOLVE HANDLER IS NOT REGISTERED!
    //
    //-----------------------------------------------------------------------------------
    //if(LINK_LONG_TYPE *_measurements = _search<LINK_LONG_TYPE> (pSeqLim, MR_TAG_MEASUREMENTS))
    //{
    //    _measurements->registerSolveHandler (fEPIStdSolveLongParam); 
    //}

    //-----------------------------------------------------------------------------------
    // register solve handler for 'delay in TR'
    //-----------------------------------------------------------------------------------
    m_DelayInTR.registerSolveHandler(rSeqLim, MR_TAG_DELAY_IN_TR, EpCommonUINS::fEPIStdSolveDoubleParam);

#ifdef SUPPORT_iPAT_a_ep2d
    //-----------------------------------------------------------------------------------
    // register standard solve handlers for PAT
    //
    // CHARM 315966: restored PAT-specific selection solve-handler, which was not
    //               correctly copied from VA21B to VA25A archive.
    //
    // CHARM 370315: Ref lines must be a multiple of the acceleration factor
    //-----------------------------------------------------------------------------------
    m_PATMode.registerSolveHandler(rSeqLim, MR_TAG_PAT_MODE, EpCommonUINS::fEPIPATModeSolveSelection);

#ifndef COMPILE_EP2D_DIFF
    m_PATMode.registerSetValueHandler(rSeqLim, MR_TAG_PAT_MODE, EpCommonUINS::fUILinkPATModeSetValue);
#endif

    m_PATRefScan.registerSolveHandler(rSeqLim, MR_TAG_PAT_REF_SCAN_MODE, EpCommonUINS::fEPIStdSolveSelection);

    m_PATAccelerationFactor.registerSolveHandler(rSeqLim, MR_TAG_PAT_ACC_PE, EpCommonUINS::fEPIStdSolveLongParam);
    m_PATAccelerationFactor.registerSetValueHandler(rSeqLim, MR_TAG_PAT_ACC_PE, EpCommonUINS::fEPIAccPESetValue);
    m_PATAccelerationFactor.registerToolTipHandler(rSeqLim, MR_TAG_PAT_ACC_PE, EpCommonUINS::fUILinkPATAccelFactorGetToolTipID);

    m_PATReferenceLines.registerSolveHandler(rSeqLim, MR_TAG_PAT_LINES_PE, EpCommonUINS::fEPIStdSolveLongParam);
    m_PATReferenceLines.registerGetLimitsHandler(rSeqLim, MR_TAG_PAT_LINES_PE, EpCommonUINS::fEPIRefLinesPEGetLimits);
    m_PATReferenceLines.registerIsAvailableHandler(rSeqLim, MR_TAG_PAT_LINES_PE, EpCommonUINS::fEPIRefLinesPEIsAvailable);

    m_AccelerationFactorSlice.registerSetValueHandler(rSeqLim, MR_TAG_SLICE_ACCEL_FACTOR, EpCommonUINS::fUILinkSliceAccelSetValue);
    m_AccelerationFactorSlice.registerSolveHandler(rSeqLim, MR_TAG_SLICE_ACCEL_FACTOR, EpCommonUINS::fUILinkSliceAccelSolve);

    m_FOVShiftFactor.registerGetLimitsHandler(rSeqLim, MR_TAG_SLICE_ACCEL_FOV_SHIFT_FACTOR, EpCommonUINS::fUILinkFOVShiftFactorGetLimits);
    m_FOVShiftFactor.registerIsAvailableHandler(rSeqLim, MR_TAG_SLICE_ACCEL_FOV_SHIFT_FACTOR, EpCommonUINS::fUILinkFOVShiftFactorIsAvailable);

    m_LocalShim.registerTryHandler(rSeqLim, MR_TAG_LOCAL_SHIM, EpCommonUINS::fUILinkLocalShimTry);

    m_FatSup.registerSolveHandler(rSeqLim, MR_TAG_FAT_WATER_CONTRAST, EpCommonUINS::fUILinkFatSuppresionSolve);

    m_CoilElements.registerElmSetValueHandler(rSeqLim, MR_TAG_COIL_ELEMENTS, EpCommonUINS::fUILinkCoilElementsSetValue);

#ifdef EP2D_MS
    m_NumberOfSlices.registerElmGetLimitsHandler( rSeqLim, MR_TAG_SLICE_GROUP_LIST, MR_TAG_SG_SIZE, EpCommonUINS::fUILinkNumberOfSlicesGetLimits );
    m_FatSupOpt.registerSolveHandler( rSeqLim, MR_TAG_FAT_SUP_OPTIMIZATION, EpCommonUINS::fEPIStdSolveSelection     );

    m_FreezeSuppressedTissue.registerSolveHandler   ( rSeqLim, MR_TAG_FREEZE_SUPPRESSED_TISSUE, EpCommonUINS::fEPIStdSolveBoolParamConflict         );
    m_FreezeSuppressedTissue.registerSetValueHandler( rSeqLim, MR_TAG_FREEZE_SUPPRESSED_TISSUE, EpCommonUINS::fUILinkFreezeSuppressedTissueSetValue );
#endif //#ifdef EP2D_MS

#endif //#ifdef SUPPORT_iPAT_a_ep2d

    m_DistortionCorrection.registerSetValueHandler    ( rSeqLim, MR_TAG_FLT_DISCOR_PROP_MODE,                        EpCommonUINS::fDistortionCorrectionSetValue    );

    m_SFC.registerSetValueHandler(rSeqLim, MR_TAG_STATIC_FIELD_CORRECTION, EpCommonUINS::fUILinkSFCSetValue);

#ifdef ZOOM_2DRF
    //-------------------------------------------------------------------------------------
    // register ZOOM solve handlers
    //-------------------------------------------------------------------------------------
    m_ExitPulse.registerSetValueHandler(rSeqLim, MR_TAG_EXCIT_PULSE, Ep2d_zoom_UINS::fUILinkExcitationPulseSetValueNew);
    m_ExitPulse.registerSolveHandler(rSeqLim, MR_TAG_EXCIT_PULSE, EpCommonUINS::fEPIStdSolveSelection);

    m_PtxVolProp.registerElmGetOptionsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY, Ep2d_zoom_UINS::fUILinkPTXVolPropGetOptionsNew);
    m_PTXVolPDim.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_PDIM, Ep2d_zoom_UINS::fUILinkPTXVolPDimGetLimitsNew);
    m_PTXVolPDim.registerElmSolveHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_PDIM, EpCommonUINS::_solveTETITR_TI);

    // refresh PTXVolume during GSP action 
    m_PTXVolRot.registerElmGetValueHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_ROTATION, Ep2d_zoom_UINS::fUILinkPTXVolRotGetValueNew);
#ifdef BOLD
    m_PTXVolPosGSP.registerElmGetValueHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POSITION_SBCS, Ep2d_zoom_UINS::fUILinkPTXVolGSPPosGetValueNew);
    m_PTXVolRotGSP.registerElmGetValueHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_ORI_AND_ROTATION, Ep2d_zoom_UINS::fUILinkPTXVolGSPRotGetValueNew);
#endif //BOLD

    // make PTXVolume fields NOT editable for type optimization
    m_PTXVolRot.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_ROTATION, Ep2d_zoom_UINS::fUILinkPTXVolRotGetLimitsNew);
    m_PTXVolRDim.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_RDIM, Ep2d_zoom_UINS::fUILinkPTXVolRDimGetLimitsNew);
    m_PTXVolSDim.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_SDIM, Ep2d_zoom_UINS::fUILinkPTXVolSDimGetLimitsNew);
#ifdef ZOOM_EXTENDED
    m_PTXVolPDim.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_SDIM, Ep2d_zoom_UINS::fUILinkPTXVolPDimGetLimitsNew);

#endif
    m_PTXVolPosSag.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POS_SAG, Ep2d_zoom_UINS::fUILinkPTXVolPosSagGetLimitsNew);
    m_PTXVolPosSagSBCS.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POS_SBCS_SAG, Ep2d_zoom_UINS::fUILinkPTXVolPosSag_SBCSGetLimitsNew);
    m_PTXVolPosCor.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POS_COR, Ep2d_zoom_UINS::fUILinkPTXVolPosCorGetLimitsNew);
    m_PTXVolPosCorSBCS.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POS_SBCS_COR, Ep2d_zoom_UINS::fUILinkPTXVolPosCor_SBCSGetLimitsNew);
    m_PTXVolPosTra.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POS_TRA, Ep2d_zoom_UINS::fUILinkPTXVolPosTraGetLimitsNew);
    m_PTXVolPosTraSBCS.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POS_SBCS_TRA, Ep2d_zoom_UINS::fUILinkPTXVolPosTra_SBCSGetLimitsNew);
    m_PTXVolOriAlpha.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_ORI_ALPHA, Ep2d_zoom_UINS::fUILinkPTXVolOriAlphaGetLimitsNew);
    m_PTXVolOriBeta.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_ORI_BETA, Ep2d_zoom_UINS::fUILinkPTXVolOriBetaGetLimitsNew);

    m_PTXVolPos.registerElmGetOptionsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POSITION2, Ep2d_zoom_UINS::fUILinkPTXVolPosGetOptionsNew);
    m_PTXVolPosSBCS.registerElmGetOptionsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_POSITION2_SBCS, Ep2d_zoom_UINS::fUILinkPTXVolPos_SBCSGetOptionsNew);
    m_PTXVolOri.registerElmGetOptionsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_ORIENTATION2, Ep2d_zoom_UINS::fUILinkPTXVolOriGetOptionsNew);
    m_PTXVolOriHistory.registerElmGetOptionsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_ORI_DESCR, Ep2d_zoom_UINS::fUILinkPTXVolOriHistoryGetOptionsNew);
#ifdef ZOOM_EXTENDED
    m_PTXVolVisibility.registerElmGetOptionsHandler(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_VISIBILITY, Ep2d_zoom_UINS::fUILinkPTXVolVisibilityGetOptionsNew);
#endif // ZOOM_EXTENDED

    // TR solve handler for AdjVolCoupling
    m_AdjVolCoupling.registerSolveHandler(rSeqLim, MR_TAG_ADJ_VOL_COUPLE_TO, EpCommonUINS::fEPIStdSolveSelection);

    // update PTXPulse fields / make them non-zero-sized with EXCITATION_ZOOMED
    m_PTXPulseArray.registerSizeHandler(rSeqLim, MR_TAG_PTX_PULSES_ARRAY, Ep2d_zoom_UINS::fUILinkPTXPulseArraySizeNew);
    m_PTXPulseTxAcc.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_PULSES_ARRAY, MR_TAG_PTX_TX_ACC, Ep2d_zoom_UINS::fUILinkPTXPulseTxAccGetLimitsNew);
    m_PTXPulseTxAcc.registerSetValueHandler(rSeqLim, MR_TAG_PTX_TX_ACC, Ep2d_zoom_UINS::fUILinkPTXPulseTxAccSetValueNew);
    m_PTXPulseTxAcc.registerElmSolveHandler(rSeqLim, MR_TAG_PTX_PULSES_ARRAY, MR_TAG_PTX_TX_ACC, EpCommonUINS::_solveTETITR_TI);
#ifdef ZOOM_EXTENDED
    m_PTXPulseTxAcc.registerIsAvailableHandler(rSeqLim, MR_TAG_PTX_TX_ACC, Ep2d_zoom_UINS::fUILinkPTXPulseTxAccIsAvailableNew);
#endif // ZOOM_EXTENDED

    // make PTXPulse fields NOT editable
    m_PTXPulseFlipAngle.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_PULSES_ARRAY, MR_TAG_PTX_FLIP_ANGLE, Ep2d_zoom_UINS::fUILinkPTXPulseFlipAngleGetLimitsNew);
    m_PTXPulsePhaseFOE.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_PULSES_ARRAY, MR_TAG_PTX_FOE_PHASE, Ep2d_zoom_UINS::fUILinkPTXPulsePhaseFoEGetLimitsNew);
    m_PTXPulsePhaseMatrix.registerElmGetLimitsHandler(rSeqLim, MR_TAG_PTX_PULSES_ARRAY, MR_TAG_PTX_TX_MATRIX_PHASE, Ep2d_zoom_UINS::fUILinkPTXPulsePhaseMatrixSizeGetLimitsNew);

    //-----------------------------------------------------------------------------------
    // register additional ZOOM solve handlers 
    //-----------------------------------------------------------------------------------
    Ep2d_zoom_UINS::fEPIRegisterZoomHandlers(rSeqLim);

    //-----------------------------------------------------------------------------------
    // for the rotated trajectory, the sidelobe distance (and therefore the pTx volume) are automatically chosen
    // therefore, the following UI elements are not required
    //-----------------------------------------------------------------------------------
    
#ifdef ZOOM_EXTENDED
    m_PTXVolRot.registerIsAvailableHandler(rSeqLim, MR_TAG_ROTATION, Ep2d_zoom_UINS::fUILinkPTXVolRotIsAvailableNew);
    m_PTXVolPDim.registerIsAvailableHandler(rSeqLim, MR_TAG_PDIM, Ep2d_zoom_UINS::fUILinkPTXVolPDimIsAvailableNew);
    m_PTXVolRDim.registerIsAvailableHandler(rSeqLim, MR_TAG_RDIM, Ep2d_zoom_UINS::fUILinkPTXVolRDimIsAvailableNew);
    m_PTXVolSDim.registerIsAvailableHandler(rSeqLim, MR_TAG_SDIM, Ep2d_zoom_UINS::fUILinkPTXVolSDimIsAvailableNew);
    m_PTXVolPos.registerIsAvailableHandler(rSeqLim, MR_TAG_POSITION2, Ep2d_zoom_UINS::fUILinkPTXVolPosIsAvailableNew);
    m_PTXVolOri.registerIsAvailableHandler(rSeqLim, MR_TAG_ORIENTATION2, Ep2d_zoom_UINS::fUILinkPTXVolOriIsAvailableNew);
    m_PTXVolVisibility.registerIsAvailableHandler(rSeqLim, MR_TAG_VOL_VISIBILITY, Ep2d_zoom_UINS::fUILinkPTXVolVisibilityIsAvailableNew);
#endif // ZOOM_EXTENDED

#endif  // ZOOM_2DRF

    //-------------------------------------------------------------------------------------
    // register modified is-available handler for "Save uncombined" (to remove checkbox)
    //-------------------------------------------------------------------------------------
#ifdef EPI_DISABLE_SAVE_UNCOMBINED

    LINK_BOOL_TYPE *_pSaveUnCo = _search< LINK_BOOL_TYPE >(rSeqLim, MR_TAG_SAVE_UNCOMBINED);

    if(_pSaveUnCo)
    {
        _pSaveUnCo->registerIsAvailableHandler(EpCommonUINS::fEPISaveUncombinedIsAvailable);
    }

#endif

#endif // #ifdef WIN32

    return ( MRI_SEQ_SEQU_NORMAL );
}

//  --------------------------------------------------------------------------
//
//  Name        : EpCommonUI::initializeUI
//
//  Description :
/// \brief        This function initializes UI functions and members
///
///               It can be executed on the measuement system, too, but is empty there.
///
///               It returns an NLS status
///

#ifdef WIN32
NLS_STATUS EpCommonUI::initializeUI(MrProt &rMrProt, SeqLim &rSeqLim)
#else
NLS_STATUS EpCommonUI::initializeUI(MrProt& /* rMrProt */, SeqLim& /* rSeqLim */)
#endif
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

#ifdef WIN32

#ifdef ZOOM_2DRF
    lStatus = Ep2d_zoom_UINS::initializeZoomUI(rMrProt, rSeqLim);
#endif  // ZOOM_2DRF

#endif // #ifdef WIN32

    return (lStatus);
}


// ------------------------------------------------------------------------------
// Function    : updateB1ControlLoopParametersInProtocol()
// ------------------------------------------------------------------------------
// Description : Updates the B1 control loop parameters in the protocol. The decision
//               if and how the B1 control loop is activated is carried out in the
//               function "useB1ControlLoop", which is called by this method.
//               
// Return      : true in case of success, false otherwise
//
// ------------------------------------------------------------------------------
#ifdef WIN32
bool EpCommonUINS::updateB1ControlLoopParametersInProtocol(MrProt &rMrProt)
{
    // local variables
    bool bSuccess               = true;

    // the B1 control loop so far does not exist for the segmented epi sequnces
#if !defined SEQUENCE_CLASS_EP_SEG && !defined SEQUENCE_CLASS_AslCsl   // only for COMPILE_EP2D_DIFF // SUPPORT_iPAT_TGSE

    float fCorrectionFactorMax  = 1.0f;
    float fPeakReserveFactor    = 0.0f;

    // get B1 control loop structure 
    MrTXSpecData &rTXSpec = rMrProt.getsTXSPEC();
    B1CorrectionParametersType &rB1CorrectionParameters = rTXSpec.getB1CorrectionParameters();

    // retrieve information if and how the control loop is used  
    rB1CorrectionParameters.setbActive(useB1ControlLoop(rMrProt, fCorrectionFactorMax, fPeakReserveFactor));
    rB1CorrectionParameters.setflCorrectionFactorMax(fCorrectionFactorMax);
    rB1CorrectionParameters.setflPeakReserveFactor(fPeakReserveFactor);
    rB1CorrectionParameters.setbValid(true); // this flag has to be set true, otherwise the values above will be set to the values from the SeqLims later on
#endif //#ifndef SEQUENCE_CLASS_EP_SEG 

    // return
    return bSuccess;
}


// ------------------------------------------------------------------------------
// handler for slice acceleration
// ------------------------------------------------------------------------------


int32_t EpCommonUINS::fUILinkSliceAccelSetValue(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    long lResult = getUI(pThis)->EpCommonUI::m_AccelerationFactorSlice.getOrigSetValueHandler()(pThis, lNewVal, lPos);

    // set recommended FOV shift factor
    // In the current implementation the UI element for the FOV shift factor is hidden thus the factor is written to the protocol directly.
    // If the UI parameter is supposed to be used again (see fUILinkFOVShiftFactorIsAvailable) use the line which is currently commented out.
    // We need to specify if this is a diffusion sequence as we get different values for 1.5T there
    MrProt	rMrProt(&pThis->prot());

    SMSProperties::setFOVShiftFactor(rMrProt, SMSProperties::getRecommendedFOVShiftFactorCumulative(pThis->prot(), getUI(pThis)->isSpinEcho()));

    //setProtocolParameter<LINK_LONG_TYPE>(pThis, MR_TAG_SLICE_ACCEL_FOV_SHIFT_FACTOR, SMSProperties::getRecommendedFOVShiftFactorCumulative(pThis->prot()), 0, false, false);

    // For acceleration factors >= 3 we recommend to use low SAR pulses
    if (getUI(pThis)->isSpinEcho())
    {
        if (lNewVal >= 3)
        {
            setProtocolParameter<LINK_SELECTION_TYPE>(
                pThis, MR_TAG_RFPULSE_TYPE, MRI_STD_RFPULSE_TYPE_LOW_SAR, 0, true, false);
        }
    }

    // Check for high channel configuration and activate channel compression in case of SMS
    const MrRxCoilSelect coilSelect(rMrProt.coilInfo().Meas().getaRxCoilSelectData()[0]);
    if (coilSelect.getNumOfUsedADCChan() > 32 && lNewVal > 1)
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_CHANNEL_MIXING_MODE, MRI_STD_MATRIX_OPT_PERFORMANCE, 0, true, false);

    return lResult;
}


unsigned EpCommonUINS::fUILinkSliceAccelSolve(LINK_LONG_TYPE* const pThis, char** arg_list, const void* pAddMem, const MrProtocolData::MrProtData * pOrigProt, int32_t lIndex)
{
    return solveSliceAccelerationSettings<LINK_LONG_TYPE>(pThis);
}


bool EpCommonUINS::fUILinkFOVShiftFactorGetLimits(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t lPos)
{
    rulVerify = LINK_LONG_TYPE::VERIFY_SCAN_ALL;

    rLimitVector.clear();
    MrProt rMrProt(&pThis->prot());

    MrProtFacade protFacade(rMrProt);
    MrLimitLong sLimInterval;
    if(protFacade.isSliceAcceleration())
    {
        //  Ice reconstruction currently supports only factors up to 4
        sLimInterval.setEqualSpaced(1, 4 /*SMSProperties::MAX_FOV_SHIFT_FACTOR*/, 1, MrLimit<int32_t>::GREEN);
        rLimitVector.push_back(sLimInterval);
    }
    else
    {
        sLimInterval.setEqualSpaced(1, 1, 1, MrLimit<int32_t>::GREEN);
        rLimitVector.push_back(sLimInterval);
    }
    return true;
}


bool EpCommonUINS::fUILinkFOVShiftFactorIsAvailable(LINK_LONG_TYPE* const pThis, int32_t lPos)
{
    // the FOV shift factor parameter is not available in product; it can be made available by setting the INI file parameter
    if (SysProperties::ReadSeqSettingGeneral("EPI_GENERAL/show_FOV_shift_factor", false, true))
    {
        return getUI(pThis)->EpCommonUI::m_FOVShiftFactor.getOrigIsAvailableHandler()(pThis, lPos);
    }

    return false;
}


bool EpCommonUINS::fUILinkNumberOfSlicesGetLimits(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t lPos)
{
    // Call base class implementation (which already considers basic slice-acceleration related restrictions)
    bool bResult = getUI(pThis)->EpCommonUI::m_NumberOfSlices.getOrigGetLimitsHandler()(pThis, rLimitVector, rulVerify, lPos);

    if(bResult)
    {
        MrProt rMrProt(&pThis->prot());
        MrProtFacade protFacade(rMrProt);

        if ( protFacade.isSliceAcceleration() )
        {
            long lMultibandFactor = SMSProperties::getMultiBandFactor( rMrProt );

            if ( ( lMultibandFactor > 1 ) && ( rMrProt.sliceSeries().mode() == SEQ::INTERLEAVED ) )
            {
                if ( ( rMrProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE ) || ( rMrProt.getsPrepPulses().getucSatRecovery() == SEQ::SATREC_SLICE_SELECTIVE ) )
                {
                    // For interleaved, slice-selective preparations, multiple slices will tyically get prepared successively.
                    // In combination with SMS, crosstalk of virtually adjacent slices might become an issue. Note that this
                    // can only happen for exactly twofold interleaved preparations.

                    // Example:  15 slices, multi-band factor 3, 2 IR-interleaves
                    //           Interleave #1:  Slice indices                  0)      1)      2)
                    //                           Corresponding slice numbers    1,6,11  3,8,13  5,10,15   => inversion of adjacent slices 5 and 6!
                    //           Interleave #2:  Slice indices                  3)      4)
                    //                           Corresponding slice numbers    2,7,12  4,9,14

                    // Possible solutions:
                    // a) Prohibit twofold interleaved preparations.              Drawback: less effective protocols.
                    // b) Add 'hidden' slices in case of need.                    Drawback: complex implementation, affecting various modules.
                    // c) Limit slice numbers to even multiples of the MB-factor. Drawback: unneccessary limitations for >twofold interleaved preparations.
                    // => For the moment, implement solution c)
                    lMultibandFactor *= 2;
                }
            }

            std::vector<MrLimitLong> rLimitVectorNew;

            for ( size_t iIndex = 0; iIndex < rLimitVector.size(); ++iIndex )
            {
                for (long lNumberOfSlices = rLimitVector[iIndex].minimum(); lNumberOfSlices <= rLimitVector[iIndex].maximum(); lNumberOfSlices += rLimitVector[iIndex].incr())
                {
                    if ((lNumberOfSlices % lMultibandFactor) == 0)
                    {
                        MrLimitLong rLonelyLimit;
                        rLonelyLimit.setLonely(lNumberOfSlices);
                        rLimitVectorNew.push_back(rLonelyLimit);
                    }
                }
            }

            rLimitVector = rLimitVectorNew;
        }
    }

    return bResult;
}


unsigned EpCommonUINS::fUILinkPATAccelFactorGetToolTipID(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t lPos)
{
    return MRI_STD_ACCELERATION_FACTOR_PE_TOOLTIP;
}


unsigned EpCommonUINS::fUILinkFatSuppresionSolve(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lPos)
{
    unsigned uResult = EpCommonUINS::fEPIStdSolveSelection(pThis, arg_list, pVoid, pOrigProt, lPos);

    if(uResult == 0)
    {
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_PAT_MODE, MRI_STD_PAT_MODE_GRAPPA, 0);
        if(pThis->sequence().prepareForBinarySearch(&pThis->prot()))
            uResult = MRI_STD_CONFIRMATION_MSG;

    }

    return uResult;
}


bool EpCommonUINS::fUILinkCoilElementsSetValue(LINK_BOOL_TYPE* const pThis, bool bDesiredState, int32_t lIndex)
{
    long lResult = getUI(pThis)->m_CoilElements.getOrigSetValueHandler()(pThis, bDesiredState, lIndex);

    MrProt rMrProt(&pThis->prot());
    const MrRxCoilSelect rCoilSelect(rMrProt.coilInfo().Meas().getaRxCoilSelectData()[0]);

    if (!IterativeDenoisingUIParameter::areSwiftBrainCoilSettingsSatisfied(rMrProt))
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_ADVANCED_RECONSTRUCTION_MODE, MRI_STD_OFF);

    if (strcmp("Body", rCoilSelect.element(0).getCoilElementID().gettCoilID().c_str()) == 0)
    {
        setProtocolParameter<LINK_LONG_TYPE>(pThis, MR_TAG_SLICE_ACCEL_FACTOR, 1, 0);
    }

    return lResult;
}

bool EpCommonUINS::fUILinkFreezeSuppressedTissueSetValue( LINK_BOOL_TYPE* const pThis, bool bDesiredState, int32_t lIndex )
{
    bool bResult = getUI( pThis )->m_FreezeSuppressedTissue.getOrigSetValueHandler()( pThis, bDesiredState, lIndex );

    if ( bDesiredState == false )
    {
        // When disabling freezing of suppressed tissue: put TI on UI raster
        int32_t lTImin = pThis->seqLimits().getTI()[0].getMin();
        int32_t lTImax = pThis->seqLimits().getTI()[0].getMax();
        int32_t lTIinc = pThis->seqLimits().getTI()[0].getInc();
        int32_t lTI    = pThis->prot().getalTI()[0];

        // NOTE: We assume something about TI-limit-handling in MrUILink!
        if ( lTI >   10000 )
        {
            lTIinc *= 10;
        }
        if ( lTI > 1000000 )
        {
            lTIinc *= 10;
        }

        int32_t lTInew = shift_to_grid_low( lTI, lTImin, lTImax, lTIinc );

        pThis->prot().getalTI()[0] = lTInew;
    }

    return bResult;
}

#ifdef SUPPORT_iPAT_a_ep2d
bool SEQ_NAMESPACE::EpCommonUINS::fEPIRefLinesPEIsAvailable(LINK_LONG_TYPE* const pThis, int32_t /*pos*/)
{
    MrProt rMrProt(&pThis->prot());

    if(rMrProt.PAT().getlAccelFactPE() < 2)
        return false;
    else
        return true;
}
#endif

unsigned SEQ_NAMESPACE::EpCommonUINS::fUILinkAdjustmentModeSolve(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lPos)
{
    unsigned uResult = EpCommonUINS::fEPIStdSolveSelection(pThis, arg_list, pVoid, pOrigProt, lPos);

    return uResult;
}

unsigned SEQ_NAMESPACE::EpCommonUINS::fUILinkAdjustmentModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    MrProt rMrProt(&pThis->prot());
    MrProtFacade protFacade(rMrProt);

    switch(uNewVal)
    {
        case MRI_STD_ADJUSTMENT_MODE_STANDARD:
            // when switching from SliceAdjust to standard 
            if (protFacade.isSliceAdj())
            {
                if (!SysProperties::isUHFSystem())
                {
                    if (rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_WHOLE_BODY
                        || rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_TUNEUP)
                    {
                        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_SHIM_MODE, MRI_STD_SHIM_MODE_STANDARD);
                    }
                }
            }
            break;

        case MRI_STD_ADJUSTMENT_MODE_SLICE_BY_SLICE:
            if (!SysProperties::isUHFSystem())
            {
                // set B0 shim to mode whole body
                if (rMrProt.getsAdjData().getuiAdjShimMode() != SEQ::ADJSHIM_ABSOLUTE)
                    setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_SHIM_MODE, MRI_STD_SHIM_MODE_WHOLE_BODY);
            }

            // disable ZOOMit
            setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_EXCIT_PULSE, MRI_STD_EXCITATION_STANDARD, 0);
            
            // turn off SMS for SliceAdjust
            if(protFacade.isSliceAcceleration())
            {
                setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_PAT_MODE, MRI_STD_PAT_MODE_GRAPPA, 0);                
            }

            // turn local shim off
            if(rMrProt.getsAdjData().getuiLocalShim() != LocalShim_OFF)
            {
                setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_LOCAL_SHIM, MRI_STD_LOCAL_SHIM_OFF, 0);
            }
            break;

        case MRI_STD_ADJUSTMENT_MODE_FAST_VIEW:  
            break;
    }

    return getUI(pThis)->m_AdjustmentMode.getOrigSetValueHandler()(pThis, uNewVal, lPos);;
}

bool SEQ_NAMESPACE::EpCommonUINS::fUILinkLocalShimTry(LINK_SELECTION_TYPE* const pThis, void* pVoid, const MrProtocolData::MrProtData* pOrig, int32_t pos)
{
    MrProt rMrProt(&pThis->prot());
    MrProtFacade protFacade(rMrProt);

    if(protFacade.isSliceAdj())
    {
        // B0 shim mode 'whole body' and coil shim are not allowed at the same time
        if(rMrProt.getsAdjData().getuiAdjShimMode() == SEQ::ADJSHIM_WHOLE_BODY && rMrProt.getsAdjData().getuiLocalShim() != LocalShim_OFF)
            return false;
    }

    // call original handler
    return getUI(pThis)->m_LocalShim.getOrigTryHandler()(pThis, pVoid, pOrig, pos);
}

unsigned SEQ_NAMESPACE::EpCommonUINS::fDistortionCorrectionSetValue(
    LINK_SELECTION_TYPE* const pThis, unsigned uNewValue, int32_t lIndex)
{
    if ((uNewValue == MRI_STD_OFF) && (pThis->prot().getucStaticFieldCorrection() != 0))
    {
        // Static field correction is enabled and user wishes to disable distortion correction
        // => Disable static field correction
        setProtocolParameter<LINK_BOOL_TYPE>(pThis, MR_TAG_STATIC_FIELD_CORRECTION, false);
    }

    return (getUI(pThis)->EpCommonUI::m_DistortionCorrection.getOrigSetValueHandler()(pThis, uNewValue, lIndex));
}

bool SEQ_NAMESPACE::EpCommonUINS::fUILinkSFCSetValue(LINK_BOOL_TYPE* const pThis, bool bDesiredState, int32_t lIndex)
{
    if (bDesiredState == true)
    {
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_SHIM_MODE, MRI_STD_SHIM_MODE_ABSOLUTE);
        setProtocolParameterElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_SLICE_GROUP_LIST, MR_TAG_SG_PE_DIR, MRI_STD_POSTERIOR_TO_ANTERIOR);
    }
    else
    {
        setProtocolParameterElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_SLICE_GROUP_LIST, MR_TAG_SG_PE_DIR, MRI_STD_ANTERIOR_TO_POSTERIOR);
    }
    
    if (pThis->prot().getsDistortionCorrFilter().getucMode() == SEQ::DISTCORR_NDIS)
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_FLT_DISCOR_PROP_MODE, MRI_STD_DISTCORR_2D);

    return getUI(pThis)->m_SFC.getOrigSetValueHandler()(pThis, bDesiredState, lIndex);

}

#endif // #ifdef WIN32
