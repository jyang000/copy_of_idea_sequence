//-----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 1999  All Rights Reserved.  Confidential
//-----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//
//    File: \n4_servers1\pkg\MrServers\MrImaging\seq\common\iPAT\iPAT.cpp
//
//  Author: NITTMA69
//          
//    Date: 2015-02-11 10:12:23 +01:00
//
//    Lang: C++
//
// Descrip: Implements common functionality for iPAT sequences.
//
//-----------------------------------------------------------------------------

#include "MrProtSrv/Domain/CoreNative/MrPat.h"
#include "MrProtSrv/Domain/CoreNative/MrFastImaging.h"
#include "MrProtSrv/Domain/CoreNative/MrFilter.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/KSpace/MrKSpace.h"
#include "MrMeasSrv/SeqIF/Sequence/sequmsg.h"
#include "MrImaging/seq/common/iPAT/iPAT.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

#ifdef WIN32
    #include    <vector>
    #include    "MrImagingFW/libSeqSysProp/SysProperties.h"
    #include    "MrProtSrv/Domain/MrProtocol/UILink/MrStdNameTags.h"
    
    // For PAT averaging / MRI_STD... macros
    #include    "MrProtSrv/Domain/MrProtocol/StdProtRes/StdProtRes.h"
#endif

#ifdef SEQUENCE_CLASS
    #include "MrMeasSrv/SeqIF/Sequence/Sequence.h"            // Sequence. Required for getPointerToUIHandlersPAT()
#endif

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif

namespace SEQ_NAMESPACE
{


// 
// class UIHandlersPAT
//
UIHandlersPAT::UIHandlersPAT()
#ifdef WIN32
        : m_pPAT_OrigSetValue_TurboFactor           (nullptr)
        , m_pPAT_OrigSetValue_Segments              (nullptr)
        , m_pPAT_OrigSetValue_BaseResolution        (nullptr)
        , m_pPAT_OrigSetValue_PATAccelFActorPE      (nullptr)
        , m_pPAT_OrigSetValue_ReadFOV               (nullptr)
        , m_pPAT_OrigSetValue_PhaseFOV              (nullptr)
        , m_pPAT_OrigSetValue_PhaseOversampling     (nullptr)
        , m_pPAT_OrigSetValue_PhaseResolution       (nullptr)
        , m_pPAT_OrigSetValue_PhasePartialFourier   (nullptr)
        , m_pPAT_OrigSetValue_PATMode               (nullptr)
        // For PATaveraging
        , m_pPAT_OrigSetValue_RefScanMode           (nullptr)
        , m_pOrigSetValue_Averages                  (nullptr)
        // - end PAT averaging
        , m_pPAT_OrigGetLimits_PATRefLinesPE        (nullptr)
        , m_lUILinkPATRefLinesOpt                   (0)
        , m_pReorder                                (nullptr)
#endif
    {};





#ifdef WIN32

// Auxiliary function: translate sequence PatRefScanMode values to UILink "uVal" values
// - see MrUILinkPat::fUILinkRefScanModeGetValue
// - \ToDo: check whether this can be moved to, say, MrUILinkBase (and be exported)
unsigned _SeqRefScanMode2UILinkVal ( SEQ::PATRefScanMode seqRefScanMode )
{
    unsigned uRet = MRI_STD_EMPTY;
    switch ( seqRefScanMode )
    {
    case SEQ::PAT_REF_SCAN_UNDEFINED:
        uRet = MRI_STD_PAT_REF_SCAN_UNDEFINED;
        break;
    case SEQ::PAT_REF_SCAN_INPLACE:
        uRet = MRI_STD_PAT_REF_SCAN_INPLACE;
        break;
    case SEQ::PAT_REF_SCAN_EXTRA:
        uRet = MRI_STD_PAT_REF_SCAN_EXTRA;
        break;
    case SEQ::PAT_REF_SCAN_EXTRA_EPI:
        uRet = MRI_STD_PAT_REF_SCAN_EXTRA_EPI;
        break;
    case SEQ::PAT_REF_SCAN_EXTRA_TSE:
        uRet = MRI_STD_PAT_REF_SCAN_EXTRA_TSE;
        break;
    case SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST:
        uRet = MRI_STD_PAT_REF_SCAN_EXTRA_GRE_FAST;
        break;
    case SEQ::PAT_REF_SCAN_PRESCAN:
        uRet = MRI_STD_PAT_REF_SCAN_PRESCAN;
        break;
    case SEQ::PAT_REF_SCAN_INTRINSIC_AVE:
        uRet = MRI_STD_PAT_REF_SCAN_INTRINSIC_AVE;
        break;
    case SEQ::PAT_REF_SCAN_INTRINSIC_REP:
        uRet = MRI_STD_PAT_REF_SCAN_INTRINSIC_REP;
        break;
    case SEQ::PAT_REF_SCAN_INTRINSIC_PHS:
        uRet = MRI_STD_PAT_REF_SCAN_INTRINSIC_PHS;
        break;
    case SEQ::PAT_REF_SCAN_INPLACE_LET:
        uRet = MRI_STD_PAT_REF_SCAN_INPLACE_LET;
        break;
    default:
        assert(0);
        uRet = MRI_STD_EMPTY;
    } // switch

    return uRet;
}

// more auxiliary functions: _extraPhOSMax and _checkPhOSAndPEFTLen
//   from MrUILinkKSpace are not accessible from the outside (they
//   are not exported). We need this functionality here to check
//   whether a change in TurboFactor will help conflicts when (auto-
//   matically) switching the PatRefScanMode, see below.
// \ToDo: check whether this can be exported from MrUILinkKSpace
//   - 20060419: still old situation (not exported), but code is still
//               the same as in \n4_servers1\pkg\MrServers\MrProtSrv\MrProtocol\UILink\Config.cpp
// \ToDo: find a way if possible to avoid duplicating UI functionality here
uint32_t _iPATextraPhOSMax(const MrUILinkBase* const pThis)
{
    static uint32_t dwExtraPhOSMax = 10;
    uint32_t dwExtraPhOSMax_SEG = pThis->seqLimits().getPhaseOSMax().getDef();

    switch (pThis->seqLimits().getPELines().getIncMode())
    {
    case SEQ::INC_TURBO_FACTOR:
    case SEQ::INC_TURBO_FACTOR_IPAT:
    case SEQ::INC_EPI_FACTOR:
    case SEQ::INC_EPI_FACTOR_IPAT:
    case SEQ::INC_GRE_SEGMENTS:
    case SEQ::INC_GRE_SEGMENTS_IPAT:
    case SEQ::INC_TGSE_FACTOR:
    case SEQ::INC_TGSE_FACTOR_IPAT:
    case SEQ::INC_SEGMENTED:
    case SEQ::INC_SEGMENTED_IPAT:
    case SEQ::INC_TWICE_EPI_FACTOR:
    case SEQ::INC_TWICE_EPI_FACTOR_IPAT:
        return dwExtraPhOSMax_SEG;
    default:
        return dwExtraPhOSMax;
    }
}

// Cannot access this in MrUiLinkKSpace - see above
// - we restrict the copy operation to the minimum, so that
//   we don't have to copy too many functions.
bool _iPATcheckPhOSAndPEFTLen(const MrUILinkBase* const pThis)
{
    const ParLim<int32_t>& rPELimits = pThis->seqLimits().getPELines();
    if(rPELimits.isAvailable())
    {
        MrProt rMrProt(&pThis->prot());
        const double dOS = rMrProt.phaseOversampling();

        // Really only worry about excessive phase oversampling
        double dAcceptedExtraOS = _iPATextraPhOSMax(pThis)/100.;
        return (dOS >= 0.)&&(dOS-rMrProt.kSpace().phaseOversamplingForDialog()) <= dAcceptedExtraOS;
    }
    return true;
}
// END HACK


// Auxiliary function for switching PAT averaging on and off
void _setPATaveraging( MrUILinkBase * const pThis, bool bPatAveragingIsPossible )
{
    // Current refScanMode
    SEQ::PATRefScanMode currRefScanMode    = pThis->prot().getsPat().getucRefScanMode();

    // Default refScanMode
    SEQ::PATRefScanMode defaultRefScanMode = pThis->seqLimits().getRefScanMode().getDef();
    // - translate this to MRI_STD... macros
    //   - see MrUILinkPat::fUILinkRefScanModeGetValue
    unsigned nDefaultRefScanID             = _SeqRefScanMode2UILinkVal(defaultRefScanMode);
    unsigned nIntrinsicAveRefScanID        = _SeqRefScanMode2UILinkVal(SEQ::PAT_REF_SCAN_INTRINSIC_AVE);

    // Is PAT_REF_SCAN_INTRINSIC_AVE a sequence limits option
    // \ToDo: if sequence limits prohibit PAT averaging, we should not be checking that
    //        every time here, but rather have a global flag (for performance reasons)
    bool bPatAveragingIsAllowedBySequence  = pThis->seqLimits().getRefScanMode().hasOption(SEQ::PAT_REF_SCAN_INTRINSIC_AVE);
    
    // Try to find the setValue handler for the PatRefScanMode
    // - still not decided what to do if the refScanMode is not available... but it should since
    //   this function is only called if the seqLimits allow INTRINSIC_AVE
    LINK_SELECTION_TYPE* pRefScanMode      = 
        _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_PAT_REF_SCAN_MODE, 0, MrUILinkBase::SEARCH_MODE::SEARCH_EDITABLE);
    
    // Also: turboFactor
    LINK_LONG_TYPE*      pTurboFactor      = 
        _search<LINK_LONG_TYPE>(pThis, MR_TAG_TURBO_FACTOR, 0, MrUILinkBase::SEARCH_MODE::SEARCH_EDITABLE);
        
    // switch on or off according to the following logic:

    // PAT_REF_SCAN_INPLACE && PATaverage possible -> switch ON
    if ((currRefScanMode == SEQ::PAT_REF_SCAN_INPLACE || currRefScanMode == SEQ::PAT_REF_SCAN_EXTRA_TSE || currRefScanMode == SEQ::PAT_REF_SCAN_EXTRA || currRefScanMode == SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST) && bPatAveragingIsPossible && bPatAveragingIsAllowedBySequence)
    {
        
        // Make that user-selectable (we may lose some occasions where
        //   PATaveraging might be switched on without the detour via INPLACE)
        if ( pRefScanMode == nullptr )
        {
            // Parameter is not available (editable)
            // - send a message
            SEQ_TRACE_ALWAYS.print("PatRefScanMode UI element not accessible (-> SEQ::PAT_REF_SCAN_INTRINSIC_AVE)!");

        } else {

            // Use setValue handler
            // - forced!
            //   -> makes, e.g., higher average values "red"
            pRefScanMode->value( nIntrinsicAveRefScanID, 0, MrUILinkBase::SET_MODE::SET_FORCED );
            // Add the turboFactor as dependent parameter
            // - if we add the refScanMode itself it will appear in the checkbox also
            //   (we might actually want that)
            // - we know that the setValue handler for the refScanMode will always try
            //   to use the turbo factor to solve
            // - only "solve" if not in shots mode - that takes care of PE line consistency itself
            MrProtocolData::MrProtData*  pMrProt = &pThis->prot();
            if( pMrProt->getsFastImaging().getucSegmentationMode() != SEQ::SEGM_MODE_DEFINE_SHOTS )
            {
                // Note: this may have already been done in pRefScanMode->value() and would then be obsolete
                //       here.
                if ( pTurboFactor != nullptr )
                    pThis->addDependentParamPtr(pTurboFactor,0);
            }

            // Also set refLines to a suitable value
            int32_t lRefLines = lPATDefOptNoRefLinesIntrinsic;
            const int32_t lPELin = pMrProt->getsKSpace().getlPhaseEncodingLines();
            if ( lRefLines > lPELin )
            {
                lRefLines = lPELin;
            }
            pMrProt->getsPat().setlRefLinesPE(lRefLines);
            
        } // if ( pRefScanMode == nullptr )
        
    }

    // PAT_REF_SCAN_INTRINSIC_AVE && PATaverage *NOT* possible -> switch to default
    if ( (currRefScanMode == SEQ::PAT_REF_SCAN_INTRINSIC_AVE) && !bPatAveragingIsPossible )
    {
        
        // see above
        if ( pRefScanMode == nullptr )
        {
            SEQ_TRACE_ALWAYS.print("PatRefScanMode UI element not accessible (-> default RefScanMode)!");

        } else {

            // Use setValue handler
            // - forced.
            pRefScanMode->value( nDefaultRefScanID, 0, MrUILinkBase::SET_MODE::SET_FORCED );
            // - any dependent parameters necessary
            //   (up to now no case found where INPLACE allows fewer settings than INTRINSIC)
            // - found some cases where switching back on the AF causes red zones.

            // - may still be necessary, because we cannot go back to low numbers of averages in some cases
            MrProtocolData::MrProtData*  pMrProt = &pThis->prot();
            if( pMrProt->getsFastImaging().getucSegmentationMode() != SEQ::SEGM_MODE_DEFINE_SHOTS )
            {
                // Note: this may have already been done in pRefScanMode->value() and would then be obsolete
                //       here.
                if ( pTurboFactor != nullptr )
                    pThis->addDependentParamPtr(pTurboFactor,0);
            }

        } // if ( pRefScanMode == nullptr )
        
    }
    
} // void _setPATaveraging(...)

#endif
//-------------------------------------------------------------------------------------
// global static variables
// 
// NOTE:
// - Due to this global instance of class UIHandlersPAT this module must be compiled 
//   with every sequence that wants to use the functionality supplied.
// - myUIHandlersPAT must not be used for 'class-type' sequences (see comments in iPAT.h)
//-------------------------------------------------------------------------------------
UIHandlersPAT myUIHandlersPAT;    // stores the orig. value handlers, pointer to ReorderInfo and optimal number of ref.lines







// ------------------------------------------------------------------------------
// Function    : fPATCheckAndUpdateReferenceLineNumber
// ------------------------------------------------------------------------------
//               
// Description : Segmented sequences using PAT only work with a restricted number
//               of inplace reference lines.
//               This function is needed to find a valid number of reference 
//               lines for the current number of phase encoding lines and number
//               of segments. Due to the fact, that reordering schemes for segmented
//               sequences can become very complex, this function relies on the
//               help of a reordering class.
//
// Return      : true , if a valid number of reference lines can be found
//               false, if not
//
// ------------------------------------------------------------------------------
bool SeqIF_UIHandlersPAT::
fPATCheckAndUpdateReferenceLineNumber(MrProtocolData::MrProtData* pMrProt, SeqLim* _pSeqLim, ReorderInfo* pRI, int32_t lInpPATRefLinesOpt)
{
#ifdef WIN32
    int32_t lPPARefLinesOrig    = pMrProt->getsPat().getlRefLinesPE();
    int32_t lPPARefLinesCurrent = lPPARefLinesOrig;
    int32_t lPPARefLinesOpt     = (lInpPATRefLinesOpt < 0 ? lPATDefOptNoRefLines : lInpPATRefLinesOpt);

    MrProt sMrProt(pMrProt);

    // nothing to do, if PAT is not active
    if (pMrProt->getsPat().getucPATMode()==SEQ::PAT_MODE_NONE)
    {
        if (lPPARefLinesOrig)
        {
            pMrProt->getsPat().setlRefLinesPE(0);
        }
        return true;
    }
    
    // can do nothing, if pRI was not set
    if (pRI==nullptr)
    {
        return false;
    }

    SeqIF_UIHandlersPAT::PATInitializeReorderMethod( pMrProt, _pSeqLim );

    if( pMrProt->getsKSpace().getucTrajectory() == SEQ::TRAJECTORY_BLADE )
    {
        //  Check reference lines but do not modify them
        return pRI->prepareCalculation(sMrProt,*_pSeqLim);
    } 

    // PATaveraging: need to check and update refLines if previous boundary conditions
    //   enforced a lower number than optimal
    if (pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_INTRINSIC_AVE)
    {

        // At this point the reorderInfo class has not been prepared yet!
        // -> worst case: run it to know boundaries (number of specified
        //    refLines is irrelevant for non-INPLACE modes
        (void) pRI->prepareCalculation(sMrProt,*_pSeqLim);

        // Set number of refLines to global default
        int32_t lRefLines = lPATDefOptNoRefLinesIntrinsic;

        // Find out how many symmetric k-space lines there are
        int32_t lMinPatLinIdx = pRI->getMinPATLinNo();
        int32_t lCtrPatLinIdx = pRI->getKSCenterLin();
        // No frontend for max pat line number!
        // - calculate it from the center line and the max. line number
        //int32_t lMaxPatLinIdx = pRI->getMaxPATLinNo();
        int32_t lLastLineIgnoringPAT = pRI->getMaxLineNumber();
        int32_t lAFLin        = pMrProt->getsPat().getlAccelFactPE();
        int32_t lMaxPatLinIdx = lCtrPatLinIdx + lAFLin*((lLastLineIgnoringPAT-lCtrPatLinIdx)/lAFLin);

        // Symmetric part of the center k-space lines
        int32_t lMaxIntrinsicRefLines = 2*(lMaxPatLinIdx-lCtrPatLinIdx+1);
        if ( lCtrPatLinIdx-lMinPatLinIdx < lMaxIntrinsicRefLines )
        {
            lMaxIntrinsicRefLines = lCtrPatLinIdx-lMinPatLinIdx;
        }

        // Reduce refLines setting if not enough symmetric lines (e.g., in partial Fourier)
        if ( lRefLines > lMaxIntrinsicRefLines )
        {
            lRefLines = lMaxIntrinsicRefLines;
        }

        // Set number of refLines in the protocol
        pMrProt->getsPat().setlRefLinesPE(lRefLines);

    }

    // set min. PPARefLines:
    // Note: there is a higher limit for EXTRA refScans in prepPost below
    // - this may be an unnecessary restriction for PatAveraging (see above), but we keep it for
    //   now because a certain minimum nuber of refLines is good for image quality
    int32_t lPPARefLinesMin = pMrProt->getsPat().getlAccelFactPE()*4;
    if(lPPARefLinesOrig < lPPARefLinesMin)
    {
        lPPARefLinesOrig = lPPARefLinesMin;
        pMrProt->getsPat().setlRefLinesPE(lPPARefLinesMin);
    }

    // nothing to do, if everything is alright already
    // - in PATaveraging mode this will most probably already succeed by now
    if (pRI->prepareCalculation(sMrProt,*_pSeqLim))
    {
        return true;
    }

    // search for lPPARefLinesOpt < lPPARefLinesCurrent < lPPARefLinesOrig
    lPPARefLinesCurrent = lPPARefLinesOrig;
    pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
    while ( ! pRI->prepareCalculation(sMrProt,*_pSeqLim) )
    {
        lPPARefLinesCurrent --;
        pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
        if(lPPARefLinesCurrent<lPPARefLinesOpt) 
        {
            lPPARefLinesCurrent = 0;
            break;
        }
    }
    if(lPPARefLinesCurrent != 0)
    {
        pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
        return true;
    }

    // search for lPPARefLinesCurrent > lPPARefLinesOrig
    lPPARefLinesCurrent = lPPARefLinesOrig;
    pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
    while ( ! pRI->prepareCalculation(sMrProt,*_pSeqLim) ) 
    {
        lPPARefLinesCurrent ++;
        pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
        if(lPPARefLinesCurrent > pRI->getMaxLineNumber()+1)
        {
            lPPARefLinesCurrent = 0;
            break;
        }
    }
    if(lPPARefLinesCurrent != 0)
    {
        pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
        return true;
    }

    // search for lPPARefLinesCurrent < lPPARefLinesOpt
    lPPARefLinesCurrent = lPPARefLinesOpt;
    pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
    while ( ! pRI->prepareCalculation(sMrProt,*_pSeqLim) )
    {
        lPPARefLinesCurrent --;
        pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
        if(lPPARefLinesCurrent<lPPARefLinesMin) 
        {
            lPPARefLinesCurrent = 0;
            break;
        }
    }
    if(lPPARefLinesCurrent != 0)
    {
        pMrProt->getsPat().setlRefLinesPE(lPPARefLinesCurrent);
        return true;
    }

#endif
    // no chance ...
    return false;

}



#ifdef WIN32
/// ------------------------------------------------------------------------------
/// Functions   : getPointerToUIHandlersPAT
/// ------------------------------------------------------------------------------
///               
/// Description : this function decides, which UIHandlersPAT object is used:
///               - 'conventional' sequences use static global myUIHandlersPAT
///               - 'class-type' sequences have to provide own instance of class UIHandlerPAT
///                 (this case selected by compiler define 'SEQUENCE_CLASS')
///               For additional information refer to comments in iPAT.h
///
/// Return      : pointer to class UIHandlersPAT
///
/// ------------------------------------------------------------------------------
template<class TYPE>
const UIHandlersPAT *getPointerToUIHandlersPAT( TYPE* pThis )
{
    // sequence class: each sequence object has to provide own UIHandlers for PAT!
    return (static_cast<SeqIF_UIHandlersPAT*>(pThis->sequence().getSeq()))->getUIHandlersPAT();
}




/// ------------------------------------------------------------------------------
/// Functions   : set value handler for segmented sequences
/// ------------------------------------------------------------------------------
///               
/// Description : provides common functionality to update number of ref.lines -
///               for segmented sequences, the number of ref.lines has to be adapted, 
///               whenever a protocol parameter changes, that influences the number of 
///               phase encoding lines to measure
///               (the UI doesn't consider iPAT when calculating lines to measure
///                -> here we (i.e. the sequence) have to ensure, that lines to measure is 
///                   an integer multiple of no.of segments)
///
/// Return      : final set value
///
/// ------------------------------------------------------------------------------
int32_t _TurboFactor_SetValue(MrUILinkLimited<int32_t>* const pThis, int32_t lNewVal, int32_t lIndex)
{
      return _PATSetValueAndUpdateRefLines(pThis, lNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigTurboFactorSetFct() );
}

int32_t _Segments_SetValue(MrUILinkLimited<int32_t>* const pThis, int32_t lNewVal, int32_t lIndex)
{
    return _PATSetValueAndUpdateRefLines(pThis, lNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigSegmentsSetFct() );
}

int32_t _BaseResolution_SetValue(MrUILinkLimited<int32_t>* const pThis, int32_t lNewVal, int32_t lIndex)
{
    return _PATSetValueAndUpdateRefLines(pThis, lNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigBaseResolutionSetFct() );
}

int32_t _PATAccelFActorPE_SetValue(MrUILinkLimited<int32_t>* const pThis, int32_t lNewVal, int32_t lIndex)
{
    // It may be necessary to use the original set-value handler(s) first, then switch PATaveraging
    // - for some reason this causes problems when the original iPAT would not allow this accelFactPE (e.g., 
    //   for very large turboFactors. Just doing it after the fact seems to work OK.
    int32_t lRet = 0;

    // PATaveraging: this may affect the settings for the setucRefScanMode(and RefLinesPE)
    // - select automatically PAT_REF_SCAN_INTRINSIC_AVE if prerequisites are met,
    //   otherwise select default (most probably PAT_REF_SCAN_INPLACE)
    // - careful if PAT_REF_SCAN_EXTRA (or anything other than INPLACE) is selected!
    // - better check whether iPAT is active at all
    if ( pThis->prot().getsPat().getucPATMode() != SEQ::PAT_MODE_NONE ) {
        
        // Conditions for PAT averaging
        bool bPATaveragingIsPossible = 
               lNewVal*pThis->prot().getsPat().getlAccelFact3D() <= pThis->prot().getlAverages();
        
        // This additionally checks whether PAT averaging is allowed by the sequence limits
        _setPATaveraging( pThis, bPATaveragingIsPossible );

    } // if ( pThis->prot().getsPat().getucPATMode() != SEQ::PAT_MODE_NONE )

    // Run original handler with new settings
    lRet = _PATSetValueAndUpdateRefLines(pThis, lNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigAccelFactorPESetFct() );
    
    // Call reference line update as before
    return lRet;
    
} // _PATAccelFActorPE_SetValue ()

double _ReadFOV_SetValue(MrUILinkLimited<double>* const pThis, double dNewVal, int32_t lIndex)
{
    return _PATSetValueAndUpdateRefLines(pThis, dNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigReadFOVSetFct() );
}

double _PhaseFOV_SetValue(MrUILinkLimited<double>* const pThis, double dNewVal, int32_t lIndex)
{
    return _PATSetValueAndUpdateRefLines(pThis, dNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigPhaseFOVSetFct() );
}

double _PhaseOversampling_SetValue(MrUILinkLimited<double>* const pThis, double dNewVal, int32_t lIndex)
{
    return _PATSetValueAndUpdateRefLines(pThis, dNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigPhaseOversamplingSetFct() );
}

double _PhaseResolution_SetValue(MrUILinkLimited<double>* const pThis, double dNewVal, int32_t lIndex)
{
    return _PATSetValueAndUpdateRefLines(pThis, dNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigPhaseResolutionSetFct() );
}

unsigned _PhasePartialFourier_SetValue(MrUILinkSelection<unsigned>* const pThis, unsigned uNewVal, int32_t lIndex)
{
    return _PATSetValueAndUpdateRefLines(pThis, uNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigPhasePartialFourierSetFct() );
}

unsigned _PATMode_SetValue(MrUILinkSelection<unsigned>* const pThis, unsigned uNewVal, int32_t lIndex)
{
    
    // Save the value of PAT mode in the original protocol
    SEQ::PATSelMode originalPATMode = pThis->prot().getsPat().getucPATMode();
    
    // Try to set the PAT mode as if there were no PAT averaging (i.e., original code)
    // - here, it seems OK to first call the original handlers:
    //   - if we are switching PAT off, then the refScanMode is irrelevant
    //   - if we are switching PAT ON, then the original PAT mode is known to be working. We can switch PATaveraging on from there
    ///    and check whether it succeeds.
    unsigned uRet = _PATSetValueAndUpdateRefLines(pThis, uNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigPATModeSetFct() );
    
    // Check whether original attempt succeeds
    if ( uRet != uNewVal ) 
    {
        // new PAT-mode could not be set - don't even try PATaveraging.
        return uRet;
    }

    // Original set attempt succeeded, now check whether this was OFF -> ON
    // - uRet is not of type SEQ::PATSelMode !
    if (    (originalPATMode               == SEQ::PAT_MODE_NONE)
         && (pThis->prot().getsPat().getucPATMode() != SEQ::PAT_MODE_NONE)   ) {
        
        // PAT has just been switched ON, check whether PAT averaging is possible
        // - see also acceleration set-value handler above

        // Condition for PAT averaging: AFtotal >= averages
        bool bPATaveragingIsPossible = 
               pThis->prot().getsPat().getlAccelFactPE()*pThis->prot().getsPat().getlAccelFact3D() <= pThis->prot().getlAverages();
        
        // This additionally checks whether PAT averaging is allowed by the sequence limits
        _setPATaveraging( pThis, bPATaveragingIsPossible );

        // Run again to do a few checks (not many, we don't have to re-calculate refLines, just check)
        uRet = _PATSetValueAndUpdateRefLines(pThis, uNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_PAT_OrigPATModeSetFct() );

    } // if PAT OFF -> ON
    
    return uRet;
    
} // _PATMode_SetValue ()

/// refScanMode overload for, e.g., PATaveraging
/// - for switching automatically, as a consequence of changes in Averages, PATmode, and AF
/// - \ToDo: overload the getOptions handler to only show the automatically chosen mode later (we probably
///          cannot get rid of it completely...). For now, it would also be nice to be able 
///          to switch PATaveraging off if desired
///   - done 20060406
unsigned _RefScanMode_SetValue(MrUILinkSelection<unsigned>* const pThis, unsigned uNewVal, int32_t lIndex)
{

    // Make a first attempt at setting the new refScanMode
    // - even that may be unnecessary...
    const UIHandlersPAT * pUIHandlersPAT = getPointerToUIHandlersPAT(pThis);
    MrUILinkSelection<unsigned>::PFctSetValue pOrigRefScanModeSetValueFct = pUIHandlersPAT->get_PAT_OrigRefScanModeSetFct();
    // Note: this may return success, but still fail internally, notably if no valid number of reference lines
    //       can be found. We need to check below.
    unsigned uRet = _PATSetValueAndUpdateRefLines(pThis, uNewVal, lIndex, pOrigRefScanModeSetValueFct );

    // Changing the refScanMode will affect the PE line and PatRefLinesPE handling! 
    // - This has different consequences depending on the TurboFactor / Shots mode
    MrProtocolData::MrProtData*  pMrProt   = &pThis->prot();
    SeqLim* _pSeqLim = &pThis->seqLimits();
    MrProt sMrProt(pMrProt);
    if( pMrProt->getsFastImaging().getucSegmentationMode() == SEQ::SEGM_MODE_DEFINE_SHOTS )
    {
        // Re-set the shots value from the protocol. This will re-calculate a consistent
        //   value for the PE lines (based on the new refScanMode, which has already been
        //   set; and which causes the internal value for number of refLines to be 0)
        if( LINK_LONG_TYPE* pShots = _search<LINK_LONG_TYPE>(pThis, MR_TAG_SHOTS) )
        {
            if( pShots->isAvailable(0) )
            {
                pShots->value(pMrProt->getsFastImaging().getlShots(),0);
            }
        }
    }
    else
    {
        // Not in shots mode here: 
        // Check whether just setting the refScanMode was enough
        // - we may have run into problems with the phase oversampling: there are no inplace refLines
        //   to always stay below the maximum hidden oversampling for all values of, say, phase resolution.
        // - \ToDo: this is somewhat ugly, since we duplicate UI functionality here. Find an alternative if
        //          possible! (Maybe using the PE lines set-value handler?)
        // - 20060419: in some cases (notably large TurboFactors) _PATSetValueAndUpdateRefLines cannot
        //             find a valid number of reference lines for the current settings (when switching back
        //             to INPLACE). The previous check for _iPATcheckPhOSAndPEFTLen only did not find it
        //             necessary to "solve" by TurboFactor below, which made it impossible to go back to
        //             1 average in one known case. We really have to check for 
        //             pRI->prepareCalculation(pMrProt,_pSeqLim) and repeat _PATSetValueAndUpdateRefLines
        //             as necessary.
        ReorderInfo * pRI = pUIHandlersPAT->get_pReorder();
        bool bReoInfoOK = pRI->prepareCalculation(sMrProt,*_pSeqLim);
        bool bPhOsOK    = _iPATcheckPhOSAndPEFTLen(pThis);
        bool bFoundOK   = bPhOsOK && bReoInfoOK && (uRet==uNewVal);
        
        if ( !bFoundOK )
        {
        
            // Do not return success if changes are necessary
            uRet = MRI_STD_PAT_REF_SCAN_UNDEFINED;

            // Otherwise: check whether we can "solve" by turbo setflFactor(and possibly echo time)
            // - Note: this is necessary to trigger "solve handling", e.g., when the Averages
            //         are changed and a subsequent TF / TE/TI/TR change is necessary. This code
            //         was duplicate to the one in _RefScanMode_Solve; it's most probably
            //         obsolete there, and has been removed.
            // - set-value handler for TF
            LINK_LONG_TYPE* pTurboFactor = 
                _search<LINK_LONG_TYPE>(pThis, MR_TAG_TURBO_FACTOR, 0, MrUILinkBase::SEARCH_MODE::SEARCH_EDITABLE);
            
            // If turboFactor is accessible, try to "solve"
            if ( pTurboFactor != nullptr )
            {
                
                // Original TurboFactor value
                int32_t lOrigTF      = pTurboFactor->value(0);
                
                // Value we will be working with below
                int32_t lCurrTF      = lOrigTF;
                
                // Flag for remembering whether a new TurboFactor value has been tried
                bool bNewTFValue  = false;
                
                // Sequence limits increment of the turboFactor for playing around
                int32_t lTFIncrement = _pSeqLim->getTurboFactor().getInc();
                int32_t lTFMin       = _pSeqLim->getTurboFactor().getMin();
                
                // Increment may not be zero (or negative)
                if ( lTFIncrement > 0 )
                {
                    
                    while ( !bFoundOK && ((lCurrTF-lTFIncrement)>=lTFMin) )
                    {
                        
                        // Try to decrease turbo factor
                        // - Q: should we increase if switching back to INPLACE fails?
                        lCurrTF -= lTFIncrement;
                        
                        // Set it
                        // - necessary to find a working value, since there is no solve handler
                        // - forced, to trigger other solve handlers
                        // - \ToDo: check whether we have to clear the dependency list if we loop over different
                        //          values here.
                        pTurboFactor->value( lCurrTF, 0, MrUILinkBase::SET_MODE::SET_FORCED );
                        bNewTFValue = true;
                        
                        // Was that it?
                        unsigned uRetTemp = _PATSetValueAndUpdateRefLines(pThis, uNewVal, lIndex, pOrigRefScanModeSetValueFct );
                        bReoInfoOK = pRI->prepareCalculation(sMrProt,*_pSeqLim);
                        bPhOsOK    = _iPATcheckPhOSAndPEFTLen(pThis);
                        bFoundOK   = bPhOsOK && bReoInfoOK && (uRetTemp==uNewVal);
                        
                    } // while
                    
                } // if ( lTFIncrement > 0 )
                
                // If we have found a working protocol, signal that to the outside
                if ( bFoundOK )
                {
                    
                    // If changing the turboFactor was the solution, make it a dependent parameter
                    if ( bNewTFValue )
                    {
                        pThis->addDependentParamPtr(pTurboFactor,0);
                        
                        // Note: there may be other parameters in turn dependent on TurboFactor, i.e., TE/TI/TR and
                        //       PPF. This has to be handled in the TurboFactor solve handler (a_tse_overloadUILink)
                        
                    }
                }
                else 
                {
                    // Make the return value 0
                    // - undefined? May be obsolete here - see above
                    uRet = MRI_STD_PAT_REF_SCAN_UNDEFINED;

                    // If we had to change the turboFactor, but did not find a working protocol, change it back
                    if ( bNewTFValue )
                    {
                        // should not need forced since that is where we came from
                        // - a more elegant solution is probably to make a copy of the original protocol and copy that back...
                        pTurboFactor->value( lOrigTF, 0 );

                    } // if ( bNewTFValue )
                    
                } // if ( bFoundOK )
                
            } // if ( pTurboFactor != nullptr )

        } // if ( !bFoundOK )

    } // if( pMrProt->getsFastImaging().getucSegmentationMode() == SEQ::SEGM_MODE_DEFINE_SHOTS )

    return uRet;
}


/// Averages setValue handler overload for PATaveraging
int32_t _Averages_SetValue(MrUILinkLimited<int32_t>* const pThis, int32_t lNewVal, int32_t lIndex)
{
    // Here, it's very debatable whether to first call the original handlers.
    // - for one critical protocol (very large turboFactor and PATaveraging active) we need the call here
    //   to get the solve handling (dependentParameter) working when reducing the averages
    int32_t lRet = 0;
    
    lRet = _PATSetValueAndUpdateRefLines(pThis, lNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_OrigAveragesSetFct() );

    MrProtocolData::MrProtData*  pMrProt = &pThis->prot();

    // PATaveraging: see _PATAccelFActorPE_SetValue above
    if ( pMrProt->getsPat().getucPATMode() != SEQ::PAT_MODE_NONE ) {

        // Conditions for PAT averaging
        bool bPATaveragingIsPossible = 
               pMrProt->getsPat().getlAccelFactPE()*pMrProt->getsPat().getlAccelFact3D() <= pMrProt->getlAverages();
        
        // This additionally checks whether PAT averaging is allowed by the sequence limits
        _setPATaveraging( pThis, bPATaveragingIsPossible );

        lRet = _PATSetValueAndUpdateRefLines(pThis, lNewVal, lIndex, (getPointerToUIHandlersPAT(pThis))->get_OrigAveragesSetFct() );
        
    } // if ( pThis->prot().getsPat().getucPATMode() != SEQ::PAT_MODE_NONE )

    // Return
    return lRet;
    
} // _Averages_SetValue ()


// ------------------------------------------------------------------------------
// Functions   : _PATRefLinesPE_GetLimits
// ------------------------------------------------------------------------------
//               
// Description : - call orig. handler
//               - for segemented protocols: set VERIFY_SCAN_ALL for ref.lines PE
//
// Return      : true, if success
//
// ------------------------------------------------------------------------------
bool _PATRefLinesPE_GetLimits (LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t lIndex)
{

    /// \ToDo: for PATaveraging, we either will clamp the RefLinesPE to a certain value (e.g., 32)
    ///        or let the user choose it freely ( = original getLimits handler ? )
    /// - done 20060406: clamp at optimum (or less if less image lines are available)

    const UIHandlersPAT * pUIHandlersPAT = getPointerToUIHandlersPAT(pThis);

    bool bRet = (*pUIHandlersPAT->get_PAT_OrigRefLinesPEGetLimitsHandler())(pThis, rLimitVector, rulVerify, lIndex);

    MrProtocolData::MrProtData*  pMrProt = &(pThis->prot());
    MrProt sMrProt(pMrProt);

    // do we need VERIFY_SCAN_ALL ?

    // in shots mode, we don't need to scan all ref.lines any more (much faster)
    if( pMrProt->getsFastImaging().getucSegmentationMode() == SEQ::SEGM_MODE_DEFINE_SHOTS )
    {
        return bRet;
    }

    // Also: for PAT averaging the number of reference lines is in principle freely selectable
    // -> no scan needed
    // - ultimately, we would like to clamp the number of refLines to a certain value, e.g., 32
    //   (for this function it would probably be best to use the current value from the protocol)
    //   - done 20060406: but this number of refLines may not be the one we want (from old
    //     protocols); either fix that here, or in updateRefLines
    if ( pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_INTRINSIC_AVE )
    {
        
        // PATaveraging: there is an optimal number of reference lines which is set automatically
        //   when this mode becomes active (see setValue handler for PAT mode)
        //   -> only option is the current value
        rLimitVector.resize(1);
        if ( rLimitVector.size() == 1 )
        {
            rLimitVector[0].setLonely(pMrProt->getsPat().getlRefLinesPE());
            return true;
        }
        else
        {
            return false;
        }
        
        //return bRet;
    }

    // VERIFY_SCAN_ALL required for segmented scanning schemes
    // (only certain ref.line numbers are allowed, because the number of ref.lines has to be calculated
    //  in a way, that the scanned number of lines (i.e. PATLines + RefLines) is a multiple of the segments;
    //  this is done in _PATSetValueAndUpdateRefLines())
    // Note: in order to get a valid number of segments for the current protocol, we need to call Reorder.prepareCalculation() first

    /// \ToDo: check whether compound boolean expressions like the one below can be assumed to have left-to-right
    ///        precedence (otherwise we might be dereferencing a nullptr pointer in certain circumstances)

    (void) (static_cast<SeqIF_UIHandlersPAT*>(pThis->sequence().getSeq()))->
        SeqIF_UIHandlersPAT::PATInitializeReorderMethod( pMrProt, &(pThis->seqLimits()) );

    if ( pUIHandlersPAT->get_pReorder() && 
         pUIHandlersPAT->get_pReorder()->prepareCalculation(sMrProt, (pThis->seqLimits())) &&
         pUIHandlersPAT->get_pReorder()->getMeasuredSegments() > 1 )
    {
            rulVerify = LINK_LONG_TYPE::VERIFY_SCAN_ALL;
    }

    return bRet;
}






    ///
    ///   Some background explanation on using PAT with segmented sequences
    ///
    ///
    ///   PELines are calculated within UILink/MrProtocolData::MrProtData to be a multiple of the segementation factor -
    ///   this is required for any segmented sequence (although of course the sequence is free in how many lines
    ///   it acquires).
    ///
    ///   However, for PAT it is not the actual PELines value that is relevant for the segmentation,
    ///   but the PAT-reduced number of acquired PELines (including inplace ref.lines).
    ///   E.g. for a matrix of 256 PELines, a PAT acceleration of x2 would lead to 128 acquired PELines
    ///   plus e.g. 12 additional inplace ref.lines: 140 (PAT-)PELines have to be a multiple of the segmentation factor!
    ///
    ///   Originally, UILink/MrProtocolData::MrKSpaceData didn't know about this and behaved as if no iPAT was selected.
    ///   The sequence resolved this issue by increasing the inplace ref.lines until also the PAT-acquired
    ///   (reduced) PELines were a multiple of the segmentation.
    ///   Disadvantages of this approach:
    ///   - PELines are increased (using hidden oversampling) by UILink without any need (only PAT-PELines need to be
    ///     on the segmentation grid!)
    ///   - a lot of UILink-code within any segmented sequence to adapt the inplace ref.lines
    ///
    ///
    ///   To improve this situation, new PELines increment modes SEQ::INC_..._IPAT modes were introduced (set by _pSeqLim->setPELines(...)).
    ///   If one of these new modes is set, the increment for PELines is just set to the PE-acceleration factor.
    ///   This requires, that the calculation procedure for PELines takes care that the *reduced* PAT-PELines 
    ///   (including inplace ref.lines) are a multiple of the segmentation setflFactor(this is done in MrProtocolData::MrKSpaceData::correctPELinesForPat()).
    ///   This enables a segmented sequence to get a valid number of PELines from 
    ///   kSpace().linesToMeasure(lLinesToMeasure).
    ///   Disadvantages:
    ///   - The conventions of the PAT k-space sampling scheme (including how inplace ref.lines are integrated)
    ///     had to be copied from the sequence (MrImaging\libSeqUtil\ReorderInfo) to MrProtocolData::MrProtData\MrKspace.
    ///     Care has to be taken, that they are kept equivalent!
    ///   - Partial fourier is NOT considered when calculating PELines (because UILink in general does not know the
    ///     real partial fourier factor, this is handled within the sequence).
    ///





// NOTE: "...SegmentedSequences..." may not be the correct name anymore if we are using
//    PATaveraging for other sequences, too, and have to register a subset of these UI
//    handlers.


bool fPATRegisterUILinkHandlersForSegmentedSequences( SeqLim &rSeqLim, 
                                                      ReorderInfo *pRI, 
                                                      SEQ::Increment uPELineIncMode, 
                                                      int32_t lInpPATRefLinesOpt, 
                                                      UIHandlersPAT *pUIHandlersPAT )
{
    return fPATRegisterUILinkHandlersForSegmentedSequences( &rSeqLim, pRI, uPELineIncMode, lInpPATRefLinesOpt, pUIHandlersPAT );
}



// ------------------------------------------------------------------------------
// Function    : fPATRegisterUILinkHandlersForSegmentedSequences
// ------------------------------------------------------------------------------
//               
// Description : Registers SetValue-handlers for all parameters that influence
//               the number of phase encoding lines. Those set value handlers
//               are typically needed for segmented sequences.
//
//               NOTE: should be called late in fSEQInit, i.e. after all UILink
//                     registration steps have been performed.
//
//               Pointers to orig. setValueHandlers and additional information (e.g. pointer to
//               ReorderInfo class, default number of Ref.Lines, etc) are stored in 
//               class UIHandlersPAT.
//
//               IMPORTANT:
//               'class-type' sequences may specify a pointer to an own object of UIHandlersPAT,
//               otherwise (i.e. pointer is nullptr) the global instance 'myUIHandlersPAT' will be used instead
//
// Return      : true , if for success
//               false, else
//
// ------------------------------------------------------------------------------
bool fPATRegisterUILinkHandlersForSegmentedSequences(SeqLim *_pSeqLim, ReorderInfo *pRI, SEQ::Increment uPELineIncMode, int32_t lInpPATRefLinesOpt, UIHandlersPAT *pUIHandlersPAT )
{

    if ( pUIHandlersPAT == nullptr )
    {
        // since calling sequence didn't provide a pointer to a UIHandlersPAT object, we use the local one
        // (note: nothing wrong with that for 'conventional' sequences. However, for 'class-type' sequences this is required,
        //        because each sequence instance needs its own set of UIHandlers to avoid interference of multiple sequence instances!
        //        Refer to comments in iPAT.h)
        pUIHandlersPAT = &myUIHandlersPAT;
    }

    pUIHandlersPAT->set_lUILinkPATRefLinesOpt ( (lInpPATRefLinesOpt < 0 ? lPATDefOptNoRefLines : lInpPATRefLinesOpt) );

    if (!pRI)
    {
        SEQ_TRACE_ERROR.print("pRI is nullptr");
        return false;
    }

    pUIHandlersPAT->set_pReorder ( pRI );

    // Remark: TAGs defined in \N4\PKG\MRSERVERS\MRPROTSRV\MRPROTOCOL\UILink\MrStdNameTags.h
    MrUILinkLimited<int32_t>       *_TurboFactor        = _search< MrUILinkLimited<int32_t>       >(_pSeqLim, MR_TAG_TURBO_FACTOR         );
    MrUILinkLimited<int32_t>       *_Segments           = _search< MrUILinkLimited<int32_t>       >(_pSeqLim, MR_TAG_SEGMENTS             );
    MrUILinkLimited<int32_t>       *_BaseResolution     = _search< MrUILinkLimited<int32_t>       >(_pSeqLim, MR_TAG_BASE_RESOLUTION      );
    MrUILinkLimited<int32_t>       *_PATAccelFActorPE   = _search< MrUILinkLimited<int32_t>       >(_pSeqLim, MR_TAG_PAT_ACC_PE           );
    MrUILinkLimited<int32_t>       *_PATRefLinesPE      = _search< MrUILinkLimited<int32_t>       >(_pSeqLim, MR_TAG_PAT_LINES_PE         );
    MrUILinkLimited<double>     *_ReadFOV            = _search< MrUILinkLimited<double>     >(_pSeqLim, MR_TAG_READOUT_FOV          );
    MrUILinkLimited<double>     *_PhaseFOV           = _search< MrUILinkLimited<double>     >(_pSeqLim, MR_TAG_PHASE_FOV            );
    MrUILinkLimited<double>     *_PhaseOversampling  = _search< MrUILinkLimited<double>     >(_pSeqLim, MR_TAG_PHASE_OVERSAMPLING   );
    MrUILinkLimited<double>     *_PhaseResolution    = _search< MrUILinkLimited<double>     >(_pSeqLim, MR_TAG_PHASE_RESOLUTION     );
    MrUILinkSelection<uint32_t> *_PhasePartialFourier = _search< MrUILinkSelection<uint32_t> >(_pSeqLim, MR_TAG_PHASE_PARTIAL_FOURIER);
    MrUILinkSelection<uint32_t> *_PATMode = _search< MrUILinkSelection<uint32_t> >(_pSeqLim, MR_TAG_PAT_MODE);
    // Specifically for PATaveraging
    MrUILinkSelection<uint32_t> *_RefScanMode = _search< MrUILinkSelection<uint32_t> >(_pSeqLim, MR_TAG_PAT_REF_SCAN_MODE);
    MrUILinkLimited<int32_t>       *_Averages           = _search< MrUILinkLimited<int32_t>       >(_pSeqLim, MR_TAG_AVERAGES             );

    if (uPELineIncMode==SEQ::INC_TURBO_FACTOR)  {
        if (_TurboFactor)     {pUIHandlersPAT->set_PAT_OrigTurboFactorSetFct          ( _TurboFactor        ->registerSetValueHandler (_TurboFactor_SetValue        ) );}
    }
    if (uPELineIncMode==SEQ::INC_GRE_SEGMENTS)  {
        if (_Segments)        {pUIHandlersPAT->set_PAT_OrigSegmentsSetFct             ( _Segments           ->registerSetValueHandler(_Segments_SetValue            ) );}
    }
    if (_BaseResolution)      {pUIHandlersPAT->set_PAT_OrigBaseResolutionSetFct       ( _BaseResolution     ->registerSetValueHandler (_BaseResolution_SetValue     ) );}
    if (_PATAccelFActorPE)    {pUIHandlersPAT->set_PAT_OrigAccelFactorPESetFct        ( _PATAccelFActorPE   ->registerSetValueHandler (_PATAccelFActorPE_SetValue   ) );}
    if (_ReadFOV)             {pUIHandlersPAT->set_PAT_OrigReadFOVSetFct              ( _ReadFOV            ->registerSetValueHandler (_ReadFOV_SetValue            ) );}
    if (_PhaseFOV)            {pUIHandlersPAT->set_PAT_OrigPhaseFOVSetFct             ( _PhaseFOV           ->registerSetValueHandler (_PhaseFOV_SetValue           ) );}
    if (_PhaseOversampling)   {pUIHandlersPAT->set_PAT_OrigPhaseOversamplingSetFct    ( _PhaseOversampling  ->registerSetValueHandler (_PhaseOversampling_SetValue  ) );}
    if (_PhaseResolution)     {pUIHandlersPAT->set_PAT_OrigPhaseResolutionSetFct      ( _PhaseResolution    ->registerSetValueHandler (_PhaseResolution_SetValue    ) );}
    if (_PhasePartialFourier) {pUIHandlersPAT->set_PAT_OrigPhasePartialFourierSetFct  ( _PhasePartialFourier->registerSetValueHandler (_PhasePartialFourier_SetValue) );}
    if (_PATMode)             {pUIHandlersPAT->set_PAT_OrigPATModeSetFct              ( _PATMode            ->registerSetValueHandler (_PATMode_SetValue            ) );}
    // Specifically for PATaveraging
    if (_RefScanMode)         {pUIHandlersPAT->set_PAT_OrigRefScanModeSetFct          ( _RefScanMode        ->registerSetValueHandler (_RefScanMode_SetValue        ) );}
    if (_Averages)            {pUIHandlersPAT->set_OrigAveragesSetFct                 ( _Averages           ->registerSetValueHandler (_Averages_SetValue           ) );}

    if (_PATRefLinesPE)       {pUIHandlersPAT->set_PAT_OrigRefLinesPEGetLimitsHandler ( _PATRefLinesPE      ->registerGetLimitsHandler(_PATRefLinesPE_GetLimits     ) );}

    return true;
}

#endif






/// ------------------------------------------------------------------------------
/// Functions   : fPATSetDefaultSeqLim
/// ------------------------------------------------------------------------------
///               
/// Description : set default SeqLim for PAT: PATMode, RefScanMode, AccelFactorPE, RefLinesPE
///               Note: 3D PAT parameters (AccelFactor3D, RefLines3D) are only set, 
///                     if the corresponding flag (second argument in interface) is 'true'
///                     (if second argument is left out, default is 'false', since currently
///                      most sequences do not support PAT^2)
///               
///               The maximum accelerations factors may be imported via registry key:
///                  SOFTWARE/Siemens/Numaris4/Config/Modality/Sequence/PAT/ 
///                     -> PATMaxAccelPE
///                     -> PATMaxAccel3D
///               If this key is not available, default values are set.
///
///
/// Return      : MRI_SEQ_SEQU_NORMAL for success
///
///               MRI_SEQ_SEQU_ERROR,
///                    else
///
///
/// ------------------------------------------------------------------------------
NLS_STATUS fPATSetDefaultSeqLim( SeqLim *_pSeqLim, bool bInitPAT3D )
{  
    // default values for max. accel. factors
    // (if not modified by registry key)
    const int32_t lMaxAccelPE_default = 8;         // default accel.factor in PE
    const int32_t lMaxAccel3D_default = 8;         // default accel.factor in 3D
    int32_t lMaxAccelPE = lMaxAccelPE_default;
    int32_t lMaxAccel3D = lMaxAccel3D_default;

    _pSeqLim->setPATMode       (SEQ::PAT_MODE_NONE,SEQ::PAT_MODE_SENSE,SEQ::PAT_MODE_GRAPPA);
    _pSeqLim->setRefScanMode   (SEQ::PAT_REF_SCAN_INPLACE                      );
    _pSeqLim->setAccelFactorPE (       1, lMaxAccelPE,           1,           2);
    //  Changed Minimum number of reference lines from 2 to 8,
    //  since fPATPrepPost enforces at least 4 times acceleration factor
    // - use lPATDefOptNoRefLines as default?
    _pSeqLim->setRefLinesPE    (       8,         128,           1,          24);

    if ( bInitPAT3D )
    {
        _pSeqLim->setAccelFactor3D (       1, lMaxAccel3D,           1,           1);
        // - use lPATDefOptNoRefLines as default?
        _pSeqLim->setRefLines3D    (       2,         128,           1,          24);
    }
    else
    {
        _pSeqLim->setAccelFactor3D (       1,           1,           1,           1);
        _pSeqLim->setRefLines3D    (       1,           1,           1,           1);
    }

    return MRI_SEQ_SEQU_NORMAL;
}

    
    


    
/// ------------------------------------------------------------------------------
/// Functions   : fPATPrepPost
/// ------------------------------------------------------------------------------
///               
/// Description : - Takes care for some iPAT related restrictions
///               - Sets pSeqExpo->setFirstFourierLine/Partition and ->setFirstRefLine/Partition 
///                 (without iPAT, these values keep their default values 0 and -1)
///
/// Return      : MRI_SEQ_SEQU_NORMAL for success
///
///               MRI_SEQ_SEQU_PAT_NO_OF_RX_CHANNELS_TOO_LOW,
///                    if number of selected coil elements is lower than acceleration factor
///
///               MRI_SEQ_SEQU_ERROR,
///                    else
///
///
/// ------------------------------------------------------------------------------
NLS_STATUS fPATPrepPost(MrProtocolData::MrProtData* pMrProt, SeqLim *_pSeqLim, SeqExpo *pSeqExpo, ReorderInfo *pReorderInfo)
{
    if (pReorderInfo == nullptr)
    {
        SEQ_TRACE_ERROR.print("fPATPrepPost() failed: function call without pointer to ReorderInfo class");
        return MRI_SEQ_SEQU_ERROR;
    }

    // check accel. PE/3D: both must not be zero in order to avoid modulo or division by zero
    if (pMrProt->getsPat().getlAccelFactPE() == 0 || pMrProt->getsPat().getlAccelFact3D() == 0)
    {
        SEQ_TRACE_ERROR.print("fPATPrepPost() failed: Accel.PE or Accel3D == 0, this must never happen!");
        return MRI_SEQ_SEQU_ERROR;
    }

    MrProt sMrProtWrapper(pMrProt);
    MrKSpace & rKSpaceWrapper = sMrProtWrapper.kSpace();


    // Peschomir: Do not execute all of this for CSENSE
    // if (pMrProt->getsPat().getucPATMode()!=SEQ::PAT_MODE_NONE)
    if (pMrProt->getsPat().getucPATMode() != SEQ::PAT_MODE_NONE && 
        pMrProt->getsPat().getucPATMode() != SEQ::PAT_MODE_CSENSE_CINE &&
        pMrProt->getsPat().getucPATMode() != SEQ::PAT_MODE_CSENSE_TOF &&
        pMrProt->getsPat().getucPATMode() != SEQ::PAT_MODE_CSENSE_SPACE)
    {
#ifdef WIN32   // registry entries not available on MCIR and thus seqLim has not been modified!
        //
        // Consistency Check:   The acceleration factor in the protocol must not be larger than the maximum
        //                      from SeqLim. This case might occur, if the maximum acceleration factor was
        //                      locally increased via registry key and the protocol is imported to
        //                      another scanner (or the registry value has been reduced).
        if ((pMrProt->getsPat().getlAccelFactPE() > _pSeqLim->getAccelFactorPE().getMax()) ||
            (pMrProt->getsPat().getlAccelFact3D() > _pSeqLim->getAccelFactor3D().getMax()))
        {
            SEQ_TRACE_ALWAYS.print("PAT accel. factor exceeds SequenceLimits; obviously this protocol was set up with a locally changed (via registry) max. accel. factor.");
            return MRI_SEQ_SEQU_ERROR;
        }
#endif

        //
        // Check some iPAT related restrictions:
        //
        // SENSE is not possible with 'save uncombined images'
        if ((pMrProt->getsPat().getucPATMode() == SEQ::PAT_MODE_SENSE) && (pMrProt->getucUncombImages() == true))
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())  {
                SEQ_TRACE_ERROR.print("fPATPrepPost() failed: iPAT mode 'SENSE' not allowed with 'save uncombined images'");
            }
            return MRI_SEQ_SEQU_ERROR;
        }
        // SENSE is not possible with other reconModes than 'Magnitude'
        if ((pMrProt->getsPat().getucPATMode() == SEQ::PAT_MODE_SENSE) &&
            (pMrProt->getucReconstructionMode() != SEQ::RECONMODE_MAGNITUDE))
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())  {
                SEQ_TRACE_ERROR.print("fPATPrepPost() failed: iPAT mode 'SENSE' only allowed with recon.mode magnitude");
            }
            return MRI_SEQ_SEQU_ERROR;
        }
        // SENSE is not possible with 'save unfiltered images'
        if ((pMrProt->getsPat().getucPATMode() == SEQ::PAT_MODE_SENSE) && (pMrProt->getsNormalizeFilter().getucSaveUnfiltered() == true))
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())  {
                SEQ_TRACE_ERROR.print("fPATPrepPost() failed: iPAT mode 'SENSE' not allowed with 'save unfiltered images'");
            }
            return MRI_SEQ_SEQU_ERROR;
        }
        // SENSE requires noise adjust
        if ((pMrProt->getsPat().getucPATMode() == SEQ::PAT_MODE_SENSE) && (pMrProt->getucEnableNoiseAdjust() == false))
        {
            if (_pSeqLim->isContextPrepForMrProtUpdate())
            {
                pMrProt->setucEnableNoiseAdjust(true);
            }
            else
            {
                if (!_pSeqLim->isContextPrepForBinarySearch())  {
                    SEQ_TRACE_ERROR.print("fPATPrepPost() failed: iPAT mode 'SENSE' requires noise adjust");
                }
                return MRI_SEQ_SEQU_ERROR;
            }
        }
        //
        if ((pMrProt->getsPat().getucPATMode() == SEQ::PAT_MODE_SENSE) && (pMrProt->getucCoilCombineMode() == SEQ::COILCOMBINE_ADAPTIVE_COMBINE))
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())  {
                SEQ_TRACE_ERROR.print("fPATPrepPost() failed: iPAT mode 'SENSE' not allowed with CoilCombineMode 'ADAPTIVE_COMBINE'");
            }
            return MRI_SEQ_SEQU_ERROR;
        }



        // see charm 00308583 
        // for GRAPPA inplace: ReferenceLines >= #AccelerationFactor * 4 
        // for GRAPPA extra: ReferenceLines >= #AccelerationFactor * 6 
        if (pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA || pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST || pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_EPI || pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_TSE)
        {
            if (pMrProt->getsPat().getlRefLinesPE() < pMrProt->getsPat().getlAccelFactPE() * 6)
            {
                if (!_pSeqLim->isContextPrepForBinarySearch())  {
                    SEQ_TRACE_ERROR.print("fPATPrepPost failed. RefLines < Accel * 6");
                }
                return MRI_SEQ_SEQU_ERROR;
            }
        }
        else
        {
            if (pMrProt->getsKSpace().getucTrajectory() == SEQ::TRAJECTORY_BLADE)
            {
                if (pMrProt->getsPat().getlRefLinesPE() < pMrProt->getsPat().getlAccelFactPE() * 2)
                {
                    if (!_pSeqLim->isContextPrepForBinarySearch())  {
                        SEQ_TRACE_ERROR.print("fPATPrepPost failed. RefLines < Accel * 2");
                    }
                    return MRI_SEQ_SEQU_ERROR;
                }
            }
            else
            {
                if (pMrProt->getsPat().getlRefLinesPE() < pMrProt->getsPat().getlAccelFactPE() * 4)
                {
                    if (!_pSeqLim->isContextPrepForBinarySearch())
                    {
                        SEQ_TRACE_ALWAYS.print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
                        SEQ_TRACE_ALWAYS.print("Ref.Lines PE = %i", pMrProt->getsPat().getlRefLinesPE());
                        SEQ_TRACE_ALWAYS.print("Accel PE     = %i", pMrProt->getsPat().getlAccelFactPE());
                        SEQ_TRACE_ERROR.print("fPATPrepPost failed. RefLines < Accel * 4");
                    }
                    return MRI_SEQ_SEQU_ERROR;
                } 
            }
        }

        if (pMrProt->getsKSpace().getucDimension() == SEQ::DIM_3 && pMrProt->getsPat().getlAccelFact3D() > 1)
        {
            if (pMrProt->getsPat().getlRefLines3D() < pMrProt->getsPat().getlAccelFact3D() * 4)
            {
                if (!_pSeqLim->isContextPrepForBinarySearch())  {
                    SEQ_TRACE_ERROR.print("fPATPrepPost failed. RefLines3D < Accel3D * 4");
                }
                return MRI_SEQ_SEQU_ERROR;
            }
        }

        
        if (pMrProt->getsPat().getucPATMode() != SEQ::PAT_MODE_CSENSE_SEMAC)
        {
            // get multi band factor
            long lMultibandfactor = pMrProt->getsSliceAcceleration().getlMultiBandFactor();
            if (lMultibandfactor < 1)
                lMultibandfactor = 1;
    
            bool bIsSliceAcceleration = false;
    
            if (lMultibandfactor > 1 && pMrProt->getsPat().getucPATMode() == SEQ::PAT_MODE_SLICE_ACCELERATION)
                bIsSliceAcceleration = true;
    
            // 2D-PAT: check total PAT > 1
            long lTotalAccelerationFactor = 0;
            if (bIsSliceAcceleration)
            {
                lTotalAccelerationFactor = pMrProt->getsPat().getlAccelFactPE() * lMultibandfactor;
            }
            else
            {
                lTotalAccelerationFactor = pMrProt->getsPat().getlAccelFactPE() * pMrProt->getsPat().getlAccelFact3D();
            }
    
            if (lTotalAccelerationFactor < 2)
            {
                if (!_pSeqLim->isContextPrepForBinarySearch())  {
                    SEQ_TRACE_ERROR.print("fPATPrepPost failed. The total acceleration factor (PE and accel. factor3D / slice) is < 2.");
                }
                return MRI_SEQ_SEQU_ERROR;
            }
        }


        //
        // 2D-PAT restrictions:  
        //

        //  - number of reconstructed partitions (PaFT length) must be multiple of accel.factor3D
        //    (accel.factor3D is always 1 for 1D-PAT) 
        //  - it must be even-numbered always for FFT
        if ((rKSpaceWrapper.partitionFTLen() % 2) || (rKSpaceWrapper.partitionFTLen() % pMrProt->getsPat().getlAccelFact3D()))
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())  {
                SEQ_TRACE_ERROR.print("fPATPrepPost failed. Number of partitions must be multiple of accel.factor3D");
            }
            return MRI_SEQ_SEQU_ERROR;
        }

        //  - slice resolution > 100% not allowed
        if (pMrProt->getsPat().getlAccelFact3D() > 1 && pMrProt->getsKSpace().getdSliceResolution() > 1)
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())  {
                SEQ_TRACE_ERROR.print("fPATPrepPost failed. Slice resolution > 100%% not allowed in combination with iPAT accel.3D");
            }
            return MRI_SEQ_SEQU_ERROR;
        }

        //  - margosian not allowed
        if (pMrProt->getsPat().getlAccelFact3D() > 1 && pSeqExpo->getPCAlgorithm() == SEQ::PC_ALGORITHM_MARGOSIAN)
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())  {
                SEQ_TRACE_ERROR.print("fPATPrepPost failed. Margosian phase correction not allowed in combination with iPAT accel.3D");
            }
            return MRI_SEQ_SEQU_ERROR;
        }

        if (pMrProt->getsPat().getucPATMode() != SEQ::PAT_MODE_CSENSE_SEMAC)
        {
            //  - number of reference partitions must be multiple of accel.factor3D
            if (pMrProt->getsPat().getlRefLines3D() % pMrProt->getsPat().getlAccelFact3D())
            {
                if (!_pSeqLim->isContextPrepForBinarySearch())  {
                    SEQ_TRACE_ERROR.print("fPATPrepPost failed. Number of ref.lines3D must be multiple of accel.factor3D");
                }
                return MRI_SEQ_SEQU_ERROR;
            }
        }

        //  PEFT-Length * Accel3D must not be larger than 2048 (CHARM 323328)
        if ((pMrProt->getsPat().getlAccelFact3D() > 1) && (pSeqExpo->getMeasuredPELines()*pMrProt->getsPat().getlAccelFact3D() > 2048))
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())  {
                SEQ_TRACE_ERROR.print("fPATPrepPost failed. pSeqExpo->getMeasuredPELines() * AccelFact3D > 2048");
            }
            return MRI_SEQ_SEQU_ERROR;
        }


        // check if Accel factor3D  <  3Dshift+1
        //-----------------------------------------------------------------------------------------------------------------------------
        if (pMrProt->getsPat().getlAccelFact3D() - 1 < pMrProt->getsPat().getlReorderingShift3D())
        {
            if (!_pSeqLim->isContextPrepForBinarySearch()) {
                SEQ_TRACE_ERROR.print("check PAT3D and 3D shift  failed. ----- PAT3D-1 > 3Dshift ");
            }
            return MRI_SEQ_SEQU_ERROR;
        }

        // PATaveraging: if REF_SCAN_MODE_INTRINSIC_AVE is active, the number of averages must be greater
        //   or equal to the acceleration factor
        // - may be obsolete for automatic PATaveraging mode
        bool bPATaveragingIsPossible = pMrProt->getsPat().getlAccelFactPE() <= pMrProt->getlAverages()
            && !pMrProt->getucCompT2Decay()  // TSE-T2Comp. not supported with pat-averaging
            && !(pMrProt->getsPat().getucPATMode() == SEQ::PAT_MODE_SLICE_ACCELERATION
                && (!pMrProt->getsSliceAcceleration().getucAverage())
                && pMrProt->getsSliceAcceleration().getlMultiBandFactor() > 1
                ); // PATave not supported with conventional SMS

        if ((pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_INTRINSIC_AVE) && !bPATaveragingIsPossible)
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("fPATPrepPost failed. In mode PAT_REF_SCAN_INTRINSIC_AVE the number of averages must not be smaller than the total acceleration factor");
            }
            return MRI_SEQ_SEQU_ERROR;
        }

#ifndef HASTE
        // PATaveraging: if possible and not selected, this is now an error (except for the TSE/Separate reference scan)
        bool bPatAveragingIsAllowedBySequence = _pSeqLim->getRefScanMode().hasOption(SEQ::PAT_REF_SCAN_INTRINSIC_AVE);
        if (bPATaveragingIsPossible && bPatAveragingIsAllowedBySequence
            && !(pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_INTRINSIC_AVE) && !(pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_TSE || pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA))
        {
            if (!_pSeqLim->isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("fPATPrepPost failed. If the number of averages is greater or equal than the total acceleration factor, the refScanMode must be PAT_REF_SCAN_INTRINSIC_AVE if allowed by the sequence");
            }
            return MRI_SEQ_SEQU_ERROR;
        }
#endif


        /// set YAPS parameters
        ///
        /// Background: The ICE program needs to be informed (via YAPS) about the k-space positions of the first actually measured line and partition
        /// (without iPAT, this value is assumed to be zero. However, with iPAT, the first line(s)/partition(s) might be 'gaps'!).
        /// Since the reference lines are stored in a separate object, the k-space positions of the first reference line/partition are required, too.
        ///

        /// Acq. scheme and YAPS/Mdh parameters for PAT
        /// (8 ref.lines, inplace reference scan mode)
        /// -------------------------------------------
        ///
        ///      x  -  PATScan
        ///      X  -  RefScan (=RefAndImaScan for inplace mode)
        ///      .  -  Gap
        ///
        ///
        ///
        ///
        ///
        ///       . x . x . X X X X X X X X x . x . x
        ///      +-------------------------------------> PE (AccelPE=2)
        ///       0 1 2 3...        ^CenterLin
        ///                 ^
        ///                 firstRefLinNo
        ///         ^
        ///         firstLinNo




        /// Acq. schemes and YAPS/Mdh parameters for PAT^2
        /// (8 ref.lines/8 ref.part, inplace reference scan mode)
        /// ------------------------------------------------------
        ///
        ///      x  -  PATScan
        ///      X  -  RefScan (=RefAndImaScan for inplace mode)
        ///      .  -  Gap
        ///
        ///
        ///      3D (Accel3D=2)
        ///      ^
        ///      |x . x . x . x . x . x . x . x . x . x
        ///      |. . . . . . . . . . . . . . . . . . .
        ///      |x . x . x . x . x . x . x . x . x . x
        ///      |. . . . . . . . . . . . . . . . . . .
        ///      |x . x . x . x . x . x . x . x . x . x
        ///      |. . . . . . X X X X X X X X . . . . .
        ///      |x . x . x . X X X X X X X X x . x . x
        ///      |. . . . . . X X X X X X X X . . . . .
        ///      |x . x . x . X X X X X X X X x . x . x <CenterLin
        ///      |. . . . . . X X X X X X X X . . . . .
        ///      |x . x . x . X X X X X X X X x . x . x
        ///    . |. . . . . . X X X X X X X X . . . . .
        ///    . |x . x . x . X X X X X X X X x . x . x       <- firstRefParNo
        ///    . |. . . . . . . . . . . . . . . . . . .
        ///    3 |x . x . x . x . x . x . x . x . x . x
        ///    2 |. . . . . . . . . . . . . . . . . . .
        ///    1 |x . x . x . x . x . x . x . x . x . x       <- firstParNo
        ///    0 |. . . . . . . . . . . . . . . . . . .
        ///      +-------------------------------------> PE (AccelPE=2)
        ///       0 1 2 3...        ^CenterLin
        ///                 ^
        ///                 firstRefLinNo
        ///       ^
        ///       firstLinNo
        ///
        ///
        ///
        ///
        ///
        ///      3D (Accel3D=2)
        ///      ^
        ///      |x x x x x x x x x x x x x x x x x x
        ///      |. . . . . . . . . . . . . . . . . .
        ///      |x x x x x x x x x x x x x x x x x x
        ///      |. . . . . . . . . . . . . . . . . .
        ///      |x x x x x x x x x x x x x x x x x x
        ///      |. . . . . X X X X X X X X . . . . .
        ///      |x x x x x X X X X X X X X x x x x x
        ///      |. . . . . X X X X X X X X . . . . .
        ///      |x x x x x X X X X X X X X x x x x x <CenterPar
        ///      |. . . . . X X X X X X X X . . . . .
        ///      |x x x x x X X X X X X X X x x x x x
        ///    . |. . . . . X X X X X X X X . . . . .      
        ///    . |x x x x x X X X X X X X X x x x x x      <- firstRefParNo
        ///    . |. . . . . . . . . . . . . . . . . .
        ///    3 |x x x x x x x x x x x x x x x x x x
        ///    2 |. . . . . . . . . . . . . . . . . .
        ///    1 |x x x x x x x x x x x x x x x x x x      <- firstParNo
        ///    0 |. . . . . . . . . . . . . . . . . .
        ///      +-------------------------------------> PE (no acceleration)
        ///       0 1 2 3...        ^CenterLin
        ///                 ^
        ///                 firstRefLinNo
        ///       ^
        ///       firstLinNo=0 (=default, if no accel. in PE)


        // set firstFourierLine/Partition
        if (pMrProt->getsPat().getlAccelFactPE() > 1)
        {
            pSeqExpo->setFirstFourierLine((int32_t)pReorderInfo->getMinPATLinNo());
        }
        else    // no accel in PE direction -> set default values
        {
            pSeqExpo->setFirstFourierLine(0);
        }

        if (pMrProt->getsPat().getlAccelFact3D() > 1)
        {
            pSeqExpo->setFirstFourierPartition((int32_t)pReorderInfo->getMinPATParNo());
        }
        else    // no accel in partition direction -> set default values
        {
            pSeqExpo->setFirstFourierPartition(0);
        }


        // set firstRefLine/Partition
        if (pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_INPLACE)
        {
            // --- INPLACE ---
            // Note: For AccelPE=1 (i.e. accel. only in partition direction), 
            //       'FirstRefLine' is still required to define the range of reference partitions in line direction.
            //       Thus we set it now (without regard of AccelFactPE)
            pSeqExpo->setFirstRefLine((int32_t)pReorderInfo->getMinPATRefLinNo());

            if (pMrProt->getsPat().getlAccelFact3D() > 1)
            {
                pSeqExpo->setFirstRefPartition((int32_t)pReorderInfo->getMinPATRefParNo());
            }
            else    // no accel in partition direction -> set default values
            {
                pSeqExpo->setFirstRefPartition(-1);
            }
        }
        else
        {
            // for SEQ::PAT_REF_SCAN_EXTRA:
            // If the ref.lines are scanned separately, the ICE program
            // assumes that they are scanned symmetric around their kSpaceCenterLine/partition
            // thus FirstRefLine/Partition is not required
            pSeqExpo->setFirstRefLine(-1);
            pSeqExpo->setFirstRefPartition(-1);
        }



        // elliptical scanning:
        // Some ref.lines may be skipped by ellipt.scanning
        int32_t lICEPara = (int32_t)pSeqExpo->getICEProgramParam(ICE_PROGRAM_PARA_MULTI_OFFLINE_CALL);
        if (pMrProt->getsPat().getucRefScanMode() == SEQ::PAT_REF_SCAN_INPLACE  &&  pMrProt->getsKSpace().getucEnableEllipticalScanning())
        {
            lICEPara |= ICE_PROGRAM_MSK_IPAT_TRIGGER_LASTSCANINMEAS;
        }
        else
        {
            lICEPara &= ~ICE_PROGRAM_MSK_IPAT_TRIGGER_LASTSCANINMEAS;
        }
        pSeqExpo->setICEProgramParam(ICE_PROGRAM_PARA_MULTI_OFFLINE_CALL, lICEPara);

    }
    else  // PATMode is 'none':
    {
        // set YAPS parameters
        // reset to default, in case they were previously changed for iPAT
        pSeqExpo->setFirstFourierLine(0);
        pSeqExpo->setFirstRefLine(-1);
        pSeqExpo->setFirstFourierPartition(0);
        pSeqExpo->setFirstRefPartition(-1);

        // elliptical scanning: reset  ICE_PROGRAM_MSK_IPAT_TRIGGER_LASTSCANINMEAS
        int32_t lICEPara = (int32_t)pSeqExpo->getICEProgramParam(ICE_PROGRAM_PARA_MULTI_OFFLINE_CALL) & ~ICE_PROGRAM_MSK_IPAT_TRIGGER_LASTSCANINMEAS;
        pSeqExpo->setICEProgramParam(ICE_PROGRAM_PARA_MULTI_OFFLINE_CALL, lICEPara);
    }

    /*
    if (! _pSeqLim->isContextPrepForBinarySearch() )  {
    {
    cout  << endl
    << "pSeqExpo->getFirstFourierLine     ():" <<      pSeqExpo->getFirstFourierLine     ()<<endl
    << "pSeqExpo->getFirstRefLine         ():" <<      pSeqExpo->getFirstRefLine         ()<<endl
    << "pSeqExpo->getFirstFourierPartition():" <<      pSeqExpo->getFirstFourierPartition()<<endl
    << "pSeqExpo->getFirstRefPartition    ():" <<      pSeqExpo->getFirstRefPartition    ()<<endl
    << "pReorderInfo->getKSCenterLin      ():" <<      pReorderInfo->getKSCenterLin      ()<<endl
    << "pReorderInfo->getKSCenterPar      ():" <<      pReorderInfo->getKSCenterPar      ()<<endl
    << endl;
    }
    */

    return MRI_SEQ_SEQU_NORMAL;
}


bool SeqIF_UIHandlersPAT::PATInitializeReorderMethod(MrProtocolData::MrProtData*  pMrProt, SeqLim* _pSeqLim)
{
#ifdef WIN32
    MrProt sMrProt( pMrProt );
    return PATInitializeReorderMethod ( sMrProt, *_pSeqLim );
#else
    return true;
#endif
}

}
