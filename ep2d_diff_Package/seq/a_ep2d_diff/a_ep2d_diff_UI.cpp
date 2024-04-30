//-----------------------------------------------------------------------------
// <copyright file="a_ep2d_diff_UI.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
//-----------------------------------------------------------------------------

// --------------------------------------------------------------------------
// General Includes
// --------------------------------------------------------------------------
#ifdef WIN32
// MrProt
#include "MrProtSrv/Domain/MrProtocol/UILink/MrStdNameTags.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkSelection.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkArray.h"
#include "MrProtSrv/Domain/MrProtocol/StdProtRes/StdProtRes.h"
#include "MrProtSrv/Domain/MrProtocol/UILink/StdRoutines.h"              // _formatFloat
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/SetProtocolParameter.h"
#include "MrProtSrv/Domain/CoreNative/MrWipMemBlock.h"
#include "MrProtSrv/Domain/CoreNative/MrNavigator.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Physiology/MrPhysiology.h"

#include "MrMeasSrv/SeqFW/SeqInfo/SeqInfoBlockGuard.h"

//#include "MrProtSrv/Domain/MrProtocol/UILink/IFix/IFixConc.h"
#include "MrProtSrv/Domain/MrProtocol/UILink/MrUILinkPhysio.h"

#include "MrProtSrv/Domain/CoreNative/MrSliceGroupData.h"

// MrProt

#include  <vector>
#endif

// MrProt facade support
#include "MrImaging/seq/common/MrProtFacade/IMrProtFacade.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

//  -------------------------------------------------------------------------- 
//  Application includes                                                       
//  -------------------------------------------------------------------------- 
#include "MrImaging/seq/a_ep2d.h"
#include "MrImaging/seq/a_ep_CommonUI.h"
#include "MrImaging/seq/a_ep2d_diff/a_ep2d_diff_UI.h"

#ifdef WIN32    
//  -----------------------------------------------------------------
//  Used Interfaces
//
#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Base.h"   // WIP type definitions
#include "MrImaging/seq/a_ep2d_diff/didi.h"

//  IMULT 
#include "MrImaging/SequenceLibraries/libPace/MODULE/MODULE_Routines.h"

//namespace for XProtocol
using namespace MED_MRES_XProtocol;


//  -----------------------------------------------------------------
//  Exported interface
//  -----------------------------------------------------------------

#ifdef BUILD_SEQU
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"  // import/export control

#endif // WIN32

using namespace SEQ_NAMESPACE;


#ifdef WIN32
using namespace UICtrl;

#ifdef BUILD_WIPParameterTool
//---------------------------------------------------------------------------
// _solveDoubleParamConflict
//---------------------------------------------------------------------------
unsigned _solveDoubleParamConflict(
    MrUILinkLimited<double>* const    _this,
    char**                            arg_list,
    const void*                       pVoid,
    const MrProtocolData::MrProtData* pOrigProt,
    int32_t                           lIndex)
{
    SEQ_TRACE_DEBUG.print("_solveDoubleParamConflict: called, _this->value(%d)=%f", lIndex, _this->value(lIndex));
    return EpCommonUINS::fEPIStdSolveDoubleParam(_this, arg_list, pVoid, pOrigProt, lIndex);
}

//---------------------------------------------------------------------------
// _solveBoolParamConflict
//---------------------------------------------------------------------------
unsigned _solveBoolParamConflict(
    MrUILinkSelection<bool>* const    _this,
    char**                            arg_list,
    const void*                       pVoid,
    const MrProtocolData::MrProtData* pOrigProt,
    int32_t                           lIndex)
{
    SEQ_TRACE_DEBUG.print("_solveBoolParamConflict: called, _this->value(%d)=", lIndex, _this->value(lIndex));
    return EpCommonUINS::fEPIStdSolveBoolParamConflict(_this, arg_list, pVoid, pOrigProt, lIndex);
}
#endif

// ===========================================================================
/// Solve handler for 'DiffusionDirections' parameter
/**
This solve handler is called to allow the number of diffusion
directions to be increased although TE is too short.
The solve strategy in this parameter conflict is to call the standard
solve handler which tries to increase TE and TR.
*/
// ===========================================================================
unsigned Ep2d_diff_UINS::fUILinkNoDiffDirectionsSolve
(
LINK_LONG_TYPE*      const _this,             /**< this pointer                                    */
char**               arg_list,                /**< receives confirmation message                   */
const void*          pAddMem,                 /**< additional memory needed for user b value array */
const MrProtocolData::MrProtData*  pOrigProt, /**< Original protocol with old rf mode              */
int32_t                 lIndex                   /**< Array index reserved                            */
)
// ===========================================================================
{
    return (fUICSolveLongParamConflict(_this, arg_list, pAddMem, pOrigProt, lIndex, NULL, NULL, NULL));
}

// ===========================================================================
/// Tooltip handler for 'DiffusionDirections' parameter
/**
In FREE mode, show user comment as tooltip
*/
// ===========================================================================
unsigned Ep2d_diff_UINS::fUILinkDiffDirectionsGetToolTipId
(
LINK_LONG_TYPE* const _this,       /**< this pointer                                    */
char*                 arg_list[],
int32_t                  lIndex       /**< Array index reserved                            */
)
{
    MrProt rMrProt(&_this->prot());

    if(rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_FREE)
    {
        static char pszToolTip[1024];

        // Prepare tooltip (take care of proper termination in case of clipping)
        strncpy(pszToolTip, getUI(_this)->m_sToolTipString[Ep2d_diff_UI::ToolTipStringDVSInfo].c_str(), 1023);
        pszToolTip[1023] = '\0';

        arg_list[0] = pszToolTip;
        return MRI_STD_STRING;
    }

    return 0;
}



// ===========================================================================
/// SetValue handler for 'Diffusion Mode' parameter
/**
Calls original SetValue handler and takes care of the following settings:
- For NOT (MDDW, QSPACE or FREE), the MultipleSeriesMode must be MULTIPLE_SERIES_EACH_MEASUREMENT.
- For FREE mode, the number of directions is set to the first vector set
found in the external definition file, reconstruction of 'diffusion weighted images'
will be switched on and all other reconstructions will be switched off.
- For MDDW mode, the number of directions is set to the default value
- For QSPACE mode, the number of weightings is set to the default value, and
- For 1-scan trace, SMS is not allowed
*/
// ===========================================================================
unsigned Ep2d_diff_UINS::fUILinkDiffusionModeSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uValue, int32_t lIndex)
{
    MrProt rMrProt(pThis->prot());

    if (rMrProt.getsDiffusion().getulMode() == SEQ::DIFFMODE_QSPACE
        || rMrProt.getsDiffusion().getulMode() == SEQ::DIFFMODE_TENSOR
        || rMrProt.getsDiffusion().getulMode() == SEQ::DIFFMODE_FREE)
    {
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_FLT_DISCOR_PROP_MODE, MRI_STD_DISTCORR_2D, 0);
    }

    // Call original handler
    getUI(pThis)->m_DiffMode.getOrigSetValueHandler()(pThis, uValue, lIndex);

    MrProtFacade protFacade(rMrProt);

    // Apply dedicated modifications
    switch(uValue)
    {
        case MRI_STD_DIFFUSION_TRACE:
            // SMS with 1-scan trace is not allowed
            if (protFacade.isSliceAcceleration())
            {
                if (rMrProt.PAT().getlAccelFactPE() > 1)
                {
                    if (!setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_PAT_MODE, MRI_STD_PAT_MODE_GRAPPA, 0))
                        return false;
                }
                else
                {
                    if (!setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_PAT_MODE, MRI_STD_PAT_MODE_NONE, 0))
                        return false;
                }
            }
            // Disable MultipleSeriesMode
            rMrProt.multipleSeriesMode(SEQ::MULTIPLE_SERIES_OFF);
            break;

        case MRI_STD_DIFFUSION_TENSOR:
            // Switch to default directions
            rMrProt.diffusion().setlDiffDirections(MDDW_DEF_DIRECTIONS);

            // turn distortion correction off
            if(!setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_FLT_DISCOR_PROP_MODE, MRI_STD_OFF, 0))
                return false;

            break;
        case MRI_STD_DIFFUSION_FREE:
        {
            // turn distortion correction off
            if (!setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_FLT_DISCOR_PROP_MODE, MRI_STD_OFF, 0))
                return false;

            // Does the protocol already contain a free diffusion vector set?
            int32_t lFreeDiffDirections = rMrProt.diffusion().getsFreeDiffusionData().getlDiffDirections();

            // No: Check whether we should support the FREE mode at all
            if(lFreeDiffDirections == 0)
            {
                // If no Didi has been registered, we cannot get any directions and bail out
                if(!getUI(pThis)->m_pDidi)
                {
                    return pThis->value(lIndex);
                }

                // If no directory containing external diffusion vector sets exists bail out
                if(!getUI(pThis)->m_pDidi->isExtDiffDirPresent())
                {
                    return pThis->value(lIndex);
                }

                // FREE mode is supported, but not directions imported so far: use fallback solution
                lFreeDiffDirections = FALLBACK_DIRECTIONS;

                // Store fallback vector set to protocol
                if(!StoreDiffusionDataToProt(pThis, lFreeDiffDirections))
                {
                    // No success
                    return pThis->value(lIndex);
                }
            }

            // Use stored vector set
            rMrProt.diffusion().setlDiffDirections(lFreeDiffDirections);
        }
            break;
        case MRI_STD_DIFFUSION_QSPACE:
            // Switch to default number of weightings
            rMrProt.diffusion().setlQSpaceSteps(QSPACE_DEF_STEPS);

            // turn distortion correction off
            if (!setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_FLT_DISCOR_PROP_MODE, MRI_STD_OFF, 0))
                return false;

            break;
        default:
            // Disable MultipleSeriesMode
            rMrProt.multipleSeriesMode(SEQ::MULTIPLE_SERIES_OFF);
            break;
    }

    // set B1 control loop active / inactive
    // this has to be done every time the diffusion mode is changed as the
    // logic to turn the B1 CL on / off depends on the diffusion mode
#ifdef WIN32
    SEQ_NAMESPACE::EpCommonUINS::updateB1ControlLoopParametersInProtocol(rMrProt);
#endif // #ifdef WIN32

    // Perform TE minimization without calling additional SetValue handlers 
    // (template mechanism requires casting of NULL pointer to desired function pointer)
    return (basicSetValue(pThis, pThis->value(lIndex), lIndex, static_cast<unsigned(*)(LINK_SELECTION_TYPE* const, unsigned, int32_t)>(NULL)));
}

// ===========================================================================
/// GetOptions handler for 'Diffusion Mode' parameter
/**
Calls original GetOptions handler and removes FREE mode from option
vector if neither the current protocol includes a user defined vector
set nor a corresponding import directory exists.
*/
// ===========================================================================
bool Ep2d_diff_UINS::fUILinkDiffusionModeGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos)
{
    MrProt rMrProt(pThis->prot());

    // Call original handler
    getUI(pThis)->m_DiffMode.getOrigGetOptionsHandler()(pThis, rOptionVector, rulVerify, pos);

    // Check wheter FREE mode should be available
    bool bFreeModeAvailable = false;

    // Does the protocol already contain a free diffusion vector set?
    if(rMrProt.diffusion().getsFreeDiffusionData().getlDiffDirections() != 0)
    {
        bFreeModeAvailable = true;
    }

    // Does a directory containing external diffusion vector sets exist?
    if(getUI(pThis)->m_pDidi)
    {
        if(getUI(pThis)->m_pDidi->isExtDiffDirPresent())
        {
            bFreeModeAvailable = true;
        }
    }

    // Remove FREE mode from available options
    if(!bFreeModeAvailable)
    {
        std::vector<unsigned>::iterator itOption_Cur = rOptionVector.begin();

        while(itOption_Cur != rOptionVector.end())
        {
            if(*itOption_Cur == MRI_STD_DIFFUSION_FREE)
            {
                rOptionVector.erase(itOption_Cur);
                break;
            }
            itOption_Cur++;
        }
    }

    return rOptionVector.size() > 0;
}

unsigned Ep2d_diff_UINS::fUILinkSolveDiffusionScheme(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pAddMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    // standard solving mechanism
    unsigned int uiResult = EpCommonUINS::fEPIStdSolveSelection(pThis, arg_list, pAddMem, pOrigProt, lIndex);

    if(uiResult != 0)
        return uiResult;


    // if not successful: reset TX acceleration to 1.0
    MrProt rMrProt(pThis->prot());
    MrProtFacade protFacade(rMrProt);

    if(protFacade.isAcceleratedZOOMit())
        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_PTX_PULSES_ARRAY, MR_TAG_PTX_TX_ACC, 1.0);

    if(pThis->sequence().prepareForBinarySearch(rMrProt))
        return (MRI_STD_CONFIRMATION_MSG);


    // no solution found
    return 0;
}

unsigned Ep2d_diff_UINS::fUILinkSolvePTXAcceleration(LINK_DOUBLE_TYPE* const pThis, char** arg_list, const void* pAddMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    // standard solving mechanism
    unsigned int uiResult = getUI(pThis)->m_PTXPulseTxAccDiff.getOrigSolveHandler() (pThis, arg_list, pAddMem, pOrigProt, lIndex);

    if(uiResult != 0)
        return uiResult;


    // if not successful: set diffusion scheme to bipolar
    setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_DIFF_SCHEME, MRI_STD_DIFFSCHEME_BIPOLAR);

    MrProt rMrProt(pThis->prot());
    if(pThis->sequence().prepareForBinarySearch(rMrProt))
        return (MRI_STD_CONFIRMATION_MSG);


    // no solution found
    return 0;
}


// ===========================================================================
/// Tooltip handler for 'Diffusion Scheme' parameter
/**
Show delta (duration of diffusion encoding gradients) and Delta (distance
of diffusion encoding gradients) if the diffusion scheme 'monopolar' is
selected.
*/
// ===========================================================================
unsigned Ep2d_diff_UINS::fUILinkDiffusionSchemeGetToolTipId(LINK_SELECTION_TYPE* const pThis, char* arg_list[], int32_t /* lIndex */)
{
    MrProt sMrProt(&pThis->prot());

    if(!pThis->sequence().prepareForBinarySearch(sMrProt))
    {
        // can't continue, should never happen
        return 0;
    }

    // Tooltip only makes sense for diffusion schemes 'monopolar' and 'STEAM'
    if((sMrProt.diffusion().getdsScheme() != SEQ::DIFFSCHEME_MONOPOLAR) &&
       (sMrProt.diffusion().getdsScheme() != SEQ::DIFFSCHEME_STEAM)
       )
    {
        return 0;
    }

    double dDiffGradDuration = getUI(pThis)->m_dToolTipParam[Ep2d_diff_UI::ToolTipParamDiffGradDuration];
    double dDiffGradSpacing = getUI(pThis)->m_dToolTipParam[Ep2d_diff_UI::ToolTipParamDiffGradSpacing];

    // Limit parameter range to [0, 99999ms]
    // (Available number of digits is restricted) 
    if((dDiffGradDuration < 0.) || (dDiffGradDuration > 99999.) ||
       (dDiffGradSpacing  < 0.) || (dDiffGradSpacing  > 99999.)
       )
    {
        return 0;
    }

    // Assemble tooltip: 'delta = XX.Xms   Delta = YY.Yms' using greek characters
    // "%1!r! = %2!r! %12!r!\t%13!r! = %14!r! %22!r!"
    arg_list[0] = (char*)MRI_STD_DELTA_LOWER_CASE;
    arg_list[1] = (char*)static_cast<int64_t>(_formatFloat(arg_list + 2, dDiffGradDuration, 1));
    arg_list[11] = (char*)MRI_STD_UNIT_MS;

    arg_list[12] = (char*)MRI_STD_DELTA_UPPER_CASE;
    arg_list[13] = (char*)static_cast<int64_t>(_formatFloat(arg_list + 14, dDiffGradSpacing, 1));
    arg_list[21] = (char*)MRI_STD_UNIT_MS;

    return MRI_STD_DIFFUSION_SCHEME_TOOLTIP;
}


// ------------------------------------------------------------------------------
// Function    : fEPIDiffSolveSelectionConflict
// ------------------------------------------------------------------------------
//
// Description : Deals with simple incompatibilities of selection parameters.
//               Calls corresponding common UI handler at the end.
//
// Return      : 0, if no solution possible
//               MRI_STD_STRING on success. The text in arg_list is then used
//               to format the confirmation message and the popup.
//
// ------------------------------------------------------------------------------
unsigned Ep2d_diff_UINS::fEPIDiffSolveSelectionConflict(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    MrProt rOrigProt(pOrigProt);
    MrProt rNewProt(pThis->prot());

    // --------------------------------------------------------------------------------------
    // Handle multiple concats constrains
    // --------------------------------------------------------------------------------------

    if(rNewProt.concatenations() > 1)
    {
        bool bMultiConcatsAllowed = false;

        // Multiple concatenations are allowed if int32_t TR triggering mode is enabled
        // and standard triggering is active

        if(SEQ_NAMESPACE::Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis))
        {
            if(pEPI->m_mySeqLoop.isLongTRTrigMode())
            {
                SEQ::PhysioSignal FirstSignal;
                SEQ::PhysioMethod FirstMethod;
                SEQ::PhysioSignal SecondSignal;
                SEQ::PhysioMethod SecondMethod;

                rNewProt.physiology().getPhysioMode(FirstSignal, FirstMethod, SecondSignal, SecondMethod);

                if(FirstMethod == SEQ::METHOD_TRIGGERING)
                {
                    bMultiConcatsAllowed = true;
                }
            }
        }

        // Multiple concatenations are allowed if navigator triggering is active
#ifdef SUPPORT_PACE
        if(rNewProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF)
        {
            bMultiConcatsAllowed = true;
        }
#endif

        if((rNewProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rNewProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
        {
            bMultiConcatsAllowed = true;
        }

        if(bMultiConcatsAllowed == false)
        {
            if(LINK_LONG_TYPE* pConc = _search<LINK_LONG_TYPE>(pThis, MR_TAG_CONCATENATIONS))
            {
                if(pConc->isAvailable(0))
                {
                    pConc->value(1, 0, MrUILinkBase::SET_MODE::SET_FORCED);
                }
            }
        }
    }

    // --------------------------------------------------------------------------------------
    // Solve conflict between slice selective inversion and non-interleaved acquisition order
    // --------------------------------------------------------------------------------------
    if(!(rNewProt.sliceSeries().mode() == SEQ::INTERLEAVED) && (rOrigProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE))
    {
        // IR is on and the user wished to select non-interleaved acquisition order: 
        // IR should be turned off
        if(!setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_INVERSION, MRI_STD_MAGN_PREP_NONE))
            return 0;

        // Calls prepareForBinarySearch by default
        if(pThis->sequence().prepareForBinarySearch(rNewProt))
        {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }
    else if(!(rOrigProt.sliceSeries().mode() == SEQ::INTERLEAVED) && (rNewProt.preparationPulses().getucInversion() == SEQ::SLICE_SELECTIVE))
    {
        // Non-interleaved acquisition order is on and the user wishes to select IR:
        // acquisition order should be switched to interleaved
        if(!setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_SERIES_MODE, MRI_STD_SERIES_INTERLEAVED))
            return 0;

        // Calls prepareForBinarySearch by default
        if(pThis->sequence().prepareForBinarySearch(rNewProt))
        {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }

    // --------------------------------------------------------------------------------------
    // Solve conflict between diffusion mode ONE_SCAN_TRACE and dynamic field correction
    // --------------------------------------------------------------------------------------
    if((rNewProt.diffusion().getulMode() == SEQ::DIFFMODE_ONE_SCAN_TRACE) && !(rOrigProt.getsDynDistCorrFilter().getucMode() == SEQ::DYN_DISTCORR_NONE))
    {
        // Dynamic distortion correction is on and the user wishes to select ONE_SCAN_TRACE: 
        // Dynamic distortion correction should be turned off
        //
        // Note: by purpose no solve handling for the opposite case is implemented -
        //       an appropriate getOptions handler is registered.
        if(!setProtocolParameter<LINK_BOOL_TYPE>(pThis, MR_TAG_FLT_DYN_DISCOR, false))
            return 0;

        // Calls prepareForBinarySearch by default
        if(pThis->sequence().prepareForBinarySearch(rNewProt))
        {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }

    // --------------------------------------------------------------------------------------
    // Solve conflict between PACE triggering and dynamic field correction
    // --------------------------------------------------------------------------------------
    if((rNewProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF) && !(rOrigProt.getsDynDistCorrFilter().getucMode() == SEQ::DYN_DISTCORR_NONE))
    {
        // Dynamic distortion correction is on and the user wishes to enable PACE: 
        // Dynamic distortion correction should be turned off
        //
        // Note: by purpose no solve handling for the opposite case is implemented -
        //       an appropriate getOptions handler is registered.
        if(!setProtocolParameter<LINK_BOOL_TYPE>(pThis, MR_TAG_FLT_DYN_DISCOR, false))
            return 0;

        // Calls prepareForBinarySearch by default
        if(pThis->sequence().prepareForBinarySearch(rNewProt))
        {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }

    // Call EPI common UI solve handler for selection conflicts, the handler however 
    // does not perform correct check of the protocol so that we need to call a prepareForBinarySearch
    // afterwards
    unsigned uReturn = (EpCommonUINS::fEPIStdSolveSelection(pThis, arg_list, pVoid, pOrigProt, lIndex));
    if(uReturn != 0)
    {
        // Calls prepareForBinarySearch by default
        if(pThis->sequence().prepareForBinarySearch(rNewProt))
        {
            return (MRI_STD_CONFIRMATION_MSG);
        }
    }

    return 0;
}

// ===========================================================================
/// GetValue handler for 'Diffusion Directions Import' parameter
/**
This getValue handler provides the file path for the diffusion
directions import dialogue if diffusion mode FREE is selected.
*/
// ===========================================================================
unsigned Ep2d_diff_UINS::fUILinkDiffDirsImportGetValue(MrUILinkString* const pThis, char* arg_list[], int32_t pos)
{
    static char pszTargetDir[256];

    // Assemble mask for file open dialogue
    std::stringstream sVectorFileMask;
    char* ptCustomerSeq = NULL;
    ptCustomerSeq = getenv(VECTORFILEENV);

    if(ptCustomerSeq == NULL)
    {
        sVectorFileMask << "ThisFileDoesNotExist";
    }
    else
    {
        // "%CustomerSeq%/DiffusionVectorSets/*.dvs"
        sVectorFileMask << ptCustomerSeq << "\\" << VECTORFILEPATH << "\\" << VECTORFILEMASK;
    }

    // Copy path (take care of proper termination in case of clipping)
    strncpy(pszTargetDir, sVectorFileMask.str().c_str(), 255);
    pszTargetDir[255] = '\0';

    arg_list[0] = pszTargetDir;
    return MRI_STD_STRING;
}

// ===========================================================================
/// SetValue handler for 'Diffusion Directions Import' parameter
/**
Require Didi to load diffusion vector sets defined by the user
in the file name given.
*/
// ===========================================================================
void Ep2d_diff_UINS::fUILinkDiffDirsImportSetValue(MrUILinkString* const pThis, const char* newVal, int32_t lPos)
{
    const char pszInternalError[] = "Internal error: Reading diffusion directions from protocol failed!";

    MrProt rMrProt(pThis->prot());

    LINK_LONG_TYPE* pDiffDirections = _search<LINK_LONG_TYPE>(pThis, MR_TAG_NUMBER_DIFF_DIRS);

    // If no Didi has been registered, we cannot get any directions and bail out:
    if(!getUI(pThis)->m_pDidi || !pDiffDirections || !pDiffDirections->isEditable(0))
    {
        // Set error status and store error message
        getUI(pThis)->m_bImportExportError = true;
        getUI(pThis)->m_sImportExportError = pszInternalError;

        return;
    }

    // Set import status in case of error
    getUI(pThis)->m_bVectorFileImported = false;

    // Fallback: If the import fails, switch to 
    // - internal fallback diffusion direction set OR
    // - previously imported direction set with same number of directions
    pDiffDirections->value(FALLBACK_DIRECTIONS, 0);

    // Assure protocol consistency: Store corresponding content
    if(!StoreDiffusionDataToProt(pThis, FALLBACK_DIRECTIONS))
    {
        // Set error status and store error message
        getUI(pThis)->m_bImportExportError = true;
        getUI(pThis)->m_sImportExportError = pszInternalError;

        // No success: set invalid filename
        getUI(pThis)->m_pDidi->setVectorFileName("");
        return;
    }

    // Provide Didi with file name selected by the user and scan for valid diffusion vector sets.
    std::string sFileName(newVal);
    if(!getUI(pThis)->m_pDidi->setVectorFileName(sFileName.substr(sFileName.find_last_of("/\\") + 1).c_str()))
    {
        // Set error status and store error message
        getUI(pThis)->m_bImportExportError = true;
        getUI(pThis)->m_sImportExportError = getUI(pThis)->m_pDidi->getErrorMessage();

        // No success: set invalid filename
        getUI(pThis)->m_pDidi->setVectorFileName("");
        return;
    }

    // Check validity of all contained diffusion vector sets
    for(int32_t lI = 1; lI <= MAX_DIRECTIONS; lI++)
    {
        if(getUI(pThis)->m_pDidi->isDirectionExternal(lI))
        {
            if(!getUI(pThis)->m_pDidi->prepExternal(lI, true))
            {
                // Set error status and store error message
                getUI(pThis)->m_bImportExportError = true;
                getUI(pThis)->m_sImportExportError = getUI(pThis)->m_pDidi->getErrorMessage();

                // No success: set invalid filename
                getUI(pThis)->m_pDidi->setVectorFileName("");
                return;
            }
        }
    }

    // Assure protocol consistency: Store corresponding (possibly just imported) content
    if(!StoreDiffusionDataToProt(pThis, FALLBACK_DIRECTIONS))
    {
        // Set error status and store error message
        getUI(pThis)->m_bImportExportError = true;
        getUI(pThis)->m_sImportExportError = pszInternalError;

        // No success: set invalid filename
        getUI(pThis)->m_pDidi->setVectorFileName("");
        return;
    }

    // Add parameter dependencies
    pThis->addDependentParamPtr(pDiffDirections, 0);

    // Success: set import status
    getUI(pThis)->m_bVectorFileImported = true;

    // * ----------------------------------------------------------------------- *
    // * Take care of TE minimzation (template not possible for MrUILinkString*) *
    // * ----------------------------------------------------------------------- *
    // Only if TOM is enabled
    if(rMrProt.TOM() == SEQ::TOM_MINIMIZE_TE)
    {
        int32_t         lTEMin           = 0;
        int32_t         lTRMin           = 0;

        // * -------------------------------------------------------------------- *
        // * Calculate the min. TE and TR for the selected parameter set          *
        // * -------------------------------------------------------------------- *
        pThis->sequence().prepareForBinarySearch(rMrProt);

        for(int32_t lI = 0; lI < rMrProt.contrasts(); lI++)
        {
            lTEMin = getUI(pThis)->m_lNeededTE;

            // * -------------------------------------------------------------------- *
            // * Set the minimum TE                                                   *
            // * -------------------------------------------------------------------- *
            (*TESetValue)(pThis, lTEMin / 1000.0, lI);
        }

        // * -------------------------------------------------------------------- *
        // * Set TR to the minimum TR value that is necessary to update TE        *
        // * -------------------------------------------------------------------- *
        lTRMin = getUI(pThis)->m_lNeededTR;

        if(lTRMin > rMrProt.tr()[0])
        {
            double dTRMin_ms = static_cast<double>(lTRMin / 1000.0);     // * TR in ms *
            fUICSetWithinLimits(pThis, dTRMin_ms, MR_TAG_TR, UIC_SET_WITHIN_LIMITS_ROUNDUP | UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, 0);
        }
    }
}

// ===========================================================================
/// IsAvailable handler for 'Diffusion Directions Import' parameter
/**
Parameter is visible if FREE mode is selected AND the directory
containing external diffusion vector sets exists.
*/
// ===========================================================================
bool Ep2d_diff_UINS::fUILinkDiffDirsImportIsAvailable(MrUILinkString* const pThis, int32_t lIndex)
{
    // If no Didi has been registered, we cannot get any directions and bail out:
    if(!getUI(pThis)->m_pDidi)
    {
        return false;
    }

    // If Didi reports that directory containing external diffusion directions does not exist: import not possible
    if(!getUI(pThis)->m_pDidi->isExtDiffDirPresent())
    {
        return false;
    }

    // Call orignal IsAvailable handler
    return getUI(pThis)->m_DiffDirsImport.getOrigIsAvailableHandler()(pThis, lIndex);
}

// ===========================================================================
/// Solve handler for 'Diffusion Directions Import' parameter
/**
Parameter is visible if FREE mode is selected AND the directory
containing external diffusion vector sets exists.
*/
// ===========================================================================
unsigned Ep2d_diff_UINS::fUILinkDiffDirsImportSolve(MrUILinkString* const pThis, char** arg_list, const void*, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    static char pszErrorMessage[256];

    // Check error status
    if(getUI(pThis)->m_bImportExportError == false)
    {
        return 0;
    }

    // Copy error message
    strcpy(pszErrorMessage, "$&OK$$");
    strncpy(pszErrorMessage + strlen("$&OK$$"), getUI(pThis)->m_sImportExportError.c_str(), 255 - strlen("$&OK$$"));
    pszErrorMessage[255] = '\0';
    arg_list[0] = pszErrorMessage;

    // Reset error status
    getUI(pThis)->m_bImportExportError = false;

    return MRI_STD_STRING;
}

// ===========================================================================
/// GetValue handler for 'Diffusion Directions Export' parameter
/**
This getValue handler provides the file path for the diffusion
directions export dialogue if diffusion mode FREE is selected.
*/
// ===========================================================================
unsigned Ep2d_diff_UINS::fUILinkDiffDirsExportGetValue(MrUILinkString* const pThis, char* arg_list[], int32_t pos)
{
    static char pszTargetDir[256];

    // Assemble mask for file open dialogue
    std::stringstream sVectorFileMask;
    char* ptCustomerSeq = NULL;
    ptCustomerSeq = getenv(VECTORFILEENV);

    if(ptCustomerSeq == NULL)
    {
        sVectorFileMask << "ThisFileDoesNotExist";
    }
    else
    {
        // "%CustomerSeq%/DiffusionVectorSets/*.dvs"
        sVectorFileMask << ptCustomerSeq << "\\" << VECTORFILEPATH << "\\" << VECTORFILEMASK;
    }

    // Copy path (take care of proper termination in case of clipping)
    strncpy(pszTargetDir, sVectorFileMask.str().c_str(), 255);
    pszTargetDir[255] = '\0';

    arg_list[0] = pszTargetDir;
    return MRI_STD_STRING;
}

// ===========================================================================
/// SetValue handler for 'Diffusion Directions Export' parameter
/**
Stores free diffusion vector set currently selected to file
specified by the user.
*/
// ===========================================================================
void Ep2d_diff_UINS::fUILinkDiffDirsExportSetValue(MrUILinkString* const pThis, const char* newVal, int32_t pos)
{
    const char pszInternalError[] = "Internal error: Storing diffusion directions to protocol failed!";

    MrProt rMrProt(pThis->prot());

    // If no Didi has been registered, we cannot get any directions and bail out:
    if(!getUI(pThis)->m_pDidi)
    {
        // Set error status and store error message
        getUI(pThis)->m_bImportExportError = true;
        getUI(pThis)->m_sImportExportError = pszInternalError;

        return;
    }

    // Prepare Didi with currently selected diffusion vector set
    if(!getUI(pThis)->m_pDidi->prepExternal(rMrProt.diffusion().getlDiffDirections(), true))
    {
        // Set error status and store error message
        getUI(pThis)->m_bImportExportError = true;
        getUI(pThis)->m_sImportExportError = getUI(pThis)->m_pDidi->getErrorMessage();

        return;
    }

    // Export diffusion vector file
    std::string sFileName(newVal);
    if(!getUI(pThis)->m_pDidi->writeToFile(sFileName.substr(sFileName.find_last_of('\\') + 1).c_str()))
    {
        // Set error status and store error message
        getUI(pThis)->m_bImportExportError = true;
        getUI(pThis)->m_sImportExportError = getUI(pThis)->m_pDidi->getErrorMessage();

        return;
    }

    return;
}

// ===========================================================================
/// IsAvailable handler for 'Diffusion Directions Export' parameter
/**
Parameter is visible if FREE mode is selected AND the directory
containing external diffusion vector sets exists.
*/
// ===========================================================================
bool Ep2d_diff_UINS::fUILinkDiffDirsExportIsAvailable(MrUILinkString* const pThis, int32_t lIndex)
{
    // If no Didi has been registered, we cannot get any directions and bail out:
    if(!getUI(pThis)->m_pDidi)
    {
        return false;
    }

    // If Didi reports that directory containing external diffusion directions does not exist: import not possible
    if(!getUI(pThis)->m_pDidi->isExtDiffDirPresent())
    {
        return false;
    }

    // Call orignal IsAvailable handler
    return getUI(pThis)->m_DiffDirsExport.getOrigIsAvailableHandler()(pThis, lIndex);
}

//================================================================================
// UI handlers for setting the content of MrFreeDiffusionData object
//================================================================================
// Note: These handlers are not linked to any UI element. They get just called
//       during protocol conversion and are responsible for copying the content
//       of the free diffusion direction protocol section (MrFreeDiffusionData
//       object) from the original protocol to the new one. While the GetValue
//       handler just returns the actual content of this object, the SetValue
//       handler also ensures a consistent setting of the number of diffusion
//       directions as shown in the UI.

//--------------------------------------------------------------------------------
// Function:        Ep2d_diff_UINS::fUIFreeDiffDataGetValueGeneric()
//
// Source:          MrImaging/seq/a_ep2d_diff/a_ep2d_diff_UI.cpp
//
// Description:     Get-value handler for 'MrFreeDiffusionData' object.
//--------------------------------------------------------------------------------
MrPtr<MrGenericDC::IValueNode> Ep2d_diff_UINS::fUIFreeDiffDataGetValueGeneric(LINK_GENERIC_TYPE* const pThis, int32_t /* lIndex */)
{
    // Create new MrFreeDiffusionData object pointer
    MrProtocolData::MrFreeDiffusionData::Pointer pFreeDiffDirsData = MrProtocolData::MrFreeDiffusionData::create();

    // Fill new MrFreeDiffusionData object with data
    pFreeDiffDirsData->copyFrom(pThis->prot().getsDiffusion().getsFreeDiffusionData());

    return MrPtr<MrGenericDC::IValueNode>(pFreeDiffDirsData.get());
}   // Resolve_UINS::fUIFreeDiffDataGetValueGeneric

//--------------------------------------------------------------------------------
// Function:        Ep2d_diff_UINS::fUIFreeDiffDataSetValueGeneric()
//
// Source:          MrImaging/seq/a_ep2d_diff/a_ep2d_diff_UI.cpp
//
// Description:     Set-value handler for 'MrFreeDiffusionData' object.
//                  Data ist written to  the protocol in any case. If the
//                  diffusion mode 'free' is already selected, the acutal 
//                  number of directions is applied to the corresponding 
//                  protocol parameter.
//--------------------------------------------------------------------------------
MrPtr<MrGenericDC::IValueNode> Ep2d_diff_UINS::fUIFreeDiffDataSetValueGeneric(LINK_GENERIC_TYPE* const pThis, const MrPtr<MrGenericDC::IValueNode>& sNewVal, int32_t lIndex)
{
    // Get MrFreeDiffusionData object from method parameter IItem object
    const MrProtocolData::MrFreeDiffusionData* pFreeDiffData = dynamic_cast<const MrProtocolData::MrFreeDiffusionData*>(sNewVal.get());

    if(NULL != pFreeDiffData)
    {
        // Store content of MrFreeDiffusionData to protocol
        pThis->prot().getsDiffusion().getsFreeDiffusionData().copyFrom(pFreeDiffData);

        LINK_SELECTION_TYPE* pDiffMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_DIFF_MODE, 0, MrUILinkBase::SEARCH_MODE::SEARCH_AVAILABLE);

        // Check wether FREE diffusion mode is selected
        if(pDiffMode && (pDiffMode->value(0) == MRI_STD_DIFFUSION_FREE))
        {
            // Set number of diffusion directions correspondingly
            // Note: Do not use the diffusion directions SetValue hander (fUILinkNumberDiffDirsSetValueNew) 
            //       here, since MrFreeDiffusionData has already been stored in the protocol. By setting the
            //       parameter using the basicSetValue handler, implicit dependencies get considered.
            LINK_LONG_TYPE* pDiffDirs = _search<LINK_LONG_TYPE>(pThis, MR_TAG_NUMBER_DIFF_DIRS, 0, MrUILinkBase::SEARCH_MODE::SEARCH_EDITABLE);
            if(pDiffDirs)
            {
                basicSetValue(pDiffDirs, pThis->prot().getsDiffusion().getsFreeDiffusionData().getlDiffDirections(), 0, getUI(pDiffDirs)->m_NumberDiffDirs.getOrigSetValueHandler());
            }
        }
    }

    return pThis->value(lIndex);
}

// ------------------------------------------------------------------------------
// Function    : fEPIDiffFltDynDistCorrGetOptions
// ------------------------------------------------------------------------------
//
// Description : Deals with inavailabilities of dynamic distortion correction
//
// Return      : true, if a list of options has been prepared
//               false otherwise
//
// ------------------------------------------------------------------------------
bool Ep2d_diff_UINS::fEPIDiffFltDynDistCorrGetOptions(LINK_BOOL_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos)
{
    bool   bRet = false;
    const  UI_ELEMENT_BOOL&  rFilterDynDistCorr = getUI(pThis)->m_FilterDynDistCorr;
    MrProt rMrProt(pThis->prot());

    if(rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_ONE_SCAN_TRACE)
    {
        // No dynamic distortion correction with diffusion mode ONE_SCAN_TRACE

        rulVerify = LINK_BOOL_TYPE::VERIFY_OFF;
        rOptionVector.clear();
        rOptionVector.push_back(false);

        return rOptionVector.size() > 0;
    }
    else
    {
        // In other cases call the standard get-options handler if it has been defined

        if(rFilterDynDistCorr.getOrigGetOptionsHandler())
        {
            bRet = (*rFilterDynDistCorr.getOrigGetOptionsHandler())(pThis, rOptionVector, rulVerify, pos);
        }
    }

    return bRet;
}


// ===========================================================================
/// GetLimits handler for 'Directions' parameter
/**
As the tensor mode provides in general a huge number of possible directions,
a GetLimits handler is necessary to provide an acceptable performance.
The reason is that the directions do not form a contiguous area. Therefore
we cannot use the binary search, must much sample each possibility
individually.

There are three different ways to specify the limits of this parameter:

For the diffusion modes 'tensor' and 'free', there is (in general) a choice
of different directions. The directions are usually provided by Didi,
the diffusion direction provider class.

If the Didi pointer has been initialized correctly, Didi will provide
the possible directions. This will result in a very efficient performance
of the UI.

If the Didi pointer is NULL (which should basically never occur), we use
the old verify-scan-all method, i.e. the framework will check all possible
directions between 1 and MAX_DIRECTIONS. This is very time-consuming.

For all other diffusion modes, there is no choice for the directions,
and the limits will be fixed to to currently selected number of directions.

\return The handler should return true when it has specified at least
one limit interval. It should return false when the parameter has no limits.
**/
// ===========================================================================
bool Ep2d_diff_UINS::fUILinkNoDiffDirectionsGetLimits
(
LINK_LONG_TYPE* const _this,              /*!< IMP: pointer to the UILink object for which the function is called */
std::vector<MrLimitLong>& rLimitVector,   /*!< EXP: vector of limit intervals                                     */
uint32_t& rulVerify,                 /*!< EXP: mode of further processing of the UILink Framework.           */
int32_t                                      /*! unused in this case                                                 */
)
{
    MrProt rMrProt(_this->prot());

    rLimitVector.resize(1);
    rulVerify = LINK_LONG_TYPE::VERIFY_SCAN_ALL;
    const ParLim<int32_t>& _seqLimits = _this->seqLimits().getDiffusionDirections();

    if((rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_TENSOR) ||
       (rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_FREE)
       )
    {
        if(getUI(_this)->m_pDidi)
        {
            MrLimitLong sLim;
            rLimitVector.clear();

            int32_t lFreeDiffDirectionsFromProt = rMrProt.diffusion().getsFreeDiffusionData().getlDiffDirections();

            // Assemble valid directions (internal or external)
            if((rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_FREE)
               && (getUI(_this)->m_bVectorFileImported == false)
               && (lFreeDiffDirectionsFromProt > 0))
            {
                // First case:
                // - free mode selected
                // - no external diffusion vector set imported so far
                // - protocol contains a free diffusion vector set
                // => Allow using the diffusion directions stored in the protocol only
                sLim.setLonely(lFreeDiffDirectionsFromProt);
                rLimitVector.push_back(sLim);
            }
            else
            {
                // Second case:
                // - mddw mode selected 
                // => allow using all internal diffusion vector sets
                // Third case:
                // - free mode selected AND external diffusion vector set imported 
                // => allow using all external diffusion vector sets (and the internal fallback vector set)
                for(int32_t lI = 1; lI <= MAX_DIRECTIONS; lI++)
                {
                    if(((rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_TENSOR) && getUI(_this)->m_pDidi->isDirectionInternal(lI))
                       || ((rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_FREE) && getUI(_this)->m_pDidi->isDirectionExternal(lI) && getUI(_this)->m_bVectorFileImported)
                       || ((rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_FREE) && (lI == FALLBACK_DIRECTIONS) && getUI(_this)->m_bVectorFileImported))
                    {
                        sLim.setLonely(lI);
                        rLimitVector.push_back(sLim);
                    }
                }
            }

            return rLimitVector.size() > 0;
        }
        else
        {
            // default MRES-UI handler
            SEQ_TRACE_ERROR.print("fUILinkNoDiffDirectionsGetLimitsNew ERROR: (static_cast<Ep2d_diff_UI*>(EpCommonUINS::getUI(_this)))->m_pDidi==NULL");
            return rLimitVector[0].setEqualSpaced(_seqLimits.getMin(),
                                                  _seqLimits.getMax(),
                                                  _seqLimits.getInc());
        }
    }
    else
    {
        return rLimitVector[0].setEqualSpaced(rMrProt.diffusion().getlDiffDirections(),
                                              rMrProt.diffusion().getlDiffDirections(),
                                              _seqLimits.getInc());
    }

}


// ===========================================================================
/// GetLimits handler for 'Weightings' parameter
/**
Q-space mode requires special handling, since only dedicated sampling
schemes are supported.
**/
// ===========================================================================
bool Ep2d_diff_UINS::fUILinkNoQSpaceStepsGetLimits
(
LINK_LONG_TYPE* const    /* pThis */,         /*!< IMP: Pointer to the UILink object for which the function is called */
std::vector<MrLimitLong> &rLimitVector,       /*!< EXP: Vector of limit intervals                                     */
uint32_t            &rulVerify,          /*!< EXP: Mode of further processing of the UILink Framework.           */
int32_t                     /* lIndex */         /*!< IMP: Index                                                         */
)
// ===========================================================================
{
    rLimitVector.resize(1);
    rulVerify = LINK_LONG_TYPE::VERIFY_SCAN_ALL;
    return rLimitVector[0].setEqualSpaced(QSPACE_MIN_STEPS, QSPACE_MAX_STEPS, 1);
}

// ------------------------------------------------------------------------------
// Function    : _bValueTry
// ------------------------------------------------------------------------------
//
// Description : calls original try-handler for b-value and dynamically builds up
//               a table which holds valid results (speeds up binary search, when
//               the same protocol - with the exception of the individual
//               b-values (!) - is tried several times)
// Return      : whatever original try-handler says
//
// ------------------------------------------------------------------------------
bool Ep2d_diff_UINS::_bValueTry(LINK_LONG_TYPE* const pThis, void* pVoid, const MrProtocolData::MrProtData*  pOrig, int32_t lIndex)
{
    MrProt                  rTempProt      = pThis->prot().clone();
    int32_t                 lBValue        = 0;
    bool                    bValidProtocol = false;
    int                     iI = 0;
    const UI_ELEMENT_LONG&  rBValue = getUI(pThis)->m_BValue;

    getUI(pThis)->m_bNeedOtherTETITR = false;

    // Protocol comparison has to disregard the actual b-values and number of diffusion weightings
    rTempProt.getsDiffusion().setlDiffWeightings(1);
    rTempProt.getsDiffusion().getalBValue()[0] = 0;

    // For the protocol validity check, the _maximum_ b-value within the list is relevant!
    for(iI = 0; iI < pThis->prot().getsDiffusion().getlDiffWeightings(); ++iI)
    {
        if(pThis->prot().getsDiffusion().getalBValue()[iI] > lBValue)
        {
            lBValue = pThis->prot().getsDiffusion().getalBValue()[iI];
        }
    }

    if(rTempProt == getUI(pThis)->m_sMemProt)
    {
        bool bFound = false;

        // Current protocol was evaluated in advance - did we already check the current b-value?
        iI = 0;

        while((iI < getUI(pThis)->m_iMemCount) && !bFound)
        {
            if(getUI(pThis)->m_asMemory[iI].lBValue == lBValue)
            {
                bFound = true;
            }
            else
            {
                ++iI;
            }
        }

        if(!bFound)
        {
            // Not yet checked
            if(rBValue.getOrigTryHandler())
            {
                bValidProtocol = (*rBValue.getOrigTryHandler())(pThis, pVoid, pOrig, lIndex);;
            }

            // Store result for current try (including values calculated by fSeqPrep)
            if(getUI(pThis)->m_iMemCount < iMaxProtMemory)
            {
                getUI(pThis)->m_asMemory[getUI(pThis)->m_iMemCount].lBValue = lBValue;
                getUI(pThis)->m_asMemory[getUI(pThis)->m_iMemCount].lNeededTE = getUI(pThis)->m_lNeededTE;
                getUI(pThis)->m_asMemory[getUI(pThis)->m_iMemCount].lNeededTR = getUI(pThis)->m_lNeededTR;
                getUI(pThis)->m_asMemory[getUI(pThis)->m_iMemCount].lNeededTI = getUI(pThis)->m_lNeededTI;
                getUI(pThis)->m_asMemory[getUI(pThis)->m_iMemCount].bValidProtocol = bValidProtocol;
                getUI(pThis)->m_asMemory[getUI(pThis)->m_iMemCount].bNeedOtherTETITR = getUI(pThis)->m_bNeedOtherTETITR;

                getUI(pThis)->m_iMemCount++;
            }

            return bValidProtocol;
        }
        else
        {
            // Already checked: recall stored values
            getUI(pThis)->m_lNeededTE = getUI(pThis)->m_asMemory[iI].lNeededTE;
            getUI(pThis)->m_lNeededTI = getUI(pThis)->m_asMemory[iI].lNeededTI;
            getUI(pThis)->m_lNeededTR = getUI(pThis)->m_asMemory[iI].lNeededTR;
            getUI(pThis)->m_bNeedOtherTETITR = getUI(pThis)->m_asMemory[iI].bNeedOtherTETITR;

            return getUI(pThis)->m_asMemory[iI].bValidProtocol;
        }
    }
    else
    {
        // Current protocol is unkown => store
        getUI(pThis)->m_sMemProt = rTempProt;

        // Initialize memory
        for(iI = 0; iI < iMaxProtMemory; ++iI)
        {
            getUI(pThis)->m_asMemory[iI].lBValue = -1;      // Invalid b-value
        }

        if(rBValue.getOrigTryHandler())
        {
            bValidProtocol = (*rBValue.getOrigTryHandler())(pThis, pVoid, pOrig, lIndex);
        }

        // Store result for current try (including values calculated by fSeqPrep)
        getUI(pThis)->m_asMemory[0].lBValue = lBValue;
        getUI(pThis)->m_asMemory[0].lNeededTE = getUI(pThis)->m_lNeededTE;
        getUI(pThis)->m_asMemory[0].lNeededTR = getUI(pThis)->m_lNeededTR;
        getUI(pThis)->m_asMemory[0].lNeededTI = getUI(pThis)->m_lNeededTI;
        getUI(pThis)->m_asMemory[0].bValidProtocol = bValidProtocol;
        getUI(pThis)->m_asMemory[0].bNeedOtherTETITR = getUI(pThis)->m_bNeedOtherTETITR;
        getUI(pThis)->m_iMemCount = 1;

        return bValidProtocol;
    }

    return true;    // Never reached
}

// ------------------------------------------------------------------------------
// Function    : _bandwidthGetLimits
// ------------------------------------------------------------------------------
//
// Description : calls original getLimits-handler for bandwidth, but sets search
//               mode to VERIFY_SCAN_ALL
// Return      : whatever original getLimits-handler says
//
// ------------------------------------------------------------------------------
bool Ep2d_diff_UINS::_bandwidthGetLimits(LINK_DOUBLE_TYPE * const pThis, std::vector< MrLimitDouble >& rLimitVector, uint32_t& rulVerify, int32_t lIndex)
{
    const UI_ELEMENT_DOUBLE&  rBandwidth = getUI(pThis)->m_Bandwidth;
    bool                      bReturn = false;

    if(rBandwidth.getOrigGetLimitsHandler())
    {
        bReturn = (*rBandwidth.getOrigGetLimitsHandler())(pThis, rLimitVector, rulVerify, lIndex);
    }

    rulVerify = LINK_DOUBLE_TYPE::VERIFY_SCAN_ALL;

    return bReturn;
}




// ===========================================================================
/// IsAvailable handler for 'Multiple Series Mode' parameter
/**
This parameter on the Contrast / Dynamic UI card specifies if all
repetitions go into one series (SEQ::MULTIPLE_SERIES_OFF) or are
saved as individual series (SEQ::MULTIPLE_SERIES_EACH_MEASUREMENT).

This option can only be used if the diffusion scans run over the
repetition loop, i.e. for diffusion mode MDDW or FREE.

\return The handler should return true if the UI parameter should
be displayed.
**/
// ===========================================================================
bool Ep2d_diff_UINS::fUILink_MultipleSeriesMode_IsAvailable
(
LINK_SELECTION_TYPE* const pThis, /*!< IMP: pointer to the UILink object for which the function is called */
int32_t                              /*!< lIndex                                                             */
)
{
    MrProt rMrProt(pThis->prot());
    MrProtFacade protFacade(rMrProt);

    return protFacade.iSDTI();
}




/*[ Function ****************************************************************\
*
* Name        : TESetValue
*
* Description : This function sets the TE time of the specified contrast
*
* Return      : New TE[lContrast] time
*
\****************************************************************************/
double Ep2d_diff_UINS::TESetValue(MrUILinkBase* const pThis, double dDesiredTE_ms, int32_t lContrast)
{
    MrProt   rMrProt(pThis->prot());
    SeqLim   &rSeqLim = pThis->seqLimits();

    int32_t     lNewTE_us = shift_to_grid_high(static_cast<int32_t>(dDesiredTE_ms * 1000 + 0.5),
                                               rSeqLim.getTE()[lContrast].getMin(),
                                               rSeqLim.getTE()[lContrast].getMax(),
                                               rSeqLim.getTE()[lContrast].getInc());

    rMrProt.te()[lContrast] = lNewTE_us;

    return (rMrProt.te()[lContrast] / 1000.0);
}

bool Ep2d_diff_UINS::StoreDiffusionDataToProt(MrUILinkBase* const pThis, int32_t lDiffDirs)
{
    MrProt rMrProt(pThis->prot());

    // If no Didi has been registered, we cannot set any directions and bail out
    if(!getUI(pThis)->m_pDidi)
    {
        return false;
    }

    if(getUI(pThis)->m_pDidi->prepExternal(lDiffDirs, true))
    {
        // Number of diffusion directions
        rMrProt.diffusion().getsFreeDiffusionData().setlDiffDirections(lDiffDirs);
        // Coordinate sytstem
        rMrProt.diffusion().getsFreeDiffusionData().setulCoordinateSystem(getUI(pThis)->m_pDidi->getCoordinateSystem());
        // Normalization
        rMrProt.diffusion().getsFreeDiffusionData().setulNormalization(MrProtocolData::DIFFDIR_NORM_NONE);
        // Comment, consisting of the file name (without path) in the first line followed by a user comment
        // Note: It is essential that the file name is stored in the protocol! Otherwise, after
        //       creating a new diffusion class instance (in SBBDiffusion.cpp) the file name gets
        //       lost (since each diffusion class instance recreates its Didi instance).
        std::string sComment;
        sComment.append(getUI(pThis)->m_pDidi->getVectorFileName());
        sComment.append("\n");
        sComment.append(getUI(pThis)->m_pDidi->getComment());
        rMrProt.diffusion().getsFreeDiffusionData().setsComment(sComment.c_str());
        // Diffusion vectors
        rMrProt.diffusion().getsFreeDiffusionData().getasDiffDirVector().resize(lDiffDirs);
        for(int32_t lI = 0; lI < lDiffDirs; ++lI)
        {
            rMrProt.diffusion().getsFreeDiffusionData().getasDiffDirVector()[lI].setdSag(getUI(pThis)->m_pDidi->getX(lI));
            rMrProt.diffusion().getsFreeDiffusionData().getasDiffDirVector()[lI].setdCor(getUI(pThis)->m_pDidi->getY(lI));
            rMrProt.diffusion().getsFreeDiffusionData().getasDiffDirVector()[lI].setdTra(getUI(pThis)->m_pDidi->getZ(lI));
        }

        return true;
    }

    return false;
}

// * -------------------------------------------------------------------------- *
// *                        LINK_LONG_TYPE set value handler                    *
// * -------------------------------------------------------------------------- *

/*[ Function ****************************************************************\
*
* Name        : fUILinkBaseResolutionSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
int32_t Ep2d_diff_UINS::fUILinkBaseResolutionSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_BaseResolution.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkPATFactorSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
int32_t Ep2d_diff_UINS::fUILinkPATFactorSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_PATAccelerationFactor.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkPATRefLinesSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
int32_t Ep2d_diff_UINS::fUILinkPATRefLinesSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_PATReferenceLines.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkNumberDiffDirsSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
int32_t Ep2d_diff_UINS::fUILinkNumberDiffDirsSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    MR_SEQUENCE_EXPORT_TLS(pThis->sequence());
    MrProt rMrProt(pThis->prot());

    // If we are in FREE mode and a valid Didi is registered: 
    // Read diffusion vector set from file and store content to protocol
    if((rMrProt.diffusion().getulMode() == SEQ::DIFFMODE_FREE) && getUI(pThis)->m_pDidi)
    {
        if(!StoreDiffusionDataToProt(pThis, lNewVal))
        {
            // Return old value
            return rMrProt.diffusion().getlDiffDirections();
        }
    }

    // Write new value to protocol and perform TE minimzation
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_NumberDiffDirs.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkBValueSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
int32_t Ep2d_diff_UINS::fUILinkBValueSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    MR_SEQUENCE_EXPORT_TLS(pThis->sequence());
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_BValue.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkBValueIsAvailable
*
* Description : This function calls the original UI Handler.
*
* Return      : New value
*
\****************************************************************************/

bool Ep2d_diff_UINS::fUILinkBValueIsAvailable(LINK_LONG_TYPE* const pThis, int32_t)
{
    MR_SEQUENCE_EXPORT_TLS(pThis->sequence());
    return getUI(pThis)->m_BValue.getOrigIsAvailableHandler();
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkBValueSizeSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/

int32_t Ep2d_diff_UINS::fUILinkBValueSizeSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_BValueSize.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkNumberSlicesSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
int32_t Ep2d_diff_UINS::fUILinkNumberSlicesSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_NumberSlices.getOrigSetValueHandler()));
}


/*[ Function ****************************************************************\
*
* Name        : fUILinkQSpaceStepsSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
int32_t Ep2d_diff_UINS::fUILinkQSpaceStepsSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_QSpaceSteps.getOrigSetValueHandler()));
}


/*[ Function ****************************************************************\
*
* Name        : fUILinkQSpaceMaxBValueSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
int32_t Ep2d_diff_UINS::fUILinkQSpaceMaxBValueSetValueNew(LINK_LONG_TYPE* const pThis, int32_t lNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, lNewVal, lPos, getUI(pThis)->m_QSpaceMaxBValue.getOrigSetValueHandler()));
}



// * -------------------------------------------------------------------------- *
// *                        LINK_DOUBLE_TYPE set value handler                  *
// * -------------------------------------------------------------------------- *

/*[ Function ****************************************************************\
*
* Name        : fUILinkBandwidthSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
double Ep2d_diff_UINS::fUILinkBandwidthSetValueNew(LINK_DOUBLE_TYPE* const pThis, double dNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, dNewVal, lPos, getUI(pThis)->m_Bandwidth.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkPhaseFOVSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
double Ep2d_diff_UINS::fUILinkPhaseFOVSetValueNew(LINK_DOUBLE_TYPE* const pThis, double dNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, dNewVal, lPos, getUI(pThis)->m_PhaseFOV.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkReadFOVSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
double Ep2d_diff_UINS::fUILinkReadFOVSetValueNew(LINK_DOUBLE_TYPE* const pThis, double dNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, dNewVal, lPos, getUI(pThis)->m_ReadFOV.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkPhaseResolutionSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
double Ep2d_diff_UINS::fUILinkPhaseResolutionSetValueNew(LINK_DOUBLE_TYPE* const pThis, double dNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, dNewVal, lPos, getUI(pThis)->m_PhaseResolution.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkSliceThicknessSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
double Ep2d_diff_UINS::fUILinkSliceThicknessSetValueNew(LINK_DOUBLE_TYPE* const pThis, double dNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, dNewVal, lPos, getUI(pThis)->m_SliceThick.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkPhaseOSSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
double Ep2d_diff_UINS::fUILinkPhaseOSSetValueNew(LINK_DOUBLE_TYPE* const pThis, double dNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, dNewVal, lPos, getUI(pThis)->m_PhaseOS.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkEchoSpacingSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
double Ep2d_diff_UINS::fUILinkEchoSpacingSetValueNew(LINK_DOUBLE_TYPE* const pThis, double dNewVal, int32_t lPos)
{
    // Perform TE minimization and call overloaded SetValue handler from base class
    return (basicSetValue(pThis, dNewVal, lPos, getUI(pThis)->m_EchoSpacing.getOrigSetValueHandler()));
}


// * -------------------------------------------------------------------------- *
// *                        LINK_SELECTION_TYPE set value handler               *
// * -------------------------------------------------------------------------- *

/*[ Function ****************************************************************\
*
* Name        : fUILinkPATModeSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkPATModeSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    // mSENSE can only be used with coil combine SOS
    if(uNewVal == MRI_STD_PAT_MODE_SENSE)
        setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_COIL_COMBINE_MODE, MRI_STD_COIL_COMBINE_SUM_OF_SQUARES, 0, true, false);

    // for SMS @ 1.5T we always switch to monopolar as this scheme does not have problems with Maxwell (concomitant fields) terms
    if(SysProperties::getNominalB0() < 2.0)
    {
        if(uNewVal == MRI_STD_PAT_MODE_SLICE_ACCEL)
            setProtocolParameter<LINK_SELECTION_TYPE>(pThis, MR_TAG_DIFF_SCHEME, MRI_STD_DIFFSCHEME_MONOPOLAR, 0, true, false);
    }

    unsigned uResult = (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_PATMode.getOrigSetValueHandler()));

    EpCommonUINS::adaptRefLinesPE(pThis, true);

    return uResult;
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkFatSuppressionSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkFatWaterContrastSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_FatSup.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkRFPulseTypeSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkRFPulseTypeSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_RFPulseType.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkFatSatModeSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkFatSatModeSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_FatSatMode.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkPhasePFSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkPhasePFSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_PhasePartF.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkTOMSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkTOMSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_TOM.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkGradModeSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkGradModeSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_GradMode.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkDiffSchemeSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkDiffSchemeSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_DiffScheme.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkQSpaceCoverageSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkQSpaceCoverageSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_QSpaceCoverage.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkQSpaceSamplingSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkQSpaceSamplingSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_QSpaceSampling.getOrigSetValueHandler()));
}

/*[ Function ****************************************************************\
*
* Name        : fUILinkPatRefScanModeSetValueNew
*
* Description : This function calls the original UI method and sets the minimum TE.
*               TR is only changed if necessary.
*
* Return      : New value
*
\****************************************************************************/
unsigned Ep2d_diff_UINS::fUILinkPatRefScanModeSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_PatRefScanMode.getOrigSetValueHandler()));
}

//  ----------------------------------------------------------------------
//
//  Name        :  fAcquisitionWindowGetLimits
//
//  Description :
//
//  Return      :  true  - if B1 control loop is used
//                 false - else
//
//  ----------------------------------------------------------------------
bool Ep2d_diff_UINS::fAcquisitionWindowGetLimits(LINK_LONG_TYPE * const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t lIndex)
{
    if (SEQ::SIGNAL_RESPIRATION == pThis->prot().getsPhysioImaging().getlSignal1() &&
        SEQ::ACQUISITION_WINDOW_MS == pThis->prot().getsPhysioImaging().getsPhysioResp().getucAcquisitionWindowSelectMode())
    {
        // "calculate" limits for Acq Window for resp triggering, in case acq wnd is provided in ms
        rulVerify = LINK_LONG_TYPE::VERIFY_BINARY_SEARCH;

        MrLimitLong limit;
 
        long high = 10000;
        long low = 10;
        long incr = 10;

        limit.setEqualSpaced(low, high, incr);
        rLimitVector.push_back(limit);

        return true;
    }
    else
    {
        //call original handler
        return getUI(pThis)->m_AcqWindow.getOrigGetLimitsHandler()(pThis, rLimitVector, rulVerify, lIndex);
    }
}

static void adaptNumberOfConcats(MrUILinkBase* const pThis)
{
    MrProt  rProt(&pThis->prot());
    MrProtFacade protFacade(rProt);

    // SMS within EP2d does currently does not support multiple concats, such no adaption of number of concats is possible
    if (!protFacade.isSliceAcceleration())
    {
        MrUILinkLimited<int32_t>*  _pConcats = _search < LINK_LONG_TYPE >(pThis, MR_TAG_CONCATENATIONS);

        if (_pConcats && _pConcats->isAvailable(0) && _pConcats->isEditable(0))
        {
            uint32_t ulVerify = 0;
            std::vector<MrLimitLong>  limits;
            _pConcats->getLimits(limits, ulVerify, 0);
            if (!limits.size())
            {
                return;
            }
            //  try maximum number of concatenations.
            int32_t lMax = limits.back().maximum();

            // create MrProt on heap instead of stack; we use MrSmartPointer to get automatic deletion on block exit
            MrProtocolData::MrProtData::Pointer tempProt = MrProtocolData::MrProtData::create();
            void* pExtra = _pConcats->storeInMemory(tempProt.get(), 0);
            _pConcats->value(lMax, 0);
            bool bResult = _pConcats->tryProt(pExtra, tempProt.get(), 0);

            //  calculate soft minimum
            if (bResult)
            {
                bResult = _pConcats->calcRestrictedLimits(limits, LINK_LONG_TYPE::CALC_SOFT_LIMITS, 0);
            }

            _pConcats->recallMemory(pExtra, tempProt.get(), 0);
            _pConcats->clearMemory(pExtra, 0);

            if (bResult)
            {
                long lMinConc = limits.front().minimum();
                _pConcats->value(lMinConc, 0);
            }
        }
    }
}


int32_t Ep2d_diff_UINS::fAcquisitionWindowSetValue(LINK_LONG_TYPE* const pThis, int32_t lNewScanWindow_ms, int32_t lIndex)
{

    //call original handler
    int32_t result = getUI(pThis)->m_AcqWindow.getOrigSetValueHandler()(pThis, lNewScanWindow_ms, lIndex);

    MrProt  rProt(&pThis->prot());
    if (rProt.getsPhysioImaging().getlSignal1() == SEQ::SIGNAL_RESPIRATION &&
        rProt.getsSliceArray().getucConcatenationsSelectModeResp() == SEQ::CONCAT_AUTOMATIC)
    {
        adaptNumberOfConcats(pThis);
    }

    return result;
}

int32_t Ep2d_diff_UINS::fAcquisitionWindowInternalSetValue(LINK_LONG_TYPE* const pThis, int32_t lNewScanWindow_ms, int32_t lIndex)
{

    //call original handler
    int32_t result = getUI(pThis)->m_AcqWindowInternal.getOrigSetValueHandler()(pThis, lNewScanWindow_ms, lIndex);

    MrProt  rProt(&pThis->prot());
    if (rProt.getsPhysioImaging().getlSignal1() == SEQ::SIGNAL_RESPIRATION &&
        rProt.getsSliceArray().getucConcatenationsSelectModeResp() == SEQ::CONCAT_AUTOMATIC)
    {
        adaptNumberOfConcats(pThis);
    }

    return result;
}


unsigned Ep2d_diff_UINS::fAcquisitionWindowSolve(LINK_LONG_TYPE* const pThis, char** arg_list, const void* pAddMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex)
{
    MrProt  rProt(&pThis->prot());
    if (rProt.getsPhysioImaging().getlSignal1() == SEQ::SIGNAL_RESPIRATION)
    {
        if (rProt.getsSliceArray().getlConc() != rProt.getsSliceArray().getlSize())
        {
            MrUILinkLimited<int32_t>*  _pConcats = _search < LINK_LONG_TYPE >(pThis, MR_TAG_CONCATENATIONS);

            if (_pConcats && _pConcats->isAvailable(0) &&_pConcats->isEditable(0))
            {
                adaptNumberOfConcats(pThis);

                if (rProt.getsSliceArray().getucConcatenationsSelectModeResp() != SEQ::CONCAT_AUTOMATIC)
                {
                    pThis->addDependentParamPtr(_pConcats, 0);
                }
            }
        }
        return MrUILinkBase::stdSolveHandler(pThis, arg_list, pAddMem, pOrigProt, lIndex);
    }
    else
    {
        //call original handler
        return getUI(pThis)->m_AcqWindow.getOrigSolveHandler()(pThis, arg_list, pAddMem, pOrigProt, lIndex);

    }
}

bool     Ep2d_diff_UINS::fSetCapture(LINK_BOOL_TYPE* const pThis, bool bSetCapture, int32_t lIndex)
{
    MrProt  sProt(&pThis->prot());  // MrProt wrapper
    bool bStatus = bSetCapture;

    if (sProt.getsPhysioImaging().getlSignal1() != SEQ::SIGNAL_RESPIRATION)
    {
        //call original handler
        return getUI(pThis)->m_SetCapture.getOrigSetValueHandler()(pThis, bSetCapture, lIndex);
    }
    else
    {
        int32_t minimumAcqWindow_ms = -1;
        LINK_LONG_TYPE* pAcqWnd = _search<LINK_LONG_TYPE>(pThis, MR_TAG_FIRST_ACQUISITION_WINDOW);
        if (pAcqWnd)
        {
            std::vector< MrLimitLong > acqWndLimitVector;
            if (sProt.getsSliceArray().getucConcatenationsSelectModeResp() == SEQ::CONCAT_AUTOMATIC)
                pAcqWnd->calcRestrictedLimits(acqWndLimitVector, LINK_LONG_TYPE::CALC_EXT_LIMITS, lIndex);
            else
                pAcqWnd->calcRestrictedLimits(acqWndLimitVector, LINK_LONG_TYPE::CALC_SOFT_LIMITS, lIndex);
            if (acqWndLimitVector.size() > 0)
                minimumAcqWindow_ms = acqWndLimitVector[0].minimum();
        }

        if (sProt.getsPhysioImaging().getsPhysioResp().getucAcquisitionWindowSelectMode() == SEQ::ACQUISITION_WINDOW_PERCENT)
            bStatus = fUILinkCaptureRRSetValue(pThis, bSetCapture, lIndex, minimumAcqWindow_ms, static_cast<short>(sProt.getsPhysioImaging().getsPhysioResp().getlAcqusitionWindowPercent()));
        else
            bStatus = fUILinkCaptureRRSetValue(pThis, bSetCapture, lIndex, minimumAcqWindow_ms, -1);

        return bStatus;
    }

}

#ifdef BUILD_WIPParameterTool
bool Ep2d_diff_UINS::fUIECCompensationSetValue(LINK_BOOL_TYPE* const pThis, bool value, int32_t lIndex)
{
    return (basicSetValue(pThis, value, lIndex, getUI(pThis)->m_EddyCurrentComp.getOrigSetValueHandler()));
}
#endif

#ifdef WIP


// ===========================================================================
/*!
\author   Michael.Zwanger@med.siemens.de

\brief This function returns the string tag for registering a UI WIP parameter.

To register a UI backend handler for WIP parameters, a tag is necessary to
identify the position of the parameter on the UI sequence special card.
Unfortunately this tag is a string, which is not related to the index of
the WIP parameter in the protocol structure. This function translates the
index of the protocol WIP mem block into a string tag. This is useful
in case of an index conflict: changing the WIP array index will automatically
change the position on the UI card.

Example for usage:
\code
if ( LINK_LONG_TYPE* pLong=_create<LINK_LONG_TYPE>(pSeqLim, fWIPString(WIP_Spoiler), WIP_Spoiler) ) {
pLong->registerGetValueHandler(_WIP_LONG_GetValue);
...
}
\endcode

*/
// ===========================================================================

bool Ep2d_diff_UINS::fWIPString
(
int   index,                     /**< Input: The index of the WIP mem block         */
char  ptWIPString[20],           /**< Output: String tag                            */
char  ptFamily[14]               /**< Input (optional): base tag. Should be one of
                                 "seq_wip"(default), "seq_res", "eva_seq_wip". */
                                 )
{
    sprintf(ptWIPString, "%s%d", ptFamily, index);
    return true;
}



// ---------------------------------------------------------------
// Implementation of UI functions for a 'selection' WIP parameter                              
// ---------------------------------------------------------------


// ===========================================================================
///  Return the text to be placed in front of selection boxes
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_SELECTION_GetLabelId(LINK_SELECTION_TYPE* const /* _this */, char* arg_list[], int32_t lIndex)
{
    static const char pszLabel0[] = "Anylabel";
    static const char pszLabelInvalid[] = "NotSupported";

    switch(lIndex)
    {
        case WIP_Anything:
            arg_list[0] = (char*)pszLabel0;
            break;
        default:
            arg_list[0] = (char*)pszLabelInvalid;
            break;
    }
    return MRI_STD_STRING;
}


// ===========================================================================
/// Selects the text of the alternatives in the selection boxes
// ===========================================================================
int Ep2d_diff_UINS::_WIP_SELECTION_Format(LINK_SELECTION_TYPE* const /* _this */, unsigned nID, char* arg_list[], int32_t lIndex)
{
    static const char pszFormat0[] = "Anyformat";

    unsigned uVal = GET_MODIFIER(nID);

    switch(lIndex)
    {
        case WIP_Anything:
            switch(uVal)
            {
                case WIP_Anything:                        /* WIP_Anyoption */
                    arg_list[0] = (char*)pszFormat0;
                    return 1;
                default:
                    arg_list[0] = NULL;
                    return 0;
            }
    }

    return 0;
}



// ===========================================================================
/// Get the value of the selection box with ID lIndex
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_SELECTION_GetValue(LINK_SELECTION_TYPE* const _this, int32_t lIndex)
{
    unsigned nRet = MRI_STD_STRING;
    SET_MODIFIER(nRet, static_cast<unsigned char>(_this->prot().getsWipMemBlock().getalFree()[lIndex]));
    return nRet;
}


// ===========================================================================
/// Get the possible options of Selection box with ID lIndex
// ===========================================================================
bool Ep2d_diff_UINS::_WIP_SELECTION_GetOptions(LINK_SELECTION_TYPE* const /* _this */, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t lIndex)
{
    switch(lIndex)
    {
        case WIP_Anything:
            rulVerify = LINK_SELECTION_TYPE::VERIFY_ON;
            rOptionVector.resize(1);
            rOptionVector[0] = MRI_STD_STRING;
            SET_MODIFIER(rOptionVector[0], WIP_Anything);   /* WIP_Anyoption */
            return true;
        default:
            break;
    }

    return false;
}


// ===========================================================================
/// Set a specified value to the selection box lIndex
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_SELECTION_SetValue(LINK_SELECTION_TYPE* const _this, unsigned nNewVal, int32_t lIndex)
{
    _this->prot().getsWipMemBlock().getalFree()[lIndex] = GET_MODIFIER(nNewVal);
    return _this->value(lIndex);
}



// ===========================================================================
/// Decide if parameter is available on card
// ===========================================================================
bool Ep2d_diff_UINS::_WIP_SELECTION_IsAvailable(LINK_SELECTION_TYPE* const /* _this */, int32_t lIndex)
{
    switch(lIndex)
    {
        case WIP_Anything:
            return true;
        default:
            return true;
    }
    return true;
}



// ===========================================================================
/// Searches for the correct tool tip text for CheckBoxes
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_SELECTION_GetToolTipId(LINK_SELECTION_TYPE* const /* _this */, char* arg_list[], int32_t lIndex)
{
    char tLine[100];
    char tToolTip[1000];

    tToolTip[0] = '\0';
    switch(lIndex)
    {
        case WIP_Anything:
            sprintf(tLine, "Anytext\n");
            strcat(tToolTip, tLine);
            arg_list[0] = tToolTip;
            return MRI_STD_STRING;
            break;

        default: break;
    }

    return 0;
}




// ---------------------------------------------------------
// Implementation of UI functions for a 'int32_t' WIP parameter                              
// ---------------------------------------------------------


// ===========================================================================
/// Defines the text in front of the box for int32_t values
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_LONG_GetLabelId(LINK_LONG_TYPE* const /* _this */, char* arg_list[], int32_t lIndex)
{
    static const char* const pszLabelInvalid = "NotSupported";

    switch(lIndex) // depending on the Index of available boxes give them a name 
    {
        case WIP_Anything:
        default:
            arg_list[0] = (char*)pszLabelInvalid;
            break;
    }
    return MRI_STD_STRING;
}


// ===========================================================================
/// Defines the text behind the box for int32_t values
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_LONG_GetUnitId(LINK_LONG_TYPE* const /* _this */, char* arg_list[], int32_t lIndex)
{
    static const char* const pszLabelInvalid = "NotSupported";

    // depending on the Index of available int32_t boxes, set the correct text behind them
    switch(lIndex)
    {
        case WIP_Anything:
        default:
            arg_list[0] = (char*)pszLabelInvalid;
            break;
    }
    return MRI_STD_STRING;
}



// ===========================================================================
/// Get the int32_t value entered in the box
// ===========================================================================
int32_t Ep2d_diff_UINS::_WIP_LONG_GetValue(LINK_LONG_TYPE* const _this, int32_t lIndex)
{
    return _this->prot().getsWipMemBlock().getalFree()[lIndex];
}



// ===========================================================================
/// Get the int32_t value entered in the box
// ===========================================================================
int32_t Ep2d_diff_UINS::_WIP_LONG_SetValue(LINK_LONG_TYPE* const _this, int32_t value, int32_t lIndex)
{
    return (_this->prot().getsWipMemBlock().getalFree()[lIndex] = value);
}



// ===========================================================================
/// Defines the valid range of the int32_t parameter
// ===========================================================================
bool Ep2d_diff_UINS::_WIP_LONG_GetLimits(LINK_LONG_TYPE* const /* _this */, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t lIndex)
{
    int32_t lMin, lMax, lInc;

    switch(lIndex)
    {
        // depending on the index of available int32_t boxes, define a range for valid values.
        case WIP_Anything:
        default:
            lMin = 1;
            lMax = 1;
            lInc = 1;
            break;
    }

    rulVerify = LINK_LONG_TYPE::VERIFY_BINARY_SEARCH;
    rLimitVector.resize(1);
    rLimitVector[0].setEqualSpaced(lMin, lMax, lInc);
    return true;
}




// ===========================================================================
/// Decide if parameter is available on card
// ===========================================================================
bool Ep2d_diff_UINS::_WIP_LONG_IsAvailable(LINK_LONG_TYPE* const /* _this */, int32_t /* lIndex */)
{
    return false;
}




// ===========================================================================
/// Searches for the correct tool tip text for int32_t parameters
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_LONG_GetToolTipId(LINK_LONG_TYPE* const /* _this */, char* arg_list[], int32_t lIndex)
{
    char tLine[100];
    char tToolTip[1000];

    tToolTip[0] = '\0';
    // MrProtocolData::MrProtData* pMrProt = &pThis->prot();
    switch(lIndex)
    {
        case WIP_Anything:
            sprintf(tLine, "Anytext\n");
            strcat(tToolTip, tLine);
            arg_list[0] = tToolTip;
            return MRI_STD_STRING;
            break;
        default:
            break;
    }

    return 0;
}



// -----------------------------------------------------------
// Implementation of UI functions for a 'double' WIP parameter                              
// -----------------------------------------------------------


// ===========================================================================
/// Defines the text in front of the box for double values
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_DOUBLE_GetLabelId(LINK_DOUBLE_TYPE* const /* _this */, char* arg_list[], int32_t lIndex)
{
    static const char* const pszLabel0 = "Anylabel";
    static const char* const pszLabelInvalid = "NotSupported";

    switch(lIndex) // depending on the Index of available boxes give them a name 
    {
        case WIP_Anything:
            arg_list[0] = (char*)pszLabel0;
            break;
        default:
            arg_list[0] = (char*)pszLabelInvalid;
            break;
    }
    return MRI_STD_STRING;
}



// ===========================================================================
/// Defines the text behind the box for double values
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_DOUBLE_GetUnitId(LINK_DOUBLE_TYPE* const /* _this */, char* arg_list[], int32_t lIndex)
{
    static const char* const pszLabel0 = "Anyunit";
    static const char* const pszLabelInvalid = "NotSupported";

    // depending on the Index of available int32_t boxes, set the correct text behind them
    switch(lIndex)
    {
        case WIP_Anything:
            arg_list[0] = (char*)pszLabel0;
            break;
        default:
            arg_list[0] = (char*)pszLabelInvalid;
            break;
    }
    return MRI_STD_STRING;
}



// ===========================================================================
/// Get the double value entered in the box
// ===========================================================================
double Ep2d_diff_UINS::_WIP_DOUBLE_GetValue(LINK_DOUBLE_TYPE* const _this, int32_t lIndex)
{
    return _this->prot().getsWipMemBlock().getadFree()[lIndex];
}


// ===========================================================================
/// Get the double value entered in the box
// ===========================================================================
double Ep2d_diff_UINS::_WIP_DOUBLE_SetValue(LINK_DOUBLE_TYPE* const _this, double value, int32_t lIndex)
{
    return (_this->prot().getsWipMemBlock().getadFree()[lIndex] = value);
}




// ===========================================================================
/// Defines the valid range of the double parameter
// ===========================================================================
bool Ep2d_diff_UINS::_WIP_DOUBLE_GetLimits(LINK_DOUBLE_TYPE* const /* _this */, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t lIndex)
{
    double dMin, dMax, dInc;

    switch(lIndex)
    {
        // depending on the index of available int32_t boxes, define a range for valid values.
        case WIP_Anything:
            dMin = -1.;
            dMax = 10.;
            dInc = 1.;
            break;
        default:
            dMin = 1.;
            dMax = 1.;
            dInc = 1.;
            break;
    }

    rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;
    rLimitVector.resize(1);
    rLimitVector[0].setEqualSpaced(dMin, dMax, dInc);
    //           ^ 1st Interval to set
    return true;
}



// ===========================================================================
/// Decide if parameter is available on card
// ===========================================================================
bool Ep2d_diff_UINS::_WIP_DOUBLE_IsAvailable(LINK_DOUBLE_TYPE* const /* _this */, int32_t lIndex)
{
    switch(lIndex)
    {
        case WIP_Anything:
            return true;
        default:
            return true;
    }
    return true;
}




// ===========================================================================
/// Searches for the correct tool tip text for double parameters 
// ===========================================================================
unsigned Ep2d_diff_UINS::_WIP_DOUBLE_GetToolTipId(LINK_DOUBLE_TYPE* const _this, char* arg_list[], int32_t lIndex)
{
    char tLine[256];
    char tToolTip[1000];

    tToolTip[0] = '\0';
    MrProtocolData::MrProtData* pMrProt = &_this->prot();

    switch(lIndex)
    {
        case WIP_Anything:
            sprintf(tLine, "Anytext\n");
            strcat(tToolTip, tLine);
            arg_list[0] = tToolTip;
            return MRI_STD_STRING;
            break;
        default: break;
    }

    return 0;
}

#endif     // of if WIP
#endif   // of #ifdef WIN32





///  \brief Constructor
///
Ep2d_diff_UI::Ep2d_diff_UI()
#ifdef WIN32
    : m_pDidi(NULL)
    , m_sImportExportError("")
    , m_bImportExportError(false)
    , m_bVectorFileImported(false)
    , m_iMemCount(0)
    , m_bThermalBalancing(false)
    , m_pFctTRGetLimits_orig(0)
    , m_pFctTRSetValue_orig(0)
    , m_pFctTESetValue_orig(0)
    , m_pFctTIGetLimits_orig(0)
    , m_pFctTISetValue_orig(0)
    , m_pFctSlicesSetValue_orig(0)
    , m_pFctConcSetValue_orig(0)
#endif
{
#ifdef WIN32
    m_dToolTipParam[ToolTipParamDiffGradDuration] = 0.;
    m_dToolTipParam[ToolTipParamDiffGradSpacing] = 0.;

    m_sToolTipString[ToolTipStringDVSInfo] = "";

#endif
    for(int i = 0; i < iMaxProtMemory; i++)
    {
        m_asMemory[i].bNeedOtherTETITR      = false;
        m_asMemory[i].bValidProtocol        = false;
        m_asMemory[i].lBValue               = -1;
        m_asMemory[i].lNeededTE             = -1;
        m_asMemory[i].lNeededTI             = -1;
        m_asMemory[i].lNeededTR             = -1;
    }

    m_isDiffusion = true;
    m_isSpinEcho  = true;
}


///  \brief Destructor
///
Ep2d_diff_UI::~Ep2d_diff_UI()
{
}


#ifdef WIN32

// ===========================================================================
// setting the thermal balancing flag
// ===========================================================================
void Ep2d_diff_UI::setThermalBalancing(bool bThermalBalancing)
{
    m_bThermalBalancing = bThermalBalancing;
}


// ===========================================================================
/// This function registers Didi to these UILink functions
/**
The UILink software relies on information about the available
diffusion directions. Therefore we must provide access from UILink
to Didi. This is accomplished by this global pointer address.
It requires that this fUILinkRegisterDidi function is called
during the initialization of the sequence.

The effect of this funtion is just to write the supplied pointer
address into the global variable (static_cast<Ep2d_diff_UI*>(EpCommonUINS::getUI(_this)))->m_pDidi.
*/
// ===========================================================================
void Ep2d_diff_UI::fUILinkRegisterDidi
(
DiffusionDirections *DidiAddress 	/*!< The address of the 'didi' actually used. A NULL pointer is
                                    also allowed to mark that no didi is actually registered.   */
                                    )
{
    m_pDidi = DidiAddress;
}

#ifdef SUPPORT_PACE
// ------------------------------------------------------------------------------
// Function    : fRespCompSetValue
// ------------------------------------------------------------------------------
//               
// Description : If PACE respiratory triggering is used the protocol parameter
//               MrProt::repetitions is not supported (i.e. must be 1) since
//               average loop must be inside of concatenation loop.
//               Therefore fRespCompSetValue multiplies the b-value specific averages
//               with the number of repetitions and set the number of repetitions
//               to one, if PACE is turned on.
// ------------------------------------------------------------------------------
unsigned Ep2d_diff_UINS::fRespCompSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    const UI_ELEMENT_SELECTION&  rRespComp = static_cast<Ep2d_diff_UI*>(getUI(pThis))->m_RespComp;

    if((rRespComp.getOrigSetValueHandler()) != 0)
    {
        uNewVal = (*(rRespComp.getOrigSetValueHandler()))(pThis, uNewVal, lPos);
    }
    return uNewVal;
}

bool Ep2d_diff_UINS::fUILinkBValueGetLimitsNew(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    MR_SEQUENCE_EXPORT_TLS(pThis->sequence());
    bool const bIsContextPrepForBinarySearch = pThis->seqLimits().isContextPrepForBinarySearch();
    ParLim<int32_t> const& rBValueLimits = pThis->seqLimits().getBValue();
    rulVerify = LINK_LONG_TYPE::VERIFY_BINARY_SEARCH;
    MrProt rMrProt(&pThis->prot());


    SEQ_NAMESPACE::Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis);

    if(!pEPI)
    {
        SEQ_TRACE_ERROR.print("fUILinkBValueGetLimitsNew ERROR: cannot get sequence pointer");
        return false;
    }
    // ================ IVIM support ==========================
    // IVIM works with monopolar diffusion mode and offers small b value increments up to 200s/mm2
    // the small b values are only allowed if the monopolar diffusion module does not use spoiler
    // gradients around the refocusing pulse


    // no IVIM
    if(rMrProt.diffusion().getdsScheme() != SEQ::DIFFSCHEME_MONOPOLAR)
    {
        rLimitVector.resize(1);
        return rLimitVector[0].setEqualSpaced(
            rBValueLimits.getMin(),
            rBValueLimits.getMax(),
            rBValueLimits.getInc()
            );
    }

    // IVIM
    rLimitVector.clear();

    // set 0 as first possible b value
    MrLimitLong rLonelyLimit;
    rLonelyLimit.setLonely(0);
    rLimitVector.push_back(rLonelyLimit);

    // set all values from smallest possible b value up to 200 with increment of 10
    MrLimitLong sVector;
    sVector.setEqualSpaced(pEPI->m_EPIKernel.getSmallestIVIMbValuePossible(bIsContextPrepForBinarySearch), pEPI->m_EPIKernel.getMaxBValueSmallIVIMIncrement(), pEPI->m_EPIKernel.getIVIMIncrement());
    rLimitVector.push_back(sVector);

    // set all values up to max value
    sVector.setEqualSpaced(pEPI->m_EPIKernel.getMaxBValueSmallIVIMIncrement() + rBValueLimits.getInc(), rBValueLimits.getMax(), rBValueLimits.getInc());
    rLimitVector.push_back(sVector);

    return true;

}

#endif //  SUPPORT_PACE

unsigned Ep2d_diff_UINS::fPhaseCorrSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    return (basicSetValue(pThis, uNewVal, lPos, getUI(pThis)->m_PhaseCorrMode.getOrigSetValueHandler()));
}


// ------------------------------------------------------------------------------
// IR Scheme
// ------------------------------------------------------------------------------
static unsigned fIRSchemeSolve(LINK_SELECTION_TYPE* const pThis, char *arg_list[], const void* pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lPos)
{
    MrProtocolData::MrProtData& rProt = pThis->prot();
    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_AUTO))
    {
        // First set IR scheme temporarily back to sequential.
        pThis->value(MRI_STD_SEQUENTIAL, lPos);
        if(rProt.getsSliceArray().getlConc() > 1)
        {
            if(rProt.getsGroupArray().getasGroup()[0].getdDistFact() < 0)
            {
                setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_SLICE_GROUP_LIST, MR_TAG_SG_DISTANCE_FACTOR, 0.0);
            }
            setProtocolParameter<LINK_LONG_TYPE>(pThis, MR_TAG_CONCATENATIONS, 1);
        }

        //  Next set TR to minimum
        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR, NULL, getUI(pThis)->m_lNeededTR / 1000.0);

        //  Back to Automatic
        pThis->value(MRI_STD_AUTO, lPos);

        if(pThis->tryProt(const_cast<void*>(pVoid), pOrigProt, lPos))
        {
            //  Confirmation text is formatted automatically ...
            return MRI_STD_CONFIRMATION_MSG;
        }
    }
    return 0;
}

//  Calculates 


unsigned fIRSchemeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos)
{
    MrProtocolData::MrProtData& rProt = pThis->prot();
    if(uNewVal == MRI_STD_SEQUENTIAL)
    {
        if(pThis->sequence().prepareForBinarySearch(&rProt))
        {
            if(SEQ_NAMESPACE::Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis))
            {
                const int iSBBScanTime_us = int(pEPI->m_mySeqLoop.getlSBBScanTime())
                    , iKernelTime_us = int(pEPI->m_mySeqLoop.getlKernelScanTime())
                    , iTI_us = rProt.getalTI()[0]
                    , iNSlices = rProt.getsSliceArray().getlSize()
                    , iNConc = std::max(int(1), int(rProt.getsSliceArray().getlConc()))
                    , iNSlcPerConc = (iNSlices + iNConc - 1) / iNConc
                    , iIRTime_us = int(const_cast<SeqBuildBlockIRsel*>(&pEPI->m_mySeqLoop.getSBBIRsel())->getDurationPerRequest())
                    , iCoolPauseWithinKernelTime_us = std::max(0, int(pEPI->m_lCoolPauseTotal - pEPI->m_lCoolPauseImplicit))
                    , iTRIncr_prefered_us = 10000
                    ;
                int iKernelOffset = 2;

                for(;iKernelOffset >= 0;--iKernelOffset)
                {
                    int iTBlock_us = 0;
                    if(fCalcTBlock(iTBlock_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                    {
                        //  Valid solution found
                        int iTR_us = iNSlcPerConc*iTBlock_us;
                        int iTRIncr_try_us = iTRIncr_prefered_us;
                        for(;iTRIncr_try_us >= GRAD_RASTER_TIME;iTRIncr_try_us /= 10)
                        {
                            //  i) Find the least common multiple of iTRIncr_try_us and iNSlcPerConc
                            const int iLCM_us = MODULE::LCM(iNSlcPerConc, iTRIncr_try_us)
                                , iTBlock_prefered_us = MODULE::IMULT_CEIL(iTBlock_us, iLCM_us / iNSlcPerConc);

                            if(fTestTBlock(iTBlock_prefered_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                            {
                                iTBlock_us = iTBlock_prefered_us;
                                iTR_us = iNSlcPerConc*iTBlock_prefered_us;
                                break;
                            }
                        }
                        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR, NULL, iTR_us / 1000.);
                        break;
                    }
                }
            }
        }
        rProt.getsPrepPulses().setucIRScheme(SEQ::IR_SCHEME_SEQUENTIAL);
    }
    else
    {
        rProt.getsPrepPulses().setucIRScheme(SEQ::IR_SCHEME_AUTO);
    }
    return pThis->value(lPos);
}

// ------------------------------------------------------------------------------
// # slices/concatenations
// ------------------------------------------------------------------------------

int32_t fSlicesSetValue(LINK_LONG_TYPE* const pThis, int32_t lNSlices_new, int32_t lPos)
{
    MrProtocolData::MrProtData& rProt = pThis->prot();

    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        //  Changing the # slices does not change TBlock
        int iNSlices_old = rProt.getsSliceArray().getlSize()
            , iNConc = std::max(int(1), int(rProt.getsSliceArray().getlConc()))
            , iNSlcPerConc_old = (iNSlices_old + iNConc - 1) / iNConc
            , iNSlcPerConc_new = (int(lNSlices_new) + iNConc - 1) / iNConc
            , iTR_us_old = rProt.getalTR()[0]
            , iTBloc_us = iTR_us_old / iNSlcPerConc_old
            , iTR_us_new = iTBloc_us*iNSlcPerConc_new
            ;
        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR, NULL, iTR_us_new / 1000.);
    }
    if(Ep2d_diff_UI* pDiffUI = dynamic_cast<Ep2d_diff_UI*>(getUI(pThis)))
    {
        if(pDiffUI->m_pFctSlicesSetValue_orig != 0)
        {
            return (*(pDiffUI->m_pFctSlicesSetValue_orig))(pThis, lNSlices_new, lPos);
        }
    }
    return pThis->value(lPos);
}

int32_t fConcSetValue(LINK_LONG_TYPE* const pThis, int32_t lNConc_new, int32_t lPos)
{
    MrProtocolData::MrProtData& rProt = pThis->prot();

    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        //  Changing the # slices does not change TBlock
        int iNSlices = rProt.getsSliceArray().getlSize()
            , iNConc_old = std::max(int(1), int(rProt.getsSliceArray().getlConc()))
            , iNSlcPerConc_old = (iNSlices + iNConc_old - 1) / iNConc_old
            , iNConc_new = std::max(int(1), int(lNConc_new))
            , iNSlcPerConc_new = (int(iNSlices) + iNConc_new - 1) / iNConc_new
            , iTR_us_old = rProt.getalTR()[0]
            , iTBloc_us = iTR_us_old / iNSlcPerConc_old
            , iTR_us_new = iTBloc_us*iNSlcPerConc_new
            ;
        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR, NULL, iTR_us_new / 1000.);
    }
    if(Ep2d_diff_UI* pDiffUI = dynamic_cast<Ep2d_diff_UI*>(getUI(pThis)))
    {
        if(pDiffUI->m_pFctConcSetValue_orig != 0)
        {
            return (*(pDiffUI->m_pFctConcSetValue_orig))(pThis, lNConc_new, lPos);
        }
    }
    return pThis->value(lPos);
}

unsigned fConcSolve(LINK_LONG_TYPE* const pThis, char* arg_list[], const void* pVoid, const MrProtocolData::MrProtData*  pOrig, int32_t pos)
{
    MrProt rMrProt(&pThis->prot());
    Ep2d_diff_UI* pDiffUI = dynamic_cast<Ep2d_diff_UI*>(getUI(pThis));

    if(!pDiffUI)
        return 0;

    unsigned uResult = (*(pDiffUI->m_pFctConcSolve_orig))(pThis, arg_list, pVoid, pOrig, pos);

    if(uResult == MRI_STD_CONFIRMATION_MSG)
        return uResult;

    if (SEQ::SIGNAL_RESPIRATION != pThis->prot().getsPhysioImaging().getlSignal1())
    {
        // set acquisition window to TR so that the # of concats can be decreased
        setProtocolParameter<LINK_DOUBLE_TYPE>(pThis, MR_TAG_ACQUISITION_WINDOW_PACE, rMrProt.getalTR()[0] / 1000.0);

        if (pThis->sequence().prepareForBinarySearch(rMrProt))
            return MRI_STD_CONFIRMATION_MSG;
    }

    return 0;
}

// ------------------------------------------------------------------------------
// TE
// ------------------------------------------------------------------------------

double fTESetValue(LINK_DOUBLE_TYPE* const pThis, double dNewVal_ms, int32_t lPos)
{
    if(Ep2d_diff_UI* pDiffUI = dynamic_cast<Ep2d_diff_UI*>(getUI(pThis)))
    {
        if(pDiffUI->m_pFctTESetValue_orig != 0)
        {
            dNewVal_ms = (*(pDiffUI->m_pFctTESetValue_orig))(pThis, dNewVal_ms, lPos);
        }
    }
    MrProtocolData::MrProtData& rProt = pThis->prot();

    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis);
        if(pEPI != 0)
        {
            //  Used as indicator that the preparation of the kernel succeeded
            pEPI->m_mySeqLoop.setlSBBScanTime(0);
            pEPI->m_mySeqLoop.setlKernelScanTime(0);
            pEPI->m_mySeqLoop.setCoolPauseWithinKernelTime_us(0);
        }
        pThis->sequence().prepareForBinarySearch(&rProt);
        if(pEPI != 0)
        {
            const int   iSBBScanTime_us = int(pEPI->m_mySeqLoop.getlSBBScanTime())
                , iKernelTime_us = int(pEPI->m_mySeqLoop.getlKernelScanTime())
                , iTI_us = rProt.getalTI()[0]
                , iNSlices = rProt.getsSliceArray().getlSize()
                , iNConc = std::max(int(1), int(rProt.getsSliceArray().getlConc()))
                , iNSlcPerConc = (iNSlices + iNConc - 1) / iNConc
                , iIRTime_us = int(const_cast<SeqBuildBlockIRsel*>(&pEPI->m_mySeqLoop.getSBBIRsel())->getDurationPerRequest())
                , iCoolPauseWithinKernelTime_us = std::max(0, int(pEPI->m_lCoolPauseTotal - pEPI->m_lCoolPauseImplicit))
                , iTRIncr_prefered_us = 10000
                ;
            if(iKernelTime_us > 0)
            {
                int iKernelOffset = rProt.getsPrepPulses().getlKernelOffset();
                int iTBlock_us = 0;
                if(fCalcTBlock(iTBlock_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                {
                    //  Valid solution found
                    int iTR_us = iNSlcPerConc*iTBlock_us;
                    const int iProtTR_us = rProt.getalTR()[0];
                    int iTRIncr_try_us = iTRIncr_prefered_us;
                    for(;iTRIncr_try_us >= GRAD_RASTER_TIME;iTRIncr_try_us /= 10)
                    {
                        //  i) Find the least common multiple of iTRIncr_try_us and iNSlcPerConc
                        const int iLCM_us = MODULE::LCM(iNSlcPerConc, iTRIncr_try_us)
                            , iTBlock_prefered_us = MODULE::IMULT_CEIL(iTBlock_us, iLCM_us / iNSlcPerConc)
                            ;
                        if(fTestTBlock(iTBlock_prefered_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                        {
                            iTBlock_us = iTBlock_prefered_us;
                            iTR_us = iNSlcPerConc*iTBlock_prefered_us;
                            break;
                        }
                    }
                    if(iTR_us != iProtTR_us)
                    {
                        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR, NULL, iTR_us / 1000.0);
                    }
                }
            }
        }
    }
    return pThis->value(lPos);
}

unsigned fTESolve(LINK_DOUBLE_TYPE* const pThis, char* arg_list[], const void* pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lPos)
{
    MrProtocolData::MrProtData& rProt = pThis->prot();

    // try common EPI solve hander first
    unsigned uResult = EpCommonUINS::_solveTETITR_TE(pThis, arg_list, pVoid, pOrigProt, lPos);
    if(uResult == MRI_STD_CONFIRMATION_MSG)
        return uResult;

    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        const int iInvTE_us = rProt.getalTE()[0], iValTE_us = pOrigProt->getalTE()[0];
        int iKernelOffset = rProt.getsPrepPulses().getlKernelOffset();
        if((iInvTE_us != iValTE_us) && (iKernelOffset > 0))
        {
            // A larger kernel Offset may work due to the larger implicit cooling pause
            if(Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis))
            {
                const int iSBBScanTime_us = int(pEPI->m_mySeqLoop.getlSBBScanTime())
                    , iKernelTime_us = int(pEPI->m_mySeqLoop.getlKernelScanTime())
                    , iTI_us = rProt.getalTI()[0]
                    , iNSlices = rProt.getsSliceArray().getlSize()
                    , iNConc = std::max(int(1), int(rProt.getsSliceArray().getlConc()))
                    , iNSlcPerConc = (iNSlices + iNConc - 1) / iNConc
                    , iIRTime_us = int(const_cast<SeqBuildBlockIRsel*>(&pEPI->m_mySeqLoop.getSBBIRsel())->getDurationPerRequest())
                    , iCoolPauseWithinKernelTime_us = std::max(0, int(pEPI->m_lCoolPauseTotal - pEPI->m_lCoolPauseImplicit))
                    , iTRIncr_prefered_us = 10000
                    ;
                if(iKernelTime_us > 0)
                {
                    for(;--iKernelOffset >= 0;)
                    {
                        int iTBlock_us = 0;
                        if(fCalcTBlock(iTBlock_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                        {
                            //  Valid solution found
                            int iTR_us = iNSlcPerConc*iTBlock_us;
                            const int iProtTR_us = rProt.getalTR()[0];
                            int iTRIncr_try_us = iTRIncr_prefered_us;
                            for(;iTRIncr_try_us >= GRAD_RASTER_TIME;iTRIncr_try_us /= 10)
                            {
                                //  i) Find the least common multiple of iTRIncr_try_us and iNSlcPerConc
                                const int iLCM_us = MODULE::LCM(iNSlcPerConc, iTRIncr_try_us)
                                    , iTBlock_prefered_us = MODULE::IMULT_CEIL(iTBlock_us, iLCM_us / iNSlcPerConc);

                                if(fTestTBlock(iTBlock_prefered_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                                {
                                    iTBlock_us = iTBlock_prefered_us;
                                    iTR_us = iNSlcPerConc*iTBlock_prefered_us;
                                    break;
                                }
                            }
                            if(iTR_us != iProtTR_us)
                            {
                                // TR is NOT added as dependent parameter as the limits are all green. The limits are not determined by
                                // the usual search but in the function "fTIGetLimits" which does not have any solving mechanism. See CHARM 469623.
                                setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR, NULL, iTR_us / 1000., 0, false);

                                if(pThis->sequence().prepareForBinarySearch(&rProt))
                                {
                                    return MRI_STD_CONFIRMATION_MSG;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

// ------------------------------------------------------------------------------
// TI
// ------------------------------------------------------------------------------

bool fTIGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t lPos)
{
    MrProtocolData::MrProtData& rProt = pThis->prot();
    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        if(pThis->sequence().prepareForBinarySearch(&rProt))
        {
            if(SEQ_NAMESPACE::Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis))
            {

                int   iSBBScanTime_us = int(pEPI->m_mySeqLoop.getlSBBScanTime())
                    , iKernelTime_us = int(pEPI->m_mySeqLoop.getlKernelScanTime())
                    , iIRTime_us = int(const_cast<SeqBuildBlockIRsel*>(&pEPI->m_mySeqLoop.getSBBIRsel())->getDurationPerRequest())
                    , iCoolPauseWithinKernelTime_us = std::max(0, int(pEPI->m_lCoolPauseTotal - pEPI->m_lCoolPauseImplicit))
                    , iKernelOffset = rProt.getsPrepPulses().getlKernelOffset()
                    ;

                const ParLim<int32_t>& rSeqLimits = pThis->seqLimits().getTI()[lPos];
                const int iIncTI_us = std::max(int(rSeqLimits.getInc()), 5000)
                    , iMinTI_us = std::max
                    (int(rSeqLimits.getMin())
                    , iIncTI_us*((iIRTime_us+iSBBScanTime_us+iIncTI_us-1)/iIncTI_us)
                    )
                    , iMaxTI_us = iIncTI_us*(int(rSeqLimits.getMax()) / iIncTI_us)
                    ;
                rLimitVector.clear();
                rLimitVector.reserve(1 + (iMaxTI_us - iMinTI_us) / iIncTI_us);

                int iTI_us = iMinTI_us;
                for(;iTI_us <= iMaxTI_us;iTI_us += iIncTI_us)
                {
                    int iTBlock_us = 0;
                    if(fCalcTBlock(iTBlock_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                    {
                        //  Valid solution found
                        MrLimitDouble sLim;
                        sLim.setLonely(iTI_us / 1000.);
                        rLimitVector.push_back(sLim);
                    }
                } //  KernelOffsetLoop
                rulVerify = LINK_DOUBLE_TYPE::VERIFY_OFF;
                return rLimitVector.size() > 0;
            }
        }
    }
    if(Ep2d_diff_UI* pDiffUI = dynamic_cast<Ep2d_diff_UI*>(getUI(pThis)))
    {
        if(pDiffUI->m_pFctTIGetLimits_orig != 0)
        {
            return (*(pDiffUI->m_pFctTIGetLimits_orig))(pThis, rLimitVector, rulVerify, lPos);
        }
    }
    return false;
}

double fTISetValue(LINK_DOUBLE_TYPE* const pThis, double dNewVal_ms, int32_t lPos)
{
    if(Ep2d_diff_UI* pDiffUI = dynamic_cast<Ep2d_diff_UI*>(getUI(pThis)))
    {
        if(pDiffUI->m_pFctTISetValue_orig != 0)
        {
            dNewVal_ms = (*(pDiffUI->m_pFctTISetValue_orig))(pThis, dNewVal_ms, lPos);
        }
    }
    MrProtocolData::MrProtData& rProt = pThis->prot();

    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis);
        if(pEPI != 0)
        {
            //  Used as indicator that the preparation of the kernel succeeded
            pEPI->m_mySeqLoop.setlSBBScanTime(0);
            pEPI->m_mySeqLoop.setlKernelScanTime(0);
            pEPI->m_mySeqLoop.setCoolPauseWithinKernelTime_us(0);
        }
        pThis->sequence().prepareForBinarySearch(&rProt);
        if(pEPI != 0)
        {
            const int iSBBScanTime_us = int(pEPI->m_mySeqLoop.getlSBBScanTime())
                , iKernelTime_us = int(pEPI->m_mySeqLoop.getlKernelScanTime())
                , iTI_us = rProt.getalTI()[0]
                , iNSlices = rProt.getsSliceArray().getlSize()
                , iNConc = std::max(int(1), int(rProt.getsSliceArray().getlConc()))
                , iNSlcPerConc = (iNSlices + iNConc - 1) / iNConc
                , iIRTime_us = int(const_cast<SeqBuildBlockIRsel*>(&pEPI->m_mySeqLoop.getSBBIRsel())->getDurationPerRequest())
                , iCoolPauseWithinKernelTime_us = std::max(0, int(pEPI->m_lCoolPauseTotal - pEPI->m_lCoolPauseImplicit))
                , iTRIncr_prefered_us = 10000
                ;
            if(iKernelTime_us > 0)
            {
                int iKernelOffset = rProt.getsPrepPulses().getlKernelOffset();
                int iTBlock_us = 0;
                if(fCalcTBlock(iTBlock_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                {
                    //  Valid solution found
                    int iTR_us = iNSlcPerConc*iTBlock_us;
                    const int iProtTR_us = rProt.getalTR()[0];
                    int iTRIncr_try_us = iTRIncr_prefered_us;
                    for(;iTRIncr_try_us >= GRAD_RASTER_TIME;iTRIncr_try_us /= 10)
                    {
                        //  i) Find the least common multiple of iTRIncr_try_us and iNSlcPerConc
                        const int iLCM_us = MODULE::LCM(iNSlcPerConc, iTRIncr_try_us)
                            , iTBlock_prefered_us = MODULE::IMULT_CEIL(iTBlock_us, iLCM_us / iNSlcPerConc);

                        if(fTestTBlock(iTBlock_prefered_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                        {
                            iTBlock_us = iTBlock_prefered_us;
                            iTR_us = iNSlcPerConc*iTBlock_prefered_us;
                            break;
                        }
                    }
                    if(iTR_us != iProtTR_us)
                    {
                        setProtocolParameterElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR, NULL, (double)(iTR_us / 1000.f));
                    }
                }
            }
        }
    }
    return pThis->value(lPos);
}
// ------------------------------------------------------------------------------
// TR
// ------------------------------------------------------------------------------

double fTRSetValue(LINK_DOUBLE_TYPE* const pThis, double dNewVal_ms, int32_t lPos)
{
    if(Ep2d_diff_UI* pDiffUI = dynamic_cast<Ep2d_diff_UI*>(getUI(pThis)))
    {
        if(pDiffUI->m_pFctTRSetValue_orig != 0)
        {
            dNewVal_ms = (*(pDiffUI->m_pFctTRSetValue_orig))(pThis, dNewVal_ms, lPos);
        }
    }
    MrProtocolData::MrProtData& rProt = pThis->prot();

    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        if(Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis))
        {
            if(pThis->sequence().prepareForBinarySearch(&rProt))
            {
                int iKernelOffset = pEPI->m_mySeqLoop.getKernelOffset();
                rProt.getsPrepPulses().setlKernelOffset(iKernelOffset);
            }
        }
    }
    return pThis->value(lPos);
}

static unsigned fTRGetPrec(LINK_DOUBLE_TYPE* const, int32_t)
{
    return 0;
}
bool fTRGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t lPos)
{
    MrProtocolData::MrProtData& rProt = pThis->prot();
    if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
    {
        if(pThis->sequence().prepareForBinarySearch(&rProt))
        {
            if(SEQ_NAMESPACE::Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis))
            {

                int   iSBBScanTime_us = int(pEPI->m_mySeqLoop.getlSBBScanTime())
                    , iKernelTime_us = int(pEPI->m_mySeqLoop.getlKernelScanTime())
                    , iTI_us = rProt.getalTI()[0]
                    , iNSlices = rProt.getsSliceArray().getlSize()
                    , iNConc = std::max(int(1), int(rProt.getsSliceArray().getlConc()))
                    , iNSlcPerConc = (iNSlices + iNConc - 1) / iNConc
                    , iIRTime_us = int(const_cast<SeqBuildBlockIRsel*>(&pEPI->m_mySeqLoop.getSBBIRsel())->getDurationPerRequest())
                    , iTRIncr_us = 10000
                    , iCoolPauseWithinKernelTime_us = std::max(0, int(pEPI->m_lCoolPauseTotal - pEPI->m_lCoolPauseImplicit))
                    ;
                int iKernelOffset = 5;
                rLimitVector.clear();
                rLimitVector.reserve(iKernelOffset + 1);

                //  First past: Check whether a particular Kernel offset can be realized at all
                //  std::pair<int,int>: first is kernel offset, second is minimum block duration in us
                std::vector< std::pair<int, int> > asFirstPass;
                asFirstPass.reserve(iKernelOffset + 1);

                //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+
                //  |IR|  |S|Kernel_n         |  |IR|  |S|Kernel_n+1       |  |IR|  |S|Kernel_n+2       |
                //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  ...
                //  <--------- TBlock ----------->
                //                               <--------- TBlock ----------->     <-> SBBSccanTime
                //  <-------> TI for KernelOffset==0 

                //  IR  ...  inversion recovery pulse
                //  S   ...  Saturation bands, chemical sat, spoilers before kernel, ...

                for(;iKernelOffset >= 0;--iKernelOffset)
                {
                    int iTBlock_us = 0;
                    if(fCalcTBlock(iTBlock_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                    {
                        //  Valid solution found
                        std::pair<int, int> sValid(iKernelOffset, iTBlock_us);
                        //  However, we prefer a TR which is on the 10 ms raster
                        int iTRIncr_try_us = iTRIncr_us;
                        for(;iTRIncr_try_us >= GRAD_RASTER_TIME;iTRIncr_try_us /= 10)
                        {
                            //  i) Find the least common multiple of iTRIncr_try_us and iNSlcPerConc
                            const int iLCM_us = MODULE::LCM(iNSlcPerConc, iTRIncr_try_us)
                                , iTBlock_prefered_us = MODULE::IMULT_CEIL(iTBlock_us, iLCM_us / iNSlcPerConc);

                            if(fTestTBlock(iTBlock_prefered_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                            {
                                sValid.second = iTBlock_prefered_us;
                                break;
                            }
                        }
                        asFirstPass.push_back(sValid);
                    }
                } //  KernelOffsetLoop
                rulVerify = LINK_DOUBLE_TYPE::VERIFY_OFF;
                const int iDeltaTR_us = MODULE::LCM(iNSlcPerConc, iTRIncr_us)
                    , iDeltaTBlock_us = iDeltaTR_us / iNSlcPerConc
                    ;
                const ParLim<int32_t>& rSeqLimits = pThis->seqLimits().getTR()[0];
                std::vector< std::pair<int, int> >::iterator it = asFirstPass.begin(), itEnd = asFirstPass.end();
                for(;it != itEnd;++it)
                {
                    MrLimitDouble sLim;
                    if((*it).first == 0)
                    {
                        // No interleaving. Hence we should be able to realize any TR
                        const int iTRMax_us = iNSlcPerConc*(*it).second + MODULE::IMULT_FLOOR(rSeqLimits.getMax() - iNSlcPerConc*(*it).second, iDeltaTR_us);
                        if((iTRMax_us - iNSlcPerConc*(*it).second) > iDeltaTR_us)
                        {
                            sLim.setEqualSpaced((iNSlcPerConc*(*it).second) / 1000., iTRMax_us / 1000., iDeltaTR_us / 1000.);
                        }
                        else
                        {
                            sLim.setLonely((iNSlcPerConc*(*it).second) / 1000.);
                        }
                        rLimitVector.push_back(sLim);
                    }
                    else
                    {
                        sLim.setLonely((iNSlcPerConc*(*it).second) / 1000.);
                        rLimitVector.push_back(sLim);
                        std::vector< std::pair<int, int> >::iterator itNext = it;
                        if(++itNext == itEnd)
                        {
                            continue;
                        }
                        //  Test additional block times 
                        int iTBlock4Test_us = (*it).second + iDeltaTBlock_us;
                        for(;((*itNext).second - iTBlock4Test_us) > iDeltaTBlock_us;)
                        {
                            if(fTestTBlock(iTBlock4Test_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us - iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, (*it).first))
                            {
                                sLim.setLonely((iNSlcPerConc*iTBlock4Test_us) / 1000.);
                                rLimitVector.push_back(sLim);
                                iTBlock4Test_us += iDeltaTBlock_us;
                            }
                            else
                            {
                                iTBlock4Test_us += GRAD_RASTER_TIME;
                            }
                        }
                    }
                }

                return rLimitVector.size() > 0;
            } // SEQ_NAMESPACE::getSeq(pThis)
        }  //  pThis->sequence().prepareForBinarySearch(&rProt)
        return false;
    }  //  IR_SCHEME_SEQUENTIAL
    
    if (SEQ::SIGNAL_RESPIRATION == pThis->prot().getsPhysioImaging().getlSignal1())
    {
        rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;
        rLimitVector.clear();

        const ParLim<int32_t>& _seqLimits = pThis->seqLimits().getTR()[lPos];

        double maxTR_ms = _seqLimits.getMax() / 1000.0;

        return _timeLimits(rLimitVector, _seqLimits.getMin() / 1000.0, maxTR_ms, _seqLimits.getInc() / 1000.0) != 0;
    }
    else
    {
        if (Ep2d_diff_UI* pDiffUI = dynamic_cast<Ep2d_diff_UI*>(getUI(pThis)))
        {
            if (pDiffUI->m_pFctTRGetLimits_orig != 0)
            {
                return (*(pDiffUI->m_pFctTRGetLimits_orig))(pThis, rLimitVector, rulVerify, lPos);
            }
        }
    }
    return false;
}

#endif //  WIN32


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

NLS_STATUS Ep2d_diff_UI::registerUI(SeqLim &rSeqLim)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    // -------------------------------------
    // Call common EPI registerUI method
    // -------------------------------------
    if((NLS_SEV & (lStatus = EpCommonUI::registerUI(rSeqLim))) == NLS_SEV)
    {
        SEQ_TRACE_ERROR.print("EpCommonUI::registerUI failed.");
        return lStatus;
    }

    // -------------------------------------
    // Register variant specific UI handlers
    // -------------------------------------

    // Note: handlers registered in EpCommonUI potentially get replaced!
    //
    // If they still are supposed to be called AFTER the variant
    // specific handler, one has to explicitly take care of
    // this (e.g. store the old handler and call it at the
    // end of the new handler).

#ifdef WIN32
    //-------------------------------------------------------------------------------------
    // Register solve-handler IR-Pulse, TI getLimit- and solve-handler
    //-------------------------------------------------------------------------------------
    if(!fEPIStdRegisterTIHandlers(rSeqLim, &EpCommonUINS::pfGetTRTIFillTimes))
    {
        SEQ_TRACE_ERROR.print("fEPIStdRegisterTIHandlers failed");
        return MRI_SEQ_SEQU_ERROR;
    }

    //-----------------------------------------------------------------------------------
    // Register diffusion directions handlers
    //-----------------------------------------------------------------------------------
    m_NumberDiffDirs.registerSolveHandler(rSeqLim, MR_TAG_NUMBER_DIFF_DIRS, Ep2d_diff_UINS::fUILinkNoDiffDirectionsSolve);
    m_NumberDiffDirs.registerGetLimitsHandler(rSeqLim, MR_TAG_NUMBER_DIFF_DIRS, Ep2d_diff_UINS::fUILinkNoDiffDirectionsGetLimits);
    m_NumberDiffDirs.registerToolTipHandler(rSeqLim, MR_TAG_NUMBER_DIFF_DIRS, Ep2d_diff_UINS::fUILinkDiffDirectionsGetToolTipId);

    //-----------------------------------------------------------------------------------
    // Register diffusion directions import / export handlers
    //-----------------------------------------------------------------------------------
    m_DiffDirsImport.registerGetValueHandler(rSeqLim, MR_TAG_DIFF_DIRS_IMPORT, Ep2d_diff_UINS::fUILinkDiffDirsImportGetValue);
    m_DiffDirsImport.registerSetValueHandler(rSeqLim, MR_TAG_DIFF_DIRS_IMPORT, Ep2d_diff_UINS::fUILinkDiffDirsImportSetValue);
    m_DiffDirsImport.registerIsAvailableHandler(rSeqLim, MR_TAG_DIFF_DIRS_IMPORT, Ep2d_diff_UINS::fUILinkDiffDirsImportIsAvailable);
    m_DiffDirsImport.registerSolveHandler(rSeqLim, MR_TAG_DIFF_DIRS_IMPORT, Ep2d_diff_UINS::fUILinkDiffDirsImportSolve);
    m_DiffDirsExport.registerGetValueHandler(rSeqLim, MR_TAG_DIFF_DIRS_EXPORT, Ep2d_diff_UINS::fUILinkDiffDirsExportGetValue);
    m_DiffDirsExport.registerSetValueHandler(rSeqLim, MR_TAG_DIFF_DIRS_EXPORT, Ep2d_diff_UINS::fUILinkDiffDirsExportSetValue);
    m_DiffDirsExport.registerIsAvailableHandler(rSeqLim, MR_TAG_DIFF_DIRS_EXPORT, Ep2d_diff_UINS::fUILinkDiffDirsExportIsAvailable);
    m_DiffDirsExport.registerSolveHandler(rSeqLim, MR_TAG_DIFF_DIRS_EXPORT, Ep2d_diff_UINS::fUILinkDiffDirsImportSolve);    // Same solve handler for import and export 

    //----------------------------------------------------------------------------------------------
    // Register dedicated handlers for setting free diffusion directions during protocol conversion
    //----------------------------------------------------------------------------------------------
    if(LINK_GENERIC_TYPE* pFreeDiffData = _create<LINK_GENERIC_TYPE>(rSeqLim, MR_TAG_FREE_DIFF_DATA))
    {
        pFreeDiffData->registerGetValueHandler(Ep2d_diff_UINS::fUIFreeDiffDataGetValueGeneric);
        pFreeDiffData->registerSetValueHandler(Ep2d_diff_UINS::fUIFreeDiffDataSetValueGeneric);
    }

    //-----------------------------------------------------------------------------------
    // Register diffusion mode handlers
    //-----------------------------------------------------------------------------------
    m_DiffMode.registerSetValueHandler(rSeqLim, MR_TAG_DIFF_MODE, Ep2d_diff_UINS::fUILinkDiffusionModeSetValueNew);
    m_DiffMode.registerSolveHandler(rSeqLim, MR_TAG_DIFF_MODE, Ep2d_diff_UINS::fEPIDiffSolveSelectionConflict);
    m_DiffMode.registerGetOptionsHandler(rSeqLim, MR_TAG_DIFF_MODE, Ep2d_diff_UINS::fUILinkDiffusionModeGetOptions);

    //-----------------------------------------------------------------------------------
    // Register diffusion scheme handlers
    //-----------------------------------------------------------------------------------
    //m_DiffScheme.registerSolveHandler           ( rSeqLim, MR_TAG_DIFF_SCHEME,      EpCommonUINS::fEPIStdSolveSelection                    );
    m_DiffScheme.registerSolveHandler(rSeqLim, MR_TAG_DIFF_SCHEME, Ep2d_diff_UINS::fUILinkSolveDiffusionScheme);
    m_DiffScheme.registerToolTipHandler(rSeqLim, MR_TAG_DIFF_SCHEME, Ep2d_diff_UINS::fUILinkDiffusionSchemeGetToolTipId);

    //-----------------------------------------------------------------------------------
    // Register pTX acceleration handlers
    //-----------------------------------------------------------------------------------
    m_PTXPulseTxAccDiff.registerElmSolveHandler(rSeqLim, MR_TAG_PTX_PULSES_ARRAY, MR_TAG_PTX_TX_ACC, Ep2d_diff_UINS::fUILinkSolvePTXAcceleration);


    //-----------------------------------------------------------------------------------
    // Register q-space related handlers
    //-----------------------------------------------------------------------------------
    m_QSpaceSteps.registerGetLimitsHandler(rSeqLim, MR_TAG_DIFF_QSPACE_STEPS, Ep2d_diff_UINS::fUILinkNoQSpaceStepsGetLimits);
    m_QSpaceSteps.registerSolveHandler(rSeqLim, MR_TAG_DIFF_QSPACE_STEPS, EpCommonUINS::fEPIStdSolveLongParam);
    m_QSpaceMaxBValue.registerSolveHandler(rSeqLim, MR_TAG_DIFF_QSPACE_MAX_B_VALUE, EpCommonUINS::fEPIStdSolveLongParam);
    m_QSpaceCoverage.registerSolveHandler(rSeqLim, MR_TAG_DIFF_QSPACE_COVERAGE, EpCommonUINS::fEPIStdSolveSelection);
    m_QSpaceSampling.registerSolveHandler(rSeqLim, MR_TAG_DIFF_QSPACE_SAMPLING_SCHEME, EpCommonUINS::fEPIStdSolveSelection);

    //-----------------------------------------------------------------------------------
    // Register dynamic distortion correction filter handlers
    //-----------------------------------------------------------------------------------
    m_FilterDynDistCorr.registerGetOptionsHandler(rSeqLim, MR_TAG_FLT_DYN_DISCOR, Ep2d_diff_UINS::fEPIDiffFltDynDistCorrGetOptions);

    //-----------------------------------------------------------------------------------
    // Register multiple series handlers
    //-----------------------------------------------------------------------------------
    m_MultipleSeries.registerIsAvailableHandler(rSeqLim, MR_TAG_MULTIPLE_SERIES, Ep2d_diff_UINS::fUILink_MultipleSeriesMode_IsAvailable);

    //-----------------------------------------------------------------------------------
    // Register inversion handlers
    //-----------------------------------------------------------------------------------
    // Note: this overrides the base class handler!
    m_Inversion.registerSolveHandler(rSeqLim, MR_TAG_INVERSION, Ep2d_diff_UINS::fEPIDiffSolveSelectionConflict);

    if(LINK_SELECTION_TYPE* pSel = _search<LINK_SELECTION_TYPE>(rSeqLim, MR_TAG_IR_SCHEME))
    {
        pSel->registerSolveHandler(fIRSchemeSolve);
        pSel->registerSetValueHandler(fIRSchemeSetValue);
    }
    if(LINK_DOUBLE_TYPE* pTR = _searchElm<LINK_DOUBLE_TYPE>(rSeqLim, MR_TAG_TR))
    {
        this->m_pFctTRGetLimits_orig = pTR->registerGetLimitsHandler(fTRGetLimits);
        this->m_pFctTRSetValue_orig = pTR->registerSetValueHandler(fTRSetValue);
    }
    if(LINK_DOUBLE_TYPE* pTE = _searchElm<LINK_DOUBLE_TYPE>(rSeqLim, MR_TAG_TE))
    {
        this->m_pFctTESetValue_orig = pTE->registerSetValueHandler(fTESetValue);
        pTE->registerSolveHandler(fTESolve);
    }
    if(LINK_DOUBLE_TYPE* pTI = _searchElm<LINK_DOUBLE_TYPE>(rSeqLim, MR_TAG_TI))
    {
        this->m_pFctTIGetLimits_orig = pTI->registerGetLimitsHandler(fTIGetLimits);
        this->m_pFctTISetValue_orig = pTI->registerSetValueHandler(fTISetValue);
    }
    if(LINK_LONG_TYPE* pLng = _searchElm<LINK_LONG_TYPE>(rSeqLim, MR_TAG_SLICE_GROUP_LIST, MR_TAG_SG_SIZE))
    {
        this->m_pFctSlicesSetValue_orig = pLng->registerSetValueHandler(fSlicesSetValue);
    }
    if(LINK_LONG_TYPE* pLng = _search<LINK_LONG_TYPE>(rSeqLim, MR_TAG_CONCATENATIONS))
    {
        this->m_pFctConcSetValue_orig   = pLng->registerSetValueHandler(fConcSetValue);
        this->m_pFctConcSolve_orig      = pLng->registerSolveHandler(fConcSolve);
    }

    //-----------------------------------------------------------------------------------
    // Register slice series mode handlers
    //-----------------------------------------------------------------------------------
    m_SliceSeriesMode.registerSolveHandler(rSeqLim, MR_TAG_SERIES_MODE, Ep2d_diff_UINS::fEPIDiffSolveSelectionConflict);

    //-----------------------------------------------------------------------------------
    // Register b value handlers
    //-----------------------------------------------------------------------------------
    m_BValue.registerElmSolveHandler(rSeqLim, MR_TAG_BVALUE, EpCommonUINS::fEPIStdSolveLongParam);
    m_BValue.registerSolveHandler(rSeqLim, MR_TAG_BVALUE, EpCommonUINS::fEPIStdSolveLongParam);

    // if thermal balancing is active, set solve handlers for the averages
    if (m_bThermalBalancing)
    {
        m_Averages.registerElmSolveHandler(rSeqLim, MR_TAG_DIFF_AVERAGES, EpCommonUINS::fEPIStdSolveLongParam);
        m_Averages.registerSolveHandler(rSeqLim, MR_TAG_DIFF_AVERAGES, EpCommonUINS::fEPIStdSolveLongParam);
    }
    else
    {
        // this try handler can only be used if thermal balancing is turned off
        m_BValue.registerTryHandler(rSeqLim, MR_TAG_BVALUE, Ep2d_diff_UINS::_bValueTry);
    }

    //-----------------------------------------------------------------------------------
    // Register averages / repetitions handlers
    //-----------------------------------------------------------------------------------
#ifdef SUPPORT_PACE
    m_RespComp.registerSetValueHandler(rSeqLim, MR_TAG_RESP_COMP, Ep2d_diff_UINS::fRespCompSetValue);
    m_RespComp.registerSolveHandler(rSeqLim, MR_TAG_RESP_COMP, Ep2d_diff_UINS::fEPIDiffSolveSelectionConflict);
#endif
    m_Measurements.registerSetValueHandler(rSeqLim, MR_TAG_MEASUREMENTS, NULL);
    m_Measurements.registerGetValueHandler(rSeqLim, MR_TAG_MEASUREMENTS, NULL);

    //-----------------------------------------------------------------------------------
    // Register bandwidth get-Limits handler
    //-----------------------------------------------------------------------------------
    char tMrTagBW[64];
    sprintf(tMrTagBW, "%s.0\0", MR_TAG_BANDWIDTH);

    m_Bandwidth.registerGetLimitsHandler(rSeqLim, tMrTagBW, Ep2d_diff_UINS::_bandwidthGetLimits);

    //-----------------------------------------------------------------------------------
    // Register RF pulse type handlers
    //-----------------------------------------------------------------------------------
    m_RFPulseType.registerSolveHandler(rSeqLim, MR_TAG_RFPULSE_TYPE, EpCommonUINS::fEPIStdSolveSelection);

    //-----------------------------------------------------------------------------------
    // Register fat saturation mode handlers
    //-----------------------------------------------------------------------------------
    m_FatSatMode.registerSolveHandler(rSeqLim, MR_TAG_FAT_SAT_MODE, EpCommonUINS::fEPIStdSolveSelection);

    //-----------------------------------------------------------------------------------
    // Remove Unwanted Options
    //-----------------------------------------------------------------------------------
    // Remove several 'save unfiltered' options. Reason: incompatible with Dynamic Field Correction (DFC)
    m_FilterNorm.registerIsAvailableHandler(rSeqLim, MR_TAG_FLT_NORM_PROP_UNFILTERED_IMAGES, EpCommonUINS::fEPIFilterBoolIsAvailable);
    m_FilterNormPreScan.registerIsAvailableHandler(rSeqLim, MR_TAG_FLT_PRESCAN_NORMALIZE_PROP_UNFILTERED_IMAGES, EpCommonUINS::fEPIFilterBoolIsAvailable);
    m_FilterNormBific.registerIsAvailableHandler(rSeqLim, MR_TAG_FLT_BIFIC_PROP_UNFILTERED_IMAGES, EpCommonUINS::fEPIFilterBoolIsAvailable);
    m_FilterDistCorr.registerIsAvailableHandler(rSeqLim, MR_TAG_FLT_DISCOR_PROP_UNFILTERED_IMAGES, EpCommonUINS::fEPIFilterBoolIsAvailable);

    //-------------------------------------------------------------------------------------
    // Register template set value handlers for automatic TE minimization
    //-------------------------------------------------------------------------------------

    // Again: consider that handlers might have already been overloaded in EpCommonUI!
    m_Bandwidth.registerSetValueHandler(rSeqLim, tMrTagBW, Ep2d_diff_UINS::fUILinkBandwidthSetValueNew);
    m_BaseResolution.registerSetValueHandler(rSeqLim, MR_TAG_BASE_RESOLUTION, Ep2d_diff_UINS::fUILinkBaseResolutionSetValueNew);
    m_PATMode.registerSetValueHandler(rSeqLim, MR_TAG_PAT_MODE, Ep2d_diff_UINS::fUILinkPATModeSetValueNew);
    m_PATAccelerationFactor.registerSetValueHandler(rSeqLim, MR_TAG_PAT_ACC_PE, Ep2d_diff_UINS::fUILinkPATFactorSetValueNew);
    m_PATReferenceLines.registerSetValueHandler(rSeqLim, MR_TAG_PAT_LINES_PE, Ep2d_diff_UINS::fUILinkPATRefLinesSetValueNew);
    m_PhaseFOV.registerSetValueHandler(rSeqLim, MR_TAG_PHASE_FOV, Ep2d_diff_UINS::fUILinkPhaseFOVSetValueNew);
    m_ReadFOV.registerSetValueHandler(rSeqLim, MR_TAG_READOUT_FOV, Ep2d_diff_UINS::fUILinkReadFOVSetValueNew);
    m_PhaseResolution.registerSetValueHandler(rSeqLim, MR_TAG_PHASE_RESOLUTION, Ep2d_diff_UINS::fUILinkPhaseResolutionSetValueNew);
    m_SliceThick.registerSetValueHandler(rSeqLim, MR_TAG_SLICE_THICKNESS, Ep2d_diff_UINS::fUILinkSliceThicknessSetValueNew);
    m_PhaseOS.registerSetValueHandler(rSeqLim, MR_TAG_PHASE_OVERSAMPLING, Ep2d_diff_UINS::fUILinkPhaseOSSetValueNew);
    m_FatSup.registerSetValueHandler(rSeqLim, MR_TAG_FAT_WATER_CONTRAST, Ep2d_diff_UINS::fUILinkFatWaterContrastSetValueNew);
    m_RFPulseType.registerSetValueHandler(rSeqLim, MR_TAG_RFPULSE_TYPE, Ep2d_diff_UINS::fUILinkRFPulseTypeSetValueNew);
    m_FatSatMode.registerSetValueHandler(rSeqLim, MR_TAG_FAT_SAT_MODE, Ep2d_diff_UINS::fUILinkFatSatModeSetValueNew);
    m_PhasePartF.registerSetValueHandler(rSeqLim, MR_TAG_PHASE_PARTIAL_FOURIER, Ep2d_diff_UINS::fUILinkPhasePFSetValueNew);
    m_TOM.registerSetValueHandler(rSeqLim, MR_TAG_TOM, Ep2d_diff_UINS::fUILinkTOMSetValueNew);
    m_GradMode.registerSetValueHandler(rSeqLim, MR_TAG_GRADIENT_MODE, Ep2d_diff_UINS::fUILinkGradModeSetValueNew);
    m_NumberDiffDirs.registerSetValueHandler(rSeqLim, MR_TAG_NUMBER_DIFF_DIRS, Ep2d_diff_UINS::fUILinkNumberDiffDirsSetValueNew);
    m_BValue.registerSetValueHandler(rSeqLim, MR_TAG_BVALUE, Ep2d_diff_UINS::fUILinkBValueSetValueNew);
    m_BValue.registerGetLimitsHandler(rSeqLim, MR_TAG_BVALUE, Ep2d_diff_UINS::fUILinkBValueGetLimitsNew);
    m_BValue.registerIsAvailableHandler(rSeqLim, MR_TAG_BVALUE, Ep2d_diff_UINS::fUILinkBValueIsAvailable);
    m_BValueSize.registerSetValueHandler(rSeqLim, MR_TAG_BVALUE_SIZE, Ep2d_diff_UINS::fUILinkBValueSizeSetValueNew);
    m_NumberSlices.registerSetValueHandler(rSeqLim, MR_TAG_SG_SIZE, Ep2d_diff_UINS::fUILinkNumberSlicesSetValueNew);
    m_DiffScheme.registerSetValueHandler(rSeqLim, MR_TAG_DIFF_SCHEME, Ep2d_diff_UINS::fUILinkDiffSchemeSetValueNew);
    m_QSpaceSteps.registerSetValueHandler(rSeqLim, MR_TAG_DIFF_QSPACE_STEPS, Ep2d_diff_UINS::fUILinkQSpaceStepsSetValueNew);
    m_QSpaceMaxBValue.registerSetValueHandler(rSeqLim, MR_TAG_DIFF_QSPACE_MAX_B_VALUE, Ep2d_diff_UINS::fUILinkQSpaceMaxBValueSetValueNew);
    m_QSpaceCoverage.registerSetValueHandler(rSeqLim, MR_TAG_DIFF_QSPACE_COVERAGE, Ep2d_diff_UINS::fUILinkQSpaceCoverageSetValueNew);
    m_QSpaceSampling.registerSetValueHandler(rSeqLim, MR_TAG_DIFF_QSPACE_SAMPLING_SCHEME, Ep2d_diff_UINS::fUILinkQSpaceSamplingSetValueNew);
    m_PatRefScanMode.registerSetValueHandler(rSeqLim, MR_TAG_PAT_REF_SCAN_MODE, Ep2d_diff_UINS::fUILinkPatRefScanModeSetValueNew);

    // m_DiffMode set value handler holds additional functionality and is registered above

    // Overloaded m_EchoSpacing set value handler has already been registered in EpCommonUI
    // => register a new handler here
    // => EpCommonUI-handler remains available via getOrigSetValueHandler
    // => Original (UILink) handler is no longer available and MUST NOT be called within EpCommonUI !!!
    // (Is there a more elegant solution that keeps track of all handlers (including the original one)???)
    m_EchoSpacing.registerSetValueHandler(rSeqLim, MR_TAG_ECHO_SPACING, Ep2d_diff_UINS::fUILinkEchoSpacingSetValueNew);

    // register acq window handlers
    m_AcqWindow.registerGetLimitsHandler(rSeqLim, MR_TAG_FIRST_ACQUISITION_WINDOW, Ep2d_diff_UINS::fAcquisitionWindowGetLimits);
    m_AcqWindow.registerSetValueHandler(rSeqLim, MR_TAG_FIRST_ACQUISITION_WINDOW, Ep2d_diff_UINS::fAcquisitionWindowSetValue);
    m_AcqWindow.registerSolveHandler(rSeqLim, MR_TAG_FIRST_ACQUISITION_WINDOW, Ep2d_diff_UINS::fAcquisitionWindowSolve);

    m_AcqWindowInternal.registerSetValueHandler(rSeqLim, MR_TAG_FIRST_ACQUISITION_WINDOW_FOR_INTERNAL, Ep2d_diff_UINS::fAcquisitionWindowInternalSetValue);

    m_SetCapture.registerSetValueHandler(rSeqLim, MR_TAG_CAPTURE_RR, Ep2d_diff_UINS::fSetCapture);

    m_PhaseCorrMode.registerSolveHandler(rSeqLim, MR_TAG_PHASE_CORRECTION, EpCommonUINS::fEPIStdSolveSelection);
    m_PhaseCorrMode.registerSetValueHandler(rSeqLim, MR_TAG_PHASE_CORRECTION, Ep2d_diff_UINS::fPhaseCorrSetValue);
    
#if defined BUILD_WIPParameterTool && defined SUPPORT_EC_COM
    m_EddyCurrentComp.registerSetValueHandler(
        rSeqLim, MR_TAG_EDDY_CURRENT_COMP, Ep2d_diff_UINS::fUIECCompensationSetValue);
    m_EddyCurrentComp.registerSolveHandler(
        rSeqLim, MR_TAG_EDDY_CURRENT_COMP, EpCommonUINS::fEPIStdSolveBoolParamConflict);
#endif

    if (SysProperties::isLowField())
    {
        m_FatSupOpt.registerSolveHandler(rSeqLim, MR_TAG_FAT_SUP_OPTIMIZATION, EpCommonUINS::fEPIStdSolveSelection);
    }

    //-----------------------------------------------------------------------------------
    // Register WIP parameter handlers 
    //-----------------------------------------------------------------------------------
#if WIP
    char ptWIPString[20] = "";

    if(Ep2d_diff_UINS::fWIPString(WIP_Anything, ptWIPString))
    {
        if(LINK_DOUBLE_TYPE* pDouble = _create< LINK_DOUBLE_TYPE >(rSeqLim, ptWIPString, WIP_Anything))
        {
            pDouble->registerGetLabelIdHandler(Ep2d_diff_UINS::_WIP_DOUBLE_GetLabelId);
            pDouble->registerGetUnitIdHandler(Ep2d_diff_UINS::_WIP_DOUBLE_GetUnitId);
            pDouble->registerGetValueHandler(Ep2d_diff_UINS::_WIP_DOUBLE_GetValue);
            pDouble->registerSetValueHandler(Ep2d_diff_UINS::_WIP_DOUBLE_SetValue);
            pDouble->registerGetLimitsHandler(Ep2d_diff_UINS::_WIP_DOUBLE_GetLimits);
            pDouble->registerGetToolTipIdHandler(Ep2d_diff_UINS::_WIP_DOUBLE_GetToolTipId);
            pDouble->registerIsAvailableHandler(Ep2d_diff_UINS::_WIP_DOUBLE_IsAvailable);
        }
    }
#endif   // of if WIP  
#endif // of #ifdef WIN32

    return MRI_SEQ_SEQU_NORMAL;
}

