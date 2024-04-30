//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep_CommonUI.h
//     Version: \main\17
//      Author: Clinical
//        Date: 2013-11-05 16:26:08 +01:00
//
//        Lang: C++
//
//
//
///  \file   a_ep_CommonUI.h
///  \brief  File containing declaration of the UI base class of
///         - ep2d_diff
///         - ep2d_perf
///         - ep2d_se
///         - ep2d_fid
///         - ep2d_bold
///         - ep2d_pace
///
///  This file contains the implementation of the class EpCommonUI.
///
//    -----------------------------------------------------------------------------



#ifndef a_ep_CommonUI_h
#define a_ep_CommonUI_h 1



//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"


#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkSelection.h"
#include "MrProtSrv/Domain/MrProtocol/UILink/MrStdNameTags.h"
#endif

#include "MrImaging/libSBB/SBBEPIReadOut.h"


#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

namespace MrProtocolData
{
    class SeqExpo;
}
class SeqLim;

namespace MrMeasSrv
{
    class ISequence;
}



typedef bool(*fEPIStdInit_pfCalculateTRTIFillTimes) (MrProtocolData::MrProtData*, SeqLim*, MrProtocolData::SeqExpo*, MrMeasSrv::ISequence*, long*, long*);

namespace SEQ_NAMESPACE
{
    //  --------------------------------------------------------------------------
    //
    //  Name        : EpCommonUI
    //
    //  Description :
    /// \brief        This class basically is a storage for the pointers to the
    ///                original setValue / getValue / solve - handlers.
    ///
    ///               The sequence registers new UI handlers, which usually do
    ///                something, then call the original UI handler, and then
    ///                do something else. To keep the information of the original
    ///                UI handlers, the EpCommonUI class stores the pointers
    ///
    ///               It also provides the method registerUI to execute the
    ///                registration of all new handlers (and the storage of
    ///                 the original pointers)
    ///
    //  --------------------------------------------------------------------------

    class EpCommonUI
    {
    public:

        //  --------------------------------------------------------------
        //
        //  Name        :  EpCommonUI::EpCommonUI
        //
        //  Description :
        /// \brief         Initialization of class members
        //
        //  Return      :
        //
        //  --------------------------------------------------------------
        EpCommonUI();


        //  --------------------------------------------------------------
        //
        //  Name        :  EpCommonUI::~EpCommonUI
        //
        //  Description :
        /// \brief         Destructor
        //
        //  Return      :
        //
        //  --------------------------------------------------------------
        virtual ~EpCommonUI();



        //  --------------------------------------------------------------------------
        //
        //  Name        : EpCommonUI::registerUI
        //
        //  Description :
        /// \brief        This function initializes the UI functions and
        ///                registers all given set / get / Solve - handlers
        ///
        ///               It can be executed on the measuement system, too, but is empty there.
        ///
        ///               On the host, it executes these steps
        ///               - Declaration of pointers to UI classes
        ///               - Registration of overloaded set value handlers
        ///
        ///               It returns an NLS status
        ///
        virtual NLS_STATUS registerUI(SeqLim &rSeqLim);

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
        virtual NLS_STATUS initializeUI(MrProt &rMrProt, SeqLim &rSeqLim);

        // ------------------------------------------------------------------------------
        // Name        : fEPIStdResetSolveHandlerControlTETITR
        // ------------------------------------------------------------------------------
        //
        // Description :
        // \brief      : resets solve handler control flag and needed TE,TR,TI.
        //
        // ------------------------------------------------------------------------------
        void fEPIStdResetSolveHandlerControlTETITR();

        // ------------------------------------------------------------------------------
        // Name        : fEPIStdUICheckNeededTETITR
        // ------------------------------------------------------------------------------
        //
        // Description :
        // \brief      : - Checks, if TE, TR, TI required by the sequence are identical
        //                 to the ones in the current protocol.
        //               - If required: Puts the required values for TE, TR, TI into the
        //                 current protocol so that it becomes consistent.
        //
        // ------------------------------------------------------------------------------
        bool fEPIStdUICheckNeededTETITR(MrProt &rMrProt, SeqLim &rSeqLim, int32_t lSEQNeedsTE, int32_t lSEQNeedsTI, int32_t lSEQNeedsTR, int32_t lTEContrastIndex = 0);

        // ------------------------------------------------------------------------------
        // Name        : fEPIStdInit
        // ------------------------------------------------------------------------------
        //
        // Description :
        // \brief      : - sets static variables of the module
        //               - specifies standard hard limits for epi sequences
        //
        // ------------------------------------------------------------------------------
        bool fEPIStdInit(SeqLim &rSeqLim, SeqBuildBlockEPIReadOut* _pEPIRO, ReorderInfo* _pREOInfo);

        // ------------------------------------------------------------------------------
        // Name        : fEPIStdRegisterTIHandlers
        // ------------------------------------------------------------------------------
        //
        // Description :
        // \brief      : registers TI handlers defined above
        //
        // ------------------------------------------------------------------------------
        bool fEPIStdRegisterTIHandlers(SeqLim &rSeqLim, fEPIStdInit_pfCalculateTRTIFillTimes _pfCalculateTRTIFill);

        // ------------------------------------------------------------------------------
        // Name        : fEPIStdRegisterEPIFactorHandlers
        // ------------------------------------------------------------------------------
        //
        // Description :
        // \brief      : registers EPI factor handlers defined above
        //
        // ------------------------------------------------------------------------------
        bool fEPIStdRegisterEPIFactorHandlers(SeqLim &rSeqLim);

        // ------------------------------------------------------------------------------
        // Name        : fEPIRegisterEP2D_DIFFHandlers
        // ------------------------------------------------------------------------------
        //
        // Description :
        // \brief      : register additional UI handlers relevant for the diffusion
        //               sequence variant
        //
        // ------------------------------------------------------------------------------
        bool fEPIRegisterEP2D_DIFFHandlers(SeqLim &rSeqLim);

        // get the flag designating a diffusion sequence
        bool isDiffusion() const
        {
            return m_isDiffusion;
        }

        // get the flag designating a spin-echo sequence
        bool isSpinEcho() const
        {
            return m_isSpinEcho;
        }

#ifdef WIN32

        //  --------------------------------------------------------------
        ///  \brief Helper class instances for UI handlers
        ///         - register new handler functions
        ///         - save pointer to original handler function
        ///         These classes exist only on the host.
        //  --------------------------------------------------------------
        UI_ELEMENT_LONG      m_BaseResolution;
        UI_ELEMENT_LONG      m_PATAccelerationFactor;
        UI_ELEMENT_LONG      m_AccelerationFactorSlice;
        UI_ELEMENT_LONG      m_FOVShiftFactor;
        UI_ELEMENT_LONG      m_PATReferenceLines;
        UI_ELEMENT_LONG      m_NumberOfSlices;
        UI_ELEMENT_DOUBLE    m_Bandwidth;
        UI_ELEMENT_DOUBLE    m_EchoSpacing;
        UI_ELEMENT_DOUBLE    m_TE;
        UI_ELEMENT_DOUBLE    m_TI;
        UI_ELEMENT_DOUBLE    m_ReadFOV;
        UI_ELEMENT_DOUBLE    m_PhaseFOV;
        UI_ELEMENT_DOUBLE    m_SliceThick;
        UI_ELEMENT_DOUBLE    m_PhaseResolution;
        UI_ELEMENT_DOUBLE    m_PhaseOS;
        UI_ELEMENT_DOUBLE    m_DelayInTR;
        UI_ELEMENT_SELECTION m_PhasePartF;
        UI_ELEMENT_SELECTION m_FatSup;
        UI_ELEMENT_SELECTION m_FatSupOpt;
        UI_ELEMENT_SELECTION m_Dimension;
        UI_ELEMENT_SELECTION m_GradMode;
        UI_ELEMENT_SELECTION m_PATMode;
        UI_ELEMENT_SELECTION m_PATRefScan;
        UI_ELEMENT_SELECTION m_Inversion;
        UI_ELEMENT_SELECTION m_AdjustmentMode;
        UI_ELEMENT_SELECTION m_LocalShim;
        UI_ELEMENT_SELECTION  m_DistortionCorrection;
        UI_ELEMENT_BOOL      m_CoilElements;
        UI_ELEMENT_BOOL      m_SFC;
        UI_ELEMENT_BOOL      m_FreezeSuppressedTissue;
        
#ifdef ZOOM_2DRF
        UI_ELEMENT_LONG      m_PTXPulsePhaseMatrix;
        UI_ELEMENT_DOUBLE    m_PTXVolRot;
        UI_ELEMENT_DOUBLE    m_PTXVolPDim;
        UI_ELEMENT_DOUBLE    m_PTXVolRDim;
        UI_ELEMENT_DOUBLE    m_PTXVolSDim;
        UI_ELEMENT_DOUBLE    m_PTXVolPosSag;
        UI_ELEMENT_DOUBLE    m_PTXVolPosSagSBCS;
        UI_ELEMENT_DOUBLE    m_PTXVolPosCor;
        UI_ELEMENT_DOUBLE    m_PTXVolPosCorSBCS;
        UI_ELEMENT_DOUBLE    m_PTXVolPosTra;
        UI_ELEMENT_DOUBLE    m_PTXVolPosTraSBCS;
        UI_ELEMENT_DOUBLE    m_PTXVolOriAlpha;
        UI_ELEMENT_DOUBLE    m_PTXVolOriBeta;
        UI_ELEMENT_DOUBLE    m_PTXPulseFlipAngle;
        UI_ELEMENT_DOUBLE    m_PTXPulseTxAcc;
        UI_ELEMENT_DOUBLE    m_PTXPulsePhaseFOE;
        UI_ELEMENT_SELECTION m_ExitPulse;
        UI_ELEMENT_SELECTION m_PtxVolProp;
        UI_ELEMENT_SELECTION m_PTXVolPos;
        UI_ELEMENT_SELECTION m_PTXVolPosSBCS;
        UI_ELEMENT_SELECTION m_PTXVolOri;
        UI_ELEMENT_SELECTION m_PTXVolOriHistory;
        UI_ELEMENT_SELECTION m_PTXPulseType;
        UI_ELEMENT_SELECTION m_AdjVolCoupling;
        UI_ELEMENT_SELECTION m_PTXPulseTrajectory;
        UI_ELEMENT_SELECTION m_PTXPulseB0Corr;
        UI_ELEMENT_ARRAY     m_PTXPulseArray;
        UI_ELEMENT_SELECTION m_PTXVolVisibility;

        // UI parameters for GSP pTX volume
        UI_ELEMENT_GRAD_DIR  m_PTXVolRotGSP;
        UINumericElement<LINK_VECTOR_TYPE> m_PTXVolPosGSP;
#endif

#endif

        int32_t m_lNeededTR;
        int32_t m_lNeededTI;
        int32_t m_lNeededTE;
        bool m_bNeedOtherTETITR;
        int32_t m_lDebug_SEQ_UILink;

        SeqBuildBlockEPIReadOut*             m_pEPIRO;
        ReorderInfo*                         m_pREOInfo;
        fEPIStdInit_pfCalculateTRTIFillTimes m_pfCalculateTRTIFill;

        //Static variables in the old code. They're now declared as member variables
        //to avoid vilation between threads.
        char m_tLine_epi_EchoSpacing_GetToolTip[100];

        // store if this is a diffusion sequence. due to the entangled behavior, not all functionality is moved to the diffusion UI
        bool m_isDiffusion{false};

        // store if this is a spin-echo sequence. due to the entangled behavior, not all functionality is moved to the ep2d_se UI
        bool m_isSpinEcho{false};
    };

#ifdef WIN32

    namespace EpCommonUINS
    {

        template<class UILinkParameterType> unsigned solveSliceAccelerationSettings(UILinkParameterType* const pThis);

#ifdef SUPPORT_IIR
        // Copy/Paste from a_tse_ui.cpp
        // Needed by TI solve handlers because SeqLoopIIR does not properly support the needed TR/TI mechanism.
        // SeqLoopIIR automatically assigns an optimum interleaving scheme (e.g. ieie), and the needed TR/TI
        // will depend on this selection. This might interfere with needed TR/TI based solve handling techniques.
        struct TE_TR_TI_HANDLER;
        unsigned fUISearchTR( MrUILinkBase* const _this, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
#endif

        // ------------------------------------------------------------------------------
        // Function    : _getTENeeded
        // ------------------------------------------------------------------------------
        //
        // Description : callback function needed for libUICtrl's standard solve handlers
        //               to get the required TE
        // Return      : TRUE , if parameter conflict can be solved with another TE
        //               FALSE, if not
        //
        // ------------------------------------------------------------------------------
        bool _getTENeeded(MrProtocolData::MrProtData*, SeqLim* rSeqLim, SeqExpo*, MrMeasSrv::ISequence*, long* alNeededTE);


        // ------------------------------------------------------------------------------
        // Function    : _getTRNeeded
        // ------------------------------------------------------------------------------
        //
        // Description : callback function needed for libUICtrl's standard solve handlers
        //               to get the required TR
        // Return      : TRUE , if parameter conflict can be solved with another TR
        //               FALSE, if not
        //
        // ------------------------------------------------------------------------------
        bool _getTRNeeded(MrProtocolData::MrProtData*, SeqLim* rSeqLim, SeqExpo*, MrMeasSrv::ISequence*, long* alNeededTR);


        // ------------------------------------------------------------------------------
        // Function    : _getTINeeded
        // ------------------------------------------------------------------------------
        //
        // Description : callback function needed for libUICtrl's standard solve handlers
        //               to get the required TI
        // Return      : TRUE , if parameter conflict can be solved with another TI
        //               FALSE, if not
        //
        // ------------------------------------------------------------------------------
        bool _getTINeeded(MrProtocolData::MrProtData*, SeqLim* rSeqLim, SeqExpo*, MrMeasSrv::ISequence*, long* alNeededTI);

        // --------------------------------------------------------------------------
        // Function    : pfGetTRTIFillTimes
        // --------------------------------------------------------------------------
        //
        // Description : Serve as a delegation of the function calculateTRTIFillTimes
        //
        // Return      : true , if success
        //               false, if not
        //
        //  --------------------------------------------------------------------------
        bool pfGetTRTIFillTimes
            (
            MrProtocolData::MrProtData  *rMrProt,
            SeqLim                      *rSeqLim,
            MrProtocolData::SeqExpo     *rSeqExpo,
            MrMeasSrv::ISequence        *pSeq,
            long                        *plNeededTI,
            long                        *plNeededTR
            );


        // *****************************************************************************
        //
        // Name        : AdaptRefLinesPE
        //
        // Description : This function checks if reference lines are within the range
        ///              If this is not the case, it adapts the no. of ref lines
        ///
        // *****************************************************************************
        void adaptRefLinesPE(MrUILinkBase* const pThis, bool bAlwaysSetToMax);


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
        unsigned fEPIStdSolveBoolParamConflict(MrUILinkSelection<bool> * const pThis, char **arg_list, const void * pMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex);


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
        unsigned fEPIStdSolveLongParam(LINK_LONG_TYPE* const pThis, char** arg_list, const void*, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);


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
        unsigned fEPIStdSolveDoubleParam(LINK_DOUBLE_TYPE* const pThis, char** arg_list, const void* pMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);


        // ------------------------------------------------------------------------------
        // Function    : fEPIStdSolveSelection
        // ------------------------------------------------------------------------------
        //
        // Description : invokes libUICtrl's standard solve handler for selction parameters
        // Return      : 0, if no solution possible
        //               MRI_STD_STRING on success. The text in arg_list is then used
        //               to format the confirmation message and the popup.
        //
        // ------------------------------------------------------------------------------
        unsigned fEPIStdSolveSelection(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);


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
        unsigned fEPIPATModeSolveSelection(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void*, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);

        // ---------------------------------------------------------------------------------
        // Function    : fUILinkPATModeSetValue
        // ---------------------------------------------------------------------------------
        // Description : Invokes standard set handler and adapts the number of PAT reference
        //               lines 
        // ---------------------------------------------------------------------------------
        unsigned fUILinkPATModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos);

        // ------------------------------------------------------------------------------
        // Function    : fEPILinesPEGetLimits
        // ------------------------------------------------------------------------------
        //
        // Description : Enforces that the limits for PAT reference lines is a multiple
        //               of the acceleration factor
        //
        // ------------------------------------------------------------------------------
        bool fEPIRefLinesPEGetLimits(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t /* pos */);
        bool fEPIRefLinesPEIsAvailable(LINK_LONG_TYPE* const pThis, int32_t /*pos*/);

        // ------------------------------------------------------------------------------
        // Function    : fEPIAccPESetValue
        // ------------------------------------------------------------------------------
        //
        // Description : Enforces that the limits for PAT reference lines is a multiple
        //               of the acceleration factor
        //
        // ------------------------------------------------------------------------------
        int32_t fEPIAccPESetValue(LINK_LONG_TYPE* const pThis, int32_t val, int32_t pos);

        void _addTimingDependencyPtr(MrUILinkBase* const pThis);

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
        unsigned _solveBaseResolution(LINK_LONG_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);



        // ------------------------------------------------------------------------------
        // Function    : _epi_BaseResolutionTry
        // ------------------------------------------------------------------------------
        //
        // Description : resets m_bNeedOtherTETITR1 and calls original try-handler for base-
        //               resolution
        // Return      : whatever original try-handler says
        //
        // ------------------------------------------------------------------------------
        static bool _epi_BaseResolutionTry(LINK_LONG_TYPE* const pThis, void* pVoid, const MrProtocolData::MrProtData*  pOrig, int32_t lIndex);


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
        unsigned _solveTETITR_TE(LINK_DOUBLE_TYPE* const pThis, char** arg_list, const void*, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);


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
        static unsigned _solveTETITR_PPF_EPIFactorConflict(MrUILinkLimited<int32_t> * const pThis, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);


        // ------------------------------------------------------------------------------
        // Function    : _epi_EchoSpacingIsAvailable
        // ------------------------------------------------------------------------------
        //
        // Description : determines, if parameter echo spacing can be seen on the UI
        // Return      : true or false
        //
        // ------------------------------------------------------------------------------
        bool _epi_EchoSpacingIsAvailable(LINK_DOUBLE_TYPE* const /*pThis*/, int32_t /*lIndex*/);


        // ------------------------------------------------------------------------------
        // Function    : _epi_EchoSpacingGetValue
        // ------------------------------------------------------------------------------
        //
        // Description : determines value of echo spacing to be shown on the UI
        // Return      : echo sapcing [ms]
        //
        // ------------------------------------------------------------------------------
        double _epi_EchoSpacingGetValue(LINK_DOUBLE_TYPE* const pThis, int32_t /*lIndex*/);


        // ------------------------------------------------------------------------------
        // Function    : _epi_EchoSpacingSetValue
        // ------------------------------------------------------------------------------
        //
        // Description : Called when user modifies echo spacing in the UI. Puts appropriate
        //               value into the protocol.
        // Return      : new value of echo spacing
        //
        // ------------------------------------------------------------------------------
        double _epi_EchoSpacingSetValue(LINK_DOUBLE_TYPE* const pThis, double dNewVal_ms, int32_t lIndex);


        // ------------------------------------------------------------------------------
        // Function    : _epi_EchoSpacingGetLimits
        // ------------------------------------------------------------------------------
        //
        // Description : Calculates limits for echo spacing. Needed, because the echo
        //               spacing has a forbidden range due to resonances of the gradient
        //               coil.
        // Return      : false, if echospacing may not be edited
        //                      or if no valid limits can be found
        //               true, if echo spacing has valid limits
        //
        // ------------------------------------------------------------------------------
        bool _epi_EchoSpacingGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector< MrLimitDouble >& rLimitVector, uint32_t& rulVerify, int32_t /*lIndex*/);


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
        unsigned _epi_EchoSpacing_GetToolTip(LINK_DOUBLE_TYPE* const pThis, char* arg_list[], int32_t);


        // ------------------------------------------------------------------------------
        // Function    : _epi_FreeEchoSpacing_GetLabelId
        // ------------------------------------------------------------------------------
        //
        // Description : determines NLS supported text shown on the UI
        // Return      : ressource-ID
        //
        // ------------------------------------------------------------------------------
        unsigned _epi_FreeEchoSpacing_GetLabelId(LINK_BOOL_TYPE* const, char* /*arg_list*/[], int32_t /*lIndex*/);


        // ------------------------------------------------------------------------------
        // Function    : _epi_FreeEchoSpacing_GetOptions
        // ------------------------------------------------------------------------------
        //
        // Description : determines options the user may select
        // Return      : true
        //
        // ------------------------------------------------------------------------------
        bool _epi_FreeEchoSpacing_GetOptions(LINK_BOOL_TYPE* const /*pThis*/, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t /*lIndex*/);


        // ------------------------------------------------------------------------------
        // Function    : _epi_FreeEchoSpacing_GetValue
        // ------------------------------------------------------------------------------
        //
        // Description : determines current setting of free echo spacing
        // Return      : true or false
        //
        // ------------------------------------------------------------------------------
        bool _epi_FreeEchoSpacing_GetValue(LINK_BOOL_TYPE* const pThis, int32_t /*lIndex*/);


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
        bool _epi_FreeEchoSpacing_SetValue(LINK_BOOL_TYPE* const pThis, bool value, int32_t /*lIndex*/);


        // ------------------------------------------------------------------------------
        // Function    : _epi_FreeEchoSpacing_Solve
        // ------------------------------------------------------------------------------
        //
        // Description : free echo spacing solve handler
        // Return      : 0, i.e. we can't help
        //
        // ------------------------------------------------------------------------------
        unsigned _epi_FreeEchoSpacing_Solve(LINK_BOOL_TYPE* const /*pThis*/, char* /*arg_list*/[], const void* /*pVoid*/, const MrProtocolData::MrProtData *  /*rMrProt*/, int32_t /*lIndex*/);

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
        static bool _TIGetLimits(LINK_DOUBLE_TYPE * const pThis, std::vector< MrLimitDouble >& rLimitVector, uint32_t& rulVerify, int32_t lIndex);


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
        static unsigned _solveTETITR_TI(MrUILinkLimited<double>* const pThis, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex);


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
        bool fEPIGradientModeGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t /*pos*/);


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
        unsigned fEPIGradientModeGetValue(LINK_SELECTION_TYPE* const pThis, int32_t /*pos*/);


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
        unsigned fEPIGradientModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned newVal, int32_t pos);

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
        bool fEPISaveUncombinedIsAvailable(LINK_BOOL_TYPE* const /* pThis */, int32_t);

        bool fEPIFltRawGetOptions(LINK_BOOL_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);

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
        unsigned fEPITriggeringSolve
            (
            LINK_SELECTION_TYPE*  const _this,   // ...
            char**                arg_list,      // receives confirmation message
            const void*       /*  pToAddMemory */,
            const MrProtocolData::MrProtData*          pOrigProt,     // Original protocol
            int32_t                  lIndex          // Array index reserved
            );

        // ------------------------------------------------------------------------------
        // Function    : fEPIFilterBoolIsAvailable()
        // ------------------------------------------------------------------------------
        //
        // Overloaded is-available handler for "Save unfiltered".
        //
        // This overloaded version is used to remove the option
        // "Save unfiltered" from various filters. It always returns false.
        // ------------------------------------------------------------------------------
        bool fEPIFilterBoolIsAvailable(LINK_BOOL_TYPE* const /*pThis*/, int32_t /*pos*/);

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
        __IMP_EXP bool updateB1ControlLoopParametersInProtocol
            (
            MrProt &rMrProt               /* reference to the protocol */
            );

        unsigned fUILinkAdjustmentModeSolve(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lPos);
        unsigned fUILinkAdjustmentModeSetValue  (LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
        bool     fUILinkLocalShimTry(LINK_SELECTION_TYPE* const pThis, void* pVoid, const MrProtocolData::MrProtData*  pOrig, int32_t pos);

        // ------------------------------------------------------------------------------
        // handler definition for slice acceleration
        // ------------------------------------------------------------------------------
        int32_t  fUILinkSliceAccelSetValue(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
        bool     fUILinkSliceAccelGetLimits(LINK_LONG_TYPE*      const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t lPos);
        unsigned fUILinkSliceAccelSolve(LINK_LONG_TYPE*      const pThis, char** arg_list, const void* pAddMem, const MrProtocolData::MrProtData * pOrigProt, int32_t lIndex);

        bool     fUILinkFOVShiftFactorGetLimits(LINK_LONG_TYPE*      const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t lPos);
        bool     fUILinkFOVShiftFactorIsAvailable(LINK_LONG_TYPE*      const pThis, int32_t lPos);

        //        int      fUILinkPatLightFormat             (LINK_BOOL_TYPE*       const pThis, bool nID, char* arg_list[], int32_t pos);

        bool     fUILinkNumberOfSlicesGetLimits(LINK_LONG_TYPE*      const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t lPos);

        unsigned fUILinkPATAccelFactorGetToolTipID(LINK_LONG_TYPE*      const pThis, char* arg_list[], int32_t lPos);

        unsigned fUILinkFatSuppresionSolve(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lPos);

        bool     fUILinkCoilElementsSetValue(LINK_BOOL_TYPE* const pThis, bool bDesiredState, int32_t lIndex);

        unsigned fDistortionCorrectionSetValue    ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );

        bool fUILinkSFCSetValue(LINK_BOOL_TYPE* const pThis, bool bDesiredState, int32_t lIndex);
    
        // ------------------------------------------------------------------------------
        //  Function    : fUILinkFreezeSuppressedTissueSetValue()
        // ------------------------------------------------------------------------------
        //  Description : Puts TI on raster when disabling freezing of suppressed tissue.
        // ------------------------------------------------------------------------------
        bool     fUILinkFreezeSuppressedTissueSetValue(LINK_BOOL_TYPE* const pThis, bool bDesiredState, int32_t lIndex);

    } // end of namespace EpCommonUINS
#endif // #ifdef WIN32
} // end of namespace SEQ_NAMESPACE

#endif

