//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\a_ep2d_diff_UI.h
//     Version: \main\18
//      Author: Clinical
//        Date: 2014-06-18 08:38:02 +02:00
//
//        Lang: C++
//
//
//
///  \file   a_ep2d_diff_UI.h
///  \brief  File containing declaraion of the UI class
///         - ep2d_diff
///
///  This file contains the implementation of the class Ep2d_diff_UI.
///
//    -----------------------------------------------------------------------------



#ifndef a_ep2d_diff_UI_h
#define a_ep2d_diff_UI_h 1



//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrImaging/seq/a_ep_CommonUI.h"
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"


#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkSelection.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkGeneric.h"
#endif

//  IMULT 
#include "MrImaging/SequenceLibraries/libPace/MODULE/MODULE_Routines.h"

//  -------------------------------------------------------------------------- *
//  Defines and typedefs                                                       *
//  -------------------------------------------------------------------------- *

// Required for b-value try handler protocol memory (see below)
const int iMaxProtMemory = 25;

namespace SEQ_NAMESPACE
{
    //  --------------------------------------------------------------------------
    //
    //  Name        : Ep2d_diff_UI
    //
    //  Description :
    /// \brief        This class basically is a storage for the pointers to the
    ///                original setValue / getValue / solve - handlers.
    ///
    ///               The sequence registers new UI handlers, which usually do
    ///                something, then call the original UI handler, and then
    ///                do something else. To keep the information of the original
    ///                UI handlers, the Ep2d_diff_UI class stores the pointers
    ///
    ///               It also provides the method registerUI to execute the
    ///                registration of all new handlers (and the storage of
    ///                 the original pointers)
    ///
    //  --------------------------------------------------------------------------

    class Ep2d_diff_UI: public EpCommonUI
    {

    public:

        //  --------------------------------------------------------------
        //
        //  Name        :  Ep2d_diff_UI::Ep2d_diff_UI
        //
        //  Description :
        /// \brief         Initialization of class members
        //
        //  Return      :
        //
        //  --------------------------------------------------------------
        Ep2d_diff_UI();


        //  --------------------------------------------------------------
        //
        //  Name        :  Ep2d_diff_UI::~Ep2d_diff_UI
        //
        //  Description :
        /// \brief         Destructor
        //
        //  Return      :
        //
        //  --------------------------------------------------------------
        virtual ~Ep2d_diff_UI();

        //  --------------------------------------------------------------------------
        //
        //  Name        : Ep2d_diff_UI::registerUI
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


#ifdef WIN32
        // ===========================================================================
        /// This function registers Didi to these UILink functions
        /**
            The UILink software relies on information about the available
            diffusion directions. Therefore we must provide access from UILink
            to Didi. This is accomplished by this global pointer address.
            It requires that this fUILinkRegisterDidi function is called
            during the initialization of the sequence.

            The effect of this funtion is just to write the supplied pointer
            address into the global variable pDidi.
            */
        virtual void fUILinkRegisterDidi(DiffusionDirections *DidiAddress);


		// set the thermal balancing flag for solve handlers
		void setThermalBalancing(bool bThermalBalancing);


        /// Pointer to the Didi diffusion directions manager actually used
        /**
        The UILink software relies on information about the available
        diffusion directions. Therefore we must provide access from UILink
        to Didi. This is accomplished by this global pointer address.
        It requires that the fUILinkRegisterDidi function is called
        during the initialization of the sequence.
        */
        DiffusionDirections *m_pDidi;

        /// Double array used to store sequence information for tooltips
        /**
        The UILink software (tooltips) relies on information about some
        diffusion encoding details. Therefore we must provide an appropriate
        interface. Here, we define a double array that is accessible from
        both, UILink and sequence.
        */
        enum eToolTipParamIndex
        {
            ToolTipParamDiffGradDuration = 0,
            ToolTipParamDiffGradSpacing
        };
        double m_dToolTipParam[2];

        /// String array used to store sequence information for tooltips
        /**
        The UILink software (tooltips) relies on information about some
        diffusion encoding details. Therefore we must provide an appropriate
        interface. Here, we define a double array that is accessible from
        both, UILink and sequence.
        */
        enum eToolTipStringIndex
        {
            ToolTipStringDVSInfo = 0
        };
        std::string m_sToolTipString[1];

        /// String used to store import / export error messages
        std::string m_sImportExportError;

        /// Import / export error status
        bool m_bImportExportError;

        /// Flag indicates that a an external diffusion vector set has been imported
        bool m_bVectorFileImported;

        //  --------------------------------------------------------------
        ///  \brief Helper class instances for UI handlers
        ///         - register new handler functions
        ///         - save pointer to original handler function
        ///         These classes exit only on the host.
        ///
        ///  The following line is an example which can be removed for
        ///  other sequences.
        //  --------------------------------------------------------------
        UI_ELEMENT_LONG      m_NumberDiffDirs;
        UI_ELEMENT_LONG      m_QSpaceSteps;
        UI_ELEMENT_LONG      m_Averages;
        UI_ELEMENT_LONG      m_Measurements;
        UI_ELEMENT_LONG      m_BValue;
        UI_ELEMENT_LONG      m_BValueSize;
        UI_ELEMENT_LONG      m_AcqWindow;
        UI_ELEMENT_LONG      m_AcqWindowInternal;
        UI_ELEMENT_BOOL      m_DiffRecon;
        UI_ELEMENT_BOOL      m_FilterDynDistCorr;
        UI_ELEMENT_BOOL      m_FilterNorm;
        UI_ELEMENT_BOOL      m_FilterNormPreScan;
        UI_ELEMENT_BOOL      m_FilterNormBific;
        UI_ELEMENT_BOOL      m_FilterDistCorr;
        UI_ELEMENT_BOOL      m_SetCapture;
        UI_ELEMENT_SELECTION m_DiffMode;
        UI_ELEMENT_SELECTION m_MultipleSeries;
        UI_ELEMENT_SELECTION m_RespComp;
        UI_ELEMENT_SELECTION m_Inversion;
        UI_ELEMENT_SELECTION m_SliceSeriesMode;
        UI_ELEMENT_SELECTION m_PatRefScanMode;
        UI_ELEMENT_DOUBLE    m_PTXPulseTxAccDiff;
        UI_ELEMENT_STRING    m_DiffDirsImport;
        UI_ELEMENT_STRING    m_DiffDirsExport;
        UI_ELEMENT_SELECTION m_PhaseCorrMode;

#ifdef BUILD_WIPParameterTool
        UI_ELEMENT_BOOL m_EddyCurrentComp;
#endif
        // The following set value handlers are required for the automatic TE minimization:
        // whenever one of these UI values is changed, the basicSetValue template handler
        // is called to set the minimum TE. Each parameter that affects the minimum TE
        // has to be considered here.

        // These are already members of the base class:
        // UI_ELEMENT_LONG      m_BaseResolution;
        // UI_ELEMENT_LONG      m_PATAccelerartionFactor;
        // UI_ELEMENT_LONG      m_PATReferenceLines;
        // UI_ELEMENT_DOUBLE    m_Bandwidth;
        // UI_ELEMENT_DOUBLE    m_EchoSpacing;
        // UI_ELEMENT_DOUBLE    m_ReadFOV;
        // UI_ELEMENT_DOUBLE    m_PhaseFOV;
        // UI_ELEMENT_DOUBLE    m_SliceThick;
        // UI_ELEMENT_DOUBLE    m_PhaseResolution;
        // UI_ELEMENT_DOUBLE    m_PhaseOS;
        // UI_ELEMENT_SELECTION m_PhasePartF;
        // UI_ELEMENT_SELECTION m_FatSup;
        // UI_ELEMENT_SELECTION m_PATMode;
        // UI_ELEMENT_SELECTION m_GradMode;
        // UI_ELEMENT_SELECTION m_PatRefScanMode;

        // These are already considered above
        // UI_ELEMENT_LONG      m_NumberDiffDirs;
        // UI_ELEMENT_LONG      m_BValue;
        // UI_ELEMENT_LONG      m_BValueSize;	
        // UI_ELEMENT_SELECTION m_DiffMode;

        // These have to be considered additionally
        UI_ELEMENT_LONG      m_NumberSlices;
        UI_ELEMENT_LONG      m_QSpaceMaxBValue;
        UI_ELEMENT_SELECTION m_QSpaceCoverage;
        UI_ELEMENT_SELECTION m_QSpaceSampling;
        UI_ELEMENT_SELECTION m_DiffScheme;
        UI_ELEMENT_SELECTION m_RFPulseType;
        UI_ELEMENT_SELECTION m_FatSatMode;
        UI_ELEMENT_SELECTION m_TOM;
        LINK_DOUBLE_TYPE::PFctGetLimits  m_pFctTRGetLimits_orig;
        LINK_DOUBLE_TYPE::PFctSetValue   m_pFctTRSetValue_orig;
        LINK_DOUBLE_TYPE::PFctSetValue   m_pFctTESetValue_orig;
        LINK_DOUBLE_TYPE::PFctGetLimits  m_pFctTIGetLimits_orig;
        LINK_DOUBLE_TYPE::PFctSetValue   m_pFctTISetValue_orig;
        LINK_LONG_TYPE::PFctSetValue     m_pFctSlicesSetValue_orig;
        LINK_LONG_TYPE::PFctSetValue     m_pFctConcSetValue_orig;
        LINK_LONG_TYPE::PFctSolve        m_pFctConcSolve_orig;
#endif

        // Memory for bValue try handler (required for UI speedup)
        // Static variables in the old code.
        typedef struct
        {
            long int lBValue;
            long int lNeededTE;
            long int lNeededTI;
            long int lNeededTR;
            bool     bValidProtocol;
            bool     bNeedOtherTETITR;
        } MyProtMemory;

        MrProt                               m_sMemProt;                    // Protocol memory for _bvalueTry
        MyProtMemory                         m_asMemory[iMaxProtMemory];    // Parameter store for _bvalueTry
        int                                  m_iMemCount;                   // Number of stored parameter sets for _bvalueTry

    protected:

	private:

		// thermal balancing flag
		bool m_bThermalBalancing;

    };

    namespace Ep2d_diff_UINS
    {

#ifdef WIN32
        // ===========================================================================
        /// Solve handler for 'DiffusionDirections' parameter
        /**
            This solve handler is called to allow the number of diffusion
            directions to be increased although TE is too short.
            The solve strategy in this parameter conflict is to call the standard
            solve handler which tries to increase TE and TR.
            */
        unsigned fUILinkNoDiffDirectionsSolve
            (
            LINK_LONG_TYPE*      const _this,     // ...
            char**               arg_list,  /**< receives confirmation message */
            const void*          pAddMem,   /**< additional memory needed for user b value array */
            const MrProtocolData::MrProtData*         pOrigProt, /**< Original protocol with old rf mode */
  int32_t                 lIndex     /**< Array index reserved */
            );

        // ===========================================================================
        /// Tooltip handler for 'DiffusionDirections' parameter
        /**
            In FREE mode, show user comment as tooltip
            */
unsigned fUILinkDiffDirectionsGetToolTipId (LINK_LONG_TYPE* const _this, char* arg_list[], int32_t lIndex );


        // ===========================================================================
        /// SetValue handler for 'Diffusion Mode' parameter
        /**
            This UILink handler is basically a copy of the original fUILinkDiffusionModeSetValue
            function found in \\n4\\pkg\\MrServers\\MrProtSrv\\MrProtocol\\UILink\\MrUILinkPerfDiff.cpp.
            The modification was necessary to implement the following actions:
            - For NOT(MDDW OR FREE), the MultipleSeriesMode must be MULTIPLE_SERIES_EACH_MEASUREMENT.
            - For FREE mode, the number of directions is set to the first vector set
            found in the external definition file, reconstruction of 'diffusion weighted images'
            will be switched on and all other reconstructions will be switched off.
            */
unsigned fUILinkDiffusionModeSetValueNew (LINK_SELECTION_TYPE* const _this, unsigned _x, int32_t pos);

        // ===========================================================================
        /// GetOptions handler for 'Diffusion Mode' parameter
        /**
            Calls original GetOptions handler and removes FREE mode from option
            vector if neither the current protocol includes a user defined vector
            set nor a corresponding import directory exists.
            */
bool fUILinkDiffusionModeGetOptions (LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);

        // ===========================================================================
        /// Tooltip handler for 'Diffusion Scheme' parameter
        /**
            Show delta (duration of diffusion encoding gradients) and Delta (distance
            of diffusion encoding gradients) if the diffusion scheme 'monopolar' is
            selected.
            */
unsigned fUILinkDiffusionSchemeGetToolTipId (LINK_SELECTION_TYPE* const _this, char* arg_list[], int32_t lIndex );

        // ===========================================================================
        /// Solve handler for 'Diffusion Scheme' parameter
        /**
            The solve handler can deactivate ZOOMit TX acceleration
            */
unsigned fUILinkSolveDiffusionScheme(LINK_SELECTION_TYPE* const pThis, char** arg_list, const void* pAddMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);

        // ===========================================================================
        /// Solve handler for 'pTX acceleration' parameter
        /**
            The solve handler can set the diffusion scheme to bipolar
            */
unsigned fUILinkSolvePTXAcceleration(LINK_DOUBLE_TYPE* const pThis, char** arg_list, const void* pAddMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);


        // ===========================================================================
        /// Generic solve handler for selection parameters
        /**
            This solve handler is called if a parameter conflict prevents
            the selection of one option. After solving diffusion variant
            specific conficts, the corresponding solve handler of the
            EPI common UI is called.
            */
unsigned fEPIDiffSolveSelectionConflict (LINK_SELECTION_TYPE* const pThis, char** arg_list, const void*, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);

        // ===========================================================================
        /// GetOptions handler for 'Dynamic Distortion Correction' parameter
        /**
            This getOptions handler disables the availability of the
            dynamic distortion correction if diffusion mode ONE_SCAN_TRACE
            is selected.
            */
bool fEPIDiffFltDynDistCorrGetOptions (LINK_BOOL_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos);

        // ===========================================================================
        /// GetValue handler for 'Diffusion Directions Import' parameter
        /**
            This getValue handler provides the file path for the diffusion
            directions import dialogue if diffusion mode FREE is selected.
            */
unsigned fUILinkDiffDirsImportGetValue( MrUILinkString* const pThis, char* arg_list[], int32_t lIndex );

        // ===========================================================================
        /// SetValue handler for 'Diffusion Directions Import' parameter
        /**
            Require Didi to load diffusion vector sets defined by the user
            in the file name given.
            */
void fUILinkDiffDirsImportSetValue( MrUILinkString* const pThis, const char* newVal, int32_t lIndex );

        // ===========================================================================
        /// IsAvailable handler for 'Diffusion Directions Import' parameter
        /**
            Parameter is visible if FREE mode is selected AND the directory
            containing external diffusion vector sets exists.
            */
bool fUILinkDiffDirsImportIsAvailable( MrUILinkString* const pThis, int32_t lIndex );

        // ===========================================================================
        /// Solve handler for 'Diffusion Directions Import' parameter
        /**
            Parameter is visible if FREE mode is selected AND the directory
            containing external diffusion vector sets exists.
            */
unsigned fUILinkDiffDirsImportSolve( MrUILinkString* const pThis, char** arg_list, const void*, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);

        // ===========================================================================
        /// GetValue handler for 'Diffusion Directions Export' parameter
        /**
            This getValue handler provides the file path for the diffusion
            directions export dialogue if diffusion mode FREE is selected.
            */
unsigned fUILinkDiffDirsExportGetValue( MrUILinkString* const pThis, char* arg_list[], int32_t lIndex );

        // ===========================================================================
        /// SetValue handler for 'Diffusion Directions Export' parameter
        /**
            Stores free diffusion vector set currently selected to file
            specified by the user.
            */
void fUILinkDiffDirsExportSetValue( MrUILinkString* const pThis, const char* newVal, int32_t lIndex );

        // ===========================================================================
        /// IsAvailable handler for 'Diffusion Directions Export' parameter
        /**
            Parameter is visible if FREE mode is selected AND the directory
            containing external diffusion vector sets exists.
            */
bool fUILinkDiffDirsExportIsAvailable( MrUILinkString* const pThis, int32_t lIndex );
       
        // ===========================================================================
        /// GetValue handler for MrFreeDiffusionData object
        /**
            Required for converting protocols using the FREE diffusion mode.
            */
MrPtr<MrGenericDC::IValueNode> fUIFreeDiffDataGetValueGeneric( LINK_GENERIC_TYPE* const pThis, int32_t lIndex );

        // ===========================================================================
        /// SetValue handler for MrFreeDiffusionData object
        /**
            Required for converting protocols using the FREE diffusion mode.
            */
MrPtr<MrGenericDC::IValueNode> fUIFreeDiffDataSetValueGeneric( LINK_GENERIC_TYPE* const pThis, const MrPtr<MrGenericDC::IValueNode>& sNewVal, int32_t lIndex );

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
        bool fUILinkNoDiffDirectionsGetLimits
            (
            LINK_LONG_TYPE* const     _this,          /*!< IMP: pointer to the UILink object for which the function is called */
            std::vector<MrLimitLong>& rLimitVector,   /*!< EXP: vector of limit intervals                                     */
    uint32_t&            rulVerify,      /*!< EXP: mode of further processing of the UILink Framework.           */
    int32_t                                      /*! unused in this case                                                 */
            );

        // ===========================================================================
        /// GetLimits handler for 'q-Space Weightings' parameter
        /**
            Q-space mode requires special handling, since only dedicated sampling
            schemes are supported.
            **/
        bool fUILinkNoQSpaceStepsGetLimits
            (
            LINK_LONG_TYPE* const     pThis,              /*!< IMP: Pointer to the UILink object for which the function is called */
            std::vector<MrLimitLong> &rLimitVector,       /*!< EXP: Vector of limit intervals                                     */
    uint32_t            &rulVerify,          /*!< EXP: Mode of further processing of the UILink Framework.           */
    int32_t                      lIndex              /*!< IMP: Index                                                         */
            );

        // ===========================================================================
        /// SetValue handler for 'Directions' parameter
        /**
            In FREE mode, the selected diffusion directions set is written to the
            free diffusion directions section of the protocol.
            */
int32_t fUILinkNoDiffDirectionsSetValue
            (
            LINK_LONG_TYPE* const     pThis,          /*!< IMP: pointer to the UILink object for which the function is called */
 int32_t                      lNewVal,        /*!< IMP: new value                                                     */
 int32_t                      lPos            /*!< IMP: position index                                                */
            );

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
bool _bValueTry (LINK_LONG_TYPE* const pThis, void* pVoid, const MrProtocolData::MrProtData*  pOrig, int32_t lIndex);

        // ------------------------------------------------------------------------------
        // Function    : _bandwidthGetLimits
        // ------------------------------------------------------------------------------
        //
        // Description : calls original getLimits-handler for bandwidth, but sets search
        //               mode to VERIFY_SCAN_ALL
        // Return      : whatever original getLimits-handler says
        //
        // ------------------------------------------------------------------------------
bool _bandwidthGetLimits (LINK_DOUBLE_TYPE * const pThis, std::vector< MrLimitDouble >& rLimitVector, uint32_t& rulVerify, int32_t lIndex);

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
unsigned fRespCompSetValue ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
#endif // #ifdef SUPPORT_PACE

unsigned fPhaseCorrSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos);

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
        bool fUILink_MultipleSeriesMode_IsAvailable
            (
            LINK_SELECTION_TYPE* const pThis /*!< IMP: pointer to the UILink object for which the function is called */,
 int32_t /*lIndex*/
            );


        // Template based set value handlers for automatic TE minimization
int32_t     fUILinkBaseResolutionSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
int32_t     fUILinkPATFactorSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
int32_t     fUILinkPATRefLinesSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
int32_t     fUILinkNumberDiffDirsSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
int32_t     fUILinkBValueSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
bool		fUILinkBValueIsAvailable(LINK_LONG_TYPE* const pThis, int32_t);
int32_t     fUILinkBValueSizeSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
int32_t     fUILinkNumberSlicesSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
int32_t     fUILinkQSpaceStepsSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
int32_t     fUILinkQSpaceMaxBValueSetValueNew(LINK_LONG_TYPE*      const pThis, int32_t     lNewVal, int32_t lPos);
double   fUILinkBandwidthSetValueNew       ( LINK_DOUBLE_TYPE*    const pThis, double   dNewVal, int32_t lPos );
double   fUILinkPhaseFOVSetValueNew        ( LINK_DOUBLE_TYPE*    const pThis, double   dNewVal, int32_t lPos );
double   fUILinkReadFOVSetValueNew         ( LINK_DOUBLE_TYPE*    const pThis, double   dNewVal, int32_t lPos );
double   fUILinkPhaseResolutionSetValueNew ( LINK_DOUBLE_TYPE*    const pThis, double   dNewVal, int32_t lPos );
double   fUILinkSliceThicknessSetValueNew  ( LINK_DOUBLE_TYPE*    const pThis, double   dNewVal, int32_t lPos );
double   fUILinkPhaseOSSetValueNew         ( LINK_DOUBLE_TYPE*    const pThis, double   dNewVal, int32_t lPos );
double   fUILinkEchoSpacingSetValueNew     ( LINK_DOUBLE_TYPE*    const pThis, double   dNewVal, int32_t lPos );
unsigned fUILinkPATModeSetValueNew         ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkFatWaterContrastSetValueNew( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkRFPulseTypeSetValueNew     ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkFatSatModeSetValueNew      ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkPhasePFSetValueNew         ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkTOMSetValueNew             ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkGradModeSetValueNew        ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkDiffSchemeSetValueNew      ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkQSpaceCoverageSetValueNew  ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkQSpaceSamplingSetValueNew  ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
unsigned fUILinkPatRefScanModeSetValueNew  ( LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, int32_t lPos );
bool     fUILinkBValueGetLimitsNew(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t pos);

bool     fAcquisitionWindowGetLimits(LINK_LONG_TYPE * const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t lIndex);
int32_t  fAcquisitionWindowSetValue(LINK_LONG_TYPE* const pThis, int32_t lNewScanWindow_ms, int32_t lIndex);
int32_t  fAcquisitionWindowInternalSetValue(LINK_LONG_TYPE* const pThis, int32_t lNewScanWindow_ms, int32_t lIndex);
unsigned fAcquisitionWindowSolve(LINK_LONG_TYPE* const pThis, char** arg_list, const void* pAddMem, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);

bool     fSetCapture(LINK_BOOL_TYPE* const pThis, bool lNewSetCapture, int32_t lIndex);



        /*[ Function ****************************************************************\
        *
        * Name        : TESetValue
        *
        * Description : This function sets the TE time of the specified contrast
        *
        * Return      : New TE[lContrast] time
        *
        \****************************************************************************/
        double TESetValue
            (
            MrUILinkBase* const pThis,              /*!< Imp: pointer to the UILink object for which the function is called */
            double              dDesiredTE_ms,      /*!< Imp: new TE value [ms]                                             */
    int32_t                lContrast           /*!< Imp: contrast ID                                                   */
            );

        /*[ Function ****************************************************************\
        *
        * Name        : StoreDiffusionDataToProt
        *
        * Description : Store user defined diffusion vector set to dedicated section
        *               of the protocol
        *
        * Return      : true in case of success, false otherwise
        *
        \****************************************************************************/
        bool StoreDiffusionDataToProt
            (
            MrUILinkBase* const pThis,              /*!< Imp: pointer to the UILink object for which the function is called */
    int32_t                lDiffDirs           /*!< Imp: number of diffusion directions                                */
            );

        // * -------------------------------------------------------------------------- *
        // *                        template set value handler                          *
        // * -------------------------------------------------------------------------- *

        /*[ Function ****************************************************************\
        *
        * Name        : basicSetValue
        *
        * Description : This function is a template set handler that can be for
        *               various UI parameters.
        *
        *               It calls the original UI set handler to set the relevant UI
        *               parameter. Afterwards TE is set to its minimum value (if desired).
        *               TE is updated only if TOM is enabled AND (currentValue and desiredValue
        *               are different OR bForcedTEUpdate is set).
        *
        * Return      : New value of the relevant UI parameter
        *
        \****************************************************************************/
        template <class ParamType, class Param, class OrigFuncType>
        Param basicSetValue(ParamType         pThis,
                            Param             desiredValue,
                      int32_t              lIndex,
                            OrigFuncType      pOriginalSetHandler,
                            bool              bForceTEUpdate = false)
        {




    int32_t         lTEMin           = 0;
    int32_t         lTRMin           = 0;
            Param        currentValue     = pThis->value(lIndex);

            MrProtocolData::MrProtData& rProt = pThis->prot();

            // * ------------------------------------------------------------------------ *
            // * Call original set function                                               *
            // * ------------------------------------------------------------------------ *
            if(pOriginalSetHandler)
            {
                desiredValue = (*pOriginalSetHandler)(pThis, desiredValue, lIndex);
            }

            MrProt  rMrProt(rProt);
            // * ------------------------------------------------------------------------ *
            // * Automatically set minimum TE                                             *
            // * ------------------------------------------------------------------------ *
            if(   // Only if TOM is enabled ...
               (rMrProt.TOM() == SEQ::TOM_MINIMIZE_TE) &&
               // ... and the value really changed OR the recalculation is enforced.
               // The latter is required if WIP_Selection SetValue handlers use this
               // template: although the value might be identical, the modifier might
               // be different.
               (!(currentValue == desiredValue) || bForceTEUpdate))
            {
                // * -------------------------------------------------------------------- *
                // * Calculate the min. TE and TR for the selected parameter set          *
                // * -------------------------------------------------------------------- *
                pThis->sequence().prepareForBinarySearch(rMrProt);

        for ( int32_t lI = 0; lI < rMrProt.contrasts(); lI++ )  
                {
                    lTEMin = getUI(pThis)->m_lNeededTE;

                    // * -------------------------------------------------------------------- *
                    // * Set the minimum TE                                                   *
                    // * -------------------------------------------------------------------- *
                    // * Set a valid TE value                                                 *
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

            if((rProt.getsPrepPulses().getucInversion() == SEQ::SLICE_SELECTIVE) && (rProt.getsPrepPulses().getucIRScheme() == SEQ::IR_SCHEME_SEQUENTIAL))
            {
                if(SEQ_NAMESPACE::Ep2d* pEPI = SEQ_NAMESPACE::getSeq(pThis))
                {
                    //  Used as indicator that the preparation of the kernel succeeded
                    pEPI->m_mySeqLoop.setlSBBScanTime(0);
                    pEPI->m_mySeqLoop.setlKernelScanTime(0);
                    pEPI->m_mySeqLoop.setCoolPauseWithinKernelTime_us(0);

                    pThis->sequence().prepareForBinarySearch(&rProt);

                    const int iSBBScanTime_us           = int(pEPI->m_mySeqLoop.getlSBBScanTime())
                        , iKernelTime_us                = int(pEPI->m_mySeqLoop.getlKernelScanTime())
                        , iTI_us                        = rProt.getalTI()[0]
                        , iNSlices                      = rProt.getsSliceArray().getlSize()
                        , iNConc                        = std::max(int(1), int(rProt.getsSliceArray().getlConc()))
                        , iNSlcPerConc                  = (iNSlices+iNConc-1)/iNConc
                        , iIRTime_us                    = int(const_cast<SeqBuildBlockIRsel*>(&pEPI->m_mySeqLoop.getSBBIRsel())->getDurationPerRequest())
                        , iCoolPauseWithinKernelTime_us = std::max(0, int(pEPI->m_lCoolPauseTotal-pEPI->m_lCoolPauseImplicit))
                        , iTRIncr_prefered_us           = 10000
                        ;
                    if(iKernelTime_us > 0)
                    {
                        int iKernelOffset = rProt.getsPrepPulses().getlKernelOffset();
                        int iTBlock_us = 0;
                        if(fCalcTBlock(iTBlock_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us-iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                        {
                            //  Valid solution found
                            int iTR_us = iNSlcPerConc*iTBlock_us;
                            const int iProtTR_us = rProt.getalTR()[0];
                            int iTRIncr_try_us = iTRIncr_prefered_us;
                            for(;iTRIncr_try_us >= GRAD_RASTER_TIME;iTRIncr_try_us/=10)
                            {
                                //  i) Find the least common multiple of iTRIncr_try_us and iNSlcPerConc
                                const int iLCM_us             = MODULE::LCM(iNSlcPerConc, iTRIncr_try_us)
                                    , iTBlock_prefered_us = MODULE::IMULT_CEIL(iTBlock_us, iLCM_us/iNSlcPerConc);

                                if(fTestTBlock(iTBlock_prefered_us, iTI_us, iIRTime_us, iSBBScanTime_us, iKernelTime_us-iCoolPauseWithinKernelTime_us, iCoolPauseWithinKernelTime_us, iKernelOffset))
                                {
                                    iTBlock_us = iTBlock_prefered_us;
                                    iTR_us     = iNSlcPerConc*iTBlock_prefered_us;
                                    break;
                                }
                            }
                            if(iTR_us != iProtTR_us)
                            {
                                if(LINK_DOUBLE_TYPE* pTR = _searchElm<LINK_DOUBLE_TYPE>(pThis, MR_TAG_TR))
                                {
                                    if(pTR->isAvailable(0))
                                    {
                                        pTR->value(iTR_us/1000., 0);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            return (desiredValue);
        }


#ifdef BUILD_WIPParameterTool
        bool fUIECCompensationSetValue(LINK_BOOL_TYPE* const pThis, bool value, int32_t lIndex);
#endif


#if WIP


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
        if ( LINK_LONG_TYPE* pLong=_create<LINK_LONG_TYPE>(rSeqLim, fWIPString(WIP_Spoiler), WIP_Spoiler) ) {
        pLong->registerGetValueHandler(_WIP_LONG_GetValue);
        ...
        }
        \endcode

        */
        // ===========================================================================
        bool fWIPString
            (
            int   index,                     /**< Input:  The index of the WIP mem block         */
            char  ptWIPString[20],           /**< Output: String tag                             */
            char  ptFamily[14] = "seq_wip"   /**< Input (optional): base tag. Should be one of
                                                  "seq_wip"(default), "seq_res", "eva_seq_wip".  */
                                                  );



        // ===========================================================================
        ///  Return the text to be placed in front of selection boxes
        unsigned _WIP_SELECTION_GetLabelId
            (
            LINK_SELECTION_TYPE* const _this,
            char*                      arg_list[],
 int32_t                       lIndex
            );


        // ===========================================================================
        /// Selects the text of the alternatives in the selection boxes
        int _WIP_SELECTION_Format
            (
            LINK_SELECTION_TYPE* const _this,
            unsigned                   nID,
            char*                      arg_list[],
 int32_t                       lIndex
            );



        // ===========================================================================
        /// Get the value of the selection box with ID lIndex
        unsigned _WIP_SELECTION_GetValue
            (
            LINK_SELECTION_TYPE* const _this,
 int32_t                       lIndex
            );



        // ===========================================================================
        /// Get the possible options of Selection box with ID lIndex
        bool _WIP_SELECTION_GetOptions
            (
            LINK_SELECTION_TYPE* const _this,
            std::vector<unsigned>&     rOptionVector,
 uint32_t&             rulVerify, 
 int32_t                       lIndex
            );




        // ===========================================================================
        /// Set a specified value to the selection box lIndex
unsigned _WIP_SELECTION_SetValue ( LINK_SELECTION_TYPE* const _this, unsigned nNewVal, int32_t lIndex );



        // ===========================================================================
        /// Decide if parameter is available on card
bool _WIP_SELECTION_IsAvailable ( LINK_SELECTION_TYPE* const _this, int32_t lIndex );



        // ===========================================================================
        /// Searches for the correct tool tip text for CheckBoxes
unsigned _WIP_SELECTION_GetToolTipId ( LINK_SELECTION_TYPE* const _this, char* arg_list[], int32_t lIndex );



        // ===========================================================================
/// Defines the text in front of the box for long values
unsigned _WIP_LONG_GetLabelId ( LINK_LONG_TYPE* const, char* arg_list[], int32_t lIndex );


        // ===========================================================================
/// Defines the text behind the box for long values
unsigned _WIP_LONG_GetUnitId ( LINK_LONG_TYPE* const, char* arg_list[], int32_t lIndex );



        // ===========================================================================
/// Get the long value entered in the box
int32_t _WIP_LONG_GetValue ( LINK_LONG_TYPE* const _this, int32_t lIndex );



        // ===========================================================================
/// Get the long value entered in the box
int32_t _WIP_LONG_SetValue ( LINK_LONG_TYPE* const _this, int32_t value, int32_t lIndex );



        // ===========================================================================
/// Defines the valid range of the long parameter
        bool _WIP_LONG_GetLimits
            (
            LINK_LONG_TYPE* const     _this,
            std::vector<MrLimitLong>& rLimitVector,
 uint32_t&            rulVerify, 
 int32_t                      lIndex
            );




        // ===========================================================================
        /// Decide if parameter is available on card
bool _WIP_LONG_IsAvailable ( LINK_LONG_TYPE* const _this, int32_t );




        // ===========================================================================
/// Searches for the correct tool tip text for long parameters
unsigned _WIP_LONG_GetToolTipId ( LINK_LONG_TYPE* const _this, char* arg_list[], int32_t lIndex );



        // ===========================================================================
        /// Defines the text in front of the box for double values
unsigned _WIP_DOUBLE_GetLabelId ( LINK_DOUBLE_TYPE* const _this, char* arg_list[], int32_t lIndex );



        // ===========================================================================
        /// Defines the text behind the box for double values
unsigned _WIP_DOUBLE_GetUnitId ( LINK_DOUBLE_TYPE* const _this, char* arg_list[], int32_t lIndex );



        // ===========================================================================
        /// Get the double value entered in the box
double _WIP_DOUBLE_GetValue ( LINK_DOUBLE_TYPE* const _this, int32_t lIndex );


        // ===========================================================================
        /// Get the double value entered in the box
double _WIP_DOUBLE_SetValue ( LINK_DOUBLE_TYPE* const _this, double value, int32_t lIndex );




        // ===========================================================================
        /// Defines the valid range of the double parameter
        bool _WIP_DOUBLE_GetLimits
            (
            LINK_DOUBLE_TYPE* const     _this,
            std::vector<MrLimitDouble>& rLimitVector,
 uint32_t&              rulVerify, 
 int32_t                        lIndex
            );



        // ===========================================================================
        /// Decide if parameter is available on card
bool _WIP_DOUBLE_IsAvailable ( LINK_DOUBLE_TYPE* const _this, int32_t lIndex );




        // ===========================================================================
        /// Searches for the correct tool tip text for double parameters
unsigned _WIP_DOUBLE_GetToolTipId ( LINK_DOUBLE_TYPE* const _This, char* arg_list[], int32_t lIndex );


#endif   // of if WIP
#endif   // of #ifdef WIN32

    } // end of namespace Ep2d_diff_UINS
} // end of namespace SEQ_NAMESPACE
#endif  // of #ifndef a_ep2d_diff_UI_h

