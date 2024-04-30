//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2019  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4\pkg\MrServers\MrImaging\seq\a_ep_SegmentationUI.h
//	 Version:
//	  Author: NEUR
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: Common UI for segmented ep2d variants
//
//	 Classes:
//
//	-----------------------------------------------------------------------------

#pragma once

//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrGlobalDefinitions/MrResult.h"                   // NLS_STATUS

#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"   // UI_ELEMENT_...
#endif


#if defined BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

#include <set>
#include <string>
#include <vector>

// Forward declarations
namespace MrProtocolData
{
    class SeqExpo;
}
class SeqLim;
class Sequence;
class MrProt;

namespace SEQ_NAMESPACE
{
    //  --------------------------------------------------------------------------
    //
    //  Name        : EpSegmentationUI
    //
    //  Description :
    /// \brief        This class basically is a storage for the pointers to the
    ///                original setValue / getValue / solve - handlers.
    ///         
    ///               All common handlers which are relevant in the context of 
    ///               segmented multi-shot EPI sequences are considered here,
    ///               including those required for an automated TE minimization.
    //  --------------------------------------------------------------------------

    class EpSegmentationUI
    {
    public:
        //  --------------------------------------------------------------
        //  Name        :  EpSegmentationUI::EpSegmentationUI
        //  Description :
        /// \brief         Constructor
        //  --------------------------------------------------------------
        EpSegmentationUI() = default;


        //  --------------------------------------------------------------
        //  Name        :  EpSegmentationUI::~EpSegmentationUI
        //  Description :
        /// \brief         Destructor
        //  --------------------------------------------------------------
        virtual ~EpSegmentationUI() = default;



        //  --------------------------------------------------------------------------
        //  Name        : EpSegmentationUI::registerUI
        //  Description : Register UI methods
        /// \brief        This function initializes the UI functions and
        ///                registers all given set / get / Solve - handlers
        //  --------------------------------------------------------------------------
        virtual NLS_STATUS registerUI (SeqLim &rSeqLim );

#ifdef WIN32

        //  --------------------------------------------------------------
        ///  \brief Helper class instances for UI handlers
        ///         - register new handler functions
        ///         - save pointer to original handler function
        ///         These classes exist only on the host.
        //  --------------------------------------------------------------
        UI_ELEMENT_BOOL       m_ManualEchoSpacing;

        UI_ELEMENT_LONG       m_BaseResolution;
        UI_ELEMENT_LONG       m_Contrasts;
#ifdef COMPILE_EP2D_DIFF
        UI_ELEMENT_LONG       m_DiffusionAverages;
#endif
        UI_ELEMENT_LONG       m_DistanceFactor;
        UI_ELEMENT_LONG       m_EPIFactor;
        UI_ELEMENT_LONG       m_Measurements;
        UI_ELEMENT_LONG       m_NoiseLevel;
        UI_ELEMENT_LONG       m_NumberOfSlices;
        UI_ELEMENT_LONG       m_PATAccelerationFactor;
        UI_ELEMENT_LONG       m_Segments;
        UI_ELEMENT_LONG       m_SliceAccelerationFactor;
        UI_ELEMENT_LONG       m_SliceAccelFOVShiftFactor;

        UI_ELEMENT_DOUBLE     m_Bandwidth;
        UI_ELEMENT_DOUBLE     m_EchoSpacing;
        UI_ELEMENT_DOUBLE     m_PhaseFOV;
        UI_ELEMENT_DOUBLE     m_PhaseOS;
        UI_ELEMENT_DOUBLE     m_PhaseResolution;
        UI_ELEMENT_DOUBLE     m_ReadFOV;
        UI_ELEMENT_DOUBLE     m_SliceDistance;
        UI_ELEMENT_DOUBLE     m_SliceThickness;
        UI_ELEMENT_DOUBLE     m_TE;

        UI_ELEMENT_SELECTION  m_CoilCombineMode;
        UI_ELEMENT_SELECTION  m_DistortionCorrection;
        UI_ELEMENT_SELECTION  m_FatSatMode;
        UI_ELEMENT_SELECTION  m_FatSup;
        UI_ELEMENT_SELECTION  m_GradMode;
        UI_ELEMENT_SELECTION  m_Inversion;
        UI_ELEMENT_SELECTION  m_PATMode;
        UI_ELEMENT_SELECTION  m_PATRefScanMode;
        UI_ELEMENT_SELECTION  m_PhaseCorrMode;
        UI_ELEMENT_SELECTION  m_PhasePF;
        UI_ELEMENT_SELECTION  m_POCS;
        UI_ELEMENT_SELECTION  m_RFPulseType;
        UI_ELEMENT_SELECTION  m_TOM;
        UI_ELEMENT_SELECTION  m_MTCMode;
        //UI_ELEMENT_SELECTION  m_WIP_ReconstructionMode;

#endif
    };

#ifdef WIN32

    namespace EpSegmentationUINS
    {
        //  --------------------------------------------------------------
        //   Overloader UI handlers
        //  --------------------------------------------------------------
        // Bool parameter handlers
        bool     fManualEchoSpacingSetValue       ( LINK_BOOL_TYPE*      const pThis, bool bNewValue,                                                                     int32_t lIndex );
        // Long parameter handlers
        int32_t  fBaseResolutionSetValue          ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
        bool     fBaseResolutionGetLimits         ( LINK_LONG_TYPE*      const pThis, std::vector<MrLimitLong>&   rLimitVector, uint32_t& rulVerify,                      int32_t lIndex );
        unsigned fBaseResolutionSolve             ( LINK_LONG_TYPE*      const pThis, char* arg_list[], const void* pAddMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
        int32_t  fContrastsSetValue               ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
#ifdef COMPILE_EP2D_DIFF
        int32_t  fDiffusionAveragesElmSetValue    ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
#endif
        int32_t  fDistanceFactorSetValue          ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
        int32_t  fEPIFactorGetValue               ( LINK_LONG_TYPE*      const pThis,                                                                                     int32_t lIndex );
        int32_t  fMeasurementsSetValue            ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
        bool     fNoiseLevelIsAvailable           ( LINK_LONG_TYPE*      const pThis,                                                                                     int32_t lIndex );
        int32_t  fNumberOfSlicesSetValue          ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
        int32_t  fPATAccelerationFactorSetValue   ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
        unsigned fPATAccelerationFactorSolve      ( LINK_LONG_TYPE*      const pThis, char* arg_list[], const void* pAddMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
        bool     fSegmentsGetLimits               ( LINK_LONG_TYPE*      const pThis, std::vector<MrLimitLong>&   rLimitVector, uint32_t& rulVerify,                      int32_t lIndex );
        int32_t  fSegmentsSetValue                ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
        unsigned fSegmentsSolve                   ( LINK_LONG_TYPE*      const pThis, char* arg_list[], const void* pAddMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
        int32_t  fSliceAccelerationFactorSetValue ( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
        int32_t  fSliceAccelFOVShiftFactorSetValue( LINK_LONG_TYPE*      const pThis, int32_t lNewValue,                                                                  int32_t lIndex );
        // Double parameter handlers
        double   fBandwidthElmSetValue            ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        double   fEchoSpacingSetValue             ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        double   fPhaseFOVSetValue                ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        unsigned fPhaseFOVSolve                   ( LINK_DOUBLE_TYPE*    const pThis, char* arg_list[], const void* pAddMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
        double   fPhaseOSSetValue                 ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        unsigned fPhaseOSSolve                    ( LINK_DOUBLE_TYPE*    const pThis, char* arg_list[], const void* pAddMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
        double   fPhaseResolutionSetValue         ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        unsigned fPhaseResolutionSolve            ( LINK_DOUBLE_TYPE*    const pThis, char* arg_list[], const void* pAddMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
        double   fReadFOVSetValue                 ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        double   fSliceDistanceElmSetValue        ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        double   fSliceThicknessSetValue          ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        double   fTEElmGetValue                   ( LINK_DOUBLE_TYPE*    const pThis,                                                                                     int32_t lIndex );
        double   fTEElmSetValue                   ( LINK_DOUBLE_TYPE*    const pThis, double dNewValue,                                                                   int32_t lIndex );
        bool     fTEElmGetLimits                  ( LINK_DOUBLE_TYPE*    const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify,                      int32_t lIndex );
        // Selection parameter handlers
        bool     fAsymmetricEchoModeGetOptions    ( LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify,                          int32_t lIndex );
        unsigned fAsymmetricEchoModeGetValue      ( LINK_SELECTION_TYPE* const pThis,                                                                                     int32_t lIndex );
        unsigned fAsymmetricEchoModeSetValue      ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        bool     fAsymmetricEchoModeIsAvailable   ( LINK_SELECTION_TYPE* const pThis,                                                                                     int32_t lIndex );
        unsigned fAsymmetricEchoModeGetToolTipId  ( LINK_SELECTION_TYPE* const pThis, char* arg_list[],                                                                   int32_t lIndex );
        unsigned fCoilCombineModeSetValue         ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fDistortionCorrectionSetValue    ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fFatSatModeSetValue              ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fFatSupSetValue                  ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fFlowCompensationElmSetValue     ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fGradModeSetValue                ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fInversionSetValue               ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fPATModeSetValue                 ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fPATModeSolve                    ( LINK_SELECTION_TYPE* const pThis, char* arg_list[], const void* pAddMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
        unsigned fPATRefScanModeSetValue          ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fPhaseCorrSetValue               ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fPhasePFSetValue                 ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fPhasePFSolve                    ( LINK_SELECTION_TYPE* const pThis, char* arg_list[], const void* pAddMem, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex );
        unsigned fPOCSSetValue                    ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fRFPulseTypeSetValue             ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fTOMSetValue                     ( LINK_SELECTION_TYPE* const pThis, unsigned uNewValue,                                                                 int32_t lIndex );
        unsigned fMTCModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewValue, int32_t lIndex);
    } // end of namespace EpSegmentationUINS
#endif // #ifdef WIN32

} // end of namespace SEQ_NAMESPACE

