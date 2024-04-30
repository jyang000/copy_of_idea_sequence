//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/X
//        File: \src\MrImaging\seq\a_tgse_asl\a_tgse_asl_UI.h
//      Author: pfeujodj
//        Date: 2018-07-03 11:44:21 +02:00
//
//        Lang: C++
//
//
//
///  \file   a_tgse_asl_UI.h
///  \brief  File containing declaration of the UI class
///         - tgse_asl
///
///  This file contains the implementation of the class Tgse_asl_UI.
///
//    -----------------------------------------------------------------------------
#pragma once

//  -------------------------------------------------------------------------- *
//  Application includes                                                       *
//  -------------------------------------------------------------------------- *
#include "MrImaging/seq/a_ep_CommonUI.h"
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"

#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkSelection.h"
#else
class LINK_LONG_TYPE;
#endif


//  -------------------------------------------------------------------------- *
//  Defines and typedefs                                                       *
//  -------------------------------------------------------------------------- *

//  -------------------------------------------------------------------------- *
//  Forward declarations                                                       *
//  -------------------------------------------------------------------------- *
namespace MrProtocolData
{
class MrProtData;
}
class SeqLim;
class SeqExpo;
class Sequence;

namespace SEQ_NAMESPACE
{

  const bool ASL_UI_bIsAdvancedMode = true;

//  --------------------------------------------------------------------------
//
//  Name        : Tgse_asl_UI
//
//  Description :
/// \brief        This class basically is a storage for the pointers to the
///                original setValue / getValue / solve - handlers.
///
///               The sequence registers new UI handlers, which usually do
///                something, then call the original UI handler, and then
///                do something else. To keep the information of the original
///                UI handlers, the Tgse_asl_UI class stores the pointers
///
///               It also provides the method registerUI to execute the
///                registration of all new handlers (and the storage of
///                 the original pointers)
///
//  --------------------------------------------------------------------------

class Tgse_asl_UI: public EpCommonUI
{
public:
    //  --------------------------------------------------------------
    //
    //  Name        :  Tgse_asl_UI::Tgse_asl_UI
    //
    //  Description :
    /// \brief         Initialization of class members
    //
    //  Return      :
    //
    //  --------------------------------------------------------------
	  Tgse_asl_UI();

    //  --------------------------------------------------------------
    //
    //  Name        :  Tgse_asl_UI::~Tgse_asl_UI
    //
    //  Description :
    /// \brief         Destructor
    //
    //  Return      :
    //
    //  --------------------------------------------------------------
    virtual ~Tgse_asl_UI();

    //  --------------------------------------------------------------------------
    //
    //  Name        : Tgse_asl_UI::registerUI
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
    //  Name        : Tgse_asl_UI::calculateKspace
    //
    //  Description :
    /// \brief        This function provides central computation for
    ///               kspace and fastimaging dependencies: it supports UI changes in
    ///               base resolution, phaseOS, sliceOS, FOV phase, phase resolution
    ///
    ///               - Number of phase encoding lines and partitions are set
    ///               - Required changes to Turbo/Epi Factors and segments are made
    ///
    ///               A call to prepareForBinarySearch after this method is typically needed
    ///
    ///               It returns an NLS status
    ///
    virtual NLS_STATUS calculateKspace(MrProtocolData::MrProtData  *rMrProt);

#ifdef SUPPORT_iPAT_TGSE
    double getPartialFourierFactor(SEQ::PartialFourierFactor pffPFF) const;

    int32_t measuredPartitions(MrProtocolData::MrProtData  &rMrProt);
    int32_t measuredLinesPE   (MrProtocolData::MrProtData  &rMrProt);

    template<typename T> NLS_STATUS changeTurboFactor(T);
    template<typename T> NLS_STATUS changeRefLines3D(T);
    NLS_STATUS changeSliceOS(LINK_LONG_TYPE* const);
    template<typename T> NLS_STATUS changeEPIFactor(T);
    template<typename T> NLS_STATUS changeRefLinesPE(T);
    NLS_STATUS changePhaseOS(LINK_LONG_TYPE* const);
#endif

    //Here we want to store some user history for soft paramters
    double m_dWantedPhaseRes;
    double m_dWantedPhaseFOV;

    static bool m_bFirst_Time;

    inline static bool getFirstTime()         { return m_bFirst_Time; }
    inline static void setFirstTime(bool t)   { m_bFirst_Time = t; }

#ifdef WIN32

    UI_ELEMENT_SELECTION m_FatSatMode;
    UI_ELEMENT_SELECTION m_AslMode;
    UI_ELEMENT_BOOL      m_SaveUncombined;
    UI_ELEMENT_BOOL      m_FilterNorm;
    UI_ELEMENT_BOOL      m_FilterPrescan;
    UI_ELEMENT_BOOL      m_FilterDiscor;
    UI_ELEMENT_BOOL      m_FilterImage;
    UI_ELEMENT_BOOL      m_FilterNormBific;
    UI_ELEMENT_DOUBLE    m_AslFlowCrush;
    UI_ELEMENT_DOUBLE    m_TE;
    UI_ELEMENT_LONG      m_Segments;
    UI_ELEMENT_LONG      m_Repetitions;
    UI_ELEMENT_LONG      m_DelayArraySize;
    UI_ELEMENT_LONG      m_InversionArraySize;
    UI_ELEMENT_DOUBLE    m_InversionTimeElm;
    UI_ELEMENT_LONG      m_epiFactor;
    UI_ELEMENT_LONG      m_TurboFactor;
    UI_ELEMENT_LONG      m_Partitions;
    UI_ELEMENT_LONG      m_BaseResolution;
    UI_ELEMENT_SELECTION m_PhasePF;
    UI_ELEMENT_DOUBLE    m_PhaseFOV;
    UI_ELEMENT_DOUBLE    m_PhaseRes;
    UI_ELEMENT_DOUBLE    m_SliceResolution;
    UI_ELEMENT_DOUBLE    m_SliceOS;
    UI_ELEMENT_DOUBLE    m_PhaseOS;
    UI_ELEMENT_BOOL      m_FreeEchoSpacing;
    UI_ELEMENT_DOUBLE    m_AslBolusDuration;
    UI_ELEMENT_LONG      m_LabelingDuration;
    UI_ELEMENT_LONG      m_PostLabelingDelayElm;
    UI_ELEMENT_SELECTION m_AslSuppression;
    UI_ELEMENT_SELECTION m_PTSat;

#endif

protected:

};


namespace Tgse_asl_UINS
{

#ifdef WIN32

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

// ------------------------------------------------------------------------------
// Function    : _UnavailableOption
// ------------------------------------------------------------------------------
//
// Overloaded is-available handlers for generic removal of UI options
//
// ------------------------------------------------------------------------------
bool _UnavailableOption(LINK_BOOL_TYPE* const /*pThis*/, int32_t /*pos*/);

bool _UnavailableOptionDouble(LINK_DOUBLE_TYPE* const /*pThis*/, int32_t /*pos*/);

bool _AvailableOptionLong(LINK_LONG_TYPE* const /*pThis*/, int32_t /*pos*/);

// ------------------------------------------------------------------------------
// Function    : 
// ------------------------------------------------------------------------------
//
// Revised handling of the measurements UI to accommodate 
// for multiple TI capabilities
//
// ------------------------------------------------------------------------------
#ifndef SUPPORT_CSL // Handlers for repetitions disabled as we use phases for encoding number of TIs
int32_t _asl_RepetitionsGetValue(LINK_LONG_TYPE* const, int32_t);
#endif

int32_t _asl_RepetitionsSetValue(LINK_LONG_TYPE* const, int32_t, int32_t);


unsigned _asl_AslLabelingDurationGetLabelId(LINK_LONG_TYPE* const /*pThis*/, char* /*arg_list*/[], int32_t /*pos*/);

unsigned _asl_AslLabelingDurationGetTooltipId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t pos);

int32_t _asl_DelayArraySizeGetValue(LINK_LONG_TYPE* const, int32_t);

int32_t _asl_DelayArraySizeSetValue(LINK_LONG_TYPE* const, int32_t, int32_t);

unsigned _asl_DelayArraySizeGetLabelId(LINK_LONG_TYPE* const /*pThis*/, char* /*arg_list*/[], int32_t /*pos*/);

unsigned _asl_DelayArraySizeGetTooltipId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t pos);

//unsigned _asl_DelayArraySizeSolve(LINK_LONG_TYPE* const, char**, const void*, const MrProtocolData::MrProtData*, int32_t);

int32_t _asl_PostLabelingDelayMaxSize(MrUILinkArray* const pThis, int32_t /*pos*/);

unsigned _asl_PostLabelingDelayElmGetLabelId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t pos);

int32_t _asl_PostLabelingDelayElmGetValue(LINK_LONG_TYPE* const pThis, int32_t pos);

int32_t _asl_PostLabelingDelayElmSetValue(LINK_LONG_TYPE* const pThis, int32_t val, int32_t pos);

bool _asl_PostLabelingDelayElmGetLimits(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t pos);

unsigned _asl_PostLabelingDelayGetTooltipId(LINK_LONG_TYPE* const pThis, char* arg_list[], int32_t pos);

int32_t _asl_InversionArraySizeGetValue(LINK_LONG_TYPE* const, int32_t);

int32_t _asl_InversionArraySizeSetValue(LINK_LONG_TYPE* const, int32_t, int32_t);

bool   _asl_InversionArraySizeGetLimits(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t lpos);

double _asl_InversionTimeElmGetValue(LINK_DOUBLE_TYPE* const pThis, int32_t lpos);

double _asl_InversionTimeElmSetValue(LINK_DOUBLE_TYPE* const pThis, double dNewVal, int32_t lpos);

bool   _asl_InversionTimeElmGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t lpos);

unsigned _asl_InversionTimeElmFlowLimitSolve(LINK_DOUBLE_TYPE* const pThis, char* arg_list[], const void*  pAddMem, const MrProtocolData::MrProtData*  pOrig, int32_t lPos);

double _asl_BolusDurationSetValue(LINK_DOUBLE_TYPE* const pThis, double val, int32_t /*pos*/);

bool   _asl_BolusDurationGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t lpos);

bool   _asl_BolusDurationIsAvailable(LINK_DOUBLE_TYPE* const pThis, int32_t lpos);

//  -----------------------------------------------------------------
//  EPI-Factor - recalculate segments for display in UI
//  -----------------------------------------------------------------

bool _asl_EPIFactorIsAvailable(LINK_LONG_TYPE* const pThis, int32_t pos);

int32_t _asl_EPIFactorSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t pos);

int32_t _asl_EPIFactorGetValue(LINK_LONG_TYPE* const pThis, int32_t pos);

bool _asl_EPIFactorGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t pos);

//  -----------------------------------------------------------------
//  Turbo-Factor - recalculate segments for display in UI
//  -----------------------------------------------------------------

int32_t _asl_TurboSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t pos);

int32_t _asl_TurboGetValue(LINK_LONG_TYPE* const pThis, int32_t pos);

bool _asl_TurboGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t pos);

//  -----------------------------------------------------------------
//  Images per Slab - recalculate turbo factor
//  -----------------------------------------------------------------
int32_t _asl_ImagesPerSlabSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t pos);

int32_t _asl_ImagesPerSlabGetValue(LINK_LONG_TYPE* const pThis, int32_t pos);

bool _asl_ImagesPerSlabGetLimits(LINK_LONG_TYPE* const pThis, std::vector< MrLimitLong >& rLimitVector, uint32_t& rulVerify, int32_t pos);

#ifndef SUPPORT_CSL
bool _asl_ImagesPerSlabTry(LINK_LONG_TYPE* const pThis, void* pVoid, const MrProtocolData::MrProtData*, int32_t);
#endif

//  -----------------------------------------------------------------
//  Oversampling - recalculate kspace
//  -----------------------------------------------------------------
double _asl_SliceOSSetValue(LINK_DOUBLE_TYPE* const pThis, double newVal, int32_t int32_t);

bool _asl_SliceOSGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t /*pos*/);

double _asl_PhaseOSSetValue(LINK_DOUBLE_TYPE* const pThis, double newVal, int32_t pos);

bool _asl_PhaseOSGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t /*pos*/);

#ifdef SUPPORT_iPAT_TGSE

//  -----------------------------------------------------------------
//  Phase FOV - recalculate kspace/epifactor
//  -----------------------------------------------------------------
double   _asl_PhaseFOVSetValue(LINK_DOUBLE_TYPE* const pThis, double newVal, int32_t pos);
bool     _asl_PhaseFOVGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t /*pos*/);

//  -----------------------------------------------------------------
//  Phase and Slice Partial Fourier - recalculate kspace/epifactor
//  -----------------------------------------------------------------
unsigned _asl_PhasePartialFourier_SetValue(MrUILinkSelection<unsigned>* const pThis, unsigned newVal, int32_t pos);
unsigned _asl_SlicePartialFourier_SetValue(MrUILinkSelection<unsigned>* const pThis, unsigned newVal, int32_t pos);

//  -----------------------------------------------------------------
//  iPAT - recalculate kspace/turbofactor/epifactor
//  -----------------------------------------------------------------
unsigned _asl_PATMode_SetValue(MrUILinkSelection<unsigned>* const pThis, unsigned newVal, int32_t pos);
int32_t     _asl_PATAccelFactorPE_SetValue(MrUILinkLimited<int32_t>* const pThis, int32_t newVal, int32_t pos);
int32_t     _asl_PATAccelFactor3D_SetValue(MrUILinkLimited<int32_t>* const pThis, int32_t newVal, int32_t pos);

bool fCoilCombineGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t lIndex);

// modified CAIPI handler
bool _asl_fUILinkPatShift3DGetLimits(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t pos);

#endif // SUPPORT_iPAT_TGSE

//  -----------------------------------------------------------------
//  BaseResolution - recalculate kspace/epifactor
//  -----------------------------------------------------------------
int32_t _asl_BaseResolutionSetValue(LINK_LONG_TYPE* const pThis, int32_t newVal, int32_t pos);

bool _asl_BaseResolutionGetLimits(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t /*pos*/);

bool _asl_BaseResolutionTry( LINK_LONG_TYPE* const pThis, void* /*pVoid*/, const MrProtocolData::MrProtData* /*pOrigProt*/, int32_t /*lIndex*/);

double _asl_PhaseResSetValue(LINK_DOUBLE_TYPE* const pThis, double newVal, int32_t pos);

bool _asl_PhaseResGetLimits(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t /*pos*/);

// ------------------------------------------------------------------------------
// Function    : _solve  ....
// ------------------------------------------------------------------------------
//
// Solve handling for lines to measure and epi factor
//
// ------------------------------------------------------------------------------
static unsigned _solve_LinesToMeasure_EPIFactor_SelectionConflict (MrUILinkSelection<unsigned> * const pThis, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);

static unsigned _solve_LinesToMeasure_EPIFactorConflict (MrUILinkLimited<int32_t> * const pThis, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);

static unsigned _solve_EPIFactor_LinesToMeasureConflict (MrUILinkLimited<int32_t> * const pThis, char **arg_list, const void *pVoid, const MrProtocolData::MrProtData*  pOrigProt, int32_t lIndex);
   
// ------------------------------------------------------------------------------
// Function    : AslMode_GetToolTip
// ------------------------------------------------------------------------------
//
// remove irrevalent tooltip
//
// ------------------------------------------------------------------------------
unsigned _asl_AslModeGetToolTipId(LINK_SELECTION_TYPE* const /*pThis*/, char * /*arg_list*/[], int32_t /*pos*/);

//unsigned _asl_AslModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned val, int32_t pos);

//unsigned _asl_AslModeSolve(LINK_SELECTION_TYPE* const pThis, char* arg_list[], const void*  pAddMem, const MrProtocolData::MrProtData*  pOrig, int32_t lPos);

//bool fUILinkPTSatSelGetOptions(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t);
//static void _filterPatSelOptions4Asl(const MrUILinkBase* const pThis, std::vector<unsigned>& rOptionVector);
//unsigned int _aslSeq2UI(SEQ::AslMode mode);
  // modified to enable both directions for PCASL
//bool _checkPTSatIsForAsl(unsigned int ptsatMode, unsigned int aslMode);

bool _asl_AslSuppressionModeIsAvailable(LINK_SELECTION_TYPE* const pThis, int32_t /*pos*/);

#endif   // of #ifdef WIN32

} // end of namespace Tgse_asl_UINS
} // end of namespace SEQ_NAMESPACE
