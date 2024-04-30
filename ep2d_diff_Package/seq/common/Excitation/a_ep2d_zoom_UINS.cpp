//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 2012  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4\pkg\MrServers\MrImaging\seq\common\Excitation\a_ep2d_zoom_UINS.cpp
//     Version: \main\20
//      Author: pfeujodj
//
//        Lang: C++
//
//
//
///  \file   a_ep2d_zoom_UINS.cpp
///  \brief  File implementing the UI class for EP2D with 2D RF pulses
///          - basic support functions: fEPIRegisterZoomHandlers(), initializeZoomUI()
///
//    -----------------------------------------------------------------------------

#ifdef WIN32

// --------------------------------------------------------------------------
// General Includes
// --------------------------------------------------------------------------
#include "MrImaging/seq/a_ep2d.h" // getSeq

#include "MrImaging/seq/common/Excitation/a_ep2d_zoom_UINS.h"
#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"
#include "MrImagingFW/libSeqSysProp/SysProperties.h" // ReadRegistry
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrMeasSrv/SeqIF/Sequence/Sequence.h"
#include "MrProtSrv/Domain/CoreNative/MrPTXData.h" // PTXVolume
#include "MrProtSrv/Domain/CoreNative/MrPat.h"     // iPAT
#include "MrProtSrv/Domain/MrProtData/MrProt/MrSliceGroup.h"

#define PRODUCT_DEF_NUM_VOL_B1SHIM 2 // max. number of volumes
#define PRODUCT_DEF_NUM_VOL_OPTIM 1  // max. number of volumes
#define PRODUCT_DEF_IND_VOL_OPTIM 0  // fixed index in "pTX Volumes" array
#define PRODUCT_DEF_IND_PTX_PULSE 0  // fixed index in "pTX Pulses" array

#define DEF_SIDELOBE_DISTANCE_mm 250

using namespace std;

// ===========================================================================
//
//      Ep2d_zoom_UINS::fEPIRegisterZoomHandlers
///
/// \brief  register UI handler for ZOOM_2DRF option
///
/// \return void
///
// ===========================================================================
void Ep2d_zoom_UINS::fEPIRegisterZoomHandlers(SeqLim& rSeqLim)
{
    if (MrUILinkArray* pArray = _search<MrUILinkArray>(rSeqLim, MR_TAG_PTX_VOLUME_ARRAY))
    {
        pArray->registerCanInsertHandler(Ep2d_zoom_UINS::fUILinkPTXVolArrayCanInsertNew);
        pArray->registerCanEraseHandler(Ep2d_zoom_UINS::fUILinkPTXVolArrayCanEraseNew);
    }

    // tooltip to show excitation RF pulse info
    UI_ELEMENT_DOUBLE m_PhaseFOV;
    m_PhaseFOV.registerToolTipHandler(rSeqLim, MR_TAG_PHASE_FOV, Ep2d_zoom_UINS::fUILinkPhaseFOV_GetToolTip);
}

// ===========================================================================
//
//      Ep2d_zoom_UINS::initializeZoomUI
///
/// \brief  initialize UI functions and members for ZOOM_2DRF option
///
/// \return void
///
// ===========================================================================
NLS_STATUS Ep2d_zoom_UINS::initializeZoomUI(MrProt& rMrProt, SeqLim& rSeqLim)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    calcPTXVolFromSliceGroup(rMrProt);
    calcPTXMPRFromSliceGroup(rMrProt);
    setPTXCalculation(rMrProt);

    // if( !calcOvSatFromSliceGroup(rMrProt, rSeqLim) )
    //    lStatus = MRI_SEQ_SEQU_ERROR;

    return (lStatus);
}

//	calculate PTXVolume by fitting it to the first slice group (that has to be in the protocol)
void Ep2d_zoom_UINS::calcPTXVolFromSliceGroup(MrProt& rMrProt)
{
    // multiple slice groups are NOT considered!
    const int32_t lSGPos = 0;
    const int32_t nIndex = PRODUCT_DEF_IND_VOL_OPTIM;

    if ((rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED) && (rMrProt.getsPTXData().getasPTXVolume().size() > nIndex)
        && (rMrProt.getsPTXData().getasPTXVolume()[nIndex].getlVolProperty() == MrProtocolData::PTxVolProp_Optimization))
    {
        Slice ptxSlice(rMrProt.getsPTXData().getasPTXVolume()[nIndex].getsSliceData());
        Slice aPtxSlice(ptxSlice);

        //-----------------------------------------------------------------
        //  set PTX.vol. to the slice group, except PhaseFOV
        //-----------------------------------------------------------------
        ptxSlice.setdThickness(ceil(rMrProt.sliceGroupList()[lSGPos].thickness() + (rMrProt.getsSliceArray().getlSize() - 1) * rMrProt.sliceGroupList()[lSGPos].distance()));
        ////ptxSlice.setdPhaseFOV   ( ceil(rMrProt.sliceGroupList()[lSGPos].phaseFOV())  ); // do NOT copy: this parameter remains flexible
        ptxSlice.setdReadoutFOV(ceil(rMrProt.sliceGroupList()[lSGPos].readoutFOV()));
        ptxSlice.setdInPlaneRot(rMrProt.sliceGroupList()[lSGPos].rotationAngle());
        aPtxSlice.position(rMrProt.sliceGroupList()[lSGPos].position());
        aPtxSlice.normal(rMrProt.sliceGroupList()[lSGPos].normal());

        // do not use this volume for Adjust framework calculations
        rMrProt.getsPTXData().getasPTXVolume()[nIndex].setlIsUsedForCalculation(MrProtocolData::PTxVolCalcMode_Off);
    }

    return;
}

//	calculate PTXMPR slice group by fitting it to the slice group and PhaseFOV of the PTXVolume
void Ep2d_zoom_UINS::calcPTXMPRFromSliceGroup(MrProt& rMrProt)
{
    const int32_t nIndex = PRODUCT_DEF_IND_VOL_OPTIM;

    if ((rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED) && (rMrProt.getsPTXData().getasPTXVolume().size() > nIndex)
        && (rMrProt.getsPTXData().getasPTXVolume()[nIndex].getlVolProperty() == MrProtocolData::PTxVolProp_Optimization))
    {
        // retrieve FoEOuterBorder from the PhaseFOV of the *first* PTXVolume with the property 'Optimization'
        double dFoEOuterBorder = rMrProt.getsPTXData().getasPTXVolume()[nIndex].getsSliceData().getdPhaseFOV(); // [mm]

        // initialize PTXMPR slice group from SliceArray and GroupArray
        rMrProt.getsPTXData().getsPTXMPRSliceArray().copyFrom(rMrProt.getsSliceArray());
        rMrProt.getsPTXData().getsPTXMPRGroupArray().copyFrom(rMrProt.getsGroupArray());

        // adapt PhaseFOV of PTXMPR slice group to the FoEOuterBorder from Optimization PTXVolume
        for (int32_t iIndSlc = 0; iIndSlc < rMrProt.getsPTXData().getsPTXMPRSliceArray().getlSize(); iIndSlc++)
        {
            rMrProt.getsPTXData().getsPTXMPRSliceArray().getasSlice()[iIndSlc].setdPhaseFOV(dFoEOuterBorder);
        }

        // make PTXMPR slice group valid: flag defines whether the original slice group or the new PTX-MPR slice group should be used
        rMrProt.getsPTXData().setuiPTXMPRSliceGroupValid(true);
    }
    else
    {
        // reset PTXMPR slice group
        rMrProt.getsPTXData().getsPTXMPRSliceArray().clear();
        rMrProt.getsPTXData().getsPTXMPRGroupArray().clear();
        rMrProt.getsPTXData().setuiPTXMPRSliceGroupValid(false);
    }

    return;
}

//	calculate Outer-volume Sats (ovSat) from slice group
bool Ep2d_zoom_UINS::calcOvSatFromSliceGroup(MrProt& rMrProt, SeqLim& rSeqLim)
{
    const double dMinThkOvSat   = 10.5;  // [mm] minimum thk: see RSat.cpp, below no OvSat will be set
    const double dMaxThkOvSat   = 110.0; // [mm] maximum thk
    int32_t      lNumberOfOvSat = 2;

    if (rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED)
    {
        if (rMrProt.getsRSatArray().getlSize() < lNumberOfOvSat)
        {
            return true; // nothing to do
        }

#ifdef DEBUG
        SEQ_TRACE_INFO.print("calculate new OvSat from SliceGroup.");
#endif

        double dRelExcitationSize = RF2DArbGenerator::Default2DExciteRelativeExcitationSize;
        double dThkFactor         = 1.0;

        if (SysProperties::getNominalBZero() > 2.8) // bIs3TeslaOrMore
        {
            dThkFactor = 1.15 * 1.25;
        }

        double dReducedFoE     = rMrProt.sliceSeries().aFront().phaseFOV() * dRelExcitationSize;
        double dFoEOuterBorder = -1; // [mm]

        // retrieve FoEOuterBorder from the PhaseFOV of the *first* PTXVolume with the property 'Optimization'
        int32_t lPTXVolumeSize = rMrProt.getsPTXData().getasPTXVolume().size();

        if (lPTXVolumeSize > 0)
        {
            for (long iVol = 0; ((iVol < lPTXVolumeSize) && (dFoEOuterBorder < 0)); iVol++)
            {
                if (rMrProt.getsPTXData().getasPTXVolume()[iVol].getlVolProperty() == MrProtocolData::PTxVolProp_Optimization)
                {
                    dFoEOuterBorder = rMrProt.getsPTXData().getasPTXVolume()[iVol].getsSliceData().getdPhaseFOV(); // [mm]
                }
            }
        }

        if (dFoEOuterBorder <= 0.0)
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_INFO.print("could not retrieve FoEOuterBorder.");
            }
        }

        // Re-calc slice thk RSat implementation: effective thickness is extended by dThkFactor
        double dOvSatNominalThickness = (dFoEOuterBorder - dReducedFoE) * 0.5;

        if (dOvSatNominalThickness * dThkFactor <= dMinThkOvSat)
        {
            if (!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_INFO.print("OvSatThickness = %.1f <= %.1f - No pulse is set.", dOvSatNominalThickness * dThkFactor, dMinThkOvSat);
            }
            return true;
        }

        // RSat is limited to dMaxThkOvSat (110 mm)
        if (dOvSatNominalThickness * dThkFactor > dMaxThkOvSat)
        {
            lNumberOfOvSat += 2;
        }
        dOvSatNominalThickness = std::min(dMaxThkOvSat, dOvSatNominalThickness * dThkFactor) / dThkFactor;

        // cover the outer volume up to the outer border,
        double dOvSatShift = (dReducedFoE + dOvSatNominalThickness) * 0.5;
        // optionally: shift in addition by the edge width (dXPhaseResolution) not to touch the inner volume excitation with the RSats
        ////DISABLED: dOvSatShift += dXPhaseResolution;

        double dOvSatShift2 = dOvSatShift + dOvSatNominalThickness * 0.8; // 20% overlap of Sats

#ifdef DEBUG
        double dXPhaseResolution = RF2DArbGenerator::Default2DExciteResolutionPhase;
        SEQ_TRACE_INFO.print(
            "FoEOuterBorder/ReducedFoE/XPhaseResolution/OvSatNominalThickness/OvSatShift = %f / %f / %f / %f / %f",
            dFoEOuterBorder,
            dReducedFoE,
            dXPhaseResolution,
            dOvSatNominalThickness,
            dOvSatShift);
#endif

        Slice _slice(rMrProt.getsSliceArray().getasSlice()[0]); // align to first slice in SliceArray

        VectorPat<double> e_phase, e_readout;
        if (!_slice.orientation(e_phase, e_readout))
        {
            SEQ_TRACE_INFO.print("calculation of slice orientation failed.");
            return false;
        }

        int32_t lShiftSign = +1; // shift in positive direction
        for (int32_t lIndexOvSat = 0; lIndexOvSat < lNumberOfOvSat; lIndexOvSat++)
        {
            double dTmpShift;
            if (lIndexOvSat < 2)
            {
                dTmpShift = dOvSatShift;
            }
            else
            {
                dTmpShift = dOvSatShift2;
            }

            // - only set OvSat for RSat's that are switch on in MrProt
            // - do not set OvSat if thk is too small
            if (rMrProt.getsRSatArray().getlSize() > lIndexOvSat)
            {
                RSat _rsat(rMrProt.getsRSatArray().getasElm()[lIndexOvSat]);

                VECTOR_PAT_DBL sNewPosition = e_phase;
                sNewPosition *= lShiftSign * dTmpShift; // multiply normal vector with shift
                sNewPosition += _slice.position();      // shift from original position

                _rsat.normal(e_phase);
                _rsat.position(sNewPosition);
                _rsat.thickness(dOvSatNominalThickness * dThkFactor);

#ifdef DEBUG
                SEQ_TRACE_INFO.print("OvSat%d - OvSatThickness/Shift = %f / %f.", lIndexOvSat, dOvSatNominalThickness * dThkFactor, lShiftSign * dTmpShift);
#endif

                // next RSat located in inverse direction
                if (lShiftSign > 0)
                    lShiftSign = -1;
                else
                    lShiftSign = +1;
            }
        }
    }

    return true;
}

//  enables PTX pulses
void Ep2d_zoom_UINS::setPTXCalculation(MrProt& rMrProt)
{
    // pointer for MrProtFacade for easier protocol queries
    MrProtFacade protFacade(rMrProt);

    const long nIndex = PRODUCT_DEF_IND_PTX_PULSE;

    // set default
    rMrProt.getsTXSPEC().getasNucleusInfo()[0].setRFPulseExcitationMode(ePulseExcitationDefault);

    // for SliceAdj with pTX we need to reactivate the property
    if (protFacade.isSliceAdj())
    {
        if (rMrProt.getsAdjData().getuiAdjSliceBySlicePtx())
        {
            rMrProt.getsTXSPEC().getasNucleusInfo()[0].setRFPulseExcitationMode(ePulseExcitationPTX);
        }
    }
    // default: no B1 nor B0 mapping
    rMrProt.getsAdjData().setuiAdjB0AcqMode(SEQ::ADJB0ACQ_NONE);
    rMrProt.getsAdjData().setuiAdjB1AcqMode(SEQ::ADJB1ACQ_NONE);

    if ((rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED) && (rMrProt.getsTXSPEC().getaPTXRFPulse().size() > nIndex))
    {
        rMrProt.getsTXSPEC().getasNucleusInfo()[0].setRFPulseExcitationMode(ePulseExcitationDynPulse1TX);

        // PulseAcceleration larger than 1.0 enables pTX pulses: B1 mapping, optional B0
        if ((rMrProt.getsTXSPEC().getaPTXRFPulse()[nIndex].getlTrajectoryType() == MrProtocolData::PTXTrajectoryType_EPI_1D)
            && (rMrProt.getsTXSPEC().getaPTXRFPulse()[nIndex].getdPulseAcceleration() > 1.0))
        {
            rMrProt.getsTXSPEC().getasNucleusInfo()[0].setRFPulseExcitationMode(ePulseExcitationPTX);
            rMrProt.getsAdjData().setuiAdjB1AcqMode(SEQ::ADJB1ACQ_STANDARD_SINGLE); // measure B1 map

            if (rMrProt.getsTXSPEC().getaPTXRFPulse()[nIndex].getlB0CorrectionType() != MrProtocolData::PTXB0CorrectionType_Off)
            {
                rMrProt.getsAdjData().setuiAdjB0AcqMode(SEQ::ADJB0ACQ_STANDARD); // measure B0 map
            }
        }

        // ExternalIni enables pTX pulses: no B0/B1 mapping necessary
        if (rMrProt.getsTXSPEC().getaPTXRFPulse()[nIndex].getlTrajectoryType() == MrProtocolData::PTXTrajectoryType_ExternalIni)
        {
            rMrProt.getsTXSPEC().getasNucleusInfo()[0].setRFPulseExcitationMode(ePulseExcitationPTX);
        }
    }
}

void Ep2d_zoom_UINS::_refreshPTXVol(MrUILinkBase* const pThis)
{
    try
    {
        MrProt rMrProt(pThis->prot());
        calcPTXVolFromSliceGroup(rMrProt);
        calcPTXMPRFromSliceGroup(rMrProt);
        setPTXCalculation(rMrProt);
    }
    catch (...)
    {
        std::cerr << "_refreshPTXVol: Caught unknown exception" << std::endl;
    }
}

//  overloaded Excitation handler: MrServers\MrProtSrv\MrProtocol\UILink\MrUILinkPrepPulses.cpp
//  -----------------------------------------------------------------
//  Excitation Pulse
//
void Ep2d_zoom_UINS::_fUILinkInsertOptimizationVolume(MrUILinkBase* const pThis)
{
    MrUILinkArray* pArray = _search<MrUILinkArray>(pThis, MR_TAG_PTX_VOLUME_ARRAY);
    if (pArray && pArray->isAvailable(0))
    {
        // insert new volume at default position
        int32_t nIndex = PRODUCT_DEF_IND_VOL_OPTIM;

        pArray->insert(nIndex, 0);

        if (pArray->isAvailable(nIndex))
        {
            // set default phaseFOV
            pThis->prot().getsPTXData().getasPTXVolume()[nIndex].getsSliceData().setdPhaseFOV(DEF_SIDELOBE_DISTANCE_mm);
            // set volume property
            pThis->prot().getsPTXData().getasPTXVolume()[nIndex].setlVolProperty(MrProtocolData::PTxVolProp_Optimization);
            // set PTxVolCalcMode off
            pThis->prot().getsPTXData().getasPTXVolume()[nIndex].setlIsUsedForCalculation(MrProtocolData::PTxVolCalcMode_Off);
        }
    }
}

unsigned Ep2d_zoom_UINS::fUILinkExcitationPulseSetValueNew(LINK_SELECTION_TYPE* const pThis, unsigned newVal, int32_t pos)
{
    // call original handler
    newVal = (*(getSeq(pThis)->getUI()->m_ExitPulse.getOrigSetValueHandler()))(pThis, newVal, pos);

    MrProt rMrProt(pThis->prot());

    Ep2d_zoom_UINS::setPTXCalculation(rMrProt);

    if (newVal == MRI_STD_EXCITATION_ZOOMED)
    {
        Ep2d_zoom_UINS::_fUILinkInsertOptimizationVolume(pThis);

#ifdef SUPPORT_iPAT_a_ep2d
        // switch iPAT off upon selection of ZOOMit (PAT still remains switchable)
        LINK_SELECTION_TYPE* pPATMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_PAT_MODE, 0, MrUILinkBase::SEARCH_MODE::SEARCH_AVAILABLE);
        if (pPATMode && (pPATMode->value(0) != MRI_STD_PAT_MODE_NONE))
        {
            pPATMode->value(MRI_STD_PAT_MODE_NONE, 0);
            pThis->addDependentParamPtr(pPATMode, 0);
        }
#endif
    }

    return newVal;
}

//  overloaded PTXVolume handler: MrServers\MrProtSrv\MrProtocol\UILink\MrUILinkPTXVolume.cpp
//  -----------------------------------------------------------------

//  /////////////////////////////////////////////////////////////////
//  Volume Property
//  /////////////////////////////////////////////////////////////////

bool Ep2d_zoom_UINS::fUILinkPTXVolPropGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos)
{
    rulVerify = LINK_BOOL_TYPE::VERIFY_ON;
    rOptionVector.clear();
    rOptionVector.reserve(3);

    LINK_SELECTION_TYPE* pB1ShimMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_ADJ_B1_SHIM_MODE);
    LINK_SELECTION_TYPE* pExitPulse  = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_EXCIT_PULSE);

    if (pExitPulse && pExitPulse->isAvailable(pos) && (pExitPulse->value(pos) == MRI_STD_EXCITATION_ZOOMED))
    {
        // allow PTXVolType Optimization option only for a single volume
        if ((pos == PRODUCT_DEF_IND_VOL_OPTIM) && pThis->seqLimits().getpTxVolProperty().hasOption(MrProtocolData::PTxVolProp_Optimization))
        {
            rOptionVector.push_back(MRI_STD_PROP_OPTIMIZATION);
        }
        else if (
            pThis->seqLimits().getpTxVolProperty().hasOption(MrProtocolData::PTxVolProp_B1Shim) && pB1ShimMode && pB1ShimMode->isAvailable(pos)
            && (pB1ShimMode->value(pos) == MRI_STD_B1_SHIM_MODE_VOL_SEL))
        {
            rOptionVector.push_back(MRI_STD_PROP_B1_SHIM);
        }
    }
    else
    {
        if (pThis->seqLimits().getpTxVolProperty().hasOption(MrProtocolData::PTxVolProp_B1Shim) && pB1ShimMode && pB1ShimMode->isAvailable(pos)
            && (pB1ShimMode->value(pos) == MRI_STD_B1_SHIM_MODE_VOL_SEL))
        {
            rOptionVector.push_back(MRI_STD_PROP_B1_SHIM);
        }
    }

    return rOptionVector.size() > 0;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolArrayCanInsertNew(MrUILinkArray* const pThis, int32_t newPos, int32_t pos)
{
    const int32_t i32NavArraySize    = pThis->prot().getsPTXData().getasPTXVolume().size();
    const int32_t i32NavArrayMaxSize = static_cast<int32_t>(pThis->prot().getsPTXData().getasPTXVolume().maxSize());

    if ((i32NavArraySize >= i32NavArrayMaxSize) || (newPos < 0) || (newPos > i32NavArraySize))
    {
        return false;
    }

    // count PTXVolTypes in current volumes
    int32_t i32ArraySize_B1Shim = 0;
    int32_t i32ArraySize_Optim  = 0;

    LINK_LONG_TYPE* pPTXVolCount = _search<LINK_LONG_TYPE>(pThis, MR_TAG_PTX_VOLUME_COUNT);
    if (pPTXVolCount && pPTXVolCount->isAvailable(0) && pPTXVolCount->value(0) > 0)
    {
        MrUILinkArray* pArray = _search<MrUILinkArray>(pThis, MR_TAG_PTX_VOLUME_ARRAY);
        if (pArray && pArray->isAvailable(0))
        {
            LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
            for (int32_t nIndex = pPTXVolCount->value(0) - 1; nIndex >= 0; nIndex--)
            {
                if (pPTXVolType && pPTXVolType->isAvailable(nIndex))
                {
                    if ((pPTXVolType->value(nIndex) == MRI_STD_PROP_B1_SHIM))
                    {
                        i32ArraySize_B1Shim++;
                    }
                    else if ((pPTXVolType->value(nIndex) == MRI_STD_PROP_OPTIMIZATION))
                    {
                        i32ArraySize_Optim++;
                    }
                }
            }
        }
    }

    if ((i32ArraySize_B1Shim < PRODUCT_DEF_NUM_VOL_B1SHIM) && pThis->seqLimits().getpTxVolProperty().hasOption(MrProtocolData::PTxVolProp_B1Shim))
    {
        LINK_SELECTION_TYPE* pB1ShimMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_ADJ_B1_SHIM_MODE);
        if (pB1ShimMode && (pB1ShimMode->isAvailable(pos)) && (pB1ShimMode->value(pos) == MRI_STD_B1_SHIM_MODE_VOL_SEL))
        {
            return true;
        }
    }

    if ((i32ArraySize_Optim < PRODUCT_DEF_NUM_VOL_OPTIM) && pThis->seqLimits().getpTxVolProperty().hasOption(MrProtocolData::PTxVolProp_Optimization))
    {
        LINK_SELECTION_TYPE* pExitPulse = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_EXCIT_PULSE);
        if (pExitPulse && (pExitPulse->isAvailable(pos)) && (pExitPulse->value(pos) == MRI_STD_EXCITATION_ZOOMED))
        {
            return true;
        }
    }

    return false;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolArrayCanEraseNew(MrUILinkArray* const pThis, int32_t doomed, int32_t pos)
{
    // do not allow to erase the optimization PTXVolume
    LINK_SELECTION_TYPE* pExitPulse = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_EXCIT_PULSE);
    if (pExitPulse && pExitPulse->isAvailable(pos) && (pExitPulse->value(pos) == MRI_STD_EXCITATION_ZOOMED) && doomed == PRODUCT_DEF_IND_VOL_OPTIM)
    {
        LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
        if (pPTXVolType && pPTXVolType->isAvailable(doomed) && (pPTXVolType->value(doomed) == MRI_STD_PROP_OPTIMIZATION))
        {
            return false;
        }
    }

    return ((doomed >= 0) && (doomed < (int32_t)pThis->prot().getsPTXData().getasPTXVolume().size()));
}

// make PTXVolume phase valid for binary search or make field not editable in case of rotated trajectory
bool Ep2d_zoom_UINS::fUILinkPTXVolPDimGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPDim.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
#ifdef ZOOM_EXTENDED
        return false;
#else
        // enable binary search
        rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;
#endif
    }

    return bRet;
}

// make PTXVolume fields NOT editable for type optimization
bool Ep2d_zoom_UINS::fUILinkPTXVolRDimGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolRDim.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolSDimGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolSDim.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolRotGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolRot.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPosSagGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPosSag.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPosSag_SBCSGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPosSagSBCS.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPosCorGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPosCor.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPosCor_SBCSGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPosCorSBCS.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPosTraGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPosTra.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPosTra_SBCSGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPosTraSBCS.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolOriAlphaGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolOriAlpha.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolOriBetaGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolOriBeta.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPosGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPos.getOrigGetOptionsHandler()))(pThis, rOptionVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPos_SBCSGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPosSBCS.getOrigGetOptionsHandler()))(pThis, rOptionVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolOriGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolOri.getOrigGetOptionsHandler()))(pThis, rOptionVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolOriHistoryGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolOriHistory.getOrigGetOptionsHandler()))(pThis, rOptionVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

#ifdef ZOOM_EXTENDED
bool Ep2d_zoom_UINS::fUILinkPTXVolVisibilityGetOptionsNew(LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, uint32_t& rulVerify, int32_t pos)
{
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolVisibility.getOrigGetOptionsHandler()))(pThis, rOptionVector, rulVerify, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}
#endif // ZOOM_EXTENDED

// update PTXVolume fields by calling calcPTXVolFromSliceGroup()
double Ep2d_zoom_UINS::fUILinkPTXVolRotGetValueNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    _refreshPTXVol(pThis);

    return (*(getSeq(pThis)->getUI()->m_PTXVolRot.getOrigGetValueHandler()))(pThis, pos);
}

double Ep2d_zoom_UINS::fUILinkPTXVolPDimGetValueNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    return (*(getSeq(pThis)->getUI()->m_PTXVolPDim.getOrigGetValueHandler()))(pThis, pos);
}

double Ep2d_zoom_UINS::fUILinkPTXVolRDimGetValueNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    return (*(getSeq(pThis)->getUI()->m_PTXVolRDim.getOrigGetValueHandler()))(pThis, pos);
}

double Ep2d_zoom_UINS::fUILinkPTXVolSDimGetValueNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    return (*(getSeq(pThis)->getUI()->m_PTXVolSDim.getOrigGetValueHandler()))(pThis, pos);
}

// update PTXVolume fields by calling _refreshPTXVol() for GSP rotation parameter
GradDirPat Ep2d_zoom_UINS::fUILinkPTXVolGSPRotGetValueNew(LINK_GRAD_DIR_TYPE* const pThis, int32_t lPos)
{
    _refreshPTXVol(pThis);

    return (*(getSeq(pThis)->getUI()->m_PTXVolRotGSP.getOrigGetValueHandler()))(pThis, lPos);
}

// update PTXVolume fields by calling _refreshPTXVol() for GSP position parameter
VectorPat<double> Ep2d_zoom_UINS::fUILinkPTXVolGSPPosGetValueNew(LINK_VECTOR_TYPE* const pThis, int32_t pos)
{
    _refreshPTXVol(pThis);

    return (*(getSeq(pThis)->getUI()->m_PTXVolPosGSP.getOrigGetValueHandler()))(pThis, pos);
}

// update PTXPulse fields / make them non-zero-sized with EXCITATION_ZOOMED
int32_t Ep2d_zoom_UINS::fUILinkPTXPulseArraySizeNew(MrUILinkArray* const pThis, int32_t pos)
{
    LINK_SELECTION_TYPE* pExitPulse = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_EXCIT_PULSE);
    if (pExitPulse && pExitPulse->isAvailable(pos) && (pExitPulse->value(pos) == MRI_STD_EXCITATION_ZOOMED))
    {
        return pThis->prot().getsTXSPEC().getaPTXRFPulse().size();
    }

    return 0;
}

bool Ep2d_zoom_UINS::fUILinkPTXPulseIsAvailableNew(MrUILinkMap* const pThis, int32_t pos)
{
    LINK_SELECTION_TYPE* pExitPulse = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_EXCIT_PULSE);
    if (pExitPulse && pExitPulse->isAvailable(pos) && (pExitPulse->value(pos) == MRI_STD_EXCITATION_ZOOMED))
    {
        return (SysProperties::isPTxSystem() && pThis->seqLimits().getPTXPulses().isAvailable() && (int32_t)pThis->prot().getsTXSPEC().getaPTXRFPulse().size() > pos);
    }

    return false;
}

double Ep2d_zoom_UINS::fUILinkPTXPulseTxAccSetValueNew(LINK_DOUBLE_TYPE* const pThis, double newVal, int32_t pos)
{
    MrProt rMrProt(pThis->prot());

    // call original handler
    newVal = (*(getSeq(pThis)->getUI()->m_PTXPulseTxAcc.getOrigSetValueHandler()))(pThis, newVal, pos);

    if ((newVal > 1.0) && rMrProt.preparationPulses().getucFatSatMode() != SEQ::FAT_SAT_STRONG)
    {
        LINK_SELECTION_TYPE* _FatSatMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_FAT_SAT_MODE);
        if (_FatSatMode && _FatSatMode->isEditable(0) && pThis->seqLimits().getFatSatMode().hasOption(SEQ::FAT_SAT_STRONG))
        {
            _FatSatMode->value(MRI_STD_STRONG, 0);
            pThis->addDependentParamPtr(_FatSatMode, 0);
        }
    }

    return newVal;
}

bool Ep2d_zoom_UINS::fUILinkPTXPulseTxAccGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXPulseTxAcc.getOrigGetLimitsHandler()))(pThis, rLimitVector, rulVerify, pos);

    // enable binary search
    rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;

    return bRet;
}

#ifdef ZOOM_EXTENDED
bool Ep2d_zoom_UINS::fUILinkPTXPulseTxAccIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    return false;
}
#endif

// make PTXPulse fields NOT editable
bool Ep2d_zoom_UINS::fUILinkPTXPulseFlipAngleGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    rulVerify = LINK_DOUBLE_TYPE::VERIFY_OFF;
    rLimitVector.clear();

    return false;
}

bool Ep2d_zoom_UINS::fUILinkPTXPulsePhaseFoEGetLimitsNew(LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    rulVerify = LINK_DOUBLE_TYPE::VERIFY_OFF;
    rLimitVector.clear();

    return false;
}

bool Ep2d_zoom_UINS::fUILinkPTXPulsePhaseMatrixSizeGetLimitsNew(LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t pos)
{
    rulVerify = LINK_DOUBLE_TYPE::VERIFY_OFF;
    rLimitVector.clear();

    return false;
}

#ifdef ZOOM_EXTENDED
// the following fields are not required for the rotated trajectory -> make them not available
bool Ep2d_zoom_UINS::fUILinkPTXVolPDimIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPDim.getOrigIsAvailableHandler()))(pThis, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
} // ZOOM_EXTENDED

bool Ep2d_zoom_UINS::fUILinkPTXVolRotIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolRot.getOrigIsAvailableHandler()))(pThis, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolRDimIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolRDim.getOrigIsAvailableHandler()))(pThis, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolSDimIsAvailableNew(LINK_DOUBLE_TYPE* const pThis, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolSDim.getOrigIsAvailableHandler()))(pThis, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolPosIsAvailableNew(LINK_SELECTION_TYPE* const pThis, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolPos.getOrigIsAvailableHandler()))(pThis, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolOriIsAvailableNew(LINK_SELECTION_TYPE* const pThis, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolOri.getOrigIsAvailableHandler()))(pThis, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        return false;
    }

    return bRet;
}

bool Ep2d_zoom_UINS::fUILinkPTXVolVisibilityIsAvailableNew(LINK_SELECTION_TYPE* const pThis, int32_t pos)
{
    // call original handler
    bool bRet = (*(getSeq(pThis)->getUI()->m_PTXVolVisibility.getOrigIsAvailableHandler()))(pThis, pos);

    LINK_SELECTION_TYPE* pPTXVolType = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_PTX_VOLUME_ARRAY, MR_TAG_VOL_PROPERTY);
    if (pPTXVolType && pPTXVolType->isAvailable(pos) && (pPTXVolType->value(pos) == MRI_STD_PROP_OPTIMIZATION))
    {
        // set volume visibility to off
        pThis->prot().getsPTXData().getasPTXVolume()[pos].setlVolVisible(SEQ::OFF);

        // hide UI element
        return false;
    }

    return bRet;
}
#endif

unsigned Ep2d_zoom_UINS::fUILinkPhaseFOV_GetToolTip(LINK_DOUBLE_TYPE* const _this, char* arg_list[], int32_t /*pos*/)
{
    MrProt      rMrProt(_this->prot());
    static char tToolTip[1000];
    tToolTip[0] = '\0';

#ifdef DEBUG
    _this->sequence().prepareForScanTimeCalculation(rMrProt);

    // determine some parameters from the *first* PTXRFPulse (index 0)
    const size_t                      uiCurVol        = 0;
    const int32_t                     lPTXRFPulseSize = rMrProt.getsTXSPEC().getaPTXRFPulse().size();
    MrProtocolData::PTXTrajectoryType eExcType        = MrProtocolData::PTXTrajectoryType_EPI_1D;
    if (lPTXRFPulseSize > 0)
    {
        eExcType = rMrProt.getsTXSPEC().getaPTXRFPulse()[uiCurVol].getlTrajectoryType();
    }

    if ((rMrProt.getsTXSPEC().getucExcitMode() == SEQ::EXCITATION_ZOOMED) && (eExcType == MrProtocolData::PTXTrajectoryType_EPI_1D))
    {
        SBB2DPtx* pSBBExcite = dynamic_cast<SBB2DPtx*>(getSeq(_this)->m_EPIKernel.getExcitationPointer());
        if (pSBBExcite != nullptr)
        {
            const std::string strTmp = pSBBExcite->getToolTipInfo();
            sprintf(tToolTip, "ZOOM_2DRF Info:\n%s", strTmp.c_str());
        }
    }
#endif

    arg_list[0] = tToolTip;
    return MRI_STD_STRING;
}

#endif // of #ifdef WIN32
