//--------------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2018  All Rights Reserved.  Confidential
//--------------------------------------------------------------------------------

// *------------------------------------------------------------------ *
// *  wip parameters for seq-special card                              *
// *    --> The enums define positions in the UI and in the proptocol   *
// *        (e.g. in WipMemBlock.adFree[])                              *  
// * ------------------------------------------------------------------ *
enum WIPUIPositions {
	POS_dCompSpoilPercentRead = 0,
	POS_dCompSpoilPercentSlice,
	POS_bCompSpoilDecay,
	POS_dECTau
};

#ifdef WIN32

#include "MrImagingFW/WIPParameterTool/WIPParameterTool.h"


//Pre-Declarations
unsigned _solveDoubleParamConflict(MrUILinkLimited<double>* const _this, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex);
unsigned _solveBoolParamConflict(MrUILinkSelection<bool>* const _this, char** arg_list, const void* pVoid, const MrProtocolData::MrProtData* pOrigProt, int32_t lIndex);


//  /////////////////////////////////////////////////////////////////
//  Register/Unregister handlers
//  /////////////////////////////////////////////////////////////////
bool initUIEddyCurrentComp(SeqLim &rSeqLim, WPT_NAMESPACE::WIPParameterTool& sWIPParamTool)
{
	//  Show special card?
	const bool bSeqSpecialCard = SysProperties::ReadSeqSettingGeneral("EDDYCURRCOMP_WIP_PARAMS_SPECIALCARD_VISIBLE", false, true);

	if (bSeqSpecialCard)
	{
		sWIPParamTool.createDoubleParameter(POS_dCompSpoilPercentRead, rSeqLim, "Comp. Read Spoil Moment", "%%", 0, 0.0, 150.0, 1.0, 100.0);
		sWIPParamTool.registerSolveHandler(POS_dCompSpoilPercentRead, _solveDoubleParamConflict);

		sWIPParamTool.createDoubleParameter(POS_dCompSpoilPercentSlice, rSeqLim, "Comp. Slice Spoil Moment", "%%", 0, 0.0, 150.0, 1.0, 100.0);
		sWIPParamTool.registerSolveHandler(POS_dCompSpoilPercentSlice, _solveDoubleParamConflict);

		sWIPParamTool.createBoolParameter(POS_bCompSpoilDecay, rSeqLim, "Comp: Consider Decay", false);
		sWIPParamTool.registerSolveHandler(POS_bCompSpoilDecay, _solveBoolParamConflict);

		sWIPParamTool.createDoubleParameter(POS_dECTau, rSeqLim, "EC Time Constant", "ms", 0, 0, 300, 1, 150);
		sWIPParamTool.registerSolveHandler(POS_dECTau, _solveDoubleParamConflict);	
	}

	return true;
}

#endif //  WIN32



