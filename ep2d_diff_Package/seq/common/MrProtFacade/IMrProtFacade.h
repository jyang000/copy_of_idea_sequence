//-----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 2012 All Rights Reserved.  Confidential.
//-----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//    File: \n4_servers1\pkg\MrServers\MrImaging\seq\common\MrProtFacade\IMrProtFacade.h
// Version: \main\7
//  Author: Simon Bauer (z002yrac), Uvo Hoelscher (z0038c7r)
//    Date: 2015-02-11 10:13:28 +01:00
//
//    Lang: C++
//
// Descrip: 
//
// Classes: IMrProtFacade
//
//-----------------------------------------------------------------------------
#ifndef IMRPROTFACADE_H
#define IMRPROTFACADE_H 1

#include "MrMeasSrv/MeasSections/MeasCoilContext.h"
class IMrProtFacade
{
public:
    virtual ~IMrProtFacade(){};
    typedef std::shared_ptr<IMrProtFacade> Pointer;

    virtual bool isTwist() = 0;
    virtual bool isIcePat() = 0;
    virtual bool isTomEquidistantWith2OrMoreEchoes() = 0;
    virtual bool isAsymmetricEcho() = 0;
    virtual bool isInterpolation() = 0;
    virtual bool isImageBasedInterpolation() = 0;
    virtual bool isSwitchOffZeroFillingInICE() = 0;
    virtual bool isBookkeepingConditionForDFC() = 0;
    virtual bool isContinuousSat() = 0;
    virtual bool isZOOMit() = 0;
    virtual bool isAcceleratedZOOMit() = 0;
    virtual bool is3DiPAT() = 0;
    virtual bool isCAIPIRINHA() = 0;
    virtual bool isReducedMotionSens() = 0;
    virtual bool isElliptical() = 0;
    virtual bool isMinTTC() = 0;
    virtual bool isUserDefinedTTC() = 0;

    virtual bool isBreastCoil() = 0;
    virtual bool isBreastApplication() = 0;
    virtual bool isAnyCartesianRadialReorderingInPhaseEncodingPlane() = 0;
    virtual bool isGRASP() = 0;
    virtual bool isGradientReversal() = 0;
    virtual bool isRadialSelfGating() = 0;
    virtual bool isIceProgramRadialVIBESelfGatingAdd() = 0;

    // DIXON
    virtual bool isDixon() = 0;
	virtual	bool isFastDixon() = 0;
	virtual bool isConventionalDixon() = 0;
    virtual bool isOppAndInTwoPointDixon()              = 0;
    virtual bool isMinAndOppTwoPointDixon()              = 0;
    virtual bool isFlexTwoPointDixon() = 0;
    virtual bool isTwoPointDixon() = 0;
    virtual bool isMultiEchoDixon() = 0;
    virtual bool isOppAndInTwoPointScreeningDixon() = 0;
    virtual bool isTwoPointDixonWithoutAsymmetricEcho() = 0;
    virtual bool isTwoPointDixonWithFirstOppInPhases() = 0;
    virtual bool isDixonT2StarOrR2StarEvaluation() = 0;
    virtual bool isDixonFatOrWaterFractionEvaluation() = 0; 
    virtual bool isAnyMultiEchoDixonEvaluationOption() = 0;
    virtual bool isDixonFatOrWaterImage() = 0;

    // Diffusion
    virtual bool isGradientReversalDiffusion() = 0;
    virtual bool isBetterSliceProfileDiffusion() = 0;
    virtual bool iSDTI() = 0;
    virtual bool isIVIM() = 0;

    // SMS
    virtual bool isSliceAcceleration() = 0;

    // SliceAdj
    virtual bool isSliceAdj() = 0;

    virtual bool isAnyFatSat() = 0;
    virtual bool isAnySat() = 0;
    virtual bool isWaterExcitation() = 0;
    virtual bool isFastWaterExcitation() = 0;
    virtual bool isAnyWaterExcitation() = 0;
    virtual bool isQuickFatSat() = 0;
    virtual bool isSPAIRFatSat() = 0;
};

#endif



