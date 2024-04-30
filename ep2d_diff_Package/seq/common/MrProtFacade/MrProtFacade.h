//-----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 2012 All Rights Reserved.  Confidential.
//-----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//    File: \n4_servers1\pkg\MrServers\MrImaging\seq\common\MrProtFacade\MrProtFacade.h
// Version: \main\7
//  Author: Simon Bauer (z002yrac), Uvo Hoelscher (z0038c7r)
//    Date: 2015-02-11 10:15:37 +01:00
//
//    Lang: C++
//
// Descrip: 
//
// Classes: MrProtFacade
//
//-----------------------------------------------------------------------------
#ifndef MRPROTFACADE_H
#define MRPROTFACADE_H 1

#include "IMrProtFacade.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"


class MrProtFacade: public IMrProtFacade
{
public:
    MrProtFacade(MrProt& rMrProt);
    MrProtFacade(const MrProt& rMrProt);
    ~MrProtFacade();

    virtual bool isTwist();
    virtual bool isIcePat();
    virtual bool isTomEquidistantWith2OrMoreEchoes();
    virtual bool isAsymmetricEcho();
    virtual bool isInterpolation();
    virtual bool isImageBasedInterpolation();
    virtual bool isSwitchOffZeroFillingInICE();
    virtual bool isBookkeepingConditionForDFC();
    virtual bool isZOOMit();
    virtual bool isAcceleratedZOOMit();
    virtual bool isContinuousSat();
    virtual bool is3DiPAT();
    virtual bool isCAIPIRINHA();
    virtual bool isReducedMotionSens();
    virtual bool isElliptical();
    virtual bool isMinTTC();
    virtual bool isUserDefinedTTC();

    virtual bool isBreastCoil();
    virtual bool isBreastApplication();
    virtual bool isAnyCartesianRadialReorderingInPhaseEncodingPlane();
    virtual bool isGRASP();
    virtual bool isGRASPDynamicPreview();
    virtual bool isGradientReversal();
    virtual bool isRadialSelfGating();
    virtual bool isIceProgramRadialVIBESelfGatingAdd();

    // DIXON
    virtual bool isDixon();
    virtual	bool isFastDixon();
    virtual bool isConventionalDixon();
    virtual bool isOppAndInTwoPointDixon();
    virtual bool isMinAndOppTwoPointDixon();
    virtual bool isFlexTwoPointDixon();
    virtual bool isTwoPointDixon();
    virtual bool isMultiEchoDixon();
    virtual bool isOppAndInTwoPointScreeningDixon();
    virtual bool isTwoPointDixonWithoutAsymmetricEcho();
    virtual bool isTwoPointDixonWithFirstOppInPhases();
    virtual bool isDixonT2StarOrR2StarEvaluation();
    virtual bool isAnyMultiEchoDixonEvaluationOption();
    virtual bool isDixonFatOrWaterFractionEvaluation();
    virtual bool isDixonFatOrWaterImage();

    // Diffusion
    virtual bool isGradientReversalDiffusion();
    virtual bool isBetterSliceProfileDiffusion();
    virtual bool iSDTI();
    virtual bool isIVIM();

    // SMS
    virtual bool isSliceAcceleration();

    // SliceAdj
    virtual bool isSliceAdj();

    // Fat water contrast
    virtual bool isAnyFatSat();
    virtual bool isAnySat();
    virtual bool isWaterExcitation();
    virtual bool isFastWaterExcitation();
    virtual bool isAnyWaterExcitation();
    virtual bool isQuickFatSat();
    virtual bool isSPAIRFatSat();

private:
    MrProtFacade();
    MrProt mrProt;
};


#endif
