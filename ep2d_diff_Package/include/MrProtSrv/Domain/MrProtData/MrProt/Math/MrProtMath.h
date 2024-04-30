//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProt\Math\MrProtMath.h
//	 Version:
//	  Author: Comp_ProBe
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: Numaris4::Measurement::CSequence::Prot::KSpace
//
//	 Classes:
//
//	-----------------------------------------------------------------------------

#pragma once

#ifdef _MSC_VER
#pragma once
#endif

#include "MrGlobalDefinitions/MrBasicTypes.h"


#undef __IMP_EXP

#ifdef BUILD_MrProt
  #define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"


class __IMP_EXP MrProtMath
{
private:
    // prevent instantation
    MrProtMath(void);
    
public:
    /// round up or down to nearest integer with a given precision
    static double roundToPrecision(double dVal, int iPrecision);
    static int32_t dbl2Lng(double x);
};

#undef __OWNER