//	---------------------------------------------------------
//	  Copyright (C) Siemens AG 1999  All Rights Reserved.
//	---------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4\pkg\MrServers\MrImaging\seq\common\MODULE\MODULE_Routines.h
//	 Version: \main\6
//	  Author: STEMAL8Q
//	    Date: 2013-06-06 14:53:34 +02:00
//
//	    Lang: C++
//
//	 Descrip: Definition of some helpers.
//
//	 Classes:
//
//	---------------------------------------------------------
//	Siemens AG Medical Solutions

#ifndef __MODULE_Routines_H
#define __MODULE_Routines_H

#include "MrGlobalDefinitions/MrResult.h"


//  ceil, floor
#include <math.h>
//	---------------------------------------------------------
//  Used interfaces.
//	---------------------------------------------------------

namespace MODULE
{
    static inline int SIGN(const double& dIn)
    {
        return dIn < 0 ? -1 : 1; 
    }

    //  Euclid's algorithm
    //  Definition: An algorithm to compute the greatest common divisor
    //  of two positive integers.
    //  It is Euclid(a,b){if (b=0) then return a; else return Euclid(b, a mod b);}.
    //  The run time complexity is O(( log a)( log b)) bit operations.
    template<class INT_TYPE>
    static INT_TYPE GCD(INT_TYPE iA, INT_TYPE iB)
    {
        if(iB == 0)
        {
            return iA;
        }
        return GCD(iB,iA%iB);
    }
    template<class INT_TYPE>
    static INT_TYPE GCD(INT_TYPE iA, INT_TYPE iB, INT_TYPE iC)
    {
        return GCD(GCD(iA,iB),iC);
    }

    // calculates smallest common divisor
    template<class INT_TYPE>
    static INT_TYPE SCD(INT_TYPE iA, INT_TYPE iB)
    {
        if(iB == 0)
        {
            return iA;
        }
        else if(iA == 0)
        {
            return iB;
        }
        else if(iA > iB)
        {
            return SCD(iA - iB,iB);
        }
        else 
        {
            return SCD(iA,iB - iA);
        }
    }

    // calculates least common multiple of iA and iB
    template<class INT_TYPE>
    inline INT_TYPE LCM(INT_TYPE iA, INT_TYPE iB)
    {
        INT_TYPE iGCD = GCD(iA,iB);
        if( iGCD == 0 )
        {
            return iA*iB;
        }
        return ((iA*iB)/iGCD);
    }

    //  Returns nearest integer to dIn which is an integer multiple of lIncr.
    //  lIncr must be greater than zero.
    template<class INT_TYPE>
    static inline INT_TYPE IMULT(double dIn, INT_TYPE lIncr)
    {
        return dIn < 0
            ? lIncr*(static_cast<INT_TYPE>(dIn-double(lIncr)/2)/lIncr)
            : lIncr*(static_cast<INT_TYPE>(dIn+double(lIncr)/2)/lIncr);
    }


    //  Returns integer multiple of lIncr which is equal to or
    //  greater than lIn/dIn. Thereby lIncr must be greater than
    //  zero.
    template<class INT_TYPE>
    static inline INT_TYPE IMULT_CEIL(INT_TYPE lIn, INT_TYPE lIncr)
    {
        if(lIn%lIncr) lIn = lIncr*(lIn < 0 ? (lIn/lIncr) : (lIn/lIncr+1));
        return lIn;
    }
    template<class INT_TYPE>
    static inline INT_TYPE IMULT_CEIL(double dIn, INT_TYPE lIncr)
    {
        return IMULT_CEIL(static_cast<INT_TYPE>(ceil(dIn)), lIncr);
    }
    

    //  Returns integer multiple of lIncr which is equal to or
    //  less than lIn. Thereby lIncr must be greater than
    //  zero.
    template<class INT_TYPE>
    static inline INT_TYPE IMULT_FLOOR(INT_TYPE lIn, INT_TYPE lIncr)
    {
        if(lIn%lIncr != 0) lIn = lIncr*(lIn < 0 ? (lIn/lIncr-1) : (lIn/lIncr));
        return lIn;
    }
    template<class INT_TYPE>
    static inline INT_TYPE IMULT_FLOOR(double dIn, INT_TYPE lIncr)
    {
        return IMULT_FLOOR(static_cast<INT_TYPE>(floor(dIn)), lIncr);
    }

    //  Returns 2^uExp
    template<class INT_TYPE>
    static inline INT_TYPE RAISE_2_TO(INT_TYPE uExp)
    {
        return 0x1 << uExp;
    }

    static inline bool NLS_FAILED(NLS_STATUS i32Status)
    {
        return (i32Status & NLS_SEV) != NLS_SUCCESS;
    }

}  //  namespace MODULE


#endif  //  __MODULE_Routines_H
