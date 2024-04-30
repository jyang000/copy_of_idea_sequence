//	---------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	---------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProtocol\UILink\StdRoutines.h
//	 Version: \main\28
//	  Author: Comp_ProBe
//	    Date: 2011-09-30 16:14:36 +02:00
//
//	    Lang: C++
//
//	 Descrip: MR::ExamUI::ExamDb::Protocol
//
//	 Classes:
//
//	---------------------------------------------------------

#pragma once

#ifndef _STDROUTINES_H
#define _STDROUTINES_H

#include <stdio.h>
#include <stdlib.h>
#include <ostream>

#include "MrProtSrv/Domain/MrProtocol/libUILink/StdRoutines.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkLimited.h"
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Math/MrProtMath.h"

#ifdef BUILD_MrUILink
  #define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

namespace ArgListTools
{
	template<typename T>
	void Insert(char* argList[], uint32_t index, T value)
	{
		
#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4312) // yeah, I'll meet the people who made me do this cast in hell.
#endif
		
		argList[index] = (char*)value;
		
#ifdef _MSC_VER
# pragma warning(pop)
#endif
		
	}
}

inline double _rad2deg(double dRadiants, double _pi = 3.141592653589793)
{
    return static_cast<double>(dRadiants/_pi * 180.);
}

inline double _deg2rad(double lDegrees, double _pi = 3.141592653589793)
{
    return lDegrees * _pi/180.0;
}

// computes greatest common divisor
inline int gcd(int a, int b)
{
    int c = 0;
    while (b != 0)
    {    
       c = a%b;
       a = b;
       b = c; 
    }       
    return abs(a);
}

// computes least common multiple
inline int lcm(int a, int b)
{
    return a*b/gcd(a,b);
}

// returns 0 if n==0
// returns 1 if n==1
// returns n if n is a power of two
// returns the next highest power of 2 to n otherwise; e.g., n==126 -> returns 128 
unsigned int __IMP_EXP ceil_pow_of_two( unsigned int n );

/*// help functions to load stdProtRes to retrieve strings
typedef struct DLLModule* MODULE_HANDLE;
int UILinkLoadStdProtResString (unsigned  nId, LPTSTR pBuffer, size_t uBufferSize);
*/
template <class T>
class MrLimit;

typedef MrLimit<double> MrLimitDouble;

unsigned __IMP_EXP _formatFloat(char* argList[], double value, int decimalPlaces);
unsigned __IMP_EXP _formatFloat(char* arg_list[], double _val, double _incr);
double   _retrieveFloat(unsigned id,char *arg_list[]);

unsigned __IMP_EXP _formatTime(char* arg_list[], double dUSec);
unsigned __IMP_EXP _formatTATime(char* arg_list[], int seconds);

//	/////////////////////////////////////////////////////////
//  Given a min-max-incr-intervall, the function appends one or
//  more MrLimit intervalls to '_out'. Each of these MrLimit
//  intervalls 'x' is a subset of the input intervall and fullfills
//  exactly one of the following additional conditions:

//  a) x is subset of [0   (_incr)    10[      [ms]   
//  b) x is subset of [10  (10*_incr) 1000[    [ms]
//  c) x is subset of [1000(100*_incr)100000[  [ms]
//  d) ...
//  The function returns the number of limit intervalls appended.
//
int __IMP_EXP _timeLimits(std::vector<MrLimitDouble>& _out, double min_ms, double max_ms, double _incr_ms);

void __IMP_EXP _addTimingDependencyPtr(MrUILinkBase* const pThis);

bool __IMP_EXP _solveTiming(MrUILinkBase* const pThis, int32_t pos);
bool _solveTimingOptimized(MrUILinkBase* const pThis, int32_t pos, bool bSolveTomTR, bool bSolveTomTE);

class MrUILinkBase;



//	/////////////////////////////////////////////////////////
//  Given the direction in patient orientated coordinates
//  the function returns the resource id of the corresponding
//  UI-Label 
//
template<class TYPE>
class VectorPat;
  
unsigned __IMP_EXP _dimLabel(const VectorPat<double>& _in);
unsigned _points(const VectorPat<double>& _in);
unsigned _pointsTo(const VectorPat<double>& _in);


//	/////////////////////////////////////////////////////////
//  Given a Slice Object the function returns the resource id
//  of the corresponding string resource and writes the 
//  insertion values to arg_list.
//
class Slice;
namespace MrProtocolData
{
    class MrSliceData;
}
unsigned _oriId(const MrProtocolData::MrSliceData& e_in, char* arg_list[]);
unsigned _oriId(const Slice& e_in, char* arg_list[]);
unsigned _oriId(const VectorPat<double>& e_norm, char* arg_list[]);


//	/////////////////////////////////////////////////////////
//  Given the position vector in patient orientated coordinates
//  the function returns the resource id of the corresponding
//  string resource and writes the insertion values to arg_list.
//                                                             
#define _SHIFT_INCREMENT 0.1
#define _SHIFT_PRECISION 1
unsigned _posId(const VectorPat<double>& e_in, char* arg_list[], int iPrecision);
unsigned _posId_SBCS(MrUILinkBase* pThis, const VectorPat<double>& e_in, char* arg_list[], int iPrecision);
VectorPat<double> get_VectorPat_PCS_to_SBCS(MrUILinkBase* pThis, const VectorPat<double>& e_in);
VectorPat<double> get_VectorPat_SBCS_to_PCS(MrUILinkBase* pThis, const VectorPat<double>& e_in);

inline double _roundToPrecision(double _val, int _precision)
{
    return MrProtMath::roundToPrecision(_val, _precision);  
}

bool __IMP_EXP getSliceGroupOri(MrUILinkBase* const pThis, 
					  VectorPat<double> &e_phase,
					  VectorPat<double> &e_read,
					  VectorPat<double> &e_slice);

/********************************************************************************		
 Function: _calcRotAngleRad
				calculates from the vectors of phase encoding and the slice orientation 
               (the slice normal vector) the inplane rotation
               have a positive sign.
 ********************************************************************************/		
double _calcRotAngleRad( double Vnormal[3], double Vphase[3],SEQ::MainOrientation lOrientation);


static inline const char* _b2s(bool b)
{
    return b ? "true" : "false";
}



#endif
