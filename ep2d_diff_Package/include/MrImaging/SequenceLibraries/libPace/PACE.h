//  -----------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//  -----------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4\pkg\MrServers\MrImaging\seq\common\libPace\PACE.h
//     Version: \main\6
//      Author: stemal8q
//        Date: 2011-10-27 13:18:15 +02:00
//
//        Lang: C++
//
//     Descrip: Declaration of some functions exported from libPACE.
//
//     Classes:
//
//  -----------------------------------------------------------------

#ifndef _PACE_h
#define _PACE_h

//  Definition of NLS_STATUS
#include "MrGlobalDefinitions/MrResult.h"


//  Definition of 
//#include "MrProtSrv/Domain/CoreNative/NavigatorDefines.h"


//  -----------------------------------------------------------------
//  Import - Export - Control
//  -----------------------------------------------------------------

#ifdef __IMP_EXP
#undef __IMP_EXP
#endif 

#ifdef LOCAL_PACE_BUILD
#define __IMP_EXP
#pragma message( "Local PACE build" )
#else  //   LOCAL_PACE_BUILD not defined
#ifdef BUILD_libPACE
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control
#endif  //  LOCAL_PACE_BUILD

class SeqLim;
namespace MrProtocolData
{
    class MrProtData;
}
class MrProt;

#include <string>

namespace PACE
{
    //  Conversion of PACE parameters between different SW versions
    __IMP_EXP NLS_STATUS fConvProt(const MrProtocolData::MrProtData& rProtSrc,MrProtocolData::MrProtData& rProtDst);

    //  Standard initialization of the PACE related SeqLim parameters.
    //  The parameters iDef, iMode1, ..., iMode13 receive the modi
    //  supported by the sequence.
    //  The mode passed to the iDef argument defines the deafult.
    //  e.g. fInit(&rSeqLim,SEQ::RESP_COMP_OFF,SEQ::RESP_COMP_BREATH_HOLD,SEQ::RESP_COMP_TRIGGER);
    //  It is strongly recommended that the function is invoked after the other SeqLim parameters have been initialized.
    __IMP_EXP NLS_STATUS fInit(SeqLim* pSeqLim, int iDef, int iMode2 = 0, int iMode3 = 0, int iMode4 = 0, int iMode5 = 0, int iMode6 = 0, int iMode7 = 0, int iMode8 = 0, int iMode9 = 0, int iMode10 = 0, int iMode11 = 0, int iMode12 = 0, int iMode13 = 0);

    //  Writes statistic data of patient with name 'pszPatientName' to file '%MEASDAT%/PACE_STAT.txt'.
    //  Overwrites any information stored during previous calls of fSetStatData. 
    __IMP_EXP NLS_STATUS fSetStatData(int32_t i32NRespIntrv,double  dRespCycle_us, const char* pszPatientName);

    //  Retrievs patient name from AccMeasTemp and then calls fSetStatData(.,.,.)
    __IMP_EXP NLS_STATUS fSetStatData(int32_t i32NRespIntrv,double  dRespCycle_us);

    //  Reads file '%MEASDAT%/PACE_STAT.txt'.
    __IMP_EXP NLS_STATUS fGetStatData(int32_t& ri32NRespIntrv,double&  rdRespCycle_us, std::string& rPatientName);

    //  First calls fGetStatData(.,.,.). Then retrieves patient name from AccMeasTemp.
    //  Suceeds only if both patient names are equal.
    __IMP_EXP NLS_STATUS fGetStatData(int32_t& ri32NRespIntrv,double&  rdRespCycle_us);


    //  First stops the measurement and then calls fSetStatData(i32NRespIntrv,dRespCycle_us).
    //  After this starts the measurement again.
    //  Should be called at the end of fSEQRun
    __IMP_EXP NLS_STATUS fFlushStatData(int32_t i32NRespIntrv, double  dRespCycle_us, int32_t i32FillBeforeHalt_us, int32_t i32FillAfterHalt_us);


    //  Returns 0, if PACE respiratory triggering is not selected.
    //  Otherwise:  it tries to retrieve the respiratory cycle of the current patient. 
    __IMP_EXP double fEstimateMeasTime_us(const MrProtocolData::MrProtData* pProt, int32_t i32TriggerHalts);


} //  namespace PACE

#ifdef __IMP_EXP_DATA
#undef __IMP_EXP_DATA
#endif

#ifdef __IMP_EXP
#undef __IMP_EXP
#endif

#endif // _PACE_h
