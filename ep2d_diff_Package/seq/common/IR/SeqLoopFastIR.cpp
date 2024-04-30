//  -----------------------------------------------------------------
//  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential
//  -----------------------------------------------------------------
//
//          Project: NUMARIS/4
//             File: \n4\pkg\MrServers\MrImaging\seq\common\IR\SeqLoopFastIR.cpp
//          Version: \main\3
//
//         Language: C++
//
//      Description: SeqLoop variant which implements a (Fast) (I)nversion-(R)ecovery scheme.
//          Structs:
//
//  -----------------------------------------------------------------

//  -----------------------------------------------------------------
//  Used interfaces
//  -----------------------------------------------------------------

//  Definition of class SeqLoopFastIR
#include "MrImaging/seq/common/IR/SeqLoopFastIR.h"

//  Definition of KERNEL_IMAGE
#include "MrMeasSrv/SeqIF/csequence.h"

//  Definition of MDH_LAST_BLADE_IN_TR
#include "MrMeasSrv/SeqIF/MDH/MdhProxy.h"

//  Definition of fSBBFillTimeRun
#include "MrImagingFW/libSBBFW/libSBBFW.h"

//  Definition of fSUVerifyRFSpoil
#include "MrImagingFW/libSeqUtilFW/libSeqUtilFW.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

#include "MrProtSrv/Domain/CoreNative/MrNavigator.h"    // * For CNavigator
#include "MrProtSrv/Domain/CoreNative/MrKSpace.h"         // * For MrKSpace
#include "MrProtSrv/Domain/CoreNative/MrPreparationPulses.h"         
#include "MrProtSrv/Domain/CoreNative/MrBlade.h"    
#include "MrProtSrv/Domain/CoreNative/MrPhysiology.h"
#include "MrProtSrv/Domain/CoreNative/MrMds.h"
#include "MrProtSrv/Domain/CoreNative/MrPat.h"
//#include "MrProtSrv/Domain/CoreNative/MrSliceGroupData.h"



// includes for old classes
#include "MrProtSrv/Domain/MrProtData/MrProt/MrSliceGroup.h"

//#include "MrProtSrv/Domain/MrProtData/MrProt/PreparationPulses/MrPreparationPulses.h"
//#include "MrProtSrv/Domain/MrProtData/MrProt/MDS/MDS.h"  
#include "MrProtSrv/Domain/MrProtData/MrProt/Physiology/MrPhysiology.h"
//#include "MrProtSrv/Domain/MrProtData/MrProt/KSpace/KSpaceData.h"

#include "MrImaging/libSBB/SeqLoopPARA.h"
// FileNameExp
#include "MrMeasSrv/MeasUtils/MeasUtils.h"


//  -----------------------------------------------------------------
//  Additional Definitions
//  -----------------------------------------------------------------

//extern SeqLoopFastIR& theSeqLoopFastIR();

namespace SEQ_NAMESPACE
{
    /*
    static void LOCATE_ERROR(const char* pszFile,int iLine)
    {
    SEQ_TRACE_ERROR.print("Error at %s(%d).",pszFile,iLine);
    }

    static void CAUGHT_EXCEPT(const char* pszFile,int iLine)
    {
    SEQ_TRACE_ERROR.print("Caught unknown exception at %s(%d).",pszFile,iLine);
    }
    */

    static bool ENFORCE_IR_SCHEME_AUTO(const MrProt &rMrProt)
    {
        return  rMrProt.preparationPulses().getucInversion() != SEQ::SLICE_SELECTIVE
            || (rMrProt.NavigatorParam().getlRespComp() != SEQ::RESP_COMP_OFF )
            || (rMrProt.getsPhysioImaging().getlMethod1() != SEQ::METHOD_NONE && rMrProt.getsPhysioImaging().getlMethod1() != SEQ::METHOD_TRIGGERING)
            || rMrProt.sliceGroupList().size() > 1
            || rMrProt.getsKSpace().getucDimension() != SEQ::DIM_2
            || rMrProt.getsPhysioImaging().getlPhases() > 1
            || rMrProt.getsPrepPulses().getucDarkBlood()
            || rMrProt.getsKSpace().getucMultiSliceMode() != SEQ::MSM_INTERLEAVED
            || rMrProt.getsMds().getulMdsModeMask() != SEQ::MDS_OFF
            ;
    }

    //  Returns integer multiple of lIncr which is equal to or
    //  greater than lIn/dIn. Thereby lIncr must be greater than
    //  zero.
    static inline int IMULT_CEIL(int iIn, int iIncr)
    {
        if(iIn%iIncr) iIn = iIncr*(iIn < 0 ? (iIn/iIncr) : (iIn/iIncr+1));
        return iIn;
    }

    static inline int IMULT_FLOOR(int iIn, int iIncr)
    {
        if(iIn%iIncr) iIn = iIncr*(iIn < 0 ? (iIn/iIncr-1) : (iIn/iIncr));
        return iIn;
    }

    static inline int IMULT_FLOOR(double dIn, int iIncr)
    {
        return IMULT_FLOOR(int(floor(dIn)), iIncr);
    }

    //  Returns nearest integer to dIn which is an integer multiple of iIncr.
    //  iIncr must be greater than zero.
    static inline int IMULT(double dIn, int iIncr)
    {
        return dIn < 0
            ? iIncr*(static_cast<int>(dIn-iIncr/2.)/iIncr)
            : iIncr*(static_cast<int>(dIn+iIncr/2.)/iIncr);
    }

    template <class _Tp>
    static inline const _Tp& loc_max(const _Tp& __a, const _Tp& __b) {
        return  __a < __b ? __b : __a;
    }

} // namespace SEQ_NAMESPACE

bool SEQ_NAMESPACE::fCalcTBlock(int& riTBlock_us, int iTI_us, int iIRTime_us, int iSBBScanTime_us, int iKernelTime_us, int iMinFillTime_us, int iKernelOffset)
{
    //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+
    //  |IR|  |S|Kernel_n         |  |IR|  |S|Kernel_n+1       |  |IR|  |S|Kernel_n+2       |
    //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  ...
    //  <--------- TBlock ----------->
    //                               <--------- TBlock ----------->     <-> SBBSccanTime
    //  <-------> TI for KernelOffset==0 

    //  IR  ...  inversion recovery pulse
    //  S   ...  Saturation bands, chemical sat, spoilers before kernel, ...

    //  TI = kernelOffset*TBloc+IRTime+FillTimePastIR+SBBScanTime
    int iFillTimePastIR_us = 0, iFillTimePastKernel_us = 0;
    riTBlock_us = iSBBScanTime_us+iKernelTime_us+iIRTime_us;
    for(;iFillTimePastIR_us >= 0;riTBlock_us+=GRAD_RASTER_TIME)
    {
        //  TI = kernelOffset*TBloc+IRTime+FillTimePastIR+SBBScanTime
        iFillTimePastIR_us  = iTI_us-iKernelOffset*riTBlock_us-iIRTime_us-iSBBScanTime_us;
        //  TBlock = iIRTime_us+FillTimePastIR+iSBBScanTime_us+iKernelTime_us+FillTimePastKernel
        iFillTimePastKernel_us = riTBlock_us-iIRTime_us-iFillTimePastIR_us-iSBBScanTime_us-iKernelTime_us;

        if( (iFillTimePastIR_us >= 0) && (iFillTimePastKernel_us >= 0) && ((iFillTimePastIR_us+iFillTimePastKernel_us) >= iMinFillTime_us) )
        {
            //  Valid solution found
            return true;
        }
    }
    return false;
}

bool SEQ_NAMESPACE::fTestTBlock(int iTBlock_us, int iTI_us, int iIRTime_us, int iSBBScanTime_us, int iKernelTime_us, int iMinFillTime_us, int iKernelOffset)
{
    //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+
    //  |IR|  |S|Kernel_n         |  |IR|  |S|Kernel_n+1       |  |IR|  |S|Kernel_n+2       |
    //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  ...
    //  <--------- TBlock ----------->
    //                               <--------- TBlock ----------->     <-> SBBSccanTime
    //  <-------> TI for KernelOffset==0 

    //  IR  ...  inversion recovery pulse
    //  S   ...  Saturation bands, chemical sat, spoilers before kernel, ...

    //  TI = kernelOffset*TBloc+IRTime+FillTimePastIR+SBBScanTime

    //  TI = kernelOffset*TBloc+IRTime+FillTimePastIR+SBBScanTime
    const int iFillTimePastIR_us  = iTI_us-iKernelOffset*iTBlock_us-iIRTime_us-iSBBScanTime_us;
    //  TBlock = iIRTime_us+FillTimePastIR+iSBBScanTime_us+iKernelTime_us+FillTimePastKernel
    const int iFillTimePastKernel_us = iTBlock_us-iIRTime_us-iFillTimePastIR_us-iSBBScanTime_us-iKernelTime_us;

    return (iFillTimePastIR_us >= 0)
        && (iFillTimePastKernel_us >= 0)
        && ((iFillTimePastIR_us+iFillTimePastKernel_us) >= iMinFillTime_us
        );
}


//  -----------------------------------------------------------------
//  Implementation of SEQ_NAMESPACE::SeqLoopFastIR::CONC out of line members
//  -----------------------------------------------------------------

int SEQ_NAMESPACE::SeqLoopFastIR::getNSweepsPerConc() const
{
    return this->m_iNSweepsPerConc;
}
void SEQ_NAMESPACE::SeqLoopFastIR::setNSweepsPerConc(int iVal)
{
    this->m_iNSweepsPerConc = iVal;
}

int  SEQ_NAMESPACE::SeqLoopFastIR::getCoolPauseWithinKernelTime_us() const
{
    return this->m_iCoolPauseWithinKernelTime_us;
}
void SEQ_NAMESPACE::SeqLoopFastIR::setCoolPauseWithinKernelTime_us(int iVal)
{
    this->m_iCoolPauseWithinKernelTime_us = iVal;
}

bool SEQ_NAMESPACE::SeqLoopFastIR::runBloc_uniform(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS asSlc[],sREADOUT* pADC)
{
    if( size_t(m_lConcatenationCounter) >= m_asConc.size() )
    {
        setNLSStatus( MRI_SBB_SBB_ERROR );
        SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
        return false;
    }
    //  Pointer to current concatenation
    std::vector<CONC>::iterator itConc = m_asConc.begin()+m_lConcatenationCounter;

    //  i) execute IR-pulse
    const long lCIR = this->m_iCExe%this->m_iModulo;
    if( (size_t(lCIR) >= (*itConc).m_aiChronPos.size()) || ((*itConc).m_aiChronPos[lCIR] == m_SlicesToMeasure) || (this->m_iCExe+this->m_iKernelOffset >= this->m_iNExe) )
    {
        //  Skip IR-pulse and execute empty event-bloc instead
        setNLSStatus( fSBBFillTimeRun(SBBIRsel.getDurationPerRequest()+this->m_iTFillPastIR_us) );
        if ( (m_NLSStatus & NLS_SEV) != NLS_SUCCESS )
        {
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }
    }
    else
    {
        //  Excute IR-pulse
        sSLICE_POS* pIRPos = &asSlc[(*itConc).m_aiChronPos[lCIR]];
        if( !SBBIRsel.run(rMrProt,rSeqLim,rSeqExpo,pIRPos) )
        {
            setNLSStatus( SBBIRsel.getNLSStatus() );
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }
#ifdef WIN32
        if( this->m_sChron.is_open() )
        {
            this->m_sChron << "IR[";
            this->m_sChron.width(2);
            this->m_sChron << pIRPos->getSliceIndex() << "] : ";
            this->m_sChron.precision(2);
            this->m_sChron.flags(std::ios::fixed);
            this->m_sChron << pIRPos->getSliceShift() << std::endl;
        }
#endif
        //  Execute fill time past IR
        if( this->m_iTFillPastIR_us > 0)
        {
            setNLSStatus( fSBBFillTimeRun(this->m_iTFillPastIR_us) );
            if ( (m_NLSStatus & NLS_SEV) != NLS_SUCCESS )
            {
                SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                return false;
            }
        }
    }

    //  ii) Execute kernels
    if( this->m_iCExe < this->m_iKernelOffset )
    {
        //  Skip Kernel and execute fill-time
        setNLSStatus( fSBBFillTimeRun(m_lEffScanTime+this->m_iTFillPastKernel_us) );
        if ( (m_NLSStatus & NLS_SEV) != NLS_SUCCESS )
        {
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }
    }
    else
    {
        //  Chronological blade counter
        //const long lCExcit = (this->m_iCExe-this->m_iKernelOffset) < m_lPreparingScans*this->m_iModulo
        //    ? (this->m_iCExe-this->m_iKernelOffset)/this->m_iModulo
        //    : (this->m_iCExe-this->m_iKernelOffset-m_lPreparingScans*this->m_iModulo)/this->m_iModulo
        //    ;
        //  Total number of blades
        //const long lNExcit    = (this->m_iNExe-this->m_iKernelOffset-m_lPreparingScans*this->m_iModulo)/this->m_iModulo;
        m_lKernelMode = (this->m_iCExe-this->m_iKernelOffset) < m_lPreparingScans*this->m_iModulo
            ? KERNEL_PREPARE
            : KERNEL_IMAGE
            ;

        MdhProxy& rMDH = pADC->getMDH();

        rMDH.setClin(0);
        rMDH.setCacq(0);
        rMDH.setCpar(0);
        rMDH.setCphs(0);
        const int iCKernel = (this->m_iCExe-this->m_iKernelOffset)%this->m_iModulo;
        //  Chronological BLADE counter
        //rMDH.setFirstScanInSlice(m_bIsFirstScanInSlice = (lCBld == 0));
        //rMDH.setLastScanInSlice (m_bIsLastScanInSlice  = (lCBld == lNBld-1));


        //const long lCKernel = (this->m_iCExe-this->m_iKernelOffset)%this->m_iModulo;

        if( size_t(iCKernel) >= (*itConc).m_aiChronPos.size() )
        {
            setNLSStatus(MRI_SBB_SBB_ERROR);
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }

        m_bFirstSliceInConcat = iCKernel==0;
        m_lSliceIndex = (*itConc).m_aiChronPos[iCKernel];
        rMDH.setCslc((unsigned short)m_lSliceIndex);

        if( m_lSliceIndex < 0 || m_lSliceIndex >= m_SlicesToMeasure )
        {
            //  Play out fill time
            setNLSStatus( fSBBFillTimeRun(m_lEffScanTime) );
            if ( (m_NLSStatus & NLS_SEV) != NLS_SUCCESS )
            {
                SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                return false;
            }
        }
        else
        {
            //  Unfortunately RF spoiling control is not performed
            //  within SeqLoop::runKernelCallsLoop
            if(rMrProt.fastImaging().getulEnableRFSpoiling() && m_PerformSATs)
            {
                sSLICE_POS* pSlcPos = &asSlc[m_lSliceIndex];
                setNLSStatus( fSUVerifyRFSpoil(
                    rMrProt,                            //  IMP: user choice parameters
                    rSeqLim,                            //  IMP: limits from fSEQInit()
                    m_bFirstSliceInConcat?0:1,          //  IMP: copy data   TRUE/FALSE
                    pSlcPos->getSliceIndex(),           //  IMP: Anatomic slice number
                    &m_dRFSpoilPhase,                      //  EXP: RF spoiling parameter
                    &m_dRFSpoilIncrement,                  //  EXP: RF spoiling parameter
                    &m_dRFSpoilPhasePrevSlice,          //  EXP: RF spoiling parameter
                    &m_dRFSpoilIncrementPrevSlice,      //  EXP: RF spoiling parameter
                    &m_dRFPrevSlicePosSag,              //  EXP: RF spoiling parameter
                    &m_dRFPrevSlicePosCor,              //  EXP: RF spoiling parameter
                    &m_dRFPrevSlicePosTra,              //  EXP: RF spoiling parameter
                    &m_dRFPrevSliceNormalSag,              //  EXP: RF spoiling parameter
                    &m_dRFPrevSliceNormalCor,              //  EXP: RF spoiling parameter
                    &m_dRFPrevSliceNormalTra              //  EXP: RF spoiling parameter
                    ) );
                if ( (m_NLSStatus & NLS_SEV) != NLS_SUCCESS )
                {
                    SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                    return false;
                }
            }
            if( m_FirstMethod == SEQ::METHOD_TRIGGERING && m_TReffective )
            {
                if(m_lKernelMode & KERNEL_PREPARE || m_bFirstSliceInConcat == false )
                {
                    rMDH.setTREffectiveEnd   (false);
                    rMDH.setTREffectiveBegin (false);
                }
                else
                {
                    if( m_bIsFirstScanInSlice )
                    {
                        m_lFlipFlop = 0;
                    }
                    if ( ++m_lFlipFlop %= 2 )
                    {
                        rMDH.setTREffectiveBegin (true );
                        rMDH.setTREffectiveEnd   (false);
                    }
                    else
                    {
                        rMDH.setTREffectiveBegin (false);
                        rMDH.setTREffectiveEnd   (true );
                    }
                }
            }
            rMDH.setCslc((unsigned short)m_lSliceIndex);
            //rMDH.setCidd(static_cast<unsigned short>(lCBld));
#ifdef WIN32
            if( this->m_sChron.is_open() )
            {
                sSLICE_POS* pSlcPos = &asSlc[m_lSliceIndex];
                this->m_sChron << "SL[";
                this->m_sChron.width(2);
                this->m_sChron << pSlcPos->getSliceIndex() << "] : ";
                this->m_sChron.precision(2);
                this->m_sChron.flags(std::ios::fixed);
                this->m_sChron << pSlcPos->getSliceShift() << std::endl;
            }
#endif
            //  Execute kernel
            if( !runKernelCallsLoop(rMrProt,rSeqLim,rSeqExpo,asSlc,pADC) )
            {
                SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                return false;
            }

            fRTSetReadoutEnable(1);
            m_bFirstSliceInConcat = false;
        }


        //  Play out fill time past kernel
        if( this->m_iTFillPastKernel_us > 0 )
        {
            setNLSStatus( fSBBFillTimeRun(this->m_iTFillPastKernel_us) );
            if( (m_NLSStatus & NLS_SEV) != NLS_SUCCESS )
            {
                SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                return false;
            }
        }
    }
    return true;
}



//  -----------------------------------------------------------------
//  Implementation of SeqLoopFastIR out-of-line member functions
//  -----------------------------------------------------------------

SEQ_NAMESPACE::SeqLoopFastIR::SeqLoopFastIR()
    : BASE_TYPE()
    , m_iCExe(0)
    , m_iNExe(0)
    , m_iModulo(1)
    , m_iKernelOffset(0)
    , m_iTFillPastIR_us(0)
    , m_iNIRPerMeas(0)
    , m_iTFillPastKernel_us(0)
    , m_iTBlock_us(0)
    , m_iCoolPauseWithinKernelTime_us(0)
    , m_iNSweepsPerConc(2)
{}

SEQ_NAMESPACE::SeqLoopFastIR::~SeqLoopFastIR()
{}

//  Maps the anatomical slice position to the chronological position
static int MAP_APOS_TO_CPOS(const SliceSeries& rSeries, int iAPos)
{
    const int iSize = rSeries.getlSize();
    if( iAPos < 0 || iAPos >= iSize )
    {
        return iSize;
    }
    const Slice* pSlice = &(rSeries.aAt(iAPos));
    int iCPos = 0;
    for(;iCPos != iSize;++iCPos)
    {
        if( &(rSeries.chronological(iCPos)) == pSlice )
        {
            return iCPos;
        }
    }
    SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
    return iSize;
}


bool SEQ_NAMESPACE::SeqLoopFastIR::calcFillTimes_uniform(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &,long lWantedTR, long, long* plNegativeFillTime)
{
    setNLSStatus(MRI_SBB_SBB_NORMAL);
    m_lTRneeded = m_lTIneeded = 0;

    //  Number of slice groups must be one
    if( rMrProt.sliceGroupList().size() != 1 || rMrProt.physiology().method(1) != SEQ::METHOD_NONE )
    {
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        setNLSStatus(MRI_SBB_SBB_ERROR);
        return false;
    }
    //  The following line is needed within limit calculation
    m_SlicesToMeasure = rMrProt.sliceSeries().getlSize();

    //  Number Conc must be less than number slices
    if( m_lConcatenations < 1 || m_SlicesToMeasure < m_lConcatenations )
    {
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        setNLSStatus (MRI_SBB_SBB_INVALID_CONCATENATIONS);
        return false;
    }

    //  Actual TR

    long lTR_us = m_lActualTR = lWantedTR < 0
        ? long(rMrProt.tr()[0])
        : lWantedTR
        ;

    //  Number of blocs executed until the slice series repeats
    this->m_iModulo = int(this->m_SlicesToMeasure+this->m_lConcatenations-1)/int(this->m_lConcatenations);

    //  In each bloc we must have time for m_lKernel+IR-pulse
    m_DummyIRTime = SBBIRsel.getDurationPerRequest();

    m_lTRMin = IMULT_CEIL(int(m_lEffScanTime+m_DummyIRTime),GRAD_RASTER_TIME)*this->m_iModulo;

    if( lTR_us < m_lTRMin )
    {
        *plNegativeFillTime = lTR_us-m_lTRMin;
        setNLSStatus(MRI_SBB_SBB_NEGATIV_TRFILL);
        m_lTRneeded = m_lTRMin;
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);

        //  Proceed with minimum TR
        lTR_us = m_lTRMin;
        if( rMrProt.preparationPulses().getucFreezeSupTissue() )
        {
            m_lTRneeded = m_lTIneeded = 0;
            return false;
        }
    }
    //  Duration of one Bloc
    this->m_iTBlock_us = IMULT_CEIL(int(lTR_us)/this->m_iModulo,GRAD_RASTER_TIME);

    //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+
    //  |IR|  |S|Kernel_n-2       |  |IR|  |S|Kernel_n-1       |  |IR|  |S|Kernel_n         |
    //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  ...
    //  <--------- TBlock ----------->
    //                               <--------- TBlock ----------->     <-> SBBSccanTime
    //  <-------------------- TI for KernelOffset==2 --------------------->


    //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+
    //  |IR|  |S|Kernel_n         |  |IR|  |S|Kernel_n+1       |  |IR|  |S|Kernel_n+2       |
    //  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  +--+  +-+-------+-+-------+  ...
    //  <--------- TBlock ----------->
    //                               <--------- TBlock ----------->     <-> SBBSccanTime
    //  <-------> TI for KernelOffset==0 

    //  IR  ...  inversion recovery pulse
    //  S   ...  Saturation bands, chemical sat, spoilers before kernel, ...

    long lTI_us = rMrProt.ti()[0];

    if( lTI_us < (m_DummyIRTime+m_lSBBScanTime) )
    {
        m_lTIneeded = m_DummyIRTime+m_lSBBScanTime;
        *plNegativeFillTime = lTI_us-m_lTIneeded;
        setNLSStatus (MRI_SBB_SBB_NEGATIV_TIFILL);
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        if( rMrProt.preparationPulses().getucFreezeSupTissue() )
        {
            m_lTRneeded = m_lTIneeded = 0;
        }
        return false;
    }
    this->m_iKernelOffset = int(lTI_us-(m_DummyIRTime+m_lSBBScanTime))/this->m_iTBlock_us;

    //  TI is defined from the start of the IR-SBB to the start of the kernel
    this->m_iTFillPastIR_us = IMULT_CEIL(int(lTI_us-m_DummyIRTime-m_lSBBScanTime-this->m_iKernelOffset*this->m_iTBlock_us),GRAD_RASTER_TIME);
    // m_lEffScanTime = Scan time required by one kernel call + Sats + Spoilers
    this->m_iTFillPastKernel_us = this->m_iTBlock_us-(this->m_iTFillPastIR_us+int(this->m_lEffScanTime)+int(this->m_DummyIRTime));
    if( this->m_iTFillPastKernel_us < 0 )
    {
        this->m_lTIneeded = IMULT_FLOOR(this->m_iKernelOffset*this->m_iTBlock_us+int(this->m_DummyIRTime)+int(this->m_lSBBScanTime)+this->m_iTBlock_us-int(m_lEffScanTime+m_DummyIRTime),GRAD_RASTER_TIME);
        *plNegativeFillTime = lTI_us-m_lTIneeded;
        setNLSStatus (MRI_SBB_SBB_NEGATIV_TIFILL);
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        if( rMrProt.preparationPulses().getucFreezeSupTissue() )
        {
            m_lTRneeded = m_lTIneeded = 0;
        }
        return false;
        //m_lTRneeded = m_lTRMin = (this->m_iTFillPastIR_us+m_lEffScanTime*m_lNKernelPerIR+m_DummyIRTime)*this->m_iKernelOffset;
        //m_lTRneeded = (this->m_iTFillPastIR_us+m_lEffScanTime*m_lNKernelPerIR+m_DummyIRTime)*this->m_iModulo;
        //*plNegativeFillTime = m_lActualTR-m_lTRneeded;
        //setNLSStatus (MRI_SBB_SBB_NEGATIV_TRFILL);
        //SEQ_TRACE_ERROR_COND(pSeqLim->isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        //return false;
    }
    m_lTRFillEnd  = 0;
    this->m_iNExe = int(m_LinesToMeasure*m_FreeLoopLength*m_AcquisitionsOuter+m_lPreparingScans)*this->m_iModulo+this->m_iKernelOffset;
    m_dMeasureTimeInSecondMeasUsec = m_dMeasureTimeInFirstMeasUsec = int(this->m_lConcatenations)*this->m_iNExe*this->m_iTBlock_us;
    if(rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA || rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST)
    {
        // duration of PAT ref scans calls = duration * slices 
        m_dMeasureTimeInFirstMeasUsec += SBBPATRefScan.getDurationPerRequest() * SBBPATRefScan.getRequestsPerMeasurement();
    }
    return m_lTIneeded == 0 && m_lTRneeded == 0;
}

bool SEQ_NAMESPACE::SeqLoopFastIR::calcFillTimesOnly (MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo& rSeqExpo, long lWantedTR, long lTimeForOsc)
{
    //  The initialization of the member variable m_lIIRScheme
    //  is repeated since the TI solve handler uses calculateTRTIFillTimes without preparation
    if( rMrProt.preparationPulses().getucInversion() != SEQ:: SLICE_SELECTIVE )
    {
        this->m_iIRScheme = SEQ::IR_SCHEME_AUTO;
    }
    else
    {
        this->m_iIRScheme = rMrProt.getsPrepPulses().getucIRScheme();
    }
    if( this->m_iIRScheme != SEQ::IR_SCHEME_SEQUENTIAL )
    {
        return BASE_TYPE::calcFillTimesOnly(rMrProt,rSeqLim,rSeqExpo,lWantedTR,lTimeForOsc);
    }

    if( ENFORCE_IR_SCHEME_AUTO(rMrProt) )
    {
        setNLSStatus (MRI_SBB_SBB_ERROR);
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        return false;
    }
    //  Used by calcFillTimesOnly. Stored for TI calculation

    //  Scan time required by one kernel call + Sats + Spoilers
    m_lEffScanTime = m_lSBBScanTime + m_lKernelScanTime - this->m_iCoolPauseWithinKernelTime_us;
    if(m_bSpoilGradAfterKernel)
    {
        m_lEffScanTime += getSpoilGradTime();
    }
    this->m_lTRFill = 0;
    long lNegativeFillTime_us = 0;
    if( !calcFillTimes_uniform(rMrProt,rSeqLim,rSeqExpo,lWantedTR,lTimeForOsc,&lNegativeFillTime_us) )
    {
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        return false;
    }
    if( (this->m_iTFillPastIR_us + this->m_iTFillPastKernel_us) < this->m_iCoolPauseWithinKernelTime_us )
    {
        this->setNLSStatus(MRI_SBB_SBB_NEGATIV_TRFILL);
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        return false;
    }
    return true;
}

void SEQ_NAMESPACE::SeqLoopFastIR::calcMeasurementTimeUsec (MrProt &rMrProt, SeqLim &rSeqLim)
{
    if( this->m_iIRScheme != SEQ::IR_SCHEME_SEQUENTIAL )
    {
        return BASE_TYPE::calcMeasurementTimeUsec(rMrProt,rSeqLim);
    }
    m_dMeasureTimeInSecondMeasUsec = m_dMeasureTimeInFirstMeasUsec = int(this->m_lConcatenations)*this->m_iNExe*this->m_iTBlock_us;
    if(rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA || rMrProt.PAT().getucRefScanMode() == SEQ::PAT_REF_SCAN_EXTRA_GRE_FAST)
    {
        // duration of PAT ref scans calls = duration * slices 
        m_dMeasureTimeInFirstMeasUsec += SBBPATRefScan.getDurationPerRequest() * SBBPATRefScan.getRequestsPerMeasurement();
    }
    m_dTrigHaltTimeInFirstMeasUsec  = m_dTrigHaltTimeInSecondMeasUsec = 0;
    m_lTotalPhysioHalts             = 0;
    m_NumberOfRelevantADCs          = 0;
    m_dPrepareTimeUsec              = 0.0;
    m_dPhaseCorTimeUsec             = 0.0;
}
void SEQ_NAMESPACE::SeqLoopFastIR::calcNoOfRelevantADCs (MrProt &rMrProt, SeqLim &rSeqLim)
{
    if( this->m_iIRScheme != SEQ::IR_SCHEME_SEQUENTIAL )
    {
        BASE_TYPE::calcNoOfRelevantADCs(rMrProt,rSeqLim);
        return;
    }
    m_lTotalPhysioHalts             = 0;
    m_NumberOfRelevantADCs          = 0;
    m_ADCCounterStep       = 0;
    std::vector<PACE_PARA::CONC>::iterator itC = this->m_pPACEPara->m_asConc.begin()
        ,                                 itE = this->m_pPACEPara->m_asConc.end()
        ;

    for(;itC != itE;++itC)
    {
        (*itC).m_lNRRx = 0;
    }
}

bool SEQ_NAMESPACE::SeqLoopFastIR::prep (MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo)
{
    if( rMrProt.preparationPulses().getucInversion() != SEQ:: SLICE_SELECTIVE )
    {
        this->m_iIRScheme = SEQ::IR_SCHEME_AUTO;
    }
    else
    {
        this->m_iIRScheme = rMrProt.getsPrepPulses().getucIRScheme();
    }
    if( this->m_iIRScheme != SEQ::IR_SCHEME_AUTO )
    {
        if( ENFORCE_IR_SCHEME_AUTO(rMrProt) )
        {
            setNLSStatus (MRI_SBB_SBB_ERROR);
            SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }
    }
    if( !BASE_TYPE::prep(rMrProt,rSeqLim,rSeqExpo) )
    {
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        return false;
    }
    if( !prepConc(rMrProt,rSeqLim,rSeqExpo) )
    {
        SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
        return false;
    }
    return true;
}

bool SEQ_NAMESPACE::SeqLoopFastIR::prepConc(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo&)
{
    //  Note: prepConc ist not a base class member! 
    if( this->m_iIRScheme == SEQ::IR_SCHEME_SEQUENTIAL )
    {
        //const SliceSeries& rSeries = pProt->sliceSeries();
        SliceSeries rSeries = rMrProt.sliceSeries();

        m_SlicesToMeasure = rSeries.getlSize();

        //  Number of slice groups must be one
        if( rMrProt.sliceGroupList().size() != 1 || m_SlicesToMeasure < m_lConcatenations )
        {
            SEQ_TRACE_ERROR_COND(rSeqLim.isContextNormal() , "Error at %s(%d).",__FILE__,__LINE__);
            setNLSStatus(MRI_SBB_SBB_ERROR);
            return false;
        }

        //  Number of blocs executed until the slice series repeats
        this->m_iEChronPos_1stSweep = this->m_iModulo = int(m_SlicesToMeasure+m_lConcatenations-1)/int(this->m_lConcatenations);

        m_asConc.resize(m_lConcatenations);
        std::vector<CONC>::iterator itConc = m_asConc.begin();

        int iAPos, iCIR;
        std::vector<int>::iterator itChronPos;


        m_lKernelRequestsInFirstMeas = 0;
        this->m_iNIRPerMeas                = 0;
        const SEQ::SeriesMode iExcitOrder = rMrProt.getsSliceArray().getucMode();

        if( this->m_iNSweepsPerConc < 1 )
        {
            this->m_iNSweepsPerConc = (this->m_lConcatenations == 1) && (iExcitOrder == SEQ::INTERLEAVED)
                ? 2
                : 1
                ;
        }
        if( m_lConcatenations == 1 )
        {
            if( iExcitOrder == SEQ::ASCENDING )
            {
                iAPos = 0;
            }
            else if( iExcitOrder == SEQ:: DESCENDING )
            {
                iAPos = int(m_SlicesToMeasure)-1;
            }
            else
            {
                iAPos = m_SlicesToMeasure%2 == 0 ? 1 : 0;
            }
            //  this->m_iModulo = SlicesPerConc //((m_SlicesToMeasure)+m_lConcatenations-1)/m_lConcatenations;
            (*itConc).m_aiChronPos.resize(this->m_iModulo);
            itChronPos = (*itConc).m_aiChronPos.begin();

            for(iCIR = 0; iCIR != this->m_iModulo;++iCIR)
            {
                if( iAPos >= m_SlicesToMeasure )
                {
                    if( iExcitOrder != SEQ::INTERLEAVED )
                    {
                        setNLSStatus (MRI_SBB_SBB_ERROR);
                        SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                        return false;
                    }
                    //  Second swap
                    iAPos = m_SlicesToMeasure%2 == 0 ? 0 : 1;
                    m_iEChronPos_1stSweep = int(itChronPos-(*itConc).m_aiChronPos.begin());
                }
                if( iAPos < m_SlicesToMeasure )
                {
                    ++this->m_iNIRPerMeas;
                }
                *itChronPos = MAP_APOS_TO_CPOS(rSeries,iAPos);
                if( *itChronPos++ < m_SlicesToMeasure )
                {
                    ++m_lKernelRequestsInFirstMeas;
                }
                switch( iExcitOrder )
                {
                case SEQ::DESCENDING:
                    --iAPos;
                    break;
                case SEQ::ASCENDING:
                    ++iAPos;
                    break;
                case SEQ::INTERLEAVED:
                    iAPos += 2;
                    break;
                default:
                    setNLSStatus (MRI_SBB_SBB_ERROR);
                    SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                    return false;
                }
            }
        }
        else
        {
            for(;itConc != m_asConc.end();++itConc)
            {
                (*itConc).m_aiChronPos.resize(this->m_iModulo);
                itChronPos = (*itConc).m_aiChronPos.begin();

                iAPos = int(itConc-m_asConc.begin());
                //  Count actual # slices
                int iNSlicesInConc = 0;
                {
                    int _iAPos = int(itConc-this->m_asConc.begin());
                    for(;_iAPos < this->m_SlicesToMeasure;_iAPos += int(this->m_lConcatenations))
                    {
                        ++iNSlicesInConc;
                    }
                }
                if( (this->m_iNSweepsPerConc == 2) && (iNSlicesInConc%2 == 0) )
                {
                    iAPos += int(this->m_lConcatenations); 
                }
                int iAPos_incr = int(this->m_lConcatenations)*this->m_iNSweepsPerConc;
                int iSweep = 0;
                for(iCIR = 0; iCIR != this->m_iModulo;++iCIR)
                {
                    if( iAPos < m_SlicesToMeasure )
                    {
                        ++this->m_iNIRPerMeas;
                    }
                    *itChronPos = MAP_APOS_TO_CPOS(rSeries,iAPos);
                    if( *itChronPos++ < m_SlicesToMeasure )
                    {
                        ++m_lKernelRequestsInFirstMeas;
                    }
                    iAPos += iAPos_incr;
                    if( (iAPos > this->m_SlicesToMeasure) && (this->m_iNSweepsPerConc > 1) && (iSweep != this->m_iNSweepsPerConc-1) )
                    {
                        ++iSweep;
                        //  Start position of first sweep
                        iAPos = int(itConc-m_asConc.begin());
                        if( (this->m_iNSweepsPerConc > 2) || (iNSlicesInConc%2 != 0) )
                        {
                            iAPos += iSweep*int(this->m_lConcatenations); 
                        }
                    }
                }
                itChronPos = itConc->m_aiChronPos.begin();
                if( ((*itConc).m_aiChronPos.size() < 1) || (*itChronPos < 0) || (*itChronPos >= rSeries.size()) )
                {
                    setNLSStatus (MRI_SBB_SBB_ERROR);
                    SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                    return false;
                }
            }
        }
        //SBBIRsel.setIdent("TFL");
        //SBBIRsel.setGSWDGradientPerformance(rMrProt,rSeqLim);
        //if( !SBBIRsel.prep(rMrProt, rSeqLim, rSeqExpo) )
        //{
        //    setNLSStatus(SBBIRsel.getNLSStatus());
        //    SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
        //    return false;
        //}
        const long lNPhsCorrScans = m_PhaseCorScans || rMrProt.phaseStabilize() ? 1 : 0;
        m_lKernelRequestsInSecondMeas = (m_lKernelRequestsInFirstMeas *= (rMrProt.averages()*m_LinesToMeasure*m_PartitionsToMeasure*m_FreeLoopLength+m_lPreparingScans+lNPhsCorrScans));
        this->m_iNIRPerMeas                                                 *= int(rMrProt.averages()*m_LinesToMeasure*m_PartitionsToMeasure*m_FreeLoopLength+m_lPreparingScans);
        m_RFInfoInSATs = int(m_lKernelRequestsInFirstMeas+m_RepetitionsToMeasure*m_lKernelRequestsInSecondMeas) *
            ( CSatFat.getRFInfoPerRequest()
            + CSatWater.getRFInfoPerRequest()
            + MSat.getRFInfoPerRequest()
            + RSat1.getRFInfoPerRequest()
            + RSat2.getRFInfoPerRequest()
            + RSat3.getRFInfoPerRequest()
            + RSat4.getRFInfoPerRequest()
            + RSat5.getRFInfoPerRequest()
            + RSat6.getRFInfoPerRequest()
            );
        return true;
    }
    return true;
}

int SEQ_NAMESPACE::SeqLoopFastIR::getTBlock_us() const
{
    if( (this->m_iIRScheme == SEQ::IR_SCHEME_SEQUENTIAL) && (this->m_iModulo > 0) )
    {
        return this->m_iTBlock_us;
    }
    return 0;
}

int SEQ_NAMESPACE::SeqLoopFastIR::getKernelOffset() const
{
    if( this->m_iIRScheme == SEQ::IR_SCHEME_SEQUENTIAL) 
    {
        return this->m_iKernelOffset;
    }
    return 0;
}

int SEQ_NAMESPACE::SeqLoopFastIR::getTRMin_us() const
{
    return int(this->m_lTRMin);
}



bool SEQ_NAMESPACE::SeqLoopFastIR::runKernel(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS asSlc[], sREADOUT* pADC)
{
    try
    {

        if( this->m_iIRScheme == SEQ::IR_SCHEME_SEQUENTIAL )
        {
            //  Chronological blade counter
            int iCExcit = (this->m_iCExe-this->m_iKernelOffset) < m_lPreparingScans*this->m_iModulo
                ? (this->m_iCExe-this->m_iKernelOffset)/this->m_iModulo
                : (this->m_iCExe-this->m_iKernelOffset-int(this->m_lPreparingScans)*this->m_iModulo)/this->m_iModulo
                ;

            //  Total number of blades
            //const long lNExcit    = (this->m_iNExe-this->m_iKernelOffset-m_lPreparingScans*this->m_iModulo)/this->m_iModulo;

            if( (this->m_iCExe-this->m_iKernelOffset) < m_lPreparingScans*this->m_iModulo )
            {
                m_lKernelMode = KERNEL_PREPARE;
                this->m_lPrepareLoopCounter      = iCExcit;
                this->m_lLineCounter             = 0;
                this->m_FreeLoopCounter          = 0;
                this->m_lOuterAcquisitionCounter = 0;

            }
            else
            {
                m_lKernelMode = KERNEL_IMAGE;
                //  Innermost supported counter
                this->m_lLineCounter             = (iCExcit%this->m_LinesToMeasure);
                this->m_FreeLoopCounter          = (iCExcit/this->m_LinesToMeasure)%this->m_FreeLoopLength;
                // outermost supported counter
                this->m_lOuterAcquisitionCounter = (iCExcit/(this->m_LinesToMeasure*this->m_FreeLoopLength))%this->m_AcquisitionsOuter;
            }

            fRTSetReadoutEnable(m_lKernelMode==KERNEL_PREPARE ? 0 : 1);
            MdhProxy& rMDH = pADC->getMDH();
            rMDH.setCacq(static_cast<unsigned short>(this->m_lOuterAcquisitionCounter));
            rMDH.setClin(static_cast<unsigned short>(this->m_lLineCounter));

            m_bIsFirstScanInSlice = (this->m_lLineCounter == 0);
            rMDH.setFirstScanInSlice(m_bIsFirstScanInSlice);

            m_bIsLastScanInSlice = (this->m_lLineCounter == this->m_LinesToMeasure - 1);
            rMDH.setLastScanInSlice (m_bIsLastScanInSlice);

            m_bIsLastScanInConcat = (this->m_iCExe == this->m_iNExe - 1);/*(lCExcit == lNExcit-1)*/
            rMDH.setLastScanInConcat( m_bIsLastScanInConcat );

            m_bIsLastScanInMeas = (this->m_iCExe == this->m_iNExe - 1)/*(lCExcit == lNExcit-1*/ && (this->m_lConcatenationCounter == this->m_lConcatenations - 1);
            rMDH.setLastScanInMeas  (m_bIsLastScanInMeas);

            if( rMDH.isLastScanInMeas() )
            {
                SEQ_TRACE_ERROR.print("Last-Scan-In-Meas %d(%d).",this->m_iCExe,this->m_iNExe);
            }
#ifdef SUPPORT_BLADE
            if( pProt->kSpace().trajectory() == SEQ::TRAJECTORY_BLADE )
            {
                const bool bLastBladeInTR = isLastBladeInTR(pProt);
                if( bLastBladeInTR )
                {
                    rMDH.addToEvalInfoMask(MDH_LAST_BLADE_IN_TR);
                    const int iNExcit    = (this->m_iNExe-this->m_iKernelOffset-int(this->m_lPreparingScans)*this->m_iModulo)/this->m_iModulo;
                    rMDH.setLastScanInConcat( m_bIsLastScanInConcat = (iCExcit == iNExcit-1) );
                    rMDH.setLastScanInMeas  ( m_bIsLastScanInMeas   = (iCExcit == iNExcit-1 && m_lConcatenationCounter == m_lConcatenations-1) );
                }
                else
                {
                    rMDH.deleteFromEvalInfoMask(MDH_LAST_BLADE_IN_TR);
                    rMDH.setLastScanInConcat( m_bIsLastScanInConcat = false);
                    rMDH.setLastScanInMeas  ( m_bIsLastScanInMeas   = false);
                }
                rMDH.setCidd(static_cast<unsigned short>(iCExcit));
            }
#endif // SUPPORT_BLADE
        }
        return BASE_TYPE::runKernel(rMrProt,rSeqLim,rSeqExpo,asSlc,pADC);
    }
    catch(...)
    {
        SEQ_TRACE_ERROR.print("Caught unknown exception at %s(%d).",__FILE__,__LINE__);
        setNLSStatus(MRI_SEQ_SEQU_ERROR);
        return false;
    }
}


bool SEQ_NAMESPACE::SeqLoopFastIR::runConcatenationKernelPre(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS* pSlcPos,sREADOUT* psADC)
{
    //  Reset loop counters
    m_FreeLoopCounter
        = m_lOuterAcquisitionCounter
        = m_lInnerAcquisitionCounter
        = m_lLineCounter
        = m_lPartitionCounter
        = m_l3dECounter
        = m_l3dEAcquisitionCounter
        = m_lOuterSliceCounter
        = m_lInnerSliceNumber
        = m_lInnerSliceCounter
        = m_lPhaseCounter
        = m_lSliceIndex
        = m_lConcIndex
        = m_KernelCallsLoopCounter
        = 0;
    return BASE_TYPE::runConcatenationKernelPre(rMrProt,rSeqLim,rSeqExpo,pSlcPos,psADC);
}

bool SEQ_NAMESPACE::SeqLoopFastIR::runConcatenationLoop(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS aSlcPos[],sREADOUT* pADC)
{
    if( this->m_iIRScheme != SEQ::IR_SCHEME_SEQUENTIAL )
    {
        return BASE_TYPE::runConcatenationLoop(rMrProt,rSeqLim,rSeqExpo,aSlcPos,pADC);
    }
    MdhProxy& rMDH = pADC->getMDH();


    for(m_lConcatenationCounter = m_lSliceOffsetConc = 0; m_lConcatenationCounter < m_lConcatenations; ++m_lConcatenationCounter)
    {
        rMDH.setCacq (0);
        rMDH.setClin (0);
        rMDH.setCpar (0);
        rMDH.setFirstScanInSlice(m_bIsFirstScanInSlice = false);
        rMDH.setLastScanInSlice (m_bIsLastScanInSlice  = false);


        // * ---------------------------------------------------------------------------- 
        // * Play the PATRefScans for the mode SEQ::IR_SCHEME_SEQUENTIAL
        // * This is necessary as the ref scans are called in the runConcatenationLoop() 
        // * method the base class
        // * ---------------------------------------------------------------------------- 

        // which loop mode is active
        if (    (m_ePATRefScanLoopMode == PATRefScanLoopMode_EACH_REPETION)
            || ((m_ePATRefScanLoopMode == PATRefScanLoopMode_ONCE) && (m_lRepetitionCounter == 0)))
        {
            // pointer to current concatenation
            std::vector<CONC>::iterator itConc = m_asConc.begin() + m_lConcatenationCounter;

            // get start and end positions of the vector 
            std::vector<int>::iterator itChronPos = (*itConc).m_aiChronPos.begin();
            std::vector<int>::iterator itChronPosE = (*itConc).m_aiChronPos.end();

            // loop over chronological positions
            for(; itChronPos != itChronPosE; ++itChronPos)
            {
                // there are cases where there are actually less slices than elements in the vector
                if( (*itChronPos) < m_SlicesToMeasure )
                {
                    // tell the ref scan SBB about the slice position
                    SBBPATRefScan.setSlice(*itChronPos);

                    // partitions loop (number of partitions to measure is provided by SBBPATRefScan itself)
                    for( int lPartitionsCounter = 0; lPartitionsCounter < SBBPATRefScan.getlPartitionsToMeasure(); lPartitionsCounter++)
                    {
                        SBBPATRefScan.setPartition(lPartitionsCounter);

                        if( !SBBPATRefScan.run(rMrProt, rSeqLim, rSeqExpo, aSlcPos) )
                        {
                            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                            setNLSStatus ( SBBPATRefScan.getNLSStatus() );
                            return ( false );
                        }
                    }
                }
            }
        }

        if( !runConc(rMrProt,rSeqLim,rSeqExpo,aSlcPos,pADC) )
        {
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }

        //  Execute delay time past concatenation
        if( m_TDFill > 0 && (m_lConcatenationCounter < m_lConcatenations-1) )
        {
            setNLSStatus( fSBBFillTimeRun(m_TDFill) );
            if( (m_NLSStatus & NLS_SEV) != NLS_SUCCESS )
            {
                SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                return false;
            }
        }
    }
    return true;
}

bool SEQ_NAMESPACE::SeqLoopFastIR::runConc(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS aSlcPos[],sREADOUT* pADC)
{
    try
    {
        if( this->m_iIRScheme == SEQ::IR_SCHEME_SEQUENTIAL )
        {
            for(this->m_iCExe = 0; this->m_iCExe != this->m_iNExe; ++this->m_iCExe)
            {
                if( !runBloc_uniform(rMrProt,rSeqLim,rSeqExpo,aSlcPos,pADC) )
                {
                    SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                    return false;
                }
            }
        }
        else
        {
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            setNLSStatus(MRI_SEQ_SEQU_ERROR);
            return false;
        }
        return true;
    }
    catch(...)
    {
        SEQ_TRACE_ERROR.print("Caught unknown exception at %s(%d).",__FILE__,__LINE__);
        setNLSStatus(MRI_SEQ_SEQU_ERROR);
        return false;
    }
}

bool SEQ_NAMESPACE::SeqLoopFastIR::check(MrProt& rMrProt,SeqLim &rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* pADC)
{
    if( this->m_iIRScheme != SEQ::IR_SCHEME_SEQUENTIAL )
    {
        return BASE_TYPE::check(rMrProt,rSeqLim,rSeqExpo,pSlcPos,pADC);
    }
    //  SeqLoop::check skips the free loop, Therefore we set the free loop length to 1 and the line/partition loop to 'linesToCheck', 'partitionsToCheck'.
    int iNExe_run = this->m_iNExe;
    this->m_iNExe = int(this->m_seqCheck.m_NoLines*1*this->m_AcquisitionsOuter+this->m_lPreparingScans)*this->m_iModulo+this->m_iKernelOffset;

    MdhProxy& rMDH = pADC->getMDH();
    for(m_lConcatenationCounter = m_lSliceOffsetConc = 0; m_lConcatenationCounter < m_lConcatenations; ++m_lConcatenationCounter)
    {
        rMDH.setCacq (0);
        rMDH.setClin (0);
        rMDH.setCpar (0);
        rMDH.setFirstScanInSlice(m_bIsFirstScanInSlice = false);
        rMDH.setLastScanInSlice (m_bIsLastScanInSlice  = false);

        if( !runConc(rMrProt,rSeqLim,rSeqExpo,pSlcPos,pADC) )
        {
            this->m_iNExe = iNExe_run;
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }

        //  Excute delay time past concatenation
        if( m_TDFill > 0 && (m_lConcatenationCounter < m_lConcatenations-1) )
        {
            setNLSStatus( fSBBFillTimeRun(m_TDFill) );
            if( (m_NLSStatus & NLS_SEV) != NLS_SUCCESS )
            {
                this->m_iNExe = iNExe_run;
                SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
                return false;
            }
        }
    }
    this->m_iNExe = iNExe_run;

    return true;
}

#ifdef WIN32
bool SEQ_NAMESPACE::SeqLoopFastIR::run_new (MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS* pSlcPos, sREADOUT* psADC)
{
    if( this->m_iIRScheme == SEQ::IR_SCHEME_SEQUENTIAL )
    {
        FileNameExp  sFileNameExp;
        //  Chronological log
        const char* pszFName = sFileNameExp.Translate("%MEASDAT%/IIR_CHRON.txt");
        this->m_sChron.open(pszFName,std::ios::out);
        if( !this->m_sChron.is_open() )
        {
            SEQ_TRACE_WARN.print("Warning at %s(%d): Failed to open file '%s'.",__FILE__,__LINE__,pszFName);
        }
    }
    bool bRet = BASE_TYPE::run_new(rMrProt,rSeqLim,rSeqExpo,pSlcPos,psADC);
    if( this->m_sChron.is_open() )
    {
        this->m_sChron.close();
    }
    return bRet;
}
#endif

bool SEQ_NAMESPACE::SeqLoopFastIR::isLastBladeInTR(const MrProt &) const
{
    if( this->m_iIRScheme == SEQ::IR_SCHEME_SEQUENTIAL )
    {
        if( size_t(m_lConcatenationCounter) >= m_asConc.size() )
        {
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }
        //  Pointer to current concatenation
        std::vector<CONC>::const_iterator itConc = m_asConc.begin()+m_lConcatenationCounter;

        const int iCKernel = (this->m_iCExe-this->m_iKernelOffset)%this->m_iModulo;

        if( iCKernel < 0 || (size_t(iCKernel) >= (*itConc).m_aiChronPos.size()) )
        {
            SEQ_TRACE_ERROR.print("Error at %s(%d).",__FILE__,__LINE__);
            return false;
        }
        //m_bFirstSliceInConcat = lCKernel==0;
        const bool bLastBladeInTR = iCKernel == int((*itConc).m_aiChronPos.size())-1
            || (*itConc).m_aiChronPos[iCKernel+1] == m_SlicesToMeasure
            ;
        return bLastBladeInTR;

    }
#ifdef SUPPORT_BLADE
    return BASE_TYPE::isLastBladeInTR(rMrProt);
#else
    const bool bLastBladeInTR = (m_lKernelMode == KERNEL_IMAGE || m_lKernelMode == KERNEL_PHASECOR)
        && (m_lOuterSliceCounter == m_SeqConcat[m_lConcatenationCounter].m_LoopsOuter-1)
        && (m_lPhaseCounter == m_PhasesToMeasure-1)
        && (m_lInnerSliceCounter == m_lInnerSliceNumber-1)
        && (m_KernelCallsLoopCounter == m_KernelCallsLoopLength-1)
        ;
    return bLastBladeInTR;
#endif
}
