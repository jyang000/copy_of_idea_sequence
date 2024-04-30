//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4\pkg\MrServers\MrImaging\seq\a_ep_seg.h
//	 Version:
//	  Author: Clinical
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: Deklarations for a_ep_seg.cpp
//
//
/// \brief  File containing declarations for the sequences 
///         - ep_seg_fid
///         - ep_seg_se
///         
/// This file contains the declaration of the class ep_seg.
/// The sequence ep_seg_fid and ep_seg_se use it to generate very nice images.
/// 
//  ***************************************************************************/
/// 

#pragma once

//-------------------------------------------------------------------------------------
// define sequence variant
//-------------------------------------------------------------------------------------
#ifdef EP_SEG_FID
#define SEQUENCE_VARIANT_DEFINED_a_ep_seg
#define NEED_a_ep_seg_fConfigureExcitationForEPIKernelFID
#define NEED_a_ep_seg_fConfigureExcitationForEPIKernelFID_REGISTERED_IN_INIT
#define NEED_a_ep_seg_calculateBasicScanTime
#define NEED_a_ep_seg_calculateTRTIFillTimes
#define NEED_a_ep_seg_calculateSatEtcScanTime
#define NEED_a_ep_seg_RFSpoiling
#define NEED_a_ep_seg_RestrictedPhasePF
#define NEED_a_ep_seg_3D
#define NEED_a_ep_seg_EarlyFIDPhaseCorrection_IN_INIT
#define NEED_a_ep_seg_SetMaxSliceThicknessForHighGain
#define NEED_a_ep_seg_FIDLikeSpoilerSettings
#define NEED_a_ep_seg_RestrictPhasePartialFourierForEPIFactor3
#define NEED_a_ep_seg_fSEQRunKernel
#endif

#ifdef EP_SEG_SE
#define SEQUENCE_VARIANT_DEFINED_a_ep_seg
#define NEED_a_ep_seg_fConfigureExcitationForEPIKernelSE
#define NEED_a_ep_seg_fConfigureExcitationForEPIKernelSE_REGISTERED_IN_INIT
#define NEED_a_ep_seg_calculateBasicScanTime
#define NEED_a_ep_seg_calculateSatEtcScanTime
#define NEED_a_ep_seg_calculateTRTIFillTimes
#define NEED_a_ep_seg_SEPlugInForEPIKernel
#define NEED_a_ep_seg_SEPlugInForEPIKernel_REGISTERED_IN_INIT
#define NEED_a_ep_seg_Inversion
#define NEED_a_ep_seg_EffectiveTRForSats
#define NEED_a_ep_seg_3D
#define NEED_a_ep_seg_SetMaxSliceThicknessForHighGain
#define NEED_a_ep_seg_SELikeSpoilerSettings
#define NEED_a_ep_seg_RestrictPhasePartialFourierForEPIFactor3
#define NEED_a_ep_seg_fSEQRunKernel
#endif



//  -------------------------------------------------------------------------- 
//  General Includes                                                           
//  -------------------------------------------------------------------------- 
#include "MrImagingFW/libSBBFW/StdSeqIF.h"
#ifdef SUPPORT_PACE
#include  "MrImaging/seq/common/libPace/SeqLoopBH.h"        
#else
#include  "MrImaging/libSBB/SEQLoop.h"
#endif
#include      "MrImaging/libSBB/SBBBinomialPulses.h"
#include      "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include      "MrImaging/seq/SystemProperties.h"
//#include      "MrImaging/seq/epi_StdUILink.h"
#include      "MrImaging/libSeqUtil/ReorderInfoEPI.h"

#ifdef NEED_a_ep_seg_RFSpoiling
#include      "MrImagingFW/libSeqUtilFW/RFSpoiling.h"
#endif


namespace SEQ_NAMESPACE
{
    class Ep_seg;
}//end of namespace SEQ_NAMESPACE


#ifdef NEED_a_ep_seg_SEPlugInForEPIKernel
#include "MrImaging/seq/a_ep2d_se/SBBRefocSE.h"
#include "MrImaging/seq/a_ep_seg_se/SBBEPIsegKernelSE.h"
#else
#include "MrImaging/seq/a_ep_seg_fid/SBBEPIsegKernelFID.h"
#endif   



//  -------------------------------------------------------------------------- 
//  Application includes                                                       
//  -------------------------------------------------------------------------- 
#include "MrImaging/seq/a_ep_seg_UI.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"


//  -------------------------------------------------------------------------- 
//  Forward declarations                                                       
//  -------------------------------------------------------------------------- 
class MrProt;
class SeqLim;
class SeqExpo;

namespace MrMeasSrv
{
    class ISequence;
}



namespace SEQ_NAMESPACE
{
    class EpCommonUI;

#ifdef WIN32
    //  ----------------------------------------------------------------------
    //
    //  Name        :  getUI
    //
    //  Description :
    /// \brief         Returns the pointer to the UI class
    ///
    //
    //  Return      :  EpCommonUI*
    //
    //  ----------------------------------------------------------------------
    Ep_segUI* getUI(MrUILinkBase* const pThis);

    //  ----------------------------------------------------------------------
    //
    //  Name        :  getUI
    //
    //  Description :
    /// \brief         Returns the pointer to the UI class
    ///
    //
    //  Return      :  EpCommonUI*
    //
    //  ----------------------------------------------------------------------
    Ep_segUI* getUI(MrMeasSrv::ISequence* const pSeq);

    //  ----------------------------------------------------------------------
    //
    //  Name        :  getSeq
    //
    //  Description :
    /// \brief         Returns the pointer to the sequence Ep2d
    //
    //  Return      :  Ep2d*
    //
    //  ----------------------------------------------------------------------
    Ep_seg* getSeq(MrUILinkBase* const pThis);

#endif // #ifdef WIN32

    //  -------------------------------------------------------------------------- 
    //
    /// \brief <b> Class definition of Ep_seg. This class is used by the sequences
    ///         - Ep_seg
    ///         - Medic1    </b>
    ///         
    /// This file contains the declaration of the class Ep_seg.
    /// The sequence Ep_seg and Medic1 use it to generate very nice images.
    /// 
    //  -------------------------------------------------------------------------- 
    class __IMP_EXP Ep_seg: public StdSeqIF
    {

    public:

        KernelCalculationLimits     m_myCalcLimits;                            // kernel calculation limits

        //  ------------------------------------------------------------------ 
        //                                                                     
        //  Name        :  Ep_seg::Ep_seg                                        
        //                                                                     
        //  Description :  
        /// \brief         Initialization of class members                     
        //                                                                     
        //  Return      :  void                                                
        //                                                                     
        //  ------------------------------------------------------------------ 
        Ep_seg();



        //  ------------------------------------------------------------------ 
        //                                                                     
        //  Name        :  Ep_seg::~Ep_seg                                       
        //                                                                     
        //  Description :  
        /// \brief         Destructor. Deletes existing MedicUI instances.                                          
        //                                                                     
        //  Return      :  void                                                
        //                                                                     
        //  ------------------------------------------------------------------ 
        virtual ~Ep_seg();


        Ep_seg(const Ep_seg& right) = delete;
        Ep_seg& operator=(const Ep_seg& right) = delete;
        Ep_seg(Ep_seg&& right)                 = delete;
        Ep_seg& operator=(Ep_seg&& right) = delete;

        //   --------------------------------------------------------------------------  
        //                                                                               
        //   Name        :  Ep_seg::initialize                                            
        //                                                                               
        //   Description :  
        ///  \brief        Initialization of the sequence                               
        ///                                                                               
        ///                On the host, the object m_pUI will actually contain sensible 
        ///                  data after Ep_seg::initialize. On the measurement system, it
        ///                  is basically an empty object behind it. 
        ///                                                                               
        //   Return      :  NLS status                                                   
        //                                                                               
        //   --------------------------------------------------------------------------  
        virtual NLSStatus initialize(SeqLim & pSeqLim);



        //  -------------------------------------------------------------------------- 
        //                                                                             
        //  Name        :  Ep_seg::prepare                                              
        //                                                                             
        //  Description :  
        /// \brief <b>     Preparation of the sequence during binary search and prior  
        ///                 to sequence execution  </b>                                    
        //                                                                             
        //  Return      :  NLS status                                                  
        //                                                                             
        //  -------------------------------------------------------------------------- 
        virtual NLSStatus prepare(MrProt &pMrProt, SeqLim &pSeqLim, SeqExpo &pSeqExpo);



#ifdef EP_SEG_FID
        // --------------------------------------------------------------------------
        //
        //  Name        :  Ep_seg::prePrepare
        //
        //  Description :
        /// \brief <b>     Check on protocol, if there are forbidden combinations of  parameters.
        ///                This speeds up the binary search, because the method is called very early in
        ///                the prepare process and therefore can save lots of lengthy calculations.
        ///                The method must not prepare any objects like pulses or SBBs.</b>
        //
        //  Return      :  NLS status
        //
        //  --------------------------------------------------------------------------
        virtual NLSStatus prePrepare(const MrProt &rMrProt, const SeqLim &rSeqLim, SeqExpo &rSeqExpo);
#endif // EP_SEG_FID



        //  -------------------------------------------------------------------------- 
        //                                                                             
        //  Name        :  Ep_seg::check                                                
        //                                                                             
        //  Description :  
        /// \brief  <b>    Check of the sequence for gradient stimulation </b>     
        ///                                                                             
        ///                This method is called by the framework prior to a 
        ///                 measurement on the host to ensure, that 
        ///                 - no gradient overflow occurs                                                                             
        ///                 - the stimulation will not exceed the threshold  
        ///
        //  Return      :  NLS status                                                  
        //                                                                             
        //  -------------------------------------------------------------------------- 
        virtual NLSStatus check(MrProt & pMrProt, SeqLim & pSeqLim, SeqExpo & pSeqExpo, SEQCheckMode *  pSEQCheckMode);


        //  -------------------------------------------------------------------------- 
        //                                                                             
        //  Name        :  Ep_seg::run                                                  
        //                                                                             
        //  Description :                                                              
        ///     \brief     Execution of the sequence                                   
        //                                                                             
        //  Return      :  NLS status                                                  
        //                                                                             
        //  -------------------------------------------------------------------------- 
        virtual NLSStatus run(MrProt & pMrProt, SeqLim & pSeqLim, SeqExpo & pSeqExpo);



        //   -------------------------------------------------------------------------- 
        //                                                                              
        //   Name        :  Ep_seg::runKernel                                            
        //                                                                              
        //   Description :  
        ///  \brief <b>     Executes the basic timing of the real-time sequence.   </b>     
        ///                                                                              
        ///                 The method runKernel plays out a sequence "Kernel",
        ///                  consisting of one or more lines in k-Space.
        ///                 It is called by SeqLoop.
        ///                                                                              
        //   Return      :  NLS status                                                  
        //                                                                              
        //   -------------------------------------------------------------------------- 
        virtual NLS_STATUS runKernel(MrProt &pMrProt, SeqLim &pSeqLim, SeqExpo &pSeqExpo, long lKernelMode, long lSlice, long lPartition, long lLine);

#ifdef SUPPORT_PACE
        // ------------------------------------------------------------------------------
        // Function    : fSEQReceive_a_ep_seg
        // ------------------------------------------------------------------------------
        //               
        // Description : Receives data from host/IR at run-time.
        //
        // Return      : A NLSStatus code.
        //
        // ------------------------------------------------------------------------------
        NLSStatus receive(SeqLim *pSeqLim, SeqExpo *pSeqExpo, const SEQData &rSEQData) override;

        /*[ Function ****************************************************************\
        *
        * Name        : fSEQConvProt
        *
        * Description : - converts old PACE protocols:
        * Return      : MRI_SEQ_SEQU_NORMAL for success
        *               MRI_SEQ_SEQU_ERROR  for error
        *
        \****************************************************************************/

#ifdef WIN32
        NLSStatus convProt(const MrProt &rMrProtSrc, MrProt &rMrProtDst);
#endif

#endif  //  SUPPORT_PACE        

        //  -------------------------------------------------------------- 
        //                                                                 
        //  Name        :  getUI                                           
        //                                                                 
        //  Description :  
        /// \brief <b>     Returns the pointer to the Ep_seg UI class  </b>     
        ///                                                                 
        ///                This method is only sensible on the host.
        ///                On the measurement system, it will return an nearly empty object.                                                 
        ///                                                                 
        //  Return      :  MedicUI*                                        
        //                                                                 
        //  -------------------------------------------------------------- 
        Ep_segUI* getUI(void) const;

#ifdef NEED_a_ep_seg_calculateBasicScanTime
        // ------------------------------------------------------------------------------
        // Function    : calculateBasicScanTime
        // ------------------------------------------------------------------------------
        //               
        // Description : Calculates basic scan time of the sequence kernel for function
        //               calculateTRTIFillTimes, i.e. the time from the beginning of the
        //               RTEB containing the excitation RF until the end of the kernel.
        //
        // Return      : long
        //
        // ------------------------------------------------------------------------------
        long calculateBasicScanTime();

#endif

#ifdef NEED_a_ep_seg_calculateSatEtcScanTime
        // ------------------------------------------------------------------------------
        // Function    : calculateSatEtcScanTime
        // ------------------------------------------------------------------------------
        //               
        // Description : Calculates the time portion of the run Kernel before the RTEB
        //               containing the excitation RF.
        //  
        // Return      : long
        //
        // ------------------------------------------------------------------------------
        long calculateSatEtcScanTime();

#endif

#ifdef NEED_a_ep_seg_calculateTRTIFillTimes
        // ------------------------------------------------------------------------------
        // Function    : calculateTRTIFillTimes
        // ------------------------------------------------------------------------------
        //               
        // Description : Calculates TI and TR needed to execute the measurement with the
        //               current protocol and the according fill times for the sequence
        //               timing.
        //
        // Return      : true (if success) or false
        //
        // ------------------------------------------------------------------------------
        bool calculateTRTIFillTimes
            (
            MrProt&  pMrProt,
            SeqLim&  pSeqLim,
            SeqExpo& pSeqExpo,
            long*    plNeededTI,
            long*    plNeededTR
            );
#endif


#ifdef EP_SEG_FID
        ReorderInfoEPI * getREOInfoPointer();
#endif // EP_SEG_FID

    protected:

        double                      m_dMaxAllowedDataRate{1000.0};                     // maximum allowed data rate for online regridding
        double                      m_dMaxSliceThicknessForHighGain{1000.0};           // controls high/low-gain-switching
        sSLICE_POS                  m_asSLC[K_NO_SLI_MAX];                     // Slice position information (rotation matrices and shifts)
#ifdef SUPPORT_PACE
        PACE::SeqLoopBH             m_mySeqLoop;                               // standard loop structure augmented by multiple breath-hold feature
#else
        SeqLoop                     m_mySeqLoop;                               // standard loop structure class
#endif
#ifdef EP_SEG_FID 
        SBBEPIsegKernelFID m_EPIKernel{nullptr}; // the EPI kernel class
#elif defined EP_SEG_SE
        SBBEPIsegKernelSE           m_EPIKernel{nullptr};                               // the EPI kernel class
#endif
        SeqBuildBlockSpoilGrad      m_SpoilBeforeEPIKernel{nullptr};                    // the spoiler before the EPI kernel
        SeqBuildBlockSpoilGrad      m_SpoilAfterEPIKernel{nullptr};                     // the spoiler after the EPI kernel
        ReorderInfoEPI              m_REOInfo;                                 // reordering data storage class

#ifdef NEED_a_ep_seg_RFSpoiling
        RFSpoiling                  m_myRFSpoil;                               // class handling RF-spoiling data
#endif


        //-------------------------------------------------------------------------------------
        // sequence hint text - text variable declaration moved from fSEQInit() - CHARM 312244
        //-------------------------------------------------------------------------------------
#ifdef SHOW_LOOP_STRUCTURE
        char m_tLoopText[200]{};
#endif
        char m_tHintText[512]{};
        //-------------------------------------------------------------------------------------
        // define sequence variant
        //-------------------------------------------------------------------------------------
        const char *ptVariant;

        //  -------------------------------------------------------------- 
        /// \brief <b> UI class for Ep_seg                                             
        ///
        ///         This class is basically empty on the measurement system
        //  -------------------------------------------------------------- 
        Ep_segUI* m_pUI{nullptr};




        //  ------------------------------------------------------------------ 
        //                                                                     
        //  Name        :  Ep_seg::createUI                                     
        //                                                                     
        //  Description :  
        /// \brief <b>     Instantiation of UI classes   </b>                      
        //                                                                     
        //  Return      :  NLS status                                          
        //                                                                     
        //  ------------------------------------------------------------------ 
        virtual NLS_STATUS createUI(SeqLim &pSeqLim);


    private:

#ifdef EP_SEG_FID
        long m_lKSpaceCenterSegment{0};
#endif // EP_SEG_FID
    };
};

#ifdef EP_SEG_FID
    inline ReorderInfoEPI * SEQ_NAMESPACE::Ep_seg::getREOInfoPointer()
    {
        return &m_REOInfo;
    }
#endif // EP_SEG_FID
