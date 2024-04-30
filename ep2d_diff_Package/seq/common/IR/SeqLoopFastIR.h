//  -----------------------------------------------------------------
//  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential
//  -----------------------------------------------------------------
//
//          Project: NUMARIS/4
//             File: \n4\pkg\MrServers\MrImaging\seq\common\IR\SeqLoopFastIR.h
//          Version: \main\2
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

//  Definition of base class
#ifdef SUPPORT_BLADE
#include "MrImaging/seq/a_tse/SeqLoopBLADE.h"
typedef SeqLoopBLADE SeqLoopFastIR_BASE_TYPE;
#elif COMPILE_EP2D_DIFF
#include "MrImaging/seq/common/SeqLoopLongTRTrig/SeqLoopLongTRTrig.h"
typedef SeqLoopLongTRTrig  SeqLoopFastIR_BASE_TYPE;
#else
#include  "MrImaging/libSBB/SEQLoop.h"
typedef SeqLoop  SeqLoopFastIR_BASE_TYPE;
#endif

#include "MrProtSrv/Domain/MrProtData/MrProt/ProtDefines.h"

//  SEQ::IRScheme
#include "MrProtSrv/Domain/CoreNative/PrepPulsesDefines.h" 

#ifndef __SeqLoopFastIR_H
#define __SeqLoopFastIR_H

#include <vector>
#ifdef WIN32
#include <fstream>
#endif
//  -----------------------------------------------------------------
//  Import - Export - Control
//  -----------------------------------------------------------------

//  -----------------------------------------------------------------
//  Definition of class SeqLoopIIR
//  -----------------------------------------------------------------
namespace SEQ_NAMESPACE
{ 
    //  Calculates the duration of a single block, if IRScheme SEQUENTIAL is used
    bool fCalcTBlock(int& riTBlock_us, int iTI_us, int iIRTime_us, int iSBBScanTime_us, int iKernelTime_us, int iMinFillTime_us, int iKernelOffset);

    //  Tests wether a particular Block time can be realized for a given TI time and kernel offset. 
    bool fTestTBlock(int iTBlock_us, int iTI_us, int iIRTime_us, int iSBBScanTime_us, int iKernelTime_us, int iMinFillTime_us, int iKernelOffset);

     //  Definition of SeqLoopFastIR
    class SeqLoopFastIR : public SeqLoopFastIR_BASE_TYPE
    {
    public:



        //  -------------------------------------------------------------
        //  Some typedef and local definitions
        //  -------------------------------------------------------------

        typedef SeqLoopFastIR            MYTYPE;
        typedef SeqLoopFastIR_BASE_TYPE  BASE_TYPE;

        //  -------------------------------------------------------------
        //  The interface of the class
        //  -------------------------------------------------------------

        //  Default Constructor
        SeqLoopFastIR();

        //  Destructor
        virtual ~SeqLoopFastIR();

        //  -------------------------------------------------------------
        //  The static member functions of the class.
        //  -------------------------------------------------------------


        //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //  Preparation of the class
        //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        virtual bool prep (MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo);


        //  Preps IR position, if LTIIIR is used.
        //  Must be invoked after SeqLoop::prep
        virtual bool prepConc(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo);

        virtual bool calcFillTimesOnly (MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo& , long lWantedTR = -1, long lTimeForOsc = -1);
        virtual void calcMeasurementTimeUsec (MrProt &rMrProt, SeqLim &rSeqLim);
        virtual void calcNoOfRelevantADCs (MrProt &rMrProt, SeqLim &rSeqLim);

        //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //  The member functions in the following block can be invoked
        //  after the class has been prepared (prep + TrTiFillTimes/calcFillTimes).
        //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        //  Returns the duration of one bloc, if selective IR pulse is used and IRScheme SEQUENTIAL is selected.
        //  Otherwise the return value is undefined.
        int getTBlock_us() const;

        //  Returns the number of IR pulse kernel pairs  between the IR pulse and the Kernel of a particular slice,
        //  if selective IR pulse is used and IRScheme SEQUENTIAL is selected
        //   //  Otherwise the return value is undefined.
        int getKernelOffset() const;


        //  Minimum TR
        int getTRMin_us() const;

        virtual bool check                    (MrProt &rMrProt,SeqLim &rSeqLim,MrProtocolData::SeqExpo &rSeqExpo,sSLICE_POS* pSlcPos,sREADOUT* psADC);

        //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //  Member function invoked at run-time by the client sequence.
        //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        virtual bool runKernel(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC);

        //  The return value is undefined, unless the member function is invoked
        //  by fSEQRunKernel and unless BLADE is used
        virtual bool isLastBladeInTR(const MrProt &rMrProt) const;


#ifdef WIN32
        virtual bool run_new (MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS* pSlcPos, sREADOUT* psADC);
#endif

        //  Number of sweeps through the stack of slices assigned to the current concatenation.
        //  Default is 2, if interleaved and Number of concatenations == 1
        //             1, otherwise
        int getNSweepsPerConc() const;
        void setNSweepsPerConc(int iVal);

        //  Relevant only if IRSCheme SEQUENTIAL is active
        //  Specifies minimum fill time between two Kernel modules without gradient switching.
        //  If specified the specified CoolPause must not be executed by the Kernel
        int  getCoolPauseWithinKernelTime_us() const;
        void setCoolPauseWithinKernelTime_us(int iVal);
        //  -------------------------------------------------------------
        //  The protected member functions of the class
        //  -------------------------------------------------------------

    protected:

        //  Actual Interleaved IR scheme
        SEQ::IRScheme  m_iIRScheme;

        //  Execution counter of current concatenation
        int   m_iCExe;

        //  Number of time a basic bloc is repeated (per concatenation)
        int   m_iNExe;

        //  Number of blocs excuted until the slice series repeats
        int   m_iModulo;

        //  Offset (in number of blocs) between a particular slice and
        //  its IR-pulse
        int   m_iKernelOffset;

        //  Fill-time inserted past IR-pulse
        int   m_iTFillPastIR_us;

        //  Total number of IR pulses per meas (needed for energy calculation)
        int   m_iNIRPerMeas;

        //  Fill time inserted past kernel
        int   m_iTFillPastKernel_us;

        //  Time between two successive inversion recovery pulses
        int   m_iTBlock_us;

        //  Cooling pause
        int   m_iCoolPauseWithinKernelTime_us;

        //  Concatenation specific information
        struct CONC
        {
            //  The length of the array is equal to or less than m_lModulo and stores the
            //  position of the IR-pulse.
            //std::vector<sSLICE_POS>  m_asIRPos;


            //  The length of the array is equal to or less than m_lNSlcPerIR*m_asIRPos.getlSize()
            //  and contains the chronlogical slice
            std::vector<int>        m_aiChronPos;
        };

        std::vector<CONC>  m_asConc;

        //  Only relevant for two sweep veriant of blockwise scheme
        //  End position of first sweep through slice stack
        int               m_iEChronPos_1stSweep;

        //  Default: 2, if interleaved & # conc == 1
        //           1, if interleaved & # conc > 1
        int                m_iNSweepsPerConc;


        //  The following member functions are usesd by TrTiFillTimes to realize
        //  a particular IR-scheme
        //  Additional input prameters are
        //
        //  m_lSBBScanTime ... Time between end of TI-Fill and true kernel (lTIMinAdd1)
        //  m_lEffScanTime ... Scan time required by one kernel call + Sats + Spoilers


        virtual bool calcFillTimes_uniform(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,long lWantedTR, long lTimeForOsc, long* plNegativeFillTime);

#ifdef WIN32
        std::ofstream m_sChron;
#endif
        virtual bool runConcatenationKernelPre(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS* pSlcPos,sREADOUT* psADC);

        virtual bool runConcatenationLoop(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS* pSlcPos,sREADOUT* psADC);

        //  executes one concatenation
        virtual bool runConc(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS* pSlcPos,sREADOUT* psADC);

        //  executes one basic bloc in mode IIR_SCHEME_UNIFORM
        virtual bool runBloc_uniform(MrProt &rMrProt,SeqLim &rSeqLim,SeqExpo &rSeqExpo,sSLICE_POS* pSlcPos,sREADOUT* psADC);
    };
} // namespace
#endif //  __SeqLoopFastIR_H
