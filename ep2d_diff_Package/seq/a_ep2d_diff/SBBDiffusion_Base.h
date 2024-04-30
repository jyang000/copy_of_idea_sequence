//-----------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2021. All Rights Reserved.
//-----------------------------------------------------------------------------

// double include protection:
#pragma once

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------
// MrProt
#include "MrProtSrv/Domain/MrProtData/MrProt/Filter/MrFilter.h"
// MrProt

#include "MrImaging/seq/a_ep2d_diff/didi.h"
#include "MrImagingFW/libBalance/GPABalance.h"   // GPA balance model
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"
#include "MrMeasSrv/SeqIF/libRT/sREADOUT.h"    // for MDH access ...
#include "MrMeasSrv/SeqIF/libRT/sSLICE_POS.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/MeasUtils/MeasMath.h"             // M_PI definition

// SMS support
#include "MrImaging/libSBB/SBBMultibandRF.h"
// SMS support

#include "MrImaging/seq/a_ep2d_diff/SBBCompGrad.h" // compensation gradients

#ifdef BUILD_SEQU
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"     // import/export control
#include "MrImagingFW/libSBBFW/SeqBuildBlock.h"
#include "MrImagingFW/libSBBFW/SBBList.h"
#include "MrImagingFW/libSBBFW/SliceAdjSupport.h"
#include "MrImaging/seq/a_ep2d_diff/DiffusionRFPulseProperties.h"
#include "MrImaging/seq/a_ep2d_diff/SequenceDebugSettings.h"

#include "MrImaging/seq/a_ep2d_diff/DiffusionOrdering.h"
#include "MrProtSrv/Domain/CoreNative/MrApplication.h"

// ---------------------------------------------------------------------------
// Defines
// ---------------------------------------------------------------------------

/// If enabled, very verbose cout's are active to debug this source code
#define DEBUG_DIFFUSION 0

/// If a lower b value is selected, a spoiler gradient will be played out.
#define SPOILER_THRESHOLD               50

// increment for IVIM scans
#define IVIM_B_VALUE_INCREMENT          10

// max value for small increments for IVIM
#define IVIM_B_VALUE_MAX_SMALL_INC      200

/// In orthogonal mode, measurements with a b value up to this will be considered isotropic.
/** changing the define above requires a redesgin: You must then export the number
of b values which are represented by a one-for-three image in ortho mode. */
#define ONE_FOR_THREE_THRESHOLD         0

/// Limits and default values
#define MDDW_DEF_DIRECTIONS             6       // MDDW:    Default number of directions
#define QSPACE_MIN_STEPS                3       // Q-SPACE: Minimum number of b-values
#define QSPACE_MAX_STEPS                5       // Q-SPACE: Maximum number of b-values
#define QSPACE_DEF_STEPS                4       // Q-SPACE: Default number of b-values

/// Number of averages for adjustment scans
/** By using values > 1, this allows to increase the precision of the adjustment scans.
Could be linked to protocol parameters (not implemented yet).
*/
#define ADJ_PREP_SCAN_AVERAGES          2

// --------------------------------------------------------------------------
// Forward declaration
// --------------------------------------------------------------------------
class SBBTestHelperDiffusion;      // SBB iTest Helper Class to Dump Class Member

namespace SEQ_NAMESPACE
{
    // ---------------------------------------------------------------------------
    // Type definitions
    // ---------------------------------------------------------------------------


    /// This enumeration provides the different sequence schemes
    /** Each sequence scheme can be directly mapped to one class implementing the
    corresponding diffusion mode.
    The first value has index 1 because this enumeration can also be used for
    a WIP parameter.
    */
    enum EnumDiffusionScheme {
        DiffusionSchemeNone = 1,
        DiffusionSchemeBipolar,
        DiffusionSchemeTrace,
        DiffusionSchemeStejskal,
        DiffusionSchemeStejskalPlus,
        DiffusionSchemeBipolarPlus,
        DiffusionSchemeSTEAM,
    };


#ifdef WIP

    /// This enumeration specifies the additional WIP parameters
    //  (and thus also the position on the Sequence/Special card)
    enum EnumDiffSpecialParameters {
        WIP_Anything = 1,
    };

#endif   // of #if WIP


    // ---------------------------------------------------------------------------
    // Function definitions
    // ---------------------------------------------------------------------------


    // ===========================================================================
    class __IMP_EXP SBBDiffusion_Base: public SeqBuildBlock
        // ===========================================================================

    {
    public:
        SBBDiffusion_Base(SBBList* pSBBList = nullptr);
        
        virtual ~SBBDiffusion_Base() = default;

        // No copy constructor, copy assignment, move constructor, or move assignment
        SBBDiffusion_Base(SBBDiffusion_Base const&) = delete;

        SBBDiffusion_Base(SBBDiffusion_Base&&) = delete;

        SBBDiffusion_Base& operator=(SBBDiffusion_Base const&) = delete;

        SBBDiffusion_Base& operator=(SBBDiffusion_Base&&) = delete;


        friend class SBBTestHelperDiffusion;      // SBB iTest Helper Class to Dump Class Member

        ///  Prepare timing and real time events
        /**
        \b Input:
        \n rMrProt, rSeqLim, rSeqExpo

        \b Output:
        \n m_lActalTE, m_RFInfoPerRequest, m_lPreEchoTimeContrib,
        m_lPostEchoTimeContrib, m_lSBBStorageTimePerRequest_us, m_bIsMagnetizationInverted,
        m_bPrepared. Real time events are prepared after calling.

        \b Return value:
        \n Success = true, failure = false

        The central part of this method consists of the TE search strategy. It identifies
        the minimum TE that is necessary to realize the desired maximum b-value of the
        actual protocol. This TE is exported indirectly by setting the values for
        m_lPreEchoTimeContrib and m_lPostEchoTimeContrib.

        Derived diffusion SBBs must provide an implementation of the pure virtual methods
        ::prepInit(), ::prepTiming() and ::prepFinish (see below):

        - ::prepInit() is called at the beginning.
        - ::prepTiming() is called within the TE search procedure
        - ::prepFinish() is called at the end

        */
        bool prepSBB(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo) override;


        ///  The pure virtual 'run'  method plays out the diffusion gradients.
        /**
        Note that setAdjustmentScan() and setLoopCounters() have to be called
        in advance to set the loop counters (member variables) to the actual values.
        */
        bool runSBB(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo, sSLICE_POS* pSLC) override = 0;


        virtual bool   setADCforDiffusionMDHentries(sREADOUT* pADC);
        /*!< This needs to be called before the run method: the correct values for the diffusion related MDH loop counters are set. */


        virtual void   setADCusTillEcho(long lADCusTillEcho);
        virtual void   setSpinPrepTimeus(long lSpinPrepTimeus);
        /// Get total number of diffusion encoding steps
        /**
        \b Input:
        \n m_vdBValues, m_lDirections, bIncludeLocalAverages

        \b Output:
        \n n.a.

        \b Return value:
        \n Number of diffusion encoding steps (considering directions, b-values and optionally local averages)
        */
        virtual long   getTotalScans(bool bIncludeLocalAverages = true);
        virtual long   getPreEchoTimeContrib();
        virtual long   getPostEchoTimeContrib();
        virtual MrProtocolData::SeqExpoRFInfo getRFInfoPerRequestMB() const;

        /// Return TE contribution per request [us]
        virtual long   getTEContributionPerRequest();
        /// Return storage time per request [us]
        virtual long   getStorageTimePerRequest();

        virtual void   getActualCounter(long  lDiffLoopCounter, long& lActualAverageCounter, long& lActualBValueCounter, long& lActualDirectionCounter);
        virtual long   getActualAverageCounter(long   lDiffLoopCounter);
        virtual long   getActualBValueCounter(long   lDiffLoopCounter);
        virtual long   getActualDirectionCounter(long   lDiffLoopCounter);
        virtual void   setMaxAmplitude(double dMaxAmplitude);
        virtual void   setMinRiseTime(double dMinRiseTime);
        virtual void   setRFPulseThicknessFactor(double dRFPulseThicknessFactor);
        virtual void   setTEArrayIndex(int    iIndex);
        virtual void   setNoiseThreshold(long   lNoiseThreshold);
        virtual bool   isMagnetizationInverted();
        virtual long   getNoOfDirections();

#ifdef QUIETDWI
        virtual void   setPhaseMomentOffset(double   dMomentum);
        virtual double getPhaseMomentOffset();
        virtual void   setReadMomentOffset(double   dMomentum);
        virtual double getReadMomentOffset();
        virtual void   setSliceMomentOffset(double   dMomentum);
        virtual double getSliceMomentOffset();
        virtual bool   prepMomentOffset(long   lavailableTime);
#endif // QUIETDWI

        /// Get total number of adjustment scans for dynamic distortion correction
        /**
        \b Input:
        \n m_eDynDistMode, m_AdjDidi

        \b Output:
        \n n.a.

        \b Return value:
        \n Number of adjustment scans
        */
        virtual long   getNoOfAdjScans();

        /// Set patient direction and position (required for coordinate transformation of b-matrices)
        virtual void   setPatPosDir(int iPatDirection, int iPatPosition);

        /// Store GPA load of gradient events outside (!) the diffusion module (e.g. EPI readout)
        virtual void setKernelGPALoad(GPABalance &sKernelBalance);

        /// Export GPA load of diffusion module gradient events
        virtual bool getGPALoad(GPABalance &sDiffBalance);

        /// Set / get m_bUseGPABalance
        virtual bool setUseGPABalance(bool bUseGPABalance);
        virtual bool getUseGPABalance();

        /// Set / get m_bMinimizeTE
        virtual void   setMinimizeTE(bool bMinimize);
        virtual bool   getMinimizeTE();

        /// Set / get m_lTRIncrement
        virtual void   setTRIncrement(long lTRIncrement);
        virtual long   getTRIncrement();

        /// Set / get m_iAdjScan
        virtual void setAdjustmentScan(int iAdjScan);
        virtual int  getAdjustmentScan();

        /// Set method for the loop counters
        /** Based on the provided information, the following internal
        loop counters / indices are set:
        m_lRepLoopCounter
        m_lDiffLoopCounter
        m_lBValueCounter
        m_lDirectionCounter
        Return value    : true (ok), false (error)
        */
        virtual bool setLoopCounters
            (
            int       iAdjustmentScan,    /**< Imp:     Adjustment scan index            */
            long      lRepLoopCounter,    /**< Imp:     SeqLoop repetition counter       */
            long      lDiffLoopCounter,   /**< Imp:     SeqLoop diffusion counter        */
            sREADOUT* pADC                /**< Imp/Exp: Mdh                              */
            );

        // Get absolute maximum gradient amplitude
        virtual double getAbsMaxAmpl() const;

        // Set absolute maximum gradient amplitude  
        virtual void setAbsMaxAmpl(double dMaxAmpl);

        // Get method for RESOLVE sequence type
        virtual bool getResolve() const;

        // Set method for RESOLVE sequence type
        virtual void setResolve(bool bValue);

        // Get pointer to DiffusionDirections object
        virtual DiffusionDirections* getDidiPointer();

        // get diffusion gradient duration in ms for display in UI
        virtual double getDiffGradDuration_ms() const;

        // get diffusion gradient spacing in ms for display in UI
        virtual double getDiffGradSpacing_ms() const;

        virtual bool getDiffusionGradientsEnabled() const { return m_bDiffusionGradientsEnabled; }
        virtual void setDiffusionGradientsEnabled(bool val) { m_bDiffusionGradientsEnabled = val; }


        /// set the run mode: either single band or multi band
        virtual void setRunMode(SliceAccelRFRunMode eRunMode);


        /// See SeqBuildBlock for detailed information: Calculate RF info
        bool calcSliceAdjSBBRFInfo(
            MrProt                                     &rMrProt,        //< IMP: points to the protocol structure.
            SeqLim                                     &rSeqLim,        //< IMP: points to the sequence limits structure.
            SeqExpo                                    &rSeqExpo,       //< IMP: points to the sequence exports structure
            const SLICEADJ::sCuboidGeometry              &sSliceAdjCuboid,  //< IMP: cuboid geometry (single)
            std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo        //< EXP: RF info
            ) override = 0;

        // determine smallest b value for IVIM
        virtual long getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch = false) = 0;

        // get IVIM increment
        long getIVIMIncrement() const;

        // get max b value for small IVI increments
        long getMaxBValueSmallIVIMIncrement() const;

        // diffusion loop counter for the highest b-value for kernel check
        long getDiffLoopCounterForHighestBValue();

        // diffusion loop counter for the scan with largest read-axis component
        long getDiffLoopCounterForHighestReadComponent(sSLICE_POS* pSlice);

        // set thermal balancing flag
        void setThermalBalancing(bool bThermalBalancing);

        // get thermal balancing flag
        bool getThermalBalancing() const;

        // set and get the scale factor between the used amplitude and max amplitude for diffusion gradients
        void                setvdScaleFactorinRun(double dScale_p, double dScale_r, double dScale_s);
        std::vector<double> getvdScaleFactorinRun() const;

        // set and get compensation flag
        void setbCompensationEnable(bool bCompensationEnable);
        bool getbCompensationEnable() const;

        // set and get duration between the end of diffusion gradient to the start of compensation gradient
        void setlPlugInToCompGradTime(long lPlugInToCompGradTime);
        long getlPlugInToCompGradTime() const;

        // set up compensation related parameters
        void setCompensationPara(bool bCompensationDecay, double dCompensationFraction, double dEddycurrentTau);

        // set parameters needed before the preparation of compensation gradients
        // e.g. min rise time, max amplitude.
        void prePrepareCompGrad(MrProt& rMrProt, SeqLim& rSeqLim, long lRampTimeOutsideSBB = 0);

        // get the full access to compensation gradient
        SeqBuildBlockCompGrad* getPointerCompGrad();

    protected:

        ///    If true, the diffusion gradient events will be played out
        /** This variable controls if the events for the diffusion gradients will
        be played out during run(). The default value is 'true'. In some
        special cases (such as PAT reference scans) it might be necessary
        to enforce a b=0 scan which can be obtained by setting this value
        temporarily to false during the execution of the sequence run kernel.
        */
        bool m_bDiffusionGradientsEnabled{true};

        /// Copy relevant UI (protocol) parameters to member variables
        /**
        \b Input:
        \n pMrProt, pSeqLim, pSeqExpo

        \b Output:
        \n member variables storing protocol information

        \b Return value:
        \n Success = true, failure = false

        This method is called at the very beginning of ::prep. It copies relevant
        protocol entries to internal member variables. This should help to
        decouple the SBB code from changes in the protocol framework.
        */
        virtual bool setParametersFromUI(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo);

        /// This pure virtual method initializes the preparation
        /**
        \b Input:
        \n Member variables set by setParametersFromUI

        \b Output:
        \n m_dRefSpoilMoment, possibly some prepared RT events (e.g. spoilers - member variables)

        \b Return value:
        \n Success = true, failure = false

        This method is called at the beginning of ::prep. It is used by derived diffusion SBBs to
        set up anything required before the actual timing calculation (::prepTiming) starts.

        Communication with the other parts of the preparation takes place via member variables.
        */
        virtual bool prepInit
            (
            MrProt                  &rMrProt,       /**< Input: protocol               */
            SeqLim                  &rSeqLim,       /**< Input: sequence limits        */
            MrProtocolData::SeqExpo &rSeqExpo       /**< Input: sequence exports       */
            ) = 0;

        /// This pure virtual method calculates the timing
        /**
        \b Input:
        \n lActualTE, member variables set by setParametersFromUI

        \b Output:
        \n m_dMaxPossibleBValue, possibly some prepared RT events (e.g. PE-axis events)

        \b Return value:
        \n Success = true, failure = false

        This method is called within the TE search strategy of ::prep. For the provided protocol,
        it has to set up a valid timing of real time events and to calculate the corresponding
        maximum obtainable b-value.

        Communication with the other parts of the preparation takes place via member variables.
        */
        virtual bool prepTiming
            (
            MrProt                  &rMrProt,       /**< Input: protocol               */
            SeqLim                  &rSeqLim,       /**< Input: sequence limits        */
            MrProtocolData::SeqExpo &rSeqExpo,      /**< Input: sequence exports       */
            long                    lActualTE       /**< Input: acutal TE [ms]         */
            ) = 0;

        /// This pure virtual method finalzes the preparation
        /**
        \b Input:
        \n dMaxRequestedBValue, member variables set by setParametersFromUI
        \n bIsContextPrepForBinarySearch, to indicate preparation context

        \b Output:
        \n m_dAmpl, m_lTRIncrement, m_RFInfoPerRequest, m_lPreEchoTimeContrib,
        m_lPostEchoTimeContrib, m_lSBBStorageTimePerRequest_us, m_bIsMagnetizationInverted, prepared RT events

        \b Return value:
        \n Success = true, failure = false

        This method is called after the timing has been calculated (::prepTiming). All real time
        events have to be prepared afterwards such that the provided maximum requested b-value
        is realized. If necessary, a TR increment can be specified.

        Communication with the other parts of the preparation takes place via member variables.
        */
        virtual bool prepFinal
            (
            double                          dMaxRequestedBValue,     /**< Input: maximum b-value        */
            bool                            bIsContextPrepForBinarySearch = false  /**< Input: is the context for binary search */
            ) = 0;


        ///  The pure virtual method calculates the current b matrix.
        /**  This method must be implemented by each derived diffusion mode and is
        to calculate the b matrix components which are obtained with the actually
        prepared gradients.

        The calculation is usually done by using the BMatrix class. The results
        are to be written to the member variables SBBDiffusion_Base::m_dBxx, m_dByy,
        m_dBzz, m_dBxy, m_dBxz, and m_dByz.
        */
        virtual void calcBMatrix() = 0;


        virtual bool prepParameters(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo);
        virtual bool prepSpoilGrad(double dMoment);
        virtual bool prepGradMoment(sGRAD_PULSE_TRAP *sGrad, double dMoment);
        virtual long getlSpoilerTotalTime() const;

        /// Prepare Dicom header information
        virtual NLS_STATUS CalculateDicomHeaderInformation
            (
            const Slice         &slice,         /**< Imp: reference to current slice           */
            sSLICE_POS*          pSLC,          /**< Imp: pointer to current slice object            */
            DiffusionDirections* pDidi,         /**< Imp: pointer to diffusion directions            */
            double               dAmpl,         /**< Imp: Diffusion gradient amplitude factor [mT/m] */
            long                 lPolarity,     /**< Imp: polarity of diffusion vector               */
            long                 lBValueCounter,    /**< Imp: b-value index                              */
            long                 lDirectionCounter, /**< Imp: diffusion direction index                  */
            sREADOUT*            pADC           /**< Exp: pointer to current readout structure       */
            );

        /// Round up TE in accordance with MrUILink
        virtual long lCalcTEOnInc
            (
            long    lTENotOnInc                 /**< Imp: TE              */
            );

        /// Get transformation from PRS(=GCS) or XYZ(=DCS) coordinate system to PCS
        virtual NLS_STATUS CalculateGradientsToPCS
            (
            sSLICE_POS*          pSLC,
            sROT_MATRIX&         rotMatrix
            );

        /// matC = matA*matB
        void MatMult(sROT_MATRIX const& matA, sROT_MATRIX const& matB, sROT_MATRIX& matC) const;

        /// matC = matA*matB'
        void MatMultTrans(sROT_MATRIX const& matA, sROT_MATRIX const& matB, sROT_MATRIX& matC) const;

        void DidiXYZ2PRS(DiffusionDirections* pDidi, sSLICE_POS* pSLC, long lDirectionCounter,
            double *dDidiP, double *dDidiR, double *dDidiS) const;

        void DidiPRS2XYZ(DiffusionDirections* pDidi, sSLICE_POS* pSLC, long lDirectionCounter,
                                 double *dDidiX, double *dDidiY, double *dDidiZ) const;

        /// Pure virtual prepare of diffusion module GPA load
        /**
        \b Input:
        \n dAmplitude

        \b Output:
        \n m_sBalanceDiff (containing gradient events)

        \b Return value:
        \n Success = true, failure = false

        m_sBalanceDiff has to be filled with all relevant gradient events of the
        diffusion encoding module. Diffusion encoding gradients are set to the
        provided amplitude. Used by the method ::CalcMaximumAmplitude.

        Two variants have to be supported:
        Single axis:  Only the gradient events on the GPABALANCE_X_AXIS have
        to be provided.
        Multi axes :  Gradient events on all axes have to be provided

        */
        virtual bool   prepGPALoadDiff(double dAmplitudeX) = 0;

        virtual bool   prepGPALoadDiff(double dAmplitudeX, double dAmplitudeY, double dAmplitudeZ) = 0;

        /// Pure virtual scaling of diffusion encoding gradient events
        /**
        \b Input:
        \n m_sBalanceDiff, dScale

        \b Output:
        \n m_sBalanceDiff (scaled gradient events)

        \b Return value:
        \n Success = true, failure = false

        All diffusion encoding gradient events within m_sBalanceDiff are scaled by
        the given factor (the same holds for other gradients that are supposed to have
        a fixed amplitude relation to the diffusion encoding gradients). Used by the
        method ::CalcMaximumAmplitude.

        Two variants have to be supported:
        Single axis:  Only the gradient events on the GPABALANCE_X_AXIS have
        to be provided.
        Multi axes :  Gradient events on all axes have to be provided

        */
        virtual bool   scaleGPALoadDiff(double dScaleX) = 0;

        virtual bool   scaleGPALoadDiff(double dScaleX, double dScaleY, double dScaleZ) = 0;

        /// Calculate maximum diffusion encoding amplitudes compatible with the current timing.
        /**
        \b Input:
        \n dUpperLimit [mT/m], m_sBalanceAcq

        \b Output:
        \n dMaxAmpl [mT/m], lTRIncrement [us]

        \b Return value:
        \n Success = true, failure = false

        Given the GPA load of the sequence and the timing of the diffusion gradients,
        the maximum amplitude for the latter is calculated based on an elaborate GPA
        supervision model. Depending on the capabilities of the GPA and the actual
        timing of the gradient events, an additional TR increment is also calculated
        which needs to be applied after each pair of diffusion - readout events.

        By limiting the actual upper limit of the gradient amplitude (input parameter),
        it is possible to change the compromise between minimum TE and minimum TR. In
        general, a reduced upper limit will yield a longer TE and a smaller TR increment.
        If the upper limit is negative, only the absolute maximum diffusion gradient
        amplitude is calculated - no TR increment will be calculated in this mode.

        The actual diffusion SBB implementation has to provide the methods
        prepGPADiff and scaleGPADiff.

        The sequence has to provide an appropriately prepared instance of GPABalance
        (by SBBDiffusion->setKernelGPALoad). It is essential that a 'worst case' gradient
        scheme on the x-axis is provided (diffusion gradients will be considered on
        the same axis).

        CalcMaximumAmplitude and CalcTRIncrement are the wrappers used
        internally to calculate only the maximum amplitude or the required
        TR increment, respectively.
        */
        virtual bool CalcMaximumAmplitude
            (
            MrProt &rMrProt,
            double &dMaxAmplitude         /**< Output: maximum possible amplitude [mT/m] */
            );

        virtual bool CalcTRIncrement
            (
            double dAmplitude,            /**< Input:  actual amplitude [mT/m]    */
            long  &lTRIncrement           /**< Output: required TR increment [us] */
            );

        virtual bool CalcMaxAmpAndTRInc
            (
            double &dMaxAmplitude,        /**< Output: maximum possible amplitude [mT/m]  */
            long   &lTRIncrement,         /**< Output: required TR increment      [us]    */
            double  dUpperLimit           /**< Input:  maximum amplitude to use   [mT/m]  */
            );

        /// Helper function that returns the sign of x (+1 or -1)
        int sign(double x);

        ///    This is the number of diffusion weightings = number of b values.
        /**   This is the same as rMrProt.diffusion().getlDiffWeightings().
        */
        int m_NoOfWeightings{1};

        ///    This is the number of diffusion directions.
        /**   The value is maintained by the SBB and should be the same as
        rMrProt.diffusion().getlDiffDirections().
        */
        long m_lDirections{1};


        ///	This is the time in microseconds used for the diffusion gradient ramps. 
        /**   All diffusion gradients used by this SBB have the same ramp time.
        Simply calculated by the product
        m_lRampTime = fSDSRoundUpGRT (m_dMinRiseTime * m_dMaxAmpl)
        */
        long m_lRampTime{1000};


        ///	This is the time in microseconds used for the selection gradient ramps. 
        /**   All slice selection gradients used by this SBB have the same ramp time.
        Simply calculated by the product
        m_lRampTimeRF = fSDSRoundUpGRT (m_dMinRiseTime * SysProperties::getGradMaxAmplFast())
        // is just used by SBBDiffusion_trace
        */
        long m_lRampTimeRF{1000};

        ///    This is the rise time in [us/(mT/m)] for all gradients used in this SBB
        /**   All gradients used by this SBB have the rise time. The value can be set
        by the method setMinRiseTime).
        */
        double m_dMinRiseTime{0.0};


        ///    This is the absolute maximum allowed gradient amplitude in [mT/m].
        /**
        It applies for all gradients used in this SBB.
        The value has to obey GPA restrictions.
        */
        double m_dMaxAmpl{0.0};


        ///    Actual maximum amplitude of the diffusion gradients for the current timing.
        /**
        If there is no GPA model available, this amplitude is identical to m_dMaxAmpl.
        If a GPA model is available, this amplitude depends on the current timing.
        The value is set during prep().
        */
        double m_dAmpl{0.0};


        ///   This defines a reference spoil moment [ms mT/m]
        /**   Most diffusion encoding schemes require additional spoiling (at lest
        for small b-values) in order to prevent interferences with undesired
        coherence pathways. A reference spoil moment (usually depending on
        the dimension of the acquired k-space region) can be stored herein.
        */
        double m_dRefSpoilMoment{0};


        ///    Duration from the beginning of the SBB till the first pi pulse in us.
        /**   The diffusion SBB contains the pi pulse of the resulting spin echo sequence.
        This value can be retrieved by the method getPreEchoTimeContrib().
        */
        long m_lPreEchoTimeContrib{0};


        ///    This member contains the duration from the last pi pulse to the end of the SBB in us.
        /**   The diffusion SBB contains the pi pulse of the resulting spin echo sequence.
        This value can be retrieved by the method getPostEchoTimeContrib().
        */
        long m_lPostEchoTimeContrib{0};


        ///    Time in microseconds from the start of the ADC till the occurence of the echo.
        /**   It is supposed that the ADC starts at the the end of the Diffusion SBB.
        This value must be specified by calling the method setADCusTillEcho()
        before calling prep() .
        */
        long m_lADCusTillEcho{0};


        ///    Time from the center of the RF excitation pulse and the start of the diffusion SBB.
        /**   Usually this are the durations of the half RF pulse plus rephasing gradient
        plus phase correction scans etc.
        This value must be specified in microseconds by calling the method setSpinPrepTimeus()
        before calling prep() .
        */
        long m_lSpinPrepTimeus{0};

        ///    Pause for decoupling the stimulation potential of diffusion and readout gradients
        /**   A pause with this duration should be inserted after the last diffusion encoding gradient
        for each diffusion module. Even if the diffusion gradients themselves do not exceed
        stimulation limits, adding the potential of the first succeeding EPI readout gradient
        might do so.
        */
        long m_lStimoDelayus{0};

        // compensation gradients used to compensated the eddy current field caused by monopolar diffusion gradients
        SeqBuildBlockCompGrad m_CompGrad{nullptr};

        ///    Actual TE required to realize the requested b-value
        /**
        This value is calculated within ::prep(). If the requested b-value cannot be
        realized with the protocol TE, the minimum required TE is calculated and stored
        herein.
        */
        long m_lActualTE{0};


        ///    Array index of the TE time to use
        /**   Some more sophisticated sequences may have several TE times.
        This variable is used to select the TE array element which is relevant
        for timing calculation. By default, the timing calculation of this SBB
        is based on the first TE time of the array, i.e. MrProt.TE()[0].
        */
        int m_iTEArrayIndex{0};


        ///    Energy applied during one execution of the run-function [Ws].
        MrProtocolData::SeqExpoRFInfo m_RFInfoPerRequestMB;

        ///	Time during which magnetization is stored along the longitudinal axis within the SBB-Run function [us].
        /**   \note Since no T2 relaxation takes place during this time, it will
        be considered as not contribution to TE (the EPI Kernel takes care
        of this).
        */
        long m_lSBBStorageTimePerRequest_us{0};

        // SBBs for refocusing pulse which supports SMS
        SBBMultibandRF m_SBB_RF_Refoc1{nullptr};
        SBBMultibandRF m_SBB_RF_Refoc2{nullptr};

        // prepares all properties of the slice accelerated refocusing SBB
        virtual bool prep_SBBMultibandRF(MrProt &rMrProt, SeqLim  &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo, SBBMultibandRF & SBBMultiband, IRF_PULSE * pBaseRF);


        ///   Correction factor for slice thickness excited by RF inversion pulse
        /**   The slice thickness for which the RF inversion pulse is calculated
        will be multiplied by this factor. Default is 1.0.
        */
        double m_dRFPulseThicknessFactor{1.0};


        ///    Noise threshold used by the ICE program for calculating the ADC maps
        /**   The value is usually set by the method SBBDiffusion_Base::setNoiseThreshold.
        This value is only exported, but not used directly in this SBB.
        */
        long m_lNoiseThreshold{0};


        ///   This flag is true, if the magnetization is inverted by the SBB Diffusion.
        /**   It has to be set by the implementation of the individual diffusion modes
        and can be read out by SBBDiffusion_Base::isMagnetizationInverted().
        Is has to be set before calling prep() (because this calls calcBValue()).
        */
        bool m_bIsMagnetizationInverted{false};

        ///  Spoiler gradient in slice direction.
        /**  The gradient will be prepared by SBBDiffusion_Base::prepSpoilGrad.
        Usually, it will  be called several times in your run function.
        */
        sGRAD_PULSE_TRAP m_DSs1{"DSs1"};

        ///  Spoiler gradient in slice direction.
        /**  The gradient will be prepared by SBBDiffusion_Base::prepSpoilGrad.
        Usually, it will  be called several times in your run function.
        */
        sGRAD_PULSE_TRAP m_DSp1{"DSp1"};

        ///  Spoiler gradient in slice direction.
        /**  The gradient will be prepared by SBBDiffusion_Base::prepSpoilGrad.
        Usually, it will  be called several times in your run function.
        */
        sGRAD_PULSE_TRAP m_DSr1{"DSr1"};

        /// Sequence gradient load
        /** This member is required for the GPA model based diffusion gradient
        amplitude calculation. It has to be provided by the sequence and
        contains all gradient pulses that are applied on the x-axis in a
        worst case situation.

        Example: for a diffusion epi sequence, all readout gradients are
        included (with alternating polarity).
        */
        GPABalance m_sBalanceAcq;

        /// Diffusion encoding gradient load
        /** This member is required for the GPA model based diffusion gradient
        amplitude calculation. It has to be provided by the actual diffusion
        encoding SBB.
        */
        GPABalance m_sBalanceDiff;

        /// Internal switch to enable / disable GPA load model
        /** Will be set to 'true' if setGPALoad is called. Default value is
        'false' which means that the diffusion module uses a predefined
        constant diffusion gradient amplitude (similar to older versions).
        A dynamic switching between 'GPA load' and 'constant amplitude'
        is not designed!
        */
        bool m_bUseGPABalance{false};

        /// Storage of required TR increment
        /** Depending on the GPA load, a TR increase might be required in
        order to guarantee that the protocol can run 'forever'. This
        value is stored by the diffusion SBB and must be considered by
        the sequence.
        */
        long m_lTRIncrement{0};

        /// Actual repetition loop counter from SeqLoop
        /** Set within setLoopCounters.
        This is the actual value of the SeqLoop repetitions counter. It
        is required to set the correct MDH entries.
        */
        long m_lRepLoopCounter{-1};

        /// Actual diffusion loop counter from SeqLoop
        /** Set within setLoopCounters.
        This is the internal diffusion scan counter which runs over the b values
        and directions. For further documentation refer to the methods
        getActualBValueCounter() and getActualDirectionCounter().
        */
        long m_lDiffLoopCounter{-1};

        /// Actual b-value
        /** Stores b-value index. Set within setLoopCounters.
        */
        long m_lBValueCounter{-1};

        /// Actual diffusion direction
        /** Stores diffusion direction index. Set within setLoopCounters
        */
        long m_lDirectionCounter{-1};

        /// Setting for actual adjustment scan
        /** For values > 0, adjustment diffusion scans will be
        played out (with dedicated directions and b-values) instead
        of imaging diffusion scans.
        */
        int m_iAdjScan{-1};

        /// TE minimization mode
        /** If set to true (usually by the sequence), the ::prep method
        does not only check whether the maximum b-value can be realized
        with the protocol TE value, but the minimum possible TE value
        is calculated.
        The sequence (UI handlers) has to take care about actually
        setting this TE value and updating the protocol.
        */
        bool m_bMinimizeTE{false};

        /// The \f$b_{xx}\f$ component of the b matrix
        double m_dBxx{0.0};

        /// The \f$b_{xy}\f$ component of the b matrix
        double m_dBxy{0.0};

        /// The \f$b_{xz}\f$ component of the b matrix
        double m_dBxz{0.0};

        /// The \f$b_{yy}\f$ component of the b matrix
        double m_dByy{0.0};

        /// The \f$b_{yz}\f$ component of the b matrix
        double m_dByz{0.0};

        /// The \f$b_{zz}\f$ component of the b matrix
        double m_dBzz{0.0};

        /// The \f$b_{zz}\f$ component of the b matrix
        double m_dBValue{0.0};

        ///   Maximum possible b value with current timing
        /**    This is the b value which will be obtained by the diffusion gradients with
        the actually calculated timing and the maximum permitted gradient amplitudes.
        */
        double m_dMaxPossibleBValue{0.0};

        /// Management of different directions is provided by the DiffusionDirections class
        DiffusionDirections m_Didi;

        /// Directions for dynamic distortion correction adjustment scans
        DiffusionDirections m_AdjDidi;

        /// Coordinate system of the external vector set (FREE mode)
        MrProtocolData::DiffDirCoordinateSystem m_eFreeCoordinateSystem{MrProtocolData::DIFFDIR_CS_XYZ};

        /// Comment of the external vector set (FREE mode)
        std::string m_strFreeUserComment;

        /// External diffusion directions (FREE mode)
        std::vector<VectorStruct> m_vFreeDiffDir;

        /// Actual patient direction (from MeasPatient)
        int m_iPatDirection{0};

        /// Actual patient position (from MeasPatient)
        int m_iPatPosition{0};

        /// Gyromagnetic ratio of current nucleus [1/T s] - set by setParametersFromUI
        double m_dGamma{0.0};

        /// Number of repetitions - set by setParametersFromUI
        long m_lRepetitions{0};

        // Maximum value in average array, not including global average 
        long m_lMaxValueInAveArray{0};

        /// RF pulse type - set by setParametersFromUI
        SEQ::RFPulseType m_eRFPulseType{SEQ::RF_NORMAL};

        /// Diffusion mode - set by setParametersFromUI
        SEQ::DiffusionMode m_eDiffusionMode{SEQ::DIFFMODE_NONE};

        /// Dynamic distortion correction mode - set by setParametersFromUI
        SEQ::DynamicDistortionCorrMode m_eDynDistMode{SEQ::DYN_DISTCORR_NONE};

        /// Number of diffusion directions used in MDDW mode - set by setParametersFromUI
        long m_lDiffusionDirectionsMDDW{0};

        /// Q-Space coverage
        SEQ::DiffQSpaceCoverageMode m_eQSpaceCoverage{SEQ::DIFF_QSPACE_COVERAGE_FULL};

        /// Q-Space sampling scheme
        SEQ::DiffQSpaceSamplingScheme m_eQSpaceSampling{SEQ::DIFF_QSPACE_SAMPLING_CARTESIAN};

        /// Q-Space maximum b-value
        double m_dQSpaceMaxBValue{0.0};

        /// Q-Space steps
        long m_lQSpaceSteps{0};

        /// b-value increment (from SeqLim) - set by setParametersFromUI
        long m_lBValueInc_Limit{0};

        /// TE increment (from SeqLim) - set by setParametersFromUI
        long m_lTEInc_Limit{0};

        /// Number of slices - set by setParametersFromUI
        long m_lSlices{0};

        /// Slice thickness - set by setParametersFromUI
        double m_dSliceThickness{0.0};

        /// Readout moment - set by setParametersFromUI
        double m_dReadoutMoment{0.0};

#ifdef QUIETDWI
        /// offset Moment for diffusion-encoding gradient in Phase direction 
        double m_dPhaseMomentOffset{0.0};

        /// offset Moment for diffusion-encoding gradient in Slice direction
        double m_dSliceMomentOffset{0.0};

        /// offset Moment for diffusion-encoding gradient in Read direction 
        double m_dReadMomentOffset{0.0};

        /// offset diffusion-encoding gradient in phase direction
        sGRAD_PULSE_TRAP m_DGoffp{"RTEIdentDGoffp"};

        /// offset diffusion-encoding gradient in slice direction
        sGRAD_PULSE_TRAP m_DGoffs{"RTEIdentDGoffs"};

        /// offset diffusion-encoding gradient in read direction
        sGRAD_PULSE_TRAP m_DGoffr{"RTEIdentDGoffr"};

#endif // QUIETDWI

        /// Vector containing diffusion weightings (b-values) - set by setParametersFromUI
        std::vector<double> m_vdBValues{0.0};

        /// Vector containing average number for different B-values  
        std::vector<long> m_vlLocalAverages;


        /// Vector containing maximum TE (from SeqLim) - set by setParametersFromUI
        std::vector<long> m_vlTEMax_Limit{0};

        /// Vector containing actual TE - set by setParametersFromUI
        std::vector<long> m_vlTE{0};

        ///   Flag to specify RESOLVE sequence type
        /**	This parameter is used to activate specific behavior
        when SBBDiffusion is used by the RESOLVE sequence
        */
        bool m_bResolve{false};

        /// Diffusion gradient duration in ms (read by sequence and passed to UI for display)
        double m_dDiffGradDuration_ms{0};

        /// Diffusion gradient spacing in ms (read by sequence and passed to UI for display)
        double m_dDiffGradSpacing_ms{0};

        /// Zoomed Excitation (ZOOMit) is used in this protocol
        bool m_bZoomedExcitationIsUsed{false};

        /// Factor to reduce maximum amplitude of diffusion gradients.
        double m_dDiffGradMaxAmpReduction{0.9};

        SliceAccelRFRunMode m_eSliceAccelRFRunMode{SINGLE_BAND};

        sREADOUT* m_pADC{nullptr};

        // RF Pulse library which holds all all information about the used RF pulses
        DiffusionRFPulseProperties m_RFPulseLibrary;

        // debug settings
        SequenceDebugSettings::SequenceDebugSettings m_debugSettings = SequenceDebugSettings::SequenceDebugSettings("USE_EPI_DEBUG_SETTINGS");

        // used in the rounding function for TE
        const double m_MaxValueForRounding{3600000000.0};

    private:

        // container for diffusion ordering information
        DiffusionOrdering m_sDiffusionOrderInfo;

        // prepare diffusion ordering
        bool prepDiffusionOrder(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo);

        // flag for thermal balancing to be set by feature toggle
        bool m_bThermalBalancing{false};

        ///   This is the total time of a spoiler gradient in us.
        /**   This variable is usually set by SBBDiffusion_Base::prepSpoilGrad.

        (In general, a second spoiler is required for refocusing. This time
        specifies only the time for one.)
        */
        long m_lSpoilerTotalTime{0};

        // duration between the end of diffusion gradient to the start of compensation gradient
        long m_lPlugInToCompGradTime{0};

        // compensation flag
        bool m_bCompensationEnable{false};

        // compensation related parameters
        bool   m_bCompensationDecay{false};
        double m_dCompensationFraction{1.0};
        double m_dEddycurrentTau{0.0};

        // the scale factor between the used amplitude and max amplitude for diffusion gradients
        std::vector<double> m_vdScaleFactor_prs{1.0, 1.0, 1.0};
    };




    // -------------------------------------------------------------
    // Inline function declarations for Class SBBDiffusion_Base
    // -------------------------------------------------------------

#ifdef QUIETDWI
    /// set offset Moment Phase
    inline void SBBDiffusion_Base::setPhaseMomentOffset(double dMomentum)
    {
        m_dPhaseMomentOffset = -dMomentum;
    }
    /// Get offset Moment Phase
    inline double SBBDiffusion_Base::getPhaseMomentOffset(void)
    {
        return m_dPhaseMomentOffset;
    }

    /// set offset Moment Read
    inline void SBBDiffusion_Base::setReadMomentOffset(double dMomentum)
    {
        m_dReadMomentOffset = -dMomentum;
    }
    /// Get offset Momentum Read
    inline double SBBDiffusion_Base::getReadMomentOffset(void)
    {
        return m_dReadMomentOffset;

    }

    /// set offset Momentum Slice
    inline void SBBDiffusion_Base::setSliceMomentOffset(double dMomentum)
    {
        m_dSliceMomentOffset = -dMomentum;
    }
    /// Get offset Momentum Slice
    inline double SBBDiffusion_Base::getSliceMomentOffset(void)
    {
        return m_dSliceMomentOffset;
    }

    //prepare offset moment
    inline bool SBBDiffusion_Base::prepMomentOffset(long lAvailableTime)
    {
        if(!m_DGoffp.prepSymmetricTOTExactMomentum(getPhaseMomentOffset(), (double)lAvailableTime)){
            setNLSStatus(m_DGoffp.getNLSStatus());
            return false;
        }
        if(!m_DGoffs.prepSymmetricTOTExactMomentum(getSliceMomentOffset(), (double)lAvailableTime)){
            setNLSStatus(m_DGoffs.getNLSStatus());
            return false;
        }
        if(!m_DGoffr.prepSymmetricTOTExactMomentum(getReadMomentOffset(), (double)lAvailableTime)){
            setNLSStatus(m_DGoffr.getNLSStatus());
            return false;
        }

        if(!m_DGoffp.check()){
            setNLSStatus(m_DGoffp.getNLSStatus());
            return false;
        }
        if(!m_DGoffs.check()){
            setNLSStatus(m_DGoffs.getNLSStatus());
            return false;
        }
        if(!m_DGoffr.check()){
            setNLSStatus(m_DGoffr.getNLSStatus());
            return false;
        }
        return true;
    }
#endif // QUIETDWI


    ///    Time in microseconds from the start of the ADC till the occurence of the echo.
    /**   This inline method just writes the specified time to  m_lADCusTillEcho.
    The time of the end of the Diffusion SBB is supposed to be identical
    to the start time of the ADC.
    The sequence programmer must specify this time before calling prep().
    */
    inline void SBBDiffusion_Base::setADCusTillEcho(long lADCusTillEcho)
    {
        if(lADCusTillEcho != m_lADCusTillEcho)
        {
            m_lADCusTillEcho = lADCusTillEcho;
            resetPrepared();
        }
    }

    ///   Time in microseconds from the center of the RF excitation pulse and the start of the diffusion  module.
    /**    Usually this is the duration of the half RF pulse plus rephasing gradient
    plus phase correction scans etc.

    This inline method just writes the specified  time to m_lSpinPrepTimeus.
    The sequence programmer must specify this time before calling prep().
    */
    inline void SBBDiffusion_Base::setSpinPrepTimeus(long lSpinPrepTimeus)
    {
        if(lSpinPrepTimeus != m_lSpinPrepTimeus)
        {
            m_lSpinPrepTimeus = lSpinPrepTimeus;
            resetPrepared();
        }
    }

    inline void SBBDiffusion_Base::setlPlugInToCompGradTime(long lPlugInToCompGradtime)
    {
        if (lPlugInToCompGradtime != m_lPlugInToCompGradTime)
        {
            m_lPlugInToCompGradTime = lPlugInToCompGradtime;
            resetPrepared();
        }
    }

    inline long SBBDiffusion_Base::getlPlugInToCompGradTime() const
    {
        return m_lPlugInToCompGradTime;
    }

    inline std::vector<double> SBBDiffusion_Base::getvdScaleFactorinRun() const
    {
        return m_vdScaleFactor_prs;
    }

    inline void SBBDiffusion_Base::setvdScaleFactorinRun(double dScale_p, double dScale_r, double dScale_s)
    {
        m_vdScaleFactor_prs = {dScale_p, dScale_r, dScale_s};
    }

    inline bool SBBDiffusion_Base::getbCompensationEnable() const
    {
        return m_bCompensationEnable;
    }

    inline void SBBDiffusion_Base::setbCompensationEnable(bool bCompensationEnable)
    {
        m_bCompensationEnable = bCompensationEnable;
        m_CompGrad.resetPrepared();
    }

    inline void SBBDiffusion_Base::setCompensationPara(
        bool bCompensationDecay, double dCompensationFraction, double dEddycurrentTau)
    {
        m_bCompensationDecay    = bCompensationDecay;
        m_dCompensationFraction = dCompensationFraction;
        m_dEddycurrentTau       = dEddycurrentTau;
    }

    inline SeqBuildBlockCompGrad* SBBDiffusion_Base::getPointerCompGrad()
    {
        return &m_CompGrad;
    }

    /// Store kernel GPA balance
    /**
    \b Input:
    sKernelBalance

    \b Output:
    n.a.

    */
    inline void SBBDiffusion_Base::setKernelGPALoad(GPABalance &sKernelBalance)
    {
        m_sBalanceAcq = sKernelBalance;
        resetPrepared();
    }

    /// Export GPA load of diffusion module gradient events
    inline bool SBBDiffusion_Base::getGPALoad(GPABalance &sDiffBalance)
    {
        if((getUseGPABalance() == false) || (!isPrepared()))
        {
            return false;
        }

        sDiffBalance = m_sBalanceDiff;
        return true;
    }


    /// Enable / disable GPA balance calculations
    /**
    \b Input:
    bUseGPABalance

    \b Output:
    true if new value has been successfully applied, false otherwise

    */
    inline bool SBBDiffusion_Base::setUseGPABalance(bool bUseGPABalance)
    {
        // Check whether acutal system supports balance models
        if(bUseGPABalance && (m_sBalanceDiff.lGetStatus() != 0))
        {
            m_bUseGPABalance = false;

            return false;
        }

        if(bUseGPABalance != m_bUseGPABalance)
        {
            m_bUseGPABalance = bUseGPABalance;
            resetPrepared();
        }

        return true;
    }

    /// Get setting for GPA balance calculations
    /**
    \b Return value:
    \n m_bUseGPABalance

    */
    inline bool SBBDiffusion_Base::getUseGPABalance()
    {
        return m_bUseGPABalance;
    }

    /// Store required TR increment (cooling pause) for TRFill calculations
    /**
    \b Input:
    \n lTRIncrement

    \b Output:
    \n n.a.

    */
    inline void SBBDiffusion_Base::setTRIncrement(long lTRIncrement)
    {
        m_lTRIncrement = lTRIncrement;

        return;
    }

    /// Receive required TR increment (cooling pause) for TRFill calculations
    /**
    \b Return value:
    \n m_lTRIncrement

    */
    inline long SBBDiffusion_Base::getTRIncrement()
    {
        return m_lTRIncrement;
    }

    /// Set TE minimization mode
    /**
    \b Input:
    \n bMinimize

    \b Output:
    \n n.a.

    */
    inline void SBBDiffusion_Base::setMinimizeTE(bool bMinimize)
    {
        if(bMinimize != m_bMinimizeTE)
        {
            m_bMinimizeTE = bMinimize;
            resetPrepared();
        }
    }

    /// Get TE minimization mode
    /**
    \b Return value:
    \n m_bMinimizeTE

    */
    inline bool SBBDiffusion_Base::getMinimizeTE()
    {
        return m_bMinimizeTE;
    }

    /// Set adjustment scan
    /**
    Default value indicates imaging scans.

    \b Input:
    \n iAdjScan

    \b Output:
    \n n.a.

    */
    inline void SBBDiffusion_Base::setAdjustmentScan(int iAdjScan = 0)
    {
        m_iAdjScan = iAdjScan;

        return;
    }

    /// Get adjustment scan
    /**
    \b Return value:
    \n m_iAdjScan

    */
    inline int SBBDiffusion_Base::getAdjustmentScan()
    {
        return m_iAdjScan;
    }




    ///    Returns the duration from the beginning of the SBB till the first pi pulse.
    /**   This inline method returns the value of the member m_lPreEchoTimeContrib.
    */
    inline long SBBDiffusion_Base::getPreEchoTimeContrib()
    {
        return m_lPreEchoTimeContrib;
    }

    ///    Return the time between the last RF pulse and the end of the SBB
    /**   The diffusion SBB contains the pi pulse of the resulting spin echo sequence.
    This methods returns the duration from the last pi pulse till the end of the
    SBB which is stored in the member variable m_lPostEchoTimeContrib.
    */
    inline long SBBDiffusion_Base::getPostEchoTimeContrib()
    {
        return m_lPostEchoTimeContrib;
    }

    ///   This method returns true, if the magnetization is inverted by the SBB Diffusion.
    /**   It returns just the value of the member SBBDiffusion_Base::m_bIsMagnetizationInverted.
    This flag is true for diffusion modes which apply an odd number of RF refocussing pulse.
    */
    inline bool SBBDiffusion_Base::isMagnetizationInverted()
    {
        return m_bIsMagnetizationInverted;
    }

    ///   Get number of directions a diffusion weighted measurement is to be performed in.
    /**   This value is the same as rMrProt.diffusion().getlDiffDirections().
    This information is required to configure the SeqUT.
    \attention It is obsolete and will be removed as soon as VA15A support has gone.
    */
    inline long SBBDiffusion_Base::getNoOfDirections()
    {
        return m_lDirections;
    }

    ///   Get number of adjustment scans for the dynamic distortion correction
    inline long SBBDiffusion_Base::getNoOfAdjScans()
    {

        // From DiffusionECCUtils.h: Maximum b-value for DFC reference images
        const double maxRefBValue     = 75.0; // [s/mm2]
        const auto  firstImageBValue = (m_vdBValues.size() > 0) ? fabs(m_vdBValues[0]) : 100.f;
        const bool  isReferenceImage = ((firstImageBValue - maxRefBValue) < std::numeric_limits<double>::epsilon());

        switch(m_eDynDistMode)
        {
            case SEQ::DYN_DISTCORR_ADJ:
                // Full set of adjustment scans with dedicated diffusion directions and weightings
                return ADJ_PREP_SCAN_AVERAGES * m_AdjDidi.getNumberOfDirections();
            case SEQ::DYN_DISTCORR_DIRECT:
                // Acquire reference images without diffusion encoding
                if ((m_eDiffusionMode == SEQ::DIFFMODE_QSPACE) || (m_eDiffusionMode == SEQ::DIFFMODE_FREE)
                    || !isReferenceImage)
                {
                    return ADJ_PREP_SCAN_AVERAGES;
                }
                else
                {
                    // If the first imaging volume gets acquired with a b-value which is used
                    // as an (undistorted) reference by the DFC algorithm, there is no need
                    // to acquire extra reference images.
                    return 0;
                }
                break;
            default:
            case SEQ::DYN_DISTCORR_NONE:
                // No adjustment scans
                return 0;
        }
    }

    ///   This method sets the thickness factor for the RF pulse
    /**   This is a simple set function for the member SBBDiffusion_Base::m_dRFPulseThicknessFactor.
    */
    inline void SBBDiffusion_Base::setRFPulseThicknessFactor(double dRFPulseThicknessFactor)
    {
        if(dRFPulseThicknessFactor != m_dRFPulseThicknessFactor)
        {
            m_dRFPulseThicknessFactor = dRFPulseThicknessFactor;
            resetPrepared();
        }
    }


    ///   Returns the energy of one execution of the SBB-run function with a multi-band pulse.
    /**   The energy has default value if m_bPrepared flag has not been set.
    */
    inline MrProtocolData::SeqExpoRFInfo SBBDiffusion_Base::getRFInfoPerRequestMB() const
    {
        if(isPrepared())
        {
            return m_RFInfoPerRequestMB;
        }
        else
        {
            return MrProtocolData::SeqExpoRFInfo();
        }
    }


    ///	Returns the TE contribution [us] of one execution of the SBB-run function
    inline long SBBDiffusion_Base::getTEContributionPerRequest()
    {
        if(isPrepared())
        {
            return getDurationPerRequest();
        }
        else
        {
            return 0;
        }
    }

    ///	Returns the storage (mixing) time [us] of one execution of the SBB-run function
    /**     This inline method returns the value of the member m_lSBBStorageTimePerRequest_us.

    The duration is 0 if m_bPrepared flag has not been set.
    */
    inline long SBBDiffusion_Base::getStorageTimePerRequest()
    {
        if(isPrepared())
        {
            return m_lSBBStorageTimePerRequest_us;
        }
        else
        {
            return 0;
        }
    }

    ///    Return the total time for spoiler gradient in us.
    inline long SBBDiffusion_Base::getlSpoilerTotalTime() const
    {
        return m_lSpoilerTotalTime;
    }

    inline void SBBDiffusion_Base::setPatPosDir(int iPatDirection, int iPatPosition)
    {
        m_iPatDirection = iPatDirection;
        m_iPatPosition  = iPatPosition;
    }

    // Get absolute maximum gradient amplitude
    inline double SBBDiffusion_Base::getAbsMaxAmpl() const { return m_dMaxAmpl; }

    ///	Set absolute maximum gradient amplitude
    inline void SBBDiffusion_Base::setAbsMaxAmpl(double dMaxAmpl)
    {
        m_dMaxAmpl = dMaxAmpl;
        resetPrepared();
    }

    ///	Get value of flag that specifies RESOLVE sequence type
    inline bool SBBDiffusion_Base::getResolve() const  { return m_bResolve; }

    ///	Set value of flag that specifies RESOLVE sequence type
    inline void SBBDiffusion_Base::setResolve(bool bValue) { m_bResolve = bValue; }

    // Get pointer to DiffusionDirections object
    inline DiffusionDirections* SBBDiffusion_Base::getDidiPointer()
    {
        return &m_Didi;
    }

    // get diffusion gradient duration in ms for display in UI
    inline double SBBDiffusion_Base::getDiffGradDuration_ms() const
    {
        return (m_dDiffGradDuration_ms);
    }

    // get diffusion gradient spacing in ms for display in UI
    inline double SBBDiffusion_Base::getDiffGradSpacing_ms() const
    {
        return (m_dDiffGradSpacing_ms);
    }

    inline void SBBDiffusion_Base::setRunMode(SliceAccelRFRunMode eRunMode)
    {
        m_eSliceAccelRFRunMode = eRunMode;
    }

}//end of namespace SEQ_NAMESPACE





