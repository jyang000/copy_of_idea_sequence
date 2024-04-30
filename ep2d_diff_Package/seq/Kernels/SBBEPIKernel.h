//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2021  All Rights Reserved.
//    -----------------------------------------------------------------------------

#pragma once

#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sSYNC.h"
#include "MrMeasSrv/SeqIF/libRT/sSYNCDATA.h"
#include "MrImaging/libSBB/SBBEPIReadOut.h"
#include "MrImaging/libSBB/SBBExcitation.h"
#include "MrImagingFW/libBalance/GPABalance.h"
#include "MrImaging/libSBB/SBBMultibandRF.h"
#include "MrImaging/libSBB/SBBBinomialPulses.h"
#include <algorithm>

#ifdef ZOOM_2DRF
#include "MrImaging/libSeqPTX/SBB2DExc.h"
#include "MrImaging/libSeqPTX/RF2DArbGenerator.h"
#include "MrImaging/libSeqPTX/SBB2DPtxPulsesIni.h"
#include "MrImaging/libSeqPTX/SBB2DExcPulsesIni.h"
#include "MrImaging/libSeqPTX/SBB2DIniPulses_ID_defines.h"
#endif //ZOOM_2DRF

#define SBBEPIKernel_EPIRO_GRAD_PERF_RO                 (SBBEPIReadOut_GRAD_PERF_RO            )
#define SBBEPIKernel_EPIRO_GRAD_PERF_BLIPS              (SBBEPIReadOut_GRAD_PERF_BLIPS         )
#define SBBEPIKernel_GRAD_PERF_GRAD_MOMENTS             (SBBEPIReadOut_NO_OF_GRAD_PERF_GROUPS+0)
#define SBBEPIKernel_NO_OF_GRAD_PERF_GROUPS             (SBBEPIReadOut_NO_OF_GRAD_PERF_GROUPS+1)

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

namespace SEQ_NAMESPACE
{

    //
    //
    //    Author(s):
    //    11/2000-02/2002: Thomas Kluge, @siemens.com, Tel. +49 9131 84 8049
    //
    //    A graphical representation of the EPI-Kernel can be found in
    //    /src/MrImaging/seq/Kernels\SBBEPIKernel.pdf.
    //
    //    This kernel can be used for single-shot or segmented EPI sequences. It is
    //    derived from the class SeqBuildBlockEPIReadOut. For excitation it uses an
    //    excitation SBB which has to be  instantiated and configured within the
    //    sequence. The kernel itself takes care that the gradient moments required
    //    for prephasing and rephasing are realized. It also takes care of TE-fill
    //    time calculations. The kernel can send sync-bits before the excitation SBB.
    //    The kernel supports echo-shifting.
    //
    //    To execute the kernel it must be supplied with a function which defines and
    //    configures the excitation SBB. The sequence has to implement a function of
    //    the following type:
    //    SeqBuildBlockExcitation* configureMyExcitationSBB
    //    (
    //      MrProtocolData::MrProtData*  pMrProt,
    //      SeqLim* pSeqLim,
    //      SeqExpo* pSeqExpo
    //    );
    //    The pointer to this function must be registered at the kernel. Preparation,
    //    execution etc. of the excitation SBB is handled by the EPI kernel itself.
    //
    //    Between excitation and the EPI read-out train there can be sent a real time
    //    event block plug-in (RTEB-plug-in). The implementation of the plug-in
    //    functionality has to be implemented in the sequence. Without RTEB-plug-in
    //    the kernel behaves like a FID-EPI sequence. Using an RTEB-plug-in for
    //    example an SE-EPI sequence can be realized as well as advanced diffusion-EPI
    //    measurements. The net sum of the gradient moments of order zero applied by
    //    the RTEB-plug-in must be zero. The plug-in must be executed in its own
    //    event-block.
    //    The kernel needs two functions to execute the RTEB-plug-in:
    //
    //    NLS_STATUS fPrepMyPlugInForEPIKernel
    //    (
    //        MrProt  *pMrProt,
    //        SeqLim  *pSeqLim,
    //        SeqExpo *pSeqExpo,
    //        long     lEPIKernelTEContributionBeforeRTEBPlugIn,
    //        long     lEPIKernelTEContributionAfterRTEBPlugIn,
    //        long    *plRTEBPlugInDurationPerRequest,
    //        long    *plRTEBPlugInTEContribution,
    //        MrProtocolData::SeqExpoRFInfo  *pRFInfoPerRequest,
    //        bool    *pRTEBPlugInInvertsMagnetization
    //    );
    //    This function has to prepare the plug-in. Also it has to supply the kernel
    //    with the minimum required information: duration of the RTEB plugged in, the
    //    TE contribution of the RTEB plugged in (this is the time from the "center
    //    concerning TE" of the plug-in to its end), the energy applied by the plug-in
    //    and a statement, if the magnetization is inverted by the plug in or not.
    //
    //    NLS_STATUS fRunMyPlugInForEPIKernel
    //    (
    //        MrProt      *pMrProt,
    //        SeqLim      *pSeqLim,
    //        SeqExpo     *pSeqExpo,
    //        sSLICE_POS  *pSLC,
    //        long         lRTEBPlugInTEFillBefore,
    //        long         lRTEBPlugInTEFillAfter
    //    )
    //    This function executes the RTEB-plug-in AND has to apply the TE-fill times.
    //
    //    An alternative prepare method is available that supports balance models:
    //
    //    NLS_STATUS fPrepMyPlugInBalanceForEPIKernel
    //    (
    //        MrProt     *pMrProt,
    //        SeqLim     *pSeqLim,
    //        SeqExpo    *pSeqExpo,
    //        long        lEPIKernelTEContributionBeforeRTEBPlugIn,
    //        long        lEPIKernelTEContributionAfterRTEBPlugIn,
    //        long       *plRTEBPlugInDurationPerRequest,
    //        long       *plRTEBPlugInTEContribution,
    //        MrProtocolData::SeqExpoRFInfo     *pRFInfoPerRequest,
    //        bool       *pRTEBPlugInInvertsMagnetization
    //        GPABalance *psEPIKernelBalanceContribution,
    //        long       *plRTEBPlugInRequiredPausePerRequest
    //    );
    //    If the usage of balance models is enabled (setUseBalanceBalanceModel), this
    //    plugin preparation variant is called. The EPI kernel gradient load is 
    //    provided by a pointer to a corresponding GPABalance instance. The
    //    plugin can use this information internally - it has to export the total 
    //    required pause (considering EPI kernel and plugin gradient events) per
    //    execution. A negative pause value indicates that the kernel should take care
    //    of the corresponding calculation completely ignoring any plugin contribution
    //    (reasonable e.g. for a very short plugin without significant gradient activity).
    //
    //    The EPI kernel supports three different types of phase correction methods,
    //    all of them using three EPI-echos to generate two phase correction data
    //    sets, one for ADCs acquired under a negative and one for  ADCs acquired
    //    under a positive read out gradient.
    //
    //    Internal phase correction scans:
    //    That means the three phase correction echos are located directly after the
    //    excitation SBB and before the RTEB-plug-in. So the phase-correction scans
    //    are acquired with each call of the kernel. So this method in general would
    //    be used only for single-shot EPI sequences.
    //
    //    Extra phase correction scan:
    //    The sequence must put the EPI kernel into the phase-correction mode. Then
    //    the EPI kernel acquires three phase correction echos instead of the imaging
    //    data, i.e. the imaging EPI read-out is replaced by the phase correction
    //    scans and appropriate fill times. The second echo has exactly the TE as the
    //    k-space center line/partition of the imaging echos. This method should be
    //    useful for segmented SE EPI sequences.
    //
    //    Extra early FID phase correction scan:
    //    The sequence must put the EPI kernel into the phase-correction mode. Then
    //    the EPI kernel acquires three phase correction echos directly after the
    //    excitation SBB. The rest of the kernel, except the final gradient moment
    //    section, is replaced by fill time. The phase correction scans therefore have
    //    the minimum possible TE. This method should be useful for segmented FID EPI
    //    sequences.
    //
    //    The user can tell the kernel to preface PE/3D-blips and/or the RO-gradients
    //    after the RTEB-plug-in directly before the imaging EPI read-out.
    //
    //    The user can rotate the frame of reference by specifying an additional
    //    phase. Note that the method setAdditionalPhase also applies the additional
    //    phase to the excitation SBB.
    //
    //      2015-02-02: Mario Zeller
    //      For multi-band imaging, the aforementioned methods take an additional
    //      parameter MrProtocolData::SeqExpoRFInfo for multi-band related energy
    //      calculations as last parameter.

    class __IMP_EXP SeqBuildBlockEPIKernel:
        public SeqBuildBlockEPIReadOut
    {
    public:
        // No public default constructor
        SeqBuildBlockEPIKernel() = delete;

        //--------------------------------------------------------------------
        //    Constructor

        SeqBuildBlockEPIKernel(
            SBBList* pSBBList
            );

        //--------------------------------------------------------------------
        //  Destructor
        virtual ~SeqBuildBlockEPIKernel() = default;

        SeqBuildBlockEPIKernel(const SeqBuildBlockEPIKernel& right) = delete;
        SeqBuildBlockEPIKernel& operator=(const SeqBuildBlockEPIKernel& right) = delete;
        SeqBuildBlockEPIKernel(SeqBuildBlockEPIKernel&& right)                 = delete;
        SeqBuildBlockEPIKernel& operator=(SeqBuildBlockEPIKernel&& right) = delete;

        /// Overloaded base class method: Calculate a single SBB-specific cuboid
        bool calcSingleSliceAdjSBBCuboid(
            MrProt                  &rMrProt,       //< IMP: points to the protocol structure.
            SeqLim                  &rSeqLim,       //< points to the sequence limits structure.
            SeqExpo                 &rSeqExpo,      //< points to the sequence exports structure.
            const sSLICE_POS*              pSLC,          //< IMP: points to the slice position structure of the current slice
            SLICEADJ::sCuboidGeometry &sSliceAdjCuboid  //< EXP: cuboid geometry (single)
            ) override;

        /// Overloaded base class method: Calculate all SBB-specific cuboids
        bool calcSliceAdjSBBCuboids(
            MrProt                               &rMrProt,         //< IMP: points to the protocol structure.
            SeqLim                               &rSeqLim,       //< points to the sequence limits structure.
            SeqExpo                              &rSeqExpo,      //< points to the sequence exports structure.
            std::vector<SLICEADJ::sCuboidGeometry> &vsSliceAdjCuboids  //< EXP: cuboid geometries (multiple elements)
            ) override;

        // The kernel does not store any energy at all but always asks the incorporated SBBs for their energy
        // and passes this information ==> we can overload the method and do nothing here
        bool updateSliceAdjLocalRFInfo(
            MrProt                                     &rMrProt,        //< IMP: points to the protocol structure.
            SeqLim                                     &rSeqLim,        //< IMP: points to the sequence limits structure.
            SeqExpo                                    &rSeqExpo        //< IMP: points to the sequence exports structure
            ) override;

        //--------------------------------------------------------------------
        //    Tells the kernel whether sync-bits should be sent
        //  before the excitation SBB or not. If so, an osc bit
        //  and an external trigger bit are prepared.
        //    Member m_lMaxSyncBitDuration is set.
        //    If the function executed with success, true is
        //  returned. The execution of the osc-bit can be disabled/enabled
        //  during run-time of the sequence using the method setDoNotSendOscBit.
        //    The execution of the external trigger bit can be
        //  disabled/enabled during run-time of the sequence using the method
        virtual bool setUseSyncBits(bool bValue, long lOscChannel = 0,     // channel for the osc bit
                                    long lOscDuration = 10,         // duration of osc bit
                                    long lOscStartTime = 0,         // start time of osc bit in event block
                                    long lExtTrigDuration = 10,     // duration of external trigger bit
                                    long lExtTrigStartTime = 0    // start time of external trigger bit in event block
                                    );

        //--------------------------------------------------------------------
        //    Tells the kernel during run-time, that the osc-bit
        //  should not be sent within the next call(s) of the run-method.
        //  If true is passed as argument, the osc bit is NOT sent.
        //    setUseSyncBits must have been called with true as 1st
        //  argument. Otherwise the flag is ignored.
        void setDoNotSendOscBit(bool bValue = true);

        //--------------------------------------------------------------------
        //    Tells the kernel during run-time, that the external
        //  trigger bit should not be sent within the next call(s)
        //  of the run-method. If true is passed as
        //    argument, the external trigger bit is NOT sent.
        //    setUseSyncBits must have been called with true as 1st
        //  argument. Otherwise the flag is ignored.
        void setDoNotSendExtTrigger(bool bValue = true);

        //--------------------------------------------------------------------
        //    Tells the kernel whether it should run internal phase
        //  correction scans or not.
        //    Deactivates early FID phase correction scans.
        virtual void setInternalPhaseCorrection(bool bValue);

        //--------------------------------------------------------------------
        //    Tells the kernel whether it should run early FID phase
        //  correction scans or not.
        //    Deactivates internal phase correction scans.
        virtual void setEarlyFIDPhaseCorrection(bool bValue);

        //--------------------------------------------------------------------
        //    Tells the kernel whether it should run full train phase
        //  correction scans or not.
        //    Deactivates early FID phase correction scans.
        virtual void setExternalEPIPhaseCorrection(bool bValue);

        //--------------------------------------------------------------------
        //  setExecuteKernelAsPhaseCorrectionScan
        //    If no internal phase correction scans are acquired,
        //  the sequence has to use this method to switch the
        //  kernel into the phase correction mode. This method
        //    also has to be used to switch the kernel back to the
        //  imaging mode. Internal phase correction must be disabled.
        virtual void setExecuteKernelAsPhaseCorrectionScan(bool bValue);

        //--------------------------------------------------------------------
        //    Set total number of contrasts (-> echo trains) to acquire.
        void setNumberOfContrasts( long lNumberOfContrasts );

        //--------------------------------------------------------------------
        //    Set number of contrasts to acquire before the RTEB-plug-in
        void setNumberOfContrastsBeforeRTEBPlugIn( long lNumberOfContrastsBeforeRTEBPlugIn );

        //--------------------------------------------------------------------
        //    Set contrast index which realizes TE
        void setTEContrastIndex( long lTEContrastIndex );

        //--------------------------------------------------------------------
        //    Get actual TE for given contrast index (valid after preparation)
        long getActualTE( size_t lContrastIndex ) const;

        //--------------------------------------------------------------------
        //    Overloaded base class method: ensure reset of TE's
        inline void resetPrepared();

        //--------------------------------------------------------------------
        //    Tell the kernel wether the RO gradient should be
        //  prephased after the RTEB-plug-in or not. If not,
        //  it is prephased before the RTEB-plug-in.
        //    A RTEB-plug-in must have been registered.
        void setPrephaseROAfterRTEBPlugIn(bool bValue);
        bool getPrephaseROAfterRTEBPlugIn () const;

        //--------------------------------------------------------------------
        //  setPrephaseBlipsAfterRTEBPlugIn
        //    Tell the kernel whether the PE/3D-blips should be
        //  prephased after the RTEB-plug-in or not. If not,
        //  they are prephased before the RTEB-plug-in.
        //    A RTEB-plug-in must have been registered.
        void setPrephaseBlipsAfterRTEBPlugIn(bool bValue);
        bool getPrephaseBlipsAfterRTEBPlugIn () const;

        //--------------------------------------------------------------------
        //    If true is passed to this function, the gradients of
        //  the gradient moment
        //    section are calculated so that they always can be sent
        //  independent of the
        //    current slice orientation. Otherwise the maximum
        //  gradient moment on one axis
        //    is used to calculate the timing of the gradient moment
        //  sections, which may
        //    lead to gradient overflows, when slices are angulated.
        void setRotationProofGradientMoments(bool bValue);

        //--------------------------------------------------------------------
        //    Advise the kernel to start the imaging EPI echo train
        //  with a negative
        //    gradient or not (i.e. positive).
        void setStartImagingReadOutWithNegativeGradient(bool bValue);
        bool getStartImagingReadOutWithNegativeGradient () const;

        //--------------------------------------------------------------------
        //    Tells the kernel the number of the echo within the
        //  imaging echo train which contains the ADC with the
        //  k-space-center line- and partition-number starting
        //    to count at zero. This number is important for the
        //  correct calculation of the TE-fill times.
        void setCenterSegment(long lCenterSegment = 0);

        //--------------------------------------------------------------------
        //    Activates or deactivates echo-shifting. If echo
        //  shifting should be activated it is important for
        //  the kernel to know the number of (lines res. partitions)
        //    counters per segment and the (line res. partition)
        //  counter within the center segment with which the k-space
        //  center is measured.For deactivation only false needs to
        //  be passed to this function.
        virtual void setUseEchoShifting(bool bValue, long lCountersPerSegment = 0, long lCounterInSegmentWithEcho = 0);

        //--------------------------------------------------------------------
        //    Tells the kernel, if the PE- and 3D-blips should be
        //  rewound after acquiring the imaging scans.
        void setRewindBlips(bool bValue);

        //--------------------------------------------------------------------
        //    Tells the kernel, if the RO-gradient should be
        //  rewound after acquiring the imaging scans.
        void setRewindRO(bool bValue);

        //--------------------------------------------------------------------
        //    Tells the kernel, if the RO-gradient should be
        //  spoiled after acquiring the imaging scans.
        void setSpoilRO(bool bValue);

        //--------------------------------------------------------------------
        //    Tells the kernel, if the GS-gradient should be prefaced
        // for the next kernel after acquiring the imaging scans.
        void setPrefaceExcitationSBB(bool bValue);

        //--------------------------------------------------------------------
        //    Can be used to specify an additional phase to the NCO
        //  during the measurement.Note that this will overwrite any
        //  additional phase specified for the excitation SBB within
        //  the sequence before calling this method!
        void setAdditionalPhase(double dValue    // the additional phase in deg
                                        ) override;

        //--------------------------------------------------------------------
        //    Can be used to specify a certain desired TE for the
        //  Kernel BEFORE prep is called. After prep was called with
        //  success, the method increaseTE must be used.
        //  If no wanted TE is specified or if it is not possible to
        //  execute the measurement with the wanted TE, the Kernel is
        //  prepared with the minimum possible TE.
        void setWantedTE(long lTE);

        //--------------------------------------------------------------------
        //    Prepares the kernel. At least an excitation SBB configuration
        //  function must have been registered. All other configuration
        //  steps have to be performed.
        bool prepSBB(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;



        // prepare the plugin; this method can be overloaded by all kernel variants
        // which use a plug-in (e.g. diffusion, SE and ASL)
        virtual bool prepPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        // init the excitation SBBs
        virtual bool initExcitation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        //--------------------------------------------------------------------
        //    Prepares the additional Gradients,e.g. flow compensated gradients
        virtual bool prepAdditionalGradients(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        //--------------------------------------------------------------------
        //    adjust TE contributions due to the applied additional Gradients,e.g. flow compensated gradients
        virtual bool adjustTEContributions(MrProt &rMrProt);

        //    adjust additional timing config due to the applied additional Gradients,e.g. flow compensated gradients
        virtual bool adjustAdditionalTiming(MrProt &rMrProt);

        //    adjust local TR fill due to the applied additional Gradients,e.g. flow compensated gradients
        virtual bool adjustLocalTRFill(MrProt &rMrProt);

        // *--------------------------------------------------------------------*
        // * run flow compensated slice refocusing gradient                     *
        virtual bool runAdditionalGradExcitation(long &lEventBlockEnd);

#ifdef ZOOM_2DRF
        // prepare the plugin; this method can be overloaded by all kernel variants
        // which use a plug-in (e.g. diffusion, SE and ASL)
        virtual bool prepZOOMitExcitation(MrProt &rMrProt,const SeqLim &rSeqLim);
        virtual bool setZOOMitParameterFromUI(MrProt &rMrProt, const SeqLim &rSeqLim);

        virtual bool configureSBB2DExc(MrProt& rMrProt, const SeqLim&);
#ifdef ZOOM_EXTENDED
        virtual bool updateZOOMitRotation(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);
#endif
#endif //ZOOM_2DRF

        // get excitation pointer
        virtual SeqBuildBlockExcitation* getExcitationPointer() { return m_pSBBExcite; };

        //--------------------------------------------------------------------
        //    Checks all gradients involved for gradient specification violations
        // in the LOGICAL coordinate system.
        bool checkGradients(MrProt &rMrProt, SeqLim &rSeqLim) override;

        //--------------------------------------------------------------------
        //	Overloaded base class method
        //  Set gradient performance:
        //  - If no GPA balance model is used: call base class method
        //  - If GPA model is activated: call base class method and change
        //    maximum amplitude afterwards
        void setDefaultGradientPerformance() override;

        //--------------------------------------------------------------------
        //	Enable use of GPA balance model for gradient pulse calculations.
        //  If enabled, the maximum usable gradient amplitude will be increased
        //  and a balance calculation ensures that the gradient events stay 
        //  within the GPA limitations IF an appropriate delay time is applied 
        //  by the sequence (see getTRIncrement).
        //
        //  If dGradMaxAmpl [mT/m] is provided, it is used as the upper limit 
        //  of the readout gradients. It is not recommended to use values above 
        //  90% of GradMaxAmplAbsolute
        //
        //  Note: do not use GPA balance model for 3D-EPI
        //  Note: 
        virtual bool setUseGPABalance(bool bUseGPABalance, double dGradMaxAmpl = -1.);

        //--------------------------------------------------------------------
        //	Check whether use of GPA balance model is enabled
        virtual bool getUseGPABalance() const;

        //--------------------------------------------------------------------
        //	Set balance values of preceding gradient events (e.g. 
        //  contrast preparation modules) if these should be considered
        virtual bool setBalanceIn(GPABalance::GPABalanceValue sBalanceIn, GPABalance::GrmsValue sBalanceRmsIn);

        //--------------------------------------------------------------------
        //	Get full access to internal balance instance (m_sSBBBalance)
        //  Result is valid only after successful ::prep
        virtual bool getBalance(GPABalance* &psBalance);

        //--------------------------------------------------------------------
        // Trigger sending of slice adjust SEQData objects during ::run()
        // Might get changed during sequence execution
        virtual void setSendSliceAdjustData( bool bSendSliceAdjustData );

        //--------------------------------------------------------------------
        //	Get TR increment required per kernel execution.  
        //  This pause has to be applied once per kernel request by
        //  the sequence.
        //  Prerequisite: calcTRIncrement has been called with success
        virtual long getTRIncrement() const;

        //--------------------------------------------------------------------
        //    Returns duration of the SBB of one call of the runSBB function.
        long getSBBDurationPerRequest() override;

        //--------------------------------------------------------------------
        //    Returns the actual energy per request. Needs to be overloaded
        //    by the individual kernel flavors.
        MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest() override;

        //--------------------------------------------------------------------
        //    Returns the actual TE of the kernel or the minimum TE required
        //  to execute the measurement with the current protocol parameters.
        long getNeededTE() const;

        //--------------------------------------------------------------------
        //    Returns the kernel duration left after the center of the echo
        //  which acquires the k-space center line/partition.
        long getDurationAfterEcho() const;

        //--------------------------------------------------------------------
        //    Increases the current TE to a certain value.
        virtual bool increaseTE(long lNewTE);

        //--------------------------------------------------------------------
        //    Sets the TR-fill the kernel should apply.
        void setTRFill(long lTRFill);

        //--------------------------------------------------------------------
        //    Executes the real time part of the SBB.
        bool runSBB(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC) override;

        //--------------------------------------------------------------------
        //  Set number of interleaves this information is currently
        // only used in the effective echo-spacing calculation
        virtual void setNumInterleaves(long lNumInterleaves);

        //--------------------------------------------------------------------
        //  Calculate: (1) An effective echo-spacing, which takes number of interleaves
        //                 and PAT acceleration factor into account
        //             (2) The bandwidth per pixel in the phase-encoding direction
        //                 FOR THE RECONSTRUCTED IMAGE
        //
        virtual bool calcEffEchoSpacingAndBWPerPixelPE(MrProt &rMrProt, long& lEffectiveEchoSpacing, double& dBandwidthPerPixelPE);

        //--------------------------------------------------------------------
        // get navigator start time
        virtual long getNavigatorStartTime() const;

        //--------------------------------------------------------------------
        //    Internal/External function to run the sync-bit and the
        //    excitation SBB. Calls runSync Bits.
        virtual bool runExcitationAndSyncBits(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

    protected:

        //--------------------------------------------------------------------
        //    Internal function to update the performance data of
        //  all gradients involved within the kernel.
        void updateGradientPerformance(SEQ::Gradients eGradMode) override;

        //--------------------------------------------------------------------
        //    Internal function to calculate the timing of the gradients
        //  realizing gradient moments between the SBBs used by the kernel.
        virtual bool calculateTimingOfGradMoments(MrProt &rMrProt, SeqLim &rSeqLim);

        //--------------------------------------------------------------------
        //	Calculate TR increment necessary to repeat the readout train
        //  forever. Used only if no plugin is specified (otherwise, the
        //  plugin has to take care of the corresponding calculation).
        //  TE fill time is considered appropriately.
        //  Prerequisite: checkBalance has been called with success OR
        //                appropriate plugin has been prepared
        virtual bool calcTRIncrement(long lTEFillTime);

        //--------------------------------------------------------------------
        //	Check whether balance allows to apply the readout train
        //  at least once
        //  Prerequisite: EPI readout has to be prepared
        virtual bool checkBalance();

        //--------------------------------------------------------------------
        //    Internal function to run one gradient moment section.
        virtual bool runGradMoments(long lSectionIndex, MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

        //--------------------------------------------------------------------
        //    update additional gradient moments, e.g.flow compensated gradients
        virtual bool updateAdditionalGradMoments(long lSectionIndex, MrProt &rMrProt, SeqLim &rSeqLim, long &lFillTimeBeforeGradients, long &lEventBlockEnd);

        //--------------------------------------------------------------------
        //    run gradient moment
        virtual bool runFinalGradMoments(long lSectionIndex, MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC, long lFillTimeBeforeGradients, long lEventBlockEnd);
        
        //--------------------------------------------------------------------
        //    Internal function to run the sync-bits.
        virtual bool runSyncBits(MrProt&, SeqLim&, SeqExpo&, sSLICE_POS*);

        //--------------------------------------------------------------------
        //    Internal function to run the internal phase correction echos.
        virtual bool runInternalPhaseCorrection(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

        //--------------------------------------------------------------------
        //    Internal function to run the early FID phase correction echos.
        virtual bool runEarlyFIDPhaseCorrection(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

        //--------------------------------------------------------------------
        //    Internal function to run the phase correction EPI read-out and the
        //    appropriate fill times, if the kernel is in the phase correction mode.
        virtual bool runEPIReadOutPhaseCorrection(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

        //--------------------------------------------------------------------
        //    Internal function to run the imaging EPI read-out.
        virtual bool runEPIReadOut(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

        //--------------------------------------------------------------------
        //    Internal function to run the RTEB-plug-in.
        virtual bool runRTEBPlugIn(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

#ifdef ZOOM_2DRF
        //--------------------------------------------------------------------
        // check if the conditions of OptPTXVolume are present (different in diffusion and others)
        virtual bool isOptPTXVolumeCondition(MrProt& rMrProt);
#endif

        //--------------------------------------------------------------------
        // Enlarge the max amplitude of Gradient Slice of ep2d_diff to allow thinner slice thickness (3mm) for low-field
        // (only in derived Diffusion kernel)
        virtual void increaseSliceSelectionGradientMaxAmplitude(MrProt& rMrProt);

        //--------------------------------------------------------------------
        // Data Members for Class Attributes

        bool m_bUseSyncBits{false};

        long m_lMaxSyncBitDuration{0};

        bool m_bExcitationSBBInternallyDefined{false};

        long m_lSBBExciteStartTime{0};

        bool m_bInternalPhaseCorrection{true};

        bool m_bExternalEPIPhaseCorrection{false};

        bool m_bEarlyFIDPhaseCorrection{false};

        bool m_bRewindBlips{false};

        bool m_bRewindRO{false};

        // Apply spoiler (rather than rewinder) along readout direction
        bool m_bSpoilRO{false};

        // Readout spoil moment factor (relative to the rewinding moment)
        double m_dSpoilROFactor{2.0};

        bool m_bPrefaceExcitationSBB{false};

        bool m_bExecuteKernelAsPhaseCorrectionScan{false};

        bool m_bPlugInAvailable{false};

        bool m_bRotationProofGradientMoments{true};

        bool m_bPrephaseROAfterRTEBPlugIn{false};

        bool m_bPrephaseBlipsAfterRTEBPlugIn{false};

        long m_lRTEBPlugInDurationPerRequest{0};

        bool m_bRTEBPlugInInvertsMagnetization{false};

        long m_lRTEBPlugInTEContribution{0};

        long m_lRTEBPlugInStorageTime{0};

        long m_lRTEBPlugInTEFillBefore{0};

        long m_lRTEBPlugInTEFillAfter{0};

        // Total number of echo-trains (at least one)
        long m_lNumberOfContrasts{1};

        // Number of echo-trains before RTEB-PlugIn
        long m_lNumberOfContrastsBeforeRTEBPlugIn{0};

        // Contrast index which realizes TE
        long m_lTEContrastIndex{0};

        // Actual TE for each of the acquired contrasts (valid after successful preparation)
        std::vector<long> m_vlActualTE;

        long m_lWantedTE{0};

        long m_lNeededTE{0};

        long m_lTRFill{0};

        long m_lDurationAfterEcho{0};

        long m_lCenterSegment{0};

        long m_lCounterInSegmentWithEcho{0};

        long m_lNumInterleaves{0};

        bool m_bEPIReadOutAppliesTEFill{false};

        bool m_bEPIReadOutAppliesTRFill{false};

        bool m_bGradMoment2AppliesTEFill{false};

        bool m_bGradMoment3AppliesTRFill{false};

        bool m_bStartImagingReadOutWithNegativeGradient{false};

        bool m_bDoNotSendOscBit{false};

        bool m_bDoNotSendExtTrigger{false};

        bool m_bUseGPABalance{false};

        double m_dGradMaxAmpl{-1.0};

        long m_lTRIncrement{0};

        GPABalance::GPABalanceValue m_sBalanceIn;

        GPABalance::GrmsValue       m_sBalanceRmsIn;

        GPABalance                  m_sSBBBalance;

        std::array<sGRAD_PULSE_TRAP, 4> m_aGradMoments;

        SeqBuildBlockExcitation *m_pSBBExcite;

        sSYNC_OSC m_OscBit;

        sSYNC_EXTTRIGGER m_ExtTrig;

        // Sequence to Ice communication
        sSYNCDATA   m_sSliceAdjustSyncData{"SliceAdjustSyncData"};
        SEQData     m_sSliceAdjustSEQData{};
        bool        m_bSendSliceAdjustData{false};

        sGRAD_PULSE m_sGrad[3];
        const long  m_lGP{0};
        const long  m_lGR{1};
        const long  m_lGS{2};
        char        m_tIdent[64];

        long m_lNavigatorStartTime{0};

        // Simultaneous multislice imaging active
        bool m_bIsSliceAcceleration{false};

        // TE contribution times
        long m_lTEContributionBeforeRTEBPlugIn{0};
        long m_lTEContributionAfterRTEBPlugIn {0};

        // PlugIn to end of echo train(start of compensation gradients)
        long m_lRTEBPlugInToCompGradTime{0};

        //-------------------------------------------------------------------------------------
        // Excitation Sequence Building Blocks
        //-------------------------------------------------------------------------------------
        SeqBuildBlockBinomialPulses m_SBBExcitationBinomial{nullptr, true, 180.0, 3, "SBBExcBinom"};
        SBBMultibandRF              m_SBBExcitation{nullptr};

#ifdef ZOOM_2DRF
        SBB2DExc                           m_2DExcite{nullptr, "Ptx2DExc"};
        RF2DArbGenerator::RF2DArbParameter m_Param; // configuration struct for SBB2DExc class

        SBB2DPtxPulsesIni m_2DExciteIniMultiCh{nullptr, PtxRfPulseIdEpiExc, "IniNch"}; // iRFPulseID = 10
        SBB2DExcPulsesIni m_2DExciteIniOneCh{nullptr, PtxRfPulseIdEpiExc, "Ini1ch"}; // iRFPulseID = 10

        // parameters for m_2DExcite method
        MrProtocolData::PTXTrajectoryType   m_ExcType         {MrProtocolData::PTXTrajectoryType_EPI_1D};
        MrProtocolData::PTXB0CorrectionType m_B0CorrectionType{MrProtocolData::PTXB0CorrectionType_Off};
        std::vector<double>                 m_vdOptPTXVolumeThickness;
        std::vector<double>                 m_vdOptPTXVolumeShift;
#ifdef ZOOM_EXTENDED
        // flag to use DebugSettings
        bool m_bUseDebugSettings{false};

        double m_dRotationAngleDegDefault   {5.0};
        double m_dRotationAngleDegLowerLimit{2.0};
        double m_dRotationAngleDegUpperLimit{8.0};
#endif // ZOOM_EXTENDED
#endif //ZOOM_2DRF
    };
}



//--------------------------------------------------------------------
//    Tells the kernel during run-time, that the osc-bit should
//  not be sent within the next call(s) of the run-method.
//  If true is passed as argument, the osc bit is NOT sent.
//    setUseSyncBits must have been called with true as 1st
//  argument. Otherwise the flag is ignored.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setDoNotSendOscBit(bool bValue)
{
    m_bDoNotSendOscBit=bValue;
}


//--------------------------------------------------------------------
//    Tells the kernel during run-time, that the external
//  trigger bit should not be sent within the next call(s) of
//  the run-method. If true is passed as argument, the external
//  trigger bit is NOT sent.
//  setUseSyncBits must have been called with true as 1st argument.
//  Otherwise the flag is ignored.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setDoNotSendExtTrigger(bool bValue)
{
    m_bDoNotSendExtTrigger=bValue;
}


//--------------------------------------------------------------------
//    Tells the kernel whether it should run internal phase correction
//  scans or not. Deactivates early FID phase correction scans.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setInternalPhaseCorrection(bool bValue)
{

    m_bInternalPhaseCorrection = bValue;
    
    if(bValue)
    {
        // Early FID phase correction is not compatible with internal phase correction
        m_bEarlyFIDPhaseCorrection = false;
    }
    resetPrepared();
}

//--------------------------------------------------------------------
// Tells the kernel whether it should run full train phase correction
// scans or not. Deactivates early FID phase correction scans.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setExternalEPIPhaseCorrection(bool bValue)
{
    m_bExternalEPIPhaseCorrection = bValue;

    if (bValue)
    {
        
        m_bExecuteExternalEPIPhaseCorrectionScan = true;
        // Early FID phase correction is not compatible with full train phase correction
        m_bEarlyFIDPhaseCorrection = false;
    }

    resetPrepared();
}


//--------------------------------------------------------------------
//    Tells the kernel whether it should run early FID phase correction
//  scans or not.Deactivates internal phase correction scans.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setEarlyFIDPhaseCorrection(bool bValue)
{
    m_bEarlyFIDPhaseCorrection = bValue;

    if(bValue)
    {
        m_bExecuteExternalEPIPhaseCorrectionScan = false;
        m_bExternalEPIPhaseCorrection = false;
        m_bInternalPhaseCorrection = false;
    }

    resetPrepared();
}


//--------------------------------------------------------------------
//    If no internal phase correction scans are acquired, the
//  sequence has to use this method to switch the kernel
//  into the phase correction mode. This method
//    also has to be used to switch the kernel back to the
//  imaging mode.Internal phase correction must be disabled.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setExecuteKernelAsPhaseCorrectionScan(bool bValue)
{
    m_bExecuteKernelAsPhaseCorrectionScan = bValue;
}

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setNumberOfContrasts( long lNumberOfContrasts )
{
    m_lNumberOfContrasts = std::max( lNumberOfContrasts, 0L );
    resetPrepared();
}

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setNumberOfContrastsBeforeRTEBPlugIn( long lNumberOfContrastsBeforeRTEBPlugIn )
{
    m_lNumberOfContrastsBeforeRTEBPlugIn = std::max( lNumberOfContrastsBeforeRTEBPlugIn, 0L );
    resetPrepared();
}

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setTEContrastIndex( long lTEContrastIndex )
{
    m_lTEContrastIndex = lTEContrastIndex;
    resetPrepared();
}

inline long SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getActualTE( size_t lContrastIndex ) const
{
    if ( m_vlActualTE.size() == 0 )
    {
        return 0;
    }
    lContrastIndex = std::max<size_t>( 0, std::min( lContrastIndex, m_vlActualTE.size() - 1 ) );
    return m_vlActualTE[lContrastIndex];
}

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::resetPrepared( void )
{
    m_vlActualTE.clear();
    SeqBuildBlockEPIReadOut::resetPrepared();
}

//--------------------------------------------------------------------
//    Tell the kernel whether the RO gradient should be prephased
//  after the RTEB-plug-in or not. If not, it is prephased before the
//  RTEB-plug-in.A RTEB-plug-in must have been registered.
inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setPrephaseROAfterRTEBPlugIn(bool bValue)
{
    m_bPrephaseROAfterRTEBPlugIn=bValue;
    resetPrepared();
}

inline bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getPrephaseROAfterRTEBPlugIn( void ) const
{
    return m_bPrephaseROAfterRTEBPlugIn;
}

//--------------------------------------------------------------------
//    Tell the kernel whether the PE/3D-blips should be prephased
//  after the RTEB-plug-in or not. If not, they are prephased before the
//  RTEB-plug-in. A RTEB-plug-in must have been registered.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setPrephaseBlipsAfterRTEBPlugIn(bool bValue)
{
    m_bPrephaseBlipsAfterRTEBPlugIn=bValue;
    resetPrepared();
}

inline bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getPrephaseBlipsAfterRTEBPlugIn() const
{
    return m_bPrephaseBlipsAfterRTEBPlugIn;
}

//--------------------------------------------------------------------
//    If true is passed to this function, the gradients of the
//  gradient moment section are calculated so that they always
//  can be sent independent of the current slice orientation.
//  Otherwise the maximum gradient moment on one axis is used
//  to calculate the timing of the gradient moment sections,
//  which may lead to gradient overflows, when slices are angulated.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setRotationProofGradientMoments(bool bValue)
{
    m_bRotationProofGradientMoments=bValue;
    resetPrepared();
}


//--------------------------------------------------------------------
//    Advise the kernel to start the imaging EPI echo train with
//  a negative gradient or not (i.e. positive).

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setStartImagingReadOutWithNegativeGradient(bool bValue)
{
    m_bStartImagingReadOutWithNegativeGradient=bValue;
}

inline bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getStartImagingReadOutWithNegativeGradient() const
{
    return m_bStartImagingReadOutWithNegativeGradient;
}

//--------------------------------------------------------------------
//    Tells the kernel the number of the echo within the imaging
//  echo train which contains the ADC with the k-space-center line-
//  and partition-number starting to count at zero. This number is
//  important for the correct calculation of the TE-fill times.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setCenterSegment(long lCenterSegment)
{
    m_lCenterSegment=lCenterSegment;
    resetPrepared();
}


//--------------------------------------------------------------------
//    Tells the kernel, if the PE- and 3D-blips should be
//  rewound after acquiring the imaging scans.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setRewindBlips(bool bValue)
{
    m_bRewindBlips=bValue;
    if ( m_bRewindRO == true )
    {
        m_bSpoilRO = false;
    }
    resetPrepared();
}


//--------------------------------------------------------------------
inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setSpoilRO(bool bValue)
{
    m_bSpoilRO = bValue;
    if ( m_bSpoilRO == true )
    {
        m_bRewindRO = false;
    }
    resetPrepared();
}
//    Tells the kernel, if the RO-gradient should be rewound
//  after acquiring the imaging scans.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setRewindRO(bool bValue)
{
    m_bRewindRO=bValue;
    resetPrepared();
}


//--------------------------------------------------------------------
//    Tells the kernel, if the GS-gradient should be prefaced
//  for the next kernel after acquiring the imaging scans.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setPrefaceExcitationSBB(bool bValue)
{
    m_bPrefaceExcitationSBB=bValue;
    resetPrepared();
}


//--------------------------------------------------------------------
//    Can be used to specify a certain desired TE for the Kernel
//  BEFORE prep is called. After prep was called with success,
//  the method increaseTE must be used. If no wanted TE is
//  specified or if it is not possible to execute the
//    measurement with the wanted TE, the Kernel is prepared
//  with the minimum possible TE.

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setWantedTE(long lTE)
{
    m_lWantedTE = std::max(0L, lTE);
    resetPrepared();
}


//--------------------------------------------------------------------
//    Returns the actual TE of the kernel or the minimum TE required
//  to execute the measurement with the current protocol parameters.

inline long SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getNeededTE() const
{
    return isPrepared() ? m_lNeededTE : 0;
}

//--------------------------------------------------------------------
//	Enable use of GPA balance model for gradient pulse calculations
inline bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setUseGPABalance(bool bUseGPABalance, double dGradMaxAmpl)
{
    // Check whether gradient system supports balance models
    if(bUseGPABalance && (m_sSBBBalance.lGetStatus() != 0))
    {
        m_bUseGPABalance = false;
        m_dGradMaxAmpl   = -1.;

        return false;
    }

    if((bUseGPABalance != m_bUseGPABalance) || (dGradMaxAmpl != m_dGradMaxAmpl))
    {
        // Mode has changed
        m_bUseGPABalance = bUseGPABalance;
        m_dGradMaxAmpl   = dGradMaxAmpl;

        // Update gradient limits
        setDefaultGradientPerformance();

        resetPrepared();
    }

    return true;
}

//--------------------------------------------------------------------
//	Check whether use of GPA balance model is enabled
inline bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getUseGPABalance() const
{
    return m_bUseGPABalance;
}

//--------------------------------------------------------------------
//	Set balance values of preceding gradient events (e.g. 
//  contrast preparation modules) if these should be considered
inline bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setBalanceIn(GPABalance::GPABalanceValue sBalanceIn, GPABalance::GrmsValue sBalanceRmsIn)
{
    m_sBalanceIn    = sBalanceIn;
    m_sBalanceRmsIn = sBalanceRmsIn;
    resetPrepared();

    return true;
}

//--------------------------------------------------------------------
//	Get full access to internal balance instance (m_sSBBBalance)
inline bool SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getBalance(GPABalance* &psBalance)
{
    if(isPrepared() && getUseGPABalance())
    {
        psBalance = &m_sSBBBalance;
        return true;
    }

    psBalance = nullptr;
    return false;
}

//--------------------------------------------------------------------
inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setSendSliceAdjustData( bool bSendSliceAdjustData )
{
    m_bSendSliceAdjustData = bSendSliceAdjustData;
}
//	Get TR increment required per kernel execution
 inline long SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getTRIncrement() const
{
    if(isPrepared() && getUseGPABalance())
    {
        return m_lTRIncrement;
    }

    return 0;
}

//--------------------------------------------------------------------
//    Returns the kernel duration left after the center of the echo
//  which acquires the k-space center line/partition.

inline long SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getDurationAfterEcho() const
{
    return isPrepared() ? m_lDurationAfterEcho : 0;
}


//--------------------------------------------------------------------
//  Sets the number of interleaves.
//  Only used in effective echo-spacing calculation.
//  NO effect on EPI readout gradient waveform

inline void SEQ_NAMESPACE::SeqBuildBlockEPIKernel::setNumInterleaves(long lNumInterleaves)
{
    m_lNumInterleaves = lNumInterleaves;
}


//--------------------------------------------------------------------
//  Return the navigator start time
inline long SEQ_NAMESPACE::SeqBuildBlockEPIKernel::getNavigatorStartTime() const
{
    return m_lNavigatorStartTime;
}
