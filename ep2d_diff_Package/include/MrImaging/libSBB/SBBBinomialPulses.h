//----------------------------------------------------------------------------------
// <copyright file="SBBBinomialPulses.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 1998-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015-2020. All Rights Reserved. Confidential.
// </copyright>
// <description>
//   SBB for adiabatic binomial pulses
// </description>
//----------------------------------------------------------------------------------

#pragma once

#ifndef SBBBinomialPulses_h
#define SBBBinomialPulses_h

#include "MrImaging/libSBB/SBBExcitation.h"
#include "MrImaging/libSBB/libSBBmsg.h"
#include "MrMeasSrv/SeqIF/csequence.h" // for CHEMICAL_SHIFT_FAT_PPM
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"

#ifdef BUILD_libSBB
#define __OWNER
#endif

#ifdef SEQLIB_STATIC
#define __STATIC
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"

#define MAX_NO_RF_PULSES 10
#define MIN_NON_SELECTIVE_RF_DURATION 100
#define MAX_SIZE_FLIPANGLEARRAY_SBBBINOMIAL 30

#define SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF (SBBExcitation_GRAD_PERF_GS) // define kept for BACKWARDS COMPATIBILITY
#define SBBBinomialPulses_GRAD_PERF_BINOMIAL (SBBExcitation_NO_OF_GRAD_PERF_GROUPS + 0)
#define SBBBinomialPulses_NO_OF_GRAD_PERF_GROUPS (SBBExcitation_NO_OF_GRAD_PERF_GROUPS + 1)

#define MIN_RF_SPACING (GRAD_RASTER_TIME) // after asking F.Rabe

#include <algorithm>
#include <array>

class SBBTestHelperBinomialPulses; // SBB iTest Helper Class to Dump Class Member

// ----------------------------------------------------------------------------
/// \ingroup grp_libsbb
/// \brief  SBB for adiabatic binomial pulses
/// \todo add documentation
///    This building block provides binomial pulses which can be used to excite the
///    on-resonant component of the spin ensemble while leaving the magnetization of
///    an off-resonant component aligned along the longitudinal axis of the
///    rotating frame and vise versa.
///    Also a non-frequency-selective excitation RF-pulse can be sent by this SBB
///    so that this SBB can replace the complete excitation part of a sequence
///    kernel. This non-frequency-selective RF-pulse pulse is realized by the class
///    this SBB is derived from.
///    A graphical representation of the SBB can be found in
///    /src/MrImaging/libSBB/SBBBinomialPulses.pdf.
///
///
///
///    The RF-spacing between the binomial RF-pulses is:
///    RFSpacing[us] = (1.0E6/360.0) * PhaseEvolutionWanted[deg] / (Chemical
///    Shift[ppm] * dNominalB0[T] * GAMMA[MHz/T])
///    Note that this RF-spacing must be rounded to the gradient raster time and
///    therefore the true phase evolution may differ slightly from the wanted one.
///
///
///
///    This SBB is dedicated to do fat suppression using water excitation or water
///    suppression using fat excitation. The SBB itself decides (within its virtual
///    method setPrepAndRunMode) which sort of pulse it has to apply depending on
///    the following protocol parameters:
///    - pMrProt->getsPrepPulses().getucFatSat()
///    - pMrProt->getsPrepPulses().getucWaterSat()
///    - pMrProt->getsTXSPEC().getucExcitMode()
///    If setPrepAndRunMode decided, that no binomial pulses should be applied,
///    normal non frequency selective excitation will be performed using the
///    methods of the excitation-RF base class.
///
///
///
///    This SBB offers the following configuration possibilities for binomial
///    pulses (a detailed explanation follows):
///    - variable off-resonant component
///    - variable number of pulses for the binomial rf-pulse train
///    - bipolar or monopolar gradients
///    - variable phase evolution between binomial rf-pulses
///    - variable additional phase
///    - variable bandwidth-time-product
///    - run with possible bandwidth-time-product
///    - optional flip angle adaption
///    - slice thickness
///    - variable non selective RF-pulse duration
///    - variable offset frequency
///    - variable RF-pulse-type
///    - variable polarity of GS-gradients
///
///    variable off-resonant component:
///    The off-resonant component is specified by its chemical shift in ppm. The
///    default is -3.3ppm.
///
///    variable number of pulses for the binomial rf-pulse train:
///    The binomial pulse train can consist of 2 up to MAX_NO_RF_PULSES (currently
///    10) rf pulses. The default is 3.
///
///    bipolar or monopolar gradients:
///    It can be chosen, if the slice selective gradients should be rewound
///    between the rf-pulses and the rf-pulses are always sent over a positive
///    gradient. In this case the SBB sends a "monopolar" binomial pulse train. The
///    alternative is to send every other rf-pulse over a negative gradient. This
///    "bipolar" pulse train allows to excite thinner slices or to specify a higher
///    slice profile quality compared to the monopolar one.
///
///    variable phase evolution between binomial rf-pulses:
///    The phase evolution of the binomial pulses is the phase the off-resonant
///    component collects during the time between two rf-pulses in degree. The
///    default is 180deg.
///
///    variable bandwidth-time-product:
///    The slice selective binomial pulses are sinc-Pulses (class sRF_PULSE_SINC).
///    The requested slice profile quality of the binomial pulses is specified
///    using the bandwidth-time-product of a sinc-pulse. The default is 2.0.
///
///    run with possible bandwidth-time-product:
///    The bandwidth time product of the sinc-pulses restricts the minimum slice
///    thickness. The client can choose to use the maximum possible slice-profile
///    quality for the requested slice thickness, if the desired quality can not be
///    achieved within in the available time. If this option is not set (which is
///    the default), the prep-function returns an error, if the slice thickness can
///    not be realized with the desired quality.
///
///    optional flip angle adaption:
///    Using binomial coefficients the flip angle selected is distributed over the
///    rf-pulses. Only for a phase evolution of 180deg the resulting flip angle for
///    the excited component is equal to the selected one. For phase evolutions
///    less than 180deg the resulting flip angle is less than the selected one. So
///    the SBB tries to increase the flip angle of each rf-pulse so that the
///    resulting angle becomes as close as possible to the desired one. The SBB
///    uses an iterative algorithm and the simple vector model with rotation
///    matrices for this adaption. Per default this algorithm is invoked.
///
///    slice thickness:
///    Thickness of the slices that should be excited in mm.
///
///    variable non selective RF-pulse duration:
///    The client can specify the duration in us of the non-slice-selective
///    rectangular rf-pulses used for frequency selective excitation. Default is
///    300us.
///
///    variable offset frequency:
///    An offset-frequency in Hz can be specified. This frequency results in an
///    additional phase applied to the binomial rf-pulses. This phase scales linear
///    with the specified frequency, the rf-spacing and the number of the rf-pulse.
///    Default is zero.
///
///    variable RF-pulse-type
///    If the user does not want to use the RF-pulse types offered by this class
///    (sinc/rect) or if the user is not satisfied with timing calculated for the
///    sinc-RF-pulses, then any pulse type derived from sRF_PULSE can be used for
///    excitation. The client has to create one instance of the RF-pulse for EACH
///    RF-pulse of the binomial rf-pulse train and do the basic configuration steps
///    for EACH pulse. Then EACH pulse must be registered at the SBB. Note that the
///    following attributes of sRF_PULSE are set by this class during preparation
///    of the SBB:
///    - ident
///    - thickness
///    - initial phase
///    - flip-angle
///    - required GS-polarity
///    The following sRF_PULSE-attributes must be set by the client:
///    - duration
///    - RF-pulse specific attributes (like bandwidth-time-product and no. of
///    samples for sRF_PULSE_SINC or family name for sRF_PULSE_EXT)
///    If the sequence wants to switch back to the usage of internal RF-pulses
///    again, then the RF-pulses must be unregistered.
///
///    variable polarity of GS-gradients
///    If the user wants to run the binomial pulses over a negative GS gradient,
///    then this can be done by calling the appropriate set-method. For bipolar
///    gradients the order of positive and negative gradients is inverted.
///
///    Flip angle array support
///    like the base class SeqBuildBlockExcitationRFPulse, this class supports the use of
///    setFlipAngleArray() and setRunIndex() to prepare an array of RF pulse amplitudes.
///    Some restrictions apply:
///    - Only 11 and 121 binomial pulses are supported for this feature.
///    - A maximum number of MAX_SIZE_FLIPANGLEARRAY_SBBBINOMIAL = 30
///     of different flip angles is allowed.
///    - If clipping occurs, the array is not chopped, but downscaled.
///

class __IMP_EXP SeqBuildBlockBinomialPulses : public SeqBuildBlockExcitationRFPulse
{
  public:
    SeqBuildBlockBinomialPulses(
        SBBList*    pSBBList,                 ///<  SBB-List the SBB should be added to
        bool        bBipolar        = false,  ///<  selects "bipolar" slice select gradients or not
        double      dPhaseEvolution = 180.0,  ///<  phase evolution of the off-resonant component between two rf-pulses in deg
        long        lNoOfRFPulses   = 3,      ///<  number of bipolar rf-pulses
        const char* ptIdent         = nullptr ///<  ident of the SBB (note that using more than one instance of this SBB
                                              ///<  requires different idents for different SBBs!)
    );

    ///     The default constructor: creates monopolar 121
    /// binomial pulses with 180deg
    ///     phase evolution and a default ident.
    SeqBuildBlockBinomialPulses(SBBList* pSBBList = nullptr); ///<  SBB-List the SBB should be added to

    virtual ~SeqBuildBlockBinomialPulses();

    SeqBuildBlockBinomialPulses(const SeqBuildBlockBinomialPulses& right) = delete;
    SeqBuildBlockBinomialPulses& operator=(const SeqBuildBlockBinomialPulses& right) = delete;
    SeqBuildBlockBinomialPulses(SeqBuildBlockBinomialPulses&& right)                 = delete;
    SeqBuildBlockBinomialPulses& operator=(SeqBuildBlockBinomialPulses&& right) = delete;

    /// Overloaded base class method: Calculate RF info for given cuboid geometry
    bool calcSliceAdjSBBRFInfo(
        MrProt&                                     rMrProt,         //< IMP: points to the protocol structure.
        SeqLim&                                     rSeqLim,         //< IMP: points to the sequence limits structure.
        SeqExpo&                                    rSeqExpo,        //< IMP: points to the sequence exports structure
        const SLICEADJ::sCuboidGeometry&            sSliceAdjCuboid, //< IMP: cuboid geometry (single)
        std::vector<MrProtocolData::SeqExpoRFInfo>& vsRFInfo         //< EXP: RF info
        ) override;

    ///     Specifies the number of binomial RF-pulses.
    void setNoOfWERFPulses(long lValue); ///<  desired number of pulses

    ///     If true is passed the binomial pulses are bipolar,
    /// i.e. the slice select
    ///     gradients of the binomial pulses are not refocused
    /// between the pulses and
    ///     every other RF-pulse is sent over a negative gradient.
    void setUseBipolarGradients(bool bValue); ///<  user selection

    ///     Specifies the off-resonant component for the binomial
    /// pulses.
    void setChemicalShift(double dValue); ///<  chemical shift of the off-resonant component in ppm

    ///     Specifies the phase the off-resonant component gathers
    /// between two rf-pulses.
    void setPhaseEvolution(double dValue); ///<  desired phase evolution in degrees

    ///     Specifies the slice profile quality for the slice
    /// selective binomial sinc
    ///     pulses.
    virtual void setBandwidthTimeProduct(double dBTProduct); ///<  desired bandwidth-time-product for the sinc pulse

    /// Specifies the TBW adaption setting for the slice
    /// selective binomial sinc pulses.
    virtual void setTBWAdaption(bool bTBWAdaption); ///<  desired TBW adaption setting for the sinc pulse

    /// If true is passed, the SBB tries to realize the
    /// desired slice profile quality. If it can't be
    /// realized within the given time, the best possible
    /// value for the bandwidth-time-product is used.
    /// If a minimum slice quality profile is required
    /// (bandwidth-time-product > 0.0) and this quality
    /// cannot be achieved, the SBB returns an error later.
    ///
    /// If false is passed, the SBBs returns an error later,
    /// if the slice profile quality can not be achieved.
    virtual void setUsePossibleBandwidthTimeProduct(
        bool   bValue,              ///<  user selection
        double dMinBTProduct = -1.0 ///<  minimum bandwidth-time-product
    );

    ///     Specifies the duration of the non selective binomial
    /// rect rf-pulses.
    virtual void setNonSelectiveWERFDuration(long lDuration); ///<  duration of rf-pulses in us

    ///     Can be used to tell the SBB to use externally defined
    /// RF-pulses for the
    ///     binomial RF-pulse-train. Note that this function must
    /// be called for EVERY
    ///     pulse of the pulse train and that different instances
    /// of the RF-pulse-class
    ///     must be registered for different RF-Pulse-numbers.
    ///     Must be called before prep.
    virtual bool registerBinomialRFPulse(
        long       lRFPulseNumber, ///<  the number of the binomial RF-pulse
        IRF_PULSE* pRFPulse        ///<  the address of the RF-pulse that should be used
    );

    ///     Tells the SBB to use the SBB-internal RF-pulses for
    /// the binomial
    ///     RF-pulse-train again.
    void unregisterBinomialRFPulses();

    ///     Activates or deactivates the flip-angle adaption
    /// algorithm which is used to
    ///     get the effective flip angle of the off-resonant
    /// component as close as
    ///     possible to the desired flip angle.
    virtual void setUseFlipAngleAdaption(bool bValue); ///<  true to activate, false to deactivate

    ///     Specifies the slice thickness for the slice selective
    /// binomial rf-pulses.
    virtual void setThickness(double dThickness); ///<  the desired thickness in mm

    ///     Can be used to tell the SBB to use negative
    /// GS-gradients.
    ///     Must be called before prep.
    void setRequiredGSPolarity(double dPolarity); ///<  must be < 0 to force the SBB to use negative GS-gradients

    ///     Specifies the offset frequency for the binomial pulses
    /// with respect to the
    ///     system frequency.
    virtual void setUseWEOffsetFrequency(
        bool   bUseIt,          ///<  true, if offset frequency should be applied, else false
        double dFrequency = 0.0 ///<  offset-frequency in Hz
    );

    ///     Tries to scale the flip angles of the binomial pulses
    /// so that the effective
    ///     flip angle of the excited component becomes as close
    /// as possible to the
    ///     desired one.
    ///     Typically called by prepBinomialRFPulses only in
    /// context "normal".
    virtual void adaptFlipAngle(MrProt& rMrProt, SeqLim& rSeqLim);

    ///     Uses the simple vector model and rotation matrices
    /// for the binomial pulses
    ///     to predict the angles phi and theta in a spherical
    /// coordinate system of the
    ///     on-resonant and the off-resonant component.
    ///     prep must have been executed with success.
    virtual void checkAnglesWE(
        SeqLim& rSeqLim,
        double* pdPhiOn    = nullptr, ///<  pointer to double for the angle phi in degree for the on-resonant component
        double* pdPhiOff   = nullptr, ///<  pointer to double for the angle phi in degree for the off-resonant component
        double* pdThetaOn  = nullptr, ///<  pointer to double for the angle theta in degree for the on-resonant component
        double* pdThetaOff = nullptr  ///<  pointer to double for the angle theta in degree for the off-resonant component
    );

    ///     Checks, if the actual slice profile quality of the
    /// binomial sinc pulses is
    ///     equal to the one desired. Returns true if so, else
    /// false. An optional
    ///     parameter can be passed to get the used value of the
    /// bandwidth-time-product.
    ///     Note that this function always returns true, if
    /// m_bSincUsePossibleBTProduct
    ///     is true.
    ///     CalculateTiming must have been executed with success.
    virtual bool checkBandwidthTimeProduct(double* pdActualBTProduct = nullptr); ///<  pointer to the double this function should assign the actual used bandwidth time product to

    ///     Compare SeqBuildBlockExcitation::checkGradients.
    ///     prep must have been executed with success.
    bool checkGradients(MrProt& rMrProt, SeqLim& rSeqLim) override;

    /// getEffectiveFactorForSliceThickness%3ACAE4E9015A
    ///     The slice select gradient may restrict the minimum for
    /// the slice thickness
    ///     due to gradient amplitude or risetime limitations.
    /// This restriction may
    ///     cause that the slice thickness can only be decreased
    /// by a factor smaller
    ///     than the one given by the calculation limits.
    /// Therefore the sequence might
    ///     waste time, if it only takes into account the
    /// calculation limits for timing
    ///     calculations instead of calling this function.
    double getEffectiveFactorForSliceThickness(MrProt& rMrProt) override;

    ///     prints debug info
    void printDebugInformation() override;

    virtual void setProperties(bool, double, long, const char*);

    ///  Is called during run time, if an flip angle array is used.
    ///   For binomial pulses, this is not supported yet, so it will return false, if
    ///   water excitation is used.
    bool setRunIndex(unsigned int uiFlipAngleIndex) override;

    ///  This inline method setFlipAngleArray first calls the base class, which stores the passed parameters.
    ///  Then, it sets the member m_dFlipAngle_deg to the maximum value of the
    ///   passed flip angle array.
    bool setFlipAngleArray(double* padFlipAngleArray, unsigned int uiNmbrOfFlipAngles) override;

  protected:
    /// This method prepares the dedicated pulses for B1 control by setting the command setRelevantForB1Ctrl(true).
    /** It overrides the method from the base class and deals with single pulses and with binomial pulses. In case
    of just one pulse this pulse is used, in case of multiple pulses (binomial mode) one of the center pulses of
    the pulse train is used.*/
    bool prepareB1ControlLoop(SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    ///     Internal function calculating the flip angle required
    /// for a certain rf-pulse.
    virtual double calcFlipAngle(long lWERFPulseNo, MrProt& rMrProt, SeqLim& rSeqLim);

    ///     Compare
    /// SeqBuildBlockExcitation::updateGradientPerformance.
    void updateGradientPerformance(SEQ::Gradients eGradMode) override;

    ///     Internal function which determines the sort of pulse
    /// the SBB will prepare
    ///     and run depending on the protocol. Sets the following
    /// attributes:
    ///     - m_bUseBinomialPulses
    ///     - m_bExciteOffResComp
    ///     - m_bNonSelectiveExcitation
    virtual void setPrepAndRunMode(MrProt& rMrProt, SeqLim& rSeqLim);

    ///     Internal function which calculates the following
    /// members:
    ///     - m_lRFSpacing
    ///     - m_dPhaseEvolution
    bool setRFSpacingAndPhaseEvolution(MrProt&, SeqLim&);

    ///     Internal function which calculates the following
    /// members:
    ///     - m_dRequiredBinomialGSMoment
    ///     - m_bRequiredBinomialGSMomentSet
    virtual bool setRequiredBinomialGSMoment(MrProt&, SeqLim&);

    ///     Internal function which calculates the following
    /// members:
    ///     - m_lFixedBinomialRFDuration
    ///     - m_bFixedBinomialRFDurationSet
    virtual bool setFixedBinomialRFDuration(MrProt&, SeqLim&);

    ///     Internal function which calculates the timing of the
    /// following members:
    ///     - m_GSa
    ///     - m_GSb
    ///     The gradient timing is optimized for bipoloar gradients.
    bool calcBipolarGradients(MrProt& rMrProt, SeqLim& rSeqLim);

    ///     Internal function which calculates the timing of the
    /// following members:
    ///     - m_GSa
    ///     - m_GSb
    ///     The gradient timing is optimized for monopolar
    /// gradients.
    bool calcMonopolarGradients(MrProt& rMrProt, SeqLim& rSeqLim);

    ///     Internal function which calculates the following
    /// members:
    ///     - m_lRFDuration
    virtual bool setBinomialRFDuration(MrProt&, SeqLim& rSeqLim);

    bool prepExternalBinomialRFPulses(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    ///     Internal function which configures the following
    /// members:
    ///     - m_aRFns[*] or m_alRFsel[*]
    ///     - m_apRF[*]
    virtual bool configBinomialRFPulses(MrProt&, SeqLim& rSeqLim, SeqExpo&);

    ///     Internal function which prepares all RF-pulses the
    /// pointers m_apRF[*] point
    ///     to. All RF-configuration steps required for binomial
    /// pulses such as setting
    ///     of flip angle or initial phase are performed before
    /// preparation of the
    ///     pulses.
    ///     Calculates m_RFInfoPerRequest for binomial pulses.
    ///     NOTE THIS:
    ///     The center RF-Pulse is prepared first. If this
    /// RF-pulse should have been
    ///     clipped (because the requested flip angle is to high),
    /// than the flip angles
    ///     of the remaining RF-pulses are scaled so that the
    /// ratio of flip angles is
    ///     correct according to the binimial coefficients.
    bool prepBinomialRFPulses(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    ///     Internal function which prepares following gradients
    /// with the correct
    ///     amplitudes for bipolar gradients:
    ///     - m_GSa
    ///     - m_GSb
    bool prepBipolarGradients(MrProt&, SeqLim& rSeqLim, SeqExpo&);

    ///     Internal function which prepares following gradients
    /// with the correct
    ///     amplitudes for monopolar gradients:
    ///     - m_GSa
    ///     - m_GSb
    bool prepMonopolarGradients(MrProt&, SeqLim& rSeqLim, SeqExpo&);

    /// setExcitationBaseClassMembersBinomialPulses%3AC9F5B4039A
    ///     Sets all members of SeqBuildBlockExcitation described
    /// in documentation for
    ///     SeqBuildBlockExcitation::prep_insertEvents.
    bool setExcitationBaseClassMembersBinomialPulses(MrProt&, SeqLim& rSeqLim);

    ///     Sets start times of the following RTEs:
    ///     - m_GSa
    ///     - m_GSb
    ///     - m_apRF[*]
    bool setBinomialRTEStartTimes(MrProt&, SeqLim&);

    ///     Calls the following internal functions in order to
    /// calculate the timing of
    ///     the SBB and to prepare all RTEs:
    ///     - setPrepAndRunMode
    ///     - SeqBuildBlockExcitationRFPulse::prep_insertEvents
    ///       and returns, if m_bUseBinomialPulses is not true
    ///     - setRFSpacingAndPhaseEvolution
    ///     - prepExternalRFPulses, if
    /// m_bUseExternalBinomialRFPulses is true
    ///     - setRequiredBinomialGSMoment
    ///     - setFixedBinomialRFDuration
    ///     - calcBipolarGradients or calcMonopolarGradients
    ///     - sets m_dMaxPossibleBinomialGSMoment
    ///     - setBinomialRFDuration
    ///     - configBinomialRFPulses, if
    /// m_bUseExternalBinomialRFPulses is false
    ///     - prepBinomialRFPulses, if
    /// m_bUseExternalBinomialRFPulses is false
    ///     - prepBipolarGradients or prepMonopolarGradients
    ///     - setExcitationBaseClassMembersBinomialPulses
    ///     - setBinomialRTEStartTimes
    ///     Compare description of
    /// SeqBuildBlockExcitation::prep_insertEvents.
    bool prep_insertEvents(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    ///     Internal function which inserts all the gradients for
    /// a bipolar gradient
    ///     scheme into the timing table.
    bool insertBipolarGradients(MrProt&, SeqLim& rSeqLim, SeqExpo&, sSLICE_POS*);

    ///     Internal function which inserts all the gradients for
    /// a monopolar gradient
    ///     scheme into the timing table.
    bool insertMonopolarGradients(MrProt&, SeqLim& rSeqLim, SeqExpo&, sSLICE_POS*);

    ///     Internal functions which inserts all the binomial
    /// RF-pulses into the timing
    ///     table. The according frequency/phase-events are also
    /// created, prepared and
    ///     executed within this function.
    bool insertBinomialRFPulses(MrProt&, SeqLim& rSeqLim, SeqExpo&, sSLICE_POS* pSLC);

    ///     Calls the following internal functions in order to
    /// create the timing of the
    ///     SBB in an open event block:
    ///     - SeqBuildBlockExcitationRFPulse::run_insertEvents
    ///       and returns, if m_bUseBinomialPulses is not true
    ///     - insertBipolarGradients or insertMonopolarGradients
    ///     - insertBinomialRFPulses
    ///     Compare description of
    /// SeqBuildBlockExcitation::run_insertEvents.
    ///     prep must have been executed with success.
    ///     checkGradients must have been called with success.
    bool run_insertEvents(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC) override;

    bool   m_bBipolar{false};
    long   m_lNoOfWERFPulses{0};
    long   m_lCenterWERF{0};
    double m_dChemicalShift{CHEMICAL_SHIFT_FAT_PPM};
    long   m_lRFSpacing{0};
    long   m_lRFDuration{0};
    bool   m_bExciteOffResComp{false};
    bool   m_bNonSelectiveExcitation{false};
    bool   m_bAdaptFlipAngle{true};
    double m_dPhaseEvolutionWanted{180};
    double m_dPhaseEvolution{0.0};
    double m_dSincBTProductWanted{2.0};
    bool   m_bSincTBWAdaption{false};
    double m_dSincBTProductMin{-1.0};
    long   m_lRectRFDuration{300};
    bool   m_bUseExternalBinomialRFPulses{false};
    double m_dRequiredBinomialGSMoment{0.0};
    bool   m_bRequiredBinomialGSMomentSet{false};
    double m_dMaxPossibleBinomialGSMoment{0.0};
    long   m_lFixedBinomialRFDuration{0};
    bool   m_bFixedBinomialRFDurationSet{false};
    double m_dThickness{-1.0};
    double m_dRequiredGSPolarity{1.0};
    bool   m_bSincUsePossibleBTProduct{false};
    bool   m_bUseWEOffsetFrequency{false};
    double m_dWEOffsetFrequency{0.0};
    bool   m_bUseBinomialPulses{false};

    ///  Member that stores the Flip Angle Array for the 1st Binomial Pulse.
    ///  Note, that only Binomial pulses with up to 3 pulses are supported for the
    ///  Flip angle array mechanism.
    std::array<double, MAX_SIZE_FLIPANGLEARRAY_SBBBINOMIAL> m_adInternalFlipAngleArray0{};

    ///  Member that stores the Flip Angle Array for the 2nd Binomial Pulse.
    std::array<double, MAX_SIZE_FLIPANGLEARRAY_SBBBINOMIAL> m_adInternalFlipAngleArray1{};

    ///  Member that stores the Flip Angle Array for the 3rd Binomial Pulse.
    std::array<double, MAX_SIZE_FLIPANGLEARRAY_SBBBINOMIAL> m_adInternalFlipAngleArray2{};

    ///  Member that stores the correction phase for the Binomial Pulses when using Flip Angle Arrays.
    std::array<double, MAX_SIZE_FLIPANGLEARRAY_SBBBINOMIAL> m_adInternalCorrPhase{};

    ///  Member that stores last set RunIndex when using Flip Angle Arrays. Set by SeqBuildBlockBinomialPulses::setRunIndex()
    unsigned int m_uiCurrentRunIndex{0};

    IRF_PULSE**      m_apRF{nullptr};
    sGRAD_PULSE_TRAP m_GSa{"GS A"};
    sGRAD_PULSE_TRAP m_GSb{"GS B"};
    sRF_PULSE_RECT*  m_aRFns{nullptr};
    sRF_PULSE_SINC*  m_aRFsel{nullptr};

    virtual void init();

  private:
    ///  sFREQ_PHASE member that controls the NCO phase and frequency at the begin of each binomial rf pulse
    std::array<sFREQ_PHASE, MAX_NO_RF_PULSES> m_BinoRFSet{};

    ///  sFREQ_PHASE member that controls the NCO phase and frequency at the end of each binomial rf pulse
    std::array<sFREQ_PHASE, MAX_NO_RF_PULSES> m_BinoRFNeg{};

    friend class SBBTestHelperBinomialPulses; // SBB iTest Helper Class to Dump Class Member
};

// Class SeqBuildBlockBinomialPulses

//    Specifies the number of binomial RF-pulses.
inline void SeqBuildBlockBinomialPulses::setNoOfWERFPulses(long lValue)
{
    if (m_lNoOfWERFPulses != lValue)
    {
        resetTimingCalculated();

        m_lNoOfWERFPulses = std::min((long)MAX_NO_RF_PULSES, std::max(2L, lValue));

        for (long lI = 0; lI < MAX_NO_RF_PULSES; lI++)
        {
            m_aRFsel[lI].setTypeUndefined();
            m_aRFns[lI].setTypeUndefined();
        }

        m_lCenterWERF = m_lNoOfWERFPulses / 2;
    }
}

//    If true is passed the binomial pulses are bipolar,
// i.e. the slice select
//    gradients of the binomial pulses are not refocused
// between the pulses and
//    every other RF-pulse is sent over a negative gradient.
inline void SeqBuildBlockBinomialPulses::setUseBipolarGradients(bool bValue)
{
    if (m_bBipolar != bValue)
    {
        resetTimingCalculated();
        m_bBipolar = bValue;
    }
}

//    Specifies the off-resonant component for the binomial
// pulses.
inline void SeqBuildBlockBinomialPulses::setChemicalShift(double dValue)
{
    if (m_dChemicalShift != dValue)
    {
        resetTimingCalculated();
        m_dChemicalShift = dValue;
    }
}

//    Specifies the phase the off-resonant component gathers
// between two rf-pulses.
inline void SeqBuildBlockBinomialPulses::setPhaseEvolution(double dValue)
{
    if (m_dPhaseEvolutionWanted != dValue)
    {
        resetTimingCalculated();

        if (dValue < 10.0 || dValue > 360.0)
            m_dPhaseEvolutionWanted = 180.0;
        else
            m_dPhaseEvolutionWanted = dValue;
    }
}

//    Specifies the slice profile quality for the slice
// selective binomial sinc
//    pulses.
inline void SeqBuildBlockBinomialPulses::setBandwidthTimeProduct(double dBTProduct)
{
    if (m_dSincBTProductWanted != dBTProduct)
    {
        resetTimingCalculated();
        m_dSincBTProductWanted = std::max(0.5, dBTProduct);
    }
}

// Specifies the TBW adaption setting for the slice
// selective binomial sinc pulses.
inline void SeqBuildBlockBinomialPulses::setTBWAdaption(bool bTBWAdaption)
{
    if (m_bSincTBWAdaption != bTBWAdaption)
    {
        resetTimingCalculated();
        m_bSincTBWAdaption = bTBWAdaption;
    }
}

// setUsePossibleBandwidthTimeProduct%36B94C520251
//
// If true is passed, the SBB tries to realize the
// desired slice profile quality. If it can't be
// realized within the given time, the best possible
// value for the bandwidth-time-product is used.
// If a minimum slice quality profile is required
// (bandwidth-time-product > 0.0) and this quality
// cannot be achieved, the SBB returns an error later.
//
// If false is passed, the SBBs returns an error later,
// if the slice profile quality can not be achieved.
inline void SeqBuildBlockBinomialPulses::setUsePossibleBandwidthTimeProduct(bool bValue, double dMinBTProduct)
{
    if ((m_bSincUsePossibleBTProduct != bValue) || (m_dSincBTProductMin != dMinBTProduct))
    {
        resetTimingCalculated();
        m_bSincUsePossibleBTProduct = bValue;
        m_dSincBTProductMin         = dMinBTProduct;
    }
}

//    Specifies the duration of the non selective binomial
// rect rf-pulses.
inline void SeqBuildBlockBinomialPulses::setNonSelectiveWERFDuration(long lDuration)
{
    m_lRectRFDuration = std::max((long)MIN_NON_SELECTIVE_RF_DURATION, lDuration);
}

//    Tells the SBB to use the SBB-internal RF-pulses for
// the binomial
//    RF-pulse-train again.
inline void SeqBuildBlockBinomialPulses::unregisterBinomialRFPulses()
{
    m_bUseExternalBinomialRFPulses = false;
    resetTimingCalculated();
}

//    Activates or deactivates the flip-angle adaption
// algorithm which is used to
//    get the effective flip angle of the off-resonant
// component as close as
//    possible to the desired flip angle.
inline void SeqBuildBlockBinomialPulses::setUseFlipAngleAdaption(bool bValue)
{
    m_bAdaptFlipAngle = bValue;
}

//    Specifies the slice thickness for the slice selective
// binomial rf-pulses.
inline void SeqBuildBlockBinomialPulses::setThickness(double dThickness)
{
    if (m_dThickness != dThickness)
    {
        resetTimingCalculated();
        if (dThickness > 0.0)
            m_dThickness = dThickness;
    }
}

//    Can be used to tell the SBB to use negative
// GS-gradients.
//    Must be called before prep.
inline void SeqBuildBlockBinomialPulses::setRequiredGSPolarity(double dPolarity)
{
    if (dPolarity * m_dRequiredGSPolarity < 0.0)
    {
        if (dPolarity < 0)
            m_dRequiredGSPolarity = -1.0;
        else
            m_dRequiredGSPolarity = 1.0;
        resetPrepared();
    }
}

//    Specifies the offset frequency for the binomial pulses
// with respect to the
//    system frequency.
inline void SeqBuildBlockBinomialPulses::setUseWEOffsetFrequency(bool bUseIt, double dFrequency)
{
    if (bUseIt)
    {
        m_bUseWEOffsetFrequency = true;
        m_dWEOffsetFrequency    = dFrequency;
    }
    else
    {
        m_bUseWEOffsetFrequency = false;
        m_dWEOffsetFrequency    = 0.0;
    }
}

//    Checks, if the actual slice profile quality of the
// binomial sinc pulses is
//    equal to the one desired. Returns true if so, else
// false. An optional
//    parameter can be passed to get the used value of the
// bandwidth-time-product.
//    Note that this function always returns true, if
// m_bSincUsePossibleBTProduct
//    is true.
//    CalculateTiming must have been executed with success.
inline bool SeqBuildBlockBinomialPulses::checkBandwidthTimeProduct(double* pdActualBTProduct)
{
    double dBTProduct = m_dSincBTProductWanted;

    if (m_bSincUsePossibleBTProduct)
    {
        dBTProduct = std::min(m_dSincBTProductWanted, m_dMaxPossibleBinomialGSMoment * m_dGamma * m_dThickness / 1.0E6);
    }

    if (pdActualBTProduct != nullptr)
        *pdActualBTProduct = dBTProduct;

    if (fabs(dBTProduct - m_dSincBTProductWanted) / m_dSincBTProductWanted > 0.005)
    {
        setNLSStatus(MRI_SBB_SBB_ERROR, "SeqBuildBlockBinomialPulses::checkBandwidthTimeProduct", "actual BT-product not equal to wanted one");
        return false;
    }
    else
        return true;
}

/// This method setFlipAngleArray first calls the base class, which stores the passed parameters.
/// Then, it sets the member m_dFlipAngle_deg to the maximum value of the
///  passed flip angle array.

inline bool SeqBuildBlockBinomialPulses::setFlipAngleArray(double* padFlipAngleArray, unsigned int uiNmbrOfFlipAngles)
{
    if (!SeqBuildBlockExcitationRFPulse::setFlipAngleArray(padFlipAngleArray, uiNmbrOfFlipAngles))
    {
        setNLSStatus(MRI_SEQ_SEQU_ERROR, "SeqBuildBlockBinomialPulses::setFlipAngleArray", "SeqBuildBlockExcitationRFPulse::setFlipAngleArray failed");
        return false;
    }

    double dMaxRequestedFlipAngle = 0.0;

    for (unsigned int ui = 0; ui < uiNmbrOfFlipAngles; ui++)
    {
        dMaxRequestedFlipAngle = std::max(dMaxRequestedFlipAngle, m_padFlipAngleArray[ui]);
    }

    m_dFlipAngle_deg = dMaxRequestedFlipAngle;

    return true;
}

#endif
