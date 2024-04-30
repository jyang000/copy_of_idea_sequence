//----------------------------------------------------------------------------------
// <copyright file="SBBEPIReadOut.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 1998-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015-2020. All Rights Reserved. Confidential.
// </copyright>
// <description>
//   SeqBuildBlockEPIReadOut implements a Sequence Building Block to execute an EPI read out train
// </description>
//----------------------------------------------------------------------------------

#pragma once

#ifndef SBBEPIReadOut_h
#define SBBEPIReadOut_h

#include "MrImaging/libSBB/libSBBmsg.h"
#include "MrImaging/libSeqUtil/CAIPIRINHABlip.h"
#include "MrImaging/libSeqUtil/SMSProperties.h"
#include "MrImaging/libSeqUtil/SUForbiddenTR.h"
#include "MrImagingFW/libSBBFW/SeqBuildBlock.h"
#include "MrImagingFW/libSeqUtilFW/ReorderInfo.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"
#include "MrMeasSrv/SeqIF/libMES/SEQData.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sREADOUT.h"
#include "MrMeasSrv/SeqIF/libRT/sROT_MATRIX.h"
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrRXSpec.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"

#include <algorithm>
#include <array>
#include <cmath>

#define SBBEPIReadOut_GRAD_PERF_RO 0
#define SBBEPIReadOut_GRAD_PERF_BLIPS 1
#define SBBEPIReadOut_NO_OF_GRAD_PERF_GROUPS 2
#define SBBEPIReadOut_PLACE_MAXCOND 5 // max. number of conditions == length of the Place vector

#ifdef BUILD_libSBB
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"


// --------------------------------------------------------------------------
// Forward declaration
// --------------------------------------------------------------------------
class SBBTestHelperEPIReadout; // SBB iTest Helper Class to Dump Class Member

/*!
  \ingroup grp_libsbb
  \todo add documentation
  \brief Sequence Building Block to execute an EPI read out train

This class provides an EPI read out train which can be used for all
sequences incorporating EPI like k-space sampling like single shot EPI,
segmented EPI or TGSE. Trapezoidal gradient pulses are used for the RO
gradient as well as for the phase encoding blips. The gradients are
calculated provide the optimum (i.e. the shortest possible) inter-echo
spacing (see calculateTiming).
A graphical representation of the SBB can be found in
/src/MrImaging/libSBB/SBBEPIReadOut.pdf.

The order of K-space sampling in 2d and 3d phase encoding direction is
controlled by the data stored within a class derived from ReorderInfo (oder
ReorderInfo itself) which must be instantiated by the client. The client has
to make this instance known to the EPI read out.

NOTE:
To correctly calculate the timing for segmented k-space sampling the client
has to tell the SBB how many lines/partitions the phase encoding blips
should be able to cover during the measurement.

This SBB offers the following configuration possibilities for the EPI echo
train (a detailed explanation follows):
- use regridding in RO-direction for faster k-space sampling
- use EPI like or TGSE like phase correction
- support echo-shifting
- use a fixed echo spacing instead of the shortest
- ignore the forbidden echo spacing range
- add RO gradients to the EPI like phase correction scan
- arbitrary echo train length
- arbitrary fill times before and after the echo train
- blind RO-gradients during the imaging scan
- start RO train with negative RO-gradient
- arbitrary NCO-phase offset
- independently mark first and/or last ADC as relevant for measurement time

use regridding in RO-direction for faster k-space sampling:
The ADC could already be switched on before the RO gradient has finished its
ramp up. If this is done, the sampled k-space does not have a constant
increment between the sampled points. This effect must be compensated by the
image reconstruction ("regridding").
As a result the echo spacing becomes shorter compared to sampling data
during the flat top of the RO gradient only. Nevertheless, if the required
RO-gradient amplitude becomes very close to the maximum allowed amplitude
for high resolution measurements with high bandwidth per pixel, it may not be
possible to use regridding any more, the SBB switches to non-regridding
automatically.

use EPI like or TGSE like phase correction:
Depending on the sequence the SBB is used for the client may want to use
different phase correction models. The supported methods require different
settings of Mdh entries during run, so the SBB must know which phase
correction method is used. If m_bTGSELikePhaseCorr is set to true, TGSE like
phase correction is used, else EPI like phase correction is used.

EPI like phase correction:
Only two different phase correction vectors are used within image
reconstruction. One for ADCs acquired with positive RO gradient and one for
ADCs acquired with negative RO gradient. The phase correction scan itself
acquires three ADCs. The first and the third ADC are averaged together to
one phase correction vector. So the Cseg-counter of the Mdh can only have
two values: 0 for positive and 1 for negative RO gradient amplitude.
NOTE:
To acquire EPI like phase correction data, the SBB must be put into the
phase correction mode. If the read out should be executed for imaging the
SBB must be put into the imaging mode again. The appropriate set methods
take care that the number of RO gradients applied is correct. They set the
ReorderInfo class into the correct phase correction mode and the MDH_
PHASCOR-flag is set correctly in the Mdh.

TGSE like phase correction:
Each echo of a multi-echo sequence has its own phase correction vector
within image reconstruction. The different phase correction vectors are
numbered with the Cseg-counter of the ADC from 0 to echos-1. In this mode
the EPI read out INCREASES the CSeg-counter already stored in the Mdh of m_
ADC when the run-function is called for each ADC sent by the actual internal
ADC-counter. This allows a TGSE sequence to correctly set the Cseg-counter
for multiple calls of the EPI read out run-function for different spin
echoes.
NOTE:
Typically the loop structure of the sequence controls the acquisition of
phase correction scans in this mode. The ReorderInfo-class must be put into
the phase correction mode by the sequence itself. Also the MDH_PHASCOR-flag
must be set correctly by the sequence.

support echo-shifting:
Suppose you have 12 PE-lines and EPI-factor 3 for a segmented EPI sequence.
If you plot the time at which a specific line is measured against the line
number you get the left plot (you may need a nice font ;-)). As the signal
amplitudes of the echos decrease with T2* similar steps can be seen in the
k-space data in PE-direction leading to artifacts. Phase errors from line to
line show similar behaviour. Applying echo-shifting leads to the right plot.
Amplitude and phase in k-space now show a linear behaviour which reduces
image artifacts. The duration of the SBB is increased by activating
echo-shifting.

           time                                 time
            ^                                    ^
            |        xxxx                        |           x
            |                                    |          x
            |                                    |         x
            |                                    |        x
without ES: |    xxxx                   with ES: |       x
            |                                    |      x
            |                                    |     x
            |                                    |    x
            |xxxx                                |   x
            |                                    |  x
            |                                    | x
            |                                    |x
            +------------> line no.              +------------> line no.

Echo shifting is activated with the method setUseEchoShifting. Before
calling the run function setCounterInSegmentForEchoShifting must be called
to tell the SBB how far the echos should be shifted within the next
execution of the run function.

use a fixed echo spacing instead of the shortest:
For special applications one does not want the echo spacing to change, when
relevant protocol parameters have been changed. This can be prevented by
advising the SBB to use a fixed echo spacing. If it is not possible to
realize the measurement with the current set of parameters within the given
echo spacing, the SBB will notify this with an error.

ignore the forbidden echo spacing range:
The SBB looks up the values getAcousticResonanceFrequency() and getAcoustic
ResonanceBandwidth() from the GC-proxy and calculates a forbidden echo
spacing range. To prevent the system from making too much noise the timing
calculation takes care that the echo spacing is not within this forbidden
range. The SBB can be advised to ignore this forbidden echo spacing range
during timing calculation.

add RO gradients to the EPI like phase correction scan:
The default EPI like phase correction scan has 3 RO gradients and 3 ADCs.
One may want to increase the number of RO gradient pulses before and after
the acquisition of the 3 phase correction echos. This is possible.

arbitrary echo train length:
The echo train for the imaging scans is equal to the number of RO gradients
switched within the imaging echo train. It can be set to any value the user
wants.

arbitrary fill times before and after the echo train:
The SBB can insert arbitrary fill times into the sequence timing table
before ("TE-fill") and after ("TR-fill") the RO gradient train.

blind RO-gradients during the imaging scan:
ADCs can be blinded at the beginning and at the end of the imaging echo
train. Note that this does not affect the number of RO gradients switched.
It reduces the number of ADCs acquired! The SBB starts acquisition with the
same line of k-space. This line is shifted to a later echo.

start RO train with negative RO-gradient:
Before every call of the run function the SBB can be asked to start with a
negative or a positive RO gradient. Note that these flags can also be set
with arguments to the methods setNextExecutionIsImaging and setNextExecution
IsPhaseCorrection.

arbitrary NCO-phase offset:
The frame of reference for the read out train can be rotated by a constant
phase which is applied to all ADCs in the echo train.

independently mark first and/or last ADC as relevant for measurement time:
Per default the first ADC put into the timing table by the run function is
marked as relevant for measurement time, if m_ADC is marked as relevant for
measurement time. For some applications it may be necessary to independently
mark the first and/or the last ADC of the echo train. A special set method
must be used for this purpose.

ignore the forbidden echo spacing ranges:
The SBB looks up the forbidden resonance frequencies using an instance
of the Forbidden class (libUICtrl) and calculates the forbidden echo
spacing ranges. To prevent the system from making too much noise the timing
calculation takes care that the echo spacing is not within any forbidden
range. The SBB can be advised to ignore all forbidden echo spacing ranges

*/

class __IMP_EXP SeqBuildBlockEPIReadOut : public SeqBuildBlock
{
  public:
    /// Calculates forbidden echo spacing ranges using
    ///  an instance of the Forbidden class.

    // No public default constructor
    SeqBuildBlockEPIReadOut() = delete;

    SeqBuildBlockEPIReadOut(SBBList* pSBBList); ///< points to the SBBList object which can be instantiated from the SBBList class. If you don't use the SBBList, nullptr can be passed.

    virtual ~SeqBuildBlockEPIReadOut() = default;

    SeqBuildBlockEPIReadOut(const SeqBuildBlockEPIReadOut& right) = delete;
    SeqBuildBlockEPIReadOut& operator=(const SeqBuildBlockEPIReadOut& right) = delete;
    SeqBuildBlockEPIReadOut(SeqBuildBlockEPIReadOut&& right)                 = delete;
    SeqBuildBlockEPIReadOut& operator=(SeqBuildBlockEPIReadOut&& right) = delete;

    /// Sets m_bPerformPhaseCorr.
    /// If EPI like phase correction is used within the same
    ///  kernel as the imaging
    /// scans, then this flag should be set to true. In this
    ///  case the time needed
    /// for the phase correction scans is taken into account
    ///  within the method get
    /// DurationPerRequest.
    virtual void setPerformPhaseCorrection(bool getalBValue);

    /// For EPI like phase correction:
    /// - set m_bExecuteAsPhaseCorrScan to true
    /// - call m_pRI->enablePhaseCorrection()
    /// - sets bStartWithNegativeROGrad
    /// Next call of the run function will execute an EPI like
    ///  phase correction scan
    /// with 3 ADCs.
    virtual void setNextExecutionIsPhaseCorrection(bool bNegativeROGrad = true); ///< if true, next run function will start with a negative RO gradient

    /// For EPI like phase correction:
    /// - set m_bExecuteAsPhaseCorrScan to false
    /// - call m_pRI->disablePhaseCorrection()
    /// - sets bStartWithNegativeROGrad
    /// Next call of the run function will execute an EPI
    ///  imaging read out train.
    virtual void setNextExecutionIsImaging(bool bNegativeROGrad = false); ///< if true, next run function will start with a negative RO gradient

    /// Sets m_bStartWithNegativeROGrad.
    /// If set to true, next run function will start with a
    ///  negative RO gradient
    virtual void setStartWithNegativeROGrad(bool getalBValue);

    /// Can be used to specify an additional phase to the NCO
    ///  during the measurement.
    virtual void setAdditionalPhase(double dValue); ///< the additional phase in deg

    /// Sets m_bUseShortestEchoSpacing to true.
    /// Method calculateTiming will provide the shortest
    ///  possible echo spacing.
    virtual void setUseShortestEchoSpacing();

    /// Sets m_lFixedEchoSpacing to the specified value and
    ///  m_bUseShortestEcho
    /// Spacing to false.
    /// Only, if argument is greater zero.
    /// Method calculateTiming will use the specified value
    ///  for the echo spacing.
    virtual void setUseFixedEchoSpacing(long lValue); ///< desired echo spacing in us

    /// Sets m_bUseRegriddingForRO.
    virtual void setUseRegriddingForRO(bool getalBValue);

    /// Sets m_dMaxRegridROAmplFactor.
    /// This value specifies the factor by which the amplitude
    ///  of the RO gradient
    /// with regridding may be greater than the amplitude of
    ///  the RO gradient without
    /// regridding (i.e. sampling on flat top completely) for
    ///  the same coverage of
    /// RO k-space.
    /// Note that the frequency range of the hardware filter
    ///  applied to each ADC
    /// depends on the amplitude of the RO gradient without
    ///  regridding. This means,
    /// if you are sampling center of k-space with a higher RO
    ///  gradient amplitude
    /// due to regridding, you are in danger of loosing data
    ///  due to this hardware
    /// filter. For a factor of 1.3 the filter proved to be
    ///  broad enough to let all
    /// required data pass through.
    void setMaxRegridROAmplFactor(double dValue);

    /// Sets m_lEchoTrainLength, i.e. the number of created
    ///  gradient echos.
    virtual void setEchoTrainLength(long NoOfEchos);

    /// Sets m_pRI.
    /// Sampling of k-space in 2d and 3d phase encoding
    ///  direction is controlled by
    /// the data stored in a ReorderInfo-object. The SBB has
    ///  to have access to this
    /// object and therefore needs a pointer to it.
    virtual bool setPointerToReorderInfo(ReorderInfo* aPointer);

    /// Sets m_lExpMaxLinesForPEBlip.
    /// This value is important for the method calculateTiming
    ///  to estimate the
    /// gradient moment required for the blip in 2d phase
    ///  encoding direction.
    virtual void setExpectedMaxLinesForPEBlip(long lMaxNoOfLines);

    /// Sets m_lExpMaxPartitionsFor3DBlip.
    /// This value is important for the method calculateTiming
    ///  to estimate the
    /// gradient moment required for the blip in 3d phase
    ///  encoding direction.
    virtual void setExpectedMaxPartitionsFor3DBlip(long lMaxNoOfPartitions);

    ///    Sets m_bTGSELikePhaseCorr to true.
    void setUseTGSELikePhaseCorrection();

    void setExecuteExternalEPIPhaseCorrectionScan(bool bValue);

    ///    Sets m_bTGSELikePhaseCorr to false.
    void setUseEPILikePhaseCorrection();

    ///    Sets m_bFlagPCforRTFeedback.
    ///    (Only effective if TGSELikePhaseCorr is disabled)
    void setFlagPCforRTFeedback(bool bFlagPCforRTFeedback);

    /// Increases the echo train length for the EPI like phase
    ///  correction.
    /// Sets m_lAddROGradientsBeforeEPIPhaseCorrScans and
    ///  m_lAddROGradientsAfter
    /// EPIPhaseCorrScans.
    /// The number of ADCs for phase correction remains
    ///  unchanged!
    void setAddROGradientsToPhaseCorrScan(long lAddROGradientsBefore, long lAddROGradientsAfter);

    /// Can be used to specify a fill time which is inserted
    ///  into the sequence
    /// timing before the first RO gradient pulse.
    /// Sets m_lLocalTEFill.
    virtual void setLocalTEFill(long lFillTime); ///< desired fill time in us

    /// Can be used to specify a fill time which is inserted
    ///  into the sequence
    /// timing after the last RO gradient pulse.
    /// Sets m_lLocalTRFill.
    virtual void setLocalTRFill(long lFillTime); ///< desired fill time in us

    /// Per default the EPIReadOut will not apply any fill
    ///  time if m_bExecuteAsPhase
    /// CorrScan is true. The EPIReadOut can be forced to
    ///  apply m_lLocalTEFill and m_
    /// lLocalTRFill even if phase correction is active, if
    ///  true is passed to this
    /// method.
    void setApplyFillTimesIfPhaseCorrection(bool getalBValue); ///< the additional phase in deg

    /// Activates or deactivates echo shifting.
    /// Sets m_lCountersPerSegmentForEchoShifting.
    /// This number controls, if echo shifting is active and
    ///  is also used for the
    /// calculation of the time increment which is multiplied
    ///  with m_lCounterIn
    /// SegmentForEchoShifting during run. It defines the
    ///  maximum allowed number for
    /// m_lCounterInSegmentForEchoShifting during run.
    void setUseEchoShifting(bool getalBValue, long lCountersPerSegment = 0); ///< number of steps for echo shifting calculation

    /// Sets m_lCounterInSegmentForEchoShifting.
    /// Must be between 0 and
    ///  m_lCountersPerSegmentForEchoShifting-1.
    void setCounterInSegmentForEchoShifting(long lCounter);

    /// Sets m_bIgnoreForbiddenEchoSpacingRange.
    /// If set to true calculateTiming will ignore the
    ///  forbidden echo spacing ranges
    virtual void setIgnoreForbiddenEchoSpacingRange(bool getalBValue);

    /// This function calculates the timing of the echo train.
    /// After this function finished with success the echo
    ///  spacing is known.
    ///
    /// The ADC is prepared depending an the bandwidth per
    ///  pixel selected by the
    /// user.
    ///
    /// All gradient calculations depend on:
    /// - the larmor constant of
    ///    pMrProt->getsTXSPEC().getasNucleusInfo()[0].gettNucleus()
    ///  - pMrProt->getsGRADSPEC().getulMode()
    /// - the available gradient performance stored in the
    ///    array which can be
    ///   accessed with the  methods setMinRiseTimes and
    ///    setMaxMagnitudes
    ///   of the base class SeqBuildBlock
    ///
    /// The RO gradient calculation depends on
    ///  - pMrProt->getsRXSPEC().effDwellTime(...)
    ///  - pMrProt->kSpace().getlBaseResolution()
    /// - pMrProt->sliceSeries().front().readoutFOV()
    /// - m_bUseRegriddingForRO
    /// - m_dMaxRegridROAmplFactor, if m_bUseRegriddingForRO
    ///    is true
    /// - m_pCalcLimits->getFactorForPixelSizeRO(...)
    ///
    /// The calculation of the phase encoding blip depends on
    /// - pMrProt->sliceSeries().front().phaseFOV()
    /// - m_lExpMaxLinesForPEBlip
    /// - m_pCalcLimits->getFactorForPhaseFOV(...)
    ///  - pMrProt->sliceSeries().front().getdThickness()
    /// - m_lExpMaxPartitionsFor3DBlip
    /// - m_pCalcLimits->getFactorForSliceThickness(...)
    ///
    /// In addition to that the function takes care that
    /// - the echo spacing is not within any forbidden
    ///    range if m_bIgnoreForbiddenEchoSpacingRange
    ///   is false
    /// - the resulting echo spacing is m_lFixedEchoSpacing, if
    ///   m_bUseShortestEchoSpacing is true
    /// - the time between two ADCs is not below the value
    ///    read with
    ///   getMinDurationBetweenReadoutAndReadout(...) from the
    ///    SysProperties
    virtual bool calculateTiming(MrProt& rMrProt, SeqLim& rSeqLim);

    /// Does final preparation steps and sets the exports which
    ///  control the
    /// regridding algorithm of the image reconstruction
    ///  program.
    /// calculateTiming must have been executed with success.
    bool prepSBB(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    /// Checks, if a gradient of the SBB violates its
    ///  performance specifications.
    /// Returns false, if at least one does.
    virtual bool checkGradients(MrProt& rMrProt, SeqLim& rSeqLim);

    /// Returns the echo-spacing in us.
    virtual long getEchoSpacing();

    /// getForbiddenEchoSpacingRangeSize
    /// Returns number of forbidden echo spacing ranges
    unsigned int getForbiddenEchoSpacingRangeSize() const;

    /// getForbiddenEchoSpacingRangeMin%39B885F00351
    /// Returns minimum of the forbidden echo spacing range iIndex in
    ///  us.
    long getForbiddenEchoSpacingRangeMin(unsigned int iIndex = 0) const;

    /// Returns maximum of the forbidden echo spacing range iIndex in
    ///  us.
    long getForbiddenEchoSpacingRangeMax(unsigned int iIndex = 0) const;

    /// Returns the duration of the SBB in us when executed as
    ///  imaging read out.
    virtual long getDurationEPIReadOutPerRequest();

    /// Returns the number of RO-gradients applied for the
    ///  EPI-like phase correction
    /// scan.
    /// This number may be higher than the number of echos
    ///  acquired for phase
    /// correction.
    long getEchoTrainLengthPhaseCorrScan() const;

    /// Returns the duration of the SBB in us when executed as
    ///  EPI-like phase
    /// correction scan.
    virtual long getDurationPhaseCorrScanPerRequest() const;

    /// Returns getDurationEPIReadOutPerRequest() +
    ///  getDurationPhaseCorrScanPer
    /// Request().
    //  Note: Overloaded base class method. SBBEPIReadout does not employ m_lSBBDurationPerRequest_us
    long getSBBDurationPerRequest() override;

    /// Returns a pointer to the internal member m_ADC. This
    ///  pointer is necessary to
    /// provide the SBB with all Mdh-entries which are not set
    ///  by the SBB's run
    /// function itself.
    virtual sREADOUT* getReadOutAddress();

    /// Returns the 2d phase encoding moment necessary to
    ///  PRE-phase the
    /// magnetization before the measurement of the first ADC.
    virtual double getPEPrePhasingMoment(MrProt& rMrProt, SeqLim&);

    /// Returns the 2d phase encoding moment necessary to
    ///  RE-phase the magnetization
    /// after the measurement of the last ADC.
    virtual double getPERePhasingMoment(MrProt& rMrProt, SeqLim&);

    /// Returns the 3d phase encoding moment necessary to
    ///  PRE-phase the
    /// magnetization before the measurement of the first ADC.
    virtual double get3DPrePhasingMoment(MrProt& rMrProt, SeqLim&);

    /// Returns the 3d phase encoding moment necessary to
    ///  RE-phase the magnetization
    /// after the measurement of the last ADC.
    virtual double get3DRePhasingMoment(MrProt& rMrProt, SeqLim&);

    /// Returns the read out moment necessary to PRE-phase the
    ///  magnetization before
    /// the measurement of the first ADC.
    virtual double getROPrePhasingMoment(MrProt&, SeqLim&);

    /// Returns a read out gradient containing the data stored
    ///  in m_GRO.
    sGRAD_PULSE_RO getGRO() const;

    /// Executes the EPI read out train.
    bool runSBB(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, sSLICE_POS* pSLC) override;

    /// run Inner Event block.
    virtual bool runInnerEvents(long lvalue);

    /// Sets boolean flags, indicating whether the first/last
    ///  ADC within the next
    /// call to the run-function of the derived type must be
    ///  flagged as valid for
    /// measurement time.
    void setADCRelevant(bool bFirst, bool bLast);

    /// Can be used to blind imaging ADCs at the beginning and
    ///  at the end of the
    /// echo train. Note that the number of RO gradients
    ///  remains identical to m_lEcho
    /// TrainLength, the number of ADCs switched becomes
    ///  smaller than the number of
    /// RO gradients switched.
    void setBlindImagingADCs(long lBlindADCsHead, long lBlindADCsTail);

    /// Operation: setMinRiseTimeScaleFactor.
    /// Sets the value of the scaling factor
    /// m_dMinRiseTimeScaleFactor, which is used
    /// to modify the gradient slew rate.
    virtual void setMinRiseTimeScaleFactor(double dValue);

    /// Additional Public Declarations

    /// PLACE (Phase Labeling for Additional Coordinate Encoding) functions:
    ///      set mode for shifting (constant or interleaved) the PEPrePhase by d number of k-space lines
    virtual void setPlaceMode(long lValue);

    /// sets shift the PEPrePhase by d number of k-space lines
    virtual void setPlaceShift(long lValue);

    /// retrieves the current value of the PLACE vector
    virtual long getPlaceCurrent();

    /// updates the PlaceVector/m_lPlaceCurrent using m_lPlaceMode/m_lPlaceShift
    virtual void updatePlace();

    /// advances the PLACE iterator of the PLACE vector
    virtual long setPlaceNext();

    /// reverse phase encoding direction
    virtual void setReversePE(bool bValue);

    /// status of PE direction
    virtual bool getIsReversePE();

    //  Enable slice-specific single-band Maxwell correction
    virtual void setMaxwellCorrSB_Slice(bool bEnable);

    //  Enable Maxwell correction for Flow Comp gradient in PE direction
    virtual void setMaxwellCorrFlowcomp_PE(bool bEnable);

    //  Enable Maxwell correction for Flow Comp gradient in PE direction
    virtual void setdAddtionalMaxwellPhase(double dAdditionalPhase);

    virtual void setRunMode(SliceAccelRFRunMode eRunMode);

    virtual bool setMultiBandFactor(long lMultiBandFactor);

    virtual bool setFOVShiftFactor(long lFOVShiftFactor);

    /// get 3D delta moment (slice acceleration or 3D-encoding) for a given protocol
    virtual double get3DDeltaMoment(MrProt& rMrProt);

  protected:
    /// When echo shifting is used a certain fill time
    ///  depending on the current
    /// (line res. partition) counter in the segment has to
    ///  be inserted before and
    /// after the EPI read-out train. This function calculates
    ///  those two fill times.
    virtual bool getEchoShiftingData(long lCounterInSegment, long* plDelay, long* plFillEnd);

    /// Is called by calculateTiming only, if the gradient
    ///  mode stored in the
    /// protocol has changed.
    /// Updates the performance of each gradient of the SBB
    ///  depending on the user
    /// selection for the gradient mode.
    void updateGradientPerformance(SEQ::Gradients eGradMode) override;

    /// Called by calculateTiming to find out the minimum time
    ///  between two ADCs.
    long minADCDistance(MrProt& rMrProt);

    /// Internal function called by the run-function to get
    ///  the number of RO
    /// gradients that should be switched with the current
    ///  execution of the SBB.
    virtual long run_getNoOfROGradients();

    /// Internal function called by the run-function to get
    ///  the number of the ADC
    /// that should be switched together with the current
    ///  RO-gradient.
    /// NOTE: Due to members set with
    ///  setAddROGradientsToPhaseCorrScan or setBlind
    /// ImagingADCs the number of ADCs switched may become
    ///  smaller than the number
    /// of RO gradients switched.
    virtual long run_getADCCounter(long lROGradientNo); /// RO gradient counter

    /// Internal function called by the run-function to get
    ///  information, wether the
    /// current ADC is the last ADC switched within the echo
    ///  train or not.
    virtual bool run_getIsLastADC(long lADCCounter);

    /// Internal function called by the run-function to set
    ///  the Mdh-entries which
    /// change from ADC to ADC.
    virtual bool run_setADCMdh(long lADCCounter, MdhProxy* pMdh, bool bIsLastADC);

    /// Operation: setDefaultGradientPerformance
    /// Sets maximum gradient amplitude using the
    /// scanner's default values.
    /// Sets the minimum gradient rise time using the
    /// scanners's default values and the scaling
    /// factor m_dMinRiseTimeScaleFactor
    virtual void setDefaultGradientPerformance();

    /// Operation: getAllowedEchoSpacing
    /// - If lCurrentEchoSpacing is outside the forbidden ranges
    ///   OR forbidden ranges shall be ignored,
    ///   the method returns the input value.
    /// - If lCurrentEchoSpacing is inside a forbidden range
    ///   AND those ranges shall be considered,
    ///   the method returns the first allowed echo spacing
    ///   that is larger than the input value.
    virtual long getAllowedEchoSpacing(long lCurrentEchoSpacing);

    bool isRunMaxwellSyncEvent();
    void setRequestMaxwellSyncEvent(bool request);
    void createAndSendSyncObject(double dEffectiveMaxwellFrequency);

    /// Data Members for Class Attributes

    bool m_bPerformPhaseCorr{true};

    bool m_bExecuteAsPhaseCorrScan{false};

    bool m_bExecuteExternalEPIPhaseCorrectionScan{false};

    bool m_bStartWithNegativeROGrad{false};

    double m_dAdditionalPhase{0.0};

    long m_lEchoTrainLength{0};

    long m_lEchoSpacing{0};

    bool m_bUseShortestEchoSpacing{true};

    long m_lFixedEchoSpacing{0};

    bool m_bUseRegriddingForRO{true};

    double m_dMaxRegridROAmplFactor{1.4};

    bool m_bTGSELikePhaseCorr{false};

    /// Indicates whether MDH_RTFEEDBACK flag will get applied to phase correction scans
    bool m_bFlagPCforRTFeedback{false};

    long m_lAddROGradientsBeforeEPIPhaseCorrScans{0};

    long m_lAddROGradientsAfterEPIPhaseCorrScans{0};

    long m_lBlindImagingADCsHead{0};

    long m_lBlindImagingADCsTail{0};

    long m_lExpMaxLinesForPEBlip{1};

    long m_lExpMaxPartitionsFor3DBlip{0};

    long m_lCountersPerSegmentForEchoShifting{0};

    long m_lCounterInSegmentForEchoShifting{0};

    long m_lLocalTEFill{0};

    long m_lLocalTRFill{0};

    bool m_bApplyFillTimesIfPhaseCorrection{false};

    long m_lShiftADC{0};

    long m_lShiftBlip{0};

    bool m_bIgnoreForbiddenEchoSpacingRange{false};

    /// Class handling forbidden echo spacings
    SUForbidden m_ForbiddenEchoSpacings;

    long m_lPlaceShift{0}; /// PLACE (Phase Labeling for Additional Coordinate Encoding)
    long m_lPlaceMode{0};
    long m_lPlaceCurrent{0};
    long m_lPlaceCurrentMin{0};
    long m_lPlaceCurrentMax{0};
    long m_lPlaceVectorPosition{-99999};
    long m_lPlaceVectorSize{1};

    bool m_isRequestMaxwellSyncEvent{false};

    std::array<long, SBBEPIReadOut_PLACE_MAXCOND> m_alPlaceVector{};

    bool m_bIsReversePE{false};

    /// Run-time information for the run function of the
    ///  derived type. If m_b
    /// ADCRelevantValid is true, the boolean flag should be
    ///  passed to the first ADC
    /// within the actual kernel call.
    bool m_bIsFirstADCRelevant{false};

    /// Run-time information for the run function of the
    ///  derived type. If m_b
    /// ADCRelevantValid is true, the boolean flag should be
    ///  passed to the last ADC
    /// within the actual kernel call.
    bool m_bIsLastADCRelevant{false};

    /// Run-time information for the run function of the
    ///  derived type. If the
    /// boolean flag is true, the ADC relevant flags have been
    ///  set immediately
    /// before the run-function of the kernel is called. In
    ///  this case the run
    /// function of the derived type is responsible to pass the
    ///  flags m_bIs
    /// First(Last)ADCRelevant to the related sREADOUT objects
    ///  and to reset the
    /// valid flag before run gives control back to the caller.
    bool m_bADCRelevantValid{false};

    sGRAD_PULSE_RO m_GRO{"m_GRO"};
    sGRAD_PULSE_PE m_Blip{"m_Blip"};
    sREADOUT       m_ADC{"m_ADC"};
    ReorderInfo*   m_pRI{nullptr};

    /// Scaling factor used to modify slew rate specified in meas perm section
    double m_dMinRiseTimeScaleFactor{1.0};

    // Slice-specific single-band Maxwell correction
    bool m_bIsMaxwellCorrSB_Slice{false};

    // The maxwell phase introduced by the flow compensation gradient in phase direction. 
    double m_dAddtionalMaxwellPhase_PE{0.0};

    // Maxwell correction for concomitant field term FC_PE
    bool m_bIsMaxwellCorrFlowcomp_PE{false};

    SliceAccelRFRunMode m_eRunMode{SINGLE_BAND};
    //CAIPIRINHABlip      m_CAIPIRINHABlip;
    long                m_lMultiBandFactor{1};
    long                m_lFOVShiftFactor{1};

    // flag used to check if echo shifting performed in SBBEPIReadout
    bool m_bEchoShifting{false};

  private:
    //double m_dPrePhasingSliceMoment{0.0}; // This variable is changed via set function from the SBBEPIKernel class

    inline bool isCenterRegionOfKSpace(long lROGradientNo) const;

    friend class SBBTestHelperEPIReadout; // SBB iTest Helper Class to Dump Class Member
};

/// Class SeqBuildBlockEPIReadOut

/// Sets m_bPerformPhaseCorr.
/// If EPI like phase correction is used within the same
///  kernel as the imaging
/// scans, then this flag should be set to true. In this
///  case the time needed
/// for the phase correction scans is taken into account
///  within the method get
/// DurationPerRequest.
inline void SeqBuildBlockEPIReadOut::setPerformPhaseCorrection(bool getalBValue)
{
    if (m_bTGSELikePhaseCorr)
        return;

    m_bPerformPhaseCorr = getalBValue;
}

/// For EPI like phase correction:
/// - set m_bExecuteAsPhaseCorrScan to true
/// - call m_pRI->enablePhaseCorrection()
/// - sets bStartWithNegativeROGrad
/// Next call of the run function will execute an EPI like
///  phase correction scan
/// with 3 ADCs.
inline void SeqBuildBlockEPIReadOut::setNextExecutionIsPhaseCorrection(bool bNegativeROGrad)
{
    if (m_bTGSELikePhaseCorr)
        return;

    m_bExecuteAsPhaseCorrScan  = true;
    m_bStartWithNegativeROGrad = bNegativeROGrad;

    m_pRI->enablePhaseCorrection();
}

/// For EPI like phase correction:
/// - set m_bExecuteAsPhaseCorrScan to false
/// - call m_pRI->disablePhaseCorrection()
/// - sets bStartWithNegativeROGrad
/// Next call of the run function will execute an EPI
///  imaging read out train.
inline void SeqBuildBlockEPIReadOut::setNextExecutionIsImaging(bool bNegativeROGrad)
{
    m_bExecuteAsPhaseCorrScan  = false;
    m_bStartWithNegativeROGrad = bNegativeROGrad;

  if ( !m_bExecuteExternalEPIPhaseCorrectionScan )
  {
    m_pRI->disablePhaseCorrection();
  }
  else
  {
      m_pRI->enablePhaseCorrection();
  }
}

/// Sets m_bStartWithNegativeROGrad.
/// If set to true, next run function will start with a
///  negative RO gradient
inline void SeqBuildBlockEPIReadOut::setStartWithNegativeROGrad(bool getalBValue)
{
    m_bStartWithNegativeROGrad = getalBValue;
}

/// Sets m_bUseShortestEchoSpacing to true.
/// Method calculateTiming will provide the shortest
///  possible echo spacing.
inline void SeqBuildBlockEPIReadOut::setUseShortestEchoSpacing()
{
    m_bUseShortestEchoSpacing = true;
    m_lFixedEchoSpacing       = 0;
}

/// Sets m_lFixedEchoSpacing to the specified value and
///  m_bUseShortestEcho
/// Spacing to false.
/// Only, if argument is greater zero.
/// Method calculateTiming will use the specified value
///  for the echo spacing.
inline void SeqBuildBlockEPIReadOut::setUseFixedEchoSpacing(long lValue)
{
    if (lValue > 0)
    {
        m_bUseShortestEchoSpacing = false;
        m_lFixedEchoSpacing       = lValue;
    }
    else
    {
        m_bUseShortestEchoSpacing = true;
        m_lFixedEchoSpacing       = 0;
    }
}

/// Sets m_bUseRegriddingForRO.
inline void SeqBuildBlockEPIReadOut::setUseRegriddingForRO(bool getalBValue)
{
    m_bUseRegriddingForRO = getalBValue;
}

/// Sets m_dMaxRegridROAmplFactor.
/// This value specifies the factor by which the amplitude
///  of the RO gradient
/// with regridding may be greater than the amplitude of
///  the RO gradient without
/// regridding (i.e. sampling on flat top completely) for
///  the same coverage of
/// RO k-space.
/// Note that the frequency range of the hardware filter
///  applied to each ADC
/// depends on the amplitude of the RO gradient without
///  regridding. This means,
/// if you are sampling center of k-space with a higher RO
///  gradient amplitude
/// due to regridding, you are in danger of loosing data
///  due to this hardware
/// filter. For a factor of 1.3 the filter proved to be
///  broad enough to let all
/// required data pass through.

inline void SeqBuildBlockEPIReadOut::setMaxRegridROAmplFactor(double dValue)
{
    m_dMaxRegridROAmplFactor = std::max(1.0, dValue);
}

/// Sets m_lEchoTrainLength, i.e. the number of created
///  gradient echos.
inline void SeqBuildBlockEPIReadOut::setEchoTrainLength(long NoOfEchos)
{
    m_lEchoTrainLength = std::max(1L, NoOfEchos);
}

/// Sets m_pRI.
/// Sampling of k-space in 2d and 3d phase encoding
///  direction is controlled by
/// the data stored in a ReorderInfo-object. The SBB has
///  to have access to this
/// object and therefore needs a pointer to it.
inline bool SeqBuildBlockEPIReadOut::setPointerToReorderInfo(ReorderInfo* aPointer)
{
    if (aPointer == nullptr)
    {
        setNLSStatus(MRI_SBB_SBB_MISSING_REORDERINFO, "SeqBuildBlockEPIReadOut::setPointerToReorderInfo", "no reordering info");
        return false;
    }
    m_pRI = aPointer;
    return true;
}

/// Sets m_lExpMaxLinesForPEBlip.
/// This value is important for the method calculateTiming
///  to estimate the
/// gradient moment required for the blip in 2d phase
///  encoding direction.
inline void SeqBuildBlockEPIReadOut::setExpectedMaxLinesForPEBlip(long lMaxNoOfLines)
{
    m_lExpMaxLinesForPEBlip = lMaxNoOfLines;
}

/// Sets m_lExpMaxPartitionsFor3DBlip.
/// This value is important for the method calculateTiming
///  to estimate the
/// gradient moment required for the blip in 3d phase
///  encoding direction.
inline void SeqBuildBlockEPIReadOut::setExpectedMaxPartitionsFor3DBlip(long lMaxNoOfPartitions)
{
    m_lExpMaxPartitionsFor3DBlip = lMaxNoOfPartitions;
}

///    Sets m_bTGSELikePhaseCorr to true.
inline void SeqBuildBlockEPIReadOut::setUseTGSELikePhaseCorrection()
{
    m_bTGSELikePhaseCorr      = true;
    m_bPerformPhaseCorr       = false;
    m_bExecuteAsPhaseCorrScan = false;
}

inline void SeqBuildBlockEPIReadOut::setExecuteExternalEPIPhaseCorrectionScan(bool bValue)
{
    m_bExecuteExternalEPIPhaseCorrectionScan = bValue;
}

inline void SeqBuildBlockEPIReadOut::setFlagPCforRTFeedback(bool bFlagPCforRTFeedback)
{
    m_bFlagPCforRTFeedback = bFlagPCforRTFeedback;
}

///    Sets m_bTGSELikePhaseCorr to false.
inline void SeqBuildBlockEPIReadOut::setUseEPILikePhaseCorrection()
{
    m_bTGSELikePhaseCorr = false;
}

/// Increases the echo train length for the EPI like phase
///  correction.
/// Sets m_lAddROGradientsBeforeEPIPhaseCorrScans and
///  m_lAddROGradientsAfter
/// EPIPhaseCorrScans.
/// The number of ADCs for phase correction remains
///  unchanged!
inline void SeqBuildBlockEPIReadOut::setAddROGradientsToPhaseCorrScan(long lAddROGradientsBefore, long lAddROGradientsAfter)
{
    m_lAddROGradientsBeforeEPIPhaseCorrScans = std::max(0L, lAddROGradientsBefore);
    m_lAddROGradientsAfterEPIPhaseCorrScans  = std::max(0L, lAddROGradientsAfter);
}

/// Per default the EPIReadOut will not apply any fill
///  time if m_bExecuteAsPhase
/// CorrScan is true. The EPIReadOut can be forced to
///  apply m_lLocalTEFill and m_
/// lLocalTRFill even if phase correction is active, if
///  true is passed to this
/// method.
inline void SeqBuildBlockEPIReadOut::setApplyFillTimesIfPhaseCorrection(bool getalBValue)
{
    m_bApplyFillTimesIfPhaseCorrection = getalBValue;
}

/// Activates or deactivates echo shifting.
/// Sets m_lCountersPerSegmentForEchoShifting.
/// This number controls, if echo shifting is active and
///  is also used for the
/// calculation of the time increment which is multiplied
///  with m_lCounterIn
/// SegmentForEchoShifting during run. It defines the
///  maximum allowed number for
/// m_lCounterInSegmentForEchoShifting during run.
inline void SeqBuildBlockEPIReadOut::setUseEchoShifting(bool getalBValue, long lCountersPerSegment)
{
    m_bEchoShifting = getalBValue;
    if (m_bEchoShifting == false || lCountersPerSegment < 2)
        m_lCountersPerSegmentForEchoShifting = 0;
    else
        m_lCountersPerSegmentForEchoShifting = lCountersPerSegment;
}

/// Sets m_lCounterInSegmentForEchoShifting.
/// Must be between 0 and
///  m_lCountersPerSegmentForEchoShifting-1.
inline void SeqBuildBlockEPIReadOut::setCounterInSegmentForEchoShifting(long lCounter)
{
    m_lCounterInSegmentForEchoShifting = std::max(0L, std::min(m_lCountersPerSegmentForEchoShifting, lCounter));
}

/// Sets m_bIgnoreForbiddenEchoSpacingRange.
/// If set to true calculateTiming will ignore the
///  forbidden ranges specified by m_ForbiddenEchoSpacings
inline void SeqBuildBlockEPIReadOut::setIgnoreForbiddenEchoSpacingRange(bool bValue)
{
    m_bIgnoreForbiddenEchoSpacingRange = bValue;
}

/// Called by calculateTiming to find out the minimum time
///  between two ADCs.
inline long SeqBuildBlockEPIReadOut::minADCDistance(MrProt& rMrProt)
{
    MrRXSpec RxSpecWrapper(rMrProt.getsRXSPEC());
    return std::lround(SysProperties::getMinDurationBetweenReadoutAndReadout(RxSpecWrapper.realDwellTime()[0] / 1000.0));
}

/// Returns the echo-spacing in us.
inline long SeqBuildBlockEPIReadOut::getEchoSpacing()
{
    return m_lEchoSpacing;
}

/// Returns number of forbidden echo spacing ranges
inline unsigned int SeqBuildBlockEPIReadOut::getForbiddenEchoSpacingRangeSize() const
{
    return m_ForbiddenEchoSpacings.getSize();
}

/// Returns minimum of the forbidden echo spacing range iIndex in
///  us.
inline long SeqBuildBlockEPIReadOut::getForbiddenEchoSpacingRangeMin(unsigned int iIndex) const
{
    if (iIndex < m_ForbiddenEchoSpacings.getSize())
    {
        /// convert from milliseconds to microseconds (factor 1000),
        /// convert from full period  to echo spacing (factor 1/2),
        /// round to gradient raster time
        return fSDSRoundDownGRT(0.5 * 1000. * m_ForbiddenEchoSpacings.getInterval()[iIndex].getMin());
    }

    return 0;
}

/// Returns maximum of the forbidden echo spacing range iIndex in
///  us.
inline long SeqBuildBlockEPIReadOut::getForbiddenEchoSpacingRangeMax(unsigned int iIndex) const
{
    if (iIndex < m_ForbiddenEchoSpacings.getSize())
    {
        /// convert from milliseconds to microseconds (factor 1000),
        /// convert from full period  to echo spacing (factor 1/2),
        /// round to gradient raster time
        return fSDSRoundUpGRT(0.5 * 1000. * m_ForbiddenEchoSpacings.getInterval()[iIndex].getMax());
    }

    return 0;
}

/// Returns the number of RO-gradients applied for the
///  EPI-like phase correction
/// scan.
/// This number may be higher than the number of echos
///  acquired for phase
/// correction.
inline long SeqBuildBlockEPIReadOut::getEchoTrainLengthPhaseCorrScan() const
{
    return 3 + m_lAddROGradientsBeforeEPIPhaseCorrScans + m_lAddROGradientsAfterEPIPhaseCorrScans;
}

/// Returns the duration of the SBB in us when executed as
///  EPI-like phase
/// correction scan.
inline long SeqBuildBlockEPIReadOut::getDurationPhaseCorrScanPerRequest() const
{
    return getEchoTrainLengthPhaseCorrScan() * m_lEchoSpacing * m_bPerformPhaseCorr;
}

/// Returns getDurationEPIReadOutPerRequest() +
///  getDurationPhaseCorrScanPer
/// Request().
inline long SeqBuildBlockEPIReadOut::getSBBDurationPerRequest()
{
    if (!isPrepared())
    {
        return 0;
    }

    if (m_bExecuteExternalEPIPhaseCorrectionScan)
        return getDurationEPIReadOutPerRequest();

    return getDurationEPIReadOutPerRequest() + getDurationPhaseCorrScanPerRequest();
}

/// Returns a pointer to the internal member m_ADC. This
///  pointer is necessary to
/// provide the SBB with all Mdh-entries which are not set
///  by the SBB's run
/// function itself.
inline sREADOUT* SeqBuildBlockEPIReadOut::getReadOutAddress()
{
    return &m_ADC;
}

/// Returns a read out gradient containing the data stored
///  in m_GRO.
inline sGRAD_PULSE_RO SeqBuildBlockEPIReadOut::getGRO() const
{
    return m_GRO;
}

/// Sets boolean flags, indicating whether the first/last
///  ADC within the next
/// call to the run-function of the derived type must be
///  flagged as valid for
/// measurement time.
inline void SeqBuildBlockEPIReadOut::setADCRelevant(bool bFirst, bool bLast)
{
    m_bIsFirstADCRelevant = bFirst;
    m_bIsLastADCRelevant  = bLast;
    m_bADCRelevantValid   = true;
}

/// Can be used to blind imaging ADCs at the beginning and
///  at the end of the
/// echo train. Note that the number of RO gradients
///  remains identical to m_lEcho
/// TrainLength, the number of ADCs switched becomes
///  smaller than the number of
/// RO gradients switched.
inline void SeqBuildBlockEPIReadOut::setBlindImagingADCs(long lBlindADCsHead, long lBlindADCsTail)
{
    if (lBlindADCsHead < 0)
        lBlindADCsHead = 0;
    if (lBlindADCsTail < 0)
        lBlindADCsTail = 0;

    if (lBlindADCsHead + lBlindADCsTail > m_lEchoTrainLength)
    {
        m_lBlindImagingADCsHead = m_lEchoTrainLength;
        m_lBlindImagingADCsTail = 0;
    }
    else
    {
        m_lBlindImagingADCsHead = lBlindADCsHead;
        m_lBlindImagingADCsTail = lBlindADCsTail;
    }
}

/// retrieves the current value of the PLACE vector
inline long SeqBuildBlockEPIReadOut::getPlaceCurrent()
{
    return (m_lPlaceCurrent);
}

/// retrieves the current value of Neg phase encoding
inline bool SeqBuildBlockEPIReadOut::getIsReversePE()
{
    return (m_bIsReversePE);
}

/// set single-band or multi-band mode for slice accelerated acquisition
inline void SeqBuildBlockEPIReadOut::setRunMode(SliceAccelRFRunMode eRunMode)
{
    m_eRunMode = eRunMode;
}

/// set multi-band factor for slice accelerated acquisition
inline bool SeqBuildBlockEPIReadOut::setMultiBandFactor(long lMultiBandFactor)
{
    if (lMultiBandFactor < 1 || lMultiBandFactor > SMSProperties::MAX_MULTIBAND_FACTOR) // To avoid passing incorrect value
    {
        SEQ_TRACE_ERROR.print("SeqBuildBlockEPIReadOut::setMultiBandFactor: ERROR - Multiband factor %li out of range [1, %li]. Setting to 1.", lMultiBandFactor, SMSProperties::MAX_MULTIBAND_FACTOR);
        m_lMultiBandFactor = 1;
        return false;
    }
    if (m_lMultiBandFactor == lMultiBandFactor)
        return true;

    m_lMultiBandFactor = lMultiBandFactor;
    resetPrepared();
    return true;
}

/// set FOV shift factor for slice accelerated acquisition
inline bool SeqBuildBlockEPIReadOut::setFOVShiftFactor(long lFOVShiftFactor)
{
    if ((lFOVShiftFactor < 1) || (lFOVShiftFactor > SMSProperties::MAX_FOV_SHIFT_FACTOR))
    {
        SEQ_TRACE_ERROR.print("SeqBuildBlockEPIReadOut::setFOVShiftFactor: ERROR - FOV shift factor %li out of range [1, %li]. Setting to 1.", lFOVShiftFactor, SMSProperties::MAX_FOV_SHIFT_FACTOR);
        m_lFOVShiftFactor = 1;
        return false;
    }
    if (m_lFOVShiftFactor == lFOVShiftFactor)
        return true;

    m_lFOVShiftFactor = lFOVShiftFactor;
    resetPrepared();
    return true;
}

inline bool SeqBuildBlockEPIReadOut::isCenterRegionOfKSpace(long lEchoNumber) const
{
    // ReorderInfo knows the segment (i.e. echo) index which acquires
    // k-space center. Note that for multi-shot acquisitions, only
    // one shot samples the actual k-space center line.
    const long lCenterEchoIndex = std::max(m_pRI->getSegmentsBeforeKSpaceCenter(), 1L);

    // Check whether the given echo number refers to the central k-space region (+/-1)
    if (labs(lCenterEchoIndex - lEchoNumber) <= 1)
    {
        return true;
    }

    return false;
}

inline void SeqBuildBlockEPIReadOut::setMaxwellCorrSB_Slice(bool bEnable)
{
    // Single-band correction gets applied on the fly - does not require any preparation
    m_bIsMaxwellCorrSB_Slice = bEnable;
}

inline void SeqBuildBlockEPIReadOut::setMaxwellCorrFlowcomp_PE(bool bEnable)
{
    // Maxwell correction for flow compensated gradient in PE direction gets applied on the fly - does not require any preparation
    m_bIsMaxwellCorrFlowcomp_PE = bEnable;
}

inline void SeqBuildBlockEPIReadOut::setdAddtionalMaxwellPhase(double dAdditionalPhase)
{
    m_dAddtionalMaxwellPhase_PE = dAdditionalPhase;
}
#endif
