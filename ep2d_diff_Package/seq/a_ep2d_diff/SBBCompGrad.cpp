//----------------------------------------------------------------------------------
// <copyright file="SBBCompGrad.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 1998-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#include "MrImaging/seq/a_ep2d_diff/SBBCompGrad.h"
#include <algorithm>

using namespace SEQ_NAMESPACE;

SeqBuildBlockCompGrad::SeqBuildBlockCompGrad(SBBList* pSBBList) : SeqBuildBlockSpoilGrad(pSBBList)
{
    setIdent("Comp Grad");

    // ----------------------------------------------------
    // Configure SBB-specific dynamic adjustment properties
    // ----------------------------------------------------
    seteSliceAdjOptimizationMode(SLICEADJ::HOLD_OPTIMIZATION); // Works on existing coherence (no RF) => hold current parameter optimization
    setsSliceAdjParametersRequestedBySBB(SLICEADJ::ADJNONE);   // Irrelevant for HOLD mode without RF => no need to consider anything
}

void SeqBuildBlockCompGrad::calcDiffCompGrad(
    double dMaxDiffAmplitude, long lDiffGradDuration, long lDiffGradSpacing, long lDiffToCompGrad)
{
    // calculate const. spoil moments
    double dMomentSumLimit = m_dMaxCompGradMonment; // limit spoil moment (this is a trade off: low limit is not so
                                                    // efficient in terms of balancing the gradients, a larger limit
                                                    // results in longer spoilers and thus eg. increases min.TR )
    long   lCompGradDuration = 0L;
    double dCompGradAmpl     = m_dMaxMagnitude;
    double dDiffAmp          = m_dCompensationFraction * dMaxDiffAmplitude;
    double dCompGradMoment   = 0.0;

    // MK: calculate spoiling moment on all axes to achieve perfect compensation at the time of the fat
    // sat pulse

    if (m_bConsiderCompensationDecay && m_dEddyCurrentTau > 0.0)
    {
        //---------------- MK: new way to calculate duration of compensation gradient, given that the amplitude is
        // fixed and known a priori --------------------------

        /*
            For calculation of compensation gradient: consider three equivalent gradient blocks:
                                                             __
                                                            /  \             ___
                                                           /    \           /   \
                                                          /comp. \         / FS  \  [ Fat sat pulse]
                       ____           ______________ ____/ grad.  \_______/ Pulse \__________________
            \  diff   /    \  diff   /
             \ grad1./      \ grad2./
              \_____/        \_____/



            Total time, amplitudes of above gradients:
            |--- total time: TD1 and TD2, Diff amplitude: D1 and D2  ---|  | TB, B  |

            Times between end of (equivalent) gradients and beginning of fat sat pulse:
                        |---------TD01--------------------------------------| (Diff1 Grad -> FS pulse)
                                        |---TD02----------------------------| (Diff2 Grad -> FS pulse)
                                                                     |-TB0- | (Comp Grad  -> FS pulse)

            At the time of the fat sat pulse, eddy currents of the three gradient blocks should
            compensate each other (tau: time constant of eddy currents):

            D1 (exp(-TD1/tau) - 1) exp(-TD01/tau) + D2 (exp(-TD2/tau) - 1) exp(-TD02/tau) = B
            (exp(-TB/tau) - 1) exp(-TB0/tau)

            Below, the above equation is solved:
            For TB to calculate the duration of the compensating gradient, assuming the amplitude B is
            fixed (using the maximum specified amplitude). This will be done for the axis (read/slice/phase) that has
            diffusion gradient.
        */
        long lCompGradDurationStartValue = 12000L; // start with initial value for the iterative calculation

        // first: calculate the duration and moment of the compensating gradient, based on the larger moment of
        // read and slice axes
        int nIter = 0; // count iterations; limit to 10 iterations
        while (abs(lCompGradDuration - lCompGradDurationStartValue) > GRAD_RASTER_TIME && nIter < 10)
        {
            if (nIter > 0)
                lCompGradDurationStartValue = lCompGradDuration;

            if (dCompGradAmpl > 0.0 && m_dEddyCurrentTau > 0.0) // avoid division by zero
            {
                double dLogArg = dDiffAmp / dCompGradAmpl * (exp(-lDiffGradDuration / m_dEddyCurrentTau) - 1)
                                     * (exp(-(lCompGradDurationStartValue + lDiffToCompGrad) / m_dEddyCurrentTau)
                                        * (1 + exp(-lDiffGradSpacing / m_dEddyCurrentTau)))
                                 + 1;

                // check for positive argument, make sure that ln(...) delivers negative result
                if (dLogArg > 0.0 && dLogArg < 1.0)
                {
                    lCompGradDuration = fSDSRoundUpGRT(-m_dEddyCurrentTau * log(dLogArg));
                }
                else
                {
                    lCompGradDuration = 0L;
                }
            }
            else
            {
                lCompGradDuration = 0L;
            }

            nIter++;

            if (lCompGradDuration == 0) // if the calculated duration is 0: stop calculation - either,there is
                                        // something wrong, or the compensating gradient is not needed
            {
                break;
            }
        }
    }
    else
    {
        double dGradMomentMaxComp = 2 * (dDiffAmp * lDiffGradDuration);
        lCompGradDuration         = fSDSRoundUpGRT(dGradMomentMaxComp / dCompGradAmpl);
    }

    // assign these values for SBB
    setCalcMode(SeqBuildBlockSpoilGrad::eSpoilMode::eFixedDurationAndMoment);

    // max Gradient Duration is actually the total time, so add the ramp down time; limit to maximum
    // specified duration
    if (lCompGradDuration + fSDSRoundUpGRT(m_dMinRiseTime * dCompGradAmpl) > m_lMaxGradientDuration)
    {
        lCompGradDuration = m_lMaxGradientDuration - fSDSRoundUpGRT(m_dMinRiseTime * dCompGradAmpl);
    }
    setAvailableTime(lCompGradDuration + fSDSRoundUpGRT(m_dMinRiseTime * dCompGradAmpl));

    // the corresponding moment is just the product of the duration and amplitude
    dCompGradMoment = lCompGradDuration * dCompGradAmpl;

    if (dCompGradMoment > dMomentSumLimit)
    {
        dCompGradMoment = dMomentSumLimit;
    }

    // the gradient is not valid if the calculated duration is less than raster time
    // then no preparation needed
    if (lCompGradDuration < GRAD_RASTER_TIME)
        m_bValid = false;
    else
        m_bValid = true;

    if (m_bValid)
    {
        setMomentsConstant(/*RO:*/ -dCompGradMoment, /*PE:*/ -dCompGradMoment, /*SS:*/ -dCompGradMoment);
    }
}

bool SeqBuildBlockCompGrad::prepGPALoad(double dAmplitudeX)
{
    // clear
    m_sBalanceCompGrad.Reset();

    // set typical frequency content of compensation module [Hz]
    // same to that for diffusion module
    m_sBalanceCompGrad.bSetFrequency(100);

    // ------------------------------------------------------------------------
    // prepare compensation gradient events for balance calculations.
    // consider GPABALANCE_X_AXIS only

    // add compensation gradients and store ID's for scaling purposes
    m_lIDEX          = m_sBalanceCompGrad.lAddGradient(0, dAmplitudeX, m_GPX->getRampUpTime(), m_GPX->getRampDownTime(), m_GPX->getFlatTopTime(), GPABALANCE_X_AXIS);

    // check error status
    if (m_sBalanceCompGrad.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}

bool SeqBuildBlockCompGrad::prepGPALoad(double dAmplitudeX, double dAmplitudeY, double dAmplitudeZ)
{
    // clear
    m_sBalanceCompGrad.Reset();

    // set typical frequency content of compensation module [Hz]
    // same to that for diffusion module
    m_sBalanceCompGrad.bSetFrequency(100);

    // ------------------------------------------------------------------------
    // prepare diffusion gradient events for balance calculations.
    // consider GPABALANCE_X_AXIS only

    // add compensation gradients and store ID's for scaling purposes
    m_lIDEX = m_sBalanceCompGrad.lAddGradient(0, dAmplitudeX, m_GPX->getRampUpTime(), m_GPX->getRampDownTime(), m_GPX->getFlatTopTime(), GPABALANCE_X_AXIS);
    m_lIDEY = m_sBalanceCompGrad.lAddGradient(0, dAmplitudeY, m_GPY->getRampUpTime(), m_GPY->getRampDownTime(), m_GPY->getFlatTopTime(), GPABALANCE_Y_AXIS);
    m_lIDEZ = m_sBalanceCompGrad.lAddGradient(0, dAmplitudeZ, m_GPZ->getRampUpTime(), m_GPZ->getRampDownTime(), m_GPZ->getFlatTopTime(), GPABALANCE_Z_AXIS);

    // check error status
    if (m_sBalanceCompGrad.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}

bool SeqBuildBlockCompGrad::scaleGPALoad(double dScaleX)
{
    const int iNumberOfEvents = 1;

    const std::array<long, iNumberOfEvents>   alID    = {m_lIDEX};
    const std::array<double, iNumberOfEvents> adScale = {dScaleX};

    // ------------------------------------------------------------------------
    // scale compensation gradient events on x-axis
    if (!m_sBalanceCompGrad.bScaleEvent(alID, adScale))
    {
        SEQ_TRACE_ERROR << "error scaling gradient events.";
        return false;
    }
    return true;
}

bool SeqBuildBlockCompGrad::scaleGPALoad(double dScaleX, double dScaleY, double dScaleZ)
{
    const int iNumberOfEvents = 3;

    const std::array<long, iNumberOfEvents>   alID    = {m_lIDEX, m_lIDEY, m_lIDEZ};
    const std::array<double, iNumberOfEvents> adScale = {dScaleX, dScaleY, dScaleZ};

    // ------------------------------------------------------------------------
    // scale compensation gradients events
    if (!m_sBalanceCompGrad.bScaleEvent(alID, adScale))
    {
        SEQ_TRACE_ERROR << "error scaling gradient events.";
        return false;
    }

    return true;
}