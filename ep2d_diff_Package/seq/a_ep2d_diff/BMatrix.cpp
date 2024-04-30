//-----------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2021. All Rights Reserved.
//-----------------------------------------------------------------------------

#include "MrImaging/seq/a_ep2d_diff/BMatrix.h"

// ---------------------------------------------------------------------------
// Definitions
// ---------------------------------------------------------------------------

// Assignment logical axes => matrix elements
#define BMATRIX_X_AXIS SEQ::AXIS_PHASE
#define BMATRIX_Y_AXIS SEQ::AXIS_READOUT
#define BMATRIX_Z_AXIS SEQ::AXIS_SLICE

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------
#include "MrImaging/seq/SeqDebug.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include <cmath>

#define PRINT0(mask, format)              \
    {                                     \
        if (m_lDebugLevel & mask)         \
        {                                 \
            SEQ_TRACE_ALWAYS.print(format); \
        }                                 \
    }

#define PRINT1(mask, format, arg1)              \
    {                                           \
        if (m_lDebugLevel & mask)               \
        {                                       \
            SEQ_TRACE_ALWAYS.print(format, arg1); \
        }                                       \
    }

#define PRINT2(mask, format, arg1, arg2)              \
    {                                                 \
        if (m_lDebugLevel & mask)                     \
        {                                             \
            SEQ_TRACE_ALWAYS.print(format, arg1, arg2); \
        }                                             \
    }

#define PRINT3(mask, format, arg1, arg2, arg3)              \
    {                                                       \
        if (m_lDebugLevel & mask)                           \
        {                                                   \
            SEQ_TRACE_ALWAYS.print(format, arg1, arg2, arg3); \
        }                                                   \
    }

#define PRINT4(mask, format, arg1, arg2, arg3, arg4)              \
    {                                                             \
        if (m_lDebugLevel & mask)                                 \
        {                                                         \
            SEQ_TRACE_ALWAYS.print(format, arg1, arg2, arg3, arg4); \
        }                                                         \
    }

#define PRINT5(mask, format, arg1, arg2, arg3, arg4, arg5)              \
    {                                                                   \
        if (m_lDebugLevel & mask)                                       \
        {                                                               \
            SEQ_TRACE_ALWAYS.print(format, arg1, arg2, arg3, arg4, arg5); \
        }                                                               \
    }

#define PRINT6(mask, format, arg1, arg2, arg3, arg4, arg5, arg6)              \
    {                                                                         \
        if (m_lDebugLevel & mask)                                             \
        {                                                                     \
            SEQ_TRACE_ALWAYS.print(format, arg1, arg2, arg3, arg4, arg5, arg6); \
        }                                                                     \
    }

#ifndef SEQ_NAMESPACE
    #error SEQ_NAMESPACE not defined
#endif
using namespace SEQ_NAMESPACE;


// ===========================================================================
///  Initialize member variables and event list.
// ===========================================================================
void BMatrix::Reset()
// ===========================================================================
{
    // Set all member variables to default values
    m_dBxx = 0.;
    m_dBxy = 0.;
    m_dBxz = 0.;
    m_dByy = 0.;
    m_dByz = 0.;
    m_dBzz = 0.;

    m_dMx  = 0.;
    m_dMy  = 0.;
    m_dMz  = 0.;

    m_bIsBMatrixValid = false;
    m_bIsBValueValid  = false;
    m_bIsSorted       = false;

    // Reset event list
    m_sEventMap.clear();

    // Reset error status
    m_lStatus = 0;

    m_sStatus = "Ok";
}

// ===========================================================================
///  Get status of module.
// ===========================================================================
long BMatrix::lGetStatus(bool bDumpMessage) const
// ===========================================================================
{
    if (bDumpMessage)
    {
        SEQ_TRACE_ALWAYS.print("Status %s", m_sStatus.c_str());
    }

    return m_lStatus;
}

// ===========================================================================
///  Calculate the b-Matrix and return the elements
// ===========================================================================
bool BMatrix::bCalcBMatrix(double &dBxx, double &dBxy, double &dBxz, double &dByy, double &dByz, double &dBzz)
// ===========================================================================
{
    // Calculate all b-matrix elements and all moments

    // Do we already have valid results?
    if (!m_bIsBMatrixValid)
    {
        if (!m_bIsBValueValid)
        {
            m_dBxx = dCalcB12(BMATRIX_X_AXIS, BMATRIX_X_AXIS, m_dMx, m_dMx);
            m_dByy = dCalcB12(BMATRIX_Y_AXIS, BMATRIX_Y_AXIS, m_dMy, m_dMy);
            m_dBzz = dCalcB12(BMATRIX_Z_AXIS, BMATRIX_Z_AXIS, m_dMz, m_dMz);
        }
        m_dBxy = dCalcB12(BMATRIX_X_AXIS, BMATRIX_Y_AXIS, m_dMx, m_dMy);
        m_dBxz = dCalcB12(BMATRIX_X_AXIS, BMATRIX_Z_AXIS, m_dMx, m_dMz);
        m_dByz = dCalcB12(BMATRIX_Y_AXIS, BMATRIX_Z_AXIS, m_dMy, m_dMz);

        m_bIsBMatrixValid = true;
        m_bIsBValueValid  = true;
    }

    dBxx = m_dBxx;
    dBxy = m_dBxy;
    dBxz = m_dBxz;
    dByy = m_dByy;
    dByz = m_dByz;
    dBzz = m_dBzz;

    PRINT3 (DEBUG_RESULT, "Moments  : M = (% 7.1f   % 7.1f   % 7.1f) mT/m * ms", m_dMx,  m_dMy,  m_dMz  );
    PRINT3 (DEBUG_RESULT, "                % 7.1f   % 7.1f   % 7.1f           ", m_dBxx, m_dBxy, m_dBxz );
    PRINT3 (DEBUG_RESULT, "B-matrix : b = (% 7.1f   % 7.1f   % 7.1f) s/mm2    ", m_dBxy, m_dByy, m_dByz );
    PRINT3 (DEBUG_RESULT, "                % 7.1f   % 7.1f   % 7.1f           ", m_dBxz, m_dByz, m_dBzz );

    // Check whether 0^{th} moment is almost zero on all axes
    if ( (fabs(m_dMx) > dCheckM) || (fabs(m_dMy) > dCheckM) || (fabs(m_dMz) > dCheckM) )
    {
        return false;
    }
    
    return true;
}


// ===========================================================================
///  Calculate the diagonal elements and return the scalar b-value
// ===========================================================================
bool BMatrix::bCalcBValue(double &dBValue)
// ===========================================================================
{
    // Calculate all b-matrix elements and all moments

    // Do we already have valid results?
    if (!m_bIsBValueValid)
    {
        m_dBxx = dCalcB12(BMATRIX_X_AXIS, BMATRIX_X_AXIS, m_dMx, m_dMx);
        m_dByy = dCalcB12(BMATRIX_Y_AXIS, BMATRIX_Y_AXIS, m_dMy, m_dMy);
        m_dBzz = dCalcB12(BMATRIX_Z_AXIS, BMATRIX_Z_AXIS, m_dMz, m_dMz);

        m_bIsBValueValid  = true;
    }

    dBValue = m_dBxx + m_dByy + m_dBzz;

    PRINT3 (DEBUG_RESULT, "Moments  : M = (% 7.1f   % 7.1f   % 7.1f) mT/m * ms", m_dMx,  m_dMy,  m_dMz  );
    PRINT3 (DEBUG_RESULT, "                % 7.1f   % 7.1f   % 7.1f           ", m_dBxx, m_dBxy, m_dBxz );
    PRINT3 (DEBUG_RESULT, "B-matrix : b = (% 7.1f   % 7.1f   % 7.1f) s/mm2    ", m_dBxy, m_dByy, m_dByz );
    PRINT3 (DEBUG_RESULT, "                % 7.1f   % 7.1f   % 7.1f           ", m_dBxz, m_dByz, m_dBzz );

    // Check whether 0^{th} moment is almost zero on all axes
    if ( (fabs(m_dMx) > dCheckM) || (fabs(m_dMy) > dCheckM) || (fabs(m_dMz) > dCheckM) )
    {
        return false;
    }
    
    return true;
}


// ===========================================================================
///  Calculate one element of the b-Matrix.
// ===========================================================================
double BMatrix::dCalcB12(SEQ::GradientAxis eAxis1, SEQ::GradientAxis eAxis2, double &dMoment1, double &dMoment2)
// ===========================================================================
{
    double      dSlope1 = 0.;
    double      dSlope2 = 0.;

    double      dGrad1  = 0.;
    double      dGrad2  = 0.;

    double      dM1     = 0.;
    double      dM2     = 0.;

    double      dB12    = 0.;

    double      dSign   = 1.;

    long        lLastTime    = m_sEventMap.begin()->lStartTime;

    double      dGamma2      = (m_dGamma * 1.e-15) * (m_dGamma * 1.e-15);

    bool        bStore       = false;       // Indicates store / restore conditions

    // Method might be time critical (e.g. within binary search) => avoid
    // access to INI-database within loops
    bool        bDebugInternal = m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/debug_BMatrix", false);

    // ===========================================================================
    // Sort list of events by time
    // ===========================================================================
    if ( !m_bIsSorted )
    {
        std::sort(m_sEventMap.begin(), m_sEventMap.end(), [](const BMatrixEvent& sEvent1, const BMatrixEvent& sEvent2) {
            
            if (sEvent1.lStartTime < sEvent2.lStartTime)
            {
                return true;
            }
            if (sEvent1.lStartTime > sEvent2.lStartTime)
            {
                return false;
            }

            // Same start time: event type determines order
            // 1. SEQ::AXIS_READOUT
            // 2. SEQ::AXIS_PHASE
            // 3. SEQ::AXIS_SLICE
            // 4. SEQ::AXIS_UNDEFINED
            // => RF events latest (required)
            // => Gradient event order doesn't matter (however, std::sort requires strict weak implementation)

            if (sEvent1.eAxis == sEvent2.eAxis)
            {
                return false;
            }

            switch (sEvent1.eAxis)
            {
                case SEQ::AXIS_READOUT:
                    return true;
                case SEQ::AXIS_PHASE:
                    if (sEvent2.eAxis == SEQ::AXIS_READOUT)
                    {
                        return false;
                    }
                    return true;
                case SEQ::AXIS_SLICE:
                    if (sEvent2.eAxis == SEQ::AXIS_UNDEFINED)
                    {
                        return true;
                    }
                    return false;
                case SEQ::AXIS_UNDEFINED:
                default:
                    return false;
            }
        });

        // Event list is sorted now
        m_bIsSorted = true;
    }

    // ===========================================================================
    // Run over the complete list of events
    // ===========================================================================
    for (auto& sEvent : m_sEventMap)
    {
        // Get start time of current event
        const long lCurrentTime = sEvent.lStartTime;

        // Calculate evolution from last event to current one
        if (lCurrentTime != lLastTime)
        {
            double dTimeStep  = static_cast<double>(lCurrentTime - lLastTime);

            if (!bStore)
            {
            // Calculate change of b-value
            // (unit: (mT/m)^2 us^3)
                dB12 += dM1 * dM2 * dTimeStep
                    + 1. / 2. * (dM1 * dGrad2 + dM2 * dGrad1)                           * std::pow(dTimeStep, 2.0)
                    + 1. / 6. * (dM1 * dSlope2 + dM2 * dSlope1 + 2. * dGrad1 * dGrad2)  * std::pow(dTimeStep, 3.0)
                    + 1. / 8. * (dGrad1 * dSlope2 + dGrad2 * dSlope1)                   * std::pow(dTimeStep, 4.0)
                    + 1. / 20. * (dSlope1 * dSlope2)                                    * std::pow(dTimeStep, 5.0);

            // Calculate change of 0th moment (time-integral of gradient)
            // (unit: mT/m us)
                dM1 += dGrad1 * dTimeStep + 1. / 2. * dSlope1 * std::pow(dTimeStep, 2.0);
                dM2 += dGrad2 * dTimeStep + 1. / 2. * dSlope2 * std::pow(dTimeStep, 2.0);
            }
            else
            {
                // We are in a 'store' state: gradients don't contribute to b-matrix or first moment

                // Calculate change of b-value
                // (unit: (mT/m)^2 us^3)
                dB12 += dM1 * dM2 * dTimeStep;
            }

            // Calculate change of gradient amplitude 
            // (unit: mT/m)
            dGrad1 += dSlope1 * dTimeStep;
            dGrad2 += dSlope2 * dTimeStep;

            if (bDebugInternal)
            {
                PRINT3 (DEBUG_INTERNAL, "%7li us: new gradients: G = (% 7.1f   % 7.1f) mT/m     ", lCurrentTime, dGrad1,     dGrad2    );
                PRINT3 (DEBUG_INTERNAL, "%7li us: new moments  : M = (% 7.1f   % 7.1f) mT/m * ms", lCurrentTime, dM1 / 1.e3, dM2 / 1.e3);
                PRINT2 (DEBUG_INTERNAL, "%7li us: new b-value  : b =  % 7.1f s/mm2              ", lCurrentTime, dB12 * dGamma2        );
            }
        }

        // Consider the effect of the current event
        //
        if (sEvent.eAxis != SEQ::AXIS_UNDEFINED)
        {
            // This is a gradient pulse
            if (sEvent.eAxis == eAxis1)
            {
                dSlope1 += dSign * sEvent.dSlope;
                
                if (bDebugInternal)
                {
                    PRINT3(
                        DEBUG_INTERNAL,
                        "%7li us: Event Gradient on 1st axis - change slope by % 7.2f mT/m/ms. New slope: % 7.2f "
                        "mT/m/ms",
                        lCurrentTime,
                        sEvent.dSlope * 1.e3,
                        dSlope1 * 1.e3);
                }
            }
            
            if (sEvent.eAxis == eAxis2)
            {
                dSlope2 += dSign * sEvent.dSlope;
                
                if (bDebugInternal)
                {
                    PRINT3(
                        DEBUG_INTERNAL,
                        "%7li us: Event Gradient on 2nd axis - change slope by % 7.2f mT/m/ms. New slope: % 7.2f "
                        "mT/m/ms",
                        lCurrentTime,
                        sEvent.dSlope * 1.e3,
                        dSlope2 * 1.e3);
                }
            }
        }
        else if (sEvent.eType == TXTYPE_REFOCUS)
        {
            // This is a refocusing pulse
            dSign *= -1.;

            dGrad1 *= -1.;
            dGrad2 *= -1.;

            dSlope1 *= -1.;
            dSlope2 *= -1.;

            if (bDebugInternal)
            {
                PRINT1 (DEBUG_INTERNAL, "%7li us: Event refocusing RFPulse - invert sign", lCurrentTime);
            }
        }
        else if (sEvent.eType == TXTYPE_EXCITATION)
        {
            // This is an excitation pulse
            dB12 = 0.;

            dM1  = 0.;
            dM2  = 0.;

            bStore = false;

            if (bDebugInternal)
            {
                PRINT1 (DEBUG_LANDMARK, "%7li us: Event excitation RFPulse - reset evolution", lCurrentTime);
            }
        }
        else if (sEvent.eType == TXTYPE_UNDEFINED)
        {
            // This is a store / restore pulse

            if (!bStore)
            {
                // Switch to 'store' state
                bStore = true;

                if (bDebugInternal)
                {
                    PRINT1 (DEBUG_LANDMARK, "%7li us: Event store RFPulse - ignore gradients", lCurrentTime);
                }
            }
            else
            {
                // Switch to 'restore' state
                bStore = false;

                // Restore pulse acts like refocusing 
                dSign *= -1.;
                
                dGrad1 *= -1.;
                dGrad2 *= -1.;
                
                dSlope1 *= -1.;
                dSlope2 *= -1.;

                if (bDebugInternal)
                {
                    PRINT1 (DEBUG_LANDMARK, "%7li us: Event restore RFPulse - consider gradients", lCurrentTime);
                }
            }

        }


        // Prepare for next event
        lLastTime = lCurrentTime;
    }

    // Export: convert gradient moments to mT/m ms
    dMoment1 = dM1 / 1.e3;
    dMoment2 = dM2 / 1.e3;

    return dB12 * dGamma2;  // Return b-matrix element in s/mm^2
}

// ===========================================================================
///  Dump the b-Matrix and 0^{th} moments.
// ===========================================================================
void BMatrix::Dump() const
// ===========================================================================
{
    SEQ_TRACE_ALWAYS.print("   ( % 7.1f )           ( % 7.1f   % 7.1f  % 7.1f )     ", 
        m_dMx, m_dBxx, m_dBxy, m_dBxz ) ;
    SEQ_TRACE_ALWAYS.print("m0=( % 7.1f ) mT/m ms b=( % 7.1f   % 7.1f  % 7.1f ) s/mm2  tr b=%7.1f s/mm2", 
        m_dMy, m_dBxy, m_dByy, m_dByz, m_dBxx + m_dByy + m_dBzz ) ;
    SEQ_TRACE_ALWAYS.print("   ( % 7.1f )           ( % 7.1f   % 7.1f  % 7.1f )     ", 
        m_dMz, m_dBxz, m_dByz, m_dBzz ) ;
}

// ===========================================================================
///  Add gradient event to list.
// ===========================================================================
bool BMatrix::bAddEvent (long lStartTime, const sGRAD_PULSE_TRAP &sGradient)
// ===========================================================================
{
    // Get axis from gradient
    SEQ::GradientAxis eSetAxis = sGradient.getAxis();

    return bAddEvent (lStartTime, sGradient, eSetAxis);
}

bool BMatrix::bAddEvent (long lStartTime, const sGRAD_PULSE_TRAP &sGradient, SEQ::GradientAxis eAxis)
// ===========================================================================
{
    return bAddGradient (lStartTime, sGradient.getAmplitude(), sGradient.getRampUpTime(), sGradient.getRampDownTime(), sGradient.getFlatTopTime(), eAxis);
}

// ===========================================================================
///  Add RF pulse event to list.
// ===========================================================================
bool BMatrix::bAddEvent (long lStartTime, const IRF_PULSE *sRFPulse)
// ===========================================================================
{
    eTXTYPE eSetTxType;

        if (sRFPulse->isTypeExcitation())
        {
            eSetTxType = TXTYPE_EXCITATION;
        }
        else if (sRFPulse->isTypeRefocussing())
        {
            eSetTxType = TXTYPE_REFOCUS;
        }
        else
        {
        // Undefined RF pulse type: error
        m_lStatus = 1;

        m_sStatus = "Error: undefined RF pulse type";

        return false;
    }

    return bAddRFPulse (lStartTime, eSetTxType);
}

bool BMatrix::bAddEvent (long lStartTime, eTXTYPE eTxType)
// ===========================================================================
{
    return bAddRFPulse (lStartTime, eTxType);
}

// ===========================================================================
///  Add RF pulse event to list.
// ===========================================================================
bool BMatrix::bAddRFPulse (long lStartTime, eTXTYPE eTxType)
// ===========================================================================
{
    if (!(eTxType == TXTYPE_EXCITATION || eTxType == TXTYPE_REFOCUS || eTxType == TXTYPE_UNDEFINED))
    {
        // Undefined RF pulse type: error
        m_lStatus = 1;

        m_sStatus = "Error: undefined RF pulse type";

        return false;
    }

    // Insert refocusing event
    m_sEventMap.emplace_back(eTxType, SEQ::AXIS_UNDEFINED, 0.0, lStartTime);

    // New event inserted: all results are invalid now
    m_bIsBMatrixValid = false;
    m_bIsBValueValid  = false;
    m_bIsSorted       = false;

    PRINT1(DEBUG_LANDMARK, "RF pulse event added at %ldus", lStartTime);

    return true;
}

// ===========================================================================
///  Add gradient event to list.
// ===========================================================================
bool BMatrix::bAddGradient (long lStartTime, double dAmplitude, long lRampUp, long lRampDown, long lFlatTop, SEQ::GradientAxis eAxis)
// ===========================================================================
{
    if ( eAxis == SEQ::AXIS_UNDEFINED )
    {
        // Undefined gradient axis: error
        m_lStatus = 1; 

        m_sStatus = "Error: invalid gradient axis";

        return false;
    }

    if ( lRampUp == 0 )
    {
        // Undefined ramp up time: error
        m_lStatus = 1;

        m_sStatus = "Error: invalid gradient ramp up time";

        return false;
    }

    if ( lRampDown == 0 )
    {
        // Undefined ramp down time: error
        m_lStatus = 1;

        m_sStatus = "Error: invalid gradient ramp down time";

        return false;
    }

    // Insert four gradient related events:
    // 1. start of ramp up
    // 2. end   of ramp up
    // 3. start of ramp down
    // 4. end   of ramp down
    m_sEventMap.emplace_back(
        TXTYPE_UNDEFINED,
        eAxis,
        +dAmplitude / static_cast<double>(lRampUp),
        lStartTime);

    m_sEventMap.emplace_back(
        TXTYPE_UNDEFINED,
        eAxis,
        -dAmplitude / static_cast<double>(lRampUp),
        lStartTime + lRampUp);

    m_sEventMap.emplace_back(
        TXTYPE_UNDEFINED,
        eAxis,
        -dAmplitude / static_cast<double>(lRampDown),
        lStartTime + lRampUp + lFlatTop);

    m_sEventMap.emplace_back(
        TXTYPE_UNDEFINED,
        eAxis,
        +dAmplitude / static_cast<double>(lRampDown),
        lStartTime + lRampUp + lFlatTop + lRampDown);

    // New event inserted: all results are invalid now
    m_bIsBMatrixValid = false;
    m_bIsBValueValid  = false;
    m_bIsSorted       = false;

    PRINT6 (DEBUG_LANDMARK, "Gradient event on axis %d added at %ldus: %7fmT/m, %ldus, %ldus, %ldus", eAxis, lStartTime, dAmplitude, lRampUp, lFlatTop, lRampDown);

    return true;
}

bool BMatrix::bAddExcitationPulse(long lStartTime)
{
    return bAddRFPulse(lStartTime, TXTYPE_EXCITATION);
}

bool BMatrix::bAddRefocussingPulse(long lStartTime)
{
    return bAddRFPulse(lStartTime, TXTYPE_REFOCUS);
}

bool BMatrix::bAddReStorePulse(long lStartTime)
{
    return bAddRFPulse(lStartTime, TXTYPE_UNDEFINED);
}
