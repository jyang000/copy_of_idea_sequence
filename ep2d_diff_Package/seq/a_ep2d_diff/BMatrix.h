//-----------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2020. All Rights Reserved.
//-----------------------------------------------------------------------------

// double include protection:
#pragma once

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------
#include <vector>                                               /**< STL vector                     */
#include "MrProtSrv/Domain/CoreNative/GradSpecDefines.h"               /**< SEQ::  definitions             */
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"          /**< RF_PULSE and TXTYPE definition */
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"        /**< GRAD_PULSE definition          */
#include "MrMeasSrv/MeasNuclei/Base/MeasNucleus.h"    /**< MeasNucleus                    */
#include "MrMeasSrv/MeasUtils/MeasMath.h"
#include "MrImaging/seq/a_ep2d_diff/SequenceDebugSettings.h"
#ifndef SEQ_NAMESPACE
#pragma message ("NOTE: SEQ_NAMESPACE not defined");
#endif
namespace SEQ_NAMESPACE
{
    // ---------------------------------------------------------------------------
    // Constants
    // ---------------------------------------------------------------------------
    const double dCheckM              = 1.e-6;                /**< Final moment is checked against this value [mT/m ms] */
    const double dDefaultGamma        = 42.5756e6*2.*M_PI;    /**< Default gyromagnetic ratio */

    // ---------------------------------------------------------------------------
    // Type definitions
    // ---------------------------------------------------------------------------

    // ===========================================================================
    /*!
    \class BMatrix

    \brief This class implements a generic b-matrix calculation based on a series
    of imaging events (trapezoidal gradients, refocusing pulses).

    The b-matrix describes the signal attenuation effect of a series of
    diffusion encoding gradients for a given diffusion tensor. It is defined
    by the following equation:

    ln A(b)/A(0) = - sum_{i=1}^3 sum_{j=1}^3 b_{ij} D_{ij}

    where D is the diffusion tensor and A describes the signal intensity
    with (b != 0) and without (b = 0) diffusion encoding.

    Calculation takes place based on a list of diffusion encoding events
    (trapezoidal gradient pulses and RF pulses), which might even overlap.
    Events are collected in an ordered list (vector based), which is the
    basis of the b-value calculation.

    All event timings must be provided in us.

    Coordinate system is determined by the gradient events: no
    coordinate transformations take place. If gradient events are
    provided in XYZ-coordinates, the b-matrix is calculated in this
    reference frame.

    INI-configurations (since this is merely a helper class, the INI-file of
    the calling instance will be used):

    - For each section, traces can be enabled by inserting the line
    TRACES=CALL,INPUT,RETURN,RESULT,LANDMARK,INTERNAL
    (or a selection of the listed severities)



    \author Thorsten.Feiweier@siemens.com

    */
    // ===========================================================================
    class BMatrix
    {
        // ------------------------------------
        // Public methods
        // ------------------------------------

    public:

        BMatrix() = default;

        virtual ~BMatrix() = default;

        BMatrix(const BMatrix&) = default;

        BMatrix& operator=(const BMatrix&) = default;

        // No move constructor or move assignment
        BMatrix(BMatrix&&) = delete;
        BMatrix& operator=(BMatrix&&) = delete;

        /// Reset the b-matrix and clear event list.
        void Reset();

        /// Set the gyromagnetic ratio.
        /**
            \b Input:
            \n dGamma

            \b Output:
            \n m_dGamma

            Since the diffusion encoding depends on the actual
            nucleus, the gyromagnetic ratio needs to be defined.
            By default, 1H is assumed (2 pi * 42.5756e6).
            */
        void setGamma(
            double dGamma                   /**< Imp: Gyromagnetic ratio [1/T s] */
            );

        /// Set the gyromagnetic ratio.
        /**
            \b Input:
            \n rNucleus

            \b Output:
            \n m_dGamma
            */
        void setGamma(
            const MeasNucleus &rNucleus     /**< Imp: Actual nucleus */
            );

        /// Return the error status of this module.
        /**
            \b Input:
            \n bDumpMessage

            \b Output:
            \n Trace of status message

            \b Return value:
            \n BMatrix error status (0 = no error)
            */
        long lGetStatus(
            bool bDumpMessage = false       /**< Imp: Switch enables trace dump of status message */
            ) const;


        /// Calculate all b-matrix elements and 0^{th} moments, return all b-matrix elements
        /**	\b Input:
            \n m_sEventMap

            \b Output:
            \n m_dBxx, m_dBxy, m_dBxz, m_dByy, m_dByz, m_dBzz,
            m_dMx, m_dMy, m_dMz

            \b Return value:
            \n true:   0th moments are zero
            \n false:  at least one 0th moment is not zero

            This is how the calculation works:

            The b-matrix elements are defined as follows:

            b_{ij} = \gamma^2 \int_0^{TE} dt' F_i(t') F_j(t')

            where

            F_i(t') = \int_0^{t'} dt G_i(t)                     | 0   <= t' <= T_1

            F_i(t') = \int_{T_1}^{T_2} dt G_i(t) - F_i(T_1)     | T_1 <  t' <= T_2

            ...

            T_1, T_2, ... denote the time points where refocusing RF
            pulses are applied. Essentially, F_i describes the 0^{th}
            gradient moment.
            (See for example Jones et al., MRM 42:515 (1999))

            The elements of the b-matrix can be calculated step-wise.
            In order to do so, the diffusion encoding time is split
            into time slices t_0 = 0, t_1, t_2, ..., t_N = TE. The t_i
            are defined by the diffusion encoding events: either a
            gradient ramp (up/down) starts or ends, or an RF pulse
            (excitation/refocusing) is applied. Within each time slice
            from t_k to t_{k+1}, the gradient time course can be described
            as a linear function (assuming trapezoid gradient pulses):

            G_i(t) = G_i(t_k) + \frac{\Delta G_i}{\Delta t} (t - t_k)

            = G_i(t_k) + S_i(t_k) (t - t_k)

            with

            S_i(t_k) = \frac{G_i(t_{k+1}) - G_i(t_k)}{t_{k+1} - t_k}

            Next, consider the first moment F_i in the interval from t_k
            to t_{k+1}:

            F_i(t) = F_i(t_k) + \frac{\Delta F_i}{\Delta t} (t - t_k) + \frac{1}{2} \frac{\Delta^2 F_i}{\Delta t^2} (t - t_k)^2

            Starting at F_i(t_k), the evolution can be described by a
            polynomial of degree 2. It turns out that

            \frac{\Delta F_i}{\Delta t}     = G_i(t_k)

            \frac{\Delta^2 F_i}{\Delta t^2} = S_i(t_k)

            so that

            F_i(t) = F_i(t_k) + G_i(t_k) (t - t_k) + \frac{1}{2} S_i(t_k) (t - t_k)^2

            If the time point t_k falls together with a refocusing pulse,

            F_i(t_k+) = - F_i(t_k-)

            For the calculation of b_{ij}, the product of F_i and F_j has to
            be considered. Again, in the interval from t_k to t_{k+1}, one gets:

            F_i(t)F_j(t) =     F_i(t_k)F_j(t_k)

            +             (F_i(t_k)G_j(t_k) + F_j(t_k)G_i(t_k))                      * (t - t_k)

            + \frac{1}{2} (F_i(t_k)S_j(t_k) + F_j(t_k)S_i(t_k) + 2 G_i(t_k)G_j(t_k)) * (t - t_k)^2

            + \frac{1}{2} (G_i(t_k)S_j(t_k) + G_j(t_k)S_i(t_k))                      * (t - t_k)^3

            + \frac{1}{4} (S_i(t_k)S_j(t_k))                                         * (t - t_k)^4

            Finally, b_{ij} is calculated in the same interval. From the product
            F_i F_j, it is obvious that b is a polynomial of degree five:

            b_{ij}(t) = \sum_{n=0}^5 \frac{1}{n!} \frac{\Delta^n b}{\Delta t^n} (t - t_k)^n

            It turns out that

            \frac{\Delta b}{\Delta t}     =              F_i(t_k)F_j(t_k)

            \frac{\Delta^2 b}{\Delta t^2} =             (F_i(t_k)G_j(t_k) + F_j(t_k)G_i(t_k))

            \frac{1}{2}  \frac{\Delta^3 b}{\Delta t^3} = \frac{1}{2} (F_i(t_k)S_j(t_k) + F_j(t_k)S_i(t_k) + 2 G_i(t_k)G_j(t_k))

            \frac{1}{6}  \frac{\Delta^4 b}{\Delta t^4} = \frac{1}{2} (G_i(t_k)S_j(t_k) + G_j(t_k)S_i(t_k))

            \frac{1}{24} \frac{\Delta^5 b}{\Delta t^5} = \frac{1}{4} (S_i(t_k)S_j(t_k))

            so that

            b_{ij}(t) =          b_{ij}(t_k)

            +                F_i(t_k)F_j(t_k)                                          * (t - t_k)

            + \frac{1}{2}   (F_i(t_k)G_j(t_k) + F_j(t_k)G_i(t_k))                      * (t - t_k)^2

            + \frac{1}{6}   (F_i(t_k)S_j(t_k) + F_j(t_k)S_i(t_k) + 2 G_i(t_k)G_j(t_k)) * (t - t_k)^3

            + \frac{1}{8}   (G_i(t_k)S_j(t_k) + G_j(t_k)S_i(t_k))                      * (t - t_k)^4

            + \frac{1}{20}  (S_i(t_k)S_j(t_k))                                         * (t - t_k)^5

            For calculating b, it is thus sufficient to know the gradient slope S
            within each interval and the positions of the refocusing pulses.
            Gradient G, moment F and b-matrix B can then be calculated successively
            from interval to interval.


            For crosscheck purposes, results can be compared with the analytical
            expression for e.g. a Stejskal scheme with trapezoidal gradient pulses.
            See e.g. Matiello et al., JMR 108:131 (1994):

            b_{ij} = \gamma^2 G_i G_j ( \delta^2 (\Delta - \delta/3) + \epsilon^3/30 - \delta \epsilon^2/6)

            with

            \delta:   effective gradient pulse duration (ramp-up + flat-top)

            \Delta:   time between start of ramp-up of the two gradients

            \epsilon: ramp-up time

            */
        bool bCalcBMatrix(
            double &dBxx,   /**< Exp: Matrix element Bxx [s/mm^2] */
            double &dBxy,   /**< Exp: Matrix element Bxy [s/mm^2] */
            double &dBxz,   /**< Exp: Matrix element Bxz [s/mm^2] */
            double &dByy,   /**< Exp: Matrix element Byy [s/mm^2] */
            double &dByz,   /**< Exp: Matrix element Byz [s/mm^2] */
            double &dBzz    /**< Exp: Matrix element Bzz [s/mm^2] */
            );

        bool bCalcBValue(
            double &dBValue     /**< Exp: Scalar b-value (Bxx + Byy + Bzz) [s/mm^2] */
            );


        /// Dump the b-matrix elements and 0^{th} moments.
        void Dump() const;

        /// Add diffusion encoding gradient to event list.
        /**
            Two variants for implicit / explicit axis definition.

            \b Input:
            \n sGradient, lStartTime, eAxis

            \b Output:
            \n new entries in m_sEventMap, m_lStatus

            \b Return value:
            \n true:  success
            \n false: no success
            */
        bool bAddEvent(
            long lStartTime,                                /**< Imp: Event start time [us] */
            const sGRAD_PULSE_TRAP &sGradient               /**< Imp: Gradient pulse with valid amplitude and timing */
            );

        bool bAddEvent(
            long lStartTime,                                /**< Imp: Event start time [us] */
            const sGRAD_PULSE_TRAP &sGradient,              /**< Imp: Gradient pulse with valid amplitude and timing */
            SEQ::GradientAxis eAxis                         /**< Imp: Gradient axis         */
            );

        /// Add RF pulse to event list. An excitation pulse will reset the magnetization evolution.
        /**
            Two variants for implicit / explicit RF pulse type definition.

            Note: TXTYPE_UNDEFINED is used to toggle a store / restore condition. This can be
            used for stimulated echo experiments. Gradients applied between two RF pulses
            of this type do not contribute to the diffusion encoding (however, the
            evolution time is still appropriately considered).

            \b Input:
            \n sRFPulse, lStartTime, eTxType

            \b Output:
            \n new entry in m_sEventMap, m_lStatus

            \b Return value:
            \n true:  success
            \n false: no success
            */
        bool bAddEvent(
            long lStartTime,                                /**< Imp: Event start time [us]  */
            const IRF_PULSE *sRFPulse                       /**< Imp: RF pulse of valid type */
            );

        bool bAddEvent(
            long lStartTime,                                /**< Imp: Event start time [us]  */
            eTXTYPE eTxType                                 /**< Imp: RF pulse type          */
            );

        /// Add diffusion encoding gradient to event list.
        /**
            \b Input:
            \n lStartTime, dAmplitude, lRampUp, lRampDownm, lFlatTop, eAxis

            \b Output:
            \n four new entry in m_sEventMap (start ramp-up, end ramp-up,
            start ramp-down, end ramp-down), m_lStatus

            \b Return value:
            \n true:  success
            \n false: no success   (SEQ::AXIS_UNDEFINED)
            */
        bool bAddGradient(
            long lStartTime,            /**< Imp: Event start time [us] */
            double dAmplitude,          /**< Imp: Gradient amplitude [mT/m] */
            long lRampUp,               /**< Imp: Ramp up time [us] */
            long lRampDown,             /**< Imp: Ramp down time [us] */
            long lFlatTop,              /**< Imp: Flat top time [us] */
            SEQ::GradientAxis eAxis     /**< Imp: Gradient axis - one of SEQ::AXIS_PHASE, SEQ::AXIS_READ, SEQ::AXIS_SLICE
                                                  (convention: use PHASE for X, READ for Y, SLICE for Z) */
                                                  );

        /// Add RF pulse to event list. An excitation pulse will reset the magnetization evolution.
        /**
            \b Input:
            \n lStartTime, eTxType

            \b Output:
            \n new entry in m_sEventMap, m_lStatus

            \b Return value:
            \n true:  success
            \n false: no success
            */
        bool bAddRFPulse(
            long lStartTime,            /**< Imp: Event start time [us] */
            eTXTYPE eTxType             /**< Imp: RF pulse type */
            );


        // RF-type-specific functions still used by non-ep2d sequences
        bool bAddExcitationPulse(
            long lStartTime             /**< Imp: Event start time [us] */
            );

        bool bAddRefocussingPulse(
            long lStartTime             /**< Imp: Event start time [us] */
            );

        bool bAddReStorePulse(
            long lStartTime             /**< Imp: Event start time [us] */
            );

    protected:

        // ---------------------------------------------------------------------------
        // Protected type definitions
        // ---------------------------------------------------------------------------

        /// This structure holds the information of each diffusion encoding event.
        struct BMatrixEvent
        {
            BMatrixEvent() = default;

            // parametrized constructor for easier creation
            BMatrixEvent(eTXTYPE eType, SEQ::GradientAxis eAxis, double dSlope, long lStartTime)
                : eType(eType), eAxis(eAxis), dSlope(dSlope), lStartTime(lStartTime)
            {}

            eTXTYPE           eType;                            /**< Excitation, refocusing or undefined (=> ReStore)  */
            SEQ::GradientAxis eAxis;                            /**< Phase, read, slice or undefined (=> RF)            */
            double            dSlope;                           /**< Gradient slope until next event [mT/m/us]          */
            long              lStartTime;                       /**< Start time [us]                                    */
        };                                        


        // ------------------------------------
        // Protected methods
        // ------------------------------------

        /// Calculate one individual element of the b-matrix and the corresponding gradient moments
        /**	\b Input:
            \n m_sEventMap, eAxis1, eAxis2

            \b Output:
            \n dMoment1, dMoment2

            \b Return value:
            \n b-matrix element B_{12} [s/mm^2]

            Note that the member variables m_dB12, m_dM1, m_dM2 are NOT touched!
            This function is used internally by the public calculation methods.
            */
        double dCalcB12(
            SEQ::GradientAxis eAxis1,   /**< Imp: 1st gradient axis */
            SEQ::GradientAxis eAxis2,   /**< Imp: 2nd gradient axis */
            double &dMoment1,           /**< Exp: 1st axis moment [mT/m ms] */
            double &dMoment2            /**< Exp: 2nd axis moment [mT/m ms] */
            );

        // ------------------------------------
        // Member variables
        // ------------------------------------
        // If there are any changes to the member variables,
        // remember to update the Reset method!

        /// Error status of this module
        long m_lStatus{0};
        /// Error status text of this module
        std::string m_sStatus{"Ok"};

        /// Sorted list of diffusion encoding events. See definition of BEventMap.
        std::vector<BMatrixEvent> m_sEventMap;

        /// The \f$b_{xx}\f$ component of the b matrix [s/mm^2]
        double m_dBxx{0.0};
        /// The \f$b_{xy}\f$ component of the b matrix [s/mm^2]
        double m_dBxy{0.0};
        /// The \f$b_{xz}\f$ component of the b matrix [s/mm^2]
        double m_dBxz{0.0};
        /// The \f$b_{yy}\f$ component of the b matrix [s/mm^2]
        double m_dByy{0.0};
        /// The \f$b_{yz}\f$ component of the b matrix [s/mm^2]
        double m_dByz{0.0};
        /// The \f$b_{zz}\f$ component of the b matrix [s/mm^2]
        double m_dBzz{0.0};

        /// The 0^{th} gradient moment along the x-direction [mT/m us]
        double m_dMx{0.0};
        /// The 0^{th} gradient moment along the y-direction [mT/m us]
        double m_dMy{0.0};
        /// The 0^{th} gradient moment along the z-direction [mT/m us]
        double m_dMz{0.0};

        /// The gyromagnetic ratio [1/T s]
        double m_dGamma{dDefaultGamma};

        /// Indicates that all elements of the b-matrix have been calculated
        bool m_bIsBMatrixValid{false};
        /// Indicates that only the diagonal elements have been calculated
        bool m_bIsBValueValid{false};

        /// Indicates that the event list is sorted (by time)
        bool m_bIsSorted{false};

        long m_lDebugLevel{0};

        SequenceDebugSettings::SequenceDebugSettings m_debugSettings = SequenceDebugSettings::SequenceDebugSettings("USE_EPI_DEBUG_SETTINGS");

    };

    inline void   BMatrix::setGamma(double dGamma)
    {
        m_dGamma = dGamma;
    }

    inline void   BMatrix::setGamma(const MeasNucleus &rNucleus)
    {
        m_dGamma = rNucleus.getLarmorConst() * 2. * M_PI * 1.e6;
    }

}//end of namespace SEQ_NAMESPACE


