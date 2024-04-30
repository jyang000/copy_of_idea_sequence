/*!
***************************************************************************
\file   SBBDiffusion_Bipolar.h

\brief  Interface of the class Diffusion_Bipolar for diffusion mode "Bipolar"

This file provides the interface information of class Diffusion_Bipolar.

<b>Archive Information:</b>
\verbatim
File-name: SBBDiffusion_Bipolar.h
Archive File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\SBBDiffusion_Bipolar.h
\endverbatim

\b Language: C++

\author PLM AW Neuro

\b Copyright: &copy; Siemens AG (http://www.siemensmedical.com/MR).
All rights reserved.
This software may be only modified as long as the original author is credited
in any subsequent revisions or modifications.
This software must not be sold or distributed as part of any commercial
software package without the written permission of the author.

***************************************************************************

\changed     20-Jan-2003; M.Zwanger; 4a21a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- file created from SBBDiffusion.h
- comments changed to doxygen format (http://www.doxygen.org)
- BOOL changed to bool
- Tensor directions controlled by class DiffusionDirections

***************************************************************************

\changed     27-Nov-2003; M.Zwanger; 4a25a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- include paths adapted for VA25A archive
- DG1 renamed in DG1p etc.
- class calcBMatrix() added
- calcTiming() has no longer MeasNucleus in parameter list

***************************************************************************
*/


/// double include protection
#ifndef SBBDiffusion_Bipolar_h
#define SBBDiffusion_Bipolar_h 1


// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Base.h"        // base class
#include "MrImaging/seq/a_ep2d_diff/BMatrix.h"             // BMatrix

#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sFREQ_PHASE.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"     // import/export control
#include "MrImagingFW/libSBBFW/SBBList.h"

#include <map>
#include <vector>
#include <algorithm>

// --------------------------------------------------------------------------
// Forward declaration
// --------------------------------------------------------------------------
class SBBTestHelperDiffusion;      // SBB iTest Helper Class to Dump Class Member


namespace SEQ_NAMESPACE
{

    // ***************************************************************************
    // class Diffusion_Bipolar
    // ***************************************************************************

    class __IMP_EXP SBBDiffusion_Bipolar: public SBBDiffusion_Base
    {

        // ------------------------------------
        // Public typedefs and definitions
        // ------------------------------------

    public:

        /// This enumeration includes all relevant event timings
        /** Timing information is required for different purposes (e.g. GPA load
        calculation, b-matrix calculation). The method getTimingInformation
        is the central instance that distributes this information.
        */
        enum EnumEventTime
        {
            DiffGrad1_Start = 1,
            DiffGrad2_Start,
            DiffGrad3_Start,
            DiffGrad4_Start,
            SpoilGrad1_Start,
            SpoilGrad2_Start,
            SpoilGrad3_Start,
            SpoilGrad4_Start,
            SliceGrad1_Start,
            SliceGrad2_Start,
            RefocRF1_Center,
            RefocRF2_Center
        };

        // ------------------------------------
        // Public methods
        // ------------------------------------

        ///	The constructor initializes the starting time of the diffusion gradients 
        ///   with 0 and the maximum possible gradient amplitudes.
        SBBDiffusion_Bipolar(SBBList* pSBBList);

        ///   This destructor does nothing.
        virtual ~SBBDiffusion_Bipolar();

        friend class SBBTestHelperDiffusion;      // SBB iTest Helper Class to Dump Class Member

        /// Actual implementation of the virtual methods of the base class
        /**
        This method calculates the reference spoil moment (m_dRefSpoilMoment) and
        prepares the spoiler gradients on one reference axis (PE-axis)
        */
        virtual bool prepInit(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo);

        /**
        This method prepares the diffusion encoding gradients (using an elaborate
        GPA model) and the refocusing pulses. The maximum possible b-value obtainable
        with the current protocol is calculated and stored (m_dMaxPossibleBValue).
        */
        virtual bool prepTiming(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo, long lActualTE);

        virtual bool prepRF(MrProt & rMrProt, SeqLim &rSeqLim, SeqExpo & rSeqExpo);

        /**
        This method prepares all spoiler gradients (based on the already prepared
        reference spoiler) and all diffusion gradients (based on the already prepared
        reference gradient). The provided b-value has to be realized. A required
        TR increment is stored (m_lTRIncrement). All exports (see base class) are set.
        */
        virtual bool prepFinal( double dMaxRequestedBValue, bool bIsContextPrepForBinarySearch = false);

        ///	Play out the run-time events for the scan actually selected.
        virtual bool runSBB(MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo, sSLICE_POS* pSLC);

        ///  Calculate the b-value which will be obtained by the current timing.
        /**  For the calculation, it is supposed that the norms of all direction
        vectors are equal. (If not, be informed that the norm of the last
        direction vector is used for calculation.)

        \b Input:
        \n m_dGamma, all diffusion encoding events

        \b Output:
        \n m_theBMatrix (calculated diagonal elements)

        \b Return value:
        \n actual b-value [s/mm^2]
        */
        virtual double calcBValue(void);

        /// Calculate the complete b-matrix
        /**
        \b Input:
        \n m_dGamma, bApplySpoiler, all diffusion encoding events

        \b Output:
        \n m_theBMatrix, m_dBxx, m_dBxy, ...
        */
        virtual void calcBMatrix(void);

        /// Prepare gradient events for b-value or b-matrix calculation, respectively
        /**
        \b Input:
        \n bBxxOnly, bApplySpoiler

        \b Output:
        \n modified m_theBMatrix; if !bBxxOnly also m_dBxx, m_dBxy, ...

        \b Return value:
        \n Success = true, failure = false
        */
        virtual bool PrepBMatrix(bool bBValueOnly = false);       /**< If argument = true, only the b-value is calculated */

        /// Actual implementation of the virtual methods of the base class
        virtual bool prepGPALoadDiff(double dAmplitudeX);

        virtual bool prepGPALoadDiff(double dAmplitudeX, double dAmplitudeY, double dAmplitudeZ);

        virtual bool scaleGPALoadDiff(double dScaleX);

        virtual bool scaleGPALoadDiff(double dScaleX, double dScaleY, double dScaleZ);

        /// Get timing information
        /**
        \b Input:
        \n eEvent

        \b Output:
        \n n.a.

        \b Return value:
        \n Requested event time (-1 in case of unknown event)

        */
        virtual long getEventTime(EnumEventTime eEvent);

        /// See SeqBuildBlock for detailed information: Calculate RF info
        //  Note: Instead of adding a complete list of cuboids to this PlugIn, we just
        //        provide the current amplitude factor (used to scale the flip angles).
        //        Usually, this factor is available from the superior instance (i.e. SBBEPIKernel).
        virtual bool calcSliceAdjSBBRFInfo(
            MrProt                                     &rMrProt,        //< IMP: points to the protocol structure.
            SeqLim                                     &rSeqLim,        //< IMP: points to the sequence limits structure.
            SeqExpo                                    &rSeqExpo,       //< IMP: points to the sequence exports structure
            const SLICEADJ::sCuboidGeometry              &sSliceAdjCuboid,  //< IMP: cuboid geometry (single)
            std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo        //< EXP: RF info
            );

        /// Get stored RF info for a certain geometry
        virtual bool getSliceAdjRFInfo(
            const SLICEADJ::sCuboidGeometry              &sSliceAdjCuboid,  //< IMP: cuboid geometry
            std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo        //< EXP: RF info
            );

        virtual MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest();

        // determine smallest b value for IVIM
        virtual long getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch = false);

    protected:

        // ------------------------------------
        // Protected methods
        // ------------------------------------

        // ------------------------------------
        // Member variables
        // ------------------------------------

        ///    This is an additional fill time inserted after the last diffusion gradient.
        /**   Usually it is zero. But as some sequence are not able to handle the fill
        times properly (lazyness, ignorance, dullness ;-), this class tries to
        recognize and fix TE fill problems using this variable.
        */
        long m_AddFillTime;


        //Static variables in the old code. They're now declared as member variables
        //to avoid violation between threads.
        double   m_dSpoilerAmplitude;
        double   m_dSpoilerFactor1;
        double   m_dSpoilerFactor2;
        bool     m_bRunOnlyOnce;             // Hook for one-time jobs   (usually executed during check() phase)
        bool     m_bDumpBMatrix;             // Hook for one-time jobs   (usually executed during check() phase)

        /// Slice thickness for which the RF pulse has been prepared
        /**  Timing analysis showed that most of the function execution time
        is spent for the preparation of the external RF pulses.
        Therefore the pulses are only prepared if necessary, i.e.
        for a changed slice thickness. This variable is used to remember
        the slice thickness of the last preparation.
        */
        double   m_dPreparedSlcThk;

        /// ID's of gradient events stored within m_sBalanceDiff
        long     m_lID1X;
        long     m_lID2X;
        long     m_lID3X;
        long     m_lID4X;

        long     m_lID1Y;
        long     m_lID2Y;
        long     m_lID3Y;
        long     m_lID4Y;

        long     m_lID1Z;
        long     m_lID2Z;
        long     m_lID3Z;
        long     m_lID4Z;


    private:

        ///   Instance of the BMatrix class
        /**     Required for the calculation of the b-matrix elements and the trace
        thereof (b-value) based on a series of diffusion encoding events.
        */
        BMatrix m_theBMatrix;

        ///   An additional fill time at the start of the SBB
        /**    Some diffusion schemes require that the SpinPrepTime and the ADCusTillEcho
        lay on a special gradient raster. If these conditions are not met, the time
        specified by this variable will be sent immediately at the start of the run
        method of this SBB. It is set by the prep methods of the SBB.
        Generally the sequence should avoid this condition by providing a proper Spin
        PrepTime and ADCusTillEcho.
        */
        long m_lSpinPrepTimeEnhancement;


        /// Refocussing RF pi pulse
        sRF_PULSE_EXT m_DRF1;
        /// Second RF pi pulse
        sRF_PULSE_EXT m_DRF2;

        /// First diffusion-encoding gradient (in phase direction)
        sGRAD_PULSE_TRAP m_DG1p;
        /// First diffusion-encoding gradient (in read direction)
        sGRAD_PULSE_TRAP m_DG1r;
        /// First diffusion-encoding gradient (in slice direction)
        sGRAD_PULSE_TRAP m_DG1s;

        /// Second diffusion-encoding gradient (in phase direction)
        sGRAD_PULSE_TRAP m_DG2p;
        /// Second diffusion-encoding gradient (in read direction)
        sGRAD_PULSE_TRAP m_DG2r;
        /// Second diffusion-encoding gradient (in slice direction)
        sGRAD_PULSE_TRAP m_DG2s;

        /// Third diffusion-encoding gradient (in phase direction)
        sGRAD_PULSE_TRAP m_DG3p;
        /// Third diffusion-encoding gradient (in read direction)
        sGRAD_PULSE_TRAP m_DG3r;
        /// Third diffusion-encoding gradient (in slice direction)
        sGRAD_PULSE_TRAP m_DG3s;

        /// Fourth diffusion-encoding gradient (in phase direction)
        sGRAD_PULSE_TRAP m_DG4p;
        /// Fourth diffusion-encoding gradient (in read direction)
        sGRAD_PULSE_TRAP m_DG4r;
        /// Fourth diffusion-encoding gradient (in slice direction)
        sGRAD_PULSE_TRAP m_DG4s;

        // map of arrays to store pointers to all diffusion gradients
        std::map< std::string, std::array<sGRAD_PULSE_TRAP*, 4> > m_maDG;


        /// Spoiler around the first refocussing pulse for b > 0 in phase direction
        /** This spoiler is played out in the event-block that contain the
        first diffusion pulse, i.e. it uses the rotation matrix of the
        diffusion gradients.
        */
        sGRAD_PULSE_TRAP m_DSpoil1P;
        /// Spoiler around the second  refocussing pulse for b > 0 in phase direction
        /** This spoiler is played out in the event-block that contains the
        second diffusion pulse, i.e. it uses the rotation matrix of the
        diffusion gradients.
        */
        sGRAD_PULSE_TRAP m_DSpoil2P;

        /// Spoiler around the first refocussing pulse for b > 0 in readout direction
        /** This spoiler is played out in the event-block that contain the
        first diffusion pulse, i.e. it uses the rotation matrix of the
        diffusion gradients.
        */
        sGRAD_PULSE_TRAP m_DSpoil1R;
        /// Spoiler around the second  refocussing pulse for b > 0 in readout direction
        /** This spoiler is played out in the event-block that contains the
        second diffusion pulse, i.e. it uses the rotation matrix of the
        diffusion gradients.
        */
        sGRAD_PULSE_TRAP m_DSpoil2R;

        /// Spoiler around the first refocussing pulse for b > 0 in slice direction
        /** This spoiler is played out in the event-block that contain the
        first diffusion pulse, i.e. it uses the rotation matrix of the
        diffusion gradients.
        */
        sGRAD_PULSE_TRAP m_DSpoil1S;
        /// Spoiler around the second  refocussing pulse for b > 0 (in slice direction)
        /** This spoiler is played out in the event-block that contains the
        second diffusion pulse, i.e. it uses the rotation matrix of the
        diffusion gradients.
        */
        sGRAD_PULSE_TRAP m_DSpoil2S;

    };

}// end of namespace SEQ_NAMESPACE
#endif // of ifndef SBBDiffusion_Bipolar_h
