/*! 
***************************************************************************
\file   SBBDiffusion_Stejskal.h

\brief  Interface of class Diffusion_Stejskal for the diffusion mode "Stejskal" 

This file provides the class Diffusion_Stejskal.
For further comments, please refer to the class description of Diffusion_Stejskal.

\author PLM AW Neuro

<b>Archive Information:</b>
\verbatim
- File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\SBBDiffusion_Stejskal.h
- Last Author: ZWANGER
- Archive Date: 2015-01-19 18:15:41 +01:00
\endverbatim

\b Language: C++

\b Copyright: &copy; Siemens AG (http://www.siemensmedical.com).
All rights reserved.   

***************************************************************************

\changed     20-Jan-2003; M.Zwanger; 4a21a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- new class 'Diffusion_Stejskal' added
- the spoiler pad moment is taken from the UI in WIP mode

***************************************************************************

\changed     27-Nov-2003; M.Zwanger; 4a25a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- include paths adapted for VA25A archive
- DG1 renamed in DG1p etc.


\changed    30-Sept-2008; k.liu, 4vb15a, 
\description 
- correct full path for *.h

***************************************************************************

\changed     01-Jan-2009; T.Feiweier; 4b15a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- complete rework: match with current SBBDiffusion_Bipolar code

***************************************************************************
*/


/// double include protection:
#ifndef SBBDiffusion_Stejskal_h
#define SBBDiffusion_Stejskal_h 1



// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Base.h"        // base class
#include "MrImaging/seq/a_ep2d_diff/BMatrix.h"             // BMatrix

#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sFREQ_PHASE.h"

#include "MrImaging/libSBB/SBBMultibandRF.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"     // import/export control
#include "MrImagingFW/libSBBFW/SBBList.h"


// --------------------------------------------------------------------------
// Forward declaration
// --------------------------------------------------------------------------
class SBBTestHelperDiffusion;      // SBB iTest Helper Class to Dump Class Member


namespace SEQ_NAMESPACE
{
    // ***************************************************************************
    // class Diffusion_Stejskal
    // ***************************************************************************

    class __IMP_EXP SBBDiffusion_Stejskal: public SBBDiffusion_Base
    {

        // ------------------------------------
        // Public typedefs and definitions
        // ------------------------------------

        friend class SBBTestHelperDiffusion;      // SBB iTest Helper Class to Dump Class Member

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
            SpoilGrad1_Start,
            SpoilGrad2_Start,
            SliceGrad1_Start,
            RefocRF1_Center,
        };


        // ------------------------------------
        // Public methods
        // ------------------------------------

        ///	The constructor initializes the starting time of the diffusion gradients 
        ///   with 0 and the maximum possible gradient amplitudes.
        SBBDiffusion_Stejskal(SBBList* pSBBList);

        ///   This destructor does nothing.
        virtual ~SBBDiffusion_Stejskal();

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

        /**
        This method prepares all spoiler gradients (based on the already prepared
        reference spoiler) and all diffusion gradients (based on the already prepared
        reference gradient). The provided b-value has to be realized. A required
        TR increment is stored (m_lTRIncrement). All exports (see base class) are set.
        */
        virtual bool prepFinal(double dMaxRequestedBValue, bool bIsContextPrepForBinarySearch = false);

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
        \n m_dGamma, all diffusion encoding events

        \b Output:
        \n m_theBMatrix, m_dBxx, m_dBxy, ...
        */
        virtual void calcBMatrix(void);

        /// Calculate the complete b-matrix
        /**
        \b Input:
        \n m_dGamma, bApplySpoiler, all diffusion encoding events

        \b Output:
        \n m_theBMatrix, m_dBxx, m_dBxy, ...
        */
        virtual void calcBMatrix(bool bApplySpoiler);             /**< Indicate whether spoiler or diffusion encoding gradients are applied */
            
        /// Prepare gradient events for b-value or b-matrix calculation, respectively
        /** 
        \b Input:
        \n bBxxOnly, bApplySpoiler

        \b Output:
        \n modified m_theBMatrix; if !bBxxOnly also m_dBxx, m_dBxy, ...

        \b Return value:
        \n Bxx or b-value (Trace of b-matrix)
        */
        virtual bool PrepBMatrix 
            (
            bool bBValueOnly,             /**< If true, only the b-value is calculated */
            bool bApplySpoiler            /**< Indicate whether spoiler or diffusion encoding gradients are applied */
            );

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


        bool prepRF(MrProt & rMrProt, SeqLim &rSeqLim, SeqExpo & rSeqExpo);


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

        // check whether the moment of the spoiler is larger than the diffusion moment
        virtual bool isDiffusionMomentSmallerThanSpoilerMoment();

        // determine smallest b value for which no spoiler around refocussing pulse(s) are used
        virtual long getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch = false);

    protected:

        // ------------------------------------
        // Protected methods
        // ------------------------------------

        // ------------------------------------
        // Member variables
        // ------------------------------------

        /// Fill time at the beginning of the SBB
        /**	This fill time (in microseconds) will be sent before 
        the first gradient will be played out in the run method. */
        long m_lPreFill;

        /// Fill time at the end of the SBB
        /**	This fill time (in microseconds) will be sent after the last gradient 
        has been played out in the run method. */
        long m_lPostFill;

        bool    m_bRunOnlyOnce;      // Hook for one-time jobs   (usually executed during check() phase)
        bool    m_bDumpBMatrix;      // Hook for one-time jobs   (usually executed during check() phase)

        /// Slice thickness for which the RF pulse has been prepared
        /**  Timing analysis showed that most of the function execution time
        is spent for the preparation of the external RF pulses.
        Therefore the pulses are only prepared if necessary, i.e.
        for a changed slice thickness. This variable is used to remember
        the slice thickness of the last preparation.
        */      
        double  m_dPreparedSlcThk;

        /// ID's of gradient events stored within m_sBalanceDiff
        long     m_lID1X;
        long     m_lID2X;

        long     m_lID1Y;
        long     m_lID2Y;

        long     m_lID1Z;
        long     m_lID2Z;



    private:

        ///   Instance of the BMatrix class
        /**     Required for the calculation of the b-matrix elements and the trace
        thereof (b-value) based on a series of diffusion encoding events.
        */
        BMatrix m_theBMatrix;

        /// Refocussing RF pi pulse
        sRF_PULSE_EXT m_DRF1;

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
    };



}// end of namespace SEQ_NAMESPACE


#endif // of ifndef SBBDiffusion_Stejskal_h
