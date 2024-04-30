/*!
***************************************************************************
\file   SBBDiffusion_Trace.h

\brief  Interface of class Diffusion_Trace.

This file provides the class Diffusion_Trace.
For further comments, please refer to the class description of Diffusion_Trace.

<b>Archive Information:</b>
\verbatim
File-name: "h:\ep2d_diff_WIP\SBBDiffusion_Trace.h"
Archive File: \n4\pkg\MrServers\MrImaging\seq\a_ep2d_diff\SBBDiffusion_Trace.h
\endverbatim

\b Language: C++

\author Michael Zwanger (alias mizwa, Michael.Zwanger@siemens.com)

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
- file created by splitting SBBDiffusion.h
- comments changed to doxygen format (http://www.doxygen.org)
- BOOL changed to bool

***************************************************************************

\changed     27-Nov-2003; M.Zwanger; 4a25a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- include paths adapted for VA25A archive

***************************************************************************
*/





/// double include protection:
#ifndef SBBDiffusion_Trace_h
#define SBBDiffusion_Trace_h 1



// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------


#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Base.h"        // base class
#include "MrImaging/seq/a_ep2d_diff/BMatrix.h"             // BMatrix

#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sFREQ_PHASE.h"
#include "MrMeasSrv/SeqIF/libRT/sSYNC.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"  // import/export control
#include "MrImagingFW/libSBBFW/SBBList.h"


namespace SEQ_NAMESPACE
{

    // ===========================================================================
    /*!
    \class  Diffusion_Trace

    \brief  This class implements a one-scan trace diffusion scheme.

    This class is derived from the SeqBuildBlockDiffusion class and provides a
    trace-weighted diffusion measurement by a single-scan, i.e. the whole
    spatial information is encoded during a single call of the run() method.
    This sequene scheme has been invented by Oliver Heid.
    This module consists of 16 diffusion gradients and one refocussing RF pulse.
    In the UI, it is connected to the Diffusion mode called "trace" (VA12) or
    "1-scan trace" (from VA15A on).

    \image  html Diffusion_Trace.gif "Timing diagram of Diffusion_Trace"

    Since it is desirable to run all three gradients simultaneously at
    an amplitude as high as possible, the diffusion gradients are applied
    in the fixed coordinate system of the magnet. The gradient pattern is
    therefore fixed for all slice orientations and does not limit the use of
    obliquely orientated slices. This is allowed due to the rotation
    invariance of the trace.

    If the b value is smaller than the SPOILER_THRESHOLD, no diffusion gradients,
    but spoiler gradients will be applied instead. The spoilers are prepared by
    calling SeqBuildBlockDiffusion::prepSpoilGrad().

    \image  html Diffusion_Trace0.gif "Timing diagram of Diffusion_Trace for b=0"



    \author Michael.Zwanger@med.siemens.de
    */

    class __IMP_EXP SBBDiffusion_Trace : public SBBDiffusion_Base
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
            DiffGrad5_Start,
            DiffGrad6_Start,
            DiffGrad7_Start,
            DiffGrad8_Start,
            DiffGrad9_Start,
            DiffGrad10_Start,
            DiffGrad11_Start,
            DiffGrad12_Start,
            DiffGrad13_Start,
            DiffGrad14_Start,
            DiffGrad15_Start,
            DiffGrad16_Start,
            SpoilGrad1_Start,
            SpoilGrad2_Start,
            SliceGrad1_Start,
            RefocRF1_Center,
        };

        ///    Class constructor.
        /**           Initializes the gradient sign arrays and maximum gradient amplitudes.
        */
        SBBDiffusion_Trace(SBBList* pSBBList);

        virtual ~SBBDiffusion_Trace();


        /// Actual implementation of the virtual methods of the base class
        /**
        This method calculates the reference spoil moment (m_dRefSpoilMoment) and
        prepares the spoiler gradients on one reference axis (PE-axis)
        */
        virtual bool prepInit   (MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo);

        /**
        This method prepares the diffusion encoding gradients (using an elaborate
        GPA model) and the refocusing pulses. The maximum possible b-value obtainable
        with the current protocol is calculated and stored (m_dMaxPossibleBValue).
        */
        virtual bool prepTiming (MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo, long lActualTE);

        /**
        This method prepares all spoiler gradients (based on the already prepared
        reference spoiler) and all diffusion gradients (based on the already prepared
        reference gradient). The provided b-value has to be realized. A required
        TR increment is stored (m_lTRIncrement). All exports (see base class) are set.
        */
        virtual bool prepFinal  ( double dMaxRequestedBValue, bool bIsContextPrepForBinarySearch = false);

        ///	Play out the run-time events for the scan actually selected.
        virtual bool runSBB  (MrProt &rMrProt, SeqLim &rSeqLim, MrProtocolData::SeqExpo &rSeqExpo, sSLICE_POS* pSLC);


        ///  Calculate the b-value which will be obtained by the current timing.
        /**   
        For the calculation, it is supposed that the norms of all direction
        vectors are equal. (If not, be informed that the norm of the last
        direction vector is used for calculation.)

        \b Input:
        \n m_dGamma, all diffusion encoding events

        \b Output:
        \n m_theBMatrix (calculated diagonal elements)

        \b Return value:
        \n actual b-value [s/mm^2]
        */
        virtual double calcBValue (void);

        /// Calculate the complete b-matrix
        /**
        \b Input:
        \n all diffusion encoding events

        \b Output:
        \n m_theBMatrix, m_dBxx, m_dBxy, ...
        */
        virtual void calcBMatrix (void);

        /// Calculate the complete b-matrix
        /**
        \b Input:
        \n m_dGamma, bApplySpoiler, all diffusion encoding events

        \b Output:
        \n m_theBMatrix, m_dBxx, m_dBxy, ...
        */
        virtual void calcBMatrix 
            (
            bool bApplySpoiler             /**< Indicate whether spoiler or diffusion encoding gradients are applied */
            );

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

        // Since this module does not use the GPA model, these methods are not required. Still,
        // they have to be implemented (pure virtual methods in base class)
        virtual bool prepGPALoadDiff  ( double dAmplitudeX                                         );

        virtual bool prepGPALoadDiff  ( double dAmplitudeX, double dAmplitudeY, double dAmplitudeZ );

        virtual bool scaleGPALoadDiff ( double dScaleX                                 );

        virtual bool scaleGPALoadDiff ( double dScaleX, double dScaleY, double dScaleZ );

        /// Get timing information
        /**
        \b Input:
        \n eEvent

        \b Output:
        \n n.a.

        \b Return value:
        \n Requested event time (-1 in case of unknown event)

        */
        virtual long getEventTime (EnumEventTime eEvent);

        /// See SeqBuildBlock for detailed information: Calculate RF info
        //  Note: Instead of adding a complete list of cuboids to this PlugIn, we just
        //        provide the current amplitude factor (used to scale the flip angles).
        //        Usually, this factor is available from the superior instance (i.e. SBBEPIKernel).
        bool calcSliceAdjSBBRFInfo(
            MrProt                                     &rMrProt,        //< IMP: points to the protocol structure.
            SeqLim                                     &rSeqLim,        //< IMP: points to the sequence limits structure.
            SeqExpo                                    &rSeqExpo,       //< IMP: points to the sequence exports structure
            const SLICEADJ::sCuboidGeometry              &sSliceAdjCuboid,  //< IMP: cuboid geometry (single)
            std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo        //< EXP: RF info
            );

        // determine smallest b value for IVIM
        virtual long getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch = false);

    private:
        // Data Members for Class Attributes

        ///   Instance of the BMatrix class
        /**     Required for the calculation of the b-matrix elements and the trace
        thereof (b-value) based on a series of diffusion encoding events.
        */
        BMatrix m_theBMatrix;

        ///    This array contains the signs (either +1 or -1) for each diffusion gradient
        /**           applied on the phase gradient axis.
        */
        int m_PhaseGradSign[16];

        ///    This array contains the signs (either +1 or -1) for each diffusion gradient
        /**           applied on the read gradient axis.
        */
        int m_ReadGradSign[16];

        ///    This array contains the signs (either +1 or -1) for each diffusion gradient
        /**           applied on the slice gradient axis.
        */
        int m_SliceGradSign[16];

        ///    This fill time (in microseconds) will be sent before
        /// the first gradient will be played out in the run method.
        long m_lPreFill;

        ///    This fill time (in microseconds) will be sent after the last gradient
        ///has been played out in the run method.
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

        // Data Members for Associations

        sRF_PULSE_EXT m_DRF;

        sFREQ_PHASE m_DFPset;

        sFREQ_PHASE m_DFPneg;

        sGRAD_PULSE_TRAP m_DGR[16];

        sGRAD_PULSE_TRAP m_DGP[16];

        sGRAD_PULSE_TRAP m_DGS[16];

        sGRAD_PULSE_TRAP m_DGSS;

    private: //## implementation
        SBBDiffusion_Trace(const SBBDiffusion_Trace &right);

        SBBDiffusion_Trace & operator=(const SBBDiffusion_Trace &right);

        // Additional Implementation Declarations

    };

}//end of namespace SEQ_NAMESPACE
#endif // of ifndef SBBDiffusion_Trace_h
