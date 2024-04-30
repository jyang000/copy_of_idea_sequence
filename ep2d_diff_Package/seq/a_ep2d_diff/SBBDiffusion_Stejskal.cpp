//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 2010  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\SBBDiffusion_Stejskal.cpp
//
//      Author: PLM AW NEUR
//
//        Lang: C++
//
//     Descrip: Implementation of the class Diffusion_Stejskal.
//
//     Classes: Diffusion_Stejskal
//
//    -----------------------------------------------------------------------------

/**
***************************************************************************

\changed     1-Sep-2002; M.Zwanger; 4a21a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- new class 'Diffusion_Stejskal' added
- Initial phase of m_DRF2 changed to 180 deg (Thanks, Juergen :-)
- shorter RF pulse for 3T (Thanks, Larry :-)

\changed    31-Okt-2003; M.Zwanger; 4a21a; CHARM: n.a.
\description
- 'm_bDiffusionGradientsEnabled' support added
- include paths adapted for VA25A archive
- wipMemBlock().adFree bugs fixed

\changed    12-Oct-2005; M.Zwanger; 4b13a; Customer wish
\description
- Effictive Amplitude based on m_Didi.getGreatestNorm()

\changed    30-Sept-2008; k.liu, 4vb15a,
\description
- correct full path for *.h


***************************************************************************
*/


// MrProt
#include "MrProtSrv/Domain/MrProtData/MrProt/Application/Application.h"
// MrProt

#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"

#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Stejskal.h"



#include "MrImagingFW/libSeqUTIF/libsequt.h"               // mSEQTest
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

#include "DiffusionRFPulseProperties.h"


#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"


#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif
using namespace SEQ_NAMESPACE;

// ---------------------------------------------------------------------------
// Type definitions
// ---------------------------------------------------------------------------


// ===========================================================================
/*!
\class Diffusion_Stejskal

\brief This class implements Stejskal-Tanner diffusion weighting
(i.e. with one spin echo and two diffusion gradients).

This class is derived from the SeqBuildBlockDiffusion class and provides a
diffusion weighting with two bipolar diffusion gradients and one refocussing
RF pulse in between.
This scheme has been suggested by Stejskal and Tanner:
"Spin Diffusion Measurement", J. Chem. Phys. 42(1), 288 (January 1965).

It can be used for diffusion measurements in read, slice and phase direction
as well as for the orthogonal and tensor mode.

\image html Diffusion_Stejskal.gif "Timing diagram of Diffusion_Stejskal"

The timing consists of the following event blocks:
-# a fill time,
-# the first diffusion gradient lobe,
-# the refocussing RF pulse (with a slice-selection gradient)
-# the second diffusion gradient and
-# a fill time.

If the b value is smaller than the SPOILER_THRESHOLD, no diffusion gradients,
but spoiler gradients will be applied instead. The spoiler gradient moment
is fixed coded in Diffusion_Stejskal::calcTiming  or can be changed on the UI
in WIP mode. The spoilers are prepared by calling SeqBuildBlockDiffusion::prepSpoilGrad().

\image html Diffusion_Stejskal0.gif "Timing diagram of Diffusion_Stejskal for b=0"


The <b>fill times</b> are calculated automatically (from the condition that
the refocussing RF pulse must be centered in the TE interval).
Actually, there is only one fill time played out. The other fill time is zero,
because this gradient is made as long as possible. If TE is longer than
necessary, the gradient amplitude will be reduced and the whole available
time will be filled. Therefore only one fill time is necessary.

The spoiler gradients and the RF pulse sent in an event block with the
<b>rotation matrix</b> of the current slice, whereas the rotation matrix of
the diffusion gradients is determined by DiffusionDirections.


INI-configurations:

- For each section, traces can be enabled by inserting the line
TRACES=CALL,INPUT,RETURN,RESULT,LANDMARK,INTERNAL
(or a selection of the listed severities)
- Scaling of spoiler amplitudes:
[DiffusionStejskal::prep]
SpoilerFactor = 3.0
- Dump of the B-matrix:
[DiffusionStejskal::run]
DumpBMatrix = 1
- Dump of the rotation matrix
[DiffusionStejskal::run]
DumpVector = 1
- Dump Ice parameters
[DiffusionStejskal::run]
DumpICE = 1
- Scaling of individual diffusion gradient amplitudes
=> not considered within b-value calculation!
=> Gmax will never be exceeded
[DiffusionStejskal::run]
DiffGradAmplXFactor = 1.0
DiffGradAmplYFactor = 1.0
DiffGradAmplZFactor = 1.0
- Overall scaling of diffusion gradient amplitudes
=> considered within b-value calculation
=> Gmax will never be exceeded
[DiffusionStejskal::CalcMaximumAmplitude]
DiffGradScale = 1.0


*/


//   ===========================================================================
///  The constructor initializes the starting time of the diffusion gradients 
///  with 0 and the maximum possible gradient amplitudes.
SBBDiffusion_Stejskal::SBBDiffusion_Stejskal(SBBList* pSBBList):
SBBDiffusion_Base(pSBBList),
m_lPreFill(0),
m_lPostFill(0),
m_bRunOnlyOnce(false),
m_bDumpBMatrix(false),
m_dPreparedSlcThk(0.0),
m_lID1X(0),
m_lID2X(0),
m_lID1Y(0),
m_lID2Y(0),
m_lID1Z(0),
m_lID2Z(0),
m_DRF1("RefocRF1"),
m_DG1p("RTEIdentDG1p"),
m_DG1r("RTEIdentDG1r"),
m_DG1s("RTEIdentDG1s"),
m_DG2p("RTEIdentDG2p"),
m_DG2r("RTEIdentDG2r"),
m_DG2s("RTEIdentDG2s")

// ===========================================================================
{
    setIdent("Diffusion_Stejskal");

    m_DG1p.setStartTime(0);
    m_DG2p.setStartTime(0);

    // Set maximum possible(!) gradient amplitudes (used during check())
    double dAmpl = SysProperties::getGradMaxAmplAbsolute();

    m_DG1p.setMaxMagnitude(dAmpl);
    m_DG1r.setMaxMagnitude(dAmpl);
    m_DG1s.setMaxMagnitude(dAmpl);
    m_DG2p.setMaxMagnitude(dAmpl);
    m_DG2r.setMaxMagnitude(dAmpl);
    m_DG2s.setMaxMagnitude(dAmpl);

    // If the gradient system does not support balance models, use a
    // conservative amplitude for the monopolar diffusion gradients.
    if(!m_sBalanceAcq.bIsBalanceModelSupported())
    {
        m_dMaxAmpl = SysProperties::getGradMaxAmplNominal();
        m_dAmpl    = m_dMaxAmpl;
    }
}



// ===========================================================================
///   This destructor does nothing.
SBBDiffusion_Stejskal::~SBBDiffusion_Stejskal()
// ===========================================================================
{
}


bool SBBDiffusion_Stejskal::calcSliceAdjSBBRFInfo(
    MrProt                                     &rMrProt,        //< IMP: points to the protocol structure.
    SeqLim                                     &rSeqLim,        //< IMP: points to the sequence limits structure.
    SeqExpo                                    &rSeqExpo,       //< IMP: points to the sequence exports structure
    const SLICEADJ::sCuboidGeometry              &sSliceAdjCuboid,  //< IMP: cuboid geometry (single)
    std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo        //< EXP: RF info
    )
{
    // nothing to do here as all energy is handled by pulse SBBs.
    return true;
}


// ===========================================================================
/// Implementation of the pure virtual base class method
// ===========================================================================
bool SBBDiffusion_Stejskal::prepInit(MrProt & /* &rMrProt */, SeqLim & /* &rSeqLim */, SeqExpo & /* &rSeqExpo */)
// ===========================================================================
{
    // -------------------------
    // Prepare spoiler gradients
    // -------------------------
    // This must be done before timing calculation. Note that only
    // the spoilers in P-direction are available within the timing
    // calculation.
    //

    // This factor determines the separation of the undesired signal pathways
    // from the desired one. A value of 1.5 is the absolute minimum that should
    // be used here.
    double dSpoilerFactor = m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/diffusion_spoil_factor", 3.0);

    m_dRefSpoilMoment     = m_dReadoutMoment * dSpoilerFactor / 1000.;

    if(! prepSpoilGrad(m_dRefSpoilMoment))
    {
        SEQ_TRACE_ERROR.print("ERROR: prepSpoilGrad() failed.");
        return false;
    }

    return true;
}

// ===========================================================================
/// Implementation of the pure virtual base class method
// ===========================================================================
bool SBBDiffusion_Stejskal::prepFinal(double dMaxRequestedBValue, bool bIsContextPrepForBinarySearch)
// ===========================================================================
{
    // ---------------------------
    // Prepare diffusion gradients
    // ---------------------------

    // m_dMaxPossibleBValue might be higher than the desired b-value. Here, 
    // the actual gradient amplitudes are scaled to the desired b-value and
    // the corresponding TR increment is calculated.
    // Note: m_dAmpl is used within fSeqRun in order to scale the diffusion
    // gradient amplitudes to the desired b-value and must not be changed
    // any more!

    // Scale Amplitude exactly to requested maximum b value
    double dAmpl = m_dAmpl * sqrt(dMaxRequestedBValue/m_dMaxPossibleBValue);

    long   lTRIncrement = 0;

    // Calculate required TR increment
    if(!CalcTRIncrement(dAmpl, lTRIncrement))
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_WARN.print("WARNING! Inconsistent GPA balance calculation!");
        }
    }

    // Update m_lTRIncrement
    setTRIncrement(lTRIncrement);

    m_DG1p.setAmplitude(dAmpl * m_Didi.getX(0));
    m_DG1p.setDuration(m_DG1p.getDuration());    // Note: timing of m_DG1p has been set in ::prepTiming
    m_DG1p.setAxis(SEQ::AXIS_PHASE);
    m_DG1p.setRampUpTime(m_DG1p.getRampUpTime());
    m_DG1p.setRampDownTime(m_DG1p.getRampDownTime());

    if(! m_DG1p.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: m_DG1p.prep() returned false ");
        }
        setNLSStatus(m_DG1p.getNLSStatus());
        return false;
    }

    if(! m_DG1p.check())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: m_DG1p.check() returned false ");
        }
        setNLSStatus(m_DG1p.getNLSStatus());
        return false;
    }

    if(m_DG1p.getTotalTime() < 300)
    {
        // CHARM 324739: Each rotation matrix must be effective for at least 300µs before the next change
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: Rotation matrix too short for DG1");
        }
        setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }

    m_DG1r.setAmplitude(dAmpl * m_Didi.getY(0));
    m_DG1r.setDuration(m_DG1p.getDuration());
    m_DG1r.setAxis(SEQ::AXIS_READOUT);
    m_DG1r.setRampUpTime(m_DG1p.getRampUpTime());
    m_DG1r.setRampDownTime(m_DG1p.getRampDownTime());
    m_DG1r.setStartTime(m_DG1p.getStartTime());

    if(! m_DG1r.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG1r.getNLSStatus());
        return false;
    }

    if(! m_DG1r.check())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG1r.getNLSStatus());
        return false;
    }

    m_DG1s.setAmplitude(dAmpl * m_Didi.getZ(0));
    m_DG1s.setDuration(m_DG1p.getDuration());
    m_DG1s.setAxis(SEQ::AXIS_SLICE);
    m_DG1s.setRampUpTime(m_DG1p.getRampUpTime());
    m_DG1s.setRampDownTime(m_DG1p.getRampDownTime());
    m_DG1s.setStartTime(m_DG1p.getStartTime());

    if(! m_DG1s.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG1s.getNLSStatus());
        return false;
    }

    if(! m_DG1s.check())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG1s.getNLSStatus());
        return false;
    }

    m_DG2p.setAmplitude(dAmpl * m_Didi.getX(0));
    m_DG2p.setDuration(m_DG2p.getDuration());
    m_DG2p.setAxis(SEQ::AXIS_PHASE);
    m_DG2p.setRampUpTime(m_DG2p.getRampUpTime());
    m_DG2p.setRampDownTime(m_DG2p.getRampDownTime());

    if(! m_DG2p.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG2p.getNLSStatus());
        return false;
    }

    if(! m_DG2p.check())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG2p.getNLSStatus());
        return false;
    }

    if(m_DG2p.getTotalTime() < 300)
    {
        // CHARM 324739: Each rotation matrix must be effective for at least 300µs before the next change
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: Rotation matrix too short for DG2");
        }
        setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }

    m_DG2r.setAmplitude(dAmpl * m_Didi.getY(0));
    m_DG2r.setDuration(m_DG2p.getDuration());
    m_DG2r.setAxis(SEQ::AXIS_READOUT);
    m_DG2r.setRampUpTime(m_DG2p.getRampUpTime());
    m_DG2r.setRampDownTime(m_DG2p.getRampDownTime());
    m_DG2r.setStartTime(m_DG2p.getStartTime());

    if(! m_DG2r.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG2r.getNLSStatus());
        return false;
    }

    if(! m_DG2r.check())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG2r.getNLSStatus());
        return false;
    }

    m_DG2s.setAmplitude(dAmpl * m_Didi.getZ(0));
    m_DG2s.setDuration(m_DG2p.getDuration());
    m_DG2s.setAxis(SEQ::AXIS_SLICE);
    m_DG2s.setRampUpTime(m_DG2p.getRampUpTime());
    m_DG2s.setRampDownTime(m_DG2p.getRampDownTime());
    m_DG2s.setStartTime(m_DG2p.getStartTime());

    if(! m_DG2s.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG2s.getNLSStatus());
        return false;
    }

    if(! m_DG2s.check())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: function returned false ");
        }
        setNLSStatus(m_DG2s.getNLSStatus());
        return false;
    }

    //----------------------------
    // Calculate Export Parameters
    //----------------------------
    m_RFInfoPerRequest      = m_SBB_RF_Refoc1.getRFInfoPerRequest();
    m_RFInfoPerRequestMB    = m_SBB_RF_Refoc1.getRFInfoPerRequestMB();


    // The times must be calculated from the events itself
    // (MrProt->TE must not be used, as SBBDurationPerRequest contains the time *requested* 
    // to run the SBB for all b values. (Sorry, don't blame me, I didn't design this)
    m_lPreFill  =  m_lActualTE/2 - m_lSpinPrepTimeus                   - m_SBB_RF_Refoc1.getDurationPerRequest()/2 - m_DG1p.getTotalTime();
    m_lPostFill =  m_lActualTE/2 - m_lADCusTillEcho  - m_lStimoDelayus - m_SBB_RF_Refoc1.getDurationPerRequest()/2 - m_DG2p.getTotalTime();

    // In the next line, note that the gradient ramps do not overlap
    m_lPreEchoTimeContrib  = m_lPreFill                    + m_DG1p.getTotalTime() + m_SBB_RF_Refoc1.getDurationPerRequest()/2;
    m_lPostEchoTimeContrib = m_lPostFill + m_lStimoDelayus + m_DG2p.getTotalTime() + m_SBB_RF_Refoc1.getDurationPerRequest()/2;

    setSBBDurationPerRequest(m_lActualTE - m_lADCusTillEcho - m_lSpinPrepTimeus);

    /*
    #ifdef DEBUG
    SEQ_TRACE_ALWAYS.print("on exit: Energy SB=%f,   Energy MB=%f,   Duration=%ld,   PreTime=%ld,   PostTime=%ld, MrProt->TE[%d]=%ld",
    m_RFInfoPerRequest.getPulseEnergyWs(), m_RFInfoPerRequestMB.getPulseEnergyWs(),
    m_lPreEchoTimeContrib, m_lPostEchoTimeContrib, m_iTEArrayIndex, m_vlTE[m_iTEArrayIndex] ) ;
    #endif
    */

    m_bIsMagnetizationInverted = true;

    // calculate diffusion gradient duration (small delta) - exported to sequence and used in UI
    // parameter not really required on MARS, but may be useful, so no WIN32 restriction
    m_dDiffGradDuration_ms = static_cast<double>(m_DG1p.getDuration()) / 1000.0;

    // calculate diffusion gradient spacing (big delta) - exported to sequence and used in UI
    // parameter not really required on MARS, but may be useful, so no WIN32 restriction  
    m_dDiffGradSpacing_ms = static_cast<double>(getEventTime(DiffGrad2_Start) - getEventTime(DiffGrad1_Start)) / 1000.0;

    // finish
    return true;
}



// ===========================================================================
///	Play out the run-time events for the scan actually selected.
/**   The diffusion-specific MDH flags (REP counter, ICEProgramParams
for diffusion direction vector) will be filled out.

In MDDW and 3-scan trace mode, the diffusion gradients are in general
specified in the XYZ magnet coordinate system. In former versions of the
SBBDiffusion, the were played out using the Unity rotation matrix.
This behavior has been changed because adaptive spoiling requires
to now the sign of the diffusion gradient with respect to the imaging
gradients. In order to optimize the gradient amplitudes, we will not
just use the slice gradient matrix from now on, but we will transform
the XYZ gradient into the PRS coordinate system.

*/
// ===========================================================================
bool SBBDiffusion_Stejskal::runSBB(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
// ===========================================================================
{
    MrProtFacade protFacade(rMrProt);

    long   lPolarity          = 1;        // sign of diffusion direction
    double dBValue            = 0.;       // actual b-value
    bool   bSpoiler           = false;    // Triggers application of spoilers instead of diffusion encoding gradients

    DiffusionDirections *pDidi = NULL;

    // Initialize error return code in case of unexpected bail-outs
    setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);

    if(! isPrepared())
    {
        SEQ_TRACE_ERROR.print("ERROR: Module is not prepared");
        return false;
    }

    if(m_pADC == NULL)
    {
        SEQ_TRACE_ERROR.print("ERROR: NULL pointer for m_pADC");
        return false;
    }

    // Check for valid loop counters
    if((m_iAdjScan < 0) || (m_lBValueCounter < 0)  || (m_lDirectionCounter < 0) || (m_lDiffLoopCounter < 0))
    {
        SEQ_TRACE_ERROR.print("ERROR: invalid loop counters (probably setLoopCounters has not been called)");
        return false;
    }

    // Set diffusion encoding properties
    if(m_iAdjScan == 0)
    {
        // Imaging scans
        pDidi = &m_Didi;
    }
    else
    {
        // Adjustment scans (for dynamic distortion correction)
        pDidi = &m_AdjDidi;
    }


    if(! m_bRunOnlyOnce)
    {
        m_bRunOnlyOnce = true;
        m_bDumpBMatrix = m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_BMatrix", false);

        // Dump diffusion direction table
        pDidi->dump();

        if(m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_diffusion_vector", false))
        {
            // Dump the rotation matrix
            SEQ_TRACE_ALWAYS.print("               = ( %4.1f  %4.1f  %4.1f  ) ",
                       pSLC->getROT_MATRIX().dMat[0][0], pSLC->getROT_MATRIX().dMat[1][0], pSLC->getROT_MATRIX().dMat[2][0]);
            SEQ_TRACE_ALWAYS.print("RotationMatrix = ( %4.1f  %4.1f  %4.1f  ) ",
                       pSLC->getROT_MATRIX().dMat[0][1], pSLC->getROT_MATRIX().dMat[1][1], pSLC->getROT_MATRIX().dMat[2][1]);
            SEQ_TRACE_ALWAYS.print("               = ( %4.1f  %4.1f  %4.1f  ) ",
                       pSLC->getROT_MATRIX().dMat[0][2], pSLC->getROT_MATRIX().dMat[1][2], pSLC->getROT_MATRIX().dMat[2][2]);
        }
    }

    // -------------------------------------
    // Prepare diffusion gradient amplitudes 
    // -------------------------------------

    if(m_iAdjScan == 0)
    {
        // b-value for imaging scan (from protocol)
        dBValue = (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE) ? m_dQSpaceMaxBValue : m_vdBValues[m_lBValueCounter];
    }
    else
    {
        switch(m_eDynDistMode)
        {
            case SEQ::DYN_DISTCORR_ADJ:
            {
                // Dedicated b-value for adjustment scan:
                // - not more than 500s/mm^2 (empirical value)
                // - not more than 50% of maximum b-value
                // - do not exceed gradient amplitude used within imaging scans

                // Get maximum b value from protocol
                double dMaxRequestedBValue = (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE) ? m_dQSpaceMaxBValue : *std::max_element(m_vdBValues.begin(), m_vdBValues.end());

                // Restrict to 50% of maximum b-value or 500s/mm^2 (whatever is smaller)
                dBValue = std::min<double>(500., dMaxRequestedBValue / 2.);

                // The calculation of m_dMaxPossibleBValue (used for amplitude scaling below)
                // considers the maximum norm of the diffusion directions. Since the norm
                // of the adjustment scan direction will be different, the following scaling
                // is required in order to acquire the desired b-value.

                dBValue *= m_Didi.getGreatestNorm() * m_Didi.getGreatestNorm() / (m_AdjDidi.getGreatestNorm() * m_AdjDidi.getGreatestNorm());

                // Restrict to maximum single axis gradient amplitude used within imaging scans
                double dMaxAmpl = m_Didi.getGreatestComponent() * m_dAmpl * sqrt(dMaxRequestedBValue / m_dMaxPossibleBValue);

                if(dMaxAmpl < m_AdjDidi.getGreatestComponent() * m_dAmpl * sqrt(dBValue / m_dMaxPossibleBValue))
                {
                    dBValue = dMaxAmpl * dMaxAmpl / (m_dAmpl * m_dAmpl * m_AdjDidi.getGreatestComponent() * m_AdjDidi.getGreatestComponent()) * m_dMaxPossibleBValue;
                }
            }
                break;
            case SEQ::DYN_DISTCORR_DIRECT:
            {
                // Dedicated b-value for undistorted reference scan
                // - not more than 50s/mm^2 (empirical value)
                // - not more than minimum b-value
                // - do not exceed gradient amplitude used within imaging scans

                // Note: the first diffusion gradient direction specified within m_AdjDidi
                //       will be used for this reference scan. The calculations below are
                //       meaningful only if the first direction owns greatest norm and component.

                // Get minimum b value (smaller than default value) from protocol
                double dMinRequestedBValue = (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE) ? 0. : *std::min_element(m_vdBValues.begin(), m_vdBValues.end());

                // Restrict to minimum b-value or 50s/mm^2 (whatever is smaller)
                dBValue = std::min<double>(50., dMinRequestedBValue);

                // The calculation of m_dMaxPossibleBValue (used for amplitude scaling below)
                // considers the maximum norm of the diffusion directions. Since the norm
                // of the adjustment scan direction will be different, the following scaling
                // is required in order to acquire the desired b-value.

                dBValue *= m_Didi.getGreatestNorm() * m_Didi.getGreatestNorm() / (m_AdjDidi.getGreatestNorm() * m_AdjDidi.getGreatestNorm());

                // Restrict to maximum single axis gradient amplitude used within imaging scans
                double dMaxAmpl = m_Didi.getGreatestComponent() * m_dAmpl * sqrt(dMinRequestedBValue / m_dMaxPossibleBValue);

                if(dMaxAmpl < m_AdjDidi.getGreatestComponent() * m_dAmpl * sqrt(dBValue / m_dMaxPossibleBValue))
                {
                    dBValue = dMaxAmpl * dMaxAmpl / (m_dAmpl * m_dAmpl * m_AdjDidi.getGreatestComponent() * m_AdjDidi.getGreatestComponent()) * m_dMaxPossibleBValue;
                }
            }
                break;
            default:
            case SEQ::DYN_DISTCORR_NONE:
                SEQ_TRACE_ERROR.print("ERROR: b-value for adjustment scans undefined");
                return false;
        }
    }

    // scale as  b \propto g^2
    double dAmpl = m_dAmpl * sqrt(dBValue/m_dMaxPossibleBValue);


    double dAmpl_p = 0.0;
    double dAmpl_r = 0.0;
    double dAmpl_s = 0.0;

#ifdef QUIETDWI
    if(!prepMomentOffset(m_lPreFill / 2)){
        return false;
    }
#endif // QUIETDWI

    /** In contrast to older versions of this module, the diffusion gradients are
    always applied in the PRS coordinate system. However, if the gradient
    directions are specified in the magnet coordinate system, the inverse
    rotation matrix is applied to the gradient in order to get the vector
    in the logical PRS coordinate system. When the gradient is then played out
    in this PRS coordinate system (i.e. after the rotation matrix
    has been applied again), the gradients are actually in the XYZ system.
    So the user will not notice this double transformation.
    */

    DidiXYZ2PRS(
        pDidi,           // Input: Diffusion direction vector information 
        pSLC,            // Input: getSliceIndex, m_sROT_MATRIX 
        m_lDirectionCounter,  // Input: Current direction vector index 
        &dAmpl_p,        // Output: Phase component of diffusion vector 
        &dAmpl_r,        // Output: Read  component of diffusion vector 
        &dAmpl_s         // Output: Slice component of diffusion vector 
        );

    if(m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_diffusion_vector", false))
    {
        SEQ_TRACE_ALWAYS.print("Vector[%ld] in Didi: ( %f / %f / %f )", m_lDirectionCounter,
                   pDidi->getX(m_lDirectionCounter), pDidi->getY(m_lDirectionCounter), pDidi->getZ(m_lDirectionCounter));
        SEQ_TRACE_ALWAYS.print("Vector[%ld] in PRS:  ( %f / %f / %f )",
                   m_lDirectionCounter, dAmpl_p, dAmpl_r, dAmpl_s);
    }

    dAmpl_p *= dAmpl * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/diff_grad_ampl_X_factor", 1.0);
    dAmpl_r *= dAmpl * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/diff_grad_ampl_Y_factor", 1.0);
    dAmpl_s *= dAmpl * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/diff_grad_ampl_Z_factor", 1.0);

    // Here we have the opportunity to invert the direction of the diffusion encoding
    // gradients. With respect to residual eddy currents, this affects image appearance
    // in the following way:
    // - If the PE-component of the diffusion vector runs parallel to the PE-blips, the
    //   eddy current fields yield a stretched image
    // - If it runs anti-parallel, eddy currents yield a compressed image
    // In general, compressed high-b-value-images lead to a bright rim in calculated ADC maps,
    // while stretched high-b-images images don't yield pronounced ADC artefacts. Thus, a 
    // parallel orientation of the two gradient directions might be preferable and is enforced
    // here. Be aware that for other tensor data (e.g. FA maps), there is no benefit of this
    // selection (but also no drawback).
    //
    // Exceptions:
    // 1. In FREE and QSPACE mode, don't change the dedicated orientations.
    // 2. For eddy current adjustment scans, don't change the predefined orientations.
    if((dAmpl_p          >  0.)
       && (m_eDiffusionMode != SEQ::DIFFMODE_FREE)
       && (m_eDiffusionMode != SEQ::DIFFMODE_QSPACE)
       && (m_iAdjScan       == 0))
    {
        dAmpl_p *= -1.;
        dAmpl_r *= -1.;
        dAmpl_s *= -1.;

        lPolarity = -1;
    }

    m_DG1p.prepAmplitude(dAmpl_p);
    m_DG1r.prepAmplitude(dAmpl_r);
    m_DG1s.prepAmplitude(dAmpl_s);
    m_DG2p.prepAmplitude(dAmpl_p);
    m_DG2r.prepAmplitude(dAmpl_r);
    m_DG2s.prepAmplitude(dAmpl_s);

    // --------------------------------------------------
    // Check whether spoiler gradients have to be applied
    // --------------------------------------------------

    // For small b-values, the implicit spoiling of the diffusion gradients will
    // not be sufficient to spoil any undesired coherence pathways. In that
    // case, spoiler gradients will be applied instead of diffusion gradients.

    // Spoiler gradients are also applied if diffusion gradients are
    // deliberately disabled (e.g. for iPAT reference scans)

    if(!m_bDiffusionGradientsEnabled || isDiffusionMomentSmallerThanSpoilerMoment())
    {
        bSpoiler = true;
        dAmpl    = 0.;
        dAmpl_p  = 0.;
        dAmpl_r  = 0.;
        dAmpl_s  = 0.;

        m_DG1p.prepAmplitude(0.);
        m_DG1r.prepAmplitude(0.);
        m_DG1s.prepAmplitude(0.);
        m_DG2p.prepAmplitude(0.);
        m_DG2r.prepAmplitude(0.);
        m_DG2s.prepAmplitude(0.);
    }
    else
    {
        bSpoiler = false;
    }

    // store the scale factors between used diffusion amplitude and the max diffusion amplitude
    if (m_dMaxAmpl > 0)
    {
        setvdScaleFactorinRun(dAmpl_p / m_dMaxAmpl, dAmpl_r / m_dMaxAmpl, dAmpl_s / m_dMaxAmpl);
    }

    // ---------------------------------------------------------------------------
    // Calculate b matrix and prepare header information
    // ---------------------------------------------------------------------------

    calcBMatrix(bSpoiler);

    if(m_bDumpBMatrix)
    {
        SEQ_TRACE_ALWAYS.print("               ( %5.2f )|   ( %5.1f )|      ( %7.1f   %7.1f  %7.1f )",
                   pDidi->getX(m_lDirectionCounter), dAmpl_p, m_dBxx, m_dBxy, m_dBxz);
        SEQ_TRACE_ALWAYS.print("g[%2ld]=%6.2f * ( %5.2f )| = ( %5.1f )|    b=( %7.1f   %7.1f  %7.1f )",
                   m_lDirectionCounter, dAmpl, pDidi->getY(m_lDirectionCounter), dAmpl_r, m_dBxy, m_dByy, m_dByz);
        SEQ_TRACE_ALWAYS.print("-              ( %5.2f )|D  ( %5.1f )|P   = ( %7.1f   %7.1f  %7.1f )",
                   pDidi->getZ(m_lDirectionCounter), dAmpl_s, m_dBxz, m_dByz, m_dBzz);
        SEQ_TRACE_ALWAYS.print("tr b[%ld]=%6.1f", m_lDirectionCounter, m_dBxx + m_dByy + m_dBzz);
    }

    /*
    double maxb = std::max(m_dBxx, std::max(m_dByy, m_dBzz));
    double dirx = m_dBxx/maxb * sign(m_dBxy) * sign(m_dBxz);
    double diry = m_dByy/maxb * sign(m_dBxy) * sign(m_dByz);
    double dirz = m_dBzz/maxb * sign(m_dBxz) * sign(m_dByz);
    SEQ_TRACE_ALWAYS.print("Estimated Directions %f / %f / %f", dirx, diry, dirz);
    */

    if(m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_gradient_moments", false))
    {
        SEQ_TRACE_ALWAYS.print("Moments X:  %8.0f  %8.0f",
                   dAmpl_p * static_cast<double>(m_DG1p.getDuration()), dAmpl_p * static_cast<double>(m_DG2p.getDuration()));
        SEQ_TRACE_ALWAYS.print("Moments Y:  %8.0f  %8.0f",
                   dAmpl_r * static_cast<double>(m_DG1r.getDuration()), dAmpl_r * static_cast<double>(m_DG2r.getDuration()));
        SEQ_TRACE_ALWAYS.print("Moments Z:  %8.0f  %8.0f",
                   dAmpl_s * static_cast<double>(m_DG1s.getDuration()), dAmpl_s * static_cast<double>(m_DG2s.getDuration()));
    }

    NLS_STATUS lStatus = CalculateDicomHeaderInformation
        (
        rMrProt.sliceSeries()[static_cast<int32_t>(pSLC->getSliceIndex())], // Input: sliceSeries()
        pSLC,            // Input: getSliceIndex, getROT_MATRIX()
        pDidi,           // Input: Diffusion direction vector information
        dAmpl,           // Input: Diffusion amplitude information
        lPolarity,       // Input: Direction sign
        m_lBValueCounter,     // Input: Current bvalue index 
        m_lDirectionCounter,  // Input: Current direction vector index 
        m_pADC             // Output: Mdh.setIceProgramPara 
        );

    if(NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        return setNLSStatus(lStatus);
    }


    // --------------------------
    // Now play out the events...
    // --------------------------

    // The PreFill time is inserted asymmetrically in order to maximize the b-value:
    // Excitation - Diffusion - Fill Time - Refocussing - Diffusion

    // Event table #1: diffusion gradient 1 and fill time
    // --------------------------------------------------
    fRTEBInit(pSLC->getROT_MATRIX());

    if(bSpoiler)
    {
        fRTEI(m_DG1p.getTotalTime() + m_lPreFill - m_DSp1.getTotalTime(), 0, 0, 0, &m_DSp1, &m_DSr1, &m_DSs1, 0);
    }
    else
    {
        runGradient(&m_DG1p);
        runGradient(&m_DG1r);
        runGradient(&m_DG1s);
    }

#ifdef QUIETDWI
    fRTEI(m_DG1p.getDuration() + m_lPreFill - m_DGoffp.getTotalTime(), 0, 0, 0, &m_DGoffp, 0, 0, 0);
    fRTEI(m_DG1r.getDuration() + m_lPreFill - m_DGoffr.getTotalTime(), 0, 0, 0, 0, &m_DGoffr, 0, 0);
    fRTEI(m_DG1s.getDuration() + m_lPreFill - m_DGoffs.getTotalTime(), 0, 0, 0, 0, 0, &m_DGoffs, 0);
#endif // QUIETDWI

    fRTEI(m_DG1p.getTotalTime() + m_lPreFill, 0, 0, 0, 0, 0, 0, 0);

    if(SeqUT.isUnitTestActive())
    {
        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'D', 0, pSLC->getSliceIndex(), 0, 0);
    }

    if(setNLSStatus(fRTEBFinish()))
    {
        SEQ_TRACE_ERROR.print("ERROR: Execution of event block 1 failed.");
        return false;
    }


    // Event table #2: inversion pulse
    // -------------------------------------

    m_SBB_RF_Refoc1.setRunMode(m_eSliceAccelRFRunMode);

    // play out SBB

    if(!m_SBB_RF_Refoc1.run(rMrProt, rSeqLim, rSeqExpo, pSLC))
    {
        SEQ_TRACE_ERROR.print("ERROR: m_SBBMultibandRFRefoc1.run failed.");
        return false;
    }


    // Event table #3: diffusion gradients / stimulation delay
    // -------------------------------------------------------

    fRTEBInit(pSLC->getROT_MATRIX());

    if(bSpoiler)
    {
        fRTEI(0, 0, 0, 0, &m_DSp1, &m_DSr1, &m_DSs1, 0);
    }
    else
    {
        runGradient(&m_DG2p);
        runGradient(&m_DG2r);
        runGradient(&m_DG2s);
    }
    fRTEI(m_DG2p.getTotalTime() + m_lStimoDelayus, 0, 0, 0, 0, 0, 0, 0);

    if(SeqUT.isUnitTestActive())
    {
        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'D', 0, pSLC->getSliceIndex(), 0, 0);
    }

    if(setNLSStatus(fRTEBFinish()))
    {
        SEQ_TRACE_ERROR.print("ERROR: Execution of event block 3 failed.");
        return false;
    }

    // ----------------------------------
    // TE fill time after the second half
    // ----------------------------------
    if(m_lPostFill)
    {
        fRTEBInit(pSLC->getROT_MATRIX());
        fRTEI(m_lPostFill, 0, 0, 0, 0, 0, 0, 0);

        if(SeqUT.isUnitTestActive())
        {
            mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'R', 0, pSLC->getSliceIndex(), 0, 0);
        }

        if(setNLSStatus(fRTEBFinish()))
        {
            SEQ_TRACE_ERROR.print("ERROR: Execution of TE fill time failed.");
            return false;
        }
    }

    setNLSStatus (MRI_SBB_SBB_NORMAL);
    return true;
}

// ===========================================================================
/// Calculate the timing of the diffusion gradients in oder achieve the highest b value.
/**
\pre The member variables used for configuration of the SBB must have been set.
*/

bool SBBDiffusion_Stejskal::prepTiming(MrProt &rMrProt, SeqLim & rSeqLim, SeqExpo &rSeqExpo, long lActualTE)
{
    MrProtFacade protFacade(rMrProt);

    // Flag: Print a message before every return statement
    bool bDebugReturn = ((!(rSeqLim.isContextPrepForBinarySearch() || rSeqLim.isContextPrepForScanTimeCalculation())) && !m_bResolve);

    long lt            = 0;

    // -------------------------------------------------
    // first we have to prepare the RF SBBs as their timing  
    // does not depend on the diffusion gradients
    // -------------------------------------------------
    if(!prepRF(rMrProt, rSeqLim, rSeqExpo))
        return false;

    //  ----------------------------------------------------
    /** Set up timing and calculate maximum possible b value */
    //  ----------------------------------------------------
    // Gradient amplitude for dummy preparation - real preparations will take place below
    double dAmpl = 5.;

    // These durations have to be set first: they are required
    // in order to calculate lt
    m_DG1p.setRampTimes(m_lRampTime);
    m_DG2p.setRampTimes(m_lRampTime);

    // lt = "delta" = Total time of one diffusion gradient
    lt = fSDSRoundUpGRT(lActualTE/2 - m_SBB_RF_Refoc1.getDurationPerRequest()/2 - std::max(m_lSpinPrepTimeus, m_lADCusTillEcho + m_lStimoDelayus));

    // Do some cross checks
    // --------------------

    // lt holds at least: 
    // - spoiler, slicegrad up, half RF pulse OR 
    // - diffgrad ramp up & down, slicegrad up, half RF pulse

    if((lt <= 0)                                               ||
       (getlSpoilerTotalTime() + m_SBB_RF_Refoc1.getDurationPerRequest()/2 > lt) ||
       (2 * m_lRampTime        + m_SBB_RF_Refoc1.getDurationPerRequest()/2 > lt))
    {
        setNLSStatus (MRI_SBB_SBB_NEGATIV_TEFILL);

        if(bDebugReturn)
        {
            SEQ_TRACE_ERROR.print("ERROR: lt = %ldus is too short or negative (spoiler time: %ldus)", lt, getlSpoilerTotalTime());
        }

        return false;
    }

    // -------------------------------------------------
    // first event block (diffusion gradient)
    // -------------------------------------------------

    m_DG1p.setAxis(SEQ::AXIS_PHASE);
    m_DG1p.setAmplitude(dAmpl);
    m_DG1p.setStartTime(0);
    m_DG1p.setDuration(lt - m_lRampTime);       // Remember: gradient pulse duration = RampUp + FlatTop


    // ---------------------------------------------------------------------------
    // 3rd event block
    // ---------------------------------------------------------------------------

    m_DG2p.setAxis(SEQ::AXIS_PHASE);
    m_DG2p.setAmplitude(dAmpl);
    m_DG2p.setStartTime(0);
    m_DG2p.setDuration(lt - m_lRampTime);


    // Calculate fill times (required for b-value calculation)
    m_lPreFill  =  lActualTE/2 - m_lSpinPrepTimeus                   - m_SBB_RF_Refoc1.getDurationPerRequest()/2 - m_DG1p.getTotalTime();
    m_lPostFill =  lActualTE/2 - m_lADCusTillEcho  - m_lStimoDelayus - m_SBB_RF_Refoc1.getDurationPerRequest()/2 - m_DG2p.getTotalTime();

    // prepare compensation gradients
    if (getbCompensationEnable())
    {
        long lDiffGradDuration = m_DG1p.getDuration();
        long lDiffGradSpacing  = getEventTime(DiffGrad2_Start) - getEventTime(DiffGrad1_Start);
        long lDiffToCompGrad   = m_lStimoDelayus + m_lPostFill + getlPlugInToCompGradTime();

        // calculate the duration of compensation gradients based on the max amplitude of diffusion gradients
        m_CompGrad.calcDiffCompGrad(m_dMaxAmpl, lDiffGradDuration, lDiffGradSpacing, lDiffToCompGrad);

        // only valid when calculated duration is more than gradient raster time
        if (m_CompGrad.getbValid())
        {
            if (!m_CompGrad.prep(rMrProt, rSeqLim, rSeqExpo))
            {
                if (bDebugReturn)
                {
                    SEQ_TRACE_DEBUG.print("m_CompGrad.prep failed.");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
            }
        }
        else
        {
            m_CompGrad.resetPrepared();
        }
    }  

    // Calculate maximum diffusion encoding gradient amplitude that is compatible with the current timing.
    // This yields the maximum possible b-value. All necessary information is read from member variables.
    if(!CalcMaximumAmplitude(rMrProt, dAmpl))
    {
        if(bDebugReturn)
        {
            SEQ_TRACE_ERROR.print("ERROR: CalcMaximumAmplitude failed.");
        }
        setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }
    SEQ_TRACE_DEBUG.print("CalcMaximumAmplitude: dAmpl: %f mT/m", dAmpl);

    // Store calculated amplitude - the required TR increment is calculated in prepFinal
    setMaxAmplitude(dAmpl);             // Update m_dAmpl

    m_DG1p.setAmplitude(dAmpl);
    m_DG2p.setAmplitude(dAmpl);

    // ------------------------------------------------------------
    /** Calculate the b value that will be obtained with this timing
    and store it in SBBDiffusion_Base::m_dMaxPossibleBValue.
    Remember that all gradients have been prepared with maximum amplitude.
    Only spoilers on the P-axis are considered in this b value estimation.
    */
    // ------------------------------------------------------------

    m_dMaxPossibleBValue = calcBValue();

    return true;
}




// ===========================================================================
double SBBDiffusion_Stejskal::calcBValue(void)
// ===========================================================================

{
    // For calculation, the ramps and durations of dg1 and dg2 are used
    // as well as the relative difference of the start times between dg2 and dg3.

    // This function is called within protocol preparation, thus no actual
    // diffusion direction is defined. At least the diffusion encoding pulses
    // in phase encoding direction (m_DG1p and m_DG2p, see calcTiming) are
    // prepared with maximum amplitude. For the calculation of the maximum
    // b-value, it is sufficient to look at the corresponding diagonal element
    // (Bxx) and multiply this b-value with the square (since the
    // b-value scales with the square of the gradient amplitude) of the 
    // maximum norm of the diffusion direction.
    //
    // e.g. orthogonal diffusion directions: norm = 1
    // e.g. room diagonals (3-scan-trace)  : norm = sqrt(3)

    double dBValue = 0.;

    // Set gamma of actual nucleus
    m_theBMatrix.setGamma(m_dGamma);

    // Prepare for b-value calculation (set up diffusion encoding events):
    //  b-value only, no spoilers
    PrepBMatrix(true, false);

    // Calculate b-value
    if(!m_theBMatrix.bCalcBValue(dBValue))
    {
        SEQ_TRACE_WARN.print("WARNING: Non-null moment");
    }

    return (dBValue * m_Didi.getGreatestNorm() * m_Didi.getGreatestNorm());

}


// ===========================================================================
/// Calculate the b matrix which will be obtained by the current timing.
/** Each component is calculated individually by summing up the moments
generated by the gradients called on the axis.

\pre For the calculation, it is supposed that all diffusion gradients have
been prepared completely.

\return The result will be returned in the member variables
SBBDiffusion_Base::m_dBxx, m_dByy, m_dBzz, m_dBxy, m_dBxz, and m_dByz.

In the DEBUG version, there is a check that the gradient moments are
zero at the end of the kernel. If this check fails, a warning message
will be generated.
*/
void SBBDiffusion_Stejskal::calcBMatrix(bool bApplySpoiler)
// ===========================================================================
{
    // Set gamma of actual nucleus
    m_theBMatrix.setGamma(m_dGamma);

    // Prepare for b-matrix calcuation (set up diffusion encoding events):
    //  complete b-matrix, spoilers on demand
    PrepBMatrix(false, bApplySpoiler);

    // Calculate b-matrix
    if(!m_theBMatrix.bCalcBMatrix(m_dBxx, m_dBxy, m_dBxz, m_dByy, m_dByz, m_dBzz))
    {
        SEQ_TRACE_WARN.print("WARNING: Non-null moment");
    }

    // Calculate b-value
    m_theBMatrix.bCalcBValue(m_dBValue);
}

void SBBDiffusion_Stejskal::calcBMatrix(void)
// ===========================================================================
{
    // Set gamma of actual nucleus
    m_theBMatrix.setGamma(m_dGamma);

    // Prepare for b-matrix calcuation (set up diffusion encoding events):
    //  complete b-matrix, no spoilers
    PrepBMatrix(false, false);

    // Calculate b-matrix with default gamma
    if(!m_theBMatrix.bCalcBMatrix(m_dBxx, m_dBxy, m_dBxz, m_dByy, m_dByz, m_dBzz))
    {
        SEQ_TRACE_WARN.print("WARNING: Non-null moment");
    }


    // Calculate b-value
    m_theBMatrix.bCalcBValue(m_dBValue);
}

bool SBBDiffusion_Stejskal::PrepBMatrix(bool bBValueOnly, bool bApplySpoiler)
// ===========================================================================
{
    // Initialize b-matrix calculation
    m_theBMatrix.Reset();

    // Fill in basic diffusion encoding events
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DG1p);       // PE-axis dephase
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DG2p);       // PE-axis rephase

    m_theBMatrix.bAddEvent(getEventTime(RefocRF1_Center), m_SBB_RF_Refoc1.getRFPulsePointer());           // First refocussing RF

    if(!bBValueOnly)
    {
        // Full b-matrix calculation
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DG1r);   // RO-axis dephase
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DG1s);   // SL-axis dephase
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DG2r);   // RO-axis rephase
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DG2s);   // SL-axis rephase

        //add versed slice gradient for accurate b-matrix calculation
        sGRAD_PULSE_TRAP* SliceGrad1 = m_SBB_RF_Refoc1.getGradientPointer(0);
        long StartTimeSliceGrad1 = getEventTime(SliceGrad1_Start);

        m_theBMatrix.bAddEvent(StartTimeSliceGrad1, *SliceGrad1);         // SL-axis first verse slice selection gradient (or the normal slice selective gradient)
        if(m_SBB_RF_Refoc1.isVerse()) // is Verse
        {
            sGRAD_PULSE_TRAP* SliceGrad2 = m_SBB_RF_Refoc1.getGradientPointer(1);
            sGRAD_PULSE_TRAP* SliceGrad3 = m_SBB_RF_Refoc1.getGradientPointer(2);
            long StartTimeSliceGrad2 = StartTimeSliceGrad1 + SliceGrad1->getDuration();
            long StartTimeSliceGrad3 = StartTimeSliceGrad2 + SliceGrad2->getDuration();

            m_theBMatrix.bAddEvent(StartTimeSliceGrad2, *SliceGrad2);         // SL-axis second verse slice selection gradient
            m_theBMatrix.bAddEvent(StartTimeSliceGrad3, *SliceGrad3);         // SL-axis third  verse slice selection gradient
        }


        if(bApplySpoiler)
        {
            // Note: if spoilers are applied instead of diffusion gradients, 
            // diffusion gradient amplitudes are set to zero (so it is still
            // ok to include them in the b-matrix calculation)
            m_theBMatrix.bAddEvent(getEventTime(SpoilGrad1_Start), m_DSp1);   // PE-axis spoiler 1
            m_theBMatrix.bAddEvent(getEventTime(SpoilGrad1_Start), m_DSr1);   // RO-axis spoiler 1
            m_theBMatrix.bAddEvent(getEventTime(SpoilGrad1_Start), m_DSs1);   // SL-axis spoiler 1
            m_theBMatrix.bAddEvent(getEventTime(SpoilGrad2_Start), m_DSp1);   // PE-axis spoiler 2
            m_theBMatrix.bAddEvent(getEventTime(SpoilGrad2_Start), m_DSr1);   // RO-axis spoiler 2
            m_theBMatrix.bAddEvent(getEventTime(SpoilGrad2_Start), m_DSs1);   // SL-axis spoiler 2
        }
    }

    if(m_theBMatrix.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}

bool SBBDiffusion_Stejskal::prepGPALoadDiff(double dAmplitudeX)
{
    // Clear
    m_sBalanceDiff.Reset();

    // Set typical frequency content of diffusion module [Hz]
    m_sBalanceDiff.bSetFrequency(100);

    // ------------------------------------------------------------------------
    // Prepare diffusion gradient events for balance calculations.
    // Consider GPABALANCE_X_AXIS only
    // See also the b-matrix calculation above (timing has to match!)

    // Add diffusion encoding gradients and store ID's for scaling purposes
    m_lID1X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad1_Start), dAmplitudeX, m_DG1p.getRampUpTime(), m_DG1p.getRampDownTime(), m_DG1p.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_lID2X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad2_Start), dAmplitudeX, m_DG2p.getRampUpTime(), m_DG2p.getRampDownTime(), m_DG2p.getFlatTopTime(), GPABALANCE_X_AXIS);

    // Check error status
    if(m_sBalanceDiff.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}

bool SBBDiffusion_Stejskal::prepGPALoadDiff(double dAmplitudeX, double dAmplitudeY, double dAmplitudeZ)
{
    // Clear
    m_sBalanceDiff.Reset();

    // Set typical frequency content of diffusion module [Hz]
    m_sBalanceDiff.bSetFrequency(100);

    // ------------------------------------------------------------------------
    // Prepare diffusion gradient events for balance calculations.
    // Consider GPABALANCE_X_AXIS only
    // See also the b-matrix calculation above (timing has to match!)

    // Add diffusion encoding gradients and store ID's for scaling purposes
    m_lID1X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad1_Start), dAmplitudeX, m_DG1p.getRampUpTime(), m_DG1p.getRampDownTime(), m_DG1p.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_lID2X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad2_Start), dAmplitudeX, m_DG2p.getRampUpTime(), m_DG2p.getRampDownTime(), m_DG2p.getFlatTopTime(), GPABALANCE_X_AXIS);

    m_lID1Y = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad1_Start), dAmplitudeY, m_DG1p.getRampUpTime(), m_DG1p.getRampDownTime(), m_DG1p.getFlatTopTime(), GPABALANCE_Y_AXIS);
    m_lID2Y = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad2_Start), dAmplitudeY, m_DG2p.getRampUpTime(), m_DG2p.getRampDownTime(), m_DG2p.getFlatTopTime(), GPABALANCE_Y_AXIS);

    m_lID1Z = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad1_Start), dAmplitudeZ, m_DG1p.getRampUpTime(), m_DG1p.getRampDownTime(), m_DG1p.getFlatTopTime(), GPABALANCE_Z_AXIS);
    m_lID2Z = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad2_Start), dAmplitudeZ, m_DG2p.getRampUpTime(), m_DG2p.getRampDownTime(), m_DG2p.getFlatTopTime(), GPABALANCE_Z_AXIS);

    // Check error status
    if(m_sBalanceDiff.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}

bool SBBDiffusion_Stejskal::scaleGPALoadDiff(double dScaleX)
{
    const int iNumberOfEvents       = 2;

    const std::array<long, iNumberOfEvents>   alID    = {m_lID1X, m_lID2X};
    const std::array<double, iNumberOfEvents> adScale = {dScaleX, dScaleX};

    // ------------------------------------------------------------------------
    // Scale diffusion encoding events on x-axis
    if (!m_sBalanceDiff.bScaleEvent(alID, adScale))
    {
        SEQ_TRACE_ERROR << "error scaling gradient events.";
        return false;
    }

    return true;
}

bool SBBDiffusion_Stejskal::scaleGPALoadDiff(double dScaleX, double dScaleY, double dScaleZ)
{
    const int iNumberOfEvents       = 6;

    const std::array<long, iNumberOfEvents>   alID    = {m_lID1X, m_lID2X, m_lID1Y, m_lID2Y, m_lID1Z, m_lID2Z};
    const std::array<double, iNumberOfEvents> adScale = {dScaleX, dScaleX, dScaleY, dScaleY, dScaleZ, dScaleZ};

    // ------------------------------------------------------------------------
    // Scale diffusion encoding events
    if (!m_sBalanceDiff.bScaleEvent(alID, adScale))
    {
        SEQ_TRACE_ERROR << "error scaling gradient events.";
        return false;
    }

    return true;
}

long SBBDiffusion_Stejskal::getEventTime(EnumEventTime eEvent)
{
    // Start of first diffusion encoding gradient is zero
    long lTime0 = 0;
    // Start time of the second diffusion encoding gradient
    long lTime1 = m_DG1p.getTotalTime() + m_lPreFill + m_SBB_RF_Refoc1.getDurationPerRequest();
    // Refocussing time
    long lTime2 = m_DG1p.getTotalTime() + m_lPreFill + m_SBB_RF_Refoc1.getDurationPerRequest()/2;
    // Start time of first spoiler (applied if b-value is below threshold)
    long lTime3 = m_DG1p.getTotalTime() + m_lPreFill - m_DSp1.getTotalTime();
    // Start time of second spoiler (applied if b-value is below threshold)
    long lTime4 = lTime1;
    // Start time of slice selection gradient
    long lTime5 = m_DG1p.getTotalTime() + m_lPreFill;

    switch(eEvent)
    {
        case DiffGrad1_Start:
            return lTime0;
        case DiffGrad2_Start:
            return lTime1;
        case SpoilGrad1_Start:
            return lTime3;
        case SpoilGrad2_Start:
            return lTime4;
        case SliceGrad1_Start:
            return lTime5;
        case RefocRF1_Center:
            return lTime2;
        default:
            SEQ_TRACE_ERROR.print("unknown event enumerator: %i", eEvent);
            return -1;
    }
}


bool SBBDiffusion_Stejskal::prepRF(MrProt & rMrProt, SeqLim &rSeqLim, SeqExpo & rSeqExpo)
{
    MrProtFacade protFacade(rMrProt);


    // -------------------------------------------
    // 2nd event block (slice-selective RF pulse)
    // -------------------------------------------

    // get properties of the pulse and pass to pulse object
    sRFPulseProperties myRFPulseProperties = m_RFPulseLibrary.getPulsePropertiesRefocusing(rMrProt);

    m_DRF1.setFamilyName(myRFPulseProperties.sFamilyName.c_str());
    m_DRF1.setTypeRefocussing();
    m_DRF1.setFlipAngle(180.0);
    m_DRF1.setInitialPhase(0.0);
    m_DRF1.setThickness(m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_DRF1.setDuration(myRFPulseProperties.lDuration_us);

    // Apply gradient reversal to refocusing pulse (just for ep2d_diff but not for RESOLVE)
    if(protFacade.isGradientReversalDiffusion() && !m_bResolve)
        m_DRF1.setRequiredGSPolarity(-1.0);
    else
        m_DRF1.setRequiredGSPolarity(+1.0);


    if((m_dPreparedSlcThk != m_dRFPulseThicknessFactor * m_dSliceThickness) || (!rSeqLim.isContextPrepForBinarySearch()))
    {
        if(! m_DRF1.prepExternal(rMrProt, rSeqExpo))
        {
            if(!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("ERROR: DRF1.prepExternal failed.");
            }
            setNLSStatus(m_DRF1.getNLSStatus());
            return false;
        }
        m_dPreparedSlcThk = m_dRFPulseThicknessFactor * m_dSliceThickness;
    }

    m_SBB_RF_Refoc1.setIdent("SBB_MB_Refoc1");
    if(!prep_SBBMultibandRF(rMrProt, rSeqLim, rSeqExpo, m_SBB_RF_Refoc1, m_DRF1))
    {
        if (!rSeqLim.isContextPrepForBinarySearch())
        {
            SEQ_TRACE_ERROR.print("ERROR: prep_SBBMultibandRF() failed.");
        }
        return false;
    }

#ifdef WIN32
    // If dRFPulseThickness is deliberately set unequal 1.0, SeqUT needs to know this.
    SeqUT.setRFThicknessInfo(&m_DRF1, m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_SBB_RF_Refoc1.setRunMode(SINGLE_BAND);
    SeqUT.setRFThicknessInfo(m_SBB_RF_Refoc1.getRFPulsePointer(), m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_SBB_RF_Refoc1.setRunMode(MULTI_BAND);
    SeqUT.setRFThicknessInfo(m_SBB_RF_Refoc1.getRFPulsePointer(), m_dRFPulseThicknessFactor * m_dSliceThickness);

    m_SBB_RF_Refoc1.setRunMode(m_eSliceAccelRFRunMode);
#endif

    return true;
}

bool SEQ_NAMESPACE::SBBDiffusion_Stejskal::getSliceAdjRFInfo(const SLICEADJ::sCuboidGeometry &sSliceAdjCuboid, std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo)
{
    return m_SBB_RF_Refoc1.getSliceAdjRFInfo(sSliceAdjCuboid, vsRFInfo);
}

MrProtocolData::SeqExpoRFInfo SEQ_NAMESPACE::SBBDiffusion_Stejskal::getRFInfoPerRequest()
{
    return m_SBB_RF_Refoc1.getRFInfoPerRequest();
}

bool SEQ_NAMESPACE::SBBDiffusion_Stejskal::isDiffusionMomentSmallerThanSpoilerMoment()
{
    return 
    ((fabs(m_DG1p.getAmplitude()) * static_cast<double>(m_DG1p.getDuration()) < m_dRefSpoilMoment * 1000.) &&
     (fabs(m_DG1r.getAmplitude()) * static_cast<double>(m_DG1r.getDuration()) < m_dRefSpoilMoment * 1000.) &&
     (fabs(m_DG1s.getAmplitude()) * static_cast<double>(m_DG1s.getDuration()) < m_dRefSpoilMoment * 1000.));
}

long SEQ_NAMESPACE::SBBDiffusion_Stejskal::getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch)
{
    long lSmallesPossibleBValue = m_lBValueInc_Limit;
    bool bResult = false;


    // step through all smaller b values and check if spoiler is used
    for(long lBValue = lSmallesPossibleBValue - IVIM_B_VALUE_INCREMENT; lBValue > 0; lBValue -= IVIM_B_VALUE_INCREMENT)
    {
        bResult = prepFinal((double)lBValue, bIsContextPrepForBinarySearch);
        if((!isDiffusionMomentSmallerThanSpoilerMoment()) && bResult)
        {
            lSmallesPossibleBValue = lBValue;
        }
        else 
        {
            return lSmallesPossibleBValue;
        }
    }

    return lSmallesPossibleBValue;
}

