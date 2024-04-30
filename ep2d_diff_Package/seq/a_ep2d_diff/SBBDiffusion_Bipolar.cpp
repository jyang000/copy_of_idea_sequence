//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 2010  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\SBBDiffusion_Bipolar.cpp
//
//      Author: PLM AW NEUR
//
//        Lang: C++
//
//     Descrip: Implementation of the class Diffusion_Bipolar.
//
//     Classes: Diffusion_Bipolar



/**
***************************************************************************

\changed     1-Sep-2002; M.Zwanger; 4a21a
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- ancient monster file split up in several files
- class 'Diffusion_Stejskal' added
- comments changed to doxygen format (http://www.doxygen.org)
- some .h includes removed
- removed member m_Images
- m_lSpoilerDuration renamed in m_lSpoilerTotalTime
- m_Norm removed
- Tensor directions controlled by class DiffusionDirections
- the spoiler duration is taken into account for b value calculation
by Diffusion_Bipolar::calcBValue
- Spoiler gradients only for b < SPOILER_THRESHOLD
- Default lambda = 30 ms
- Initial phase of m_DRF2 changed to 180 deg (Thanks, Juergen :-)
- shorter RF pulse for 3T (Thanks, Larry :-)

\changed     21-Oct-2003; M.Zwanger; 4b13a
\requirement CHARM 324742, CHARM 324739,
N4_elh_DTI_Measurement_b_matrix_DICOM_header
\description
- CHARM 324742: Maxwell crossterm compensation for AC44
- 'm_bDiffusionGradientsEnabled' support added
- compiler define 'EP2D_DIFF_VA15A' no longer supported
- include paths adapted for VA25A archive
- getsWiPMemBlock().getadFree() bug fixed
- CHARM 324739: Each rotation matrix must be effective for at least 300us
- DG1 renamed in DG1p etc.
- class calcBMatrix() added
- pi pulses have now a phase of 90 and 270 deg resp.
- calcTiming() has no longer MeasNucleus in parameter list
- Spoiler schema completely revised

\changed    07-Sep-2005; M.Zwanger; 4b13a; CHARM: n.a.
\description
- MDH FreeParameters are used to export b matrix
- Nominal b value and direction index exported in MDH
- MDH CSet counter no longer used for MDDH mode

\changed    29-Sep-2005; M.Zwanger; 4b13a; CHARM: n.a.
\description
- The nominal b value for the image text is taken from the b value table;
for FREE mode the real value is used.
- The REP MDH counter runs to (diff_scans * repetitions) to handle more
than 1 average.

\changed    12-Oct-2005; M.Zwanger; 4b13a; Customer wish
\description
- Effective Amplitude based on m_Didi.getGreatestNorm()

\changed    124-Nov-2005; M.Zwanger; 4b13a; Bug fix
\description
- Apply rotation matrix to b matrix

\changed    9-Dec-2005; S. Huwer; 4b13a; CHARM: 344712
\description
- CalculateGradientsToPCS, MatMult, MatMultTrans

\changed    9-Dec-2005; S. Huwer; 4b13a; CHARM: 344712
\description
- rotMatrix transforms from GCS/DCS to PCS

\changed    13-Nov-2006; T.Feiweier; WIP
\description
- EC compensation time constant lambda configurable within ini-File

\changed    5-Sep-2007; T.Feiweier; Development
- Separate class for b-matrix calculations
- GPA-model based optimization of diffusion gradient amplitudes

***************************************************************************
*/
// MrProt
#include "MrProtSrv/Domain/MrProtData/MrProt/Application/Application.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Filter/MrFilter.h"
// MrProt

#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"

#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Bipolar.h"

#include "MrImaging/seq/SeqDebug.h"
#include "MrImagingFW/libSeqSysProp/SysProperties.h"

#include <math.h>                                          // ln
#include <assert.h>                                        // assert
#include <algorithm>                                       // std::max_element
#include <utility>										   // std::make_pair

#include "MrImaging/libSBB/libSBBmsg.h"          // SBB_ message codes
#include "MrMeasSrv/SeqIF/libRT/libRT.h"         // for RTEBInit,...
#include "MrImagingFW/libSeqUTIF/libsequt.h"               // mSEQTest
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrImaging/libSeqUtil/libSeqUtil.h"     // for runGradient, printGradient
#include "MrMeasSrv/SeqIF/libRT/sROT_MATRIX.h"   // sROT_MATRIX

// DiffusionRFPulseProperties
#include "MrImaging/seq/a_ep2d_diff/DiffusionRFPulseProperties.h"
// DiffusionRFPulseProperties

// MRProtFacade
#include "MrImaging/seq/common/MrProtFacade/MrProtFacade.h"
// MRProtFacade

// SMS support
#include "MrImaging/libSBB/SBBMultibandRF.h"
// SMS support

#include "MrImaging/seq/common/MaxwellCorrection/MaxwellCorr.h"

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif
using namespace SEQ_NAMESPACE;

// ---------------------------------------------------------------------------
// Type definitions
// ---------------------------------------------------------------------------

// Empirical factor: provide desired slice thickness for the twice refocused spin echo
const double dRFPulseEmpiricalFactor = 1.15;

// ===========================================================================
/*!
\class Diffusion_Bipolar

\brief This class implements diffusion weighting with a double spin echo and
four (bipolar) diffusion gradients.

\author Michael.Zwanger@med.siemens.de

This class is derived from the SBBDiffusion_Base class and provides a
diffusion weighting with four bipolar diffusion gradients and two refocussing
RF pulses in between.
This scheme has been suggested by Oliver Heid: "Eddy-Current-Nulled Diffusion
Weighting", Proc. ISMRM 2000, p. 799 (<A HREF="../ISMRM_2000_799_Heid_Bipolar.pdf">[pdf]</A>).

Its advantage against the standard Stejskal-Tanner scheme
(class Diffusion_Stejskal) is its ability to compensate eddy currents.

It can be used for diffusion measurements in read, slice and phase direction
as well as for the orthogonal and tensor mode.

\image html Diffusion_Bipolar.gif "Timing diagram of Diffusion_Biploar"

The timing consists of the following event blocks:
-# a fill time,
-# the first diffusion gradient lobe,
-# the refocussing RF pulse (with a slice-selection gradient)
-# the second and third diffusion gradient,
-# the second refocussing RF pulse (with a slice-selection gradient)
-# the fourth diffusion gradient, and
-# a fill time.

This diffusion module that generates diffusion-weighted images with
sensitivity in a single direction uses a "dual bipolar" design to reduce
the effects of eddy currents. In a standard "STEJSKAL-TANNER" implementation
of diffusion weighting, two equal-length and same-polarity diffusion gradient
lobes straddle the single 180-degree refocusing RF pulse.
These large amplitude diffusion gradients can generate eddy currents,
to which EPI acquisitions can be particularly sensitive. The two equal sign lobes can
be split into two bipolar pairs by the application of a second 180-degree RF pulse. The
resulting diffusion-sensitizing gradients consist of 4 lobes: the first and fourth lobes have
equal length and opposite polarity and the second and third also have equal length and
opposite polarity. Using gradient pulses of different duration or amplitude, an exact eddy
current compensation can be achieved for eddy currents with a finite half-life time.
It makes sense to chose the eddy current delay constant (referred to as "lambda")
in a way that the eddy currents are suppressed at the center position of the echo.

The ratio of the pulse durations is calculated depending on the current
protocol parameters and must meet the following conditions:
- magnetization is rephased if a 180 pulse is placed between the first and
second as well as the third and fourth gradient pulse,
- all stimulated echoes are dephased (therefore we need the spoilers),
- the spin echo will appear at 10% of the total duration of the pulse train
after the end of the SBB (neglecting the RF pulses and the gradient ramps),
- eddy currents with a decay time \f$\lambda\f$ are completely compensated.
Between the two pi pulses (or, to be more precise, between the second and
third gradient), the SBB inserts the time TE/2 for correct echo adjustment.

In older versions of the code, two different types of spoilers (i.e. spoilers
applied for all b-values and spoilers only applied for low b-values) were used.
Since the implementation of the b matrix calculation only one type of spoiler
is used which is applied in all cases.
The gradient events are named ::m_DSpoil1S, ::m_DSpoil1R and ::m_DSpoil1P
(spoiler pad for the 1st RF pulse) and ::m_DSpoil2S, ::m_DSpoil2R and
::m_DSpoil2P (for the 2nd RF pulse). (The spoiler gradients provided by
the base class are not used in this module.)


\image html Diffusion_Bipolar0.gif "Timing diagram of Diffusion_Bipolar for b=0"


**/


//   ===========================================================================
///  The constructor initializes the starting time of the diffusion gradients 
///  with 0 and the maximum possible gradient amplitudes.
SBBDiffusion_Bipolar::SBBDiffusion_Bipolar(SBBList* pSBBList):
SBBDiffusion_Base(pSBBList),
m_AddFillTime(0),
m_dSpoilerAmplitude(0),
m_dSpoilerFactor1(0),
m_dSpoilerFactor2(0),
m_bRunOnlyOnce(false),
m_bDumpBMatrix(false),
m_dPreparedSlcThk(0.0),
m_lID1X(0),
m_lID2X(0),
m_lID3X(0),
m_lID4X(0),
m_lID1Y(0),
m_lID2Y(0),
m_lID3Y(0),
m_lID4Y(0),
m_lID1Z(0),
m_lID2Z(0),
m_lID3Z(0),
m_lID4Z(0),
m_lSpinPrepTimeEnhancement(0),
m_DRF1("RefocRF1"),
m_DRF2("RefocRF2"),
m_DG1p("RTEIdentDG1p"),
m_DG1r("RTEIdentDG1r"),
m_DG1s("RTEIdentDG1s"),
m_DG2p("RTEIdentDG2p"),
m_DG2r("RTEIdentDG2r"),
m_DG2s("RTEIdentDG2s"),
m_DG3p("RTEIdentDG3p"),
m_DG3r("RTEIdentDG3r"),
m_DG3s("RTEIdentDG3s"),
m_DG4p("RTEIdentDG4p"),
m_DG4r("RTEIdentDG4r"),
m_DG4s("RTEIdentDG4s"),
m_DSpoil1P("RTEIdDSpP"),
m_DSpoil2P("RTEIdDSp2P"),
m_DSpoil1R("RTEIdDSpR"),
m_DSpoil2R("RTEIdDSp2R"),
m_DSpoil1S("RTEIdDSpS"),
m_DSpoil2S("RTEIdDSp2S")

// ===========================================================================
{
    setIdent("Diffusion_Bipolar");

    m_maDG.insert(std::pair<std::string, std::array<sGRAD_PULSE_TRAP*, 4> >("p", { &m_DG1p, &m_DG2p, &m_DG3p, &m_DG4p }));
    m_maDG.insert(std::pair<std::string, std::array<sGRAD_PULSE_TRAP*, 4> >("r", { &m_DG1r, &m_DG2r, &m_DG3r, &m_DG4r }));
    m_maDG.insert(std::pair<std::string, std::array<sGRAD_PULSE_TRAP*, 4> >("s", { &m_DG1s, &m_DG2s, &m_DG3s, &m_DG4s }));


    m_DG1p.setStartTime(0);
    m_DG2p.setStartTime(0);
    m_DG3p.setStartTime(0);
    m_DG4p.setStartTime(0);

    // Set maximum possible(!) gradient amplitudes (used during check())
    double dAmpl = SysProperties::getGradMaxAmplAbsolute();

    for (auto& elem : m_maDG)
    {
        std::for_each(std::begin(elem.second), std::end(elem.second), [=](sGRAD_PULSE_TRAP* pgrad) {pgrad->setMaxMagnitude(dAmpl); });
    }

    m_DSpoil1P.setMaxMagnitude(dAmpl);
    m_DSpoil1R.setMaxMagnitude(dAmpl);
    m_DSpoil1S.setMaxMagnitude(dAmpl);
    m_DSpoil2P.setMaxMagnitude(dAmpl);
    m_DSpoil2R.setMaxMagnitude(dAmpl);
    m_DSpoil2S.setMaxMagnitude(dAmpl);

    m_DSpoil1S.setAxis(SEQ::AXIS_SLICE);
    m_DSpoil2S.setAxis(SEQ::AXIS_SLICE);
    m_DSpoil1P.setAxis(SEQ::AXIS_PHASE);
    m_DSpoil2P.setAxis(SEQ::AXIS_PHASE);
    m_DSpoil1R.setAxis(SEQ::AXIS_READOUT);
    m_DSpoil2R.setAxis(SEQ::AXIS_READOUT);

    // The bipolar SBB uses a Maxwell correction and has to take care of including 
    // a margin for the additional moments into the max amplitude. A margin of 0.1mT/m
    // is enough to cover all possible Maxwell corrections.
    m_dMaxAmpl -= 0.1;

    setRFPulseThicknessFactor(dRFPulseEmpiricalFactor);
}



// ===========================================================================
///   This destructor does nothing.
SBBDiffusion_Bipolar::~SBBDiffusion_Bipolar()
// ===========================================================================
{
}


bool SBBDiffusion_Bipolar::calcSliceAdjSBBRFInfo(
    MrProt                                     &rMrProt,        //< IMP: points to the protocol structure.
    SeqLim                                     &rSeqLim,        //< IMP: points to the sequence limits structure.
    SeqExpo                                    &rSeqExpo,       //< IMP: points to the sequence exports structure
    const SLICEADJ::sCuboidGeometry              &sSliceAdjCuboid,  //< IMP: cuboid geometry (single)
    std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo        //< EXP: RF info
    )
{
   // nothing to do here as all energy is handled by pulse SBBs
    return true;
}


// ===========================================================================
/// Implementation of the pure virtual base class method
// ===========================================================================
bool SBBDiffusion_Bipolar::prepInit(MrProt & /* &rMrProt */, SeqLim & /* &rSeqLim */, SeqExpo & /* &rSeqExpo */)
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
    // from the desired one. A value of 3.0 is the absolute minimum that should
    // be used here.
    double dSpoilerFactor = m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/diffusion_spoil_factor", 3.0);
    m_dRefSpoilMoment = m_dReadoutMoment * dSpoilerFactor / 1000;

    m_DSpoil1P.setStartTime(0);
    if(!prepGradMoment(&m_DSpoil1P, m_dRefSpoilMoment))
    {
        SEQ_TRACE_ERROR.print("ERROR: prepGradMoment() failed.");
        return false;
    }

    m_DSpoil2P = m_DSpoil1P;

    return true;
}

// ===========================================================================
/// Implementation of the pure virtual base class method
// ===========================================================================
bool SBBDiffusion_Bipolar::prepFinal(double dMaxRequestedBValue, bool bIsContextPrepForBinarySearch)
// ===========================================================================
{
    // ---------------------------
    // Prepare spoiler gradients
    // ---------------------------

    // m_DSpoil1P is the reference spoiler which has already been prepared in ::prepInit

    m_DSpoil1R = m_DSpoil1P;
    m_DSpoil1R.setAxis(SEQ::AXIS_READOUT);   // axis must be set again!

    if(! m_DSpoil1R.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: DSpoil1R.prep failed.");
        }
        setNLSStatus(m_DSpoil1R.getNLSStatus());
        return false;
    }

    m_DSpoil1S = m_DSpoil1P;
    m_DSpoil1S.setAxis(SEQ::AXIS_SLICE);

    if(! m_DSpoil1S.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: DSpoil1S.prep failed.");
        }
        setNLSStatus(m_DSpoil1S.getNLSStatus());
        return false;
    }

    m_DSpoil2R = m_DSpoil2P;
    m_DSpoil2R.setAxis(SEQ::AXIS_READOUT);

    if(! m_DSpoil2R.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: DSpoil2R.prep failed.");
        }
        setNLSStatus(m_DSpoil2R.getNLSStatus());
        return false;
    }

    m_DSpoil2S = m_DSpoil2P;
    m_DSpoil2S.setAxis(SEQ::AXIS_SLICE);

    if(! m_DSpoil2S.prep())
    {
        if(!bIsContextPrepForBinarySearch)
        {
            SEQ_TRACE_ERROR.print("ERROR: DSpoil2S.prep failed.");
        }
        setNLSStatus(m_DSpoil2S.getNLSStatus());
        return false;
    }


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


    // prepare and check all diffusion gradients
    for (auto& elem : m_maDG)
    {
        int index = 0;
        for (auto pgrad : elem.second)
        {
            if (!elem.first.compare("p"))
            {
                pgrad->setAmplitude(dAmpl * m_Didi.getX(0));
                pgrad->setAxis(SEQ::AXIS_PHASE);
            }
            if (!elem.first.compare("r"))
            {
                pgrad->setAmplitude(dAmpl * m_Didi.getY(0));
                pgrad->setAxis(SEQ::AXIS_READOUT);
            }
            if (!elem.first.compare("s"))
            {
                pgrad->setAmplitude(dAmpl * m_Didi.getZ(0));
                pgrad->setAxis(SEQ::AXIS_SLICE);
            }

            if ((index % 2) == 1)
            {
                pgrad->setAmplitude(pgrad->getAmplitude() * (-1.0));
            }

            pgrad->setDuration(m_maDG["p"][index]->getDuration());
            pgrad->setRampUpTime(m_maDG["p"][index]->getRampUpTime());
            pgrad->setRampDownTime(m_maDG["p"][index]->getRampDownTime());
            pgrad->setStartTime(m_maDG["p"][index]->getStartTime());
            ++index;

            if (!pgrad->prep())
            {
                if (!bIsContextPrepForBinarySearch)
                {
                    SEQ_TRACE_ERROR.print("ERROR: diffusion gradient preparation failed ");
                }
                setNLSStatus(m_DG1p.getNLSStatus());
                return false;
            }

            if (!pgrad->check())
            {
                if (!bIsContextPrepForBinarySearch)
                {
                    SEQ_TRACE_ERROR.print("ERROR: diffusion gradient check failed ");
                }
                setNLSStatus(m_DG1p.getNLSStatus());
                return false;
            }

            if (pgrad->getTotalTime() < 300)
            {
                if (!bIsContextPrepForBinarySearch)
                {
                    SEQ_TRACE_ERROR.print("ERROR: Rotation matrix too short for diffusion gradient");
                }
                setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                return false;
            }
        }
    }





    //----------------------------
    // Calculate Export Parameters
    //----------------------------
    m_RFInfoPerRequest    = m_SBB_RF_Refoc1.getRFInfoPerRequest() + m_SBB_RF_Refoc2.getRFInfoPerRequest();
    m_RFInfoPerRequestMB  = m_SBB_RF_Refoc1.getRFInfoPerRequestMB()  + m_SBB_RF_Refoc2.getRFInfoPerRequestMB();

    // The times must be calculated from the events itself
    // (MrProt->TE must not be used, as SBBDurationPerRequest contains the time *requested* 
    // to run the SBB for all b values. (Sorry, don't blame me, I didn't design this)
    m_lPreEchoTimeContrib  = m_lSpinPrepTimeEnhancement + m_DG1p.getTotalTime() + m_DSpoil1P.getTotalTime() + m_SBB_RF_Refoc1.getDurationPerRequest()/2;
    m_lPostEchoTimeContrib = m_lStimoDelayus            + m_DG4p.getTotalTime() + m_DSpoil2P.getTotalTime() + m_SBB_RF_Refoc2.getDurationPerRequest()/2;

    setSBBDurationPerRequest(m_lActualTE - m_lADCusTillEcho - m_lSpinPrepTimeus);
    /*
#ifdef DEBUG
    SEQ_TRACE_ALWAYS.print("on exit: Energy SB=%f,   Energy MB=%f,   Duration=%ld,   PreTime=%ld,   PostTime=%ld, MrProt->TE[%d]=%ld",
    m_RFInfoPerRequest.getPulseEnergyWs(), m_RFInfoPerRequestMB.getPulseEnergyWs(),
    m_lPreEchoTimeContrib, m_lPostEchoTimeContrib, m_iTEArrayIndex, m_vlTE[m_iTEArrayIndex] ) ;
    #endif
    */


    m_bIsMagnetizationInverted = false;

    return true;
}



// ===========================================================================
///	Play out the run-time events for the scan actually selected.
/**   The diffusion-specific MDH flags (REP counter, ICEProgramParams
for diffusion direction vector) will be filled out.

In MDDW and 3-scan trace mode, the diffusion gradients are in general
specified in the XYZ magnet coordinate system. In former versions of the
SBBDiffusion, the were played out using the Unity rotation matrix.
This behaviour has been changed because adaptive spoiling requires
to now the sign of the diffusion gradient with respect to the imaging
gradients. In order to optimize the gradient amplitudes, we will not
just use the slice gradient matrix from now on, but we will transform
the XYZ gradient into the PRS coordinate system.

*/
bool SBBDiffusion_Bipolar::runSBB(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC)
// ===========================================================================
{
    MrProtFacade protFacade(rMrProt);

    MaxwellCorrection sMaxwellCorrection{};

    long   lPolarity     = 1;         // sign of diffusion direction
    double dBValue       = 0.;        // actual b-value

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

        double da = m_DSpoil1S.getAmplitude();
        m_dSpoilerAmplitude = (da < 0 ? -da : da);

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
    const double dAmpl = m_dAmpl * sqrt(dBValue/m_dMaxPossibleBValue);

    double dAmpl_p = 0.0;
    double dAmpl_r = 0.0;
    double dAmpl_s = 0.0;

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
        pSLC,            // Input: getSliceIndex, getROT_MATRIX()
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

    m_DG1p.prepAmplitude(dAmpl_p);
    m_DG1r.prepAmplitude(dAmpl_r);
    m_DG1s.prepAmplitude(dAmpl_s);
    m_DG2p.prepAmplitude(-dAmpl_p);
    m_DG2r.prepAmplitude(-dAmpl_r);
    m_DG2s.prepAmplitude(-dAmpl_s);
    m_DG3p.prepAmplitude(dAmpl_p);
    m_DG3r.prepAmplitude(dAmpl_r);
    m_DG3s.prepAmplitude(dAmpl_s);
    m_DG4p.prepAmplitude(-dAmpl_p);
    m_DG4r.prepAmplitude(-dAmpl_r);
    m_DG4s.prepAmplitude(-dAmpl_s);




    // -----------------------------------
    // Prepare spoiler gradient amplitudes 
    // -----------------------------------

    if(m_DG1p.getDuration() - m_DG4p.getDuration() < 0)
    {
        m_dSpoilerFactor1 = -0.42;
        m_dSpoilerFactor2 = -1.0;
    }
    else
    {
        m_dSpoilerFactor1 = +1.0;
        m_dSpoilerFactor2 = +0.42;
    }

    m_DSpoil1P.prepAmplitude(sign(dAmpl_p) * m_dSpoilerFactor1 * m_dSpoilerAmplitude * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/spoiler1_ampl_X", 1.0));
    m_DSpoil1R.prepAmplitude(sign(dAmpl_r) * m_dSpoilerFactor2 * m_dSpoilerAmplitude * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/spoiler1_ampl_Y", 1.0));
    m_DSpoil1S.prepAmplitude(-sign(dAmpl_s) * m_dSpoilerFactor1 * m_dSpoilerAmplitude * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/spoiler1_ampl_Z", 1.0));

    m_DSpoil2P.prepAmplitude(sign(dAmpl_p) * m_dSpoilerFactor2 * m_dSpoilerAmplitude * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/spoiler2_ampl_X", 1.0));
    m_DSpoil2R.prepAmplitude(sign(dAmpl_r) * m_dSpoilerFactor1 * m_dSpoilerAmplitude * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/spoiler2_ampl_Y", 1.0));
    m_DSpoil2S.prepAmplitude(-sign(dAmpl_s) * m_dSpoilerFactor2 * m_dSpoilerAmplitude * m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/spoiler2_ampl_Z", 1.0));


    if(m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_spoiler_info", false))
    {
        SEQ_TRACE_ALWAYS.print("AmplDG1=%5.1f, DeltaMom1/4: %ld, MomSp1P=%7.0f, MomSp2P=%7.0f",
                   m_DG1p.getAmplitude(),
                   m_DG1p.getDuration() - m_DG4p.getDuration(),
                   m_DSpoil1P.getAmplitude() * static_cast<double> (m_DSpoil1P.getDuration()),
                   m_DSpoil2P.getAmplitude() * static_cast<double> (m_DSpoil2P.getDuration()));
    }


    // ---------------------------------------------------------------------------
    // Calculate b matrix and prepare header information
    // ---------------------------------------------------------------------------

    calcBMatrix();

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


    if(m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_diffusion_gradient_value", false))
    {
        SEQ_TRACE_ALWAYS.print("Moments X:  %8.0f  %8.0f  %8.0f   %8.0f  | Delta=%8.0f",
                   dAmpl_p * static_cast<double>(m_DG1p.getDuration()), -dAmpl_p * static_cast<double>(m_DG2p.getDuration()),
                   dAmpl_p * static_cast<double>(m_DG3p.getDuration()), -dAmpl_p * static_cast<double>(m_DG4p.getDuration()),
                   dAmpl_p * static_cast<double>((labs(m_DG1p.getDuration())-labs(m_DG4p.getDuration()))));
        SEQ_TRACE_ALWAYS.print("Moments Y:  %8.0f  %8.0f  %8.0f   %8.0f  | Delta=%8.0f",
                   dAmpl_r * static_cast<double>(m_DG1r.getDuration()), -dAmpl_r * static_cast<double>(m_DG2r.getDuration()),
                   dAmpl_r * static_cast<double>(m_DG3r.getDuration()), -dAmpl_r * static_cast<double>(m_DG4r.getDuration()),
                   dAmpl_r * static_cast<double>((labs(m_DG1r.getDuration())-labs(m_DG4r.getDuration()))));
        SEQ_TRACE_ALWAYS.print("Moments Z:  %8.0f  %8.0f  %8.0f   %8.0f  | Delta=%8.0f",
                   dAmpl_s * static_cast<double>(m_DG1s.getDuration()), -dAmpl_s * static_cast<double>(m_DG2s.getDuration()),
                   dAmpl_s * static_cast<double>(m_DG3s.getDuration()), -dAmpl_s * static_cast<double>(m_DG4s.getDuration()),
                   dAmpl_s * static_cast<double>((labs(m_DG1s.getDuration())-labs(m_DG4s.getDuration()))));
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
        m_pADC           // Output: Mdh.setIceProgramPara 
        );

    if(NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        return setNLSStatus(lStatus);
    }


    // --------------------------
    // Maxwell Corrections
    // --------------------------

    // Note: The correction gets limited to a) the concomittant fields of the 
    //       (dominating) readout and b) to the induced local slice-gradient.
    //       In principle, one could also consider the induced gradients along
    //       the read- and phase-directions. However, this would imply to adapt
    //       the geometrical corrections applied in the Maxwell-correction
    //       Ice-functor correspondingly.
    // Note: Once the 'dynamic-adjustment'-framework is available, one could
    //       integrate all 0th- and 1st-order Maxwell corrections (which are
    //       currently distributed in various sequence modules).
    //---------------------------------------------------------------------------

    double dMaxwellCorrectionTermP = 0.0;
    double dMaxwellCorrectionTermR = 0.0;
    double dMaxwellCorrectionTermS = 0.0;

    // Unit test cannot handle Maxwell correction only apply correction if UT is not active
    // Also, slice position is not unique in SMS, also disable the correction in this case
    if (!(SeqUT.isUnitTestActive() || protFacade.isSliceAcceleration()))
    {
        try
        {
            sMaxwellCorrection.setParameters(pSLC, &m_DG1p, &m_DG1r, &m_DG1s);

            // the correction is calculated once per EPI readout and used for all diffusion gradients
            // with different timing in the bipolar scheme, so ramp correction is not taken into account
            sMaxwellCorrection.calcMaxwellCorrection(false);
        }
        catch (const std::exception& e)
        {
            setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
            return false;
        }

        dMaxwellCorrectionTermP = sMaxwellCorrection.getMaxwellGradAmplP();
        dMaxwellCorrectionTermR = sMaxwellCorrection.getMaxwellGradAmplR();
        dMaxwellCorrectionTermS = sMaxwellCorrection.getMaxwellGradAmplS();

        if (m_debugSettings.getDefaultSetting<bool>("EP2D_DIFF/dump_diffusion_maxwell_correction", false))
        {
            SEQ_TRACE_ALWAYS.print("Maxwell compensation: SlicePos P=%fmm, R=%fmm, S=%fmm", pSLC->getSliceOffCenterPE(), pSLC->getSliceOffCenterRO(), pSLC->getSliceShift());
            SEQ_TRACE_ALWAYS.print("Maxwell compensation: DiffGrad P=%fmT/m, R=%fmT/m, S=%fmT/m", dAmpl_p, dAmpl_r, dAmpl_s);
            SEQ_TRACE_ALWAYS.print("Maxwell compensation: P=%fmT/m, R=%fmT/m, S=%fmT/m", dMaxwellCorrectionTermP, dMaxwellCorrectionTermR, dMaxwellCorrectionTermS);
        }
    }

#ifdef QUIETDWI
        if (!prepMomentOffset(5000)) {
            return false;
        }
#endif // QUIETDWI


        // Inverting the sign of all diffusion gradients does not change the sign of
        // the corresponding Maxwell correction terms.

        m_DG1p.prepAmplitude(dAmpl_p - dMaxwellCorrectionTermP);
        m_DG1r.prepAmplitude(dAmpl_r - dMaxwellCorrectionTermR);
        m_DG1s.prepAmplitude(dAmpl_s - dMaxwellCorrectionTermS);
        m_DG2p.prepAmplitude(-dAmpl_p - dMaxwellCorrectionTermP);
        m_DG2r.prepAmplitude(-dAmpl_r - dMaxwellCorrectionTermR);
        m_DG2s.prepAmplitude(-dAmpl_s - dMaxwellCorrectionTermS);
        m_DG3p.prepAmplitude(dAmpl_p - dMaxwellCorrectionTermP);
        m_DG3r.prepAmplitude(dAmpl_r - dMaxwellCorrectionTermR);
        m_DG3s.prepAmplitude(dAmpl_s - dMaxwellCorrectionTermS);
        m_DG4p.prepAmplitude(-dAmpl_p - dMaxwellCorrectionTermP);
        m_DG4r.prepAmplitude(-dAmpl_r - dMaxwellCorrectionTermR);
        m_DG4s.prepAmplitude(-dAmpl_s - dMaxwellCorrectionTermS);


    // --------------------------
    // Now play out the events...
    // --------------------------

    if(m_lSpinPrepTimeEnhancement)
    {
        // Event table #0: Insert fill time before SBB
        // -------------------------------------------
        fRTEBInit(pSLC->getROT_MATRIX());
        fRTEI(m_lSpinPrepTimeEnhancement, 0, 0, 0, 0, 0, 0, 0);

        if(SeqUT.isUnitTestActive())
        {
            mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'C', 0, pSLC->getSliceIndex(), 0, 0);
        }

        if(setNLSStatus(fRTEBFinish()))
        {
            SEQ_TRACE_ERROR.print("ERROR: Execution of event block 0 failed.");
            return false;
        }
    }


    // Event table #1: diffusion gradient 1 + spoiler
    // ----------------------------------------------
    fRTEBInit(pSLC->getROT_MATRIX());

    if(m_bDiffusionGradientsEnabled)
    {
        runGradient(&m_DG1p);
        runGradient(&m_DG1r);
        runGradient(&m_DG1s);
    }

    fRTEI(m_DG1p.getTotalTime(), 0, 0, 0, &m_DSpoil1P, &m_DSpoil1R, &m_DSpoil1S, 0);
    fRTEI(m_DG1p.getTotalTime()+m_DSpoil1P.getTotalTime(), 0, 0, 0, 0, 0, 0, 0);

    if(SeqUT.isUnitTestActive())
    {
        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'D', 0, pSLC->getSliceIndex(), 0, 0);
    }

    if(setNLSStatus(fRTEBFinish()))
    {
        SEQ_TRACE_ERROR.print("ERROR: Execution of event block 1 failed.");
        return false;
    }


    // Event table #2: first inversion pulse
    // -------------------------------------
    m_SBB_RF_Refoc1.setRunMode(m_eSliceAccelRFRunMode);

    // play out SBB

    if(!m_SBB_RF_Refoc1.run(rMrProt, rSeqLim, rSeqExpo, pSLC))
    {
        SEQ_TRACE_ERROR.print("ERROR: m_SBBMultibandRFRefoc1.run failed.");
        return false;
    }

    // Event table #3: spoiler / low b spoiler resp. diffusion gradients / spoiler
    // -----------------------------------

    fRTEBInit(pSLC->getROT_MATRIX());

    fRTEI(0, 0, 0, 0, &m_DSpoil1P, &m_DSpoil1R, &m_DSpoil1S, 0);

#ifdef QUIETDWI
    // Events for additional offset diffusion gradient in qDWI 
    fRTEI(m_DG3p.getStartTime() - m_DGoffp.getTotalTime(), 0, 0, 0, &m_DGoffp, 0, 0, 0);
    fRTEI(m_DG3r.getStartTime() - m_DGoffr.getTotalTime(), 0, 0, 0, 0, &m_DGoffr, 0, 0);
    fRTEI(m_DG3s.getStartTime() - m_DGoffs.getTotalTime(), 0, 0, 0, 0, 0, &m_DGoffs, 0);
#endif // QUIETDWI

    if(m_bDiffusionGradientsEnabled)
    {
        runGradient(&m_DG2p);
        runGradient(&m_DG2r);
        runGradient(&m_DG2s);
        runGradient(&m_DG3p);
        runGradient(&m_DG3r);
        runGradient(&m_DG3s);
    }

    fRTEI(m_DG3p.getStartTime()+m_DG3p.getTotalTime(), 0, 0, 0, &m_DSpoil2P, &m_DSpoil2R, &m_DSpoil2S, 0);
    fRTEI(m_DSpoil1P.getTotalTime() + m_DG2p.getTotalTime() + m_DG3p.getTotalTime()+m_DSpoil2P.getTotalTime(), 0, 0, 0, 0, 0, 0, 0);

    if(SeqUT.isUnitTestActive())
    {
        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'D', 0, pSLC->getSliceIndex(), 0, 0);
    }

    if(setNLSStatus(fRTEBFinish()))
    {
        SEQ_TRACE_ERROR.print("ERROR: Execution of event block 3 failed.");
        return false;
    }


    // Event table #4: second inversion pulse
    // --------------------------------------
    m_SBB_RF_Refoc2.setRunMode(m_eSliceAccelRFRunMode);

    // play out SBB

    if(!m_SBB_RF_Refoc2.run(rMrProt, rSeqLim, rSeqExpo, pSLC))
    {
        SEQ_TRACE_ERROR.print("ERROR: m_SBBMultibandRFRefoc2.run failed.");
        return false;
    }



    // Event table #5: spoiler / low b spoiler resp. diffusion gradient 4 / stimulation delay
    // --------------------------------------------------------------------------------------
    fRTEBInit(pSLC->getROT_MATRIX());

    fRTEI(0, 0, 0, 0, &m_DSpoil2P, &m_DSpoil2R, &m_DSpoil2S, 0);

    if(m_bDiffusionGradientsEnabled)
    {
        runGradient(&m_DG4p);
        runGradient(&m_DG4r);
        runGradient(&m_DG4s);
    }

    fRTEI(m_DSpoil2P.getTotalTime() + m_DG4p.getTotalTime() + m_lStimoDelayus, 0, 0, 0, 0, 0, 0, 0);

    if(SeqUT.isUnitTestActive())
    {
        mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunKernel, 'D', 0, pSLC->getSliceIndex(), 0, 0);
    }

    if(setNLSStatus(fRTEBFinish()))
    {
        SEQ_TRACE_ERROR.print("ERROR: Execution of event block 5 failed.");
        return false;
    }
    // SEQ_TRACE_ALWAYS.print("EventBlock DG4 duration=%d", m_DG4p.getTotalTime() );


    setNLSStatus (MRI_SBB_SBB_NORMAL);
    return true;
}



// ===========================================================================
/// Calculate the timing of the diffusion gradients in oder achieve the highest b value.
/**
\pre The spoilers m_DSpoil1P and m_DSpoil2P must already have been prepared.
\pre The member variables used for configuration of the SBB must have been set.
*/
bool SBBDiffusion_Bipolar::prepTiming(MrProt  &rMrProt, SeqLim & rSeqLim, SeqExpo &rSeqExpo, long lActualTE)
// ===========================================================================
{
    MrProtFacade protFacade(rMrProt);

    // Flag: Print a message before every return statement
    bool bDebugReturn = ((!(rSeqLim.isContextPrepForBinarySearch() || rSeqLim.isContextPrepForScanTimeCalculation())) && !m_bResolve);

    // Duration of gradient no 1 / 2 / 3 / 4 :
    std::array<double, 4> adt;


    // -------------------------------------------------
    // first we have to prepare the RF SBBs as their timing  
    // does not depend on the diffusion gradients
    // -------------------------------------------------
    if(!prepRF(rMrProt, rSeqLim, rSeqExpo))
        return false;


    // -------------------------
    /** Calculate gradient ratios */
    // -------------------------
    {
        // Algorithm from Oliver Heid: Technical Note XIV
        // or O. Heid: Eddy-Current-Nulled Diffusion Weighting, Proc. ISMRM 2000, p. 799

        // Physics plus hardware raster constrictions require that ts has a 
        // 2*GRAD_RASTER_TIME raster. This may be achieved by an enhanced fill time:
        m_lSpinPrepTimeEnhancement = (m_lADCusTillEcho + m_lSpinPrepTimeus + m_lStimoDelayus) % (2*GRAD_RASTER_TIME);

        // dts := time between end of SBB and echo center
        // Due to the upper tests, rounding is no longer necessary here: ts is on double gradient raster
        double dts     = static_cast<double>(m_lADCusTillEcho + m_lSpinPrepTimeus + m_lSpinPrepTimeEnhancement + m_lStimoDelayus);
        // dt := time available for the diffusion encoding SBB
        // dt is on double gradient raster:
        double dt      = static_cast<double>(lActualTE)- dts;

        double dTimeConstant = m_debugSettings.getDefaultSetting<double>("EP2D_DIFF/eddy_current_time_constant", 30000.0);
        double dlambda = 0.0;

        if(dTimeConstant > 0)
            dlambda = 1.0 / dTimeConstant;
        else
            dlambda = 1.0 / 30000.0;

        adt[0] = (1. / dlambda) * log((2. + exp(dlambda * dt / 2.)  + exp(-dlambda * dt / 2.)) /
                                   (2. * exp(dlambda * dts / 2.) + 2. * exp(-dlambda * dt / 2.)));

        adt[0] = fSDSDoubleRoundUp(0.0, m_MaxValueForRounding, adt[0], (double)GRAD_RASTER_TIME);
        adt[1] = (dt        / 2. - adt[0]);
        adt[2] = (dts       / 2. + adt[0]);
        adt[3] = ((dt - dts) / 2. - adt[0]);


        // Do some cross checks
        // --------------------

        // adt[0] holds at least: Spoiler, slicegrad up, half rf pulse, trapezoid diff grad
        if(adt[0] < (m_DSpoil1P.getTotalTime() + m_SBB_RF_Refoc1.getGradientPointer()->getFlatTopTime()/2 + 2 * m_lRampTime + m_SBB_RF_Refoc1.getGradientPointer()->getRampUpTime()))
        {
            if(bDebugReturn)
            {
                SEQ_TRACE_ERROR.print("ERROR: adt[0]=%f too short", adt[0]);
            }
            setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);
            return false;
        }

        // adt[3] holds at least: grad up, grad down, slicegrad up, half rf pulse, trapezoid diff grad
        if(adt[3] < (m_DSpoil2P.getTotalTime() + m_SBB_RF_Refoc2.getGradientPointer()->getFlatTopTime()/2 + 2 * m_lRampTime + m_SBB_RF_Refoc2.getGradientPointer()->getRampUpTime()))
        {
            if(bDebugReturn)
            {
                SEQ_TRACE_ERROR.print("ERROR: adt[3]=%f too short", adt[3]);
            }
            setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);
            return false;
        }

        // The sum of the four blocks must be the total SBB time
        if(fGSLAlmEqual(std::accumulate(std::begin(adt),std::end(adt),0.0), dt) == false)
        {
            if(bDebugReturn)
            {
                SEQ_TRACE_ERROR.print("ERROR: Sum of ratios %f not total time t %f", std::accumulate(std::begin(adt), std::end(adt), 0.0), dt);
            }
            setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);
            return false;
        }

#ifdef DEBUG
        bool bWorkAsDebug = true;
#else
        bool bWorkAsDebug = SeqUT.isUnitTestActive();
#endif
        if(bWorkAsDebug)
        {
            // some lines above we have assured that adt[0] is on double gradient raster.
            // This should enforce that the other values are  also on gradient raster.
            // Therefore the following checks are not really necessary.
            for (auto dti : adt)
            {
                if (static_cast<long>(dti) % GRAD_RASTER_TIME)
                {
                    if (bDebugReturn)
                    {
                        SEQ_TRACE_ERROR.print("ERROR: dt[]=%f not on gradient raster", dti);
                    }
                    setNLSStatus(MRI_SBB_SBB_DIFFUSION_ERROR);
                    return false;
                }
            }
        } // of DEBUG      


    } // END of calculation gradient ratios

    //  ----------------------------------------------------
    /** Set up timing and calculate maximum possible b value */
    //  ----------------------------------------------------
    // Gradient amplitude for dummy preparation - real preparations will take place below
    double dAmpl = 5.;

    int index = 0;
    for (auto pgrad : m_maDG["p"])
    {
        pgrad->setRampTimes(m_lRampTime);
        pgrad->setAxis(SEQ::AXIS_PHASE);
        
        if ((index % 2) == 1)
        {
            pgrad->setAmplitude(-dAmpl);
        }
        else
        {
            pgrad->setAmplitude(dAmpl);
        }

        // duration and start time is different for all of them
        long lDuration = 0;
        long lStartTime = 0;

        switch (index)
        {
        case 0:
            lDuration = static_cast<long>(adt[index]) - m_lRampTime - m_SBB_RF_Refoc1.getDurationPerRequest() / 2 - m_DSpoil1P.getTotalTime();
            lStartTime = 0;
            break;
        case 1:
            lDuration = static_cast<long>(adt[index]) - m_SBB_RF_Refoc1.getDurationPerRequest() / 2 - m_DSpoil1P.getTotalTime() - m_DG2p.getRampDownTime();
            lStartTime = m_DSpoil1P.getTotalTime();
            break;
        case 2:
            lDuration = static_cast<long>(adt[index]) - m_SBB_RF_Refoc2.getDurationPerRequest() / 2 - m_DSpoil2P.getTotalTime() - m_DG3p.getRampDownTime();
            lStartTime = m_DG2p.getStartTime() + m_DG2p.getTotalTime();
            break;
        case 3:
            lDuration = static_cast<long>(adt[index]) - m_SBB_RF_Refoc2.getDurationPerRequest() / 2 - m_DSpoil2P.getTotalTime() - m_lRampTime;
            lStartTime = m_DSpoil2P.getTotalTime();
            break;
        default:
            return false;
        }

        if (lDuration < 0)
            return false;

        pgrad->setDuration(lDuration);
        pgrad->setStartTime(lStartTime);

        if (pgrad->getDuration() < pgrad->getRampUpTime())
            return false;

        ++index;
    }


    // ------------
    // Final checks and prep
    // ------------

    /** Now everything has been calculated.
    There remain some calculations which will be done in run()
    to calulate event start times.
    We now check all these calculations to be sure they are ok.
    */

    if(m_DG1p.getTotalTime() - m_DSp1.getTotalTime() < 0)
    {
        if(bDebugReturn)
        {
            SEQ_TRACE_ERROR.print("ERROR: negative event start time ");
        }
        setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
    }

    if(m_DG2p.getTotalTime() + m_DG3p.getTotalTime() - m_DSp1.getTotalTime() < 0)
    {
        if(bDebugReturn)
        {
            SEQ_TRACE_ERROR.print("ERROR: negative event start time: %ld + %ld - %ld < 0 ",
                        m_DG2p.getTotalTime(), m_DG3p.getTotalTime(), m_DSp1.getTotalTime());
        }
        setNLSStatus (MRI_SBB_SBB_DIFFUSION_ERROR);
        return false;
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


    // final prepare of phase axis diffusion gradients
    index = 0;
    for (auto pgrad : m_maDG["p"])
    {
        if ((index % 2) == 1)
        {
            pgrad->setAmplitude(-dAmpl);
        }
        else
        {
            pgrad->setAmplitude(dAmpl);
        }

        if (!pgrad->prep())
        {
            if (bDebugReturn)
            {
                SEQ_TRACE_ERROR.print("ERROR: diffusion gradient preparation failed.");
            }
            setNLSStatus(pgrad->getNLSStatus());
            return false;
        }
        if (!pgrad->check())
        {
            if (bDebugReturn)
            {
                SEQ_TRACE_ERROR.print("ERROR: diffusion gradient check failed.");
            }
            setNLSStatus(pgrad->getNLSStatus());
            return false;
        }

        ++index;
    }


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
double SBBDiffusion_Bipolar::calcBValue(void)
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

    // Prepare for b-value calcuation (set up diffusion encoding events):
    //  b-value only
    PrepBMatrix(true);

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
void SBBDiffusion_Bipolar::calcBMatrix(void)
// ===========================================================================
{
    // Set gamma of actual nucleus
    m_theBMatrix.setGamma(m_dGamma);

    // Prepare for b-matrix calcuation (set up diffusion encoding events):
    //  complete b-matrix
    PrepBMatrix();

    // Calculate b-matrix
    if(!m_theBMatrix.bCalcBMatrix(m_dBxx, m_dBxy, m_dBxz, m_dByy, m_dByz, m_dBzz))
    {
        SEQ_TRACE_WARN.print("WARNING: Non-null moment");
    }

    // Calculate b-value
    m_theBMatrix.bCalcBValue(m_dBValue);
}

bool SBBDiffusion_Bipolar::PrepBMatrix(bool bBValueOnly)
// ===========================================================================
{
    // Initialize b-matrix calculation
    m_theBMatrix.Reset();

    // Fill in basic diffusion encoding events
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DG1p);           // PE-axis first diffusion encoding gradient
    m_theBMatrix.bAddEvent(getEventTime(SpoilGrad1_Start), m_DSpoil1P);       // PE-axis first spoiling gradient

    m_theBMatrix.bAddEvent(getEventTime(RefocRF1_Center), m_SBB_RF_Refoc1.getRFPulsePointer());           // First refocussing RF

    m_theBMatrix.bAddEvent(getEventTime(SpoilGrad2_Start), m_DSpoil1P);       // PE-axis second spoiling gradient
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DG2p);           // PE-axis second diffusion encoding gradient
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad3_Start), m_DG3p);           // PE-axis third diffusion encoding gradient
    m_theBMatrix.bAddEvent(getEventTime(SpoilGrad3_Start), m_DSpoil2P);       // PE-axis third spoiling gradient

    m_theBMatrix.bAddEvent(getEventTime(RefocRF2_Center), m_SBB_RF_Refoc2.getRFPulsePointer());           // Second refocussing RF

    m_theBMatrix.bAddEvent(getEventTime(SpoilGrad4_Start), m_DSpoil2P);       // PE-axis fourth spoiling gradient
    m_theBMatrix.bAddEvent(getEventTime(DiffGrad4_Start), m_DG4p);           // PE-axis fourth diffusion encoding gradient

    if(!bBValueOnly)
    {
        // Full b-matrix calculation

        m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DG1r);           // RO-axis first diffusion encoding gradient
        m_theBMatrix.bAddEvent(getEventTime(SpoilGrad1_Start), m_DSpoil1R);       // RO-axis first spoiling gradient
        m_theBMatrix.bAddEvent(getEventTime(SpoilGrad2_Start), m_DSpoil1R);       // RO-axis second spoiling gradient
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DG2r);           // RO-axis second diffusion encoding gradient
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad3_Start), m_DG3r);           // RO-axis third diffusion encoding gradient
        m_theBMatrix.bAddEvent(getEventTime(SpoilGrad3_Start), m_DSpoil2R);       // RO-axis third spoiling gradient
        m_theBMatrix.bAddEvent(getEventTime(SpoilGrad4_Start), m_DSpoil2R);       // RO-axis fourth spoiling gradient
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad4_Start), m_DG4r);           // RO-axis fourth diffusion encoding gradient

        m_theBMatrix.bAddEvent(getEventTime(DiffGrad1_Start), m_DG1s);           // SL-axis first diffusion encoding gradient
        m_theBMatrix.bAddEvent(getEventTime(SpoilGrad1_Start), m_DSpoil1S);       // SL-axis first spoiling gradient

        //add versed slice gradient for accurate b-matrix calculation
        sGRAD_PULSE_TRAP* SliceGrad1_1 = m_SBB_RF_Refoc1.getGradientPointer(0);
        long StartTimeSliceGrad1_1 = getEventTime(SliceGrad1_Start);

        m_theBMatrix.bAddEvent(StartTimeSliceGrad1_1, *SliceGrad1_1);           // SL-axis first  verse slice selection gradient
        if(m_SBB_RF_Refoc1.isVerse()) // is Verse
        {
            sGRAD_PULSE_TRAP* SliceGrad1_2 = m_SBB_RF_Refoc1.getGradientPointer(1);
            sGRAD_PULSE_TRAP* SliceGrad1_3 = m_SBB_RF_Refoc1.getGradientPointer(2);
            long StartTimeSliceGrad1_2 = StartTimeSliceGrad1_1 + SliceGrad1_1->getDuration();
            long StartTimeSliceGrad1_3 = StartTimeSliceGrad1_2 + SliceGrad1_2->getDuration();

            m_theBMatrix.bAddEvent(StartTimeSliceGrad1_2, *SliceGrad1_2);           // SL-axis second verse slice selection gradient
            m_theBMatrix.bAddEvent(StartTimeSliceGrad1_3, *SliceGrad1_3);           // SL-axis third  verse slice selection gradient
        }

        m_theBMatrix.bAddEvent(getEventTime(SpoilGrad2_Start), m_DSpoil1S);       // SL-axis second spoiling gradient
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad2_Start), m_DG2s);           // SL-axis second diffusion encoding gradient
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad3_Start), m_DG3s);           // SL-axis third diffusion encoding gradient
        m_theBMatrix.bAddEvent(getEventTime(SpoilGrad3_Start), m_DSpoil2S);       // SL-axis third spoiling gradient


        //add versed slice gradient for accurate b-matrix calculation
        sGRAD_PULSE_TRAP* SliceGrad2_1 = m_SBB_RF_Refoc2.getGradientPointer(0);
        long StartTimeSliceGrad2_1 = getEventTime(SliceGrad2_Start);

        m_theBMatrix.bAddEvent(StartTimeSliceGrad2_1, *SliceGrad2_1);           // SL-axis first  verse slice selection gradient
        if(m_SBB_RF_Refoc2.isVerse()) // is Verse
        {
            sGRAD_PULSE_TRAP* SliceGrad2_2 = m_SBB_RF_Refoc2.getGradientPointer(1);
            sGRAD_PULSE_TRAP* SliceGrad2_3 = m_SBB_RF_Refoc2.getGradientPointer(2);
            long StartTimeSliceGrad2_2 = StartTimeSliceGrad2_1 + SliceGrad2_1->getDuration();
            long StartTimeSliceGrad2_3 = StartTimeSliceGrad2_2 + SliceGrad2_2->getDuration();

            m_theBMatrix.bAddEvent(StartTimeSliceGrad2_2, *SliceGrad2_2);           // SL-axis second verse slice selection gradient
            m_theBMatrix.bAddEvent(StartTimeSliceGrad2_3, *SliceGrad2_3);           // SL-axis third  verse slice selection gradient
        }


        m_theBMatrix.bAddEvent(getEventTime(SpoilGrad4_Start), m_DSpoil2S);       // SL-axis fourth spoiling gradient
        m_theBMatrix.bAddEvent(getEventTime(DiffGrad4_Start), m_DG4s);           // SL-axis fourth diffusion encoding gradient
    }

    if(m_theBMatrix.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}


bool SBBDiffusion_Bipolar::prepGPALoadDiff(double dAmplitudeX)
{
    double dSign1 = 1.;
    double dSign2 = 1.;

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
    m_lID2X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad2_Start), -dAmplitudeX, m_DG2p.getRampUpTime(), m_DG2p.getRampDownTime(), m_DG2p.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_lID3X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad3_Start), dAmplitudeX, m_DG3p.getRampUpTime(), m_DG3p.getRampDownTime(), m_DG3p.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_lID4X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad4_Start), -dAmplitudeX, m_DG4p.getRampUpTime(), m_DG4p.getRampDownTime(), m_DG4p.getFlatTopTime(), GPABALANCE_X_AXIS);

    // Add spoiler gradients (ensure worst case polarity)
    if(m_DSpoil1P.getAmplitude() < 0.)
    {
        // First and second spoiler get same sign as second (long) diffusion encoding gradient
        dSign1 =  1.;
    }
    else
    {
        dSign1 = -1.;
    }

    if(m_DSpoil2P.getAmplitude() < 0.)
    {
        // Third and fourth spoiler get same sign as third (long) diffusion encoding gradient
        dSign2 = -1.;
    }
    else
    {
        dSign2 =  1.;
    }

    // Spoiler gradient amplitudes are constant => no need to store the ID's
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad1_Start), dSign1 * m_DSpoil1P.getAmplitude(), m_DSpoil1P.getRampUpTime(), m_DSpoil1P.getRampDownTime(), m_DSpoil1P.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad2_Start), dSign1 * m_DSpoil1P.getAmplitude(), m_DSpoil1P.getRampUpTime(), m_DSpoil1P.getRampDownTime(), m_DSpoil1P.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad3_Start), dSign2 * m_DSpoil2P.getAmplitude(), m_DSpoil2P.getRampUpTime(), m_DSpoil2P.getRampDownTime(), m_DSpoil2P.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad4_Start), dSign2 * m_DSpoil2P.getAmplitude(), m_DSpoil2P.getRampUpTime(), m_DSpoil2P.getRampDownTime(), m_DSpoil2P.getFlatTopTime(), GPABALANCE_X_AXIS);

    // Check error status
    if(m_sBalanceDiff.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}

bool SBBDiffusion_Bipolar::prepGPALoadDiff(double dAmplitudeX, double dAmplitudeY, double dAmplitudeZ)
{
    double dSign1 = 1.;
    double dSign2 = 1.;

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
    m_lID2X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad2_Start), -dAmplitudeX, m_DG2p.getRampUpTime(), m_DG2p.getRampDownTime(), m_DG2p.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_lID3X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad3_Start), dAmplitudeX, m_DG3p.getRampUpTime(), m_DG3p.getRampDownTime(), m_DG3p.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_lID4X = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad4_Start), -dAmplitudeX, m_DG4p.getRampUpTime(), m_DG4p.getRampDownTime(), m_DG4p.getFlatTopTime(), GPABALANCE_X_AXIS);

    m_lID1Y = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad1_Start), dAmplitudeY, m_DG1p.getRampUpTime(), m_DG1p.getRampDownTime(), m_DG1p.getFlatTopTime(), GPABALANCE_Y_AXIS);
    m_lID2Y = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad2_Start), -dAmplitudeY, m_DG2p.getRampUpTime(), m_DG2p.getRampDownTime(), m_DG2p.getFlatTopTime(), GPABALANCE_Y_AXIS);
    m_lID3Y = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad3_Start), dAmplitudeY, m_DG3p.getRampUpTime(), m_DG3p.getRampDownTime(), m_DG3p.getFlatTopTime(), GPABALANCE_Y_AXIS);
    m_lID4Y = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad4_Start), -dAmplitudeY, m_DG4p.getRampUpTime(), m_DG4p.getRampDownTime(), m_DG4p.getFlatTopTime(), GPABALANCE_Y_AXIS);

    m_lID1Z = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad1_Start), dAmplitudeZ, m_DG1p.getRampUpTime(), m_DG1p.getRampDownTime(), m_DG1p.getFlatTopTime(), GPABALANCE_Z_AXIS);
    m_lID2Z = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad2_Start), -dAmplitudeZ, m_DG2p.getRampUpTime(), m_DG2p.getRampDownTime(), m_DG2p.getFlatTopTime(), GPABALANCE_Z_AXIS);
    m_lID3Z = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad3_Start), dAmplitudeZ, m_DG3p.getRampUpTime(), m_DG3p.getRampDownTime(), m_DG3p.getFlatTopTime(), GPABALANCE_Z_AXIS);
    m_lID4Z = m_sBalanceDiff.lAddGradient(getEventTime(DiffGrad4_Start), -dAmplitudeZ, m_DG4p.getRampUpTime(), m_DG4p.getRampDownTime(), m_DG4p.getFlatTopTime(), GPABALANCE_Z_AXIS);

    // Add spoiler gradients (ensure worst case polarity)
    if(m_DSpoil1P.getAmplitude() < 0.)
    {
        // First and second spoiler get same sign as second (long) diffusion encoding gradient
        dSign1 =  1.;
    }
    else
    {
        dSign1 = -1.;
    }

    if(m_DSpoil2P.getAmplitude() < 0.)
    {
        // Third and fourth spoiler get same sign as third (long) diffusion encoding gradient
        dSign2 = -1.;
    }
    else
    {
        dSign2 =  1.;
    }

    // Spoiler gradient amplitudes are constant => no need to store the ID's
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad1_Start), dSign1 * m_DSpoil1P.getAmplitude(), m_DSpoil1P.getRampUpTime(), m_DSpoil1P.getRampDownTime(), m_DSpoil1P.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad2_Start), dSign1 * m_DSpoil1P.getAmplitude(), m_DSpoil1P.getRampUpTime(), m_DSpoil1P.getRampDownTime(), m_DSpoil1P.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad3_Start), dSign2 * m_DSpoil2P.getAmplitude(), m_DSpoil2P.getRampUpTime(), m_DSpoil2P.getRampDownTime(), m_DSpoil2P.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad4_Start), dSign2 * m_DSpoil2P.getAmplitude(), m_DSpoil2P.getRampUpTime(), m_DSpoil2P.getRampDownTime(), m_DSpoil2P.getFlatTopTime(), GPABALANCE_X_AXIS);

    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad1_Start), dSign1 * m_DSpoil1P.getAmplitude(), m_DSpoil1P.getRampUpTime(), m_DSpoil1P.getRampDownTime(), m_DSpoil1P.getFlatTopTime(), GPABALANCE_X_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad2_Start), dSign1 * m_DSpoil1P.getAmplitude(), m_DSpoil1P.getRampUpTime(), m_DSpoil1P.getRampDownTime(), m_DSpoil1P.getFlatTopTime(), GPABALANCE_Y_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad3_Start), dSign2 * m_DSpoil2P.getAmplitude(), m_DSpoil2P.getRampUpTime(), m_DSpoil2P.getRampDownTime(), m_DSpoil2P.getFlatTopTime(), GPABALANCE_Y_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad4_Start), dSign2 * m_DSpoil2P.getAmplitude(), m_DSpoil2P.getRampUpTime(), m_DSpoil2P.getRampDownTime(), m_DSpoil2P.getFlatTopTime(), GPABALANCE_Y_AXIS);

    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad1_Start), dSign1 * m_DSpoil1P.getAmplitude(), m_DSpoil1P.getRampUpTime(), m_DSpoil1P.getRampDownTime(), m_DSpoil1P.getFlatTopTime(), GPABALANCE_Z_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad2_Start), dSign1 * m_DSpoil1P.getAmplitude(), m_DSpoil1P.getRampUpTime(), m_DSpoil1P.getRampDownTime(), m_DSpoil1P.getFlatTopTime(), GPABALANCE_Z_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad3_Start), dSign2 * m_DSpoil2P.getAmplitude(), m_DSpoil2P.getRampUpTime(), m_DSpoil2P.getRampDownTime(), m_DSpoil2P.getFlatTopTime(), GPABALANCE_Z_AXIS);
    m_sBalanceDiff.lAddGradient(getEventTime(SpoilGrad4_Start), dSign2 * m_DSpoil2P.getAmplitude(), m_DSpoil2P.getRampUpTime(), m_DSpoil2P.getRampDownTime(), m_DSpoil2P.getFlatTopTime(), GPABALANCE_Z_AXIS);

    // Check error status
    if(m_sBalanceDiff.lGetStatus() != 0)
    {
        SEQ_TRACE_ERROR.print("error adding gradient events.");
        return false;
    }

    return true;
}

bool SBBDiffusion_Bipolar::scaleGPALoadDiff(double dScaleX)
{
    const int iNumberOfEvents       = 4;

    const std::array<long, iNumberOfEvents>   alID    = {m_lID1X, m_lID2X, m_lID3X, m_lID4X};
    const std::array<double, iNumberOfEvents> adScale = {dScaleX, dScaleX, dScaleX, dScaleX};

    // ------------------------------------------------------------------------
    // Scale diffusion encoding events on x-axis
    if (!m_sBalanceDiff.bScaleEvent(alID, adScale))
    {
        SEQ_TRACE_ERROR <<"error scaling gradient events.";
        return false;
    }

    return true;
}

bool SBBDiffusion_Bipolar::scaleGPALoadDiff(double dScaleX, double dScaleY, double dScaleZ)
{
    const int iNumberOfEvents       = 12;

    const std::array<long, iNumberOfEvents>   alID    = {m_lID1X, m_lID2X, m_lID3X, m_lID4X, m_lID1Y, m_lID2Y, m_lID3Y, m_lID4Y, m_lID1Z, m_lID2Z, m_lID3Z, m_lID4Z};
    const std::array<double, iNumberOfEvents> adScale = {dScaleX, dScaleX, dScaleX, dScaleX, dScaleY, dScaleY, dScaleY, dScaleY, dScaleZ, dScaleZ, dScaleZ, dScaleZ};

    // ------------------------------------------------------------------------
    // Scale diffusion encoding events
    if (!m_sBalanceDiff.bScaleEvent(alID, adScale))
    {
        SEQ_TRACE_ERROR << "error scaling gradient events.";
        return false;
    }

    return true;
}

long SBBDiffusion_Bipolar::getEventTime(EnumEventTime eEvent)
{
    // Start of first diffusion encoding gradient is zero
    long lTime0 = 0;
    // Start of first spoiling gradient
    long lTime1 = m_DG1p.getTotalTime();
    // Start of first slice selection gradient 
    long lTime2 = lTime1 + m_DSpoil1P.getTotalTime();
    // First refocussing time
    long lTime3 = lTime2 + m_SBB_RF_Refoc1.getDurationPerRequest()/2;
    // Start of second spoiling gradient
    long lTime4 = lTime2 + m_SBB_RF_Refoc1.getDurationPerRequest();
    // Start of second diffusion encoding gradient
    long lTime5 = lTime4 + m_DSpoil1P.getTotalTime();
    // Start of third diffusion encoding gradient
    long lTime6 = lTime5 + (m_DG3p.getStartTime() - m_DG2p.getStartTime());
    // Start of third spoiling gradient
    long lTime7 = lTime6 + m_DG3p.getTotalTime();
    // Start of second slice selection gradient
    long lTime8 = lTime7 + m_DSpoil2P.getTotalTime();
    // Second refocussing time
    long lTime9 = lTime8 + m_SBB_RF_Refoc2.getDurationPerRequest()/2;
    // Start of fourth spoiling gradient
    long lTime10 = lTime8 + m_SBB_RF_Refoc2.getDurationPerRequest();
    // Start of fourth diffusion encoding gradient
    long lTime11 = lTime10 + m_DSpoil2P.getTotalTime();

    switch(eEvent)
    {
        case DiffGrad1_Start:
            return lTime0;
        case DiffGrad2_Start:
            return lTime5;
        case DiffGrad3_Start:
            return lTime6;
        case DiffGrad4_Start:
            return lTime11;
        case SpoilGrad1_Start:
            return lTime1;
        case SpoilGrad2_Start:
            return lTime4;
        case SpoilGrad3_Start:
            return lTime7;
        case SpoilGrad4_Start:
            return lTime10;
        case SliceGrad1_Start:
            return lTime2;
        case SliceGrad2_Start:
            return lTime8;
        case RefocRF1_Center:
            return lTime3;
        case RefocRF2_Center:
            return lTime9;
        default:
            SEQ_TRACE_ERROR.print("unknown event enumerator: %i", eEvent);
            return -1;
    }
}

bool SBBDiffusion_Bipolar::prepRF(MrProt & rMrProt, SeqLim &rSeqLim, SeqExpo & rSeqExpo)
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


    // ---------------------------------------------------------------------------
    // 4th event block (slice-selective RF pulse)
    // ---------------------------------------------------------------------------

    m_DRF2.setFamilyName(myRFPulseProperties.sFamilyName.c_str());
    m_DRF2.setTypeRefocussing();
    m_DRF2.setFlipAngle(180.0);
    m_DRF2.setInitialPhase(180.0);
    m_DRF2.setThickness(m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_DRF2.setDuration(myRFPulseProperties.lDuration_us);

    if((m_dPreparedSlcThk != m_dRFPulseThicknessFactor * m_dSliceThickness) || (!rSeqLim.isContextPrepForBinarySearch()))
    {
        if(! m_DRF2.prepExternal(rMrProt, rSeqExpo))
        {
            if(!rSeqLim.isContextPrepForBinarySearch())
            {
                SEQ_TRACE_ERROR.print("ERROR: DRF2.prepExternal failed.");
            }
            setNLSStatus(m_DRF2.getNLSStatus());
            return false;
        }
        m_dPreparedSlcThk = m_dRFPulseThicknessFactor * m_dSliceThickness;
    }

    m_SBB_RF_Refoc2.setIdent("SBB_MB_Refoc2");
    if(!prep_SBBMultibandRF(rMrProt, rSeqLim, rSeqExpo, m_SBB_RF_Refoc2, m_DRF2))
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
    SeqUT.setRFThicknessInfo(&m_DRF2, m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_SBB_RF_Refoc1.setRunMode(SINGLE_BAND);
    SeqUT.setRFThicknessInfo(m_SBB_RF_Refoc1.getRFPulsePointer(), m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_SBB_RF_Refoc1.setRunMode(MULTI_BAND);
    SeqUT.setRFThicknessInfo(m_SBB_RF_Refoc1.getRFPulsePointer(), m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_SBB_RF_Refoc2.setRunMode(SINGLE_BAND);
    SeqUT.setRFThicknessInfo(m_SBB_RF_Refoc2.getRFPulsePointer(), m_dRFPulseThicknessFactor * m_dSliceThickness);
    m_SBB_RF_Refoc2.setRunMode(MULTI_BAND);
    SeqUT.setRFThicknessInfo(m_SBB_RF_Refoc2.getRFPulsePointer(), m_dRFPulseThicknessFactor * m_dSliceThickness);

    m_SBB_RF_Refoc1.setRunMode(m_eSliceAccelRFRunMode);
    m_SBB_RF_Refoc2.setRunMode(m_eSliceAccelRFRunMode);
#endif

    return true;
}

bool SEQ_NAMESPACE::SBBDiffusion_Bipolar::getSliceAdjRFInfo(const SLICEADJ::sCuboidGeometry &sSliceAdjCuboid, std::vector<MrProtocolData::SeqExpoRFInfo> &vsRFInfo)
{
    std::vector<MrProtocolData::SeqExpoRFInfo> helpEnergy;

    if(!m_SBB_RF_Refoc1.getSliceAdjRFInfo(sSliceAdjCuboid, vsRFInfo))
        return false;

    if(!m_SBB_RF_Refoc2.getSliceAdjRFInfo(sSliceAdjCuboid, helpEnergy))
        return false;

    // Add contribution from second RF pulse
    for(size_t iI = 0; iI < helpEnergy.size(); ++iI)
        vsRFInfo[iI] += helpEnergy[iI];

    return true;
}

MrProtocolData::SeqExpoRFInfo SEQ_NAMESPACE::SBBDiffusion_Bipolar::getRFInfoPerRequest()
{
    return m_SBB_RF_Refoc1.getRFInfoPerRequest() + m_SBB_RF_Refoc2.getRFInfoPerRequest();
}

long SEQ_NAMESPACE::SBBDiffusion_Bipolar::getSmallestIVIMbValuePossible(bool bIsContextPrepForBinarySearch)
{
    return m_lBValueInc_Limit;
}



