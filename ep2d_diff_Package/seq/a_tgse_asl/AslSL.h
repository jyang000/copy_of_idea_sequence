//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2013  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/X
//	    File: \src\MrImaging\seq\a_tgse_asl\AslSL.h
//	  Author: pfeujodj
//	    Date: 2017-10-16 09:30:03 +02:00
//
//	    Lang: C++
//
//	 Descrip: 3D ASL sequence with CompositeSeqLoop support
//            derived from MrImaging\libSL\StdSL.h with extended structures: ASLLoop, M0Scans
//
//	-----------------------------------------------------------------------------

#pragma once

#include "MrImagingFW/libCSL/RootNode.h"
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

class MrProt;
class SeqLim;
class ASLLoop;          // SUPPORT_CSL

namespace MrProtocolData { class SeqExpo; }

namespace SL_NAME_SPACE  
{

// Forward declaration  
class M0Scans;          // SUPPORT_CSL
class ASLLoop;          // SUPPORT_CSL

class NoiseMeas;
class MeasurementLoop;
class TokTokTok;
class ConcatLoop;
class BreathHoldSync;
class InnerSliceLoop;
class MeasRepDelay;
class ASLPATRefScan;    // SUPPORT_CSL
class PreparingScans;
class PhaseCorrScans;
class ConstLoop;
class kSpaceLoop;
class StdIRASL;         // ECG_TRIGGERING
class InterleavedIR;
class MSat;
class RFSpoiledLeaf;
class CSat;
class RSat;
class SpoilGrad;
class ConcatFillTime;

// ----------------------------------------------------------------------------
/// \ingroup grp_libSL 
///
/// \brief This class implements the standard asl (SUPPORT_CSL) loop structure used for MR imaging 
/// sequences in N4 .
///
/// The standard Composite SeqLoop supports the following features:
/// - <b>Introduction</b> The execution of the \c TokTokTok can be controlled 
///   using the key \c SLK::e_PERF_TOK_TOK_TOK. Allowed values: \n 
///   - \c SL::ALWAYS   Perform TokTokTok prior to each measurement. \c TokTokTok
///                     will be located within the \c MeasurementLoop.
///   - \c SL::ONLY_FIRST_REPETITION    Perform \c TokTokTok prior to the first 
///                     measurement, only. \c TokTokTok will be located before the 
///                     \c MeasurementLoop.
///   - \c SL::NEVER      Do not perform TokTokTok
/// - <b>Preparing scans</b> The excution of preparing sans can be controlled by the
///   key \c SLK::e_PERF_PREP_SCAN. The duration in us of the preparing scans is 
///   given by the key \c SL::l_PREPARING_TIME. Allowed values are:
///   - \c SL::ALWAYS   Perform preparing scans prior to each measurement
///   - \c SL::ONLY_FIRST_REPETITION  Perform preparing scans prior to the first 
///     measurement, only
///   - \c SL::NEVER    Do not perform preparing scans
/// - <b>External reference scans for iPAT</b> The standard Composite SeqLoop supports 
///   external reference scans for iPAT. These scans are enabled / disabled by the 
///   UI corresponding parameters \c rMrProt.getsPat().getucPATMode() and 
///   \c rMrProt.getsPat().getucRefScanMode(). In addition the execution can be 
///   controlled by the key \c SLK::e_PERF_PAT_REF_SCAN. The following values are 
///   supported:
///   - \c SL::ALWAYS     Perform PAT reference scan prior to each measurement
///   - \c SL::ONLY_FIRST_REPETITION  Perform PAT reference scan prior to the first 
///     measurement, only
///   - \c SL::NEVER      Do not perform PAT reference scan
/// - <b>Long and short term averaging</b> The averaging mode is controlled by the UI
///   parameter \c rMrProt.kSpace().averagingMode(). For long term averaging the k-space
///   sampling loop is nested within the averaging loop. For short term averaging the
///   averaging loop is nested within the k-space sampling loop.
/// - <b>Physiological syncronization</b> Two modes are supported which differ in the 
///   position of the phyiologic trigger. The modes are controlled by the key 
///   \c SLK::b_SINGLE_SHOT_TRIGGER. If set to \c true the trigger is executed once per 
///   slice. This mode is intended for single shot measurements. If set to \c false the 
///   trigger is executed prior to the phases loop. This mode is intended for standard 
///   measurements.
/// - <b>Multi slice mode</b> Multi slice modes interleaved, single shot as well as 
///   sequenctial are supported.
///   <b>Note:</b> In contrast to previous SeqLoop implementations there is no concatenation 
///   loop serving as a slice loop for sequential multi slice.
/// - <b>Standard inversion recovery</b>
/// - <b>Dark blood preparation</b> Dark blood preparation is enabled / disabled by the 
///   UI parameter \c rMrProt.preparationPulses().getucDarkBlood(). \n
///   <b>Note:</b> In case of dark blood preparation the TR fill time will be inserted 
///   between the dark blood preparation pulse and the \c SBB_Pre_MP mount point. 
///
/// The class is marked as root node.
/// 
/// <b>Loop Structure</b> \n
///
/// Interleaved multi slice: \n
/// <pre>
/// |--------------------------------------------------------------------  StdSL   ..........................................
/// |
/// | |--------------------------------------------------------U---------  Begin: MeasurementLoop 0 ... 0   .................
/// | |
/// | | |----------------------------------------|                         ..................................................
/// | | | TokTokTok                              |
/// | | |----------------------------------------|
/// | |
/// | | |------------------------------------------------------U---------  Begin: ConcatLoop 0 ... 0   ......................
/// | | |
/// | | | |--------------------------------------------------------------  Begin: PreparingScans 0 ... 0   ..................
/// | | | |
/// | | | | |------------------------------------------------------------  StdIR   ..........................................
/// | | | | |
/// | | | | | |----------------------------------------------------------  Begin: PhaseLoop 0 ... 0   .......................
/// | | | | | |
/// | | | | | | |--------------------------------------------------------  Begin: InnerSliceLoop 0 ... 0   ..................
/// | | | | | | |
/// | | | | | | | |------------------------------------------------------  SBB_Pre_MP   .....................................
/// | | | | | | | |
/// | | | | | | | |------------------------------------------------------
/// | | | | | | |
/// | | | | | | | |------------------------------------------------------  Kernel_MP   ......................................
/// | | | | | | | |
/// | | | | | | | | |----------------------------------------|             ..................................................
/// | | | | | | | | | MedicKernel                            |
/// | | | | | | | | |----------------------------------------|
/// | | | | | | | |
/// | | | | | | | |------------------------------------------------------
/// | | | | | | |
/// | | | | | | | |----------------------------------------|               ..................................................
/// | | | | | | | | TRFill                                 |
/// | | | | | | | |----------------------------------------|
/// | | | | | | |
/// | | | | | | | |------------------------------------------------------  SBB_Post_MP   ....................................
/// | | | | | | | |
/// | | | | | | | |------------------------------------------------------
/// | | | | | | |
/// | | | | | | |--------------------------------------------------------
/// | | | | | |
/// | | | | | | |----------------------------------------|                 ..................................................
/// | | | | | | | TRFillEnd                              |
/// | | | | | | |----------------------------------------|
/// | | | | | |
/// | | | | | |----------------------------------------------------------
/// | | | | |
/// | | | | |------------------------------------------------------------
/// | | | |
/// | | | |--------------------------------------------------------------
/// | | |
/// | | | |--------------------------------------------------------------  PhaseCorScans   ..................................
/// | | | |
/// | | | |--------------------------------------------------------------
/// | | |
/// | | | |--------------------------------------------------------------  kSpaceLoop   .....................................
/// | | | |
/// | | | | |------------------------------------------------------------  Begin: ScanLoop 0 ... 255   ......................
/// | | | | |
/// | | | | | |----------------------------------------------------------  Begin: AcquisitionLoop 0 ... 0   .................
/// | | | | | |
/// | | | | | | |--------------------------------------------------------  StdIR   ..........................................
/// | | | | | | |
/// | | | | | | | |------------------------------------------------------  Begin: PhaseLoop 0 ... 0   .......................
/// | | | | | | | |
/// | | | | | | | | |----------------------------------------------------  Begin: InnerSliceLoop 0 ... 0   ..................
/// | | | | | | | | |
/// | | | | | | | | | |--------------------------------------------------  SBB_Pre_MP   .....................................
/// | | | | | | | | | |
/// | | | | | | | | | |--------------------------------------------------
/// | | | | | | | | |
/// | | | | | | | | | |--------------------------------------------------  Kernel_MP   ......................................
/// | | | | | | | | | |
/// | | | | | | | | | | |----------------------------------------|         ..................................................
/// | | | | | | | | | | | MedicKernel                            |
/// | | | | | | | | | | |----------------------------------------|
/// | | | | | | | | | |
/// | | | | | | | | | |--------------------------------------------------
/// | | | | | | | | |
/// | | | | | | | | | |----------------------------------------|           ..................................................
/// | | | | | | | | | | TRFill                                 |
/// | | | | | | | | | |----------------------------------------|
/// | | | | | | | | |
/// | | | | | | | | | |--------------------------------------------------  SBB_Post_MP   ....................................
/// | | | | | | | | | |
/// | | | | | | | | | |--------------------------------------------------
/// | | | | | | | | |
/// | | | | | | | | |----------------------------------------------------
/// | | | | | | | |
/// | | | | | | | | |----------------------------------------|             ..................................................
/// | | | | | | | | | TRFillEnd                              |
/// | | | | | | | | |----------------------------------------|
/// | | | | | | | |
/// | | | | | | | |------------------------------------------------------
/// | | | | | | |
/// | | | | | | |--------------------------------------------------------
/// | | | | | |
/// | | | | | |----------------------------------------------------------
/// | | | | |
/// | | | | |------------------------------------------------------------
/// | | | |
/// | | | |--------------------------------------------------------------
/// | | |
/// | | |----------------------------------------------------------------
/// | |
/// | | |----------------------------------------|                         ..................................................
/// | | | MeasRepDelay                           |
/// | | |----------------------------------------|
/// | |
/// | |------------------------------------------------------------------
/// |
/// |--------------------------------------------------------------------
/// </pre>
///
/// Sequential multi slice: \n 
///
/// <pre>
/// |--------------------------------------------------------------------  StdSL   ..........................................
/// |
/// | |--------------------------------------------------------U---------  Begin: MeasurementLoop 0 ... 0   .................
/// | |
/// | | |----------------------------------------|                         ..................................................
/// | | | TokTokTok                              |
/// | | |----------------------------------------|
/// | |
/// | | |----------------------------------------------------------------  Begin: InnerSliceLoop 0 ... 0   ..................
/// | | |
/// | | | |--------------------------------------------------------------  Begin: PreparingScans 0 ... 0   ..................
/// | | | |
/// | | | | |------------------------------------------------------------  StdIR   ..........................................
/// | | | | |
/// | | | | | |----------------------------------------------------------  Begin: PhaseLoop 0 ... 0   .......................
/// | | | | | |
/// | | | | | | |----------------------------------------|                 ..................................................
/// | | | | | | | TRFillEnd                              |
/// | | | | | | |----------------------------------------|
/// | | | | | |
/// | | | | | | |--------------------------------------------------------  SBB_Pre_MP   .....................................
/// | | | | | | |
/// | | | | | | |--------------------------------------------------------
/// | | | | | |
/// | | | | | | |--------------------------------------------------------  Kernel_MP   ......................................
/// | | | | | | |
/// | | | | | | | |----------------------------------------|               ..................................................
/// | | | | | | | | MedicKernel                            |
/// | | | | | | | |----------------------------------------|
/// | | | | | | |
/// | | | | | | |--------------------------------------------------------
/// | | | | | |
/// | | | | | | |----------------------------------------|                 ..................................................
/// | | | | | | | TRFill                                 |
/// | | | | | | |----------------------------------------|
/// | | | | | |
/// | | | | | | |--------------------------------------------------------  SBB_Post_MP   ....................................
/// | | | | | | |
/// | | | | | | |--------------------------------------------------------
/// | | | | | |
/// | | | | | |----------------------------------------------------------
/// | | | | |
/// | | | | |------------------------------------------------------------
/// | | | |
/// | | | |--------------------------------------------------------------
/// | | |
/// | | | |--------------------------------------------------------------  PhaseCorScans   ..................................
/// | | | |
/// | | | |--------------------------------------------------------------
/// | | |
/// | | | |--------------------------------------------------------------  kSpaceLoop   .....................................
/// | | | |
/// | | | | |------------------------------------------------------------  Begin: ScanLoop 0 ... 255   ......................
/// | | | | |
/// | | | | | |----------------------------------------------------------  Begin: AcquisitionLoop 0 ... 0   .................
/// | | | | | |
/// | | | | | | |--------------------------------------------------------  StdIR   ..........................................
/// | | | | | | |
/// | | | | | | | |------------------------------------------------------  Begin: PhaseLoop 0 ... 0   .......................
/// | | | | | | | |
/// | | | | | | | | |----------------------------------------|             ..................................................
/// | | | | | | | | | TRFillEnd                              |
/// | | | | | | | | |----------------------------------------|
/// | | | | | | | |
/// | | | | | | | | |----------------------------------------------------  SBB_Pre_MP   .....................................
/// | | | | | | | | |
/// | | | | | | | | |----------------------------------------------------
/// | | | | | | | |
/// | | | | | | | | |----------------------------------------------------  Kernel_MP   ......................................
/// | | | | | | | | |
/// | | | | | | | | | |----------------------------------------|           ..................................................
/// | | | | | | | | | | MedicKernel                            |
/// | | | | | | | | | |----------------------------------------|
/// | | | | | | | | |
/// | | | | | | | | |----------------------------------------------------
/// | | | | | | | |
/// | | | | | | | | |----------------------------------------|             ..................................................
/// | | | | | | | | | TRFill                                 |
/// | | | | | | | | |----------------------------------------|
/// | | | | | | | |
/// | | | | | | | | |----------------------------------------------------  SBB_Post_MP   ....................................
/// | | | | | | | | |
/// | | | | | | | | |----------------------------------------------------
/// | | | | | | | |
/// | | | | | | | |------------------------------------------------------
/// | | | | | | |
/// | | | | | | |--------------------------------------------------------
/// | | | | | |
/// | | | | | |----------------------------------------------------------
/// | | | | |
/// | | | | |------------------------------------------------------------
/// | | | |
/// | | | |--------------------------------------------------------------
/// | | |
/// | | |----------------------------------------------------------------
/// | |
/// | | |----------------------------------------|                         ..................................................
/// | | | MeasRepDelay                           |
/// | | |----------------------------------------|
/// | |
/// | |------------------------------------------------------------------
/// |
/// |--------------------------------------------------------------------
/// </pre>
///
/// \par Default Values
/// - \c SLK::e_PERF_TOK_TOK_TOK        Controls the execution of the TokTokTok. \n
///                                     Allowed values: \n 
///                                     - \c SL::ALWAYS     Perform TokTokTok prior 
///                                                         to each measurement
///                                     - \c SL::ONLY_FIRST_REPETITION  Perform 
///                                                         TokTokTok prior to the 
///                                                         first measurement, only
///                                     - \c SL::NEVER      Do not perform TokTokTok
/// - \c SLK::b_SINGLE_SHOT_TRIGGER     Flag determining the position of the trigger
///                                     - \c true   intended for single shot sequences.
///                                                 The trigger is executed once per 
///                                                 slice
///                                     - \c false  default, intended for none single 
///                                                 shot sequences. The trigger is 
///                                                 executed prior to the phases loop.
/// - \c SLK::b_INTERLEAVED_IR          Flag determining whether Interleaved IR or 
///                                     standard IR has to be used
/// - \c SLK::b_PERF_RF_SPOILED_SATS    Determines whether saturation pulses should use
///                                     rf-spoiling
/// - \c SLK::b_SPOIL_BEFORE_KERNEL     Determines whether a spoiling gradient should be 
///                                     apllied prior to the kernel
/// - \c SLK::b_PERF_SATS_FOR_CHECK     Determines whether saturation pulses should be
///                                     used during sequence check
///
/// \par Exported Values
/// - \c SLK::b_PERF_TR_FILL            Determines whether an OscBit event has to be 
///                                     executed
/// - \c SLK::b_SPOIL_AFTER_KERNEL      Determines whether a spoiling gradient should be 
///                                     apllied after to the kernel
/// - \c SLK::b_PERF_NOISE_MEAS         Controls the execution of the noise measurement 
///                                     that is used for iPAT image reconstruction
/// - \c SLK::l_ACQ_LOOP_LENGTH         Number of acquisitions, taken from the protocol
/// - \c SLK::l_MEASUREMENT_LOOP_LENGTH Number of measuremnts, value is retrieved 
///                                     from the protocol
/// - \c SLK::l_PHASE_LOOP_LENGTH       Number of phases, taken from the protocol
/// - \c SLK::l_CONCAT_LOOP_COUNTER     Current concatenation counter
///
/// \par Imported Values
/// - \c SLK::b_PERF_NOISE_MEAS         Controls the execution of the noise measurement 
///                                     that is used for iPAT image reconstruction
/// - \c SLK::b_INTERLEAVED_IR          Flag determining whether Interleaved IR or 
///                                     standard IR has to be used
/// - \c SLK::b_PERF_SATS_FOR_CHECK     Determines whether saturation pulses should be
///                                     used during sequence check
/// - \c SLK::e_PERF_TOK_TOK_TOK        Controls the execution of the TokTokTok. \n
///                                     Allowed values: \n 
///                                     - \c SL::ALWAYS     Perform TokTokTOk prior 
///                                                         to each measurement
///                                     - \c SL::ONLY_FIRST_REPETITION  Perform 
///                                                         TokTokTok prior to the 
///                                                         first measurement, only
///                                     - \c SL::NEVER      Do not perform TokTokTok
/// - \c SLK::e_PERF_PAT_REF_SCAN       Controls the execution of the PAT reference 
///                                     scan. \n
///                                     Allowed values: \n 
///                                     - \c SL::ALWAYS     Perform PAT reference scan 
///                                                         prior to each measurement
///                                     - \c SL::ONLY_FIRST_REPETITION  Perform 
///                                                         PAT reference scan prior 
///                                                         to the first measurement, 
///                                                         only
///                                     - \c SL::NEVER      Do not perform PAT reference 
///                                                         scan
/// - \c SLK::e_PERF_PREP_SCAN          Controls the execution of the preparing scans. \n
///                                     Allowed values: \n 
///                                     - \c SL::ALWAYS     Perform preparing scans 
///                                                         prior to each measurement
///                                     - \c SL::ONLY_FIRST_REPETITION  Perform 
///                                                         preparing scans prior 
///                                                         to the first measurement, 
///                                                         only
///                                     - \c SL::NEVER      Do not perform preparing
///                                                         scans
///
// ----------------------------------------------------------------------------
class  AslSL : public RootNode 
{
    public:


        // * ------------------------------------------------------------------ *
        // * Maximum number of RSats                                            *
        // * ------------------------------------------------------------------ *
        enum  { lMaxNoSat = 8 };


        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::AslSL
        //
        // Description :
        ///              \brief     Creates all composite elements that are required
        ///                         to build the standard Composite SeqLoop
        ///
        ///                         - m_pNoiseMeas              Noise measurement used for
        ///                                                     iPAT reconstruction
        ///                         - m_pMeasurementLoop        Measurement loop
        ///                         - m_pTokTokTok              Introduction, series of three
        ///                                                     gradient pulses
        ///                         - m_pConcatLoop             Concatenation loop
        ///                         - m_pInnerSliceLoop         Slice loop
        ///                         - m_pMeasRepDelay           Delay between consecutive 
        ///                                                     measurementrts
        ///                         - m_pPATRefScan             PAT external reference scans
        ///                         - m_pPreparingScans         Preparing scans to drive 
        ///                                                     magnetization onto steady-state
        ///                         - m_pPhaseCorScans          Phase correction scans
        ///                         - m_pAcquisitionLoop        Acquisition loop
        ///                         - m_pKSpaceLoop             Either nested line and partiton
        ///                                                     loops or a single k-space scan 
        ///                                                     loop
        ///                         - m_pStdIR                  Module for standard IR 
        ///                                                     measurements
        ///                         - m_pInterleavedIR          Module for interleaved IR 
        ///                                                     measurements
        ///                         - m_pMSat                   Magnetization transfer pulse
        ///                         - m_pCSatFat                Fat saturation pulse
        ///                         - m_pRFS_CSatFat            Fat saturation pulse using 
        ///                                                     rf spoiling
        ///                         - m_pCSatWater              Water saturation pulse
        ///                         - m_pRFS_CSatWater          Water saturation pulse using
        ///                                                     rf spoiling
        ///                         - m_pRSat                   Array of regional saturation
        ///                                                     pulses
        ///                         - m_pRFS_RSat               Array of regional saturation
        ///                                                     pulses that use rf spoiling
        ///                         - m_pSpoilerBeforeKernel    Spoiler element that is 
        ///                                                     applied before the SBBs /
        ///                                                     kernel
        ///                         - m_pSpoilerAfterKernel     Spoiler element that is 
        ///                                                     applied after the SBBs /
        ///                                                     kernel
        //
        // Parameters  :
        ///               \param    sIdent              String identifier of the leaf
        ///               \param    pMediator           Pointer to Mediator class 
        //
        // Return      :
        ///              \return    void
        // 
        // ------------------------------------------------------------------------------
        AslSL (const std::string& sIdent, Mediator* pMediator);

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::AslSL
        //
        // Description :
        ///              \brief     Creates all composite elements that are required
        ///                         to build the standard Composite SeqLoop
        ///
        ///                         - m_pNoiseMeas              Noise measurement used for
        ///                                                     iPAT reconstruction
        ///                         - m_pMeasurementLoop        Measurement loop
        ///                         - m_pTokTokTok              Introduction, series of three
        ///                                                     gradient pulses
        ///                         - m_pConcatLoop             Concatenation loop
        ///                         - m_pInnerSliceLoop         Slice loop
        ///                         - m_pMeasRepDelay           Delay between consecutive 
        ///                                                     measurementrts
        ///                         - m_pPATRefScan             PAT external reference scans
        ///                         - m_pPreparingScans         Preparing scans to drive 
        ///                                                     magnetization onto steady-state
        ///                         - m_pPhaseCorScans          Phase correction scans
        ///                         - m_pAcquisitionLoop        Acquisition loop
        ///                         - m_pKSpaceLoop             Either nested line and partiton
        ///                                                     loops or a single k-space scan 
        ///                                                     loop
        ///                         - m_pStdIR                  Module for standard IR 
        ///                                                     measurements
        ///                         - m_pInterleavedIR          Module for interleaved IR 
        ///                                                     measurements
        ///                         - m_pMSat                   Magnetization transfer pulse
        ///                         - m_pCSatFat                Fat saturation pulse
        ///                         - m_pRFS_CSatFat            Fat saturation pulse using 
        ///                                                     rf spoiling
        ///                         - m_pCSatWater              Water saturation pulse
        ///                         - m_pRFS_CSatWater          Water saturation pulse using
        ///                                                     rf spoiling
        ///                         - m_pRSat                   Array of regional saturation
        ///                                                     pulses
        ///                         - m_pRFS_RSat               Array of regional saturation
        ///                                                     pulses that use rf spoiling
        ///                         - m_pSpoilerBeforeKernel    Spoiler element that is 
        ///                                                     applied before the SBBs /
        ///                                                     kernel
        ///                         - m_pSpoilerAfterKernel     Spoiler element that is 
        ///                                                     applied after the SBBs /
        ///                                                     kernel
        //
        // Parameters  :
        ///              \param     sScope      String identifier of the scope of this element
        ///               \param    sIdent      String identifier of the leaf
        ///               \param    pMediator   Pointer to Mediator class 
        //
        // Return      :
        ///              \return    void
        // 
        // ------------------------------------------------------------------------------
        AslSL (const std::string& sScope, const std::string& sIdent, Mediator* pMediator);

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::~AslSL
        //
        // Description :
        ///              \brief     Deletes all elements that have been created by the
        ///                         constructor
        //
        // Parameters  :
        ///              \param     -
        //
        // Return      :
        ///              \return    void
        // 
        // ------------------------------------------------------------------------------
        virtual ~AslSL ();

        void createObjects(const std::string& sScope, const std::string& sIdent, Mediator* pMediator);

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::registerDefaultValues
        //
        // Description :
        ///              \brief     Registers the default values of the standard loop
        //
        // Parameters  :
        ///              \param     rMrProt     Reference to the protocol class
        //
        // Return      :
        ///              \return    
        ///                         \arg \c     true  : successful execution \n
        ///                         \arg \c     false : error occurred
        // 
        // ------------------------------------------------------------------------------
        virtual bool registerDefaultValues (const MrProt& rMrProt);

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::configure
        //
        // Description :
        ///              \brief     Configuration of AslSL
        ///
        ///                         This function registers the following values:
        ///
        ///                             \li \c SLK::b_PERF_TR_FILL
        ///                             \li \c SLK::b_SPOIL_AFTER_KERNEL
        ///                             \li \c SLK::b_PERF_NOISE_MEAS
        ///                             \li \c SLK::l_CONCAT_LOOP_LENGTH
        ///                             \li \c SLK::l_ACQ_LOOP_LENGTH
        ///                             \li \c SLK::l_MEASUREMENT_LOOP_LENGTH
        ///                             \li \c SLK::l_PHASE_LOOP_LENGTH
        ///
        ///                         and expands the standard tree.
        //
        // Parameters  :
        ///              \param     rMrProt     Reference to the protocol class
        ///              \param     rSeqLim     Reference to the sequence limits class
        //
        // Return      :
        ///              \return    \arg \c     true  : successful execution \n
        ///                         \arg \c     false : error occurred
        // 
        // ------------------------------------------------------------------------------
        virtual bool configure (const MrProt& rMrProt, const SeqLim& rSeqLim);

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::expand
        //
        // Description :
        ///              \brief     Expands the loop standard ASL Composite SeqLoop structure.
        //
        // Parameters  :
        ///              \param     rMrProt     Reference to the protocol class
        //
        // Return      :
        ///              \return    \arg \c     true  : loop structure successfully 
        ///                                             expanded \n
        ///                         \arg \c     false : error occurred
        // 
        // ------------------------------------------------------------------------------
        virtual IComponent* expand (const MrProt& rMrProt);

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::check
        //
        // Description :
        ///              \brief     Executes the loop structure for the GSWD check and 
        ///                         maximum gradient amplitude check
        ///
        ///                         The scan mode of all elements is stored and set to 
        ///                         scan mode CHECK before execution. Afterwars the 
        ///                         previous scan mode is rstored.
        //
        // Parameters  :
        ///              \param     rMrProt     Reference to the protocol structure
        ///              \param     rSeqLim     Reference to the sequence limits structure
        ///              \param     rSeqExpo    Reference to the sequence exports structure
        ///              \param     pSLC        Pointer to the rotation matrix and slice
        ///                                     position  information
        //
        // Return      :
        ///              \return   
        ///                         \arg \c     true  : successful execution \n
        ///                         \arg \c     false : error occurred
        // 
        // ------------------------------------------------------------------------------
        virtual bool check (MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo, sSLICE_POS* pSLC);

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::getbPerformTRFill
        //
        // Description :
        ///              \brief     Returns whether TR fill times will be executed by
        ///                         the composite seqloop or by the kernel.
        //
        // Parameters  :
        ///              \param     -
        //
        // Return      :
        ///              \return    \arg \c     true  : TR fill time will be performed 
        ///                                             StdSL \n
        ///                         \arg \c     false : TR fill time has to be performed
        ///                                             by the kernel element
        // 
        // ------------------------------------------------------------------------------
        bool getbPerformTRFill (void);

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::getDurationOfMeasurement
        //
        // Description :
        ///              \brief     Returns the duration of the specified measurement
        ///                         in us
        //
        // Parameters  :
        ///              \param     lCounter    Index of the measurement of interest
        //
        // Return      :
        ///              \return    Duration of the specified measurement in us
        // 
        // ------------------------------------------------------------------------------
        virtual int64_t getDurationOfMeasurement (int32_t lCounter);

        // ------------------------------------------------------------------------------
        //
        // Name        :            StdSL::calculateTotalMeasurementTime
        //
        // Description :
        ///              \brief     Returns the total duration and considers special timing 
        ///                         
        //
        // Parameters  :
        ///              \param     rMrProt     protocol
        //
        // Return      :
        ///              \return    total measurement time to be inserted in SeqExpo
        // 
        // ------------------------------------------------------------------------------
        virtual int64_t calculateTotalMeasurementTime(MrProt &rMrProt);

        // ------------------------------------------------------------------------------
        //
        // Name        :            StdSL::getTimeToFirstImagAcq
        //
        // Description :
        ///              \brief     Relevant for TimeToCentre calclulations. Returns the
        ///                         time spend in nodes (e.g. TokTokTok, noiseMeas, etc.)
        ///                         before the first imaging acquisition.
        //
        //
        // Return      :
        ///              \return    Time to first k-space acquisition
        // 
        // ------------------------------------------------------------------------------
        virtual int64_t getTimeToFirstImagAcq();

        // ------------------------------------------------------------------------------
        //
        // Name        :            AslSL::addASLComponents (SUPPORT_CSL)
        //
        // Description :
        ///              \brief     Convenient method, that adds the Nodes for the ASL preparation as well as the readout to the loop structure.
        //
        // Parameters  :
        ///              \param     preparation    Pointer to the preparation LeafAdapter
        ///              \param     readout    Pointer to the readout LeafAdapter
        //
        // Return      :
        ///              \return    Duration of the specified measurement in us
        // 
        // ------------------------------------------------------------------------------
        void addASLComponents   (IComponent* preparation, IComponent* readout);
        void removeASLComponents();

    protected:
        // --------------------------------------------------------------------
        /// M0Scan loop (SUPPORT_CSL)
        // --------------------------------------------------------------------
        M0Scans* m_pM0Scans;

        // --------------------------------------------------------------------
        /// NoiseMeas
        // --------------------------------------------------------------------
        NoiseMeas* m_pNoiseMeas;

        // --------------------------------------------------------------------
        /// Measurement loop
        // --------------------------------------------------------------------
        MeasurementLoop*  m_pMeasurementLoop;

        // --------------------------------------------------------------------
        /// TokTokTok
        // --------------------------------------------------------------------
        TokTokTok* m_pTokTokTok;

        // --------------------------------------------------------------------
        /// Concatenation loop
        // --------------------------------------------------------------------
        ConcatLoop* m_pConcatLoop;

        // --------------------------------------------------------------------
        /// Multi breathhold; synchronization with breathing commands
        // --------------------------------------------------------------------
        BreathHoldSync* m_pBreathHoldSync;

        // --------------------------------------------------------------------
        /// Slice loop used for sequential multislice imaging
        // --------------------------------------------------------------------
        InnerSliceLoop* m_pInnerSliceLoop;

        // --------------------------------------------------------------------
        /// MeasRepDelay
        // --------------------------------------------------------------------
        MeasRepDelay* m_pMeasRepDelay;

        // --------------------------------------------------------------------
        /// External PAT reference scans
        // --------------------------------------------------------------------
        ASLPATRefScan* m_pPATRefScan;   // SUPPORT_CSL

        // --------------------------------------------------------------------
        /// Preparation scans
        // --------------------------------------------------------------------
        PreparingScans* m_pPreparingScans;

        // --------------------------------------------------------------------
        /// Phase correction scans
        // --------------------------------------------------------------------
        PhaseCorrScans* m_pPhaseCorrScans;

        // --------------------------------------------------------------------
        /// Acquisition loop
        // --------------------------------------------------------------------
        ConstLoop* m_pAcquisitionLoop;

        // --------------------------------------------------------------------
        /// Scan loop
        /// This composite expands either to a pair lines/partition loops or
        /// to a single 3D k-space loop
        // --------------------------------------------------------------------
        kSpaceLoop* m_pKSpaceLoop;

        // --------------------------------------------------------------------
        /// Composite for standard inversion
        // --------------------------------------------------------------------
        StdIRASL* m_pStdIR;   // ECG_TRIGGERING

        // --------------------------------------------------------------------
        /// Composite for interleaved inversion
        // --------------------------------------------------------------------
        InterleavedIR* m_pInterleavedIR;

        // --------------------------------------------------------------------
        /// Composite for looping over AslLabelState(s) (SUPPORT_CSL)
        // --------------------------------------------------------------------
        ASLLoop* m_pAslLoop;

        // --------------------------------------------------------------------
        /// Composite for looping over segments. (SUPPORT_CSL)
        // --------------------------------------------------------------------
        ConstLoop * m_pSegmentsLoop;

        // --------------------------------------------------------------------
        /// Magnetization transfer pulse
        // --------------------------------------------------------------------
        MSat* m_pMSat;

        // --------------------------------------------------------------------
        /// Fat saturation pulse
        // --------------------------------------------------------------------
        CSat* m_pCSatFat;

        // --------------------------------------------------------------------
        /// Fat saturation pulse using rf spoiling
        // --------------------------------------------------------------------
        RFSpoiledLeaf* m_pRFS_CSatFat;

        // --------------------------------------------------------------------
        /// Water saturation pulse
        // --------------------------------------------------------------------
        CSat* m_pCSatWater;

        // --------------------------------------------------------------------
        /// Water saturation pulse using rf spoiling
        // --------------------------------------------------------------------
        RFSpoiledLeaf* m_pRFS_CSatWater;

        // --------------------------------------------------------------------
        /// RSat
        // --------------------------------------------------------------------
        std::vector<RSat*> m_RSat;

        // --------------------------------------------------------------------
        /// RSat using rf spoiling
        // --------------------------------------------------------------------
        std::vector<RFSpoiledLeaf*> m_RFS_RSat;

        // --------------------------------------------------------------------
        /// Spoiler element applied before the kernel
        // --------------------------------------------------------------------
        SpoilGrad* m_pSpoilerBeforeKernel;

        // --------------------------------------------------------------------
        /// Spoiler element applied after the kernel
        // --------------------------------------------------------------------
        SpoilGrad* m_pSpoilerAfterKernel;

        // --------------------------------------------------------------------
        /// fill time between concatenations
        // --------------------------------------------------------------------
        ConcatFillTime* m_pConcatFillTime;

    private:
        // --------------------------------------------------------------------
        // Copy constructor not implemented
        // --------------------------------------------------------------------
        AslSL (const AslSL& right);

        // --------------------------------------------------------------------
        // Assignment operator not implemented
        // --------------------------------------------------------------------
        AslSL & operator=(const AslSL &right);
};

} // namespace SL_NAME_SPACE


