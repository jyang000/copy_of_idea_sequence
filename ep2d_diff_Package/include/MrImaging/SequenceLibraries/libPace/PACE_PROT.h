
//  -----------------------------------------------------------------
//  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential
//  -----------------------------------------------------------------
//
//          Project: NUMARIS/4
//             File: \n4_servers1\pkg\MrServers\MrImaging\seq\common\libPace\PACE_PROT.h
//          Version: \main\27
//
//         Language: C++
//
//      Description: Wrapper for PACE parameter access.
//          Structs: 
//
//  -----------------------------------------------------------------

#ifndef __PACE_PROT_H
#define __PACE_PROT_H

//  -----------------------------------------------------------------
//  Used Interfaces

//  Definition of int32_t
#include "MrGlobalDefinitions/MrBasicTypes.h"
//  Definition of ICE_PROGRAM_PARA_LEN
#include "MrProtSrv/Domain/CoreNative/MeasAqcDefs.h"

//  Forward declarations
#ifdef USE_XPROTOCOL
class XProtocol;
#else
namespace MrProtocolData
{
    class MrProtData;
}
#endif // USE_XPROTOCOL

//  Definition of enum PulseMode
#include "MrProtSrv/Domain/CoreNative/NavigatorDefines.h"

class MdhProxy;
class Navigator;

#ifdef __IMP_EXP
#undef __IMP_EXP
#endif 

#ifdef LOCAL_PACE_BUILD
#define __IMP_EXP
#pragma message( "Local PACE build" )
#else  //   LOCAL_PACE_BUILD not defined
#ifdef BUILD_libPACE
#define __OWNER
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control
#else 
#ifdef BUILD_libSL
#undef __OWNER
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control
#else
#define __IMP_EXP 
#endif
#endif

#endif  //  LOCAL_PACE_BUILD


namespace PACE
{
#ifdef USE_XPROTOCOL
    typedef XProtocol PROT_TYPE;
#else
    typedef MrProtocolData::MrProtData PROT_TYPE;
#endif        
        //  ScoutMode
        enum Send
        {
            SEND_ROI_PROJ_SUM_LIN     = 0x0
        ,   SEND_NAVIGATOR_IMAGE      = 0x1
        ,   SEND_ROI_PROJ_SUM_COL     = 0x2
        //  Overlay data only
        ,   SEND_PHYSIO_SIGNAL        = 0x4
        ,   SEND_NOTHING              = 0x8

        };
        
        //  Motion direction
        enum MotionDir
        {
            MOTION_DIR_COL            = 0x0
        ,   MOTION_DIR_INPLANE        = 0x1
        };
        
        //  Motion detection algorithm
        enum MoDeAlgo
        {
            MODE_ALGO_PATTERN_MATCH   = 0x0
        ,   MODE_ALGO_GAUSS_NEWTON    = 0x1   
        ,   MODE_ALGO_PATTERN_CMP     = 0x2
        };
        
        enum PaceTriggerAlgo
        {
            // The standard trigger algo uses a acceptance window in which breathing position is detected
            PACE_TRIGGER_ALGO_STANDARD = 0x00

            // this algo triggers only on a level.
        ,   PACE_TRIGGER_ALGO_TRIGGER_LEVEL = 0x01
        };

        enum AcqWindowDuration
        {
            // the default settings: ACQ window is its actual duration
            ACQ_DURATION_ACTUAL = 0x00

            // the ACQ window duration is prolonged synthetically in order to center the TE in the window
        ,   ACQ_DURATION_CENTER_TE = 0x01
        };


        enum DiffMeas
        {
            DIFF_MEAS_NCC             = 0x0
        ,   DIFF_MEAS_SSD             = 0x1                
        ,   DIFF_MEAS_ZNCC            = 0x2
        };

        enum Extract
        {
            EXTRACT_MAG               = 0x0
        ,   EXTRACT_PHS               = 0x1
        };
        
        //  Window function
        enum Filter
        {
            FILTER_OFF                = 0x0
        ,   FILTER_MEAN               = 0x1
        ,   FILTER_GAUSS              = 0x2
        };
        
        //  Simulation mode
        enum SimMode
        {
            SIM_MODE_OFF              = 0
        ,   SIM_MODE_SIN_PI_2         = 1
        ,   SIM_MODE_SIN_PI_4         = 2
        ,   SIM_MODE_SIN_PI_8         = 3
        ,   SIM_MODE_SIN_PI_16        = 4
        };
        
        enum ExcitMode
        {
            EXCIT_MODE_2D_PENCIL      = SEQ::EXCIT_MODE_2D_PENCIL
        ,   EXCIT_MODE_GRE            = SEQ::EXCIT_MODE_GRE
        ,   EXCIT_MODE_SLICE          = SEQ::EXCIT_MODE_SLICE
        };
        
        
        enum Seed
        {
            SEED_MOST_PROBABLE                = 0x0
                
            //  Left-most,Posterior-most, or head-most
        ,   SEED_LPH_MOST                     = 0x1
        
            //  Right-most,Anterior-most, or feet-most
        ,   SEED_RAF_MOST                     = 0x2
        
            //  Percental value between LPH-most and RAF-most.
        ,   SEED_PERCENTAL                    = 0x4
        };
        
        enum Reorder
        {
            REORDER_LIN_ASC                   = 0x0
        ,   REORDER_INV_CENTRIC               = 0x1
        };
        
        
        enum TriggerPoint
        {
            //  The respiratory cycle
            TRIGGER_POINT_MASK_RESP_CYCLE                      = 0xf0
                
            //  Expiration
        ,   TRIGGER_POINT_MASK_RESP_CYCLE_EXPIRATION           = 0x00
                
            //  Inspiration
        ,   TRIGGER_POINT_MASK_RESP_CYCLE_INSPIRATION          = 0x10

            //  Tries to find the quieter phase
        ,   TRIGGER_POINT_MASK_RESP_CYCLE_AUTO                 = 0x20
            
            //  Trigger condition
        ,   TRIGGER_POINT_MASK_TRIGGER_COND                    = 0x0f
            
            //  Absolute diaphragm position
        ,   TRIGGER_POINT_MASK_TRIGGER_COND_ABS_DIAPHRAGM_POS  = 0x00
            
            //  Relative diaphragm position in respiratory cycle
        ,   TRIGGER_POINT_MASK_TRIGGER_COND_REL_DIAPHRAGM_POS  = 0x01
            
            //  Delay after end inspiration/expiration
        ,   TRIGGER_POINT_MASK_TRIGGER_COND_DELAY_AFTER_END_IE = 0x02
                
        };
        
        enum DblTriggerMode
        {
            DBLTRIGGER_MODE_SIMULTANEOUSLY                      = 0x0
        ,   DBLTRIGGER_MODE_CARDIAC_AFTER_RESP                  = 0x1
        };
        
        enum QRS_DETECTION
        {
            //  Note: Enumeration name 'QRSDetection' is conform
            //  with naming convention but conflicts with method
            //  QRSDetection.
            
            //  Standard trigger of PMU software
            QRS_DETECTION_STANDARD               = 0x0
            //  "Engelse and Zeelenberg" 
        ,   QRS_DETECTION_DF1                    = 0x1
        };
        
        enum DVAMode
        {
            //  (Reacq)uire (R)ight-(A)nterior-(F)eet-most
            DVA_MODE_REACQ_RAF_MOST              = 0x0
                
            //  (Reacq)uire (L)eft-(P)osterior-(H)ead-most
        ,   DVA_MODE_REACQ_LPH_MOST              = 0x1
                
            //  (Reacq)uire (M)aximum (A)bsolute (D)eviate 
            //  from (Median) diaphragm position
        ,   DVA_MODE_REACQ_MAD_MEDIAN            = 0x2
                
        };
        
        enum GatingAlgorithm
        {
            GATING_ALGORITHM_PAWS_B2                     = 0x00
        ,   GATING_ALGORITHM_PAWS_B3                     = 0x01
        //  (A)utomatic (M)ultiple (W)indow (S)election
        ,   GATING_ALGORITHM_AMWS                        = 0x20
        //  (G)ating (A)fter (T)rigger
        ,   GATING_ALGORITHM_GAT                         = 0x40
        ,
        };
        
        enum TOM
        {
            TOM_MINIMIZE_TE                      = 0x00
        ,   TOM_IN_PHASE                         = 0x01
        ,   TOM_OPPOSED_PHASE                    = 0x02
        ,   TOM_OPPOSED_IN_PHASE                 = 0x04
        ,   TOM_IN_OPPOSED_PHASE                 = 0x10
        ,   TOM_OFF                              = 0x20
        };
        
        enum ExtraSwitch
        {
            EXTRA_SWITCH_MASK_SHOW_CORR_SLICE_DATA              = 0x0000001
        ,   EXTRA_SWITCH_MASK_STORE_NAVIGATOR_IMA               = 0x0000002
        ,   EXTRA_SWITCH_MASK_LOG_FILE                          = 0x0000004
        ,   EXTRA_SWITCH_MASK_VERBOSE                           = 0x0000008
        ,   EXTRA_SWITCH_MASK_REC_ENTIRE_IMA                    = 0x0000100
        ,   EXTRA_SWITCH_MASK_SERIES_INTERL_DESCENDING          = 0x0000200
        ,   EXTRA_SWITCH_MASK_SHOW_OUTER_NAV                    = 0x0000400
        ,   EXTRA_SWITCH_MASK_PTTHRESHOLD_ADAPTIVE              = 0x0000800
        ,   EXTRA_SWITCH_MASK_CLX_COMBINE                       = 0x0001000
        ,   EXTRA_SWITCH_MASK_CLX_COMB_WITH_PHASE_CORR          = 0x0002000
        ,   EXTRA_SWITCH_MASK_MANUAL_TRIGGER_THRESHOLD          = 0x0040000
        ,   EXTRA_SWITCH_MASK_SEPARATE_STAT_PHASE               = 0x0080000
        ,   EXTRA_SWITCH_MASK_TRIGGER_ADAPTIVE_LEVEL            = 0x0100000        
        ,   EXTRA_SWITCH_MASK_TRIGGER_ADAPTIVE_OFF              = 0x0200000
        ,   EXTRA_SWITCH_MASK_CONTINUOUS_MODE                   = 0x0400000
        ,   EXTRA_SWITCH_MASK_INVERT_SIGNAL                     = 0x0800000
        ,   EXTRA_SWITCH_MASK_CONTINUOUS_SAT                    = 0x1000000
        ,   EXTRA_SWITCH_MASK_STEADY_STATE_PREP_ON              = 0x2000000
        ,   EXTRA_SWITCH_MASK_NOISE_ADJ_SCAN                    = 0x4000000
        ,   EXTRA_SWITCH_MASK_NOISE_DECORRELATION_OFF           = 0x8000000
        };

        enum SelectAcqWindow
        {
            SELECT_ACQ_WINDOW_AUTOMATIC   = SEQ::SELECT_ACQ_AUTOMATIC
        ,   SELECT_ACQ_WINDOW_MANUAL      = SEQ::SELECT_ACQ_MANUAL
        };

        enum PositionNavigator
        {
            POSITION_NAVIGATOR_AUTOMATIC  = SEQ::POSITION_NAVIGATOR_AUTOMATIC
        ,   POSITION_NAVIGATOR_MANUAL     = SEQ::POSITION_NAVIGATOR_MANUAL
        };

        //  
        enum
        {
            //  Free int32_t parameters            
            FREELNG_INDX_B0_SCOUT_BASE_RES             = 0,
            FREELNG_SBIT_B0_SCOUT_BASE_RES             = 0,
            FREELNG_INDX_GRE_SCOUT_BASE_RES            = 0,
            FREELNG_SBIT_GRE_SCOUT_BASE_RES            = 16,
            FREELNG_INDX_NAV_SCOUT_BASE_RES            = 1,
            FREELNG_SBIT_NAV_SCOUT_BASE_RES            = 0,
            FREELNG_INDX_NAV_SCOUT_PREP_SCANS          = 1,
            FREELNG_SBIT_NAV_SCOUT_PREP_SCANS          = 16,


            FREELNG_INDX_GRE_SCOUT_DWELL_TIME          = 2,
            FREELNG_INDX_SLEEP_INTERV                  = 3,
            FREELNG_INDX_POLL_INTERV                   = 4,
            FREELNG_INDX_EXTRA_SWITCH                  = 5,
            FREELNG_INDX_GRE_SCOUT_TRECOVERY_BEFORE_SCOUT = 6,
            FREELNG_INDX_B0_SCOUT_DWELL_TIME           = 7,
            //  16 bit each
            FREELNG_INDX_B0_SCOUT_NECO                 = 8,
            FREELNG_SBIT_B0_SCOUT_NECO                 = 0,

            FREELNG_INDX_GRE_SCOUT_NECO                = 8,
            FREELNG_SBIT_GRE_SCOUT_NECO                = 16,

            FREELNG_INDX_B0_SCOUT_TRECOVERY_BEFORE_SCOUT = 9,
            FREELNG_INDX_B0_SCOUT_PREP_SCANS             = 10,
            FREELNG_SBIT_B0_SCOUT_PREP_SCANS             = 0,
            FREELNG_INDX_GRE_SCOUT_PREP_SCANS            = 10,
            FREELNG_SBIT_GRE_SCOUT_PREP_SCANS            = 16,
            
            FREELNG_INDX_SIM_MODE                      = 11,
            FREELNG_SBIT_SIM_MODE                      = 16,
            
            FREELNG_INDX_MODE_ALGO                     = 12,
            FREELNG_SBIT_MODE_ALGO                     = 0,
            FREELNG_INDX_ACQ_DURATION                  = 12,
            FREELNG_SBIT_ACQ_DURATION                  = 24,
            FREELNG_INDX_DIFF_MEAS                     = 12,
            FREELNG_SBIT_DIFF_MEAS                     = 8,
            FREELNG_INDX_EXTRACT                       = 12,
            FREELNG_SBIT_EXTRACT                       = 16,
            
            FREELNG_INDX_SEND                          = 13,
            FREELNG_SBIT_SEND                          = 0,
            FREELNG_INDX_SEND_CHA                      = 13,
            FREELNG_SBIT_SEND_CHA                      = 8,
            FREELNG_INDX_SEARCH_SEEDS                  = 13,
            FREELNG_SBIT_SEARCH_SEEDS                  = 16,
            
            FREELNG_INDX_COARSEST_RES_COL              = 14,
            FREELNG_SBIT_COARSEST_RES_COL              = 0,
            FREELNG_INDX_FINEST_RES_COL                = 14,
            FREELNG_SBIT_FINEST_RES_COL                = 8,
            
            FREELNG_INDX_COARSEST_RES_LIN              = 14,
            FREELNG_SBIT_COARSEST_RES_LIN              = 16,
            FREELNG_INDX_FINEST_RES_LIN                = 14,
            FREELNG_SBIT_FINEST_RES_LIN                = 24,
            
            FREELNG_INDX_STATISTIC_SCANS               = 15,
            FREELNG_INDX_NAV_SCOUT_DWELL_TIME          = 16,

            FREELNG_INDX_B0_SCOUT_PE_LIN               = 17,
            FREELNG_SBIT_B0_SCOUT_PE_LIN               = 0,
            FREELNG_INDX_B0_SCOUT_PEFT_LEN             = 17,
            FREELNG_SBIT_B0_SCOUT_PEFT_LEN             = 16,

            
            FREELNG_INDX_B0_SCOUT_TRIGGER_POINT        = 18,
            FREELNG_SBIT_B0_SCOUT_TRIGGER_POINT        = 0,            
            
            FREELNG_INDX_GRE_SCOUT_TRIGGER_POINT       = 18,
            FREELNG_SBIT_GRE_SCOUT_TRIGGER_POINT       = 8,

            FREELNG_INDX_NAV_SCOUT_TRIGGER_POINT       = 18,
            FREELNG_SBIT_NAV_SCOUT_TRIGGER_POINT       = 16,
            
            FREELNG_INDX_DBLTRIGGER_MODE               = 18,
            FREELNG_SBIT_DBLTRIGGER_MODE               = 24,

            FREELNG_INDX_GATING_ALGORITHM              = 19,
            FREELNG_SBIT_GATING_ALGORITHM              = 0,

            FREELNG_INDX_GATING_EXPLOIT_RESP_PHASE     = 19,
            FREELNG_SBIT_GATING_EXPLOIT_RESP_PHASE     = 8,

            
            FREELNG_INDX_GRE_SCOUT_PE_LIN              = 20,
            FREELNG_SBIT_GRE_SCOUT_PE_LIN              = 0,
            FREELNG_INDX_GRE_SCOUT_PEFT_LEN            = 20,
            FREELNG_SBIT_GRE_SCOUT_PEFT_LEN            = 16,

            FREELNG_INDX_NAV_SCOUT_TRAVERSE            = 21,
            FREELNG_SBIT_NAV_SCOUT_TRAVERSE            = 0,

            FREELNG_INDX_NAV_SCOUT_ROTATION             = 21,
            FREELNG_SBIT_NAV_SCOUT_ROTATION             = 8,

            FREELNG_INDX_NAV_SCOUT_TRAJECTORY           = 21,
            FREELNG_SBIT_NAV_SCOUT_TRAJECTORY           = 16,

            FREELNG_INDX_NAV_SCOUT_PROFILE              = 21,
            FREELNG_SBIT_NAV_SCOUT_PROFILE              = 24,

            FREELNG_INDX_NAV_SCOUT_TRECOVERY_BEFORE_SCOUT = 22,
            FREELNG_INDX_TRIGGER_PULSE                  = 23,
            
            FREELNG_INDX_GRE_SCOUT_NRF                  = 24,
            FREELNG_SBIT_GRE_SCOUT_NRF                  = 0,

            FREELNG_INDX_B0_SCOUT_NRF                   = 24,
            FREELNG_SBIT_B0_SCOUT_NRF                   = 16,

            
            //  Reference signal passed to ICE
            FREELNG_INDX_REF_SIGNAL                    = 25,

            FREELNG_INDX_STATISTIC_DURATION            = 26,

            FREELNG_INDX_B0_SCOUT_TOM                  = 27,
            FREELNG_SBIT_B0_SCOUT_TOM                  = 0,

            FREELNG_INDX_B0_SCOUT_FLAGS                = 27,
            FREELNG_SBIT_B0_SCOUT_FLAGS                = 8,
            FREELNG_MASK_B0_SCOUT_VEL_COMP             = 0x1,
            FREELNG_MASK_B0_SCOUT_BIPOLAR              = 0x2,
            FREELNG_MASK_B0_SCOUT_SHORTEST_DWELL_TIME  = 0x4,
            FREELNG_MASK_B0_SCOUT_AUTOCOILSELECT       = 0x10,
            FREELNG_MASK_B0_SCOUT_RF_SPOILING          = 0x20,


            FREELNG_INDX_GRE_SCOUT_TOM                 = 27,
            FREELNG_SBIT_GRE_SCOUT_TOM                 = 16,

            FREELNG_INDX_GRE_SCOUT_FLAGS               = 27,
            FREELNG_SBIT_GRE_SCOUT_FLAGS               = 24,
            FREELNG_MASK_GRE_SCOUT_VEL_COMP            = 0x1,
            FREELNG_MASK_GRE_SCOUT_BIPOLAR             = 0x2,
            FREELNG_MASK_GRE_SCOUT_SHORTEST_DWELL_TIME = 0x4,
            FREELNG_MASK_GRE_SCOUT_AUTOCOILSELECT      = 0x10,
            FREELNG_MASK_GRE_SCOUT_RF_SPOILING         = 0x20,

            FREELNG_INDX_PTTHRESHOLD_FOR_PTTHRESHOLD_ADAPTIVE = 28, 
            FREELNG_INDX_IDLE_TIME_FOR_TRIGGER_ADAPTIVE = 29, 


            FREELNG_INDX_B0_SCOUT_RF_DURATION          = 31,
            FREELNG_INDX_GRE_SCOUT_RF_DURATION         = 32,

            FREELNG_INDX_TRECOVERY_CONTINUOUS_MODE     = 33,

            FREELNG_INDX_B0_SCOUT_SCOUT_TR             = 34,

            FREELNG_INDX_NAV_SCOUT_SCOUT_TR            = 35,


            //  enum { MaxSizeOf_alFree = 36 };

            //  Free double parameters
            FREEDBL_INDX_GRE_SCOUT_FLIP_ANGLE          = 0,
            FREEDBL_INDX_NAV_SCOUT_FLIP_ANGLE          = 1,
            FREEDBL_INDX_B0_SCOUT_FLIP_ANGLE           = 2,
            FREEDBL_INDX_PHASE_OS                      = 3,
            FREEDBL_INDX_TRIGGER_LOCK_TIME             = 4,
            FREEDBL_INDX_MAX_SIM_TRA                   = 5,
            FREEDBL_INDX_PT_THRESHOLD                  = 6,
            FREEDBL_INDX_SENSITIVITY                   = 7,
            FREEDBL_INDX_IDLE_TIME                     = 8,
            FREEDBL_INDX_PT_THRESHOLD_LEVEL             = 9,
            FREEDBL_INDX_FFT_SCALE                     = 11,
            FREEDBL_INDX_TRIGGER_THRESHOLD             = 12,
            FREEDBL_INDX_MIN_SEND_INTRV                = 13,
            FREEDBL_INDX_GRE_SCOUT_SPOILER_PHASE_DISPERSION = 14,
            FREEDBL_INDX_PT_THRESHOLD_REL              = 15,
            FREEDBL_INDX_ACCEPT_WND_REL                = 16,

            FREEDBL_INDX_ACQ_WINDOW_REL                = 17,

            FREEDBL_INDX_B0_SCOUT_SPOILER_PHASE_DISPERSION = 18,

            FREEDBL_INDEX_STEADY_STATE_PREP_T1             = 19,

            FREEDBL_INDX_B0_SCOUT_BWTP                     = 20,
            FREEDBL_INDX_GRE_SCOUT_BWTP                    = 21,
            FREEDBL_INDX_NAV_SCOUT_RESOLUTION              = 22,
            FREEDBL_INDX_SIM_DEVIATION                     = 23,

            
            //  enum { MaxSizeOf_adFree = 24 };
            
            //  Free SeqExpo parameters
            //  Echo time of the first navigator echo
            SEQ_EXPO_FB_PARA_TE1_us                            = 0,
            //  Echo spacing, if navigator acquires multiple echoes
            SEQ_EXPO_FB_PARA_DELTA_TE_us                       = 1,
            
            SEQ_EXPO_FB_PARA_OFFSET                            = 4,
            
            //  Receives enum SEQ::GradientAxis that specifies the readout (COL)-axis
            //  of the first navigator object.
            SEQ_EXPO_FB_PARA_INDX_COL_DIR                      = 0+SEQ_EXPO_FB_PARA_OFFSET,
            
            //  Receives enum SEQ::GradientAxis that specifies the phase-encoding (LIN)-axis
            //  of the first navigator object.
            SEQ_EXPO_FB_PARA_INDX_LIN_DIR                      = 1+SEQ_EXPO_FB_PARA_OFFSET,
            
            
            //  Receives the number N of feedback channels used.
            SEQ_EXPO_FB_PARA_INDX_FBCHA_SIZE                   = 2+SEQ_EXPO_FB_PARA_OFFSET,
            
            //  At most three channel types:
            
            //  Configures FB channel 0 
            SEQ_EXPO_FB_PARA_INDX_CONF_FBCHA0                  = 3+SEQ_EXPO_FB_PARA_OFFSET,
            
            //  Configures FB channel 1, if N > 1.
            SEQ_EXPO_FB_PARA_INDX_CONF_FBCHA1                  = 4+SEQ_EXPO_FB_PARA_OFFSET,
            
            //  Configures FB channel 2, ..., N-1, if N > 2.
            SEQ_EXPO_FB_PARA_INDX_CONF_FBCHA2                  = 5+SEQ_EXPO_FB_PARA_OFFSET,
            
            //  Configures (F)unction (BLOCK) used to (REC)onstruct (SCOUT) (IMA)ges.
            SEQ_EXPO_FB_PARA_SBIT_FBLOCK_REC                   = 0,
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_REC_DEF               = 0x1,
            //SEQ_EXPO_FB_PARA_MASK_FBLOCK_REC_EXT               = 0x2,
            //SEQ_EXPO_FB_PARA_MASK_FBLOCK_REC_ECG               = 0x4,
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_REC_B0S               = 0x8,
            
            //  Configures function block used to (FILTER) the reconstructed scout images.
            SEQ_EXPO_FB_PARA_SBIT_FBLOCK_FILTER                = 8,
            //  Standard image filter (e.g. GAUSS)
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_FILTER_IMA            = 0x1,
            //  (E)dge (ENHANCE)ement
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_ENHANCE               = 0x2,
            
            //  Configures function block used for (MO)tion (DE)tection.
            SEQ_EXPO_FB_PARA_SBIT_FBLOCK_MODE                  = 16,
            //  (P)attern (MATCH) algorithm.
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_MODE_PMATCH           = 0x1,
            
            
            //  Configures function block used for (EVAL)uation.
            SEQ_EXPO_FB_PARA_SBIT_FBLOCK_EVAL                  = 24,
            //  (TRIGGER) algorithm
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_EVAL_TRIGGER          = 0x1,
            //  (R)espiratory (O)rdered (P)hase (E)ncoding
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_EVAL_ROPE             = 0x2,
            
            //  (G)ating (A)fter (T)rigger
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_EVAL_GAT              = 0x4,
            
            //  (P)ure (F)ollowing
            SEQ_EXPO_FB_PARA_MASK_FBLOCK_EVAL_FOLLOW           = 0x10,
            
            
            //  Channel specific data in the range [._BEGIN,._END[
            SEQ_EXPO_FB_PARA_INDX_CPARA_BEGIN                  = 6+SEQ_EXPO_FB_PARA_OFFSET,
            SEQ_EXPO_FB_PARA_INDX_CPARA_END                    = ICE_PROGRAM_PARA_LEN,
            
            //  ID of the channel specific data. 
            //  Statistic time window
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_STAT_TIME_WINDOW  = 0x1,
            
            //  Maximum number of reacquire phases.
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_MAX_REACQ_PHASES  = 0x2,
            
            //  Receives the total number of respiartory trigger halts
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_TRIGGER_HALTS     = 0x3,
            
            //  Receives the trigger delay [ms]
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_TRIGGER_DELAY_MS  = 0x4,
            
            //  Receives acquisition period [ms]
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_ACQ_PERIOD_MS     = 0x5,
            
            //  Actual acquisition (quite) period
            //  = AcqPeriod_ms for simple triggered measurement
            //  = AcqPeriod_ms + AcqPeriod_RR * AverageRR for a double
            //    triggered measurement.
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_ACQ_PERIOD_RR     = 0x6,
            
            //  Receives # profiles with different priority
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_NPRIO             = 0x7,
            
            //  Receives position of profile with highest priority
            //  within [0,NPRIO[
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_POS_PRIO_ZERO     = 0x8,
            
            //  Receives the number of profiles (lines times partitions)
            //  acquired after each navigator
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_PROFILES_PER_NAV  = 0x9,
            
            //  Time between adjacent physiological signal points. 
            //  If not set Scout TR is used.
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_SAMPL_INTRV_MS       = 0xa,

            //  Only if gating algorithm GAT is used:
            //  Number of respiratory phases per trigger
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_GAT_NRESP_PHASE    = 0xb,

            //  The respiratory phase with (h)ighest (p)riority.
            SEQ_EXPO_FB_PARA_MASK_CPARA_ID_GAT_CRESP_PHASE_HP = 0xc
                
    };
    
    //  -------------------------------------------------------------
    //  IScout static member functions
    //  -------------------------------------------------------------
    //  Set/retrieve function bloc from int32_t variable that is used to
    //  configure the FB channel.
    int32_t&         FBLOCK_REC(int32_t& rlDst, unsigned char ucMask);
    unsigned char FBLOCK_REC(int32_t lSrc);
    
    int32_t&         FBLOCK_FILTER(int32_t& rlDst, unsigned char ucMask);
    unsigned char FBLOCK_FILTER(int32_t lSrc);
    
    int32_t&         FBLOCK_MODE(int32_t& rlDst, unsigned char ucMask);
    unsigned char FBLOCK_MODE(int32_t lSrc);
    
    int32_t&         FBLOCK_EVAL(int32_t& rlDst, unsigned char ucMask);
    unsigned char FBLOCK_EVAL(int32_t lSrc);
    
    //  Set channel specific parameter.
    //
    //  Layout:  |      24|      16|       8|       0|
    //           +--------+--------+--------+--------+
    //           |  CHA   |   ID   |    V A L U E    |
    //           +--------+--------+-----------------+
    //  Range:     0...255  0...255         0...65535
    //
    bool PACK_CPARA(int32_t& rlDst, int32_t lCha, int32_t lID, int32_t lVal);
    __IMP_EXP void UNPACK_CPARA(int32_t& rlCha, int32_t& rlID, int32_t& lVal, int32_t lSrc);
    

    

/*    
    //  The static member function returns the FTLen nearest to dFTLen that
    //  avoids the ChirpZ algorithm.
    static inline int32_t findFTLen(double dFTLen);
    
    //  The static member function returns an integer that is the multiple of
    //  nMult closest to dIn.
    static inline int32_t multipleOf(double dIn, int32_t nMult);
    
    //  Masks bit flag of parameter pProt->navigatorPara().alFree[FREELNG_INDX_EXTRA_SWITCH]
    static inline bool getExtraFlag(const PROT_TYPE* pProt, uint32_t ulMask);
    
    //  Sets bit ON/OFF masked by ulMask in pProt->navigatorPara().alFree[FREELNG_INDX_EXTRA_SWITCH]
    static inline bool setExtraFlag(PROT_TYPE* pProt, uint32_t ulMask, bool bNewVal);
    
    //  The first static member function masks the bit flags ulMask of the free int32_t
    //  parameter with index FREELNG_INDX_VALID_FLAG
    //  The second static member function sets bit ON/OFF masked by ulMask.
    static inline bool getValidFlag(const PROT_TYPE* pProt, int32_t i32Mask);
    static inline bool setValidFlag(PROT_TYPE* pProt, uint32_t ulMask, bool bNewVal);
    
    
    //  Returns the 8 bits from position iStart
    uint8_t getBits(const PROT_TYPE* pProt, int32_t iPos, int32_t iStart);
       */
 
    
    
    
    //  -------------------------------------------------------------
    //  Note: If one of the following static member functions is needed,
    //        the file "Measurement/Sequence/Prot/MrProt.h" must precede
    //        IScout.h in the clients include hierachy.
    //  -------------------------------------------------------------
    
    //static inline uint8_t version(const PROT_TYPE* pProt);
    //static inline Version setVersion(PROT_TYPE* pProt);
    
    //static inline uint8_t context(const PROT_TYPE* pProt);
    //static inline uint8_t context(PROT_TYPE* pProt, uint8_t  ucNewVal);
  
    //  The first static member function returns the
    //  'base resolution' of the navigator image.
    //  The second static member function is the associated set function.
    int16_t baseRes(const PROT_TYPE* pProt);
    int16_t baseRes(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t B0ScoutBaseRes(const PROT_TYPE* pProt);
    int16_t B0ScoutBaseRes(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t GREScoutBaseRes(const PROT_TYPE* pProt);
    int16_t GREScoutBaseRes(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t NavScoutBaseRes(const PROT_TYPE* pProt);
    int16_t NavScoutBaseRes(PROT_TYPE* pProt, int16_t i16NewVal);


    //  The static member function returns the number of image lines.
    int32_t imageLines(const PROT_TYPE* pProt);
    
    //  The static member function returns the number of image columns.
    int32_t imageColumns(const PROT_TYPE* pProt);

    
    //  The first static member function returns the # of phase encoding
    //  lines of the navigator image. The second static member function
    //  is the associated set function.
    int16_t B0ScoutNPELin(const PROT_TYPE* pProt);
    int16_t B0ScoutNPELin(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t GREScoutNPELin(const PROT_TYPE* pProt);
    int16_t GREScoutNPELin(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t NPELin(const PROT_TYPE* pProt);
    int16_t NPELin(PROT_TYPE* pProt, int16_t i16NewVal);
    
    //  The first static member function returns the number of 
    //  preparing scans. The second static member function
    //  is the associated set function.
    int16_t GREScoutPrepScans(const PROT_TYPE* pProt);
    int16_t GREScoutPrepScans(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t B0ScoutPrepScans(const PROT_TYPE* pProt);
    int16_t B0ScoutPrepScans(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t NavScoutPrepScans(const PROT_TYPE* pProt);
    int16_t NavScoutPrepScans(PROT_TYPE* pProt, int16_t i16NewVal);
    
    //  The static member function returns the OS-factor in phase encoding
    //  direction.
    double phaseOS(const PROT_TYPE* pProt);
    double phaseOS(PROT_TYPE* pProt, double dNewVal);

    //  The first function returns the phase dispersion of the spoiler gradient
    //  at the end of the scout sequence along the slice selection direction 
    //  in units of PI
    //  The second function is the associated set-function.
    double GREScoutSpoilerPhaseDispersion_pi(const PROT_TYPE* pProt);
    double GREScoutSpoilerPhaseDispersion_pi(PROT_TYPE* pProt, double dNewVal);

    double B0ScoutSpoilerPhaseDispersion_pi(const PROT_TYPE* pProt);
    double B0ScoutSpoilerPhaseDispersion_pi(PROT_TYPE* pProt, double dNewVal);

    
    
    //  Initializes navigator object from XProtocol
    //  Fails if index lPos is invalid.
    //  The static member calculates the size and the upper left corner of the RoI.
    //  The unit is Pixel.
    bool calcImaWindow(const PROT_TYPE* pProt, int32_t& rlStartCol, int32_t& rlStartLin, int32_t& rlSizeCol, int32_t& rlSizeLin, uint32_t u32Nav, uint32_t u32Wnd);

    bool navAt     (Navigator& rDst, const PROT_TYPE* pProt, int32_t i32Pos);
    bool navPhFoVAt(double& rdPhFoV, const PROT_TYPE* pProt, int32_t i32Pos);
    bool navRdFoVAt(double& rdRdFoV, const PROT_TYPE* pProt, int32_t i32Pos);
    bool navThickAt(double& rdThick, const PROT_TYPE* pProt, int32_t i32Pos);
    
    
    //  The static member function calculates the image pixel size in
    //  column and lines direction and writes the result to
    //  rdPixelSizeCol_mm/double& rdPixelSizeLin_mm respectively.
    //static inline bool calcPixelSize(const PROT_TYPE* pProt, double& rdPixelSizeCol_mm, double& rdPixelSizeLin_mm, int32_t lNavPos);
    
    //  The static member function returns the PEFTLen.
    int16_t B0ScoutPEFTLen(const PROT_TYPE* pProt);
    int16_t B0ScoutPEFTLen(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t GREScoutPEFTLen(const PROT_TYPE* pProt);
    int16_t GREScoutPEFTLen(PROT_TYPE* pProt, int16_t i16NewVal);

    int16_t PEFTLen(const PROT_TYPE* pProt);
    int16_t PEFTLen(PROT_TYPE* pProt, int16_t i16NewVal);
    
    //  The static member function returns the real (including oversampling) dwell time in [ns]
    //  used for the readout of the GREScout.
    int32_t GREScoutDwellTime_ns(const PROT_TYPE* pProt);
    int32_t GREScoutDwellTime_ns(PROT_TYPE* pProt, int32_t lNewVal_ns);
    
    //  The static member function returns the real (including oversampling) dwell time in [ns]
    //  used for the readout of the B0-Scout.
    int32_t B0ScoutDwellTime_ns(const PROT_TYPE* pProt);
    int32_t B0ScoutDwellTime_ns(PROT_TYPE* pProt, int32_t lNewVal_ns);

    //  The static member function returns tpencil beam.
    int32_t NavScoutDwellTime_ns(const PROT_TYPE* pProt);
    int32_t NavScoutDwellTime_ns(PROT_TYPE* pProt, int32_t lNewVal_ns);


    //  The static member function returns/sets the duartion of the excitation
    //  pulse in us. The value must be an integer multiple of 20 us.
    int32_t GREScoutRFDuration_us(const PROT_TYPE* pProt);
    int32_t GREScoutRFDuration_us(PROT_TYPE* pProt, int32_t i32NewVal_us);

    int32_t B0ScoutRFDuration_us(const PROT_TYPE* pProt);
    int32_t B0ScoutRFDuration_us(PROT_TYPE* pProt, int32_t i32NewVal_us);

  
    //  The static member function returns/sets the flip angle
    //  of the excitation pulse in [deg].
    double GREScoutFlipAngle_deg(const PROT_TYPE* pProt);
    double GREScoutFlipAngle_deg(PROT_TYPE* pProt, double dNewVal_deg);

    double NavScoutFlipAngle_deg(const PROT_TYPE* pProt);
    double NavScoutFlipAngle_deg(PROT_TYPE* pProt, double dNewVal_deg);

    double B0ScoutFlipAngle_deg(const PROT_TYPE* pProt);
    double B0ScoutFlipAngle_deg(PROT_TYPE* pProt, double dNewVal_deg);

    //  The static member function returns/sets the bandwidth-time
    //  product of the excitation pulse.
    double GREScoutBWTP(const PROT_TYPE* pProt, double dLarmorConst_MHz_T);
    double GREScoutBWTP(PROT_TYPE* pProt, double dLarmorConst_MHz_T, double dNewVal);

    double B0ScoutBWTP(const PROT_TYPE* pProt, double dLarmorConst_MHz_T);
    double B0ScoutBWTP(PROT_TYPE* pProt, double dLarmorConst_MHz_T, double dNewVal);


    //  If liver dome scout is used:
    //  The static member function returns an enumeration which defined whether
    //  magnitude (EXTRACT_MAG) or phase images (EXTRACT_PHS) are calculated from
    //  the navigator data.
    uint8_t extract(const PROT_TYPE* pProt);
    uint8_t extract(PROT_TYPE* pProt, uint8_t u8NewVal);

    //  The function returns the search window 'u32Pos' in column direction. 
    //  The second static member function is the associated set function.
    double maxTraCol_mm(const PROT_TYPE* pProt);
    double maxTraCol_mm(PROT_TYPE* pProt, double dNewVal_mm);

    //  The function returns the search window 'u32Pos' in line direction. 
    //  The second static member function is the associated set function.
    double maxTraLin_mm (const PROT_TYPE* pProt);
    double maxTraLin_mm (PROT_TYPE* pProt, double dNewVal_mm);

    //  The function returns maxTraCol(Lin)_mm(pProt,lPos) / (Pixel Size Col(Lin)[mm])
    int32_t maxTraCol_pix(const PROT_TYPE* pProt);
    int32_t maxTraLin_pix(const PROT_TYPE* pProt);


    //  The static member function returns the difference measure used
    //  for the pattern match algorithm.
    uint8_t PACETriggerAlgo(const PROT_TYPE* pProt);
    uint8_t PACETriggerAlgo(PROT_TYPE* pProt, unsigned char ucNewVal);
    
    //  GRE-Scout only:
    //  The first static member function returns a boolean value that
    //  indicates whether the PT-threshold is adapted to the breathing pattern
    //  of the patient. In this case the  PTThreshold_rel must be specified.
    //  The second static member function is the related set function.
    bool PTThresholdAdaptive(const PROT_TYPE* pProt);
    bool PTThresholdAdaptive(PROT_TYPE* pProt, bool bNewVal);


    //  PTThreshold for adaptive
    //  Return value devided by 10 specifies the factor by which the median respiratory cycle
    //  is multiplied to calculate the idle time (time since last trigger event) after which the
    //  trigger level is adapted
    int32_t idleTimeForTriggerAdaptive_x10(const PROT_TYPE* pProt);
    int32_t idleTimeForTriggerAdaptive_x10(PROT_TYPE* pProt, int32_t iNewVal_x10);

    //  The static member function returns the difference measure used
    //  for the pattern match algorithm.
    uint8_t ACQDuration(const PROT_TYPE* pProt);
    AcqWindowDuration ACQDuration(PROT_TYPE* pProt, unsigned char ucNewVal);
    
    //  The static member function returns the difference measure used
    //  for the pattern match algorithm.
    uint8_t diffMeas(const PROT_TYPE* pProt);
    DiffMeas diffMeas(PROT_TYPE* pProt, unsigned char ucNewVal);
    
    //  Parameters for the coarse to fine minimization of the pattern match
    //  algorithm:
    
    //  The finest resolution in column direction is 2ª [pix]. Thereby a is
    //  the return value of the static member function.
    char finestResCol(const PROT_TYPE* pProt);
    char finestResCol(PROT_TYPE* pProt, char cNewVal);
    
    
    //  The coarsest resolution in column direction is 2ª [pix]. Thereby a is
    //  the return value of the static member function.
    char coarsestResCol(const PROT_TYPE* pProt);
    char coarsestResCol(PROT_TYPE* pProt, char cNewVal);
    
    //  The finest resolution in lines direction is 2ª [pix]. Thereby a is
    //  the return value of the static member function.
    char finestResLin(const PROT_TYPE* pProt);
    char finestResLin(PROT_TYPE* pProt, char cNewVal);
    
    //  The coarsest resolution in lines direction is 2ª [pix]. Thereby a is
    //  the return value of the static member function.
    char coarsestResLin(const PROT_TYPE* pProt);
    char coarsestResLin(PROT_TYPE* pProt, char cNewVal);
    
    //  The member function returns the final resolution in column
    //  direction in units of image pixels;
    double finalResCol_pix(const PROT_TYPE* pProt);
    
    
    //  The static member function returns the number of search seeds evaluated
    //  by the coarse-to-fine pattern match algorithm.
    int16_t searchSeeds(const PROT_TYPE* pProt);
    short searchSeeds(PROT_TYPE* pProt, short shNewVal);

    //  The static member function returns an enum type that specifies the image
    //  to be send.
    unsigned char send(const PROT_TYPE* pProt);
    unsigned char send(PROT_TYPE* pProt, unsigned char ucNewVal);
    
    //  The static member function returns the scout channel that sends
    //  images to the Online display.
    unsigned char sendCha(const PROT_TYPE* pProt);
    unsigned char sendCha(PROT_TYPE* pProt, unsigned char ucNewVal);
        
    //  The static member function returns a boolean value indicating whether the
    //  detected motion is stored on disc.
    bool logFile(const PROT_TYPE* pProt);
    bool logFile(PROT_TYPE* pProt, bool bNewVal);
    
    //  The static member function returns a boolean value indicating whether the
    //  outer navigator is displayed in the GSP.
    //static inline bool showOuterNav(const PROT_TYPE* pProt);
    //static inline bool showOuterNav(PROT_TYPE* pProt, bool bNewVal);
    
    //  If the static member function returns true diagnostic run-time information,
    //  is written to the MPCU-trace-logger.
    bool verbose(const PROT_TYPE* pProt);
    bool verbose(PROT_TYPE* pProt, bool bNewVal);

    //  The static member function returns a boolean value indicating whether it
    //  is required to reconstruct the entire image (even when this is not necessary.
    bool recEntireIma(const PROT_TYPE* pProt);
    bool recEntireIma(PROT_TYPE* pProt, bool bNewVal);

    //  The static member function returns a boolean value indicating whether
    //  the first navigator object is visualized in GSP.
    bool SHOW_OUTER_NAV(const PROT_TYPE* pProt);
    bool SHOW_OUTER_NAV(PROT_TYPE* pProt, bool bNewVal);

    
    //  If the static member function returns true, the raw data of the individual
    //  channels are added prior to RoFT. If the return value is false, separate
    //  images for each channel are calculated and added afterwords:
    //  Ima_comb(x,y) ~ sqrt( sum_c FFTScaleFact_c * IMA_c(x,y,c)^2 )
    bool clxCombine(const PROT_TYPE* pProt);
    bool clxCombine(PROT_TYPE* pProt, bool bNewVal);

    
    //  Like clxCombine. In addition an individual phase factor for each channel
    //  is determined in an initial phase correction scan. clxCombine and 
    //  clxCombWithPhaseCorr are mutually exclusive.
    bool clxCombWithPhaseCorr(const PROT_TYPE* pProt);
    bool clxCombWithPhaseCorr(PROT_TYPE* pProt, bool bNewVal);
    
    
    //  After the scout scan has been send to the DSP's the sequence stops feeding 
    //  the DSP's with meaningful data as int32_t as the feedback is pending.
    //  It feeds them with a series of synchronization events, each followed
    //  by a pause time. Whenever the GSU sends the echo of the synchronization event
    //  it is checked whether the feedback is still pending. If so another WAKEUP/Pause
    //  interval is inserted, otherwise the sequence breaks out of the loop. 
    //  The return value of the static member function determines the time of a single
    //  pause interval in [us]. The duration of the first pause interval is returned by
    //  sleepIntervAfterScout_us. To retrieve the duration off all following intervals
    //  pollIntervAfterScout_us must be used. If one of the pause 
    //  intervals is choosen too short a FIFO-underflow error may occur. A int32_t pause time
    //  worsen the time resolution. 
    int32_t sleepIntervAfterScout_us(const PROT_TYPE* pProt);
    int32_t pollIntervAfterScout_us (const PROT_TYPE* pProt);
    int32_t sleepIntervAfterScout_us(PROT_TYPE* pProt, int32_t i32NewVal_us);
    int32_t pollIntervAfterScout_us (PROT_TYPE* pProt, int32_t i32NewVal_us);
    
    //  To measure the position of the diaphragm a navigator, perpendicular to the image
    //  plane, is placed through the dome of the hemidiaphragm. The measured respiratory
    //  displacement of the diaphragm is used to adapt the slice position.
    //  Since the anatomy under study moves in general less than the diaphragm, an empirical
    //  correction factor must be used for the determination of the corrected slice position.
    double corrFactCol(const PROT_TYPE* pProt);
    double corrFactLin(const PROT_TYPE* pProt);
    double corrFactCol(PROT_TYPE* pProt, double dNewVal);
    double corrFactLin(PROT_TYPE* pProt, double dNewVal);

    //  In the simulation mode the navigator image is acquired at a position that is shifted
    //  compared with the position stored in the protocol. This simulates motion.
    //  The static member function returns an enum that describes the motion.
    uint8_t simMode(const PROT_TYPE* pProt);
    uint8_t simMode(PROT_TYPE* pProt, uint8_t  ui8NewVal);
    
    //  The static member function returns the maximum simulated translation [mm].
    double maxSimTra_mm(const PROT_TYPE* pProt);
    double maxSimTra_mm(PROT_TYPE* pProt, double dNewVal_mm);
    
    //  The static member function returns the deviation of the breathing period for the simulation cycle.
    double SimCycleDev(const PROT_TYPE* pProt);
    double SimCycleDev(PROT_TYPE* pProt, double dNewVal);
   
    /* The static member function returns the factor <sensitivity> to control triggering in case of breathing deviations
    Short cycle mode :
    The time between successive triggers is limited by the protocol TR to avoid that
        SAR is exceeded during online tracing.This can shift the trigger towards the start of
        inspiration when the respiratory cycle is in the order of the protocol TR.
        To avoid this the WIP now skips every second cycle if the following condition is
    fulfilled :
    Median respiratory cycle - <sensitivity> × standard deviation < TR
    */
    double SensitivityForShortCycleMode(const PROT_TYPE* pProt);
    double SensitivityForShortCycleMode(PROT_TYPE* pProt, double dNewVal);
    
    //  PTThreshold for adaptive
    //  Return value specifies the factor by which the median respiratory cycle
    //  is multiplied to calculate the idle time (time since last trigger event) after which the
    //  trigger level is adapted
    double IdleTimeForTriggerAdaptive(const PROT_TYPE* pProt);
    double IdleTimeForTriggerAdaptive(PROT_TYPE* pProt, double dNewVal);

    //  The static member function returns an enum that specifies the desired
    //  navigator. Supported values are SEQ::EXCIT_MODE_2D_PENCIL, SEQ::EXCIT_MODE_GRE, SEQ::EXCIT_MODE_SLICE
    int32_t  excitMode(const PROT_TYPE* pProt);
    int32_t  excitMode(PROT_TYPE* pProt, int32_t i32NewVal);

    //  The static member function returns a boolean value, that indicates whether
    //  the profile images are stored as an extra series in the patient's data base.
    bool storeProfileIma(const PROT_TYPE* pProt);
    bool storeProfileIma(PROT_TYPE* pProt, bool bNewVal);

    //  The static member function returns an enum that specifies the method
    //  used to suppress artifacts caused by respiratory motion.
    __IMP_EXP int32_t respComp(const PROT_TYPE* pProt);
    int32_t respComp(PROT_TYPE* pProt, int32_t i32NewVal);
    

    //  Returns the repetition time of the navigator during scout mode/trigger phase [us].
    int32_t GREScoutScoutTR_us(const PROT_TYPE* pProt);
    int32_t GREScoutScoutTR_us(PROT_TYPE* pProt, int32_t i32NewVal_us);
    int32_t B0ScoutScoutTR_us(const PROT_TYPE* pProt);
    int32_t B0ScoutScoutTR_us(PROT_TYPE* pProt, int32_t i32NewVal_us);
    int32_t NavScoutScoutTR_us(const PROT_TYPE* pProt);
    int32_t NavScoutScoutTR_us(PROT_TYPE* pProt, int32_t i32NewVal_us);
    __IMP_EXP int32_t scoutTR_us(const PROT_TYPE* pProt);
    int32_t scoutTR_us(PROT_TYPE* pProt, int32_t lNewVal_us);

    //  The function returns a boolean value, that indicates whether
    //  the sequence runs in scout mode (no real image are acquired) or not.
    __IMP_EXP bool scoutMode(const PROT_TYPE* pProt);
    bool scoutMode(PROT_TYPE* pProt, bool bNewVal);
        
    //  Returns the number of measurements of a scout protocol. 
    __IMP_EXP int32_t scoutDuration_us(const PROT_TYPE* pProt);
    int32_t scoutDuration_us(PROT_TYPE* pProt, int32_t lNewVal_us);
    
    //  The static member function returns the minimum recovery time between
    //  an imaging scan and the application of the RF-pulse. This time delay
    //  is necessary to avoid pattern mismatch due to saturation effects.
    int32_t B0ScoutTRecoveryBeforeScout_us(const PROT_TYPE* pProt);
    int32_t B0ScoutTRecoveryBeforeScout_us(PROT_TYPE* pProt, int32_t i32NewVal_us);
    int32_t GREScoutTRecoveryBeforeScout_us(const PROT_TYPE* pProt);
    int32_t GREScoutTRecoveryBeforeScout_us(PROT_TYPE* pProt, int32_t i32NewVal_us);
    int32_t NavScoutTRecoveryBeforeScout_us(const PROT_TYPE* pProt);
    int32_t NavScoutTRecoveryBeforeScout_us(PROT_TYPE* pProt, int32_t i32NewVal_us);


    //  Recovery time, if continuous mode is used
    __IMP_EXP int32_t TRecoveryContinuousMode_us(const PROT_TYPE* pProt);
    int32_t TRecoveryContinuousMode_us(PROT_TYPE* pProt, int32_t i32NewVal_us);

    
    //  The static member function returns the acceptace window in COL/LIN
    //  direction [mm].
    double   acceptanceWindowCol_mm  (const PROT_TYPE* pProt);
    double   acceptanceWindowLin_mm  (const PROT_TYPE* pProt);
    int32_t acceptanceWindowCol_pix (const PROT_TYPE* pProt);
    int32_t acceptanceWindowLin_pix (const PROT_TYPE* pProt);
    double   dAcceptanceWindowCol_pix(const PROT_TYPE* pProt);
    double   dAcceptanceWindowLin_pix(const PROT_TYPE* pProt);
    
    
    double acceptanceWindowCol_mm(PROT_TYPE* pProt, double dNewVal_mm);
    double acceptanceWindowLin_mm(PROT_TYPE* pProt, double dNewVal_mm);


    //  If B0-Scout is used: the acceptance window as fractional value 
    //  between the minimum and maximum value of the physiological signal
    //  found during the statistic phase.
    double acceptanceWindow_rel(const PROT_TYPE* pProt);
    double acceptanceWindow_rel(PROT_TYPE* pProt, double dNewVal_frac);
    
    //  The static member function returns the duration of the statistic
    //  phase in units of the respiratory cycles.
    int32_t statisticScans(const PROT_TYPE* pProt);
    int32_t statisticScans(PROT_TYPE* pProt, int32_t lNewVal);

    //  If B0-Scout is used: The duration of an additional learning phase
    //  needed to determine the acceptance window and the PTThreshold.
    int32_t statisticDuration_us(const PROT_TYPE* pProt);
    int32_t statisticDuration_us(PROT_TYPE* pProt, int32_t i32NewVal);

        
    //  The first static member function returns a boolean value that
    //  indicates whether the trigger threshold is set manually or
    //  automatically. In the manual mode the trigger threshold is
    //  determined by the return value of triggerThreshold, in the
    //  automatic mode the optimal trigger threshold is determined
    //  retrospectively at the end of the learning phase of the
    //  gating algorithm.
    //  The second static member function is the related set function.
    bool manualTriggerThreshold(const PROT_TYPE* pProt);
    bool manualTriggerThreshold(PROT_TYPE* pProt, bool bNewVal);
    
    //  If the manual trigger mode is used the first static member
    //  function returns the trigger threshold in [0,1].
    //  The second static member function is the realated set function.
    double triggerThreshold(const PROT_TYPE* pProt);
    double triggerThreshold(PROT_TYPE* pProt, double dNewVal_us);
    
    //  The first static member function returns an enum that determines
    //  the physical quantity that determines the trigger point. The second
    //  static member function is the related set function.
    unsigned char B0ScoutTriggerPoint(const PROT_TYPE* pProt);
    unsigned char B0ScoutTriggerPoint(PROT_TYPE* pProt, unsigned char ucNewVal);
    unsigned char GREScoutTriggerPoint(const PROT_TYPE* pProt);
    unsigned char GREScoutTriggerPoint(PROT_TYPE* pProt, unsigned char ucNewVal);
    unsigned char NavScoutTriggerPoint(const PROT_TYPE* pProt);
    unsigned char NavScoutTriggerPoint(PROT_TYPE* pProt, unsigned char ucNewVal);
    unsigned char triggerPoint(const PROT_TYPE* pProt);
    unsigned char triggerPoint(PROT_TYPE* pProt, unsigned char ucNewVal);

    
    //  The first two static member function returns an enum that determines
    //  the way how double triggering is implemented. The third
    //  static member function is the related set function.
    unsigned char dblTriggerMode(const PROT_TYPE* pProt);
    unsigned char dblTriggerMode(PROT_TYPE* pProt, unsigned char ucNewVal);
        
    //  The first static member function returns a boolean value that
    //  indicates whether the trigger algorithm recalculates statistic data
    //  and adapts trigger conditions past the learning phase.
    //  The second static member function is the related set function.
    bool triggerAdaptive(const PROT_TYPE* pProt);
    bool triggerAdaptive(PROT_TYPE* pProt, bool bNewVal);

    bool triggerAdaptive_Level(const PROT_TYPE* pProt);
    bool triggerAdaptive_Level(PROT_TYPE* pProt, bool bNewVal);


    //  The first static member function returns a boolean value that
    //  indicates whether the detected physiological signal should be inverted.
    //  The second static member function is the related set function.
    bool invertSignal(const PROT_TYPE* pProt);
    bool invertSignal(PROT_TYPE* pProt, bool bNewVal);
    
    //  The first static member functions return the threshold used
    //  by the peak/trough detection algorithm in [mm].
    //  The second static member function returns the threshold used
    //  by the peak/trough detection algorithm in units of pixel.
    //  The third static member function is the related set function.
    double PTThreshold_mm (const PROT_TYPE* pProt);
    double PTThreshold_pix(const PROT_TYPE* pProt);
    double PTThreshold_mm(PROT_TYPE* pProt, double dNewVal_mm);

    //  If B0-Scout is used: the PTThreshold as fractional value 
    //  between the minimum and maximum value of the physiological signal
    //  found during the statistic phase.
    double PTThreshold_rel (const PROT_TYPE* pProt);
    double PTThreshold_rel(PROT_TYPE* pProt, double dNewVal_mm);

    // if new trigger algo "Level" is used, a unit-less value is used
    double PTThreshold_Level(const PROT_TYPE* pProt);
    double PTThreshold_Level(PROT_TYPE* pProt, double dNewVal);

    // if new trigger algo "Level" is used, in automatic positioning mode a minimum trigger Level will be read from protocol
    double MinimumTriggerLevel(const PROT_TYPE* pProt);
    double MinimumTriggerLevel(PROT_TYPE* pProt, double dNewVal);

    // if new trigger algo "Level" is used, in manual positioning mode a trigger Level will be read from protocol
    double ManualTriggerLevel(const PROT_TYPE* pProt);

    //  Channel independent scaling factor applide to complex image space images before sum-of-squares
    double FFTScale(const PROT_TYPE* pProt);
    double FFTScale(PROT_TYPE* pProt, double dNewVal);
    
    //  Trigger lock time in units of the median respiratory
    //  cycle.
    //  Trigger-lock-time*<Resp Cycle> is the minimum time
    //  interval between adjacent trigger events.
    double triggerLockTime(const PROT_TYPE* pProt);
    double triggerLockTime(PROT_TYPE* pProt, double dNewVal);
        
    //  The first static member function returns the number of respiratory
    //  intervals to be skipped between two adjacent trigger events.
    //  The second static member function is the related set function.
    int32_t triggerPulse(const PROT_TYPE* pProt);
    int32_t triggerPulse(PROT_TYPE* pProt, int32_t i32NewVal);
    
    //  The first static member functions return the minimum temporal
    //  distance between two profile images send to the host.
    //  The second static member function is the realated set function.
    double minSendIntrv_us(const PROT_TYPE* pProt);
    double minSendIntrv_us(PROT_TYPE* pProt, double dNewVal_us);

    //  In continuous mode the navigator is also played out during the
    //  imaging phase of the sequence.
    __IMP_EXP bool continuousMode(const PROT_TYPE* pProt);
    bool continuousMode(PROT_TYPE* pProt, bool bNewVal);

    //  Continuous fat sat mode means that the fat saturation pulse is also
    //  played out during the trigger phase.
    __IMP_EXP bool continuousSat(const PROT_TYPE* pProt);
    __IMP_EXP bool continuousSat(PROT_TYPE* pProt, bool bNewVal);

    //  Interleaved from Head to Feet id main orientation is tra
    __IMP_EXP bool seriesInterlDescending(const PROT_TYPE* pProt);
    bool seriesInterlDescending(PROT_TYPE* pProt, bool bNewVal);

    //  Exceute the very first navigator as noise adjust scan.
    bool noiseAdjScan(const PROT_TYPE* pProt);
    bool noiseAdjScan(PROT_TYPE* pProt, bool bNewVal);

    //  Perform noise decorrelation (using system functor, i.e. independent of noiseAdjScan).
    bool noiseDecorrelation(const PROT_TYPE* pProt);
    bool noiseDecorrelation(PROT_TYPE* pProt, bool bNewVal);


    //  Perform steady state preparation (i.e. modified flip-angle) at the 
    //  beginning of a series of navigator sequence with equal time delay.
    bool steadyStatePrep(const PROT_TYPE* pProt);
    bool steadyStatePrep(PROT_TYPE* pProt, bool bNewVal);

    //  T1 relaxation time of tissue component used for the calcualtion
    //  of steady state preparation module.
    double steadyStatePrepT1_ms(const PROT_TYPE* pProt);
    double steadyStatePrepT1_ms(PROT_TYPE* pProt, double dNewVal);

    //  TOM_MINIMIZE_TE or TOM_IN_PHASE or TOM_OPPOSED_PHASE
    uint8_t B0ScoutTOM(const PROT_TYPE* pProt);
    uint8_t B0ScoutTOM(PROT_TYPE* pProt, uint8_t u8NewVal);

    uint8_t GREScoutTOM(const PROT_TYPE* pProt);
    uint8_t GREScoutTOM(PROT_TYPE* pProt, uint8_t u8NewVal);

    //  Echo time of the phase scout in the case B0ScoutTOM(.) == TOM_OFF
    int32_t B0ScoutTE_us(const PROT_TYPE* pProt);
    int32_t B0ScoutTE_us(PROT_TYPE* pProt, int32_t i32NewVal);

    //  Echo time of the liver dome scout in the case B0ScoutTOM(.) == TOM_OFF
    int32_t GREScoutTE_us(const PROT_TYPE* pProt);
    int32_t GREScoutTE_us(PROT_TYPE* pProt, int32_t i32NewVal);

    //  Returns number of echoes for B0-Scout/GRE-Scout
    int16_t B0ScoutNEco(const PROT_TYPE* pProt);
    int16_t B0ScoutNEco(PROT_TYPE* pProt, int16_t i16NewVal);
    int16_t GREScoutNEco(const PROT_TYPE* pProt);
    int16_t GREScoutNEco(PROT_TYPE* pProt, int16_t i16NewVal);

    //  Returns number of excitations pulses for B0-Scout
    //  Return value greater than 1 means that a binominal excitation pulse
    //  for spectral selective water esxcitation is used.
    int16_t  B0ScoutNRF(const PROT_TYPE* pProt);
    int16_t  B0ScoutNRF(PROT_TYPE* pProt, int16_t  i16NewVal);
    int16_t  GREScoutNRF(const PROT_TYPE* pProt);
    int16_t  GREScoutNRF(PROT_TYPE* pProt, int16_t  i16NewVal);

    //  Returns boolean flag, indicating whether velocity compensation is used
    bool B0ScoutVelComp(const PROT_TYPE* pProt);
    bool B0ScoutVelComp(PROT_TYPE* pProt, bool bNewVal);
    bool GREScoutVelComp(const PROT_TYPE* pProt);
    bool GREScoutVelComp(PROT_TYPE* pProt, bool bNewVal);

    //  If number of echoes is greater than one:
    //  Returns boolean flag, indicating whether bipolar (EPI like) readout gradient
    //  scheme should be used.
    bool B0ScoutBipolar(const PROT_TYPE* pProt);
    bool B0ScoutBipolar(PROT_TYPE* pProt, bool bNewVal);
    bool GREScoutBipolar(const PROT_TYPE* pProt);
    bool GREScoutBipolar(PROT_TYPE* pProt, bool bNewVal);

    //  Returns boolean flag, indicating whether the readout time should be minimized.
    bool B0ScoutShortestDwellTime(const PROT_TYPE* pProt);
    bool B0ScoutShortestDwellTime(PROT_TYPE* pProt, bool bNewVal);
    bool GREScoutShortestDwellTime(const PROT_TYPE* pProt);
    bool GREScoutShortestDwellTime(PROT_TYPE* pProt, bool bNewVal);



    //  If true, uses separate coil select for navigator.
    bool B0ScoutAutoCoilSelect(const PROT_TYPE* pProt);
    bool B0ScoutAutoCoilSelect(PROT_TYPE* pProt, bool bNewVal);
    bool GREScoutAutoCoilSelect(const PROT_TYPE* pProt);
    bool GREScoutAutoCoilSelect(PROT_TYPE* pProt, bool bNewVal);

    //  Rf spoiling for navigator scan
    bool B0ScoutRFSpoiling(const PROT_TYPE* pProt);
    bool B0ScoutRFSpoiling(PROT_TYPE* pProt, bool bNewVal);
    bool GREScoutRFSpoiling(const PROT_TYPE* pProt);
    bool GREScoutRFSpoiling(PROT_TYPE* pProt, bool bNewVal);

    //  SELECT_ACQ_WINDOW_AUTOMATIC  ...  Acquisition window is defined as fractional value of the respiratory cycle 
    //  SELECT_ACQ_WINDOW_MANUAL     ...  Acquisition window is defined as absolute value 
    int32_t selectAcqWindow(const PROT_TYPE* pProt);
    int32_t selectAcqWindow(PROT_TYPE* pProt, int32_t i32NewVal);

    //  Acquisition window in us
    int32_t acqWindow_us(const PROT_TYPE* pProt);
    int32_t acqWindow_us(PROT_TYPE* pProt, int32_t i32NewVal);

    //  Acquisition window as fractional value of the respiratory cycle.
    __IMP_EXP double acqWindow_rel(const PROT_TYPE* pProt);
    __IMP_EXP double acqWindow_rel(PROT_TYPE* pProt, double dNewVal);

    //  POSITION_NAVIGATOR_AUTOMATIC  ...  Navigator is adapted relative to the slice stack.
    //  POSITION_NAVIGATOR_MANUAL     ...  Navigator is set manually
    int32_t positionNavigator(const PROT_TYPE* pProt);
    int32_t positionNavigator(PROT_TYPE* pProt, int32_t i32NewVal);

    //  Returns enum of type SEQ::PhysioSignal. Default SIGNAL_NONE.
    int32_t refSignal(const PROT_TYPE* pProt);
    int32_t refSignal(PROT_TYPE* pProt, int32_t i32NewVal);

    //  The first static member function returns an enumeration, which 
    //  determines the k-space traverse speed.
    //  The returned value is either
    //  TRAVERSE_CONST_ANGULAR_RATE ... Gradient limits of outmost ring determines the duration of all rings
    //  TRAVERSE_CONST_SLEW_RATE    ... Duration of each ring is limited by the gardient ampkitude/slew rate.
    //  The second static member function is the associated set-function.
    uint8_t NavScoutTraverse(const PROT_TYPE* pProt);
    uint8_t NavScoutTraverse(PROT_TYPE* pProt, uint8_t u8NewVal);

    //  The first static member function returns an enumeration, which 
    //  determines the k-space rotation direction, in the case of a concentric
    //  ring trajectory.
    //  ROTATION_COUNTERCLOCKWISE  ... clockwise
    //  ROTATION_CLOCKWISE         ... counterclockwise
    //  ROTATION_INTERLEAVED       ... counterclockwise-clockwise-counterclockwise-...
    uint8_t NavScoutRotation(const PROT_TYPE* pProt);
    uint8_t NavScoutRotation(PROT_TYPE* pProt, uint8_t u8NewVal);

    //  Excitation k-space trajectory
    //  TRAJECTORY_CR  ...  concentric ring trajectory
    //  TRAJECTORY_EPI ...  EPI like trajectory.
    //  The second static member function is the associated set-function.
    uint8_t NavScoutTrajectory(const PROT_TYPE* pProt);
    uint8_t NavScoutTrajectory(PROT_TYPE* pProt, uint8_t u8NewVal);

    //  The returned value is either
    //  PROFILE_SINC   ... constant k-space weighting
    //  PROFILE_GAUSS  ... gaussian k-space weighting
    //  The second static member function is the associated set-function.
    uint8_t NavScoutProfile(const PROT_TYPE* pProt);
    uint8_t NavScoutProfile(PROT_TYPE* pProt, uint8_t u8NewVal);

    //  The dimensionless return value times the resolution dr of the 2D Rf-pulse is equal to the diameter of the excitation profile:
    //  dr = Diameter/Resolution, where Resolution is the return value of the function.
    //  The resolution dr determines the k-space radius: k_max = 1/(2*dr)
    double NavScoutResolution(const PROT_TYPE* pProt);
    double NavScoutResolution(PROT_TYPE* pProt, double dNewVal);

    //  Gating algorithm
    uint8_t GatingAlgorithm(const PROT_TYPE* pProt);
    uint8_t GatingAlgorithm(PROT_TYPE* pProt, uint8_t u8NewVal);
    
    // Boolean value, indicating that the same diaphragm position maps to two different
    //  states, if the first is acquired during inspiration and the second during expiration.
    bool GatingExploitRespPhase(const PROT_TYPE* pProt);
    bool GatingExploitRespPhase(PROT_TYPE* pProt, bool bNewVal);


};

#endif // __PACE_PROT_H
