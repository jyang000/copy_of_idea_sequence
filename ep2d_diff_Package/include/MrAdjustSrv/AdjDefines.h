// -----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 2009  All Rights Reserved.  Confidential
// -----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//    File: \n4_servers1\pkg\MrServers\MrAdjustSrv\AdjDefines.h
// Version: \main\110
//  Author: cc_adjust
//    Date: 2013-11-14 15:08:44 +01:00
//
// -----------------------------------------------------------------------------

#ifndef ADJDEFINES_H_INCLUDED
#define ADJDEFINES_H_INCLUDED

// -----------------------------------------------------------------------------
// AdjData
// -----------------------------------------------------------------------------

#define ADJ_FALSE 0
#define ADJ_TRUE  1

// -----------------------------------------------------------------------------

#define ADJ_RXGAIN_LOW  0
#define ADJ_RXGAIN_HIGH 1

// -----------------------------------------------------------------------------
// AdjType
// -----------------------------------------------------------------------------

#define ADJ_KIND_ALL       -1

// -----------------------------------------------------------------------------

#define ADJ_MODE_DEFAULT       0

// AdjTuneLC

// AdjFre
#define ADJ_MODE_STEAM10       0
#define ADJ_MODE_FID400        1
#define ADJ_MODE_FID10         2
#define ADJ_MODE_FID           3
#define ADJ_MODE_XFRECALC     10

// AdjTra
#define ADJ_MODE_DOUBLEECHO    0
#define ADJ_MODE_DOUBLEECHO1D  1
#define ADJ_MODE_XTRACALC     10

// AdjDico
#define ADJ_MODE_NORMAL        0
#define ADJ_MODE_GAINVAR       1
#define ADJ_MODE_CONTINUOUS    2

// AdjFieldMap
#define ADJ_MODE_STANDARD      0
#define ADJ_MODE_ADVANCED      1
#define ADJ_MODE_MULTIPLE      2

// AdjShim
#define ADJ_MODE_TUNEUP        0
#define ADJ_MODE_CALC3DSHIM   10
#define ADJ_MODE_INTERACTIVE  20

// AdjWatSup

// AdjRFMap

// AdjB1Shim

// AdjCoilCorr

// AdjCoilSens

// AdjCoilPos
#define ADJ_MODE_STATIC        0
#define ADJ_MODE_MOVING       10

// AdjMDS

// AdjMDSScout

#define ADJ_MODE_MAX          99

// -----------------------------------------------------------------------------

#define ADJ_TYPE_MAXSTRLENGTH 128

// -----------------------------------------------------------------------------
// AdjStatus
// -----------------------------------------------------------------------------

#define ADJ_STATUS_INIT           0
#define ADJ_STATUS_MEMORY        11
#define ADJ_STATUS_EXECUTION     12
#define ADJ_STATUS_ABORTED      101
#define ADJ_STATUS_USERSTOP     102
#define ADJ_STATUS_HALTED       111
#define ADJ_STATUS_CALCULATED   112
#define ADJ_STATUS_PROTOCOL     201
#define ADJ_STATUS_SECTION      202
#define ADJ_STATUS_MANUAL       203
#define ADJ_STATUS_SYSTEM       204
#define ADJ_STATUS_MDS          205
#define ADJ_STATUS_CONVERGED    211

#define ADJ_STATUS_MAX          999

#define ADJ_STATUS_SIMULATED  65536

// -----------------------------------------------------------------------------

#define ADJ_STATUS_MAXSTRLENGTH 128

// -----------------------------------------------------------------------------
// AdjContext
// -----------------------------------------------------------------------------

#define ADJ_CONTEXT_MAXSTRLENGTH      256

// -----------------------------------------------------------------------------
// AdjResult
// -----------------------------------------------------------------------------

#define ADJ_RESULT_MAXSTRLENGTH       256

// -----------------------------------------------------------------------------
// AdjQuery
// -----------------------------------------------------------------------------

#define ADJ_QUERY_MAXSTRLENGTH     256

// -----------------------------------------------------------------------------
// AdjVector
// -----------------------------------------------------------------------------

#define ADJ_VECTOR_INTERMEDIATE 0
#define ADJ_VECTOR_MANDATORY    1

// -----------------------------------------------------------------------------

#define ADJ_VECTOR_MAXLENGTH        64

// -----------------------------------------------------------------------------
// AdjError
// -----------------------------------------------------------------------------

#define ADJ_ERROR_ONSTACK        1
#define ADJ_ERROR_TODISPLAY      2
#define ADJ_ERROR_STACKTODISPLAY 4

// -----------------------------------------------------------------------------

#define ADJ_ERROR_MAXSTRLENGTH 128

// -----------------------------------------------------------------------------
// AdjInfo
// -----------------------------------------------------------------------------

#define ADJ_INFO_TUNELC    1
#define ADJ_INFO_FRE       2
#define ADJ_INFO_XFRE      3
#define ADJ_INFO_TRA       4
#define ADJ_INFO_DICO      5
#define ADJ_INFO_FIELDMAP  6
#define ADJ_INFO_3DSHIM    7
#define ADJ_INFO_INTSHIM   8
#define ADJ_INFO_WATSUP    9
#define ADJ_INFO_RFMAP    10
#define ADJ_INFO_B1SHIM   11
#define ADJ_INFO_COILCORR 12
#define ADJ_INFO_COILSENS 13
#define ADJ_INFO_COILPOS  14
#define ADJ_INFO_MDS      15
#define ADJ_INFO_MDSSCOUT 16

// -----------------------------------------------------------------------------

#define ADJ_INFO_START  1
#define ADJ_INFO_IMAGES 2
#define ADJ_INFO_STOP   3

// -----------------------------------------------------------------------------

#define ADJ_INFO_MAXSTRLENGTH 128

// -----------------------------------------------------------------------------
// AdjSpec
// -----------------------------------------------------------------------------

#define ADJ_SPEC_MAXSTRLENGTH     256

// -----------------------------------------------------------------------------
// AdjPlot
// -----------------------------------------------------------------------------

#define ADJ_PLOT_COMBINED_ADC -1

// -----------------------------------------------------------------------------
// AdjServerIFImpl
// -----------------------------------------------------------------------------

#define ADJ_SERVERIFIMPL_CONNECTION_TIMEOUT 10000 // [ms]
#define ADJ_SERVERIFIMPL_LOOKUP_TIMEOUT     10000 // [ms]
#define ADJ_SERVERIFIMPL_ADJUST_TIMEOUT        -1 // [ms] no timeout
#define ADJ_SERVERIFIMPL_MESSAGE_TIMEOUT    10000 // [ms]

// -----------------------------------------------------------------------------
// AdjControlCache
// -----------------------------------------------------------------------------

#define ADJ_CONTROLCACHE_MAXLENGTH 128

// -----------------------------------------------------------------------------
// AdjCommand
// -----------------------------------------------------------------------------

#define ADJ_COMMAND_CONNECT                0
#define ADJ_COMMAND_LOOKUP                 1
#define ADJ_COMMAND_PREPARE_VECTOR         2
#define ADJ_COMMAND_PREPARE_CONTEXT_RESULT 3
#define ADJ_COMMAND_ADJUST                 4
#define ADJ_COMMAND_MEASURE                5
#define ADJ_COMMAND_INVALIDATE_KIND        6
#define ADJ_COMMAND_INVALIDATE_UID         7
#define ADJ_COMMAND_VALIDATE               8

// -----------------------------------------------------------------------------
// AdjSchedulerList
// -----------------------------------------------------------------------------

#define ADJ_SCHEDULERLIST_MAXLENGTH     32
#define ADJ_SCHEDULERLIST_MAXDEPENDENCY  8

// -----------------------------------------------------------------------------
// AdjCoilSelect
// -----------------------------------------------------------------------------

#define ADJ_COILSELECT_MAXLENGTH    256
#define ADJ_COILSELECT_MAXSTRLENGTH 128

// -----------------------------------------------------------------------------
// AdjMemoryCore
// -----------------------------------------------------------------------------

#define ADJ_MEMORYCORE_MAXLENGTH 256

// -----------------------------------------------------------------------------
// AdjCom
// -----------------------------------------------------------------------------

#define ADJ_COM_MAXLENGTH 128

// -----------------------------------------------------------------------------
// AdjTuneLC
// -----------------------------------------------------------------------------

#define ADJ_TUNELC_PLUGS                       2 // []
#define ADJ_TUNELC_CHANNELS_PER_PLUG           2 // []
#define ADJ_TUNELC_MAX_CHANNELS                4 // []
#define ADJ_TUNELC_FREQUENCY_STEPS           512 // []
#define ADJ_TUNELC_ALGORITHM_5_STEP            1 // []
#define ADJ_TUNELC_ALGORITHM_N_STEP            2 // []
#define ADJ_TUNELC_ALGORITHM_PHASECORRECTION   3 // []

// -----------------------------------------------------------------------------
// AdjFre
// -----------------------------------------------------------------------------

#define ADJ_XFRE_OVERSAMPLING_FACTOR 2 // []

// -----------------------------------------------------------------------------
// AdjTra
// -----------------------------------------------------------------------------

#define ADJ_TRA_DURATION_TIME_REF_PULSE 1000 // [us]

#define ADJ_TRA_USE_0D              "UseAdjTra0D"
#define ADJ_TRA_USE_1D              "UseAdjTra1D"
#define ADJ_TRA_USE_2D              "UseAdjTra2D"
#define ADJ_TRA_USE_2D_OFFCENTER    "UseAdjTraOC"
#define ADJ_TRA_IGNORE_2D_OFFCENTER "IgnAdjTraOC"

#define ADJ_TRA_MODE_0D 0
#define ADJ_TRA_MODE_1D 1
#define ADJ_TRA_MODE_2D 2

// -----------------------------------------------------------------------------
// AdjDico
// -----------------------------------------------------------------------------

#define ADJ_DICO_DURATION_TIME_RF_PULSE         1500 // [us]
#define ADJ_DICO_DURATION_TIME_RF_PULSE_SHORT    200 // [us]
#define ADJ_DICO_DURATION_TIME_READOUT          2560 // [us]
#define ADJ_DICO_DURATION_TIME_EVALUATION       1000 // [us]
#define ADJ_DICO_DURATION_TIME_EVALUATION_SHORT  120 // [us]


#define ADJ_DICO_MAX_NUM_OF_PICKUP_PROBES  2  // []
#define ADJ_DICO_MAX_NUM_OF_CHARA_SAMPLES 10  // []
#define ADJ_DICO_DICO_PULSE                0  // []
#define ADJ_DICO_PICKUP_PULSE              1  // []
#define ADJ_DICO_SSB1_FORWARD              0  // []
#define ADJ_DICO_SSB1_REFLECTED            1  // []
#define ADJ_DICO_TAS_DUMMY_REFLECTED       2  // []
#define ADJ_DICO_PICKUP_PROBE1             8  // []
#define ADJ_DICO_PICKUP_PROBE2             9  // []
#define ADJ_DICO_RCCS_PICKUP_PROBE1        0  // []
#define ADJ_DICO_RCCS_PICKUP_PROBE2        1  // []
#define ADJ_DICO_SSB2_FORWARD              16 // []
#define ADJ_DICO_SSB2_REFLECTED            17 // []


#define ADJ_DICO_MAX_VALUES            32 // []
#define ADJ_DICO_ARBITRARY_AMPLITUDE    0 // []
#define ADJ_DICO_ZERO_AMPLITUDE         1 // []
#define ADJ_DICO_TRAREF_AMPLITUDE       2 // []
#define ADJ_DICO_CALIBRATION_AMPLITUDE  3 // []
#define ADJ_DICO_MAX_AMPLITUDE          4 // []
#define ADJ_DICO_BC_SIGNAL              0 // []
#define ADJ_DICO_BC0_SIGNAL             1 // []
#define ADJ_DICO_BC1_SIGNAL             2 // []
#define ADJ_DICO_LC_SIGNAL              3 // []

// -----------------------------------------------------------------------------
// AdjFieldMap
// -----------------------------------------------------------------------------

#define ADJ_FIELDMAP_COMBINED_BY_RX_GROUPS 0 // fieldmap data are combined by rx group ids
#define ADJ_FIELDMAP_COMBINED              1 // fieldmap data are combined to 1 channel
#define ADJ_FIELDMAP_UNCOMBINED            2 // fieldmap data are uncombined

// -----------------------------------------------------------------------------
// AdjShim
// -----------------------------------------------------------------------------

#define ADJ_SHIM_MAX_CHANNELS               15 // []
#define ADJ_LOCAL_SHIM_MAX_CHANNELS          8 // []
#define ADJ_HIGHER_SHIM_MAX_CHANNELS         8 // []    
#define ADJ_INTSHIM_PHYSIO_RESOLUTION      256 // []
#define ADJ_INTSHIM_INFINITE_ITERATIONS 999999 // []

// -----------------------------------------------------------------------------
// AdjWatSup
// -----------------------------------------------------------------------------
#define ADJ_WATSUP_NUMBER_OF_PRESCANS         6 
#define ADJ_WATSUP_NUMBER_OF_CORRECTION_SCANS 11
#define ADJ_WATSUP_CORRECTION_FACTOR          1.0
#define ADJ_WATSUP_CORRECTION_DELTA           0.1

// -----------------------------------------------------------------------------
// AdjRFMap
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// AdjB1Shim
// -----------------------------------------------------------------------------

#define ADJ_B1SHIM_MAX_CHANNELS 32 // []

// -----------------------------------------------------------------------------
// AdjCoil
// -----------------------------------------------------------------------------

#define ADJ_RXCOILS_NOISE_SCANS_NAME  "AdjCoilCorr_NoiseScans_"
#define ADJ_RXCOILS_DECORR_MATRIX     "AdjCoilCorr_DecorrMatrix"
#define ADJ_RXCOILS_COILCORR_RAW_DATA "AdjCoilCorr_RawData"

#define ADJ_RXCOILS_COILSENS_RES_32  "PSNR_32"  // resolution 32x32x64
#define ADJ_RXCOILS_COILSENS_RES_32A "PSNR_32A" // resolution 32x32x64 with iPAT
#define ADJ_RXCOILS_COILSENS_RES_64A "PSNR_64A" // resolution 64x64x64 with iPAT
#define ADJ_RXCOILS_COILSENS_RES_96A "PSNR_96A" // resolution 96x96x96 with iPAT

// -----------------------------------------------------------------------------
// AdjCoilPos
// -----------------------------------------------------------------------------

#define ADJ_COILPOS_CDT_POSTERIOR_PARAMCARD_Y_POS 50 // []

// -----------------------------------------------------------------------------
// AdjDynAdjust
// -----------------------------------------------------------------------------

#define ADJ_DYNADJUST_MAX_VOLUMES 512 // []

#define ADJ_DYNADJUST_FLAG_ISLIMITEDTOVOI  0x02 // Indicates whether the cuboid specifies a saturation region               (default: no )
#define ADJ_DYNADJUST_FLAG_ISSATURATION    0x04 // Indicates whether the cuboid should get limited to the VOI               (default: yes)
#define ADJ_DYNADJUST_FLAG_ISSATURATION_EX 0x08 // Indicates whether saturation regions should get excluded from the cuboid (default: yes)

// -----------------------------------------------------------------------------
// AdjMDS
// -----------------------------------------------------------------------------

#define ADJ_MDS_MAXLENGTH     100 // []
#define ADJ_MDS_ADJUSTMENTS     5 // []
#define ADJ_MDS_FRE_INDEX       0 // []
#define ADJ_MDS_TRA_INDEX       1 // []
#define ADJ_MDS_DICO_INDEX      2 // []
#define ADJ_MDS_COILPOS_INDEX   3 // []

// -----------------------------------------------------------------------------
// AdjMDSScout
// -----------------------------------------------------------------------------

#define ADJ_MDSSCOUT_IMG_NAME "AdjMDSScout_Scout_Image_"

// -----------------------------------------------------------------------------
// AdjEvaluation
// -----------------------------------------------------------------------------

#define ADJ_EVALUATION_COMBINED_ADC -1
#define ADJ_EVALUATION_NO_ADC       -2

// -----------------------------------------------------------------------------
// Signal limits
// -----------------------------------------------------------------------------
#define ADJ_SIGNAL_TOO_HIGH_LIMIT   32767.0

// -----------------------------------------------------------------------------
// AdjSpecialMode
// -----------------------------------------------------------------------------

#define ADJ_SPECIALMODE_TRA_DOUBLE_ECHO                  1 // Use DoubleEcho 0D
#define ADJ_SPECIALMODE_TRA_DOUBLE_ECHO_1D               2 // Use DoubleEcho 1D
#define ADJ_SPECIALMODE_TRA_B1MAP_2D                     4 // Use B1Map 2D
#define ADJ_SPECIALMODE_FRE_ANALYSE_COMBINED             8 // Analyse all RX chennels in combined mode
#define ADJ_SPECIALMODE_FIELDMAP_USE_TRA_REF_AMPL       16 // Use AdjTra result for fieldmap measurement
#define ADJ_SPECIALMODE_SHIM_SKIP_MGPA_CHECK            32 // mGPA only: skip gradoffset and gradsensitivity check
#define ADJ_SPECIALMODE_SHIM_WRITE_B0MAP_ALWAYS         64 // Always activate B0MapAcq
#define ADJ_SPECIALMODE_COILSENS_10_AVERAGES           128 // Used by SW test for loudness measurement
#define ADJ_SPECIALMODE_MDS_ENABLE_TRA_FEEDBACK        256 // Apply AdjTra results during MDS adjustment
#define ADJ_SPECIALMODE_MDS_DISABLE_DICO_FEEDBACK      512 // Don't apply AdjDico results during MDS adjustment
#define ADJ_SPECIALMODE_MDS_ENABLE_IMAGING            1024 // Enable MDS Sequences (otherwise MDS sequences are blocked in AdjSafety)
#define ADJ_SPECIALMODE_MDSSCOUT_BODY_REGIONS         2048 // Show body regions menu (special parameter card)#define ADJ_SPECIALMODE_MDSSCOUT_OFFCENTER            8192 // Measure offcenter, if a part of the scan range can't be measured isocenter
#define ADJ_SPECIALMODE_MDSSCOUT_OFFCENTER            4096 // Measure offcenter, if a part of the scan range can't be measured isocenter
#define ADJ_SPECIALMODE_ADJ_EXEC_USE_PHYSSERVER       8192 // Execute calculations on physserver
#define ADJ_SPECIALMODE_EXAMDB_FROM_CUSTOMER         16384 // Use adjustment protocols from customer tree
#define ADJ_SPECIALMODE_ADJSEQ_FROM_CUSTOMER         32768 // Load adjustment sequences form customer directory
#define ADJ_SPECIALMODE_SERVERIF_DISABLE_CACHE       65536 // Disable clientside cache
#define ADJ_SPECIALMODE_SAFETY_DEBUG_MODE           131072 // Enhanced debug output
#define ADJ_SPECIALMODE_COILSELECT_DISABLE_CACHE    262144 // Disable coilselect cache (AdjServer)
#define ADJ_SPECIALMODE_NO_CALLSTACK_IN_ADJTEMP     524288 // Deadlock detection
#define ADJ_SPECIALMODE_PHYS_DEBUG_MODE            1048576 // Enhanced debug output
#define ADJ_SPECIALMODE_COILPOS_NORMALIZE_BC       2097152 // Measure body coil for normalization
                                                  
#define ADJ_SPECIALMODE_LAST_ENTRY                 2097152 // Last entry

// -----------------------------------------------------------------------------

#endif
