//Conversion factors
//------------------
// Seconds to microseconds
#define SEC2USEC 1000000.0

//Default Values and Constants
#define NMAXMEGS                             2 // * Maximum number of MEG instances *
#define MEG_DEFAULT_AXIS        WIP_Axis_Slice // Default MEG axis
#define MEG_DEFAULT_HZ                    60.0 // Default Driver frequency
#define MEG_DEFAULT_PERIOD_US            16667 // Default Driver period
#define EXT_TRIGGER_DURATION_US             20 // Necessary for the Resoundant driver
#ifdef MRE_EPI
    #define MRE_TR                          50 // default TE for the MRE sequence
    #define MEG_DEFAULT_EPI_FACTOR           1 // Sequence is a GRE for the beginning
    #define MEG_DEFAULT_DUMMYSCAN_DURATION 501 // ms dummy scan duration 

    #define MEASDATA_MCIR_DIRNAME "/opt/med/SimMeasData/MRE.ini"
    #define MEASDATA_HOST_INIFILENAME "C:\\MedCom\\MCIR\\med\\SimMeasData\\MRE.ini"

#endif

//Switches
#define SOLITARY_MEG 1       // Allow no overlap with MEG and other prephasing gradients   

//Hard limits for MEG frequency box
#define MEG_FREQUENCY_HARD_LOW_HZ 30.0
#define MEG_FREQUENCY_HARD_HIGH_HZ 300.0

//Macros for WIP specific activity
//---------------------------------
// Retrieve MEG frequency from Prot
#define GET_MEG_FREQUENCY rMrProt.wipMemBlock().getadFree()[WIP_MEGFrequency]
// Retrieve MEG axis from Prot
#define GET_MEG_AXIS rMrProt.wipMemBlock().getalFree()[WIP_MEGGradientAxis_Box]
// Retrieve MEG period
#define GET_MEG_PERIOD SEC2USEC/GET_MEG_FREQUENCY

