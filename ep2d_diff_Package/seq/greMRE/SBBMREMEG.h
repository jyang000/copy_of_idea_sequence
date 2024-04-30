//----------------------------------------------------------------------------------
// <copyright file="SBBMREMEG.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 2006-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015-2021. All Rights Reserved. Confidential.
// </copyright>
// <description>
//   Motion Encoding Gradient (MEG)
//   The MEG gradients consist of three lobes, MEG1, MEG2, and MEG3
//   The middle gradient lobe has the opposite polarities of the other two
//   gradient lobes. The middle gradient has the maximum amplitude, and the
//   other two gradients amplitudes are calculated to make the zero moment
//   zero. These two gradients have the same amplitude.
// </description>
//----------------------------------------------------------------------------------

#pragma once

//  -------------------------------------------------------------------------- 
//  Includes                                                                   
//  -------------------------------------------------------------------------- 

#include "MrImagingFW/libSBBFW/SeqBuildBlock.h"
#include "MrImaging/seq/greMRE/a_gre_mre_def.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"
#include "MrImaging/seq/greMRE/MREProperties.h"


class SBBMREMEG : public SeqBuildBlock
{

    public:

        //  ---------------------------------------------- 
        //  Constructor                                 
        //  ---------------------------------------------- 

        SBBMREMEG( SBBList *pSBBList);

        SBBMREMEG();

        //  ----------------------------------------------- 
        //  Destructor                                       
        //  ----------------------------------------------- 

        virtual ~SBBMREMEG() = default;

        SBBMREMEG(const SBBMREMEG& right) = delete;
        SBBMREMEG& operator=(const SBBMREMEG& right) = delete;
        SBBMREMEG(SBBMREMEG&& right)                 = delete;
        SBBMREMEG& operator=(SBBMREMEG&& right) = delete;

        //  ----------------------------------------------- 
        //  Specifier of the MEG axis direction          
        //  ----------------------------------------------- 

        enum eGradientAxis  { None, Slice, Read, Phase, None_With_Pulse };             
        

       //  ----------------------------------------------- 
        //  Specifier of the MEG mode          
        //  ----------------------------------------------- 

        enum eMEGMode       { Standard,                 ///< 3 MEG Gradients, full duration for all gradients (see definition of MEG_DEFAULT_PERIOD_US). 
                              Two_Gradients,            ///< only 2 MEG gradients, duration of MEG reduced to 2/3 of MEG_DEFAULT_PERIOD_US.
                              ReduceMEGDuration,        ///< all 3 MEG gradient durations are reduced. Reduction factor from MRE.ini file.
                              ReduceOuterLobDuration,   ///< reducing the outer lobs in duration, Reduction factor from MRE.ini file. center MEG stays unchanged
                              ReduceOuterLobs};         ///<  not used at the moment.z

        //  ------------------------------------------------ 
        //  Specifier of the polarity of the MEG              
        //  ------------------------------------------------ 

        enum eMEGPolarity   { Negative = -1, Positive = 1 };              

        //  --------------------------------------------------------------------
        //                                                                      
        //  Name        :  SBBMREMEG::prep                        
        //                                                                      
        //  Description :  Preparation of the MEG gradient          
        //                                                                      
        //  Parameters  :  pMrProt, pSeqLim, pSeqExpo                                                
        //                                                                      
        //  Return      :  bool                                                 
        //                                                                      
        //  ------------------------------------------------------------------- 

        virtual bool prep (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        //  --------------------------------------------------------------------
        //                                                                      
        //  Name        :  SBBMREMEG::check                        
        //                                                                      
        //  Description :  Check of the MEG gradient          
        //                                                                      
        //  Parameters  :  pSeqLim,                                                
        //                                                                      
        //  Return      :  bool                                                 
        //                                                                      
        //  ------------------------------------------------------------------- 

        virtual bool check (SeqLim &rSeqLim);

        //  --------------------------------------------------------------------
        //                                                                      
        //  Name        :  SBBMREMEG::run                        
        //                                                                      
        //  Description :  Run method for the MEG gradient       
        //                                                                      
        //  Parameters  :  pMrProt, pSeqLim, pSeqExpo, pSLC                                                
        //                                                                      
        //  Return      :  bool                                                 
        //                                                                      
        //  ------------------------------------------------------------------- 


        virtual bool run (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS *pSLC);

        //  --------------------------------------------------------------------
        //                                                                      
        //  Name        :  SBBMREMEG::updateMEGGradPerf                        
        //                                                                      
        //  Description :  Set up specific gradient performance specs for the MEG       
        //                                                                      
        //  Parameters  :  eGradientMode                                                
        //                                                                      
        //  Return      :  bool                                                 
        //                                                                      
        //  ------------------------------------------------------------------- 


        virtual bool updateMEGGradPerf (SEQ::Gradients eGradientMode);


        //  --------------------------------------------------------------------
        //                                                                      
        //  Name        :  SBBMREMEG::setMEGGradientAxis                        
        //                                                                      
        //  Description :  Specifies the axis of the MEG to be applied          
        //                                                                      
        //  Parameters  :  eAxis                                                
        //                                                                      
        //  Return      :  void                                                 
        //                                                                      
        //  ------------------------------------------------------------------- 
  
        void setMEGGradientAxis (eGradientAxis eAxis);

        //  ------------------------------------------------------------------- 
        //                                                                     
        //  Name        :  SBBMREMEG::getMEGGradientAxis                       
        //                                                                     
        //  Description :  Returns the gradient axis of the MEG,               
        //                 Slice, Read, Phase, None                            
        //                                                                     
        //  Return      :  eGradientAxis                                       
        //                                                                     
        //  ------------------------------------------------------------------ 

        eGradientAxis getMEGGradientAxis ();


        //  --------------------------------------------------------------------
        //                                                                      
        //  Name        :  SBBMREMEG::setMEGMode                        
        //                                                                      
        //  Description :  Sets the MEG timing mode          
        //                                                                      
        //  Parameters  :  eMEGMode                                                
        //                                                                      
        //  Return      :  void                                                 
        //                                                                      
        //  ------------------------------------------------------------------- 

        void setMEGMode (eMEGMode eMegMode);


        //  --------------------------------------------------------------------
        //                                                                      
        //  Name        :  SBBMREMEG::getMEGMode                        
        //                                                                      
        //  Description :  gets the MEG timing mode          
        //                                                                      
        //  Return      :  eMEGMode                                                 
        //                                                                      
        //  ------------------------------------------------------------------- 

        eMEGMode getMEGMode();

        // -------------------------------------------------------------------
        //                                                                    
        //  Name        :  SBBMREMEG::setMEGPolarity                          
        //                                                                    
        //  Description :  Specifies the polarity of the MEG Gradient         
        //                                                                    
        //  Parameters  :  ePolarity                                          
        //                                                                    
        //  Return      :  void                                               
        //                                                                    
        // -------------------------------------------------------------------

        void setMEGPolarity (eMEGPolarity ePolarity);


        //  ------------------------------------------------------------------- 
        //                                                                     
        //  Name        :  SBBMREMEG::getMEGPolarity                           
        //                                                                     
        //  Description :  Returns the polarity of the MEG gradient,                      
        //                 Slice, Read, Phase, None                            
        //                                                                     
        //  Return      :  eMEGPolarity                                        
        //                                                                     
        //  ------------------------------------------------------------------ 

        eMEGPolarity getMEGPolarity ();

        //  --------------------------------------------------------------------
        // 
        //  Name        :  SBBMREMEG::setMEGFrequency
        // 
        //  Description :  Specifies the frequency of the MEG to be applied
        // 
        //  Parameters  :  dFrequency
        // 
        //  Return      :  void
        // 
        //  -------------------------------------------------------------------

        bool setMEGFrequency (double dFrequency);

        //  -------------------------------------------------------------------
        // 
        //  Name        :  SBBMREMEG::getMEGFrequency
        // 
        //  Description :  Return the ferquency of the MEG Gradient.
        // 
        //  Parameters  :  none
        // 
        //  Return      :  double
        // --------------------------------------------------------------------
     
        double getMEGFrequency ();

        //  --------------------------------------------------------------------
        // 
        //  Name        :  SBBMREMEG::setMEGAmplitude
        // 
        //  Description :  Specifies the max amplitude of the MEG to be applied in mT/m
        // 
        //  Parameters  :  dAmplitude
        // 
        //  Return      :  void
        // 
        //  -------------------------------------------------------------------

        bool setMEGAmplitude(double dAmplitude);

        //  --------------------------------------------------------------------
        // 
        //  Name        :  SBBMREMEG::setMinRiseTime
        // 
        //  Description :  Specifies the max Rise Time of the MEG to be applied in 
        // 
        //  Parameters  :  dRiseTime
        // 
        //  Return      :  void
        // 
        //  -------------------------------------------------------------------

        bool setMinRiseTime(double dRiseTime);

        // --------------------------------------------------------------------
        // 
        //  Name        :  SBBMREMEG::setlStartTime
        // 
        //  Description :  Specifies the start time of the MEG Gradient.
        // 
        //  Parameters  :  none
        // 
        //  Return      :  long
        // 
        // --------------------------------------------------------------------

        void setlStartTime (long lTime);

        // --------------------------------------------------------------------
        // 
        //  Name        :  SBBMREMEG::getlStartTime
        // 
        //  Description :  Return the start time of the MEG Gradient.
        // 
        //  Parameters  :  none
        // 
        //  Return      :  long
        // 
        // -------------------------------------------------------------------

        long getlStartTime ();

        // --------------------------------------------------------------------
        // 
        //  Name        :  SBBMREMEG::getTriggerStep
        // 
        //  Description :  Specifies the start time of the MEG Gradient.
        // 
        //  Parameters  :  none
        // 
        //  Return      :  long
        // 
        // --------------------------------------------------------------------

        long getTriggerStep();


        void setReductionFactor(double ReductFactor);

        double getReductionFactor();

        void setMREProperties(MREProperties);
        MREProperties getMREProperties();

    protected:

        //  ------------------------------------------------------------------ 
        //  Specifies whether the MEG axis should be slice, read, phase or     
        //  none                                                               
        //  ------------------------------------------------------------------ 

        eGradientAxis m_eGradientAxis{Slice};  

        //  -----------------------------------------------------------------
        //  Specifies whether the MEG polarity should have positive or         
        //  negative gradient polarity.                                        
        //  ------------------------------------------------------------------ 

        eMEGPolarity m_eMEGPolarity{Positive};


        //  -----------------------------------------------------------------
        //  Specifies the MEG mode (for improving timing of the sequence and
        //  shortening TE of the sequence                                        
        //  ------------------------------------------------------------------ 

        eMEGMode m_eMEGMode;

        // -----------------------------------------------------
        // Frequency of the MEG
        // -----------------------------------------------------

        double m_dMEGFrequency;

        // -----------------------------------------------------
        // Amplitude of the MEG
        // -----------------------------------------------------

        double m_dMEGAmplitude{-1.0};

        // -----------------------------------------------------
        // Rise Time of the MEG
        // -----------------------------------------------------

        double m_dMinRiseTime{MEG_RISETIME_HARD_DEFAULT};

        // -----------------------------------------------------
        // MEGs might have to be reduced in timing due to low 
        // signal intensity (T2* loss) 
        // Reduction may have different representations/forms
        // this factor represents the reduction corresponding
        // to the chosen reduction mode 
        //  mode 'standard'                 factor is always '1'
        //  mode 'Two_Gradients'            factor describes (orig. MEG(3 pulses))/(two pulses)
        //  mode 'ReduceMEGDuration'        rate of MEG/orig. MEG duration
        //  mode 'ReduceOuterLobDuration'   rate of (MEG pulse 3)/(orig. MEG 3)
        //  mode 'ReduceOuterLobs'          not used at the moment
        // TODO: Zeichnung
        // -----------------------------------------------------

        double m_dReductionFactor{1.0};

        // this factor represents the reduction of MEGs Amplitude
        double m_dAmpReductionFactor{MEG_AMPLITUDE_REDUCTION_FACTOR};

        // ------------------------------------------------------
        // Period of the MEG
        // ------------------------------------------------------

        long m_lMEGPeriod;

        // ------------------------------------------------------
        // Number of trigger steps
        // ------------------------------------------------------

        long m_lNTriggerSteps{MEG_DEFAULT_TRIGGER_STEPS};

        // ------------------------------------------------------
        // Start time of the MEG 
        // ------------------------------------------------------

        long m_lStartTime{-1};


        bool m_bTwoGrad{false};
        bool m_bReduceAll{false};
        bool m_bReduceLobeDuration{false};
        bool m_KeepOuterGradientMoment{false};

        // ------------------------------------------------------
        // The first gradient lobe of the MEG
        // ------------------------------------------------------
        
        sGRAD_PULSE_TRAP m_MEG1;

        // ------------------------------------------------------
        // The second gradient lobe of the MEG
        // ------------------------------------------------------

        sGRAD_PULSE_TRAP m_MEG2;

        // ------------------------------------------------------
        //  The third gradient lobe of the MEG
        // ------------------------------------------------------
        
        sGRAD_PULSE_TRAP m_MEG3;

        // ------------------------------------------------------
        //  MRE Properties of the MEG
        // ------------------------------------------------------
        MREProperties m_MEGMREProperties;

};

// ------------------------------------------------------------
// Set the MEG gradient axis for the SBB
// ------------------------------------------------------------

inline void SBBMREMEG::setMEGGradientAxis (eGradientAxis eAxis)
{
    m_eGradientAxis = eAxis;
}

//-------------------------------------------------------------
// Set the MEG polarity for the SBB
//-------------------------------------------------------------

inline void SBBMREMEG::setMEGPolarity (eMEGPolarity ePolarity)
{
    m_eMEGPolarity = ePolarity;
}

//-------------------------------------------------------------
// Set the MEG frequency for the SBB
//-------------------------------------------------------------

inline bool SBBMREMEG::setMEGFrequency (double dFrequency)
{
    if (dFrequency == 0) {
        return false;
    } else {
        m_dMEGFrequency = dFrequency;
        m_lMEGPeriod = (long) (SEC2USEC/dFrequency);  //MEG period is dependent
    }
    return true;
}

//-------------------------------------------------------------
// Set the start time for the MEG gradient
//-------------------------------------------------------------

inline void SBBMREMEG::setlStartTime (long lTime)
{
    m_lStartTime = lTime;
}

//-----------------------------------------------------------------
// Retrieve the start time for the MEG gradient as known by the SBB
//-----------------------------------------------------------------

inline long SBBMREMEG::getlStartTime ()
{
    return(m_lStartTime);
}

//-----------------------------------------------------------------
// Retrieve the MEG gradient axis as known by the SBB
//-----------------------------------------------------------------

inline enum SBBMREMEG::eGradientAxis SBBMREMEG::getMEGGradientAxis ()
{  
   return (m_eGradientAxis);
}

//------------------------------------------------------------------
// Retrieve the MEG gradient polarity as known by the SBB
//------------------------------------------------------------------

inline enum SBBMREMEG::eMEGPolarity SBBMREMEG::getMEGPolarity ()
{
   return (m_eMEGPolarity);
}

//-----------------------------------------------------------------
// Retrieve the MEG gradient frequency as known by the SBB
//-----------------------------------------------------------------

inline double SBBMREMEG::getMEGFrequency ()
{
   return (m_dMEGFrequency);
}

//------------------------------------------------------------------
// Set the amplitude for the MEG gradient
//------------------------------------------------------------------

inline bool SBBMREMEG::setMEGAmplitude(double dAmplitude)
{
    m_dMEGAmplitude = dAmplitude;
    return (true);
}

//------------------------------------------------------------------
// Set the Rise Time for the MEG gradient
//------------------------------------------------------------------

inline bool SBBMREMEG::setMinRiseTime(double dRiseTime)
{
    m_dMinRiseTime = dRiseTime;
    return (true);
}

// ----------------------------------------------------------------
// Retrieve trigger step size from the SBB
// ----------------------------------------------------------------

inline long SBBMREMEG::getTriggerStep()
{
   return (fSDSRoundUpGRT(m_lMEGPeriod/m_lNTriggerSteps));
}

// ----------------------------------------------------------------
// Sets MEG mode 
// ----------------------------------------------------------------

inline void SBBMREMEG::setMEGMode (eMEGMode eMegMode)
{
    m_eMEGMode = eMegMode;
}

// ----------------------------------------------------------------
// returns the currently set MEG mode
// ----------------------------------------------------------------

inline enum SBBMREMEG::eMEGMode SBBMREMEG::getMEGMode()
{
    return (m_eMEGMode);
}


inline void SBBMREMEG::setReductionFactor(double ReductFactor)
{
    m_dReductionFactor = ReductFactor;
};


inline double SBBMREMEG::getReductionFactor()
{
    return (m_dReductionFactor);
} 

inline void SBBMREMEG::setMREProperties(MREProperties MREProps)
{
    m_MEGMREProperties = MREProps;
    setMEGFrequency(m_MEGMREProperties.getDriverSettings().getFrequency());
    setReductionFactor(m_MEGMREProperties.getFractionalFactor());
    m_dAmpReductionFactor = m_MEGMREProperties.getAmpReductionFactor();
    // set some member values
};

inline MREProperties SBBMREMEG::getMREProperties()
{
    return m_MEGMREProperties;
};
