//-----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 2012  All Rights Reserved.  Confidential
//-----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//
//    File: \n4_servers1\pkg\MrServers\MrImaging\libSeqPTX\RF2DArbGenerator.h
//
//  Author: Josef Pfeuffer  
//
//    Lang: C++
//
// Descrip: Container for pTX RF Pulse
//      
//-----------------------------------------------------------------------------


#ifndef RF2DArbGenerator_h
#define RF2DArbGenerator_h 1


#include <vector>
#include "MrMeasSrv/SeqIF/libRT/libRTDefines.h" // for IRF_MAX_ENVELOPE_ELEMENTS (max. rf sample points)
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"    // sRF_PULSE
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"

class sSLICE_POS;
class C2DGradientsBlippedEPI;

#define   DEF_FILTER_COEFF  0.15 

//-----------------------------------------------------------------------------
// import/export control
//-----------------------------------------------------------------------------
#ifdef __IMP_EXP
#undef __IMP_EXP
#endif 

#ifdef LOCAL_PTX_BUILD
#define __IMP_EXP
#pragma message( "Local PTX build" )
#else  //   LOCAL_PTX_BUILD not defined
#ifdef BUILD_libSeqPTX
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control
#endif  //  LOCAL_PTX_BUILD


class  __IMP_EXP  RF2DArbGenerator
{
public:
	// hardcoded settings for product version 
	static const long DefaultSideLobeDistance; // [mm]
	static const double Default2DExciteRelativeExcitationSize;
	static const double Default2DExciteResolutionPhase; // [mm]
	static const double Default2DExciteFilterCoefficient;
	static const double Default2DExciteGradientMaxAmplitudeValue;      // [%] rel. scaling 
	static const double Default2DExciteSliceBWT;
	static const double Default2DExciteRotationAngleDeg;

	struct RFFLYBACK   // 
	{
		enum RFFlyBackType
		{
			PSEUDO = 1,    // pseudo flyback, by applying double FoE
			FULL = 2     // every odd excitation line it truely *nulled*
		};
	};

	//-------------------------------------------------------------
	//  configuration object
	struct RF2DArbParameter
	{
		RF2DArbParameter() = default;

		bool   bIsX2drfDumpMode{ false };
		long   lSidelobeDistance{ DefaultSideLobeDistance };	// recalculated with UI input
		double dFoERelativeSize{ Default2DExciteRelativeExcitationSize };	// recalculated with UI input
		double dRelExcitationSize{ Default2DExciteRelativeExcitationSize };
		double dXReadResolution{ 0.0 };	// recalculated: RESOLUTION_FREQUENCY_mm
		double dXPhaseResolution{ Default2DExciteResolutionPhase };
		double dXReadBWT{ Default2DExciteSliceBWT };
		double dPulseAcceleration{1.0};
		double dPulsePartialFourierFactor{1.0};
		double dFlipAngle{90.0};	// recalculated with MrProt input 
		double dRotationAngleDeg{ Default2DExciteRotationAngleDeg };
		double dFilterCoeff1{0.0};	// !! do not filter in the slice select direction (XRead)
		double dFilterCoeff2{ Default2DExciteFilterCoefficient };
		RFFLYBACK::RFFlyBackType eFlyBackMode{ RFFLYBACK::PSEUDO };
		double dAdditionalPhase{0.0};	// offset by adding a phase to each sample of the rf envelope
		bool   bRampSampling{ true };
		double dGradMaxAmpl{ Default2DExciteGradientMaxAmplitudeValue };
		double dRequiredGSPolarity{+1.0};
	};

public:
	RF2DArbGenerator();
	virtual ~RF2DArbGenerator();

	RF2DArbGenerator(const RF2DArbGenerator& right) = delete;
	RF2DArbGenerator & operator=(const RF2DArbGenerator& right) = delete;
	RF2DArbGenerator(RF2DArbGenerator&& right) = delete;
	RF2DArbGenerator & operator=(RF2DArbGenerator&& right) = delete;

	// only access to member: memory is not to be deleted from outside!
	virtual sRF_PULSE_ARB* getArbitraryRFPulse() const;

	// only access to member: memory is not to be deleted from outside!
	virtual C2DGradientsBlippedEPI* getTrajectory() const;

	// Sets all generally used parameters like sidelobedistance,	resolution
	bool initialize(const RF2DArbParameter &sParam);

	// generate / prepare the rf pulse
	//    initialize() is necessary before!
	virtual bool prepare(MrProt &rMrProt, const SeqLim &rSeqLim, SeqExpo &rSeqExpo);

	void dumpClassInfo() const;

	void setPrepared();
	void resetPrepared();
	bool isPrepared() const;

	//	set context to state, which enables call by calculatePTX()
	void setContextExportForCalculatePtx(bool bContextExportForCalculatePtx);

	// Get functions to retrieve (unnormalized) RF and gradient information of the 2D pulse 
	void getRFAmplitude(std::vector<float>& returnVector) const;
	void getRFPhase(std::vector<float>& returnVector) const;
	void getRFDensity(std::vector<int>  & returnVector) const;
	void getGrad1(std::vector<float>& returnVector) const;
	void getGrad2(std::vector<float>& returnVector) const;
	long getRFRasterTime() const;

	///  Can be used to tell the SBB to use negative GS-gradients. Must be called before prep.
	virtual void setRequiredGSPolarity(double dPolarity);    ///<  must be < 0 to force the SBB to use negative GS-gradients
	virtual double getRequiredGSPolarity() const;

	virtual std::string getToolTipInfo() const;       // returns the info string characterizing the 2DRF pulse parameters

	// supported access to SeqBuildBlock type parameters
	long getDurationPerRequest() const;
	long getTEContribution() const;
	MrProtocolData::SeqExpoRFInfo getRFInfoPerRequest() const;

protected:
	// dump into INI file using internal class INI writer
	virtual bool dumpINIFile() const;

	//	calculates the rf pulse for the 2d excitation
	bool calculateRF(MrProt &rMrProt, sSLICE_POS* pSLC = nullptr);

	//	calculates the dwell time of the rf pulse. Can be 1, 2, 5 or 10us
	bool calcRFRasterTime();


private:
	// this is the central ARB pulse of this class
	sRF_PULSE_ARB m_ArbitraryRFPulse{ "ArbRf2DExc" };

	// this is the central gradient of the TX trajectory initialized by setTrajectory()
	C2DGradientsBlippedEPI *	  m_pTrajectory;
	//-------------------------------------------------------------------------------------------------------------------------------------------

	// supported SeqBuildBlock type parameters
	long    m_lSBBDurationPerRequest_us{ 0 };       ///    Time needed for one execution of the SBB-Run function
	long    m_lTEContribution{ 0 };       // TE contribution of RF and Grad as joint object (relative to SBBDurationPerRequest_us)
	MrProtocolData::SeqExpoRFInfo m_RFInfoPerRequest{};   /// detailed Energy applied during one execution of the run function
	double  m_dFlipAngle_deg{ 90 };
	bool    m_bIsPrepared{ false };

	//	'Complex' array that holds the rf data.	Uses maximal available points for an rf pulse structure.
	sSample m_asSample[IRF_MAX_ENVELOPE_ELEMENTS]{};

	//	Effective amplitude integral. Necessary to calculate the rf power.
	double m_dEffAmplInt{ 0.0 };
	bool m_bGradientsOK{ false };

	bool    m_bIsClassInfoDumpDone{ false };
	bool    m_bIsDumpMode{ false };	// no dump to INI if not explicitely set

	//	Dwell time of the rf pulse, can be 1, 2, 5, 10 or 20.
	long    m_lRFRasterTime{ 0 };
	double  m_dGradMaxAmpl{ 100.0 };	// [%] factor to scale Exc. gradient
	double  m_dRelExcitationSize{ 1.0 };	// factor to scale phaseFOV
	RFFLYBACK::RFFlyBackType m_eFlyBackMode{ RFFLYBACK::PSEUDO };       // default: PSEUDO, if FULL: nulling of RF every odd Exc. line

	double m_dFilterCoefficient{ DEF_FILTER_COEFF };	//	Filter coefficient [is between 0 and 1]
	double m_dGaussExponent{ 0.0 };						//	gauss exponent

	// Factor to invert GS-gradient polarity (+1.0/-1.0)
	double m_dRequiredGSPolarity{ +1.0 };

	std::vector<float> m_vfRFAmplitude{};     // Holds the complete RF amplitude information for the pulse 
	std::vector<float> m_vfRFPhase{};         // Holds the complete RF phase information for the pulse 
	std::vector<float> m_vfRFDensity{};       // Holds the complete RF density information for the pulse 
	std::vector<int>   m_viRFDensity{};       // Holds the complete (discrete) RF density information for the pulse 
	std::vector<float> m_vfGrad1{};           // Holds the complete readout gradient information for the pulse 
	std::vector<float> m_vfGrad2{};           // Holds the complete phase   gradient information for the pulse 

	bool m_bContextExportForCalculatePtx{ false };

	const double getSinc(double dK, double dSize) const;
};


inline sRF_PULSE_ARB* RF2DArbGenerator::getArbitraryRFPulse() const
{
	return &m_ArbitraryRFPulse;
}


inline C2DGradientsBlippedEPI* RF2DArbGenerator::getTrajectory() const
{
	return m_pTrajectory;
}


inline void RF2DArbGenerator::setPrepared()
{
	m_bIsPrepared = true;
}


inline void RF2DArbGenerator::resetPrepared()
{
	m_bIsPrepared = false;
}


inline bool RF2DArbGenerator::isPrepared() const
{
	return m_bIsPrepared;
}

inline void RF2DArbGenerator::setContextExportForCalculatePtx(bool bContextExportForCalculatePtx)
{
	m_bContextExportForCalculatePtx = bContextExportForCalculatePtx;
}


inline void RF2DArbGenerator::getRFAmplitude(std::vector<float>& returnVector) const
{
	returnVector = m_vfRFAmplitude;
}


inline void RF2DArbGenerator::getRFPhase(std::vector<float>& returnVector) const
{
	returnVector = m_vfRFPhase;
}


inline void RF2DArbGenerator::getRFDensity(std::vector<int>& returnVector) const
{
	returnVector = m_viRFDensity;   // !! Discrete values
}


inline void RF2DArbGenerator::getGrad1(std::vector<float>& returnVector) const
{
	returnVector = m_vfGrad1;
}


inline void RF2DArbGenerator::getGrad2(std::vector<float>& returnVector) const
{
	returnVector = m_vfGrad2;
}

//returns the RF dwell time
inline long RF2DArbGenerator::getRFRasterTime() const
{
	return (m_lRFRasterTime);
}

inline void RF2DArbGenerator::setRequiredGSPolarity(double dPolarity)
{
	if (dPolarity*m_dRequiredGSPolarity < 0.0)
	{
		if (dPolarity < 0) m_dRequiredGSPolarity = -1.0;
		else             m_dRequiredGSPolarity = 1.0;
		resetPrepared();
	}
}

inline double RF2DArbGenerator::getRequiredGSPolarity() const
{
	return m_dRequiredGSPolarity;
}

//    Returns the duration (in us) of one execution of the SBB run(.) function.
//    The duration is 0 if the m_bPrepared flag has not been set.
inline long RF2DArbGenerator::getDurationPerRequest() const
{
	if (isPrepared())
	{
		return m_lSBBDurationPerRequest_us;
	}
	else
	{
		return 0;
	}
}


inline long RF2DArbGenerator::getTEContribution() const
{
	if (isPrepared())
	{
		return m_lTEContribution;
	}
	else
	{
		return 0;
	}
}


//    Returns by value a structure with the complete energy information of one
//    execution of the SBB run(.) function.
inline MrProtocolData::SeqExpoRFInfo RF2DArbGenerator::getRFInfoPerRequest() const
{
	if (isPrepared())
	{
		return m_RFInfoPerRequest;
	}
	else
	{
		return MrProtocolData::SeqExpoRFInfo();
	}
}



#endif  // RF2DArbGenerator_h
