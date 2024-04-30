//-----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 2010  All Rights Reserved.  Confidential
//-----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//
//    File: \n4_servers1\pkg\MrServers\MrImaging\libSeqPTX\C2DGradientsBlippedEPI.h
//
//  Author: pfeujodj  
//
//    Lang: C++
//
// Descrip: SBB for 2DRF Excitation
//
//-----------------------------------------------------------------------------


#ifndef C2DGradientsBlippedEPI_h
#define C2DGradientsBlippedEPI_h 1

#include "MrImagingFW/libSBBFW/SeqBuildBlock.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"	// sGRAD_PULSE
#include "MrImaging/libSeqUtil/SUForbiddenTR.h" // SUForbidden
#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"      // for fSDSRoundDownGRT, ...
#include "MrImaging/libSBB/libSBBmsg.h"         // for MRI_SBB_SBB_ERROR etc.
#include "MrProtSrv/Domain/MrProtData/MrProt/MeasParameter/MrTXSpec.h"  // MrProt Logic

//#define SBB2DGradients_MAX_NO_OF_GRAD_CONCATS   11

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

class __IMP_EXP C2DGradientsBlippedEPI : public SeqBuildBlock
{

public:
	static const long SBB2DGradientsPerfBlip;
	static const long SBB2DGradientsPerfLine;
	static const long SBB2DGradientsPerfRefocus;

	C2DGradientsBlippedEPI(SBBList* pSBBList);

	virtual ~C2DGradientsBlippedEPI() = default;

	// No copy and move operations
	C2DGradientsBlippedEPI(const C2DGradientsBlippedEPI& right) = delete;
	C2DGradientsBlippedEPI & operator=(const C2DGradientsBlippedEPI& right) = delete;
	C2DGradientsBlippedEPI(C2DGradientsBlippedEPI&& right) = delete;
	C2DGradientsBlippedEPI & operator=(C2DGradientsBlippedEPI&& right) = delete;

	//	prepares trajectory and (if requested) refocussing gradients
	bool prep(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;

	//	runs trajectory's gradients and (if requested) refocussing gradients.
	//	RTEBInit(.) and RTEBFinish(.) must be called before and after run(.). 
	// (N.B.	this is done in SeqBuildBlock2DExcitation...)
	bool run(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC = nullptr) override;

	//	get gradient amplitude in [mT/m]   for time lTime in [us] from start of	trajectory.
	void getG(long lTime, double& pdG1, double& pdG2);

	//	returns the factor to get from gradient moments (G) to K-space positions.
	double getMomentToKFactor();

	//	gets the start position in K-space
	void getKStart(double& dK1, double& dK2);

	//	gets the end position in K-s-pace (also used to determine	the refocussing moment)
	void getKEnd(double& dK1, double& dK2);

	//	set distance of side peaks in [mm]
	void setSidelobeDistance(int nDistance);

	//	set resolution in direction X1 in [mm]
	void setReadResolution(double dReadResolution);

	//	set resolution in direction X2 in [mm]
	void setPhaseResolution(double dPhaseResolution);

	// set rotation angle in rad
	void setRotationAngleRad(double dRotationAngle);
	double getRotationAngleRad();

	//	returns a correction factor, depending upon the absolute time ltime	from the start of the trajectory.  
	// This factor may be used to suppress the rf e.g. during switching of blips, or refocusser.
	double getDensityCorrection(long lTime, double dK1, double dK2);

	long getStartTimeInEventBlock();
	void setStartTimeInEventBlock(long lStartTime);
	//void setStartTimeInEventBlockNoResetPrepared(long lStartTime);

	//	returns the echo time contribution of trajectory
	long getTrajectoryTEContribution();

	//	Set gyromagnetic ratio and the MomentToKFactor. Necessary for upward compatibility because the nucleus can soon be selected in the UI.
	void setGamma(MrProt &rMrProt);

	//	returns m_bSelfRefocussed
	bool isSelfRefocussing();

	//	get gradient amplitude in mT/m for time lTime in us from start of refocussing gradient
	void getRefocussingGrad(long lTime, double& pdG1, double& pdG2);

	long getRefocussingTotalTime();

	// If true rf is played out during gradient ramping
	void setRampSampling(bool bFlag);

	long getTargetExcitationEchoSpacing();

	long getCurrentExcitationEchoSpacing();

	long getBlippedEPINumberOfKLines();
	long getBlippedEPINumberOfKLinesNotAccel();

	void setPulseAcceleration(double dCurPulseAcceleration);
	double getPulseAcceleration() const;

	// an asymmetric rf pulse can be generated: default PFFactor is 8/8
	void setPulsePartialFourierFactor(double dPulsePFFactor);

	// retrieve actual asymmetry of rf pulse
	double getAsymmetry() const;

	void set2DGradMaxAmpl(double dGradMaxAmpl);     // [%] rel. scaling

	//	sets m_bSelfRefocussed. If true, the refocussing gradients are played out by the run(.) function and getRefocussingMoments(.) 
	// returns 0.
	void setSelfRefocussing(bool bFlag);

	//	returns the time need by the refocussing gradients in [us].
	//	If m_bSelfRefocussed =true, 0 is returned.
	long getRefocussingDuration();

	//	returns the time need by the trajectory without refocussing gradients in [us]
	long getTrajectoryDuration();

	///  Can be used to tell the SBB to use negative GS-gradients.
	///     Must be called before prep.
	virtual void setRequiredGSPolarity(double dPolarity);    ///<  must be < 0 to force the SBB to use negative GS-gradients

	// access to gradient objects
	sGRAD_PULSE getXLine();
	sGRAD_PULSE getXBlip();

	//	updates trajectory's gradient amplitudes and risetimes.
	void updateGradientPerformanceTrajectory(SEQ::Gradients eGradMode);


protected:
	// for asymmetric trajectories: NumberOfKLines = FirstHalf + 1 CenterLine + SecondHalf
	int getNumberOfKLines();

private:
	sGRAD_PULSE m_XLine{"XLine"};
	sGRAD_PULSE m_XBlip{"XBlip"};

	// gradient objects required for rotated trajectory
	sGRAD_PULSE m_XLinePhase{"XLinePhase"};
	sGRAD_PULSE m_XLineSlice{"XLineSlice"};

	sGRAD_PULSE m_XBlipPhase{"XBlipPhase"};
	sGRAD_PULSE m_XBlipSlice{"XBlipSlice"};

	sGRAD_PULSE m_rK1{"rK1"};
	sGRAD_PULSE m_rK2{"rK2"};

	//	distance of side wiggles in [mm]
	int m_nSidelobeDist{ 0 };

	//	Resolution in read direction in [mm]
	double m_dResolutionRead{ 1. };

	//	Resolution in phase direction in [mm]
	double m_dResolutionPhase{ 1. };
	double m_dRotationAngleRad{ 0. };

	double m_dK1Start{ 0. };
	double m_dK1End{ 0. };
	double m_dK2Start{ 0. };
	double m_dK2End{ 0. };

	//	start time in event block in [us]
	long m_lStartTimeInRTEB{ 0 };

	//	gyromagnetic ratio, necessary for upward compatibility because the nucleus can soon be selected in the UI
	//	and therefore the access to the gamma through pMRProt->...	is performed only once and stored locally.
	double m_dGamma{ 0. };

	bool m_bSelfRefocussed{ false };

	//	factor to get from gradient switching  to K-space position change
	double m_dMomentToKFactor{ 0. };

	// Factor to invert GS-gradient polarity (+1.0/-1.0)
	double m_dRequiredGSPolarity{ 1.0 };

	double m_dGradMaxAmpl{ 100.0 };

	long m_lTargetExcitationEchoSpacing{ 0 };
	long m_lCurrentExcitationEchoSpacing{ 0 };    // in case of forbidden ranges, Current can be different from target

	long m_lBlippedEPINumberOfKLines{ 0 };
	long m_lBlippedEPINumberOfKLinesNotAccel{ 0 };

	double m_dPulseAcceleration{ 1.0 };

	double m_dPulsePartialFourierFactor{ 1.0 };      // input value, ie. planned asymmetry
	double m_dAsymmetry{ 0.0 };          // calculated value

	// for asymmetric trajectories: NumberOfKLines = FirstHalf + 1 CenterLine + SecondHalf
	int m_nNumberOfKLinesFirstHalf{ 0 };
	int m_nNumberOfKLinesSecondHalf{ 0 };
	bool m_bRampSampling{ true };

	// Class handling forbidden echo spacings
	SUForbidden m_ForbiddenEchoSpacings{};

	// array [2,N] of ForbiddenEchoSpacing bands: minimum/maximum esp
	std::vector<std::vector<int> > m_vi2DEspMinMax{};

	// initializes forbidden bands from meas perm and system specific settings
	void initArrayForbiddenEchoSpacing();

	void getKSpacePosition(const double& dIn1, const double& dIn2, double& dOut1, double& dOut2);

};


inline sGRAD_PULSE C2DGradientsBlippedEPI::getXLine()
{
	return m_XLine;
}


inline sGRAD_PULSE C2DGradientsBlippedEPI::getXBlip()
{
	return m_XBlip;
}


//	returns the factor to get from gradient moments (G) to K-space positions.
inline double C2DGradientsBlippedEPI::getMomentToKFactor()
{
	return m_dMomentToKFactor;
}

//	set distance of side peaks in [mm]
inline void C2DGradientsBlippedEPI::setSidelobeDistance(int nDistance)
{
	m_nSidelobeDist = nDistance;
}

//	set resolution in direction X1 in [mm]
inline void C2DGradientsBlippedEPI::setReadResolution(double dReadResolution)
{
	m_dResolutionRead = dReadResolution;
}

//	set resolution in direction X2 in [mm]
inline void C2DGradientsBlippedEPI::setPhaseResolution(double dPhaseResolution)
{
	m_dResolutionPhase = dPhaseResolution;
}

//	set rotation angle in rad
inline void C2DGradientsBlippedEPI::setRotationAngleRad(double dRotationAngle)
{
	m_dRotationAngleRad = dRotationAngle;
}


inline double C2DGradientsBlippedEPI::getRotationAngleRad()
{
	return m_dRotationAngleRad;
}


inline long C2DGradientsBlippedEPI::getStartTimeInEventBlock()
{
	return m_lStartTimeInRTEB;
}


inline void C2DGradientsBlippedEPI::setStartTimeInEventBlock(long lStartTime)
{
	m_lStartTimeInRTEB = lStartTime;
}

//	returns the echo time contribution of trajectory
inline long C2DGradientsBlippedEPI::getTrajectoryTEContribution()
{
	return static_cast<long>(floor((static_cast<double>(m_nNumberOfKLinesSecondHalf) + 0.5)*static_cast<double>(m_XLine.getTotalTime())) + getRefocussingDuration());
}


//	Set gyromagnetic ratio and the MomentToKFactor. Necessary for upward compatibility because the nucleus can soon be selected in  the UI.
inline void C2DGradientsBlippedEPI::setGamma(MrProt &rMrProt)
{
	MeasNucleus myMeasNucleus(rMrProt.txSpec().nucleusInfoArray()[0].gettNucleus().c_str()); // get the selected nucleus
	m_dGamma = myMeasNucleus.getLarmorConst() / 1000.; //  in 1/(us T) 

	//------------------------
	// set gamma2K
	//------------------------
	m_dMomentToKFactor = m_dGamma * 2.* M_PI; // factor to get from G to K
}


inline bool C2DGradientsBlippedEPI::isSelfRefocussing()
{
	return (m_bSelfRefocussed);
}


inline long C2DGradientsBlippedEPI::getRefocussingTotalTime()
{
	// assumption: m_rK1 / m_rK2 have same duration: is secured in prepRefocussingGradients()
	return (m_rK1.getTotalTime());
}


//	sets m_bSelfRefocussed. If true, the refocussing gradients are played out by the run(.) function and getRefocussingMoments(.) returns 0.
inline void C2DGradientsBlippedEPI::setSelfRefocussing(bool bFlag)
{
	m_bSelfRefocussed = bFlag;
}


inline long C2DGradientsBlippedEPI::getTargetExcitationEchoSpacing()
{
	return m_lTargetExcitationEchoSpacing;
}


inline long C2DGradientsBlippedEPI::getCurrentExcitationEchoSpacing()
{
	return m_lCurrentExcitationEchoSpacing;
}


inline long C2DGradientsBlippedEPI::getBlippedEPINumberOfKLines()
{
	return m_lBlippedEPINumberOfKLines;
}


inline long C2DGradientsBlippedEPI::getBlippedEPINumberOfKLinesNotAccel()
{
	return m_lBlippedEPINumberOfKLinesNotAccel;
}


inline void C2DGradientsBlippedEPI::setPulseAcceleration(double dCurPulseAcceleration)
{
	// set trajectory under/oversampling factor
	// must be not zero!
	m_dPulseAcceleration = std::max(0.001, dCurPulseAcceleration);
}


inline double C2DGradientsBlippedEPI::getPulseAcceleration() const
{
	return m_dPulseAcceleration;
}


// an asymmetric rf pulse can be generated: default PFFactor is 8/8
inline void C2DGradientsBlippedEPI::setPulsePartialFourierFactor(double dPulsePFFactor)
{
	// set pulse PF factor: must be between 0.5 and 1.0!
	m_dPulsePartialFourierFactor = std::max(0.5, std::min(1.0, dPulsePFFactor));
	//m_dPulsePartialFourierFactor = std::clamp(dPulsePFFactor, 0.5, 1.0);
}


// retrieve actual asymmetry of rf pulse: successful prep() necessary
inline double C2DGradientsBlippedEPI::getAsymmetry() const
{
	return m_dAsymmetry;
}


inline void C2DGradientsBlippedEPI::setRampSampling(bool bFlag)
{
	m_bRampSampling = bFlag;
}


inline void C2DGradientsBlippedEPI::set2DGradMaxAmpl(double dGradMaxAmpl)
{
	m_dGradMaxAmpl = dGradMaxAmpl;
}


//	returns the time need by the refocussing gradients in [us].	If m_bSelfRefocussed =true, 0 is returned.
inline long C2DGradientsBlippedEPI::getRefocussingDuration()
{
	if (isSelfRefocussing())
	{
		return (maximum(m_rK1.getTotalTime(), m_rK2.getTotalTime()));
	}
	else
	{
		return (0);
	}
}


//	returns the time need by the trajectory without refocussing gradients in [us]
inline long C2DGradientsBlippedEPI::getTrajectoryDuration()
{
	return m_lBlippedEPINumberOfKLines * m_XLine.getTotalTime();
}


inline void C2DGradientsBlippedEPI::setRequiredGSPolarity(double dPolarity)
{
	if (dPolarity*m_dRequiredGSPolarity < 0.0)
	{
		if (dPolarity < 0)
			m_dRequiredGSPolarity = -1.0;
		else
			m_dRequiredGSPolarity = 1.0;

		resetPrepared();
	}
}

inline int C2DGradientsBlippedEPI::getNumberOfKLines()
{
	return m_nNumberOfKLinesFirstHalf + 1 + m_nNumberOfKLinesSecondHalf;
}

#endif
