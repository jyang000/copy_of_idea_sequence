//----------------------------------------------------------------------------------
// <copyright file="SBBCompGrad.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 1998-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#ifndef SBBCompGrad_h
#define SBBCompGrad_h

#include "MrImaging/libSBB/SBBSpoilGrad.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrImagingFW/libBalance/GPABalance.h"

namespace SEQ_NAMESPACE
{
class SeqBuildBlockCompGrad : public SeqBuildBlockSpoilGrad
{
  public:
    // no public default constructor
    SeqBuildBlockCompGrad() = delete;

    SeqBuildBlockCompGrad(SBBList* pSBBList); ///< points to the SBBList object which can be instantiated from the
                                               ///< SBBList class. If you don't use the SBBList, nullptr can be passed.

    virtual ~SeqBuildBlockCompGrad() = default;

    // no copy and move constructor and assignment operator
    SeqBuildBlockCompGrad(const SeqBuildBlockCompGrad& right) = delete;
    SeqBuildBlockCompGrad& operator=(const SeqBuildBlockCompGrad& right) = delete;
    SeqBuildBlockCompGrad(SeqBuildBlockCompGrad&& right)                 = delete;
    SeqBuildBlockCompGrad& operator=(SeqBuildBlockCompGrad&& right) = delete;

    void scaleAmplitude(double dScaleAmpl_x, double dScaleAmpl_y, double dScaleAmpl_z);
    void setCompensationPara(bool bDecayEnable, double dCompFraction, double dTau);
    void setCompGradParaLimit(long lMaxDuration, double dMaxMoment);
    void setdMaxAmplitude(double dMaxAmplitue);
    double getdMaxAmplitude();
    void   setdMinRiseTime(double dMinRiseTime);

    void setbValid(bool bValid);
    bool getbValid();

    //--------------------------------------------------------------------
    //  get full access to internal balance instance (m_sSBBBalance)
    //  tesult is valid only after successful ::prep
    GPABalance getGPALoad();
    void       setGPALoad(GPABalance& sGPALoad);

    bool setbUseGPABalance(bool bUseGPABalance);
    bool getbUseGPABalance() const;

    // GPA load calculation
    bool prepGPALoad(double dAmplitudeX);

    bool prepGPALoad(double dAmplitudeX, double dAmplitudeY, double dAmplitudeZ);

    bool scaleGPALoad(double dScaleX);

    bool scaleGPALoad(double dScaleX, double dScaleY, double dScaleZ);

    // calculate the duration of compensation gradients based on the max amplitude of diffusion gradients
    void calcDiffCompGrad(double dMaxDiffAmplitude, long lDiffGradDuration, long lDiffGradSpacing, long lDiffToCompGrad);

   private:
    bool   m_bConsiderCompensationDecay{false};
    double m_dCompensationFraction{1.0};
    double m_dEddyCurrentTau{150.0};
    double m_dMaxMagnitude{10.0};
    double m_dMinRiseTime{20.0};
    double m_dMaxCompGradMonment{1000000.0};
    bool   m_bValid{false};
    bool   m_bUseGPABalance{false};
    GPABalance m_sBalanceCompGrad;

    /// ID's of gradient events stored within m_sBalanceCompGrad
    long m_lIDEX{0};
    long m_lIDEY{0};
    long m_lIDEZ{0};
};

inline void SeqBuildBlockCompGrad::setCompensationPara(bool bDecayEnable, double dCompFraction, double dTau)
{
    m_bConsiderCompensationDecay  = bDecayEnable;
    m_dCompensationFraction      = dCompFraction;
    m_dEddyCurrentTau     = dTau;
}

inline void SeqBuildBlockCompGrad::setCompGradParaLimit(long lMaxTotalTime, double dMaxMoment)
{
    setMaxGradientDuration(lMaxTotalTime);
    m_dMaxCompGradMonment  = dMaxMoment;
}

inline void SeqBuildBlockCompGrad::setbValid(bool bValid)
{
    m_bValid = bValid;
}

inline bool SeqBuildBlockCompGrad::getbValid()
{
    return m_bValid;
}

inline void SeqBuildBlockCompGrad::setdMaxAmplitude(double dMaxAmplitue)
{
    m_dMaxMagnitude = dMaxAmplitue;
}

inline double SeqBuildBlockCompGrad::getdMaxAmplitude()
{
    return m_dMaxMagnitude;
}

inline void SeqBuildBlockCompGrad::setdMinRiseTime(double dMinRiseTime)
{
    m_dMinRiseTime = dMinRiseTime;
}

inline bool SeqBuildBlockCompGrad::setbUseGPABalance(bool bUseGPABalance)
{
    // check whether actual system supports balance models
    if (bUseGPABalance && (m_sBalanceCompGrad.lGetStatus() != 0))
    {
        m_bUseGPABalance = false;

        return false;
    }

    if (bUseGPABalance != m_bUseGPABalance)
    {
        m_bUseGPABalance = bUseGPABalance;
        resetPrepared();
    }
    return true;
}

//--------------------------------------------------------------------
//  check whether use of GPA balance model is enabled
inline bool SeqBuildBlockCompGrad::getbUseGPABalance() const
{
    return m_bUseGPABalance;
}

//--------------------------------------------------------------------
//  get full access to internal balance instance (m_sSBBBalance)
inline GPABalance SeqBuildBlockCompGrad::getGPALoad()
{
    return m_sBalanceCompGrad;
}

//--------------------------------------------------------------------
// set GPA load of compensation gradient events
inline void SeqBuildBlockCompGrad::setGPALoad(GPABalance& sKernelBalance)
{
    m_sBalanceCompGrad = sKernelBalance;
}

//--------------------------------------------------------------------
// reset the amplitude of compensation gradient events
inline void SeqBuildBlockCompGrad::scaleAmplitude(double dScaleAmpl_x, double dScaleAmpl_y, double dScaleAmpl_z)
{
    if (fabs(dScaleAmpl_x) < 1)
    {
        m_GPX->setAmplitude(dScaleAmpl_x * m_dMaxMagnitude);
    }

    if (fabs(dScaleAmpl_y) < 1)
    {
        m_GPY->setAmplitude(dScaleAmpl_y * m_dMaxMagnitude);
    }

    if (fabs(dScaleAmpl_z) < 1)
    {
        m_GPZ->setAmplitude(dScaleAmpl_z * m_dMaxMagnitude);
    }
}

} // namespace SEQ_NAMESPACE

#endif // #ifndef SBBSpoilGrad_h
