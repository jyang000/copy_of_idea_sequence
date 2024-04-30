//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens Healthcare GmbH 2016  All Rights Reserved.
//	-----------------------------------------------------------------------------

#pragma once

#include "MrMeasSrv/MeasNuclei/Base/MeasNucleus.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sSLICE_POS.h"
#include "MrProtSrv/Domain/CoreNative/GradSpecDefines.h"

#include <algorithm>
#include <numeric>
#include <vector>

class MaxwellCorrection
{
  public:
    // constructor
    MaxwellCorrection() = default;

    // default destructor
    virtual ~MaxwellCorrection() = default;

    MaxwellCorrection(const MaxwellCorrection& right) = delete;
    MaxwellCorrection& operator=(const MaxwellCorrection& right) = delete;
    MaxwellCorrection(MaxwellCorrection&& right)                 = delete;
    MaxwellCorrection& operator=(MaxwellCorrection&& right) = delete;

    // setting parameters from slice position and single gradient
    void setParameters(sSLICE_POS const* pSLC, sGRAD_PULSE const& sGradient, SEQ::GradientAxis const& eAxis);

    // setting any number of gradients, with the missing ones substituted by 0 in the argument (like in fRTEI)
    void setParameters(
        sSLICE_POS const*  pSLC,
        sGRAD_PULSE const* pGradientP,
        sGRAD_PULSE const* pGradientR,
        sGRAD_PULSE const* pGradientS);

    // set parameters independent of gradients
    void setBasicParameters(sSLICE_POS const* pSLC);

    // calculate the maxwell terms of the input gradient
    void calcMaxwellCorrection(bool bRampCorrection = true);

    // prepare the correction gradient of the required axis
    void prepMaxwellCorrectionGradient(sGRAD_PULSE& sCorrectionGradient, SEQ::GradientAxis const& eAxis) const;

    // calc frequency and phase modifications due to maxwell terms and additional gradient, depending on the larmor freq
    // (B0 and nucleus)
    void calcFrequencyAndPhase(
        MeasNucleus const& Nucleus, sGRAD_PULSE const& sCorrectionGradient, SEQ::GradientAxis const& eAxis);

    // getting the phase increments
    void getPhaseIncrements(double& dMaxwellPhaseIncrement, double& dMaxwellTotalPhaseIncrement) const;

    // get resulting correction gradient amplitudes, separately for each axis
    double getMaxwellGradAmplP() const;

    double getMaxwellGradAmplR() const;

    double getMaxwellGradAmplS() const;

  private:
    // system-dependent variables
    double m_dAlpha{0.0};

    double m_dPsiX{0.0};

    double m_dPsiY{0.0};

    double m_dNominalB0{0.0};

    // constant and linear terms
    double m_dMaxwellB0{0.0};

    // input gradient amplitudes in PRS and XYS
    std::vector<double> m_vdGradAmpl{0.0, 0.0, 0.0};

    std::vector<double> m_vdGradAmplXYZ{0.0, 0.0, 0.0};

    // input gradient other properties
    bool m_bIsGradientPrepared{false};

    long m_lRampUpTime{0};

    long m_lRampDownTime{0};

    long m_lFlatTopTime{0};

    long m_lDuration{0};

    // input bool for ramp correction
    bool m_bRampCorrection{false};

    // input position data, in PRS and XYZ
    std::vector<double> m_vdRefPosition{0.0, 0.0, 0.0};

    std::vector<double> m_vdRefPositionXYZ{0.0, 0.0, 0.0};

    // input rotation matrix
    std::vector<std::vector<double>> m_vvdRotMatrix;

    // inverse of rotation matrix
    std::vector<std::vector<double>> m_vvdRotMatrixInv;

    // resulting maxwell correction terms
    std::vector<double> m_vdMaxwellTermsXYZ{0.0, 0.0, 0.0}; // X, Y, Z, [mT/m]

    // resulting gradient amplitude in PRS
    std::vector<double> m_vdMaxwellGradAmpl{0.0, 0.0, 0.0};

    // output for NCO and ICE
    double m_dMaxwellPhaseIncrement{0.0};

    double m_dMaxwellTotalPhaseIncrement{0.0};
};

// function for matrix-vector multiplication
template<class T>
std::vector<T> MatrixVectorMultiply(std::vector<std::vector<T>> const& vMatrix, std::vector<T> const& vVector)
{
    std::vector<T> vRes{vVector};

    std::transform(std::begin(vMatrix), std::end(vMatrix), std::begin(vRes), [vVector](std::vector<T> row) {
        return std::inner_product(std::begin(row), std::end(row), std::begin(vVector), 0.0);
    });

    return vRes;
}

// function to copy from 2D C array to vector of vectors.
// this would be more elegant initializing the index variable within the lambda capture, but the linux compiler cannot
// handle it currently.

std::vector<std::vector<double>> FromRotMatrixToVectorOfVectors(sROT_MATRIX const& sRotMatrix);

// function to get index from axis enum
long indexFromGradAxisEnum(SEQ::GradientAxis const& eAxis);
