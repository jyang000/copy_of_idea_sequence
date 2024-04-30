//----------------------------------------------------------------------------------
// <copyright file="PaceFeedback.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 1998-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015-2021. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#include "MrImaging/seq/Kernels/SBBEPIKernel.h"
#include "MrMeasSrv/SeqIF/libRT/SEQSemaphore.h"

#include <array>

// define maximum number of measurements (repet.s+1)
#define MAX_NO_MEAS 4096
#ifndef ASL
#define TWAKEUP 20000 // Wakeuptime for realtime stuff
#else
#define TWAKEUP 70000 // Increased Wakeuptime for support for FOCI pulse update
#endif

namespace SEQ_NAMESPACE
{
class PaceFeedback
{
  public:
    PaceFeedback() = default;

    virtual ~PaceFeedback() = default;

    PaceFeedback(const PaceFeedback& right) = delete;
    PaceFeedback& operator=(const PaceFeedback& right) = delete;
    PaceFeedback(PaceFeedback&& right)                 = delete;
    PaceFeedback& operator=(PaceFeedback&& right) = delete;

    void Init();

    NLS_STATUS Prep(MrProt& rMrProt, sSLICE_POS* asSLC);

    // use feedback parameters for next meas
    NLS_STATUS SyncAndIncorporateFeedback(
        MrProt&                 rMrProt,
        SeqLim&                 rSeqLim,
        SeqExpo&                rSeqExpo,
        sSLICE_POS*             asSLC,
        SeqBuildBlockEPIKernel* pEPIKernel,
        long                    lKernelMode);

    // provide pointer to semaphore
    SEQSemaphore* GetFeedbackSemaphore();

    // sets current FB Data variables
    void SetCurrFBData(long lCurrFBRepetNo, double pdCurrFBRotMatLog[3][3], double pdCurrFBTransLog[3]);

  protected:
    void MatMulVec3D(double pdResultVec[3], double pdMatrix[3][3], double pdVector[3]);

    void MatMulMat3D(double pdResult[3][3], double pdIn1[3][3], double pdIn2[3][3]);

    void MatCopy3D(double pdOut[3][3], double pdIn[3][3]);

    void TransposeMat3D(double pdResult[3][3], double pdIn[3][3]);

    void VecAddVec3D(double pdResultVec[3], double pdInVec1[3], double pdInVec2[3]);

    bool m_bInitPerformed{false};

    bool m_bPrepPerformed{false};

    std::array<bool, MAX_NO_MEAS> m_bFeedbackPerformedOnRepetNo;

    bool m_bNewFeedbackOccurred{false};

    long m_lLastRepetWithFeedback{-1};

    long m_lLastRepetCalledWithSync{-1};

    long m_lLastKernelCalledWithSync{-1};

    SEQSemaphore m_sFeedbackSemaphore; // used to synchronise run and receive

    double m_pdCurrFBTransLog[3]; // data from last feedback in logical coordinates (PRS)

    double m_pdCurrFBRotMatLog[3][3];

    long m_lCurrFBRepetNo{-1};

    double m_pdToTPaceRotMatLog[3][3]; // total PACE transformation in logical coordinates

    double m_pdToTPaceTransLog[3];

    std::array<sSLICE_POS, K_NO_SLI_MAX> m_asOrigSLC; // original slice pos information on measurement start

    sSYNC_WAKEUP m_sWakeUp;
};

inline SEQSemaphore* PaceFeedback::GetFeedbackSemaphore()
{
    return &m_sFeedbackSemaphore;
}

inline void PaceFeedback::MatMulVec3D(double pdResultVec[3], double pdMatrix[3][3], double pdVector[3])
{
    pdResultVec[0] = pdMatrix[0][0] * pdVector[0] + pdMatrix[0][1] * pdVector[1] + pdMatrix[0][2] * pdVector[2];
    pdResultVec[1] = pdMatrix[1][0] * pdVector[0] + pdMatrix[1][1] * pdVector[1] + pdMatrix[1][2] * pdVector[2];
    pdResultVec[2] = pdMatrix[2][0] * pdVector[0] + pdMatrix[2][1] * pdVector[1] + pdMatrix[2][2] * pdVector[2];
}

inline void PaceFeedback::MatMulMat3D(double pdResult[3][3], double pdIn1[3][3], double pdIn2[3][3])
{
    pdResult[0][0] = pdIn1[0][0] * pdIn2[0][0] + pdIn1[0][1] * pdIn2[1][0] + pdIn1[0][2] * pdIn2[2][0];
    pdResult[0][1] = pdIn1[0][0] * pdIn2[0][1] + pdIn1[0][1] * pdIn2[1][1] + pdIn1[0][2] * pdIn2[2][1];
    pdResult[0][2] = pdIn1[0][0] * pdIn2[0][2] + pdIn1[0][1] * pdIn2[1][2] + pdIn1[0][2] * pdIn2[2][2];

    pdResult[1][0] = pdIn1[1][0] * pdIn2[0][0] + pdIn1[1][1] * pdIn2[1][0] + pdIn1[1][2] * pdIn2[2][0];
    pdResult[1][1] = pdIn1[1][0] * pdIn2[0][1] + pdIn1[1][1] * pdIn2[1][1] + pdIn1[1][2] * pdIn2[2][1];
    pdResult[1][2] = pdIn1[1][0] * pdIn2[0][2] + pdIn1[1][1] * pdIn2[1][2] + pdIn1[1][2] * pdIn2[2][2];

    pdResult[2][0] = pdIn1[2][0] * pdIn2[0][0] + pdIn1[2][1] * pdIn2[1][0] + pdIn1[2][2] * pdIn2[2][0];
    pdResult[2][1] = pdIn1[2][0] * pdIn2[0][1] + pdIn1[2][1] * pdIn2[1][1] + pdIn1[2][2] * pdIn2[2][1];
    pdResult[2][2] = pdIn1[2][0] * pdIn2[0][2] + pdIn1[2][1] * pdIn2[1][2] + pdIn1[2][2] * pdIn2[2][2];
}

inline void PaceFeedback::MatCopy3D(double pdOut[3][3], double pdIn[3][3])
{
    pdOut[0][0] = pdIn[0][0];
    pdOut[0][1] = pdIn[0][1];
    pdOut[0][2] = pdIn[0][2];

    pdOut[1][0] = pdIn[1][0];
    pdOut[1][1] = pdIn[1][1];
    pdOut[1][2] = pdIn[1][2];

    pdOut[2][0] = pdIn[2][0];
    pdOut[2][1] = pdIn[2][1];
    pdOut[2][2] = pdIn[2][2];
}

inline void PaceFeedback::TransposeMat3D(double pdResult[3][3], double pdIn[3][3])
{
    pdResult[0][0] = pdIn[0][0];
    pdResult[0][1] = pdIn[1][0];
    pdResult[0][2] = pdIn[2][0];

    pdResult[1][0] = pdIn[0][1];
    pdResult[1][1] = pdIn[1][1];
    pdResult[1][2] = pdIn[2][1];

    pdResult[2][0] = pdIn[0][2];
    pdResult[2][1] = pdIn[1][2];
    pdResult[2][2] = pdIn[2][2];
}

inline void PaceFeedback::VecAddVec3D(double pdResultVec[3], double pdInVec1[3], double pdInVec2[3])
{
    pdResultVec[0] = pdInVec1[0] + pdInVec2[0];
    pdResultVec[1] = pdInVec1[1] + pdInVec2[1];
    pdResultVec[2] = pdInVec1[2] + pdInVec2[2];
}

} // namespace SEQ_NAMESPACE
