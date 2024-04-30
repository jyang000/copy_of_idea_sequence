//----------------------------------------------------------------------------------
// <copyright file="calcSPAIRTime.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 2006-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#include "MrImaging/seq/common/calcSPAIRTime/calcSPAIRTime.h"

#include "MrImaging/libSeqUtil/SMSProperties.h"
#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"
#include "MrProtSrv/Domain/CoreNative/MrNavigator.h"

calcSPAIRTime::calcSPAIRTime()
{
}

bool calcSPAIRTime::calcSPAIRTimeEPI(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, long lScanTimeBasic, long& lScanTimeSatsEtc, SeqBuildBlockOptfsPrep* pSBBOptfsPrep, SeqBuildBlockOptfs* pSBBOptfs)
{
    // If the calculate button is pressed then the TI value should be calculated automatically
    // For non navigator sequences then the repeat time can be calculated from the TR time, the slices and concatenations
    // For navigator sequences the repeat time must be calc. from the min TR as the TR time has no real time
    double dOptTI         = 1.0;
    double dT1            = 290.0; // Need to query the field strength here and set accordingly
    bool   bNavFlag       = false; // a Nav_flag is needed because the timing will be different
    bool   bExtend_calc   = false;
    long   lIntialTrguess = 1; // this will be the TR value
    double dTrepeat       = 1.0;

    long lSpoilerGradTime = pSBBOptfs->getSpoilerDuration();
    long lRFDuration = pSBBOptfs->getRFDuration();
    long lDwellTime     = 0;

    if (SysProperties::isLowField())
    {
        dT1 = 210.0;
    }
    else if (SysProperties::getNominalBZero() < 2.5)
    {
        dT1 = 230.0; // we are operating at 1.5T therefore we need a reduced T1 value for the fat
    }

    long lMinTR = 0; // lMinTR is the minimum possible TR.  When we have PACE this will be the repeat time
    if (SysProperties::isLowField()) // for Skewed RF pulse, the RF duration is long, the duration of the OPTFS SBB can
                                     // not be ignored.
    {
        lMinTR = ((lScanTimeBasic + lScanTimeSatsEtc + pSBBOptfs->OptfsGetDuration()) / 1000 + 1) * 1000;
    }
    else
    {
        lMinTR = ((lScanTimeBasic + lScanTimeSatsEtc) / 1000 + 1) * 1000;
    }

    const SEQ::RspCompMode currRspCompMode = (rMrProt.getsNavigatorPara().getlRespComp());

    if ((currRspCompMode == SEQ::RESP_COMP_TRIGGER) || (currRspCompMode == SEQ::RESP_COMP_TRIGGER_AND_FOLLOW))
    {
        bNavFlag = true; // and also check to see if we have PACE triggering
    }

    if (bNavFlag == false) // then the TR in the protocol is "real" and can be used for calc the repeat time
    {
        lIntialTrguess = rMrProt.getalTR()[0] / 1000;
        dTrepeat
            = (double)((lIntialTrguess / SMSProperties::getNReducedSlices(rMrProt)) * rMrProt.getsSliceArray().getlConc()); // but for all other variants it is a fn of slices and concats
    }
    else // when we have pace triggering then the protocol TR is meaningless
    {
        dTrepeat = (double)(fSDSRoundUpGRT(lMinTR)) / 1000.0;
    }

    double dNew_repeat = dTrepeat; // we will use this for our first guess of what the values will be

    // Calculate our intial guess for what the TI should be
    // this will be valid if minTR + guessTI > getalTR()[0].
    double dGuess_optTI = 0.;
    double dNewMinTR    = 0.;
    if (SysProperties::isLowField()) // for Skewed RF pulse, the RF duration is long, the duration of the OPTFS SBB can
                                     // not be ignored.
    {
        // Fitted model of SPAIR TI   TI =  A*(1-exp(-B*pow(Trepeat/T1,C) ) )
        // For skewed RF pulse(FA = 240, dur = 69120us, offset = -40Hz) in EPI at low-field:  A = 107.31, B = 0.51,
        // C= 1.22 TI was defined as the time interval from the end of the skewed spair to the start of the excitation
        // pulse.
        const double dTemp = -0.51 * pow((dNew_repeat / dT1), 1.22);
        dGuess_optTI       = 107.31 * (1.0 - exp(dTemp));

        dNewMinTR = ((double)(lMinTR - lSpoilerGradTime - pSBBOptfs->getDwellTime() - lScanTimeSatsEtc)) + (dGuess_optTI * 1000.);
    }
    else
    {
        dGuess_optTI
            = 1.0 * dT1
              * (0.693 - log(1.0 + exp(-(dNew_repeat - 0.0) / dT1))); // this is based on a formulae modelled from the
                                                                      // steady state behaviour and modified by the
                                                                      // slice profile of the pulse used
        dNewMinTR
            = ((double)(lMinTR - lScanTimeSatsEtc))
              + (dGuess_optTI
                 * 1000.); // If a lower field is to be used then the T1 value in this formulae should be modified also.
    }

    if (((dNewMinTR * SMSProperties::getNReducedSlices(rMrProt)) / rMrProt.getsSliceArray().getlConc()) > rMrProt.getalTR()[0])
    {
        bExtend_calc = true; // in order to insert the spir pulse we need to increase TR, which in turn will lengthen the spir pulse ....
    }
    else
    {
        bExtend_calc = false; // we have space to insert the new pulse and all is happy in the world
    }

    //  here we iteratively calculate what our IT will be. This is because if we don not have enough space when we add
    //  our pulse we increase the repeat time Increasing the repeat time in turn increases the necessary dwell time and
    //  so on.

    if ((bExtend_calc == true) || (bNavFlag == true))
    {
        for (long lI = 0; lI < 15; lI++)
        {
            double dLocal_optTI = 0;

            if (SysProperties::isLowField())
            {
                double dTemp = -0.51 * pow((dNew_repeat / dT1), 1.22);
                dLocal_optTI = 107.31 * (1.0 - exp(dTemp));

                dNew_repeat = (double)(lMinTR - lSpoilerGradTime - pSBBOptfs->getDwellTime() - lScanTimeSatsEtc) / 1000. + dLocal_optTI;
            }
            else
            {
                dLocal_optTI
                    = 1.0 * dT1 * (0.693 - log(1.0 + exp(-(dNew_repeat - 0.0) / dT1))); //  calculate the required TI
                dNew_repeat = ((double)lMinTR) / 1000. + dLocal_optTI;
            }
            const double dResidual = dLocal_optTI - dOptTI;
            if (dResidual < 0.001)
            {
                lI = 20;
            }
            dOptTI = dLocal_optTI;
        }
    }
    else
    {
        dOptTI = dGuess_optTI;
    }

    if (bNavFlag == true) // here we tell SeqLoop to set the flip-angle in the Tickle module
    {
        if (!pSBBOptfsPrep->adaptFlipAngle(dTrepeat, rMrProt, rSeqLim, rSeqExpo))
        {
            SEQ_TRACE_ERROR.print("The adapt angle function for tickle pulse has failed");
            return false;
        }
    }

    // For the dwell time calculation, the SliceAdj transit RTEB must not be considered

    if (SysProperties::isLowField())
    {
        // For Skewed SPAIR, TI was defined as the time interval from the end of the skewed spair to the start of the
        // excitation pulse.
        lDwellTime = std::lround(1000.0 * dOptTI) - lSpoilerGradTime - lScanTimeSatsEtc;
    }
    else
    {
        lDwellTime = std::lround(1000.0 * dOptTI) - (lRFDuration + lSpoilerGradTime)
                       - lScanTimeSatsEtc; // we calc the time needed to be added to the SBBOPTFS in order to achieve
                                           // the desired null time
    }

    SEQ_TRACE_DEBUG.print("new_repeat_ms = %lf ms, optTI_ms = %lf ms, dwelltime = %ld us", dNew_repeat, dOptTI, lDwellTime);

    if (lDwellTime < 0)
    {
        lDwellTime = 10; // this needs to accommodate any Reg SATS (or others) between the OPtfs and the imaging module
    }
    lDwellTime = fSDSRoundUpGRT(lDwellTime);
    // For the total time calculation, the SliceAdj transit RTEB has to be considered
    lScanTimeSatsEtc = lScanTimeSatsEtc + pSBBOptfs->OptfsGetDuration(); // and then update the scan time for all the sats

    pSBBOptfs->setDwellTime(lDwellTime);

    return true;
}

