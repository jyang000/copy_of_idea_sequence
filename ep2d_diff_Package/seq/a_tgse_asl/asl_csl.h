//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2013  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/X
//      File: \src\MrImaging\seq\a_tgse_asl\asl_csl.h
//    Author: pfeujodj
//      Date: 2018-08-07 12:22:22 +02:00
//
//	    Lang: C++
//
//	 Descrip: 3D ASL sequence with CompositeSeqLoop support
//
//	-----------------------------------------------------------------------------
#pragma once

#include "MrImagingFW/libSBBFW/StdSeqIF.h"
#include "MrGlobalDefinitions/MrResult.h"

#include "MrImaging/seq/SystemProperties.h"
#include "MrImagingFW/libSeqSysProp/SysProperties.h"

#include "MrImaging/libSL/StdSL.h"
#include "MrImaging/libSL/StdSL_ID.h"
#include "MrImaging/libSL/StdMediator.h"
#include "MrImagingFW/libCSL/RFSpoiledLeaf.h"

#include "MrImagingFW/libKSpace/SamplingParams.h"
#include "MrImagingFW/libKSpace/KSpaceException.h"
#include "MrImagingFW/libKSpace/StdCheckSampling.h"

#include "MrMeasSrv/SeqIF/Sequence/sequmsg.h"
#include "MrMeasSrv/SeqIF/libRT/sSLICE_POS.h"

#ifdef WIN32
#include "MrProtSrv/Domain/MrProtocol/StdProtRes/StdProtRes.h"
#include "MrMeasSrv/SeqIF/Sequence/Sequence.h"
#include <vector>
#endif

#include "MrImaging/seq/a_tgse_asl/tgse_asl_csl.h"


#ifdef BUILD_SEQU
    #define __OWNER
#endif

// The following include is necessary for the DLL generation
#include "MrGlobalDefinitions/ImpExpCtrl.h"

#include "MrImaging/seq/a_tgse_asl/a_tgse_asl_UI.h"
#include "MrImaging/seq/a_tgse_asl/ASLDefines.h"
   
//  --------------------------------------------------------------------------
//  Forward declarations
//  --------------------------------------------------------------------------
namespace MrProtocolData
{
class MrProtData;
}

class SeqLim;
class SeqExpo;
class Sequence;

namespace SEQ_NAMESPACE
{
        class     Tgse_asl_UI;         // Forward declaration
        typedef   Tgse_asl_UI EpUI;    // Variant specific UI hook
		
#ifdef WIN32
    //  ----------------------------------------------------------------------
    //
    //  Name        :  getUI
    //
    //  Description :
    /// \brief         Returns the pointer to the UI class
    ///
    //
    //  Return      :  EpCommonUI*
    //
    //  ----------------------------------------------------------------------
    EpUI* getUI (MrUILinkBase* const pThis);

    //  ----------------------------------------------------------------------
    //
    //  Name        :  getUI
    //
    //  Description :
    /// \brief         Returns the pointer to the UI class
    ///
    //
    //  Return      :  EpCommonUI*
    //
    //  ----------------------------------------------------------------------
    EpUI* getUI(MrMeasSrv::ISequence* const pSeq);

#endif // WIN32

	/**
	 * @brief The AslCsl sequence class is intended as a demo implementation for teaching purposes. Therefore, all
	 * functions are implemented in a "minimalistic" fashion.
	 */
class __IMP_EXP AslCsl : public StdSeqIF
{
public:
    /**
     * @brief Standard constructor.
     *
     * Creates a AslCsl object with all class member objects being initialized.
     */
    AslCsl(void);

    /**
     * @brief Destructor.
     *
     * The AslCsl sequence is destructed. A potentially associated AslCslUI object will
     * be deleted.
     */
    virtual ~AslCsl();

    // --------------------------------------------------------------------------------------------
    ///  /brief    This method performs several initialization steps.
    ///
    ///            The following information is stored in the provided SeqLim object:
    ///            - General information about this sequence and the used version
    ///            - Hardware requirements
    ///            - Hard limits for UI parameters
    ///
    ///            Furthermore, the FlashUI instance is created.
    ///
    ///  /param    rSeqLim  Object to store sequence information, limits and requirements
    ///
    ///  /return   NLSStatus (SEQU__NORMAL if no error occurred)
    // --------------------------------------------------------------------------------------------
    virtual NLSStatus initialize(SeqLim& rSeqLim);

    // --------------------------------------------------------------------------------------------
    ///  /brief    This method prepares the sequence and its members.
    ///
    ///            In general, this method will check that allowed parameter limits (i.e. soft
    ///            limits) are not exceeded. Depending on the "prepare context", the method will
    ///            perform additional tasks.
    ///
    ///            The "normal context" is the final preparation prior to an actual measurement
    ///            (really executed or simulated). In this case, the method performs all required
    ///            calculations and adjustments.
    ///
    ///            Additionally, there is the "context for binary search". This context is used by
    ///            the user interface (UI) to determine parameter soft limits. The binary search
    ///            is a systematic trial-and-error method which requires a fast feedback method.
    ///            Therefore, the prepare method is performance optimized in the "context for
    ///            binary search".
    ///
    ///            Finally, there is the "context for protocol update". This context is used to set
    ///            set or modify protocol values (while the protocol should never be modified in
    ///            other contexts). This context is required for example for setting WIP parameter
    ///            default values. Another example can be found in sequences which always run at
    ///            minimum TR.
    ///
    ///            The mechanism to obtain the "context for protocol update" is as follows. If the
    ///            prepare method returns SEQU_ERROR, the prepare method will be called once again
    ///            with the same protocol in "context for protocol update" to allow protocol
    ///            modifications.
    ///
    ///  /param    rMrProt   Measurement protocol
    ///  /param    rSeqLim   Limits and requirements for this sequence
    ///  /param    rSeqExpo  Object to store sequence export information
    ///
    ///  /return   NLSStatus (e.g. MRI_SBB_SBB_NORMAL or MRI_SBB_SBB_ERROR)
    // --------------------------------------------------------------------------------------------
    virtual NLSStatus prepare(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    // --------------------------------------------------------------------------------------------
    ///  /brief    This method checks the protocol for gradient stimulation (maximum gradient
    ///            amplitude and GSWD look ahead).
    ///
    ///            This method will be called prior to a measurement to ensure that no gradient
    ///            overflow will occur and stimulation of nerves due to switching gradients does
    ///            not exceed the allowed limits.
    ///
    ///  /param    rMrProt        Measurement protocol
    ///  /param    rSeqLim        Limits and requirements for this sequence
    ///  /param    rSeqExpo       Object to store sequence export information
    ///  /param    pSEQCheckMode  Object to tell which tests need to be performed
    ///
    ///  /return   NLSStatus (SEQU__NORMAL if no error occurred)
    // --------------------------------------------------------------------------------------------
    virtual NLSStatus check(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, SEQCheckMode* pSEQCheckMode);

    // --------------------------------------------------------------------------------------------
    ///  /brief    This method executes the sequence with the given measurement protocol.
    ///
    ///  /param    rMrProt   Measurement protocol
    ///  /param    rSeqLim   Limits and requirements for this sequence
    ///  /param    rSeqExpo  Object to store sequence export information
    ///
    ///  /return   NLSStatus (SEQU__NORMAL if no error occurred)
    // --------------------------------------------------------------------------------------------
    virtual NLSStatus run(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo);

    // --------------------------------------------------------------------------------------------
    ///  /brief    Obsolete method that is required for backward compatibility.
    // --------------------------------------------------------------------------------------------
    virtual NLS_STATUS runKernel(MrProt&, SeqLim&, SeqExpo&, long, long, long, long) { return MRI_SBB_SBB_NORMAL; }

    Tgse_asl_UI* getUI() const   { return m_pUI; };

    virtual NLS_STATUS createUI(SeqLim &rSeqLim);

protected:
	  TGSEAslCsl m_tgseAslCsl;
   
    //  --------------------------------------------------------------
    /// /brief <b> UI class for Templ
    ///
    ///         This class is basically empty on the measurement system
    //  --------------------------------------------------------------
    Tgse_asl_UI* m_pUI;

private:
    AslCsl (const AslCsl &right);
    AslCsl & operator=(const AslCsl &right);
};
}  //SEQ_NAMESPACE

