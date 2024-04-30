//-----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 2012  All Rights Reserved.  Confidential
//-----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//
//    File: \n4\pkg\MrServers\MrImaging\libSeqPTX\SBB2DExc.h
//
//  Author: pfeujodj  
//
//    Lang: C++
//
// Descrip: SBB base class for 2DRF Excitation
//
//-----------------------------------------------------------------------------


#ifndef SBB2DExc_h
#define SBB2DExc_h 1


#include "MrImaging/libSeqPTX/SBB2DPtx.h"     // base class
#include "MrImaging/libSeqPTX/RF2DArbGenerator.h"


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


class  __IMP_EXP  SBB2DExc : public SBB2DPtx
{
public:

	SBB2DExc(SBBList* pSBBList = nullptr, const char* ptIdent = "Exc2D");
	virtual ~SBB2DExc() = default;

	SBB2DExc(const SBB2DExc &right) = delete;
	SBB2DExc & operator=(const SBB2DExc &right) = delete;

	SBB2DExc(SBB2DExc&& right) = delete;
	SBB2DExc & operator=(SBB2DExc&& right) = delete;


	/// Implements prep-functionality required by base class SeqBuildBlockExcitation.
	bool prep_insertEvents(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;

	/// Implements run-functionality required by base class SeqBuildBlockExcitation.
	bool run_insertEvents(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC) override;

	// Implementation of the call to the pTX pulse calculation (sRF_PULSE_PTX events)
	//  - will be called automatically from the sequence method SeqIF::calculatePTX() after initialize() and before final prepare()
	//  - SeqIF::calculatePTX() will be called only once by MesSer - similar as SeqIF::check() method
	bool calculatePTX(MrProt &rMrProt, const SeqLim &rSeqLim) override;

	// only access to generator member: memory is not to be deleted from outside!
	RF2DArbGenerator* getRFArbGenerator();

	bool setRFParams_Epi(const RF2DArbGenerator::RF2DArbParameter & sParam);

	virtual bool isCompatibilityMode() const;

protected:

	// necessary call to allocate m_pRFGenerator with class RF2DPtxGenerator_Internal
	bool allocateRFGenerator() override;

	// allocate and initalize RF class
	bool initializeRFGenerator(MrProt &rMrProt, const SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;

	// prepares RF class and assigns m_pRFRef
	bool prepareAndAssignRF(MrProt &rMrProt, const SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;

private:

	/// excitation RF pulse containing a sRF_PULSE_ARB and a C2DGradients_Internal object
	RF2DArbGenerator m_RF2DArbGenerator{};

};


#endif    // SBB2DExc_h
