//----------------------------------------------------------------------------------
// <copyright file="SeqLoopMultiband.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens AG, 2013-2015. All Rights Reserved. Confidential.
//   Copyright (C) Siemens Healthcare GmbH, 2015. All Rights Reserved. Confidential.
// </copyright>
// <description>
//   This file contains SeqLoopMultiBand functionality
//
//   For multi-band imaging:
//   - runOuterLoop handles switching between single-band (preparation) and multi-band scans
//   - Slices indices are re-calculated to reflect the grouping of multi-band slices
// </description>
//----------------------------------------------------------------------------------

#pragma once

#ifndef SeqLoopMultiBand_h
#define SeqLoopMultiBand_h

#include "MrImaging/libSeqUtil/SMSProperties.h"
#include "MrProtSrv/Domain/CoreNative/MrPat.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MDS/MDS.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

//------------------------------------------------------------
// Debug
//------------------------------------------------------------
#ifdef DEBUG_ORIGIN
#undef DEBUG_ORIGIN
#endif
#define DEBUG_ORIGIN DEBUG_SEQLOOP

//------------------------------------------------------------
// before refactoring to a template, SeqLoopMultiband used to
// derive from:
// SeqLoopFastIR (e.g. ep2d_diff), or
// SeqLoopLongTRTrig (non-diff EPI, Resolve)
// IIR (some versions of TSE)
// SeqLoop
//
// Now this class is refactored to be a template to enable instance-specific
// inheritance hierarchy.
//------------------------------------------------------------
#include "MrImaging/libSBB/SEQLoop.h"

namespace SEQ_NAMESPACE
{
//------------------------------------------------------------
// SeqLoopMultiBand class
//------------------------------------------------------------
template <class SeqLoop_BASE_TYPE>
class __IMP_EXP_EXPORT_DECL SeqLoopMultiBand : public SeqLoop_BASE_TYPE
{
  public:
    // constructor
    SeqLoopMultiBand() = default;

    // destructor
    virtual ~SeqLoopMultiBand() = default;

    // overloaded SeqLoop_BASE_TYPE functions
    bool prep(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;

    // Setter for multi-band factor
    virtual bool setMultibandFactor(long lMultibandFactor);

    // Returns TRFillEnd from SeqConcat
    virtual long getlTRFillEndInConcat(long) const;

    // Sets multi-band/single runmode
    virtual bool setRunMode(SliceAccelRFRunMode eRunMode);

    // Gets multi-band/single runmode
    virtual SliceAccelRFRunMode getRunMode() const;

    // Gets inner slice counter
    virtual long getInnerSliceCounter() const;

    // Gets inner slice number
    virtual long getInnerSliceNumber() const;

    inline long getMultibandFactor() const
    {
        return m_lMultibandFactor;
    }

    inline void setUseMultiBandOnlyForIR(bool bVal)
    {
		SeqLoop_BASE_TYPE::SBBIRsel.setUseMultiBandOnly(bVal);
    }

  protected:
    // this function calculates and returns the slice index from on a set of looping parameters specified in the interface
    bool mapLoopCounterToSliceIndex(MrProt& rMrProt) override;

    long mapLoopCounterToSliceIndex(
        MrProt&                    rMrProt,
        const SeqConcat*           pSeqConcat,
        SpecialSliceInterleaveMode eSpecialSliceInterleaveMode,
        long                       lInnerSliceCounter,
        long                       lOuterSliceCounter,
        long                       lConcatenationCounter,
        long                       lConcatenations,
        long                       lSlicesToMeasure,
        long                       lSliceOffset,
        long                       lTotalNumberOfOuterLoops) override;

    // Helper function: Copies SeqConcat depending on runmode
    bool copySeqConcat(SeqConcat* pDestination, SeqConcat* pSource);

    // multi-band factor equals number of RF bands
    long m_lMultibandFactor{1};

    // single or multi-band
    SliceAccelRFRunMode m_eRunMode{MULTI_BAND};

    // Bookkeeping in runOuterSliceLoop
    long m_lStoredSlicesToMeasure{-1};

    // SeqConcat for single-band mode
    SeqConcat m_SeqConcatSingleBand[K_NO_SLI_MAX];

    // SeqConcat for multi-band mode
    SeqConcat m_SeqConcatMultiBand[K_NO_SLI_MAX];
};


} // namespace SEQ_NAMESPACE




namespace SEQ_NAMESPACE
{

	//-----------------------------------------------------------------------
	// function:    prep()
	//
	// description: Overloaded SeqLoop function.
	//              Modify number of total slices to be measured (m_SlicesToMeasure)
	//              based on slice acceleration factor
	//-----------------------------------------------------------------------
	template <class SeqLoop_BASE_TYPE>
	bool SeqLoopMultiBand<SeqLoop_BASE_TYPE>::prep(MrProt& rMrProt, SeqLim& rSeqLim, MrProtocolData::SeqExpo& rSeqExpo)
	{
		// Reset stored number of slices for runOuterSliceLoop
		m_lStoredSlicesToMeasure = -1;

		const_cast<SeqBuildBlockIRsel*>(&this->getSBBIRsel())->setSMSProperties(rMrProt);
		const_cast<SeqBuildBlockDB*>(&this->getSBBDB())->setSMSProperties(rMrProt);

		// Reduce slices to adequately prepare the SeqLoop base class
		const long lPreviousSliceArraySize = SMSProperties::getNSlices(rMrProt);
		rMrProt.getsSliceArray().setlSize((int32_t)SMSProperties::getNReducedSlices(rMrProt));

		// Call base class with reduced slice number
		bool bRet = SeqLoop_BASE_TYPE::prep(rMrProt, rSeqLim, rSeqExpo);

		// Reset slice array size
		rMrProt.getsSliceArray().setlSize((int32_t)lPreviousSliceArraySize);

		return bRet;
	}

	//-----------------------------------------------------------------------
	// function:    mapLoopCounterToSliceIndex ()
	//
	// description: Overloaded SeqLoop function.
	//              Maps loop counter to slice index, adapted to simultaneous multi slice imaging.
	//              Special interleave modes like multi concat currently not supported.
	//----------------------------------------------------------------------
	template <class SeqLoop_BASE_TYPE>
	bool SeqLoopMultiBand<SeqLoop_BASE_TYPE>::mapLoopCounterToSliceIndex(MrProt& rMrProt)
	{
		SeqLoop_BASE_TYPE::m_lSliceIndex = mapLoopCounterToSliceIndex(
			rMrProt,
			SeqLoop_BASE_TYPE::m_SeqConcat,
			SeqLoop_BASE_TYPE::m_eSpecialSliceInterleaveMode,
			SeqLoop_BASE_TYPE::m_lInnerSliceCounter,
			SeqLoop_BASE_TYPE::m_lOuterSliceCounter,
			SeqLoop_BASE_TYPE::m_lConcatenationCounter,
			SeqLoop_BASE_TYPE::m_lConcatenations,
			SeqLoop_BASE_TYPE::m_SlicesToMeasure,
			SeqLoop_BASE_TYPE::m_lSliceOffset,
			SeqLoop_BASE_TYPE::m_lTotalNumberOfOuterLoops);
		return true;
	}

	//-----------------------------------------------------------------------
	// function:    mapLoopCounterToSliceIndex ()
	//
	// description: Overloaded SeqLoop function.
	//              Maps loop counter to slice index, adapted to simultaneous multi slice imaging.
	//              Special interleave modes like multi concat currently not supported.
	//-----------------------------------------------------------------------
	template <class SeqLoop_BASE_TYPE>
	long SeqLoopMultiBand<SeqLoop_BASE_TYPE>::mapLoopCounterToSliceIndex(
		MrProt&                    rMrProt,
		const SeqConcat*           pSeqConcat,
		SpecialSliceInterleaveMode eSpecialSliceInterleaveMode,
		long                       lInnerSliceCounter,
		long                       lOuterSliceCounter,
		long                       lConcatenationCounter,
		long                       lConcatenations,
		long                       lSlicesToMeasure,
		long                       lSliceOffset,
		long                       lTotalNumberOfOuterLoops)
	{
		//-----------------------------------------------------------------------
		// Call base SeqLoop mapLoopCounterToSliceIndex function
		//-----------------------------------------------------------------------
		long lSliceIndex = SeqLoop_BASE_TYPE::mapLoopCounterToSliceIndex(
			rMrProt, pSeqConcat, eSpecialSliceInterleaveMode, lInnerSliceCounter, lOuterSliceCounter, lConcatenationCounter, lConcatenations, lSlicesToMeasure, lSliceOffset, lTotalNumberOfOuterLoops);

		//-----------------------------------------------------------------------
		// All single-band situations
		//-----------------------------------------------------------------------
		if (m_eRunMode == SINGLE_BAND || !SMSProperties::isSMS(rMrProt))
		{
			// Not multi-band, use standard routine
			return lSliceIndex;
		}

		//-----------------------------------------------------------------------
		// Rest of the code is for multi-band, i.e. m_lMultibandFactor > 1
		//-----------------------------------------------------------------------

		//-----------------------------------------------------------------------
		// Special TSE case: SMS, interleaved in BH and one slice band only
		//-----------------------------------------------------------------------
		if (m_lMultibandFactor > 1 && eSpecialSliceInterleaveMode == SpecialInterleaveMode_INTERLEAVED_IN_BH && lSlicesToMeasure == 1)
		{
			lSliceIndex = m_lMultibandFactor - 1;

			return lSliceIndex;
		}

		//-----------------------------------------------------------------------
		// Only one slice group
		//-----------------------------------------------------------------------
		if (lSlicesToMeasure == 1)
		{
			// For multi-band scan
			if (m_lMultibandFactor % 2)
			{
				// m_lMultibandFactor is odd
				if (rMrProt.sliceSeries().mode() == SEQ::DESCENDING) // DESCENDING
				{
					// For descending slice order, the "first" slice is the highest slice in location
					lSliceIndex = m_lMultibandFactor - 1;
				}
				else
				{
					lSliceIndex = 0;
				}
			}
			else
			{
				// m_lMultibandFactor is even
				if (rMrProt.sliceSeries().mode() == SEQ::INTERLEAVED) // INTERLEAVED
				{
					lSliceIndex = m_lMultibandFactor / 2;
				}
				else
				{
					if (rMrProt.sliceSeries().mode() == SEQ::DESCENDING) // DESCENDING
					{
						// For descending slice order, the "first" slice is the highest slice in location
						lSliceIndex = m_lMultibandFactor - 1;
					}
					else
					{
						lSliceIndex = 0;
					}
				}
			}

			return lSliceIndex;
		}

		//-----------------------------------------------------------------------
		// Slice mapping for SMS mode
		//-----------------------------------------------------------------------
		lSliceIndex = SMSProperties::mapLoopCounterToSliceIndex(rMrProt, lConcatenations, lSlicesToMeasure, lSliceIndex, m_lMultibandFactor);

		return lSliceIndex;
	}

	template <class SeqLoop_BASE_TYPE>
	bool SeqLoopMultiBand<SeqLoop_BASE_TYPE>::setRunMode(SliceAccelRFRunMode eRunMode)
	{
		m_eRunMode = eRunMode;
		switch (m_eRunMode)
		{
		case SINGLE_BAND:
			copySeqConcat(SeqLoop_BASE_TYPE::m_SeqConcat, m_SeqConcatSingleBand);
			break;
		case MULTI_BAND:
			copySeqConcat(SeqLoop_BASE_TYPE::m_SeqConcat, m_SeqConcatMultiBand);
			break;
		default:
			// Throw error if run mode is unknown
			SEQ_TRACE_ERROR.print("ERROR: RunMode not supported.");
			return false;
		}

		return true;
	}

	template <class SeqLoop_BASE_TYPE>
	long SeqLoopMultiBand<SeqLoop_BASE_TYPE>::getlTRFillEndInConcat(long lIndex) const
	{
		return SeqLoop_BASE_TYPE::m_SeqConcat[lIndex].m_TRFillEnd;
	}

	template <class SeqLoop_BASE_TYPE>
	bool SeqLoopMultiBand<SeqLoop_BASE_TYPE>::copySeqConcat(SeqConcat* pDestination, SeqConcat* pSource)
	{
		if (pDestination != nullptr && pSource != nullptr)
		{
			memcpy((void*)pDestination, pSource, sizeof(SeqConcat) * K_NO_SLI_MAX);
			return true;
		}

		// Throw error
		SEQ_TRACE_ERROR.print("ERROR: copySeqConcat failed.");
		return false;
	}

	template <class SeqLoop_BASE_TYPE>
	SliceAccelRFRunMode SeqLoopMultiBand<SeqLoop_BASE_TYPE>::getRunMode() const
	{
		return m_eRunMode;
	}

	template <class SeqLoop_BASE_TYPE>
	long SeqLoopMultiBand<SeqLoop_BASE_TYPE>::getInnerSliceCounter() const
	{
		return SeqLoop_BASE_TYPE::m_lInnerSliceCounter;
	}

	template <class SeqLoop_BASE_TYPE>
	long SeqLoopMultiBand<SeqLoop_BASE_TYPE>::getInnerSliceNumber() const
	{
		return SeqLoop_BASE_TYPE::m_lInnerSliceNumber;
	}

	template <class SeqLoop_BASE_TYPE>
	bool SeqLoopMultiBand<SeqLoop_BASE_TYPE>::setMultibandFactor(long lMultibandFactor)
	{
		// Check range
		if (lMultibandFactor < 1 || lMultibandFactor > SMSProperties::MAX_MULTIBAND_FACTOR)
		{
			// Everything else out of range
			m_lMultibandFactor = 1;
			SEQ_TRACE_ERROR.print("ERROR: lMultibandFactor = %ld is not in the expected range [1, %ld].", lMultibandFactor, SMSProperties::MAX_MULTIBAND_FACTOR);
			return false;
		}

		m_lMultibandFactor = lMultibandFactor;

		return true;
	}

}; // namespace SEQ_NAMESPACE




#endif
