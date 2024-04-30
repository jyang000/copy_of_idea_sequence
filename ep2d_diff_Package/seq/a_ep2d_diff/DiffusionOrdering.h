//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2019  All Rights Reserved. Confidential.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/X
//        File: \src\MrImaging\seq\a_ep2d_diff\DiffusionOrdering.h
//     Version: \main\1
//      Author: Kettinger, Adam (SHS DI MR DL AR ACQ)
//        Date: 2019-03-06 14:44:00 +01:00
//
//        Lang: C++
//
//     Descrip: Classes for defining the diffusion reordering
//
//       Class: DiffusionOrdering
//
//    -----------------------------------------------------------------------------

// Double include protection:
#pragma once

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------
#include <vector>                                           // std::vector
#include <list>                                             // std::list
#include "MrProtSrv/Domain/CoreNative/ApplDefines.h"		// for diffusion type enum
#include "MrImaging/seq/a_ep2d_diff/didi.h"					// for diffusion directions
#include "MrMeasSrv/SeqIF/libRT/sSLICE_POS.h"

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif
namespace SEQ_NAMESPACE
{

    // ===========================================================================
    /*!
    \struct BValueOccurence

    \brief container to store number of directions and averages for a b-value occurrence
    */
    // ===========================================================================
    struct BValueOccurence
    {
        long lNumberOfDirections;
        long lNumberOfAverages;

        BValueOccurence() :
            lNumberOfDirections(0),
            lNumberOfAverages(0) {};

        BValueOccurence(long lDirections, long lAverages) :
            lNumberOfDirections(lDirections),
            lNumberOfAverages(lAverages) {};

        bool operator==(const BValueOccurence& rhs) const
        {
            const bool bIsIdentical =
                (lNumberOfDirections == rhs.lNumberOfDirections) &&
                (lNumberOfAverages == rhs.lNumberOfAverages);
            return bIsIdentical;
        }
    };


    // ===========================================================================
    /*!
    \struct DiffusionLoopInfo

    \brief container to store information for each diffusion meas
        */
    // ===========================================================================
    struct DiffusionLoopInfo
    {
        long lDirectionIndex;
        long lBValueIndex;
        long lAverageIndex;
        bool bReverseDiffusionDirection;
        bool bReversePhaseEncodingDirection;

        DiffusionLoopInfo() :
            lDirectionIndex(0),
            lBValueIndex(0),
            lAverageIndex(0),
            bReverseDiffusionDirection(false),
            bReversePhaseEncodingDirection(false) {};

        DiffusionLoopInfo(long lDirection, long lBValue, long lAverage, bool bReverseDiffusion, bool bReversePhase) :
            lDirectionIndex(lDirection),
            lBValueIndex(lBValue),
            lAverageIndex(lAverage),
            bReverseDiffusionDirection(bReverseDiffusion),
            bReversePhaseEncodingDirection(bReversePhase) {};

        bool operator==(const DiffusionLoopInfo& rhs) const
        {
            const bool bIsIdentical =
                (lDirectionIndex == rhs.lDirectionIndex) &&
                (lBValueIndex == rhs.lBValueIndex) &&
                (lAverageIndex == rhs.lAverageIndex) &&
                (bReverseDiffusionDirection == rhs.bReverseDiffusionDirection) &&
                (bReversePhaseEncodingDirection == rhs.bReversePhaseEncodingDirection);
            return bIsIdentical;
        }
    };


    // ===========================================================================
    /*!
    \class DiffusionOrdering

    \brief Class for storing and calculating diffusion order information

    */
    // ===========================================================================
    class DiffusionOrdering
    {

    public:

        // default Constructor
        DiffusionOrdering();

        // parametrized constructor to set the member variables:
        DiffusionOrdering(
            SEQ::DiffusionMode	eDiffusionMode,
            std::vector<double> vdBValues,
            std::vector<long>	vlLocalAverages,
            long				lDirections,
            bool				bThermalBalancing
        );

        // Virtual destructor
        ~DiffusionOrdering();

        // set functions for member variables:
        void setDiffusionMode(SEQ::DiffusionMode eDiffusionMode);

        void setBValues(const std::vector<double>& vdBValues);

        void setLocalAverages(const std::vector<long>& vlLocalAverages);

        void setDirections(long lDirections);
        
        void setThermalBalancingReorder(bool bThermalBalancing);
        
        // Get average index from the single loop counter
        long getAverageIndex(long lDiffLoopCounter);

        // Get b-value index from the single loop counter
        long getBValueIndex(long lDiffLoopCounter);

        // Get direction index from the single loop counter
        long getDirectionIndex(long	lDiffLoopCounter);
        
        // Get thermal balancing flag
        bool getThermalBalancingReorder();

        // get pointer to the used direction indices
        std::vector<long>* getUsedDirectionIndices();

        // get number of diffusion indices (number of measurements)
        long getNumberOfDiffLoopIndices();

        // Function for preparing everything:
        bool prepAll(
            bool bFastPrep,
            DiffusionDirections* pDidi
        );

        // virtual function for getting the prepared flag:
        bool isFullyPrepared();

        // getting the diffusion counter for the first occurrence of the given b-value
        long getFirstDiffLoopCounterForBValueIndex(long lBValueIndex);

        // getting the diffusion counter for the last occurrence of the given b-value
        long getLastDiffLoopCounterForBValueIndex(long lBValueIndex);

        // getting the diffusion counter for the first occurrence of the highest b-value.
        // (It will be called in the runKernel)
        long getDiffLoopCounterForHighestBValue();

        // getting the diffusion counter for the scan with the largest read-axis component
        long getDiffLoopCounterForHighestReadComponent(sSLICE_POS* pSlice, DiffusionDirections* pDidi);


    protected:

        // prepare the original order, i.e. create a container containing all measurement, taking into account b-values, averages, directions:
        bool prepOriginalOrder();

        // function for diffusion reordering
        bool prepTBBValueReorder();

        // function for reordering the diffusion directions
        bool prepTBDirectionReorder(
            bool bFastPrep,
            DiffusionDirections* pDidi
        );

    private:

        SEQ::DiffusionMode		m_eDiffusionMode{ SEQ::DIFFMODE_NONE };
        
        std::vector<double>		m_vdBValues;
        
        std::vector<long>		m_vlLocalAverages;
        
        long					m_lDirections{ 0 };
        
        bool					m_bThermalBalancingReorder{ false };

        std::vector<long>		m_vlOriginalDirectionIndices;

        std::vector<long>		m_vlBalancedDirectionIndices;

        std::vector<long>*		m_pvlUsedDirectionIndices{ nullptr };
        
        bool					m_bOriginalOrderPrepared{ false };
        
        bool					m_bTBOrderPrepared{ false };
        
        bool					m_bFullyPrepared{ false };

        bool					m_bParametersUnchanged{ false };

        bool					m_bParamsForDirectionReorderUnchanged{ false };
        
        std::vector<DiffusionLoopInfo> m_vsDiffusionLoopReordering;

        std::vector<BValueOccurence> m_vsBValueInfoArray;

        // storing a reordered list in the member loop vector
        bool storeDiffLoopList(std::list<DiffusionLoopInfo> lsDiffLoopInfos);

        // checking if supplied index is within boundaries
        bool checkIndex(long lDiffLoopCounter);
 
    }; // end of declaration of the interface class DiffusionOrdering



}//end of namespace SEQ_NAMESPACE



