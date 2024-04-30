//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens Healthcare GmbH 2019  All Rights Reserved. Confidential.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/X
//        File: \src\MrImaging\seq\a_ep2d_diff\DiffusionOrdering.cpp
//     Version: \main\1
//      Author: Kettinger, Adam (SHS DI MR DL AR ACQ)
//        Date: 2019-03-06 14:44:00 +01:00
//
//        Lang: C++
//
//     Descrip: Classes for defining the diffusion ordering in thermal balancing
//
//       Class: DiffusionOrdering
//
//    -----------------------------------------------------------------------------

#include <algorithm>											// std::max_element, std::find_if
#include <numeric>												// std::iota
#include <stdexcept>											// out-of-range exceptions		
#include "MrImaging/seq/a_ep2d_diff/DiffusionOrdering.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "DiffusionOrdering.h"
#include "MrMeasSrv/SeqIF/libRT/sROT_MATRIX.h"

namespace SEQ_NAMESPACE
{

    // default constructor
    DiffusionOrdering::DiffusionOrdering()
    {
    }


    // parametrized constructor to set member variables:
    DiffusionOrdering::DiffusionOrdering(
        SEQ::DiffusionMode	eDiffusionMode,
        std::vector<double> vdBValues,
        std::vector<long>	vlLocalAverages,
        long				lDirections,
        bool				bThermalBalancingReorder
    ) :
        m_eDiffusionMode(eDiffusionMode),
        m_vdBValues(vdBValues),
        m_vlLocalAverages(vlLocalAverages),
        m_lDirections(lDirections),
        m_bThermalBalancingReorder(bThermalBalancingReorder),
        m_pvlUsedDirectionIndices(nullptr),
        m_bOriginalOrderPrepared(false),
        m_bTBOrderPrepared(false),
        m_bFullyPrepared(false),
        m_bParametersUnchanged(false),
        m_bParamsForDirectionReorderUnchanged(false)
    {
    }


    // destructor does nothing
    DiffusionOrdering::~DiffusionOrdering()
    {
    }


    // set functions for member variables:
    void DiffusionOrdering::setDiffusionMode(SEQ::DiffusionMode eDiffusionMode)
    {
        if (m_eDiffusionMode != eDiffusionMode)
        {
            m_bParametersUnchanged = false;
            m_bParamsForDirectionReorderUnchanged = false;
        }

        m_eDiffusionMode = eDiffusionMode;
    }


    void DiffusionOrdering::setBValues(const std::vector<double>& vdBValues)
    {
        if (m_vdBValues != vdBValues)
            m_bParametersUnchanged = false;

        m_vdBValues = vdBValues;
    }


    void DiffusionOrdering::setLocalAverages(const std::vector<long>& vlLocalAverages)
    {
        if (m_vlLocalAverages != vlLocalAverages)
            m_bParametersUnchanged = false;

        m_vlLocalAverages = vlLocalAverages;
    }


    void DiffusionOrdering::setDirections(long lDirections)
    {
        if (m_lDirections != lDirections)
        {
            m_bParametersUnchanged = false;
            m_bParamsForDirectionReorderUnchanged = false;
        }

        m_lDirections = lDirections;
    }


    void DiffusionOrdering::setThermalBalancingReorder(bool bThermalBalancingReorder)
    {
        if (m_bThermalBalancingReorder != bThermalBalancingReorder)
        {
            m_bParametersUnchanged = false;
            m_bParamsForDirectionReorderUnchanged = false;
        }

        m_bThermalBalancingReorder = bThermalBalancingReorder;
    }


    // ==========================
    // checking if index is valid
    bool DiffusionOrdering::checkIndex(long lDiffLoopCounter)
    {
        if (!m_bFullyPrepared)
        {
            SEQ_TRACE_ERROR << "DiffusionOrdering::checkIndex: ERROR - DiffusionReorderInfo is not prepared!";
            return false;
        }
        if (lDiffLoopCounter < 0 || lDiffLoopCounter >= static_cast<long>(m_vsDiffusionLoopReordering.size()))
        {
            SEQ_TRACE_ERROR << "DiffusionOrdering::checkIndex: ERROR - lDiffLoopCounter " << lDiffLoopCounter << " exceeds reorder table size " << static_cast<long>(m_vsDiffusionLoopReordering.size());
            return false;
        }

        return true;
    }


    // ==============================================
    // Get average index from the single loop counter
    long DiffusionOrdering::getAverageIndex(long lDiffLoopCounter)
    {
        if (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
            return 0;

        if (!checkIndex(lDiffLoopCounter))
        {
            SEQ_TRACE_ERROR << "DiffusionOrdering::getAverageIndex: ERROR - Invalid lDiffLoopCounter index " << lDiffLoopCounter;
            throw std::out_of_range("ERROR: invalid lDiffLoopCounter");
        }

        return m_vsDiffusionLoopReordering[lDiffLoopCounter].lAverageIndex;
    }

    
    // ==============================================
    // Get b-value index from the single loop counter
    long DiffusionOrdering::getBValueIndex(long	lDiffLoopCounter)
    {
        if (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
            return 0;

        if (!checkIndex(lDiffLoopCounter))
        {
            SEQ_TRACE_ERROR << "DiffusionOrdering::getBValueIndex: ERROR - Invalid lDiffLoopCounter index " << lDiffLoopCounter;
            throw std::out_of_range("ERROR: invalid lDiffLoopCounter");
        }
        
        return m_vsDiffusionLoopReordering[lDiffLoopCounter].lBValueIndex;
    }


    // ================================================
    // Get direction index from the single loop counter
    long DiffusionOrdering::getDirectionIndex(long lDiffLoopCounter)
    {
        if (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
            return lDiffLoopCounter;

        if (!checkIndex(lDiffLoopCounter))
        {
            SEQ_TRACE_ERROR << "DiffusionOrdering::getDirectionIndex: ERROR - Invalid lDiffLoopCounter index " << lDiffLoopCounter;
            throw std::out_of_range("ERROR: invalid lDiffLoopCounter");;
        }
        
        return m_vsDiffusionLoopReordering[lDiffLoopCounter].lDirectionIndex;
    }
    
    
    // Get thermal balancing reorder flag
    bool DiffusionOrdering::getThermalBalancingReorder()
    {
        return m_bThermalBalancingReorder;
    }
    

    // get fully prepared flag
    bool DiffusionOrdering::isFullyPrepared()
    {
        return m_bFullyPrepared;
    }


    // ============================================================================
    // getting the diffusion counter for the first occurrence of the given b-value
    long DiffusionOrdering::getFirstDiffLoopCounterForBValueIndex(long lBValueIndex)
    {
        // find the first occurrence of the given b-value index:
        auto it = std::find_if(std::begin(m_vsDiffusionLoopReordering), std::end(m_vsDiffusionLoopReordering),
            [lBValueIndex](const DiffusionLoopInfo& sDiffInfo) { return sDiffInfo.lBValueIndex == lBValueIndex; });

        // if we did not find any, something went wrong
        if (it == std::end(m_vsDiffusionLoopReordering))
        {
            SEQ_TRACE_ERROR << "DiffusionOrdering::getFirstDiffLoopCounterForBValueIndex: ERROR - did not find any occurrence of b-value index " << lBValueIndex;
            return 0;
        }
        else
        {
            return static_cast<long>(std::distance(std::begin(m_vsDiffusionLoopReordering), it));
        }
    }


    // ============================================================================
    // getting the diffusion counter for the last occurrence of the given b-value
    long DiffusionOrdering::getLastDiffLoopCounterForBValueIndex(long lBValueIndex)
    {
        // find the first occurrence of the given b-value index:
        auto it = std::find_if(std::rbegin(m_vsDiffusionLoopReordering), std::rend(m_vsDiffusionLoopReordering),
            [lBValueIndex](const DiffusionLoopInfo& sDiffInfo) { return sDiffInfo.lBValueIndex == lBValueIndex; });

        // if we did not find any, something went wrong, give back the last element
        if (it == std::rend(m_vsDiffusionLoopReordering))
        {
            SEQ_TRACE_ERROR << "DiffusionOrdering::getLastDiffLoopCounterForBValueIndex: ERROR - did not find any occurrence of b-value index " << lBValueIndex;
            return (static_cast<long>(m_vsDiffusionLoopReordering.size()) - 1);
        }
        else
        {
            return static_cast<long>(std::distance(it, std::rend(m_vsDiffusionLoopReordering)) - 1);
        }
    }


    // ======================================================
    // getting the diffusion counter for the last occurrence
    // of the highest b-value for kernel check
    long DiffusionOrdering::getDiffLoopCounterForHighestBValue()
    {
        // for q-space, return the last direction index
        if (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
            return (m_lDirections - 1);

        return getLastDiffLoopCounterForBValueIndex(static_cast<long>(m_vdBValues.size()) - 1);
    }


    // ======================================================
    // getting the diffusion counter for the scan with the
    // largest Read-axis component
    long DiffusionOrdering::getDiffLoopCounterForHighestReadComponent(sSLICE_POS* pSlice, DiffusionDirections* pDidi)
    {
        // for q-space, return the last direction index
        if (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
            return (m_lDirections - 1);

        auto rot_matrix = pSlice->getROT_MATRIX();

        // find the index with the biggest readout component for GSWD calculation
        auto it = std::max_element(std::begin(m_vsDiffusionLoopReordering), std::end(m_vsDiffusionLoopReordering),
            [&](const DiffusionLoopInfo& sDiffInfo1, const DiffusionLoopInfo& sDiffInfo2)
            { 
                double  dRcomp1, dRcomp2;

                if (pDidi->getCoordinateSystem() == MrProtocolData::DIFFDIR_CS_PRS)
                {
                    // If current vector set is already specified in PRS, we are already done:
                    dRcomp1 = pDidi->getY(sDiffInfo1.lDirectionIndex);
                    dRcomp2 = pDidi->getY(sDiffInfo2.lDirectionIndex);
                }
                else
                {
                    // if vector is specified in XYZ, we need to transform it to PRS.
                    // similar calculations are implemented in SBBDiffusion_Base::DidiXYZ2PRS() but
                    // 1) we do not see that here, and
                    // 2) we do not need to calculate all components, only read
                    std::vector<double> const Vec1{ pDidi->getX(sDiffInfo1.lDirectionIndex),  pDidi->getY(sDiffInfo1.lDirectionIndex),  pDidi->getZ(sDiffInfo1.lDirectionIndex) };
                    std::vector<double> const Vec2{ pDidi->getX(sDiffInfo2.lDirectionIndex),  pDidi->getY(sDiffInfo2.lDirectionIndex),  pDidi->getZ(sDiffInfo2.lDirectionIndex) };

                    // we need one row of the inverse matrix, i.e. one column of the original rotation matrix (transpose = inverse)
                    std::vector<double> InvRotMatrixRead{ rot_matrix.dMat[0][1] , rot_matrix.dMat[1][1] , rot_matrix.dMat[2][1] };

                    dRcomp1 = std::inner_product(std::begin(InvRotMatrixRead), std::end(InvRotMatrixRead), std::begin(Vec1), 0.0);
                    dRcomp2 = std::inner_product(std::begin(InvRotMatrixRead), std::end(InvRotMatrixRead), std::begin(Vec2), 0.0);
                }

                double const dBValue1 = m_vdBValues[sDiffInfo1.lBValueIndex];
                double const dBValue2 = m_vdBValues[sDiffInfo2.lBValueIndex];

                // larger actual Read-axis component
                return (abs(sqrt(dBValue1) * dRcomp1) < abs(sqrt(dBValue2) * dRcomp2));
            });

        return static_cast<long>(std::distance(std::begin(m_vsDiffusionLoopReordering), it));
    }

    // =================================================
    // getting the pointer to the used direction indices
    std::vector<long>* DiffusionOrdering::getUsedDirectionIndices()
    {
        return m_pvlUsedDirectionIndices;
    }

    // ==================================================
    // getting the total number of diffusion loop indices
    long DiffusionOrdering::getNumberOfDiffLoopIndices()
    {
        return static_cast<long>(m_vsDiffusionLoopReordering.size());
    }


    // ==================================================
    // storing a reordered direction list in the member loop vector
    bool DiffusionOrdering::storeDiffLoopList(std::list<DiffusionLoopInfo> lsDiffLoopInfos)
    {
        // Consistency check
        if (m_vsDiffusionLoopReordering.size() != lsDiffLoopInfos.size())
        {
            // Something went wrong ...
            SEQ_TRACE_ERROR << "DiffusionOrdering::storeDiffLoopList: ERROR - Thermal balancing direction reorder failed!";
            return false;
        }

        // Store reordering information
        std::list<DiffusionLoopInfo>::iterator itDiffLoopList = lsDiffLoopInfos.begin();

        for (size_t iI = 0; iI < lsDiffLoopInfos.size(); ++iI)
        {
            m_vsDiffusionLoopReordering[iI] = *itDiffLoopList;
            ++itDiffLoopList;
        }

        return true;
    }


    // ===============================================================================
    // prepare the original order, i.e. create a container containing all measurements,
    // taking into account b-values, averages, directions:
    bool DiffusionOrdering::prepOriginalOrder()
    {

        m_bOriginalOrderPrepared = false;
        m_bTBOrderPrepared = false;
        m_bFullyPrepared = false;

        //setting up the BValueInfoArray:
        m_vsBValueInfoArray.clear();

        for (long lBValueIndex = 0; lBValueIndex < static_cast<long>(m_vdBValues.size()); lBValueIndex++)
        {
            const long lTempDirections = m_vdBValues[lBValueIndex] <= 0 ? 1 : m_lDirections;
            BValueOccurence NewBValue(lTempDirections, m_vlLocalAverages[lBValueIndex]);
            m_vsBValueInfoArray.push_back(NewBValue);
        }

        // Initialization
        m_vsDiffusionLoopReordering.clear();

        // Find maximum number of averages
        auto maxIt = std::max_element(m_vsBValueInfoArray.begin(), m_vsBValueInfoArray.end(), [](BValueOccurence const& b1, BValueOccurence const& b2)
            {return b1.lNumberOfAverages < b2.lNumberOfAverages; });

        long lMaxNumberOfAverages = maxIt->lNumberOfAverages;

        // Assemble all required combinations: direction loop innermost, then b-values, then averages.
        for (long lAverageIndex = 0; lAverageIndex < lMaxNumberOfAverages; lAverageIndex++)
        {
            for (long lBValueIndex = 0; lBValueIndex < static_cast<long>(m_vsBValueInfoArray.size()); lBValueIndex++)
            {
                if (m_vsBValueInfoArray[lBValueIndex].lNumberOfAverages >= (lAverageIndex + 1))
                {
                    for (long lDirectionIndex = 0; lDirectionIndex < m_vsBValueInfoArray[lBValueIndex].lNumberOfDirections; lDirectionIndex++)
                    {
                        DiffusionLoopInfo sInfo((*m_pvlUsedDirectionIndices)[lDirectionIndex], lBValueIndex, lAverageIndex, false, false);

                        m_vsDiffusionLoopReordering.push_back(sInfo);
                    }
                }
            }
        }

        m_bOriginalOrderPrepared = true;

        if (!m_bThermalBalancingReorder)
        {
            m_bFullyPrepared = true;
        }

        return true;

    }


    // ==========================
    // prepare b-value reordering
    bool DiffusionOrdering::prepTBBValueReorder()
    {

        // if thermal balancing is not required, do nothing else
        if (!m_bThermalBalancingReorder)
        {   
            m_bTBOrderPrepared = false;
            m_bFullyPrepared = true;

            return true;
        }

        m_bTBOrderPrepared = false;
        m_bFullyPrepared = false;

        // first, delete the first occurrence of the first b-value index, because we explicitly want this to be at the very beginning. store it separately for later.
        long lFirstBValueFirstIndex = getFirstDiffLoopCounterForBValueIndex(0);
        DiffusionLoopInfo sFirstBValueMeas = m_vsDiffusionLoopReordering[lFirstBValueFirstIndex];
        m_vsDiffusionLoopReordering.erase(m_vsDiffusionLoopReordering.begin() + lFirstBValueFirstIndex);

        // Find b-value index with maximum abundance (considering averages * directions)
        auto maxIt = std::max_element(m_vsBValueInfoArray.begin(), m_vsBValueInfoArray.end(), [](BValueOccurence const& b1, BValueOccurence const& b2)
            {return ( b1.lNumberOfAverages * b1.lNumberOfDirections ) < ( b2.lNumberOfAverages * b2.lNumberOfDirections ); });

        long lMaxBValueAbundanceIndex = static_cast<long>(std::distance(m_vsBValueInfoArray.begin(), maxIt));

        // Populate list with elements referring to b-value index with maximum abundance
        std::list<DiffusionLoopInfo> lsDiffLoopInfos;
                
        for (size_t iI = 0; iI < m_vsDiffusionLoopReordering.size(); ++iI)
        {
            if (m_vsDiffusionLoopReordering[iI].lBValueIndex == lMaxBValueAbundanceIndex)
            {
                lsDiffLoopInfos.push_back(m_vsDiffusionLoopReordering[iI]);
            }
        }

        // Loop over all b-values
        for (long lBValueIndex = 0; lBValueIndex < static_cast<long>(m_vsBValueInfoArray.size()); ++lBValueIndex)
        {
            // Skip b-value with maximum abundance (already considered above)
            if (lBValueIndex != lMaxBValueAbundanceIndex)
            {
                const long  lBValueAbundance = m_vsBValueInfoArray[lBValueIndex].lNumberOfDirections * m_vsBValueInfoArray[lBValueIndex].lNumberOfAverages;

                if (lBValueAbundance != 0)
                {
                    // Required distance (number of list entries) between elements with current b-value index
                    const float fInsertAfterEach = static_cast<float>(lsDiffLoopInfos.size()) / static_cast<float>(lBValueAbundance);

                    // Set iterator to start of list
                    std::list<DiffusionLoopInfo>::iterator itDiffLoopList = lsDiffLoopInfos.begin();

                    // Initialize distance counter
                    float fDistance = 0.;

                    // before inserting the first found element, move fInsertAfterEach/2 to achieve a more balanced positioning
                    do
                    {
                        ++itDiffLoopList;
                        fDistance += 1.f;
                    } while ((fDistance <= fInsertAfterEach / 2.0f) && (itDiffLoopList != lsDiffLoopInfos.end()));

                    // reset distance counter
                    fDistance = 0.;

                    // Insert elements referring to current b-value with appropriate frequency
                    for (size_t iI = 0; iI < m_vsDiffusionLoopReordering.size(); ++iI)
                    {

                        if (m_vsDiffusionLoopReordering[iI].lBValueIndex == lBValueIndex)
                        {
                            // Insert element before current position
                            lsDiffLoopInfos.insert(itDiffLoopList, m_vsDiffusionLoopReordering[iI]);

                            // Advance iterator to next position (skip appropriate number of elements)
                            do
                            {
                                if (itDiffLoopList != lsDiffLoopInfos.end())
                                {
                                    ++itDiffLoopList;
                                }
                                else
                                {
                                    itDiffLoopList = lsDiffLoopInfos.begin();
                                }
                                fDistance += 1.f;
                            } while ((fDistance < fInsertAfterEach) && (itDiffLoopList != lsDiffLoopInfos.end()));

                            // Initialize distance counter for next element, considering the remainder
                            fDistance -= fInsertAfterEach;
                        }
                    }
                }
            }
        }

        // store reordering information
        storeDiffLoopList(lsDiffLoopInfos);

        // adding the first b-value to the beginning:
        m_vsDiffusionLoopReordering.insert(m_vsDiffusionLoopReordering.begin(), sFirstBValueMeas);

        m_bTBOrderPrepared = true;
        m_bFullyPrepared = true;

        return true;
    }


    // ============================
    // prepare direction reordering
    bool DiffusionOrdering::prepTBDirectionReorder(bool bFastPrep, DiffusionDirections* pDidi)
    {
        // if fast prepare is required, and if everything is prepared, do nothing:
        if (bFastPrep && m_bFullyPrepared && m_bParamsForDirectionReorderUnchanged)
        {
            return true;
        }

        m_bOriginalOrderPrepared = false;
        m_bFullyPrepared = false;
        m_bTBOrderPrepared = false;

        // initializing original direction index set
        m_vlOriginalDirectionIndices.resize(m_lDirections);
        std::iota(m_vlOriginalDirectionIndices.begin(), m_vlOriginalDirectionIndices.end(), 0);

        // without thermal balancing, or with a low number of directions, we are done
        if (!m_bThermalBalancingReorder || m_lDirections < 6)
        {
            m_pvlUsedDirectionIndices = &m_vlOriginalDirectionIndices;
            return true;
        }

        // initializing balanced direction index set
        m_vlBalancedDirectionIndices.clear();

        // initializing auxiliary direction index set:
        std::list<long> llTempDirectionIndices(m_lDirections);
        std::iota(llTempDirectionIndices.begin(), llTempDirectionIndices.end(), 0);

        // initialize balanced vector set
        std::list<VectorStruct> lsTempDirectionOrder;

        // fill the new sets with the original directions:
        for (long lDir = 0; lDir < m_lDirections; lDir++)
        {
            VectorStruct sTemp = { pDidi->getX(lDir), pDidi->getY(lDir), pDidi->getZ(lDir) };
            lsTempDirectionOrder.push_back(sTemp);
        }

        VectorStruct sLastVector = { 0., 0., 0. };
        VectorStruct sLastButOneVector = { 0., 0., 0. };

        // Loop over all vectors
        for (long iJ = 0; iJ < m_lDirections; ++iJ)
        {
            // Iterator pointing to the best match
            std::list<VectorStruct>::iterator itBalancedDirectionOrder;
            std::list<long>::iterator itBalancedDirectionIndices;

            if (iJ == 0)
            {
                // Assign first vector from balanced order
            itBalancedDirectionOrder = lsTempDirectionOrder.begin();
            itBalancedDirectionIndices = llTempDirectionIndices.begin();
            }
            else if (iJ == 1)
            {
                // Find 'most orthogonal' vector in list
                double                            dMinValue = 0.;
                std::list<VectorStruct>::iterator itFind = lsTempDirectionOrder.begin();
                std::list<long>::iterator itIndexFind = llTempDirectionIndices.begin();

                while (itFind != lsTempDirectionOrder.end())
                {
                    // Orthogonality criterion: minimum absolute dot product
                    double dValue = fabs(VectorStruct::DotProduct(sLastVector, *itFind));

                    if ((itFind == lsTempDirectionOrder.begin()) || (dValue < dMinValue))
                    {
                        dMinValue = dValue;
                        itBalancedDirectionOrder = itFind;
                        itBalancedDirectionIndices = itIndexFind;
                    }

                    ++itFind;
                    ++itIndexFind;
                }
            }
            else
            {
                // Calculate cross product of the last two vectors
                VectorStruct sOrthogonalDirection = VectorStruct::CrossProduct(sLastButOneVector, sLastVector);

                // Find 'most parallel' vector in list
                double                            dMaxValue = 0.;
                std::list<VectorStruct>::iterator itFind = lsTempDirectionOrder.begin();
                std::list<long>::iterator itIndexFind = llTempDirectionIndices.begin();

                while (itFind != lsTempDirectionOrder.end())
                {
                    // Most parallel: maximum dot product
                    double dValue = VectorStruct::DotProduct(sOrthogonalDirection, *itFind);

                    if ((itFind == lsTempDirectionOrder.begin()) || (fabs(dValue) > fabs(dMaxValue)))
                    {
                        dMaxValue = dValue;
                        itBalancedDirectionOrder = itFind;
                        itBalancedDirectionIndices = itIndexFind;
                    }

                    ++itFind;
                    ++itIndexFind;
                }
            }

            // store the found diffusion index in array:
            m_vlBalancedDirectionIndices.push_back(*itBalancedDirectionIndices);

            // Remember last but one vector
            sLastButOneVector = sLastVector;

            // Remember last vector
            sLastVector = *itBalancedDirectionOrder;

            // Delete element from temporary storage
            lsTempDirectionOrder.erase(itBalancedDirectionOrder);
            llTempDirectionIndices.erase(itBalancedDirectionIndices);
        }

        // storing balanced indices
        m_pvlUsedDirectionIndices = &m_vlBalancedDirectionIndices;
        
        return true;
    }


    // ==================
    // prepare everything
    bool DiffusionOrdering::prepAll(bool bFastPrep, DiffusionDirections* pDidi)
    {
        // set methods are called before this prepAll, so if any parameter is changed, this will not be true.
        if (bFastPrep && m_bParametersUnchanged && m_bFullyPrepared)
            return true;

        // no need for q-space and free mode
        if (m_eDiffusionMode == SEQ::DIFFMODE_QSPACE)
        {
            m_bFullyPrepared = true;
            return true;
        }

        // we do not want any reordering in free mode
        if (m_eDiffusionMode == SEQ::DIFFMODE_FREE)
        {
            setThermalBalancingReorder(false);
        }

        prepTBDirectionReorder(bFastPrep, pDidi);
        prepOriginalOrder();
        prepTBBValueReorder();

        if (m_bFullyPrepared)
        {
            m_bParametersUnchanged = true; //if the next time this is called in the context PrepForBinarySearch, nothing needs to be done.
            m_bParamsForDirectionReorderUnchanged = true;
        }
        return m_bFullyPrepared;
    }

}   // end namespace SEQ_NAMESPACE

