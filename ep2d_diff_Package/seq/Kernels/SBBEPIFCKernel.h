
//    -----------------------------------------------------------------------------
//      Copyright (C) Siemens AG 1998  All Rights Reserved.
//    -----------------------------------------------------------------------------
//
//     Project: NUMARIS/4
//        File: \n4_servers1\pkg\MrServers\MrImaging\seq\Kernels\SBBEPIKernel.h
//     Version: \main\12
//      Author: CC_SEQUENCES CC_FMRI STEMAL8Q LANDWIPD
//        Date: 2015-02-03 11:51:39 +01:00
//
//        Lang: C++
//
//     Descrip: MR::MrServers::MrImaging::seq::Kernels
//
//     Classes:
//
//    -----------------------------------------------------------------------------

#pragma once

#include "MrImaging/seq/Kernels/SBBEPIKernel.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

namespace SEQ_NAMESPACE
{
    class __IMP_EXP SeqBuildBlockEPIFCKernel :
        public SeqBuildBlockEPIKernel
    {
    public:

        // No public default constructor
        SeqBuildBlockEPIFCKernel() = delete;

        SeqBuildBlockEPIFCKernel(
            SBBList* pSBBList
        );

        virtual ~SeqBuildBlockEPIFCKernel() = default;

        SeqBuildBlockEPIFCKernel(const SeqBuildBlockEPIFCKernel& right) = delete;
        SeqBuildBlockEPIFCKernel& operator=(const SeqBuildBlockEPIFCKernel& right) = delete;
        SeqBuildBlockEPIFCKernel(SeqBuildBlockEPIFCKernel&& right)                 = delete;
        SeqBuildBlockEPIFCKernel& operator=(SeqBuildBlockEPIFCKernel&& right) = delete;

        //--------------------------------------------------------------------
        //    Prepares the additional Gradients,e.g. flow compensated gradients
        bool prepAdditionalGradients(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo) override;

        //--------------------------------------------------------------------
        //    adjust TEContributions due to the applied additional Gradients,e.g. flow compensated gradients
        bool adjustTEContributions(MrProt &rMrProt) override;


        //    adjust additional timing config due to the applied additional Gradients,e.g. flow compensated gradients
        bool adjustAdditionalTiming(MrProt &rMrProt) override;

        //    adjust local TR fill due to the applied additional Gradients,e.g. flow compensated gradients
        bool adjustLocalTRFill(MrProt &rMrProt) override;

        // *--------------------------------------------------------------------*
        // * run flow compensated slice refocusing gradient                     *
        // * ------------------------------------------------------------------ * 
        bool runAdditionalGradExcitation(long &lEventBlockEnd) override;

        // *------------------------------------------------------------------- *
        // * Activates or deactivates echo-shifting. If echo  shifting should   *
        // * be activated it is important for the kernel to know the number of	*
        // * the number of (lines res. partitions) counters per segment and the	*
        // * (line res. partition) counter within the center segment with which *
        // * the k-space center is measured.For deactivation only false needs to*
        // * be passed to this function.                                        *
        // * ------------------------------------------------------------------ * 
        virtual void setUseEchoShifting(bool bValue, bool bFlowcomp = false, long lCountersPerSegment = 0, long lCounterInSegmentWithEcho = 0);

        //* --------------------------------------------------------------------*
        //* Prepares the flow compensation gradients in phase encoding direction*
        //* Only the central echo per shot will be fully compensated.	        *
        //* --------------------------------------------------------------------* 
        virtual bool prepFC_PE(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        //* --------------------------------------------------------------------*
        //* Prepares the flow compensation gradients in partition encoding      *
        //* direction Only the central echo per shot will be fully compensated .*
        //* --------------------------------------------------------------------* 
        virtual bool prepFC_PA(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        // * ------------------------------------------------------------------ *
        // * Prepare flow compensated readout gradients                         *
        // * Only the odd echoes in EPI readout will be fully compensated       *
        // * ------------------------------------------------------------------ * 
        virtual bool prepFC_RO(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);

        // * ------------------------------------------------------------------ *
        // * Prepare flow compensated slice refocusing gradient                 *
        // * ------------------------------------------------------------------ * 
        virtual bool prepFC_SS(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);


        //flow compensation configuration functors

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::setbFlowCompRead          *
        // *                                                                    *
        // * Description :  Disables or enables flow compensation in readout    *
        // *                direction                                           *
        // *                                                                    *
        // * Parameter   :  bValue                                              *
        // *                    true  : enables flow compensation               *
        // *                    false : disables flow compensation              *
        // *                                                                    *
        // * Return      :  void                                                *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        void setbFlowCompRead(bool bValue);

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::getbFlowCompRead          *
        // *                                                                    *
        // * Description :  Returns the selection for flow compensation in      *
        // *                readout direction                                   *
        // *                                                                    *
        // * Parameter   :  void                                                *
        // *                                                                    *
        // * Return      :  bool                                                *
        // *                    true  : flow compensation enabled               *
        // *                    false : flow compensation disabled              *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        bool getbFlowCompRead() const;

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::setbFlowCompPhase         *
        // *                                                                    *
        // * Description :  Disables or enables flow compensation in phase      *
        // *                encoding direction                                  *
        // *                                                                    *
        // * Parameter   :  bValue                                              *
        // *                    true  : enables flow compensation               *
        // *                    false : disables flow compensation              *
        // *                                                                    *
        // * Return      :  void                                                *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        void setbFlowCompPhase(bool bValue);

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::getbFlowCompPhase         *
        // *                                                                    *
        // * Description :  Returns the selection for flow compensation in      *
        // *                phase encoding direction                            *
        // *                                                                    *
        // * Parameter   :  void                                                *
        // *                                                                    *
        // * Return      :  bool                                                *
        // *                    true  : flow compensation enabled               *
        // *                    false : flow compensation disabled              *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        bool getbFlowCompPhase() const;

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::setbFlowCompSlice         *
        // *                                                                    *
        // * Description :  Disables or enables flow compensation in slice      *
        // *                selection direction. This function enables /        *
        // *                disables flow compensation of the slice selection   *
        // *                gradient and flow compensation of the partition     *
        // *                encoding simultaneously.                            *
        // *                                                                    *
        // * Parameter   :  bValue                                              *
        // *                    true  : enables flow compensation               *
        // *                    false : disables flow compensation              *
        // *                                                                    *
        // * Return      :  void                                                *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        void setbFlowCompSlice(bool bValue);

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::getbFlowCompSlice         *
        // *                                                                    *
        // * Description :  Returns the status of the flow compensation in      *
        // *                slice selection direction. Both gradients (slice    *
        // *                selection gradient and partiition encoding table)   *
        // *                considered simultaneously.                          *
        // *                                                                    *
        // * Parameter   :  void                                                *
        // *                                                                    *
        // * Return      :  bool                                                *
        // *                    true  : if flow compensation is enabled for     *
        // *                            the slice selection gradient AND for    *
        // *                            the partition encoding table            *
        // *                    false : if flow compensation is disabled for at *
        // *                            least one gradient (slice selection or  *
        // *                            partition encoding table)               *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        bool getbFlowCompSlice() const;

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::setbFlowCompSliceSelection*
        // *                                                                    *
        // * Description :  Disables or enables flow compensation of slice      *
        // *                selection gradient. MrProtocolData::MrFlowData      *
        // *                compensation of the partition encoding table is not *
        // *                effected by this function.                          *
        // *                                                                    *
        // * Parameter   :  bValue                                              *
        // *                    true  : enables flow compensation               *
        // *                    false : disables flow compensation              *
        // *                                                                    *
        // * Return      :  void                                                *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        void setbFlowCompSliceSelection(bool bValue);

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::getbFlowCompSliceSelection*
        // *                                                                    *
        // * Description :  Returns the selection for flow compensation of      *
        // *                the slice selection gradient.                       *
        // *                MrProtocolData::MrFlowData compensation of the      *
        // *                partition encoding table is not effected by this    *
        // *                function.                                           *
        // *                                                                    *
        // * Parameter   :  void                                                *
        // *                                                                    *
        // * Return      :  bool                                                *
        // *                    true  : if flow compensation is enabled for     *
        // *                            the slice selection gradient            *
        // *                    false : if flow compensation is disabled for    *
        // *                            slice selection)                        *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        bool getbFlowCompSliceSelection() const;

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::setbFlowCompPartitionEncode*
        // *                                                                    *
        // * Description :  Disables or enables flow compensation of the        *
        // *                partition encoding table. MrProtocolData::MrFlowData*
        // *                compensation of the slice selection gradient is not *
        // *                effected by this function.                          *
        // *                                                                    *
        // * Parameter   :  bValue                                              *
        // *                    true  : enables flow compensation               *
        // *                    false : disables flow compensation              *
        // *                                                                    *
        // * Return      :  void                                                *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        void setbFlowCompPartitionEncode(bool bValue);

        // * ------------------------------------------------------------------ *
        // *                                                                    *
        // * Name        :  SeqBuildBlockEPIFCKernel::getbFlowCompPartitionEncode*
        // *                                                                    *
        // * Description :  Returns the selection for flow compensation of      *
        // *                the partition encoding table.                       *
        // *                MrProtocolData::MrFlowDatacompensation of the slice *
        // *                selection gradient is not effected by this function.*
        // *                                                                    *
        // * Parameter   :  void                                                *
        // *                                                                    *
        // * Return      :  bool                                                *
        // *                    true  : if flow compensation is enabled         *
        // *                            the partition encoding table            *
        // *                    false : if flow compensation is disabled        *
        // *                            partition encoding table)               *
        // *                                                                    *
        // * ------------------------------------------------------------------ *
        bool getbFlowCompPartitionEncode() const;

    protected:

        //--------------------------------------------------------------------
        //    update flow compensated gradient moment
        bool updateAdditionalGradMoments(long lSectionIndex, MrProt &rMrProt, SeqLim &rSeqLim, long &lFillTimeBeforeGradients, long &lEventBlockEnd) override;

        //--------------------------------------------------------------------
        //    run gradient moment
        bool runFinalGradMoments(long lSectionIndex, MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC, long lFillTimeBeforeGradients, long lEventBlockEnd) override;

        // * --------------------------------------------------------------------
        // *   Calculate the maxwell term from flow compensation gradient in phase direction
        // * ------------------------------------------------------------------ *
        bool calculateMaxwellTerm(MrProt &rMrProt, sSLICE_POS* pSLC);

        // * ------------------------------------------------------------------ *
        // * Perform flow compensation in read direction                        *
        // * ------------------------------------------------------------------ *
        bool m_bFlowCompRead{false};

        // * ------------------------------------------------------------------ *
        // * Perform flow compensation in phase encoding direction              *
        // * ------------------------------------------------------------------ *
        bool m_bFlowCompPhase{false};

        // * ------------------------------------------------------------------ *
        // * Perform flow compensation of the slice selection gradient          *
        // * ------------------------------------------------------------------ *
        bool m_bFlowCompSliceSelect{false};

        // * ------------------------------------------------------------------ *
        // * Perform flow compensation of the partition encoding table          *
        // * ------------------------------------------------------------------ *
        bool m_bFlowCompPartitionEncode{false};

        //  -------------------------------------------------------------------
        /// Gradient raster time: 10us
        //  -------------------------------------------------------------------
        enum { eGradRasterTime = 10 };

        // * ------------------------------------------------------------------ *
        // * Gradient pulse used for slice refocusing with flow compensation in *
        // * slice select direction                                             *
        // * ------------------------------------------------------------------ *
        sGRAD_PULSE_TRAP m_sP_3D_FC{ "GS refoc FC" };

        // * ------------------------------------------------------------------ *
        // * Gradient pulse used for slice refocusing in slice select           *
        // * direction                                                          *
        // * ------------------------------------------------------------------ *
        sGRAD_PULSE_TRAP m_sP_3D{ "GS refoc" };

        // * ------------------------------------------------------------------ *
        // * Gradient pulse used for a bipolar flow compensation in phase       *
        // * encoding direction                                                 *
        // * ------------------------------------------------------------------ *
        sGRAD_PULSE_TRAP m_sPE_FC1{ "PE FC1" };
        sGRAD_PULSE_TRAP m_sPE_FC2{ "PE FC2" };

        // * ------------------------------------------------------------------ *
        // * Gradient pulse used for flow compensation in readout direction. It *
        // * is applied prior to the readout prephasing gradient sP_ROP.        *
        // * sP_ROP_FC is applied with zero duration and amplitude if no flow   *
        // * compensation is selected for readout.                              *
        // * ------------------------------------------------------------------ *
        sGRAD_PULSE_TRAP m_sP_ROP_FC{ "GR prephase FC" };

        // * ------------------------------------------------------------------ *
        // * Gradient pulse used for prephasing in readout direction            *
        // * ------------------------------------------------------------------ *
        sGRAD_PULSE_TRAP m_sP_ROP{ "GR prephase" };

        // * ------------------------------------------------------------------ *
        // * Gradient pulse used for a bipolar flow compensation in partition   *
        // * encoding direction                                                 *
        // * ------------------------------------------------------------------ *
        sGRAD_PULSE_TRAP m_sPA_FC1{ "PA FC1" };
        sGRAD_PULSE_TRAP m_sPA_FC2{ "PA FC2" };

        long   m_lEchoshiftDelay_end{0};
        long   m_lEchoshiftDelay{0};
        double m_dBlipsMoment{0.0};
    };
}


inline void SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::setbFlowCompRead(bool bValue)
{
    m_bFlowCompRead = bValue;
}

inline bool SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::getbFlowCompRead() const
{
    return m_bFlowCompRead;
}

inline void SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::setbFlowCompPhase(bool bValue)
{
    m_bFlowCompPhase = bValue;
}

inline bool SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::getbFlowCompPhase() const
{
    return m_bFlowCompPhase;
}

inline void SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::setbFlowCompSlice(bool bValue)
{
    m_bFlowCompSliceSelect = bValue;
    m_bFlowCompPartitionEncode = bValue;
}

inline bool SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::getbFlowCompSlice() const
{
    return m_bFlowCompSliceSelect && m_bFlowCompPartitionEncode;
}

inline void SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::setbFlowCompSliceSelection(bool bValue)
{
    m_bFlowCompSliceSelect = bValue;
}

inline bool SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::getbFlowCompSliceSelection() const
{
    return m_bFlowCompSliceSelect;
}

inline void SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::setbFlowCompPartitionEncode(bool bValue)
{
    m_bFlowCompPartitionEncode = bValue;
}

inline bool SEQ_NAMESPACE::SeqBuildBlockEPIFCKernel::getbFlowCompPartitionEncode() const
{
    return m_bFlowCompPartitionEncode;
}
