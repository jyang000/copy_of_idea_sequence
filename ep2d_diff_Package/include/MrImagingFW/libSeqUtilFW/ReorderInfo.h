//  -----------------------------------------------------------------------------
//    Copyright (C) Siemens Healthcare GmbH 2015  All Rights Reserved.
//  -----------------------------------------------------------------------------
//
//   Project: NUMARIS/4
//      File: \src\MrImagingFW\libSeqUtilFW\ReorderInfo.h
//   Version: \main\40
//    Author: KUEHN
//      Date: 2013-05-17 18:37:52 +02:00
//
//      Lang: C++
//
//   Descrip: MR::Measurement::Sequence::libSeqUtilFW
//
//   Classes:
//
//  -----------------------------------------------------------------------------

/*!
\file ReorderInfo.h
\brief Base class for reordering schemes
*/

#pragma once

#include "MrImagingFW/libSeqUtilFW/IReorderInfo.h"

#include "MrGlobalDefinitions/MrResult.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrAPI/ITrace.h"

#include <algorithm>

// MrProt KSpace Wrapper
#include "MrProtSrv/Domain/MrProtData/MrProt/KSpace/MrKSpace.h"

// SeqLim
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"

// libSeqUtil
#include "MrImagingFW/libSeqUtilFW/libSeqUtilFW.h"
#include "MrImaging/libSeqUtil/libSeqUtilmsg.h"  // for SU_* messages


typedef struct
{
    unsigned short uiCLin;     // line                    number = 0, ...,65535
    unsigned short uiCPar;     // partition               number = 0, ...,65535
    unsigned char  uiCEco;     // echo (or anything else) number = 0, ...,255
    unsigned char  uiPEFT : 1; // flag to indicate last measured line of current partition and echo
    unsigned char  uiPAFT : 1; // flag to indicate last measured partition of current line and echo
    unsigned char  uiPost : 1; // post shared flag
    unsigned char  uiLast : 1; // last shared flag
    unsigned char  uiFree1: 1; // free bit 1
    unsigned char  uiFree2: 1; // free bit 2
} RawLineData;


enum ContourMode { Ascending, Descending };


#define LIN_IN_PAR        (1)
#define PAR_IN_LIN        (2)
#define SQUARE_SPIRAL     (3)
#define LINEAR_ASCENDING  (1) //fSUReorder_LINEAR_ASCENDING   (1)
#define LINEAR_DESCENDING (2) //fSUReorder_LINEAR_DESCENDING  (2)
#define CENTRIC_DOWN      (3) //fSUReorder_CENTRIC_DOWN       (3)
#define CENTRIC_UP        (4) //fSUReorder_CENTRIC_UP         (4)

#define AbortReorderInfo(A)                                                      \
    {                                                                            \
        if (rSeqLim.isContextNormal() || rSeqLim.isContextPrepForMrProtUpdate()) \
        {                                                                        \
            setNLSStatus(MRI_SUT_SU_ERROR, MRAPI_PRETTY_FUNCTION, A);            \
        }                                                                        \
        else                                                                     \
        {                                                                        \
            setNLSStatus(MRI_SUT_SU_ERROR);                                      \
        }                                                                        \
        return false;                                                            \
    }

//-----------------------------------------------------------------------------
// import/export control
//-----------------------------------------------------------------------------
#ifdef BUILD_libSeqUtilFW
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h" // import/export control

//  ======================================================================
//
//  Class : ReorderInfo
//  Author: Thomas Kluge (tk) for Siemens AG Med MRIA-Seq; (09131)/848049
//          (1998-2001)
//
/// \brief This is the base class for all reordering classes.
/// Supports iPAT.
///
/// It contains the important functions
///
/// - prepareCalculation
/// For basic calculations concerning max./min. line/partition-numbers, center
/// line/partition-number etc..
///
/// - reorder3dE,
/// A reordering function supporting all kinds of non-segmented 2d/3d-reordering
/// schemes especially square-spiral reordering (for "care bolus"
/// 3d-measuremenst) and the elliptical-scanning technique.
//
//  ======================================================================
class __IMP_EXP ReorderInfo : public IReorderInfo
{

public:
    ReorderInfo (long lMaxNoOfReorderIndices = -1, long lCalcBufferSize = -1);

    virtual ~ReorderInfo();


    /// Gives the NLS status. Should be called after a member function returned with error.
    NLS_STATUS getNLSStatus ();

    long getNoOfReorderIndices ();

    long getMaxNoOfReorderIndices ();

    virtual long getTotalNoOfReorderIndices ();

    void setReorderIndexOffset (long ReorderIndexOffset);

    void increaseReorderIndexOffset (long ReorderIndexOffsetIncrement);

    //  *****************************************************
    //
    ///  MaxLineNumber = lines to measure -1
    ///  Note: MaxLineNumber does NOT consider iPAT
    ///  (with iPAT, less lines are actually measured)
    //
    //  *****************************************************
    long getMaxLineNumber ();

    long getLinNoCenterZero (long DeltaReorderIndex = 0);

    long getMinLinNoCenterZero ();

    long getMaxLinNoCenterZero ();

    long getLinesToMeasure ();

    long getKSCenterLin ();

    bool isKSCenterLin (long DeltaReorderIndex = 0);

    //  *****************************************************
    //
    ///  MaxPartitionNumber = partitions to measure -1
    ///  Note: MaxPartitionNumber does NOT consider iPAT (applied in partition direction)
    ///  (with iPATin partition direction, less partitions are actually measured)
    //
    //  *****************************************************
    long getMaxPartitionNumber ();

    virtual long getParNoCenterZero (long DeltaReorderIndex = 0);

    virtual long getMinParNoCenterZero ();

    virtual long getMaxParNoCenterZero ();

    long getPartitionsToMeasure ();

    virtual long getLinNo (long DeltaReorderIndex = 0);

    virtual long getParNo (long DeltaReorderIndex = 0);

    virtual long getKSCenterPar ();

    virtual bool isKSCenterPar (long DeltaReorderIndex = 0);

    virtual long getEcoNo (long DeltaReorderIndex = 0);

    virtual bool isPhaseFT (long DeltaReorderIndex = 0);

    bool existsPhaseFTFlag ();

    void deletePhaseFTFlags ();

    virtual bool isPartitionFT (long DeltaReorderIndex = 0);

    bool existsPartitionFTFlag ();

    void deletePartitionFTFlags ();

    bool isPartitionFTBeforePhaseFT ();

    bool isLastAcquisition();
    void setIsLastAcquisition (bool bValue);

    SEQ::SharedDimension getSharedDimension ();

    virtual bool isPostShared (long DeltaReorderIndex = 0);

    virtual bool isLastShared (long DeltaReorderIndex = 0);

    virtual bool isFreeBit1 (long DeltaReorderIndex = 0);
    virtual bool isFreeBit2 (long DeltaReorderIndex = 0);

    virtual bool isFirstScanInSlice (long DeltaReorderIndex = 0);

    virtual bool isLastScanInSlice (long DeltaReorderIndex = 0);

    //  *****************************************************
    //
    ///  Set up some basic parameters and check for limitations, i.e.
    ///
    ///  m_lMaxLineNumber
    ///  m_lKSpaceCenterLineNumber
    ///  m_lMaxPartitionNumber
    ///  m_lKSpaceCenterPartitionNumber
    ///  m_bEllipticalScanning
    ///  m_dEllipseLineRadius
    ///  m_dEllipsePartitionRadius
    ///
    /// If iPAT is selected, MaxLine/PartNumber is adapted in way,
    /// that MaxLine/PartNumber does not coincide with a 'gap'
    /// (assumption for iPAT: kSpaceCenterLine/Part is not a 'gap')
    //
    //  *****************************************************
    virtual bool prepareCalculation (MrProt&, SeqLim &rSeqLim);

    //  *****************************************************
    //
    ///  Provides a generic method to calculate common
    ///  2D- and 3D-reorderings. It supports elliptical scanning,
    ///  iPAT (=PPA) with and without inplace reference lines/
    ///  partitions, several schemes (par_in_lin, lin_in_par,
    ///  square spiral).
    //
    //  *****************************************************
    virtual bool reorder3dE (MrProt &rMrProt, SeqLim &rSeqLim);

    //  *****************************************************
    //
    ///  Set the following reordering mode parameters
    ///  BasicScanningScheme: e.g. LIN_IN_PAR, PAR_IN_LIN, SQUARE_SPIRAL
    ///  InnerLoopDirection : e.g. LINEAR_ASCENDING, LINEAR_DESCENDING,
    ///                            CENTRIC_DOWN, CENTRIC_UP
    ///  OuterLoopDirection : e.g. LINEAR_ASCENDING, ...
    //
    //  *****************************************************
    void setReorderMode (long lBasicScanningScheme, long lInnerLoopDirection = LINEAR_ASCENDING, long lOuterLoopDirection = LINEAR_ASCENDING);

    long getBasicScanningScheme (void) const;

    long getEchoTrainLength ();

    long getLinesPerSegment ();

    long getMeasuredSegments ();

    //  *****************************************************
    ///
    ///  Calculates the number of raw lines that have to
    ///  be acquired to reach the k-space center
    ///
    //  *****************************************************
    long getRawLinesToKSpaceCenter (SeqLim &rSeqLim);

    //  *****************************************************
    ///
    ///  If PhaseCorrectionMode is active, the method
    ///  getLin/Par() will always return KSpaceCenterLin/
    ///  Par
    ///
    //  *****************************************************
    void enablePhaseCorrection ();

    void disablePhaseCorrection ();

    void enableCheckMode ();

    void disableCheckMode ();

    void setAllowEllipticalScanningAndPartialFourier (bool bPhase, bool bSlice);

    void setOmitLowerPartOfKSpaceIfPF (bool bDoItInSliceDirection, bool bDoItInPhaseDirection);

    void usePrivatePartialFourierFactors (double dPhasePFF = -1.0, double dSlicePFF = -1.0);

    void doNotUsePrivatPartialFourierFactors ();

    double getPhasePartialFourierFactor ();

    SEQ::PartialFourierFactor getPhasePartialFourierFactor_x_8 ();

    double getSlicePartialFourierFactor ();

    SEQ::PartialFourierFactor getSlicePartialFourierFactor_x_8 ();

    //  *****************************************************
    //  *
    ///   printData() will dump all reordering parameters
    ///   via TRACEPUT
    //  *
    //  *****************************************************
    virtual void printData (const char* tText);

    void enableContourCheckMode (ContourMode eMode, long lSkipLine = 1);

    void disableContourCheckMode ();

    long getSegmentsBeforeKSpaceCenter ();

    bool isPATActive ();

    long getPATAccelerationFactorPE ();

    long getPATAccelerationFactor3D ();

    //  *****************************************************
    //
    ///  reduced number of lines(partitions) to measure with PAT
    ///  (including inplace reference lines if Ref.Scan Mode 'Inplace' is selected )
    ///  requires MaxLineNumber and KSpaceCenter members are set
    ///  may be used even if PATMode='NONE', assuming that accel.factor=1 in this case -> will give same result as getLinesToMeasure()
    //
    //  *****************************************************
    virtual long getPATLinesToMeasure ();
    virtual long getPATPartitionsToMeasure ();

    long getNoOfPATRefLines ();

    long getNoOfPATRefPartitions ();

    bool isPATInplaceRefScans ();

    //  *****************************************************
    ///
    ///  acquireForPAT() decides whether the given Lin/Par
    ///  shall be actually scanned or skipped for a
    ///  PAT scan. Inplace reference lines/partions are
    ///  are considered
    ///
    //  *****************************************************
    virtual bool acquireForPAT (long lLin, long lPar, long lCenterLin, long lCenterPar);

    bool acquireForPAT (long lLin, long lPar);

    //  *****************************************************
    //
    ///  isPATScan() decides, whether the given Lin/Par
    ///  is a scanned line/part (inplace reference
    ///  lines/part. are NOT considered) or a 'gap'
    ///  for a PAT scan.
    //
    //  *****************************************************
    virtual bool isPATScan (long lLin, long lPar, long lCenterLin, long lCenterPar);

    bool isPATScan (long lLin, long lPar);

    virtual bool isPATScan (long DeltaReorderIndex = 0);

    //  *****************************************************
    //
    ///  isPATRefScan decides, whether the given Lin/Par
    ///  is an inplace reference line/part. for a PAT scan.
    //
    //  *****************************************************
    virtual bool isPATRefScan (long lLin, long lPar, long lCenterLin, long lCenterPar);

    bool isPATRefScan (long lLin, long lPar);

    virtual bool isPATRefScan (long DeltaReorderIndex = 0);

    //  *****************************************************
    ///
    ///  isPATRefScanInGap() decides, whether the given
    ///  Lin/Par is
    ///    an inplace reference line/part.
    ///  AND
    ///    is located on a 'gap' (i.e. would not be
    ///    acquired without inplace reference lines/part)
    ///
    //  *****************************************************
    virtual bool isPATRefScanInGap (long lLin, long lPar, long lCenterLin, long lCenterPar);

    bool isPATRefScanInGap (long lLin, long lPar);

    virtual bool isPATRefScanInGap (long DeltaReorderIndex = 0);

    ///
    ///      x  -  PATLines   (without ref.lines)
    ///      R  -  Ref.lines
    ///      r  -  Ref.lines in gap (ref.scan mode 'inplace')
    ///      .  -  Gap
    ///
    ///
    ///      Example: AccelPE = 2, RefLines = 8, RefScanMode = 'inplace':
    ///
    ///                 R R R R R R R R
    ///       . x . x . x r x r x r x r x . x . x
    ///      +-------------------------------------> PE (AccelPE=2)
    ///       0 1 2 3...        ^CenterLin
    ///                 ^
    ///                 firstRefLinNo (:= CenterLin - RefLines/2)
    ///         ^
    ///         firstLinNo
    ///
    ///      -> total measured lines for iPAT = PATLines + PATRefLinesInGap
    ///

    long getMinPATRefLinNo ();

    long getMaxPATRefLinNo ();

    long getNoOfPATRefLinesInLowerKSpaceInGap ();

    long getNoOfPATRefLinesInUpperKSpaceInGap ();

    long getNoOfPATRefLinesInGap ();

    long getNoOfPATLinesInLowerKSpace ();

    long getNoOfPATLinesInUpperKSpace ();

    //  *****************************************************
    ///
    ///  returns the number of actually scanned lines
    ///  for PAT (*without* any inplace ref.lines!!!)
    ///
    //  *****************************************************
    long getNoOfPATLines ();

    long getMinPATRefParNo ();

    long getMaxPATRefParNo ();

    long getNoOfPATRefPartitionsInLowerKSpaceInGap ();

    long getNoOfPATRefPartitionsInUpperKSpaceInGap ();

    long getNoOfPATRefPartitionsInGap ();

    long getNoOfPATPartitionsInLowerKSpace ();

    long getNoOfPATPartitionsInUpperKSpace ();

    //  *****************************************************
    ///
    ///  returns the number of actually scanned partitions
    ///  for PAT (*without* any inplace ref.partitions!!!)
    ///
    //  *****************************************************
    long getNoOfPATPartitions ();

    long getNoOfPATRefScansInGap ();

    //  *****************************************************
    ///
    ///   Default definition of PATRefAndImaScan
    ///   (as used in Mdh flags for image reconstruction):
    ///   - is a reference scan
    ///   - is on the PAT grid (i.e. not an 'additional'
    ///     scan in a gap)
    ///
    //  *****************************************************
    virtual bool isPATRefAndImaScan (long DeltaReorderIndex = 0);

    long getPhaseEncodingLines ();

    long getPartitions ();

    // Additional Public Declarations

    //  *****************************************************
    //
    ///  Returns the first actually measured line in k-space.
    ///  E.g. for iPAT it is possible that line=0 is a 'gap' and thus is NOT measured;
    ///  this is in consequence of the rule, that the k-space center line is always measured (i.e. may not be a 'gap').
    ///  From this rule also follows  MinPATLineNo = m_lKSpaceCenterLineNumber % m_lPATAccelerationFactorPE.
    ///
    ///  Note:
    ///  The ICE program expects, that the first lines in k-space are scanned,
    ///  otherwise it needs to be informed via Seq.Exports(YAPS) about the first measured line.
    ///  This might happen e.g. for iPAT, if k-space starts with one or more 'gaps' (i.e. not measured lines).
    ///
    ///  PrepareCalculation() is required before calling this method.
    //
    //  *****************************************************
    virtual long getMinPATLinNo ();

    //  *****************************************************
    //
    ///  Returns the first actually measured partition in k-space.
    ///  E.g. for iPAT it is possible that partition=0 is a 'gap' and thus is NOT measured;
    ///  this is in consequence of the rule, that the k-space center partition is always measured (i.e. must not be a 'gap').
    ///  From this rule also follows  MinPATPartNo = m_lKSpaceCenterPartitionNumber % m_lPATAccelerationFactor3D.
    ///
    ///  Note:
    ///  The ICE program expects, that the first partitions in k-space are scanned,
    ///  otherwise it needs to be informed via Seq.Exports(YAPS) about the first measured partition.
    ///  This might happen e.g. for iPAT, if k-space starts with one or more gaps (i.e. not measured partitions).
    ///
    ///  PrepareCalculation() is required before calling this method.
    //
    //  *****************************************************
    virtual long getMinPATParNo ();

    bool getOmitLowerPartOfKSpaceIfSlicePF () const;

    bool getOmitLowerPartOfKSpaceIfPhasePF () const;

    /// For PATaverage, TSENSE etc.: line offset handling
    void setLineOffset ( long offset ) {
        m_lLineOffset = offset;
        // might want to add a check feature: but go through all stored reorder values
        //   does take long
    }
    /// For PATaverage, TSENSE etc.: line offset handling
    long getLineOffset ( ) {
        return m_lLineOffset;
    }
    /// For PATaverage, TSENSE etc.: partition offset handling
    void setPartitionOffset ( long offset ) {
        m_lPartitionOffset = offset;
        // might want to add a check feature: but go through all stored reorder values
        //   does take long
    }
    /// For PATaverage, TSENSE etc.: partition offset handling
    long getPartitionOffset ( ) {
        return m_lPartitionOffset;
    }
    // PATaverage END



    // ------------------------------------------------------------------------
    /// Supported reorder methods
    // ------------------------------------------------------------------------
    // use this:
    enum {
        Method_undefined = 1, Method_3dE = 2
    };



    // ------------------------------------------------------------------------
    /// Specification of the prepareCalculation (...) method
    // ------------------------------------------------------------------------
    void seteReorderMethod ( long lMethod);

    // ------------------------------------------------------------------------
    /// Returns the selected prepareCalculation (...) method
    // ------------------------------------------------------------------------
    long geteReorderMethod (void) const;

  protected:

      bool setNLSStatus (NLS_STATUS NLSStatus);

      bool setNLSStatus (NLS_STATUS NLSStatus, const char* moduleName, const char* ptAdditionalText = NULL);

      virtual void resetBasic ();

      virtual void resetAll ();

      bool setNoOfReorderIndices (long lNumber);

      void setMaxLineNumber (long lNewMaxLine);

      void setMaxPartitionNumber (long lNewMaxPartition);

      void setLinNo (long DeltaReorderIndex, long lLinNo);

      void setParNo (long DeltaReorderIndex, long lParNo);

      void setEcoNo (long DeltaReorderIndex, long lEcoNo);

      void setPhaseFT (long DeltaReorderIndex, bool bValue);

      void setPartitionFT (long DeltaReorderIndex, bool bValue);

      void setPostShared (long DeltaReorderIndex, bool bValue);

      void setLastShared (long DeltaReorderIndex, bool bValue);

      void setFreeBit1 (long DeltaReorderIndex, bool bValue);
      void setFreeBit2 (long DeltaReorderIndex, bool bValue);

      virtual bool setAllFFTFlags ();

      virtual bool addLinPar3dE (long lLin, long lPar);

      virtual long getValidDataIndex (long DeltaReorderIndex);

      void setBuf (long index, long Value);

      long getBuf (long index);

      void resetBuf ();

      long* getCalcBufferAddress ();

      long getCalcBufferSize ();

      virtual bool isPATScanOnRegularGrid (long lLin, long lPar, long lCenterLin, long lCenterPar);

      virtual bool isPatScanOnCaipirinhaGrid (long lLin, long lPar, long lCenterLin, long lCenterPar);

      virtual long getMinPATParNoOnCaipirinhaGrid();

      virtual long getMinPATParNoOnRegularGrid();



      /// For PATaverage, TSENSE etc.: line offset handling
      long m_lLineOffset;
      /// For PATaverage, TSENSE etc.: partition offset handling
      long m_lPartitionOffset;
      // PATaverage END

      long m_lKSpaceCenterLineNumber;

      bool m_bOmitLowerPartOfKSpaceIfPhasePF;

      long m_lKSpaceCenterPartitionNumber;

      bool m_bOmitLowerPartOfKSpaceIfSlicePF;

      long m_lEchoTrainLength;

      long m_lLinesPerSegment;

      long m_lMeasuredSegments;

      long m_lMaxEchoNumber;

      long m_lBasicScanningScheme;

      long m_lInnerLoopDirection;

      long m_lOuterLoopDirection;

      bool m_bEllipticalScanning;

      double m_dEllipseLineRadius;

      double m_dEllipsePartitionRadius;

      bool m_bAllowEllipticalAndSlicePF;

      bool m_bAllowEllipticalAndPhasePF;

      double m_dPrivateSlicePFFactor;

      double m_dPrivatePhasePFFactor;

      SEQ::SharedDimension m_eSharedDimension;

      long m_lReorderIndexOffset;

      //  *****************************************************
      ///
      ///  Number of raw lines that have to be acquired to
      ///  reach the k-space center
      ///
      //  *****************************************************
      long m_lRawLinesToKSpaceCenter;

      bool m_bPhaseCorrectionActive;

      bool m_bCheckModeActive;

      bool m_bIsLastAcquisition;

      bool m_bContourCheckModeActive;

      ContourMode m_eContourMode;

      long m_lContourModeSkipLine;

      long m_lSegmentsBeforeKSpaceCenter;

      bool m_bPATActive;
      
      SEQ::PATSelMode m_ePATMode;

      long m_lPATAccelerationFactorPE;

      long m_lPATAccelerationFactor3D;

      long m_lPATReorderingShift3D;

      long m_lNoOfPATRefLines;

      long m_lNoOfPATRefPartitions;

      bool m_bPATInplaceRefScans;

      long m_lPhaseEncodingLines;

      long m_lPartitions;

      bool m_bIgnorePATInAddLinPar3dE;


      // ------------------------------------------------------------------------
      /// Specification of the prepareCalculation (...) method
      // ------------------------------------------------------------------------
      long m_lReorderMethod;

      //  *****************************************************
      //
      ///  MaxLineNumber = lines to measure -1
      ///  Note: MaxLineNumber does NOT consider iPAT
      ///  (with iPAT, less lines are actually measured)
      //
      //  *****************************************************
      long m_lMaxLineNumber;

      //  *****************************************************
      //
      ///  MaxPartitionNumber = partitions to measure -1
      ///  Note: MaxPartitionNumber does NOT consider iPAT
      ///  (with iPAT, less partitions are actually measured)
      //
      //  *****************************************************
      long m_lMaxPartitionNumber;

      //  *****************************************************
      //
      /// the array m_aAllInfo[MaxNoOfReoIndices] keeps the
      /// 'chronologically' sorted lines, partitions,
      /// an additional counter and some flags (PeFT, PaFT,
      /// post/last sharing flag) in a bitfield
      /// (see definition of RawLineData).
      //
      //  *****************************************************
      RawLineData* m_aAllInfo;

      long* m_aCalcBuffer;

      long m_lNoOfReorderIndices;

      long m_lMaxNoOfReorderIndices;

      long m_lCalcBufferSize;

      NLS_STATUS m_NLSStatus;



  private:
      ReorderInfo(const ReorderInfo &right);

      ReorderInfo & operator=(const ReorderInfo &right);


};


// Class ReorderInfo


inline bool ReorderInfo::setNLSStatus (NLS_STATUS NLSStatus)
{
    m_NLSStatus = NLSStatus;

    if (NLS_SEVERITY(m_NLSStatus) == NLS_SUCCESS) return false;
    else                                          return true;
}

inline bool ReorderInfo::setNLSStatus (NLS_STATUS NLSStatus, const char* _ptModule, const char* ptAdditionalText)
{
    m_NLSStatus = NLSStatus;

    if (NLS_SEVERITY(m_NLSStatus) == NLS_SUCCESS) return false;
    else
    {
        if (ptAdditionalText==NULL)
        {
            SEQ_TRACE_INFO.print("%s: 0x%lx", _ptModule, m_NLSStatus);
        }
        else
        {
            SEQ_TRACE_INFO.print("%s %s: 0x%lx", _ptModule, ptAdditionalText, m_NLSStatus);
        }
        return true;
    }
}

/// Gives the NLS status. Should be called after a member
///  function returned with
/// error.
inline NLS_STATUS ReorderInfo::getNLSStatus ()
{
    return m_NLSStatus;
}

inline bool ReorderInfo::setNoOfReorderIndices (long lNumber)
{
    if (lNumber>m_lMaxNoOfReorderIndices)
    {
        setNLSStatus(MRI_SUT_SU_NOT_ENOUGH_MEMORY, "ReorderInfo::setNoOfReorderIndices", "Reordering Index out of boundaries");
        m_lNoOfReorderIndices = 0;
        return false;
    }
    else
    {
        m_lNoOfReorderIndices = lNumber;
        return true;
    }
}

inline long ReorderInfo::getNoOfReorderIndices ()
{
    return m_lNoOfReorderIndices;
}

inline long ReorderInfo::getMaxNoOfReorderIndices ()
{
    return m_lMaxNoOfReorderIndices;
}

inline long ReorderInfo::getTotalNoOfReorderIndices ()
{
    return m_lNoOfReorderIndices;
}

inline void ReorderInfo::setReorderIndexOffset (long ReorderIndexOffset)
{
    m_lReorderIndexOffset = std::max(0L, std::min(ReorderIndexOffset, m_lNoOfReorderIndices-1L));
}

inline void ReorderInfo::increaseReorderIndexOffset (long ReorderIndexOffsetIncrement)
{
    m_lReorderIndexOffset = std::max(0L, std::min(m_lReorderIndexOffset+ReorderIndexOffsetIncrement, m_lNoOfReorderIndices-1L));
}

inline long ReorderInfo::getMaxLineNumber ()
{
    return m_lMaxLineNumber;
}

inline void ReorderInfo::setMaxLineNumber (long lNewMaxLine)
{
    if (m_bOmitLowerPartOfKSpaceIfPhasePF)
    {
        m_lKSpaceCenterLineNumber += lNewMaxLine-m_lMaxLineNumber;
    }

    m_lMaxLineNumber = lNewMaxLine;
}

inline long ReorderInfo::getLinNoCenterZero (long DeltaReorderIndex)
{
    return getLinNo(DeltaReorderIndex)-m_lKSpaceCenterLineNumber;
}

inline long ReorderInfo::getMinLinNoCenterZero ()
{
    return -m_lKSpaceCenterLineNumber;
}

inline long ReorderInfo::getMaxLinNoCenterZero ()
{
    return m_lMaxLineNumber-m_lKSpaceCenterLineNumber;
}

inline long ReorderInfo::getLinesToMeasure ()
{
    return m_lMaxLineNumber+1;
}

inline long ReorderInfo::getKSCenterLin ()
{
    return m_lKSpaceCenterLineNumber;
}

inline bool ReorderInfo::isKSCenterLin (long DeltaReorderIndex)
{
    if (getLinNo(DeltaReorderIndex)==m_lKSpaceCenterLineNumber) return true;
    else                                                        return false;
}

inline long ReorderInfo::getMaxPartitionNumber ()
{
    return m_lMaxPartitionNumber;
}

inline void ReorderInfo::setMaxPartitionNumber (long lNewMaxPartition)
{
    if (m_bOmitLowerPartOfKSpaceIfSlicePF)
    {
        m_lKSpaceCenterPartitionNumber += lNewMaxPartition-m_lMaxPartitionNumber;
    }

    m_lMaxPartitionNumber = lNewMaxPartition;
}

inline long ReorderInfo::getParNoCenterZero (long DeltaReorderIndex)
{
    return getParNo(DeltaReorderIndex)-m_lKSpaceCenterPartitionNumber;
}

inline long ReorderInfo::getMinParNoCenterZero ()
{
    return -m_lKSpaceCenterPartitionNumber;
}

inline long ReorderInfo::getMaxParNoCenterZero ()
{
    return m_lMaxPartitionNumber-m_lKSpaceCenterPartitionNumber;
}

inline long ReorderInfo::getPartitionsToMeasure ()
{
    return m_lMaxPartitionNumber+1;
}

inline void ReorderInfo::setLinNo (long DeltaReorderIndex, long lLinNo)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiCLin = (unsigned short) std::max(0L, std::min(m_lMaxLineNumber, lLinNo));
}

inline long ReorderInfo::getLinNo (long DeltaReorderIndex)
{
    if (m_bPhaseCorrectionActive)
    {
        return m_lKSpaceCenterLineNumber;
    }
    else if (m_bCheckModeActive)
    {
        if(m_bEllipticalScanning)
        {

            switch (getValidDataIndex(DeltaReorderIndex)%4)
            {
            default:
            case  0: return m_lMaxLineNumber;
            case  1: return                0;
            case  2: return m_lKSpaceCenterLineNumber;
            case  3: return m_lKSpaceCenterLineNumber;
            }
        }
        else
        {
            switch (getValidDataIndex(DeltaReorderIndex)%4)
            {
            default:
            case  0: return m_lMaxLineNumber;
            case  1: return                0;
            case  2: return m_lMaxLineNumber;
            case  3: return                0;
            }
        }
    }
    else if ( m_bContourCheckModeActive )
    {


        if ( m_lMaxLineNumber > m_lMaxPartitionNumber )
        {

            if ( m_eContourMode == Ascending )
            {
                return ( (DeltaReorderIndex * m_lContourModeSkipLine) % getLinesToMeasure() );
            }
            else
            {
                return ( m_lMaxLineNumber - (DeltaReorderIndex * m_lContourModeSkipLine) % getLinesToMeasure() );
            }

        }
        else   // m_lMaxPartitionNumber > m_lMaxLineNumber *
        {

            long lLin, lPar;
            long lActualPar;


            if ( m_eContourMode == Ascending )
            {
                lActualPar = (DeltaReorderIndex * m_lContourModeSkipLine) % getPartitionsToMeasure();
            }
            else
            {
                lActualPar = m_lMaxPartitionNumber - (DeltaReorderIndex * m_lContourModeSkipLine) % getPartitionsToMeasure();
            }


            if ( (DeltaReorderIndex * m_lContourModeSkipLine) % (2 * getPartitionsToMeasure()) <= m_lMaxPartitionNumber )
            {

                long lMinLin = m_lMaxLineNumber;

                for ( long lI=0; lI < m_lNoOfReorderIndices; lI++ )
                {
                    lLin = (long) m_aAllInfo[lI].uiCLin;
                    lPar = (long) m_aAllInfo[lI].uiCPar;

                    if ( (lLin < lMinLin) && (lPar == lActualPar) )  {  lMinLin = lLin;  }
                }

                return (lMinLin);

            }
            else
            {

                long lMaxLin = 0;

                for ( long lI=0; lI < m_lNoOfReorderIndices; lI++ )
                {
                    lLin = (long) m_aAllInfo[lI].uiCLin;
                    lPar = (long) m_aAllInfo[lI].uiCPar;

                    if ( (lLin > lMaxLin) && (lPar == lActualPar) )  {  lMaxLin = lLin;  }
                }

                return (lMaxLin);

            }

        }

    }
    else
    {
        // PATaveraging: add line offset
        return (long) m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiCLin + getLineOffset();
    }
}

inline void ReorderInfo::setParNo (long DeltaReorderIndex, long lParNo)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiCPar = (unsigned short)std::max(0L, std::min(m_lMaxPartitionNumber, lParNo));
}

inline long ReorderInfo::getParNo (long DeltaReorderIndex)
{
    if (m_bPhaseCorrectionActive)
    {
        return m_lKSpaceCenterPartitionNumber;
    }
    else if (m_bCheckModeActive)
    {
        if(m_bEllipticalScanning )
        {
            switch (getValidDataIndex(DeltaReorderIndex)%4)
            {
            default:
            case  0: return m_lKSpaceCenterPartitionNumber;
            case  1: return m_lKSpaceCenterPartitionNumber;
            case  2: return                     0;
            case  3: return m_lMaxPartitionNumber;
            }
        }
        else
        {
            switch (getValidDataIndex(DeltaReorderIndex)%4)
            {
            default:
            case  0: return m_lMaxPartitionNumber;
            case  1: return                     0;
            case  2: return                     0;
            case  3: return m_lMaxPartitionNumber;
            }
        }
    }
    else if ( m_bContourCheckModeActive )
    {

        if ( m_lMaxLineNumber > m_lMaxPartitionNumber )
        {

            long lLin, lPar;
            long lActualLin;

            if ( m_eContourMode == Ascending )
            {
                lActualLin = (DeltaReorderIndex * m_lContourModeSkipLine) % getLinesToMeasure();
            }
            else
            {
                lActualLin = m_lMaxLineNumber - (DeltaReorderIndex * m_lContourModeSkipLine) % getLinesToMeasure();
            }



            if ( (DeltaReorderIndex * m_lContourModeSkipLine) % (2 * getLinesToMeasure()) <= m_lMaxLineNumber )
            {
                long lMinPar = m_lMaxPartitionNumber;

                for ( long lI=0; lI < m_lNoOfReorderIndices; lI++ )
                {
                    lLin = (long) m_aAllInfo[lI].uiCLin;
                    lPar = (long) m_aAllInfo[lI].uiCPar;

                    if ( (lPar < lMinPar) && (lLin == lActualLin) )  {  lMinPar = lPar;  }
                }

                return (lMinPar);
            }
            else
            {

                long lMaxPar = 0;

                for ( long lI=0; lI < m_lNoOfReorderIndices; lI++ )
                {
                    lLin = (long) m_aAllInfo[lI].uiCLin;
                    lPar = (long) m_aAllInfo[lI].uiCPar;

                    if ( (lPar > lMaxPar) && (lLin == lActualLin) )  {  lMaxPar = lPar;  }
                }

                return (lMaxPar);
            }


        }
        else   // m_lMaxPartitionNumber > m_lMaxLineNumber *
        {

            if ( m_eContourMode == Ascending )
            {
                return ( (DeltaReorderIndex * m_lContourModeSkipLine) % getPartitionsToMeasure() );
            }
            else
            {
                return ( m_lMaxPartitionNumber - (DeltaReorderIndex * m_lContourModeSkipLine) % getPartitionsToMeasure() );
            }

        }

    }
    else
    {
        // PATaveraging: add partition offset
        return (long) m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiCPar + getPartitionOffset();
    }
}

inline long ReorderInfo::getKSCenterPar ()
{
    return m_lKSpaceCenterPartitionNumber;
}

inline bool ReorderInfo::isKSCenterPar (long DeltaReorderIndex)
{
    if (getParNo(DeltaReorderIndex)==m_lKSpaceCenterPartitionNumber) return true;
    else                                                             return false;
}

inline void ReorderInfo::setEcoNo (long DeltaReorderIndex, long lEcoNo)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiCEco = (unsigned char) std::max(0L, std::min(m_lMaxEchoNumber, lEcoNo));
}

inline long ReorderInfo::getEcoNo (long DeltaReorderIndex)
{
    return (long) m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiCEco;
}

inline void ReorderInfo::setPhaseFT (long DeltaReorderIndex, bool bValue)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiPEFT = bValue;
}

inline bool ReorderInfo::isPhaseFT (long DeltaReorderIndex)
{
    return (m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiPEFT && m_bIsLastAcquisition && !m_bPhaseCorrectionActive);
}

inline bool ReorderInfo::existsPhaseFTFlag ()
{
    long i = 0;

    while ( ( i < m_lNoOfReorderIndices ) &&  !m_aAllInfo[i].uiPEFT ) 
    {
        i++;
    }

    return !( i == m_lNoOfReorderIndices );
}

inline void ReorderInfo::deletePhaseFTFlags ()
{
    for (long lI=0; lI<m_lNoOfReorderIndices; lI++)
    {
        m_aAllInfo[lI].uiPEFT = 0;
    }
}

inline void ReorderInfo::setPartitionFT (long DeltaReorderIndex, bool bValue)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiPAFT = bValue;
}

inline bool ReorderInfo::isPartitionFT (long DeltaReorderIndex)
{
    return (m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiPAFT && m_bIsLastAcquisition && !m_bPhaseCorrectionActive);
}

inline bool ReorderInfo::existsPartitionFTFlag ()
{
    long i = 0;

    while ( ( i < m_lNoOfReorderIndices ) && !m_aAllInfo[i].uiPAFT ) 
    {
        i++;
    }

    return !( i == m_lNoOfReorderIndices );
}

inline void ReorderInfo::deletePartitionFTFlags ()
{
    for (long lI=0; lI<m_lNoOfReorderIndices; lI++)
    {
        m_aAllInfo[lI].uiPAFT = 0;
    }
}

inline bool ReorderInfo::isPartitionFTBeforePhaseFT ()
{
    long i = 0, a = -1;

    while (i<m_lNoOfReorderIndices && a==-1)
    {
        if (m_aAllInfo[i].uiPAFT) a = 0;

        if (m_aAllInfo[i].uiPEFT) a = 1;

        i++;
    }

    return (!a);
}

inline bool ReorderInfo::isLastAcquisition()
{
    return m_bIsLastAcquisition;
}

inline void ReorderInfo::setIsLastAcquisition (bool bValue)
{
    m_bIsLastAcquisition=bValue;
}

inline SEQ::SharedDimension ReorderInfo::getSharedDimension ()
{
    return m_eSharedDimension;
}

inline void ReorderInfo::setPostShared (long DeltaReorderIndex, bool bValue)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiPost = bValue;
}

inline bool ReorderInfo::isPostShared (long DeltaReorderIndex)
{
    return (m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiPost && !m_bPhaseCorrectionActive);
}

inline void ReorderInfo::setLastShared (long DeltaReorderIndex, bool bValue)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiLast = bValue;
}

inline bool ReorderInfo::isLastShared (long DeltaReorderIndex)
{
    return (m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiLast && !m_bPhaseCorrectionActive);
}

inline bool ReorderInfo::isFirstScanInSlice (long DeltaReorderIndex)
{
    return (getValidDataIndex(DeltaReorderIndex)==0);
}

inline bool ReorderInfo::isLastScanInSlice (long DeltaReorderIndex)
{
    return ((getValidDataIndex(DeltaReorderIndex)==(m_lNoOfReorderIndices - 1)) && !m_bPhaseCorrectionActive);
}

inline void ReorderInfo::setFreeBit1 (long DeltaReorderIndex, bool bValue)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiFree1 = bValue;
}

inline bool ReorderInfo::isFreeBit1 (long DeltaReorderIndex)
{
    return (m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiFree1 && !m_bPhaseCorrectionActive);
}

inline void ReorderInfo::setFreeBit2 (long DeltaReorderIndex, bool bValue)
{
    m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiFree2 = bValue;
}

inline bool ReorderInfo::isFreeBit2 (long DeltaReorderIndex)
{
    return (m_aAllInfo[getValidDataIndex(DeltaReorderIndex)].uiFree2 && !m_bPhaseCorrectionActive);
}

inline void ReorderInfo::setReorderMode (long lBasicScanningScheme, long lInnerLoopDirection, long lOuterLoopDirection)
{
    m_lBasicScanningScheme = lBasicScanningScheme;
    m_lInnerLoopDirection  = lInnerLoopDirection ;
    m_lOuterLoopDirection  = lOuterLoopDirection ;
}

inline long ReorderInfo::getBasicScanningScheme (void) const
{
    return ( m_lBasicScanningScheme );
}

inline long ReorderInfo::getEchoTrainLength ()
{
    return m_lEchoTrainLength;
}

inline long ReorderInfo::getLinesPerSegment ()
{
    return m_lLinesPerSegment;
}

inline long ReorderInfo::getMeasuredSegments ()
{
    return m_lMeasuredSegments;
}

inline void ReorderInfo::enablePhaseCorrection ()
{
    m_bPhaseCorrectionActive = true;
}

inline void ReorderInfo::disablePhaseCorrection ()
{
    m_bPhaseCorrectionActive = false;
}

inline void ReorderInfo::enableCheckMode ()
{
    m_bCheckModeActive=true;
}

inline void ReorderInfo::disableCheckMode ()
{
    m_bCheckModeActive=false;
}

inline void ReorderInfo::setAllowEllipticalScanningAndPartialFourier (bool bPhase, bool bSlice)
{
    m_bAllowEllipticalAndPhasePF = bPhase;
    m_bAllowEllipticalAndSlicePF = bSlice;
}

inline void ReorderInfo::setOmitLowerPartOfKSpaceIfPF (bool bDoItInSliceDirection, bool bDoItInPhaseDirection)
{
    m_bOmitLowerPartOfKSpaceIfSlicePF = bDoItInSliceDirection;
    m_bOmitLowerPartOfKSpaceIfPhasePF = bDoItInPhaseDirection;
}

inline void ReorderInfo::usePrivatePartialFourierFactors (double dPhasePFF, double dSlicePFF)
{
    if (dPhasePFF<0.5 || dPhasePFF>1.0) m_dPrivatePhasePFFactor =-1.0;
    else m_dPrivatePhasePFFactor                                = dPhasePFF;

    if (dSlicePFF<0.5 || dSlicePFF>1.0) m_dPrivateSlicePFFactor =-1.0;
    else m_dPrivateSlicePFFactor                                = dSlicePFF;
}

inline void ReorderInfo::doNotUsePrivatPartialFourierFactors ()
{
    m_dPrivatePhasePFFactor =-1.0;
    m_dPrivateSlicePFFactor =-1.0;
}

inline double ReorderInfo::getPhasePartialFourierFactor ()
{
    return (double)getLinesToMeasure()/(double)m_lPhaseEncodingLines;
}

inline SEQ::PartialFourierFactor ReorderInfo::getPhasePartialFourierFactor_x_8 ()
{
    double dA = getPhasePartialFourierFactor()/0.125 + 0.5;
    long lA   = std::min(8L, std::max(4L, ((long)dA)));

    switch (lA)
    {
    default:
    case  8: return SEQ::PF_OFF;
    case  7: return SEQ::PF_7_8;
    case  6: return SEQ::PF_6_8;
    case  5: return SEQ::PF_5_8;
    case  4: return SEQ::PF_HALF;
    }
}

inline double ReorderInfo::getSlicePartialFourierFactor ()
{
    return (double)getPartitionsToMeasure()/(double)m_lPartitions;
}

inline SEQ::PartialFourierFactor ReorderInfo::getSlicePartialFourierFactor_x_8 ()
{
    double dA = getSlicePartialFourierFactor()/0.125 + 0.5;
    long lA   = std::min(8L, std::max(4L, ((long)dA)));

    switch (lA)
    {
    default:
    case  8: return SEQ::PF_OFF;
    case  7: return SEQ::PF_7_8;
    case  6: return SEQ::PF_6_8;
    case  5: return SEQ::PF_5_8;
    case  4: return SEQ::PF_HALF;
    }
}

inline long ReorderInfo::getValidDataIndex (long DeltaReorderIndex)
{
    return std::max(0L, std::min(m_lReorderIndexOffset+DeltaReorderIndex, m_lNoOfReorderIndices-1L));
}

inline void ReorderInfo::setBuf (long index, long Value)
{
    m_aCalcBuffer[std::max(0L, std::min(index, m_lCalcBufferSize-1L))]=Value;
}

inline long ReorderInfo::getBuf (long index)
{
    return m_aCalcBuffer[std::max(0L, std::min(index, m_lCalcBufferSize-1L))];
}

inline void ReorderInfo::resetBuf ()
{
    for (int i=0; i<m_lCalcBufferSize; i++) m_aCalcBuffer[i] = 0;
}

inline long* ReorderInfo::getCalcBufferAddress ()
{
    return m_aCalcBuffer;
}

inline long ReorderInfo::getCalcBufferSize ()
{
    return m_lCalcBufferSize;
}

inline void ReorderInfo::enableContourCheckMode (ContourMode eMode, long lSkipLine)
{
    m_bContourCheckModeActive = true;
    m_eContourMode            = eMode;
    m_lContourModeSkipLine    = lSkipLine;
}

inline void ReorderInfo::disableContourCheckMode ()
{
    m_bContourCheckModeActive = false;
}

inline long ReorderInfo::getSegmentsBeforeKSpaceCenter ()
{
    return ( m_lSegmentsBeforeKSpaceCenter );
}

inline bool ReorderInfo::isPATActive ()
{
    return m_bPATActive;
}

inline long ReorderInfo::getPATAccelerationFactorPE ()
{
    return m_lPATAccelerationFactorPE;
}

inline long ReorderInfo::getPATAccelerationFactor3D ()
{
    return m_lPATAccelerationFactor3D;
}

inline long ReorderInfo::getNoOfPATRefLines ()
{
    return m_lNoOfPATRefLines;
}

inline long ReorderInfo::getNoOfPATRefPartitions ()
{
    return m_lNoOfPATRefPartitions;
}

inline bool ReorderInfo::isPATInplaceRefScans ()
{
    return m_bPATInplaceRefScans;
}

inline bool ReorderInfo::acquireForPAT (long lLin, long lPar)
{
    return acquireForPAT(lLin,lPar,m_lKSpaceCenterLineNumber,m_lKSpaceCenterPartitionNumber);
}

inline bool ReorderInfo::isPATScan (long lLin, long lPar)
{
    return isPATScan(lLin,lPar,m_lKSpaceCenterLineNumber,m_lKSpaceCenterPartitionNumber);
}

inline bool ReorderInfo::isPATScan (long DeltaReorderIndex)
{
    return isPATScan(getLinNo(DeltaReorderIndex),getParNo(DeltaReorderIndex),m_lKSpaceCenterLineNumber,m_lKSpaceCenterPartitionNumber);
}

inline bool ReorderInfo::isPATRefScan (long lLin, long lPar)
{
    return isPATRefScan(lLin,lPar,m_lKSpaceCenterLineNumber,m_lKSpaceCenterPartitionNumber);
}

inline bool ReorderInfo::isPATRefScan (long DeltaReorderIndex)
{
    return isPATRefScan(getLinNo(DeltaReorderIndex),getParNo(DeltaReorderIndex),m_lKSpaceCenterLineNumber,m_lKSpaceCenterPartitionNumber);
}

inline bool ReorderInfo::isPATRefScanInGap (long lLin, long lPar)
{
    return isPATRefScanInGap(lLin,lPar,m_lKSpaceCenterLineNumber,m_lKSpaceCenterPartitionNumber);
}

inline bool ReorderInfo::isPATRefScanInGap (long DeltaReorderIndex)
{
    return isPATRefScanInGap(getLinNo(DeltaReorderIndex),getParNo(DeltaReorderIndex),m_lKSpaceCenterLineNumber,m_lKSpaceCenterPartitionNumber);
}

inline long ReorderInfo::getMinPATRefLinNo ()
{
    return m_lKSpaceCenterLineNumber - m_lNoOfPATRefLines/2;
}

inline long ReorderInfo::getMaxPATRefLinNo ()
{
    return getMinPATRefLinNo() + m_lNoOfPATRefLines - 1;
}

inline long ReorderInfo::getNoOfPATRefLinesInLowerKSpaceInGap ()
{
    long lRefLinesInLowerKSpace  = m_lKSpaceCenterLineNumber-getMinPATRefLinNo();
    long lPATLinesWithinRefLines = 0;

    if (m_bPATInplaceRefScans)  // ref. lines in gaps are only used in Inplace mode
    {
        if( m_lPATAccelerationFactorPE != 0 )  {
            lPATLinesWithinRefLines = lRefLinesInLowerKSpace/m_lPATAccelerationFactorPE;
        } else {
            lPATLinesWithinRefLines = lRefLinesInLowerKSpace;  // Accel.==0 should not happen! just to avoid division by zero
        }
        return lRefLinesInLowerKSpace-lPATLinesWithinRefLines;
    }
    return 0;  // no Inplace mode -> no lines in gaps
}

inline long ReorderInfo::getNoOfPATRefLinesInUpperKSpaceInGap ()
{
    long lRefLinesInUpperKSpace  = getMaxPATRefLinNo()-m_lKSpaceCenterLineNumber;
    long lPATLinesWithinRefLines = 0;

    if (m_bPATInplaceRefScans)  // ref. lines in gaps are only used in Inplace mode
    {
        if( m_lPATAccelerationFactorPE != 0 )  {
            lPATLinesWithinRefLines = lRefLinesInUpperKSpace/m_lPATAccelerationFactorPE;
        } else {
            lPATLinesWithinRefLines = lRefLinesInUpperKSpace;  // Accel.==0 should not happen! just to avoid division by zero
        }
        return lRefLinesInUpperKSpace-lPATLinesWithinRefLines;
    }
    return 0;  // no Inplace mode -> no lines in gaps
}

inline long ReorderInfo::getNoOfPATRefLinesInGap ()
{
    return getNoOfPATRefLinesInLowerKSpaceInGap()+getNoOfPATRefLinesInUpperKSpaceInGap();
}

inline long ReorderInfo::getNoOfPATLinesInLowerKSpace ()
{
    if( m_lPATAccelerationFactorPE != 0 )  {
        return m_lKSpaceCenterLineNumber/m_lPATAccelerationFactorPE;
    } else {
        return m_lKSpaceCenterLineNumber;  // Accel.==0 should not happen! just to avoid division by zero
    }
}

inline long ReorderInfo::getNoOfPATLinesInUpperKSpace ()
{
    if( m_lPATAccelerationFactorPE != 0 )  {
        return (getMaxLineNumber()-m_lKSpaceCenterLineNumber) / m_lPATAccelerationFactorPE;
    } else {
        return (getMaxLineNumber()-m_lKSpaceCenterLineNumber); // Accel.==0 should not happen! just to avoid division by zero
    }
}

inline long ReorderInfo::getNoOfPATLines ()
{
    return getNoOfPATLinesInLowerKSpace() + 1 + getNoOfPATLinesInUpperKSpace();
}

inline long ReorderInfo::getMinPATRefParNo ()
{
    return m_lKSpaceCenterPartitionNumber - m_lNoOfPATRefPartitions/2;
}

inline long ReorderInfo::getMaxPATRefParNo ()
{
    return getMinPATRefParNo() + m_lNoOfPATRefPartitions - 1;
}

inline long ReorderInfo::getNoOfPATRefPartitionsInLowerKSpaceInGap ()
{
    long lRefPartitionsInLowerKSpace  = m_lKSpaceCenterPartitionNumber-getMinPATRefParNo();
    long lPATPartitionsWithinRefLines = 0;

    if (m_bPATInplaceRefScans)  // ref. partitions in gaps are only used in Inplace mode
    {
        if( m_lPATAccelerationFactor3D != 0 )  {
            lPATPartitionsWithinRefLines = lRefPartitionsInLowerKSpace/m_lPATAccelerationFactor3D;
        } else {
            lPATPartitionsWithinRefLines = lRefPartitionsInLowerKSpace;  // Accel.==0 should not happen! just to avoid division by zero
        }
        return lRefPartitionsInLowerKSpace-lPATPartitionsWithinRefLines;
    }
    return 0;  // no Inplace mode -> no partitions in gaps
}

inline long ReorderInfo::getNoOfPATRefPartitionsInUpperKSpaceInGap ()
{
    long lRefPartitionsInUpperKSpace  = getMaxPATRefParNo()-m_lKSpaceCenterPartitionNumber;
    long lPATPartitionsWithinRefLines = 0;

    if (m_bPATInplaceRefScans)  // ref. partitions in gaps are only used in Inplace mode
    {
        if( m_lPATAccelerationFactor3D != 0 )  {
            lPATPartitionsWithinRefLines = lRefPartitionsInUpperKSpace/m_lPATAccelerationFactor3D;
        } else {
            lPATPartitionsWithinRefLines = lRefPartitionsInUpperKSpace;  // Accel.==0 should not happen! just to avoid division by zero
        }
        return lRefPartitionsInUpperKSpace-lPATPartitionsWithinRefLines;
    }
    return 0;  // no Inplace mode -> no partitions in gaps
}

inline long ReorderInfo::getNoOfPATRefPartitionsInGap ()
{
    return getNoOfPATRefPartitionsInLowerKSpaceInGap()+getNoOfPATRefPartitionsInUpperKSpaceInGap();
}

inline long ReorderInfo::getNoOfPATPartitionsInLowerKSpace ()
{
    if( m_lPATAccelerationFactor3D != 0 )  {
        return m_lKSpaceCenterPartitionNumber/m_lPATAccelerationFactor3D;
    } else {
        return m_lKSpaceCenterPartitionNumber;    // Accel.==0 should not happen! just to avoid division by zero
    }
}

inline long ReorderInfo::getNoOfPATPartitionsInUpperKSpace ()
{
    if( m_lPATAccelerationFactor3D != 0 )  {
        return (getMaxPartitionNumber()-m_lKSpaceCenterPartitionNumber) / m_lPATAccelerationFactor3D;
    } else {
        return (getMaxPartitionNumber()-m_lKSpaceCenterPartitionNumber);// Accel.==0 should not happen! just to avoid division by zero
    }
}

inline long ReorderInfo::getNoOfPATPartitions ()
{
    return getNoOfPATPartitionsInLowerKSpace() + 1 + getNoOfPATPartitionsInUpperKSpace();
}

inline long ReorderInfo::getPhaseEncodingLines ()
{
    return m_lPhaseEncodingLines;
}

inline long ReorderInfo::getPartitions ()
{
    return m_lPartitions;
}

inline bool ReorderInfo::getOmitLowerPartOfKSpaceIfSlicePF () const
{
    return (m_bOmitLowerPartOfKSpaceIfSlicePF ? true : false);
}

inline bool ReorderInfo::getOmitLowerPartOfKSpaceIfPhasePF () const
{
    return (m_bOmitLowerPartOfKSpaceIfPhasePF ? true : false);
}



inline void ReorderInfo::seteReorderMethod (long lMethod)
{
    m_lReorderMethod = lMethod;
}



inline long ReorderInfo::geteReorderMethod (void) const
{
    return (m_lReorderMethod);
}
