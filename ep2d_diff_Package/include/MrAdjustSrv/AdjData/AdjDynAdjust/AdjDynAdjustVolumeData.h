// -----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 2015  All Rights Reserved.  Confidential
// -----------------------------------------------------------------------------
//
// Project: NUMARIS/X
//    File: 
// Version: 
//  Author: 
//    Date: 
//
//    Lang: C++
//
// Descrip:
//
// Classes:
//
// -----------------------------------------------------------------------------

#ifndef ADJDYNADJUSTVOLUMEDATA_H_INCLUDED
#define ADJDYNADJUSTVOLUMEDATA_H_INCLUDED

// -----------------------------------------------------------------------------
// includes
// -----------------------------------------------------------------------------
#include "MrAdjustSrv/AdjDefines.h"
#include "MrGlobalDefinitions/MrBasicTypes.h"
#include "MrMeasSrv/MeasUtils/Complex.h"
#include "MrMeasSrv/MeasRealTimeTypes.h"


//-----------------------------------------------------------------------------
// C type data container class which corresponds to the dynamic
// data container class MrProtocolData::VectorPatDbl.
//-----------------------------------------------------------------------------
class AdjDynAdjustVectorPatDbl
{
public:
  double m_dSag; // Sagittal
  double m_dCor; // Coronal
  double m_dTra; // Transversal
};

//-----------------------------------------------------------------------------
// C type data container class which corresponds to the dynamic
// data container class MrProtocolData::MrSliceData.
//-----------------------------------------------------------------------------
class AdjDynAdjustSliceData
{
public:
  AdjDynAdjustVectorPatDbl m_sPosition; // Position vector
	AdjDynAdjustVectorPatDbl m_sNormal;   // Normal vector
	double m_dThickness;                  // Thickness
	double m_dPhaseFOV;                   // Phase FoV
	double m_dReadoutFOV;                 // Readout FoV
	double m_dInPlaneRot;                 // Inplane rotation
  MeasSrv::rt_string  m_sDataRole;
  MeasSrv::rt_string  m_sSliceID;

};

//-----------------------------------------------------------------------------
// C type data container class which corresponds to the dynamic
// data container class MrProtocolData::MrDynamicAdjustVolumeData.
//-----------------------------------------------------------------------------
class AdjDynAdjustVolumeData
{
public:
  AdjDynAdjustSliceData m_sCuboid;                   // the cuboid -->
	int32_t m_lFlags;                                  // bitmask providing additional information about the properties of this adjustment subvolume  -->
	double  m_dWeighting;                              // relative weight of this adjustment subvolume  -->
  int32_t m_lFrequencyHz;                            // f0 frequency [Hz]
  bool    m_bFrequencyIsRelative;                    // f0 frequency is relative (additiv)
  double  m_dRFAmplitudeScaleFactor;                 // RF pulse ampltitude scale factor
  bool    m_bRFAmplitudeScaleFactorIsRelative;       // RF pulse ampltitude scale factor is relative (multiplicative)
  double  m_dGradOffsetX_mTm;                        // gradient offset x  (absolute) [mT/m]  -->
	double  m_dGradOffsetY_mTm;                        // gradient offset x  (absolute) [mT/m]  -->
	double  m_dGradOffsetZ_mTm;                        // gradient offset x  (absolute) [mT/m]  -->
  bool    m_bGradOffsetsIsRelative;                  // gradient offsets are relative (additiv)
  size_t  m_iTxScaleFactorSize;                      // number of entries in m_aTxScaleFactor
  Complex m_aTxScaleFactor[ADJ_B1SHIM_MAX_CHANNELS]; // TX scale factors (indexed by TX channel number) -->
  bool    m_bTxScaleFactorIsRelative;                // TX scale factors (multiplicative)
};

#endif
