//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 2014  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/X
//	    File: \src\MrImaging\seq\a_tgse_asl\AslSL_ID.h
//	  Author: pfeujodj
//	    Date: 2015-02-19 11:28:51 +01:00
//
//	    Lang: C++
//
//	 Descrip: 3D ASL sequence with CompositeSeqLoop support
//
//	-----------------------------------------------------------------------------

#pragma once

#include "MrImagingFW/libCSL/Key.h"
#include "MrImagingFW/libCSL/StdTypes.h"
#include "MrImagingFW/libCSL/ConcInfo.h"


// Forward declaration
class MrProt;


namespace SL_NAME_SPACE 
{
  class AslKeys 
  {  
  public:
    static const Key<int32_t>  l_SEGMENTS_COUNTER;
    static const Key<int32_t>  l_SEGMENTS_LENGTH;

    static const Key<int32_t>  l_ASL_LOOP_LENGTH;
    static const Key<int32_t>  l_ASL_LOOP_COUNTER;
    // --------------------------------------------------------------------
    /// Defines if first scan is M0Scan
    // --------------------------------------------------------------------
    static const Key<bool>     b_ENABLE_FIRST_PREPSCAN_AS_M0SCAN;
    // --------------------------------------------------------------------
    /// Defines if first M0Scan should be performed.
    // --------------------------------------------------------------------    
    const static Key <bool>    b_PERF_M0_SCAN;
    // --------------------------------------------------------------------
    /// Defines if loop is currently in M0Scan
    // --------------------------------------------------------------------
    const static Key <bool>    b_IS_M0_SCAN;
    // --------------------------------------------------------------------
    /// Defines if protocol should be overwritten for M0 measurement
    // --------------------------------------------------------------------
    const static Key <bool>    b_OVERWRITE_M0_PROTOCOL;
    // --------------------------------------------------------------------
    /// New protocol to be used for M0 measurement
    // --------------------------------------------------------------------
    const static Key <MrProt*> p_M0_PROTOCOL;
  
  };

} // SL_NAME_SPACE