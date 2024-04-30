//  -----------------------------------------------------------------------------
//    Copyright (C) Siemens AG 2013  All Rights Reserved.
//  -----------------------------------------------------------------------------
//
//   Project: NUMARIS/X
//      File: \src\MrImaging\seq\a_tgse_asl\ASLDefines.h
//    Author: pfeujodj
//      Date: 2018-09-12 08:55:49 +02:00
//
//      Lang: C++
//
//   Descrip: A collection of definitions and tools that are required to store and retrieve ASL data from
//            the User Parameters in the various IDEA/ICE storage objects (MrProt, SEQExpo etc)
//            !! to be shared between SEQ and ICE
//
//  -----------------------------------------------------------------------------
#pragma once

#include "MrImagingFW/libSBBFW/StdSeqIF.h"

//
// ASL definitions
//
namespace SEQ_ASL
{
  //
  // mode of prescan to be performed
  //
  enum ASL_PRESCAN 
  {
      PRESCAN_NONE          = 1
    , PRESCAN_FULL          = 2     // inline full M0 scan
    , PRESCAN_SEPARATE_SCAN = 3     // sep. scan with RF of ASL module switched off
  };
  
  //
  // ASL data types using in SEQ_ASL::ASLBinaryFlagEncoding
  //
  enum AslImageType
  {
      IMAGE_CONTROL     = 0x01  // control image (0001)
    , IMAGE_LABEL       = 0x02  // label image   (0010)
    , IMAGE_PRESCAN     = 0x08  // prescan       (1000)
    , IMAGE_SUBTRACTION = 0x20  // subtraction image (0010 0000)
  };

  const long PRESCAN_TR_DEFAULT_MS = 5000;

  //
  // Export locations in SeqExpo (have from USER_0 to USER_9 to fill)
  //
  // NOTE: EpiNav uses ICE_PROGRAM_PARA_USER_0-5 (MrImaging\seq\common\MocoFramework\Navigators\SBBEPINav3D.cpp)
  const long  ASL_EXPO_CBF_BAT_POS      = ICE_PROGRAM_PARA_USER_7;  //[17]  Map delivered (BAT and multi-TI quant)
  const long  ASL_EXPO_NUMOF_TI_POS     = ICE_PROGRAM_PARA_USER_8;  //[18]
  const long  ASL_EXPO_NUMOF_ENC_POS    = ICE_PROGRAM_PARA_USER_9;  //[19]  //number of encodings (for ASL this will be 2: label/control)
  
  //
  // ASL encoding type storage location in MDH.IceProgramPara()
  enum MDH_STORAGE_POSITION
  {
    MDH_ICEPROGPARA_POS_ASL_TYPE = 0,   //position in MDH where scan is tagged according ASL_IMAGE_TYPE: CONTROL,LABEL,PRESCAN etc
    MDH_ICEPROGPARA_POS_CSET = 11,      //position in MDH where SET counter is stored
    MDH_ICEPROGPARA_POS_CPHS = 12,      //position in MDH where PHS counter is stored
    MDH_ICEPROGPARA_POS_CREP = 13,      //position in MDH where REP counter is stored
    MDH_ICEPROGPARA_POS_CIDA = 14       //position in MDH where IDA counter is stored
  };

  // Flag for relCBF and BAT map delivery
  enum CBF_BAT_STATE
  {
      CBF_BAT_NONE = 0x00,  //00 00, Neither CBF nor BAT
      CBF_ONLY     = 0x02,  //00 10, Only CBF 
      CBF_BAT_ALL  = 0x10   //10 10, Both CBF and BAT
  };

  // class to handle bitwise operations
  // make it a template, as we don't know what we are going to store.
  // it is also possible that the mask is bigger than the array we'd like to store. (we are most likely storing in uint16 SHORT ).
  // it is possible that our mask needs to be bigger than that. E.G, in hadamard.
  template<class storeType, class flagType>
  class ASLBinaryFlagEncoding 
  {
  public:
      // constructor
      ASLBinaryFlagEncoding(storeType* pBinArray=NULL) 
        : m_pDataArray(NULL) 
      {
        // how many data blocks do we need to encode over?
        m_sNumBlocks = sizeof(flagType)/sizeof(storeType);

        // how much of a bitshift do we need to encode over multiple blocks?
        m_sBitShift = sizeof(storeType)*CHAR_BIT;

        // init data array
        m_pDataArray = new storeType[m_sNumBlocks];

        if (pBinArray != NULL) 
        {
          for(int i=0;i<m_sNumBlocks;i++) 
          {
            m_pDataArray[i] = pBinArray[i];
          }
        }
        else
        {
          reset();    // init with zero
        }
      }

      // destructor
      ~ASLBinaryFlagEncoding( void ) 
      {
        if( m_pDataArray != NULL )
        {
          delete[] m_pDataArray;  
        }
        m_pDataArray = NULL;
      }

      bool addFlag(flagType tFlag) 
      {
        if( m_pDataArray == NULL ) 
        {
          return false;
        }

        // check that flag is already set.if it is set, then do nothing.
        if( !isFlagSet(tFlag) ) 
        {
          for(short i=0;i<m_sNumBlocks;i++) 
          {
            m_pDataArray[i] += static_cast<storeType>(tFlag >> m_sBitShift*i);
          }
        }

        return true;
      }

      bool removeFlag(flagType tFlag) 
      {
        if( m_pDataArray == NULL ) 
        {
          return false;
        }

        // check that flag is not set. if it is set, then remove
        if( isFlagSet(tFlag) ) 
        {
          for(short i=0;i<m_sNumBlocks;i++) 
          {
            m_pDataArray[i] -= static_cast<storeType>(tFlag >> m_sBitShift*i);
          }
        }

        return true;
      }

      bool isFlagSet(flagType tFlag) 
      {
        if( m_pDataArray == NULL ) 
        {
          return false;
        }

        bool bRetValue = false;
        flagType tInternalFlag = static_cast<flagType>(0);

        for(short i=0;i<m_sNumBlocks;i++) 
        {
          tInternalFlag += static_cast<flagType>(m_pDataArray[i] << m_sBitShift*i);
        }

        // perform logical & to determine if flag is present
        if( (tInternalFlag & tFlag) == tFlag) 
        {
          bRetValue = true;
        }

        return bRetValue;
      }

      //  override operator to return the encoded data array
      storeType operator[](uint16_t idx) 
      {
        if( (m_pDataArray == NULL) || (idx >= m_sNumBlocks) ) 
        {
          return static_cast<storeType>(0);
        }

        return m_pDataArray[idx];
      }


      bool reset( void ) 
      {
        if( m_pDataArray == NULL ) 
        {
          return false;
        }

        for(int i=0;i<m_sNumBlocks;i++) 
        {
          m_pDataArray[i] = static_cast<storeType>(0);
        }

        return true;
      }

  protected:
      storeType* m_pDataArray;
      short m_sNumBlocks;
      short m_sBitShift;
  };

}  //namespace: SEQ_ASL
