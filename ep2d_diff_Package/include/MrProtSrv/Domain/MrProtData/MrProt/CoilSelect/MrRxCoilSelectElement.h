//## begin module.cm preserve=no
//	---------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	---------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProt\CoilSelect\MrRxCoilSelectElement.h
//	 Version:
//	  Author: HEDER
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: MR::Measurement::CSequence::CoilSelect
//
//	 Classes:
//
//	---------------------------------------------------------

#pragma once

//-----------------------------------------------------------------------------
//  Includes
//-----------------------------------------------------------------------------
#include "MrProtSrv/Domain/CoreNative/MrCoilSelectData.h"

//-----------------------------------------------------------------------------
// Import/Export control
//-----------------------------------------------------------------------------
#undef __IMP_EXP

#ifdef BUILD_MrProt
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"


//-----------------------------------------------------------------------------
//  class prototype to define the friend state
//-----------------------------------------------------------------------------
class MrRxCoilSelect;

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4275) // non dll-interface class. Can be ignored if the whole code base is build with the same compiler and all modules share the same runtime
#endif

class __IMP_EXP MrRxCoilSelectElement: public MrProtocolData::MrRxCoilSelectElementDataDelegate
{
    typedef MrProtocolData::MrRxCoilSelectElementDataDelegate BasicImplementation;
  public:
      MrRxCoilSelectElement();
      MrRxCoilSelectElement(const MrRxCoilSelectElement& rSource);
      MrRxCoilSelectElement& operator=(const MrRxCoilSelectElement& rSource);
      explicit MrRxCoilSelectElement(MrProtocolData::MrRxCoilSelectElementData& data);
      MrRxCoilSelectElement& operator= (const MrProtocolData::MrRxCoilSelectElementData& rhs);
	  explicit MrRxCoilSelectElement(MrProtocolData::MrRxCoilSelectElementData::Pointer data);
      virtual ~MrRxCoilSelectElement(void);

    operator MrProtocolData::MrRxCoilSelectElementData*() { return m_pData.get(); }
    operator const MrProtocolData::MrRxCoilSelectElementData*() const { return m_pData.get(); }
  public:
      //	Returns the coilElementID struct, which identifies the
      //	coil element
      const MrProtocolData::MrRxCoilSelectElementData& getCoilElementData () const;
      MrProtocolData::MrRxCoilSelectElementData& getCoilElementData();

      bool isSelected () const;

      //	Switches a coil element to ON.
      void select ();

      //	Switches a coil element to OFF.
      void deselect ();

      //	Returns the coilElementID struct, which object the
      //	coil element
      const MrProtocolData::MrCoilElementIDData& getCoilElementID () const;

      //	Gets a reference to the coilElementID object.
      MrProtocolData::MrCoilElementIDData& getCoilElementID ();

      //	Returns the coilElementID struct, which object the
      //	coil element
      void setCoilElementID (MrProtocolData::MrCoilElementIDData&);

      //	Assings a coil element to an ADC channel.
      //	(If there is more than one coil element assigned to
      //	the same channel, then the ADC-signals are combined).
      void setADCChannel (const int32_t ADCChannel);

      //	Returns the number of the ADC channel, which is
      //	connected to the given coil element.
      int32_t getADCChannel () const;

      //	Assings a coil element to an RX channel.
      void setRxChannel (const int32_t RxChannel);

      //	Returns the number of the RX channel, to which the
      //	given coil element is connected.
      int32_t getRxChannel () const;

      //	The element given by coilElementID  is initialized with
      //	ADCChannel (nothing else is done). To deselect an
      //	element, just call init (id, 0).
      bool set(const MrProtocolData::MrCoilElementIDData &id, const int32_t ADCChannel);

      //	Assignment operator.
      
  private:

      friend class MrRxCoilSelect;
};

#undef __OWNER

#ifdef _MSC_VER
# pragma warning(pop)
#endif