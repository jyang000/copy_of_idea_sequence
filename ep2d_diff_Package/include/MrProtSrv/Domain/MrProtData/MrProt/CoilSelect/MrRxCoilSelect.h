//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4_servers1\pkg\MrServers\MrProtSrv\MrProt\CoilSelect\MrRxCoilSelect.h
//	 Version: \main\6
//	  Author: heumthte
//	    Date: 2011-05-09 13:27:38 +02:00
//
//	    Lang: C++
//
//	 Descrip: MR::Measurement::CSequence::CoilSelect
//
//	 Classes:
//
//	-----------------------------------------------------------------------------

#ifndef MrRxCoilSelect_h
#define MrRxCoilSelect_h 1

#ifdef _MSC_VER
#pragma once
#endif

#include "MrProtSrv/Domain/CoreNative/MrCoilSelectData.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/CoilSelect/MrRxCoilSelectElement.h"

class MeasNucleus;

//-----------------------------------------------------------------------------
// Import/Export control
//-----------------------------------------------------------------------------
#undef __IMP_EXP

#ifdef BUILD_MrProt
  #define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4275) // non dll-interface class. Can be ignored if the whole code base is build with the same compiler and all modules share the same runtime
#endif

class __IMP_EXP MrRxCoilSelect: public MrProtocolData::MrRxCoilSelectDataDelegate
{
    typedef MrProtocolData::MrRxCoilSelectDataDelegate BasicImplementation;
    friend class CoilSelectManipulator;

  public:
    MrRxCoilSelect();
    MrRxCoilSelect(const MrRxCoilSelect& rSource);
    MrRxCoilSelect(const MrProtocolData::MrRxCoilSelectData& data);
    MrRxCoilSelect(MrProtocolData::MrRxCoilSelectData::Pointer pData);
    virtual ~MrRxCoilSelect(void);

    MrProtocolData::MrRxCoilSelectData& getCoilSelectData();
    const MrProtocolData::MrRxCoilSelectData& getCoilSelectData() const;

    //	Returns const MrRxCoilSelectElement object.
    //	0 <= index < getTotalNumberOfElements().
    const MrRxCoilSelectElement element (const int32_t index) const;
    //	Returns MrRxCoilSelectElement object.
    //	0 <= index < getTotalNumberOfElements().
    MrRxCoilSelectElement element (const int32_t index);

    //	Returns MrRxCoilSelectElement object.
    //	0 <= index < getTotalNumberOfEle  ments().
    MrRxCoilSelectElement operator [] (const int32_t index);

    //	Assignment operator.
    MrRxCoilSelect& operator = (const MrRxCoilSelect& rhs);
    MrRxCoilSelect& operator = (const MrProtocolData::MrRxCoilSelectData& rhs);

    bool operator == (const MrRxCoilSelect &rhs) const;
    bool operator == (const MrProtocolData::MrRxCoilSelectData &rhs) const;

    bool operator != (const MrRxCoilSelect &rhs) const;
    bool operator != (const MrProtocolData::MrRxCoilSelectData &rhs) const;

    //Stream-Operator
     friend __IMP_EXP std::ostream& operator<< (std::ostream& rStream, const MrRxCoilSelect& rMrRxCoilSelect);

    //	Search for a free Coil-Select element and return the
    //  index. If there is no free element -1 will be returned.
    int32_t findFreeSelectionElement () const;
    //	Search for a Coil-Selection with a gived CoilElementID.
    //
    //	If found, the index will be returned else -1.
    int32_t findSelectionElementByID (const MrProtocolData::MrCoilElementIDData &id) const;

    //	Vergleicht zwei CoilSelects.
    //	Es wird vor dem Vergleich aber temporaer jweils ein
    //  "buildUniformedSelect" darauf angewendet.
    bool cmpUniformed (const MrRxCoilSelect &rCoilSelect) const;

    bool isEmpty () const;

    //	get number of used ADC Channels
    //  keep this int not int32_t to be compatible with yapsinit.cpp
    int getNumOfUsedADCChan  () const;

    /// Get the number of application channels.
    /// The BC combine mode is consided properly.
    /// If only LC elements are selected, the # of channels is equal with getNumOfUsedADCChan()
    uint32_t getNoOfApplChannels () const;

    //	is given ADC Channel used?
    bool isADCChanUsed (int32_t iADCChan) const;

    //	Returns number of selected elements
    int32_t getNumberOfSelectedElements () const;

    /**  Returns the global R-Factor.
     *   0 : System-Default (e.g. = 1 for 8-channel, = 3 for  32-channel)
     *  >0: other global R-Factor
     */
    int32_t getGlobalRFactor() const;

    /**  Sets the global R-Factor.
     *  0 : System-Default (e.g. = 1 for 8-channel, = 3 for
     *  32-channel)
     *  >0: other global R-Factor
     */
    void setGlobalRFactor(int32_t iGlobalRFactor);

     /**
      * Actually used R Factor
      */
    int32_t getUsedRFactor() const;

    /**
     * @sa getUsedRFactor
     */
    void setUsedRFactor(int32_t iUsedRFactor);

    //  Checks if the flag for the IndividualRFactors is set.
    bool isIndividualRFactors() const;

    //  Sets the flag for the IndividualRFactors.
    //  This flag indicates whether the Coil-Select uses
    //  individuel R-Factors (which are set in the expert-UI)
    //  for the elements, or not.
    void setIndividualRFactors(bool bIndividualRFactors);

    //	Returns reference to array of FFT Scale factors.
    const MrGenericDC::IVectorParameter< MrProtocolData::MrFFTScaleData >& FFTScale(void) const;
    //	Returns reference to array of FFT Scale factors.
    MrGenericDC::IVectorParameter< MrProtocolData::MrFFTScaleData >& FFTScale(void);
    // reset fft scale
    void resetFFTScale(void);

    //	The element given by   CoilElementID will be "connected"
    //  to ADCChannel and the state will be changed to "selected".
    //
    //	If the element was already connected, the ADCChannel
    //  will be replaced by the new one.
    bool add (const MrProtocolData::MrCoilElementIDData &id, int32_t lADCChannel = 0);

    //	The element given by   CoilElementID will be "connected"
    //  to ADCChannel and the state will be changed to "selected".
    //
    //	If the element was already connected, the ADCChannel
    //  will be replaced by the new one.
    bool add (const MrProtocolData::MrCoilElementIDData &id, int32_t lADCChannel, uint32_t ulTimeStamp);

    //	The element given by   CoilElementID will be "connected"
    //  to RxChannel and the state will be changed to "selected".
    //
    //	A Coil-Selection selected by an CoilElementID will be removed.
    void remove (const MrProtocolData::MrCoilElementIDData &id);

    //  Return true, if the nucleus in the CoilSelect not empty.
    bool isValid () const;

    //  Return true, if the given nucleus equal the nucleus
    //  in the CoilSelect.
    bool isValid(MeasNucleus& rNucleus) const;
    //  Set the given nucleus into the CoilSelect.
    void setNucleus(MeasNucleus& rNucleus);

    //	Deselects all elements:
    void deselectAll (void );

    //	Returns the max.  number of selected coil elements
    //  supported by Mr Coilselect. (This is identical to
    //  max.  number of connected coil elements.)
    int32_t getTotalNumberOfElements (void ) const;

    //	Setzt den Speicher aller nicht selektierten Elemente auf '\0'.
    void clearAllUnselectedElements ();

    //	Erzeugt eine einheitliche Selektion.
    //
    //	Folgende Dinge werden getan:
    //	1. Eine einheitliche Sortierung
    //	    aller CoilSelectElements
    //	2. Neuvergabe aller ADCChannels
    //	    (Das Kombinieren der Elemente
    //	    wird aber nicht veraendert)
    //
    //	ACHTUNG: Das erzeugte CoilSelect DARF NICHT direkt
    //  fuer eine Messung benutzt werden
    //	(die erzeugten Verschaltungen koennen verboten sein).
    void buildUniformedSelect ();

  private:
    //	shift selection element form index indSrc to indDest
    bool shift (uint32_t indSrc, uint32_t indDest);
};

#undef __OWNER

#ifdef _MSC_VER
# pragma warning(pop)
#endif

#endif
