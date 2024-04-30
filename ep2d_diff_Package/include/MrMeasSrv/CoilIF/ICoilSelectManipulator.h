// -----------------------------------------------------------------------------
// Copyright (C) Siemens Healthcare GmbH 2015 -2018 All Rights Reserved.
// Restricted
// -----------------------------------------------------------------------------

#pragma once

#ifndef ICoilSelectManipulator_h
#define ICoilSelectManipulator_h

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include <stdexcept>
#include "MrVista/Parc/Reflection/IObject.h"
#include "MrVista/Parc/Reflection/ModuleManager.h"
#include "MrMeasSrv/CoilIF/CoilSelectManipulatorStatus.h"
#include "MrProtSrv/Domain/CoreNative/MrCoilSelectData.h"              // wg. BCCMode
#include <iostream>
#include <set>

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------
class MeasDumpEnv;
class MeasCoilContext;
class AccCoilElement;
class MeasNucleus;

namespace MrProtocolData
{
  class MrProtData;
  class MrRxCoilSelectData;
  class MrTxCoilSelectData;
  class MrCoilPlugsData;
}


namespace
{
  struct ElementRange
  {
    int          iCoilElementIndex;
    uint32_t     uiMuxCh;
    uint32_t     uiIFindex;
    std::string  sElementUUID;
    // Possible connections Mux -> Rx
    std::set<int> possibleRxConn;

    bool operator== (const ElementRange &rhs) const
    {
      return ( iCoilElementIndex == rhs.iCoilElementIndex &&
               uiMuxCh           == rhs.uiMuxCh           &&
               uiIFindex         == rhs.uiIFindex         &&
               sElementUUID      == rhs.sElementUUID      &&
               possibleRxConn    == rhs.possibleRxConn);
    }
  };
}
 
struct sortElementRange
{
    bool operator()(const ElementRange& left, const ElementRange& right)
    {
        bool bResult;
        if (left.possibleRxConn.size() != right.possibleRxConn.size())
        {
            bResult = (left.possibleRxConn.size() < right.possibleRxConn.size());
        }
        else if (left.uiMuxCh != right.uiMuxCh)
        {
            bResult = (left.uiMuxCh < right.uiMuxCh);
        }
        else
        {
            bResult = (left.sElementUUID < right.sElementUUID);
        }

        return bResult;
    };
};


/// Class CoilSelectManipulator
class ICoilSelectManipulator
  : public Parc::Interface
{
  public:

    DECLARE_PARC_INTERFACE( ICoilSelectManipulator );
    
    /// Status bit constants
    enum
    {
      COILSEL_OK  = CoilSelectManipulatorStatus::COILSEL_OK,                    
      COILSEL_CHANGED = CoilSelectManipulatorStatus::COILSEL_CHANGED,               
      COILSEL_EMPTY = CoilSelectManipulatorStatus::COILSEL_EMPTY,                 
      ELEM_DESELECTED = CoilSelectManipulatorStatus::ELEM_DESELECTED,               
      ELEM_NOT_CONNECTED = CoilSelectManipulatorStatus::ELEM_NOT_CONNECTED,            
      NO_RX_CHAN_AVAIL = CoilSelectManipulatorStatus::NO_RX_CHAN_AVAIL,              
      RX_CHAN_INVAL = CoilSelectManipulatorStatus::RX_CHAN_INVAL,                 
      MUX_CHAN_INVAL = CoilSelectManipulatorStatus::MUX_CHAN_INVAL,                
      NO_TX_ELEM = CoilSelectManipulatorStatus::NO_TX_ELEM,                    
      CONNECTED_ELEM_CHANGED = CoilSelectManipulatorStatus::CONNECTED_ELEM_CHANGED,        
      COMP_ELEM_EXCHANGED = CoilSelectManipulatorStatus::COMP_ELEM_EXCHANGED,           
      UNNEEDED_ELEM_COMBINATION = CoilSelectManipulatorStatus::UNNEEDED_ELEM_COMBINATION,     
      ELEM_IS_NOT_SELECTED = CoilSelectManipulatorStatus::ELEM_IS_NOT_SELECTED,          
      NUCLEUS_NOT_SUPPORTED_BY_ELEM = CoilSelectManipulatorStatus::NUCLEUS_NOT_SUPPORTED_BY_ELEM, 
      NUCLEUS_DOES_NOT_MATCH = CoilSelectManipulatorStatus::NUCLEUS_DOES_NOT_MATCH,        
      PAT_ELEM_DESELECTED = CoilSelectManipulatorStatus::PAT_ELEM_DESELECTED,           
      CONNECTED_ELEM_INVAL = CoilSelectManipulatorStatus::CONNECTED_ELEM_INVAL,          
      ELEM_COULD_NOT_BE_SELECTED = CoilSelectManipulatorStatus::ELEM_COULD_NOT_BE_SELECTED,    
      ELEM_COULD_NOT_BE_DESELECTED = CoilSelectManipulatorStatus::ELEM_COULD_NOT_BE_DESELECTED,  
      BUILD_REQUIRED  = CoilSelectManipulatorStatus::BUILD_REQUIRED,              
      CSM_INTERNAL_ERROR = CoilSelectManipulatorStatus::CSM_INTERNAL_ERROR, 
      COILENV_INVALID = CoilSelectManipulatorStatus::COILENV_INVALID,               
      NOT_ENOUGH_RX_CHAN_AVAIL = CoilSelectManipulatorStatus::NOT_ENOUGH_RX_CHAN_AVAIL,      
      RX_CHAN_ALREADY_ASSIGNED = CoilSelectManipulatorStatus::RX_CHAN_ALREADY_ASSIGNED,      
      MUX_CONSTRAINT_CONFLICT = CoilSelectManipulatorStatus::MUX_CONSTRAINT_CONFLICT,       
      TX_CHAN_ALREADY_ASSIGNED = CoilSelectManipulatorStatus::TX_CHAN_ALREADY_ASSIGNED,
      BCCOMBINEMATRIX_COULD_NOT_SET = CoilSelectManipulatorStatus::BCCOMBINEMATRIX_COULD_NOT_SET,
      TX_CHAN_INVAL  = CoilSelectManipulatorStatus::TX_CHAN_INVAL,                 
      RX_CHAN_CHANGED  = CoilSelectManipulatorStatus::RX_CHAN_CHANGED,               
      TX_CHAN_CHANGED = CoilSelectManipulatorStatus::TX_CHAN_CHANGED,
      LC_ELEM_DESELECTED = CoilSelectManipulatorStatus::LC_ELEM_DESELECTED,
    };

    /// static create method
    static ICoilSelectManipulator::Pointer create();

    /// Init the class. One of these methods _MUST_ be called
    /// imediately after the default constructor
    virtual void init (MrProtocolData::MrCoilSelectData &rMrCoilSelectData) = 0;
    virtual void init (MrProtocolData::MrCoilSelectData &rMrCoilSelectData,
                       const MeasCoilContext            &rCoilContext) = 0;

    /// Fill the CoilSelectManipulator with all
    /// relevant information from the protocol
    virtual void setup (const MrProtocolData::MrProtData* pProt) = 0;

    /// setCoilContext
    virtual void setCoilContext (const MeasCoilContext& rCoilContext) = 0;

    /// getCoilContext
    virtual const MeasCoilContext& getCoilContext () const = 0;

    /// Set the nucleus in coilselect with given index
    /// without modifying the CoilSelect.
    /// The nucleus will be used with the manipulation.
    virtual bool setRxNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus) = 0;

    /// Set the TX nucleus in coilselect with given index
    /// without modifying the CoilSelect.
    /// The nucleus will be used with the manipulation.
    virtual bool setTxNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus) = 0;

    // Set the Tx and Rx nucleus in the coilselect with the given index
    // without modifying the CoilSelect.
    // The nucleus will be used with the manipulation.
    virtual bool setNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus) = 0;

    /// Set if combination allowed without modifying the
    /// current coilselect.
    virtual void setCombinationAllowed(bool bAllowed) = 0;

    /// Set the number ADC channels without modifying the coilselect.
    virtual void setNumberOfADCChannels(int iNoOfADCChs) = 0;

    /// Set, if PAT is used without modifying the current
    /// coilselect.
    ///
    /// This flag will be used to dermine the system default
    /// R-Factor.
    virtual void setPAT(bool bPAT) = 0;
    
    /// Get the PAT acceleration the given coil element.
    /// The element must be selected in coilselect with given index.
    virtual void getPATAcceleration(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement,
      float&          rflAccelerationX,
      float&          rflAccelerationY,
      float&          rflAccelerationZ) const = 0;
    
    /// Insert Body Coil in coilselect with given index,
    /// if there is a body coil.
    /// bAllElements is false (default): Insert only the first
    /// body coil element (for normal handling enough)
    /// bAllElements is true: Insert all body coil elements
    /// (Required for tuning, service-SW)
    virtual CoilSelectManipulatorStatus insertBodyCoil(
      unsigned int            uiIndex = 0,
      MrProtocolData::BCCMode eBCCMode = MrProtocolData::BCCM_UNDEFINED) = 0;

    /// Insert all Pickup coil elements in coilselect with given index.
    virtual CoilSelectManipulatorStatus insertPickupCoil(unsigned int uiIndex) = 0;

    /// Selects the given coil element with the specified
    /// R-Factor in coilselect with given index.
    ///
    /// If the element is already selected then R-Factor will
    /// be adapted if necessary.
    ///
    /// If a R-Factor of 0 is used or no R-Factor has been
    /// specified the R-Factor will be automaticaly choiced.
    ///
    /// Return the status of the modification.
    virtual CoilSelectManipulatorStatus selectRxElement(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement,
      int             iRFactorElement = 0) = 0;

    /// Select the given TX coil element in the TX coil select
    /// with index uiIndex.
    ///
    /// Return the status of the modification.
    virtual CoilSelectManipulatorStatus selectTxElement(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement) = 0;

    /// Deselects all elements in coilselect with given index
    /// which belongs to the PAT-Group of the given coil element.
    ///
    /// Return the status of the modification.
    virtual CoilSelectManipulatorStatus deselectRxElement(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement) = 0;

    /// Deselects TX element from TX coil select with index uiIndex
    ///
    /// Return the status of the modification.
    virtual CoilSelectManipulatorStatus deselectTxElement(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement) = 0;

    // Set SuppressMandatoryHandling in coilselect with given index
    // without modifying the CoilSelect.
    virtual bool setRxSuppressMandatoryProperties(unsigned int uiIndex, bool bSuppressMandatoryProperties) = 0;

    // Set ignore BC/LC excluding handling in coilselect with given index
    // without modifying the CoilSelect.
    virtual bool setRxIgnoreBCLCExcluding(unsigned int uiIndex, bool bIgnoreBCLCExcluding) = 0;

    // Set suppress coil element exclusive handling in coilselect with given index
    // without modifying the CoilSelect.
    virtual bool setRxSuppressExclusiveProperties(unsigned int uiIndex, bool bSuppressExclusiveProperties) = 0;

    /// Returns true, if the given coil element is selected in
    /// coilselect with given index.
    virtual bool isRxSelected(unsigned int uiIndex, AccCoilElement* pAccCoilElement) const = 0;

    /// Determine a possible CoilPlug object (PlugIDs) from a
    /// given coilSelect.
    ///
    /// All PlugIDs in the CoilPlug result should
    /// provide all coil elements which are
    /// selected in the coilSelect.
    virtual bool determineCoilPlug (
      MrProtocolData::MrCoilPlugsData&          rMrCoilPlugData,
      const MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData) = 0;

    /// Determine a possible CoilPlug object (PlugIDs) from the
    /// coilselect with given index.
    /// All PlugIDs in the CoilPlug result should provide all
    /// coil elements which are selected in the coilSelect.
    virtual bool determineCoilPlug(unsigned int uiIndex) = 0;

    /// Convert the coilselect with given index, determine coil plugs
    /// and make it valid.
    virtual CoilSelectManipulatorStatus adaptCoilPlugsAndSelect(unsigned int uiIndex) = 0;
    
    /// Return true, if no CoilPlug determine required.
    virtual bool isCoilPlugValid(const MrProtocolData::MrCoilPlugsData& rMrCoilPlugData) const = 0;

    /// If necessary, the coilselect with given index will be modified
    /// to create a valid RX and TX coilselect.
    ///
    /// If it is not possible to create a valid coilselect,
    /// then the original coilselect will not be modified.
    ///
    /// Return the status of the modification.
    virtual CoilSelectManipulatorStatus buildValid(unsigned int uiIndex) = 0;

    /// Deselect the given coil element and all other elements
    /// with the same PAT-group.

    /// If necessary, the coilselect with given index will be modified
    /// to create a valid coilselect.
    ///
    /// If it is not possible to create a valid coilselect,
    /// then the original coilselect will not be modified.
    ///
    /// Return the status of the modification.
    virtual CoilSelectManipulatorStatus buildValidRx(unsigned int uiIndex) = 0;

    /// Checks, if the coilselect with given index is valid. The
    /// coilselect will never be changed. Only a check will be
    /// performed.
    ///
    /// Return the status of the check.
    virtual CoilSelectManipulatorStatus check(unsigned int uiIndex) = 0;

    /// Checks, if the RX coilselect with given index is valid. The
    /// coilselect will never be changed. Only a check will be
    /// performed.
    ///
    /// Return the status of the check.
    virtual CoilSelectManipulatorStatus checkRx(unsigned int uiIndex) = 0;

    /// Build a valid TX select
    ///
    /// Return the status of the modification.
    virtual CoilSelectManipulatorStatus buildValidTx(unsigned int uiIndex) = 0;

    /// Checks, if the TX coilselect with given index is valid. The
    /// coilselect will never be changed. Only a check will be
    /// performed.
    ///
    /// Return the status of the check.
    virtual CoilSelectManipulatorStatus checkTx(unsigned int uiIndex) = 0;

    /// Get RX coil select at uiIndex 
    virtual MrProtocolData::MrRxCoilSelectData* getRxCoilSelect(unsigned int uiIndex) = 0;

    /// Todo (HJe): Comment???
    virtual bool refreshRxSelect (
      unsigned int                        uiIndex,
      MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData) = 0;

    
    /// Gets the maximum usable R-Factor for the given coil
    /// element.
    ///
    /// It is not required, that the element is selected. It
    /// must only be connected.
    ///
    /// The value depends on the maximum supported R-Factor of
    /// the element and the remaining number of receiver channels.
    virtual int getMaxRFactor(
      MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData,
      AccCoilElement*                     pAccCoilElement) const = 0;

    /// Get the maximum R factor for selected coil elements in coilselect
    /// with given index.
    virtual int getMaxRFactor(unsigned int uiIndex) const = 0;

    /// Sets the R-Factor for the given selected coil element
    /// in coilselect with given index.
    ///
    /// If the element is not selected nothing will be changged.
    ///
    /// Return the status of the modification.
    virtual CoilSelectManipulatorStatus setRFactor(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement,
      int             iRFactorElement) = 0;

    /// Calls the internal method getRFactor() for getting the
    /// R factor for the given coil element.
    /// The element must be selected in coilselect with given index.
    virtual int getRFactor(unsigned int uiIndex, AccCoilElement* pAccCoilElement) const = 0;

    /// Set the global R-Factor in coilselect with given index
    /// and adjust the current selection of the elements with
    /// this new factor.
    ///
    /// If the input factor is 0, the system default R-Factor
    /// will be used.
    virtual CoilSelectManipulatorStatus setGlobalRFactor(
      unsigned int uiIndex,
      int          iRFactor) = 0;

    /// Get the global R factor in coilselect with given index.
    virtual int getGlobalRFactor(unsigned int uiIndex) const = 0;

    /// Get the system maximum R factor.
    ///
    /// If the automatic R factor adaption disabled ->
    /// The value is dependent from the number of ADCs and if
    /// PAT reqired.
    ///
    /// Number of ADCs <= 8: R-factor = 1
    /// Number of ADCs >  8: R-factor = 3
    ///
    /// If the automatic R factor adaption enabled ->
    /// The value is always 3
    virtual int getSystemMaximumRFactor() const = 0;

    /// Returns the state of the IndividualRFactors in the given coilSelect.
    virtual bool isIndividualRFactors(
      const MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData) const = 0;

    /// Returns the state of the IndividualRFactors in coilselect
    /// with given index
    virtual bool isIndividualRFactors(unsigned int uiIndex) const = 0;

    /// Return the mode of RFactor in coilselect with given index.
    /// Following strings will be returned
    /// - "V:" -> Varius RFactors    (isIndividualRFactors() is true)
    /// - "C:" -> CP-Mode      (RFactor 0)
    /// - "D:" -> Dual-Mode    (RFactor 1)
    /// - "T:" -> Tribble-Mode (RFactor 2)
    ///
    /// If the character is in lower case, the RFactor has been
    /// specified as GlobalRFactor in the Protocol, if upper
    /// it is the system default RFactor.
    virtual const char * getModeDesc(unsigned int uiIndex) const = 0;
      
    ///------------------------------------------------------------------------
    /// misc. dump methods      
    ///------------------------------------------------------------------------

    /// Dump all selected coil elements.
    virtual void dump() = 0;

    /// Dump all selected coil elements in coilselect with given index.
    virtual void dump(unsigned int uiIndex) = 0;

    /// Dump all selected coil elements in the given
    /// MeasDumpEnv.
    virtual void dump(MeasDumpEnv& rMeasDumpEnv) = 0;

    /// Dump all selected coil elements in coilselect with given index
    /// to the given MeasDumpEnv.
    virtual void dump(unsigned int uiIndex, MeasDumpEnv& rMeasDumpEnv) = 0;

    /// Dump the plugs
    virtual void dumpPlugs(MeasDumpEnv& rMeasDumpEnv) const = 0;

    /// Dump all rx connections of the coilselect with given index
    virtual void dumpRxConn(std::ostream &rStream, unsigned int uiIndex) const = 0;

    /// Test methode to call many methodes of the CoilSelectManipulator
    /// interactive about a SelectionMenu
    virtual void edit() = 0;
    
    /// Get the internal debug level
    /// If the internal not set, get the debug level from the
    /// measperm
    virtual int getDebugLevel(void ) const = 0;

    /// Set the internal debug level
    /// The debug level controls the trace outputs
    virtual void setDebugLevel(int iDebugLevel) = 0;
    
    /// Check the RxChannels/ADCs of the RX coilselect with given index.
    /// If all Rx channels are valid, do nothing.
    /// If one or more Rx channels are invalid, determine a new one and set in into  
    /// the Rx coilselect
    /// Return the status of the check.
    virtual CoilSelectManipulatorStatus checkAndFixRxChannels(unsigned int uiIndex) = 0;
    
    /// Check the TxChannels of the TX coilselect with given index.
    /// If all Tx channels are valid, do nothing.
    /// If one or more Tx channels are invalid, determine a new one and set in into  
    /// the Tx coilselect
    /// Return the status of the check.
    virtual CoilSelectManipulatorStatus checkAndFixTxChannels(unsigned int uiIndex) = 0;


    /// unit test support
    virtual bool reloadRxConnections() = 0;
    
    /// unit test support
    /// Enable or Disable the overlap check
    virtual void  setOverlapCheck(bool bOverlapCheck) = 0;

    /// unit test support
    /// Set the scan position in z direction (for test)
    virtual void  setScanPosition(int iPosition) = 0;
};


inline  ICoilSelectManipulator::Pointer ICoilSelectManipulator::create()
{
  Pointer ptr;
  (void) ptr.CreateObject("CoilSelectManipulatorImpl@CoilManipulator");

  return ptr;
}



#endif
