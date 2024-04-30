// -----------------------------------------------------------------------------
// Copyright (C) Siemens Healthcare GmbH 2015 -2018 All Rights Reserved.
// Restricted
// -----------------------------------------------------------------------------


#pragma once


#ifndef CoilSelectManipulator_h
#define CoilSelectManipulator_h 1

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include "MrMeasSrv/CoilIF/ICoilSelectManipulator.h"

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------
class MeasDumpEnv;
class MeasCoilContext;
class AccCoilElement;
class CoilPlugs;
class MeasNucleus;
class MrRxCoilSelect;

namespace MrProtocolData
{
  class MrProtData;
  class MrRxCoilSelectData;
  class MrCoilPlugsData;
  class MrCoilSelectData;
}

/// Class CoilSelectManipulator
/// PImpl to CoilSelectManipulatorImpl
class CoilSelectManipulator
{

  public:

       /// Constructor
    CoilSelectManipulator (MrProtocolData::MrCoilSelectData&  rMrCoilSelectData);
    CoilSelectManipulator (MrProtocolData::MrCoilSelectData&  rMrCoilSelectData,
                           const MeasCoilContext&             rCoilContext);

    /// Destructor (generated)
    virtual ~CoilSelectManipulator();

    /// allows reuse of coil select manipulator
    void init(MrProtocolData::MrCoilSelectData&  rMrCoilSelectData,
              const MeasCoilContext&             rCoilContext);
    
    /// Fill the CoilSelectManipulator with all
    /// relevant information from the protocol
    void setup (const MrProtocolData::MrProtData* pProt);

    /// setCoilContext
    void setCoilContext (const MeasCoilContext& rCoilContext);

    /// getCoilContext
    const MeasCoilContext& getCoilContext () const;

    /// Set the nucleus in coilselect with given index
    /// without modifying the CoilSelect.
    /// The nucleus will be used with the manipulation.
    bool setRxNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus);

    /// Set the TX nucleus in coilselect with given index
    /// without modifying the CoilSelect.
    /// The nucleus will be used with the manipulation.
    bool setTxNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus);

    // Set the Tx and Rx nucleus in the coilselect with the given index
    // without modifying the CoilSelect.
    // The nucleus will be used with the manipulation.
    bool setNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus);

    /// Set if combination allowed without modifying the
    /// current coilselect.
    void setCombinationAllowed(bool bAllowed);

    /// Set the number ADC channels without modifying the coilselect.
    void setNumberOfADCChannels(int iNoOfADCChs);

    /// Set, if PAT is used without modifying the current
    /// coilselect.
    ///
    /// This flag will be used to dermine the system default
    /// R-Factor.
    void setPAT(bool bPAT);
    
    /// Get the PAT acceleration the given coil element.
    /// The element must be selected in coilselect with given index.
    void getPATAcceleration(unsigned int    uiIndex,
                            AccCoilElement* pAccCoilElement,
                            float&          rflAccelerationX,
                            float&          rflAccelerationY,
                            float&          rflAccelerationZ) const;
    
    /// Insert Body Coil in coilselect with given index,
    /// if there is a body coil.
    /// bAllElements is false (default): Insert only the first
    /// body coil element (for normal handling enough)
    /// bAllElements is true: Insert all body coil elements
    /// (Required for tuning, service-SW)
    CoilSelectManipulatorStatus insertBodyCoil(
            unsigned int            uiIndex = 0,
            MrProtocolData::BCCMode eBCCMode = MrProtocolData::BCCM_UNDEFINED);

    /// Insert all Pickup coil elements in coilselect with given index.
    CoilSelectManipulatorStatus insertPickupCoil(unsigned int uiIndex);

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
    CoilSelectManipulatorStatus selectRxElement(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement,
      int             iRFactorElement = 0);

    /// Select the given TX coil element in the TX coil select
    /// with index uiIndex.
    ///
    /// Return the status of the modification.
    CoilSelectManipulatorStatus selectTxElement(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement);

    /// Deselects all elements in coilselect with given index
    /// which belongs to the PAT-Group of the given coil element.
    ///
    /// Return the status of the modification.
    CoilSelectManipulatorStatus deselectRxElement(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement);

    /// Deselects TX element from TX coil select with index uiIndex
    ///
    /// Return the status of the modification.
    CoilSelectManipulatorStatus deselectTxElement(
      unsigned int    uiIndex,
      AccCoilElement* pAccCoilElement);
    
    // Set SuppressMandatoryHandling in coilselect with given index
    // without modifying the CoilSelect.
    bool setRxSuppressMandatoryProperties(unsigned int uiIndex, bool bSuppressMandatoryProperties);

    // Set ignore BC/LC excluding handling in coilselect with given index
    // without modifying the CoilSelect.
    bool setRxIgnoreBCLCExcluding(unsigned int uiIndex, bool bIgnoreBCLCExcluding);

    // Set suppress coil element exclusive handling in coilselect with given index
    // without modifying the CoilSelect.
    bool setRxSuppressExclusiveProperties(unsigned int uiIndex, bool bSuppressExclusiveProperties);
    
    /// Returns true, if the given coil element is selected in
    /// coilselect with given index.
    bool isRxSelected(unsigned int uiIndex, AccCoilElement* pAccCoilElement) const;

    /// Determine a possible CoilPlug object (PlugIDs) from a
    /// given coilSelect.
    ///
    /// All PlugIDs in the CoilPlug result should
    /// provide all coil elements which are
    /// selected in the coilSelect.
    bool determineCoilPlug (MrProtocolData::MrCoilPlugsData&          rMrCoilPlugData,
                            const MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData);

    /// Determine a possible CoilPlug object (PlugIDs) from the
    /// coilselect with given index.
    /// All PlugIDs in the CoilPlug result should provide all
    /// coil elements which are selected in the coilSelect.
    bool determineCoilPlug(unsigned int uiIndex);

    /// Convert the coilselect with given index, determine coil plugs
    /// and make it valid.
    CoilSelectManipulatorStatus adaptCoilPlugsAndSelect(unsigned int uiIndex);
    
    /// Return true, if no CoilPlug determine required.
    bool isCoilPlugValid(const MrProtocolData::MrCoilPlugsData& rMrCoilPlugData) const;


    /// If necessary, the coilselect with given index will be modified
    /// to create a valid RX and TX coilselect.
    ///
    /// If it is not possible to create a valid coilselect,
    /// then the original coilselect will not be modified.
    ///
    /// Return the status of the modification.
    CoilSelectManipulatorStatus buildValid(unsigned int uiIndex);

    /// Deselect the given coil element and all other elements
    /// with the same PAT-group.

    /// If necessary, the coilselect with given index will be modified
    /// to create a valid coilselect.
    ///
    /// If it is not possible to create a valid coilselect,
    /// then the original coilselect will not be modified.
    ///
    /// Return the status of the modification.
    CoilSelectManipulatorStatus buildValidRx(unsigned int uiIndex);

    /// Checks, if the coilselect with given index is valid. The
    /// coilselect will never be changed. Only a check will be
    /// performed.
    ///
    /// Return the status of the check.
    CoilSelectManipulatorStatus check(unsigned int uiIndex);

    /// Checks, if the RX coilselect with given index is valid. The
    /// coilselect will never be changed. Only a check will be
    /// performed.
    ///
    /// Return the status of the check.
    CoilSelectManipulatorStatus checkRx(unsigned int uiIndex);

    /// Build a valid TX select
    ///
    /// Return the status of the modification.
    CoilSelectManipulatorStatus buildValidTx(unsigned int uiIndex);

    /// Checks, if the TX coilselect with given index is valid. The
    /// coilselect will never be changed. Only a check will be
    /// performed.
    ///
    /// Return the status of the check.
    CoilSelectManipulatorStatus checkTx(unsigned int uiIndex);

    /// Get RX coil select at uiIndex 
    MrProtocolData::MrRxCoilSelectData* getRxCoilSelect(unsigned int uiIndex);

    /// Todo (HJe): Comment???
    bool refreshRxSelect (unsigned int                        uiIndex,
                          MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData);

    
    /// Gets the maximum usable R-Factor for the given coil
    /// element.
    ///
    /// It is not required, that the element is selected. It
    /// must only be connected.
    ///
    /// The value depends on the maximum supported R-Factor of
    /// the element and the remaining number of receiver channels.
    int getMaxRFactor(MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData,
                      AccCoilElement*                     pAccCoilElement) const;

    /// Get the maximum R factor for selected coil elements in coilselect
    /// with given index.
    int getMaxRFactor(unsigned int uiIndex) const;

    /// Sets the R-Factor for the given selected coil element
    /// in coilselect with given index.
    ///
    /// If the element is not selected nothing will be changged.
    ///
    /// Return the status of the modification.
    CoilSelectManipulatorStatus setRFactor(unsigned int    uiIndex,
                                           AccCoilElement* pAccCoilElement,
                                           int             iRFactorElement);

    /// Calls the internal method getRFactor() for getting the
    /// R factor for the given coil element.
    /// The element must be selected in coilselect with given index.
    int getRFactor(unsigned int uiIndex, AccCoilElement* pAccCoilElement) const;

    /// Set the global R-Factor in coilselect with given index
    /// and adjust the current selection of the elements with
    /// this new factor.
    ///
    /// If the input factor is 0, the system default R-Factor
    /// will be used.
    CoilSelectManipulatorStatus setGlobalRFactor(unsigned int uiIndex,
                                                 int          iRFactor);

    /// Get the global R factor in coilselect with given index.
    int getGlobalRFactor(unsigned int uiIndex) const;

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
    int getSystemMaximumRFactor() const;

    /// Returns the state of the IndividualRFactors in the given coilSelect.
    bool isIndividualRFactors(const MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData) const;

    /// Returns the state of the IndividualRFactors in coilselect
    /// with given index
    bool isIndividualRFactors(unsigned int uiIndex) const;

    /// Return the mode of RFactor in coilselect with given index.
    /// Following strings will be returned
    /// - "V:" -> Varius RFactors    (isIndividualRFactors() is
    /// true)
    /// - "C:" -> CP-Mode      (RFactor 0)
    /// - "D:" -> Dual-Mode    (RFactor 1)
    /// - "T:" -> Tribble-Mode (RFactor 2)
    ///
    /// If the character is in lower case, the RFactor has been
    /// specified as GlobalRFactor in the Protocol, if upper
    /// it is the system default RFactor.
    const char * getModeDesc(unsigned int uiIndex) const;
      
    ///------------------------------------------------------------------------
    /// misc. dump methods      
    ///------------------------------------------------------------------------

    /// Dump all selected coil elements.
    void dump();

    /// Dump all selected coil elements in coilselect with given index.
    void dump(unsigned int uiIndex);

    /// Dump all selected coil elements in the given
    /// MeasDumpEnv.
    void dump(MeasDumpEnv& rMeasDumpEnv);

    /// Dump all selected coil elements in coilselect with given index
    /// to the given MeasDumpEnv.
    void dump(unsigned int uiIndex, MeasDumpEnv& rMeasDumpEnv);

    /// Dump the plugs
    void dumpPlugs(MeasDumpEnv& rMeasDumpEnv) const;

    /// Dump all rx connections of the coilselect with given index
    void dumpRxConn(std::ostream &rStream, unsigned int uiIndex) const;

    /// Test methode to call many methodes of the CoilSelectManipulator
    /// interactive about a SelectionMenu
    void edit();
    
    /// Get the internal debug level
    /// If the internal not set, get the debug level from the
    /// measperm
    int getDebugLevel(void ) const;

    /// Set the internal debug level
    /// The debug level controls the trace outputs
    void setDebugLevel(int iDebugLevel);
    
    /// Check the RxChannels/ADCs of the RX coilselect with given index.
    /// If all Rx channels are valid, do nothing.
    /// If one or more Rx channels are invalid, determine a new one and set in into  
    /// the Rx coilselect
    /// Return the status of the check.
    CoilSelectManipulatorStatus checkAndFixRxChannels(unsigned int uiIndex);
    
    /// Check the TxChannels of the TX coilselect with given index.
    /// If all Tx channels are valid, do nothing.
    /// If one or more Tx channels are invalid, determine a new one and set in into  
    /// the Tx coilselect
    /// Return the status of the check.
    CoilSelectManipulatorStatus checkAndFixTxChannels(unsigned int uiIndex);

  private:
    
    /// Pointer to implementation 
    /// ICoilSelectManipulator::Pointer  m_pImpl;
    ICoilSelectManipulator::Pointer  m_pImpl;
    
    /// Not implemented
    CoilSelectManipulator           (const CoilSelectManipulator &right);
    CoilSelectManipulator &operator=(const CoilSelectManipulator &right);

  friend class CoilSelectManipulatorUT;
};


/// Constructor
inline CoilSelectManipulator::CoilSelectManipulator(
  MrProtocolData::MrCoilSelectData&  rMrCoilSelectData
)
{
  m_pImpl = ICoilSelectManipulator::create();
  m_pImpl->init(rMrCoilSelectData);
}


/// Constructor
inline CoilSelectManipulator::CoilSelectManipulator(
  MrProtocolData::MrCoilSelectData&      rMrCoilSelectData,
  const MeasCoilContext& rCoilContext
)
{
  m_pImpl = ICoilSelectManipulator::create();
  m_pImpl->init(rMrCoilSelectData, rCoilContext);
}

/// Destructor
inline CoilSelectManipulator::~CoilSelectManipulator()
{
}


/// allows reuse of coil select manipulator
inline void CoilSelectManipulator::init(MrProtocolData::MrCoilSelectData&  rMrCoilSelectData,
                                        const MeasCoilContext&             rCoilContext)
{
  m_pImpl->init(rMrCoilSelectData, rCoilContext);
}


/// Fill the CoilSelectManipulator with all
/// relevant information from the protocol
inline void CoilSelectManipulator::setup (const MrProtocolData::MrProtData* pProt)
{
  m_pImpl->setup(pProt);
}

/// setCoilContext
inline void CoilSelectManipulator::setCoilContext (const MeasCoilContext& rCoilContext)
{
  m_pImpl->setCoilContext(rCoilContext);
}

/// getCoilContext
inline const MeasCoilContext& CoilSelectManipulator::getCoilContext () const
{
  return (m_pImpl->getCoilContext());
}

/// Set the nucleus in coilselect with given index
/// without modifying the CoilSelect.
/// The nucleus will be used with the manipulation.
inline bool CoilSelectManipulator::setRxNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus)
{
  return (m_pImpl->setRxNucleus(uiIndex, rNucleus));
}

/// Set the TX nucleus in coilselect with given index
/// without modifying the CoilSelect.
/// The nucleus will be used with the manipulation.
inline bool CoilSelectManipulator::setTxNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus)
{
  return (m_pImpl->setTxNucleus(uiIndex, rNucleus));
}

// Set the Tx and Rx nucleus in the coilselect with the given index
// without modifying the CoilSelect.
// The nucleus will be used with the manipulation.
inline bool CoilSelectManipulator::setNucleus(unsigned int uiIndex, const MeasNucleus& rNucleus)
{
  return (m_pImpl->setNucleus(uiIndex, rNucleus));
}

/// Set if combination allowed without modifying the
/// current coilselect.
inline void CoilSelectManipulator::setCombinationAllowed(bool bAllowed)
{
  m_pImpl->setCombinationAllowed(bAllowed);
}

/// Set the number ADC channels without modifying the coilselect.
inline void CoilSelectManipulator::setNumberOfADCChannels(int iNoOfADCChs)
{
  m_pImpl->setNumberOfADCChannels(iNoOfADCChs);
}

/// Set, if PAT is used without modifying the current
/// coilselect.
///
/// This flag will be used to dermine the system default
/// R-Factor.
inline void CoilSelectManipulator::setPAT(bool bPAT)
{
  m_pImpl->setPAT(bPAT);
}

/// Get the PAT acceleration the given coil element.
/// The element must be selected in coilselect with given index.
inline void CoilSelectManipulator::getPATAcceleration(
  unsigned int    uiIndex,
  AccCoilElement* pAccCoilElement,
  float&          rflAccelerationX,
  float&          rflAccelerationY,
  float&          rflAccelerationZ
) const
{
  m_pImpl->getPATAcceleration(uiIndex,
                              pAccCoilElement,
                              rflAccelerationX,
                              rflAccelerationY,
                              rflAccelerationZ);
}

/// Insert Body Coil in coilselect with given index,
/// if there is a body coil.
/// bAllElements is false (default): Insert only the first
/// body coil element (for normal handling enough)
/// bAllElements is true: Insert all body coil elements
/// (Required for tuning, service-SW)
inline CoilSelectManipulatorStatus CoilSelectManipulator::insertBodyCoil(
  unsigned int uiIndex,
  MrProtocolData::BCCMode      eBCCMode
)
{
  return (m_pImpl->insertBodyCoil(uiIndex, eBCCMode));
}

/// Insert all Pickup coil elements in coilselect with given index.
inline CoilSelectManipulatorStatus CoilSelectManipulator::insertPickupCoil(
  unsigned int uiIndex
)
{
  return (m_pImpl->insertPickupCoil(uiIndex));
}

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
inline CoilSelectManipulatorStatus CoilSelectManipulator::selectRxElement(
  unsigned int    uiIndex,
  AccCoilElement* pAccCoilElement,
  int             iRFactorElement
)
{
  return (m_pImpl->selectRxElement(uiIndex, pAccCoilElement, iRFactorElement));
} 


/// Select the given TX coil element in the TX coil select
/// with index uiIndex.
///
/// Return the status of the modification.
inline CoilSelectManipulatorStatus CoilSelectManipulator::selectTxElement(
  unsigned int    uiIndex,
  AccCoilElement* pAccCoilElement
)
{
  return (m_pImpl->selectTxElement(uiIndex, pAccCoilElement));
} 


/// Deselects all elements in coilselect with given index
/// which belongs to the PAT-Group of the given coil element.
///
/// Return the status of the modification.
inline CoilSelectManipulatorStatus CoilSelectManipulator::deselectRxElement(
  unsigned int    uiIndex,
  AccCoilElement* pAccCoilElement
)
{
  return (m_pImpl->deselectRxElement(uiIndex, pAccCoilElement));
}


/// Deselects TX element from TX coil select with index uiIndex
///
/// Return the status of the modification.
inline CoilSelectManipulatorStatus CoilSelectManipulator::deselectTxElement(
  unsigned int    uiIndex,
  AccCoilElement* pAccCoilElement
)
{
  return (m_pImpl->deselectTxElement(uiIndex, pAccCoilElement));
}



// Set SuppressMandatoryHandling in coilselect with given index
// without modifying the CoilSelect.
inline bool CoilSelectManipulator::setRxSuppressMandatoryProperties(
  unsigned int uiIndex,
  bool         bSuppressMandatoryProperties
)
{
  return (m_pImpl->setRxSuppressMandatoryProperties(uiIndex, bSuppressMandatoryProperties));
}


/// Returns true, if the given coil element is selected in
/// coilselect with given index.
inline bool CoilSelectManipulator::isRxSelected(
  unsigned int    uiIndex, 
  AccCoilElement* pAccCoilElement
) const
{
  return (m_pImpl->isRxSelected(uiIndex, pAccCoilElement));
}


/// Determine a possible CoilPlug object (PlugIDs) from a
/// given coilSelect.
///
/// All PlugIDs in the CoilPlug result should
/// provide all coil elements which are
/// selected in the coilSelect.
inline bool CoilSelectManipulator::determineCoilPlug(
  MrProtocolData::MrCoilPlugsData&          rMrCoilPlugData,
  const MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData
)
{
  return (m_pImpl->determineCoilPlug(rMrCoilPlugData, rMrRxCoilSelectData));
}

/// Determine a possible CoilPlug object (PlugIDs) from the
/// coilselect with given index.
/// All PlugIDs in the CoilPlug result should provide all
/// coil elements which are selected in the coilSelect.
inline bool CoilSelectManipulator::determineCoilPlug(unsigned int uiIndex)
{
  return (m_pImpl->determineCoilPlug(uiIndex));
}

/// Convert the coilselect with given index, determine coil plugs
/// and make it valid.
inline CoilSelectManipulatorStatus CoilSelectManipulator::adaptCoilPlugsAndSelect(unsigned int uiIndex)
{
  return (m_pImpl->adaptCoilPlugsAndSelect(uiIndex));
}


/// Return true, if no CoilPlug determine required.
inline bool CoilSelectManipulator::isCoilPlugValid(const MrProtocolData::MrCoilPlugsData& rMrCoilPlugData) const
{
  return (m_pImpl->isCoilPlugValid(rMrCoilPlugData));
}

/// If necessary, the coilselect with given index will be modified
/// to create a valid RX and TX coilselect.
///
/// If it is not possible to create a valid coilselect,
/// then the original coilselect will not be modified.
///
/// Return the status of the modification.
inline CoilSelectManipulatorStatus CoilSelectManipulator::buildValid(unsigned int uiIndex)
{
  return (m_pImpl->buildValid(uiIndex));
}


/// Deselect the given coil element and all other elements
/// with the same PAT-group.

/// If necessary, the coilselect with given index will be modified
/// to create a valid coilselect.
///
/// If it is not possible to create a valid coilselect,
/// then the original coilselect will not be modified.
///
/// Return the status of the modification.
inline CoilSelectManipulatorStatus CoilSelectManipulator::buildValidRx(unsigned int uiIndex)
{
  return (m_pImpl->buildValidRx(uiIndex));
}

/// Checks, if the coilselect with given index is valid. The
/// coilselect will never be changed. Only a check will be
/// performed.
///
/// Return the status of the check.
inline CoilSelectManipulatorStatus CoilSelectManipulator::check(unsigned int uiIndex)
{
  return (m_pImpl->check(uiIndex));
}

/// Checks, if the RX coilselect with given index is valid. The
/// coilselect will never be changed. Only a check will be
/// performed.
///
/// Return the status of the check.
inline CoilSelectManipulatorStatus CoilSelectManipulator::checkRx(unsigned int uiIndex)
{
  return (m_pImpl->checkRx(uiIndex));
}

/// Build a valid TX select
///
/// Return the status of the modification.
inline CoilSelectManipulatorStatus CoilSelectManipulator::buildValidTx(unsigned int uiIndex)
{
  return (m_pImpl->buildValidTx(uiIndex));
}

/// Checks, if the TX coilselect with given index is valid. The
/// coilselect will never be changed. Only a check will be
/// performed.
///
/// Return the status of the check.
inline CoilSelectManipulatorStatus CoilSelectManipulator::checkTx(unsigned int uiIndex)
{
  return (m_pImpl->checkTx(uiIndex));
}

/// Get RX coil select at uiIndex 
inline MrProtocolData::MrRxCoilSelectData* CoilSelectManipulator::getRxCoilSelect(unsigned int uiIndex)
{
  return (m_pImpl->getRxCoilSelect(uiIndex));
}

/// Todo (HJe): Comment???
inline bool CoilSelectManipulator::refreshRxSelect(
  unsigned int        uiIndex,
  MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData
)
{
  return (m_pImpl->refreshRxSelect(uiIndex, rMrRxCoilSelectData));
}


/// Gets the maximum usable R-Factor for the given coil
/// element.
///
/// It is not required, that the element is selected. It
/// must only be connected.
///
/// The value depends on the maximum supported R-Factor of
/// the element and the remaining number of receiver channels.
inline int CoilSelectManipulator::getMaxRFactor(
  MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData,
  AccCoilElement*     pAccCoilElement
) const
{
  return (m_pImpl->getMaxRFactor(rMrRxCoilSelectData, pAccCoilElement));
}


/// Get the maximum R factor for selected coil elements in coilselect
/// with given index.
inline int CoilSelectManipulator::getMaxRFactor(unsigned int uiIndex) const
{
  return (m_pImpl->getMaxRFactor(uiIndex));
}

/// Sets the R-Factor for the given selected coil element
/// in coilselect with given index.
///
/// If the element is not selected nothing will be changged.
///
/// Return the status of the modification.
inline CoilSelectManipulatorStatus CoilSelectManipulator::setRFactor(
  unsigned int    uiIndex,
  AccCoilElement* pAccCoilElement,
  int             iRFactorElement
)
{
  return (m_pImpl->setRFactor(uiIndex, pAccCoilElement, iRFactorElement));
}

/// Calls the internal method getRFactor() for getting the
/// R factor for the given coil element.
/// The element must be selected in coilselect with given index.
inline int CoilSelectManipulator::getRFactor(
  unsigned int    uiIndex, 
  AccCoilElement* pAccCoilElement
) const
{
  return (m_pImpl->getRFactor(uiIndex, pAccCoilElement));
}

/// Set the global R-Factor in coilselect with given index
/// and adjust the current selection of the elements with
/// this new factor.
///
/// If the input factor is 0, the system default R-Factor
/// will be used.
inline CoilSelectManipulatorStatus CoilSelectManipulator::setGlobalRFactor(
  unsigned int uiIndex,
  int          iRFactor
)
{
  return (m_pImpl->setGlobalRFactor(uiIndex, iRFactor));
}

/// Get the global R factor in coilselect with given index.
inline int CoilSelectManipulator::getGlobalRFactor(unsigned int uiIndex) const
{
  return (m_pImpl->getGlobalRFactor(uiIndex));
}

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
inline int CoilSelectManipulator::getSystemMaximumRFactor() const
{
  return (m_pImpl->getSystemMaximumRFactor());
}

/// Returns the state of the IndividualRFactors in the given coilSelect.
inline bool CoilSelectManipulator::isIndividualRFactors(const MrProtocolData::MrRxCoilSelectData& rMrRxCoilSelectData) const
{
  return (m_pImpl->isIndividualRFactors(rMrRxCoilSelectData));
}

/// Returns the state of the IndividualRFactors in coilselect
/// with given index
inline bool CoilSelectManipulator::isIndividualRFactors(unsigned int uiIndex) const
{
  return (m_pImpl->isIndividualRFactors(uiIndex));
}

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
inline const char * CoilSelectManipulator::getModeDesc(unsigned int uiIndex) const
{
  return (m_pImpl->getModeDesc(uiIndex));
}

///------------------------------------------------------------------------
/// misc. dump methods      
///------------------------------------------------------------------------

/// Dump all selected coil elements.
inline void CoilSelectManipulator::dump()
{
  m_pImpl->dump();
}

/// Dump all selected coil elements in coilselect with given index.
inline void CoilSelectManipulator::dump(unsigned int uiIndex)
{
  m_pImpl->dump(uiIndex);
}

/// Dump all selected coil elements in the given
/// MeasDumpEnv.
inline void CoilSelectManipulator::dump(MeasDumpEnv& rMeasDumpEnv)
{
  m_pImpl->dump(rMeasDumpEnv);
}

/// Dump all selected coil elements in coilselect with given index
/// to the given MeasDumpEnv.
inline void CoilSelectManipulator::dump(
  unsigned int uiIndex, 
  MeasDumpEnv& rMeasDumpEnv
)
{
  m_pImpl->dump(uiIndex, rMeasDumpEnv);
}

/// Dump the plugs
inline void CoilSelectManipulator::dumpPlugs(MeasDumpEnv& rMeasDumpEnv) const
{
  m_pImpl->dumpPlugs(rMeasDumpEnv);
}

/// Dump all rx connections of the coilselect with given index
inline void CoilSelectManipulator::dumpRxConn(std::ostream &rStream, unsigned int uiIndex) const
{
  m_pImpl->dumpRxConn(rStream, uiIndex);
}

/// Test methode to call many methodes of the CoilSelectManipulator
/// interactive about a SelectionMenu
inline void CoilSelectManipulator::edit()
{
  m_pImpl->edit();
}

/// Get the internal debug level
/// If the internal not set, get the debug level from the
/// measperm
inline int CoilSelectManipulator::getDebugLevel(void ) const
{
  return (m_pImpl->getDebugLevel());
}

/// Set the internal debug level
/// The debug level controls the trace outputs
inline void CoilSelectManipulator::setDebugLevel(int iDebugLevel)
{
  m_pImpl->setDebugLevel(iDebugLevel);
}

/// Check the RxChannels/ADCs of the RX coilselect with given index.
/// If all Rx channels are valid, do nothing.
/// If one or more Rx channels are invalid, determine a new one and set in into  
/// the Rx coilselect
/// Return the status of the check.
inline CoilSelectManipulatorStatus CoilSelectManipulator::checkAndFixRxChannels(unsigned int uiIndex)
{
  return m_pImpl->checkAndFixRxChannels(uiIndex);
}

/// Check the TxChannels of the TX coilselect with given index.
/// If all Tx channels are valid, do nothing.
/// If one or more Tx channels are invalid, determine a new one and set in into  
/// the Tx coilselect
/// Return the status of the check.
inline CoilSelectManipulatorStatus CoilSelectManipulator::checkAndFixTxChannels(unsigned int uiIndex)
{
  return m_pImpl->checkAndFixTxChannels(uiIndex);
}


#endif
