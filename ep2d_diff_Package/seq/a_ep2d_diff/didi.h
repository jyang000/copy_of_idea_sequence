/*! 
***************************************************************************
\file   didi.h

\brief  Interface of the class 'DiffusionDirections'
        (which controls the direction of diffusion gradients)

<b>Archive Information:</b>
\verbatim
   File-name: "h:\ep2d_diff_WIP\didi.h"
  Time-stamp: <12-October-2005 19:16 mizwa@EH402A9C>
Archive File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d_diff\didi.h
\endverbatim

\b Language: C++

\author PLM AW NEUR

\b Copyright: &copy; Siemens AG (http://www.siemensmedical.com/MR).
All rights reserved.   
This software may be only modified as long as the original author is credited 
in any subsequent revisions or modifications.
This software must not be sold or distributed as part of any commercial 
software package without the written permission of the author. 

***************************************************************************

\changed     30-Sep-2002; M.Zwanger; 4a21a
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description 
    - first version

\changed     20-Apr-2005; M.Zwanger; 4b13a
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description 
    - Definition of NormalisationMode enum changed due to name conflicts

\changed     12-Oct-2005; M.Zwanger; 4b13a
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description 
    - getGreatestNorm() and m_dGreatestNorm added

***************************************************************************
*/



#ifndef didi_h
/// multiple include protection
#define didi_h 1

#include <math.h>                                       // sqrt
#include "MrProtSrv/Domain/CoreNative/MrApplication.h"         // MrProtocolData::DiffDirCoordinateSystem

#ifdef BUILD_SEQU
  #define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"  // import/export control

/// Maximum number of directions (consistent with MrProt data structure)
/**
   This parameter determines the size of the look-up table.
   Changing its value is not critical generally, but affects both
   the memory consumption and the performance of the UI (incremental
   search over all directions).
*/
#define MAX_DIRECTIONS 10000

/// Default number of directions for FREE mode
/**
    This number is used as the default if the user selects the FREE diffusion
    mode. A corresponding diffusion vector set has to be provided within the
    ::prep method in order to guarantee that it's always possible to use this 
    number of diffusion directions (even if no corresponding vector set exists
    in the external file).
*/
#define FALLBACK_DIRECTIONS 6

/// These defines are used to assemble the diffusion vector set path information
/** "%CustomerSeq%\DiffusionVectorSets\*.dvs"
*/
#define VECTORFILEENV      "CustomerSeq"
#define VECTORFILEPATH     "DiffusionVectorSets"
#define VECTORFILEMASK     "*.dvs"


namespace SEQ_NAMESPACE
{

/// This structure holds the coordinates of a 3D vector.
struct VectorStruct 
{
    // Coordinates
    double dx;
    double dy;
    double dz;

    // Default constructor (does nothing)
    VectorStruct() : dx(0), dy(0), dz(0) {};
    // Alternative constructor (set coordinates)
    VectorStruct( double dVx, double dVy, double dVz ) : dx( dVx ), dy( dVy ), dz( dVz ) {};

	// Copy constructor
	VectorStruct(const VectorStruct& sOther) :
		dx(sOther.dx),
		dy(sOther.dy),
		dz(sOther.dz)
	{
	}

	// Assignment
	VectorStruct& operator=(const VectorStruct &rhs)
	{
		if (this != &rhs)
		{
			dx = rhs.dx;
			dy = rhs.dy;
			dz = rhs.dz;
		}
		return *this;
	}

    // Comparator used by std::sort
    inline bool operator<( const VectorStruct& sComp ) const
    {
        const double dVal1 =       dx *       dx +       dy *       dy +       dz *       dz;
        const double dVal2 = sComp.dx * sComp.dx + sComp.dy * sComp.dy + sComp.dz * sComp.dz;

        if ( dVal1 != dVal2 )
        {
            // Different vector length: sort by length
            return ( dVal1 < dVal2 );
        }
        else if ( dx != sComp.dx )
        {
            // Same length, different x-coordinate: sort by x-coordinate
            return ( dx < sComp.dx );
        }
        else if ( dy != sComp.dy )
        {
            // Same length and x-coordinate, different y-coordinate: sort by y-coordinate
            return ( dy < sComp.dy );
        }
        else if ( dz != sComp.dz )
        {
            // Same length and x- and y-coordinate, different z-coordinate: sort by z-coordinate
            return ( dz < sComp.dz );
        }

        // Identical vectors: ensure strict weak ordering
        return false;
    }

	// Dot (inner) product
	static double DotProduct(const VectorStruct &sVector1, const VectorStruct &sVector2)
	{
		return (sVector1.dx * sVector2.dx + sVector1.dy * sVector2.dy + sVector1.dz * sVector2.dz);
	}

	// Cross product
	static VectorStruct CrossProduct(const VectorStruct &sVector1, const VectorStruct &sVector2)
	{
		VectorStruct sVector;

		sVector.dx = sVector1.dy * sVector2.dz - sVector1.dz * sVector2.dy;
		sVector.dy = sVector1.dz * sVector2.dx - sVector1.dx * sVector2.dz;
		sVector.dz = sVector1.dx * sVector2.dy - sVector1.dy * sVector2.dx;

		return sVector;
	}
};


class __IMP_EXP DiffusionDirections
{
  
public:
    DiffusionDirections();  
    virtual ~DiffusionDirections();

    // ===========================================================================
    /// Dump the current vector table and the associated control information
    /**
        This function dumps all information about the current vector set,
        i.e. the number of directions, the requested coordinate system
        and all individual vectors.
    */
    virtual void dump ( void );

    // ===========================================================================
    /// Prepare a new internal diffusion vector set
    /**
        A new internal (=built-in) vector set will be generated which holds 
        the vectors for the specifed number of directions. The memory will be
        released and allocated dynamically.

        \return true if a vector set could be prepared successfully.

        \post   m_dGreatestNorm and m_dGreatestComponent will we set.
    */
    virtual bool prepInternal
        ( 
        long lDirections,                    /**< Imp: Number of diffusion directions                                            */ 
        char cFlag,                          /**< Imp: Modifier (if more than one set with the same number of directions exists) */
        bool bIsContextPrepForBinarySearch   /**< Imp: Suppress warnings                                                         */
        );

#ifdef WIN32
    // ===========================================================================
    /// Prepare a new external diffusion vector set
    /**
        A new vector set will be generated which holds the vectors 
        for the specifed number of directions. The memory will be
        released and allocated dynamically.

        If a vector set for the requested number of directions has been found in 
        the external vetor file during the scan at the instantiation of didi, 
        this external vector set will be loaded from file now.

        As fallback (e.g. if the name of the external vector file is not yet 
        known, an internal set of FALLBACK_DIRECTIONS diffusion directions is
        provided. Calling the method with this number will always yield a 
        valid vector set.

        This method is available on the host only.

        \return true if a vector set could be prepared successfully.

        \post m_dGreatestNorm and m_dGreatestComponent will we set.
    */
    virtual bool prepExternal
        ( 
        long lDirections,                    /**< Imp: Number of diffusion directions */          
        bool bIsContextPrepForBinarySearch   /**< Imp: Suppress warnings              */
        );
#endif // #ifdef WIN32

    // ===========================================================================
    /// Prepare a diffusion vector set using explicit parameters
    /**
        A new vector set will be generated which holds the vectors 
        for the specifed number of directions. The memory will be
        released and allocated dynamically.

        \return true if a vector set could be prepared successfully.

        \post m_dGreatestNorm and m_dGreatestComponent will we set.
    */
    virtual bool prepExplicit
        ( 
        long                                     lDirections,                   /**< Imp: Number of diffusion directions */  
        MrProtocolData::DiffDirCoordinateSystem  eCoordinateSystem,             /**< Imp: Coordinates system             */
        const std::string                       &strComment,                    /**< Imp: User comment (for tooltip)     */
        const std::vector<VectorStruct>         &vDiffDir,                      /**< Imp: Array of diffusion vectors     */
        bool                                     bIsContextPrepForBinarySearch  /**< Imp: Suppress warnings              */ 
        );

    // ===========================================================================
    /// Prepare a diffusion vector set covering a defined q-space region
    /**
        A new vector set will be generated which holds the vectors 
        covering a certain region in q-space region according to
        the input parameters. Usually, the vectors will exhibit
        different lengths (i.e. generate different diffusion weightings).

        \return Actual number of diffusion directions if successful,
                0 otherwise

        \post m_dGreatestNorm and m_dGreatestComponent will we set.
    */
    virtual long prepQSpace
        ( 
        SEQ::DiffQSpaceCoverageMode   eQSpaceCoverage,               /**< Imp: Q-space coverage (full / half)        */  
        SEQ::DiffQSpaceSamplingScheme eQSpaceSampling,               /**< Imp: Q-space sampling (cartesian / radial) */
        long                          lQSpaceSteps,                  /**< Imp: Number steps to cover in Q-space      */
        bool                          bIsContextPrepForBinarySearch  /**< Imp: Suppress warnings                     */ 
        );

    // ===========================================================================
    /// Get coordinate system of the current vector set (inline code)
    /**
        Valid after a successful preparation.

        \return Current coordinate system (prs or xyz).
    */
    virtual MrProtocolData::DiffDirCoordinateSystem getCoordinateSystem ( void ); 

    // ===========================================================================
    /// Get total number of directions of the current vector set (inline code)
    /**
        Valid after a successful preparation.

        \return Current number of diffusion directions.
    */
    virtual long getNumberOfDirections ( void );

    // ===========================================================================
    /// Get number of external directions sets (inline code)
    /**
        Valid after successfully setting the external vector filename.

        \return Number of direction sets in current external file.
    */
    virtual long getNumberOfExternalDirectionSets ( void );

    // ===========================================================================
    /// Set q-space hemisphere for partial coverage
    /**
        For partial q-space coverage (e.g. half), this parameter determines
        the area that should preferably get scanned completely.

        E.g.: sHemisphere = { 0.0, -1.0, 0.0 } => scan negative qy-coordinates,
                                                  omit positive qy-coordinates

        \return true if the provided vector was valid,
                false otherwise (default hemisphere is set in this case)
    */
    virtual bool setQSpaceHemisphere ( VectorStruct sHemisphere );

    // ===========================================================================
    /// Check availability of internal direction set (inline code)
    /**
        \return true if an internal direction set with the given number of directions exists.
    */
    virtual bool isDirectionInternal 
        ( 
            long lDirections                 /**< Imp: Number of diffusion directions */ 
        );

    // ===========================================================================
    /// Check availability of external direction set (inline code)
    /**
        Useable after successfully setting the external vector filename.

        \return true if an external direction set with the given number of directions exists.
    */
    virtual bool isDirectionExternal 
        ( 
            long lDirections                 /**< Imp: Number of diffusion directions */  
        );

    // ===========================================================================
    /// Check existence of directory containing external direction sets (inline code)
    /**
        On the scanner, this will always return false.

        \return true if the path defined internally exists.
    */
    virtual bool isExtDiffDirPresent( void );

    // ===========================================================================
    /// Get norm of the specifed gradient vector
    /** 
        Valid after a successful preparation.

        The used norm is \f$\sqrt{ x^2 + y^2 + z^2}\f$.

        \return Norm of the specified gradient vector. If the index is out of bounds, 
                the norm of the last vector in the current vector set will be returned.
    */
    virtual double getNorm 
        ( 
            long lIndex                     /**< Imp: Diffusion vector index */  
        ); 

    // ===========================================================================
    /// Get greatest norm of current vector set (inline code)
    /**
        Valid after a successful preparation.

        \return Greatest norm within current vector set.
    */
    virtual double getGreatestNorm ( void ); 

    // ===========================================================================
    /// Get greatest individual component (absolute value) of current vector set (inline code)
    /**
        Valid after a successful preparation.

        \return Greatest individual component within current vector set.
    */
    virtual double getGreatestComponent ( void );

    // ===========================================================================
    /// Get x component of specified gradient vector
    /** 
        Valid after a successful preparation.
        This inline function does not(!) check if 'lIndex' is out of range!

        \return x component of specified gradient vector.
    */
    virtual double getX 
        ( 
            long lIndex                     /**< Imp: Diffusion vector index */   
        );
    // ===========================================================================
    /// Get y component of specified gradient vector
    /** 
        Valid after a successful preparation.
        This inline function does not(!) check if 'lIndex' is out of range!

        \return y component of specified gradient vector.
    */
    virtual double getY 
        ( 
            long lIndex                     /**< Imp: Diffusion vector index */    
        );

    // ===========================================================================
    /// Get z component of specified gradient vector
    /** 
        Valid after a successful preparation.
        This inline function does not(!) check if 'lIndex' is out of range!

        \return z component of specified gradient vector.
    */
    virtual double getZ 
        ( 
            long lIndex                     /**< Imp: Diffusion vector index */     
        );

    // ===========================================================================
    /// Set diffusion vector filename
    /**
        The specified filename will be written into m_sVectorFileName. This file
        will get scanned for valid diffusion vector set definitions.

        On the scanner, this method is empty and always returns true.

        \return true if successful
    */
    virtual bool setVectorFileName
        ( 
            const char* pszVectorFileName   /**< Imp: Diffusion vector filename (excluding path information) */      
        );

    // ===========================================================================
    /// Get the current diffusion vector filename
    /**
        Valid after setting a vector filename.

        \return Current vector filename (excluding path information)
    */
    virtual std::string getVectorFileName( void );

    // ===========================================================================
    /// Get error message
    /** 
        Valid after setting the external vector filename without success.

        If reading a diffusion vector set fails, this provides a more
        specific error message.

        \return Error message.
    */
    virtual std::string getErrorMessage  ( void );

    // ===========================================================================
    /// Get comment for current diffusion vector set
    /**
        Can be used e.g. for a tooltip.

        Valid after a successful preparation.

        \return Comment for current diffusion vector set.
    */
    virtual std::string getComment       ( void );

#ifdef WIN32
    // ===========================================================================
    /// Writes the current free vector set to an external file.
    /**
        The free diffusion vector set currently selected will be written to
        the file %CustomerSeq%/pszVectorFileName.

        This method is available on the host only.

        \return true if the vector set has been saved successfully.

        \pre    The environment variable $CustomerSeq must be available.
    */
    virtual bool writeToFile  ( const char* pszVectorFileName );
#endif  // #ifdef WIN32

protected:

    // ===========================================================================
    /// Calculate greatest norm and greatest component
    /**
        Should be called after successful generation of an (internal or
        external) diffusion vector set.
    */
    virtual void calcGreatestNorm ( void );  

    // ===========================================================================
    /// Normalize a complete vector set
    /**
        The provided vector set will be normalized in situ
        as requested by the mode parameter. The used norm is \f$\sqrt{ x^2 + y^2 + z^2}\f$.

        \return true if everything was ok, false if a NULL vector in the vector set
                could not be normalized.
    */
    virtual bool normalize 
      ( 
        VectorStruct                         sVector[],         /**< Imp/Exp: Array containing diffusion vectors */
        const long                           lDirections,       /**< Imp:     Number of diffusion directions     */
        MrProtocolData::DiffDirNormalization eMode              /**< Imp:     Algorithm used for normalization   */
       );

    // ===========================================================================
    /// Helper function to set a gradient vector (inline code)
    /** 
        The specified coordinates will be written into m_sVector[lIndex].
        This inline function does not(!) check if 'lIndex' is out of range!
    */
    virtual void setVector 
        ( 
            long lIndex,            /**< Imp: Diffusion direction index    */ 
            double dx,              /**< Imp: Diffusion vector x-component */
            double dy,              /**< Imp: Diffusion vector x-component */
            double dz               /**< Imp: Diffusion vector x-component */
        );

    // ===========================================================================
    /// Initialize vector lookup table
    /** 
        Marks valid internal (host and scanner) and external (host only)
        diffusion vector sets.

        \return true if successful.
    */
    virtual bool initLookupTable ( void );

#ifdef WIN32
    // ===========================================================================
    /// Reads the specified vector set from an external file.
    /**
        The file %CustomerSeq%/m_sVectorFileName will be opened and the requested 
        vector set will be searched. If the section defining the requested number 
        of directions has been found, it will be read in into a dynamically 
        allocated temporary memory. If there are no errors, m_lDirections
        will be set to its new value and the temporary buffer will be copied
        into the m_sVector buffer.

        A special meaning has the value 0 which makes this method search for 
        all defined direction sections.

        This method is available on the host only.

        \return true if the requested vector set has been found
                and loaded successfully.

        \pre    The environment variable $CustomerSeq must be available.
    */
    virtual bool readFromFile 
        ( 
            long lDirections        /**< Imp: Vector set to be loaded (characterized by the number of directions). */
        );
#endif  // #ifdef WIN32

    /// Name of the external vector file
    std::string m_sVectorFileName;

    /// Contains error message if loading of diffusion vector set fails
    std::string m_sErrorMessage;

    /// Contains user defined comment for current diffusion vector set
    std::string m_sUserComment;

    /// Identifier of this class (used for tracing)
    const char* m_tIdent;
  
    /// Severity code of debug messages
    long m_lDebugLevel;

    /// Number of directions of current vector set
    long  m_lDirections;

    /// Number of external directions sets
    long  m_lExternalDirectionSets;

    /// Directional vector indicating prefered q-space hemisphere to scan
    VectorStruct m_sQSpaceHemisphere;

    /// Flag indicates that a directory containing external diffusion vector sets is present
    bool   m_bExtDiffDirPresent;

    /// Coordinate system of the current vector set
    MrProtocolData::DiffDirCoordinateSystem m_eCoordinateSystem;

    /// The greatest norm of all vectors
    double m_dGreatestNorm;

    /// The greatest component within all vectors
    double m_dGreatestComponent;

    /// Pointer to the current vector array
    VectorStruct *m_sVector;


    /// This look-up table stores which directions are available.
    /** 
        The encoding is as follows: If a direction i is not available,
        the value is 0. If a direction is available in the external definition
        file, the value LookupTable[i] will be 1. If the direction is
        built-in, the value is -1.

        The table is filled when calling the constructor.
        It will be updated if prep() fails due to a syntax error in 
        an individual section.      
    */ 
    // Index is +1 because we want to say "m_cLookupTable[MAX_DIRECTIONS]"
    signed char m_cLookupTable[MAX_DIRECTIONS+1];

private:

    /// \name Hidden Interfaces
    //@{
    /// copy constructor
    DiffusionDirections(const DiffusionDirections &right);
    /// assignment operator
    DiffusionDirections & operator=(const DiffusionDirections &right);
    //@}

};



// -----------------------------
// Inline function declarations:
// -----------------------------

/// Get coordinate system of the current vector set (inline code)
inline MrProtocolData::DiffDirCoordinateSystem DiffusionDirections::getCoordinateSystem ( void )  
{ 
    return m_eCoordinateSystem; 
}

/// Get total number of directions of the current vector set (inline code)
inline long DiffusionDirections::getNumberOfDirections ( void ) 
{ 
    return m_lDirections; 
}

/// Get number of external directions sets (inline code)
inline long DiffusionDirections::getNumberOfExternalDirectionSets ( void )  
{ 
    return m_lExternalDirectionSets; 
}

/// Check availability of internal direction set (inline code)
inline bool DiffusionDirections::isDirectionInternal ( long lIndex )
{ 
    return ( m_cLookupTable[lIndex] == -1 ); 
}

/// Check availability of external direction set (inline code)
inline bool DiffusionDirections::isDirectionExternal ( long lIndex )
{ 
    return ( m_cLookupTable[lIndex] == +1 ); 
}

/// Check existence of directory containing external direction sets (inline code)
inline bool DiffusionDirections::isExtDiffDirPresent( void ) 
{ 
    return m_bExtDiffDirPresent; 
}

/// Get greatest norm of current vector set (inline code)
inline double DiffusionDirections::getGreatestNorm ( void )  
{ 
    return m_dGreatestNorm; 
}

/// Get greatest individual component (absolute value) of current vector set (inline code)
inline double DiffusionDirections::getGreatestComponent ( void )  
{ 
    return m_dGreatestComponent;
}

/// Get x component of specified gradient vector
inline double DiffusionDirections::getX ( long lIndex ) 
{ 
    return m_sVector[lIndex].dx;
}

/// Get y component of specified gradient vector
inline double DiffusionDirections::getY ( long lIndex )
{ 
    return m_sVector[lIndex].dy;
}

/// Get z component of specified gradient vector
inline double DiffusionDirections::getZ ( long lIndex ) 
{ 
    return m_sVector[lIndex].dz;
}

/// Helper function to set a gradient vector (inline code)
inline void DiffusionDirections::setVector ( long lIndex, double dx, double dy, double dz ) 
{
    m_sVector[lIndex].dx = dx;  
    m_sVector[lIndex].dy = dy;   
    m_sVector[lIndex].dz = dz;
}

/// Set diffusion vector filename
inline bool DiffusionDirections::setVectorFileName 
(
#ifndef WIN32
    const char* /* pszVectorFileName */ 
) 
{
    return true;
}
#else
    const char* pszVectorFileName
) 
{
    // Assign new filename
    m_sVectorFileName.assign( pszVectorFileName );

    // Update m_cLookupTable
    return initLookupTable();
}
#endif  // #ifndef WIN32

/// Get the current diffusion vector filename
inline std::string DiffusionDirections::getVectorFileName  ( void )
{
    return m_sVectorFileName;
}

/// Get error message
inline std::string DiffusionDirections::getErrorMessage  ( void )
{
    return m_sErrorMessage;
}

/// Get comment for current diffusion vector set
inline std::string DiffusionDirections::getComment  ( void )
{
    return m_sUserComment;
}

}//end of namespace SEQ_NAMESPACE
#endif     // of ifndef didi_h
