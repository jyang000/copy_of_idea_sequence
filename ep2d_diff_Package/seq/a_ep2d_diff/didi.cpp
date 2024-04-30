/**
***************************************************************************
\file   didi.cpp

\brief  Implementation of the class 'DiffusionDirections'
(which controls the direction of diffusion gradients)

<b>Archive Information:</b>
\verbatim
File-name: "h:\ep2d_diff_WIP\didi.cpp"
Time-stamp: <12-October-2005 19:13 mizwa@EH402A9C>
Archive File: \n4\pkg\MrServers\MrImaging\seq\a_ep2d_diff\didi.cpp
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

\changed     30-Sep-2002; M.Zwanger; 4a21a; CHARM: n.a.
\requirement N4_elh_DIFF_DiffusionTensorImaging
\description
- First version created


\changed    07-Sep-2005; M.Zwanger; 4b13a; CHARM: n.a.
\description
- Vector set for 256 directions added (MRM-Yellow request)
- The old vector set for 12 directions has been replaced by an octahedron.


\changed    06-Oct-2005; M.Zwanger; 4b13a; CHARM: n.a.
\description
- Avoid exception if $CustomSeq has not been set


\changed    12-Oct-2005; M.Zwanger; 4b13a; Customer wish
\description
- The vector file format now also supports the definition of null vector
lines at any arbitrary position for all normalisation modes.


\changed    09-Dec-2005; S.Huwer; 4b13a; Bug Fix
\description
- Calculate m_GreatestNorm in prep() also for external tables

***************************************************************************
*/


#include <algorithm>
#include <cstdio>
#include <ctime>
#define _USE_MATH_DEFINES
#include <cmath>

#include "MrImaging/seq/a_ep2d_diff/didi.h"
#include "MrImaging/seq/SeqDebug.h"           // mPrintTrace
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

namespace SEQ_NAMESPACE
{

    /// Internal Epsilon parameter
    const double dEpsilon = 1e-6;

    /// Internal maximum length of user comment
    const long   lMaxCommentLength = 255;

    /// Test if two doubles have the same value
    /**
        This function tests if the two doubles passed as parameters are identical
        with respect to a epsilon defined within this function.

        This function is a re-implementation of fGSLAlmEqual coded in libGSL.
        It has been reimplemented in order to avoid about 850 .h files to be
        included.
        */
    bool fAlmostEqual(double da, double db)
    {
        if((da - dEpsilon < db) && (da + dEpsilon > db))
        {
            return true;
        }

        return false;
    }

}//end of namespace SEQ_NAMESPACE


#ifndef SEQ_NAMESPACE
#pragma message ("NOTE: SEQ_NAMESPACE not defined");
#endif
using namespace SEQ_NAMESPACE;





// ***************************************************************************
// class DiffusionDirections
// ***************************************************************************

// ===========================================================================
/*!
\class DiffusionDirections

\author PLM AW NEUR

\brief This class controls the direction of the diffusion gradient vectors.

Internally a dynamic array is created, which holds the
actual number of diffusion vectors. This array can be accessed by the (protected)
pointer ::m_sVector.

An individual vector can be retrieved by addressing its index.
For this purpose the inline functions ::getX, ::getY and ::getZ
are used.

Each vector is represented as a triple of numbers.
The member ::m_eCoordinateSystem specifies which coordinate systems
these triples refer to.

There are two different ways to build up a direction table:
- Additional vector tables can be loaded from an external input file
- Built-in directions are "hard-coded" in the ::prepInternal method.
They are described in detail in section \ref Builtin of this manual.

The external definition file is searched for. This way, it is possible for
the user to supersede the internal definition. The syntax of this external
file is described in section "\ref VectorFile".

To speed up things a little bit, available directions are cached in an
array (\ref m_cLookupTable). For all directions not found prep() will
immediately return false.

\b Environment:

\e DEBUG_DiffusionDirections: This variable controls the debug level.
Its value must be in the range 0..255 and will be stored in the member
variable m_lDebugLevel. Default is 0.

*/
// ===========================================================================


// ===========================================================================
/// The constructor initializes some member variables.
/** In detail, the contructor initializes the m_lDebugLevel member from the registry,
    initializes m_sVectorFileName and checks the availability of a directory
    containing external diffusion vector sets.
    Available internal directions will be marked in the m_cLookupTable array.
    */
DiffusionDirections::DiffusionDirections()
// ===========================================================================
: m_sVectorFileName(""),
m_sErrorMessage(""),
m_sUserComment(""),
m_tIdent("DiffusionDirections"),
m_lDebugLevel(0),
m_lDirections(1),
m_lExternalDirectionSets(0),
m_sQSpaceHemisphere(0.0, 0.0, 1.0),                       // Magnitude 1.0 required. Selected default: omit negative q-space coordinates
m_bExtDiffDirPresent(false),
m_eCoordinateSystem(MrProtocolData::DIFFDIR_CS_XYZ),
m_dGreatestNorm(0.0),
m_dGreatestComponent(0.0),
m_sVector(NULL)
{
    if(getenv("DEBUG_DiffusionDirections"))
    {
        m_lDebugLevel = atoi(getenv("DEBUG_DiffusionDirections"));
    }

    // dummy initialization to allow delete operation in prep()
    m_sVector = new VectorStruct[1];
    setVector(0, 1.0, 0, 0);

#ifdef WIN32
    // --------------------------------------------------------------------------------
    // Check if directory containing external diffusion vector set (.dvs) files exists
    // --------------------------------------------------------------------------------
    std::stringstream sVectorFile;
    char* ptCustomerSeq = NULL;
    ptCustomerSeq = getenv(VECTORFILEENV);

    if(ptCustomerSeq == NULL)
    {
        sVectorFile << "ThisFileDoesNotExist";
    }
    else
    {
        // "%CustomerSeq%/DiffusionVectorSets"
        sVectorFile << ptCustomerSeq << "/" << VECTORFILEPATH;
    }

    struct stat sStatus;
    m_bExtDiffDirPresent = (stat(sVectorFile.str().c_str(), &sStatus) == 0);
#endif  // #ifdef WIN32

    // ---------------------------
    // Initialize m_cLookupTable
    // ---------------------------
    initLookupTable();
}


// ===========================================================================
/// The destructor releases the allocated memory of the vector array.
DiffusionDirections::~DiffusionDirections()
// ===========================================================================

{
#ifdef DEBUG_DIFFUSION  
    cout <<  m_tIdent << "::Destructor" << endl;
#endif // #ifdef DEBUG_DIFFUSION
    delete[] m_sVector;
}


// ===========================================================================
/// Dump the current vector table and the associated control information
/**
    This function dumps all information about the current vector set,
    i.e. the number of directions, the requested coordinate system
    and all individual vectors.
    */
void DiffusionDirections::dump(void)
// ===========================================================================

{
    long lIndex = 0;

    if(m_eCoordinateSystem == MrProtocolData::DIFFDIR_CS_XYZ)
    {
        SEQ_TRACE_INFO.print("%ld vectors specified in XYZ coordinates", m_lDirections);
    }
    else
    {
        SEQ_TRACE_INFO.print("%ld vectors specified in PRS coordinates", m_lDirections);
    }

    for(lIndex = 0; lIndex < m_lDirections; lIndex++)
    {
        SEQ_TRACE_INFO.print("Vector[%ld] = ( %6.3f, %6.3f, %6.3f )     Norm = %6.3f",
                   lIndex, m_sVector[lIndex].dx, m_sVector[lIndex].dy, m_sVector[lIndex].dz, getNorm(lIndex));
    }

}


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
bool DiffusionDirections::setQSpaceHemisphere(VectorStruct sHemisphere)
{
    m_sQSpaceHemisphere = sHemisphere;

    // Verify that at least one component is non-zero
    if(fAlmostEqual(m_sQSpaceHemisphere.dx, 0.0) && fAlmostEqual(m_sQSpaceHemisphere.dy, 0.0) && fAlmostEqual(m_sQSpaceHemisphere.dz, 0.0))
    {
        // Fallback solution: ensure proper functioning of prepQSpace
        m_sQSpaceHemisphere.dx = 0.;
        m_sQSpaceHemisphere.dy = 0.;
        m_sQSpaceHemisphere.dz = 1.;
        return false;
    }

    // Normalization
    double dMagn = sqrt(m_sQSpaceHemisphere.dx * m_sQSpaceHemisphere.dx + m_sQSpaceHemisphere.dy * m_sQSpaceHemisphere.dy + m_sQSpaceHemisphere.dz * m_sQSpaceHemisphere.dz);

    // Store hemisphere vector
    m_sQSpaceHemisphere.dx = m_sQSpaceHemisphere.dx / dMagn;
    m_sQSpaceHemisphere.dy = m_sQSpaceHemisphere.dy / dMagn;
    m_sQSpaceHemisphere.dz = m_sQSpaceHemisphere.dz / dMagn;

    return true;
}

// ===========================================================================
/// Get norm of the specifed gradient vector
/**
    Valid after a successful preparation.

    The used norm is \f$\sqrt{ x^2 + y^2 + z^2}\f$.

    \return Norm of the specified gradient vector. If the index is out of bounds,
    the norm of the last vector in the current vector set will be returned.
    */
double DiffusionDirections::getNorm
(
long lIndex                     /**< Imp: Diffusion vector index */
)
// ===========================================================================

{
    if(lIndex >= m_lDirections)
    {
        SEQ_TRACE_INFO.print("WARNING: Index %ld out of bounds %ld", lIndex, m_lDirections);
        lIndex = m_lDirections - 1;
    }

    return (sqrt(m_sVector[lIndex].dx * m_sVector[lIndex].dx +
        m_sVector[lIndex].dy * m_sVector[lIndex].dy +
        m_sVector[lIndex].dz * m_sVector[lIndex].dz));
}


// ===========================================================================
/// Initialize vector lookup table
/**
    Marks valid internal (host and scanner) and external (host only)
    diffusion vector sets.

    Convention used:
    0:  no direction set for specified number of directions
    -1:  built-in direction set
    +1:  user defined direction set

    \return true if successful.
    */
bool DiffusionDirections::initLookupTable(void)
// ===========================================================================
{
    for(int iLoop = 0; iLoop <= MAX_DIRECTIONS; iLoop++)
    {
        m_cLookupTable[iLoop] = 0;
    }

    // Flag built-in directions
    m_cLookupTable[0]   = -1;             // Adjustment mode (dynamic distortion correction)
    m_cLookupTable[1]   = -1;
    m_cLookupTable[3]   = -1;
    m_cLookupTable[4]   = -1;
    m_cLookupTable[6]   = -1;
    m_cLookupTable[10]  = -1;
    m_cLookupTable[12]  = -1;
    m_cLookupTable[20]  = -1;
    m_cLookupTable[30]  = -1;
    m_cLookupTable[64]  = -1;
    m_cLookupTable[256] = -1;

#ifdef WIN32
    // Scan external vector file for user defined directions
    if(m_bExtDiffDirPresent)
    {
        // Flag built-in fallback directions
        m_cLookupTable[FALLBACK_DIRECTIONS] = -1;

        return readFromFile(0);
    }
#endif // #ifdef WIN32

    return true;
}

// ===========================================================================
/// Prepare a new internal diffusion vector set
/**
    A new internal (=built-in) vector set will be generated which holds
    the vectors for the specifed number of directions. The memory will be
    released and allocated dynamically.

    If no internal table is found (i.e. an invalid number of lDirections),
    the function will return false.

    When fSeqPrep() runs in "binary search" mode for the "Directions" parameter,
    this is actually no binary search, but a point-by-point query.
    (There is no contingeous area of valid parameter values.)
    Therefore the performance of this function is quiet crucial for
    the UI behaviour.

    \return true if a vector set could be prepared successfully.

    \post   m_dGreatestNorm and m_dGreatestComponent will we set.
    */
bool DiffusionDirections::prepInternal
(
long lDirections,                   /**< Imp: Number of diffusion directions                                            */
char cFlag,                         /**< Imp: Modifier (if more than one set with the same number of directions exists) */
bool bIsContextPrepForBinarySearch  /**< Imp: Suppress warnings                                                         */
)
{
    if(lDirections < 0)
    {
        if(!bIsContextPrepForBinarySearch || (m_lDebugLevel & DEBUG_RETURN))
        {
            SEQ_TRACE_INFO.print("ERROR: lDirections=%ld is invalid.", lDirections);
        }
        return false;
    }

    if(m_cLookupTable[lDirections] == 0)
    {
        // A vector set for this direction number does not exist
        if(!bIsContextPrepForBinarySearch || (m_lDebugLevel & DEBUG_RETURN))
        {
            SEQ_TRACE_INFO.print("ERROR: lDirections=%ld is invalid.", lDirections);
        }
        return false;
    }

    switch(lDirections)
    {
        /** \note Take care that directions added to this function are also flagged
        in m_cLookupTable when executing the contructor.
        */
        case 0:
            // Adjustment mode
            m_lDirections       = 9;          // +x, 0, -x, +y, 0, -y, +z, 0, -z
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            m_sUserComment      = "SIEMENS\nInternal adjustment vector set";
            setVector(0, 1.0, 0.0, 0.0);
            setVector(1, 0.0, 0.0, 0.0);
            setVector(2, -1.0, 0.0, 0.0);
            setVector(3, 0.0, 1.0, 0.0);
            setVector(4, 0.0, 0.0, 0.0);
            setVector(5, 0.0, -1.0, 0.0);
            setVector(6, 0.0, 0.0, 1.0);
            setVector(7, 0.0, 0.0, 0.0);
            setVector(8, 0.0, 0.0, -1.0);
            break;

        case 1:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_PRS;
            switch(cFlag)
            {
                case 'P':
                    setVector(0, 1.0, 0, 0);
                    m_sUserComment    = "SIEMENS\nInternal 'phase' direction set";
                    break;
                case 'R':
                    setVector(0, 0, 1.0, 0);
                    m_sUserComment    = "SIEMENS\nInternal 'read' direction set";
                    break;
                case 'S':
                    setVector(0, 0, 0, 1.0);
                    m_sUserComment    = "SIEMENS\nInternal 'slice' direction set";
                    break;
                case 'D':
                    setVector(0, 1.0, 1.0, 1.0);
                    m_sUserComment    = "SIEMENS\nInternal 'diagonal' direction set";
                    m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
                    break;
                default:
                    if(! bIsContextPrepForBinarySearch ||  (m_lDebugLevel & DEBUG_RETURN))
                    {
                        SEQ_TRACE_INFO.print("ERROR: Invalid control flag '%c' in sequence code.", cFlag);
                    }
                    return false;
            }
            break;

        case 3:
            m_lDirections = lDirections;
            delete[] m_sVector;
            m_sVector     = new VectorStruct[m_lDirections];
            switch(cFlag)
            {
                case 'O':
                    m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_PRS;
                    m_sUserComment      = "SIEMENS\nInternal 'orthogonal' direction set";
                    setVector(0, 1.0, 0.0, 0.0);
                    setVector(1, 0.0, 1.0, 0.0);
                    setVector(2, 0.0, 0.0, 1.0);
                    break;
                case 'T':
                    m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
                    m_sUserComment      = "SIEMENS\nInternal '3-scan-trace' direction set";
                    setVector(0, 1.0, 1.0, -0.5);
                    setVector(1, 1.0, -0.5, 1.0);
                    setVector(2, -0.5, 1.0, 1.0);
                    break;
                default:
                    if(! bIsContextPrepForBinarySearch ||  (m_lDebugLevel & DEBUG_RETURN))
                    {
                        SEQ_TRACE_INFO.print("ERROR: Invalid control flag '%c' in sequence code.", cFlag);
                    }
                    return false;
            }
            break;

        case 4:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            switch(cFlag)
            {
                case 'T':
                    m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
                    m_sUserComment      = "SIEMENS\nInternal '4-scan-trace' direction set";
                    setVector(0, 1.0, 1.0, -1.0);
                    setVector(1, 1.0, -1.0, 1.0);
                    setVector(2, -1.0, 1.0, 1.0);
                    setVector(3, -1.0, -1.0, -1.0);
                    break;
                default:
                    if(! bIsContextPrepForBinarySearch ||  (m_lDebugLevel & DEBUG_RETURN))
                    {
                        SEQ_TRACE_INFO.print("ERROR: Invalid control flag '%c' in sequence code.", cFlag);
                    }
                    return false;
            }
            break;

        case 6:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            m_sUserComment      = "SIEMENS\nInternal 6-directions set";
            // (Note: do not modify the following comments - doxygen greps for them)
            // vector set for 6 directions (normalised to maximum)
            setVector(0, 1.0, 0.0, 1.0);
            setVector(1, -1.0, 0.0, 1.0);
            setVector(2, 0.0, 1.0, 1.0);
            setVector(3, 0.0, 1.0, -1.0);
            setVector(4, 1.0, 1.0, 0.0);
            setVector(5, -1.0, 1.0, 0.0);
            break;

        case 10:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            m_sUserComment      = "SIEMENS\nInternal 10-directions set";
            // vector set for 10 directions (normalised to maximum)
            setVector(0, 0.000000, 0.809017, 0.618034);
            setVector(1, 0.000000, 0.190983, 1.000000);
            setVector(2, -0.587785, 0.809017, 0.190983);
            setVector(3, -0.951057, 0.190983, 0.309017);
            setVector(4, -0.363271, 0.809017, -0.500000);
            setVector(5, -0.587785, 0.190983, -0.809017);
            setVector(6, 0.363271, 0.809017, -0.500000);
            setVector(7, 0.587785, 0.190983, -0.809017);
            setVector(8, 0.587785, 0.809017, 0.190983);
            setVector(9, 0.951057, 0.190983, 0.309017);
            break;

        case 12:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            m_sUserComment      = "SIEMENS\nInternal 12-directions set";
            // vector set for 12 directions (normalised to maximum)
            setVector(0, 1.000000, 0.414250, -0.414250);
            setVector(1, 1.000000, -0.414250, -0.414250);
            setVector(2, 1.000000, -0.414250, 0.414250);
            setVector(3, 1.000000, 0.414250, 0.414250);
            setVector(4, 0.414250, 0.414250, 1.000000);
            setVector(5, 0.414250, 1.000000, 0.414250);
            setVector(6, 0.414250, 1.000000, -0.414250);
            setVector(7, 0.414250, 0.414250, -1.000000);
            setVector(8, 0.414250, -0.414250, -1.000000);
            setVector(9, 0.414250, -1.000000, -0.414250);
            setVector(10, 0.414250, -1.000000, 0.414250);
            setVector(11, 0.414250, -0.414250, 1.000000);
            break;

        case 20:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            m_sUserComment      = "SIEMENS\nInternal 20-directions set";
            // vector set for 20 directions (normalised to maximum)
            setVector(0, 1.000000, 0.000000, 0.000000);
            setVector(1, 0.000000, 1.000000, 0.000000);
            setVector(2, -0.031984, 0.799591, 0.599693);
            setVector(3, 0.856706, 0.493831, -0.148949);
            setVector(4, 0.834429, 0.309159, 0.456234);
            setVector(5, 0.834429, -0.309159, 0.456234);
            setVector(6, 0.856706, -0.493831, -0.148949);
            setVector(7, 0.822228, 0.000000, -0.569158);
            setVector(8, 0.550834, 0.425872, -0.717784);
            setVector(9, 0.468173, 0.834308, -0.291108);
            setVector(10, 0.515933, 0.808894, 0.281963);
            setVector(11, 0.391890, 0.515855, 0.761785);
            setVector(12, 0.478151, 0.000000, 0.878278);
            setVector(13, 0.391890, -0.515855, 0.761785);
            setVector(14, 0.515933, -0.808894, 0.281963);
            setVector(15, 0.468173, -0.834308, -0.291108);
            setVector(16, 0.550834, -0.425872, -0.717784);
            setVector(17, 0.111012, -0.264029, -0.958105);
            setVector(18, 0.111012, 0.264029, -0.958105);
            setVector(19, 0.031984, 0.799591, -0.599693);
            break;

        case 30:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            m_sUserComment      = "SIEMENS\nInternal 30-directions set";
            // vector set for 30 directions (normalised to maximum)
            setVector(0, -0.208098, 0.525514, 0.850005);
            setVector(1, 0.202387, 0.526131, 0.851002);
            setVector(2, 0.409956, 0.175267, 0.918257);
            setVector(3, -0.412630, 0.742620, 0.565889);
            setVector(4, -0.207127, 0.959492, 0.280092);
            setVector(5, -0.872713, 0.525505, 0.064764);
            setVector(6, -0.746815, 0.526129, 0.455449);
            setVector(7, -0.415238, 0.175473, 0.915841);
            setVector(8, -0.746636, 0.175268, 0.673642);
            setVector(9, -0.665701, 0.742619, -0.217574);
            setVector(10, -0.330391, 0.959489, -0.110458);
            setVector(11, -0.331275, 0.525513, -0.809983);
            setVector(12, -0.663936, 0.526130, -0.569521);
            setVector(13, -0.999332, 0.175474, -0.111904);
            setVector(14, -0.871398, 0.175267, -0.501922);
            setVector(15, 0.001214, 0.742616, -0.700356);
            setVector(16, 0.002949, 0.959483, -0.348370);
            setVector(17, 0.667975, 0.525509, -0.565356);
            setVector(18, 0.336490, 0.526126, -0.807431);
            setVector(19, 0.202383, -0.175470, 0.985002);
            setVector(20, 0.208094, 0.175265, -0.983848);
            setVector(21, 0.666452, 0.742619, -0.215262);
            setVector(22, 0.332212, 0.959489, -0.104850);
            setVector(23, 0.205064, 0.958364, 0.285421);
            setVector(24, 0.412630, 0.742620, 0.565889);
            setVector(25, 0.746093, 0.175315, 0.674232);
            setVector(26, 0.744110, 0.525505, 0.460568);
            setVector(27, 0.871894, 0.526125, 0.070507);
            setVector(28, 0.874264, 0.175471, -0.496841);
            setVector(29, 1.000000, 0.175267, -0.106112);
            break;


        case 64:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            m_sUserComment      = "SIEMENS\nInternal 64-directions set";
            // vector set for 64 directions (normalised to maximum)
            setVector(0, 1.000000, 0.000000, 0.000000);
            setVector(1, 0.000000, 1.000000, 0.000000);
            setVector(2, -0.026007, 0.649170, 0.760199);
            setVector(3, 0.591136, -0.766176, 0.252058);
            setVector(4, -0.236071, -0.524158, 0.818247);
            setVector(5, -0.893021, -0.259006, 0.368008);
            setVector(6, 0.796184, 0.129030, 0.591137);
            setVector(7, 0.233964, 0.929855, 0.283956);
            setVector(8, 0.935686, 0.139953, 0.323891);
            setVector(9, 0.505827, -0.844710, -0.174940);
            setVector(10, 0.346220, -0.847539, -0.402256);
            setVector(11, 0.456968, -0.630956, -0.626956);
            setVector(12, -0.486997, -0.388997, 0.781995);
            setVector(13, -0.617845, 0.672831, 0.406898);
            setVector(14, -0.576984, -0.104997, -0.809978);
            setVector(15, -0.826695, -0.520808, 0.212921);
            setVector(16, 0.893712, -0.039987, -0.446856);
            setVector(17, 0.290101, -0.541189, -0.789276);
            setVector(18, 0.115951, -0.962591, -0.244896);
            setVector(19, -0.800182, 0.403092, -0.444101);
            setVector(20, 0.513981, 0.839970, 0.173994);
            setVector(21, -0.788548, 0.152912, -0.595659);
            setVector(22, 0.949280, -0.233069, 0.211062);
            setVector(23, 0.232964, 0.782880, 0.576911);
            setVector(24, -0.020999, -0.187990, -0.981946);
            setVector(25, 0.216932, -0.955701, 0.198938);
            setVector(26, 0.774003, -0.604002, 0.190001);
            setVector(27, -0.160928, 0.355840, 0.920587);
            setVector(28, -0.147035, 0.731173, -0.666158);
            setVector(29, 0.888141, 0.417066, 0.193031);
            setVector(30, -0.561971, 0.231988, -0.793959);
            setVector(31, -0.380809, 0.142928, 0.913541);
            setVector(32, -0.306000, -0.199000, -0.931001);
            setVector(33, -0.332086, -0.130034, 0.934243);
            setVector(34, -0.963226, -0.265062, 0.044010);
            setVector(35, -0.959501, 0.205107, 0.193101);
            setVector(36, 0.452965, -0.888932, 0.067995);
            setVector(37, -0.773133, 0.628108, 0.088015);
            setVector(38, 0.709082, 0.408047, 0.575066);
            setVector(39, -0.692769, 0.023992, 0.720760);
            setVector(40, 0.681659, 0.528735, -0.505747);
            setVector(41, -0.141995, -0.724976, 0.673978);
            setVector(42, -0.740168, 0.388088, 0.549125);
            setVector(43, -0.103006, 0.822044, 0.560030);
            setVector(44, 0.584037, -0.596038, 0.551035);
            setVector(45, -0.088008, -0.335031, 0.938088);
            setVector(46, -0.552263, -0.792377, 0.259123);
            setVector(47, 0.838158, -0.458086, -0.296056);
            setVector(48, 0.362995, -0.560993, 0.743990);
            setVector(49, -0.184062, 0.392133, -0.901306);
            setVector(50, -0.720938, -0.692941, 0.008999);
            setVector(51, 0.433101, 0.682159, -0.589137);
            setVector(52, 0.502114, 0.690157, 0.521119);
            setVector(53, -0.170944, -0.508833, -0.843722);
            setVector(54, 0.462968, 0.422971, 0.778946);
            setVector(55, 0.385030, -0.809064, 0.444035);
            setVector(56, -0.713102, -0.247035, 0.656094);
            setVector(57, 0.259923, 0.884737, -0.386885);
            setVector(58, 0.001000, 0.077002, -0.997030);
            setVector(59, 0.037002, -0.902057, 0.430027);
            setVector(60, 0.570320, -0.303170, -0.763428);
            setVector(61, -0.282105, 0.145054, -0.948354);
            setVector(62, 0.721098, 0.608082, 0.332045);
            setVector(63, 0.266985, 0.959945, -0.084995);
            break;

        case 256:
            m_lDirections       = lDirections;
            delete[] m_sVector;
            m_sVector           = new VectorStruct[m_lDirections];
            m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            m_sUserComment      = "SIEMENS\nInternal 256-directions set";
            // vector set for 256 directions (normalised to maximum)
            setVector(0, 0.059010, -0.013002, 0.998173);
            setVector(1, 0.070008, 0.157018, 0.985111);
            setVector(2, -0.078038, 0.083041, 0.993486);
            setVector(3, -0.092987, -0.086988, 0.991860);
            setVector(4, 0.046987, -0.178952, 0.982735);
            setVector(5, 0.199944, -0.108970, 0.973729);
            setVector(6, 0.211021, 0.052005, 0.976097);
            setVector(7, 0.210999, 0.212999, 0.953997);
            setVector(8, -0.050982, 0.250913, 0.966666);
            setVector(9, -0.208056, 0.186050, 0.960259);
            setVector(10, -0.355044, 0.120015, 0.927114);
            setVector(11, -0.225050, 0.016004, 0.974216);
            setVector(12, -0.239998, -0.156999, 0.957994);
            setVector(13, -0.107055, -0.253129, 0.961491);
            setVector(14, 0.031988, -0.337873, 0.940648);
            setVector(15, 0.191053, -0.268075, 0.944264);
            setVector(16, 0.333011, -0.202007, 0.921030);
            setVector(17, 0.350043, -0.034004, 0.936116);
            setVector(18, 0.353962, 0.131986, 0.925900);
            setVector(19, 0.321087, 0.303082, 0.897243);
            setVector(20, 0.124035, 0.337095, 0.933264);
            setVector(21, -0.028007, 0.395092, 0.918214);
            setVector(22, -0.189921, 0.343856, 0.919616);
            setVector(23, -0.340929, 0.281942, 0.896814);
            setVector(24, -0.483028, 0.214012, 0.849048);
            setVector(25, -0.492017, 0.052002, 0.869031);
            setVector(26, -0.364886, -0.049984, 0.929710);
            setVector(27, -0.376988, -0.226993, 0.897972);
            setVector(28, -0.251019, -0.327024, 0.911068);
            setVector(29, -0.116041, -0.415148, 0.902323);
            setVector(30, 0.027007, -0.490130, 0.871231);
            setVector(31, 0.176089, -0.413208, 0.893450);
            setVector(32, 0.323985, -0.380982, 0.865960);
            setVector(33, 0.451026, -0.272016, 0.850049);
            setVector(34, 0.477203, -0.106045, 0.872371);
            setVector(35, 0.487906, 0.059988, 0.870832);
            setVector(36, 0.468001, 0.231001, 0.853003);
            setVector(37, 0.426147, 0.394136, 0.814282);
            setVector(38, 0.266963, 0.442939, 0.855883);
            setVector(39, 0.120952, 0.485807, 0.865657);
            setVector(40, -0.024989, 0.534766, 0.844631);
            setVector(41, -0.176029, 0.491082, 0.853142);
            setVector(42, -0.323126, 0.437171, 0.839328);
            setVector(43, -0.463978, 0.371982, 0.803961);
            setVector(44, -0.613900, 0.139977, 0.776874);
            setVector(45, -0.617975, -0.018999, 0.785968);
            setVector(46, -0.493054, -0.117013, 0.862094);
            setVector(47, -0.503790, -0.292878, 0.812661);
            setVector(48, -0.383865, -0.394861, 0.834705);
            setVector(49, -0.253926, -0.486859, 0.835757);
            setVector(50, -0.113952, -0.565762, 0.816657);
            setVector(51, 0.037997, -0.625950, 0.778937);
            setVector(52, 0.182983, -0.543949, 0.818924);
            setVector(53, 0.328974, -0.543956, 0.771938);
            setVector(54, 0.446988, -0.437988, 0.779979);
            setVector(55, 0.565718, -0.320840, 0.759621);
            setVector(56, 0.595011, -0.163003, 0.787015);
            setVector(57, 0.611304, -0.002001, 0.791393);
            setVector(58, 0.600103, 0.163028, 0.783134);
            setVector(59, 0.568061, 0.325035, 0.756081);
            setVector(60, 0.512003, 0.482003, 0.711004);
            setVector(61, 0.362192, 0.543287, 0.757401);
            setVector(62, 0.211049, 0.594139, 0.776181);
            setVector(63, 0.058967, 0.640642, 0.765572);
            setVector(64, -0.119057, 0.637302, 0.761361);
            setVector(65, -0.274102, 0.590219, 0.759282);
            setVector(66, -0.421992, 0.527990, 0.736986);
            setVector(67, -0.560010, 0.454008, 0.693012);
            setVector(68, -0.598017, 0.298008, 0.744021);
            setVector(69, -0.730015, 0.044001, 0.682014);
            setVector(70, -0.718931, -0.123988, 0.683934);
            setVector(71, -0.607904, -0.181971, 0.772878);
            setVector(72, -0.643781, -0.311894, 0.698762);
            setVector(73, -0.524667, -0.433725, 0.732535);
            setVector(74, -0.395115, -0.537156, 0.745217);
            setVector(75, -0.255897, -0.624749, 0.737703);
            setVector(76, -0.106010, -0.694064, 0.712066);
            setVector(77, 0.051982, -0.740746, 0.669771);
            setVector(78, 0.198017, -0.665057, 0.720062);
            setVector(79, 0.338997, -0.682994, 0.646994);
            setVector(80, 0.455869, -0.583832, 0.671807);
            setVector(81, 0.567210, -0.477177, 0.671249);
            setVector(82, 0.670100, -0.360054, 0.649097);
            setVector(83, 0.704211, -0.210063, 0.678203);
            setVector(84, 0.723245, -0.057019, 0.688233);
            setVector(85, 0.720289, 0.101041, 0.686276);
            setVector(86, 0.696042, 0.258015, 0.670040);
            setVector(87, 0.650084, 0.411053, 0.639082);
            setVector(88, 0.582844, 0.557851, 0.590842);
            setVector(89, 0.441061, 0.629087, 0.640089);
            setVector(90, 0.292091, 0.685213, 0.667208);
            setVector(91, 0.136935, 0.732653, 0.666684);
            setVector(92, -0.035006, 0.744130, 0.667117);
            setVector(93, -0.207001, 0.724004, 0.658004);
            setVector(94, -0.361776, 0.671585, 0.646600);
            setVector(95, -0.506971, 0.602966, 0.615965);
            setVector(96, -0.641109, 0.517088, 0.567097);
            setVector(97, -0.690897, 0.366945, 0.622907);
            setVector(98, -0.718746, 0.206927, 0.663766);
            setVector(99, -0.821800, 0.110973, 0.558864);
            setVector(100, -0.820917, -0.048995, 0.568942);
            setVector(101, -0.787722, -0.210926, 0.578796);
            setVector(102, -0.733198, -0.353096, 0.581157);
            setVector(103, -0.632761, -0.470822, 0.614768);
            setVector(104, -0.516033, -0.577036, 0.633040);
            setVector(105, -0.386042, -0.669073, 0.635069);
            setVector(106, -0.244127, -0.745387, 0.620322);
            setVector(107, -0.092982, -0.802848, 0.588888);
            setVector(108, 0.061971, -0.837608, 0.542746);
            setVector(109, 0.207049, -0.774182, 0.598141);
            setVector(110, 0.355042, -0.784093, 0.509061);
            setVector(111, 0.480125, -0.691179, 0.540140);
            setVector(112, 0.596018, -0.585017, 0.550016);
            setVector(113, 0.703421, -0.465279, 0.537322);
            setVector(114, 0.781234, -0.304091, 0.545163);
            setVector(115, 0.812165, -0.145029, 0.565115);
            setVector(116, 0.821641, 0.016993, 0.569751);
            setVector(117, 0.808928, 0.175984, 0.560950);
            setVector(118, 0.774726, 0.330883, 0.538810);
            setVector(119, 0.718831, 0.477888, 0.504881);
            setVector(120, 0.642730, 0.615741, 0.455808);
            setVector(121, 0.509022, 0.697030, 0.505021);
            setVector(122, 0.366176, 0.758364, 0.539259);
            setVector(123, 0.213039, 0.808148, 0.549100);
            setVector(124, 0.045984, 0.830708, 0.554805);
            setVector(125, -0.123006, 0.826039, 0.550026);
            setVector(126, -0.288131, 0.793361, 0.536244);
            setVector(127, -0.439982, 0.734970, 0.515979);
            setVector(128, -0.579687, 0.658644, 0.479741);
            setVector(129, -0.705004, 0.567003, 0.426002);
            setVector(130, -0.764589, 0.424772, 0.484740);
            setVector(131, -0.800709, 0.269902, 0.534806);
            setVector(132, -0.892958, 0.161992, 0.419980);
            setVector(133, -0.900927, 0.002000, 0.433965);
            setVector(134, -0.879831, -0.158969, 0.447914);
            setVector(135, -0.833019, -0.319007, 0.452010);
            setVector(136, -0.747071, -0.470045, 0.470045);
            setVector(137, -0.637748, -0.587767, 0.497803);
            setVector(138, -0.513898, -0.689863, 0.509899);
            setVector(139, -0.377037, -0.775076, 0.507050);
            setVector(140, -0.230976, -0.842913, 0.485950);
            setVector(141, -0.081983, -0.889811, 0.448905);
            setVector(142, 0.074014, -0.912169, 0.403075);
            setVector(143, 0.219092, -0.860361, 0.460193);
            setVector(144, 0.362034, -0.861080, 0.357033);
            setVector(145, 0.490041, -0.778065, 0.393033);
            setVector(146, 0.608025, -0.678028, 0.413017);
            setVector(147, 0.718996, -0.558997, 0.412998);
            setVector(148, 0.810079, -0.408040, 0.421041);
            setVector(149, 0.870877, -0.239966, 0.428940);
            setVector(150, 0.894955, -0.073996, 0.439978);
            setVector(151, 0.896096, 0.088009, 0.435047);
            setVector(152, 0.874684, 0.243912, 0.418849);
            setVector(153, 0.831860, 0.391934, 0.392934);
            setVector(154, 0.768978, 0.530985, 0.355990);
            setVector(155, 0.687097, 0.657093, 0.310044);
            setVector(156, 0.564702, 0.745606, 0.353813);
            setVector(157, 0.432080, 0.812149, 0.392072);
            setVector(158, 0.286963, 0.863887, 0.413946);
            setVector(159, 0.126973, 0.895811, 0.425910);
            setVector(160, -0.039977, 0.903483, 0.426756);
            setVector(161, -0.206046, 0.885199, 0.417094);
            setVector(162, -0.364892, 0.841750, 0.397882);
            setVector(163, -0.510950, 0.775924, 0.369964);
            setVector(164, -0.639907, 0.695898, 0.325952);
            setVector(165, -0.753755, 0.597806, 0.272911);
            setVector(166, -0.818055, 0.471032, 0.330022);
            setVector(167, -0.862916, 0.318969, 0.391962);
            setVector(168, -0.938282, 0.216065, 0.270081);
            setVector(169, -0.956045, 0.055003, 0.288014);
            setVector(170, -0.946681, -0.106964, 0.303898);
            setVector(171, -0.909812, -0.267945, 0.316935);
            setVector(172, -0.839417, -0.429702, 0.332769);
            setVector(173, -0.741292, -0.573226, 0.349138);
            setVector(174, -0.624959, -0.687955, 0.368976);
            setVector(175, -0.493966, -0.784946, 0.373974);
            setVector(176, -0.353977, -0.861944, 0.362977);
            setVector(177, -0.211958, -0.918819, 0.332934);
            setVector(178, -0.066965, -0.953504, 0.293847);
            setVector(179, 0.084956, -0.963501, 0.253868);
            setVector(180, 0.228069, -0.924279, 0.306092);
            setVector(181, 0.333060, -0.925167, 0.182033);
            setVector(182, 0.466929, -0.853870, 0.229965);
            setVector(183, 0.595385, -0.760492, 0.259168);
            setVector(184, 0.709740, -0.647763, 0.276899);
            setVector(185, 0.814356, -0.503220, 0.289126);
            setVector(186, 0.892778, -0.337916, 0.297926);
            setVector(187, 0.938954, -0.162992, 0.302985);
            setVector(188, 0.954223, 0.004001, 0.299070);
            setVector(189, 0.944039, 0.166007, 0.285012);
            setVector(190, 0.909805, 0.321931, 0.261944);
            setVector(191, 0.852595, 0.469777, 0.228891);
            setVector(192, 0.771212, 0.608168, 0.188052);
            setVector(193, 0.644272, 0.742313, 0.184078);
            setVector(194, 0.512029, 0.830048, 0.221013);
            setVector(195, 0.371982, 0.891957, 0.256988);
            setVector(196, 0.218003, 0.934014, 0.283004);
            setVector(197, 0.147000, 0.978001, 0.148000);
            setVector(198, 0.050005, 0.955101, 0.292031);
            setVector(199, -0.120044, 0.950347, 0.287105);
            setVector(200, -0.287133, 0.918424, 0.272126);
            setVector(201, -0.444757, 0.861529, 0.244866);
            setVector(202, -0.582982, 0.786976, 0.201994);
            setVector(203, -0.707118, 0.690116, 0.154026);
            setVector(204, -0.823874, 0.547916, 0.144978);
            setVector(205, -0.895869, 0.373945, 0.239965);
            setVector(206, -0.956936, 0.266982, 0.113992);
            setVector(207, -0.984969, 0.108997, 0.133996);
            setVector(208, -0.987109, -0.053006, 0.151017);
            setVector(209, -0.962464, -0.213103, 0.168081);
            setVector(210, -0.907966, -0.372986, 0.190993);
            setVector(211, -0.822992, -0.528995, 0.206998);
            setVector(212, -0.712235, -0.666220, 0.221073);
            setVector(213, -0.584746, -0.777662, 0.230900);
            setVector(214, -0.445210, -0.866409, 0.226107);
            setVector(215, -0.302832, -0.932484, 0.196891);
            setVector(216, -0.155030, -0.975188, 0.158031);
            setVector(217, 0.003001, -0.992396, 0.123049);
            setVector(218, -0.105043, 0.994407, 0.011005);
            setVector(219, 0.184970, -0.973842, 0.131979);
            setVector(220, 0.279071, -0.960245, 0.007002);
            setVector(221, 0.426883, -0.902753, 0.052986);
            setVector(222, 0.560256, -0.822376, 0.099045);
            setVector(223, 0.687169, -0.716176, 0.122030);
            setVector(224, 0.794057, -0.589042, 0.150011);
            setVector(225, 0.887491, -0.430238, 0.165091);
            setVector(226, 0.951729, -0.254927, 0.170951);
            setVector(227, 0.983436, -0.085038, 0.160071);
            setVector(228, 0.985899, 0.079992, 0.146985);
            setVector(229, 0.961934, 0.241983, 0.126991);
            setVector(230, 0.912643, 0.396845, 0.097962);
            setVector(231, 0.837257, 0.543167, 0.063019);
            setVector(232, 0.723919, 0.687923, 0.051994);
            setVector(233, 0.588089, 0.807122, 0.052008);
            setVector(234, 0.450862, 0.888729, 0.082975);
            setVector(235, 0.304931, 0.944787, 0.119973);
            setVector(236, -0.233066, -0.972275, 0.019005);
            setVector(237, 0.073019, 0.997258, 0.012003);
            setVector(238, -0.026001, 0.988038, 0.152006);
            setVector(239, -0.200088, 0.969424, 0.142062);
            setVector(240, -0.365992, 0.922979, 0.118997);
            setVector(241, -0.510695, 0.856488, 0.074955);
            setVector(242, -0.640958, 0.766950, 0.030998);
            setVector(243, -0.764977, 0.643980, 0.010000);
            setVector(244, 0.857996, -0.512998, 0.026000);
            setVector(245, -0.901024, 0.423011, 0.096003);
            setVector(246, 0.938499, -0.343182, 0.038020);
            setVector(247, 0.983769, -0.177958, 0.022995);
            setVector(248, 0.999898, -0.012999, 0.005999);
            setVector(249, -0.988633, -0.149944, 0.010996);
            setVector(250, -0.949788, -0.310931, 0.034992);
            setVector(251, -0.883406, -0.464213, 0.064029);
            setVector(252, -0.783923, -0.615940, 0.077992);
            setVector(253, -0.660104, -0.746118, 0.087014);
            setVector(254, -0.521822, -0.848711, 0.085971);
            setVector(255, -0.380935, -0.922844, 0.056990);
            break;

        default:
            if(! bIsContextPrepForBinarySearch ||  (m_lDebugLevel & DEBUG_RETURN))
            {
                SEQ_TRACE_INFO.print("ERROR: %ld directions not supported, or error in vector definition file", lDirections);
            }
            return false;
    } // of switch

    calcGreatestNorm();

    return true;
}

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

    If no external table is found (i.e. an invalid number of lDirections),
    the function will return false.

    As fallback (e.g. if the name of the external vector file is not yet
    known, an internal set of FALLBACK_DIRECTIONS diffusion directions is
    provided. Calling the method with this number will always yield a
    valid vector set.

    When fSeqPrep() runs in "binary search" mode for the "Directions" parameter,
    this is actually no binary search, but a point-by-point query.
    (There is no contingeous area of valid parameter values.)
    Therefore the performance of this function is quiet crucial for
    the UI behaviour.

    This method is available on the host only.

    \return true if a vector set could be prepared successfully.

    \post m_dGreatestNorm and m_dGreatestComponent will we set.
    */
bool DiffusionDirections::prepExternal
(
long lDirections,                   /**< Imp: Number of diffusion directions */
bool bIsContextPrepForBinarySearch  /**< Imp: Suppress warnings              */
)
{
    if((lDirections < 0) || (m_cLookupTable[lDirections] == 0))
    {
        // A vector set for this direction number does not exist
        if(!bIsContextPrepForBinarySearch || (m_lDebugLevel & DEBUG_RETURN))
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "Invalid number of directions: " << lDirections << ".";
            m_sErrorMessage = sErrorMsg.str();

            SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str());
        }
        return false;
    }

    // Check the external file - a customer wants to use his own definition
    if(m_cLookupTable[lDirections] == +1)
    {
        if(!readFromFile(lDirections))
        {
            // The requested vector set was not found or it is not good
            // Note: m_sErrorMessage gets set by readFromFile()

            SEQ_TRACE_INFO.print("Error interpreting section for %ld directions: '%s'.", lDirections, m_sErrorMessage.c_str());
            return false;
        }
    }
    else if(lDirections == FALLBACK_DIRECTIONS)
    {
        m_lDirections       = lDirections;
        delete[] m_sVector;
        m_sVector           = new VectorStruct[m_lDirections];
        m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
        m_sUserComment      = "SIEMENS\nInternal 6-directions set";
        // (Note: do not modify the following comments - doxygen greps for them)
        // vector set for 6 directions (normalised to maximum)
        setVector(0, 1.0, 0.0, 1.0);
        setVector(1, -1.0, 0.0, 1.0);
        setVector(2, 0.0, 1.0, 1.0);
        setVector(3, 0.0, 1.0, -1.0);
        setVector(4, 1.0, 1.0, 0.0);
        setVector(5, -1.0, 1.0, 0.0);
    }
    else
    {
        if(!bIsContextPrepForBinarySearch || (m_lDebugLevel & DEBUG_RETURN))
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "'" << lDirections << "' directions not defined in external file.";
            m_sErrorMessage = sErrorMsg.str();

            SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str());
        }
        return false;
    }

    calcGreatestNorm();

    return true;
}
#endif // #ifdef WIN32


// ===========================================================================
/// Prepare a diffusion vector set using explicit parameters
/**
    A new vector set will be generated which holds the vectors
    for the specifed number of directions. The memory will be
    released and allocated dynamically.

    When fSeqPrep() runs in "binary search" mode for the "Directions" parameter,
    this is actually no binary search, but a point-by-point query.
    (There is no contingeous area of valid parameter values.)
    Therefore the performance of this function is quiet crucial for
    the UI behaviour.

    \return true if a vector set could be prepared successfully.

    \post m_dGreatestNorm and m_dGreatestComponent will we set.
    */
bool DiffusionDirections::prepExplicit
(
long                                     lDirections,                   /**< Imp: Number of diffusion directions */
MrProtocolData::DiffDirCoordinateSystem  eCoordinateSystem,             /**< Imp: Coordinates system             */
const std::string                       &strComment,                    /**< Imp: User comment (for tooltip)     */
const std::vector<VectorStruct>         &vDiffDir,                      /**< Imp: Array of diffusion vectors     */
bool                                     bIsContextPrepForBinarySearch  /**< Imp: Suppress warnings              */
)
{
    if(lDirections <= 0)
    {
        if(!bIsContextPrepForBinarySearch || (m_lDebugLevel & DEBUG_RETURN))
        {
            SEQ_TRACE_INFO.print("ERROR: lDirections=%ld is invalid.", lDirections);
        }
        return false;
    }

    m_lDirections       = lDirections;
    delete[] m_sVector;
    m_sVector           = new VectorStruct[m_lDirections];
    m_eCoordinateSystem = eCoordinateSystem;
    m_sUserComment      = strComment;

    for(long lI = 0; lI < m_lDirections; ++lI)
    {
        setVector(lI, vDiffDir[lI].dx, vDiffDir[lI].dy, vDiffDir[lI].dz);
    }

    calcGreatestNorm();

    return true;
}

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
long DiffusionDirections::prepQSpace
(
SEQ::DiffQSpaceCoverageMode   eQSpaceCoverage,               /**< Imp: Q-space coverage (full / half)        */
SEQ::DiffQSpaceSamplingScheme eQSpaceSampling,               /**< Imp: Q-space sampling (cartesian / radial) */
long                          lQSpaceSteps,                  /**< Imp: Number steps to cover in Q-space      */
bool                          bIsContextPrepForBinarySearch  /**< Imp: Suppress warnings                     */
)
{
    // Check input parameters
    if(lQSpaceSteps <= 0)
    {
        if(!bIsContextPrepForBinarySearch || (m_lDebugLevel & DEBUG_RETURN))
        {
            SEQ_TRACE_INFO.print("ERROR: lQSpaceSteps=%ld is invalid.", lQSpaceSteps);
        }
        return 0;
    }

    // Temporary vector storage
    std::vector<VectorStruct> vsTempVector;

    // Factor for normalizing vectors to unit sphere
    double dQSpaceSteps = static_cast<double>(lQSpaceSteps);

    switch(eQSpaceSampling)
    {
        case SEQ::DIFF_QSPACE_SAMPLING_CARTESIAN:
        {
            // Pattern contains all q-space coordinates on a cartesian
            // grid within a unit sphere.

            // Step 1.:     x x x x x       . x x x . 
            //              x x x x x       x x x x x
            //              x x x x x   =>  x x x x x
            //              x x x x x       x x x x x
            //              x x x x x       . x x x .

            // Reserve memory (estimated value, to avoid re-allocation)
            vsTempVector.reserve(8 * lQSpaceSteps * lQSpaceSteps * lQSpaceSteps);

            // Loop over selected q-space region
            for(long lQx = -lQSpaceSteps; lQx <= lQSpaceSteps; ++lQx)
            {
                double dQx = static_cast<double>(lQx) / dQSpaceSteps;

                for(long lQy = -lQSpaceSteps; lQy <= lQSpaceSteps; ++lQy)
                {
                    double dQy = static_cast<double>(lQy) / dQSpaceSteps;

                    for(long lQz = -lQSpaceSteps; lQz <= lQSpaceSteps; ++lQz)
                    {
                        if((lQx * lQx + lQy * lQy + lQz * lQz) <= lQSpaceSteps * lQSpaceSteps)
                        {
                            double dQz = static_cast<double>(lQz) / dQSpaceSteps;

                            // Copy vector to temporary storage
                            VectorStruct sVector(dQx, dQy, dQz);
                            vsTempVector.push_back(sVector);
                        }
                    }
                }
            }
            break;
        }
        case SEQ::DIFF_QSPACE_SAMPLING_RADIAL:
            // Not implemented so far
        default:
            // Configuration error: trace always
            SEQ_TRACE_INFO.print("ERROR: eQSpaceSampling=%d is invalid.", eQSpaceSampling);
            return 0;
    }   // of switch

    // Sort vectors according to distance from q-space center (uses dedicated VectorStruct operator '<')
    // Note: An alternative sorting schemes could toggle between low and high b-values. 
    //       - Advantage 'toggle-sort': Temporal distribution of thermal load
    //       - Advantage 'q-space-out': Consecutive filling of complete q-space spheres permits interim calculations
    std::sort(vsTempVector.begin(), vsTempVector.end());

    switch(eQSpaceCoverage)
    {
        case SEQ::DIFF_QSPACE_COVERAGE_FULL:
            // Do nothing
            break;
        case SEQ::DIFF_QSPACE_COVERAGE_HALF:
        {
            // Scan over all diffusion vectors and keep only those residing in the desired hemisphere
            for(std::vector<VectorStruct>::iterator itTempVector = vsTempVector.begin(); itTempVector != vsTempVector.end(); /* Iterator gets advanced in loop */)
            {
                VectorStruct sDiffVect((*itTempVector).dx, (*itTempVector).dy, (*itTempVector).dz);     // For convenience     
                VectorStruct sCrossProduct(0., 0., 0.);

                double dNonZeroCoord1 = 0.;
                double dNonZeroCoord2 = 0.;

                if(fAlmostEqual(sDiffVect.dx, 0.0) && fAlmostEqual(sDiffVect.dy, 0.0) && fAlmostEqual(sDiffVect.dz, 0.0))
                {
                    // Keep q-space center, continue with next vector
                    itTempVector++;
                }
                else
                {
                    // Consider angle phi between current diffusion vector D and (normalized) hemisphere vector H:
                    //      D * H / |D| = cos( phi )
                    //      angle       = asin( D * H / |D| ) = pi/2 - acos( D * H / |D| ) = pi/2 - phi
                    // =>   angle       > 0    if D and H are parallel
                    //      angle       = 0    if D and H are perpendicular
                    //      angle       < 0    if D and H are antiparallel
                    //
                    // Step 2.:     . x x x .       . x x x . 
                    //              x x x x x       x x x x x
                    //              x x x x x   =>  x x x x x
                    //              x x x x x       . . . . .
                    //              . x x x .       . . . . .
                    double dNorm  =               sqrt(sDiffVect.dx *           sDiffVect.dx + sDiffVect.dy *           sDiffVect.dy + sDiffVect.dz *           sDiffVect.dz);
                    double dAngle = 180. / M_PI * asin((sDiffVect.dx * m_sQSpaceHemisphere.dx + sDiffVect.dy * m_sQSpaceHemisphere.dy + sDiffVect.dz * m_sQSpaceHemisphere.dz) / dNorm);

                    // Consider cross product: this is required to identify the vectors in the plane
                    // intersecting the two hemisphere that actually need to be scanned (only one half
                    // of these vectors are actually required)
                    // 
                    // Step 3.:     . x x x .       . x x x . 
                    //              x x x x x       x x x x x
                    //              x x x x x   =>  x x x . .
                    //              . . . . .       . . . . .
                    //              . . . . .       . . . . .
                    sCrossProduct.dx = m_sQSpaceHemisphere.dy * sDiffVect.dz - m_sQSpaceHemisphere.dz * sDiffVect.dy;
                    sCrossProduct.dy = m_sQSpaceHemisphere.dz * sDiffVect.dx - m_sQSpaceHemisphere.dx * sDiffVect.dz;
                    sCrossProduct.dz = m_sQSpaceHemisphere.dx * sDiffVect.dy - m_sQSpaceHemisphere.dy * sDiffVect.dx;

                    if(fAlmostEqual(sCrossProduct.dx, 0.))
                    {
                        dNonZeroCoord1 = sCrossProduct.dy;
                        dNonZeroCoord2 = sCrossProduct.dz;
                    }
                    else if(fAlmostEqual(sCrossProduct.dy, 0.))
                    {
                        dNonZeroCoord1 = sCrossProduct.dx;
                        dNonZeroCoord2 = sCrossProduct.dz;
                    }
                    else
                    {
                        dNonZeroCoord1 = sCrossProduct.dx;
                        dNonZeroCoord2 = sCrossProduct.dy;
                    }

                    // Do we have to keep the current diffusion vector?
                    if((!fAlmostEqual(dAngle, 0.)                                       && (dAngle         <  0.))
                       || (fAlmostEqual(dAngle, 0.)                                       && (dNonZeroCoord1 <  0.))
                       || (fAlmostEqual(dAngle, 0.) && fAlmostEqual(dNonZeroCoord1, 0.) && (dNonZeroCoord2 <= 0.)))
                    {
                        // Remove current vector, set iterator to next item
                        itTempVector = vsTempVector.erase(itTempVector);
                    }
                    else
                    {
                        // Advance iterator
                        itTempVector++;
                    }
                }
            }
        }
        break;
        default:
            // Configuration error: trace always
            SEQ_TRACE_INFO.print("ERROR: eQSpaceCoverage=%d is invalid.", eQSpaceSampling);
            return 0;
    }

    // Delete old content
    delete[] m_sVector;

    // Apply new vector set
    m_lDirections       = static_cast<long>(vsTempVector.size());
    m_sVector           = new VectorStruct[m_lDirections];
    m_eCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
    m_sUserComment      = "SIEMENS\nInternal q-space sampling scheme";

    for(long lI = 0; lI < m_lDirections; ++lI)
    {
        setVector(lI, vsTempVector[lI].dx, vsTempVector[lI].dy, vsTempVector[lI].dz);
    }

    calcGreatestNorm();

    return m_lDirections;
}

// ===========================================================================
/// Calculate greatest norm and greatest component
/**
    Should be called after successful generation of an (internal or
    external) diffusion vector set.
    */
void DiffusionDirections::calcGreatestNorm(void)
{
    // Calculate greatest norm and greatest component
    m_dGreatestNorm       = 0.0;
    m_dGreatestComponent  = 0.0;
    double dTempNorm      = 0.0;
    double dTempComponent = 0.0;

    for(long lIndex=0; lIndex < m_lDirections; lIndex++)
    {
        dTempNorm      = getNorm(lIndex);
        dTempComponent = std::max(fabs(m_sVector[lIndex].dx), std::max(fabs(m_sVector[lIndex].dy), fabs(m_sVector[lIndex].dz)));

        if(dTempNorm > m_dGreatestNorm)
        {
            m_dGreatestNorm = dTempNorm;
        }

        if(dTempComponent > m_dGreatestComponent)
        {
            m_dGreatestComponent = dTempComponent;
        }
    }

    return;
}


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
bool DiffusionDirections::readFromFile
(
long lDirections        /**< Imp: Vector set to be loaded (characterized by the number of directions). */
)
// ===========================================================================
{
    // some macros for debugging:
#define DEBUG_readFromFile(a) // {a;}

    long lLineNumber    = 0;
    long lSectionsFound = 0;
    bool bSuccess       = true;
    char line[256];
    char ptCommand[256], ptParam[80], ptIndex[80], ptValue[256];
    int  ii;
    long lIndex         = 0;
    bool bGate          = false;
    // Dummy initializations
    std::string                             sNewComment("");
    MrProtocolData::DiffDirCoordinateSystem eNewCoordinateSystem  = MrProtocolData::DIFFDIR_CS_XYZ;
    MrProtocolData::DiffDirNormalization    eNewNormalisationMode = MrProtocolData::DIFFDIR_NORM_UNITY;

    // The strategy is to write first only a temporary vector buffer.
    // Then it will be checked and finally copied to m_sVector
    VectorStruct *sVectorBuffer = NULL;

    DEBUG_readFromFile(SEQ_TRACE_INFO.print("CALL DiffusionDirections::readFromFile (long lDirections=%ld)\n", lDirections));

    // Assemble filename
    std::stringstream sVectorFile;

    char* ptCustomerSeq = NULL;
    ptCustomerSeq = getenv(VECTORFILEENV);

    if(ptCustomerSeq == NULL)
    {
        std::stringstream sErrorMsg;
        sErrorMsg << "Environment does not exist: " << VECTORFILEENV;
        m_sErrorMessage = sErrorMsg.str();

        DEBUG_readFromFile(SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str()));
        return false;
    }
    else
    {
        // e.g. "%CustomerSeq%\\DiffusionVectorSets\\DiffusionVectors.dvs"
        sVectorFile << ptCustomerSeq << "\\" << VECTORFILEPATH << "\\" << m_sVectorFileName;
    }

    DEBUG_readFromFile(SEQ_TRACE_INFO.print("m_sVectorFileName.c_str(): %s", sVectorFile.str().c_str()));

    // open vector file
    FILE* pInputFile = NULL;
    if((pInputFile = fopen(sVectorFile.str().c_str(), "r")) == NULL)
    {
        std::stringstream sErrorMsg;
        sErrorMsg << "Cannot open diffusion vector set: " << sVectorFile.str();
        m_sErrorMessage = sErrorMsg.str();

        DEBUG_readFromFile(SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str()));
        return false;
    }


    if(lDirections)
    {
        sVectorBuffer = new VectorStruct[lDirections];
        memset(sVectorBuffer, 0, lDirections * sizeof(VectorStruct));
        // Remember: From now on no 'return's until the next 'delete'
    }

    while(fgets(line, 240, pInputFile) != NULL)
    {
        // Replace '\n' by '\0'
        for(int iMyIndex = 0; iMyIndex < 256; iMyIndex++)
        {
            if(line[iMyIndex] == '\n')
            {
                line[iMyIndex] = '\0';
            }
        }

        lLineNumber++;

        // map to lower case
        int i = 0;
        while(line[i])
        {
            line[i] = char(tolower(line[i]));  // cast due to warning C4244: '=' : conversion from 'int' to 'char'
            i++;
        }

        // strip leading blanks / tabs
        char *pline;
        pline = line;
        while((pline[0] == ' ') || (pline[0] == '\t'))
        {
            pline++;
        }

        // skip empty lines
        if(pline[0] == 0) continue;

        // Ignore comment lines
        if(pline[0] == '#') continue;

        // start of a new section:
        if(strstr(line, "directions="))
        {
            long lDir = 0;
            if(sscanf(line, "[directions=%ld]", &lDir) != 1)
            {
                std::stringstream sErrorMsg;
                sErrorMsg << "Invalid section specified:  '" << line << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
                m_sErrorMessage = sErrorMsg.str();

                bSuccess = false;
                continue;
            }

            // Are we in section scan mode? So just scan...
            if(lDirections == 0)
            {
                DEBUG_readFromFile(SEQ_TRACE_INFO.print("New Section Direction = %li", lDir));
                lSectionsFound++;
                m_cLookupTable[lDir] = 1;
                continue;
            }

            // Have we found the requested section? Then open the gate...
            if(lDirections == lDir)
            {
                DEBUG_readFromFile(SEQ_TRACE_INFO.print("Target found: %li", lDirections));
                bGate = true;
                lSectionsFound++;
                continue;
            }

            // Is the gate already open and a new sections starts? Then we are done now...
            if(bGate) break;

            // We are still looking for our section? Let's continue searching ...
            continue;
        } // end of  if ( strstr(line, "directions=") )


        if(! bGate) continue;


        // cout << "Using " << line << endl;


        // Parse the line 
        // --------------
        strcpy(ptCommand, line);
        strcpy(ptIndex, "");

        // Get parameter name index
        // Note: lenght of the input line must not exceed size of int
        ii = static_cast<int>(strchr(ptCommand, ']') - strchr(ptCommand, '['));
        if(ii)
        {
            ii--;
            strncpy(ptIndex, strchr(ptCommand, '[') + 1, ii);
        }
        ptIndex[ii] = '\0';
        lIndex = atol(ptIndex);
        if((lIndex < 0) || (lIndex >= lDirections))
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "Index out of range in command:  '" << ptCommand << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
            m_sErrorMessage = sErrorMsg.str();

            bSuccess = false;
            continue;
        }


        // Get parameter name
        ii = 0;
        while(isalnum(ptCommand[ii]))  // true for letter or digit
        {
            *(ptCommand + ii) = char(tolower(*(ptCommand + ii)));
            ii++;
        }
        if(!strncpy(ptParam, ptCommand, ii))
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "No parameter specified in command:  '" << ptCommand << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
            m_sErrorMessage = sErrorMsg.str();

            bSuccess = false;
            continue;
        }
        ptParam[ii] = '\0';


        // Get parameter value
        if(strchr(ptCommand, '=') == NULL)
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "No value assignment in command:  '" << ptCommand << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
            m_sErrorMessage = sErrorMsg.str();

            bSuccess = false;
            continue;
        }

    {
        char * charptr;
        charptr = strchr(ptCommand, '=') + 1;
        while(isspace(*charptr))
        {
            // remove leading blanks
            charptr++;
        }
        strcpy(ptValue, charptr);

        while(strlen(ptValue) && isspace(ptValue[strlen(ptValue) - 1]))
        {
            // remove trailing blanks
            ptValue[strlen(ptValue)-1] = '\0';
        }

        if(strlen(ptValue) == 0)
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "No value assignment in command:  '" << ptCommand << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
            m_sErrorMessage = sErrorMsg.str();

            bSuccess = false;
            continue;
        }
    }

    DEBUG_readFromFile(SEQ_TRACE_INFO.print(": parsed Param <%s> [%s] = %s", ptParam, ptIndex, ptValue));


    // Execute command of current line
    // ---------------

    if(strstr(ptParam, "coordinatesystem"))
    {
        // command: CoordinateSystem
        if(strstr(ptValue, "xyz"))
        {
            // cout << "KS = xyz" << endl;
            eNewCoordinateSystem = MrProtocolData::DIFFDIR_CS_XYZ;
            continue;
        }
        else if(strstr(ptValue, "prs"))
        {
            // cout << "KS = prs" << endl;
            eNewCoordinateSystem = MrProtocolData::DIFFDIR_CS_PRS;
            continue;
        }
        else
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "Invalid coordinate system in command:  '" << ptCommand << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
            m_sErrorMessage = sErrorMsg.str();

            bSuccess = false;
            continue;
        }
    }

    if(strstr(ptParam, "normalisation") || strstr(ptParam, "normalization"))
    {
        // command: normalisation
        if(strstr(ptValue, "none"))
        {
            eNewNormalisationMode = MrProtocolData::DIFFDIR_NORM_NONE;
            continue;
        }
        else if(strstr(ptValue, "unity"))
        {
            eNewNormalisationMode = MrProtocolData::DIFFDIR_NORM_UNITY;
            continue;
        }
        else if(strstr(ptValue, "maximum"))
        {
            eNewNormalisationMode = MrProtocolData::DIFFDIR_NORM_MAX;
            continue;
        }
        else
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "Invalid normalisation mode in command:  '" << ptCommand << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
            m_sErrorMessage = sErrorMsg.str();

            bSuccess = false;
            continue;
        }
    }
    else if(strstr(ptParam, "vector"))
    {
        // command: vector
        double dx, dy, dz;
        int status = sscanf(ptValue, "(%lf,%lf,%lf)", &dx, &dy, &dz);
        DEBUG_readFromFile(SEQ_TRACE_INFO.print("%s -->    dx = %f    dy = %f    dz = %f    index = %li", ptValue, dx, dy, dz, lIndex));
        if(status != 3)
        {
            std::stringstream sErrorMsg;
            sErrorMsg << "Invalid vector definition:  '" << line << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
            m_sErrorMessage = sErrorMsg.str();

            bSuccess = false;
            continue;
        }
        sVectorBuffer[lIndex].dx = dx;
        sVectorBuffer[lIndex].dy = dy;
        sVectorBuffer[lIndex].dz = dz;
        continue;
    }
    else if(strstr(ptParam, "comment"))
    {
        // command: comment
        sNewComment.assign(ptValue, lMaxCommentLength);
        continue;
    }
    else
    {
        // error: no valid command
        std::stringstream sErrorMsg;
        sErrorMsg << "Invalid command :  '" << line << "' \n(" << sVectorFile.str() << ", line " << lLineNumber << ")";
        m_sErrorMessage = sErrorMsg.str();

        bSuccess = false;
        continue;
    }

    } // end of while ( !input.eof() )

    fclose(pInputFile);


    if(lDirections == 0)
    {
        // we were just scanning for available vector sets
        m_lExternalDirectionSets = lSectionsFound;

        // return false if no valid vector set has been found
        if(!bSuccess || (m_lExternalDirectionSets == 0))
        {
            // error: no valid direction sets found
            if(m_lExternalDirectionSets == 0)
            {
                std::stringstream sErrorMsg;
                sErrorMsg << "No valid sections found in " << sVectorFile.str();
                m_sErrorMessage = sErrorMsg.str();
            }
            SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str());
            delete[] sVectorBuffer;
            return false;
        }

        SEQ_TRACE_INFO.print("%ld sections found in %s", m_lExternalDirectionSets, sVectorFile.str().c_str());
        return true;
    }

    if(!lSectionsFound)
    {
        // error: no valid direction set found
        std::stringstream sErrorMsg;
        sErrorMsg << "No section for " << lDirections << " directions found in " << sVectorFile.str();
        m_sErrorMessage = sErrorMsg.str();

        if(m_lDebugLevel & DEBUG_RETURN)
        {
            SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str());
        }
        delete[] sVectorBuffer;
        return false;
    }

    if(! bSuccess)
    {
        m_cLookupTable[lDirections] = 0;
        delete[] sVectorBuffer;
        return false;
    }

    // Perform normalization
    if(!normalize(sVectorBuffer, lDirections, eNewNormalisationMode))
    {
        // error: normalization failed
        std::stringstream sErrorMsg;
        sErrorMsg << "Normalization of vectors failed (null vector? Not all " << lDirections << " vectors defined?)";
        m_sErrorMessage = sErrorMsg.str();

        SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str());

        m_cLookupTable[lDirections] = 0;
        delete[] sVectorBuffer;
        return false;
    }

    // Final checks
    if(eNewCoordinateSystem == MrProtocolData::DIFFDIR_CS_PRS)
    {
        // For diffusion vector sets in PRS coordinates, no vector
        // norm must exceed a value off 1
        for(long lI = 0; lI < lDirections; ++lI)
        {
            if(sqrt(sVectorBuffer[lI].dx * sVectorBuffer[lI].dx +
                sVectorBuffer[lI].dy * sVectorBuffer[lI].dy +
                sVectorBuffer[lI].dz * sVectorBuffer[lI].dz) > 1. + dEpsilon)
            {
                bSuccess = false;
                lIndex   = lI;
                break;
            }
        }
    }
    else
    {
        // For diffusion vector sets in XYZ coordinates, no vector
        // component must exceed a value of 1
        for(long lI = 0; lI < lDirections; ++lI)
        {
            if((fabs(sVectorBuffer[lI].dx) > 1. + dEpsilon) ||
               (fabs(sVectorBuffer[lI].dy) > 1. + dEpsilon) ||
               (fabs(sVectorBuffer[lI].dz) > 1. + dEpsilon))
            {
                bSuccess = false;
                lIndex   = lI;
                break;
            }
        }
    }

    if(!bSuccess)
    {
        // error: limits exceeded
        std::stringstream sErrorMsg;
        sErrorMsg << "Vector index [" << lIndex << "/" << lDirections <<"] exceeds allowed range \n(prs: length <= 1, xyz: each component <= 1)";
        m_sErrorMessage = sErrorMsg.str();

        SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str());

        m_cLookupTable[lDirections] = 0;
        delete[] sVectorBuffer;
        return false;
    }

    // Vector set appears to be valid - copy content to internal structures
    m_lDirections       = lDirections;
    delete[] m_sVector;
    m_sVector           = new VectorStruct[lDirections];
    memcpy(m_sVector, sVectorBuffer, lDirections * sizeof(VectorStruct));
    m_eCoordinateSystem = eNewCoordinateSystem;
    m_sUserComment      = sNewComment;
    delete[] sVectorBuffer;

    return true;
}
#endif // #ifdef WIN32


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
bool DiffusionDirections::writeToFile
(
const char* pszVectorFileName /*!< File name (without path) */
)
// ===========================================================================
{
    // some macros for debugging:
#define DEBUG_writeToFile(a) // {a;}

    std::stringstream sVectorFile;

    if(m_cLookupTable[m_lDirections] == -1)
    {
        m_sErrorMessage = "It is not allowed to dump SIEMENS internal directions sets!";

        SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str());

        return false;
    }

    // Assemble filename
    char* ptCustomerSeq = NULL;
    ptCustomerSeq = getenv(VECTORFILEENV);

    if(ptCustomerSeq == NULL)
    {
        std::stringstream sErrorMsg;
        sErrorMsg << "Environment does not exist: " << VECTORFILEENV;
        m_sErrorMessage = sErrorMsg.str();

        DEBUG_writeToFile(SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str()));
        return false;
    }
    else
    {
        // e.g. "%CustomerSeq%\DiffusionVectorSets\DiffusionVectors.dvs"
        sVectorFile << ptCustomerSeq << "\\" << VECTORFILEPATH << "\\" << pszVectorFileName;
    }

    DEBUG_writeToFile(SEQ_TRACE_INFO.print("m_sVectorFileName.c_str(): %s", sVectorFile.str().c_str()));

    // open vector file
    FILE* pOutputFile = NULL;
    if((pOutputFile = fopen(sVectorFile.str().c_str(), "w")) == NULL)
    {
        std::stringstream sErrorMsg;
        sErrorMsg << "Cannot open diffusion vector set: " << sVectorFile.str();
        m_sErrorMessage = sErrorMsg.str();

        DEBUG_writeToFile(SEQ_TRACE_INFO.print("ERROR: %s", m_sErrorMessage.c_str()));
        return false;
    }

    // Get time and date
    time_t ltime;
    time(&ltime);

    // Write header
    fprintf(pOutputFile, "# -----------------------------------------------------------------------------\n");
    fprintf(pOutputFile, "#        Copyright (C) SIEMENS Healthcare GmbH 2016  All Rights Reserved.      \n");
    fprintf(pOutputFile, "# -----------------------------------------------------------------------------\n");
    fprintf(pOutputFile, "#                       \n");
    fprintf(pOutputFile, "#  Project: NUMARIS/X   \n");
    fprintf(pOutputFile, "#     File: %s          \n", sVectorFile.str().c_str());
    fprintf(pOutputFile, "#     Date: %s          \n", asctime(localtime(&ltime)));
    fprintf(pOutputFile, "#                       \n");
    fprintf(pOutputFile, "#  Descrip: External vector file for SBBDiffusion\n");
    fprintf(pOutputFile, "# -----------------------------------------------------------------------------\n\n");
    fprintf(pOutputFile, "[directions=%li]  \n", m_lDirections);

    // Write coordinate system
    if(m_eCoordinateSystem == MrProtocolData::DIFFDIR_CS_XYZ)
    {
        fprintf(pOutputFile, "CoordinateSystem = xyz\n");
    }
    else
    {
        fprintf(pOutputFile, "CoordinateSystem = prs\n");
    }

    // Write normalisation directive
    fprintf(pOutputFile, "Normalisation = none\n");

    // Write user comment (if there is one)
    const char *pszUserComment = m_sUserComment.c_str();
    if(strlen(pszUserComment) > 0)
    {
        fprintf(pOutputFile, "Comment = %s\n", pszUserComment);
    }

    // Write directions
    for(long lI = 0; lI < m_lDirections; ++lI)
    {
        fprintf(pOutputFile, "Vector[%li] = ( %f, %f, %f )\n", lI, m_sVector[lI].dx, m_sVector[lI].dy, m_sVector[lI].dz);
    }

    fclose(pOutputFile);

    return true;
}
#endif // #ifdef WIN32


// ===========================================================================
/// Normalize a complete vector set
/**
    The provided vector set will be normalized in situ
    as requested by the mode parameter. The used norm is \f$\sqrt{ x^2 + y^2 + z^2}\f$.

    \return true if everything was ok, false if a NULL vector in the vector set
    could not be normalized.
    */
bool DiffusionDirections::normalize
(
VectorStruct                         sVector[],         /**< Imp/Exp: Array containing diffusion vectors */
const long                           lDirections,       /**< Imp:     Number of diffusion directions     */
MrProtocolData::DiffDirNormalization eMode              /**< Imp:     Algorithm used for normalization   */
)
// ===========================================================================
{
    // macro for debugging:
#define DEBUG_normalize(a) // {a;}

    if(eMode == MrProtocolData::DIFFDIR_NORM_NONE)
    {
        DEBUG_normalize(cout << "Normalisation=none:" << endl; dump(););
        return true;
    }

    // Normalize each vector to one
    long   li                = 0;      // loop index
    double dNorm             = 0.0;    // norm of current vector
    double dNormVec0         = -1.0;   // expected norm of vectors (-1 means 'invalid')
    double dLargestComponent = 0.0;    // largest component in vector set
    for(li=0; li < lDirections; li++)
    {
        // Calculate norm
        dNorm = sqrt(sVector[li].dx * sVector[li].dx +
                     sVector[li].dy * sVector[li].dy +
                     sVector[li].dz * sVector[li].dz);

        // Check for null vector
        if(fAlmostEqual(dNorm, 0.0))
        {
            DEBUG_normalize(cout << "Null Vector" << endl;);
            sVector[li].dx = sVector[li].dy = sVector[li].dz = 0.0;
            return false;
        }

        // Check for common norm
        if(dNormVec0 < 0.)
        {
            dNormVec0 = dNorm;
        }
        else
        {
            if(! fAlmostEqual(dNorm, dNormVec0))
            {
                SEQ_TRACE_INFO.print("WARNING: Norm(set=%ld,dir=%ld)=%f differs from Norm(set=%ld,dir=0)=%f",
                           lDirections, li, dNorm, lDirections, dNormVec0);
            }
        }

        // Normalize vector
        sVector[li].dx /= dNorm;
        sVector[li].dy /= dNorm;
        sVector[li].dz /= dNorm;

        // Get larget component
        dLargestComponent = std::max(std::max(dLargestComponent, fabs(sVector[li].dx)),
                                     std::max(fabs(sVector[li].dy), fabs(sVector[li].dz)));
    }

    if(eMode == MrProtocolData::DIFFDIR_NORM_UNITY)
    {
        DEBUG_normalize(cout << "Normalisation=Unity:" << endl; dump(););
        return true;
    }


    DEBUG_normalize(cout << "LargestComponent=" << dLargestComponent << endl;);
    double dNormFactor = 1. / dLargestComponent;
    for(li = 0; li < lDirections; li++)
    {
        sVector[li].dx *= dNormFactor;
        sVector[li].dy *= dNormFactor;
        sVector[li].dz *= dNormFactor;
    }
    DEBUG_normalize(cout << "Normalisation=Maximum:" << endl; dump(););

    return true;
}
