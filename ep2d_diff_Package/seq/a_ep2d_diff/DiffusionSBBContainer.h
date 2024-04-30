//-----------------------------------------------------------------------------
// <copyright file="DiffusionSBBContainer.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
//-----------------------------------------------------------------------------

/*
\class DiffusionSBBContainer

\author Michael.Zwanger@med.siemens.de, Uvo.Hoelscher@siemens.com

\brief This class is the interface for the sequence to access all
diffusion functionality.

This is a class factory, which creates an instance of the requested diffusion
class (which is derived from the SBBDiffusion_Base base class).
For example, if the MrProt protocol requests a trace diffusion weighting,
the class factory creates an instance of the Diffusion_Trace class and
returns a pointer to this object.

Within the sequence, this class is used as the common interface
for all access to SBBDiffusion functionality.This is achieved
by overloading the -> operator.

*/

#pragma once

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"
#include "MrImaging/seq/a_ep2d_diff/SBBDiffusion_Base.h"

#ifdef BUILD_SEQU
#define __OWNER
#endif

#include "MrGlobalDefinitions/ImpExpCtrl.h"     // import/export control



// ---------------------------------------------------------------------------
// Type definitions
// ---------------------------------------------------------------------------


/// This enumeration provides the different sequence schemes
/** Each sequence scheme can be directly mapped to one class implementing the
corresponding diffusion mode.
The first value has index 1 because this enumeration can also be used for
a WIP parameter.
*/
enum EnumDiffusionScheme {
    DiffusionSchemeNone = 1,
    DiffusionSchemeBipolar,
    DiffusionSchemeTrace,
    DiffusionSchemeStejskal,
    DiffusionSchemeStejskalPlus,
    DiffusionSchemeBipolarPlus,
    DiffusionSchemeSTEAM,
};


namespace SEQ_NAMESPACE
{

    // ===========================================================================
    class __IMP_EXP DiffusionSBBContainer
        // ===========================================================================
    {

    public:
        DiffusionSBBContainer();
        virtual ~DiffusionSBBContainer();

        bool create(MrProt &rMrProt);
        
        bool getDiffusionSeqLims(SeqLim &rSeqLim);

        // This getter returns the pointer to the real SBB
        SBBDiffusion_Base * operator -> ();
        SBBDiffusion_Base* getpSBB();

        /// This is the pointer to the current diffusion scheme class
        SBBDiffusion_Base * m_pSBBDiffusion_Base{nullptr};

    private:

        void free();

        ///    This identifier is used by the trace messages.
        const char* m_tIdent{"DiffusionSBBContainer"};

        ///    This is the diffusion encoding sequence scheme currently used.
        EnumDiffusionScheme m_eDiffusionScheme{DiffusionSchemeNone};

        /// Actual patient direction (from MeasPatient)
        int m_iPatDirection{0};

        /// Actual patient position (from MeasPatient)
        int m_iPatPosition{0};
    };

} //end of namespace SEQ_NAMESPACE
