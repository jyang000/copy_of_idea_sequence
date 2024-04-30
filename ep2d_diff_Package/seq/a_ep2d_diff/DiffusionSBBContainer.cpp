//-----------------------------------------------------------------------------
// <copyright file="DiffusionSBBContainer.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
//-----------------------------------------------------------------------------

#include "DiffusionSBBContainer.h"
#include "MrMeasSrv/MeasPatient/MeasPatient.h"
#include "SBBDiffusion_Trace.h"
#include "SBBDiffusion_Bipolar.h"
#include "SBBDiffusion_Stejskal.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"

namespace SEQ_NAMESPACE
{
    // ===========================================================================
    ///  The constructor of this class does nothing.
    // ===========================================================================
    DiffusionSBBContainer::DiffusionSBBContainer()
    {
        // get registered patient position and direction
        MeasPatient theMeasPatient;
        theMeasPatient.getDirection(&m_iPatDirection);
        theMeasPatient.getPosition(&m_iPatPosition);

        // create dummy diffusion SBB
        m_pSBBDiffusion_Base = new SBBDiffusion_Stejskal(nullptr);
    }


    // ===========================================================================
    ///  The destructor frees allocated memory.
    // ===========================================================================
    DiffusionSBBContainer::~DiffusionSBBContainer()
        // ===========================================================================

    {
        free();
    }



    // ===========================================================================
    ///	Creates a pointer to the diffusion mode selected in the protocol.
    /**
    \note If the Diffusion_Bipolar (or Diffusion_Stejskal resp.) class is selected,
    the member attribute Diffusion_Bipolar::Mode will also be initialized
    (due to a friend relation).
    */
    // ===========================================================================
    bool DiffusionSBBContainer::create(MrProt &rMrProt)
        // ===========================================================================
    {
        SEQ::DiffusionMode   Mode;                                      // the diffusion mode as requested in MrProt
        SEQ::DiffusionScheme Scheme;                                    // the diffusion scheme as requested in MrProt

        EnumDiffusionScheme eDiffusionScheme = DiffusionSchemeBipolar;  //the requested sequence scheme with dummy initialization

        Mode   = rMrProt.diffusion().getulMode();
        Scheme = rMrProt.diffusion().getdsScheme();

        switch(Mode)
        {
            case SEQ::DIFFMODE_ONE_SCAN_TRACE:
                eDiffusionScheme = DiffusionSchemeTrace;
                break;
            case SEQ::DIFFMODE_TENSOR:
            case SEQ::DIFFMODE_QSPACE:
            case SEQ::DIFFMODE_FREE:
            case SEQ::DIFFMODE_THREE_SCAN_TRACE:
            case SEQ::DIFFMODE_FOUR_SCAN_TRACE:
            case SEQ::DIFFMODE_ORTHOGONAL:
            case SEQ::DIFFMODE_SLICE:
            case SEQ::DIFFMODE_READ:
            case SEQ::DIFFMODE_PHASE:
            case SEQ::DIFFMODE_DIAGONAL:
                switch(Scheme)
                {
                    case SEQ::DIFFSCHEME_BIPOLAR:
                        eDiffusionScheme = DiffusionSchemeBipolar;
                        break;
                    case SEQ::DIFFSCHEME_MONOPOLAR:
                        eDiffusionScheme = DiffusionSchemeStejskal;
                        break;
#ifdef WIP
                    case SEQ::DIFFSCHEME_BIPOLAR_PLUS:
                        eDiffusionScheme = DiffusionSchemeBipolar;
                        break;
                    case SEQ::DIFFSCHEME_MONOPOLAR_PLUS:
                        eDiffusionScheme = DiffusionSchemeStejskal;
                        break;
                    case SEQ::DIFFSCHEME_STEAM:
                        eDiffusionScheme = DiffusionSchemeStejskal;
                        break;
#endif //   #ifdef WIP
                    default:
                        SEQ_TRACE_ERROR.print("ERROR: Unknown scheme %d", int(Scheme));
                        return false;
                }
                break;
            default:
                SEQ_TRACE_ERROR.print("ERROR: Unknown mode %d", int(Mode));
                return false;
        }

        if(m_pSBBDiffusion_Base &&  eDiffusionScheme == m_eDiffusionScheme)
        {
            // The current diffusion class is still valid.
            return true;
        }

        // free old class
        free();

        // create new diffusion class :
        switch(eDiffusionScheme)
        {
            case DiffusionSchemeTrace:
                m_pSBBDiffusion_Base = new SBBDiffusion_Trace(nullptr);
                break;
            case DiffusionSchemeBipolar:
                m_pSBBDiffusion_Base = new SBBDiffusion_Bipolar(nullptr);
                break;
            case DiffusionSchemeStejskal:
                m_pSBBDiffusion_Base = new SBBDiffusion_Stejskal(nullptr);
                break;
#ifdef WIP
            case DiffusionSchemeBipolarPlus:
                m_pSBBDiffusion_Base = new Diffusion_BipolarPlus();
                break;
            case DiffusionSchemeStejskalPlus:
                m_pSBBDiffusion_Base = new Diffusion_StejskalPlus();
                break;
            case DiffusionSchemeSTEAM:
                m_pSBBDiffusion_Base = new Diffusion_STEAM();
                break;
#endif //   #ifdef WIP
            default:
                SEQ_TRACE_ERROR.print("ERROR: Unknown DiffusionScheme %d", static_cast<int>(Scheme));
                return false;
        }

        m_eDiffusionScheme = eDiffusionScheme;

        // Pass patient position and direction to new diffusion mode
        m_pSBBDiffusion_Base->setPatPosDir(m_iPatDirection, m_iPatPosition);

        // Pass current diffusion vector set filename to new diffusion mode:
        // Extract first line from comment (which holds the filename only)
        std::string sFileName(rMrProt.diffusion().getsFreeDiffusionData().getsComment());
        m_pSBBDiffusion_Base->getDidiPointer()->setVectorFileName(sFileName.substr(0, sFileName.find_first_of('\n')).c_str());

        return true;
    }


    // ===========================================================================
    ///  Destroy the current diffusion mode class and free memory.
    void DiffusionSBBContainer::free()
        // ===========================================================================
    {

        if(m_pSBBDiffusion_Base)
        {
            delete m_pSBBDiffusion_Base;
            m_pSBBDiffusion_Base = nullptr;
        }
    }


    // ===========================================================================
    ///	The -> operator is overloaded to access easily the objects instanciated by this SBBDiffusion class factory.
    /**	Without having to care about the diffusion mode actually active, you can
    just write in your sequence code e.g.:
    \code
    if (! Diff->prep (rMrProt, rSeqLim, rSeqExpo))
    return( Diff->getNLSStatus());
    \endcode

    This operator can only be used after it has been initialized by
    SBBDiffusion::create(). Violations will result in an unhandled
    exception.
    */


    // ===========================================================================
    SBBDiffusion_Base * DiffusionSBBContainer::operator -> ()
        // ===========================================================================
    {
        assert(m_pSBBDiffusion_Base);
        /* If you have been trapped by this assertion, then there might be an error
        like this in your calling sequence:
        You are accessing any Diff->... method prior to having called Diff.create().
        */
        return (m_pSBBDiffusion_Base);
    }


    // ===========================================================================
    ///	Set all diffusion related seqLims.
    /*
    The diffusion mode FREE is only available if the external vector file
    VECTORFILENAME could be found.
    */
    // ===========================================================================
    bool DiffusionSBBContainer::getDiffusionSeqLims(SeqLim &rSeqLim)
        // ===========================================================================
    {
        rSeqLim.setNoOfDiffWeightings(1, K_B_VALUE_MAX_SIZE, 1, 1); 

        rSeqLim.setBValue(0, 10000, SPOILER_THRESHOLD, 0);

        // ATTENTION: The increment should be the same as SPOILER_THRESHOLD.
        // Why?
        // For trace mode, no trace diffusion gradient train is applied, but only the spoiler
        // if the b value < SPOILER_THRESHOLD.
        // For bipolar mode, the spoiler will be sent additionally to the bipolar train
        // in this case.

        rSeqLim.setDiffusionDirections(0, MAX_DIRECTIONS, 1, 1);
        rSeqLim.setDiffusionScheme(SEQ::DIFFSCHEME_BIPOLAR, SEQ::DIFFSCHEME_MONOPOLAR);

        // Availability of FREE diffusion mode is determined within GetOption handler
        rSeqLim.setDiffusionMode(SEQ::DIFFMODE_READ,
                                 SEQ::DIFFMODE_TENSOR, SEQ::DIFFMODE_QSPACE, SEQ::DIFFMODE_FREE,
                                 SEQ::DIFFMODE_ORTHOGONAL, SEQ::DIFFMODE_THREE_SCAN_TRACE, SEQ::DIFFMODE_FOUR_SCAN_TRACE,
                                 SEQ::DIFFMODE_READ, SEQ::DIFFMODE_SLICE, SEQ::DIFFMODE_PHASE,
                                 SEQ::DIFFMODE_DIAGONAL, SEQ::DIFFMODE_ONE_SCAN_TRACE);

        // Limit max. noise level to 1000 (CHARM 308923):
        rSeqLim.setNoiseLevel(0, 1000, 1, 40);

        // Q-space related parameters
        rSeqLim.setDiffQSpaceCoverage(SEQ::DIFF_QSPACE_COVERAGE_HALF, SEQ::DIFF_QSPACE_COVERAGE_FULL);
        rSeqLim.setDiffQSpaceSampling(SEQ::DIFF_QSPACE_SAMPLING_CARTESIAN);

        // Do not display the following UI parameters:
        rSeqLim.getDiffQSpaceSampling().setDisplayMode(SEQ::DM_OFF);

        return true;
    }

    SBBDiffusion_Base* DiffusionSBBContainer::getpSBB()
    {
        // this pointer always have to point to a non nullptr pointer
        assert(m_pSBBDiffusion_Base);
        return (m_pSBBDiffusion_Base);
    }

} //NAMESPACE