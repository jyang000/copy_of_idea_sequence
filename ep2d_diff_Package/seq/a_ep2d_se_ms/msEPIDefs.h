#pragma once

#include <string>
#include <sstream>
#include <cstdlib>      // std::getenv()
#include <cstdint>
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/CoreNative/MrWipMemBlock.h"
//#ifndef SEQ_NAMESPACE // included from ICE
//#define SEQ_NAMESPACE WIP_NAMESPACE_DECL
//#endif

namespace SEQ_NAMESPACE
{


    // ---------------------------------------------------------------------------
    // Definitions
    // ---------------------------------------------------------------------------

// Ice program parameter indices (for MDH)
#define ICEPROGRAMPARA_SHOT_INDEX      ( MDH_NUMBEROFICEPROGRAMPARA - 1 )
// Currently (XA20A), up to 24 parameters are supported (see MDH_NUMBEROFICEPROGRAMPARA, \MrServers\MrMeasSrv\SeqIF\MDH\mdh.h)

// Additional Ice program parameters
#define WIP_ICE_PARA_EXPORT_TRAININGDATA ICE_PROGRAM_PARA_USER_3    // Reserved  (enables training data export in offline processing)
#define WIP_ICE_PARA_EXPORT_IMAGE        ICE_PROGRAM_PARA_USER_4    // Reserved  (enables image         export in offline processing)

// @DEBUG
#define WIP_ICE_PARA_EXPORT_SNR_PARAM   ICE_PROGRAM_PARA_USER_5    // Reserved  (enables export of SNR parameters)
#define WIP_ICE_PARA_USE_AUTOREG        ICE_PROGRAM_PARA_USER_6    // Reserved  (enables Auto-Regularization)


    /// Helper function to identify network reconstruction modes
    inline bool isDLReconMode(const MrProtData& rMrProt)
    {
        return ((rMrProt.getMrAdvancedReconstruction().getlAdvancedReconstructionMode() == AdvancedReconstructionMode::AdvancedReconstructionMode_IterativeDenoising)
            && (rMrProt.getMrAdvancedReconstruction().getlDenoisingMethod() == DenoisingMethod::DeepResolve_FastBrain));
    }

    inline bool isDLReconMode(const MrProt& rMrProt)
    {
        return ((rMrProt.getMrAdvancedReconstruction().getlAdvancedReconstructionMode() == AdvancedReconstructionMode::AdvancedReconstructionMode_IterativeDenoising)
            && (rMrProt.getMrAdvancedReconstruction().getlDenoisingMethod() == DenoisingMethod::DeepResolve_FastBrain));
    }


} //end namespace SEQ_NAMESPACE

