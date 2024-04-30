//--------------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2018  All Rights Reserved.  Confidential
//--------------------------------------------------------------------------------
#pragma once
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"

namespace IterativeDenoisingUIParameter
{
    bool areSwiftBrainCoilSettingsSatisfied(const MrProt& rMrProt);
    void initialize(SeqLim& rSeqLim);
    bool ProtocolOkForIterativeDenoising(const MrProt& rMrProt, const SeqLim& rSeqLim);
    void activateNoZeroFillingAndIcePATIfNeeded(const MrProt& rMrProt, SeqExpo& rSeqExpo);
    bool SeqExpoOkForIterativeDenoising(const MrProt& rMrProt, const SeqExpo& rSeqExpo, const SeqLim& rSeqLim);
}
