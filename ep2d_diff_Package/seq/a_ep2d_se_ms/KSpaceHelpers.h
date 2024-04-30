// 	---------------------------------------------------------
// 	  Copyright (C) Siemens Healthcare GmbH 2020  All Rights Reserved.
// 	---------------------------------------------------------

#pragma once

// MrProtocolData::MrProtData
#include "MrProtSrv/Domain/CoreNative/MrPat.h"
#include "MrProtSrv/Domain/CoreNative/MrFastImaging.h"
#include "MrProtSrv/Domain/CoreNative/MrCoilSelectData.h"
#include "MrProtSrv/Domain/CoreNative/MrFilter.h"
// MrProtocolData::MrProtData

// MrProt KSpace Wrapper
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrCoilInfo.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/CoilSelect/MrRxCoilSelect.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/KSpace/MrKSpace.h"

#include "MrImaging/seq/common/iPAT/iPAT.h"


namespace KSpace
{

class __IMP_EXP KSpaceHelpers
{
    private:
    
    public:
        KSpaceHelpers() = default;
        virtual ~KSpaceHelpers() =default;    

    
    static bool checkPhOSAndPEFTLen(const MrProt& rMrProt, const SeqLim& rSeqLim);

#ifdef WIN32
     static void _setPELines(MrUILinkBase* const pThis);
#endif
                             

private:
    KSpaceHelpers(const KSpaceHelpers&);
    KSpaceHelpers& operator=(const KSpaceHelpers&);
};


}


