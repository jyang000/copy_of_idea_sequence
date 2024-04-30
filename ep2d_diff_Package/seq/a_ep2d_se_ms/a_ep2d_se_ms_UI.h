//----------------------------------------------------------------------------------
// <copyright file="a_ep2d_se_ms_UI.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2021. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#include "MrImaging/seq/a_ep2d_se/a_ep2d_se_UI.h"
#include "MrImaging/seq/common/msEPI/a_ep_SegmentationUI.h"
#include "MrImaging/seq/common/msEPI/a_NetworkUI.h"


namespace SEQ_NAMESPACE
{

class Ep2d_se_ms_UI
    : public Ep2d_se_UI
    , public EpSegmentationUI
    , public NetworkUI
{
  public:
    Ep2d_se_ms_UI() = default;

    virtual ~Ep2d_se_ms_UI() = default;
    
    NLS_STATUS initializeUI(MrProt& rMrProt, SeqLim& rSeqLim) override;

    NLS_STATUS registerUI(SeqLim& rSeqLim) override;

#ifdef WIN32
    UI_ELEMENT_SELECTION m_RFPulseType_se;
    UI_ELEMENT_SELECTION m_FatSatMode_se;
    UI_ELEMENT_SELECTION m_FlowAttenuation_se;

    UI_ELEMENT_LONG m_Contrasts_se;
#endif
};

}//end of namespace SEQ_NAMESPACE
