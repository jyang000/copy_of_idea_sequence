//----------------------------------------------------------------------------------
// <copyright file="a_ep2d_fid_ms_UI.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2021. All Rights Reserved. Confidential.
// </copyright>
//----------------------------------------------------------------------------------

#pragma once

#include "MrImaging/seq/a_ep2d_fid/a_ep2d_fid_UI.h"
#include "MrImaging/seq/common/msEPI/a_ep_SegmentationUI.h"
#include "MrImaging/seq/common/msEPI/a_NetworkUI.h"

namespace SEQ_NAMESPACE
{
class Ep2d_fid_ms_UI
    : public Ep2d_fid_UI
    , public EpSegmentationUI
    , public NetworkUI
{
  public:
    Ep2d_fid_ms_UI() = default;

    virtual ~Ep2d_fid_ms_UI() = default;

    NLS_STATUS initializeUI(MrProt& rMrProt, SeqLim& rSeqLim) override;

    NLS_STATUS registerUI(SeqLim& rSeqLim) override;

#ifdef WIN32
    UI_ELEMENT_SELECTION m_RFPulseType_fid;
    UI_ELEMENT_SELECTION m_FatSatMode_fid;

    UI_ELEMENT_BOOL m_FilterNorm_fid;
    UI_ELEMENT_BOOL m_FilterNormPreScan_fid;
    UI_ELEMENT_BOOL m_FilterNormBific_fid;
    UI_ELEMENT_BOOL m_FilterDistCorr_fid;
    UI_ELEMENT_BOOL m_FilterImage_fid;
#endif
};

} // end of namespace SEQ_NAMESPACE
