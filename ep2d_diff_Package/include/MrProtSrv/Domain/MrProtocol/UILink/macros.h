//  ---------------------------------------------------------
//    Copyright (C) Siemens AG 2004  All Rights Reserved.
//  ---------------------------------------------------------
//
//   Project: NUMARIS/4
//      File: \n4\pkg\MrServers\MrProtSrv\MrProtocol\UILink\macros.h 
//   Version: 
//    Author: Comp_ProBe 
//      Date: n.a.
//
//      Lang: C++
//
//   Descrip: macros
//
//   Classes: 
//
//  ---------------------------------------------------------

#ifndef _UILINK_MACROS
#define _UILINK_MACROS

// Some macros to support integration of new parameters

#define DECLARE_STD_UILINK_HANDLERS_NUM(ParamName, ParamType, LinkType, LimitType) \
    unsigned    fUILink##ParamName##GetLabelId    (LinkType* const pThis, char* arg_list[], int32_t pos);\
    unsigned    fUILink##ParamName##GetUnitId     (LinkType* const pThis, char* arg_list[], int32_t pos);\
    unsigned    fUILink##ParamName##GetToolTipId  (LinkType* const pThis, char* arg_list[], int32_t pos);\
    bool        fUILink##ParamName##IsAvailable   (LinkType* const pThis, int32_t pos);\
    bool        fUILink##ParamName##IsVisible     (MrUILinkBase* const pThis, int32_t pos);\
    ParamType   fUILink##ParamName##GetValue      (LinkType* const pThis, int32_t pos);\
    ParamType   fUILink##ParamName##SetValue      (LinkType* const pThis, ParamType val, int32_t pos);\
    bool        fUILink##ParamName##GetLimits     (LinkType* const pThis, std::vector<LimitType>&, uint32_t& rulVerify, int32_t pos);

#define DEFINE_EMPTY_STD_UILINK_HANDLERS_NUM(ParamName, ParamType, LinkType, LimitType) \
    unsigned    fUILink##ParamName##GetLabelId    (LinkType* const, char* [], int32_t)         { return MRI_STD_EMPTY;}\
    unsigned    fUILink##ParamName##GetUnitId     (LinkType* const, char* [], int32_t)         { return MRI_STD_EMPTY;}\
    unsigned    fUILink##ParamName##GetToolTipId  (LinkType* const, char* [], int32_t)         { return MRI_STD_EMPTY;}\
    bool        fUILink##ParamName##IsAvailable   (LinkType* const, int32_t)                   { return true;}\
    bool        fUILink##ParamName##IsVisible     (MrUILinkBase* const, int32_t)                   { return true;}\
    ParamType   fUILink##ParamName##GetValue      (LinkType* const, int32_t)                   { return ParamType();}\
    ParamType   fUILink##ParamName##SetValue      (LinkType* const, ParamType, int32_t)        { return ParamType();}\
    bool        fUILink##ParamName##GetLimits     (LinkType* const, std::vector<LimitType>& , unsigned int32_t&, int32_t) { return false;}

#define DECLARE_STD_UILINK_HANDLERS_SEL(ParamName, ParamType, LinkType, LimitType) \
    unsigned    fUILink##ParamName##GetLabelId    (LinkType* const pThis, char* arg_list[], int32_t pos);\
    unsigned    fUILink##ParamName##GetUnitId     (LinkType* const pThis, char* arg_list[], int32_t pos);\
    unsigned    fUILink##ParamName##GetToolTipId  (LinkType* const pThis, char* arg_list[], int32_t pos);\
    bool        fUILink##ParamName##IsAvailable   (LinkType* const pThis, int32_t pos);\
    bool        fUILink##ParamName##IsVisible     (MrUILinkBase* const pThis, int32_t pos);\
    ParamType   fUILink##ParamName##GetValue      (LinkType* const pThis, int32_t pos);\
    ParamType   fUILink##ParamName##SetValue      (LinkType* const pThis, ParamType val, int32_t pos);\
    bool        fUILink##ParamName##GetOptions    (LinkType* const pThis, std::vector<LimitType>&, uint32_t& rulVerify, int32_t pos);

#define DEFINE_EMPTY_STD_UILINK_HANDLERS_SEL(ParamName, ParamType, LinkType, LimitType) \
    unsigned    fUILink##ParamName##GetLabelId    (LinkType* const, char* [], int32_t)         { return MRI_STD_EMPTY;}\
    unsigned    fUILink##ParamName##GetUnitId     (LinkType* const, char* [], int32_t)         { return MRI_STD_EMPTY;}\
    unsigned    fUILink##ParamName##GetToolTipId  (LinkType* const, char* [], int32_t)         { return MRI_STD_EMPTY;}\
    bool        fUILink##ParamName##IsAvailable   (LinkType* const, int32_t)                   { return true;}\
    bool        fUILink##ParamName##IsVisible     (MrUILinkBase* const, int32_t)                   { return true;}\
    ParamType   fUILink##ParamName##GetValue      (LinkType* const, int32_t)                   { return MRI_STD_EMPTY;}\
    ParamType   fUILink##ParamName##SetValue      (LinkType* const, ParamType, int32_t)        { return MRI_STD_EMPTY;}\
    bool        fUILink##ParamName##GetOptions    (LinkType* const, std::vector<LimitType>& , uint32_t&, int32_t) { return false;}

#define DECLARE_STD_UILINK_HANDLERS_DOUBLE(     ParamName ) DECLARE_STD_UILINK_HANDLERS_NUM(ParamName, double,  LINK_DOUBLE_TYPE,   MrLimitDouble)
#define DECLARE_STD_UILINK_HANDLERS_LONG(       ParamName ) DECLARE_STD_UILINK_HANDLERS_NUM(ParamName, int32_t, LINK_LONG_TYPE,     MrLimitLong)
#define DECLARE_STD_UILINK_HANDLERS_BOOL(       ParamName ) DECLARE_STD_UILINK_HANDLERS_SEL(ParamName, bool,    LINK_BOOL_TYPE,     unsigned)
#define DECLARE_STD_UILINK_HANDLERS_SELECTION(  ParamName ) DECLARE_STD_UILINK_HANDLERS_SEL(ParamName, unsigned,LINK_SELECTION_TYPE,unsigned)


#define DEFINE_EMPTY_STD_UILINK_HANDLERS_DOUBLE(    ParamName ) DEFINE_EMPTY_STD_UILINK_HANDLERS_NUM(ParamName, double,  LINK_DOUBLE_TYPE,   MrLimitDouble)
#define DEFINE_EMPTY_STD_UILINK_HANDLERS_LONG(      ParamName ) DEFINE_EMPTY_STD_UILINK_HANDLERS_NUM(ParamName, int32_t, LINK_LONG_TYPE,     MrLimitLong)
#define DEFINE_EMPTY_STD_UILINK_HANDLERS_BOOL(      ParamName ) DEFINE_EMPTY_STD_UILINK_HANDLERS_SEL(ParamName, bool,    LINK_BOOL_TYPE,     unsigned)
#define DEFINE_EMPTY_STD_UILINK_HANDLERS_SELECTION( ParamName ) DEFINE_EMPTY_STD_UILINK_HANDLERS_SEL(ParamName, unsigned,LINK_SELECTION_TYPE,unsigned)

#define REGISTER_STD_UILINK_HANDLERS_NUM(ParamName, LinkType, Tag, pSeqLim, container, pos) \
        if( LinkType* p = _create<LinkType>(pSeqLim , Tag, container, pos) )\
        {\
            p->registerGetLabelIdHandler(fUILink##ParamName##GetLabelId);\
            p->registerGetToolTipIdHandler(fUILink##ParamName##GetToolTipId);\
            p->registerGetUnitIdHandler(fUILink##ParamName##GetUnitId);\
            p->registerIsAvailableHandler(fUILink##ParamName##IsAvailable);\
            p->registerGetValueHandler(fUILink##ParamName##GetValue);\
            p->registerSetValueHandler(fUILink##ParamName##SetValue);\
            p->registerGetLimitsHandler(fUILink##ParamName##GetLimits);\
        }

#define REGISTER_STD_UILINK_HANDLERS_SEL(ParamName, LinkType, Tag, pSeqLim, container, pos) \
        if( LinkType* p = _create<LinkType>(pSeqLim , Tag, container, pos) )\
        {\
            p->registerGetLabelIdHandler(fUILink##ParamName##GetLabelId);\
            p->registerGetToolTipIdHandler(fUILink##ParamName##GetToolTipId);\
            p->registerGetUnitIdHandler(fUILink##ParamName##GetUnitId);\
            p->registerIsAvailableHandler(fUILink##ParamName##IsAvailable);\
            p->registerGetValueHandler(fUILink##ParamName##GetValue);\
            p->registerSetValueHandler(fUILink##ParamName##SetValue);\
            p->registerGetOptionsHandler(fUILink##ParamName##GetOptions);\
        }

#define REGISTER_STD_UILINK_HANDLERS_DOUBLE(    ParamName, Tag, pSeqLim, container, pos)    REGISTER_STD_UILINK_HANDLERS_NUM(ParamName, LINK_DOUBLE_TYPE,   Tag, pSeqLim, container, pos)
#define REGISTER_STD_UILINK_HANDLERS_LONG(      ParamName, Tag, pSeqLim, container, pos)    REGISTER_STD_UILINK_HANDLERS_NUM(ParamName, LINK_LONG_TYPE,     Tag, pSeqLim, container, pos)
#define REGISTER_STD_UILINK_HANDLERS_BOOL(      ParamName, Tag, pSeqLim, container, pos)    REGISTER_STD_UILINK_HANDLERS_SEL(ParamName, LINK_BOOL_TYPE,     Tag, pSeqLim, container, pos)
#define REGISTER_STD_UILINK_HANDLERS_SELECTION( ParamName, Tag, pSeqLim, container, pos)    REGISTER_STD_UILINK_HANDLERS_SEL(ParamName, LINK_SELECTION_TYPE,Tag, pSeqLim, container, pos)

#define TRACE_UILINK_PARAM_INFO(_pThis_, _stream_) \
    if (!_pThis_->isWithinBinarySearch()) \
		    { \
        UTRACE_STREAM(Always, 0x01) <<  _stream_; \
		    } 
#endif