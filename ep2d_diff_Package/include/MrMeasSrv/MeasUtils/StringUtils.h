//-----------------------------------------------------------------------------
// Copyright (C) Siemens Healthcare GmbH 2019 - 2020
// All Rights Reserved.  Restricted
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// include control
//-----------------------------------------------------------------------------

#pragma once

#ifndef MeasSrv_StringUtil_h
#define MeasSrv_StringUtil_h 1

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <locale>
#include <memory>
#include "MrGlobalDefinitions/MrStringUtil.h"
#include "MrMeasSrv/MeasCompilerDefs.h"


//-----------------------------------------------------------------------------
// Import/export control
//-----------------------------------------------------------------------------
#ifdef BUILD_MeasUtils
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"


//*****************************************************************************
// exports
//*****************************************************************************

namespace MeasSrv
{

  /// get default locale for encoding strings within meassrv
  static inline std::locale getDefaultLocale()
  {
#ifdef _WIN32
#ifdef _MSC_VER
#if (_MSC_VER >= 1700)
    return std::locale("en-US");
#else //(_MSC_VER >= 1700)
    return std::locale("us");
#endif //(_MSC_VER >= 1700)
#else //_MSC_VER
#error "needs porting"
#endif //_MSC_VER
#else //_WIN32
    return std::locale();
#endif //_WIN32
  }


  //-----------------------------------------------------------------------------
  //  Prototypes of strutils
  //-----------------------------------------------------------------------------
#ifndef WIN32
  __IMP_EXP int stricmp
  (
    const char * ptStr1,                  // Str 1
    const char * ptStr2                   // Str 2
  );

  __IMP_EXP int strnicmp
  (
    const char * ptStr1,                  // Str 1
    const char * ptStr2,                  // Str 2
    size_t       count                    // maximal number of characters to compare
  );

#endif

  //-----------------------------------------------------------------------------
  // Convert all characters of a string to the corresponding uppercase version
  // (The original pointer will always be returned)
  // keep this inline to avoid linker deps
  //-----------------------------------------------------------------------------
  static inline char * upperstr(char * ptStr)
  {
    assert(ptStr);

    //std::locale loc("en-US");
    std::locale loc = getDefaultLocale();
    for (char* pt = ptStr; *pt; pt++)
    {
      *pt = std::toupper(*pt, loc);
    }

    return ptStr;
  }

  //-----------------------------------------------------------------------------
  // Convert all characters of a string to the corresponding lowercase version
  // (The original pointer will always be returned)
  // keep this inline to avoid linker deps
  //-----------------------------------------------------------------------------
  static inline char * lowerstr(char * ptStr)
  {
    assert(ptStr);

    //std::locale loc("en-US");
    std::locale loc = getDefaultLocale();
    for (char* pt = ptStr; *pt; pt++)
    {
      *pt = std::tolower(*pt, loc);
    }

    return ptStr;
  }


  /// split strings into tokens
  __IMP_EXP std::vector<std::string> tokenize(const std::string& input, char separator);


  //-----------------------------------------------------------------------------
  //  Convert all characters of a string to the corresponding uppercase version
  //  (The original pointer will always be returned)
  //-----------------------------------------------------------------------------
  __IMP_EXP std::string& upperstr(std::string& str);

  //-----------------------------------------------------------------------------
  //  Convert all characters of a string to the corresponding lowercase version
  //  (The original pointer will always be returned)
  //-----------------------------------------------------------------------------
  __IMP_EXP std::string& lowerstr(std::string& str);

  ///
  /// printf() like std::string formatting
  ///
  template<typename ... Args>
  inline std::string formatString(const std::string& format, Args ... args)
  {
#ifdef WIN32
    size_t size = _snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
#else
    size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
#endif
    auto buf = std::make_unique<char[]>(size);
#ifdef WIN32
    _snprintf(buf.get(), size, format.c_str(), args ...);
#else
    snprintf(buf.get(), size, format.c_str(), args ...);
#endif
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
  }

  /// copies up to dstsize - 1 characters from the string src to dst, NUL-terminating the result if dstsize is not 0.
  /// If the src and dst strings overlap, the behavior is undefined.
  /// If the return value is >= dstsize, the output string has been truncated. It is the caller's responsibility to handle this.
  /// this method is inline since MakeSiemensSequence must not link againat any dyn. lib

  inline size_t strlcpy(char *dst, const char *src, size_t dstsize)
  {
    return ::syngo::MR::strlcpy(dst, src, dstsize);
  }


  /// appends string src to the end of dst. It will append at most dstsize - strlen(dst) - 1 characters. 
  /// It will then NUL-terminate, unless dstsize is 0 or the original dst string was longer than dstsize 
  /// (in practice this should not happen as it means that either dstsize is incorrect or that dst is not a proper string).
  // If the src and dst strings overlap, the behavior is undefined.
  __IMP_EXP size_t strlcat(char *dst, const char *src, size_t dstsize);

  /// unicode conversion stuff
  __IMP_EXP std::string unicode2utf8(const std::wstring& rStr);

  /// unicode conversion stuff
  __IMP_EXP std::wstring utf82unicode(const std::string& rStr);

  /// left trim whitespace with default locale (us)
  __IMP_EXP std::string ltrim(const std::string &in);

  /// rigt trim whitespace with default locale (us)
  __IMP_EXP std::string rtrim(const std::string &in);

  /// left and right trim whitespace with default locale (us)
  __IMP_EXP std::string trim(const std::string &in);

  /// as long as we do not have std::hash<std::string_view> (C++17), 
  /// provide own hash implementation
  __IMP_EXP std::uint32_t hash_pjw(const std::string&) MEAS_NO_THROW;

  /// as long as we do not have std::hash<std::string_view> (C++17), 
  /// provide own hash implementation
  __IMP_EXP std::uint32_t hash_pjw(const char* pstr) MEAS_NO_THROW;
  
}


using MeasSrv::upperstr;
using MeasSrv::lowerstr;

#ifndef WIN32
using MeasSrv::stricmp;
using MeasSrv::strnicmp;
#endif

#endif // MeasSrv_StringUtil_h


