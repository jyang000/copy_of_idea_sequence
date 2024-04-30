//-----------------------------------------------------------------------------
// Copyright (C) Siemens Healthcare GmbH 2019
// All Rights Reserved.  Restricted
//-----------------------------------------------------------------------------

#pragma once

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <string>
#include <list>
#include <memory>


#include "MrMeasSrv/MeasUtils/StringUtils.h" // compatibility with old clients

#include "MrGlobalDefinitions/MrResult.h"

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

#ifdef _WIN32
  /// set the minimum working set size (avoid swapping) of own process
  __IMP_EXP bool setMinWorkingSetSize(size_t minSizeMB);

  __IMP_EXP bool setProcessPriority(const std::string& sPriority);
#endif

  __IMP_EXP bool isRealTimeContext();


  // retrieve internal name of posix / xenomai objects (for reflection lib)
  __IMP_EXP void getMeasUtilsRealtimeObjectName(const char* className, char* buffer, size_t buffer_size);

  /// get brief copyright message for printouts (untranslated)
  /// uses the generated rc file on windows and a hardcoded string in linux
  __IMP_EXP std::string getCopyright();

  /// get product name for printouts (untranslated)
  /// uses the generated rc file on windows and a hardcoded default in linux
  __IMP_EXP std::string getProductName();

  /// get product version for printouts (untranslated)
  /// uses the generated rc file on windows and an empty string in linux
  __IMP_EXP std::string getProductVersion();

  /// get a part of the version resource (untranslated)
  /// uses the generated rc file on windows and an empty string in linux
  __IMP_EXP std::string queryVersionResource(std::string pcSubBlock);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// File Name Expansion
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

class __IMP_EXP FileNameExp
{
public:

  //---------------------------------------------------------------------------
  // constructor/destructor
  //---------------------------------------------------------------------------

  FileNameExp (void);
  ~FileNameExp (void);

  // Translate method. This method leaves the input string untouched
  // and reserves dynamically memory for the output (which is released
  // by the destructor). On error, a NULL pointer is returned.
  const char *Translate (const char *srcfile);

private:

  /// internal helper
  void DoReplace(size_t start, size_t end, std::string& rs);

  std::list<std::string> m_lResults;

  // Disable copy constructor and assignment - make no sense
  FileNameExp (FileNameExp& src);
  FileNameExp& operator= (FileNameExp& src);
};


#define MeasUtilsGlobalGuard() MeasUtilsModuleGuard __measUtils_guard_global__;

// This class will be used to protected
// global initialisation in MeasUtils
// from multithreading race conditions.
// It merely falls back to
// Parc::GlobalGuard which is
// using the same mutex.
// This behaviour is by design, because
// different locks will result in a higher
// chance for deadlocks...
class __IMP_EXP MeasUtilsModuleGuard
{
public:
  MeasUtilsModuleGuard();
  ~MeasUtilsModuleGuard();
};
