//---------------------------------------------------------------------------------
//  Copyright (C) Siemens Healthcare GmbH 2017.  All Rights Reserved. Confidential.
//---------------------------------------------------------------------------------
//
// Project: N/X
//    File: include\MrNFramework\MrFeatureToggle\Toggle.h
//  Author: Nestler, Thomas (HC DI MR R&D MC SFI)
//    Date: n.a.
//
//    Lang: C++
//
//---------------------------------------------------------------------------------


#pragma once

#include <memory>

#ifdef BUILD_MrFeatureToggle
#define __OWNER
#endif
#include <MrGlobalDefinitions/ImpExpCtrl.h>


namespace MrFeatureToggle
{
    class CToggleAPICtxt;

    /// \brief Toggle-instances are used as client-API to access the state of a toggle (ON or OFF)
    ///              

    /*
        usage example:

        //-------------------------------------------------------------------------------------------------------------------------
        // file f1.h:
        //-------------------------------------------------------------------------------------------------------------------------

        #define TOGGLE_OPTIMIZED_ALGORITHM_F1 "EnableOptimizedAlgorithm_f1"
        int f1();

        //-------------------------------------------------------------------------------------------------------------------------
        // file main.cpp:
        //-------------------------------------------------------------------------------------------------------------------------

        #include "f1.h"
        #include <MrNFramework/MrFeatureToggle/Toggle.h>
        #include <iostream>

        int main2(int argc, char* argv[])
        {
            try
            {
                // InitGuard to optimize lifecycle-/runtime-behaviour (optional)
                MrFeatureToggle::InitGuard pinFeatureToggleFramework; // can throw, e.g. if toggle configuration not found
                return f1();
            }
            catch (const std::exception& e)
            {
                std::cerr << "exception occurred: " << e.what() << std::endl;
                return 1;
            }
        }

        //-------------------------------------------------------------------------------------------------------------------------
        // file f1.cpp:
        //-------------------------------------------------------------------------------------------------------------------------

        #include <MrNFramework/MrFeatureToggle/Toggle.h>

        int f1()
        {
            MrFeatureToggle::Toggle f1Alogo(TOGGLE_OPTIMIZED_ALGORITHM_F1); // can throw, e.g. if this toggle not configured 
            if (f1Alogo.isOn())
            {
                // calculate with optimized algorithm
                return 3;
            }
            else
            {
                // calculate with default algorithm
                return 2;
            }
        }

        //-------------------------------------------------------------------------------------------------------------------------
    */

    class __IMP_EXP Toggle
    {
    public:
        Toggle(const char* psToggleName);
        Toggle(const Toggle&) = default;
        virtual ~Toggle();

        Toggle& operator=(const Toggle& r) = default;

        // IToggle
        bool isOn() const;
        const char* name() const;

    private:
        std::shared_ptr<CToggleAPICtxt> m_spCtx;
    };

    /// \brief Deprecated: Scoped Initialization of the toggle-library
    class InitGuard
    {
    public:
        InitGuard() {}
        InitGuard(const InitGuard&) = default;
        ~InitGuard() {}
        InitGuard& operator=(const InitGuard&) = default;
    };
}