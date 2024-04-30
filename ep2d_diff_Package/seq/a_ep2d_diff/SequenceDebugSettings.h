/*! 
***************************************************************************
\file   SequenceDebugSettings.h 

\brief  Class for storing all debug Settings

\author Uvo Hoelscher, Ralf Kartaeusch

\b Language: C++

\b Copyright: &copy; Siemens Healthineers 2018
All rights reserved.   

***************************************************************************

*/


#ifndef SequenceDebugSettings_h
#define SequenceDebugSettings_h 


// system support
#include <boost/spirit/home/support/detail/hold_any.hpp>
#include <boost/any.hpp>
#include <map>
#include <string>
#include <list>
#include <memory>
#include "MrCommon/UConfig/ConfigGuard.h"


#include "MrMeasSrv/MeasUtils/MeasUtils.h"
#include "MrCommon/UConfig/KeyValueAccess.h"
#include "MrCommon/UConfig/CompConfig.h"
#include "MrCommon/UConfig/ExtPtr.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrImagingFW/libSeqSysProp/SysProperties.h"
// system support



namespace SequenceDebugSettings
{
	using namespace UConfig;
	using namespace std;

	class SequenceDebugSettings
	{
	public:

		SequenceDebugSettings(string sParamNameUseDebugSettingsGeneral)
		{
			m_sParamNameUseDebugSettings = sParamNameUseDebugSettingsGeneral;
			// read settings
			if (ReadXMLSettings()) {
				readDebugSettings();
			}
			else {
				SEQ_TRACE_ALWAYS << "Error while reading SequenceDebugSettings.";
				m_bUseDebugSettings = false;
			}
		}
		
		

		template< typename T>
		T getDefaultSetting(std::string settingName, T defaultValue){
			if (m_bUseDebugSettings)
			{
				
				// real all settings from config file
				using namespace UConfig;
				ICompConfig* pConfig = m_pGuard->getCompConfig();
				IKeyValueAccess* pRoot = ExtPtr<IKeyValueAccess>(pConfig).get();

				if (pRoot->isAvailable(settingName.c_str()))
				{
                    T value = KeyValueAccess::get(pRoot, settingName.c_str(), defaultValue);
					SEQ_TRACE_ALWAYS << "Read Setting " << settingName << " with value: "<< value;
					return value;
				}
			
			}
			return defaultValue;
		}

		bool areDebugSettingsUsed()
		{
			return m_bUseDebugSettings;
		}

	protected:



		
		bool ReadXMLSettings()
		{
			FileNameExp aFileNameExp;
			std::string path = aFileNameExp.Translate("%MRIPRODUCT%/MrImaging/SequenceSettings_General.xml");

			using namespace UConfig;

			m_pGuard = std::make_shared<CConfigGuard>(nullptr, path.c_str());

			ICompConfig* pConfig = m_pGuard->getCompConfig();
			if (pConfig == nullptr)
				return false;
			pConfig->loadData(false); // do not create file if not exists
			IKeyValueAccess* pRoot = ExtPtr<IKeyValueAccess>(pConfig).get();
			pRoot->setAllowDefaultValues(true); // prevent Exception Trace
			return true;
		}

		void readDebugSettings()
		{


			using namespace UConfig;
			ICompConfig* pConfig = m_pGuard->getCompConfig();
			IKeyValueAccess* pRoot = ExtPtr<IKeyValueAccess>(pConfig).get();

			if (pRoot->isAvailable(m_sParamNameUseDebugSettings.c_str()))
			{
                m_bUseDebugSettings = KeyValueAccess::get(pRoot, m_sParamNameUseDebugSettings.c_str(), false);
			}
			else {
				m_bUseDebugSettings = false;
			}
		}


		


		// shall the debug settings be used in general? This is toggled by an dedicated bool
		bool m_bUseDebugSettings;
		
		// general parameter for switching debug settings on and off
		std::string m_sParamNameUseDebugSettings;

		std::shared_ptr<UConfig::CConfigGuard> m_pGuard;

	private:
		SequenceDebugSettings() {};

	};
}

#endif 