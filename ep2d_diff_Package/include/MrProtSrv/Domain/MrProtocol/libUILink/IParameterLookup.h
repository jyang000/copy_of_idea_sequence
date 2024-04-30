#ifndef IParameterLookup_h
#define IParameterLookup_h

#include <string>

class MrUILinkBase;
template<typename>
class MrUILinkSelection;

#ifdef BUILD_libUILink
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

class __IMP_EXP IParameterLookup
{
protected:
	IParameterLookup() = default;
	IParameterLookup(const IParameterLookup&) = default;
	IParameterLookup(IParameterLookup&&) = default;
	IParameterLookup& operator=(const IParameterLookup&) = default;
	IParameterLookup& operator=(IParameterLookup&&) = default;

public:
	virtual ~IParameterLookup() = default;
	virtual MrUILinkBase* Search(std::string nameTag) const = 0;
	virtual MrUILinkBase* SearchMember(std::string parent, std::string member) const = 0;
};

namespace zelda
{
	template<typename T>
	T* Editable(T* parameter, int position)
	{
		if (nullptr == parameter
			|| false == parameter->isEditable(position))
		{
			return nullptr;
		}

		return parameter;
	}

	template<typename T>
	T* Available(T* parameter, int position)
	{
		if (nullptr == parameter
			|| false == parameter->isAvailable(position))
		{
			return nullptr;
		}

		return parameter;
	}

	template<typename T>
	T* Type(MrUILinkBase* parameter)
	{
		return dynamic_cast<T*>(parameter);
	}
	
	extern template __IMP_EXP MrUILinkSelection<unsigned>* Type<MrUILinkSelection<unsigned>>(MrUILinkBase* parameter);
}

#endif