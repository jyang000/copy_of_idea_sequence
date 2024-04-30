#pragma once

#include "MrProtSrv/Domain/MrProtocol/libUILink/UILinkBase.h"
#include "MrProtSrv/Domain/MrProtocol/libUILink/IParameterLookup.h"

namespace UICtrl
{

    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  setProtocolParameter                                    *
    // *                                                                        *
    // * Description : The function takes the nametag and a value as argument   *
    // *               and sets this value into the protocol using the          *
    // *               set handler (with SET_UNFORCED, aka. without calling     *
    // *               solve handler)                                           *
    // *                                                                        *
    // * Return      : true for success, else false                             *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    template <class UILinkParameterType, typename valueType>
    bool setProtocolParameter(
        MrUILinkBase*       pThis,
        const char*         cNametag,
        const valueType     desiredValue,
        long                lIndex = 0,
        bool                bAddDependency = true,
        bool                bCheckIfDesiredValueIsRealized = true
        )
    {
        valueType               oldValue;
        valueType               realizedValue;
        UILinkParameterType*    pParameterInProtocol;

        // get pointer to parameter and check for is editable
		pParameterInProtocol = zelda::Editable(zelda::Type<UILinkParameterType>(pThis->GetParameterLookup().Search(cNametag)), lIndex);

        // check if the pointer is valid
        if (!pParameterInProtocol)
            return false;

        // get old value
        oldValue = pParameterInProtocol->value(lIndex);

        // set value
		realizedValue = pParameterInProtocol->value(desiredValue, lIndex, MrUILinkBase::SET_MODE::SET_UNFORCED);

        // add dependency
        if (bAddDependency && (realizedValue != oldValue))
            pThis->addDependentParamPtr(pParameterInProtocol, lIndex);

        // check value (maybe refine for double tolerance)
        if (bCheckIfDesiredValueIsRealized && (realizedValue != desiredValue))
            return false;

        // success
        return true;
    }

    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Name        :  setProtocolParameterElm                                 *
    // *                                                                        *
    // * Description : The function takes the container, the element and a      *
    // *               value as argument (sometimes container is enough)        *
    // *               and sets this value into the protocol using the          *
    // *               set handler (with SET_UNFORCED, aka. without calling     *
    // *               solve handler)                                           *
    // *                                                                        *
    // * Return      : true for success, else false                             *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    template <class UILinkParameterType, typename valueType>
    bool setProtocolParameterElm(
        MrUILinkBase*       pThis,
        const char*         cContainer,
        const char*         cElement,
        const valueType     desiredValue,
        long                lIndex = 0,
        bool                bAddDependency = true,
        bool                bCheckIfDesiredValueIsRealized = true
        )
    {
        valueType               oldValue;
        valueType               realizedValue;
        UILinkParameterType*    pParameterInProtocol;

        // get pointer to parameter and check for is editable
        pParameterInProtocol = _searchElm< UILinkParameterType >(pThis, cContainer, cElement);

        // check if the pointer is valid
        if(!pParameterInProtocol)
            return false;

        // get old value
        oldValue = pParameterInProtocol->value(lIndex);

        // set value
        realizedValue = pParameterInProtocol->value(desiredValue, lIndex, MrUILinkBase::SET_MODE::SET_UNFORCED);

        // add dependency
        if(bAddDependency && (realizedValue != oldValue))
            pThis->addDependentParamPtr(pParameterInProtocol, lIndex);

        // check value (maybe refine for double tolerance)
        if(bCheckIfDesiredValueIsRealized && (realizedValue != desiredValue))
            return false;

        // success
        return true;
    }

}
