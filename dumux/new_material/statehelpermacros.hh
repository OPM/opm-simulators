/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: andreas.lauser _at_ iws.uni-stuttgart.de                         *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Convenience macros to ease writting state classes.
 */
#ifndef STATE_HELPER_MACROS_HH
#define STATE_HELPER_MACROS_HH

//////////////////////////////////////////////////////////////////
// Internal macros
//////////////////////////////////////////////////////////////////

#define __PROPERTY_SETTER(Type, Name, SetterName)   \
    public:                                         \
    void SetterName(Type val)                       \
    { Name##_ = val; }

#define __PROPERTY_GETTER(Type, Name, GetterName)   \
    public:                                         \
    Type GetterName() const                         \
    { return Name##_; }

#define __PROPERTY_MEMBER(Type, Name)           \
    protected:                                  \
    Type Name##_;                               \
public:

#define __PROPERTY_PROXY_SETTER(FwdObj, Type, Name, SetterName) \
    public:                                                     \
    void SetterName(Type val)                                   \
    { (FwdObj)->SetterName(val); }

#define __PROPERTY_PROXY_GETTER(FwdObj, Type, Name, GetterName) \
    public:                                                     \
    Type GetterName() const                                     \
    { return (FwdObj)->Name(); }


//////////////////////////////////////////////////////////////////
// Macros for properties which are not mutable (i.e. that
// the getter returns a const reference)
//////////////////////////////////////////////////////////////////

#define PROPERTY(Type, Name, SetterName)                \
    __PROPERTY_SETTER(const Type &, Name, SetterName)   \
    __PROPERTY_GETTER(const Type &, Name, Name)         \
    __PROPERTY_MEMBER(Type, Name)

#define PARAMETER(Type, Name)                   \
    __PROPERTY_GETTER(const Type &, Name, Name) \
    __PROPERTY_MEMBER(Type, Name)

#define CONST_PARAMETER(Type, Name, Value)      \
    public:                                     \
    Type Name() const                           \
    { return Type(Value); }

#define PTR_PROPERTY(Type, Name, SetterName)    \
    __PROPERTY_SETTER(Type *, Name, SetterName) \
    __PROPERTY_GETTER(const Type *, Name, Name) \
    __PROPERTY_MEMBER(Type*, Name)

#define PTR_PARAMETER(Type, Name)               \
    __PROPERTY_GETTER(const Type *, Name, Name) \
    __PROPERTY_MEMBER(Type*, Name)


//////////////////////////////////////////////////////////////////
// Macros for properties which are mutable (i.e. that
// the getter returns a non-const reference)
//////////////////////////////////////////////////////////////////

#define MUTABLE_PROPERTY(Type, Name, SetterName)        \
    __PROPERTY_SETTER(const Type &, Name, SetterName)   \
    __PROPERTY_GETTER(Type &, Name, Name)               \
    __PROPERTY_MEMBER(mutable Type, Name)

#define MUTABLE_PARAMETER(Type, Name)           \
    __PROPERTY_GETTER(Type &, Name, Name)       \
    __PROPERTY_MEMBER(mutable Type, Name)

#define MUTABLE_PTR_PROPERTY(Type, Name, SetterName)    \
    __PROPERTY_SETTER(Type *, Name, SetterName)         \
    __PROPERTY_GETTER(Type *, Name, Name)               \
    __PROPERTY_MEMBER(mutable Type*, Name)

#define MUTABLE_PTR_PARAMETER(Type, Name)       \
    __PROPERTY_GETTER(Type *, Name, Name)       \
    __PROPERTY_MEMBER(mutable Type*, Name)


//////////////////////////////////////////////////////////////////
// Macros for properties non mutable proxy properties
// (proxy properties don't store the data themselves but
// forward all calls to a "forward object")
//////////////////////////////////////////////////////////////////

#define PROXY_PROPERTY(FwdObj, Type, Name, SetterName)              \
    __PROPERTY_PROXY_SETTER(FwdObj, const Type &, Name, SetterName) \
    __PROPERTY_PROXY_GETTER(FwdObj, const Type &, Name, Name)

#define PROXY_PARAMETER(FwdObj, Type, Name)                     \
    __PROPERTY_PROXY_GETTER(FwdObj, const Type &, Name, Name)

#define PROXY_PTR_PROPERTY(FwdObj, Type, Name, SetterName)      \
    __PROPERTY_PROXY_SETTER(FwdObj, Type *, Name, SetterName)   \
    __PROPERTY_PROXY_GETTER(FwdObj, const Type *, Name, Name)

#define PROXY_PTR_PARAMETER(FwdObj, Type, Name)                 \
    __PROPERTY_PROXY_GETTER(FwdObj, const Type *, Name, Name)

//////////////////////////////////////////////////////////////////
// Macros for properties mutable proxy properties
//////////////////////////////////////////////////////////////////

#define MUTABLE_PROXY_PROPERTY(FwdObj, Type, Name, SetterName)      \
    __PROPERTY_PROXY_SETTER(FwdObj, const Type &, Name, SetterName) \
    __PROPERTY_PROXY_GETTER(FwdObj, Type &, Name, Name)

#define MUTABLE_PROXY_PARAMETER(FwdObj, Type, Name)             \
    __PROPERTY_PROXY_GETTER(FwdObj, const Type &, Name, Name)

#define MUTABLE_PROXY_PTR_PROPERTY(FwdObj, Type, Name, SetterName)  \
    __PROPERTY_PROXY_SETTER(FwdObj, Type *, Name, SetterName)       \
    __PROPERTY_PROXY_GETTER(FwdObj, Type *, Name, Name)

#define MUTABLE_PROXY_PTR_PARAMETER(FwdObj, Type, Name) \
    __PROPERTY_PROXY_GETTER(FwdObj, Type *, Name, Name)

#endif // STATE_HELPER_MACROS_HH
