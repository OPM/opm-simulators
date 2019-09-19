// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief Provides the magic behind the eWoms property system.
 *
 * Properties allow to associate arbitrary data types to
 * identifiers. A property is always defined on a pair (\c TypeTag,
 * \c PropertyTag) where \c TypeTag is the identifier for the object the
 * property is defined for and \c PropertyTag is an unique identifier of
 * the property.
 *
 * Type tags are hierarchic and inherit properties defined on their
 * ancesters. At each level, properties defined on lower levels can be
 * overwritten or even made undefined.
 *
 * Properties may use other properties for the respective type tag and
 * these properties can also be defined on an arbitrary level of the
 * hierarchy. The only restriction on this is that cycles are not
 * allowed when defining properties.
 */
#ifndef EWOMS_PROPERTIES_HH
#define EWOMS_PROPERTIES_HH

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/common/classname.hh>

#include <type_traits> // required for 'is_base_of<A, B>'
#include <map>
#include <set>
#include <list>
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>
#include <tuple>
#include <memory>

//! \cond SKIP_THIS

namespace Opm {
namespace Properties {

#define EWOMS_GET_HEAD_(Arg1, ...) Arg1

#if !defined NO_PROPERTY_INTROSPECTION

//! Internal macro which is only required if the property introspection is enabled
#define PROP_INFO_(EffTypeTagName, PropKind, PropTagName, ...)          \
    template <>                                                         \
    struct PropertyInfo<TTAG(EffTypeTagName), PTAG(PropTagName)>        \
    {                                                                   \
        static int init() {                                             \
            PropertyRegistryKey key(                                    \
                /*effTypeTagName=*/ Dune::className<TTAG(EffTypeTagName)>(), \
                /*kind=*/PropKind,                                      \
                /*name=*/propertyName(),                                \
                /*value=*/#__VA_ARGS__,                                 \
                /*file=*/__FILE__,                                      \
                /*line=*/__LINE__);                                     \
            PropertyRegistry::addKey(key);                              \
            return 0;                                                   \
        }                                                               \
        static std::string propertyName() { return #PropTagName; }      \
    };                                                                  \
    namespace fooPropInfo_ ## EffTypeTagName {                          \
    static const int foo_ ## PropTagName OPM_UNUSED  =                  \
        PropertyInfo<TTAG(EffTypeTagName), PTAG(PropTagName)>::init();  \
    }


//! Internal macro which is only required if the property introspection is enabled
#define TTAG_INFO_(TagName, ...)                                        \
    template <>                                                         \
    struct TypeTagInfo<TTAG(TagName)>                                   \
    {                                                                   \
        static int init() {                                             \
            TypeTagRegistry::addAllChildren<TTAG(TagName), __VA_ARGS__>(); \
            return 0;                                                   \
        }                                                               \
    };                                                                  \
    static const int fooTypeTagInfo_ ## TagName OPM_UNUSED =            \
        TypeTagInfo<TagName>::init();

//! Internal macro which is only required if the property introspection is enabled
#define SPLICE_INFO_(SpliceName, ...)                                   \
    template <>                                                         \
    struct SpliceInfo<TTAG(SpliceName)>                                 \
    {                                                                   \
        static int init() {                                             \
            TypeTagRegistry::addAllSplices<TTAG(SpliceName), __VA_ARGS__>(); \
            return 0;                                                   \
        }                                                               \
    };                                                                  \
    static const int fooSpliceInfo_ ## SpliceName OPM_UNUSED =          \
        SpliceInfo<TTAG(SpliceName)>::init();

#else
//! Don't do anything if introspection is disabled
#define PROP_INFO_(EffTypeTagName, PropKind, PropTagName, ...)
#define TTAG_INFO_(EffTypeTagName, ...)
#define SPLICE_INFO_(EffTypeTagName, ...)
#endif

// some macros for simplification

//! \endcond

/*!
 * \ingroup Properties
 * \brief Indicates that property definitions follow
 */
#define BEGIN_PROPERTIES namespace Opm { namespace Properties {

/*!
 * \ingroup Properties
 * \brief Indicates that all properties have been specified (for now)
 */
#define END_PROPERTIES }}

/*!
 * \ingroup Properties
 * \brief Convert a type tag name to a type
 *
 * The main advantage of the type of a \c TypeTag is that it can be
 * passed as a template argument.
 */
#define TTAG(TypeTagName) Opm::Properties::TTag::TypeTagName

/*!
 * \ingroup Properties
 * \brief Makes a type out of a property tag name
 *
 * Again property type names can be passed as template argument. This
 * is rarely needed, though.
 */
#define PTAG(PropTagName) Opm::Properties::PTag::PropTagName

/*!
 * \ingroup Properties
 * \brief Define a new type tag.
 *
 * A type tag can inherit the properties defined on up to five parent
 * type tags. Examples:
 *
 * \code
 * // The type tag doesn't inherit any properties from other type tags
 * NEW_TYPE_TAG(FooTypeTag);
 *
 * // BarTypeTag inherits all properties from FooTypeTag
 * NEW_TYPE_TAG(BarTypeTag, INHERITS_FROM(FooTypeTag));
 *
 * // FooBarTypeTag inherits the properties of FooTypeTag as well as
 * // those of BarTypeTag. Properties defined on BarTypeTag have
 * // preceedence over those defined for FooTypeTag:
 * NEW_TYPE_TAG(FooBarTypeTag, INHERITS_FROM(FooTypeTag, BarTypeTag));
 * \endcode
 */
#define NEW_TYPE_TAG(...)                           \
    namespace TTag {                                \
    struct EWOMS_GET_HEAD_(__VA_ARGS__, dontcare)   \
        : public TypeTag<__VA_ARGS__>               \
    { };                                            \
    TTAG_INFO_(__VA_ARGS__, void)                   \
    }                                               \
    extern int semicolonHack_

/*!
 * \ingroup Properties
 * \brief Define splices for a given type tag.
 *
 * Splices can be seen as children which can be overridden lower in
 * the hierarchy.  It can thus be seen as a "deferred inheritance"
 * mechanism. Example:
 *
 * \code
 * // First, define type tags for two different linear solvers:
 * // BiCGStab and SuperLU. The first needs the "MaxIterations"
 * // property, the second defines the "UsePivoting" property.
 * NEW_TYPE_TAG(BiCGStabSolver);
 * NEW_PROP_TAG(MaxIterations);
 * SET_INT_PROP(BiCGStabSolver, MaxIterations, 100);
 *
 * NEW_TYPE_TAG(SuperLUSolver);
 * NEW_PROP_TAG(UsePivoting);
 * SET_BOOL_PROP(SuperLUSolver, UsePivoting, true);
 *
 * // The model type tag defines the splice 'LinearSolver' and sets it
 * // to the 'BiCGStabSolver' type tag.
 * NEW_TYPE_TAG(ModelTypeTag);
 * NEW_PROP_TAG(LinearSolver);
 * SET_SPLICES(ModelTypeTag, LinearSolver);
 * SET_TAG_PROP(ModelTypeTag, LinearSolver, BiCGStabSolver);
 *
 * // The problem type tag is derived from the model type tag, but uses
 * // the SuperLU solver. Since this is done using a splice, all properties
 * // defined for the "SuperLUSolver" are inherited and the ones for the
 * // BiCGStabSolver type tag are undefined
 * NEW_TYPE_TAG(ProblemTypeTag, INHERITS_FROM(ModelTypeTag));
 * SET_TAG_PROP(ProblemTypeTag, LinearSolver, SuperLUSolver);
 * \endcode
 */
#define SET_SPLICES(TypeTagName, ...)                               \
    namespace PTag {                                                \
    template<>                                                      \
    struct Splices<TTAG(TypeTagName)>                               \
    {                                                               \
        typedef RevertedTuple<__VA_ARGS__>::type tuple;             \
    };                                                              \
    SPLICE_INFO_(TypeTagName, __VA_ARGS__)                    \
    }                                                               \
    extern int semicolonHack_

/*!
 * \ingroup Properties
 * \brief Syntactic sugar for NEW_TYPE_TAG.
 *
 * See the documentation for NEW_TYPE_TAG.
 */
#define INHERITS_FROM(...) __VA_ARGS__

/*!
// * \ingroup Properties
 * \brief Define a property tag.
 *
 * A property tag is the unique identifier for a property. It may only
 * be declared once in your program. There is also no hierarchy of
 * property tags as for type tags.
 *
 * Examples:
 *
 * \code
 * NEW_PROP_TAG(blubbPropTag);
 * NEW_PROP_TAG(blabbPropTag);
 * \endcode
 */
#define NEW_PROP_TAG(PTagName)                      \
    namespace PTag {                                \
    struct PTagName; } extern int semicolonHack_

//! \cond SKIP_THIS
#define SET_PROP_(EffTypeTagName, PropKind, PropTagName, ...)   \
    template <class TypeTag>                                    \
    struct Property<TypeTag,                                    \
                    TTAG(EffTypeTagName),                       \
                    PTAG(PropTagName)>;                         \
    PROP_INFO_(EffTypeTagName,                                  \
               /*kind=*/PropKind,                               \
               PropTagName,                                     \
               /*value=*/__VA_ARGS__)                           \
    template <class TypeTag>                                    \
    struct Property<TypeTag,                                    \
                    TTAG(EffTypeTagName),                       \
                    PTAG(PropTagName) >
//! \endcond

/*!
 * \ingroup Properties
 * \brief Set a property for a specific type tag.
 *
 * After this macro, you must to specify a complete body of a class
 * template, including the trailing semicolon. If you need to retrieve
 * another property within the class body, you can use \c TypeTag as the
 * argument for the type tag for the \c GET_PROP macro.
 *
 * Example:
 *
 * \code
 * SET_PROP(FooTypeTag, blubbPropTag)
 * {
 *    static int value = 10;
 *    static int calculate(int arg)
 *    { calculateInternal_(arg); }
 *
 * private:
 *    // retrieve the blabbProp property for the TypeTag the
 *    // property is defined on. Note that blabbProb does not need to
 *    // be defined on FooTypeTag, but can also be defined for some
 *    // derived type tag.
 *    typedef typename GET_PROP(TypeTag, blabbProp) blabb;
 *
 *    static int calculateInternal_(int arg)
 *    { return arg * blabb::value; };
 * \endcode
 * };
 */
#define SET_PROP(EffTypeTagName, PropTagName)   \
    template <class TypeTag>                    \
    struct Property<TypeTag,                    \
                    TTAG(EffTypeTagName),       \
                    PTAG(PropTagName)>;         \
    PROP_INFO_(EffTypeTagName,                  \
               /*kind=*/"opaque",               \
               PropTagName,                     \
               /*value=*/"<opaque>")            \
    template <class TypeTag>                    \
    struct Property<TypeTag,                    \
                    TTAG(EffTypeTagName),       \
                    PTAG(PropTagName) >

/*!
 * \ingroup Properties
 * \brief Explicitly unset a property for a type tag.
 *
 * This means that the property will not be inherited from the type
 * tag's parents.
 *
 * Example:
 *
 * \code
 * // make the blabbPropTag property undefined for the BarTypeTag.
 * UNSET_PROP(BarTypeTag, blabbPropTag);
 * \endcode
 */
#define UNSET_PROP(EffTypeTagName, PropTagName) \
    template <>                                 \
    struct PropertyUnset<TTAG(EffTypeTagName),  \
                         PTAG(PropTagName) >;   \
    PROP_INFO_(EffTypeTagName,                  \
               /*kind=*/"withdraw",             \
               PropTagName,                     \
               /*value=*/<none>)                \
    template <>                                 \
    struct PropertyUnset<TTAG(EffTypeTagName),  \
                         PTAG(PropTagName) >    \
            : public PropertyExplicitlyUnset    \
    {}

/*!
 * \ingroup Properties
 * \brief Set a property to a simple constant integer value.
 *
 * The constant can be accessed by the \c value attribute.
 */
#define SET_INT_PROP(EffTypeTagName, PropTagName, /*Value*/...) \
    SET_PROP_(EffTypeTagName,                                   \
              /*kind=*/"int   ",                                \
              PropTagName,                                      \
              /*value=*/__VA_ARGS__)                            \
    {                                                           \
        typedef int type;                                       \
        static const int value = __VA_ARGS__;                   \
    }

/*!
 * \ingroup Properties
 * \brief Set a property to a simple constant boolean value.
 *
 * The constant can be accessed by the \c value attribute.
 */
#define SET_BOOL_PROP(EffTypeTagName, PropTagName, /*Value*/...)    \
    SET_PROP_(EffTypeTagName,                                       \
              /*kind=*/"bool  ",                                    \
              PropTagName,                                          \
              /*value=*/__VA_ARGS__)                                \
    {                                                               \
        typedef bool type;                                          \
        static const bool value = __VA_ARGS__;                      \
    }

/*!
 * \ingroup Properties
 * \brief Set a property which defines a type.
 *
 * The type can be accessed by the \c type attribute.
 */
#define SET_TYPE_PROP(EffTypeTagName, PropTagName, /*Value*/...)  \
    SET_PROP_(EffTypeTagName,                                     \
              /*kind=*/"type  ",                                  \
              PropTagName,                                        \
              /*value=*/__VA_ARGS__)                              \
    {                                                             \
        typedef __VA_ARGS__ type;                                 \
    }

/*!
 * \ingroup Properties
 * \brief This macro provides the boiler plate code to declare a static constant member
 *        outside of a property definition.
 *
 * This macro is fairly low-level and there should rarely be a need to use it.
 */
#define PROP_STATIC_CONST_MEMBER_DEFINITION_PREFIX_(EffTypeTagName, PropTagName) \
    template <class TypeTag>                                            \
    const typename Property<TypeTag, TTAG(EffTypeTagName), PTAG(PropTagName)>::type \
    Property<TypeTag, TTAG(EffTypeTagName), PTAG(PropTagName)>

/*!
 * \ingroup Properties
 * \brief This macro provides the boiler plate code to declare a static member outside of
 *        a property definition.
 *
 * This macro is fairly low-level and there should rarely be a need to use it.
 */
#define PROP_STATIC_MEMBER_DEFINITION_PREFIX_(EffTypeTagName, PropTagName) \
    template <class TypeTag>                                            \
    typename Property<TypeTag, TTAG(EffTypeTagName), PTAG(PropTagName)>::type \
    Property<TypeTag, TTAG(EffTypeTagName), PTAG(PropTagName)>

/*!
 * \ingroup Properties
 * \brief Set a property to a simple constant scalar value.
 *
 * The constant can be accessed by the \c value attribute. In order to
 * use this macro, the property tag \c Scalar needs to be defined for
 * the type tag.
 */
#define SET_SCALAR_PROP(EffTypeTagName, PropTagName, ...)               \
    SET_PROP_(EffTypeTagName,                                           \
              /*kind=*/"scalar",                                        \
              PropTagName,                                              \
              /*value=*/__VA_ARGS__)                                    \
    {                                                                   \
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;         \
    public:                                                             \
        typedef Scalar type;                                            \
        static const Scalar value;                                      \
    };                                                                  \
    PROP_STATIC_CONST_MEMBER_DEFINITION_PREFIX_(EffTypeTagName, PropTagName)::value(__VA_ARGS__)

/*!
 * \ingroup Properties
 * \brief Set a property to a simple constant string value.
 *
 * The constant can be accessed by the \c value attribute and is of
 * type <tt>std::string</tt>.
 */
#define SET_STRING_PROP(EffTypeTagName, PropTagName, ...)               \
    SET_PROP_(EffTypeTagName,                                           \
              /*kind=*/"string",                                        \
              PropTagName,                                              \
              /*value=*/__VA_ARGS__)                                    \
    {                                                                   \
    public:                                                             \
        typedef std::string type;                                       \
        static const std::string value;                                 \
    };                                                                  \
    PROP_STATIC_CONST_MEMBER_DEFINITION_PREFIX_(EffTypeTagName, PropTagName)::value(__VA_ARGS__)

/*!
 * \ingroup Properties
 * \brief Define a property containing a type tag.
 *
 * This is convenient for splices.
 */
#define SET_TAG_PROP(EffTypeTagName, PropTagName, ValueTypeTagName) \
    SET_PROP_(EffTypeTagName,                                     \
              /*kind=*/"tag   ",                                  \
              PropTagName,                                        \
              /*value=*/TTAG(ValueTypeTagName))                   \
    {                                                             \
        typedef TTAG(ValueTypeTagName) type;                      \
    }


/*!
 * \ingroup Properties
 * \brief Retrieve a property for a type tag.
 *
 * If you use \c GET_PROP within a template and want to refer to some
 * type (including the property itself), \c GET_PROP must be preceeded by
 * the '\c typename' keyword.
 */
#define GET_PROP(TypeTag, PropTagName) \
    ::Opm::Properties::GetProperty<TypeTag, PTAG(PropTagName)>::p
//!\cond SKIP_THIS
#define GET_PROP_(TypeTag, PropTag) \
    ::Opm::Properties::GetProperty<TypeTag, PropTag>::p
//!\endcond

/*!
 * \ingroup Properties
 * \brief Access the \c value attribute of a property for a type tag.
 *
 * This is just for convenience and equivalent to
 * <tt>GET_PROP(TypeTag, PropTag)::</tt><tt>value</tt> .  If the property doesn't
 * have an attribute named \c value, this yields a compiler error.
 */
#define GET_PROP_VALUE(TypeTag, PropTagName)                            \
    ::Opm::Properties::GetProperty<TypeTag, PTAG(PropTagName)>::p::value
//!\cond SKIP_THIS
#define GET_PROP_VALUE_(TypeTag, PropTag)                               \
    ::Opm::Properties::GetProperty<TypeTag, PropTag>::p::value
//!\endcond

/*!
 * \ingroup Properties
 * \brief Access the \c type attribute of a property for a type tag.
 *
 * This is just for convenience and equivalent to
 * <tt>GET_PROP(TypeTag, PropTag)::</tt><tt>type</tt>.  If the property doesn't
 * have an attribute named \c type, this yields a compiler error. Also,
 * if you use this macro within a template, it must be preceeded by
 * the \c typename keyword.
 */
#define GET_PROP_TYPE(TypeTag, PropTagName) \
    ::Opm::Properties::GetProperty<TypeTag, PTAG(PropTagName)>::p::type
//!\cond SKIP_THIS
#define GET_PROP_TYPE_(TypeTag, PropTag) \
    ::Opm::Properties::GetProperty<TypeTag, PropTag>::p::type
//!\endcond

#if !defined NO_PROPERTY_INTROSPECTION
/*!
 * \ingroup Properties
 * \brief Return a human readable diagnostic message how exactly a
 *        property was defined.
 *
 * This is only enabled if the \c NO_PROPERTY_INTROSPECTION macro is not
 * defined.
 *
 * Example:
 *
 * \code
 * int main()
 * {
 *    std::cout << PROP_DIAGNOSTIC(FooBarTypeTag, blabbPropTag) << "\n";
 * };
 * \endcode
 */
#define PROP_DIAGNOSTIC(TypeTag, PropTagName) \
    ::Opm::Properties::getDiagnostic<TypeTag>(#PropTagName)

#else
#define PROP_DIAGNOSTIC(TypeTag, PropTagName) "Property introspection disabled by macro NO_PROPERTY_INTROSPECTION."
#endif

//////////////////////////////////////////////
// some serious template kung fu. Don't look at it too closely, it
// might damage your brain!
//////////////////////////////////////////////

//! \cond SKIP_THIS

namespace PTag {}
namespace TTag {}

#if !defined NO_PROPERTY_INTROSPECTION

namespace TTag
{
template <class EffTypeTag>
struct TypeTagInfo
{};
}

namespace PTag
{
template <class EffTypeTagName>
struct SpliceInfo
{};
}

template <class EffTypeTagName, class PropTagName>
struct PropertyInfo
{};

class PropertyRegistryKey
{
public:
    PropertyRegistryKey()
    {}

    PropertyRegistryKey(const std::string& effTypeTagName,
                        const std::string& propertyKind,
                        const std::string& propertyName,
                        const std::string& propertyValue,
                        const std::string& fileDefined,
                        int lineDefined)
        : effTypeTagName_(effTypeTagName)
        , propertyKind_(propertyKind)
        , propertyName_(propertyName)
        , propertyValue_(propertyValue)
        , fileDefined_(fileDefined)
        , lineDefined_(lineDefined)
    { }

    // copy constructor
    PropertyRegistryKey(const PropertyRegistryKey&) = default;
    PropertyRegistryKey& operator=(const PropertyRegistryKey&) = default;

    const std::string& effTypeTagName() const
    { return effTypeTagName_; }
    const std::string& propertyKind() const
    { return propertyKind_; }
    const std::string& propertyName() const
    { return propertyName_; }
    const std::string& propertyValue() const
    { return propertyValue_; }
    const std::string& fileDefined() const
    { return fileDefined_; }
    int lineDefined() const
    { return lineDefined_; }

private:
    std::string effTypeTagName_;
    std::string propertyKind_;
    std::string propertyName_;
    std::string propertyValue_;
    std::string fileDefined_;
    int lineDefined_;
};


template <class TypeTag, class PropertyTag>
struct GetProperty;

class TypeTagRegistry
{
public:
    struct SpliceRegistryEntry {
        SpliceRegistryEntry(const std::string& name)
        { name_ = name; }

        std::string propertyName() const
        { return name_; }

    private:
        std::string name_;
    };

    typedef std::list<std::unique_ptr<SpliceRegistryEntry> > SpliceList;
    typedef std::map<std::string, SpliceList> SpliceListMap;

    typedef std::list<std::string> ChildrenList;
    typedef std::map<std::string, ChildrenList> ChildrenListMap;

    // this method adds the children of a type tag to the registry if the registry does
    // not yet contain an entry for that type tag
    template <class TypeTag, typename ... Children>
    static void addAllChildren()
    {
        std::string typeTagName = Dune::className<TypeTag>();

        if (children_().find(typeTagName) == children_().end())
            addChildren<TypeTag, Children...>();
    }


    // end of recursion. the last argument is not a child, but 'void'
    // which is required for the macro magic...
    template <class TypeTag, class DummyChild>
    static void addChildren()
    {}

    // the last argument is not a child, but 'void' which is required
    // for the macro magic...
    template <class TypeTag, class Child1, class Child2, typename ... RemainingChildren>
    static void addChildren()
    {
        std::string typeTagName = Dune::className<TypeTag>();

        children_()[typeTagName].emplace_front(Dune::className<Child1>());

        addChildren<TypeTag, Child2, RemainingChildren...>();
    }

    // this method adds the splices of a type tag to the registry if the registry does
    // not yet contain an entry for that type tag
    template <class TypeTag, typename ... Splices>
    static void addAllSplices()
    {
        std::string typeTagName = Dune::className<TypeTag>();

        if (splices_().find(typeTagName) == splices_().end())
            addSplices<TypeTag, Splices...>();
    }


    // end of recursion. the last argument is not a child, but 'void'
    // which is required for the macro magic...
    template <class TypeTag>
    static void addSplices()
    { }

    // the last argument is not a child, but 'void' which is required
    // for the macro magic...
    template <class TypeTag, class Splice1, typename ... RemainingSplices>
    static void addSplices()
    {
        const std::string& typeTagName = Dune::className<TypeTag>();
        const std::string& propName =
            PropertyInfo<typename GetProperty<TypeTag, Splice1>::template GetEffectiveTypeTag_<TypeTag>::type, Splice1>::propertyName();

        splices_()[typeTagName].emplace_front(new SpliceRegistryEntry(propName));

        addSplices<TypeTag, RemainingSplices...>();
    }

    static const SpliceList& splices(const std::string& typeTagName)
    { return splices_()[typeTagName]; }

    static const ChildrenList& children(const std::string& typeTagName)
    { return children_()[typeTagName]; }

private:
    static SpliceListMap& splices_()
    {
        static SpliceListMap data;
        return data;
    }
    static ChildrenListMap& children_()
    {
        static ChildrenListMap data;
        return data;
    }
};

class PropertyRegistry
{
    typedef Opm::Properties::TypeTagRegistry TypeTagRegistry;

public:
    typedef std::map<std::string, PropertyRegistryKey> KeyList;
    typedef std::map<std::string, KeyList> KeyListMap;

    static void addKey(const PropertyRegistryKey& key)
    {
        keys_()[key.effTypeTagName()][key.propertyName()] = key;
    }

    static const std::string& getSpliceTypeTagName(const std::string& typeTagName,
                                                   const std::string& propertyName)
    {
        const auto& keyIt = keys_().find(typeTagName);
        const auto& keyEndIt = keys_().end();
        if (keyIt == keyEndIt)
            throw std::runtime_error("Unknown type tag key '"+typeTagName+"'");

        // check whether the propery is directly defined for the type tag currently
        // checked, ...
        const auto& propIt = keyIt->second.find(propertyName);
        const auto& propEndIt = keyIt->second.end();
        if (propIt != propEndIt)
            return propIt->second.propertyValue();

        // ..., if not, check all splices,  ...
        typedef typename TypeTagRegistry::SpliceList SpliceList;
        const SpliceList& splices = TypeTagRegistry::splices(typeTagName);
        SpliceList::const_iterator spliceIt = splices.begin();
        for (; spliceIt != splices.end(); ++spliceIt) {
            const auto& spliceTypeTagName =
                PropertyRegistry::getSpliceTypeTagName(typeTagName,
                                                       (*spliceIt)->propertyName());

            if (spliceTypeTagName == "")
                continue;

            const auto& tmp = getSpliceTypeTagName(spliceTypeTagName, propertyName);
            if (tmp != "")
                return tmp;
        }

        // .. if still not, check all normal children.
        typedef TypeTagRegistry::ChildrenList ChildrenList;
        const ChildrenList& children = TypeTagRegistry::children(typeTagName);
        ChildrenList::const_iterator ttagIt = children.begin();
        for (; ttagIt != children.end(); ++ttagIt) {
            const auto& tmp = getSpliceTypeTagName(*ttagIt, propertyName);
            if (tmp != "")
                return tmp;
        }

        // if the property was not defined on a given type tag, return
        // the empty string.
        static std::string tmp("");
        return tmp;
    }

    static const PropertyRegistryKey& getKey(const std::string& effTypeTagName,
                                             const std::string& propertyName)
    { return keys_()[effTypeTagName][propertyName]; }

    static const KeyList& getKeys(const std::string& effTypeTagName)
    { return keys_()[effTypeTagName]; }

private:
    static KeyListMap& keys_()
    {
        static KeyListMap data;
        return data;
    }
};

#endif // !defined NO_PROPERTY_INTROSPECTION

struct PropertyUndefined { };
class PropertyExplicitlyUnset {};

template <class RealTypeTag,
          class EffectiveTypeTag,
          class PropertyTag>
struct Property : public PropertyUndefined
{};

template <class EffectiveTypeTag,
          class PropertyTag>
struct PropertyUnset : public PropertyUndefined
{};

template <class Tree, class PropertyTag>
struct propertyExplicitlyUnset
{
    const static bool value =
        std::is_base_of<PropertyExplicitlyUnset,
                        PropertyUnset<typename Tree::SelfType,
                                      PropertyTag>
                        >::value;
};

template <class Tree, class PropertyTag>
class propertyExplicitlyUnsetOnTree
{
    static const bool explicitlyUnset = propertyExplicitlyUnset<Tree, PropertyTag>::value;

    template <class ChildTuple>
    struct unsetOnAllChildren
    { static const bool value = true; };

    template <class Child, class ... RemainingChildren>
    struct unsetOnAllChildren<std::tuple<Child, RemainingChildren...> >
    { static const bool value =
            propertyExplicitlyUnsetOnTree<Child, PropertyTag>::value
            && unsetOnAllChildren<std::tuple<RemainingChildren...> >::value; };

public:
    static const bool value =
        (explicitlyUnset || (!Tree::isLeaf && unsetOnAllChildren<typename Tree::ChildrenTuple>::value));
};

template <class PropertyTag>
struct propertyExplicitlyUnsetOnTree<void, PropertyTag>
{
    const static bool value = std::true_type::value;
};

template <class RealTypeTag, class Tree, class PropertyTag>
struct propertyDefinedOnSelf
{
    const static bool value =
        ! std::is_base_of<PropertyUndefined,
                          Property<RealTypeTag,
                                   typename Tree::SelfType,
                                   PropertyTag> >::value;
};

// template class to revert the order or a std::tuple's arguments. This is required to
// make the properties of children defined on the right overwrite the properties of the
// children on the left. See https://sydius.me/2011/07/reverse-tuple-in-c/
template<typename... Args>
class RevertedTuple
{
private:
    template<unsigned int N, typename... All>
    struct RevertedTupleOuter
    {
        template<typename Head, typename... Tail>
        struct RevertedTupleInner: RevertedTupleOuter<N-1, Head, All...>::template RevertedTupleInner<Tail...> { };
    };

    template<typename... All>
    struct RevertedTupleOuter<0, All...>
    {
        template<typename... Tail>
        struct RevertedTupleInner {
            typedef std::tuple<All...> type;
        };
    };

public:
    typedef typename RevertedTupleOuter<sizeof...(Args)>::template RevertedTupleInner<Args...>::type type;
};

template <class SelfT,
          typename ... Children>
class TypeTag
{
public:
    typedef SelfT SelfType;

    typedef typename RevertedTuple<Children...>::type ChildrenTuple;
    static const bool isLeaf = std::is_same<ChildrenTuple, std::tuple<> >::value;
};

namespace PTag {
// this class needs to be located in the PTag namespace for reasons
// you don't really want to know...
template <class TypeTag>
struct Splices
{
    typedef typename std::tuple<> tuple;
};
} // namespace PTag

template <class TypeTag, class PropertyTag>
struct GetProperty
{
    // find the type tag for which the property is defined
    template <class CurTree, bool directlyDefined =
                  propertyDefinedOnSelf<TypeTag,
                                        CurTree,
                                        PropertyTag>::value>
    struct GetEffectiveTypeTag_
    { typedef typename CurTree::SelfType type; };

    template <class ...Elements>
    struct SearchTypeTagList_;

    template <class EffectiveTypeTag, class ...Elements>
    struct SearchTypeTagList_FirstThenRemaining_;

    template <class ...SpliceList>
    struct SearchSpliceList_;

    template <class EffectiveTypeTag, class ...Splices>
    struct SearchSpliceList_FirstThenRemaining_;

    // find the first type tag in a tuple for which the property is
    // defined
    template <class TypeTagTuple>
    struct SearchTypeTagTuple_
    { typedef void type; };

    template <class ...TypeTagList>
    struct SearchTypeTagTuple_<std::tuple<TypeTagList...> >
    { typedef typename SearchTypeTagList_<TypeTagList...>::type type; };

    template <class ...Elements>
    struct SearchTypeTagList_
    { typedef void type; };

    template <class FirstElement, class ...RemainingElements>
    struct SearchTypeTagList_<FirstElement, RemainingElements...>
    {
        typedef typename SearchTypeTagList_FirstThenRemaining_<
            typename GetEffectiveTypeTag_<FirstElement>::type,
            RemainingElements...>::type type;
    };

    template <class EffectiveTypeTag, class ...Elements>
    struct SearchTypeTagList_FirstThenRemaining_
    { typedef EffectiveTypeTag type; };

    template <class ...RemainingElements>
    struct SearchTypeTagList_FirstThenRemaining_<void, RemainingElements...>
    { typedef typename SearchTypeTagList_<RemainingElements...>::type type; };

    // find the first type tag in a tuple of splices for which the
    // property is defined
    template <class SpliceTuple>
    struct SearchSpliceTuple_
    { typedef void type; };

    template <class ...SpliceList>
    struct SearchSpliceTuple_<std::tuple<SpliceList...> >
    { typedef typename SearchSpliceList_<SpliceList...>::type type; };

    template <class ...SpliceList>
    struct SearchSpliceList_
    { typedef void type; };

    template <class FirstSplice, class ...RemainingSplices>
    struct SearchSpliceList_<FirstSplice, RemainingSplices...>
    {
        typedef typename SearchSpliceList_FirstThenRemaining_<
            typename GetEffectiveTypeTag_<typename GetProperty<TypeTag, FirstSplice>::p::type>::type,
            RemainingSplices...>::type type;
    };

    template <class EffectiveTypeTag, class ...Splices>
    struct SearchSpliceList_FirstThenRemaining_
    { typedef EffectiveTypeTag type; };

    template <class ...RemainingSplices>
    struct SearchSpliceList_FirstThenRemaining_<void, RemainingSplices...>
    { typedef typename SearchSpliceList_<RemainingSplices...>::type type; };

    // find the splice or the child type tag for which the property is defined
    template <class CurTree,
              class SpliceTypeTag = typename SearchSpliceTuple_< typename PTag::Splices<CurTree>::tuple >::type >
    struct SearchSplicesThenChildren_
    { typedef SpliceTypeTag type; };

    template <class CurTree>
    struct SearchSplicesThenChildren_<CurTree, void>
    { typedef typename SearchTypeTagTuple_<typename CurTree::ChildrenTuple>::type type; };

    template <class CurTree>
    struct GetEffectiveTypeTag_<CurTree, /*directlyDefined = */false>
    { typedef typename SearchSplicesThenChildren_<CurTree>::type type; };

public:
    typedef Property<TypeTag, typename GetEffectiveTypeTag_<TypeTag>::type, PropertyTag> p;
};

#if !defined NO_PROPERTY_INTROSPECTION
inline int myReplaceAll_(std::string& s,
                         const std::string& pattern,
                         const std::string& replacement);
inline int myReplaceAll_(std::string& s,
                         const std::string& pattern,
                         const std::string& replacement)
{
    size_t pos;
    int i = 0;
    while ( (pos = s.find(pattern)) != s.npos) {
        s.replace(pos, pattern.size(), replacement);
        ++i;
    };
    return i;
}

inline std::string canonicalTypeTagNameToName_(const std::string& canonicalName);
inline std::string canonicalTypeTagNameToName_(const std::string& canonicalName)
{
    std::string result(canonicalName);
    myReplaceAll_(result, "Opm::Properties::TTag::", "TTAG(");
    myReplaceAll_(result, "::", "");
    result += ")";
    return result;
}

inline bool getDiagnostic_(const std::string& typeTagName,
                           const std::string& propTagName,
                           std::string& result,
                           const std::string indent)
{
    const PropertyRegistryKey *key = 0;

    const typename PropertyRegistry::KeyList& keys =
        PropertyRegistry::getKeys(typeTagName);
    typename PropertyRegistry::KeyList::const_iterator it = keys.begin();
    for (; it != keys.end(); ++it) {
        if (it->second.propertyName() == propTagName) {
            key = &it->second;
            break;
        };
    }

    if (key) {
        std::ostringstream oss;
        oss << indent
            << key->propertyKind() << " "
            << key->propertyName() << " defined on '"
            << canonicalTypeTagNameToName_(key->effTypeTagName()) << "' at "
            << key->fileDefined() << ":" << key->lineDefined() << "\n";
        result = oss.str();
        return true;
    }

    // print properties defined on children
    typedef typename TypeTagRegistry::ChildrenList ChildrenList;
    const ChildrenList& children = TypeTagRegistry::children(typeTagName);
    ChildrenList::const_iterator ttagIt = children.begin();
    std::string newIndent = indent + "  ";
    for (; ttagIt != children.end(); ++ttagIt) {
        if (getDiagnostic_(*ttagIt, propTagName, result, newIndent)) {
            result.insert(0, indent + "Inherited from " + canonicalTypeTagNameToName_(typeTagName) + "\n");
            return true;
        }
    }

    return false;
}

template <class TypeTag>
const std::string getDiagnostic(std::string propTagName)
{
    std::string result;

    std::string TypeTagName(Dune::className<TypeTag>());

    propTagName.replace(0, strlen("PTag("), "");
    auto n = propTagName.length();
    propTagName.replace(n - 1, 1, "");
    //TypeTagName.replace(0, strlen("Opm::Properties::TTag::"), "");

    return result;
}

inline void print_(const std::string& rootTypeTagName,
                   const std::string& curTypeTagName,
                   const std::string& splicePropName,
                   std::ostream& os,
                   const std::string indent,
                   std::set<std::string>& printedProperties)
{
    if (indent == "") {
        os << indent << "###########\n";
        os << indent << "# Properties\n";
        os << indent << "###########\n";
        os << indent << "Properties for " << canonicalTypeTagNameToName_(curTypeTagName) << ":";
    }
    else if (splicePropName != "")
        os << indent << "Inherited from splice " << splicePropName << " (set to " << canonicalTypeTagNameToName_(curTypeTagName) << "):";
    else
        os << indent << "Inherited from " << canonicalTypeTagNameToName_(curTypeTagName) << ":";
    const typename PropertyRegistry::KeyList& keys =
        PropertyRegistry::getKeys(curTypeTagName);
    typename PropertyRegistry::KeyList::const_iterator it = keys.begin();
    bool somethingPrinted = false;
    for (; it != keys.end(); ++it) {
        const PropertyRegistryKey& key = it->second;
        if (printedProperties.count(key.propertyName()) > 0)
            continue; // property already printed
        if (!somethingPrinted) {
            os << "\n";
            somethingPrinted = true;
        }
        os << indent << "  "
           << key.propertyKind() << " " << key.propertyName();
        if (key.propertyKind() != "opaque") {
            std::string s(key.propertyValue());
            myReplaceAll_(s, "typename ", "");
            if (myReplaceAll_(s, "::Opm::Properties::TTag::", "TTAG("))
                s += ')';
            myReplaceAll_(s, "::Opm::Properties::PTag::", "");
            myReplaceAll_(s, "::Opm::Properties::GetProperty<", "GET_PROP(");
            myReplaceAll_(s, ">::p::", ")::");
            myReplaceAll_(s, "GET_PROP(TypeTag, Scalar)::type", "Scalar");

            os << " = '" << s << "'";
        }
        os << " defined at " << key.fileDefined()
           << ":" << key.lineDefined()
           << "\n";
        printedProperties.insert(key.propertyName());
    };
    if (!somethingPrinted)
        os << " (none)\n";
    // print properties defined on splices or children
    std::string newIndent = indent + "  ";

    // first, iterate over the splices, ...
    typedef TypeTagRegistry::SpliceList SpliceList;
    const SpliceList& splices = TypeTagRegistry::splices(curTypeTagName);
    SpliceList::const_iterator spliceIt = splices.begin();
    for (; spliceIt != splices.end(); ++ spliceIt) {
        const auto& spliceTypeTagName = PropertyRegistry::getSpliceTypeTagName(rootTypeTagName,
                                                                               (*spliceIt)->propertyName());
        print_(rootTypeTagName, spliceTypeTagName, (*spliceIt)->propertyName(), os, newIndent, printedProperties);
    }

    // ... then, over the children
    typedef typename TypeTagRegistry::ChildrenList ChildrenList;
    const ChildrenList& children = TypeTagRegistry::children(curTypeTagName);
    ChildrenList::const_iterator ttagIt = children.begin();
    for (; ttagIt != children.end(); ++ttagIt) {
        print_(rootTypeTagName, *ttagIt, /*splicePropName=*/"", os, newIndent, printedProperties);
    }
}

template <class TypeTag>
void printValues(std::ostream& os = std::cout)
{
    std::set<std::string> printedProps;
    print_(Dune::className<TypeTag>(), Dune::className<TypeTag>(), /*splicePropertyName=*/"", os, /*indent=*/"", printedProps);
}
#else // !defined NO_PROPERTY_INTROSPECTION
template <class TypeTag>
void printValues(std::ostream& os = std::cout)
{
    os <<
        "The eWoms property system was compiled with the macro\n"
        "NO_PROPERTY_INTROSPECTION defined.\n"
        "No diagnostic messages this time, sorry.\n";
}

template <class TypeTag>
const std::string getDiagnostic(std::string propTagName)
{
    std::string result;
    result =
        "The eWoms property system was compiled with the macro\n"
        "NO_PROPERTY_INTROSPECTION defined.\n"
        "No diagnostic messages this time, sorry.\n";
    return result;
}

#endif // !defined NO_PROPERTY_INTROSPECTION

//! \endcond

} // namespace Properties
} // namespace Opm

#endif
