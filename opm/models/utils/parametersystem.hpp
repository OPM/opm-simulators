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
 *
 * \brief This file provides the infrastructure to retrieve run-time parameters
 *
 * Internally, runtime parameters are implemented using
 * Dune::ParameterTree with the default value taken from the parameter
 * definition.
 */
#ifndef OPM_PARAMETER_SYSTEM_HPP
#define OPM_PARAMETER_SYSTEM_HPP

#if HAVE_QUAD
#include <opm/material/common/quad.hpp>
#endif // HAVE_QUAD

#include <dune/common/classname.hh>
#include <dune/common/parametertree.hh>

#include <charconv>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>

#include <unistd.h>
#include <sys/ioctl.h>

namespace Opm::Parameters {

namespace detail {

template <typename, class = void>
struct has_name : public std::false_type {};

template <typename T>
struct has_name<T, std::void_t<decltype(std::declval<T>().name)>>
: public std::true_type {};

//! get the name data member of a parameter
template<class Parameter>
auto getParamName()
{
    if constexpr (has_name<Parameter>::value) {
        return Parameter::name;
    } else {
        std::string paramName = Dune::className<Parameter>();
        paramName.replace(0, std::strlen("Opm::Parameters::"), "");
        const auto pos = paramName.find_first_of('<');
        if (pos != std::string::npos) {
            paramName.erase(pos);
        }
        return paramName;
    }
}

//! \brief Private implementation.
void Hide_(const std::string& paramName);

}

struct ParamInfo
{
    std::string paramName;
    std::string paramTypeName;
    std::string typeTagName;
    std::string usageString;
    std::string defaultValue;
    bool isHidden;

    bool operator==(const ParamInfo& other) const;
};

/*!
 * \ingroup Parameter
 *
 * \brief Retrieve a runtime parameter.
 *
 * The default value is specified in the parameter struct.
 *
 * Example:
 *
 * \code
 * // Retrieves value UpwindWeight, default
 * // is taken from the property UpwindWeight
 * ::Opm::Parameters::get<::Opm::Parameters::UpwindWeight>();
 * \endcode
 */
template <class Param>
auto Get(bool errorIfNotRegistered = true);

/*!
 * \ingroup Parameter
 *
 * \brief Set a runtime parameter.
 *
 * Override the default value specified.
 *
 * Example:
 *
 * \code
 * // Set the value UpwindWeight
 * ::Opm::Parameters::Set<::Opm::Parameters::UpwindWeight>(3.0);
 * \endcode
 */
template <class Param>
auto SetDefault(decltype(Param::value) new_value);

class ParamRegFinalizerBase_
{
public:
    virtual ~ParamRegFinalizerBase_()
    {}
    virtual void retrieve() = 0;
};

template <class Param>
class ParamRegFinalizer_ : public ParamRegFinalizerBase_
{
public:
    void retrieve() override
    {
        // retrieve the parameter once to make sure that its value does
        // not contain a syntax error.
        std::ignore = Get<Param>(/*errorIfNotRegistered=*/true);
    }
};

struct MetaData
{
    using type = Dune::ParameterTree;

    static Dune::ParameterTree& tree()
    { return *storage_().tree; }

    static std::map<std::string, ParamInfo>& mutableRegistry()
    { return storage_().registry; }

    static const std::map<std::string, ParamInfo>& registry()
    { return storage_().registry; }

    static std::list<std::unique_ptr<ParamRegFinalizerBase_>> &registrationFinalizers()
    { return storage_().finalizers; }

    static bool& registrationOpen()
    { return storage_().registrationOpen; }

    static void clear()
    {
        storage_().tree = std::make_unique<Dune::ParameterTree>();
        storage_().finalizers.clear();
        storage_().registrationOpen = true;
        storage_().registry.clear();
    }

private:
    // this is not pretty, but handling these attributes as static variables inside
    // member functions of the ParameterMetaData property class triggers a bug in clang
    // 3.5's address sanitizer which causes these variables to be initialized multiple
    // times...
    struct Storage_
    {
        Storage_()
        {
            tree = std::make_unique<Dune::ParameterTree>();
            registrationOpen = true;
        }

        std::unique_ptr<Dune::ParameterTree> tree;
        std::map<std::string, ParamInfo> registry;
        std::list<std::unique_ptr<ParamRegFinalizerBase_>> finalizers;
        bool registrationOpen;
    };

    static Storage_& storage_()
    {
        static Storage_ obj;
        return obj;
    }
};

// function prototype declarations
void printParamUsage(std::ostream& os, const ParamInfo& paramInfo);

std::string breakLines(const std::string& msg,
                       int indentWidth,
                       int maxWidth);

int getTtyWidth();

void getFlattenedKeyList(std::list<std::string>& dest,
                         const Dune::ParameterTree& tree,
                         const std::string& prefix = "");

// print the values of a list of parameters
void printParamList(std::ostream& os,
                    const std::list<std::string>& keyList,
                    bool printDefaults = false);

//! \endcond

/*!
 * \ingroup Parameter
 * \brief Print a usage message for all run-time parameters.
 *
 * \param helpPreamble The string that is printed after the error message and before the
 *                     list of parameters.
 * \param errorMsg The error message to be printed, if any
 * \param os The \c std::ostream which should be used.
 */
void printUsage(const std::string& helpPreamble,
                const std::string& errorMsg = "",
                std::ostream& os = std::cerr,
                const bool showAll = false);

/// \cond 0
int noPositionalParameters_(std::function<void(const std::string&, const std::string&)>,
                            std::set<std::string>&,
                            std::string& errorMsg,
                            int,
                            const char** argv,
                            int paramIdx,
                            int);

/// \endcond

//! \brief Callback function for command line parsing.
using PositionalArgumentCallback = std::function<int(std::function<void(const std::string&,
                                                                        const std::string&)>,
                                                     std::set<std::string>&,
                                                     std::string&,
                                                     int,
                                                     const char**,
                                                     int,
                                                     int)>;
/*!
 * \ingroup Parameter
 * \brief Parse the parameters provided on the command line.
 *
 * This function does some basic syntax checks.
 *
 * \param argc The number of parameters passed by the operating system to the
 *             main() function
 * \param argv The array of strings passed by the operating system to the main()
 *             function
 * \param helpPreamble If non-empty, the --help and -h parameters will be recognized and
 *                     the content of the string will be printed before the list of
 *                     command line parameters
 * \return Empty string if everything worked out. Otherwise the thing that could
 *         not be read.
 */
std::string
parseCommandLineOptions(int argc,
                        const char **argv,
                        const std::string& helpPreamble = "",
                        const PositionalArgumentCallback& posArgCallback = noPositionalParameters_);

/*!
 * \ingroup Parameter
 * \brief Read the parameters from an INI-style file.
 *
 * This function does some basic syntax checks.
 */
void parseParameterFile(const std::string& fileName, bool overwrite = true);

/*!
 * \ingroup Parameter
 * \brief Print values of the run-time parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 */
void printValues(std::ostream& os = std::cout);

/*!
 * \ingroup Parameter
 * \brief Print the list of unused run-time parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 *
 * \return true if something was printed
 */
bool printUnused(std::ostream& os = std::cout);

template <class Param>
auto Get(bool errorIfNotRegistered)
{
    const std::string paramName = detail::getParamName<Param>();
    if (errorIfNotRegistered) {
        if (MetaData::registrationOpen())
            throw std::runtime_error("Parameters can only retrieved after _all_ of them have "
                                     "been registered.");

        if (MetaData::registry().find(paramName) == MetaData::registry().end()) {
            throw std::runtime_error("Accessing parameter " + paramName
                                     +" without prior registration is not allowed.");
        }
    }

    using ParamType = std::conditional_t<std::is_same_v<decltype(Param::value),
                                                        const char* const>, std::string,
                                         std::remove_const_t<decltype(Param::value)>>;
    ParamType defaultValue = Param::value;

    const std::string& defVal = MetaData::mutableRegistry()[paramName].defaultValue;
    if constexpr (std::is_same_v<ParamType, std::string>) {
        defaultValue = defVal;
    }
    else if constexpr (std::is_same_v<ParamType, bool>) {
        defaultValue = defVal == "1";
    }
#if HAVE_QUAD
    else if constexpr (std::is_same_v<ParamType, quad>) {
        defaultValue = std::strtold(defVal.data(), nullptr);
    }
#endif
#if !HAVE_FLOATING_POINT_FROM_CHARS
    else if constexpr (std::is_floating_point_v<ParamType>) {
        defaultValue = std::strtod(defVal.c_str(), nullptr);
    }
#endif // !HAVE_FLOATING_POINT_FROM_CHARS
    else {
        std::from_chars(defVal.data(), defVal.data() + defVal.size(), defaultValue);
    }

    // prefix the parameter name by the model's GroupName. E.g. If
    // the model specifies its group name to be 'Stokes', in an
    // INI file this would result in something like:
    //
    // [Stokes]
    // NewtonWriteConvergence = true
    // retrieve actual parameter from the parameter tree
    return MetaData::tree().template get<ParamType>(paramName, defaultValue);
}

template <class Param>
auto SetDefault(decltype(Param::value) new_value)
{
    const std::string paramName = detail::getParamName<Param>();
    if (MetaData::registry().find(paramName) == MetaData::registry().end()) {
        throw std::runtime_error("Accessing parameter " + paramName +
                                 " without prior registration is not allowed.");
    }
    std::ostringstream oss;
    oss << new_value;
    MetaData::mutableRegistry()[paramName].defaultValue = oss.str();
}

/*!
 * \brief Retrieves the lists of parameters specified at runtime and their values.
 *
 * The two arguments besides the TypeTag are assumed to be STL containers which store
 * std::pair<std::string, std::string>.
 */
template <class Container>
void getLists(Container& usedParams, Container& unusedParams)
{
    usedParams.clear();
    unusedParams.clear();

    if (MetaData::registrationOpen()) {
        throw std::runtime_error("Parameter lists can only retrieved after _all_ of them have "
                                 "been registered.");
    }

    // get all parameter keys
    std::list<std::string> allKeysList;
    getFlattenedKeyList(allKeysList, MetaData::tree());

    for (const auto& key : allKeysList) {
        if (MetaData::registry().find(key) == MetaData::registry().end()) {
            // key was not registered
            unusedParams.emplace_back(key, MetaData::tree()[key]);
        }
        else {
            // key was registered
            usedParams.emplace_back(key, MetaData::tree()[key]);
        }
    }
}

/*!
 *  \brief Reset parameter system.
 */
void reset();

/*!
 * \brief Returns true if a parameter has been specified at runtime, false
 *        otherwise.
 *
 * If the parameter in question has not been registered, this throws an exception.
 */
template <class Param>
bool IsSet(bool errorIfNotRegistered = true)
{
    const std::string paramName = detail::getParamName<Param>();

    if (errorIfNotRegistered) {
        if (MetaData::registrationOpen()) {
            throw std::runtime_error("Parameters can only checked after _all_ of them have "
                                     "been registered.");
        }

        if (MetaData::registry().find(paramName) == MetaData::registry().end())
            throw std::runtime_error("Accessing parameter " + std::string(paramName) +
                                     " without prior registration is not allowed.");
    }

    // check whether the parameter is in the parameter tree
    return MetaData::tree().hasKey(paramName);
}

/*!
 * \ingroup Parameter
 *
 * \brief Register a run-time parameter.
 *
 * In OPM, parameters can only be used after they have been
 * registered.
 *
 * Example:
 *
 * \code
 * // Registers a run-time parameter "UpwindWeight"
 *    and the description "Relative weight of the upwind node."
 * Register<UpwindWeight>("Relative weight of the upwind node.");
 * \endcode
 */
template <class Param>
void Register(const char* usageString)
{
    const std::string paramName = detail::getParamName<Param>();
    if (!MetaData::registrationOpen()) {
        throw std::logic_error("Parameter registration was already closed before "
                               "the parameter '" + paramName + "' was registered.");
    }

    const auto defaultValue = Param::value;
    using ParamType = std::conditional_t<std::is_same_v<decltype(defaultValue),
                                                        const char* const>, std::string,
                                         std::remove_const_t<decltype(defaultValue)>>;
    MetaData::registrationFinalizers().push_back(
        std::make_unique<ParamRegFinalizer_<Param>>());

    ParamInfo paramInfo;
    paramInfo.paramName = paramName;
    paramInfo.paramTypeName = Dune::className<ParamType>();
    paramInfo.usageString = usageString;
    std::ostringstream oss;
    oss << defaultValue;
    paramInfo.defaultValue = oss.str();
    paramInfo.isHidden = false;
    if (MetaData::registry().find(paramName) != MetaData::registry().end()) {
        // allow to register a parameter twice, but only if the
        // parameter name, type and usage string are exactly the same.
        if (MetaData::registry().at(paramName) == paramInfo) {
            return;
        }
        throw std::logic_error("Parameter " + paramName
                               +" registered twice with non-matching characteristics.");
    }

    MetaData::mutableRegistry()[paramName] = paramInfo;
}

/*!
 * \brief Indicate that a given parameter should not be mentioned in the help message
 *
 * This allows to deal with unused parameters
 */
template <class Param>
void Hide()
{
    detail::Hide_(detail::getParamName<Param>());
}

/*!
 * \brief Indicate that all parameters are registered for a given type tag.
 *
 * If registerParam is called after the invocation of
 * \c endParamRegistration, a <tt>std::logic_error</tt> exception
 * will be thrown.
 */
void endRegistration();
//! \endcond

} // namespace Opm::Parameters

#endif // OPM_PARAMETER_SYSTEM_HPP
