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

#include <dune/common/classname.hh>

#include <cstring>
#include <functional>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

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
template<class ParamType>
ParamType Get_(const std::string& paramName, ParamType defaultValue,
               bool errorIfNotRegistered);

//! \brief Private implementation.
void Hide_(const std::string& paramName);

//! \brief Private implementation.
bool IsSet_(const std::string& paramName, bool errorIfNotRegistered);

//! \brief Private implementation.
void Register_(const std::string& paramName,
               const std::string& paramTypeName,
               const std::string& defaultValue,
               const char* usageString);

//! \brief Private implementation.
void SetDefault_(const std::string& paramName,
                 const std::string& paramValue);

}

//! \endcond

/*!
 * \ingroup Parameter
 * \brief Print a usage message for all run-time parameters.
 *
 * \param helpPreamble The string that is printed after the error message and before the
 *                     list of parameters.
 * \param os The \c std::ostream which should be used.
 * \param errorMsg The error message to be printed, if any
 * \param showAll True to show all parameters
 */
void printUsage(const std::string& helpPreamble,
                std::ostream& os,
                const std::string& errorMsg = "",
                const bool showAll = false);

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
 * \param posArgCallback Callback for positional arguments
 * \param helpPreamble If non-empty, the --help and -h parameters will be recognized and
 *                     the content of the string will be printed before the list of
 *                     command line parameters
 * \return Empty string if everything worked out. Otherwise the thing that could
 *         not be read.
 */
std::string
parseCommandLineOptions(int argc,
                        const char **argv,
                        const PositionalArgumentCallback& posArgCallback,
                        const std::string& helpPreamble = "");

/*!
 * \ingroup Parameter
 * \brief Read the parameters from an INI-style file.
 *
 * This function does some basic syntax checks.
 */
bool parseParameterFile(const std::string& fileName, bool overwrite = true);

/*!
 * \ingroup Parameter
 * \brief Print values of the run-time parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 */
void printValues(std::ostream& os);

/*!
 * \ingroup Parameter
 * \brief Print the list of unused run-time parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 *
 * \return true if something was printed
 */
bool printUnused(std::ostream& os);

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
auto Get(bool errorIfNotRegistered = true)
{
    using ParamType = std::conditional_t<std::is_same_v<decltype(Param::value),
                                                        const char* const>, std::string,
                                         std::remove_const_t<decltype(Param::value)>>;
    ParamType defaultValue = Param::value;
    return detail::Get_(detail::getParamName<Param>(),
                        defaultValue, errorIfNotRegistered);
}

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
void SetDefault(decltype(Param::value) new_value)
{
    const std::string paramName = detail::getParamName<Param>();
    std::ostringstream oss;
    oss << new_value;
    detail::SetDefault_(paramName, oss.str());
}

/*!
 * \brief A struct holding the key-value pair for a parameter.
 */
struct Parameter
{
    Parameter(const std::string& k, const std::string& v)
        : key(k), value(v)
    {}

    friend std::ostream& operator<<(std::ostream& os, const Parameter& param)
    {
        os << param.key << "=\"" << param.value << '"';
        return os;
    }

    bool operator==(const Parameter& setting) const
    {
        return setting.key == key
            && setting.value == value;
    }

    bool operator !=(const Parameter& setting) const
    {
        return !(*this == setting);
    }

    std::string key, value;
};

/*!
 * \brief Retrieves the lists of parameters specified at runtime and their values.
 *
 * The two arguments besides the TypeTag are assumed to be STL containers which store
 * std::pair<std::string, std::string>.
 */
void getLists(std::vector<Parameter>& usedParams,
              std::vector<Parameter>& unusedParams);

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
    return detail::IsSet_(detail::getParamName<Param>(), errorIfNotRegistered);
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
    const auto defaultValue = Param::value;
    using ParamType = std::conditional_t<std::is_same_v<decltype(defaultValue),
                                                        const char* const>, std::string,
                                         std::remove_const_t<decltype(defaultValue)>>;

    std::ostringstream oss;
    oss << defaultValue;
    detail::Register_(paramName, Dune::className<ParamType>(), oss.str(), usageString);
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
 * \brief Query whether parameter registration is open or not.
 * \return True if registration is open, false otherwise
 */
bool IsRegistrationOpen();

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
