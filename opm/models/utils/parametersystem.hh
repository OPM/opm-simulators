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
 * Dune::ParameterTree with the default value taken from the property
 * system.
 */
#ifndef EWOMS_PARAMETER_SYSTEM_HH
#define EWOMS_PARAMETER_SYSTEM_HH

#include <opm/models/utils/propertysystem.hh>

#include <opm/material/common/Unused.hpp>

#if HAVE_QUAD
#include <opm/material/common/quad.hpp>
#endif // HAVE_QUAD

#include <dune/common/classname.hh>
#include <dune/common/parametertree.hh>

#include <map>
#include <set>
#include <list>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include <unistd.h>
#include <sys/ioctl.h>

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
 * // Registers a run-time parameter "UpwindWeight" which has type
 * // "Scalar" and the description "Relative weight of the upwind node."
 * EWOMS_REGISTER_PARAM(TypeTag, Scalar, UpwindWeight,
 *                      "Relative weight of the upwind node.");
 * \endcode
 */
#define EWOMS_REGISTER_PARAM(TypeTag, ParamType, ParamName, Description)       \
    ::Opm::Parameters::registerParam<TypeTag, ParamType>( \
        #ParamName, #ParamName, getPropValue<TypeTag, Properties::ParamName>(), Description)

/*!
 * \ingroup Parameter
 *
 * \brief Indicate that a given parameter should not be mentioned in the help message
 *
 * This allows to deal with unused parameters
 */
#define EWOMS_HIDE_PARAM(TypeTag, ParamName)                \
    ::Opm::Parameters::hideParam<TypeTag>(#ParamName, getPropValue<TypeTag, Properties::ParamName>())

/*!
 * \ingroup Parameter
 *
 * \brief Indicate that all parameters are registered for a given type tag.
 *
 * If \c EWOMS_REGISTER_PARAM is called after the invokation of
 * \c END_PARAM_REGISTRATION, a <tt>std::logic_error</tt> exception
 * will be thrown.
 */
#define EWOMS_END_PARAM_REGISTRATION(TypeTag)                                  \
    ::Opm::Parameters::endParamRegistration<TypeTag>()

/*!
 * \ingroup Parameter
 *
 * \brief Retrieve a runtime parameter.
 *
 * The default value is specified via the property system.
 *
 * Example:
 *
 * \code
 * // Retrieves scalar value UpwindWeight, default
 * // is taken from the property UpwindWeight
 * EWOMS_GET_PARAM(TypeTag, Scalar, UpwindWeight);
 * \endcode
 */
#define EWOMS_GET_PARAM(TypeTag, ParamType, ParamName)                         \
    (::Opm::Parameters::get<TypeTag, ParamType>(#ParamName, #ParamName, \
                                                getPropValue<TypeTag, Properties::ParamName>()))

//!\cond SKIP_THIS
#define EWOMS_GET_PARAM_(TypeTag, ParamType, ParamName)                 \
    (::Opm::Parameters::get<TypeTag, ParamType>(#ParamName, #ParamName, \
                                                getPropValue<TypeTag, Properties::ParamName>(), \
                                                /*errorIfNotRegistered=*/false))

/*!
 * \ingroup Parameter
 *
 * \brief Retrieves the lists of parameters specified at runtime and their values.
 *
 * The two arguments besides the TypeTag are assumed to be STL containers which store
 * std::pair<std::string, std::string>.
 */
#define EWOMS_GET_PARAM_LISTS(TypeTag, UsedParamList, UnusedParamList)    \
    (::Opm::Parameters::getLists<TypeTag>(UsedParamList, UnusedParamList))

//!\cond SKIP_THIS
#define EWOMS_RESET_PARAMS_(TypeTag)            \
    (::Opm::Parameters::reset<TypeTag>())

/*!
 * \ingroup Parameter
 *
 * \brief Returns true if a parameter has been specified at runtime, false
 *        otherwise.
 *
 * If the parameter in question has not been registered, this throws an exception.
 */
#define EWOMS_PARAM_IS_SET(TypeTag, ParamType, ParamName)               \
    (::Opm::Parameters::isSet<TypeTag, ParamType>(#ParamName, #ParamName))

namespace Opm {
namespace Parameters {

struct ParamInfo
{
    std::string paramName;
    std::string paramTypeName;
    std::string typeTagName;
    std::string propertyName;
    std::string usageString;
    std::string compileTimeValue;
    bool isHidden;

    bool operator==(const ParamInfo& other) const
    {
        return other.paramName == paramName
               && other.paramTypeName == paramTypeName
               && other.typeTagName == typeTagName
               && other.propertyName == propertyName
               && other.usageString == usageString
               && other.compileTimeValue == compileTimeValue;
    }
};

// forward declaration
template <class TypeTag, class ParamType, class PropTag>
const ParamType get(const char *propTagName,
                    const char *paramName,
                    bool errorIfNotRegistered = true);
template <class TypeTag, class ParamType>
const ParamType get(const char *propTagName,
                    const char *paramName,
                    const ParamType& defaultValue,
                    bool errorIfNotRegistered = true);

class ParamRegFinalizerBase_
{
public:
    virtual ~ParamRegFinalizerBase_()
    {}
    virtual void retrieve() = 0;
};

template <class TypeTag, class ParamType>
class ParamRegFinalizer_ : public ParamRegFinalizerBase_
{
public:
    ParamRegFinalizer_(const std::string& paramName, const ParamType& defaultValue)
        : paramName_(paramName)
        , defaultValue_(defaultValue)
    {}

    virtual void retrieve() override
    {
        // retrieve the parameter once to make sure that its value does
        // not contain a syntax error.
        ParamType __attribute__((unused)) dummy =
            get<TypeTag, ParamType>(/*propTagName=*/paramName_.data(),
                                    paramName_.data(),
                                    defaultValue_,
                                    /*errorIfNotRegistered=*/true);
    }

private:
    std::string paramName_;
    ParamType defaultValue_;
};
} // namespace Parameters

} // namespace Opm

namespace Opm::Properties {

namespace TTag {

// type tag which is supposed to spliced in or inherited from if the
// parameter system is to be used
struct ParameterSystem {};

} // namespace TTag

template<class TypeTag, class MyTypeTag>
struct ParameterMetaData { using type = UndefinedProperty; };


//! Set the ParameterMetaData property
template<class TypeTag>
struct ParameterMetaData<TypeTag, TTag::ParameterSystem>
{
    using type = Dune::ParameterTree;

    static Dune::ParameterTree& tree()
    { return *storage_().tree; }

    static std::map<std::string, ::Opm::Parameters::ParamInfo>& mutableRegistry()
    { return storage_().registry; }

    static const std::map<std::string, ::Opm::Parameters::ParamInfo>& registry()
    { return storage_().registry; }

    static std::list<std::unique_ptr<::Opm::Parameters::ParamRegFinalizerBase_> > &registrationFinalizers()
    { return storage_().finalizers; }

    static bool& registrationOpen()
    { return storage_().registrationOpen; }

    static void clear()
    {
        storage_().tree.reset(new Dune::ParameterTree());
        storage_().finalizers.clear();
        storage_().registrationOpen = true;
        storage_().registry.clear();
    }

private:
    // this is not pretty, but handling these attributes as static variables inside
    // member functions of the ParameterMetaData property class triggers a bug in clang
    // 3.5's address sanitizer which causes these variables to be initialized multiple
    // times...
    struct Storage_ {
        Storage_()
        {
            tree.reset(new Dune::ParameterTree());
            registrationOpen = true;
        }

        std::unique_ptr<Dune::ParameterTree> tree;
        std::map<std::string, ::Opm::Parameters::ParamInfo> registry;
        std::list<std::unique_ptr<::Opm::Parameters::ParamRegFinalizerBase_> > finalizers;
        bool registrationOpen;
    };
    static Storage_& storage_() {
        static Storage_ obj;
        return obj;
    }
};


} // namespace Opm::Properties

namespace Opm {

namespace Parameters {
// function prototype declarations
void printParamUsage_(std::ostream& os, const ParamInfo& paramInfo);
void getFlattenedKeyList_(std::list<std::string>& dest,
                          const Dune::ParameterTree& tree,
                          const std::string& prefix = "");

inline std::string breakLines_(const std::string& msg,
                               int indentWidth,
                               int maxWidth)
{
    std::string result;
    int startInPos = 0;
    int inPos = 0;
    int lastBreakPos = 0;
    int ttyPos = 0;
    for (; inPos < int(msg.size()); ++ inPos, ++ ttyPos) {
        if (msg[inPos] == '\n') {
            result += msg.substr(startInPos, inPos - startInPos + 1);
            startInPos = inPos + 1;
            lastBreakPos = startInPos + 1;

            // we need to use -1 here because ttyPos is incremented after the loop body
            ttyPos = -1;
            continue;
        }

        if (std::isspace(msg[inPos]))
            lastBreakPos = inPos;

        if (ttyPos >= maxWidth) {
            if (lastBreakPos > startInPos) {
                result += msg.substr(startInPos, lastBreakPos - startInPos);
                startInPos = lastBreakPos + 1;
                lastBreakPos = startInPos;
                inPos = startInPos;
            }
            else {
                result += msg.substr(startInPos, inPos - startInPos);
                startInPos = inPos;
                lastBreakPos = startInPos;
                inPos = startInPos;
            }

            result += "\n";
            for (int i = 0; i < indentWidth; ++i)
                result += " ";
            ttyPos = indentWidth;
        }
    }

    result += msg.substr(startInPos);

    return result;
}

inline int getTtyWidth_()
{
    int ttyWidth = 10*1000; // effectively do not break lines at all.
    if (isatty(STDOUT_FILENO)) {
#if defined TIOCGWINSZ
        // This is a bit too linux specific, IMO. let's do it anyway
        struct winsize ttySize;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &ttySize);
        ttyWidth = std::max<int>(80, ttySize.ws_col);
#else
        // default for systems that do not implement the TIOCGWINSZ ioctl
        ttyWidth = 100;
#endif
    }

    return ttyWidth;
}

inline void printParamUsage_(std::ostream& os, const ParamInfo& paramInfo)
{
    std::string paramMessage, paramType, paramDescription;

    int ttyWidth = getTtyWidth_();

    // convert the CamelCase name to a command line --parameter-name.
    std::string cmdLineName = "-";
    const std::string camelCaseName = paramInfo.paramName;
    for (unsigned i = 0; i < camelCaseName.size(); ++i) {
        if (isupper(camelCaseName[i]))
            cmdLineName += "-";
        cmdLineName += static_cast<char>(std::tolower(camelCaseName[i]));
    }

    // assemble the printed output
    paramMessage = "    ";
    paramMessage += cmdLineName;

    // add the =VALUE_TYPE part
    bool isString = false;
    if (paramInfo.paramTypeName == Dune::className<std::string>()
        || paramInfo.paramTypeName == "const char *")
    {
        paramMessage += "=STRING";
        isString = true;
    }
    else if (paramInfo.paramTypeName == Dune::className<float>()
             || paramInfo.paramTypeName == Dune::className<double>()
             || paramInfo.paramTypeName == Dune::className<long double>()
#if HAVE_QUAD
             || paramInfo.paramTypeName == Dune::className<quad>()
#endif // HAVE_QUAD
        )
        paramMessage += "=SCALAR";
    else if (paramInfo.paramTypeName == Dune::className<int>()
             || paramInfo.paramTypeName == Dune::className<unsigned int>()
             || paramInfo.paramTypeName == Dune::className<short>()
             || paramInfo.paramTypeName == Dune::className<unsigned short>())
        paramMessage += "=INTEGER";
    else if (paramInfo.paramTypeName == Dune::className<bool>())
        paramMessage += "=BOOLEAN";
    else if (paramInfo.paramTypeName.empty()) {
        // the parameter is a flag. Do nothing!
    }
    else {
        // unknown type
        paramMessage += "=VALUE";
    }

    // fill up the up help string to the 50th character
    paramMessage += "  ";
    while (paramMessage.size() < 50)
        paramMessage += " ";


    // append the parameter usage string.
    paramMessage += paramInfo.usageString;

    // add the default value
    if (!paramInfo.paramTypeName.empty()) {
        if (paramMessage.back() != '.')
            paramMessage += '.';
        paramMessage += " Default: ";
        if (paramInfo.paramTypeName == "bool") {
            if (paramInfo.compileTimeValue == "0")
                paramMessage += "false";
            else
                paramMessage += "true";
        }
        else if (isString) {
            paramMessage += "\"";
            paramMessage += paramInfo.compileTimeValue;
            paramMessage += "\"";
        }
        else
            paramMessage += paramInfo.compileTimeValue;
    }

    paramMessage = breakLines_(paramMessage, /*indent=*/52, ttyWidth);
    paramMessage += "\n";

    // print everything
    os << paramMessage;
}

inline void getFlattenedKeyList_(std::list<std::string>& dest,
                                 const Dune::ParameterTree& tree,
                                 const std::string& prefix)
{
    // add the keys of the current sub-structure
    auto keyIt = tree.getValueKeys().begin();
    const auto& keyEndIt = tree.getValueKeys().end();
    for (; keyIt != keyEndIt; ++keyIt) {
        std::string newKey(prefix);
        newKey += *keyIt;
        dest.push_back(newKey);
    }

    // recursively add all substructure keys
    auto subStructIt = tree.getSubKeys().begin();
    const auto& subStructEndIt = tree.getSubKeys().end();
    for (; subStructIt != subStructEndIt; ++subStructIt) {
        std::string newPrefix(prefix);
        newPrefix += *subStructIt;
        newPrefix += ".";

        getFlattenedKeyList_(dest, tree.sub(*subStructIt), newPrefix);
    }
}

// print the values of a list of parameters
template <class TypeTag>
void printParamList_(std::ostream& os, const std::list<std::string>& keyList, bool printDefaults = false)
{
    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;

    const Dune::ParameterTree& tree = ParamsMeta::tree();

    auto keyIt = keyList.begin();
    const auto& keyEndIt = keyList.end();
    for (; keyIt != keyEndIt; ++keyIt) {
        const auto& paramInfo = ParamsMeta::registry().at(*keyIt);
        const std::string& defaultValue = paramInfo.compileTimeValue;
        std::string value = defaultValue;
        if (tree.hasKey(*keyIt))
            value = tree.get(*keyIt, "");
        os << *keyIt << "=\"" << value << "\"";
        if (printDefaults)
            os << " # default: \"" << defaultValue << "\"";
        os << "\n";
    }
}

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
template <class TypeTag>
void printUsage(const std::string& helpPreamble,
                const std::string& errorMsg = "",
                std::ostream& os = std::cerr,
                const bool showAll = false)
{
    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;

    if (errorMsg != "") {
        os << errorMsg << "\n"
           << "\n";
    }

    os << breakLines_(helpPreamble, /*indent=*/2, /*maxWidth=*/getTtyWidth_());
    os << "\n";

    os << "Recognized options:\n";

    if (!helpPreamble.empty()) {
        ParamInfo pInfo;
        pInfo.paramName = "h,--help";
        pInfo.usageString = "Print this help message and exit";
        printParamUsage_(os, pInfo);
        pInfo.paramName = "-help-all";
        pInfo.usageString = "Print all parameters, including obsolete, hidden and deprecated ones.";
        printParamUsage_(os, pInfo);
    }

    auto paramIt = ParamsMeta::registry().begin();
    const auto& paramEndIt = ParamsMeta::registry().end();
    for (; paramIt != paramEndIt; ++paramIt) {
        if (showAll || !paramIt->second.isHidden)
            printParamUsage_(os, paramIt->second);
    }
}

/// \cond 0
inline int noPositionalParameters_(std::set<std::string>&,
                                   std::string& errorMsg,
                                   int,
                                   const char** argv,
                                   int paramIdx,
                                   int)
{
    errorMsg = std::string("Illegal parameter \"")+argv[paramIdx]+"\".";
    return 0;
}

/// \endcond


inline void removeLeadingSpace_(std::string& s)
{
    unsigned i;
    for (i = 0; i < s.size(); ++ i)
        if (!std::isspace(s[i]))
            break;
    s = s.substr(i);
}

inline std::string transformKey_(const std::string& s,
                                 bool capitalizeFirstLetter = true,
                                 const std::string& errorPrefix = "")
{
    std::string result;

    if (s.empty())
        throw std::runtime_error(errorPrefix+"Empty parameter names are invalid");

    if (!std::isalpha(s[0]))
        throw std::runtime_error(errorPrefix+"Parameter name '" + s + "' is invalid: First character must be a letter");

    if (capitalizeFirstLetter)
        result += static_cast<char>(std::toupper(s[0]));
    else
        result += s[0];

    for (unsigned i = 1; i < s.size(); ++i) {
        if (s[i] == '-') {
            ++ i;
            if (s.size() <= i || !std::isalpha(s[i]))
                throw std::runtime_error(errorPrefix+"Invalid parameter name '" + s + "'");
            result += static_cast<char>(std::toupper(s[i]));
        }
        else if (!std::isalnum(s[i]))
            throw std::runtime_error(errorPrefix+"Invalid parameter name '" + s + "'");
        else
            result += s[i];
    }

    return result;
}

inline std::string parseKey_(std::string& s)
{
    unsigned i;
    for (i = 0; i < s.size(); ++ i)
        if (std::isspace(s[i]) || s[i] == '=')
            break;

    std::string ret = s.substr(0, i);
    s = s.substr(i);
    return ret;
}

// parse a quoted string
inline std::string parseQuotedValue_(std::string& s, const std::string& errorPrefix)
{
    if (s.empty() || s[0] != '"')
        throw std::runtime_error(errorPrefix+"Expected quoted string");

    std::string result;
    unsigned i = 1;
    for (; i < s.size(); ++i) {
        // handle escape characters
        if (s[i] == '\\') {
            ++ i;
            if (s.size() <= i)
                throw std::runtime_error(errorPrefix+"Unexpected end of quoted string");

            if (s[i] == 'n')
                result += '\n';
            else if (s[i] == 'r')
                result += '\r';
            else if (s[i] == 't')
                result += '\t';
            else if (s[i] == '"')
                result += '"';
            else if (s[i] == '\\')
                result += '\\';
            else
                throw std::runtime_error(errorPrefix+"Unknown escape character '\\" + s[i] + "'");
        }
        else if (s[i] == '"')
            break;
        else
            result += s[i];
    }

    s = s.substr(i+1);
    return result;
}

inline std::string parseUnquotedValue_(std::string& s, const std::string&)
{
    unsigned i;
    for (i = 0; i < s.size(); ++ i)
        if (std::isspace(s[i]))
            break;

    std::string ret = s.substr(0, i);
    s = s.substr(i);
    return ret;
}

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
template <class TypeTag, class PositionalArgumentCallback>
std::string parseCommandLineOptions(int argc,
                                    const char **argv,
                                    const std::string& helpPreamble = "",
                                    const PositionalArgumentCallback& posArgCallback = noPositionalParameters_)
{
    Dune::ParameterTree& paramTree = GetProp<TypeTag, Properties::ParameterMetaData>::tree();

    // handle the "--help" parameter
    if (!helpPreamble.empty()) {
        for (int i = 1; i < argc; ++i) {
            if (std::string("-h") == argv[i]
                || std::string("--help") == argv[i]) {
                printUsage<TypeTag>(helpPreamble, /*errorMsg=*/"", std::cout);
                return "Help called";
            }
            if (std::string("--help-all") == argv[i]) {
                printUsage<TypeTag>(helpPreamble, /*errorMsg=*/"", std::cout, true);
                return "Help called";
            }
        }
    }

    std::set<std::string> seenKeys;
    int numPositionalParams = 0;
    for (int i = 1; i < argc; ++i) {
        // All non-positional command line options need to start with '-'
        if (strlen(argv[i]) < 4
            || argv[i][0] != '-'
            || argv[i][1] != '-')
        {
            std::string errorMsg;
            int numHandled = posArgCallback(seenKeys, errorMsg, argc, argv, i, numPositionalParams);

            if (numHandled < 1) {
                std::ostringstream oss;

                if (!helpPreamble.empty())
                    printUsage<TypeTag>(helpPreamble, errorMsg, std::cerr);

                return errorMsg;
            }
            else {
                ++ numPositionalParams;
                i += numHandled - 1;
                continue;
            }
        }

        std::string paramName, paramValue;

        // read a --my-opt=abc option. This gets transformed
        // into the parameter "MyOpt" with the value being
        // "abc"

        // There is nothing after the '-'
        if (argv[i][2] == 0 || !std::isalpha(argv[i][2])) {
            std::ostringstream oss;
            oss << "Parameter name of argument " << i
                << " ('" << argv[i] << "') "
                << "is invalid because it does not start with a letter.";

            if (!helpPreamble.empty())
                printUsage<TypeTag>(helpPreamble, oss.str(), std::cerr);

            return oss.str();
        }

        // copy everything after the "--" into a separate string
        std::string s(argv[i] + 2);

        // parse argument
        paramName = transformKey_(parseKey_(s), /*capitalizeFirst=*/true);
        if (seenKeys.count(paramName) > 0) {
            std::string msg =
                std::string("Parameter '")+paramName+"' specified multiple times as a "
                "command line parameter";

            if (!helpPreamble.empty())
                printUsage<TypeTag>(helpPreamble, msg, std::cerr);
            return msg;
        }
        seenKeys.insert(paramName);

        if (s.empty() || s[0] != '=') {
            std::string msg =
                std::string("Parameter '")+paramName+"' is missing a value. "
                +" Please use "+argv[i]+"=value.";

            if (!helpPreamble.empty())
                printUsage<TypeTag>(helpPreamble, msg, std::cerr);
            return msg;
        }

        paramValue = s.substr(1);

        // Put the key=value pair into the parameter tree
        paramTree[paramName] = paramValue;
    }
    return "";
}

/*!
 * \ingroup Parameter
 * \brief Read the parameters from an INI-style file.
 *
 * This function does some basic syntax checks.
 */
template <class TypeTag>
void parseParameterFile(const std::string& fileName, bool overwrite = true)
{
    Dune::ParameterTree& paramTree = GetProp<TypeTag, Properties::ParameterMetaData>::tree();

    std::set<std::string> seenKeys;
    std::ifstream ifs(fileName);
    unsigned curLineNum = 0;
    while (ifs) {
        // string and file processing in c++ is quite blunt!
        std::string curLine;
        std::getline(ifs, curLine);
        curLineNum += 1;
        std::string errorPrefix = fileName+":"+std::to_string(curLineNum)+": ";

        // strip leading white space
        removeLeadingSpace_(curLine);

        // ignore empty and comment lines
        if (curLine.empty() || curLine[0] == '#' || curLine[0] == ';')
            continue;

        // TODO (?): support for parameter groups.

        // find the "key" of the key=value pair
        std::string key = parseKey_(curLine);
        std::string canonicalKey = transformKey_(key, /*capitalizeFirst=*/true, errorPrefix);

        if (seenKeys.count(canonicalKey) > 0)
            throw std::runtime_error(errorPrefix+"Parameter '"+canonicalKey+"' seen multiple times in the same file");
        seenKeys.insert(canonicalKey);

        // deal with the equals sign
        removeLeadingSpace_(curLine);
        if (curLine.empty() || curLine[0] != '=')
            std::runtime_error(errorPrefix+"Syntax error, expecting 'key=value'");

        curLine = curLine.substr(1);
        removeLeadingSpace_(curLine);

        if (curLine.empty() || curLine[0] == '#' || curLine[0] == ';')
            std::runtime_error(errorPrefix+"Syntax error, expecting 'key=value'");

        // get the value
        std::string value;
        if (curLine[0] == '"')
            value = parseQuotedValue_(curLine, errorPrefix);
        else
            value = parseUnquotedValue_(curLine, errorPrefix);

        // ignore trailing comments
        removeLeadingSpace_(curLine);
        if (!curLine.empty() && curLine[0] != '#' && curLine[0] != ';')
            std::runtime_error(errorPrefix+"Syntax error, expecting 'key=value'");

        // all went well, add the parameter to the database object
        if (overwrite || !paramTree.hasKey(canonicalKey))
            paramTree[canonicalKey] = value;
    }
}

/*!
 * \ingroup Parameter
 * \brief Print values of the run-time parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 */
template <class TypeTag>
void printValues(std::ostream& os = std::cout)
{
    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;

    const Dune::ParameterTree& tree = ParamsMeta::tree();

    std::list<std::string> runTimeAllKeyList;
    std::list<std::string> runTimeKeyList;
    std::list<std::string> unknownKeyList;

    getFlattenedKeyList_(runTimeAllKeyList, tree);
    auto keyIt = runTimeAllKeyList.begin();
    const auto& keyEndIt = runTimeAllKeyList.end();
    for (; keyIt != keyEndIt; ++keyIt) {
        if (ParamsMeta::registry().find(*keyIt) == ParamsMeta::registry().end()) {
            // key was not registered by the program!
            unknownKeyList.push_back(*keyIt);
        }
        else {
            // the key was specified at run-time
            runTimeKeyList.push_back(*keyIt);
        }
    }

    // loop over all registered parameters
    std::list<std::string> compileTimeKeyList;
    auto paramInfoIt = ParamsMeta::registry().begin();
    const auto& paramInfoEndIt = ParamsMeta::registry().end();
    for (; paramInfoIt != paramInfoEndIt; ++paramInfoIt) {
        // check whether the key was specified at run-time
        const auto& keyName = paramInfoIt->first;
        if (tree.hasKey(keyName))
            continue;
        else
            compileTimeKeyList.push_back(keyName);
    }

    // report the values of all registered (and unregistered)
    // parameters
    if (runTimeKeyList.size() > 0) {
        os << "# [known parameters which were specified at run-time]\n";
        printParamList_<TypeTag>(os, runTimeKeyList, /*printDefaults=*/true);
    }

    if (compileTimeKeyList.size() > 0) {
        os << "# [parameters which were specified at compile-time]\n";
        printParamList_<TypeTag>(os, compileTimeKeyList, /*printDefaults=*/false);
    }

    if (unknownKeyList.size() > 0) {
        os << "# [unused run-time specified parameters]\n";
        auto unusedKeyIt = unknownKeyList.begin();
        const auto& unusedKeyEndIt = unknownKeyList.end();
        for (; unusedKeyIt != unusedKeyEndIt; ++unusedKeyIt) {
            os << *unusedKeyIt << "=\"" << tree.get(*unusedKeyIt, "") << "\"\n" << std::flush;
        }
    }
}

/*!
 * \ingroup Parameter
 * \brief Print the list of unused run-time parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 *
 * \return true if something was printed
 */
template <class TypeTag>
bool printUnused(std::ostream& os = std::cout)
{
    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;

    const Dune::ParameterTree& tree = ParamsMeta::tree();
    std::list<std::string> runTimeAllKeyList;
    std::list<std::string> unknownKeyList;

    getFlattenedKeyList_(runTimeAllKeyList, tree);
    auto keyIt = runTimeAllKeyList.begin();
    const auto& keyEndIt = runTimeAllKeyList.end();
    for (; keyIt != keyEndIt; ++keyIt) {
        if (ParamsMeta::registry().find(*keyIt) == ParamsMeta::registry().end()) {
            // key was not registered by the program!
            unknownKeyList.push_back(*keyIt);
        }
    }

    if (unknownKeyList.size() > 0) {
        os << "# [unused run-time specified parameters]\n";
        auto unusedKeyIt = unknownKeyList.begin();
        const auto& unusedKeyEndIt = unknownKeyList.end();
        for (; unusedKeyIt != unusedKeyEndIt; ++unusedKeyIt) {
            os << *unusedKeyIt << "=\"" << tree.get(*unusedKeyIt, "") << "\"\n" << std::flush;
        }
        return true;
    }
    return false;
}

//! \cond SKIP_THIS
template <class TypeTag>
class Param
{
    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;

public:
    template <class ParamType, class PropTag>
    static ParamType get(const char *propTagName,
                    const char *paramName,
                    bool errorIfNotRegistered = true)
    {
        return retrieve_<ParamType>(propTagName, paramName, getPropValue<TypeTag, PropTag>(), errorIfNotRegistered);
    }
    
    template <class ParamType>
    static ParamType get(const char *propTagName,
                         const char *paramName,
                         const ParamType& defaultValue,
                         bool errorIfNotRegistered = true)
    {
        return retrieve_<ParamType>(propTagName, paramName, defaultValue, errorIfNotRegistered);
    }
    
    static void clear()
    {
        ParamsMeta::clear();
    }

    template <class ParamType>
    static bool isSet(const char *propTagName OPM_OPTIM_UNUSED,
                      const char *paramName OPM_OPTIM_UNUSED,
                      bool errorIfNotRegistered = true)
    {

#ifndef NDEBUG
        // make sure that the parameter is used consistently. since
        // this is potentially quite expensive, it is only done if
        // debugging code is not explicitly turned off.
        check_(Dune::className<ParamType>(), propTagName, paramName);
#endif

        if (errorIfNotRegistered) {
            if (ParamsMeta::registrationOpen())
                throw std::runtime_error("Parameters can only checked after _all_ of them have "
                                         "been registered.");

            if (ParamsMeta::registry().find(paramName) == ParamsMeta::registry().end())
                throw std::runtime_error("Accessing parameter "+std::string(paramName)
                                         +" without prior registration is not allowed.");
        }

        std::string canonicalName(paramName);

        // check whether the parameter is in the parameter tree
        return ParamsMeta::tree().hasKey(canonicalName);
    }


private:
    struct Blubb
    {
        std::string propertyName;
        std::string paramTypeName;
        std::string groupName;

        Blubb& operator=(const Blubb& b)
        {
            propertyName = b.propertyName;
            paramTypeName = b.paramTypeName;
            groupName = b.groupName;
            return *this;
        }
    };

    static void check_(const std::string& paramTypeName,
                       const std::string& propertyName,
                       const char *paramName)
    {
        using StaticData = std::unordered_map<std::string, Blubb>;
        static StaticData staticData;

        typename StaticData::iterator it = staticData.find(paramName);
        Blubb *b;
        if (it == staticData.end()) {
            Blubb a;
            a.propertyName = propertyName;
            a.paramTypeName = paramTypeName;
            staticData[paramName] = a;
            b = &staticData[paramName];
        }
        else
            b = &(it->second);

        if (b->propertyName != propertyName) {
            throw std::logic_error("GET_*_PARAM for parameter '"+std::string(paramName)
                                   +"' called for at least two different properties ('"
                                   +b->propertyName+"' and '"+propertyName+"')");
        }

        if (b->paramTypeName != paramTypeName) {
            throw std::logic_error("GET_*_PARAM for parameter '"+std::string(paramName)
                                   +"' called with at least two different types ("
                                   +b->paramTypeName+" and "+paramTypeName+")");
        }
    }

    template <class ParamType>
    static ParamType retrieve_(const char OPM_OPTIM_UNUSED *propTagName,
                               const char *paramName,
                               const ParamType& defaultValue,
                               bool errorIfNotRegistered = true)
    {
#ifndef NDEBUG
        // make sure that the parameter is used consistently. since
        // this is potentially quite expensive, it is only done if
        // debugging code is not explicitly turned off.
        check_(Dune::className<ParamType>(), propTagName, paramName);
#endif

        if (errorIfNotRegistered) {
            if (ParamsMeta::registrationOpen())
                throw std::runtime_error("Parameters can only retieved after _all_ of them have "
                                         "been registered.");

            if (ParamsMeta::registry().find(paramName) == ParamsMeta::registry().end())
                throw std::runtime_error("Accessing parameter "+std::string(paramName)
                                         +" without prior registration is not allowed.");
        }

        // prefix the parameter name by the model's GroupName. E.g. If
        // the model specifies its group name to be 'Stokes', in an
        // INI file this would result in something like:
        //
        // [Stokes]
        // NewtonWriteConvergence = true
        std::string canonicalName(paramName);

        // retrieve actual parameter from the parameter tree
        return ParamsMeta::tree().template get<ParamType>(canonicalName, defaultValue);
    }
};

template <class TypeTag, class ParamType, class PropTag>
const ParamType get(const char *propTagName, const char *paramName, bool errorIfNotRegistered)
{
    return Param<TypeTag>::template get<ParamType, PropTag>(propTagName,
                                                            paramName,
                                                            errorIfNotRegistered);
}

template <class TypeTag, class ParamType>
const ParamType get(const char *propTagName, const char *paramName, const ParamType& defaultValue, bool errorIfNotRegistered)
{
    return Param<TypeTag>::template get<ParamType>(propTagName, paramName, defaultValue, errorIfNotRegistered);
}

template <class TypeTag, class Container>
void getLists(Container& usedParams, Container& unusedParams)
{
    usedParams.clear();
    unusedParams.clear();

    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;
    if (ParamsMeta::registrationOpen())
        throw std::runtime_error("Parameter lists can only retieved after _all_ of them have "
                                 "been registered.");

    // get all parameter keys
    std::list<std::string> allKeysList;
    const auto& paramTree = ParamsMeta::tree();
    getFlattenedKeyList_(allKeysList, paramTree);

    for (const auto& key : allKeysList) {
        if (ParamsMeta::registry().find(key) == ParamsMeta::registry().end()) {
            // key was not registered
            unusedParams.emplace_back(key, paramTree[key]);
        }
        else {
            // key was registered
            usedParams.emplace_back(key, paramTree[key]);
        }
    }
}

template <class TypeTag>
void reset()
{
    return Param<TypeTag>::clear();
}

template <class TypeTag, class ParamType>
bool isSet(const char *propTagName, const char *paramName, bool errorIfNotRegistered = true)
{
    return Param<TypeTag>::template isSet<ParamType>(propTagName,
                                                     paramName,
                                                     errorIfNotRegistered);
}

template <class TypeTag, class ParamType>
void registerParam(const char *paramName, const char *propertyName, const ParamType& defaultValue, const char *usageString)
{
    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;
    if (!ParamsMeta::registrationOpen())
        throw std::logic_error("Parameter registration was already closed before "
                               "the parameter '"+std::string(paramName)+"' was registered.");

    ParamsMeta::registrationFinalizers().emplace_back(
        new ParamRegFinalizer_<TypeTag, ParamType>(paramName, defaultValue));

    ParamInfo paramInfo;
    paramInfo.paramName = paramName;
    paramInfo.paramTypeName = Dune::className<ParamType>();
    std::string tmp = Dune::className<TypeTag>();
    tmp.replace(0, strlen("Opm::Properties::TTag::"), "");
    paramInfo.propertyName = propertyName;
    paramInfo.usageString = usageString;
    std::ostringstream oss;
    oss << defaultValue;
    paramInfo.compileTimeValue = oss.str();
    paramInfo.isHidden = false;
    if (ParamsMeta::registry().find(paramName) != ParamsMeta::registry().end()) {
        // allow to register a parameter twice, but only if the
        // parameter name, type and usage string are exactly the same.
        if (ParamsMeta::registry().at(paramName) == paramInfo)
            return;
        throw std::logic_error("Parameter "+std::string(paramName)
                               +" registered twice with non-matching characteristics.");
    }

    ParamsMeta::mutableRegistry()[paramName] = paramInfo;
}

template <class TypeTag, class ParamType>
void hideParam(const char *paramName, const ParamType&)
{
    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;
    if (!ParamsMeta::registrationOpen())
        throw std::logic_error("Parameter '"+std::string(paramName)+"' declared as hidden"
                               " when parameter registration was already closed.");

    auto paramInfoIt = ParamsMeta::mutableRegistry().find(paramName);
    if (paramInfoIt == ParamsMeta::mutableRegistry().end())
        throw std::logic_error("Tried to declare unknown parameter '"
                               +std::string(paramName)+"' hidden.");

    auto& paramInfo = paramInfoIt->second;
    paramInfo.isHidden = true;
}

template <class TypeTag>
void endParamRegistration()
{
    using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;
    if (!ParamsMeta::registrationOpen())
        throw std::logic_error("Parameter registration was already closed. It is only possible "
                               "to close it once.");

    ParamsMeta::registrationOpen() = false;

    // loop over all parameters and retrieve their values to make sure
    // that there is no syntax error
    auto pIt = ParamsMeta::registrationFinalizers().begin();
    const auto& pEndIt = ParamsMeta::registrationFinalizers().end();
    for (; pIt != pEndIt; ++pIt)
        (*pIt)->retrieve();
    ParamsMeta::registrationFinalizers().clear();
}
//! \endcond

} // namespace Parameters
} // namespace Opm

#endif
