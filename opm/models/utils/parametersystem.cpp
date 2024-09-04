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
#include <config.h>
#include <opm/models/utils/parametersystem.hpp>

#if HAVE_QUAD
#include <opm/material/common/quad.hpp>
#endif

namespace {

void getFlattenedKeyList(std::list<std::string>& dest,
                         const Dune::ParameterTree& tree,
                         const std::string& prefix = "")
{
    // add the keys of the current sub-structure
    for (const auto& valueKey : tree.getValueKeys()) {
        std::string newKey(prefix + valueKey);
        dest.push_back(newKey);
    }

    // recursively add all substructure keys
    for (const auto& subKey : tree.getSubKeys()) {
        std::string newPrefix(prefix + subKey + '.');
        getFlattenedKeyList(dest, tree.sub(subKey), newPrefix);
    }
}

std::string parseKey(std::string& s)
{
    unsigned i;
    for (i = 0; i < s.size(); ++ i)
        if (std::isspace(s[i]) || s[i] == '=')
            break;

    std::string ret = s.substr(0, i);
    s = s.substr(i);
    return ret;
}

std::string parseQuotedValue(std::string& s, const std::string& errorPrefix)
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

std::string parseUnquotedValue(std::string& s, const std::string&)
{
    unsigned i;
    for (i = 0; i < s.size(); ++ i)
        if (std::isspace(s[i]))
            break;

    std::string ret = s.substr(0, i);
    s = s.substr(i);
    return ret;
}

void printParamList(std::ostream& os,
                    const std::list<std::string>& keyList,
                    bool printDefaults = false)
{
    const Dune::ParameterTree& tree = Opm::Parameters::MetaData::tree();

    for (const auto& key : keyList) {
        const auto& paramInfo = Opm::Parameters::MetaData::registry().at(key);
        const std::string& defaultValue = paramInfo.defaultValue;
        std::string value = defaultValue;
        if (tree.hasKey(key))
            value = tree.get(key, "");
        os << key << "=\"" << value << "\"";
        if (printDefaults)
            os << " # default: \"" << defaultValue << "\"";
        os << "\n";
    }
}

void printParamUsage(std::ostream& os,
                     const Opm::Parameters::ParamInfo& paramInfo)
{
    std::string paramMessage, paramType, paramDescription;

    int ttyWidth = Opm::Parameters::getTtyWidth();

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
            if (paramInfo.defaultValue == "0")
                paramMessage += "false";
            else
                paramMessage += "true";
        }
        else if (isString) {
            paramMessage += "\"";
            paramMessage += paramInfo.defaultValue;
            paramMessage += "\"";
        }
        else
            paramMessage += paramInfo.defaultValue;
    }

    paramMessage = Opm::Parameters::breakLines(paramMessage, /*indent=*/52, ttyWidth);
    paramMessage += "\n";

    // print everything
    os << paramMessage;
}

void removeLeadingSpace(std::string& s)
{
    unsigned i;
    for (i = 0; i < s.size(); ++ i)
        if (!std::isspace(s[i]))
            break;
    s = s.substr(i);
}

std::string transformKey(const std::string& s,
                         bool capitalizeFirstLetter,
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

} // anonymous namespace

namespace Opm::Parameters {

namespace detail {

template<class ParamType>
ParamType Get_(const std::string& paramName, ParamType defaultValue,
               bool errorIfNotRegistered)
{
    if (errorIfNotRegistered) {
        if (MetaData::registrationOpen())
            throw std::runtime_error("Parameters can only retrieved after _all_ of them have "
                                     "been registered.");

        if (MetaData::registry().find(paramName) == MetaData::registry().end()) {
            throw std::runtime_error("Accessing parameter " + paramName
                                     +" without prior registration is not allowed.");
        }
    }

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

void Hide_(const std::string& paramName)
{
    if (!MetaData::registrationOpen()) {
        throw std::logic_error("Parameter '" + paramName + "' declared as hidden"
                               " when parameter registration was already closed.");
    }

    auto paramInfoIt = MetaData::mutableRegistry().find(paramName);
    if (paramInfoIt == MetaData::mutableRegistry().end()) {
        throw std::logic_error("Tried to declare unknown parameter '"
                               + paramName + "' hidden.");
    }

    auto& paramInfo = paramInfoIt->second;
    paramInfo.isHidden = true;
}

bool IsSet_(const std::string& paramName, bool errorIfNotRegistered)
{
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

void Register_(const std::string& paramName,
               const std::string& paramTypeName,
               const std::string& defaultValue,
               const char* usageString)
{
    if (!MetaData::registrationOpen()) {
        throw std::logic_error("Parameter registration was already closed before "
                               "the parameter '" + paramName + "' was registered.");
    }

    ParamInfo paramInfo;
    paramInfo.paramName = paramName;
    paramInfo.paramTypeName = paramTypeName;
    paramInfo.usageString = usageString;
    paramInfo.defaultValue = defaultValue;
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

void SetDefault_(const std::string& paramName,
                 const std::string& paramValue)
{
    if (MetaData::registry().find(paramName) == MetaData::registry().end()) {
        throw std::runtime_error("Accessing parameter " + paramName +
                                 " without prior registration is not allowed.");
    }
    MetaData::mutableRegistry()[paramName].defaultValue = paramValue;
}

}

bool ParamInfo::operator==(const ParamInfo& other) const
{
    return other.paramName == paramName
        && other.paramTypeName == paramTypeName
        && other.typeTagName == typeTagName
        && other.usageString == usageString;
}

std::string breakLines(const std::string& msg,
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

void reset()
{
    MetaData::clear();
}

bool IsRegistrationOpen()
{
    return MetaData::registrationOpen();
}

void endRegistration()
{
    if (!MetaData::registrationOpen()) {
        throw std::logic_error("Parameter registration was already closed. It is only possible "
                               "to close it once.");
    }

    MetaData::registrationOpen() = false;
}

void getLists(std::vector<Parameter>& usedParams,
              std::vector<Parameter>& unusedParams)
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

void printUsage(const std::string& helpPreamble,
                const std::string& errorMsg,
                std::ostream& os,
                const bool showAll)
{
    if (!errorMsg.empty()) {
        os << errorMsg << "\n\n";
    }

    os << breakLines(helpPreamble, /*indent=*/2, /*maxWidth=*/getTtyWidth());
    os << "\n";

    os << "Recognized options:\n";

    if (!helpPreamble.empty()) {
        ParamInfo pInfo;
        pInfo.paramName = "h,--help";
        pInfo.usageString = "Print this help message and exit";
        printParamUsage(os, pInfo);
        pInfo.paramName = "-help-all";
        pInfo.usageString = "Print all parameters, including obsolete, hidden and deprecated ones.";
        printParamUsage(os, pInfo);
    }

    for (const auto& param : MetaData::registry()) {
        if (showAll || !param.second.isHidden)
            printParamUsage(os, param.second);
    }
}

void parseParameterFile(const std::string& fileName, bool overwrite)
{
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
        removeLeadingSpace(curLine);

        // ignore empty and comment lines
        if (curLine.empty() || curLine[0] == '#' || curLine[0] == ';')
            continue;

        // TODO (?): support for parameter groups.

        // find the "key" of the key=value pair
        std::string key = parseKey(curLine);
        std::string canonicalKey = transformKey(key, /*capitalizeFirst=*/true, errorPrefix);

        if (seenKeys.count(canonicalKey) > 0)
            throw std::runtime_error(errorPrefix+"Parameter '"+canonicalKey+"' seen multiple times in the same file");
        seenKeys.insert(canonicalKey);

        // deal with the equals sign
        removeLeadingSpace(curLine);
        if (curLine.empty() || curLine[0] != '=')
            std::runtime_error(errorPrefix+"Syntax error, expecting 'key=value'");

        curLine = curLine.substr(1);
        removeLeadingSpace(curLine);

        if (curLine.empty() || curLine[0] == '#' || curLine[0] == ';')
            std::runtime_error(errorPrefix+"Syntax error, expecting 'key=value'");

        // get the value
        std::string value;
        if (curLine[0] == '"')
            value = parseQuotedValue(curLine, errorPrefix);
        else
            value = parseUnquotedValue(curLine, errorPrefix);

        // ignore trailing comments
        removeLeadingSpace(curLine);
        if (!curLine.empty() && curLine[0] != '#' && curLine[0] != ';')
            std::runtime_error(errorPrefix+"Syntax error, expecting 'key=value'");

        // all went well, add the parameter to the database object
        if (overwrite || !MetaData::tree().hasKey(canonicalKey)) {
            MetaData::tree()[canonicalKey] = value;
        }
    }
}

std::string parseCommandLineOptions(int argc,
                                    const char **argv,
                                    const PositionalArgumentCallback& posArgCallback,
                                    const std::string& helpPreamble)
{
    // handle the "--help" parameter
    if (!helpPreamble.empty()) {
        for (int i = 1; i < argc; ++i) {
            if (std::string("-h") == argv[i]
                || std::string("--help") == argv[i]) {
                printUsage(helpPreamble, /*errorMsg=*/"", std::cout);
                return "Help called";
            }
            if (std::string("--help-all") == argv[i]) {
                printUsage(helpPreamble, /*errorMsg=*/"", std::cout, true);
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
            int numHandled = posArgCallback([](const std::string& k, const std::string& v)
                                            {
                                                MetaData::tree()[k] = v;
                                            }, seenKeys, errorMsg,
                                            argc, argv, i, numPositionalParams);

            if (numHandled < 1) {
                std::ostringstream oss;

                if (!helpPreamble.empty())
                    printUsage(helpPreamble, errorMsg, std::cerr);

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
                printUsage(helpPreamble, oss.str(), std::cerr);

            return oss.str();
        }

        // copy everything after the "--" into a separate string
        std::string s(argv[i] + 2);

        // parse argument
        paramName = transformKey(parseKey(s), /*capitalizeFirst=*/true);
        if (seenKeys.count(paramName) > 0) {
            std::string msg =
                std::string("Parameter '")+paramName+"' specified multiple times as a "
                "command line parameter";

            if (!helpPreamble.empty())
                printUsage(helpPreamble, msg, std::cerr);
            return msg;
        }
        seenKeys.insert(paramName);

        if (s.empty() || s[0] != '=') {
            std::string msg =
                std::string("Parameter '")+paramName+"' is missing a value. "
                +" Please use "+argv[i]+"=value.";

            if (!helpPreamble.empty())
                printUsage(helpPreamble, msg, std::cerr);
            return msg;
        }

        paramValue = s.substr(1);

        // Put the key=value pair into the parameter tree
        MetaData::tree()[paramName] = paramValue;
    }
    return "";
}

void printValues(std::ostream& os)
{
    std::list<std::string> runTimeAllKeyList;
    std::list<std::string> runTimeKeyList;
    std::list<std::string> unknownKeyList;

    getFlattenedKeyList(runTimeAllKeyList, MetaData::tree());
    for (const auto& key : runTimeAllKeyList) {
        if (MetaData::registry().find(key) == MetaData::registry().end()) {
            // key was not registered by the program!
            unknownKeyList.push_back(key);
        }
        else {
            // the key was specified at run-time
            runTimeKeyList.push_back(key);
        }
    }

    // loop over all registered parameters
    std::list<std::string> compileTimeKeyList;
    for (const auto& reg : MetaData::registry()) {
        // check whether the key was specified at run-time
        if (MetaData::tree().hasKey(reg.first)) {
            continue;
        } else  {
            compileTimeKeyList.push_back(reg.first);
        }
    }

    // report the values of all registered (and unregistered)
    // parameters
    if (runTimeKeyList.size() > 0) {
        os << "# [known parameters which were specified at run-time]\n";
        printParamList(os, runTimeKeyList, /*printDefaults=*/true);
    }

    if (compileTimeKeyList.size() > 0) {
        os << "# [parameters which were specified at compile-time]\n";
        printParamList(os, compileTimeKeyList, /*printDefaults=*/false);
    }

    if (unknownKeyList.size() > 0) {
        os << "# [unused run-time specified parameters]\n";
        for (const auto& unused : unknownKeyList) {
            os << unused << "=\"" << MetaData::tree().get(unused, "") << "\"\n" << std::flush;
        }
    }
}

bool printUnused(std::ostream& os)
{
    std::list<std::string> runTimeAllKeyList;
    std::list<std::string> unknownKeyList;

    getFlattenedKeyList(runTimeAllKeyList, MetaData::tree());
    for (const auto& key : runTimeAllKeyList) {
        if (MetaData::registry().find(key) == MetaData::registry().end()) {
            // key was not registered by the program!
            unknownKeyList.push_back(key);
        }
    }

    if (unknownKeyList.size() > 0) {
        os << "# [unused run-time specified parameters]\n";
        for (const auto& unused : unknownKeyList) {
            os << unused << "=\""
               << MetaData::tree().get(unused, "") << "\"\n" << std::flush;
        }
        return true;
    }
    return false;
}

int getTtyWidth()
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

namespace detail {

template bool Get_(const std::string&, bool, bool);
template double Get_(const std::string&, double, bool);
template float Get_(const std::string&, float, bool);
template int Get_(const std::string&, int, bool);
template long Get_(const std::string&, long, bool);
template std::string Get_(const std::string&, std::string, bool);
template unsigned Get_(const std::string&, unsigned, bool);

#if HAVE_QUAD
template quad Get_(const std::string&, quad, bool);
#endif

}

} // namespace Opm::Parameters
