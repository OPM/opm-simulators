/*
  Copyright 2019, 2020 SINTEF Digital, Mathematics and Cybernetics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <external/tomlc17/tomlc17.h>

#include <boost/property_tree/json_parser.hpp>

#include <algorithm>
#include <charconv>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include <fmt/format.h>

namespace {

using PTree = boost::property_tree::ptree;

/// RAII wrapper for a tomlc17 parse result.
class TomlResult
{
public:
    explicit TomlResult(toml_result_t result) : result_(result) {}
    ~TomlResult() { toml_free(result_); }

    TomlResult(const TomlResult&) = delete;
    TomlResult& operator=(const TomlResult&) = delete;

    const toml_result_t& get() const { return result_; }

private:
    toml_result_t result_;
};

std::string tomlTimestampString(const toml_datum_t& datum)
{
    const auto& ts = datum.u.ts;
    const auto date = fmt::format("{:04}-{:02}-{:02}", ts.year, ts.month, ts.day);
    auto time = fmt::format("{:02}:{:02}:{:02}", ts.hour, ts.minute, ts.second);
    if (ts.usec != 0) {
        time += fmt::format(".{:06}", ts.usec);
    }

    switch (datum.type) {
    case TOML_DATE:
        return date;
    case TOML_TIME:
        return time;
    case TOML_DATETIME:
        return date + 'T' + time;
    case TOML_DATETIMETZ: {
        if (ts.tz == 0) {
            return date + 'T' + time + 'Z';
        }
        const int absTz = ts.tz < 0 ? -ts.tz : ts.tz;
        return fmt::format("{}T{}{}{:02}:{:02}", date, time,
                           ts.tz < 0 ? '-' : '+', absTz / 60, absTz % 60);
    }
    default:
        return {};
    }
}

/// Convert parsed TOML data into a Boost property tree, using the same
/// conventions as Boost's JSON parser: all values are stored as strings
/// and array elements are stored as children with empty keys.
PTree tomlToPtree(const toml_datum_t& datum)
{
    auto node = PTree{};

    switch (datum.type) {
    case TOML_TABLE:
        for (int i = 0; i < datum.u.tab.size; ++i) {
            node.push_back({ std::string(datum.u.tab.key[i], datum.u.tab.len[i]),
                             tomlToPtree(datum.u.tab.value[i]) });
        }
        break;

    case TOML_ARRAY:
        for (int i = 0; i < datum.u.arr.size; ++i) {
            node.push_back({ std::string{}, tomlToPtree(datum.u.arr.elem[i]) });
        }
        break;

    case TOML_STRING:
        node.data() = std::string(datum.u.str.ptr, datum.u.str.len);
        break;

    case TOML_INT64:
        node.data() = std::to_string(datum.u.int64);
        break;

    case TOML_FP64:
        node.data() = fmt::format("{}", datum.u.fp64);
        break;

    case TOML_BOOLEAN:
        node.data() = datum.u.boolean ? "true" : "false";
        break;

    case TOML_DATE:
    case TOML_TIME:
    case TOML_DATETIME:
    case TOML_DATETIMETZ:
        node.data() = tomlTimestampString(datum);
        break;

    case TOML_UNKNOWN:
    default:
        break;
    }

    return node;
}

PTree parseToml(const toml_result_t& result, std::string_view source)
{
    if (! result.ok) {
        throw std::runtime_error {
            fmt::format("Failed to parse TOML {}: {}", source, result.errmsg)
        };
    }

    return tomlToPtree(result.toptab);
}

PTree readTomlFile(const std::string& fileName)
{
    const auto result = TomlResult { toml_parse_file_ex(fileName.c_str()) };
    return parseToml(result.get(), fmt::format("file '{}'", fileName));
}

PTree readTomlString(const std::string& toml)
{
    if (toml.size() > static_cast<std::size_t>(INT_MAX)) {
        throw std::invalid_argument { "TOML string is too large" };
    }

    const auto result = TomlResult {
        toml_parse(toml.data(), static_cast<int>(toml.size()))
    };
    return parseToml(result.get(), "string");
}

// ---------------------------------------------------------------------------
// TOML output
// ---------------------------------------------------------------------------

bool isArrayNode(const PTree& node)
{
    return !node.empty() &&
        std::all_of(node.begin(), node.end(),
                    [](const auto& child) { return child.first.empty(); });
}

bool isTableNode(const PTree& node)
{
    return !node.empty() && !isArrayNode(node);
}

bool isBareKey(std::string_view key)
{
    return !key.empty() &&
        std::all_of(key.begin(), key.end(), [](const char c)
        {
            return ((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z')) ||
                   ((c >= '0') && (c <= '9')) || (c == '_') || (c == '-');
        });
}

std::string quotedTomlString(std::string_view s)
{
    std::string out = "\"";
    for (const char c : s) {
        switch (c) {
        case '"':  out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\b': out += "\\b";  break;
        case '\t': out += "\\t";  break;
        case '\n': out += "\\n";  break;
        case '\f': out += "\\f";  break;
        case '\r': out += "\\r";  break;
        default:
            if ((static_cast<unsigned char>(c) < 0x20) || (c == 0x7f)) {
                out += fmt::format("\\u{:04X}", static_cast<unsigned char>(c));
            }
            else {
                out += c;
            }
        }
    }
    out += '"';
    return out;
}

std::string tomlKey(std::string_view key)
{
    return isBareKey(key) ? std::string{key} : quotedTomlString(key);
}

bool isDigit(const char c)
{
    return (c >= '0') && (c <= '9');
}

std::optional<std::string> tomlIntegerValue(std::string_view s)
{
    auto digits = s;
    if (!digits.empty() && ((digits.front() == '-') || (digits.front() == '+'))) {
        digits.remove_prefix(1);
    }

    // TOML does not allow leading zeros in integers.
    if (digits.empty() || !std::all_of(digits.begin(), digits.end(), isDigit) ||
        ((digits.size() > 1) && (digits.front() == '0')))
    {
        return std::nullopt;
    }

    std::int64_t value{};
    const auto* first = s.data() + (s.front() == '+' ? 1 : 0);
    const auto [ptr, ec] = std::from_chars(first, s.data() + s.size(), value);
    if ((ec != std::errc{}) || (ptr != s.data() + s.size())) {
        return std::nullopt;
    }

    return std::to_string(value);
}

std::optional<std::string> tomlFloatValue(std::string_view s)
{
    auto mantissa = s;
    if (!mantissa.empty() && (mantissa.front() == '-')) {
        mantissa.remove_prefix(1);
    }

    // Only consider strings that unambiguously represent finite floating
    // point numbers, e.g., "0.5", "-1.5e-3", or "1e10".
    if (mantissa.empty() || !isDigit(mantissa.front()) ||
        (mantissa.find_first_of(".eE") == std::string_view::npos))
    {
        return std::nullopt;
    }

    double value{};
    const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
    if ((ec != std::errc{}) || (ptr != s.data() + s.size())) {
        return std::nullopt;
    }

    // Shortest representation that round-trips, ensuring a TOML float.
    auto out = fmt::format("{}", value);
    if (out.find_first_of(".eE") == std::string::npos) {
        out += ".0";
    }

    return out;
}

std::string tomlScalarValue(const std::string& s)
{
    if ((s == "true") || (s == "false")) {
        return s;
    }

    if (auto i = tomlIntegerValue(s); i.has_value()) {
        return *i;
    }

    if (auto f = tomlFloatValue(s); f.has_value()) {
        return *f;
    }

    return quotedTomlString(s);
}

std::string tomlInlineValue(const PTree& node);

std::string tomlInlineTable(const PTree& node)
{
    std::string out = "{ ";
    bool first = true;
    for (const auto& [key, child] : node) {
        if (!first) {
            out += ", ";
        }
        out += tomlKey(key) + " = " + tomlInlineValue(child);
        first = false;
    }
    out += " }";
    return out;
}

std::string tomlInlineValue(const PTree& node)
{
    if (isArrayNode(node)) {
        std::string out = "[";
        bool first = true;
        for (const auto& elem : node) {
            out += first ? "" : ", ";
            out += tomlInlineValue(elem.second);
            first = false;
        }
        out += "]";
        return out;
    }

    if (isTableNode(node)) {
        return tomlInlineTable(node);
    }

    return tomlScalarValue(node.data());
}

void writeTomlComment(std::ostream& os, const std::string& comment)
{
    if (comment.empty()) {
        return;
    }

    auto lines = std::istringstream { comment };
    for (std::string line; std::getline(lines, line); ) {
        os << (line.empty() ? "#" : "# " + line) << '\n';
    }
}

std::string joinKey(const std::string& parent, const std::string& key)
{
    return parent.empty() ? key : parent + '.' + key;
}

void writeTomlTable(std::ostream& os,
                    const PTree& table,
                    const std::string& path,
                    const std::string& header,
                    const Opm::PropertyTree::TomlCommentProvider& comments,
                    bool& needsSeparator)
{
    // Key/value pairs (including arrays) must precede sub tables, as they
    // would otherwise be attributed to the last sub table.
    for (const auto& [key, child] : table) {
        if (isTableNode(child)) {
            continue;
        }

        const auto fullKey = joinKey(path, key);
        const auto comment = comments ? comments(fullKey) : std::string{};
        if (!comment.empty() && needsSeparator) {
            os << '\n';
        }
        writeTomlComment(os, comment);
        os << tomlKey(key) << " = " << tomlInlineValue(child) << '\n';
        needsSeparator = true;
    }

    for (const auto& [key, child] : table) {
        if (!isTableNode(child)) {
            continue;
        }

        const auto fullKey = joinKey(path, key);
        const auto childHeader = header.empty()
            ? tomlKey(key) : header + '.' + tomlKey(key);

        if (needsSeparator) {
            os << '\n';
        }
        writeTomlComment(os, comments ? comments(fullKey) : std::string{});
        os << '[' << childHeader << "]\n";
        needsSeparator = false;

        writeTomlTable(os, child, fullKey, childHeader, comments, needsSeparator);
        needsSeparator = true;
    }
}

bool endsWith(std::string_view s, std::string_view suffix)
{
    return (s.size() > suffix.size()) &&
        (s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0);
}

} // Anonymous namespace

namespace Opm {

PropertyTree::PropertyTree()
    : tree_(std::make_unique<boost::property_tree::ptree>())
{
}

PropertyTree::PropertyTree(const PropertyTree& tree)
    : tree_(std::make_unique<boost::property_tree::ptree>(*tree.tree_))
{
}

PropertyTree::PropertyTree(const std::string& file)
    : tree_(std::make_unique<boost::property_tree::ptree>())
{
    if (isTomlFileName(file)) {
        *tree_ = readTomlFile(file);
    }
    else {
        boost::property_tree::read_json(file, *tree_);
    }
}

PropertyTree PropertyTree::fromJsonString(const std::string& json)
{
    auto tree = boost::property_tree::ptree{};
    auto jsonStream = std::istringstream{json};
    boost::property_tree::read_json(jsonStream, tree);
    return PropertyTree(tree);
}

PropertyTree PropertyTree::fromTomlString(const std::string& toml)
{
    return PropertyTree(readTomlString(toml));
}

bool PropertyTree::isTomlFileName(std::string_view fileName)
{
    return endsWith(fileName, ".toml");
}

bool PropertyTree::isConfigFileName(std::string_view fileName)
{
    return isTomlFileName(fileName) || endsWith(fileName, ".json");
}

PropertyTree::PropertyTree(const boost::property_tree::ptree& tree)
    : tree_(std::make_unique<boost::property_tree::ptree>(tree))
{
}

PropertyTree::~PropertyTree() = default;

template<class T>
T PropertyTree::get(const std::string& key) const
{
    return tree_->get<T>(key);
}

template<class T>
T PropertyTree::get(const std::string& key, const T& defValue) const
{
    return tree_->get<T>(key, defValue);
}

template<class T>
void PropertyTree::put(const std::string& key, const T& value)
{
    tree_->put(key,value);
}

void PropertyTree::write_json(std::ostream &os, bool pretty) const
{
    boost::property_tree::write_json(os, *tree_, pretty);
}

void PropertyTree::write_toml(std::ostream& os,
                              const TomlCommentProvider& comments) const
{
    bool needsSeparator = false;
    writeTomlTable(os, *tree_, std::string{}, std::string{}, comments, needsSeparator);
}

PropertyTree
PropertyTree::get_child(const std::string& key) const
{
  auto pt = tree_->get_child(key);

  return PropertyTree(pt);
}

std::optional<PropertyTree>
PropertyTree::get_child_optional(const std::string& key) const
{
  auto pt = tree_->get_child_optional(key);
  if (!pt)
      return std::nullopt;

  return PropertyTree(pt.get());
}

std::vector<std::string> PropertyTree::get_child_keys() const
{
    std::vector<std::string> keys;
    for (auto it = tree_->begin(); it != tree_->end(); ++it) {
        keys.push_back(it->first);
    }
    return keys;
}

template <typename T>
std::optional<std::vector<T>>
PropertyTree::get_child_items_as_vector(const std::string& child) const
{
    auto items = std::optional<std::vector<T>>{};

    auto subTree = this->tree_->get_child_optional(child);
    if (! subTree) {
        return items;
    }

    items.emplace();
    // compatibility with older platforms where code is built as c++-17
#if __cplusplus>= 202002L
    std::ranges::transform(*subTree,
#else
    std::transform(subTree->begin(), subTree->end(),
#endif
                           std::back_inserter(*items),
                           [](const auto& childItem)
                           { return childItem.second.template get_value<T>(); });

    return items;
}

PropertyTree& PropertyTree::operator=(const PropertyTree& tree)
{
  tree_ = std::make_unique<boost::property_tree::ptree>(*tree.tree_);
  return *this;
}

template void PropertyTree::put(const std::string& key, const std::string& value);
template void PropertyTree::put(const std::string& key, const double& value);
template void PropertyTree::put(const std::string& key, const float& value);
template void PropertyTree::put(const std::string& key, const int& value);
template void PropertyTree::put(const std::string& key, const std::size_t& value);
template void PropertyTree::put(const std::string& key, const bool& value);

template std::string PropertyTree::get(const std::string& key) const;
template double PropertyTree::get(const std::string& key) const;
template float PropertyTree::get(const std::string& key) const;
template int PropertyTree::get(const std::string& key) const;
template std::size_t PropertyTree::get(const std::string& key) const;
template bool PropertyTree::get(const std::string& key) const;

template std::string PropertyTree::get(const std::string& key, const std::string& defValue) const;
template double PropertyTree::get(const std::string& key, const double& defValue) const;
template float PropertyTree::get(const std::string& key, const float& defValue) const;
template int PropertyTree::get(const std::string& key, const int& defValue) const;
template std::size_t PropertyTree::get(const std::string& key, const std::size_t& defValue) const;
template bool PropertyTree::get(const std::string& key, const bool& defValue) const;

template std::optional<std::vector<int>>
PropertyTree::get_child_items_as_vector(const std::string& child) const;
template std::optional<std::vector<double>>
PropertyTree::get_child_items_as_vector(const std::string& child) const;

} // namespace Opm
