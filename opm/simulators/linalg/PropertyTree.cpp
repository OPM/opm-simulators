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

#include <boost/property_tree/json_parser.hpp>

namespace Opm
{

PropertyTree::PropertyTree()
    : tree_(std::make_unique<boost::property_tree::ptree>())
{
}

PropertyTree::PropertyTree(const PropertyTree& tree)
    : tree_(std::make_unique<boost::property_tree::ptree>(*tree.tree_))
{
}

PropertyTree::PropertyTree(const std::string& jsonFile)
    : tree_(std::make_unique<boost::property_tree::ptree>())
{
    boost::property_tree::read_json(jsonFile, *tree_);
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

PropertyTree& PropertyTree::operator=(const PropertyTree& tree)
{
  tree_ = std::make_unique<boost::property_tree::ptree>(*tree.tree_);
  return *this;
}

template std::string PropertyTree::get<std::string>(const std::string& key) const;
template std::string PropertyTree::get<std::string>(const std::string& key, const std::string& defValue) const;
template double PropertyTree::get<double>(const std::string& key) const;
template double PropertyTree::get<double>(const std::string& key, const double& defValue) const;
template int PropertyTree::get<int>(const std::string& key) const;
template int PropertyTree::get<int>(const std::string& key, const int& defValue) const;
template size_t PropertyTree::get<size_t>(const std::string& key) const;
template size_t PropertyTree::get<size_t>(const std::string& key, const size_t& defValue) const;
template bool PropertyTree::get<bool>(const std::string& key) const;
template bool PropertyTree::get<bool>(const std::string& key, const bool& defValue) const;

template void PropertyTree::put<std::string>(const std::string& key, const std::string& value);
template void PropertyTree::put<double>(const std::string& key, const double& value);
template void PropertyTree::put<int>(const std::string& key, const int& value);


} // namespace Opm
