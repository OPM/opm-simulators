/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_PROPERTYTREE_HEADER_INCLUDED
#define OPM_PROPERTYTREE_HEADER_INCLUDED

#include <iosfwd>
#include <memory>
#include <optional>

namespace boost {
namespace property_tree {
    template<class T1, class T2, class T3> class basic_ptree;
    using ptree = basic_ptree<std::string,std::string,std::less<std::string>>;
}
}

namespace Opm
{

class PropertyTree {
public:
    PropertyTree();
    PropertyTree(const std::string& jsonFile);
    PropertyTree(const PropertyTree& tree);
    ~PropertyTree();

    template<class T>
    void put(const std::string& key, const T& data);

    template<class T>
    T get(const std::string& key) const;

    template<class T>
    T get(const std::string& key, const T& defValue) const;

    PropertyTree get_child(const std::string& key) const;

    std::optional<PropertyTree> get_child_optional(const std::string& key) const;

    PropertyTree& operator=(const PropertyTree& tree);

    void write_json(std::ostream& os, bool pretty) const;

protected:
    PropertyTree(const boost::property_tree::ptree& tree);

    std::unique_ptr<boost::property_tree::ptree> tree_;
};


} // namespace Opm

#endif // OPM_PROPERTYTREE_HEADER_INCLUDED
