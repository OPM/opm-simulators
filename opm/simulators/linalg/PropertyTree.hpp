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

#include <functional>
#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace boost::property_tree {
    template<class T1, class T2, class T3> class basic_ptree;
    using ptree = basic_ptree<std::string, std::string, std::less<std::string>>;
} // namespace boost::property_tree

namespace Opm {

/// Hierarchical collection of key/value pairs
class PropertyTree
{
public:
    /// Default constructor.
    ///
    /// Should typically be populated in put() before use.
    PropertyTree();

    /// Constructor
    ///
    /// Loads a property tree from an external source expected to be a text
    /// file in JSON.
    ///
    /// \param[in] jsonFile Name of file containing external property tree,
    /// linearised into JSON format.
    explicit PropertyTree(const std::string& jsonFile);

    /// Copy constructor.
    ///
    /// \param[in] tree Source object.
    PropertyTree(const PropertyTree& tree);

    /// Destructor.
    ~PropertyTree();

    /// Assignment operator
    ///
    /// \param[in] tree Source object.
    ///
    /// \return \code *this \endcode.
    PropertyTree& operator=(const PropertyTree& tree);

    /// Insert key/value pair into property tree
    ///
    /// \tparam T Value type
    ///
    /// \param[in] key Property key.  Expected to be in hierarchical
    /// notation for subtrees--i.e., using periods ('.') to separate
    /// hierarchy levels.
    ///
    /// \param[in] data Property value corresponding to \p key.
    template<class T>
    void put(const std::string& key, const T& data);

    /// Retrieve property value given hierarchical property key.
    ///
    /// \tparam T Value type
    ///
    /// \param[in] key Property key.  Expected to be in hierarchical
    /// notation for subtrees--i.e., using periods ('.') to separate
    /// hierarchy levels.
    ///
    /// \return Copy of internal property value for \p key.
    template<class T>
    T get(const std::string& key) const;

    /// Retrieve property value given hierarchical property key.
    ///
    /// \tparam T Value type
    ///
    /// \param[in] key Property key.  Expected to be in hierarchical
    /// notation for subtrees--i.e., using periods ('.') to separate
    /// hierarchy levels.
    ///
    /// \param[in] defValue Default value for when \p key is not in the
    /// property tree.
    ///
    /// \return Copy of internal property value for \p key, or a copy of \p
    /// defValue if the \p key is not in the property tree.
    template<class T>
    T get(const std::string& key, const T& defValue) const;

    /// Retrieve copy of sub tree rooted at node.
    ///
    /// Throws an exception if no sub tree exists at given root.
    ///
    /// \param[in] key Property key.  Expected to be in hierarchical
    /// notation for subtrees--i.e., using periods ('.') to separate
    /// hierarchy levels.
    ///
    /// \return Copy of property sub tree rooted at \p key.
    PropertyTree get_child(const std::string& key) const;

    /// Retrieve copy of sub tree rooted at node.
    ///
    /// \param[in] key Property key.  Expected to be in hierarchical
    /// notation for subtrees--i.e., using periods ('.') to separate
    /// hierarchy levels.
    ///
    /// \return Copy of property sub tree rooted at \p key.  Nullopt if no
    /// sub tree exists that is rooted at \p key.
    std::optional<PropertyTree>
    get_child_optional(const std::string& key) const;

    /// Retrieve node items as linearised vector.
    ///
    /// Assumes that the node's child is an array type of homongeneous
    /// elements.
    ///
    /// \tparam T Array element type.
    ///
    /// \param[in] child Property key.  Expected to be in hierarchical
    /// notation for subtrees--i.e., using periods ('.') to separate
    /// hierarchy levels.
    ///
    /// \return Array of property values.  Nullopt if no node named by \p
    /// child exists.
    template <typename T>
    std::optional<std::vector<T>>
    get_child_items_as_vector(const std::string& child) const;

    /// Emit a textual representation of the property tree in JSON form
    ///
    /// \param[in,out] os Output stream.  Typically a stream opened on a
    /// file.
    ///
    /// \param[in] pretty Whether or not to pretty-print the JSON
    /// output--i.e., whether or not to insert new lines and spaces for
    /// human readability.
    void write_json(std::ostream& os, bool pretty) const;

protected:
    /// Converting constructor.
    ///
    /// Forms a property tree object from a Boost ptree.
    ///
    /// \param[in] tree Source object represented as a Boost ptree.
    PropertyTree(const boost::property_tree::ptree& tree);

    /// Internal representation of the property tree.
    std::unique_ptr<boost::property_tree::ptree> tree_;
};

} // namespace Opm

#endif // OPM_PROPERTYTREE_HEADER_INCLUDED
