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
 * \brief Defines a type tags and some fundamental properties all models.
 */
#ifndef EWOMS_BASIC_PROPERTIES_HH
#define EWOMS_BASIC_PROPERTIES_HH

#include <dune/common/parametertree.hh>

#include <opm/models/utils/basicparameters.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#endif

namespace Opm {

template<class TypeTag> class DgfVanguard;

}

namespace Opm::Properties {

///////////////////////////////////
// Type tag definitions:
//
// NumericModel
// |
// +-> ImplicitModel
///////////////////////////////////

// Create new type tags
namespace TTag {

//! Type tag for all models.
struct NumericModel {};

//! Type tag for all fully coupled models.
struct ImplicitModel { using InheritsFrom = std::tuple<NumericModel>; };

} // end namespace TTag

///////////////////////////////////
// Property names which are always available:
//
// Scalar
///////////////////////////////////

//! Property to specify the type of scalar values.
template<class TypeTag, class MyTypeTag>
struct Scalar { using type = UndefinedProperty; };

//! Number of equations in the system of PDEs
template<class TypeTag, class MyTypeTag>
struct NumEq { using type = UndefinedProperty; };

//! Property which provides a Dune::ParameterTree.
template<class TypeTag, class MyTypeTag>
struct ParameterTree { using type = UndefinedProperty; };

//! The type of the model
template<class TypeTag, class MyTypeTag>
struct Model { using type = UndefinedProperty; };

//! Property which defines the group that is queried for parameters by default
template<class TypeTag, class MyTypeTag>
struct ModelParameterGroup { using type = UndefinedProperty; };

//! Property which provides a Vanguard (manages grids)
template<class TypeTag, class MyTypeTag>
struct Vanguard { using type = UndefinedProperty; };

//! The type of the DUNE grid
template<class TypeTag, class MyTypeTag>
struct Grid { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct GridView { using type = UndefinedProperty; };

#if HAVE_DUNE_FEM
template<class TypeTag, class MyTypeTag>
struct GridPart { using type = UndefinedProperty; };
#endif

//! level of the grid view
template<class TypeTag, class MyTypeTag>
struct GridViewLevel { using type = UndefinedProperty; };

//! Manages the simulation time
template<class TypeTag, class MyTypeTag>
struct Simulator { using type = UndefinedProperty; };

/*!
 * \brief The class which marks the border indices associated with the
 *        degrees of freedom on a process boundary.
 *
 * This is required for the algebraic overlap stuff.
 */
template<class TypeTag, class MyTypeTag>
struct BorderListCreator { using type = UndefinedProperty; };

///////////////////////////////////
// Values for the properties
///////////////////////////////////

//! Set the default type of scalar values to double
template<class TypeTag>
struct Scalar<TypeTag, TTag::NumericModel> { using type = double; };

//! Set the ParameterTree property
template<class TypeTag>
struct ParameterTree<TypeTag, TTag::NumericModel>
{
    using type = Dune::ParameterTree;

    static Dune::ParameterTree& tree()
    {
        static Dune::ParameterTree obj_;
        return obj_;
    }
};

//! use the global group as default for the model's parameter group
template<class TypeTag>
struct ModelParameterGroup<TypeTag, TTag::NumericModel>
{ static constexpr auto value = ""; };

#if HAVE_DUNE_FEM
template<class TypeTag>
struct GridPart<TypeTag, TTag::NumericModel>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using type = Dune::Fem::AdaptiveLeafGridPart<Grid>;
};

template<class TypeTag>
struct GridView<TypeTag, TTag::NumericModel>
{ using type = typename GetPropType<TypeTag, Properties::GridPart>::GridViewType; };
#else
//! Use the leaf grid view by default.
//!
//! Except for spatial refinement, there is rarly a reason to use
//! anything else...
template<class TypeTag>
struct GridView<TypeTag, TTag::NumericModel>
{ using type = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView; };
#endif

template<class TypeTag>
struct Vanguard<TypeTag, TTag::NumericModel>
{ using type = Opm::DgfVanguard<TypeTag>; };

} // namespace Opm::Properties

#endif
