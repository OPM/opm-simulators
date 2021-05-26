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
 * \copydoc Opm::VcfvDiscretization
 */
#ifndef EWOMS_VCFV_DISCRETIZATION_HH
#define EWOMS_VCFV_DISCRETIZATION_HH

#include <opm/material/densead/Math.hpp>

#include "vcfvproperties.hh"
#include "vcfvstencil.hh"
#include "p1fegradientcalculator.hh"
#include "vcfvgridcommhandlefactory.hh"
#include "vcfvbaseoutputmodule.hh"

#include <opm/simulators/linalg/vertexborderlistfromgrid.hh>
#include <opm/models/discretization/common/fvbasediscretization.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#endif

namespace Opm {
template <class TypeTag>
class VcfvDiscretization;

} // namespace Opm

namespace Opm::Properties {

//! Set the stencil
template<class TypeTag>
struct Stencil<TypeTag, TTag::VcfvDiscretization>
{
private:
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using CoordScalar = typename GridView::ctype;

public:
    using type = VcfvStencil<CoordScalar, GridView>;
};

//! Mapper for the degrees of freedoms.
template<class TypeTag>
struct DofMapper<TypeTag, TTag::VcfvDiscretization> { using type = GetPropType<TypeTag, Properties::VertexMapper>; };

//! The concrete class which manages the spatial discretization
template<class TypeTag>
struct Discretization<TypeTag, TTag::VcfvDiscretization> { using type = VcfvDiscretization<TypeTag>; };

//! The base class for the output modules (decides whether to write
//! element or vertex based fields)
template<class TypeTag>
struct DiscBaseOutputModule<TypeTag, TTag::VcfvDiscretization>
{ using type = VcfvBaseOutputModule<TypeTag>; };

//! Calculates the gradient of any quantity given the index of a flux approximation point
template<class TypeTag>
struct GradientCalculator<TypeTag, TTag::VcfvDiscretization>
{ using type = P1FeGradientCalculator<TypeTag>; };

//! The class to create grid communication handles
template<class TypeTag>
struct GridCommHandleFactory<TypeTag, TTag::VcfvDiscretization>
{ using type = VcfvGridCommHandleFactory<TypeTag>; };

//! Use two-point gradients by default for the vertex centered finite volume scheme.
template<class TypeTag>
struct UseP1FiniteElementGradients<TypeTag, TTag::VcfvDiscretization> { static constexpr bool value = false; };

#if HAVE_DUNE_FEM
//! Set the DiscreteFunctionSpace
template<class TypeTag>
struct DiscreteFunctionSpace<TypeTag, TTag::VcfvDiscretization>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>  ;
    using GridPart = GetPropType<TypeTag, Properties::GridPart>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using FunctionSpace = Dune::Fem::FunctionSpace<typename GridPart::GridType::ctype,
                                                   Scalar,
                                                   GridPart::GridType::dimensionworld,
                                                   numEq>;
public:
    // Lagrange discrete function space with unknowns at the cell vertices
    using type = Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, 1 >;
};
#endif

//! Set the border list creator for vertices
template<class TypeTag>
struct BorderListCreator<TypeTag, TTag::VcfvDiscretization>
{ private:
    using VertexMapper = GetPropType<TypeTag, Properties::VertexMapper>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
public:
    using type = Linear::VertexBorderListFromGrid<GridView, VertexMapper>;
};

//! For the vertex centered finite volume method, ghost and overlap elements must _not_
//! be assembled to avoid accounting twice for the fluxes over the process boundary faces
//! of the local process' grid partition
template<class TypeTag>
struct LinearizeNonLocalElements<TypeTag, TTag::VcfvDiscretization> { static constexpr bool value = false; };


} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup VcfvDiscretization
 *
 * \brief The base class for the vertex centered finite volume discretization scheme.
 */
template<class TypeTag>
class VcfvDiscretization : public FvBaseDiscretization<TypeTag>
{
    using ParentType = FvBaseDiscretization<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Model>;
    using DofMapper = GetPropType<TypeTag, Properties::DofMapper>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    enum { dim = GridView::dimension };

public:
    VcfvDiscretization(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Returns a string of discretization's human-readable name
     */
    static std::string discretizationName()
    { return "vcfv"; }

    /*!
     * \brief Returns the number of global degrees of freedom (DOFs) due to the grid
     */
    size_t numGridDof() const
    { return static_cast<size_t>(this->gridView_.size(/*codim=*/dim)); }

    /*!
     * \brief Mapper to convert the Dune entities of the
     *        discretization's degrees of freedoms are to indices.
     */
    const DofMapper& dofMapper() const
    { return this->vertexMapper(); }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter& res)
    { res.template serializeEntities</*codim=*/dim>(asImp_(), this->gridView_); }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        res.template deserializeEntities</*codim=*/dim>(asImp_(), this->gridView_);
        this->solution(/*timeIdx=*/1) = this->solution(/*timeIdx=*/0);
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
} // namespace Opm

#endif
