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

BEGIN_PROPERTIES

//! Set the stencil
SET_PROP(VcfvDiscretization, Stencil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::ctype CoordScalar;

public:
    typedef Opm::VcfvStencil<CoordScalar, GridView> type;
};

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(VcfvDiscretization, DofMapper, typename GET_PROP_TYPE(TypeTag, VertexMapper));

//! The concrete class which manages the spatial discretization
SET_TYPE_PROP(VcfvDiscretization, Discretization, Opm::VcfvDiscretization<TypeTag>);

//! The base class for the output modules (decides whether to write
//! element or vertex based fields)
SET_TYPE_PROP(VcfvDiscretization, DiscBaseOutputModule,
              Opm::VcfvBaseOutputModule<TypeTag>);

//! Calculates the gradient of any quantity given the index of a flux approximation point
SET_TYPE_PROP(VcfvDiscretization, GradientCalculator,
              Opm::P1FeGradientCalculator<TypeTag>);

//! The class to create grid communication handles
SET_TYPE_PROP(VcfvDiscretization, GridCommHandleFactory,
              Opm::VcfvGridCommHandleFactory<TypeTag>);

//! Use two-point gradients by default for the vertex centered finite volume scheme.
SET_BOOL_PROP(VcfvDiscretization, UseP1FiniteElementGradients, false);

#if HAVE_DUNE_FEM
//! Set the DiscreteFunctionSpace
SET_PROP(VcfvDiscretization, DiscreteFunctionSpace)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridPart) GridPart;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::Fem::FunctionSpace<typename GridPart::GridType::ctype,
                                     Scalar,
                                     GridPart::GridType::dimensionworld,
                                     numEq> FunctionSpace;
public:
    // Lagrange discrete function space with unknowns at the cell vertices
    typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, 1 > type;
};
#endif

//! Set the border list creator for vertices
SET_PROP(VcfvDiscretization, BorderListCreator)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef Opm::Linear::VertexBorderListFromGrid<GridView, VertexMapper> type;
};

//! For the vertex centered finite volume method, ghost and overlap elements must _not_
//! be assembled to avoid accounting twice for the fluxes over the process boundary faces
//! of the local process' grid partition
SET_BOOL_PROP(VcfvDiscretization, LinearizeNonLocalElements, false);


END_PROPERTIES

namespace Opm {

/*!
 * \ingroup VcfvDiscretization
 *
 * \brief The base class for the vertex centered finite volume discretization scheme.
 */
template<class TypeTag>
class VcfvDiscretization : public FvBaseDiscretization<TypeTag>
{
    typedef FvBaseDiscretization<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

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
