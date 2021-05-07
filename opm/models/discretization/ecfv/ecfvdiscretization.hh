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
 * \copydoc Opm::EcfvDiscretization
 */
#ifndef EWOMS_ECFV_DISCRETIZATION_HH
#define EWOMS_ECFV_DISCRETIZATION_HH

#include <opm/material/densead/Math.hpp>

#include "ecfvproperties.hh"
#include "ecfvstencil.hh"
#include "ecfvgridcommhandlefactory.hh"
#include "ecfvbaseoutputmodule.hh"

#include <opm/simulators/linalg/elementborderlistfromgrid.hh>
#include <opm/models/discretization/common/fvbasediscretization.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/finitevolume.hh>
#endif

namespace Opm {
template <class TypeTag>
class EcfvDiscretization;
}

namespace Opm::Properties {

//! Set the stencil
template<class TypeTag>
struct Stencil<TypeTag, TTag::EcfvDiscretization>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

public:
    using type = EcfvStencil<Scalar, GridView>;
};

//! Mapper for the degrees of freedoms.
template<class TypeTag>
struct DofMapper<TypeTag, TTag::EcfvDiscretization> { using type = GetPropType<TypeTag, Properties::ElementMapper>; };

//! The concrete class which manages the spatial discretization
template<class TypeTag>
struct Discretization<TypeTag, TTag::EcfvDiscretization> { using type = EcfvDiscretization<TypeTag>; };

//! The base class for the output modules (decides whether to write
//! element or vertex based fields)
template<class TypeTag>
struct DiscBaseOutputModule<TypeTag, TTag::EcfvDiscretization>
{ using type = EcfvBaseOutputModule<TypeTag>; };

//! The class to create grid communication handles
template<class TypeTag>
struct GridCommHandleFactory<TypeTag, TTag::EcfvDiscretization>
{ using type = EcfvGridCommHandleFactory<TypeTag>; };

#if HAVE_DUNE_FEM
//! Set the DiscreteFunctionSpace
template<class TypeTag>
struct DiscreteFunctionSpace<TypeTag, TTag::EcfvDiscretization>
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
    using type = Dune::Fem::FiniteVolumeSpace< FunctionSpace, GridPart, 0 >;
};
#endif

//! Set the border list creator for to the one of an element based
//! method
template<class TypeTag>
struct BorderListCreator<TypeTag, TTag::EcfvDiscretization>
{ private:
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
public:
    using type = Linear::ElementBorderListFromGrid<GridView, ElementMapper>;
};

//! For the element centered finite volume method, ghost and overlap elements must be
//! assembled to calculate the fluxes over the process boundary faces of the local
//! process' grid partition
template<class TypeTag>
struct LinearizeNonLocalElements<TypeTag, TTag::EcfvDiscretization> { static constexpr bool value = true; };

//! locking is not required for the element centered finite volume method because race
//! conditions cannot occur since each matrix/vector entry is written exactly once
template<class TypeTag>
struct UseLinearizationLock<TypeTag, TTag::EcfvDiscretization> { static constexpr bool value = false; };

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup EcfvDiscretization
 *
 * \brief The base class for the element-centered finite-volume discretization scheme.
 */
template<class TypeTag>
class EcfvDiscretization : public FvBaseDiscretization<TypeTag>
{
    using ParentType = FvBaseDiscretization<TypeTag>;

    using Implementation = GetPropType<TypeTag, Properties::Model>;
    using DofMapper = GetPropType<TypeTag, Properties::DofMapper>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

public:
    EcfvDiscretization(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Returns a string of discretization's human-readable name
     */
    static std::string discretizationName()
    { return "ecfv"; }

    /*!
     * \brief Returns the number of global degrees of freedom (DOFs) due to the grid
     */
    size_t numGridDof() const
    { return static_cast<size_t>(this->gridView_.size(/*codim=*/0)); }

    /*!
     * \brief Mapper to convert the Dune entities of the
     *        discretization's degrees of freedoms are to indices.
     */
    const DofMapper& dofMapper() const
    { return this->elementMapper(); }

    /*!
     * \brief Syncronize the values of the primary variables on the
     *        degrees of freedom that overlap with the neighboring
     *        processes.
     *
     * For the Element Centered Finite Volume discretization, this
     * method retrieves the primary variables corresponding to
     * overlap/ghost elements from their respective master process.
     */
    void syncOverlap()
    {
        // syncronize the solution on the ghost and overlap elements
        using GhostSyncHandle = GridCommHandleGhostSync<PrimaryVariables,
                                                        SolutionVector,
                                                        DofMapper,
                                                        /*commCodim=*/0>;

        auto ghostSync = GhostSyncHandle(this->solution(/*timeIdx=*/0),
                                         asImp_().dofMapper());
        this->gridView().communicate(ghostSync,
                                     Dune::InteriorBorder_All_Interface,
                                     Dune::ForwardCommunication);
    }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter& res)
    { res.template serializeEntities</*codim=*/0>(asImp_(), this->gridView_); }

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
        res.template deserializeEntities</*codim=*/0>(asImp_(), this->gridView_);
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
