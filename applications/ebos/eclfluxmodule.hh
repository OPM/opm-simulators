/*
  Copyright (C) 2014 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \brief This file contains the flux module which is used for ECL problems
 *        two-point flux approximation
 *
 * This is used by the ECL blackoil simulator
 */
#ifndef EWOMS_ECL_FLUX_MODULE_HH
#define EWOMS_ECL_FLUX_MODULE_HH

#include <ewoms/disc/common/fvbaseproperties.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(MaterialLaw);
}}

namespace Ewoms {
template <class TypeTag>
class EclTransIntensiveQuantities;

template <class TypeTag>
class EclTransExtensiveQuantities;

template <class TypeTag>
class EclTransBaseProblem;

/*!
 * \ingroup EclTransmissibility
 * \brief Specifies a velocity module which uses the transmissibilities.
 */
template <class TypeTag>
struct EclTransVelocityModule
{
    typedef EclTransIntensiveQuantities<TypeTag> VelocityIntensiveQuantities;
    typedef EclTransExtensiveQuantities<TypeTag> VelocityExtensiveQuantities;
    typedef EclTransBaseProblem<TypeTag> VelocityBaseProblem;

    /*!
     * \brief Register all run-time parameters for the velocity module.
     */
    static void registerParameters()
    { }
};

/*!
 * \ingroup EclTransmissibility
 * \brief Provides the defaults for the parameters required by the
 *        transmissibility based volume flux calculation.
 */
template <class TypeTag>
class EclTransBaseProblem
{ };

/*!
 * \ingroup EclTransmissibility
 * \brief Provides the intensive quantities for the Darcy velocity module
 */
template <class TypeTag>
class EclTransIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
protected:
    void update_(const ElementContext &elemCtx, int dofIdx, int timeIdx)
    { }
};

/*!
 * \ingroup EclTransmissibility
 * \brief Provides the ECL "velocity module"
 */
template <class TypeTag>
class EclTransExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \brief Returns transmissibility for a given sub-control volume face.
     */
    Scalar transmissibility() const
    { return trans_; }

    /*!
     * \brief Return the intrinsic permeability tensor at a face [m^2]
     */
    const DimMatrix& intrinsicPermeability() const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The ECL transmissibility module does not provide an explicit intrinsic permeability");
    }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase at the
     *        face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const DimVector& potentialGrad(int phaseIdx) const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The ECL transmissibility module does not provide explicit potential gradients");
    }

    /*!
     * \brief Return the filter velocity of a fluid phase at the
     *        face's integration point [m/s]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const DimVector& filterVelocity(int phaseIdx) const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The ECL transmissibility module does not provide explicit filter velocities");
    }

    /*!
     * \brief Return the volume flux of a fluid phase at the face's integration point
     *        \f$[m^3/s / m^2]\f$
     *
     * This is the fluid volume of a phase per second and per square meter of face
     * area.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar volumeFlux(int phaseIdx) const
    { return - pressureDifferential_[phaseIdx]*mobility_[phaseIdx] * trans_/faceArea_; }

protected:
    /*!
     * \brief Returns the local index of the degree of freedom in which is
     *        in upstream direction.
     *
     * i.e., the DOF which exhibits a higher effective pressure for
     * the given phase.
     */
    int upstreamIndex_(int phaseIdx) const
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return (pressureDifferential_[phaseIdx] >= 0)?exteriorDofIdx_:interiorDofIdx_;
    }

    /*!
     * \brief Returns the local index of the degree of freedom in which is
     *        in downstream direction.
     *
     * i.e., the DOF which exhibits a lower effective pressure for the
     * given phase.
     */
    int downstreamIndex_(int phaseIdx) const
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return (pressureDifferential_[phaseIdx] >= 0)?interiorDofIdx_:exteriorDofIdx_;
    }

    /*!
     * \brief Update the required gradients for interior faces
     */
    void calculateGradients_(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        Valgrind::SetUndefined(*this);

        const auto& problem = elemCtx.problem();
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        interiorDofIdx_ = scvf.interiorIndex();
        exteriorDofIdx_ = scvf.exteriorIndex();
        assert(interiorDofIdx_ != exteriorDofIdx_);

        trans_ = problem.transmissibility(stencil.globalSpaceIndex(interiorDofIdx_),
                                          stencil.globalSpaceIndex(exteriorDofIdx_));
        faceArea_ = scvf.area();

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        const auto &intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx_, timeIdx);
        const auto &intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx_, timeIdx);

        Scalar zIn = elemCtx.pos(interiorDofIdx_, timeIdx)[dimWorld - 1];
        Scalar zEx = elemCtx.pos(exteriorDofIdx_, timeIdx)[dimWorld - 1];
        Scalar zFace = scvf.integrationPos()[dimWorld - 1];

        Scalar distZIn = zIn - zFace;
        Scalar distZEx = zEx - zFace;

        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            // calculate the hydrostatic pressures at the face's integration point
            Scalar rhoIn = intQuantsIn.fluidState().density(phaseIdx);
            Scalar rhoEx = intQuantsEx.fluidState().density(phaseIdx);

            Scalar pressureInterior = intQuantsIn.fluidState().pressure(phaseIdx);
            Scalar pressureExterior = intQuantsEx.fluidState().pressure(phaseIdx);

            pressureInterior += - rhoIn*(g*distZIn);
            pressureExterior += - rhoEx*(g*distZEx);

            pressureDifferential_[phaseIdx] = pressureExterior - pressureInterior;

            const auto& up = elemCtx.intensiveQuantities(upstreamIndex_(phaseIdx), timeIdx);
            mobility_[phaseIdx] = up.mobility(phaseIdx);
        }
    }

    /*!
     * \brief Update the velocities for all fluid phases on the interior faces of the context
     */
    void calculateVelocities_(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    { }

    // the local indices of the interior and exterior degrees of freedom
    int interiorDofIdx_;
    int exteriorDofIdx_;

    // transmissibility [m^3 s]
    Scalar trans_;

    // the area of the face between the DOFs [m^2]
    Scalar faceArea_;

    // the mobility of all phases [1 / (Pa s)]
    Scalar mobility_[numPhases];

    // the difference in effective pressure between the two degrees of
    // freedom [Pa]
    Scalar pressureDifferential_[numPhases];
};

} // namespace Ewoms

#endif
