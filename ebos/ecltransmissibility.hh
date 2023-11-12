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
 * \copydoc Opm::EclTransmissibility
 */
#ifndef EWOMS_ECL_TRANSMISSIBILITY_HH
#define EWOMS_ECL_TRANSMISSIBILITY_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/LookUpData.hh>


#include <array>
#include <functional>
#include <map>
#include <cstdint>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace Opm {

class KeywordLocation;
class EclipseState;
struct NNCdata;
class TransMult;


template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
class EclTransmissibility {
    // Grid and world dimension
    enum { dimWorld = GridView::dimensionworld };
public:

    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

    EclTransmissibility(const EclipseState& eclState,
                        const GridView& gridView,
                        const CartesianIndexMapper& cartMapper,
                        const Grid& grid,
                        std::function<std::array<double,dimWorld>(int)> centroids,
                        bool enableEnergy,
                        bool enableDiffusivity);

    /*!
     * \brief Return the permeability for an element.
     */
    const DimMatrix& permeability(unsigned elemIdx) const
    { return permeability_[elemIdx]; }

    /*!
     * \brief Return the transmissibility for the intersection between two elements.
     */
    Scalar transmissibility(unsigned elemIdx1, unsigned elemIdx2) const;

    /*!
     * \brief Return the transmissibility for a given boundary segment.
     */
    Scalar transmissibilityBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const;

    /*!
     * \brief Return the thermal "half transmissibility" for the intersection between two
     *        elements.
     *
     * The "half transmissibility" features all sub-expressions of the "thermal
     * transmissibility" which can be precomputed, i.e. they are not dependent on the
     * current solution:
     *
     * H_t = A * (n*d)/(d*d);
     *
     * where A is the area of the intersection between the inside and outside elements, n
     * is the outer unit normal, and d is the distance between the center of the inside
     * cell and the center of the intersection.
     */
    Scalar thermalHalfTrans(unsigned insideElemIdx, unsigned outsideElemIdx) const;

    Scalar thermalHalfTransBoundary(unsigned insideElemIdx, unsigned boundaryFaceIdx) const;

    /*!
     * \brief Return the diffusivity for the intersection between two elements.
     */
    Scalar diffusivity(unsigned elemIdx1, unsigned elemIdx2) const;

    /*!
     * \brief Actually compute the transmissibility over a face as a pre-compute step.
     *
     * This code actually uses the direction specific "centroids" of
     * each element. These "centroids" are _not_ the identical
     * barycenter of the element, but the middle of the centers of the
     * faces of the logical Cartesian cells, i.e., the centers of the
     * faces of the reference elements. We do things this way because
     * the barycenter of the element can be located outside of the
     * element for sufficiently "ugly" (i.e., thin and "non-flat")
     * elements which in turn leads to quite wrong
     * permeabilities. This approach is probably not always correct
     * either but at least it seems to be much better.
     */
    void finishInit(const std::function<unsigned int(unsigned int)>& map = {})
    {
        this->update(true, map, /*applyNncMultRegT = */ true);
    }

    /*!
     * \brief Compute all transmissibilities
     *
     * \param[in] global Whether or not to call \c update() on all
     *   processes.  Also, this updates the "thermal half
     *   transmissibilities" if energy is enabled.
     *
     * \param[in] map Undocumented.
     *
     * \param[in] applyNncMultRegT Whether or not to apply regional
     *   multipliers to explicit NNCs.  Explicit NNCs are those entered
     *   directly in the input data, e.g., through the NNC/EDITNNC/EDITNNCR
     *   keywords, or the result of generating connections to or within
     *   numerical aquifers.  Default value: \c false, meaning do not apply
     *   regional multipliers to explicit NNCs.
     */
    void update(bool global, const std::function<unsigned int(unsigned int)>& map = {}, bool applyNncMultRegT = false);

protected:
    void updateFromEclState_(bool global);

    void removeSmallNonCartesianTransmissibilities_();

    /// \brief Apply the Multipliers for the case PINCH(4)==TOPBOT
    ///
    /// \param pinchTop Whether PINCH(5) is TOP, otherwise ALL is assumed.
    void applyAllZMultipliers_(Scalar& trans,
                               unsigned insideFaceIdx,
                               unsigned outsideFaceIdx,
                               unsigned insideCartElemIdx,
                               unsigned outsideCartElemIdx,
                               const TransMult& transMult,
                               const std::array<int, dimWorld>& cartDims,
                               bool pinchTop);

    /// \brief Creates TRANS{XYZ} arrays for modification by FieldProps data
    ///
    /// \param is_tran Whether TRAN{XYZ} will be modified. If entry is false the array will be empty
    /// \returns an array of vector (TRANX, TRANY, TRANZ}
    std::array<std::vector<double>,3>
    createTransmissibilityArrays_(const std::array<bool,3>& is_tran);

    /// \brief overwrites calculated transmissibilities
    ///
    /// \param is_tran Whether TRAN{XYZ} have been modified.
    /// \param trans Arrays with modified transmissibilities TRAN{XYZ}
    void resetTransmissibilityFromArrays_(const std::array<bool,3>& is_tran,
                                          const std::array<std::vector<double>,3>& trans);

    template <class Intersection>
    void computeFaceProperties(const Intersection& intersection,
                               const int,
                               const int,
                               const int,
                               const int,
                               DimVector& faceCenterInside,
                               DimVector& faceCenterOutside,
                               DimVector& faceAreaNormal,
                               /*isCpGrid=*/std::false_type) const;

    template <class Intersection>
    void computeFaceProperties(const Intersection& intersection,
                               const int insideElemIdx,
                               const int insideFaceIdx,
                               const int outsideElemIdx,
                               const int outsideFaceIdx,
                               DimVector& faceCenterInside,
                               DimVector& faceCenterOutside,
                               DimVector& faceAreaNormal,
                               /*isCpGrid=*/std::true_type) const;

    /*
     * \brief Applies additional transmissibilities specified via NNC keyword.
     *
     * Applies only those NNC that are actually represented by the grid. These may
     * NNCs due to faults or NNCs that are actually neighbours. In both cases that
     * specified transmissibilities (scaled by EDITNNC) will be added to the already
     * existing models.
     *
     * \param cartesianToCompressed Vector containing the compressed index (or -1 for inactive
     *                              cells) as the element at the cartesian index.
     * \return Nothing.
     */
    void applyNncToGridTrans_(const std::unordered_map<std::size_t,int>& cartesianToCompressed);

    /// \brief Multiplies the grid transmissibilities according to EDITNNC.
    void applyEditNncToGridTrans_(const std::unordered_map<std::size_t,int>& globalToLocal);

    /// \brief Resets the grid transmissibilities according to EDITNNCR.
    void applyEditNncrToGridTrans_(const std::unordered_map<std::size_t,int>& globalToLocal);

    void applyNncMultreg_(const std::unordered_map<std::size_t,int>& globalToLocal);

    void applyEditNncToGridTransHelper_(const std::unordered_map<std::size_t,int>& globalToLocal,
                                        const std::string& keyword, const std::vector<NNCdata>& nncs,
                                        const std::function<KeywordLocation(const NNCdata&)>& getLocation,
                                        const std::function<void(double&, const double&)>& apply);

    void extractPermeability_();

    void extractPermeability_(const std::function<unsigned int(unsigned int)>& map);

    void extractPorosity_();

    void computeHalfTrans_(Scalar& halfTrans,
                           const DimVector& areaNormal,
                           int faceIdx, // in the reference element that contains the intersection
                           const DimVector& distance,
                           const DimMatrix& perm) const;

    void computeHalfDiffusivity_(Scalar& halfDiff,
                                 const DimVector& areaNormal,
                                 const DimVector& distance,
                                 const Scalar& poro) const;

    DimVector distanceVector_(const DimVector& center,
                              int faceIdx, // in the reference element that contains the intersection
                              unsigned elemIdx,
                              const std::array<std::vector<DimVector>, dimWorld>& axisCentroids) const;

    void applyMultipliers_(Scalar& trans,
                           unsigned faceIdx,
                           unsigned cartElemIdx,
                           const TransMult& transMult) const;

    void applyNtg_(Scalar& trans,
                   unsigned faceIdx,
                   unsigned elemIdx,
                   const std::vector<double>& ntg) const;

    std::vector<DimMatrix> permeability_;
    std::vector<Scalar> porosity_;
    std::unordered_map<std::uint64_t, Scalar> trans_;
    const EclipseState& eclState_;
    const GridView& gridView_;
    const CartesianIndexMapper& cartMapper_;
    const Grid& grid_;
    std::function<std::array<double,dimWorld>(int)> centroids_;
    Scalar transmissibilityThreshold_;
    std::map<std::pair<unsigned, unsigned>, Scalar> transBoundary_;
    std::map<std::pair<unsigned, unsigned>, Scalar> thermalHalfTransBoundary_;
    bool enableEnergy_;
    bool enableDiffusivity_;
    std::unordered_map<std::uint64_t, Scalar> thermalHalfTrans_; //NB this is based on direction map size is ca 2*trans_ (diffusivity_)
    std::unordered_map<std::uint64_t, Scalar> diffusivity_;

    const LookUpData<Grid,GridView> lookUpData_;
    const LookUpCartesianData<Grid,GridView> lookUpCartesianData_;
};

} // namespace Opm

#endif
