/*
  Copyright 2024, SINTEF Digital

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

#ifndef OPM_EMPTY_MODEL_HPP
#define OPM_EMPTY_MODEL_HPP

// this is an empty model that having a lot of empty interfaces.
// it is use for the development when some facility class are not ready

#include <opm/output/data/Aquifer.hpp>

#include <opm/models/discretization/common/baseauxiliarymodule.hh>

namespace Opm {

template<typename TypeTag>
class EmptyModel : public BaseAuxiliaryModule<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    explicit EmptyModel(Simulator& /*simulator*/)
    {
    }

    void init(){}
    template<class Something>
    void init(Something /*A*/){}
    void prepareTracerBatches(){};
    using NeighborSet = std::set<unsigned>;
    void linearize(SparseMatrixAdapter& /*matrix*/, GlobalEqVector& /*residual*/) override {}
    unsigned numDofs() const override { return 0; }
    void addNeighbors(std::vector<NeighborSet>& /*neighbors*/) const override {}
    void initialSolutionApplied(){};
    template <class Restarter>
    void serialize(Restarter& /*res*/){};

    template <class Restarter>
    void deserialize(Restarter& /*res*/){};

    void beginEpisode(){};
    void beginTimeStep(){};
    void beginIteration(){};
    // add the water rate due to aquifers to the source term.
    template<class RateVector, class Context>
    void addToSource(RateVector& /*rates*/, const Context& /*context*/,
                     unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const {}
    template<class RateVector>
    void addToSource(RateVector& /*rates*/, unsigned /*globalSpaceIdx*/,
                     unsigned /*timeIdx*/) const {}
    void endIteration()const{};
    void endTimeStep(){};
    void endEpisode(){};
    void applyInitial() override {}
    auto aquiferData() const {
        return data::Aquifers{};
    }
};

} // end of namespace Opm

#endif // OPM_EMPTY_MODEL_HPP
