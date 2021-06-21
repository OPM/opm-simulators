/*
  Copyright 2021 Equinor AS.

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
#ifndef OPM_MATERIAL_TYPES_HPP
#define OPM_MATERIAL_TYPES_HPP

#include <opm/material/common/IntervalTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <opm/material/densead/Evaluation1.hpp>
#include <opm/material/densead/Evaluation2.hpp>
#include <opm/material/densead/Evaluation3.hpp>
#include <opm/material/densead/Evaluation4.hpp>
#include <opm/material/densead/Evaluation5.hpp>
#include <opm/material/densead/Evaluation6.hpp>
#include <opm/material/densead/Evaluation7.hpp>
#include <opm/material/densead/Evaluation8.hpp>
#include <opm/material/densead/Evaluation9.hpp>
#include <opm/material/densead/Evaluation10.hpp>
#include <opm/material/densead/Evaluation11.hpp>
#include <opm/material/densead/Evaluation12.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/blackoilpvt/BrineCo2Pvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/Co2GasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityBrinePvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DeadOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/GasPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/GasPvtThermal.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/OilPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/OilPvtThermal.hpp>
#include <opm/material/fluidsystems/blackoilpvt/SolventPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtThermal.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>

#include <opm/material/constraintsolvers/NcpFlash.hpp>

namespace Opm {

extern template class IntervalTabulated2DFunction<double>;
extern template class Tabulated1DFunction<double>;
extern template class UniformTabulated2DFunction<double>;
extern template class UniformXTabulated2DFunction<double>;

extern template class DenseAd::Evaluation<double,1>;
extern template class DenseAd::Evaluation<double,2>;
extern template class DenseAd::Evaluation<double,3>;
extern template class DenseAd::Evaluation<double,4>;
extern template class DenseAd::Evaluation<double,5>;
extern template class DenseAd::Evaluation<double,6>;
extern template class DenseAd::Evaluation<double,7>;
extern template class DenseAd::Evaluation<double,8>;
extern template class DenseAd::Evaluation<double,9>;
extern template class DenseAd::Evaluation<double,10>;
extern template class DenseAd::Evaluation<double,11>;
extern template class DenseAd::Evaluation<double,12>;

extern template class DenseAd::Evaluation<double,-1,4u>;
extern template class DenseAd::Evaluation<double,-1,5u>;
extern template class DenseAd::Evaluation<double,-1,6u>;
extern template class DenseAd::Evaluation<double,-1,7u>;
extern template class DenseAd::Evaluation<double,-1,8u>;
extern template class DenseAd::Evaluation<double,-1,9u>;
extern template class DenseAd::Evaluation<double,-1,10u>;

extern template class BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>;

extern template class BrineCo2Pvt<double>;
extern template class Co2GasPvt<double>;
extern template class ConstantCompressibilityBrinePvt<double>;
extern template class ConstantCompressibilityOilPvt<double>;
extern template class ConstantCompressibilityWaterPvt<double>;
extern template class DeadOilPvt<double>;
extern template class DryGasPvt<double>;
extern template class GasPvtThermal<double>;
extern template class LiveOilPvt<double>;
extern template class OilPvtThermal<double>;
extern template class SolventPvt<double>;
extern template class WaterPvtThermal<double>;
extern template class WetGasPvt<double>;

extern template class GasPvtMultiplexer<double, false>;
extern template class GasPvtMultiplexer<double, true>;
extern template class OilPvtMultiplexer<double, false>;
extern template class OilPvtMultiplexer<double, true>;
extern template class WaterPvtMultiplexer<double, false>;
extern template class WaterPvtMultiplexer<double, true>;

extern template class EclEpsScalingPoints<double>;

extern template class NcpFlash<double, BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>>;

}

#endif
