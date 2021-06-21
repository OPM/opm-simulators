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

#include <config.h>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>

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

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <opm/material/constraintsolvers/NcpFlash.hpp>

namespace Opm {

template class Tabulated1DFunction<double>;
template class UniformTabulated2DFunction<double>;
template class UniformXTabulated2DFunction<double>;

template class DenseAd::Evaluation<double,1>;
template class DenseAd::Evaluation<double,2>;
template class DenseAd::Evaluation<double,3>;
template class DenseAd::Evaluation<double,4>;
template class DenseAd::Evaluation<double,5>;
template class DenseAd::Evaluation<double,6>;
template class DenseAd::Evaluation<double,7>;
template class DenseAd::Evaluation<double,8>;
template class DenseAd::Evaluation<double,9>;
template class DenseAd::Evaluation<double,10>;
template class DenseAd::Evaluation<double,11>;
template class DenseAd::Evaluation<double,12>;

template class DenseAd::Evaluation<double,-1,4u>;
template class DenseAd::Evaluation<double,-1,5u>;
template class DenseAd::Evaluation<double,-1,6u>;
template class DenseAd::Evaluation<double,-1,7u>;
template class DenseAd::Evaluation<double,-1,8u>;
template class DenseAd::Evaluation<double,-1,9u>;
template class DenseAd::Evaluation<double,-1,10u>;

template class BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>;

template class BrineCo2Pvt<double>;
template class Co2GasPvt<double>;
template class ConstantCompressibilityBrinePvt<double>;
template class ConstantCompressibilityOilPvt<double>;
template class ConstantCompressibilityWaterPvt<double>;
template class DeadOilPvt<double>;
template class DryGasPvt<double>;
template class GasPvtThermal<double>;
template class LiveOilPvt<double>;
template class OilPvtThermal<double>;
template class SolventPvt<double>;
template class WaterPvtThermal<double>;
template class WetGasPvt<double>;

template class GasPvtMultiplexer<double, false>;
template class GasPvtMultiplexer<double, true>;
template class OilPvtMultiplexer<double, false>;
template class OilPvtMultiplexer<double, true>;
template class WaterPvtMultiplexer<double, false>;
template class WaterPvtMultiplexer<double, true>;

template class EclEpsScalingPoints<double>;

template class NcpFlash<double, BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>>;

}
