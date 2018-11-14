/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_AUTODIFF_VFPINJPROPERTIESLEGACY_HPP_
#define OPM_AUTODIFF_VFPINJPROPERTIESLEGACY_HPP_

#include <opm/autodiff/VFPHelpersLegacy.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>

#include <vector>
#include <map>



namespace Opm {

template <class Scalar>
class AutoDiffBlock;

class VFPInjPropertiesLegacy : public VFPInjProperties {
public:
    typedef AutoDiffBlock<double> ADB;

    /**
     * Empty constructor
     */
    VFPInjPropertiesLegacy() {}

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_table A *single* VFPINJ table
     */
    explicit VFPInjPropertiesLegacy(const VFPInjTable* inj_table) :
      VFPInjProperties(inj_table) {}

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_tables A map of different VFPINJ tables.
     */
    explicit VFPInjPropertiesLegacy(const InjTable& inj_tables) :
      VFPInjProperties(inj_tables) {}

    using VFPInjProperties::bhp;
    /**
     * Linear interpolation of bhp as function of the input parameters.
     * @param table_id Table number to use
     * @param wells Wells structure with information about wells in qs
     * @param qs Flow quantities
     * @param thp Tubing head pressure
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    ADB bhp(const std::vector<int>& table_id,
            const Wells& wells,
            const ADB& qs,
            const ADB& thp) const;

    /**
     * Linear interpolation of bhp as a function of the input parameters given as ADBs
     * Each entry corresponds typically to one well.
     * @param table_id Table number to use. A negative entry (e.g., -1)
     *                 will indicate that no table is used, and the corresponding
     *                 BHP will be calculated as a constant -1e100.
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table, for each entry in the
     * input ADB objects.
     */
    ADB bhp(const std::vector<int>& table_id,
            const ADB& aqua,
            const ADB& liquid,
            const ADB& vapour,
            const ADB& thp) const;
};


} //namespace



#endif /* OPM_AUTODIFF_VFPINJPROPERTIES_HPP_ */
