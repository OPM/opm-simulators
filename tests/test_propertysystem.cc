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
 * \brief This file tests the properties system.
 *
 * We define a few type tags and property tags, then we attach values
 * to (TypeTag, PropertyTag) tuples and finally we use them in the
 * main function and print some diagnostic messages.
 */
#include "config.h"

#include <opm/models/utils/propertysystem.hh>

#include <iostream>

namespace Opm::Properties {

///////////////////
// Define some hierarchy of type tags:
//
//  Vehicle -- CompactCar -- Sedan -_
//     \                            \.
//      \                            +- Pickup ---_
//       \                          /              \.
//        +-- Truck ---------------^                \.
//         \                                         \.
//          +- Tank ----------------------------------+- HummerH1
///////////////////
namespace TTag {
struct Vehicle {};

struct CompactCar { using InheritsFrom = std::tuple<Vehicle>; };
struct Truck { using InheritsFrom = std::tuple<Vehicle>; };
struct Tank { using InheritsFrom = std::tuple<Vehicle>; };

struct Sedan { using InheritsFrom = std::tuple<CompactCar>; };
struct Pickup { using InheritsFrom = std::tuple<Truck, Sedan>; };

struct HummerH1 { using InheritsFrom = std::tuple<Pickup, Tank, Sedan>; };
} // end namespace TTag

///////////////////
// Define the property tags:
// TopSpeed, NumSeats, CanonCaliber, GasUsage, AutomaticTransmission, Payload
///////////////////
template<class TypeTag, class MyTypeTag>
struct TopSpeed { using type = UndefinedProperty; }; // [km/h]
template<class TypeTag, class MyTypeTag>
struct NumSeats { using type = UndefinedProperty; }; // []
template<class TypeTag, class MyTypeTag>
struct CanonCaliber { using type = UndefinedProperty; }; // [mm]
template<class TypeTag, class MyTypeTag>
struct GasUsage { using type = UndefinedProperty; }; // [l/100km]
template<class TypeTag, class MyTypeTag>
struct AutomaticTransmission { using type = UndefinedProperty; }; // true/false
template<class TypeTag, class MyTypeTag>
struct Payload { using type = UndefinedProperty; }; // [t]

///////////////////
// Make the AutomaticTransmission default to false
template<class TypeTag>
struct AutomaticTransmission<TypeTag, TTag::Vehicle> { static constexpr bool value = false; };

///////////////////
// Define some values for the properties on the type tags:
//
// (CompactCar, TopSpeed) = GasUsage*35
// (CompactCar, NumSeats) = 5
// (CompactCar, GasUsage) = 4
//
// (Truck, TopSpeed) = 100
// (Truck, NumSeats) = 2
// (Truck, GasUsage) = 12
// (Truck, Payload) = 35
//
// (Tank, TopSpeed) = 60
// (Tank, GasUsage) = 65
// (Tank, CanonCaliber) = 120
//
// (Sedan, GasUsage) = 7
// (Sedan, AutomaticTransmission) = true
//
// (Pickup, TopSpeed) = 120
// (Pickup, Payload) = 5
//
// (HummmerH1, TopSpeed) = (Pickup, TopSpeed)
///////////////////

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::CompactCar> { static constexpr int value = getPropValue<TypeTag, GasUsage>() * 30; };
template<class TypeTag>
struct NumSeats<TypeTag, TTag::CompactCar> { static constexpr int value = 5; };
template<class TypeTag>
struct GasUsage<TypeTag, TTag::CompactCar> { static constexpr int value = 4; };

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::Truck> { static constexpr int value = 100; };
template<class TypeTag>
struct NumSeats<TypeTag, TTag::Truck> { static constexpr int value = 2; };
template<class TypeTag>
struct GasUsage<TypeTag, TTag::Truck> { static constexpr int value = 12; };
template<class TypeTag>
struct Payload<TypeTag, TTag::Truck> { static constexpr int value = 35; };

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::Tank> { static constexpr int value = 60; };
template<class TypeTag>
struct GasUsage<TypeTag, TTag::Tank> { static constexpr int value = 65; };
template<class TypeTag>
struct CanonCaliber<TypeTag, TTag::Tank> { static constexpr int value = 120; };

template<class TypeTag>
struct GasUsage<TypeTag, TTag::Sedan> { static constexpr int value = 7; };
template<class TypeTag>
struct AutomaticTransmission<TypeTag, TTag::Sedan> { static constexpr bool value = true; };

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::Pickup> { static constexpr int value = 120; };
template<class TypeTag>
struct Payload<TypeTag, TTag::Pickup> { static constexpr int value = 5; };

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::HummerH1> { static constexpr int value = getPropValue<TTag::Pickup, TopSpeed>(); };

} // namespace Opm::Properties


int main()
{
    using namespace Opm;
    using namespace Opm::Properties;

    // print all properties for all type tags
    std::cout << "---------------------------------------\n";
    std::cout << "-- Property values\n";
    std::cout << "---------------------------------------\n";

    std::cout << "---------- Values for CompactCar ----------\n";

    std::cout << "(CompactCar, TopSpeed) = " << getPropValue<TTag::CompactCar, TopSpeed>() << "\n";
    std::cout << "(CompactCar, NumSeats) = " << getPropValue<TTag::CompactCar, NumSeats>() << "\n";
    std::cout << "(CompactCar, GasUsage) = " << getPropValue<TTag::CompactCar, GasUsage>() << "\n";
    std::cout << "(CompactCar, AutomaticTransmission) = " << getPropValue<TTag::CompactCar, AutomaticTransmission>() << "\n";

    std::cout << "---------- Values for Truck ----------\n";

    std::cout << "(Truck, TopSpeed) = " << getPropValue<TTag::Truck, TopSpeed>() << "\n";
    std::cout << "(Truck, NumSeats) = " << getPropValue<TTag::Truck, NumSeats>() << "\n";
    std::cout << "(Truck, GasUsage) = " << getPropValue<TTag::Truck, GasUsage>() << "\n";
    std::cout << "(Truck, Payload) = " << getPropValue<TTag::Truck, Payload>() << "\n";
    std::cout << "(Truck, AutomaticTransmission) = " << getPropValue<TTag::Truck, AutomaticTransmission>() << "\n";

    std::cout << "---------- Values for Tank ----------\n";

    std::cout << "(Tank, TopSpeed) = " << getPropValue<TTag::Tank, TopSpeed>() << "\n";
    std::cout << "(Tank, GasUsage) = " << getPropValue<TTag::Tank, GasUsage>() << "\n";
    std::cout << "(Tank, AutomaticTransmission) = " << getPropValue<TTag::Tank, AutomaticTransmission>() << "\n";
    std::cout << "(Tank, CanonCaliber) = " << getPropValue<TTag::Tank, CanonCaliber>() << "\n";

    std::cout << "---------- Values for Sedan ----------\n";

    std::cout << "(Sedan, TopSpeed) = " << getPropValue<TTag::Sedan, TopSpeed>() << "\n";
    std::cout << "(Sedan, NumSeats) = " << getPropValue<TTag::Sedan, NumSeats>() << "\n";
    std::cout << "(Sedan, GasUsage) = " << getPropValue<TTag::Sedan, GasUsage>() << "\n";
    std::cout << "(Sedan, AutomaticTransmission) = " << getPropValue<TTag::Sedan, AutomaticTransmission>() << "\n";

    std::cout << "---------- Values for Pickup ----------\n";
    std::cout << "(Pickup, TopSpeed) = " << getPropValue<TTag::Pickup, TopSpeed>() << "\n";
    std::cout << "(Pickup, NumSeats) = " << getPropValue<TTag::Pickup, NumSeats>() << "\n";
    std::cout << "(Pickup, GasUsage) = " << getPropValue<TTag::Pickup, GasUsage>() << "\n";
    std::cout << "(Pickup, Payload) = " << getPropValue<TTag::Pickup, Payload>() << "\n";
    std::cout << "(Pickup, AutomaticTransmission) = " << getPropValue<TTag::Pickup, AutomaticTransmission>() << "\n";

    std::cout << "---------- Values for HummerH1 ----------\n";
    std::cout << "(HummerH1, TopSpeed) = " << getPropValue<TTag::HummerH1, TopSpeed>() << "\n";
    std::cout << "(HummerH1, NumSeats) = " << getPropValue<TTag::HummerH1, NumSeats>() << "\n";
    std::cout << "(HummerH1, GasUsage) = " << getPropValue<TTag::HummerH1, GasUsage>() << "\n";
    std::cout << "(HummerH1, Payload) = " << getPropValue<TTag::HummerH1, Payload>() << "\n";
    std::cout << "(HummerH1, AutomaticTransmission) = " << getPropValue<TTag::HummerH1, AutomaticTransmission>() << "\n";

    return 0;
}
