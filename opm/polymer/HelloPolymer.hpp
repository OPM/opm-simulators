/*===========================================================================
//
// File: HelloPolymer.hpp
//
// Created: 2012-02-01 16:15:27+0100
//
// Authors: Knut-Andreas Lie      <Knut-Andreas.Lie@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Xavier Raynaud        <Xavier.Raynaud@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_HELLOPOLYMER_HPP_HEADER
#define OPM_HELLOPOLYMER_HPP_HEADER

#include <string>

namespace Opm {
    namespace Polymer {
        class Hello {
        public:
            Hello(const ::std::string& dsgn = ::std::string("world"))
                : designee_(dsgn)
            {}

            template <class Ostream>
            friend
            Ostream&
            operator<< (Ostream& os, const Hello& h);

        private:
            ::std::string designee_;
        };

        template <class Ostream>
        Ostream &
        operator<<(Ostream& os, const Hello& h) {
            os << "Polymer::Hello, " << h.designee_;
            return os;
        }
    }  // namespace Polymer
}      // namespace Opm
#endif  /* OPM_HELLOPOLYMER_HPP_HEADER */
