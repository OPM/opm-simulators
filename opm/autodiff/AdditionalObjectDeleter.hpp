 /*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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
#ifndef OPM_ADDITIONALOBJECTDELETER_HEADER_INCLUDED
#define OPM_ADDITIONALOBJECTDELETER_HEADER_INCLUDED
namespace Opm
{
//! \brief A custom deleter that will delete an additional object, too.
//!
//! In dune-istl most parallel preconditioners hold a reference to
//! a sequential preconditioner.
//! In CPRPreconditioner and NewtonIterationBlackoilInterleaved we use unique_ptr
//! for the memory management.
//! Ergo we need to construct the sequential preconditioner with new and
//! make sure that it gets deleted together with the enclosing parallel
//! preconditioner. Therefore this deleter stores a pointer to it and deletes
//! it during destruction.
//! \tparam The type of the additional object to be deleted.
template<class T>
class AdditionalObjectDeleter
{
public:
    //! \brief empty constructor.
    AdditionalObjectDeleter()
        : additional_object_()
    {}
    //! \brief Constructor taking the object that needs to deleted.
    AdditionalObjectDeleter(T& additional_object)
    : additional_object_(&additional_object){}
    //! \brief Delete an object and the additional one.
    template<class T1>
    void operator()(T1* pt)
    {
        delete pt;
        delete additional_object_;
    }
private:
    T* additional_object_;
};

}
#endif
