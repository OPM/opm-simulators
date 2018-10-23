/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Equinor.

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

#ifndef OPM_CONVERGENCESTATUS_HEADER_INCLUDED
#define OPM_CONVERGENCESTATUS_HEADER_INCLUDED

#include <cassert>
#include <string>
#include <vector>

namespace Opm
{

    /// Represents the convergence status of the whole simulator, to
    /// make it possible to query and store the reasons for
    /// convergence failures.
    class ConvergenceStatus
    {
    public:

        // ----------- Types -----------

        enum Status { AllGood            = 0,
                      ReservoirFailed    = 1 << 0,
                      WellFailed         = 1 << 1 };
        struct ReservoirFailure
        {
            enum struct Type { Mb, Cnv };
            Type type;
            int phase;
            int cell_index;
        };
        struct WellFailure
        {
            enum struct Type { Mb, CtrlBHP, CtrlTHP, CtrlRate };
            Type type;
            int phase;
            std::string well_name;
        };

        // ----------- Mutating member functions -----------

        ConvergenceStatus()
            : status_{AllGood}
            , res_failures_{}
            , well_failures_{}
        {
        }

        void clear()
        {
            status_ = AllGood;
            res_failures_.clear();
            well_failures_.clear();
        }

        void setReservoirFailed(const ReservoirFailure& rf)
        {
            status_ = static_cast<Status>(status_ | ReservoirFailed);
            res_failures_.push_back(rf);
        }

        void setWellFailed(const WellFailure& wf)
        {
            status_ = static_cast<Status>(status_ | WellFailed);
            well_failures_.push_back(wf);
        }

        ConvergenceStatus& operator+=(const ConvergenceStatus& other)
        {
            status_ = static_cast<Status>(status_ | other.status_);
            res_failures_.insert(res_failures_.end(), other.res_failures_.begin(), other.res_failures_.end());
            well_failures_.insert(well_failures_.end(), other.well_failures_.begin(), other.well_failures_.end());
            assert(reservoirFailed() != res_failures_.empty());
            assert(wellFailed() != well_failures_.empty());
            return *this;
        }

        // ----------- Const member functions (queries) -----------

        bool converged() const
        {
            return status_ == AllGood;
        }

        bool reservoirFailed() const
        {
            return status_ & ReservoirFailed;
        }

        bool wellFailed() const
        {
            return status_ & WellFailed;
        }

        const std::vector<ReservoirFailure>& reservoirFailures() const
        {
            return res_failures_;
        }

        const std::vector<WellFailure>& wellFailures() const
        {
            return well_failures_;
        }

    private:

        // ----------- Member variables -----------
        Status status_;
        std::vector<ReservoirFailure> res_failures_;
        std::vector<WellFailure> well_failures_;
    };

} // namespace Opm

#endif // OPM_CONVERGENCESTATUS_HEADER_INCLUDED
