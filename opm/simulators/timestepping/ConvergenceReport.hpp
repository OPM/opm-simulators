/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018, 2024 Equinor.

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

#ifndef OPM_CONVERGENCEREPORT_HEADER_INCLUDED
#define OPM_CONVERGENCEREPORT_HEADER_INCLUDED

#include <algorithm>
#include <cassert>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{

    /// Represents the convergence status of the whole simulator, to
    /// make it possible to query and store the reasons for
    /// convergence failures.
    class ConvergenceReport
    {
    public:

        // ----------- Types -----------

        enum Status {
            AllGood         = 0,
            ReservoirFailed = 1 << 0,
            WellFailed      = 1 << 1,
        };
        // More severe problems should have higher numbers
        enum struct Severity {
            None       = 0,
            Normal     = 1,
            ConvergenceMonitorFailure = 2,
            TooLarge   = 3,
            NotANumber = 4,
        };

        struct PenaltyCard {
            int nonConverged{0};
            int distanceDecay{0};
            int largeWellResiduals{0};

            int total() const {
                return nonConverged + distanceDecay + largeWellResiduals;
            }

            void reset()
            {
                nonConverged = 0;
                distanceDecay = 0;
                largeWellResiduals = 0;
            }

            PenaltyCard& operator+=(const PenaltyCard& other) {
                nonConverged += other.nonConverged;
                distanceDecay += other.distanceDecay;
                largeWellResiduals += other.largeWellResiduals;
                return *this;
            }

            template <typename Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(nonConverged);
                serializer(distanceDecay);
                serializer(largeWellResiduals);
            }
        };

        using CnvPvSplit = std::pair<
            std::vector<double>,
            std::vector<int>>;

        class ReservoirFailure
        {
        public:
            enum struct Type { Invalid, MassBalance, Cnv, ConvergenceMonitorFailure };

            // Default constructor needed for object serialisation.  Don't
            // use this for anything else.
            ReservoirFailure() = default;

            ReservoirFailure(Type t, Severity s, int phase)
                : type_(t), severity_(s), phase_(phase)
            {}

            Type type() const { return type_; }
            Severity severity() const { return severity_; }
            int phase() const { return phase_; }

            template <typename Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(this->type_);
                serializer(this->severity_);
                serializer(this->phase_);
            }

        private:
            // Note to maintainers: If you change this list of data members,
            // then please update serializeOp() accordingly.
            Type type_ { Type::Invalid };
            Severity severity_ { Severity::None };
            int phase_ { -1 };
        };

        class ReservoirConvergenceMetric
        {
        public:
            // Default constructor needed for object serialisation.  Don't
            // use this for anything else.
            ReservoirConvergenceMetric() = default;

            ReservoirConvergenceMetric(ReservoirFailure::Type t, int phase, double value, double tolerance)
                : type_(t), phase_(phase), value_(value), tolerance_(tolerance)
            {}

            ReservoirFailure::Type type() const { return type_; }
            int phase() const { return phase_; }
            double value() const { return value_; }
            double tolerance() const { return tolerance_; }

            template <typename Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(this->type_);
                serializer(this->phase_);
                serializer(this->value_);
                serializer(this->tolerance_);
            }

        private:
            // Note to maintainers: If you change this list of data members,
            // then please update serializeOp() accordingly.
            ReservoirFailure::Type type_ { ReservoirFailure::Type::Invalid };
            int phase_ { -1 };
            double value_ { 0.0 };
            double tolerance_ { 0.0 };
        };

        class WellFailure
        {
        public:
            enum struct Type {
                Invalid,
                MassBalance,
                Pressure,
                ControlBHP,
                ControlTHP,
                ControlRate,
                Unsolvable,
                WrongFlowDirection,
            };

            // Default constructor needed for object serialisation.  Don't
            // use this for anything else.
            WellFailure() = default;

            WellFailure(Type t, Severity s, int phase, const std::string& well_name)
                : type_(t), severity_(s), phase_(phase), well_name_(well_name)
            {}

            Type type() const { return type_; }
            Severity severity() const { return severity_; }
            int phase() const { return phase_; }
            const std::string& wellName() const { return well_name_; }

            template <typename Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(this->type_);
                serializer(this->severity_);
                serializer(this->phase_);
                serializer(this->well_name_);
            }

        private:
            // Note to maintainers: If you change this list of data members,
            // then please update serializeOp() accordingly.
            Type type_ { Type::Invalid };
            Severity severity_ { Severity::None };
            int phase_ { -1 };
            std::string well_name_ {};
        };

        class WellConvergenceMetric
        {
        public:
            // Default constructor needed for object serialisation.  Don't
            // use this for anything else.
            WellConvergenceMetric() = default;

            WellConvergenceMetric(WellFailure::Type t, Severity s, int phase, double value, const std::string& well_name)
                : type_(t), severity_(s), phase_(phase), value_(value), well_name_(well_name)
            {}

            WellFailure::Type type() const { return type_; }
            Severity severity() const { return severity_; }
            int phase() const { return phase_; }
            double value() const { return value_; }
            const std::string& wellName() const { return well_name_; }

            template <typename Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(this->type_);
                serializer(this->severity_);
                serializer(this->phase_);
                serializer(this->value_);
                serializer(this->well_name_);
            }

        private:
            // Note to maintainers: If you change this list of data members,
            // then please update serializeOp() accordingly.
            WellFailure::Type type_ { WellFailure::Type::Invalid };
            Severity severity_ { Severity::None };
            int phase_ { -1 };
            double value_ { 0.0 };
            std::string well_name_ {};
        };

        // ----------- Mutating member functions -----------

        ConvergenceReport()
            : ConvergenceReport{0.0}
        {}

        explicit ConvergenceReport(const double reportTime)
            : reportTime_{reportTime}
            , status_{AllGood}
            , res_failures_{}
            , well_failures_{}
            , wellGroupTargetsViolated_(false)
            , network_needs_more_balancing_force_another_newton_iteration_(false)
        {}

        void clear()
        {
            status_ = AllGood;
            res_failures_.clear();
            well_failures_.clear();
            wellGroupTargetsViolated_ = false;
            network_needs_more_balancing_force_another_newton_iteration_ = false;
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

        template <typename... Args>
        void setReservoirConvergenceMetric(Args&&... args)
        {
            this->res_convergence_.emplace_back(std::forward<Args>(args)...);
        }

        template <typename... Args>
        void setWellConvergenceMetric(Args&&... args)
        {
            this->well_convergence_.emplace_back(std::forward<Args>(args)...);
        }

        void setWellGroupTargetsViolated(const bool wellGroupTargetsViolated)
        {
            wellGroupTargetsViolated_ = wellGroupTargetsViolated;
        }

        void setNetworkNotYetBalancedForceAnotherNewtonIteration(const bool network_needs_more_balancing_force_another_newton_iteration)
        {
            network_needs_more_balancing_force_another_newton_iteration_ = network_needs_more_balancing_force_another_newton_iteration;
        }

        void setCnvPoreVolSplit(const CnvPvSplit& cnvPvSplit,
                                const double eligiblePoreVolume)
        {
            this->cnvPvSplit_ = cnvPvSplit;
            this->eligiblePoreVolume_ = eligiblePoreVolume;
        }

        ConvergenceReport& operator+=(const ConvergenceReport& other)
        {
            reportTime_ = std::max(reportTime_, other.reportTime_);
            status_ = static_cast<Status>(status_ | other.status_);
            res_failures_.insert(res_failures_.end(), other.res_failures_.begin(), other.res_failures_.end());
            well_failures_.insert(well_failures_.end(), other.well_failures_.begin(), other.well_failures_.end());
            res_convergence_.insert(res_convergence_.end(), other.res_convergence_.begin(), other.res_convergence_.end());
            well_convergence_.insert(well_convergence_.end(), other.well_convergence_.begin(), other.well_convergence_.end());
            assert(reservoirFailed() != res_failures_.empty());
            assert(wellFailed() != well_failures_.empty());
            wellGroupTargetsViolated_ = (wellGroupTargetsViolated_ || other.wellGroupTargetsViolated_);
            network_needs_more_balancing_force_another_newton_iteration_ = (network_needs_more_balancing_force_another_newton_iteration_
                || other.network_needs_more_balancing_force_another_newton_iteration_);

            // Note regarding the CNV pore-volume split: We depend on the
            // fact that the quantities have already been aggregated across
            // all MPI ranks--see the implementation of member function
            // BlackoilModel::getReservoirConvergence() for details--and are
            // therefore equal on all ranks.  Consequently, we simply assign
            // 'other's values here, if it is non-empty.  Empty splits
            // typically come from well contributions.
            if (! other.cnvPvSplit_.first.empty()) {
                this->cnvPvSplit_ = other.cnvPvSplit_;
                this->eligiblePoreVolume_ = other.eligiblePoreVolume_;
            }

            return *this;
        }

        // ----------- Const member functions (queries) -----------

        double reportTime() const
        {
            return reportTime_;
        }

        double eligiblePoreVolume() const
        {
            return this->eligiblePoreVolume_;
        }

        const CnvPvSplit& cnvPvSplit() const
        {
            return this->cnvPvSplit_;
        }

        bool converged() const
        {
            return (status_ == AllGood)
                && !wellGroupTargetsViolated_
                && !network_needs_more_balancing_force_another_newton_iteration_;
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

        const std::vector<ReservoirConvergenceMetric>& reservoirConvergence() const
        {
            return res_convergence_;
        }

        const std::vector<WellFailure>& wellFailures() const
        {
            return well_failures_;
        }

        const std::vector<WellConvergenceMetric>& wellConvergence() const
        {
            return well_convergence_;
        }

        const PenaltyCard& getPenaltyCard() const
        {
            return penaltyCard_;
        }

        void addNonConvergedPenalty()
        {
            penaltyCard_.nonConverged++;
        }

        void addDistanceDecayPenalty()
        {
            penaltyCard_.distanceDecay++;
        }

        void addLargeWellResidualsPenalty()
        {
            penaltyCard_.largeWellResiduals++;
        }

        Severity severityOfWorstFailure() const
        {
            // A function to get the worst of two severities.
            auto smax = [](Severity s1, Severity s2) {
                return s1 < s2 ? s2 : s1;
            };
            auto s = Severity::None;
            for (const auto& f : res_failures_) {
                s = smax(s, f.severity());
            }
            for (const auto& f : well_failures_) {
                s = smax(s, f.severity());
            }
            return s;
        }

        template <typename Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->reportTime_);
            serializer(this->status_);
            serializer(this->res_failures_);
            serializer(this->well_failures_);
            serializer(this->res_convergence_);
            serializer(this->well_convergence_);
            serializer(this->wellGroupTargetsViolated_);
            serializer(this->network_needs_more_balancing_force_another_newton_iteration_);
            serializer(this->cnvPvSplit_);
            serializer(this->eligiblePoreVolume_);
            serializer(this->penaltyCard_);
        }

    private:
        // ----------- Member variables -----------
        // Note to maintainers: If you change this list of data members,
        // then please update serializeOp() accordingly.
        double reportTime_;
        Status status_;
        std::vector<ReservoirFailure> res_failures_;
        std::vector<WellFailure> well_failures_;
        std::vector<ReservoirConvergenceMetric> res_convergence_;
        std::vector<WellConvergenceMetric> well_convergence_;
        bool wellGroupTargetsViolated_;
        bool network_needs_more_balancing_force_another_newton_iteration_;
        CnvPvSplit cnvPvSplit_{};
        double eligiblePoreVolume_{};
        PenaltyCard penaltyCard_;
    };

    struct StepReport
    {
        int report_step;
        int current_step;
        std::vector<ConvergenceReport> report;
    };

    std::string to_string(const ConvergenceReport::ReservoirFailure::Type t);

    std::string to_string(const ConvergenceReport::Severity s);

    std::string to_string(const ConvergenceReport::WellFailure::Type t);

    std::string to_string(const ConvergenceReport::WellFailure& wf);

    std::string to_string(const ConvergenceReport::PenaltyCard& pc);



} // namespace Opm

#endif // OPM_CONVERGENCEREPORT_HEADER_INCLUDED
