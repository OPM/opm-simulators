/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_PRECONDITIONERWITHUPDATE_HEADER_INCLUDED
#define OPM_PRECONDITIONERWITHUPDATE_HEADER_INCLUDED

#include <dune/istl/preconditioner.hh>
#include <memory>

namespace Dune
{

/// Interface class adding the update() method to the preconditioner interface.
template <class X, class Y>
class PreconditionerWithUpdate : public Preconditioner<X, Y>
{
public:
    virtual void update() = 0;

    // Force derived classes to define if preconditioner has perfect update
    virtual bool hasPerfectUpdate() const = 0;
};

template <class OriginalPreconditioner>
class DummyUpdatePreconditioner : public PreconditionerWithUpdate<typename OriginalPreconditioner::domain_type,
                                                                  typename OriginalPreconditioner::range_type>
{
public:
    template <class... Args>
    explicit DummyUpdatePreconditioner(Args&&... args)
        : orig_precond_(std::forward<Args>(args)...)
    {
    }

    using X = typename OriginalPreconditioner::domain_type;
    using Y = typename OriginalPreconditioner::range_type;

    virtual void pre(X& x, Y& b) override
    {
        orig_precond_.pre(x, b);
    }

    virtual void apply(X& v, const Y& d) override
    {
        orig_precond_.apply(v, d);
    }

    virtual void post(X& x) override
    {
        orig_precond_.post(x);
    }

    virtual SolverCategory::Category category() const override
    {
        return orig_precond_.category();
    }

    // The update() function does nothing for a wrapped preconditioner.
    virtual void update() override
    {
    }

    virtual bool hasPerfectUpdate() const override {
        return true;
    }

private:
    OriginalPreconditioner orig_precond_;
};

template <class OriginalPreconditioner, class... Args>
std::shared_ptr<DummyUpdatePreconditioner<OriginalPreconditioner>>
getDummyUpdateWrapper(Args&&... args)
{
    return std::make_shared<DummyUpdatePreconditioner<OriginalPreconditioner>>(std::forward<Args>(args)...);
}

/// @brief Interface class ensuring make function is overriden
/// @tparam OriginalPreconditioner - An arbitrary Preconditioner type
template <class OriginalPreconditioner>
struct GeneralPreconditionerMaker {
    virtual std::unique_ptr<
        Preconditioner<typename OriginalPreconditioner::domain_type, typename OriginalPreconditioner::range_type>>
    make() = 0;
    virtual ~GeneralPreconditionerMaker() = default;
};

/// @brief Struct implementing a make function which creates a preconditioner
/// Can create arbitrary preconditioners using parameter packs as template arguments
/// @tparam OriginalPreconditioner - An arbitrary preconditioner type
/// @tparam ...Args - All arguments needed to construct the preconditioner of choice
template <class OriginalPreconditioner, class... Args>
struct PreconditionerMaker : public GeneralPreconditionerMaker<OriginalPreconditioner> {
    using GenericPreconditioner = Preconditioner<typename OriginalPreconditioner::domain_type, typename OriginalPreconditioner::range_type>;

    explicit PreconditionerMaker(Args&&... args)
        : args_(args...)
    {
    }

    std::unique_ptr<GenericPreconditioner>
    make() override
    {
        // return std::unique_ptr<GenericPreconditioner> {new auto(std::make_from_tuple<OriginalPreconditioner>(args_))};
        return std::apply(
            [](auto&&... args) {
                return std::make_unique<OriginalPreconditioner>(std::forward<Args>(args)...);
            }, args_);
    }

    std::tuple<Args...> args_;
};

/// @brief Wrapper class of preconditioners that should be reconstructed on update
/// @tparam OriginalPreconditioner - Preconditioner of your choice
template <class OriginalPreconditioner>
class RebuildOnUpdatePreconditioner : public PreconditionerWithUpdate<typename OriginalPreconditioner::domain_type,
                                                                      typename OriginalPreconditioner::range_type>
{
public:
    template<class... Args>
    explicit RebuildOnUpdatePreconditioner(Args... args)
        : preconditioner_maker_(std::make_unique<PreconditionerMaker<OriginalPreconditioner, Args...>>(std::forward<Args>(args)...))
    {
        update();
    }

    using X = typename OriginalPreconditioner::domain_type;
    using Y = typename OriginalPreconditioner::range_type;

    void pre(X& x, Y& b) override
    {
        orig_precond_->pre(x, b);
    }

    void apply(X& v, const Y& d) override
    {
        orig_precond_->apply(v, d);
    }

    void post(X& x) override
    {
        orig_precond_->post(x);
    }

    SolverCategory::Category category() const override
    {
        return orig_precond_->category();
    }

    // Rebuild the preconditioner on update
    void update() override
    {
        orig_precond_ = preconditioner_maker_->make();
    }

    bool hasPerfectUpdate() const override
    {
        return true;
    }

private:
    using AbstractMakerType = GeneralPreconditionerMaker<OriginalPreconditioner>;
    using GenericPreconditioner = Preconditioner<typename OriginalPreconditioner::domain_type, typename OriginalPreconditioner::range_type>;

    std::unique_ptr<AbstractMakerType> preconditioner_maker_;
    std::unique_ptr<GenericPreconditioner> orig_precond_;
};

/// @brief Wrapper function creating and return a shared pointer to a preconditioner which is reconstructed on update
/// @tparam OriginalPreconditioner - Preconditioner of your choice
/// @tparam ...Args - Types of the arguments needed in the preconditioner constructed
/// @param ...args - Arguments needed to construct the preconditioner of choice
/// @return Shared pointer to preconditioner which has an update function that reconstrcuts the preconditioner
template <class OriginalPreconditioner, class... Args>
auto
getRebuildOnUpdateWrapper(Args... args)
{
    return std::make_shared<RebuildOnUpdatePreconditioner<OriginalPreconditioner>>(
        std::forward<Args>(args)...);
}

} // namespace Dune

#endif // OPM_PRECONDITIONERWITHUPDATE_HEADER_INCLUDED
