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
};

template <class OriginalPreconditioner>
class DummyUpdatePreconditioner : public PreconditionerWithUpdate<typename OriginalPreconditioner::domain_type,
                                                                  typename OriginalPreconditioner::range_type>
{
public:
    template <class... Args>
    DummyUpdatePreconditioner(Args&&... args)
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

private:
    OriginalPreconditioner orig_precond_;
};

template <class OriginalPreconditioner>
struct AbstractPreconditionerWithUpdateMaker
{
    virtual std::unique_ptr<Preconditioner<typename OriginalPreconditioner::domain_type, typename OriginalPreconditioner::range_type>> make() = 0;
    // virtual ~AbstractMaker() {};
};

template <class OriginalPreconditioner, class... Args>
struct PreconditionerMaker : public AbstractPreconditionerWithUpdateMaker<OriginalPreconditioner>
{
    PreconditionerMaker(Args&&... args)
        : args_(args...)
    {
    }
    std::unique_ptr<Preconditioner<typename OriginalPreconditioner::domain_type, typename OriginalPreconditioner::range_type>> make() override
    {
        return std::unique_ptr<Preconditioner<typename OriginalPreconditioner::domain_type, typename OriginalPreconditioner::range_type>>{new auto(std::make_from_tuple<OriginalPreconditioner>(args_))};
    }
    std::tuple<Args...> args_;
};

// struct Contain
// {
//     Contain()
//     {
//     }

//     template <class Product, class... Args>
//     void fill(Args&&... args)
//     {
//         using MakerType = Maker<Product, Args...>;
//         maker = std::make_unique<MakerType>(std::forward<Args>(args)...);
//         recreate();
//     }
//     void recreate()
//     {
//         product = maker->make();
//     }
//     std::unique_ptr<AbstractMaker> maker;
//     std::unique_ptr<AbstractProduct> product;
// };

template <class OriginalPreconditioner, class... Args>
std::shared_ptr<DummyUpdatePreconditioner<OriginalPreconditioner>>
getDummyUpdateWrapper(Args&&... args)
{
    return std::make_shared<DummyUpdatePreconditioner<OriginalPreconditioner>>(std::forward<Args>(args)...);
}

template <class OriginalPreconditioner, class... Args>
class RebuildOnUpdatePreconditioner : public PreconditionerWithUpdate<typename OriginalPreconditioner::domain_type,
                                                                      typename OriginalPreconditioner::range_type>
{
public:
    RebuildOnUpdatePreconditioner(Args... args)
        : preconditioner_params_(args...),
        preconditioner_maker_(std::make_unique<ConcreteMakerType>(std::forward<Args>(args)...))
    {
        update();
    }

    using X = typename OriginalPreconditioner::domain_type;
    using Y = typename OriginalPreconditioner::range_type;

    virtual void pre(X& x, Y& b) override
    {
        orig_precond_->pre(x, b);
    }

    virtual void apply(X& v, const Y& d) override
    {
        orig_precond_->apply(v, d);
    }

    virtual void post(X& x) override
    {
        orig_precond_->post(x);
    }

    virtual SolverCategory::Category category() const override
    {
        return orig_precond_->category();
    }

    // Rebuild the preconditioner on update
    void update() override
    {
        // orig_precond_ = std::make_unique<OriginalPreconditioner>(mat_, n_, w_, resort_);
        orig_precond_ = preconditioner_maker_->make();
    }

private:
    using AbstractMakerType = AbstractPreconditionerWithUpdateMaker<OriginalPreconditioner>;
    using ConcreteMakerType = PreconditionerMaker<OriginalPreconditioner, Args...>;

    std::tuple<Args...> preconditioner_params_;
    std::unique_ptr<AbstractMakerType> preconditioner_maker_;
    // std::unique_ptr<ConcreteMakerType> preconditioner_maker_;
    // std::unique_ptr<OriginalPreconditioner> orig_precond_;
    std::unique_ptr<Preconditioner<typename OriginalPreconditioner::domain_type,typename OriginalPreconditioner::range_type>> orig_precond_;
};

template <class OriginalPreconditioner, class... Args>
std::shared_ptr<RebuildOnUpdatePreconditioner<OriginalPreconditioner, Args...>>
getRebuildOnUpdateWrapper(Args... args)
{
    return std::make_shared<RebuildOnUpdatePreconditioner<OriginalPreconditioner, Args...>>(std::forward<Args>(args)...);
}

// template <class OriginalPreconditioner, class Matrix>
// std::shared_ptr<RebuildOnUpdatePreconditioner<OriginalPreconditioner, Matrix>>
// getRebuildOnUpdateWrapper(const Matrix &mat, const int n, const double w, const bool resort)
// {
//     return std::make_shared<RebuildOnUpdatePreconditioner<OriginalPreconditioner, Matrix>>(mat, n, w, resort);
// }

} // namespace Dune

#endif // OPM_PRECONDITIONERWITHUPDATE_HEADER_INCLUDED
