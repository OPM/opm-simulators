#ifndef OPM_EXTRASMOOTHERS_HPP
#define OPM_EXTRASMOOTHERS_HPP

#include "DILU.hpp"

namespace Dune
{
    template <class M, class X, class Y>
    class MultithreadDILU;

namespace Amg
{
    /**
     * @brief Policy for the construction of the MultithreadDILU smoother
     */
    template <class M, class X, class Y>
    struct ConstructionTraits<MultithreadDILU<M, X, Y>> {
    using Arguments = DefaultConstructionArgs<MultithreadDILU<M, X, Y>>;

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
        static inline std::shared_ptr<MultithreadDILU<M, X, Y>> construct(Arguments& args) {
            return std::make_shared<MultithreadDILU<M, X, Y>>(args.getMatrix());
        }

#else
        static inline MultithreadDILU<M, X, Y>* construct(Arguments& args) {
            return new MultithreadDILU<M, X, Y>(args.getMatrix());
        }

        static void deconstruct(MultithreadDILU<M, X, Y>* dilu) {
            delete dilu;
        }
#endif
    };

} // namespace Amg
} // namespace Dune
#endif // OPM_EXTRASMOOTHERS_HPP