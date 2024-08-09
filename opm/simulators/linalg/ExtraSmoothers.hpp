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

        static inline std::shared_ptr<MultithreadDILU<M, X, Y>> construct(Arguments& args) {
            return std::make_shared<MultithreadDILU<M, X, Y>>(args.getMatrix());
        }

    };

} // namespace Amg
} // namespace Dune
#endif // OPM_EXTRASMOOTHERS_HPP
