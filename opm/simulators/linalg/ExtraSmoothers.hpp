#ifndef OPM_EXTRASMOOTHERS_HPP
#define OPM_EXTRASMOOTHERS_HPP

#include "DILU.hpp"

namespace Dune
{
    template <class M, class X, class Y>
    class SeqDilu;

namespace Amg
{
    /**
     * @brief Policy for the construction of the SeqDilu smoother
     */
    template <class M, class X, class Y>
    struct ConstructionTraits<SeqDilu<M, X, Y>> {
    using Arguments = DefaultConstructionArgs<SeqDilu<M, X, Y>>;

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
        static inline std::shared_ptr<SeqDilu<M, X, Y>> construct(Arguments& args) {
            return std::make_shared<SeqDilu<M, X, Y>>(args.getMatrix());
        }

#else
        static inline SeqDilu<M, X, Y>* construct(Arguments& args) {
            return new SeqDilu<M, X, Y>(args.getMatrix());
        }

        static void deconstruct(SeqDilu<M, X, Y>* dilu) {
            delete dilu;
        }
#endif
    };

} // namespace Amg
} // namespace Dune
#endif // OPM_EXTRASMOOTHERS_HPP