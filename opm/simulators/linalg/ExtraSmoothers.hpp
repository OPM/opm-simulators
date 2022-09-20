#ifndef OPM_EXTRASMOOTHERS_HPP
#define OPM_EXTRASMOOTHERS_HPP

#include "ExtraPrecontioners.hpp"

namespace Dune
{
namespace Amg
{

template <class T>
class ConstructionTraits;

/**
 * @brief Policy for the construction of the SeqJacNew smoother
 */
template <class M, class X, class Y, int l>
struct ConstructionTraits<SeqJacNew<M, X, Y, l>> {
    typedef DefaultConstructionArgs<SeqJacNew<M, X, Y, l>> Arguments;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    static inline std::shared_ptr<SeqJacNew<M, X, Y, l>> construct(Arguments& args)
    {
        return std::make_shared<SeqJacNew<M, X, Y, l>>(args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor);
    }
#else
    static inline SeqJacNew<M, X, Y, l>* construct(Arguments& args)
    {
        return new SeqJacNew<M, X, Y, l>(args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor);
    }

    static void deconstruct(SeqJacNew<M, X, Y, l>* jac)
    {
        delete jac;
    }

#endif
};

/**
 * @brief Policy for the construction of the SeqSpai0 smoother
 */
template <class M, class X, class Y, int l>
struct ConstructionTraits<SeqSpai0<M, X, Y, l>> {
typedef DefaultConstructionArgs<SeqSpai0<M, X, Y, l>> Arguments;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)    
    static inline std::shared_ptr<SeqSpai0<M, X, Y, l>> construct(Arguments& args)
{ 
    return std::make_shared<SeqSpai0<M, X, Y, l>>(args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor);
}

#else
static inline SeqSpai0<M, X, Y, l>* construct(Arguments& args)
{
    return new SeqSpai0<M, X, Y, l>(args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor);
}

static void deconstruct(SeqSpai0<M, X, Y, l>* jac)
{
    delete jac;
}
#endif
};

}
}
#endif // OPM_EXTRASMOOTHERS_HPP
