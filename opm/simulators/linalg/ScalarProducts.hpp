#ifndef OPM_SCALAR_PRODUCTS_HEADER_INCLUDED
#define OPM_SCALAR_PRODUCTS_HEADER_INCLUDED

namespace Dune
{

template<class Vector, class Comm>
class GhostLastScalarProduct : public ScalarProduct<Vector>
{
    public:

    static constexpr auto block_size = Vector::block_type::dimension;

    /*!
    * \param com The communication object for syncing overlap and copy
    * data points.
    * \param cat parallel solver category (nonoverlapping or overlapping)
    */
    GhostLastScalarProduct (std::shared_ptr<const Comm> com, SolverCategory::Category cat)
        : _communication(com), _category(cat)
    {
        count_ = getLocalCount();
        int verify = verifyLocalCount();
        if (count_ != verify) OPM_THROW(std::runtime_error, "Inconsistent local node count!!\n");
    }

    /*!
    * \param com The communication object for syncing overlap and copy
    * data points.
    * \param cat parallel solver category (nonoverlapping or overlapping)
    * \note if you use this constructor you have to make sure com stays alive
    */
    GhostLastScalarProduct (const Comm& com, SolverCategory::Category cat)
        : GhostLastScalarProduct(stackobject_to_shared_ptr(com), cat)
    {}

    /*! \brief Dot product of two vectors.
    It is assumed that the vectors are consistent on the interior+border
    partition.
    */
    virtual double dot (const Vector& vx, const Vector& vy) const override
    {

        // access underlying data
        double const *x = &vx[0][0];
        double const *y = &vy[0][0];

        // total array length
        int NN = block_size*count_;

        // unroll loop in multiples of 8
        int n=NN/8;
        int N=8*n;
        double agg[8];
        for(int i=0;i<8;i++) agg[i]=0.0;
        for(int i=0;i<N;i+=8) for(int j=0;j<8;j++) agg[j]+=x[i+j]*y[i+j];
        for(int j=0;j<4;j++) agg[j]+=agg[j+4];
        for(int j=0;j<2;j++) agg[j]+=agg[j+2];
        for(int j=0;j<1;j++) agg[j]+=agg[j+1];

        // trailing end
        for(int j=N;j<NN;j++) agg[0]+=x[j]*y[j];

        // Global summation
        auto cc = _communication->communicator();
        double result = cc.sum(agg[0]);
        return result;
    }

    /*! \brief Norm of a right-hand side vector.
    The vector must be consistent on the interior+border partition
    */
    virtual double norm (const Vector& x) const override
    {
        return sqrt(dot(x,x));
    }

    //! Category of the scalar product (see SolverCategory::Category)
    virtual SolverCategory::Category category() const override
    {
        return _category;
    }

    private:
    std::shared_ptr<const Comm> _communication;
    SolverCategory::Category _category;
    int count_;

    int getLocalCount() const
    {
        int count = 0;
        // Loop over index set
        auto indexSet = _communication->indexSet();
        for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {
            if (idx->local().attribute()==1) count++; // count non-local indices
        }
        return count;
    }

    int verifyLocalCount() const
    {
        auto indexSet = _communication->indexSet();

        size_t is = 0;
        // Loop over index set
        for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {
            //Only take "owner" indices
            if (idx->local().attribute()==1) {
                //get local index
                auto loc = idx->local().local();
                // if loc is higher than "old interior size", update it
                if (loc > is) {
                    is = loc;
                }
            }
        }
        return is + 1; //size is plus 1 since we start at 0
    }

};



//! @brief Optimized sequential scalar product.
//!
//! @tparam Vector block-vector class with data stored as contiguous double array
template<class Vector>
class SeqOptmizedProduct : public Dune::SeqScalarProduct<Vector>
{
public:

    // extract block size
    static constexpr auto block_size = Vector::block_type::dimension;

    // Compute the dot product <vx, vy>
    virtual double dot(const Vector& vx, const Vector& vy) const override
    {
        // access underlying data
        double const *x = &vx[0][0];
        double const *y = &vy[0][0];

        // total array length
        int NN = block_size*vx.N();

        // unroll loop in multiples of 8
        int n=NN/8;
        int N=8*n;
        double agg[8];
        for(int i=0;i<8;i++) agg[i]=0.0;
        for(int i=0;i<N;i+=8) for(int j=0;j<8;j++) agg[j]+=x[i+j]*y[i+j];
        for(int j=0;j<4;j++) agg[j]+=agg[j+4];
        for(int j=0;j<2;j++) agg[j]+=agg[j+2];
        for(int j=0;j<1;j++) agg[j]+=agg[j+1];

        // trailing end
        for(int j=N;j<NN;j++) agg[0]+=x[j]*y[j];

        return agg[0];
    }

    // Compute the norm ||x||
    virtual double norm(const Vector& x) const override {
        return std::sqrt(this->dot(x, x));
    }
};

}

#endif //OPM_SCALAR_PRODUCTS_HEADER_INCLUDED

