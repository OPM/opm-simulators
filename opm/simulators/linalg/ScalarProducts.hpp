#ifndef OPM_SCALAR_PRODUCTS_HEADER_INCLUDED
#define OPM_SCALAR_PRODUCTS_HEADER_INCLUDED

namespace Dune
{

/// A parallel scalar product that takes advantage of the fact that all
/// elements associated with ghost cells are located at the end of the
/// vector. This allows us to ignore the block structure of the vector
/// and eliminate the use of a mask to exclude ghost entries from being
/// included in the scalar product
template<class Vector, class Comm>
class GhostLastScalarProduct : public ScalarProduct<Vector>
{
    public:

    ///Exctract block size from vector type
    static constexpr auto block_size = Vector::block_type::dimension;

    /*! \brief constructor
    * \param com The communication object for syncing overlap and copy
    * data points.
    * \param cat parallel solver category (nonoverlapping or overlapping)
    */
    GhostLastScalarProduct (std::shared_ptr<const Comm> com, SolverCategory::Category cat)
        : _communication(com), _category(cat)
    {
        count_ = getLocalCount();        // number or local cells
        int verify = verifyLocalCount(); // redundant check on numbef of local cells
        if (count_ != verify) OPM_THROW(std::runtime_error, "Inconsistent local node count!!\n");
    }

    /*! \brief constructor
    * \param com The communication object for syncing overlap and copy
    * data points.
    * \param cat parallel solver category (nonoverlapping or overlapping)
    * \note if you use this constructor you have to make sure com stays alive
    */
    GhostLastScalarProduct (const Comm& com, SolverCategory::Category cat)
        : GhostLastScalarProduct(stackobject_to_shared_ptr(com), cat)
    {}

    /*! \brief Dot product of two vectors.
    * \param vx first input vector
    * \param vy second input vector
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

        // loop-peeling of trailing end
        for(int j=N;j<NN;j++) agg[0]+=x[j]*y[j];

        // Global summation
        auto cc = _communication->communicator();
        double result = cc.sum(agg[0]);
        return result;
    }

    /*! \brief Vector L2-norm.
    * \param vx input vector
    */
    virtual double norm (const Vector& vx) const override
    {
        return sqrt(dot(vx,vx));
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

    /*! \brief Count number of local cells.
    */
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

    /*! \brief Infer number of local cells from largest local index.
    */
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



/// A sequential scalar product that ignores block structure of the vector
/// to facilitate well-known optimization techniques
template<class Vector>
class SeqOptmizedProduct : public Dune::SeqScalarProduct<Vector>
{
public:

    // extract block size
    static constexpr auto block_size = Vector::block_type::dimension;

    /*! \brief Dot product of two vectors.
    * \param vx first input vector
    * \param vy second input vector
    */
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

        // loop-peeling of trailing end
        for(int j=N;j<NN;j++) agg[0]+=x[j]*y[j];

        return agg[0];
    }

    /*! \brief Vector L2-norm.
    * \param vx input vector
    */
    virtual double norm(const Vector& vx) const override {
        return std::sqrt(this->dot(vx, vx));
    }
};

}

#endif //OPM_SCALAR_PRODUCTS_HEADER_INCLUDED

