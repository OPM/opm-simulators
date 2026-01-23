#pragma once
namespace Dune
{
    class SeqComm
    {
    public:
        using size_type = std::size_t;
        // SolverCategory::Category category () const {
        //     return category_;
        // }

        // const Communication<MPI_Comm>& communicator() const
        // {
        //   return cc;
        // }
        static constexpr size_type size()
        {
            return 0;
        }
        template <typename T>
        void project(T & /*x*/) const
        {
            // No operation for sequential communicator
        }
        template <typename T1, typename T2>
        void dot(const T1 &x, const T1 &y, T2 &result) const
        {
            result = x.dot(y);
        }
        template <typename T>
        double norm(const T &x) const
        {
            return x.two_norm();
        };
        template <typename T>
        void copyOwnerToAll(const T &x,T &y) const
        {
            y = x;
        }
    };

  //template <typename ColCom> Communication<MPI_Comm>
    class JacComm
    {
    public:
        using size_type = std::size_t;
        // SolverCategory::Category category () const {
        //     return category_;
        // }

        // const Communication<MPI_Comm>& communicator() const
        // {
        //   return cc;
        // }
    public:
      //JacComm(const ColCom& colcom): colcom_(colcom){}
      JacComm(): colcom_(MPI_COMM_WORLD){}
        static constexpr size_type size()
        {
            return 0;
        }
        template <typename T>
        void project(T & /*x*/) const
        {
            // No operation for sequential communicator
        }
      
        template <typename T1, typename T2>
        void dot(const T1 &x, const T1 &y, T2 &result) const
        {
            result = x.dot(y);
            result = colcom_.sum(result);
        }
      
        template <typename T>
        double norm(const T &x) const
        {
            //double result = x.two_norm();
            double result = x.dot(x);
            result = colcom_.sum(result);
            result = std::sqrt(result);
            return result;
        }
      
        template <typename T>
        void copyOwnerToAll(const T &x,T &y) const
        {
            y = x;
        }
       private:
       Communication<MPI_Comm> colcom_;
    };
  
    template <typename... Args>
    class MultiCommunicator
        : public std::tuple<Args...>
    {
        /** \brief Helper type */
        typedef std::tuple<Args...> TupleType;
        typedef MultiCommunicator<Args...> type;
        using field_type = double;

    public:
        using std::tuple<Args...>::tuple;
        using size_type = std::size_t;

        static constexpr size_type size()
        {
            return sizeof...(Args);
        }
        
        template <size_type index>
        typename std::tuple_element<index, TupleType>::type &
        operator[]([[maybe_unused]] const std::integral_constant<size_type, index> indexVariable)
        {
            return std::get<index>(*this);
        }
        template <size_type index>
        const typename std::tuple_element<index, TupleType>::type &
        operator[]([[maybe_unused]] const std::integral_constant<size_type, index> indexVariable) const
        {
            return std::get<index>(*this);
        }

        template <typename T>
        void project(T &x) const
        {
            using namespace Dune::Hybrid;
            //auto size = index_constant<sizeof...(Args)>();
            forEach(integralRange(Hybrid::size(*this)), [&](auto &&i)
                    { (*this)[i].project(x[i]); });
        }
        template <typename T1, typename T2>
        void dot(const T1 &x, const T1 &y, T2& result) const
        {
            result = field_type(0);
            using namespace Dune::Hybrid;
            forEach(integralRange(Hybrid::size(*this)), [&](auto &&i)
            { double result_tmp = 0;
              (*this)[i].dot(x[i],y[i],result_tmp);
              result += result_tmp;
              //std::cout << " Dot partial result " << i << " r " << result <<std::endl;
            } );
        }
        template <typename T>
        field_type norm(const T &x) const
        {
            using namespace Dune::Hybrid;
            // return accumulate(integralRange(Hybrid::size(*this)), field_type(0), [&](auto &&a, auto &&i)
            //                   { return a + (*this)[i].norm(x[i]); });

            field_type result(0);
            using namespace Dune::Hybrid;
            forEach(integralRange(Hybrid::size(*this)), [&](auto &&i)
            { 
                //double result_tmp = (*this)[i].norm(x[i]);
                double result_tmp = 0.0;
                (*this)[i].dot(x[i], x[i], result_tmp);
                result += result_tmp;
              //std::cout << " Dot partial result " << i << " r " << result <<std::endl;
            } );
            return std::sqrt(result);
        }
        template <typename T>
        void copyOwnerToAll(const T &x,T &y) const
        {
            using namespace Dune::Hybrid;
            forEach(integralRange(Hybrid::size(*this)), [&](auto &&i)
                    { (*this)[i].copyOwnerToAll(x[i], y[i]); });
        }
        // Explicit template instantiations
    };
}
