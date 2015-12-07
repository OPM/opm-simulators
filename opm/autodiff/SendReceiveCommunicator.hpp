/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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
#ifndef OPM_ASYNCHRONOUSCOMMUNICATOR_HEADER_INCLUDED
#define OPM_ASYNCHRONOUSCOMMUNICATOR_HEADER_INCLUDED

#if HAVE_MPI
// MPI header
#include <mpi.h>
#include <vector>
#include <map>
#include <functional>
#include <dune/common/unused.hh>
#include <dune/common/parallel/interface.hh>
#include <dune/common/parallel/mpitraits.hh>

namespace std
{
inline ostream& operator<<(ostream& os,
                           const Dune::Interface::InformationMap& interface)
{
    typedef Dune::Interface::InformationMap InfoMap;
    typedef InfoMap::const_iterator Iter;
    for(Iter i=interface.begin(), end = interface.end();
        i!=end; ++i)
    {
        os<<i->first<<": [ source=[";
        for(std::size_t j=0; j < i->second.first.size(); ++j)
            os<<i->second.first[j]<<" ";
        os<<"] size="<<i->second.first.size()<<", target=[";
        for(std::size_t j=0; j < i->second.second.size(); ++j)
            os<<i->second.second[j]<<" ";
        os<<"] size="<<i->second.second.size()<<"\n";
    }
    return os;
}
} // end namespace std

namespace Opm
{
namespace Detail
{
/**
 * @brief A message buffer to be used for communication
 * @tparam T The type of data that the buffer will hold.
 */
template<class T, class Allocator=std::allocator<T> >
class MessageBuffer
{
public:
    /**
     * @brief Constructs a message.
     * @param size The number of elements that buffer should hold,
     */
    explicit MessageBuffer(int size)
        : buffer_(size), position_(0)
    {}
    /**
     * @brief Default constructor.
     *
     * Use resize to allocate the memory.
     */
    explicit MessageBuffer()
        : buffer_(), position_(0)
    {}
    /**
     * @brief Copy constructor.
     * @param o The instance to copy.
     */
    explicit MessageBuffer(const MessageBuffer& o)
        : buffer_(o.buffer_), position_(o.position_)
    {
    }
    /** @brief Destructor. */
    ~MessageBuffer()
    {}
    /**
     * @brief Change the size of the buffer.
     *
     * Will throw away allocated buffer together with the entries
     * and allocate new one.
     * @param size The number elements to hold.
     */
    void resize(std::size_t size)
    {
        buffer_.resize(size);
        position_=0;
    }
    /**
     * @brief Free the memory used for the buffer.
     */
    void free()
    {
        std::vector<T>().swap(buffer_);
        position_=0;
    }

    /**
     * @brief Write an item to the buffer.
     * @param data The data item to write.
     */
    void write(const T& data)
    {
        buffer_[position_++]=data;
    }

    /**
     * @brief Reads a data item from the buffer
     * @param[out] data Reference to where to store the read data.
     */
    void read(T& data)
    {
        data=buffer_[position_++];
    }

    /**
     * @brief Reset the buffer.
     *
     * On return the buffer will be positioned at the start again.
     */
    void reset()
    {
        position_=0;
    }

    /**
     * @brief Test whether the whole buffer was read.
     * @return True if we read or wrot until the end of the buffer.
     */
    bool finished()
    {
        return position_==buffer_.size();
    }

    /**
     * @brief Tests whether the buffer has enough space left to read/write data.
     * @param notItems The number of items to read or write.
     * @return True if there is enough space for noItems items.
     */
    bool hasSpaceForItems(int noItems)
    {
        return position_+noItems<=buffer_.size();
    }
    /**
     * @brief Get the size of the buffer.
     * @return The number of elements the buffer can hold.
     */
    std::size_t size() const
    {
        return buffer_.size();
    }
    /**
     * @brief Converts the buffer to a C array.
     * @return The underlying C array.
     */
    operator T*()
    {
        return buffer_.data();
    }

private:
    /**
     * @brief Pointer to the current insertion point of the buffer.
     */
    std::vector<T> buffer_;
    /**
     * @brief The current position in the buffer.
     */
    std::size_t position_;
};


} // end namespace Detail

/**
 * \brief A communicator where one can do a next neighbour communication
 * in two steps (send, and receive).
 *
 * One can  receive data in some parts of the code and send values at
 * different parts in the code. Note that the send and receive methods will
 * still block until the communciation finished.
 */
class SendReceiveCommunicator
{
public:
    typedef Dune::Interface::InformationMap InformationMap;

    /**
     * \brief Constructor
     * \param communicator The MPI communicator.
     * \param interface The communication interface. \see Interface in dune-common.
     */
    SendReceiveCommunicator(MPI_Comm communicator,
                            const InformationMap& interfaces)
        : interfaces_(&interfaces), communicator_(MPI_COMM_NULL), tag_(424242)
    {
        MPI_Comm_dup(communicator, &communicator_);
        std::size_t max_interface_index_=0;
        for( auto& proc_infos_pair : *interfaces_)
        {
            for(std::size_t i = 0; i < proc_infos_pair.second.first.size(); ++i)
                max_interface_index_ = std::max(max_interface_index_,
                                                proc_infos_pair.second.first[i]);
            for(std::size_t i = 0; i < proc_infos_pair.second.second.size(); ++i)
                max_interface_index_ = std::max(max_interface_index_,
                                                proc_infos_pair.second.second[i]);
        }
        size_handle_vector_size_ = max_interface_index_ + 1;
    }

    ~SendReceiveCommunicator()
    {
        MPI_Comm_free(&communicator_);
    }

    /**
     * \brief Send data according to the interface.
     * \param handle A data handle describing, gathering and scattering
     *               the data.
     * \tparam Datahandle The type of the data handle following the following
     *                    interface:
     template<class Datahandle>
     * \code{.cpp}
     * // returns whether the number of data items per entry is fixed
     * bool fixedsize() const;
     * // get the number of data items for an entry with index i
     * std::size_t size(std::size_t i) const;
     * // gather the data at index i
     * template<class MessageBuffer>
     * void gather(MessageBuffer& buf, std::size_t  i);
     * // scatter the n data items to index i
     * template<class MessageBuffer>
     * void scatter(MessageBuffer& buf, std::size_t i, std::size_t n);
     * \endcode
     */
    template<class Datahandle>
    void sendData(Datahandle& handle);

    /**
     * \brief Receive data according to the interface.
     * \see recvData
     */
    template<class Datahandle>
    void receiveData(Datahandle& handle);

private:
    template<class Datahandle>
    void sendFixedSize(Datahandle& handle);

    template<class Datahandle>
    std::vector<MPI_Request>
    sendFixedSizeWithoutWaiting(Datahandle& handle,
                                std::vector<Detail::MessageBuffer<typename Datahandle::DataType> >& buffers);

    template<class Datahandle, class R>
    void receiveFixedSize(Datahandle& handle,
                          R& postprocess_received_size);

    template<class Datahandle>
    void sendVariableSize(Datahandle& handle);

    template<class Datahandle>
    void receiveVariableSize(Datahandle& handle);

    template<class Datahandle, class R, class T>
    void unpackAndProcessData(Datahandle& handle,
                              std::vector<MPI_Request>& requests,
                              std::vector<Detail::MessageBuffer<typename Datahandle::DataType> >& buffers,
                              R& postprocess_received_size,
                              const T& sizes);

    /**
     * @brief description of the interface.
     *
     * This is a map of the neighboring process number to a pair of local index lists.
     * The first is a list of indices to gather data for sending from and the second is a list of
     * indices to scatter received data to during forward.
     */
    const InformationMap* interfaces_;
    /**
     * @brief The communicator.
     *
     * This is a cloned communicator to ensure there are no interferences.
     */
    MPI_Comm communicator_;

    int tag_;

    /// \brief The vector size needed for the SizeHandle.
    ///
    /// Needs to be indexable with the maximum index in ther interface
    std::size_t size_handle_vector_size_;

    /**
     * \brief Functor for calculating buffer sizes and receive data for variable sized case.
     *
     * In this case we need to receive the size for each data item first.
     * Once on size message is received this functor is used to send the data.
     */
    template<class Type>
    class PostProcessReceivedSizes
    {
    public:
        PostProcessReceivedSizes(const InformationMap& interfaces,
                                 std::vector<MPI_Request>& requests,
                                 std::vector<Detail::MessageBuffer<Type> >& buffers,
                                 const std::vector<std::size_t>& sizes,
                                 int tag, MPI_Comm communicator)
            : interfaces_(interfaces), requests_(requests), buffers_(buffers),
              sizes_(sizes), tag_(tag), communicator_(communicator)
        {}
        void operator()(int index);
    private:
        const InformationMap& interfaces_;
        std::vector<MPI_Request>& requests_;
        std::vector<Detail::MessageBuffer<Type> >& buffers_;
        const std::vector<std::size_t>& sizes_;
        int tag_;
        MPI_Comm communicator_;
    };

    /**
     * \brief A Datahandle for communicating the size per index.
     */
    template<class Datahandle>
    struct SizeDatahandle
    {
        typedef std::size_t DataType;

        SizeDatahandle(const Datahandle& handle)
            : handle_(&handle)
        {}
        SizeDatahandle(std::size_t max_index)
            :sizes_(max_index, std::numeric_limits<std::size_t>::max())
        {}

        bool fixedsize() const
        {
            return true;
        }
        std::size_t size(std::size_t) const
        {
            return 1;
        }
        template<class B>
        void gather(B& buffer, std::size_t i)
        {
            buffer.write(handle_->size(i));
        }
        template<class B>
        void scatter(B& buffer, std::size_t i, std::size_t n)
        {
            DUNE_UNUSED_PARAMETER(n);
            assert(n == 1);
            assert(sizes_[i] == std::numeric_limits<std::size_t>::max());
            buffer.read(sizes_[i]);
        }
        std::vector<std::size_t> sizes_;
        const Datahandle* handle_;
    };

    struct SizeContainer
    {
        explicit SizeContainer(std::size_t size)
            : size_(size)
        {}
        std::size_t operator[](std::size_t) const
        { return size_; }
        const std::size_t size_;
    };
};

template<class Datahandle>
void
SendReceiveCommunicator::sendFixedSize(Datahandle& handle)
{
    std::vector<Detail::MessageBuffer<typename Datahandle::DataType> >
        buffers(interfaces_->size());
    std::vector<MPI_Request> requests =
        sendFixedSizeWithoutWaiting(handle, buffers);
    // We need to wait for the requests to finish.
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

template<class Datahandle>
std::vector<MPI_Request>
SendReceiveCommunicator
::sendFixedSizeWithoutWaiting(Datahandle& handle,
                              std::vector<Detail::MessageBuffer<typename Datahandle::DataType> >& buffers)
{
    std::vector<MPI_Request> requests(interfaces_->size(), MPI_REQUEST_NULL);
    auto current_request = requests.begin();
    auto current_buffer  = buffers.begin();

    for( const auto& infpair : *interfaces_)
    {
        auto& interface_list = infpair.second.first;
        if ( interface_list.size() )
        {
            std::size_t buffer_size = interface_list.size() * handle.size(interface_list[0]);
            current_buffer->resize(buffer_size);
            for(std::size_t i=0; i < interface_list.size(); ++i)
            {
                handle.gather(*current_buffer, interface_list[i]);
            }
            assert(current_buffer->finished());
            MPI_Issend(*current_buffer,
                       interface_list.size() * handle.size(interface_list[0]),
                       Dune::MPITraits<typename Datahandle::DataType>::getType(),
                       infpair.first, tag_, communicator_, &(*current_request));
        }
        ++current_request;
        ++current_buffer;
    }
    return requests;
}

template<class Datahandle>
void
SendReceiveCommunicator::sendVariableSize(Datahandle& handle)
{
    SizeDatahandle<Datahandle> size_handle(handle);
    typedef typename SizeDatahandle<Datahandle>::DataType SizeDataType;
    std::vector<Detail::MessageBuffer<SizeDataType> >
        size_buffers(interfaces_->size());
    std::vector<MPI_Request> size_requests =
        sendFixedSizeWithoutWaiting(size_handle, size_buffers);
    std::vector<MPI_Request> requests(size_requests.size(), MPI_REQUEST_NULL);
    std::vector<Detail::MessageBuffer<typename Datahandle::DataType> > buffers;
    std::vector<int> finished(size_requests.size(), MPI_UNDEFINED);

    while ( true )
    {
        // Receive the size data.
        int no_finished = MPI_UNDEFINED;
        MPI_Waitsome(size_requests.size(), size_requests.data(), &no_finished,
                     finished.data(), MPI_STATUSES_IGNORE);
        if ( no_finished == MPI_UNDEFINED )
        {
            // everything has arrived already
            break;
        }
        for(int i=0; i < no_finished; ++i)
        {
            size_buffers[finished[i]].free(); // free some space.
            // For each size message received send the actual data.
            auto infpair = interfaces_->begin();
            std::advance(infpair, finished[i]);
            auto& interface_list = infpair->second.first;
            assert( interface_list.size() );
            std::size_t buffer_size = 0;
            for(std::size_t j=0; j < interface_list.size(); ++j)
            {
                buffer_size += handle.size(interface_list[j]);
            }

            buffers.emplace_back(buffer_size);
            auto& buffer = buffers.back();
            for(std::size_t j=0; j < interface_list.size(); ++j)
            {
                handle.gather(buffer, interface_list[j]);
            }

            assert(buffer.finished());
            MPI_Issend(buffer, buffer_size,
                       Dune::MPITraits<typename Datahandle::DataType>::getType(),
                       infpair->first, tag_, communicator_, &(requests[finished[i]]));
        }
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}



template<class Datahandle, class R>
void SendReceiveCommunicator::receiveFixedSize(Datahandle& handle,
                                               R& postprocess_received_size)
{
    std::vector<Detail::MessageBuffer<typename Datahandle::DataType> >
        buffers(interfaces_->size());
    std::vector<MPI_Request> requests(interfaces_->size(), MPI_REQUEST_NULL);
    auto current_request = requests.begin();
    auto current_buffer  = buffers.begin();
    std::size_t size = 0;
    for( const auto& infpair : *interfaces_)
    {
        auto& interface_list = infpair.second.second;
        if ( interface_list.size() )
        {
            current_buffer->resize(interface_list.size());
            size = handle.size(interface_list[0]);
            MPI_Irecv(*current_buffer,
                      interface_list.size() * size,
                      Dune::MPITraits<typename Datahandle::DataType>::getType(),
                      infpair.first, tag_, communicator_, &(*current_request));
        }
        ++current_request;
        ++current_buffer;
    }
    unpackAndProcessData(handle, requests, buffers, postprocess_received_size,
                         SizeContainer(size));
}
template<class Datahandle, class R, class T>
void
SendReceiveCommunicator
::unpackAndProcessData(Datahandle& handle,
                       std::vector<MPI_Request>& requests,
                       std::vector<Detail::MessageBuffer<typename Datahandle::DataType> >& buffers,
                       R& postprocess_received_size,
                       const T& sizes)
{
    std::vector<int> finished(requests.size(), MPI_UNDEFINED);

    while( true )
    {
        int no_finished = MPI_UNDEFINED;
        // Wait until the sends have completed
        MPI_Waitsome(requests.size(), requests.data(), &no_finished, finished.data(),
                     MPI_STATUSES_IGNORE);
        if( no_finished == MPI_UNDEFINED )
        {
            // Everything arrived
            break;
        }
        for(int i=0; i < no_finished; ++i)
        {
            // Unpack the received data
            auto& current_buffer = buffers[finished[i]];
            auto infpair = interfaces_->begin();
            std::advance(infpair, finished[i]);
            auto& interface_list = infpair->second.second;
            assert ( interface_list.size() );
            for(std::size_t j=0; j < interface_list.size(); ++j)
            {
                // unpack data
                handle.scatter(current_buffer, interface_list[j],
                               sizes[interface_list[j]]);
            }
            assert(current_buffer.finished());
            current_buffer.free();
            // Postprocess. In the case this was the size of the variable
            // sized data, then we now send the actual data.
            postprocess_received_size(finished[i]);
        }
    }
}

template<class Datahandle>
void
SendReceiveCommunicator::receiveVariableSize(Datahandle& handle)
{
    SizeDatahandle<Datahandle> size_handle(size_handle_vector_size_);
    std::vector<std::size_t>& size_to_receive = size_handle.sizes_;
    std::vector<Detail::MessageBuffer<typename Datahandle::DataType> >
        variable_buffers(interfaces_->size());
    std::vector<MPI_Request> variable_requests(interfaces_->size(), MPI_REQUEST_NULL);
    PostProcessReceivedSizes<typename Datahandle::DataType>
        post_process_received(*interfaces_, variable_requests,
                              variable_buffers, size_to_receive,
                              tag_, communicator_);
    receiveFixedSize(size_handle, post_process_received);

    auto dummy = [] (int) {};
    unpackAndProcessData(handle, variable_requests, variable_buffers, dummy,
                         size_to_receive);
}
template<class Type>
void
SendReceiveCommunicator::PostProcessReceivedSizes<Type>::operator()(int index)
{
    auto infpair = interfaces_.begin();
    std::advance(infpair, index);
    auto& interface_list = infpair->second.second;
    assert ( interface_list.size() );
    // Calculate buffer size
    std::size_t buffer_size = 0;
    for(std::size_t j=0; j < interface_list.size(); ++j)
    {
        buffer_size += sizes_[interface_list[j]];
    }
    auto& buffer = buffers_[index];
    buffer.resize(buffer_size);
    MPI_Irecv(buffer, buffer_size,
              Dune::MPITraits<Type>::getType(),
              infpair->first, tag_, communicator_, &(requests_[index]));
}

template<class Datahandle>
void
SendReceiveCommunicator::sendData(Datahandle& handle)
{
    if( handle.fixedsize() )
    {
        sendFixedSize(handle);
    }
    else
    {
        sendVariableSize(handle); // Takes care of everything including waiting.
    }
}

template<class Datahandle>
void
SendReceiveCommunicator::receiveData(Datahandle& handle)
{
    if( handle.fixedsize() )
    {
        auto null_processor = [](int){};
        receiveFixedSize(handle, null_processor);
    }
    else
    {
        receiveVariableSize(handle); // Takes care of everything including waiting.
    }
}
} // end namespace Opm

#endif
#endif
