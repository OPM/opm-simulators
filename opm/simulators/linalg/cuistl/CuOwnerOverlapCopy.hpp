/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_CUISTL_CUOWNEROVERLAPCOPY_HPP
#define OPM_CUISTL_CUOWNEROVERLAPCOPY_HPP
#include <dune/istl/owneroverlapcopy.hh>
#include <memory>
#include <mutex>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>

namespace Opm::cuistl
{
/**
 * @brief CUDA compatiable variant of Dune::OwnerOverlapCopyCommunication
 *
 * This class can essentially be seen as an adapter around Dune::OwnerOverlapCopyCommunication, and should work as
 * a Dune::OwnerOverlapCopyCommunication on CuVectors
 *
 *
 * @note This currently only has the functionality to parallelize the linear solve.
 *
 * @tparam field_type should be a field_type supported by CuVector (double, float)
 * @tparam block_size the block size used (this is relevant for say figuring out the correct indices)
 * @tparam OwnerOverlapCopyCommunicationType should mimic Dune::OwnerOverlapCopyCommunication.
 */
template <class field_type, int block_size, class OwnerOverlapCopyCommunicationType>
class CuOwnerOverlapCopy
{
public:
    using X = CuVector<field_type>;

    CuOwnerOverlapCopy(const OwnerOverlapCopyCommunicationType& cpuOwnerOverlapCopy)
        : m_cpuOwnerOverlapCopy(cpuOwnerOverlapCopy)
    {
    }
    /**
     * @brief dot will carry out the dot product between x and y on the owned indices, then sum up the result across MPI
     * processes.
     * @param[out] output result will be stored here
     *
     * @note This uses the same interface as its DUNE equivalent.
     */
    void dot(const X& x, const X& y, field_type& output) const
    {
        std::call_once(m_initializedIndices, [&]() { initIndexSet(); });

        const auto dotAtRank = x.dot(y, *m_indicesOwner);
        output = m_cpuOwnerOverlapCopy.communicator().sum(dotAtRank);
    }

    /**
     * @brief norm computes the l^2-norm of x across processes.
     *
     * This will compute the dot product of x with itself on owned indices, then
     * sum the result across process and return the square root of the sum.
     */
    field_type norm(const X& x) const
    {
        auto xDotX = field_type(0);
        this->dot(x, x, xDotX);

        using std::sqrt;
        return sqrt(xDotX);
    }

    /**
     * @brief project will project x to the owned subspace
     *
     * For each component i which is not owned, x_i will be set to 0
     *
     * @param[inout] x the vector to project
     */
    void project(X& x) const
    {
        std::call_once(m_initializedIndices, [&]() { initIndexSet(); });
        x.setZeroAtIndexSet(*m_indicesCopy);
    }

    /**
     * @brief copyOwnerToAll will copy source to the CPU, then call OwnerOverlapCopyCommunicationType::copyOwnerToAll on
     * the copied data, and copy the result back to the GPU
     * @param[in] source
     * @param[out] dest
     */
    void copyOwnerToAll_orig(const X& source, X& dest) const
    {
        // TODO: [perf] Can we reduce copying from the GPU here?
        // TODO: [perf] Maybe create a global buffer instead?
        auto sourceAsDuneVector = source.template asDuneBlockVector<block_size>();
        auto destAsDuneVector = dest.template asDuneBlockVector<block_size>();
        m_cpuOwnerOverlapCopy.copyOwnerToAll(sourceAsDuneVector, destAsDuneVector);
        dest.copyFromHost(destAsDuneVector);
    }

    // Georgs new code intended to use GPU direct
    void copyOwnerToAll(const X& source, X& dest) const
    {

        printf("\n\nGPU DIRECT CODE IS RUN\n\n");
            printf("Compile time check:\n");
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
    printf("This MPI library has CUDA-aware support.\n", MPIX_CUDA_AWARE_SUPPORT);
#elif defined(MPIX_CUDA_AWARE_SUPPORT) && !MPIX_CUDA_AWARE_SUPPORT
    printf("This MPI library does not have CUDA-aware support.\n");
#else
    printf("This MPI library cannot determine if there is CUDA-aware support.\n");
#endif /* MPIX_CUDA_AWARE_SUPPORT */
 
    printf("Run time check:\n");
#if defined(MPIX_CUDA_AWARE_SUPPORT)
    if (1 == MPIX_Query_cuda_support()) {
        printf("This MPI library has CUDA-aware support.\n");
    } else {
        printf("This MPI library does not have CUDA-aware support.\n");
    }
#else /* !defined(MPIX_CUDA_AWARE_SUPPORT) */
    printf("This MPI library cannot determine if there is CUDA-aware support.\n");
#endif /* MPIX_CUDA_AWARE_SUPPORT */


        assert(&source == &dest); // In this context, source == dest!!!
        std::call_once(m_initializedIndices, [&]() { initIndexSet(); });

        int rank;
        MPI_Comm_rank(m_cpuOwnerOverlapCopy.communicator(), &rank);
        dest.prepareSendBuf(*m_GPUSendBuf, *m_commpair_indicesOwner);

        // Start MPI stuff here...
        // Note: This has been taken from DUNE's parallel/communicator.hh
        MPI_Request* sendRequests = new MPI_Request[messageInformation_.size()];
        MPI_Request* recvRequests = new MPI_Request[messageInformation_.size()];
        size_t numberOfRealRecvRequests = 0;

        typedef typename InformationMap::const_iterator const_iterator;
        const const_iterator end = messageInformation_.end();
        size_t i=0;
        int* processMap = new int[messageInformation_.size()];
        for(const_iterator info = messageInformation_.begin(); info != end; ++info, ++i) {
            processMap[i]=info->first;
            if(info->second.second.size_) {
                MPI_Irecv(m_GPURecvBuf->data()+info->second.second.start_,
                          info->second.second.size_,
                          MPI_BYTE,
                          info->first,
                          commTag_,
                          m_cpuOwnerOverlapCopy.communicator(),
                          &recvRequests[i]);
                numberOfRealRecvRequests += 1;
            } else {
                recvRequests[i]=MPI_REQUEST_NULL;
            }
        }

        i=0;
        for(const_iterator info = messageInformation_.begin(); info != end; ++info, ++i) {
            if(info->second.first.size_) {
                MPI_Issend(m_GPUSendBuf->data()+info->second.first.start_,
                           info->second.first.size_,
                           MPI_BYTE,
                           info->first,
                           commTag_,
                           m_cpuOwnerOverlapCopy.communicator(),
                           &sendRequests[i]);
            } else {
                sendRequests[i]=MPI_REQUEST_NULL;
            }
        }
        i=0;
        int finished = MPI_UNDEFINED;
        MPI_Status status;
        for(i=0; i< numberOfRealRecvRequests; i++) {
            status.MPI_ERROR=MPI_SUCCESS;
            MPI_Waitany(messageInformation_.size(), recvRequests, &finished, &status);

            if(status.MPI_ERROR!=MPI_SUCCESS) {
                std::cerr<< rank << ": MPI_Error occurred while receiving message from "<< processMap[finished] << std::endl;
                assert(false);
            }
        }
        MPI_Status recvStatus;
        for(i=0; i< messageInformation_.size(); i++) {
            if(MPI_SUCCESS!=MPI_Wait(&sendRequests[i], &recvStatus)) {
                std::cerr << rank << ": MPI_Error occurred while sending message to " << processMap[finished] << std::endl;
                assert(false);
            }
        }
        delete[] processMap;
        delete[] sendRequests;
        delete[] recvRequests;
        // ...End of MPI stuff

        dest.syncFromRecvBuf(*m_GPURecvBuf, *m_commpair_indicesCopy);
    }



private:
    const OwnerOverlapCopyCommunicationType& m_cpuOwnerOverlapCopy;

    // Used to call the initIndexSet. Note that this is kind of a
    // premature optimization, in the sense that we could just initialize these indices
    // always, but they are not always used.
    mutable std::once_flag m_initializedIndices;
    mutable std::unique_ptr<CuVector<int>> m_indicesCopy;
    mutable std::unique_ptr<CuVector<int>> m_indicesOwner;

        mutable std::unique_ptr<CuVector<int>> m_commpair_indicesCopy;
    mutable std::unique_ptr<CuVector<int>> m_commpair_indicesOwner;
    mutable std::unique_ptr<CuVector<field_type>> m_GPUSendBuf;
    mutable std::unique_ptr<CuVector<field_type>> m_GPURecvBuf;

    struct MessageInformation
    {
        MessageInformation() : start_(0), size_(0) {}
        MessageInformation(size_t start, size_t size) : start_(start), size_(size) {}
        size_t start_; // offset in elements of "field_type"
        size_t size_;  // size in bytes
    };

    typedef std::map<int,std::pair<MessageInformation,MessageInformation> > InformationMap;
    mutable InformationMap messageInformation_;
    typedef std::map<int,std::pair<std::vector<int>,std::vector<int> > > IM;
    mutable IM m_im;

    constexpr static int commTag_ = 0; // So says DUNE

    void buildCommPairIdxs() const
    {
        int rank;
        MPI_Comm_rank(m_cpuOwnerOverlapCopy.communicator(), &rank);
        auto &ri = m_cpuOwnerOverlapCopy.remoteIndices();
        auto end = ri.end();
        std::vector<int> commpair_indicesCopyOnCPU;
        std::vector<int> commpair_indicesOwnerCPU;

        for(auto process = ri.begin(); process != end; ++process) {
            int size = 0;
            m_im[process->first] = std::pair(std::vector<int>(), std::vector<int>());
            for(int send = 0; send < 2; ++send) {
                auto remoteEnd = send ? process->second.first->end()
                                      : process->second.second->end();
                auto remote = send ? process->second.first->begin()
                                   : process->second.second->begin();

                while(remote != remoteEnd) {
                    if (send ? (remote->localIndexPair().local().attribute() == 1)
                             : (remote->attribute() == 1)) {
                        ++size;
                        if (send) {
                            m_im[process->first].first.push_back(remote->localIndexPair().local().local()); 
                        } else {
                            m_im[process->first].second.push_back(remote->localIndexPair().local().local());
                        }
                    }
                    ++remote;
                }
            }
        }

        int sendBufIdx = 0;
        int recvBufIdx = 0;
        for (auto it = m_im.begin(); it != m_im.end(); it++) {
            int noSend = it->second.first.size();
            int noRecv = it->second.second.size();

            if (noSend + noRecv > 0) {
                messageInformation_.insert(
                        std::make_pair(it->first,
                                       std::make_pair(MessageInformation(
                                                        sendBufIdx * block_size,
                                                        noSend * block_size * sizeof(field_type)),
                                                      MessageInformation(
                                                        recvBufIdx * block_size,
                                                        noRecv * block_size * sizeof(field_type)))));

                for(int x = 0; x < noSend; x++) {
                    for(int bs = 0; bs < block_size; bs++) {
                        commpair_indicesOwnerCPU.push_back(it->second.first[x] * block_size + bs);
                    }
                }
                for(int x = 0; x < noRecv; x++) {
                    for(int bs = 0; bs < block_size; bs++) {
                        commpair_indicesCopyOnCPU.push_back(it->second.second[x] * block_size + bs);
                    }
                }
                sendBufIdx += noSend;
                recvBufIdx += noRecv;
            }
        }

        m_commpair_indicesCopy = std::make_unique<CuVector<int>>(commpair_indicesCopyOnCPU);
        m_commpair_indicesOwner = std::make_unique<CuVector<int>>(commpair_indicesOwnerCPU);

        m_GPUSendBuf = std::make_unique<CuVector<field_type>>(sendBufIdx * block_size);
        m_GPURecvBuf = std::make_unique<CuVector<field_type>>(recvBufIdx * block_size);
    }

    void initIndexSet() const
    {
        // We need indices that we we will use in the project, dot and norm calls.
        // TODO: [premature perf] Can this be run once per instance? Or do we need to rebuild every time?
        const auto& pis = m_cpuOwnerOverlapCopy.indexSet();
        std::vector<int> indicesCopyOnCPU;
        std::vector<int> indicesOwnerCPU;
        for (const auto& index : pis) {
            if (index.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::copy) {
                for (int component = 0; component < block_size; ++component) {
                    indicesCopyOnCPU.push_back(index.local().local() * block_size + component);
                }
            }

            if (index.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                for (int component = 0; component < block_size; ++component) {
                    indicesOwnerCPU.push_back(index.local().local() * block_size + component);
                }
            }
        }

        m_indicesCopy = std::make_unique<CuVector<int>>(indicesCopyOnCPU);
        m_indicesOwner = std::make_unique<CuVector<int>>(indicesOwnerCPU);

        buildCommPairIdxs();
    }
};
} // namespace Opm::cuistl
#endif
