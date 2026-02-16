/*
  Copyright 2022-2023 SINTEF AS, 2025 Equinor ASA

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
#ifndef OPM_GPUISTL_GPUAWAREMPISENDER_HPP
#define OPM_GPUISTL_GPUAWAREMPISENDER_HPP

#include <dune/istl/owneroverlapcopy.hh>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSender.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <mpi.h>

#include <memory>
#include <vector>

namespace Opm::gpuistl
{

/**
 * @brief Derived class of GPUSender that handles MPI made with CUDA aware MPI
 * The copOwnerToAll function uses MPI calls refering to data that resides on the GPU in order
 * to send it directly to other GPUs, skipping the staging step on the CPU
 * @tparam field_type is float or double
 * @tparam block_size is the blocksize of the blockelements in the matrix
 * @tparam OwnerOverlapCopyCommunicationType is typically a Dune::LinearOperator::communication_type
 */
template <class field_type, int block_size, class OwnerOverlapCopyCommunicationType>
class GPUAwareMPISender : public GPUSender<field_type, OwnerOverlapCopyCommunicationType>
{
public:
    using X = GpuVector<field_type>;

    explicit GPUAwareMPISender(const OwnerOverlapCopyCommunicationType& cpuOwnerOverlapCopy)
        : GPUSender<field_type, OwnerOverlapCopyCommunicationType>(cpuOwnerOverlapCopy)
    {
    }

    void copyOwnerToAll(const X& source, X& dest) const override
    {
        OPM_ERROR_IF(&source != &dest,
                     "The provided GpuVectors' address did not match"); // In this context, source == dest!!!
        std::call_once(this->m_initializedIndices, [&]() { initIndexSet(); });

        int rank = this->m_cpuOwnerOverlapCopy.communicator().rank();
        dest.prepareSendBuf(*m_GPUSendBuf, *m_commpairIndicesOwner);

        // Start MPI stuff here...
        // Note: This has been taken from DUNE's parallel/communicator.hh
        std::vector<MPI_Request> sendRequests(m_messageInformation.size());
        std::vector<MPI_Request> recvRequests(m_messageInformation.size());
        std::vector<int> processMap(m_messageInformation.size());
        size_t numberOfRealRecvRequests = 0;

        using const_iterator = typename InformationMap::const_iterator;
        const const_iterator end = m_messageInformation.end();

        {
            size_t i = 0;
            for (const_iterator info = m_messageInformation.begin(); info != end; ++info, ++i) {
                processMap[i] = info->first;
                if (info->second.second.m_size) {
                    MPI_Irecv(m_GPURecvBuf->data() + info->second.second.m_start,
                              detail::to_int(info->second.second.m_size),
                              MPI_BYTE,
                              info->first,
                              m_commTag,
                              this->m_cpuOwnerOverlapCopy.communicator(),
                              &recvRequests[i]);
                    numberOfRealRecvRequests += 1;
                } else {
                    recvRequests[i] = MPI_REQUEST_NULL;
                }
            }
        }

        {
            size_t i = 0;
            for (const_iterator info = m_messageInformation.begin(); info != end; ++info, ++i) {
                if (info->second.first.m_size) {
                    MPI_Issend(m_GPUSendBuf->data() + info->second.first.m_start,
                               detail::to_int(info->second.first.m_size),
                               MPI_BYTE,
                               info->first,
                               m_commTag,
                               this->m_cpuOwnerOverlapCopy.communicator(),
                               &sendRequests[i]);
                } else {
                    sendRequests[i] = MPI_REQUEST_NULL;
                }
            }
        }
        int finished = MPI_UNDEFINED;
        MPI_Status status;
        for (size_t i = 0; i < numberOfRealRecvRequests; i++) {
            status.MPI_ERROR = MPI_SUCCESS;
            MPI_Waitany(m_messageInformation.size(), recvRequests.data(), &finished, &status);

            if (status.MPI_ERROR != MPI_SUCCESS) {
                OPM_THROW(std::runtime_error,
                          fmt::format(fmt::runtime("MPI_Error occurred while rank {} received a message from rank {}"),
                                      rank,
                                      processMap[finished]));
            }
        }
        MPI_Status recvStatus;
        for (size_t i = 0; i < m_messageInformation.size(); i++) {
            if (MPI_SUCCESS != MPI_Wait(&sendRequests[i], &recvStatus)) {
                OPM_THROW(std::runtime_error,
                          fmt::format(fmt::runtime("MPI_Error occurred while rank {} sent a message from rank {}"),
                                      rank,
                                      processMap[finished]));
            }
        }
        // ...End of MPI stuff

        dest.syncFromRecvBuf(*m_GPURecvBuf, *m_commpairIndicesCopy);
    }

private:
    mutable std::unique_ptr<GpuVector<int>> m_commpairIndicesCopy;
    mutable std::unique_ptr<GpuVector<int>> m_commpairIndicesOwner;
    mutable std::unique_ptr<GpuVector<field_type>> m_GPUSendBuf;
    mutable std::unique_ptr<GpuVector<field_type>> m_GPURecvBuf;

    struct MessageInformation {
        MessageInformation()
            : m_start(0)
            , m_size(0)
        {
        }
        MessageInformation(size_t start, size_t size)
            : m_start(start)
            , m_size(size)
        {
        }
        size_t m_start; // offset in elements of "field_type"
        size_t m_size; // size in bytes
    };

    using InformationMap = std::map<int, std::pair<MessageInformation, MessageInformation>>;
    mutable InformationMap m_messageInformation;
    using IM = std::map<int, std::pair<std::vector<int>, std::vector<int>>>;
    mutable IM m_im;

    constexpr static int m_commTag = 0; // So says DUNE

    void buildCommPairIdxs() const
    {
        const auto& ri = this->m_cpuOwnerOverlapCopy.remoteIndices();
        std::vector<int> commpairIndicesCopyOnCPU;
        std::vector<int> commpairIndicesOwnerCPU;

        for (const auto& process : ri) {
            m_im[process.first] = std::pair(std::vector<int>(), std::vector<int>());
            for (int send = 0; send < 2; ++send) {
                auto remoteEnd = send ? process.second.first->end() : process.second.second->end();
                auto remote = send ? process.second.first->begin() : process.second.second->begin();

                while (remote != remoteEnd) {
                    if (send ? (remote->localIndexPair().local().attribute() == 1) : (remote->attribute() == 1)) {
                        if (send) {
                            m_im[process.first].first.push_back(remote->localIndexPair().local().local());
                        } else {
                            m_im[process.first].second.push_back(remote->localIndexPair().local().local());
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
                m_messageInformation.insert(std::make_pair(
                    it->first,
                    std::make_pair(
                        MessageInformation(sendBufIdx * block_size, noSend * block_size * sizeof(field_type)),
                        MessageInformation(recvBufIdx * block_size, noRecv * block_size * sizeof(field_type)))));

                for (int x = 0; x < noSend; x++) {
                    for (int bs = 0; bs < block_size; bs++) {
                        commpairIndicesOwnerCPU.push_back(it->second.first[x] * block_size + bs);
                    }
                }
                for (int x = 0; x < noRecv; x++) {
                    for (int bs = 0; bs < block_size; bs++) {
                        commpairIndicesCopyOnCPU.push_back(it->second.second[x] * block_size + bs);
                    }
                }
                sendBufIdx += noSend;
                recvBufIdx += noRecv;
            }
        }

        m_commpairIndicesCopy = std::make_unique<GpuVector<int>>(commpairIndicesCopyOnCPU);
        m_commpairIndicesOwner = std::make_unique<GpuVector<int>>(commpairIndicesOwnerCPU);

        m_GPUSendBuf = std::make_unique<GpuVector<field_type>>(sendBufIdx * block_size);
        m_GPURecvBuf = std::make_unique<GpuVector<field_type>>(recvBufIdx * block_size);
    }

    void initIndexSet() const override
    {
        // We need indices that we we will use in the project, dot and norm calls.
        // TODO: [premature perf] Can this be run once per instance? Or do we need to rebuild every time?
        const auto& pis = this->m_cpuOwnerOverlapCopy.indexSet();
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

        this->m_indicesCopy = std::make_unique<GpuVector<int>>(indicesCopyOnCPU);
        this->m_indicesOwner = std::make_unique<GpuVector<int>>(indicesOwnerCPU);

        buildCommPairIdxs();
    }
};

} // namespace Opm::gpuistl
#endif // OPM_GPUISTL_GPUAWAREMPISENDER_HPP
