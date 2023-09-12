/*
  Copyright 2023 Equinor ASA

  This file is part of the Open Porous Media Project (OPM).

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

#include <config.h>

#include <opm/simulators/utils/ParallelNLDDPartitioningZoltan.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/CSRGraphFromCoordinates.hpp>

// Note: The build system guarantees that we're only built in configurations
// that have both MPI and Zoltan.  Zoltan.h redefines 'HAVE_MPI', so let it
// do so, but restore a simple definition afterwards.
#undef HAVE_MPI
#include <zoltan.h>
#undef HAVE_MPI
#define HAVE_MPI

#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

    /// Assign Zoltan control parameters.
    ///
    /// Built-in defaults possibly overriden by user-selected values.
    ///
    /// \param[in] params User-selected Zoltan control parameters.
    ///
    /// \param[in,out] zz Opaque Zoltan context.
    void setZoltanParameters(const Opm::ParallelNLDDPartitioningZoltan::ZoltanParamMap& params,
                             Zoltan_Struct*                                             zz)
    {
        Zoltan_Set_Param(zz, "DEBUG_LEVEL",             "0");
        Zoltan_Set_Param(zz, "LB_METHOD",               "GRAPH");
        Zoltan_Set_Param(zz, "LB_APPROACH",             "PARTITION");
        Zoltan_Set_Param(zz, "NUM_GID_ENTRIES",         "1");
        Zoltan_Set_Param(zz, "NUM_LID_ENTRIES",         "1");
        Zoltan_Set_Param(zz, "RETURN_LISTS",            "EXPORT"); // Self to "other"
        Zoltan_Set_Param(zz, "CHECK_GRAPH",             "2");
        Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM",         "0");
        Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM",          "0");
        Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", "0.35"); // 0-remove all, 1-remove none

        for (const auto& [ param, value ] : params) {
            Zoltan_Set_Param(zz, param.c_str(), value.c_str());
        }
    }

    /// Create consecutive numbering of reachable vertices
    class EnumerateSeenVertices
    {
    public:
        /// Constructor.
        ///
        /// \tparam Edge Type for representing a directed edge between a
        ///   pair of vertices.  Must support \code std::get<>() \endcode
        ///   protocol.  Typically a \code std::pair<> \endcode of integral
        ///   types.
        ///
        /// \param[in] numVertices Maximum number of vertices in graph.
        ///
        /// \param[in] edges Edge representation of connectivity graph.
        template <typename Edge>
        explicit EnumerateSeenVertices(const std::size_t numVertices,
                                       const std::vector<Edge>& edges)
            : index_(numVertices, -1)
        {
            auto seen = std::vector<bool>(numVertices, false);

            for (const auto& [v1, v2] : edges) {
                seen[v1] = true;
                seen[v2] = true;
            }

            for (auto vertex = 0*numVertices; vertex < numVertices; ++vertex) {
                if (seen[vertex]) {
                    this->index_[vertex] = this->numSeenVertices_++;
                }
            }
        }

        /// Retrieve number of reachable vertices.  Less than or equal to \p
        /// numVertices.
        std::size_t numVertices() const
        {
            return this->numSeenVertices_;
        }

        /// Retrieve reachable index of vertex.
        ///
        /// \param[in] vertex Vertex ID in original numbering
        ///
        /// \return Reachable index of \p vertex.  -1 if \p vertex is not
        ///   reachable.
        int operator[](const std::size_t vertex) const
        {
            return this->index_[vertex];
        }

    private:
        /// Indices of reachable vertices.
        std::vector<int> index_{};

        /// Number of reachable vertices.  One more than the maximum value
        /// in \c index_.
        int numSeenVertices_{0};
    };

    /// Unstructured graph of reachable vertices
    class VertexGraph
    {
    public:
        /// Constructor.
        ///
        /// \tparam Edge Type for representing a directed edge between a
        ///   pair of vertices.  Must support \code std::get<>() \endcode
        ///   protocol.  Typically a \code std::pair<> \endcode of integral
        ///   types.
        ///
        /// \tparam GlobalCellID Callback type for mapping a vertex ID to a
        ///   globally unique ID.  Typically the same as
        ///   \code ParallelNLDDPartitioningZoltan::GlobalCellID \endcode.
        ///
        /// \param[in] myRank Current rank in MPI communicator.  Needed by
        ///   Zoltan.
        ///
        /// \param[in] numVertices Maximum number of vertices in graph.
        ///
        /// \param[in] edges Edge representation of connectivity graph.
        ///
        /// \param[in] vertexId Enumeration of reachable vertices.
        ///
        /// \param[in] globalCell Callback for mapping (local) vertex IDs to
        ///   globally unique vertex IDs.
        template <typename Edge, typename GlobalCellID>
        explicit VertexGraph(const int                    myRank,
                             const std::size_t            numVertices,
                             const std::vector<Edge>&     edges,
                             const EnumerateSeenVertices& vertexId,
                             GlobalCellID&&               globalCell)
            : myRank_ { myRank }
        {
            // Form undirected connectivity graph.
            for (const auto& [v1, v2] : edges) {
                this->graph_.addConnection(vertexId[v1], vertexId[v2]);
                this->graph_.addConnection(vertexId[v2], vertexId[v1]);
            }

            this->graph_.compress(vertexId.numVertices());

            // Form local-to-global vertex ID mapping for reachable vertices.
            this->globalCell_.resize(vertexId.numVertices());
            for (auto vertex = 0*numVertices; vertex < numVertices; ++vertex) {
                if (const auto localIx = vertexId[vertex]; localIx >= 0) {
                    this->globalCell_[localIx] = globalCell(vertex);
                }
            }
        }

        /// Retrive my rank in current MPI communicator.
        ///
        /// Needed by Zoltan.
        int myRank() const
        {
            return this->myRank_;
        }

        /// Retrieve number of vertices in connectivity graph.
        int numVertices() const
        {
            return static_cast<int>(this->graph_.numVertices());
        }

        /// Retrieve globally unqique ID of reachable vertex.
        ///
        /// \param[in] localCell Index of locally reachable cell/vertex.
        int globalId(const int localCell) const
        {
            return this->globalCell_[localCell];
        }

        /// Get read-only access to start pointers of graph's CSR representation.
        ///
        /// \return Reference to start pointers (IA array).
        decltype(auto) startPointers() const
        {
            return this->graph_.startPointers();
        }

        /// Get read-only access to column indices of graph's CSR representation.
        ///
        /// \return Reference to column indices (JA array).
        decltype(auto) columnIndices() const
        {
            return this->graph_.columnIndices();
        }

    private:
        // VertexID = int, TrackCompressedIdx = false
        using Backend = Opm::utility::CSRGraphFromCoordinates<>;

        /// My rank in current MPI communicator.
        int myRank_{};

        /// Globally unique vertex ID of each locally reachable vertex.
        std::vector<int> globalCell_{};

        /// Vertex connectivity graph.
        Backend graph_{};
    };

// Use C linkage for Zoltan interface/query functions.  Ensures maximum compatibility.
extern "C" {
    /// Compute number of vertices in connectivity graph.
    ///
    /// Callback for Zoltan_Set_Num_Obj_Fn.
    ///
    /// \param[in] graphPtr Opaque user data.  Assumed to point to a
    ///   \c VertexGraph instance.
    ///
    /// \param[out] ierr Error code for Zoltan consumption.  Single \c int.
    int numVertices(void* graphPtr, int* ierr)
    {
        *ierr = ZOLTAN_OK;
        return static_cast<VertexGraph*>(graphPtr)->numVertices();
    }

    /// Compute number of neighbours for each vertex in connectivity graph.
    ///
    /// Callback for Zoltan_Set_Num_Edges_Multi_Fn.
    ///
    /// \param[in] graphPtr Opaque user data.  Assumed to point to a
    ///   \c VertexGraph instance.
    ///
    /// \param[in] numCells Number of cells/vertices/objects for which to
    ///   compute number of neighbours.
    ///
    /// \param[in] localID Local IDs of those cells/vertices/objects for
    ///   which to compute the respective number of neighbours.
    ///
    /// \param[in,out] numEdges Number of neighbours (in/out edges) for each
    ///   cell/vertex/object.  Size equal to \code numVertices(graphPtr)
    ///   \endcode.  Populated by this function.  Allocated by Zoltan.
    ///
    /// \param[out] ierr Error code for Zoltan consumption.  Single \c int.
    void numEdges(void*            graphPtr,
                  const int     /* sizeGID */,
                  const int     /* sizeLID */,
                  const int        numCells,
                  ZOLTAN_ID_PTR /* globalID */,
                  ZOLTAN_ID_PTR    localID,
                  int*             numEdges,
                  int*             ierr)
    {
        const auto& ia = static_cast<VertexGraph*>(graphPtr)->startPointers();

        for (auto cell = 0*numCells; cell < numCells; ++cell) {
            const auto localCell = localID[cell];
            numEdges[cell] = ia[localCell + 1] - ia[localCell + 0];
        }

        *ierr = ZOLTAN_OK;
    }

    /// Identify objects/vertices in connectivity graph.
    ///
    /// Callback for Zoltan_Set_Obj_List_Fn.
    ///
    /// \param[in] graphPtr Opaque user data.  Assumed to point to a
    ///   \c VertexGraph instance.
    ///
    /// \param[in] numElmsPerGid Number of ZOLTAN_ID_TYPE elements needed to
    ///   describe a single globally unique object/vertex ID.  Must be one
    ///   (1) in this implementation.
    ///
    /// \param[in] numElmsPerLid Number of ZOLTAN_ID_TYPE elements needed to
    ///   describe a single local object/vertex ID.  Must be one (1) in this
    ///   implementation.
    ///
    /// \param[in,out] globalIds Globally unique object/vertex IDs.  Size
    ///   equal to \code numElmsPerGid * numVertices(graphPtr) \endcode.
    ///   Populated by this function.  Allocated by Zoltan.
    ///
    /// \param[in,out] localIds Local object/vertex IDs.  Size equal to
    ///   \code numElmsPerLid * numVertices(graphPtr) \endcode.  Populated
    ///   by this function.  Allocated by Zoltan.
    ///
    /// \param[out] ierr Error code for Zoltan consumption.  Single \c int.
    void vertexList(void*            graphPtr,
                    const int        numElmsPerGid,
                    const int        numElmsPerLid,
                    ZOLTAN_ID_PTR    globalIds,
                    ZOLTAN_ID_PTR    localIds,
                    const int     /* wgtDim */,
                    float*        /* objWgts */,
                    int*             ierr)
    {
        if ((numElmsPerGid != numElmsPerLid) || (numElmsPerLid != 1)) {
            *ierr = ZOLTAN_FATAL;
        }

        const auto* graph = static_cast<VertexGraph*>(graphPtr);

        if (graph->numVertices() == 0) {
            *ierr = ZOLTAN_OK;
            return;
        }

        std::iota(&localIds[0], &localIds[0] + graph->numVertices(), ZOLTAN_ID_TYPE{0});
        std::transform(&localIds[0],
                       &localIds[0] + graph->numVertices(),
                       &globalIds[0],
                       [graph](const int localCell) -> ZOLTAN_ID_TYPE
                       {
                           return graph->globalId(localCell);
                       });

        *ierr = ZOLTAN_OK;
    }

    /// Identify neighbours in connectivity graph.
    ///
    /// Callback for Zoltan_Set_Edge_List_Multi_Fn
    ///
    /// \param[in] graphPtr Opaque user data.  Assumed to point to a
    ///   \c VertexGraph instance.
    ///
    /// \param[in] numElmsPerGid Number of ZOLTAN_ID_TYPE elements needed to
    ///   describe a single globally unique object/vertex ID.  Must be one
    ///   (1) in this implementation.
    ///
    /// \param[in] numElmsPerLid Number of ZOLTAN_ID_TYPE elements needed to
    ///   describe a single local object/vertex ID.  Must be one (1) in this
    ///   implementation.
    ///
    /// \param[in] numCells Number of cells/vertices/objects for which to
    ///   identify the neighbouring cells/vertices/objects.
    ///
    /// \param[in] localIds Local object/vertex IDs.  Size equal to \code
    ///   numElmsPerLid * numVertices(graphPtr) \endcode.
    ///
    /// \param[in,out] neighbourGid Globally unique object ID of each
    ///   neighbour cell/vertex/object.  Populated by this function.
    ///   Allocated by Zoltan.
    ///
    /// \param[in,out] neighbourProc Owner of each neighbouring
    ///   cell/vertex/object.  Populated by this function.  Allocated by
    ///   Zoltan.
    ///
    /// \param[out] ierr Error code for Zoltan consumption.  Single \c int.
    void edgeList(void*            graphPtr,
                  const int        numElmsPerGid,
                  const int        numElmsPerLid,
                  const int        numCells,
                  ZOLTAN_ID_PTR /* globalIds */,
                  ZOLTAN_ID_PTR    localIds,
                  int*          /* numEdges */,
                  ZOLTAN_ID_PTR    neighbourGid,
                  int*             neighbourProc,
                  int           /* weightDim */,
                  float*        /* edgeWeights */,
                  int*             ierr)
    {
        const auto* graph = static_cast<VertexGraph*>(graphPtr);

        if ((numElmsPerGid != numElmsPerLid) ||
            (numElmsPerLid != 1) ||
            (numCells != graph->numVertices()))
        {
            *ierr = ZOLTAN_FATAL;
        }

        if (graph->numVertices() == 0) {
            *ierr = ZOLTAN_OK;
            return;
        }

        const auto& ia = graph->startPointers();
        const auto& ja = graph->columnIndices();

        for (auto cell = 0*numCells, ix = 0; cell < numCells; ++cell) {
            const auto localCell = localIds[cell];
            for (auto neighIx = ia[localCell], end = ia[localCell + 1];
                 neighIx != end; ++neighIx, ++ix)
            {
                neighbourGid [ix] = static_cast<ZOLTAN_ID_TYPE>(graph->globalId(ja[neighIx]));
                neighbourProc[ix] = graph->myRank();
            }
        }
    }
} // extern "C"

    /// Partition VertexGraph objects using Zoltan graph partitioning software.
    class Partitioner
    {
    public:
        /// Constructor.
        ///
        /// \param[in] comm MPI communicator.
        ///
        /// \param[in] params Control parameters for Zoltan graph partiitioning procedure.
        explicit Partitioner(const Opm::Parallel::Communication                         comm,
                             const Opm::ParallelNLDDPartitioningZoltan::ZoltanParamMap& params);

        /// Destructor.
        ///
        /// Destroys Zoltan context.
        ~Partitioner();

        /// Partition VertexGraph instance.
        ///
        /// \param[in] graphPtr Graph being partitioned.  Assumed to point
        ///   to a \c VertexGraph instance that is fully populated by caller.
        ///
        /// \param[in] numCells Number of vertices in graph pointed to by \p
        ///   graphPtr.  Used to size the resulting partition vector.
        ///
        /// \return Partition vector.  Non-negative block IDs for each of
        ///   the \p numCells cells/vertices/objects in \p graphPtr.
        std::vector<int> operator()(void* graphPtr, const std::size_t numCells);

    private:
        /// Helper RAII type to manage partition ("part") memory allocated
        /// by Zoltan.
        struct ZoltanPart
        {
            /// Globally unique cell/vertex/object IDs.
            ZOLTAN_ID_PTR globalId{nullptr};

            /// Local cell/vertex/object IDs.
            ZOLTAN_ID_PTR localId{nullptr};

            /// Owning process for each cell/vertex/object.
            int* process{nullptr};

            /// Block/domain to which each cell/vertex/object is assigned.
            int* block{nullptr};

            /// Destructor.
            ///
            /// Releases partition ("part") memory acquired by Zoltan.
            ~ZoltanPart();
        };

        /// Zoltan library context.
        Zoltan_Struct* zz_{nullptr};
    };

    Partitioner::ZoltanPart::~ZoltanPart()
    {
        Zoltan_LB_Free_Part(&this->globalId,
                            &this->localId,
                            &this->process,
                            &this->block);
    }

    Partitioner::Partitioner(const Opm::Parallel::Communication                         comm,
                             const Opm::ParallelNLDDPartitioningZoltan::ZoltanParamMap& params)
    {
        const auto argc   = 0;
        char*      argv[] = { nullptr };
        auto       ver    = 0.0f;

        const auto rc = Zoltan_Initialize(argc, argv, &ver);
        if (rc != ZOLTAN_OK) {
            OPM_THROW(std::runtime_error, "Unable to Initialise Zoltan");
        }

        this->zz_ = Zoltan_Create(comm);
        setZoltanParameters(params, this->zz_);
    }

    Partitioner::~Partitioner()
    {
        Zoltan_Destroy(&this->zz_);
    }

    std::vector<int>
    Partitioner::operator()(void* graphPtr, const std::size_t numCells)
    {
        Zoltan_Set_Num_Obj_Fn        (this->zz_, &numVertices, graphPtr);
        Zoltan_Set_Num_Edges_Multi_Fn(this->zz_, &numEdges   , graphPtr);

        Zoltan_Set_Obj_List_Fn       (this->zz_, &vertexList, graphPtr);
        Zoltan_Set_Edge_List_Multi_Fn(this->zz_, &edgeList  , graphPtr);

        auto send = ZoltanPart{}; // Objects to export/send to "other" processes/parts
        auto recv = ZoltanPart{}; // Objects to import/receive from "other" processes/parts

        auto partitionChanged = 0;
        auto numElmPerGid = 1;
        auto numElmPerLid = 1;
        auto numRecv = 0;
        auto numSend = 0;
        const auto rc = Zoltan_LB_Partition
            (this->zz_,
             &partitionChanged,
             &numElmPerGid, &numElmPerLid,
             &numRecv, &recv.globalId, &recv.localId, &recv.process, &recv.block,
             &numSend, &send.globalId, &send.localId, &send.process, &send.block);

        if (rc != ZOLTAN_OK) {
            OPM_THROW(std::runtime_error, "Failed to Partition Domain Graph");
        }

        auto blocks = std::vector<int>(numCells, 0);
        for (auto e = 0*numSend; e < numSend; ++e) {
            blocks[send.localId[e]] = send.block[e];
        }

        return blocks;
    }

    void forceSameDomain(const std::vector<int>& cells,
                         std::vector<int>&       parts)
    {
        if (cells.empty()) { return; }

        const auto first = parts[cells.front()];
        for (auto& cell : cells) {
            parts[cell] = first;
        }
    }

    std::vector<int> compressPartitionIDs(std::vector<int>&& parts0)
    {
        auto partition = std::move(parts0);

        if (! partition.empty()) {
            auto mmPos = std::minmax_element(partition.begin(), partition.end());
            auto seen = std::vector<bool>(*mmPos.second - *mmPos.first + 1, false);
            for (const auto& domain : partition) {
                seen[domain - *mmPos.first] = domain >= 0;
            }

            auto num_domains = 0;
            auto compressed = std::vector<int>(*mmPos.second - *mmPos.first + 1, -1);
            for (auto i = 0*compressed.size(); i < compressed.size(); ++i) {
                if (seen[i]) {
                    compressed[i] = num_domains++;
                }
            }

            for (auto& domain : partition) {
                if (domain >= 0) {
                    domain = compressed[domain - *mmPos.first];
                }
            }
        }

        return partition;
    }

} // Anonymous namespace

std::vector<int>
Opm::ParallelNLDDPartitioningZoltan::partitionElements(const ZoltanParamMap& params) const
{
    const auto vertexId = EnumerateSeenVertices { this->numElements_, this->conns_ };

    auto graph = VertexGraph {
        this->comm_.rank(), this->numElements_,
        this->conns_, vertexId, this->globalCell_
    };

    const auto partsForReachableCells = Partitioner {
        this->comm_, params
    }(static_cast<void*>(&graph), graph.numVertices());

    // Map reachable cells back to full cell numbering.
    auto parts = std::vector<int>(this->numElements_, -1);
    for (auto elm = 0*this->numElements_; elm < this->numElements_; ++elm) {
        if (const auto reachableElmIx = vertexId[elm]; reachableElmIx >= 0) {
            parts[elm] = partsForReachableCells[reachableElmIx];
        }
    }

    for (const auto& cells : this->sameDomain_) {
        ::forceSameDomain(cells, parts);
    }

    return compressPartitionIDs(std::move(parts));
}
