// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2023 Inria, Bretagne–Atlantique Research Center
  
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::DamarisWriter
 */
#ifndef EWOMS_DAMARIS_WRITER_HH
#define EWOMS_DAMARIS_WRITER_HH

#include <dune/grid/common/partitionset.hh>

#include <ebos/collecttoiorank.hh>
#include <ebos/eclbasevanguard.hh>
#include <ebos/eclgenericwriter.hh>
#include <ebos/ecloutputblackoilmodule.hh>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <limits>
#include <stdexcept>
#include <string>

#include <fmt/format.h>

#include <Damaris.h>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutput {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutputCollective {
    using type = UndefinedProperty;
};
} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass them to Damaris server processes.
 *
 * Currently only passing through PRESSURE, GLOBAL_CELL_INDEX and MPI_RANK information.
 * This clss will be enhanced to pass through the 3D mesh information to Damaris to enable
 * in situ visualization via Paraview or Ascent. And developed so that variables specified 
 * through the Eclipse input deck will be available to Damaris.
 */
 
 
template <class TypeTag>
class DamarisWriter : public EclGenericWriter<GetPropType<TypeTag, Properties::Grid>,
                                          GetPropType<TypeTag, Properties::EquilGrid>,
                                          GetPropType<TypeTag, Properties::GridView>,
                                          GetPropType<TypeTag, Properties::ElementMapper>,
                                          GetPropType<TypeTag, Properties::Scalar>>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    
    using BaseType = EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>;
    
public:
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableDamarisOutputCollective,
                             "Write output via Damaris using parallel HDF5 to get single file per timestep instead of one per Damaris core.");
    }

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    DamarisWriter(Simulator& simulator)
        : BaseType(simulator.vanguard().schedule(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().summaryConfig(),
                   simulator.vanguard().grid(),
                   simulator.vanguard().grid().comm().rank() == 0 ? &simulator.vanguard().equilGrid() : nullptr,
                   simulator.vanguard().gridView(),
                   simulator.vanguard().cartesianIndexMapper(),
                   simulator.vanguard().grid().comm().rank() == 0 ? &simulator.vanguard().equilCartesianIndexMapper() : nullptr,
                   false, false)
        , simulator_(simulator)
    {
        this->damarisUpdate_ = true ;

        rank_ = simulator_.vanguard().grid().comm().rank() ;
        nranks_ = simulator_.vanguard().grid().comm().size();
        
        const auto& gridView = simulator_.gridView();
        const auto& interior_elements = elements(gridView, Dune::Partitions::interior);
        // Get the size of the unique vector elements (excludes the shared 'ghost' elements)
        numElements_ = std::distance(interior_elements.begin(), interior_elements.end());

        this->elements_rank_offsets_.resize(nranks_) ;
        this->damarisOutputModule_ = std::make_unique<EclOutputBlackOilModule<TypeTag>>(simulator, this->collectToIORank_);
    }

    ~DamarisWriter()
    { }

    /*!
     * \brief Writes localCellData through to Damaris servers. Sets up the unstructured mesh which is passed to Damaris.
     */
    void writeOutput(data::Solution& localCellData , bool isSubStep)
    {
        OPM_TIMEBLOCK(writeOutput);
        const int reportStepNum = simulator_.episodeIndex() + 1;

        // added this as localCellData was not being written
        if (!isSubStep)
            this->damarisOutputModule_->invalidateLocalData() ;  
        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->damarisOutputModule_->outputErrorLog(simulator_.gridView().comm());

        // The damarisWriter is not outputing well or aquifer data (yet)
        auto localWellData = simulator_.problem().wellModel().wellData(); // data::Well

        if (! isSubStep) 
        {
            if (localCellData.size() == 0) {
                this->damarisOutputModule_->assignToSolution(localCellData);
            }

            // add cell data to perforations for Rft output
            this->damarisOutputModule_->addRftDataToWells(localWellData, reportStepNum);
            
            // On first call and if the mesh and variable size change then set damarisUpdate_ to true
            if (damarisUpdate_ == true) {
                // Sets the damaris parameter values "n_elements_local" and "n_elements_total" 
                // which define sizes of the Damaris variables, per-rank and globally (over all ranks).
                // Also sets the offsets to where a ranks array data sits within the global array. 
                // This is usefull for HDF5 output and for defining distributed arrays in Dask.
                this->setupDamarisWritingPars(simulator_.vanguard().grid().comm(), numElements_, elements_rank_offsets_);
                
                // sets data for non-time-varying variables MPI_RANK and GLOBAL_CELL_INDEX
                this->setGlobalIndexForDamaris() ; 
                
                // Currently by default we assume static grid (unchanging through the simulation)
                // Set damarisUpdate_ to true if we want to update the geometry to sent to Damaris 
                this->damarisUpdate_ = false; 
            }

            if (this->damarisOutputModule_->getPRESSURE_ptr() != nullptr) 
            {
                int64_t temp_int64_t[1];
                temp_int64_t[0] = static_cast<int64_t>(this->elements_rank_offsets_[rank_]);
                dam_err_ = damaris_set_position("PRESSURE", temp_int64_t);
                if (dam_err_ != DAMARIS_OK && rank_ == 0) {
                    OpmLog::error(fmt::format("ERORR: damariswriter::writeOutput()       : ( rank:{}) damaris_set_position(PRESSURE, ...), Damaris Error: {}  ",  rank_, damaris_error_string(dam_err_) ));
                }

                dam_err_ = damaris_write("PRESSURE", (void*)this->damarisOutputModule_->getPRESSURE_ptr());
                if (dam_err_ != DAMARIS_OK) {
                   OpmLog::error(fmt::format("ERORR: damariswriter::writeOutput()       : ( rank:{}) damaris_write(PRESSURE, ...), Damaris Error: {}  ",  rank_, damaris_error_string(dam_err_) ));
                }

                dam_err_ =  damaris_end_iteration();
                if (dam_err_ != DAMARIS_OK) {
                    OpmLog::error(fmt::format("ERORR: damariswriter::writeOutput()       : ( rank:{}) damaris_end_iteration(), Damaris Error: {}  ",  rank_, damaris_error_string(dam_err_) ));
                }
            }
         } // end of ! isSubstep
    }

private:
    int dam_err_ ;
    int rank_  ;       
    int nranks_ ;
    int numElements_ ;  ///<  size of the unique vector elements
    
    Simulator& simulator_;
    std::unique_ptr<EclOutputBlackOilModule<TypeTag>> damarisOutputModule_;
    std::vector<unsigned long long> elements_rank_offsets_ ;
    bool damarisUpdate_ = false;  ///< Whenever this is true writeOutput() will set up Damaris mesh information and offsets of model fields

    static bool enableDamarisOutput_()
    { 
        return EWOMS_GET_PARAM(TypeTag, bool, EnableDamarisOutput); 
    }

    void setGlobalIndexForDamaris () 
    {
        // GLOBAL_CELL_INDEX is used to reorder variable data when writing to disk 
        // This is enabled using select-file="GLOBAL_CELL_INDEX" in the <variable> XML tag
        if ( this->collectToIORank_.isParallel() ){
            const std::vector<int>& local_to_global =  this->collectToIORank_.localIdxToGlobalIdxMapping(); 
            dam_err_ = damaris_write("GLOBAL_CELL_INDEX", local_to_global.data());
        } else {
            std::vector<int> local_to_global_filled ;
            local_to_global_filled.resize(this->numElements_) ;
            for (int i = 0 ; i < this->numElements_ ; i++)
            {
                local_to_global_filled[i] = i ;
            }
            dam_err_ = damaris_write("GLOBAL_CELL_INDEX", local_to_global_filled.data());
        }

        if (dam_err_ != DAMARIS_OK) {
            OpmLog::error(fmt::format("ERORR: damariswriter::writeOutput()       :"
                                     "( rank:{}) damaris_write(GLOBAL_CELL_INDEX, ...), Damaris Error: {}  ",  
                                     rank_, damaris_error_string(dam_err_) ));
        }

        std::vector<int> mpiRank(this->numElements_, rank_ ) ;
        dam_err_ = damaris_write("MPI_RANK", mpiRank.data() ) ;
        if (dam_err_ != DAMARIS_OK) {
           OpmLog::error(fmt::format("ERORR: damariswriter::writeOutput()       :"
                                     " ( rank:{}) damaris_write(MPI_RANK, ...), Damaris Error: {}  ",
                                     rank_, damaris_error_string(dam_err_) ));
        }
    }

    void setupDamarisWritingPars(Parallel::Communication comm, const int n_elements_local_grid, std::vector<unsigned long long>& elements_rank_offsets)
    {
        // one for each rank -- to be gathered from each client rank
        std::vector<unsigned long long> elements_rank_sizes(nranks_); 
        // n_elements_local_grid should be the full model size
        const unsigned long long n_elements_local = n_elements_local_grid;

        // This gets the n_elements_local from all ranks and copies them to a std::vector of all the values on all ranks
        // (elements_rank_sizes[]).
        comm.allgather(&n_elements_local, 1, elements_rank_sizes.data());
        elements_rank_offsets[0] = 0ULL;
        // This scan makes the offsets to the start of each ranks grid section if each local grid data was concatenated (in
        // rank order)
        for (int t1 = 1; t1 < nranks_; t1++) {
            elements_rank_offsets[t1] = elements_rank_offsets[t1 - 1] + elements_rank_sizes[t1 - 1];
        }

        // find the global/total size
        unsigned long long n_elements_global_max = elements_rank_offsets[nranks_ - 1];
        n_elements_global_max += elements_rank_sizes[nranks_ - 1]; // add the last ranks size to the already accumulated offset values

        if (rank_ == 0) {
            OpmLog::debug(fmt::format("In setupDamarisWritingPars(): n_elements_global_max = {}", n_elements_global_max));
        }

        // Set the paramater so that the Damaris servers can allocate the correct amount of memory for the variabe
        // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
        // ToDo: Do we need to check that local ranks are 0 based ?
        int temp_int = static_cast<int>(elements_rank_sizes[rank_]);
        dam_err_ = damaris_parameter_set("n_elements_local", &temp_int, sizeof(int));
        if (dam_err_ != DAMARIS_OK && rank_ == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_parameter_set(\"n_elements_local\", &temp_int, sizeof(int));");
        }
        // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
        // ToDo: Do we need to check that n_elements_global_max will fit in a C int type (INT_MAX)
        temp_int = static_cast<int>(n_elements_global_max);
        dam_err_ = damaris_parameter_set("n_elements_total", &temp_int, sizeof(int));
        if (dam_err_ != DAMARIS_OK && rank_ == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_parameter_set(\"n_elements_total\", &temp_int, sizeof(int));");
        }

        // Use damaris_set_position to set the offset in the global size of the array.
        // This is used so that output functionality (e.g. HDF5Store) knows global offsets of the data of the ranks
        int64_t temp_int64_t[1];
        temp_int64_t[0] = static_cast<int64_t>(elements_rank_offsets[rank_]);
        dam_err_ = damaris_set_position("PRESSURE", temp_int64_t);
        if (dam_err_ != DAMARIS_OK && rank_ == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_set_position(\"PRESSURE\", temp_int64_t);");
        }
        dam_err_ = damaris_set_position("GLOBAL_CELL_INDEX", temp_int64_t);
        if (dam_err_ != DAMARIS_OK && rank_ == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_set_position(\"GLOBAL_CELL_INDEX\", temp_int64_t);");
        }
        dam_err_ = damaris_set_position("MPI_RANK", temp_int64_t);
        if (dam_err_ != DAMARIS_OK && rank_ == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_set_position(\"MPI_RANK\", temp_int64_t);");
        }
    }

    void prepareLocalCellData(const bool isSubStep,
                              const int  reportStepNum)
    {
        OPM_TIMEBLOCK(prepareLocalCellData);
        if (damarisOutputModule_->localDataValid()) {
            return;
        }

        const auto& gridView = simulator_.vanguard().gridView();
        const int numElements = gridView.size(/*codim=*/0);
        const bool log = this->collectToIORank_.isIORank();

        damarisOutputModule_->allocBuffers(numElements, reportStepNum,
                                      isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
        OPM_TIMEBLOCK(prepareCellBasedData);
        for (const auto& elem : elements(gridView)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            damarisOutputModule_->processElement(elemCtx);
        }
        }
        if(!simulator_.model().linearizer().getFlowsInfo().empty()){
            OPM_TIMEBLOCK(prepareFlowsData);
            for (const auto& elem : elements(gridView)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                damarisOutputModule_->processElementFlows(elemCtx);
            }
        }
        {
        OPM_TIMEBLOCK(prepareBlockData);
        for (const auto& elem : elements(gridView)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            damarisOutputModule_->processElementBlockData(elemCtx);
        }
        }
        {
        OPM_TIMEBLOCK(prepareFluidInPlace);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int dofIdx=0; dofIdx < numElements; ++dofIdx){
                const auto& intQuants = *(simulator_.model().cachedIntensiveQuantities(dofIdx, /*timeIdx=*/0));
                const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);
                damarisOutputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
        }
        }
        damarisOutputModule_->validateLocalData();
        OPM_END_PARALLEL_TRY_CATCH("DamarisWriter::prepareLocalCellData() failed: ", simulator_.vanguard().grid().comm());
    }
};

} // namespace Opm

#endif
