#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/grid/GridHelpers.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/grid/utility/compressedToCartesian.hpp>
#include <opm/core/props/rock/RockFromDeck.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellConnections.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <iterator>
#include <numeric>

namespace WellsManagerDetail
{


    namespace ProductionControl
    {
        enum Mode { ORAT, WRAT, GRAT,
                    LRAT, CRAT, RESV,
                    BHP , THP , GRUP };
    /*
        namespace Details {
            std::map<std::string, Mode>
            init_mode_map();
        } // namespace Details
    */
        Mode mode(const std::string& control);


        Mode mode(Opm::WellProducer::ControlModeEnum controlMode);
    } // namespace ProductionControl


    namespace InjectionControl
    {
        enum Mode { RATE, RESV, BHP,
                    THP, GRUP };
    /*
        namespace Details {
            std::map<std::string, Mode>
            init_mode_map();
        } // namespace Details
    */
        Mode mode(const std::string& control);

        Mode mode(Opm::WellInjector::ControlModeEnum controlMode);

    } // namespace InjectionControl


template <int dim, class C2F, class FC>
std::array<double, dim>
getCubeDim(const C2F& c2f,
           FC         begin_face_centroids,
           int        cell)
{
    std::array< std::vector<double>, dim > X;
    {
        const std::vector<double>::size_type
            nf = std::distance(c2f[cell].begin(),
                               c2f[cell].end  ());

        for (int d = 0; d < dim; ++d) {
            X[d].reserve(nf);
        }
    }

    typedef typename C2F::row_type::const_iterator FI;

    for (FI f = c2f[cell].begin(), e = c2f[cell].end(); f != e; ++f) {
        using Opm::UgGridHelpers::increment;
        using Opm::UgGridHelpers::getCoordinate;

        const FC& fc = increment(begin_face_centroids, *f, dim);

        for (int d = 0; d < dim; ++d) {
            X[d].push_back(getCoordinate(fc, d));
        }
    }

    std::array<double, dim> cube;
    for (int d = 0; d < dim; ++d) {
        typedef std::vector<double>::iterator VI;
        typedef std::pair<VI,VI>              PVI;

        const PVI m = std::minmax_element(X[d].begin(), X[d].end());

        cube[d] = *m.second - *m.first;
    }

    return cube;
}
} // end namespace WellsManagerDetail

namespace Opm
{
template<class C2F, class FC, class NTG>
void WellsManager::createWellsFromSpecs(const std::vector<Well2>& wells, size_t timeStep,
                                        const C2F& /* c2f */,
                                        const int* cart_dims,
                                        FC /* begin_face_centroids */,
                                        int dimensions,
                                        std::vector<double>& /* dz */,
                                        std::vector<std::string>& well_names,
                                        std::vector<WellData>& well_data,
                                        std::map<std::string, int>& well_names_to_index,
                                        const PhaseUsage& phaseUsage,
                                        const std::map<int,int>& cartesian_to_compressed,
                                        const double* /* permeability */,
                                        const NTG& /* ntg */,
                                        std::vector<int>& wells_on_proc,
                                        const std::unordered_set<std::string>& ignored_wells)
{
    if (dimensions != 3) {
        OPM_THROW(std::domain_error,
                  "WellsManager::createWellsFromSpecs() only "
                  "supported in three space dimensions");
    }

    std::vector<std::vector<PerfData> > wellperf_data;
    wellperf_data.resize(wells.size());
    wells_on_proc.resize(wells.size(), 1);

    // The well index on the current process.
    // Note that some wells are deactivated as they live on the interior
    // domain of another proccess. Therefore this might different from
    // the index of the well according to the eclipse state
    int active_well_index = 0;
    for (auto wellIter= wells.begin(); wellIter != wells.end(); ++wellIter) {
        const auto& well = (*wellIter);

        if (well.getStatus() == WellCommon::SHUT) {
            continue;
        }

        if ( ignored_wells.find(well.name()) != ignored_wells.end() ) {
            wells_on_proc[ wellIter - wells.begin() ] = 0;
            continue;
        }

        {   // COMPDAT handling
            // shut completions and open ones stored in this process will have 1 others 0.

            for(const auto& completion : well.getConnections()) {
                if (completion.state() == WellCompletion::OPEN) {
                    const int i = completion.getI();
                    const int j = completion.getJ();
                    const int k = completion.getK();

                    const int* cpgdim = cart_dims;
                    const int cart_grid_indx = i + cpgdim[0]*(j + cpgdim[1]*k);
                    const std::map<int, int>::const_iterator cgit = cartesian_to_compressed.find(cart_grid_indx);
                    if (cgit == cartesian_to_compressed.end()) {
                        const std::string msg = ("Cell with i,j,k indices " + std::to_string(i) + " " + std::to_string(j)
                                    + " " + std::to_string(k)  +  " not found in grid (well = " +  well.name() + ").");
                        OPM_THROW(std::runtime_error, msg);
                    }
                    else
                    {
                        PerfData pd;
                        pd.cell = cgit->second;
                        pd.well_index = completion.CF() * completion.wellPi();
                        pd.satnumid = completion.satTableId();

                        wellperf_data[active_well_index].push_back(pd);
                    }
                } else {
                    if (completion.state() != WellCompletion::SHUT) {
                        OPM_THROW(std::runtime_error, "Completion state: " << WellCompletion::StateEnum2String( completion.state() ) << " not handled");
                    }
                }
            }
        }

        if (wellperf_data[active_well_index].empty()) {
            const std::string msg = " there is no perforations associated with the well "
                                  + well.name() + ", the well is ignored for the report step "
                                  + std::to_string(timeStep);
            OpmLog::warning(msg);
            wells_on_proc[wellIter - wells.begin()] = 0;
            continue;
        }
        {   // WELSPECS handling
            well_names_to_index[well.name()] = active_well_index;
            well_names.push_back(well.name());
            {
                WellData wd;
                wd.reference_bhp_depth = well.getRefDepth( );
                wd.welspecsline = -1;
                if (well.isInjector( ))
                    wd.type = INJECTOR;
                else
                    wd.type = PRODUCER;

                wd.allowCrossFlow = well.getAllowCrossFlow();
                well_data.push_back(wd);
            }
        }

        active_well_index++;
    }
    // Set up reference depths that were defaulted. Count perfs.

    const int num_wells = well_data.size();

    int num_perfs = 0;
    assert (dimensions == 3);
    for (int w = 0; w < num_wells; ++w) {
        num_perfs += wellperf_data[w].size();
    }
    // Create the well data structures.
    struct Wells* w = create_wells(phaseUsage.num_phases, num_wells, num_perfs);

    if (!w) {
        OPM_THROW(std::runtime_error, "Failed creating Wells struct.");
    }

    std::swap( w, w_ );
    destroy_wells( w );

    // Add wells.
    for (int w = 0; w < num_wells; ++w) {
        const int           w_num_perf = wellperf_data[w].size();
        std::vector<int>    perf_cells  (w_num_perf);
        std::vector<double> perf_prodind(w_num_perf);
        std::vector<int> perf_satnumid(w_num_perf);

        for (int perf = 0; perf < w_num_perf; ++perf) {
            perf_cells  [perf] = wellperf_data[w][perf].cell;
            perf_prodind[perf] = wellperf_data[w][perf].well_index;
            perf_satnumid[perf] = wellperf_data[w][perf].satnumid;
        }

        const double* comp_frac = NULL;

        // We initialize all wells with a null component fraction,
        // and must (for injection wells) overwrite it later.
        const int ok =
            add_well(well_data[w].type,
                     well_data[w].reference_bhp_depth,
                     w_num_perf,
                     comp_frac,
                     perf_cells.data(),
                     perf_prodind.data(),
                     perf_satnumid.data(),
                     well_names[w].c_str(),
                     well_data[w].allowCrossFlow,
                     w_);

        if (!ok) {
            OPM_THROW(std::runtime_error,
                      "Failed adding well "
                      << well_names[w]
                      << " to Wells data structure.");
        }
    }
}

template <class C2F, class FC>
WellsManager::
WellsManager(const Opm::EclipseState& eclipseState,
             const Opm::Schedule& schedule,
             const Opm::SummaryState& summaryState,
             const size_t                    timeStep,
             int                             number_of_cells,
             const int*                      global_cell,
             const int*                      cart_dims,
             int                             dimensions,
             const C2F&                      cell_to_faces,
             FC                              begin_face_centroids,
             bool                            is_parallel_run,
             const std::unordered_set<std::string>&    deactivated_wells)
    : w_(create_wells(0,0,0)), is_parallel_run_(is_parallel_run)
{
  init(eclipseState, schedule, timeStep, number_of_cells, global_cell,
         cart_dims, dimensions,
         cell_to_faces, begin_face_centroids, deactivated_wells);
}

/// Construct wells from deck.
template <class C2F, class FC>
void
WellsManager::init(const Opm::EclipseState& eclipseState,
                   const Opm::Schedule& schedule,
                   const size_t                    timeStep,
                   int                             number_of_cells,
                   const int*                      global_cell,
                   const int*                      cart_dims,
                   int                             dimensions,
                   const C2F&                      cell_to_faces,
                   FC                              begin_face_centroids,
                   const std::unordered_set<std::string>&    deactivated_wells)
{
    if (dimensions != 3) {
        OPM_THROW(std::runtime_error,
                  "We cannot initialize wells from a deck unless "
                  "the corresponding grid is 3-dimensional.");
    }

    if (schedule.numWells() == 0) {
        OPM_MESSAGE("No wells specified in Schedule section, "
                    "initializing no wells");
        return;
    }

    std::map<int,int> cartesian_to_compressed;
    setupCompressedToCartesian(global_cell, number_of_cells,
                               cartesian_to_compressed);

    // Obtain phase usage data.
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);

    // These data structures will be filled in this constructor,
    // then used to initialize the Wells struct.
    std::vector<std::string> well_names;
    std::vector<WellData> well_data;


    // For easy lookup:
    std::map<std::string, int> well_names_to_index;

    const auto wells = schedule.getWells2(timeStep);
    std::vector<int> wells_on_proc;

    well_names.reserve(wells.size());
    well_data.reserve(wells.size());

    typedef GridPropertyAccess::ArrayPolicy::ExtractFromDeck<double> DoubleArray;
    typedef GridPropertyAccess::Compressed<DoubleArray, GridPropertyAccess::Tag::NTG> NTGArray;

    DoubleArray ntg_glob(eclipseState, "NTG", 1.0);
    NTGArray    ntg(ntg_glob, global_cell);

    const auto& eclGrid = eclipseState.getInputGrid();

    // use cell thickness (dz) from eclGrid
    // dz overwrites values calculated by WellDetails::getCubeDim
    std::vector<double> dz(number_of_cells);
    {
        std::vector<int> gc = compressedToCartesian(number_of_cells, global_cell);
        for (int cell = 0; cell < number_of_cells; ++cell) {
            dz[cell] = eclGrid.getCellThicknes(gc[cell]);
        }
    }

    std::vector<double> interleavedPerm;
    RockFromDeck::extractInterleavedPermeability(eclipseState,
                                                 number_of_cells,
                                                 global_cell,
                                                 cart_dims,
                                                 0.0,
                                                 interleavedPerm);

    createWellsFromSpecs(wells, timeStep, cell_to_faces,
                         cart_dims,
                         begin_face_centroids,
                         dimensions,
                         dz,
                         well_names, well_data, well_names_to_index,
                         pu, cartesian_to_compressed, interleavedPerm.data(), ntg,
                         wells_on_proc, deactivated_wells);

    setupWellControls(wells, timeStep, well_names, pu, wells_on_proc);

    {
        const auto& fieldGroup = schedule.getGroup( "FIELD" );
        well_collection_.addField(fieldGroup, timeStep, pu);

        const auto& grouptree = schedule.getGroupTree( timeStep );
        std::vector< std::string > group_stack = { "FIELD" };

        do {
            auto parent = group_stack.back();
            group_stack.pop_back();
            const auto& children = grouptree.children( parent );
            group_stack.insert( group_stack.end(), children.begin(), children.end() );

            for( const auto& child : children ) {
                well_collection_.addGroup( schedule.getGroup( child ), parent, timeStep, pu );
            }

        } while( !group_stack.empty() );
    }

    for (size_t i = 0; i < wells_on_proc.size(); ++i) {
        // wells_on_proc is a vector of flag to indicate whether a well is on the process
        if (wells_on_proc[i]) {
            well_collection_.addWell(wells[i], timeStep, pu);
        }
    }

    well_collection_.setWellsPointer(w_);

    if (well_collection_.groupControlActive()) {
        // here does not consider the well potentials related guide rate setting
        setupGuideRates(wells, timeStep, well_data, well_names_to_index);
    }

    // Debug output.
#define EXTRA_OUTPUT
#ifdef EXTRA_OUTPUT
    /*
      std::cout << "\t WELL DATA" << std::endl;
      for(int i = 0; i< num_wells; ++i) {
      std::cout << i << ": " << well_data[i].type << "  "
      << well_data[i].control << "  " << well_data[i].target
      << std::endl;
      }

      std::cout << "\n\t PERF DATA" << std::endl;
      for(int i=0; i< int(wellperf_data.size()); ++i) {
      for(int j=0; j< int(wellperf_data[i].size()); ++j) {
      std::cout << i << ": " << wellperf_data[i][j].cell << "  "
      << wellperf_data[i][j].well_index << std::endl;
      }
      }
    */
#endif
}

} // end namespace Opm
