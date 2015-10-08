#include <opm/core/utility/Units.hpp>
#include <opm/core/grid/GridHelpers.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/compressedToCartesian.hpp>

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

double computeWellIndex(const double radius,
                        const std::array<double, 3>& cubical,
                        const double* cell_permeability,
                        const double skin_factor,
                        const Opm::WellCompletion::DirectionEnum direction,
                        const double ntg);

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
void WellsManager::createWellsFromSpecs(std::vector<WellConstPtr>& wells, size_t timeStep,
                                        const C2F& c2f,
                                        const int* cart_dims,
                                        FC begin_face_centroids,
                                        int dimensions,
                                        std::vector<double>& dz,
                                        std::vector<std::string>& well_names,
                                        std::vector<WellData>& well_data,
                                        std::map<std::string, int>& well_names_to_index,
                                        const PhaseUsage& phaseUsage,
                                        const std::map<int,int>& cartesian_to_compressed,
                                        const double* permeability,
                                        const NTG& ntg,
                                        std::vector<int>& wells_on_proc)
{
    if (dimensions != 3) {
        OPM_THROW(std::domain_error,
                  "WellsManager::createWellsFromSpecs() only "
                  "supported in three space dimensions");
    }

    std::vector<std::vector<PerfData> > wellperf_data;
    wellperf_data.resize(wells.size());
    wells_on_proc.resize(wells.size(), 1);

    int well_index = 0;
    for (auto wellIter= wells.begin(); wellIter != wells.end(); ++wellIter) {
        WellConstPtr well = (*wellIter);

        if (well->getStatus(timeStep) == WellCommon::SHUT) {
            continue;
        }

        {   // COMPDAT handling
            CompletionSetConstPtr completionSet = well->getCompletions(timeStep);
            // shut completions and open ones stored in this process will have 1 others 0.
            std::vector<std::size_t> completion_on_proc(completionSet->size(), 1);
            std::size_t shut_completions_number = 0;
            for (size_t c=0; c<completionSet->size(); c++) {
                CompletionConstPtr completion = completionSet->get(c);
                if (completion->getState() == WellCompletion::OPEN) {
                    int i = completion->getI();
                    int j = completion->getJ();
                    int k = completion->getK();

                    const int* cpgdim = cart_dims;
                    int cart_grid_indx = i + cpgdim[0]*(j + cpgdim[1]*k);
                    std::map<int, int>::const_iterator cgit = cartesian_to_compressed.find(cart_grid_indx);
                    if (cgit == cartesian_to_compressed.end()) {
                        if ( is_parallel_run_ )
                        {
                            completion_on_proc[c]=0;
                            continue;
                        }
                        else
                        {
                            OPM_MESSAGE("****Warning: Cell with i,j,k indices " << i << ' ' << j << ' '
                                      << k << " not found in grid. The completion will be igored (well = "
                                      << well->name() << ')');
                        }
                    }
                    else
                    {
                        int cell = cgit->second;
                        PerfData pd;
                        pd.cell = cell;
                        {
                            const Value<double>& transmissibilityFactor = completion->getConnectionTransmissibilityFactorAsValueObject();
                            const double wellPi = completion ->getWellPi();
                            if (transmissibilityFactor.hasValue()) {
                                pd.well_index = transmissibilityFactor.getValue();
                            } else {
                                double radius = 0.5*completion->getDiameter();
                                if (radius <= 0.0) {
                                    radius = 0.5*unit::feet;
                                    OPM_MESSAGE("**** Warning: Well bore internal radius set to " << radius);
                                }

                                std::array<double, 3> cubical =
                                    WellsManagerDetail::getCubeDim<3>(c2f, begin_face_centroids, cell);

                                // overwrite dz values calculated in getCubeDim.
                                if (dz.size() > 0) {
                                    cubical[2] = dz[cell];
                                }

                                const double* cell_perm = &permeability[dimensions*dimensions*cell];
                                pd.well_index =
                                    WellsManagerDetail::computeWellIndex(radius, cubical, cell_perm,
                                                                         completion->getSkinFactor(),
                                                                         completion->getDirection(),
                                                                         ntg[cell]);
                            }
                            pd.well_index *= wellPi;
                        }
                        wellperf_data[well_index].push_back(pd);
                    }
                } else {
                    ++shut_completions_number;
                    if (completion->getState() != WellCompletion::SHUT) {
                        OPM_THROW(std::runtime_error, "Completion state: " << WellCompletion::StateEnum2String( completion->getState() ) << " not handled");
                    }
                }
            }
            if ( is_parallel_run_ )
            {
                // sum_completions_on_proc includes completions
                // that are shut
                std::size_t sum_completions_on_proc = std::accumulate(completion_on_proc.begin(),
                                                                      completion_on_proc.end(),0);
                // Set wells that are not on this processor to SHUT.
                // A well is not here if only shut completions are found.
                if ( sum_completions_on_proc == shut_completions_number )
                {
                    // Mark well as not existent on this process
                    wells_on_proc[wellIter-wells.begin()] = 0;
                    continue;
                }
                else
                {
                    // Check that the complete well is on this process
                    if ( sum_completions_on_proc < completionSet->size() )
                    {
                        std::cout<< "Well "<< well->name() << " semms not be in "
                                 << "completely in the disjoint partition of "
                                 << "process deactivating here." << std::endl;
                        // Mark well as not existent on this process
                        wells_on_proc[wellIter-wells.begin()] = 0;
                        continue;
                    }
                }
            }
        }
        {   // WELSPECS handling
            well_names_to_index[well->name()] = well_index;
            well_names.push_back(well->name());
            {
                WellData wd;
                wd.reference_bhp_depth = well->getRefDepth();
                wd.welspecsline = -1;
                if (well->isInjector( timeStep ))
                    wd.type = INJECTOR;
                else
                    wd.type = PRODUCER;
                well_data.push_back(wd);
            }
        }

        well_index++;
    }
    // Set up reference depths that were defaulted. Count perfs.

    const int num_wells = well_data.size();

    int num_perfs = 0;
    assert (dimensions == 3);
    for (int w = 0; w < num_wells; ++w) {
        num_perfs += wellperf_data[w].size();
    }

    // Create the well data structures.
    w_ = create_wells(phaseUsage.num_phases, num_wells, num_perfs);
    if (!w_) {
        OPM_THROW(std::runtime_error, "Failed creating Wells struct.");
    }


    // Add wells.
    for (int w = 0; w < num_wells; ++w) {
        const int           w_num_perf = wellperf_data[w].size();
        std::vector<int>    perf_cells  (w_num_perf);
        std::vector<double> perf_prodind(w_num_perf);

        for (int perf = 0; perf < w_num_perf; ++perf) {
            perf_cells  [perf] = wellperf_data[w][perf].cell;
            perf_prodind[perf] = wellperf_data[w][perf].well_index;
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
                     well_names[w].c_str(),
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
WellsManager(const Opm::EclipseStateConstPtr eclipseState,
             const size_t                    timeStep,
             int                             number_of_cells,
             const int*                      global_cell,
             const int*                      cart_dims,
             int                             dimensions,
             const C2F&                      cell_to_faces,
             FC                              begin_face_centroids,
             const double*                   permeability,
             bool                            is_parallel_run)
    : w_(0), is_parallel_run_(is_parallel_run)
{
    init(eclipseState, timeStep, number_of_cells, global_cell,
         cart_dims, dimensions,
         cell_to_faces, begin_face_centroids, permeability);
}

/// Construct wells from deck.
template <class C2F, class FC>
void
WellsManager::init(const Opm::EclipseStateConstPtr eclipseState,
                   const size_t                    timeStep,
                   int                             number_of_cells,
                   const int*                      global_cell,
                   const int*                      cart_dims,
                   int                             dimensions,
                   const C2F&                      cell_to_faces,
                   FC                              begin_face_centroids,
                   const double*                   permeability)
{
    if (dimensions != 3) {
        OPM_THROW(std::runtime_error,
                  "We cannot initialize wells from a deck unless "
                  "the corresponding grid is 3-dimensional.");
    }

    if (eclipseState->getSchedule()->numWells() == 0) {
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

    ScheduleConstPtr          schedule = eclipseState->getSchedule();
    std::vector<WellConstPtr> wells    = schedule->getWells(timeStep);
    std::vector<int>          wells_on_proc;

    well_names.reserve(wells.size());
    well_data.reserve(wells.size());

    typedef GridPropertyAccess::ArrayPolicy::ExtractFromDeck<double> DoubleArray;
    typedef GridPropertyAccess::Compressed<DoubleArray, GridPropertyAccess::Tag::NTG> NTGArray;

    DoubleArray ntg_glob(eclipseState, "NTG", 1.0);
    NTGArray    ntg(ntg_glob, global_cell);

    EclipseGridConstPtr eclGrid = eclipseState->getEclipseGrid();

    // use cell thickness (dz) from eclGrid
    // dz overwrites values calculated by WellDetails::getCubeDim
    std::vector<double> dz(number_of_cells);
    {
        std::vector<int> gc = compressedToCartesian(number_of_cells, global_cell);
        for (int cell = 0; cell < number_of_cells; ++cell) {
            dz[cell] = eclGrid->getCellThicknes(gc[cell]);
        }
    }

    createWellsFromSpecs(wells, timeStep, cell_to_faces,
                         cart_dims,
                         begin_face_centroids,
                         dimensions,
                         dz,
                         well_names, well_data, well_names_to_index,
                         pu, cartesian_to_compressed, permeability, ntg,
                         wells_on_proc);

    setupWellControls(wells, timeStep, well_names, pu, wells_on_proc);

    {
        GroupTreeNodeConstPtr fieldNode =
            schedule->getGroupTree(timeStep)->getNode("FIELD");

        GroupConstPtr fieldGroup =
            schedule->getGroup(fieldNode->name());

        well_collection_.addField(fieldGroup, timeStep, pu);
        addChildGroups(fieldNode, schedule, timeStep, pu);
    }

    for (auto w = wells.begin(), e = wells.end(); w != e; ++w) {
        well_collection_.addWell(*w, timeStep, pu);
    }

    well_collection_.setWellsPointer(w_);
    well_collection_.applyGroupControls();

    setupGuideRates(wells, timeStep, well_data, well_names_to_index);

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
