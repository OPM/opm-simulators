/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

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


#ifndef OPM_STANDARDWELLSDENSE_HEADER_INCLUDED
#define OPM_STANDARDWELLSDENSE_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>
#include <tuple>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <opm/core/wells.h>
#include <opm/core/wells/DynamicListEconLimited.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterleaved.cpp>


#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

namespace Opm {


        /// Class for handling the standard well model.
        template<typename FluidSystem, typename BlackoilIndices>
        class StandardWellsDense {
        public:
            struct WellOps {
                WellOps(const Wells* wells)
                  : w2p(),
                    p2w(),
                    well_cells()
                {
                    if( wells )
                    {
                        w2p = Eigen::SparseMatrix<double>(wells->well_connpos[ wells->number_of_wells ], wells->number_of_wells);
                        p2w = Eigen::SparseMatrix<double>(wells->number_of_wells, wells->well_connpos[ wells->number_of_wells ]);

                        const int        nw   = wells->number_of_wells;
                        const int* const wpos = wells->well_connpos;

                        typedef Eigen::Triplet<double> Tri;

                        std::vector<Tri> scatter, gather;
                        scatter.reserve(wpos[nw]);
                        gather .reserve(wpos[nw]);

                        for (int w = 0, i = 0; w < nw; ++w) {
                            for (; i < wpos[ w + 1 ]; ++i) {
                                scatter.push_back(Tri(i, w, 1.0));
                                gather .push_back(Tri(w, i, 1.0));
                            }
                        }

                        w2p.setFromTriplets(scatter.begin(), scatter.end());
                        p2w.setFromTriplets(gather .begin(), gather .end());

                        well_cells.assign(wells->well_cells, wells->well_cells + wells->well_connpos[wells->number_of_wells]);
                    }
                }
                Eigen::SparseMatrix<double> w2p;              // well -> perf (scatter)
                Eigen::SparseMatrix<double> p2w;              // perf -> well (gather)
                std::vector<int> well_cells;                  // the set of perforated cells
            };

            // ---------      Types      ---------
            using ADB = AutoDiffBlock<double>;

            typedef DenseAd::Evaluation<double, /*size=*/6> EvalWell;
            typedef WellStateFullyImplicitBlackoil WellState;
            typedef BlackoilModelParameters ModelParameters;

            //typedef AutoDiffBlock<double> ADB;
            using Vector = ADB::V;
            using V = ADB::V;
            typedef ADB::M M;
            typedef double Scalar;
            typedef Dune::FieldVector<Scalar, 3    >       VectorBlockType;
            typedef Dune::FieldMatrix<Scalar, 3, 3 >      MatrixBlockType;
            typedef Dune::BCRSMatrix <MatrixBlockType>      Mat;
            typedef Dune::BlockVector<VectorBlockType>      BVector;


            // copied from BlackoilModelBase
            // should put to somewhere better
            using DataBlock =  Eigen::Array<double,
                                            Eigen::Dynamic,
                                            Eigen::Dynamic,
                                            Eigen::RowMajor>;
            // ---------  Public methods  ---------
            StandardWellsDense(const Wells* wells_arg,
                               const ModelParameters& param,
                               const bool terminal_output,
                               const std::vector<double>& pv)
                : wells_active_(wells_arg!=nullptr)
                , wells_(wells_arg)
                , wops_(wells_arg)
                , fluid_(nullptr)
                , active_(nullptr)
                , phase_condition_(nullptr)
                , vfp_properties_(nullptr)
                , well_perforation_densities_(Vector())
                , well_perforation_pressure_diffs_(Vector())
                , store_well_perforation_fluxes_(false)
                , wellVariables_(wells_arg->number_of_wells * wells_arg->number_of_phases)
                , F0_(wells_arg->number_of_wells * wells_arg->number_of_phases)
                , param_(param)
                , terminal_output_(terminal_output)
                , pv_(pv)
                //, resWell(wells_->number_of_wells)
                //, duneD( Mat(2, 2, 3, 0.4, Mat::implicit))
                //, duneB( Mat(300, 2, 6, 0.4, Mat::implicit))
                //, duneC( Mat(2, 300, 6, 0.4, Mat::implicit))
              {
              }


            template <typename Simulator>
            IterationReport assemble(const Simulator& ebosSimulator,
                                     const int iterationIdx,
                                     const double dt,
                                     WellState& well_state,
                                     LinearisedBlackoilResidual& residual) {
                resetWellControlFromState(well_state);
                updateWellControls(well_state);
                // Set the primary variables for the wells
                setWellVariables(well_state);

                if (iterationIdx == 0) {
                    computeWellConnectionPressures(ebosSimulator, well_state);
                    computeAccumWells();
                }

                IterationReport iter_report = {false, false, 0, 0};
                if ( ! wellsActive() ) {
                    return iter_report;
                }
                //wellModel().extractWellPerfProperties(state, rq_, mob_perfcells, b_perfcells);

                if (param_.solve_welleq_initially_ && iterationIdx == 0) {
                    // solve the well equations as a pre-processing step
                    iter_report = solveWellEq(ebosSimulator, dt, well_state);
                }

                std::vector<ADB> cq_s;
                int nw = wells_->number_of_wells;
                int nc = numCells();
                Mat duneD (nw, nw, 9, 0.4, Mat::implicit);
                Mat duneB (nc, nw, 9, 0.4, Mat::implicit);
                Mat duneC (nw, nc, 9, 0.4, Mat::implicit);
                Mat duneA (nc, nc, 9, 0.4, Mat::implicit);
                BVector rhs(nc);

                computeWellFluxDense(ebosSimulator, cq_s, 4, duneA, duneB, duneC, duneD, rhs);
                //std::cout << cq_s[0] << std::endl;
                //std::cout << cq_s[1] << std::endl;
                //std::cout << cq_s[2] << std::endl;


                updatePerfPhaseRatesAndPressures(cq_s, well_state);
                addWellContributionToMassBalanceEq(cq_s,residual);

                addWellFluxEq(cq_s, dt, 4, residual, duneD);
                //std::cout << residual.well_flux_eq << std::endl;
                duneA.compress();
                duneB.compress();
                duneC.compress();
                duneD.compress();

                //std::cout << "B" << std::endl;
                //print(duneB);
                //std::cout << "C" << std::endl;
                //print(duneC);
                //print(duneD);


                V resWellEigen = residual.well_flux_eq.value();
                const int np = numPhases();
                BVector resWell(nw);
                for (int i = 0; i < nw; ++i){
                    for( int p = 0; p < np; ++p ) {
                        int idx = i + p * nw;
                        resWell[i][flowPhaseToEbosCompIdx(p)] = resWellEigen(idx);
                    }
                }

                resWell_ = resWell;
                rhs_ = rhs;
                duneB_ = duneB;
                duneC_ = duneC;
                localInvert(duneD);
                invDuneD_ = duneD;
                duneA_ = duneA;


                if (param_.compute_well_potentials_) {
                    //wellModel().computeWellPotentials(mob_perfcells, b_perfcells, state0, well_state);
                }
                return iter_report;
            }
            void localInvert(Mat& istlA) const {
                for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
                    for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                        //std::cout << (*col) << std::endl;
                        (*col).invert();
                    }
                }
            }

            void print(Mat& istlA) const {
                for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
                    for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                        std::cout << (*col) << std::endl;
                    }
                }
            }
            void addRhs(BVector& x, Mat& jac) const {
                assert(x.size() == rhs.size());
                x += rhs_;
                jac += duneA_;
            }

            void apply( Mat& A,
                        BVector& res) const {

                Mat BmultinvD;
                Mat duneA;

                Dune::matMultMat(BmultinvD, duneB_ , invDuneD_);
                Dune::matMultMat(duneA, BmultinvD, duneC_);
                A -= duneA;
                BmultinvD.mmv(resWell_, res);
            }

            void recoverVariable(const BVector& x, BVector& xw) const {
                BVector resWell = resWell_;
                duneC_.mmv(x, resWell);
                invDuneD_.mv(resWell, xw);
            }



            void addDuneMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor> eigen, Mat& dune, int p1, int p2) {
                //Mat dune( eigen.rows(), eigen.cols(), eigen.nonZeros(), 0.4,Mat::implicit) ;
                const int* ia = eigen.outerIndexPtr();
                const int* ja = eigen.innerIndexPtr();
                for (int row = 0; row < eigen.rows(); ++row) {
                    dune.entry(ia[row],ja[row])[p1][p2] = eigen.coeff(ia[row] ,ja[row]);
                }
            }

            int flowPhaseToEbosCompIdx( const int phaseIdx ) const
            {
                const int phaseToComp[ 3 ] = { FluidSystem::waterCompIdx, FluidSystem::oilCompIdx, FluidSystem::gasCompIdx };
                return phaseToComp[ phaseIdx ];
            }

            int flowToEbosPvIdx( const int flowPv ) const
            {
                const int flowToEbos[ 3 ] = {
                                              BlackoilIndices::pressureSwitchIdx,
                                              BlackoilIndices::waterSaturationIdx,
                                              BlackoilIndices::compositionSwitchIdx
                                            };
                return flowToEbos[ flowPv ];
            }

            int ebosCompToFlowPhaseIdx( const int compIdx ) const
            {
                const int compToPhase[ 3 ] = { Oil, Water, Gas };
                return compToPhase[ compIdx ];
            }



            void init(const BlackoilPropsAdInterface* fluid_arg,
                      const std::vector<bool>* active_arg,
                      const std::vector<PhasePresence>* pc_arg,
                      const VFPProperties*  vfp_properties_arg,
                      const double gravity_arg,
                      const Vector& depth_arg)
            {
                fluid_ = fluid_arg;
                active_ = active_arg;
                phase_condition_ = pc_arg;
                vfp_properties_ = vfp_properties_arg;
                gravity_ = gravity_arg;
                perf_cell_depth_ = subset(depth_arg, wellOps().well_cells); 
            }

            const WellOps& wellOps() const
            {
                return wops_;
            }


            int numPhases() const { return wells().number_of_phases; }

            int numCells() const { return pv_.size();}

            template<class WellState>
            void resetWellControlFromState(WellState xw) {
                const int        nw   = wells_->number_of_wells;
                for (int w = 0; w < nw; ++w) {
                    WellControls* wc = wells_->ctrls[w];
                    well_controls_set_current( wc, xw.currentControls()[w]);
                }
            }

            const Wells& wells() const
            {
                assert(wells_ != 0);
                return *(wells_);
            }

            const Wells* wellsPointer() const
            {
                return wells_;
            }

            /// return true if wells are available in the reservoir
            bool wellsActive() const
            {
                return wells_active_;
            }
            void setWellsActive(const bool wells_active)
            {
                wells_active_ = wells_active;
            }
            /// return true if wells are available on this process
            bool localWellsActive() const
            {
                return wells_ ? (wells_->number_of_wells > 0 ) : false;
            }

            int numWellVars() const
            {
                if ( !localWellsActive() )
                {
                    return 0;
                }

                // For each well, we have a bhp variable, and one flux per phase.
                const int nw = wells().number_of_wells;
                return (numPhases() + 1) * nw;
            }

            /// Density of each well perforation
            Vector& wellPerforationDensities() // mutable version kept for BlackoilMultisegmentModel
            {
                return well_perforation_densities_;
            }

            const Vector& wellPerforationDensities() const     {
                return well_perforation_densities_;
            }

            /// Diff to bhp for each well perforation.
            Vector& wellPerforationPressureDiffs()     { // mutable version kept for BlackoilMultisegmentModel
                return well_perforation_pressure_diffs_;
            }
            const Vector& wellPerforationPressureDiffs() const
            {
                return well_perforation_pressure_diffs_;
            }

            typedef DenseAd::Evaluation<double, /*size=*/3> Eval;
            EvalWell extendEval(Eval in) const {
                EvalWell out = 0.0;
                out.value = in.value;
                for(int i = 0;i<3;++i) {
                    out.derivatives[i] = in.derivatives[flowToEbosPvIdx(i)];
                }
                return out;
            }

            void
            setWellVariables(const WellState& xw) {
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
                for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                    for (int w = 0; w < nw; ++w) {
                        wellVariables_[w + nw*phaseIdx] = 0.0;
                        wellVariables_[w + nw*phaseIdx].value = xw.wellSolutions()[w + nw* phaseIdx];
                        wellVariables_[w + nw*phaseIdx].derivatives[np + phaseIdx] = 1.0;
                    }
                }
            }

            void print(EvalWell in) const {
                std::cout << in.value << std::endl;
                for (int i = 0; i < in.derivatives.size(); ++i) {
                    std::cout << in.derivatives[i] << std::endl;
                }
            }

            void
            computeAccumWells() {
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
                for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                    for (int w = 0; w < nw; ++w) {
                        F0_[w + nw * phaseIdx] = wellVolumeFraction(w,phaseIdx).value;
                    }
                }
            }

            EvalWell extractDenseAD(const ADB& data, int i, int j) const
            {
                EvalWell output = 0.0;
                output.value = data.value()[i];
                const int np = wells().number_of_phases;
                const std::vector<Opm::AutoDiffMatrix>& jac = data.derivative();
                //std::cout << jac.size() << std::endl;
                int numblocs = jac.size();
                for (int b = 0; b < numblocs; ++b) {
                    if (b < np) { // don't copy well blocks)
                        //std::cout << jac[b].coeff(i,j) << std::endl;
                        output.derivatives[b] = jac[b].coeff(i,j);
                    }
                }
                return output;
            }

            EvalWell extractDenseADWell(const ADB& data, int i) const
            {
                EvalWell output = 0.0;
                output.value = data.value()[i];
                const int nw = wells().number_of_wells;
                const int np = wells().number_of_phases;
                const std::vector<Opm::AutoDiffMatrix>& jac = data.derivative();
                //std::cout << jac.size() << std::endl;
                int numblocs = jac.size();
                for (int b = 0; b < np; ++b) {
                    output.derivatives[b+np] = jac[numblocs-1].coeff(i, b*nw + i);
                }
                return output;
            }

            const ADB convertToADB(const std::vector<EvalWell>& local, const std::vector<int>& well_cells, const int nc, const std::vector<int>& well_id, const int nw, const int numVars) const
        {
            typedef typename ADB::M  M;
            const int nLocal = local.size();
            typename ADB::V value( nLocal );
            //const int numVars = 5;
            const int np = wells().number_of_phases;

            std::vector<Eigen::SparseMatrix<double>> mat(np, Eigen::SparseMatrix<double>(nLocal,nc));
            Eigen::SparseMatrix<double> matFlux(nLocal,np*nw);
            Eigen::SparseMatrix<double> matBHP(nLocal,nw);

            for( int i=0; i<nLocal; ++i )
            {
                value[ i ] = local[ i ].value;
                for( int d=0; d<np; ++d ) {
                    //std::cout << i << " " <<d << " "<<local[i].derivatives[d] << std::endl;
                    mat[d].insert(i, well_cells[i]) = local[i].derivatives[d];
                }

                for (int phase = 0; phase < np; ++phase) {
                    //std::cout << "well: "<< i << " " << phase << " " << local[i].derivatives[np + phase] << std::endl;

                    matFlux.insert(i, nw*phase + well_id[i]) = local[i].derivatives[np + phase];
                }
                //matBHP.insert(i,well_id[i]) = local[i].derivatives[2*np];
            }

            std::vector< M > jacs( numVars );
            if (numVars == 4) {
                for( int d=0; d<np; ++d ) {
                    //Eigen::DiagonalMatrix<double>(deri[d]);
                    jacs[ d ] = M(mat[d]);
                }

                jacs[3] = M(matFlux);
                //jacs[4] = M(matBHP);
            }
            else if (numVars == 1) {
                jacs[0] = M(matFlux);
                //jacs[1] = M(matBHP);
            }
            //std::cout << numVars << std::endl;

            return ADB::function( std::move( value ), std::move( jacs ));
        }


            template <class WellState>
            void updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                                  WellState& xw) const
            {
                if ( !localWellsActive() )
                {
                    // If there are no wells in the subdomain of the proces then
                    // cq_s has zero size and will cause a segmentation fault below.
                    return;
                }

                // Update the perforation phase rates (used to calculate the pressure drop in the wellbore).
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
                const int nperf = wells().well_connpos[nw];

                V cq = superset(cq_s[0].value(), Span(nperf, np, 0), nperf*np);
                for (int phase = 1; phase < np; ++phase) {
                    cq += superset(cq_s[phase].value(), Span(nperf, np, phase), nperf*np);
                }
                xw.perfPhaseRates().assign(cq.data(), cq.data() + nperf*np);

                // Update the perforation pressures.
                const V& cdp = wellPerforationPressureDiffs();
                for (int w = 0; w < nw; ++w  ) {
                    for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                        xw.perfPress()[perf] = cdp[perf] + xw.bhp()[w];
                    }
                }

            }

            void
            addWellFluxEq(std::vector<ADB> cq_s,
                          const double dt,
                          const int numBlocks,
                          LinearisedBlackoilResidual& residual,
                          Mat& duneD)
            {

                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;

                double volume = 0.002831684659200; // 0.1 cu ft;
                //std::vector<ADB> F = wellVolumeFractions(state);
                //std::cout << F0_[0] << std::endl;
                //std::cout << F[0] << std::endl;
                //std::cout << "før Ebos" <<residual_.well_flux_eq << std::endl;
                ADB qs = ADB::constant(ADB::V::Zero(np*nw));
                for (int p = 0; p < np; ++p) {

                    std::vector<EvalWell> res_vec(nw);
                    for (int w = 0; w < nw; ++w) {

                        EvalWell res = (wellVolumeFraction(w, p) - F0_[w + nw*p]) * volume / dt;
                        res += getQs(w, p);
                        //for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                        //    res -= cq_s[perf*np + p];
                        //}
                        res_vec[w] = res;
                        for (int i = 0; i < np; ++i) {
                            duneD.entry(w,w)[flowPhaseToEbosCompIdx(p)][flowToEbosPvIdx(i)] += res.derivatives[i+3];
                        }
                    }

                    ADB tmp = convertToADBWell(res_vec, numBlocks);
                    qs += superset(tmp,Span(nw,1,p*nw), nw*np);
                }

                //wellModel().convertToADB(res_vec, well_cells, nc, well_id, nw, numBlocks);
                //ADB qs = state.qs;
                for (int phase = 0; phase < np; ++phase) {
                    qs -= superset(wellOps().p2w * cq_s[phase], Span(nw, 1, phase*nw), nw*np);
                    //qs += superset((F[phase]-F0_[phase]) * vol_dt, Span(nw,1,phase*nw), nw*np);
                }

                residual.well_flux_eq = qs;
                //std::cout << "etter Ebos" << residual_.well_flux_eq << std::endl;

            }
            const AutoDiffBlock<double> convertToADBWell(const std::vector<EvalWell>& local, const int numVars) const
            {
                typedef typename ADB::M  M;
                const int nLocal = local.size();
                typename ADB::V value( nLocal );
                //const int numVars = 5;
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;

                Eigen::SparseMatrix<double> matFlux(nLocal,np*nw);

                for( int i=0; i<nLocal; ++i )
                {
                    value[ i ] = local[ i ].value;
                    for (int phase = 0; phase < np; ++phase) {
                        matFlux.insert(i, nw*phase + i) = local[i].derivatives[np + phase];
                    }
                }

                std::vector< M > jacs( numVars );
                if (numVars == 4) {
                    for( int d=0; d<np; ++d ) {
                        //Eigen::DiagonalMatrix<double>(deri[d]);
                        //jacs[ d ] = M(mat[d]);
                    }

                    jacs[3] = M(matFlux);
                    //jacs[4] = M(matBHP);
                }
                else if (numVars == 1) {
                    jacs[0] = M(matFlux);
                    //jacs[1] = M(matBHP);
                }
                //std::cout << numVars << std::endl;

                return ADB::function( std::move( value ), std::move( jacs ));
            }


            template<typename Simulator>
            void
            computeWellFluxDense(const Simulator& ebosSimulator,
                                 std::vector<ADB>& cq_s,
                                 const int numBlocks,
                                 Mat& duneA,
                                 Mat& duneB,
                                 Mat& duneC,
                                 Mat& duneD,
                                 BVector& rhs) const
            {
                if( ! localWellsActive() ) return ;
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
                const int nperf = wells().well_connpos[nw];
                const Opm::PhaseUsage& pu = fluid_->phaseUsage();
                V Tw = Eigen::Map<const V>(wells().WI, nperf);
                const std::vector<int>& well_cells = wellOps().well_cells;
                std::vector<int> well_id(nperf);
                std::vector<std::vector<EvalWell>> cq_s_dense(np, std::vector<EvalWell>(nperf,0.0));


                // pressure diffs computed already (once per step, not changing per iteration)
                const V& cdp = wellPerforationPressureDiffs();

                for (int w = 0; w < nw; ++w) {

                    EvalWell bhp = getBhp(w);

                    // TODO: fix for 2-phase case
                    std::vector<EvalWell> cmix_s(np,0.0);
                    for (int phase = 0; phase < np; ++phase) {
                        //int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                        cmix_s[phase] = wellVolumeFraction(w,phase);
                    }

                    for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                        const int cell_idx = well_cells[perf];
                        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                        const auto& fs = intQuants.fluidState();
                        well_id[perf] = w;
                        EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
                        EvalWell rs = extendEval(fs.Rs());
                        EvalWell rv = extendEval(fs.Rv());
                        std::vector<EvalWell> b_perfcells_dense(np, 0.0);
                        std::vector<EvalWell> mob_perfcells_dense(np, 0.0);
                        for (int phase = 0; phase < np; ++phase) {
                            int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                            b_perfcells_dense[phase] = extendEval(fs.invB(ebosPhaseIdx));
                            mob_perfcells_dense[phase] = extendEval(intQuants.mobility(ebosPhaseIdx));
                        }

                        // Pressure drawdown (also used to determine direction of flow)
                        EvalWell well_pressure = bhp + cdp[perf];
                        EvalWell drawdown = pressure - well_pressure;

                        // injection perforations
                        if ( drawdown.value > 0 )  {

                            //Do nothing if crossflow is not allowed
                            if (!wells().allow_cf[w] && wells().type[w] == INJECTOR)
                                continue;
                            // compute phase volumetric rates at standard conditions
                            std::vector<EvalWell> cq_ps(np, 0.0);
                            for (int phase = 0; phase < np; ++phase) {
                                const EvalWell cq_p = - Tw[perf] * (mob_perfcells_dense[phase] * drawdown);
                                cq_ps[phase] = b_perfcells_dense[phase] * cq_p;
                            }


                            if ((*active_)[Oil] && (*active_)[Gas]) {
                                const int oilpos = pu.phase_pos[Oil];
                                const int gaspos = pu.phase_pos[Gas];
                                const EvalWell cq_psOil = cq_ps[oilpos];
                                const EvalWell cq_psGas = cq_ps[gaspos];
                                cq_ps[gaspos] += rs * cq_psOil;
                                cq_ps[oilpos] += rv * cq_psGas;
                            }


                            // map to ADB
                            for (int phase = 0; phase < np; ++phase) {
                                cq_s_dense[phase][perf] = cq_ps[phase];

                            }

                        } else {
                            //Do nothing if crossflow is not allowed
                            if (!wells().allow_cf[w] && wells().type[w] == PRODUCER)
                                continue;

                            // Using total mobilities
                            EvalWell total_mob_dense = mob_perfcells_dense[0];
                            for (int phase = 1; phase < np; ++phase) {
                                total_mob_dense += mob_perfcells_dense[phase];
                            }
                            // injection perforations total volume rates
                            const EvalWell cqt_i = - Tw[perf] * (total_mob_dense * drawdown);

                            // compute volume ratio between connection at standard conditions
                            EvalWell volumeRatio = 0.0;
                            if ((*active_)[Water]) {
                                const int watpos = pu.phase_pos[Water];
                                volumeRatio += cmix_s[watpos] / b_perfcells_dense[watpos];
                            }

                            if ((*active_)[Oil] && (*active_)[Gas]) {
                                EvalWell well_temperature = extendEval(fs.temperature(FluidSystem::oilPhaseIdx));
                                EvalWell rsSatEval = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), well_temperature, well_pressure);
                                EvalWell rvSatEval = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), well_temperature, well_pressure);

                                const int oilpos = pu.phase_pos[Oil];
                                const int gaspos = pu.phase_pos[Gas];
                                EvalWell rvPerf = 0.0;
                                if (cmix_s[gaspos] > 0)
                                    rvPerf = cmix_s[oilpos] / cmix_s[gaspos];

                                if (rvPerf.value > rvSatEval.value) {
                                    rvPerf = rvSatEval;
                                    //rvPerf.value = rvSatEval.value;
                                }

                                EvalWell rsPerf = 0.0;
                                if (cmix_s[oilpos] > 0)
                                    rsPerf = cmix_s[gaspos] / cmix_s[oilpos];

                                if (rsPerf.value > rsSatEval.value) {
                                    //rsPerf = 0.0;
                                    rsPerf= rsSatEval;
                                }

                                // Incorporate RS/RV factors if both oil and gas active
                                const EvalWell d = 1.0 - rvPerf * rsPerf;

                                const EvalWell tmp_oil = (cmix_s[oilpos] - rvPerf * cmix_s[gaspos]) / d;
                                //std::cout << "tmp_oil " <<tmp_oil << std::endl;
                                volumeRatio += tmp_oil / b_perfcells_dense[oilpos];

                                const EvalWell tmp_gas = (cmix_s[gaspos] - rsPerf * cmix_s[oilpos]) / d;
                                //std::cout << "tmp_gas " <<tmp_gas << std::endl;
                                volumeRatio += tmp_gas / b_perfcells_dense[gaspos];
                            }
                            else {
                                if ((*active_)[Oil]) {
                                    const int oilpos = pu.phase_pos[Oil];
                                    volumeRatio += cmix_s[oilpos] / b_perfcells_dense[oilpos];
                                }
                                if ((*active_)[Gas]) {
                                    const int gaspos = pu.phase_pos[Gas];
                                    volumeRatio += cmix_s[gaspos] / b_perfcells_dense[gaspos];
                                }
                            }
                            // injecting connections total volumerates at standard conditions
                            EvalWell cqt_is = cqt_i/volumeRatio;
                            //std::cout << "volrat " << volumeRatio << " " << volrat_perf_[perf] << std::endl;
                            for (int phase = 0; phase < np; ++phase) {
                                cq_s_dense[phase][perf] = cmix_s[phase] * cqt_is; // * b_perfcells_dense[phase];                                
                            }
                        }

                        for (int p1 = 0; p1 < np; ++p1) {
                            EvalWell tmp = cq_s_dense[p1][perf];
                            for (int p2 = 0; p2 < np; ++p2) {
                                duneC.entry(w, cell_idx)[flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= tmp.derivatives[p2];
                                duneB.entry(cell_idx, w)[flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= tmp.derivatives[p2+3];
                                duneA.entry(cell_idx, cell_idx)[flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= tmp.derivatives[p2];
                                duneD.entry(w, w)[flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= tmp.derivatives[p2+3];
                            }
                            rhs[cell_idx][flowPhaseToEbosCompIdx(p1)] -= tmp.value;
                        }
                    }

                }
                cq_s.resize(np, ADB::null());
                for (int phase = 0; phase < np; ++phase) {
                    cq_s[phase] = convertToADB(cq_s_dense[phase], well_cells, numCells(), well_id, nw, numBlocks);
                }


            }

            template <typename Simulator>
            IterationReport solveWellEq(Simulator& ebosSimulator,
                                        const double dt,
                                        WellState& well_state)
            {
                const int np = wells().number_of_phases;
                std::vector<ADB> cq_s(np, ADB::null());
                WellState well_state0 = well_state;

                LinearisedBlackoilResidual residual( { std::vector<ADB>(3, ADB::null()),
                                    ADB::null(),
                                    ADB::null(),
                                    { 1.1169, 1.0031, 0.0031 }, // the default magic numbers
                                    false } );

                int it  = 0;
                bool converged;
                do {
                    // bhp and Q for the wells
                    int nw = wells().number_of_wells;
                    int nc = numCells();
                    Mat duneDlo( nw, nw, 9, 0.4,Mat::implicit);
                    Mat duneBlo( nc, nw, 9, 0.4, Mat::implicit);
                    Mat duneClo( nw, nc, 9, 0.4, Mat::implicit);
                    Mat duneAlo( nc, nc, 9, 0.4, Mat::implicit);
                    BVector rhslo(nc);
                    computeWellFluxDense(ebosSimulator, cq_s, 1, duneAlo, duneBlo, duneClo, duneDlo, rhslo);
                    updatePerfPhaseRatesAndPressures(cq_s, well_state);
                    addWellFluxEq(cq_s, dt, 1, residual, duneDlo);
                    V resWellEigen = residual.well_flux_eq.value();
                    //std::cout << "resWellEigen " << resWellEigen << std::endl;
                    BVector resWell(nw);
                    for (int i = 0; i < nw; ++i){
                        for( int p = 0; p < np; ++p ) {
                            int idx = i + p * nw;
                            resWell[i][flowPhaseToEbosCompIdx(p)] = resWellEigen(idx);
                        }
                    }
                    duneDlo.compress();
                    //print(duneDlo);
                    localInvert(duneDlo);

                    resWell_ = resWell;
                    converged = getWellConvergence(ebosSimulator, it);

                    if (converged) {
                        break;
                    }

                    ++it;
                    if( localWellsActive() )
                    {
//                        std::vector<ADB> eqs;
//                        eqs.reserve(1);
//                        eqs.push_back(residual.well_flux_eq);
//                        //eqs.push_back(residual_.well_eq);
//                        ADB total_residual = vertcatCollapseJacs(eqs);
//                        const std::vector<M>& Jn = total_residual.derivative();
//                        typedef Eigen::SparseMatrix<double> Sp;
//                        Sp Jn0;
//                        Jn[0].toSparse(Jn0);
//                        std::cout << Jn0 << std::endl;
//                        const Eigen::SparseLU< Sp > solver(Jn0);
//                        ADB::V total_residual_v = total_residual.value();
//                        std::cout << "tot res " <<total_residual_v << std::endl;
//                        const Eigen::VectorXd& dx = solver.solve(total_residual_v.matrix());
                        BVector dx_new (nw);
                        duneDlo.mv(resWell_, dx_new);
//                        std::cout << "hei" << std::endl;
//                        Sp eye(nw*np,nw*np);
//                        eye.setIdentity();

//                        Sp invD = solver.solve(eye);

//                        std::cout << invD << std::endl;
//                        print(duneDlo);

                        V dx_new_eigen(np*nw);
                        for( int p=0; p<np; ++p) {

                            for (int w = 0; w < nw; ++w) {
                                int idx = w + nw*p;

                                dx_new_eigen(idx) = dx_new[w][flowPhaseToEbosCompIdx(p)];
                            }
                        }
//                        std::cout << "new " <<dx_new_eigen << std::endl;
//                        std::cout << "old " <<dx << std::endl;


                        assert(dx.size() == total_residual_v.size());
                        updateWellState(dx_new_eigen.array(), well_state);
                        updateWellControls(well_state);
                        setWellVariables(well_state);
                    }
                } while (it < 15);

                if (!converged) {
                    well_state = well_state0;
                }

                const bool failed = false; // Not needed in this method.
                const int linear_iters = 0; // Not needed in this method
                return IterationReport{failed, converged, linear_iters, it};
            }

            std::vector<double> residual() {

                const int np = numPhases();
                const int nw = wells().number_of_wells;
                std::vector<double> res(np*nw);
                for( int p=0; p<np; ++p) {
                    for (int i = 0; i < nw; ++i) {
                        int idx = i + nw*p;
                        res[idx] = resWell_[i][flowPhaseToEbosCompIdx(p)];
                    }
                }
                return res;
            }


            template <typename Simulator>
            bool getWellConvergence(Simulator& ebosSimulator,
                                    const int iteration)
            {
                const int np = numPhases();
                const int nc = numCells();
                const double tol_wells = param_.tolerance_wells_;
                const double maxResidualAllowed = param_.max_residual_allowed_;

                std::vector<double> R_sum(np);
                std::vector<double> B_avg(np);
                std::vector<double> maxCoeff(np);
                std::vector<double> maxNormWell(np);
                Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> B(nc, np);
                Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> R(nc, np);
                Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> tempV(nc, np);

                for ( int idx = 0; idx < np; ++idx )
                {
                    V b(nc);
                    for (int cell_idx = 0; cell_idx < nc; ++cell_idx) {
                        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                        const auto& fs = intQuants.fluidState();

                        int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(idx);

                        b[cell_idx] = 1 / fs.invB(ebosPhaseIdx).value;
                    }
                    B.col(idx) = b;
                }

                detail::convergenceReduction(B, tempV, R, R_sum, maxCoeff, B_avg, maxNormWell, nc, np, pv_, residual());

                std::vector<double> well_flux_residual(np);
                bool converged_Well = true;
                // Finish computation
                for ( int idx = 0; idx < np; ++idx )
                {
                    well_flux_residual[idx] = B_avg[idx] * maxNormWell[idx];
                    converged_Well = converged_Well && (well_flux_residual[idx] < tol_wells);
                }

                // if one of the residuals is NaN, throw exception, so that the solver can be restarted
                for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                    const auto& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));

                    if (std::isnan(well_flux_residual[phaseIdx])) {
                        OPM_THROW(Opm::NumericalProblem, "NaN residual for phase " << phaseName);
                    }
                    if (well_flux_residual[phaseIdx] > maxResidualAllowed) {
                        OPM_THROW(Opm::NumericalProblem, "Too large residual for phase " << phaseName);
                    }
                }

                if ( terminal_output_ )
                {
                    // Only rank 0 does print to std::cout
                    if (iteration == 0) {
                        std::string msg;
                        msg = "Iter";
                        for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                            const std::string& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));
                            msg += "  W-FLUX(" + phaseName + ")";
                        }
                        OpmLog::note(msg);
                    }
                    std::ostringstream ss;
                    const std::streamsize oprec = ss.precision(3);
                    const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                    ss << std::setw(4) << iteration;
                    for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                        ss << std::setw(11) << well_flux_residual[phaseIdx];
                    }
                    ss.precision(oprec);
                    ss.flags(oflags);
                    OpmLog::note(ss.str());
                }
                return converged_Well;
            }



            template<typename Simulator>
            void
            computeWellConnectionPressures(const Simulator& ebosSimulator,
                                           const WellState& xw)
            {
                if( ! localWellsActive() ) return ;
                // 1. Compute properties required by computeConnectionPressureDelta().
                //    Note that some of the complexity of this part is due to the function
                //    taking std::vector<double> arguments, and not Eigen objects.
                std::vector<double> b_perf;
                std::vector<double> rsmax_perf;
                std::vector<double> rvmax_perf;
                std::vector<double> surf_dens_perf;
                computePropertiesForWellConnectionPressures(ebosSimulator, xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

                const V& pdepth = perf_cell_depth_;
                const int nperf = wells().well_connpos[wells().number_of_wells];
                const std::vector<double> depth_perf(pdepth.data(), pdepth.data() + nperf);

                computeWellConnectionDensitesPressures(xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf, depth_perf, gravity_);

            }

            template<typename Simulator, class WellState>
            void
            computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                        const WellState& xw,
                                                        std::vector<double>& b_perf,
                                                        std::vector<double>& rsmax_perf,
                                                        std::vector<double>& rvmax_perf,
                                                        std::vector<double>& surf_dens_perf)
            {
                const int nperf = wells().well_connpos[wells().number_of_wells];
                const int nw = wells().number_of_wells;
                const std::vector<int>& well_cells = wellOps().well_cells;
                const PhaseUsage& pu = fluid_->phaseUsage();
                const int np = fluid_->numPhases();
                b_perf.resize(nperf*np);
                rsmax_perf.resize(nperf);
                rvmax_perf.resize(nperf);

                std::vector<PhasePresence> perf_cond(nperf);
                for (int perf = 0; perf < nperf; ++perf) {
                    perf_cond[perf] = (*phase_condition_)[well_cells[perf]];
                }

                // Compute the average pressure in each well block
                const V perf_press = Eigen::Map<const V>(xw.perfPress().data(), nperf);
                //V avg_press = perf_press*0;
                for (int w = 0; w < nw; ++w) {
                    for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                        const int cell_idx = well_cells[perf];
                        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                        const auto& fs = intQuants.fluidState();

                        const double p_above = perf == wells().well_connpos[w] ? xw.bhp()[w] : perf_press[perf - 1];
                        const double p_avg = (perf_press[perf] + p_above)/2;
                        double temperature = fs.temperature(FluidSystem::oilPhaseIdx).value;

                        if (pu.phase_used[BlackoilPhases::Aqua]) {
                            b_perf[ pu.phase_pos[BlackoilPhases::Aqua] + perf * pu.num_phases] =
                                    FluidSystem::waterPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                        }


                        if (pu.phase_used[BlackoilPhases::Vapour]) {
                            int gaspos = pu.phase_pos[BlackoilPhases::Vapour] + perf * pu.num_phases;
                            if (perf_cond[perf].hasFreeOil()) {
                                b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                            }
                            else {
                                double rv = fs.Rv().value;
                                b_perf[gaspos] = FluidSystem::gasPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, rv);
                            }
                        }

                        if (pu.phase_used[BlackoilPhases::Liquid]) {
                            int oilpos = pu.phase_pos[BlackoilPhases::Liquid] + perf * pu.num_phases;
                            if (perf_cond[perf].hasFreeGas()) {
                                b_perf[oilpos] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                            }
                            else {
                                double rs = fs.Rs().value;
                                b_perf[oilpos] = FluidSystem::oilPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, rs);
                            }
                        }

                        if (pu.phase_used[BlackoilPhases::Liquid] && pu.phase_used[BlackoilPhases::Vapour]) {
                            rsmax_perf[perf] = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), temperature, p_avg);
                            rvmax_perf[perf] = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), temperature, p_avg);
                        }

                    }
                }

                // Surface density.
                // The compute density segment wants the surface densities as
                // an np * number of wells cells array
                V rho = superset(fluid_->surfaceDensity(0 , well_cells), Span(nperf, pu.num_phases, 0), nperf*pu.num_phases);
                for (int phase = 1; phase < pu.num_phases; ++phase) {
                    rho += superset(fluid_->surfaceDensity(phase , well_cells), Span(nperf, pu.num_phases, phase), nperf*pu.num_phases);
                }
                surf_dens_perf.assign(rho.data(), rho.data() + nperf * pu.num_phases);

            }


            void
            addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s,
                                               LinearisedBlackoilResidual& residual)
            {
                if ( !localWellsActive() )
                {
                    // If there are no wells in the subdomain of the proces then
                    // cq_s has zero size and will cause a segmentation fault below.
                    return;
                }

                // Add well contributions to mass balance equations
                const int np = numPhases();
                for (int phase = 0; phase < np; ++phase) {
                    residual.material_balance_eq[phase] -= superset(cq_s[phase], wellOps().well_cells, numCells());
                }
            }

            template <class WellState>
            void updateWellState(const Vector& dwells,
                                 WellState& well_state)
            {
                if( localWellsActive() )
                {
                    const int np = wells().number_of_phases;
                    const int nw = wells().number_of_wells;

                    // Extract parts of dwells corresponding to each part.
                    int varstart = 0;
                    const Vector dxvar_well = subset(dwells, Span(np*nw, 1, varstart));
                    //const Vector dqs = subset(dwells, Span(np*nw, 1, varstart));
                    varstart += dxvar_well.size();
                    //const Vector dbhp = subset(dwells, Span(nw, 1, varstart));
                    //varstart += dbhp.size();
                    assert(varstart == dwells.size());

                    // Qs update.
                    // Since we need to update the wellrates, that are ordered by wells,
                    // from dqs which are ordered by phase, the simplest is to compute
                    // dwr, which is the data from dqs but ordered by wells.
                    const Vector xvar_well_old = Eigen::Map<const Vector>(&well_state.wellSolutions()[0], nw*np);

                    double dFLimit = 0.2;
                    double dBHPLimit = 2;
                    double dTotalRateLimit = 0.5;
                    //std::cout << "dxvar_well "<<dxvar_well << std::endl;


                    for (int w = 0; w < nw; ++w) {
                        const WellControls* wc = wells().ctrls[w];
                        // The current control in the well state overrides
                        // the current control set in the Wells struct, which
                        // is instead treated as a default.
                        const int current = well_state.currentControls()[w];
                        const double target_rate = well_controls_iget_target(wc, current);
                        const double* distr = well_controls_iget_distr(wc, current);
                        std::vector<double> F(np,0.0);
                        const int sign2 = dxvar_well[nw + w] > 0 ? 1: -1;
                        const double dx2_limited = sign2 * std::min(std::abs(dxvar_well[nw + w]),dFLimit);
                        well_state.wellSolutions()[nw + w] = xvar_well_old[nw + w] - dx2_limited;
                        const int sign3 = dxvar_well[2*nw + w] > 0 ? 1: -1;
                        const double dx3_limited = sign3 * std::min(std::abs(dxvar_well[2*nw + w]),dFLimit);
                        well_state.wellSolutions()[2*nw + w] = xvar_well_old[2*nw + w] - dx3_limited;
                        F[Water] = well_state.wellSolutions()[nw + w];
                        F[Gas] = well_state.wellSolutions()[2*nw + w];
                        F[Oil] = 1.0 - F[Water] - F[Gas];

        //                const double dFw = dxvar_well[nw + w];
        //                const double dFg = dxvar_well[nw*2 + w];
        //                const double dFo = - dFw - dFg;
        //                //std::cout << w << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;
        //                double step = dFLimit / std::max(std::abs(dFw),std::max(std::abs(dFg),std::abs(dFo))); //)) / dFLimit;
        //                step = std::min(step, 1.0);
        //                //std::cout << step << std::endl;
        //                F[Water] = xvar_well_old[nw + w] - step*dFw;
        //                F[Gas] = xvar_well_old[2*nw + w] - step*dFg;
        //                F[Oil] = (1.0 - xvar_well_old[2*nw + w] - xvar_well_old[nw + w]) - step * dFo;

                        if (F[Water] < 0.0) {
                            F[Gas] /= (1.0 - F[Water]);
                            F[Oil] /= (1.0 - F[Water]);
                            F[Water] = 0.0;
                        }
                        if (F[Gas] < 0.0) {
                            F[Water] /= (1.0 - F[Gas]);
                            F[Oil] /= (1.0 - F[Gas]);
                            F[Gas] = 0.0;
                        }
                        if (F[Oil] < 0.0) {
                            F[Water] /= (1.0 - F[Oil]);
                            F[Gas] /= (1.0 - F[Oil]);
                            F[Oil] = 0.0;
                        }
                        well_state.wellSolutions()[nw + w] = F[Water];
                        well_state.wellSolutions()[2*nw + w] = F[Gas];

                        //std::cout << wells().name[w] << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;


                        std::vector<double> g = {1,1,0.01};

                        if (well_controls_iget_type(wc, current) == RESERVOIR_RATE) {
                            for (int p = 0; p < np; ++p) {
                                F[p] /= distr[p];
                            }
                        } else {
                            for (int p = 0; p < np; ++p) {
                                F[p] /= g[p];
                            }
                        }
                        //std::cout << w << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;






        //                const double dFw = dxvar_well[nw + w];
        //                const double dFg = dxvar_well[nw*2 + w];
        //                const double dFo = - dFw - dFg;
                            //std::cout << w << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;
        //                double step = dFLimit / std::max(std::abs(dFw),std::max(std::abs(dFg),std::abs(dFo))); //)) / dFLimit;
        //                step = std::min(step, 1.0);
        //                std::cout << step << std::endl;
        //                F[Water] = xvar_well_old[nw + w] - step*dFw;
        //                F[Gas] = xvar_well_old[2*nw + w] - step*dFg;
        //                F[Oil] = (1.0 - xvar_well_old[2*nw + w] - xvar_well_old[nw + w]) - step * dFo;
        //                double sumF = F[Water]+F[Gas]+F[Oil];
        //                F[Water] /= sumF;
        //                F[Gas] /= sumF;
        //                F[Oil] /= sumF;
        //                well_state.wellSolutions()[nw + w] = F[Water];
        //                well_state.wellSolutions()[2 * nw + w] = F[Gas];

                        switch (well_controls_iget_type(wc, current)) {
                        case BHP:
                        {
                            //const int sign1 = dxvar_well[w] > 0 ? 1: -1;
                            //const double dx1_limited = sign1 * std::min(std::abs(dxvar_well[w]),std::abs(xvar_well_old[w])*dTotalRateLimit);
                            well_state.wellSolutions()[w] = xvar_well_old[w] - dxvar_well[w];

                            switch (wells().type[w]) {
                            case INJECTOR:
                                for (int p = 0; p < np; ++p) {
                                    const double comp_frac = wells().comp_frac[np*w + p];
                                    //if (comp_frac > 0) {
                                    well_state.wellRates()[w*np + p] = comp_frac * well_state.wellSolutions()[w];
                                    //}

                                }
                                break;
                            case PRODUCER:
                                for (int p = 0; p < np; ++p) {
                                    well_state.wellRates()[w*np + p] = well_state.wellSolutions()[w] * F[p];
                                }
                                break;
                            }
                        }
                            break;
                        case SURFACE_RATE:
                        {
                            const int sign1 = dxvar_well[w] > 0 ? 1: -1;
                            const double dx1_limited = sign1 * std::min(std::abs(dxvar_well[w]),std::abs(xvar_well_old[w])*dBHPLimit);
                            well_state.wellSolutions()[w] = xvar_well_old[w] - dx1_limited;
                            //const int sign = (dxvar_well1[w] < 0) ? -1 : 1;
                            //well_state.bhp()[w] -= sign * std::min( std::abs(dxvar_well1[w]), std::abs(well_state.bhp()[w])*dpmaxrel) ;
                            well_state.bhp()[w] = well_state.wellSolutions()[w];

                            if (wells().type[w]==PRODUCER) {

                                double F_target = 0.0;
                                for (int p = 0; p < np; ++p) {
                                    F_target += wells().comp_frac[np*w + p] * F[p];
                                }
                                for (int p = 0; p < np; ++p) {
                                    //std::cout << F[p] << std::endl;
                                    well_state.wellRates()[np*w + p] = F[p] * target_rate /F_target;
                                }
                            } else {

                                for (int p = 0; p < np; ++p) {
                                    //std::cout << wells().comp_frac[np*w + p] << " " <<distr[p] << std::endl;
                                    well_state.wellRates()[w*np + p] = wells().comp_frac[np*w + p] * target_rate;
                                }
                            }


                        }
                            break;
                        case RESERVOIR_RATE: {
                            const int sign1 = dxvar_well[w] > 0 ? 1: -1;
                            const double dx1_limited = sign1 * std::min(std::abs(dxvar_well[w]),std::abs(xvar_well_old[w])*dBHPLimit);
                            well_state.wellSolutions()[w] = xvar_well_old[w] - dx1_limited;
                            //const int sign = (dxvar_well1[w] < 0) ? -1 : 1;
                            //well_state.bhp()[w] -= sign * std::min( std::abs(dxvar_well1[w]), std::abs(well_state.bhp()[w])*dpmaxrel) ;
                            well_state.bhp()[w] = well_state.wellSolutions()[w];
                            for (int p = 0; p < np; ++p) {
                                well_state.wellRates()[np*w + p] = F[p] * target_rate;
                            }
                        }
                            break;

                        }
                    }





                    const Opm::PhaseUsage& pu = fluid_->phaseUsage();
                    //Loop over all wells
        #pragma omp parallel for schedule(static)
                    for (int w = 0; w < nw; ++w) {
                        const WellControls* wc = wells().ctrls[w];
                        const int nwc = well_controls_get_num(wc);
                        //Loop over all controls until we find a THP control
                        //that specifies what we need...
                        //Will only update THP for wells with THP control
                        for (int ctrl_index=0; ctrl_index < nwc; ++ctrl_index) {
                            if (well_controls_iget_type(wc, ctrl_index) == THP) {
                                double aqua = 0.0;
                                double liquid = 0.0;
                                double vapour = 0.0;

                                if ((*active_)[ Water ]) {
                                    aqua = well_state.wellRates()[w*np + pu.phase_pos[ Water ] ];
                                }
                                if ((*active_)[ Oil ]) {
                                    liquid = well_state.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                                }
                                if ((*active_)[ Gas ]) {
                                    vapour = well_state.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                                }

                                double alq = well_controls_iget_alq(wc, ctrl_index);
                                int table_id = well_controls_iget_vfp(wc, ctrl_index);

                                const WellType& well_type = wells().type[w];
                                if (well_type == INJECTOR) {
                                    double dp = wellhelpers::computeHydrostaticCorrection(
                                            wells(), w, vfp_properties_->getInj()->getTable(table_id)->getDatumDepth(),
                                            wellPerforationDensities(), gravity_);

                                    well_state.thp()[w] = vfp_properties_->getInj()->thp(table_id, aqua, liquid, vapour, well_state.bhp()[w] + dp);
                                }
                                else if (well_type == PRODUCER) {
                                    double dp = wellhelpers::computeHydrostaticCorrection(
                                            wells(), w, vfp_properties_->getProd()->getTable(table_id)->getDatumDepth(),
                                            wellPerforationDensities(), gravity_);

                                    well_state.thp()[w] = vfp_properties_->getProd()->thp(table_id, aqua, liquid, vapour, well_state.bhp()[w] + dp, alq);
                                }
                                else {
                                    OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
                                }

                                //Assume only one THP control specified for each well
                                break;
                            }
                        }
                    }
                }
            }

            template <class WellState>
            void updateWellControls(WellState& xw)
            {
                if( !localWellsActive() ) return ;

                std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };
                // Find, for each well, if any constraints are broken. If so,
                // switch control to first broken constraint.
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
        #pragma omp parallel for schedule(dynamic)
                for (int w = 0; w < nw; ++w) {
                    WellControls* wc = wells().ctrls[w];
                    // The current control in the well state overrides
                    // the current control set in the Wells struct, which
                    // is instead treated as a default.
                    int current = xw.currentControls()[w];
                    // Loop over all controls except the current one, and also
                    // skip any RESERVOIR_RATE controls, since we cannot
                    // handle those.
                    const int nwc = well_controls_get_num(wc);
                    int ctrl_index = 0;
                    for (; ctrl_index < nwc; ++ctrl_index) {
                        if (ctrl_index == current) {
                            // This is the currently used control, so it is
                            // used as an equation. So this is not used as an
                            // inequality constraint, and therefore skipped.
                            continue;
                        }
                        if (wellhelpers::constraintBroken(
                                xw.bhp(), xw.thp(), xw.wellRates(),
                                w, np, wells().type[w], wc, ctrl_index)) {
                            // ctrl_index will be the index of the broken constraint after the loop.
                            break;
                        }
                    }
                    if (ctrl_index != nwc) {
                        // Constraint number ctrl_index was broken, switch to it.
                        // We disregard terminal_ouput here as with it only messages
                        // for wells on one process will be printed.
                        std::ostringstream ss;
                        ss << "Switching control mode for well " << wells().name[w]
                              << " from " << modestring[well_controls_iget_type(wc, current)]
                              << " to " << modestring[well_controls_iget_type(wc, ctrl_index)] << std::endl;
                        OpmLog::info(ss.str());
                        xw.currentControls()[w] = ctrl_index;
                        current = xw.currentControls()[w];
                        well_controls_set_current( wc, current);



                        // Updating well state and primary variables if constraint is broken

                        // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
                        const double target = well_controls_iget_target(wc, current);
                        const double* distr = well_controls_iget_distr(wc, current);
                        switch (well_controls_iget_type(wc, current)) {
                        case BHP:
                            xw.bhp()[w] = target;
                            break;

                        case THP: {
                            double aqua = 0.0;
                            double liquid = 0.0;
                            double vapour = 0.0;

                            const Opm::PhaseUsage& pu = fluid_->phaseUsage();

                            if ((*active_)[ Water ]) {
                                aqua = xw.wellRates()[w*np + pu.phase_pos[ Water ] ];
                            }
                            if ((*active_)[ Oil ]) {
                                liquid = xw.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                            }
                            if ((*active_)[ Gas ]) {
                                vapour = xw.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                            }

                            const int vfp        = well_controls_iget_vfp(wc, current);
                            const double& thp    = well_controls_iget_target(wc, current);
                            const double& alq    = well_controls_iget_alq(wc, current);

                            //Set *BHP* target by calculating bhp from THP
                            const WellType& well_type = wells().type[w];

                            if (well_type == INJECTOR) {
                                double dp = wellhelpers::computeHydrostaticCorrection(
                                            wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                            wellPerforationDensities(), gravity_);

                                xw.bhp()[w] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                            }
                            else if (well_type == PRODUCER) {
                                double dp = wellhelpers::computeHydrostaticCorrection(
                                            wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                            wellPerforationDensities(), gravity_);

                                xw.bhp()[w] = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                            }
                            else {
                                OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                            }
                            break;
                        }

                        case RESERVOIR_RATE:
                            // No direct change to any observable quantity at
                            // surface condition.  In this case, use existing
                            // flow rates as initial conditions as reservoir
                            // rate acts only in aggregate.
                            break;

                        case SURFACE_RATE:
                            // assign target value as initial guess for injectors and
                            // single phase producers (orat, grat, wrat)
                            const WellType& well_type = wells().type[w];
                            if (well_type == INJECTOR) {
                                for (int phase = 0; phase < np; ++phase) {
                                    const double& compi = wells().comp_frac[np * w + phase];
                                    //if (compi > 0.0) {
                                        xw.wellRates()[np*w + phase] = target * compi;
                                    //}
                                }
                            } else if (well_type == PRODUCER) {

                                // only set target as initial rates for single phase
                                // producers. (orat, grat and wrat, and not lrat)
                                // lrat will result in numPhasesWithTargetsUnderThisControl == 2
                                int numPhasesWithTargetsUnderThisControl = 0;
                                for (int phase = 0; phase < np; ++phase) {
                                    if (distr[phase] > 0.0) {
                                        numPhasesWithTargetsUnderThisControl += 1;
                                    }
                                }
                                for (int phase = 0; phase < np; ++phase) {
                                    if (distr[phase] > 0.0 && numPhasesWithTargetsUnderThisControl < 2 ) {
                                        xw.wellRates()[np*w + phase] = target * distr[phase];
                                    }
                                }
                            } else {
                                OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                            }


                            break;
                        }
                        std::vector<double> g = {1,1,0.01};
                        if (well_controls_iget_type(wc, current) == RESERVOIR_RATE) {
                            const double* distr = well_controls_iget_distr(wc, current);
                            for (int phase = 0; phase < np; ++phase) {
                                g[phase] = distr[phase];
                            }
                        }
                        switch (well_controls_iget_type(wc, current)) {
                        case BHP:
                        {
                            const WellType& well_type = wells().type[w];
                            xw.wellSolutions()[w] = 0.0;
                            if (well_type == INJECTOR) {
                                for (int p = 0; p < np; ++p)  {
                                    xw.wellSolutions()[w] += xw.wellRates()[np*w + p] * wells().comp_frac[np*w + p];
                                }
                            } else {
                                for (int p = 0; p < np; ++p)  {
                                    xw.wellSolutions()[w] += g[p] * xw.wellRates()[np*w + p];
                                }
                            }
                        }
                            break;


                        case RESERVOIR_RATE: // Intentional fall-through
                        case SURFACE_RATE:
                        {
                            xw.wellSolutions()[w] = xw.bhp()[w];
                        }
                            break;
                        }

                        double tot_well_rate = 0.0;
                        for (int p = 0; p < np; ++p)  {
                            tot_well_rate += g[p] * xw.wellRates()[np*w + p];
                        }
                        if(std::abs(tot_well_rate) > 0) {
                            xw.wellSolutions()[nw + w] = g[Water] * xw.wellRates()[np*w + Water] / tot_well_rate; //wells->comp_frac[np*w + Water]; // Water;
                            xw.wellSolutions()[2*nw + w] = g[Gas] * xw.wellRates()[np*w + Gas] / tot_well_rate ; //wells->comp_frac[np*w + Gas]; //Gas
                        } else {
                            //xw.wellSolutions()[nw + w] =  wells().comp_frac[np*w + Water];
                            //xw.wellSolutions()[2 * nw + w] =  wells().comp_frac[np*w + Gas];
                        }
                    }
                }
            }






            template <typename Simulator, class WellState>
            void computeWellConnectionPressures(const Simulator& ebosSimulator,
                                                const WellState& xw)
            {
                if( ! localWellsActive() ) return ;
                // 1. Compute properties required by computeConnectionPressureDelta().
                //    Note that some of the complexity of this part is due to the function
                //    taking std::vector<double> arguments, and not Eigen objects.
                std::vector<double> b_perf;
                std::vector<double> rsmax_perf;
                std::vector<double> rvmax_perf;
                std::vector<double> surf_dens_perf;
                computePropertiesForWellConnectionPressures(ebosSimulator, xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

                const Vector& pdepth = perf_cell_depth_;
                const int nperf = wells().well_connpos[wells().number_of_wells];
                const std::vector<double> depth_perf(pdepth.data(), pdepth.data() + nperf);

                computeWellConnectionDensitesPressures(xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf, depth_perf, gravity_);

            }

            // state0 is non-constant, while it will not be used outside of the function
            template <class SolutionState, class WellState>
            void
            computeWellPotentials(const std::vector<ADB>& mob_perfcells,
                                  const std::vector<ADB>& b_perfcells,
                                  SolutionState& state0,
                                  WellState& well_state)
            {
                const int nw = wells().number_of_wells;
                const int np = wells().number_of_phases;
                const Opm::PhaseUsage& pu = fluid_->phaseUsage();

                Vector bhps = Vector::Zero(nw);
                for (int w = 0; w < nw; ++w) {
                    const WellControls* ctrl = wells().ctrls[w];
                    const int nwc = well_controls_get_num(ctrl);
                    //Loop over all controls until we find a BHP control
                    //or a THP control that specifies what we need.
                    //Pick the value that gives the most restrictive flow
                    for (int ctrl_index=0; ctrl_index < nwc; ++ctrl_index) {

                        if (well_controls_iget_type(ctrl, ctrl_index) == BHP) {
                            bhps[w] = well_controls_iget_target(ctrl, ctrl_index);
                        }

                        if (well_controls_iget_type(ctrl, ctrl_index) == THP) {
                            double aqua = 0.0;
                            double liquid = 0.0;
                            double vapour = 0.0;

                            if ((*active_)[ Water ]) {
                                aqua = well_state.wellRates()[w*np + pu.phase_pos[ Water ] ];
                            }
                            if ((*active_)[ Oil ]) {
                                liquid = well_state.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                            }
                            if ((*active_)[ Gas ]) {
                                vapour = well_state.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                            }

                            const int vfp        = well_controls_iget_vfp(ctrl, ctrl_index);
                            const double& thp    = well_controls_iget_target(ctrl, ctrl_index);
                            const double& alq    = well_controls_iget_alq(ctrl, ctrl_index);

                            //Set *BHP* target by calculating bhp from THP
                            const WellType& well_type = wells().type[w];

                            if (well_type == INJECTOR) {
                                double dp = wellhelpers::computeHydrostaticCorrection(
                                            wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                            wellPerforationDensities(), gravity_);
                                const double bhp = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                                // apply the strictest of the bhp controlls i.e. smallest bhp for injectors
                                if ( bhp < bhps[w]) {
                                    bhps[w] = bhp;
                                }
                            }
                            else if (well_type == PRODUCER) {
                                double dp = wellhelpers::computeHydrostaticCorrection(
                                            wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                            wellPerforationDensities(), gravity_);

                                const double bhp = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                                // apply the strictest of the bhp controlls i.e. largest bhp for producers
                                if ( bhp > bhps[w]) {
                                    bhps[w] = bhp;
                                }
                            }
                            else {
                                OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                            }
                        }
                    }

                }

                // use bhp limit from control
                state0.bhp = ADB::constant(bhps);

                // compute well potentials
                Vector aliveWells;
                std::vector<ADB> well_potentials;
                computeWellFlux(state0, mob_perfcells,  b_perfcells, aliveWells, well_potentials);

                // store well potentials in the well state
                // transform to a single vector instead of separate vectors pr phase
                const int nperf = wells().well_connpos[nw];
                Vector cq = superset(well_potentials[0].value(), Span(nperf, np, 0), nperf*np);
                for (int phase = 1; phase < np; ++phase) {
                    cq += superset(well_potentials[phase].value(), Span(nperf, np, phase), nperf*np);
                }
                well_state.wellPotentials().assign(cq.data(), cq.data() + nperf*np);
            }


            /// If set, computeWellFlux() will additionally store the
            /// total reservoir volume perforation fluxes.
            void setStoreWellPerforationFluxesFlag(const bool store_fluxes)
            {
                store_well_perforation_fluxes_ = store_fluxes;
            }

            /// Retrieves the stored fluxes. It is an error to call this
            /// unless setStoreWellPerforationFluxesFlag(true) has been
            /// called.
            const Vector& getStoredWellPerforationFluxes() const
            {
                assert(store_well_perforation_fluxes_);
                return well_perforation_fluxes_;
            }

            int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
            {
                const int flowToEbos[ 3 ] = { FluidSystem::waterPhaseIdx, FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx };
                return flowToEbos[ phaseIdx ];
            }

            /// upate the dynamic lists related to economic limits
            template<class WellState>
            void
            updateListEconLimited(ScheduleConstPtr schedule,
                                  const int current_step,
                                  const Wells* wells_struct,
                                  const WellState& well_state,
                                  DynamicListEconLimited& list_econ_limited) const
            {
                const int nw = wells_struct->number_of_wells;

                for (int w = 0; w < nw; ++w) {
                    // flag to check if the mim oil/gas rate limit is violated
                    bool rate_limit_violated = false;
                    const std::string& well_name = wells_struct->name[w];
                    const Well* well_ecl = schedule->getWell(well_name);
                    const WellEconProductionLimits& econ_production_limits = well_ecl->getEconProductionLimits(current_step);

                    // economic limits only apply for production wells.
                    if (wells_struct->type[w] != PRODUCER) {
                        continue;
                    }

                    // if no limit is effective here, then continue to the next well
                    if ( !econ_production_limits.onAnyEffectiveLimit() ) {
                        continue;
                    }
                    // for the moment, we only handle rate limits, not handling potential limits
                    // the potential limits should not be difficult to add
                    const WellEcon::QuantityLimitEnum& quantity_limit = econ_production_limits.quantityLimit();
                    if (quantity_limit == WellEcon::POTN) {
                        const std::string msg = std::string("POTN limit for well ") + well_name + std::string(" is not supported for the moment. \n")
                                              + std::string("All the limits will be evaluated based on RATE. ");
                        OpmLog::warning("NOT_SUPPORTING_POTN", msg);
                    }

                    const WellMapType& well_map = well_state.wellMap();
                    const typename WellMapType::const_iterator i_well = well_map.find(well_name);
                    assert(i_well != well_map.end()); // should always be found?
                    const WellMapEntryType& map_entry = i_well->second;
                    const int well_number = map_entry[0];

                    if (econ_production_limits.onAnyRateLimit()) {
                        rate_limit_violated = checkRateEconLimits(econ_production_limits, well_state, well_number);
                    }

                    if (rate_limit_violated) {
                        if (econ_production_limits.endRun()) {
                            const std::string warning_message = std::string("ending run after well closed due to economic limits is not supported yet \n")
                                                              + std::string("the program will keep running after ") + well_name + std::string(" is closed");
                            OpmLog::warning("NOT_SUPPORTING_ENDRUN", warning_message);
                        }

                        if (econ_production_limits.validFollowonWell()) {
                            OpmLog::warning("NOT_SUPPORTING_FOLLOWONWELL", "opening following on well after well closed is not supported yet");
                        }

                        if (well_ecl->getAutomaticShutIn()) {
                            list_econ_limited.addShutWell(well_name);
                            const std::string msg = std::string("well ") + well_name + std::string(" will be shut in due to economic limit");
                            OpmLog::info(msg);
                        } else {
                            list_econ_limited.addStoppedWell(well_name);
                            const std::string msg = std::string("well ") + well_name + std::string(" will be stopped due to economic limit");
                            OpmLog::info(msg);
                        }
                        // the well is closed, not need to check other limits
                        continue;
                    }

                    // checking for ratio related limits, mostly all kinds of ratio.
                    bool ratio_limits_violated = false;
                    RatioCheckTuple ratio_check_return;

                    if (econ_production_limits.onAnyRatioLimit()) {
                        ratio_check_return = checkRatioEconLimits(econ_production_limits, well_state, map_entry);
                        ratio_limits_violated = std::get<0>(ratio_check_return);
                    }

                    if (ratio_limits_violated) {
                        const bool last_connection = std::get<1>(ratio_check_return);
                        const int worst_offending_connection = std::get<2>(ratio_check_return);

                        const int perf_start = map_entry[1];

                        assert((worst_offending_connection >= 0) && (worst_offending_connection <  map_entry[2]));

                        const int cell_worst_offending_connection = wells_struct->well_cells[perf_start + worst_offending_connection];
                        list_econ_limited.addClosedConnectionsForWell(well_name, cell_worst_offending_connection);
                        const std::string msg = std::string("Connection ") + std::to_string(worst_offending_connection) + std::string(" for well ")
                                              + well_name + std::string(" will be closed due to economic limit");
                        OpmLog::info(msg);

                        if (last_connection) {
                            list_econ_limited.addShutWell(well_name);
                            const std::string msg2 = well_name + std::string(" will be shut due to the last connection closed");
                            OpmLog::info(msg2);
                        }
                    }

                }
            }

            template <class WellState>
            void computeWellConnectionDensitesPressures(const WellState& xw,
                                                        const std::vector<double>& b_perf,
                                                        const std::vector<double>& rsmax_perf,
                                                        const std::vector<double>& rvmax_perf,
                                                        const std::vector<double>& surf_dens_perf,
                                                        const std::vector<double>& depth_perf,
                                                        const double grav) {
                // Compute densities
                std::vector<double> cd =
                        WellDensitySegmented::computeConnectionDensities(
                                wells(), xw, fluid_->phaseUsage(),
                                b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

                const int nperf = wells().well_connpos[wells().number_of_wells];

                // Compute pressure deltas
                std::vector<double> cdp =
                        WellDensitySegmented::computeConnectionPressureDelta(
                                wells(), depth_perf, cd, grav);

                // Store the results
                well_perforation_densities_ = Eigen::Map<const Vector>(cd.data(), nperf);
                well_perforation_pressure_diffs_ = Eigen::Map<const Vector>(cdp.data(), nperf);

            }

        protected:
            bool wells_active_;
            const Wells*   wells_;
            const WellOps  wops_;
            ModelParameters param_;
            bool terminal_output_;

            const BlackoilPropsAdInterface* fluid_;
            const std::vector<bool>*  active_;
            const std::vector<PhasePresence>*  phase_condition_;
            const VFPProperties* vfp_properties_;
            double gravity_;
            // the depth of the all the cell centers
            // for standard Wells, it the same with the perforation depth
            Vector perf_cell_depth_;

            Vector well_perforation_densities_;
            Vector well_perforation_pressure_diffs_;

            bool store_well_perforation_fluxes_;
            Vector well_perforation_fluxes_;

            std::vector<EvalWell> wellVariables_;
            std::vector<double> F0_;
            const std::vector<double>& pv_;
            Mat duneA_;
            Mat duneB_;
            Mat duneC_;
            Mat invDuneD_;
            BVector resWell_;
            BVector rhs_;

            // protected methods
            EvalWell getBhp(const int wellIdx) const {
                const WellControls* wc = wells().ctrls[wellIdx];
                if (well_controls_get_current_type(wc) == BHP) {
                    EvalWell bhp = 0.0;
                    const double target_rate = well_controls_get_current_target(wc);
                    bhp.value = target_rate;
                    return bhp;
                }
                return wellVariables_[wellIdx];
            }

            EvalWell getQs(const int wellIdx, const int phaseIdx) const {
                EvalWell qs = 0.0;
                const WellControls* wc = wells().ctrls[wellIdx];
                const int np = fluid_->numPhases();
                const double target_rate = well_controls_get_current_target(wc);

                if (wells().type[wellIdx] == INJECTOR) {
                    const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];
                    if (comp_frac == 0.0)
                        return qs;

                    if (well_controls_get_current_type(wc) == BHP) {
                        return wellVariables_[wellIdx];
                    }
                    qs.value = target_rate;
                    return qs;
                }

                // Producers
                if (well_controls_get_current_type(wc) == BHP) {
                    return wellVariables_[wellIdx] * wellVolumeFractionScaled(wellIdx,phaseIdx);
                }
                if (well_controls_get_current_type(wc) == SURFACE_RATE) {
                    const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];

                    if (comp_frac == 1.0) {
                        qs.value = target_rate;
                        return qs;
                    }
                    int currentControlIdx = 0;
                    for (int i = 0; i < np; ++i) {
                        currentControlIdx += wells().comp_frac[np*wellIdx + i] * i;
                    }

                    if (wellVolumeFractionScaled(wellIdx,currentControlIdx) == 0) {
                        return qs;
                    }
                    return (target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx) / wellVolumeFractionScaled(wellIdx,currentControlIdx));
                }
                // ReservoirRate
                return target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx);
            }

            EvalWell wellVolumeFraction(const int wellIdx, const int phaseIdx) const {
                assert(fluid_.numPhases() == 3);
                const int nw = wells().number_of_wells;
                if (phaseIdx == Water) {
                   return wellVariables_[nw + wellIdx];
                }

                if (phaseIdx == Gas) {
                   return wellVariables_[2*nw + wellIdx];
                }

                // Oil
                return 1.0 - wellVariables_[nw + wellIdx] - wellVariables_[2 * nw + wellIdx];
            }

            EvalWell wellVolumeFractionScaled(const int wellIdx, const int phaseIdx) const {
                const WellControls* wc = wells().ctrls[wellIdx];
                if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
                    const double* distr = well_controls_get_current_distr(wc);
                    return wellVolumeFraction(wellIdx, phaseIdx) / distr[phaseIdx];
                }
                std::vector<double> g = {1,1,0.01};
                return (wellVolumeFraction(wellIdx, phaseIdx) / g[phaseIdx]);
            }



            template <class WellState>
            bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                     const WellState& well_state,
                                     const int well_number) const
            {
                const Opm::PhaseUsage& pu = fluid_->phaseUsage();
                const int np = well_state.numPhases();

                if (econ_production_limits.onMinOilRate()) {
                    assert((*active_)[Oil]);
                    const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
                    const double min_oil_rate = econ_production_limits.minOilRate();
                    if (std::abs(oil_rate) < min_oil_rate) {
                        return true;
                    }
                }

                if (econ_production_limits.onMinGasRate() ) {
                    assert((*active_)[Gas]);
                    const double gas_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Gas ] ];
                    const double min_gas_rate = econ_production_limits.minGasRate();
                    if (std::abs(gas_rate) < min_gas_rate) {
                        return true;
                    }
                }

                if (econ_production_limits.onMinLiquidRate() ) {
                    assert((*active_)[Oil]);
                    assert((*active_)[Water]);
                    const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
                    const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
                    const double liquid_rate = oil_rate + water_rate;
                    const double min_liquid_rate = econ_production_limits.minLiquidRate();
                    if (std::abs(liquid_rate) < min_liquid_rate) {
                        return true;
                    }
                }

                if (econ_production_limits.onMinReservoirFluidRate()) {
                    OpmLog::warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
                }

                return false;
            }


            using WellMapType = typename WellState::WellMapType;
            using WellMapEntryType = typename WellState::mapentry_t;

            // a tuple type for ratio limit check.
            // first value indicates whether ratio limit is violated, when the ratio limit is not violated, the following three
            // values should not be used.
            // second value indicates whehter there is only one connection left.
            // third value indicates the indx of the worst-offending connection.
            // the last value indicates the extent of the violation for the worst-offending connection, which is defined by
            // the ratio of the actual value to the value of the violated limit.
            using RatioCheckTuple = std::tuple<bool, bool, int, double>;

            enum ConnectionIndex {
                INVALIDCONNECTION = -10000
            };


            template <class WellState>
            RatioCheckTuple checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                                 const WellState& well_state,
                                                 const WellMapEntryType& map_entry) const
            {
                // TODO: not sure how to define the worst-offending connection when more than one
                //       ratio related limit is violated.
                //       The defintion used here is that we define the violation extent based on the
                //       ratio between the value and the corresponding limit.
                //       For each violated limit, we decide the worst-offending connection separately.
                //       Among the worst-offending connections, we use the one has the biggest violation
                //       extent.

                bool any_limit_violated = false;
                bool last_connection = false;
                int worst_offending_connection = INVALIDCONNECTION;
                double violation_extent = -1.0;

                if (econ_production_limits.onMaxWaterCut()) {
                    const RatioCheckTuple water_cut_return = checkMaxWaterCutLimit(econ_production_limits, well_state, map_entry);
                    bool water_cut_violated = std::get<0>(water_cut_return);
                    if (water_cut_violated) {
                        any_limit_violated = true;
                        const double violation_extent_water_cut = std::get<3>(water_cut_return);
                        if (violation_extent_water_cut > violation_extent) {
                            violation_extent = violation_extent_water_cut;
                            worst_offending_connection = std::get<2>(water_cut_return);
                            last_connection = std::get<1>(water_cut_return);
                        }
                    }
                }

                if (econ_production_limits.onMaxGasOilRatio()) {
                    OpmLog::warning("NOT_SUPPORTING_MAX_GOR", "the support for max Gas-Oil ratio is not implemented yet!");
                }

                if (econ_production_limits.onMaxWaterGasRatio()) {
                    OpmLog::warning("NOT_SUPPORTING_MAX_WGR", "the support for max Water-Gas ratio is not implemented yet!");
                }

                if (econ_production_limits.onMaxGasLiquidRatio()) {
                    OpmLog::warning("NOT_SUPPORTING_MAX_GLR", "the support for max Gas-Liquid ratio is not implemented yet!");
                }

                if (any_limit_violated) {
                    assert(worst_offending_connection >=0);
                    assert(violation_extent > 1.);
                }

                return std::make_tuple(any_limit_violated, last_connection, worst_offending_connection, violation_extent);
            }

            template <class WellState>
            RatioCheckTuple checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                                                  const WellState& well_state,
                                                  const WellMapEntryType& map_entry) const
            {
                bool water_cut_limit_violated = false;
                int worst_offending_connection = INVALIDCONNECTION;
                bool last_connection = false;
                double violation_extent = -1.0;

                const int np = well_state.numPhases();
                const Opm::PhaseUsage& pu = fluid_->phaseUsage();
                const int well_number = map_entry[0];

                assert((*active_)[Oil]);
                assert((*active_)[Water]);

                const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
                const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
                const double liquid_rate = oil_rate + water_rate;
                double water_cut;
                if (std::abs(liquid_rate) != 0.) {
                    water_cut = water_rate / liquid_rate;
                } else {
                    water_cut = 0.0;
                }

                const double max_water_cut_limit = econ_production_limits.maxWaterCut();
                if (water_cut > max_water_cut_limit) {
                    water_cut_limit_violated = true;
                }

                if (water_cut_limit_violated) {
                    // need to handle the worst_offending_connection
                    const int perf_start = map_entry[1];
                    const int perf_number = map_entry[2];

                    std::vector<double> water_cut_perf(perf_number);
                    for (int perf = 0; perf < perf_number; ++perf) {
                        const int i_perf = perf_start + perf;
                        const double oil_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Oil ] ];
                        const double water_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Water ] ];
                        const double liquid_perf_rate = oil_perf_rate + water_perf_rate;
                        if (std::abs(liquid_perf_rate) != 0.) {
                            water_cut_perf[perf] = water_perf_rate / liquid_perf_rate;
                        } else {
                            water_cut_perf[perf] = 0.;
                        }
                    }

                    last_connection = (perf_number == 1);
                    if (last_connection) {
                        worst_offending_connection = 0;
                        violation_extent = water_cut_perf[0] / max_water_cut_limit;
                        return std::make_tuple(water_cut_limit_violated, last_connection, worst_offending_connection, violation_extent);
                    }

                    double max_water_cut_perf = 0.;
                    for (int perf = 0; perf < perf_number; ++perf) {
                        if (water_cut_perf[perf] > max_water_cut_perf) {
                            worst_offending_connection = perf;
                            max_water_cut_perf = water_cut_perf[perf];
                        }
                    }

                    assert(max_water_cut_perf != 0.);
                    assert((worst_offending_connection >= 0) && (worst_offending_connection < perf_number));

                    violation_extent = max_water_cut_perf / max_water_cut_limit;
                }

                return std::make_tuple(water_cut_limit_violated, last_connection, worst_offending_connection, violation_extent);
            }

        };


} // namespace Opm

#include "StandardWellsDense_impl.hpp"

#endif
