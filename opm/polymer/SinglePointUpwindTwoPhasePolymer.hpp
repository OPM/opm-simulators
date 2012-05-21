/*===========================================================================
//
// File: SinglePointUpwindTwoPhase.hpp
//
// Created: 2011-09-28 14:21:34+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

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

#ifndef OPM_SINGLEPOINTUPWINDTWOPHASEPOLYMER_HPP_HEADER
#define OPM_SINGLEPOINTUPWINDTWOPHASEPOLYMER_HPP_HEADER

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <vector>
#include <iostream>

namespace Opm {
    namespace polymer_reorder {
        class ModelParameterStorage {
        public:
            ModelParameterStorage(int nc, int totconn)
                : drho_(0.0), rockdensity_(0.0), mob_(0), 
                  dmobds_(0), dmobwatdc_(0), mc_(0),
                  dmcdc_(0), porevol_(0), porosity_(0), dg_(0), sw_(0), c_(0), cmax_(0),
                  ds_(0), dsc_(0), dcads_(0), dcadsdc_(0), pc_(0), dpc_(0), 
                  trans_(0), data_()
            {
                size_t alloc_sz;

                alloc_sz  = 2 * nc;        // mob_
                alloc_sz += 2 * nc;        // dmobds_
                alloc_sz +=  nc;           // dmobwatdc_
                alloc_sz +=  nc;           // mc_
                alloc_sz +=  nc;           // dmcdc_
                alloc_sz += 1 * nc;        // porevol_
                alloc_sz += 1 * nc;        // porosity_
                alloc_sz += 1 * totconn;   // dg_
                alloc_sz += 1 * nc;        // sw_
                alloc_sz += 1 * nc;        // c_
                alloc_sz += 1 * nc;        // cmax_
                alloc_sz += 1 * nc;        // ds_
                alloc_sz += 1 * nc;        // dsc_
                alloc_sz += 1 * nc;        // dcads_
                alloc_sz += 1 * nc;        // dcadsdc_
                alloc_sz += 1 * nc;        // pc_
                alloc_sz += 1 * nc;        // dpc_
                alloc_sz += 1 * totconn;   // trans_
                data_.resize(alloc_sz);

                mob_           = &data_[0];
                dmobds_        = mob_            + (2 * nc     );
                dmobwatdc_     = dmobds_         + (2 * nc     );
                mc_            = dmobwatdc_      + (1 * nc     );
                dmcdc_         = mc_             + (1 * nc     );
                porevol_       = dmcdc_          + (1 * nc     );
                porosity_      = porevol_        + (1 * nc     );
                dg_            = porosity_       + (1 * nc     );
                sw_            = dg_             + (1 * totconn);
                c_             = sw_             + (1 * nc     );
                cmax_          = c_              + (1 * nc     );
                ds_            = cmax_           + (1 * nc     );
                dsc_           = ds_             + (1 * nc     );
                dcads_         = dsc_            + (1 * nc     );
                dcadsdc_       = dcads_          + (1 * nc     );
                pc_            = dcadsdc_        + (1 * nc     );
                dpc_           = pc_             + (1 * nc     );
                trans_         = dpc_            + (1 * nc     );
            }

            double&       drho   ()            { return drho_            ; }
            double        drho   ()      const { return drho_            ; }

            double&       rockdensity()            { return rockdensity_      ; }
            double        rockdensity()    const { return rockdensity_  ; }

            double*       mob    (int cell)       { return mob_  + (2*cell + 0); }
            const double* mob    (int cell) const { return mob_  + (2*cell + 0); }

            double*       dmobds    (int cell)       { return dmobds_  + (2*cell + 0); }
            const double* dmobds    (int cell) const { return dmobds_  + (2*cell + 0); }

            double&       dmobwatdc    (int cell)       { return dmobwatdc_[cell]; }
            double        dmobwatdc    (int cell) const { return dmobwatdc_[cell]; }

            double&       mc    (int cell)       { return mc_[cell]; }
            double        mc    (int cell) const { return mc_[cell]; }

            double&       dmcdc    (int cell)       { return dmcdc_[cell]; }
            double        dmcdc    (int cell) const { return dmcdc_[cell]; }

            double*       porevol()            { return porevol_      ; }
            double        porevol(int cell) const { return porevol_[cell]   ; }

            double*       porosity()            { return porosity_      ; }
            double        porosity(int cell) const { return porosity_[cell]   ; }


            double&       dg(int i)            { return dg_[i]        ; }
            double        dg(int i)      const { return dg_[i]        ; }

            double&       sw(int cell)            { return sw_[cell]        ; }
            double        sw(int cell)      const { return sw_[cell]        ; }

            double&       c(int cell)            { return c_[cell]        ; }
            double        c(int cell)      const { return c_[cell]        ; }

            double&       cmax(int cell)            { return cmax_[cell]        ; }
            double        cmax(int cell)      const { return cmax_[cell]        ; }

            double&       ds(int cell)            { return ds_[cell]        ; }
            double        ds(int cell)      const { return ds_[cell]        ; }

            double&       dsc(int cell)            { return dsc_[cell]        ; }
            double        dsc(int cell)      const { return dsc_[cell]        ; }

            double&       dcads(int cell)            { return dcads_[cell]        ; }
            double        dcads(int cell)      const { return dcads_[cell]        ; }

            double&       dcadsdc(int cell)            { return dcadsdc_[cell]        ; }
            double        dcadsdc(int cell)      const { return dcadsdc_[cell]        ; }

            double&       pc(int cell)            { return pc_[cell]        ; }
            double        pc(int cell)      const { return pc_[cell]        ; }

            double&       dpc(int cell)           { return dpc_[cell]       ; }
            double        dpc(int cell)     const { return dpc_[cell]       ; }

            double&       trans(int f)         { return trans_[f]     ; }
            double        trans(int f)   const { return trans_[f]     ; }

        private:
            double  drho_        ;
            double  rockdensity_  ;
            double *mob_         ;
            double *dmobds_      ;
            double *dmobwatdc_   ;
            double *mc_          ;
            double *dmcdc_       ;
            double *porevol_     ;
            double *porosity_    ;
            double *dg_          ;
            double *sw_          ;
            double *c_           ;
            double *cmax_           ;
            double *ds_          ;
            double *dsc_         ;
            double *dcads_       ; // difference of cads to compute residual
            double *dcadsdc_     ; // derivative of cads
            double *pc_          ;
            double *dpc_         ;
            double *trans_       ;

            std::vector<double> data_;
        };
    }


    template <class TwophaseFluidPolymer>
    class SinglePointUpwindTwoPhasePolymer {
    public:
        template <class Grid>
        SinglePointUpwindTwoPhasePolymer(const TwophaseFluidPolymer&  fluid    ,
					 const Grid&                  g        ,
					 const std::vector<double>&   porevol  ,
					 const double*                grav  = 0,
					 const bool                   guess_previous = true)
            : fluid_   (fluid)                              ,
              gravity_ (grav)                               ,
              f2hf_    (2 * g.number_of_faces, -1)          ,
              store_   (g.number_of_cells,
			g.cell_facepos[ g.number_of_cells ]),
	      init_step_use_previous_sol_(guess_previous)   ,
	      sat_tol_  (1e-5)
        {

            if (gravity_) {
                store_.drho() = fluid_.density(0) - fluid_.density(1);
            }

            for (int c = 0, i = 0; c < g.number_of_cells; ++c) {
                for (; i < g.cell_facepos[c + 1]; ++i) {
                    const int f = g.cell_faces[i];
                    const int p = 1 - (g.face_cells[2*f + 0] == c);
                    f2hf_[2*f + p] = i;
                }
            }

            std::copy(porevol.begin(), porevol.end(), store_.porevol());
            const double* poro = fluid.porosity();
            std::copy(poro, poro + g.number_of_cells, store_.porosity());
            store_.rockdensity() = fluid.rockdensity();
        }

        void makefhfQPeriodic(  const std::vector<int>& p_faces,const std::vector<int>& hf_faces,
                                const std::vector<int>& nb_faces)
        {
            std::vector<int> nbhf(hf_faces.size());
            for(unsigned int i=0; i<p_faces.size(); ++i){
                int nbf = nb_faces[i];
                if(f2hf_[2*nbf] == -1){
                    nbhf[i] = f2hf_[2*nbf+1];
                }else{
                    assert(f2hf_[2*nbf+1]==-1);
                    nbhf[i] = f2hf_[2*nbf];
                }
            }
            for(unsigned int i=0; i<p_faces.size(); ++i){

                int f = p_faces[i];
                int hf = hf_faces[i];
                bool changed=false;

                if(f2hf_[2*f] == hf){
                    assert(f2hf_[2*f+1]==-1);
                }else{
                    assert(f2hf_[2*f]==-1);
                    f2hf_[2*f]=nbhf[i];
                    changed=true;
                }
                if(!changed){
                    if(f2hf_[2*f+1]== hf){
                        assert(f2hf_[2*f]==-1);
                    }else{
                        assert(f2hf_[2*f+1]==-1);
                        f2hf_[2*f+1]=nbhf[i];
                        changed=true;
                    }
                }
                assert(changed);
            }
        }

        // -----------------------------------------------------------------
        // System assembly innards
        // -----------------------------------------------------------------

        enum { DofPerCell = 1 };

        void
        initResidual(const int c, double* Fs, double* Fc) const {
            (void) c;       // Suppress 'unused' warning
            *Fs = 0.0;
	    *Fc = 0.0;
        }

        template <class ReservoirState,
                  class Grid          >
        void
        fluxConnection(const ReservoirState& state  ,
                       const Grid&           g      ,
                       const double          dt     ,
                       const int             cell   ,
                       const int             f      ,
                       double* F                    , // F[0] = s-residual, F[1] = c-residual
                       double* dFd1                 , //Jacobi matrix for residual with respect to variables in cell
                       double* dFd2                   //Jacobi matrix for residual with respect to variables in OTHER cell
                                                      //dFd1[0]= d(F[0])/d(s1), dFd1[1]= d(F[0])/d(c1), dFd1[2]= d(F[1])/d(s1), dFd1[3]= d(F[1])/d(c1),
                                                      //dFd2[0]= d(F[0])/d(s2), dFd2[1]= d(F[0])/d(c2), dFd2[2]= d(F[1])/d(s2), dFd2[3]= d(F[1])/d(c2).

		       ) const {

            const int *n = g.face_cells + (2 * f);
            double dflux = state.faceflux()[f];
            double gflux = gravityFlux(f);
            double pcflux, dpcflux[2];
            capFlux(f, n, pcflux, dpcflux);
            gflux += pcflux;

            int    pix[2];
            double m[2], dmds[2], dmobwatdc;
            double mc, dmcdc;
            upwindMobility(dflux, gflux, n, pix, m, dmds, dmobwatdc, mc, dmcdc);

            assert ((m[0] >= 0.0) && (m[1] >= 0.0));

            double mt = m[0] + m[1];
            assert (mt >= 0.0);

            double sgn  = 2.0*(n[0] == cell) - 1.0;
            dflux      *= sgn;
            gflux      *= sgn;


            double       f1 = m[0] / mt;
            const double v1 = dflux + m[1]*gflux;

            // Assemble residual contributions
            F[0] += dt * f1 * v1;
            F[1] += dt * mc * f1 * v1;

            // Assemble Jacobian (J1 <-> cell, J2 <-> other)
            double *dFsds[2];
            double *dFsdc[2];
            double *dFcds[2];
            double *dFcdc[2];
            if (n[0] == cell) {
                dFsds[0] = &dFd1[0]; dFsds[1] = &dFd2[0];
                dFsdc[0] = &dFd1[1]; dFsdc[1] = &dFd2[1];
                dFcds[0] = &dFd1[2]; dFcds[1] = &dFd2[2];
                dFcdc[0] = &dFd1[3]; dFcdc[1] = &dFd2[3];
                // sign is positive
                dFd1[0] += sgn*dt * f1            * dpcflux[0] * m[1];
                dFd2[0] += sgn*dt * f1            * dpcflux[1] * m[1];
                dFd1[2] += sgn*dt * f1 * mc       * dpcflux[0] * m[1];
                dFd2[2] += sgn*dt * f1 * mc       * dpcflux[1] * m[1];
		// We assume that the capillary pressure is independent of the polymer concentration.
                // Hence, no more contributions.
            } else {
                dFsds[0] = &dFd2[0]; dFsds[1] = &dFd1[0];
                dFsdc[0] = &dFd2[1]; dFsdc[1] = &dFd1[1];
                dFcds[0] = &dFd2[2]; dFcds[1] = &dFd1[2];
                dFcdc[0] = &dFd2[3]; dFcdc[1] = &dFd1[3];
                // sign is negative
                dFd1[0] += sgn*dt * f1            * dpcflux[1] * m[1];
                dFd2[0] += sgn*dt * f1            * dpcflux[0] * m[1];
                dFd1[2] += sgn*dt * f1 * mc       * dpcflux[1] * m[1];
                dFd2[2] += sgn*dt * f1 * mc       * dpcflux[0] * m[1];
		// We assume that the capillary pressure is independent of the polymer concentration.
                // Hence, no more contributions.
            }

            // dFs/dm_1 \cdot dm_1/ds
            *dFsds[ pix[0] ] += dt * (1 - f1) / mt * v1      * dmds[0];
            // dFc/dm_1 \cdot dm_1/ds
            *dFcds[ pix[0] ] += dt * (1 - f1) / mt * v1 * mc * dmds[0];

            // dFs/dm_2 \cdot dm_2/ds
            *dFsds[ pix[1] ] -= dt * f1       / mt * v1    *      dmds[1];
            *dFsds[ pix[1] ] += dt * f1            * gflux *      dmds[1];
            // dFc/dm_2 \cdot dm_2/ds
            *dFcds[ pix[1] ] -= dt * f1       / mt * v1    * mc * dmds[1];
            *dFcds[ pix[1] ] += dt * f1            * gflux * mc * dmds[1];

            // dFs/dm_1 \cdot dm_1/dc
            *dFsdc[ pix[0] ] += dt * (1 - f1) / mt * v1      * dmobwatdc;
            // dFc/dm_1 \cdot dm_1/dc
            *dFcdc[ pix[0] ] += dt * (1 - f1) / mt * v1 * mc * dmobwatdc;
            *dFcdc[ pix[0] ] += dt * f1 * v1 * dmcdc;                 // Polymer is only carried by water.
        }

        template <class Grid>
        void
        accumulation(const Grid& g,
                     const int   cell,
                     double*     F,    // Residual vector,
                     double*     dF    // Jacobian, same convention as for fluxConnection.
                     ) const {
            (void) g;

            const double pv   = store_.porevol(cell);
            const double dps  = fluid_.deadporespace();
            const double rhor = fluid_.rockdensity();
            const double poro = store_.porosity(cell);

            F[0]  += pv * store_.ds(cell);
            F[1]  += pv * (1 - dps) * store_.dsc(cell) + rhor*(1 - poro)/poro*pv*store_.dcads(cell);
            dF[0] += pv;
            dF[1] += 0.;
            dF[2] += pv * (1 - dps) * store_.c(cell);
            dF[3] += pv * (1 - dps) * store_.sw(cell) + rhor*(1 - poro)/poro*pv*store_.dcadsdc(cell);
        }

        template <class Grid       ,
                  class SourceTerms>
        void
        sourceTerms(const Grid&        g  ,
                    const SourceTerms* src,
                    const int          i  ,
                    const double       dt ,
                    double*            J  ,
                    double*            F  ) const {

            (void) g;

            double dflux = -src->flux[i]; // ->flux[] is rate of *inflow*

            if (dflux < 0) {
                // src -> cell, affects residual only.
                *F += dt * dflux * src->saturation[2*i + 0];
            } else {
                // cell -> src
                const int     cell  = src->cell[i];
                const double* m  = store_.mob (cell);
                const double* dm = store_.dmobds(cell);

                const double  mt = m[0] + m[1];

                assert (! ((m[0] < 0) || (m[1] < 0)));
                assert (mt > 0);

                const double f  = m[0] / mt;
                const double df = ((1 - f)*dm[0] - f*dm[1]) / mt;

                *F += dt * dflux *  f;
                *J += dt * dflux * df;
            }
        }
        template <class Grid>
        void
        initGravityTrans(const Grid&  g    ,
                         const std::vector<double> &  htrans) {

            assert (htrans.size() ==
                    static_cast<std::vector<double>::size_type>(g.cell_facepos[ g.number_of_cells ]));

            for (int f = 0; f < g.number_of_faces; ++f) {
                store_.trans(f) = 0.0;
            }

            for (int c = 0, i = 0; c < g.number_of_cells; ++c) {
                for (; i < g.cell_facepos[c + 1]; ++i) {
                    int f = g.cell_faces[i];

                    assert (htrans[i] > 0.0);

                    store_.trans(f) += 1.0 / htrans[i];
                }
            }

            for (int f = 0; f < g.number_of_faces; ++f) {
                store_.trans(f) = 1.0 / store_.trans(f);
            }

            if (gravity_) {
                this->computeStaticGravity(g);
            }
        }

        // -----------------------------------------------------------------
        // Newton control
        // -----------------------------------------------------------------

        template <class ReservoirState,
                  class Grid          ,
                  class JacobianSystem>
        void
        initStep(const ReservoirState& state,
                 const Grid&           g    ,
                 JacobianSystem&       sys  ) {

            (void) state;       // Suppress 'unused' warning.

            typename JacobianSystem::vector_type& x =
                sys.vector().writableSolution();

            assert (x.size() == (::std::size_t) (2*g.number_of_cells));

	    if (init_step_use_previous_sol_) {
		std::fill(x.begin(), x.end(), 0.0);
            } else {
                std::fill(x.begin(), x.end(), 0.0);
		const std::vector<double>& s = state.saturation();
		for (int cell = 0, ncell = g.number_of_cells; cell < ncell; ++cell) {
		    // Impose s=0.5 at next time level as an NR initial value.
		    x[2*cell + 0] = 0.5 - s[2*cell + 0];
		}
	    }
        }

        template <class ReservoirState,
                  class Grid          ,
                  class JacobianSystem>
        bool
        initIteration(const ReservoirState& state,
                      const Grid&           g    ,
                      JacobianSystem&       sys) {

            double s[2];
            double mob[2];
            double dmobds[4];
            double dmobwatdc;
            double c, cmax;
            double mc, dmcdc;
            double pc, dpc;

            const typename JacobianSystem::vector_type& x =
                sys.vector().solution();
            const ::std::vector<double>& sat = state.saturation();
            const ::std::vector<double>& cpoly = state.concentration();
            const ::std::vector<double>& cmaxpoly = state.maxconcentration();

            bool in_range = true;
            for (int cell = 0; cell < g.number_of_cells; ++cell) {
                // Store wat-sat, sat-change, cpoly, (sat * cpoly)-change for accumulation().
                store_.ds(cell) = x[2*cell + 0];
                s[0] = sat[cell*2 + 0] + x[2*cell + 0];
                c = cpoly[cell] + x[2*cell + 1];
                store_.sw(cell) = s[0];
                store_.c(cell) = c;
                cmax = std::max(c, cmaxpoly[cell]);
                store_.cmax(cell) = cmax;
                store_.dsc(cell) = s[0]*c - sat[cell*2 + 0]*cpoly[cell];
                double dcadsdc;
                double cads;
                fluid_.adsorption(cpoly[cell], cmaxpoly[cell], cads, dcadsdc);
                store_.dcads(cell) =  -cads;
                fluid_.adsorption(c, cmax, cads, dcadsdc);
                store_.dcads(cell) +=  cads;
                store_.dcadsdc(cell) = dcadsdc;
                double s_min = fluid_.s_min(cell);
                double s_max = fluid_.s_max(cell);

                if ( s[0] < (s_min - sat_tol_) || s[0] > (s_max + sat_tol_) ) {
                    // if (s[0] < s_min){
		    // 	std::cout << "Warning: s out of range, s-s_min = " << s_min-s[0] << std::endl;
                    // }
                    // if (s[0] > s_max){
		    // 	std::cout << "Warning: s out of range, s-s_max = " << s[0]-s_max << std::endl;
                    // }
                    in_range = false; //line search fails
                }
                s[0] = std::max(s_min, s[0]);
                s[0] = std::min(s_max, s[0]);
                s[1] = 1 - s[0];

                fluid_.mobility(cell, s, c, cmax, mob, dmobds, dmobwatdc);
                fluid_.computeMc(c, mc, dmcdc);
                fluid_.pc(cell, s, pc, dpc);

                store_.mob (cell)[0]   =  mob [0];
                store_.mob (cell)[1]   =  mob [1];
                store_.dmobds(cell)[0] =  dmobds[0*2 + 0];
                store_.dmobds(cell)[1] = -dmobds[1*2 + 1];
                store_.dmobwatdc(cell) =  dmobwatdc;
                store_.mc(cell)        = mc;
                store_.dmcdc(cell)     = dmcdc;
                store_.pc(cell)        = pc;
                store_.dpc(cell)       = dpc;
            }
	    if (!in_range) {
		std::cout << "Warning: initIteration() - s was clamped in some cells.\n";
	    }
            return in_range;
        }

        template <class ReservoirState,
                  class Grid          ,
                  class NewtonIterate >
        void
        finishIteration(const ReservoirState& state,
                        const Grid&           g    ,
                        NewtonIterate&        it   ) {
            // Nothing to do at end of iteration in this model.
            (void) state;  (void) g;  (void) it;
            typedef typename NewtonIterate::vector_type vector_t;
        }

        template <class Grid          ,
                  class SolutionVector,
                  class ReservoirState>
        void
        finishStep(const Grid&           g    ,
                   const SolutionVector& x    ,
                   ReservoirState&       state) {

            double *s    = &state.saturation()[0*2 + 0];
            double *c    = &state.concentration()[0*1 + 0];
            double *cmax = &state.maxconcentration()[0*1 + 0];

            for (int cell = 0; cell < g.number_of_cells; ++cell, s += 2, c += 1, cmax +=1) {
                s[0] += x[2*cell + 0];
                c[0] += x[2*cell + 1];
                cmax[0] = std::max(c[0], cmax[0]);
                double s_min = fluid_.s_min(cell);
                double s_max = fluid_.s_max(cell);
                assert(s[0] >= s_min - sat_tol_);
                assert(s[0] <= s_max + sat_tol_);
                s[0] = std::max(s_min, s[0]);
                s[0] = std::min(s_max, s[0]);
                s[1]  = 1.0 - s[0];
            }
        }

    private:
        void
        upwindMobility(const double dflux,
                       const double gflux,
                       const int*   n    ,
                       int*         pix  ,
                       double*      m    ,
                       double*      dmds ,
                       double&      dmobwatdc ,
                       double&      mc,
                       double&      dmcdc) const {
            bool equal_sign = ( (! (dflux < 0)) && (! (gflux < 0)) ) ||
                ( (! (dflux > 0)) && (! (gflux > 0)) );

            if (equal_sign) {

                if (! (dflux < 0) && ! (gflux < 0)) { pix[0] = 0; }
                else                                { pix[0] = 1; }

                m[0] = store_.mob(n[ pix[0] ]) [ 0 ];
                mc = store_.mc(n[ pix[0] ]);

                if (! (dflux - m[0]*gflux < 0))     { pix[1] = 0; }
                else                                { pix[1] = 1; }

                m[1] = store_.mob(n[ pix[1] ]) [ 1 ];

            } else {

                if (! (dflux < 0) && ! (gflux > 0)) { pix[1] = 0; }
                else                                { pix[1] = 1; }

                m[1] = store_.mob(n[ pix[1] ]) [ 1 ];

                if (dflux + m[1]*gflux > 0)         { pix[0] = 0; }
                else                                { pix[0] = 1; }

                m[0] = store_.mob(n[ pix[0] ]) [ 0 ];
                mc = store_.mc(n[ pix[0] ]);
            }

            dmds[0]   = store_.dmobds(n[ pix[0] ]) [ 0 ];
            dmds[1]   = store_.dmobds(n[ pix[1] ]) [ 1 ];
            dmobwatdc = store_.dmobwatdc(n[ pix[0] ]);
            dmcdc     = store_.dmcdc(n[ pix[0] ]);
        }

        template <class Grid>
        void
        computeStaticGravity(const Grid& g) {

            const int d = g.dimensions;

            for (int c = 0, i = 0; c < g.number_of_cells; ++c) {
                const double* cc = g.cell_centroids + (c * d);

                for (; i < g.cell_facepos[c + 1]; ++i) {
                    const int     f  = g.cell_faces[i];
                    const double* fc = g.face_centroids + (f * d);

                    double dg = 0.0;
                    for (int j = 0; j < d; ++j) {
                        dg += gravity_[j] * (fc[j] - cc[j]);
                    }

                    store_.dg(i) = store_.trans(f) * dg;
                }
            }
        }

        double
        gravityFlux(const int f) const {
            double gflux;

            if (gravity_) {
                int i1 = f2hf_[2*f + 0];
                int i2 = f2hf_[2*f + 1];

                assert ((i1 >= 0) && (i2 >= 0));

                gflux  = store_.dg(i1) - store_.dg(i2);
                gflux *= store_.drho();
            } else {
                gflux = 0.0;
            }

            return gflux;
        }
        void
        capFlux(const int f,const int* n,double& pcflux, double* dpcflux) const {
            //double capflux;
            int i1 = n[0];
            int i2 = n[1];
            assert ((i1 >= 0) && (i2 >= 0));
            //double sgn=-1.0;
            pcflux  = store_.trans(f)*(store_.pc(i2) - store_.pc(i1));
            dpcflux[0]  = -store_.trans(f)*store_.dpc(i1);
            dpcflux[1]  = store_.trans(f)*store_.dpc(i2);
        }

        TwophaseFluidPolymer          fluid_  ;
        const double*                 gravity_;
        std::vector<int>              f2hf_   ;
        polymer_reorder::ModelParameterStorage store_  ;
	bool init_step_use_previous_sol_;
	double sat_tol_;
    };
}
#endif  /* OPM_SINGLEPOINTUPWINDTWOPHASE_HPP_HEADER */
