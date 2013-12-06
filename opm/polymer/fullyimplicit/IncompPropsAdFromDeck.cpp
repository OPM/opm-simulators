/**/
#include <opm/polymer/fullyimplicit/IncompPropsAdFromDeck.hpp>
#include <opm/polymer/fullyimplicit/AutoDiffHelpers.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{
    IncompPropsAdFromDeck::IncompPropsAdFromDeck(const EclipseGridParser& deck,
                                                 const UnstructuredGrid& grid)
    {
        rock_.init(deck, grid);
        pvt_.init(deck);
        satprops_.init(deck, grid, 200);
        if (pvt_.numPhases() != satprops_.numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck() - "
                  "Inconsistent number of phases in pvt data (" << pvt_.numPhases()
                  << ") and saturation-dependent function data (" << satprops_.numPhases() << ").");
        }
    }
    IncompPropsAdFromDeck::~IncompPropsAdFromDeck()
    {
    }
    // rock interface
    int IncompPropsAdFromDeck::numDimensions() const
    {
        return rock_.numDimensions();
    }
    int IncompPropsAdFromDeck::numCells() const
    {
        return rock_.numCells();
    }

    const double* IncompPropsAdFromDeck::porosity() const
    {
        return rock_.porosity();
    }
    const double* IncompPropsAdFromDeck::permeability() const
    {
        return rock_.permeability();
    }

    // fluid interface
    int IncompPropsAdFromDeck::numPhases() const
    {
        return pvt_.numPhases();
    }


    const double* IncompPropsAdFromDeck::viscosity() const
    {
        return pvt_.viscosity();
    }


    const double* IncompPropsAdFromDeck::density() const
    {
        return pvt_.reservoirDensities();
    }


    const double* IncompPropsAdFromDeck::surfaceDensity() const
    {
        return pvt_.surfaceDensities();
    }


    typedef IncompPropsAdFromDeck::ADB ADB;
    typedef IncompPropsAdFromDeck::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;


    std::vector<V> 
    IncompPropsAdFromDeck::relperm(const V& sw,
                                   const V& so,
                                   const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        assert(sw.size() == n && so.size() == n);
        s_all.col(0) = sw;
        s_all.col(1) = so;
        Block kr(n, np);
        satprops_.relperm(n, s_all.data(), cells.data(), kr.data(), 0);
        std::vector<V> relperms;
        relperms.reserve(np);
        for (int phase = 0; phase < np; ++phase) {
            relperms.emplace_back(kr.col(phase));
        }
        return relperms;
    }



    std::vector<ADB> 
    IncompPropsAdFromDeck::relperm(const ADB& sw,
                                   const ADB& so,
                                   const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        assert(sw.size() == n && so.size() == n);
        s_all.col(0) = sw.value();
        s_all.col(1) = so.value();
        Block kr(n, np);
        Block dkr(n, np*np);
        satprops_.relperm(n, s_all.data(), cells.data(), kr.data(), dkr.data());
        const int num_blocks = so.numBlocks();
        std::vector<ADB> relperms;
        relperms.reserve(np);
        typedef const ADB* ADBPtr;
        ADBPtr s[2] = { &sw, &so };
        for (int phase1 = 0; phase1 < np; ++phase1) {
            const int phase1_pos = phase1;
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = ADB::M(n, s[phase1]->derivative()[block].cols());
            }
            for (int phase2 = 0; phase2 < np; ++phase2) {
                const int phase2_pos = phase2;
                // Assemble dkr1/ds2.
                const int column = phase1_pos + np*phase2_pos; // Recall: Fortran ordering from props_.relperm()
                ADB::M dkr1_ds2_diag = spdiag(dkr.col(column));
                for (int block = 0; block < num_blocks; ++block) {
                    jacs[block] += dkr1_ds2_diag * s[phase2]->derivative()[block];
                }
            }
                relperms.emplace_back(ADB::function(kr.col(phase1_pos), jacs));
        }
        return relperms;
    }
} //namespace Opm

