/*
  Copyright 2026 SINTEF Digital

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

#ifndef OPM_HYPRE_NUMERIC_RAP_HPP
#define OPM_HYPRE_NUMERIC_RAP_HPP

#include <_hypre_parcsr_ls.h>

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

namespace Opm::linalg::HypreInterface
{

/// A_H = P^T A P on fixed patterns: build() records the products once (P and its off-process
/// rows are fixed then); apply() does the multiply-adds and one exchange of remote coarse rows.
class NumericRap
{
public:
    /// Returns false (on every process) if a product falls outside A_H's pattern.
    bool build(hypre_ParCSRMatrix* A, hypre_ParCSRMatrix* P, hypre_ParCSRMatrix* AH)
    {
        *this = NumericRap{};
        comm_ = hypre_ParCSRMatrixComm(A);
        int nprocs = 1;
        MPI_Comm_size(comm_, &nprocs);

        hypre_CSRMatrix* Ad = hypre_ParCSRMatrixDiag(A);
        hypre_CSRMatrix* Ao = hypre_ParCSRMatrixOffd(A);
        hypre_CSRMatrix* Pd = hypre_ParCSRMatrixDiag(P);
        hypre_CSRMatrix* Po = hypre_ParCSRMatrixOffd(P);
        const HYPRE_Int n = hypre_CSRMatrixNumRows(Ad);
        nnzAd_ = hypre_CSRMatrixNumNonzeros(Ad);
        nnzPd_ = hypre_CSRMatrixNumNonzeros(Pd);
        nnzPo_ = hypre_CSRMatrixNumNonzeros(Po);
        const HYPRE_BigInt firstP = hypre_ParCSRMatrixFirstColDiag(P);
        const HYPRE_BigInt* cmapP = hypre_ParCSRMatrixColMapOffd(P);
        const HYPRE_BigInt ng = hypre_ParCSRMatrixGlobalNumCols(P);

        // rows of P for A's off-process columns, in the order of A's col_map_offd
        std::vector<HYPRE_Int> extI(1, 0);
        std::vector<HYPRE_BigInt> extJ;
        if (nprocs > 1 && hypre_CSRMatrixNumCols(Ao) > 0) {
            hypre_CSRMatrix* Pext = hypre_ParCSRMatrixExtractBExt(P, A, 1);
            const HYPRE_Int nr = hypre_CSRMatrixNumRows(Pext);
            const HYPRE_Int* pi = hypre_CSRMatrixI(Pext);
            extI.assign(pi, pi + nr + 1);
            extJ.assign(hypre_CSRMatrixBigJ(Pext), hypre_CSRMatrixBigJ(Pext) + pi[nr]);
            pext_.assign(hypre_CSRMatrixData(Pext), hypre_CSRMatrixData(Pext) + pi[nr]);
            hypre_CSRMatrixDestroy(Pext);
        }

        // Dense markers over the global coarse columns (coarse levels are small):
        // stamp[J] == row means pos[J] is valid for that row.
        std::vector<HYPRE_BigInt> stamp(ng, -1);
        std::vector<HYPRE_Int> pos(ng, 0);

        // A P with global coarse columns; value sources index [Ad | Ao] and [Pd | Po | Pext]
        const HYPRE_Int *Adi = hypre_CSRMatrixI(Ad), *Adj = hypre_CSRMatrixJ(Ad);
        const HYPRE_Int *Aoi = hypre_CSRMatrixI(Ao), *Aoj = hypre_CSRMatrixJ(Ao);
        const HYPRE_Int *Pdi = hypre_CSRMatrixI(Pd), *Pdj = hypre_CSRMatrixJ(Pd);
        const HYPRE_Int *Poi = hypre_CSRMatrixI(Po), *Poj = hypre_CSRMatrixJ(Po);
        std::vector<HYPRE_Int> APi(n + 1, 0);
        std::vector<HYPRE_BigInt> APJ;
        APJ.reserve(2 * static_cast<std::size_t>(nnzAd_));
        ap_a_.reserve(4 * static_cast<std::size_t>(nnzAd_));
        ap_p_.reserve(ap_a_.capacity());
        ap_dst_.reserve(ap_a_.capacity());
        auto addAP = [&](HYPRE_Int k, HYPRE_Int a, HYPRE_Int p, HYPRE_BigInt J) {
            if (stamp[J] != k) {
                stamp[J] = k;
                pos[J] = static_cast<HYPRE_Int>(APJ.size());
                APJ.push_back(J);
            }
            ap_a_.push_back(a);
            ap_p_.push_back(p);
            ap_dst_.push_back(pos[J]);
        };
        for (HYPRE_Int k = 0; k < n; ++k) {
            APi[k] = static_cast<HYPRE_Int>(APJ.size());
            for (HYPRE_Int a = Adi[k]; a < Adi[k + 1]; ++a) {
                const HYPRE_Int l = Adj[a];
                for (HYPRE_Int q = Pdi[l]; q < Pdi[l + 1]; ++q) {
                    addAP(k, a, q, firstP + Pdj[q]);
                }
                for (HYPRE_Int q = Poi[l]; q < Poi[l + 1]; ++q) {
                    addAP(k, a, nnzPd_ + q, cmapP[Poj[q]]);
                }
            }
            if (Aoi != nullptr) {
                for (HYPRE_Int a = Aoi[k]; a < Aoi[k + 1]; ++a) {
                    const HYPRE_Int e = Aoj[a];
                    for (HYPRE_Int q = extI[e]; q < extI[e + 1]; ++q) {
                        addAP(k, nnzAd_ + a, nnzPd_ + nnzPo_ + q, extJ[q]);
                    }
                }
            }
        }
        APi[n] = static_cast<HYPRE_Int>(APJ.size());
        ap_.assign(APJ.size(), 0.0);

        // all P^T (A P) products as flat records: target row I, P source p, AP slot b
        std::vector<HYPRE_BigInt> recI;
        std::vector<HYPRE_Int> recP, recB;
        for (HYPRE_Int k = 0; k < n; ++k) {
            const HYPRE_Int nb = APi[k + 1] - APi[k];
            for (HYPRE_Int q = Pdi[k]; q < Pdi[k + 1]; ++q) {
                for (HYPRE_Int b = APi[k]; b < APi[k] + nb; ++b) {
                    recI.push_back(firstP + Pdj[q]);
                    recP.push_back(q);
                    recB.push_back(b);
                }
            }
            if (Poi != nullptr) {
                for (HYPRE_Int q = Poi[k]; q < Poi[k + 1]; ++q) {
                    for (HYPRE_Int b = APi[k]; b < APi[k] + nb; ++b) {
                        recI.push_back(cmapP[Poj[q]]);
                        recP.push_back(nnzPd_ + q);
                        recB.push_back(b);
                    }
                }
            }
        }

        // slots of A_H's local rows: diag first, then offd (offset nnzHd)
        hypre_CSRMatrix* Hd = hypre_ParCSRMatrixDiag(AH);
        hypre_CSRMatrix* Ho = hypre_ParCSRMatrixOffd(AH);
        nnzHd_ = hypre_CSRMatrixNumNonzeros(Hd);
        const HYPRE_Int nH = hypre_CSRMatrixNumRows(Hd);
        const HYPRE_BigInt firstH = hypre_ParCSRMatrixFirstRowIndex(AH);
        const HYPRE_BigInt firstHc = hypre_ParCSRMatrixFirstColDiag(AH);
        const HYPRE_BigInt* cmapH = hypre_ParCSRMatrixColMapOffd(AH);
        const HYPRE_Int* Hdi = hypre_CSRMatrixI(Hd);
        const HYPRE_Int* Hdj = hypre_CSRMatrixJ(Hd);
        const HYPRE_Int* Hoi = hypre_CSRMatrixI(Ho);
        const HYPRE_Int* Hoj = hypre_CSRMatrixJ(Ho);
        auto stampRow = [&](HYPRE_Int i) {
            const HYPRE_BigInt tag = static_cast<HYPRE_BigInt>(n) + i; // above the fine-row stamps of A P
            for (HYPRE_Int h = Hdi[i]; h < Hdi[i + 1]; ++h) {
                stamp[firstHc + Hdj[h]] = tag;
                pos[firstHc + Hdj[h]] = h;
            }
            if (Hoi != nullptr) {
                for (HYPRE_Int h = Hoi[i]; h < Hoi[i + 1]; ++h) {
                    stamp[cmapH[Hoj[h]]] = tag;
                    pos[cmapH[Hoj[h]]] = nnzHd_ + h;
                }
            }
            return tag;
        };
        // resolve (I, J) pairs bucketed by local row with a counting sort
        bool ok = true;
        auto resolve = [&](const std::vector<HYPRE_BigInt>& I, auto&& colOf, std::vector<HYPRE_Int>& dst) {
            dst.assign(I.size(), -1);
            std::vector<HYPRE_Int> cnt(nH + 1, 0);
            for (const auto i : I) {
                if (i < firstH || i >= firstH + nH) {
                    ok = false;
                    return;
                }
                ++cnt[i - firstH + 1];
            }
            std::partial_sum(cnt.begin(), cnt.end(), cnt.begin());
            std::vector<std::size_t> order(I.size());
            std::vector<HYPRE_Int> fill(cnt.begin(), cnt.end() - 1);
            for (std::size_t t = 0; t < I.size(); ++t) {
                order[fill[I[t] - firstH]++] = t;
            }
            for (HYPRE_Int i = 0; i < nH; ++i) {
                const HYPRE_BigInt tag = stampRow(i);
                for (HYPRE_Int u = cnt[i]; u < cnt[i + 1]; ++u) {
                    const std::size_t t = order[u];
                    const HYPRE_BigInt J = colOf(t);
                    if (stamp[J] != tag) {
                        ok = false;
                        return;
                    }
                    dst[t] = pos[J];
                }
            }
        };

        // owner of each coarse row
        std::vector<HYPRE_BigInt> starts(nprocs + 1);
        MPI_Allgather(&firstP, 1, HYPRE_MPI_BIG_INT, starts.data(), 1, HYPRE_MPI_BIG_INT, comm_);
        starts[nprocs] = ng;
        auto ownerOf = [&](HYPRE_BigInt I) {
            return static_cast<int>(std::upper_bound(starts.begin(), starts.end(), I) - starts.begin()) - 1;
        };

        // split records into local and remote
        std::vector<HYPRE_BigInt> locI;
        std::vector<std::size_t> remote;
        for (std::size_t t = 0; t < recI.size(); ++t) {
            if (recI[t] >= firstH && recI[t] < firstH + nH) {
                locI.push_back(recI[t]);
                r_p_.push_back(recP[t]);
                r_ap_.push_back(recB[t]);
            } else {
                remote.push_back(t);
            }
        }
        {
            std::vector<HYPRE_Int> locB(r_ap_);
            resolve(locI, [&](std::size_t t) { return APJ[locB[t]]; }, r_dst_);
        }

        // remote: one send slot per distinct (owner, I, J), in owner order
        std::sort(remote.begin(), remote.end(), [&](std::size_t x, std::size_t y) {
            const int ox = ownerOf(recI[x]), oy = ownerOf(recI[y]);
            if (ox != oy) return ox < oy;
            if (recI[x] != recI[y]) return recI[x] < recI[y];
            return APJ[recB[x]] < APJ[recB[y]];
        });
        sendCounts_.assign(nprocs, 0);
        sendDispls_.assign(nprocs, 0);
        std::vector<HYPRE_BigInt> skeys;
        HYPRE_Int nslots = 0;
        for (std::size_t u = 0; u < remote.size(); ++u) {
            const std::size_t t = remote[u];
            const bool fresh = u == 0 || recI[remote[u - 1]] != recI[t] || APJ[recB[remote[u - 1]]] != APJ[recB[t]];
            if (fresh) {
                ++sendCounts_[ownerOf(recI[t])];
                skeys.push_back(recI[t]);
                skeys.push_back(APJ[recB[t]]);
                ++nslots;
            }
            s_p_.push_back(recP[t]);
            s_ap_.push_back(recB[t]);
            s_dst_.push_back(nslots - 1);
        }
        for (int o = 1; o < nprocs; ++o) {
            sendDispls_[o] = sendDispls_[o - 1] + sendCounts_[o - 1];
        }
        send_.assign(nslots, 0.0);

        recvCounts_.assign(nprocs, 0);
        recvDispls_.assign(nprocs, 0);
        if (nprocs > 1) {
            MPI_Alltoall(sendCounts_.data(), 1, MPI_INT, recvCounts_.data(), 1, MPI_INT, comm_);
            for (int o = 1; o < nprocs; ++o) {
                recvDispls_[o] = recvDispls_[o - 1] + recvCounts_[o - 1];
            }
            const int nrecv = recvDispls_[nprocs - 1] + recvCounts_[nprocs - 1];
            std::vector<HYPRE_BigInt> rkeys(2 * static_cast<std::size_t>(nrecv));
            std::vector<int> sc2(nprocs), sd2(nprocs), rc2(nprocs), rd2(nprocs);
            for (int o = 0; o < nprocs; ++o) {
                sc2[o] = 2 * sendCounts_[o];
                sd2[o] = 2 * sendDispls_[o];
                rc2[o] = 2 * recvCounts_[o];
                rd2[o] = 2 * recvDispls_[o];
            }
            MPI_Alltoallv(skeys.data(), sc2.data(), sd2.data(), HYPRE_MPI_BIG_INT,
                          rkeys.data(), rc2.data(), rd2.data(), HYPRE_MPI_BIG_INT, comm_);
            recv_.assign(nrecv, 0.0);
            std::vector<HYPRE_BigInt> rI(nrecv);
            for (int t = 0; t < nrecv; ++t) {
                rI[t] = rkeys[2 * static_cast<std::size_t>(t)];
            }
            resolve(rI, [&](std::size_t t) { return rkeys[2 * t + 1]; }, recv_dst_);
        }
        int okLocal = ok ? 1 : 0;
        int okAll = okLocal;
        MPI_Allreduce(&okLocal, &okAll, 1, MPI_INT, MPI_MIN, comm_);
        valid_ = okAll == 1;
        return valid_;
    }

    bool valid() const { return valid_; }

    void apply(hypre_ParCSRMatrix* A, hypre_ParCSRMatrix* P, hypre_ParCSRMatrix* AH)
    {
        const HYPRE_Real* ad = hypre_CSRMatrixData(hypre_ParCSRMatrixDiag(A));
        const HYPRE_Real* ao = hypre_CSRMatrixData(hypre_ParCSRMatrixOffd(A));
        const HYPRE_Real* pd = hypre_CSRMatrixData(hypre_ParCSRMatrixDiag(P));
        const HYPRE_Real* po = hypre_CSRMatrixData(hypre_ParCSRMatrixOffd(P));
        auto aval = [&](HYPRE_Int i) { return i < nnzAd_ ? ad[i] : ao[i - nnzAd_]; };
        auto pval = [&](HYPRE_Int i) {
            return i < nnzPd_ ? pd[i] : (i < nnzPd_ + nnzPo_ ? po[i - nnzPd_] : pext_[i - nnzPd_ - nnzPo_]);
        };
        std::fill(ap_.begin(), ap_.end(), 0.0);
        for (std::size_t t = 0; t < ap_dst_.size(); ++t) {
            ap_[ap_dst_[t]] += aval(ap_a_[t]) * pval(ap_p_[t]);
        }
        hypre_CSRMatrix* Hd = hypre_ParCSRMatrixDiag(AH);
        hypre_CSRMatrix* Ho = hypre_ParCSRMatrixOffd(AH);
        HYPRE_Real* hd = hypre_CSRMatrixData(Hd);
        HYPRE_Real* ho = hypre_CSRMatrixData(Ho);
        auto hslot = [&](HYPRE_Int h) -> HYPRE_Real& { return h < nnzHd_ ? hd[h] : ho[h - nnzHd_]; };
        std::fill(hd, hd + nnzHd_, 0.0);
        if (ho != nullptr) {
            std::fill(ho, ho + hypre_CSRMatrixNumNonzeros(Ho), 0.0);
        }
        for (std::size_t t = 0; t < r_dst_.size(); ++t) {
            hslot(r_dst_[t]) += pval(r_p_[t]) * ap_[r_ap_[t]];
        }
        if (recvCounts_.size() > 1) {
            std::fill(send_.begin(), send_.end(), 0.0);
            for (std::size_t t = 0; t < s_dst_.size(); ++t) {
                send_[s_dst_[t]] += pval(s_p_[t]) * ap_[s_ap_[t]];
            }
            MPI_Alltoallv(send_.data(), sendCounts_.data(), sendDispls_.data(), HYPRE_MPI_REAL,
                          recv_.data(), recvCounts_.data(), recvDispls_.data(), HYPRE_MPI_REAL, comm_);
            for (std::size_t t = 0; t < recv_dst_.size(); ++t) {
                hslot(recv_dst_[t]) += recv_[t];
            }
        }
    }

private:
    MPI_Comm comm_{};
    bool valid_ = false;
    HYPRE_Int nnzAd_ = 0, nnzPd_ = 0, nnzPo_ = 0, nnzHd_ = 0;
    std::vector<HYPRE_Real> pext_, ap_, send_, recv_;
    std::vector<HYPRE_Int> ap_a_, ap_p_, ap_dst_;   // AP[dst] += A[a] * P[p]
    std::vector<HYPRE_Int> r_p_, r_ap_, r_dst_;     // local A_H slot += P[p] * AP[ap]
    std::vector<HYPRE_Int> s_p_, s_ap_, s_dst_;     // send slot += P[p] * AP[ap]
    std::vector<HYPRE_Int> recv_dst_;               // local A_H slot of each received value
    std::vector<int> sendCounts_, sendDispls_, recvCounts_, recvDispls_;
};

} // namespace Opm::linalg::HypreInterface

#endif // OPM_HYPRE_NUMERIC_RAP_HPP
