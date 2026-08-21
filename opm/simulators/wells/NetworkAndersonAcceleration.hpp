/*
  Copyright 2026 Equinor ASA.

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

#ifndef OPM_NETWORK_ANDERSON_ACCELERATION_HPP
#define OPM_NETWORK_ANDERSON_ACCELERATION_HPP

#include <cmath>
#include <cstddef>
#include <deque>
#include <vector>

namespace Opm {

/// Anderson acceleration of the network fixed-point iteration x <- G(x), where x
/// holds the node pressures of one network and G(x) is the pressure the network
/// gives for the rates the wells produce when x is applied as their THP.
///
/// The per-node update in BlackoilWellModelNetworkGeneric treats each node on its
/// own; on a tree where several leaves share an interior node that is wrong, and
/// the nodes fight each other. Anderson uses the last few (x, G(x)) pairs of the
/// whole vector, so it sees that coupling without needing any derivative.
///
/// Deliberately self-contained: no OPM dependencies, no state outside this class,
/// and a single call site. It is off by default.
template<typename Scalar>
class NetworkAndersonAccelerator
{
public:
    explicit NetworkAndersonAccelerator(const std::size_t depth = 4)
        : depth_(depth)
    {}

    void setDepth(const std::size_t depth)
    {
        depth_ = (depth > 0) ? depth : 1;
    }

    /// Next iterate from the applied pressures x and the network's answer gx.
    /// Falls back to returning gx (plain fixed-point) while there is no history,
    /// if the vector size changed, or if the least-squares problem is degenerate.
    std::vector<Scalar> next(const std::vector<Scalar>& x, const std::vector<Scalar>& gx)
    {
        const std::size_t n = x.size();
        if (n == 0 || gx.size() != n) {
            this->clear();
            return gx;
        }
        if (!x_.empty() && x_.back().size() != n) {
            this->clear();   // the node set changed
        }

        std::vector<Scalar> f(n);
        for (std::size_t i = 0; i < n; ++i) {
            f[i] = gx[i] - x[i];
        }

        std::vector<Scalar> next_x = gx;   // plain fixed-point step
        const std::size_t m = x_.size();   // number of stored differences
        if (m > 0) {
            // dF[j] = f_j - f_{j-1}, dX[j] = x_j - x_{j-1}, with the newest pair
            // formed against the incoming (x, f).
            std::vector<std::vector<Scalar>> dF(m), dX(m);
            for (std::size_t j = 0; j < m; ++j) {
                dF[j].resize(n);
                dX[j].resize(n);
                const auto& xj = x_[j];
                const auto& fj = f_[j];
                const auto& xn = (j + 1 < m) ? x_[j + 1] : x;
                const auto& fn = (j + 1 < m) ? f_[j + 1] : f;
                for (std::size_t i = 0; i < n; ++i) {
                    dF[j][i] = fn[i] - fj[i];
                    dX[j][i] = xn[i] - xj[i];
                }
            }
            // Regularised normal equations (dF^T dF + lambda I) gamma = dF^T f.
            std::vector<std::vector<Scalar>> A(m, std::vector<Scalar>(m + 1, Scalar{0}));
            Scalar trace{0};
            for (std::size_t a = 0; a < m; ++a) {
                for (std::size_t b = 0; b < m; ++b) {
                    Scalar s{0};
                    for (std::size_t i = 0; i < n; ++i) {
                        s += dF[a][i] * dF[b][i];
                    }
                    A[a][b] = s;
                    if (a == b) {
                        trace += s;
                    }
                }
                Scalar s{0};
                for (std::size_t i = 0; i < n; ++i) {
                    s += dF[a][i] * f[i];
                }
                A[a][m] = s;
            }
            const Scalar lambda = (trace > Scalar{0}) ? Scalar{1e-10} * trace / static_cast<Scalar>(m)
                                                      : Scalar{0};
            for (std::size_t a = 0; a < m; ++a) {
                A[a][a] += lambda;
            }
            std::vector<Scalar> gamma;
            if (solve(A, m, gamma)) {
                // x_{k+1} = G(x_k) - sum_j gamma_j (dX_j + dF_j)
                for (std::size_t j = 0; j < m; ++j) {
                    for (std::size_t i = 0; i < n; ++i) {
                        next_x[i] -= gamma[j] * (dX[j][i] + dF[j][i]);
                    }
                }
                for (const auto v : next_x) {
                    if (!std::isfinite(v)) {
                        next_x = gx;    // give up on this step, keep the history
                        break;
                    }
                }
            }
        }

        x_.push_back(x);
        f_.push_back(f);
        while (x_.size() > depth_) {
            x_.pop_front();
            f_.pop_front();
        }
        return next_x;
    }

    void clear()
    {
        x_.clear();
        f_.clear();
    }

private:
    /// Gaussian elimination with partial pivoting on the m x (m+1) augmented system.
    static bool solve(std::vector<std::vector<Scalar>>& A, const std::size_t m,
                      std::vector<Scalar>& out)
    {
        for (std::size_t c = 0; c < m; ++c) {
            std::size_t piv = c;
            for (std::size_t r = c + 1; r < m; ++r) {
                if (std::abs(A[r][c]) > std::abs(A[piv][c])) {
                    piv = r;
                }
            }
            if (!(std::abs(A[piv][c]) > Scalar{0})) {
                return false;
            }
            std::swap(A[c], A[piv]);
            for (std::size_t r = c + 1; r < m; ++r) {
                const Scalar w = A[r][c] / A[c][c];
                for (std::size_t k = c; k <= m; ++k) {
                    A[r][k] -= w * A[c][k];
                }
            }
        }
        out.assign(m, Scalar{0});
        for (std::size_t ri = 0; ri < m; ++ri) {
            const std::size_t r = m - 1 - ri;
            Scalar s = A[r][m];
            for (std::size_t k = r + 1; k < m; ++k) {
                s -= A[r][k] * out[k];
            }
            out[r] = s / A[r][r];
        }
        for (const auto v : out) {
            if (!std::isfinite(v)) {
                return false;
            }
        }
        return true;
    }

    std::size_t depth_;
    std::deque<std::vector<Scalar>> x_;
    std::deque<std::vector<Scalar>> f_;
};

} // namespace Opm

#endif // OPM_NETWORK_ANDERSON_ACCELERATION_HPP
