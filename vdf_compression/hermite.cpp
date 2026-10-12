/*
 * Adaptive orthonormal Hermite transform of a dense 3D velocity cube.
 * Port of faiser/data_processing/vdf_tools.py; see hermite.h for conventions.
 */
#include "hermite.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace ASTERIX {
namespace hermite {

namespace {

// Orthonormal basis phi_n(v) = H_n(x) exp(-x^2/2) / sqrt(2^n n! sqrt(pi) vth),
// x = (v-u)/vth. Returns flat [order][vlen] array.
std::vector<double> hermite_basis(const std::vector<double>& v_ax, int order, double vth, double u) {
   const std::size_t vlen = v_ax.size();
   std::vector<double> H(static_cast<std::size_t>(order) * vlen, 0.0);
   if (order <= 0) return H;
   for (std::size_t i = 0; i < vlen; ++i) {
      const double x = (v_ax[i] - u) / vth;
      H[i] = std::exp(-0.5 * x * x);
      if (order > 1) H[vlen + i] = 2.0 * x * H[i];
   }
   for (int n = 2; n < order; ++n) {
      for (std::size_t i = 0; i < vlen; ++i) {
         const double x = (v_ax[i] - u) / vth;
         H[n * vlen + i] = 2.0 * x * H[(n - 1) * vlen + i] - 2.0 * (n - 1) * H[(n - 2) * vlen + i];
      }
   }
   // norm_n^2 = 2^n n! sqrt(pi) vth, built iteratively (norm_n^2 = 2n norm_{n-1}^2)
   double norm2 = std::sqrt(std::numbers::pi) * vth;
   for (int n = 0; n < order; ++n) {
      if (n > 0) norm2 *= 2.0 * n;
      const double inv = 1.0 / std::sqrt(norm2);
      for (std::size_t i = 0; i < vlen; ++i) H[n * vlen + i] *= inv;
   }
   return H;
}

// Cell-centre axis d of a grid.
std::vector<double> grid_axis(const Grid& g, int d) {
   std::vector<double> ax(g.shape[d]);
   for (int i = 0; i < g.shape[d]; ++i) ax[i] = g.v0[d] + i * g.dv[d];
   return ax;
}

std::size_t n_cells(const Grid& g) {
   return static_cast<std::size_t>(g.shape[0]) * g.shape[1] * g.shape[2];
}

bool valid(const Grid& g) { return g.shape[0] > 0 && g.shape[1] > 0 && g.shape[2] > 0; }

Grid cubic_grid(double vlim, int vlen) {
   const double dv = 2.0 * vlim / vlen;
   const double v0 = -vlim + dv / 2.0;
   return Grid{{v0, v0, v0}, {dv, dv, dv}, {vlen, vlen, vlen}};
}

} // namespace

std::size_t n_coeffs(int order) {
   if (order <= 0) return 0;
   const std::size_t o = static_cast<std::size_t>(order);
   return o * (o + 1) * (o + 2) / 6;
}

std::vector<double> velocity_axis(double vlim, int vlen) {
   const double dv = 2.0 * vlim / vlen;
   std::vector<double> ax(vlen);
   for (int i = 0; i < vlen; ++i) ax[i] = -vlim + dv / 2.0 + i * dv;
   return ax;
}

std::array<double, 3> drift_velocity(const std::vector<double>& cube, const Grid& grid) {
   if (!valid(grid) || cube.size() != n_cells(grid)) return {0.0, 0.0, 0.0};
   const double dv3 = grid.dv[0] * grid.dv[1] * grid.dv[2];
   const auto axx = grid_axis(grid, 0), axy = grid_axis(grid, 1), axz = grid_axis(grid, 2);
   const std::size_t nx = grid.shape[0], ny = grid.shape[1], nz = grid.shape[2];
   double sum = 0.0, sx = 0.0, sy = 0.0, sz = 0.0;
   for (std::size_t iz = 0; iz < nz; ++iz)
      for (std::size_t iy = 0; iy < ny; ++iy)
         for (std::size_t ix = 0; ix < nx; ++ix) {
            const double f = cube[(iz * ny + iy) * nx + ix];
            sum += f;
            sx += f * axx[ix];
            sy += f * axy[iy];
            sz += f * axz[iz];
         }
   const double n = sum * dv3;
   if (n == 0.0) return {0.0, 0.0, 0.0};
   return {sx * dv3 / n, sy * dv3 / n, sz * dv3 / n};
}

double thermal_velocity(const std::vector<double>& cube, const Grid& grid, const std::array<double, 3>& u) {
   if (!valid(grid) || cube.size() != n_cells(grid)) return 0.0;
   const auto axx = grid_axis(grid, 0), axy = grid_axis(grid, 1), axz = grid_axis(grid, 2);
   const std::size_t nx = grid.shape[0], ny = grid.shape[1], nz = grid.shape[2];
   // 1-D axis projections (sums over the other two axes)
   std::vector<double> px(nx, 0.0), py(ny, 0.0), pz(nz, 0.0);
   double sum = 0.0;
   for (std::size_t iz = 0; iz < nz; ++iz)
      for (std::size_t iy = 0; iy < ny; ++iy)
         for (std::size_t ix = 0; ix < nx; ++ix) {
            const double f = cube[(iz * ny + iy) * nx + ix];
            sum += f;
            px[ix] += f;
            py[iy] += f;
            pz[iz] += f;
         }
   if (!(sum > 0.0)) return 0.0;
   double pxx = 0.0, pyy = 0.0, pzz = 0.0;
   for (std::size_t i = 0; i < nx; ++i) pxx += px[i] * (axx[i] - u[0]) * (axx[i] - u[0]);
   for (std::size_t i = 0; i < ny; ++i) pyy += py[i] * (axy[i] - u[1]) * (axy[i] - u[1]);
   for (std::size_t i = 0; i < nz; ++i) pzz += pz[i] * (axz[i] - u[2]) * (axz[i] - u[2]);
   // P_dd / n: the cell volume cancels
   return std::sqrt((pxx / sum + pyy / sum + pzz / sum) / 3.0);
}

Spectrum adaptive_transform(const std::vector<double>& cube, const Grid& grid, double vth,
                            const std::array<double, 3>& u, int max_order, double tolerance) {
   Spectrum out;
   if (max_order <= 0 || !valid(grid) || cube.size() != n_cells(grid)) return out;
   if (!(vth > 0.0)) return out;

   const std::size_t nx = grid.shape[0], ny = grid.shape[1], nz = grid.shape[2];
   const double dv3 = grid.dv[0] * grid.dv[1] * grid.dv[2];

   double P_total = 0.0;
   for (double f : cube) P_total += f * f;
   P_total *= dv3;
   if (!(P_total > 0.0)) return out; // empty VDF: nothing to transform

   const auto Hx = hermite_basis(grid_axis(grid, 0), max_order, vth, u[0]);
   const auto Hy = hermite_basis(grid_axis(grid, 1), max_order, vth, u[1]);
   const auto Hz = hermite_basis(grid_axis(grid, 2), max_order, vth, u[2]);

   // Separable contraction with lazily filled caches:
   //   Px[n][z,y]  = sum_x cube[z,y,x] Hx[n,x]      O(N^3), once per n
   //   Pxy[n,m][z] = sum_y Px[n][z,y] Hy[m,y]       O(N^2), once per (n,m)
   //   c[l,m,n]    = sum_z Pxy[n,m][z] Hz[l,z] dv^3 O(N)
   const std::size_t O = max_order;
   std::vector<std::vector<double>> Px(O);
   std::vector<std::vector<double>> Pxy(O * O);

   auto get_Px = [&](int n) -> const std::vector<double>& {
      auto& p = Px[n];
      if (p.empty()) {
         p.assign(ny * nz, 0.0);
         const double* h = &Hx[n * nx];
         for (std::size_t zy = 0; zy < ny * nz; ++zy) {
            const double* row = &cube[zy * nx];
            double acc = 0.0;
            for (std::size_t ix = 0; ix < nx; ++ix) acc += row[ix] * h[ix];
            p[zy] = acc;
         }
      }
      return p;
   };
   auto get_Pxy = [&](int n, int m) -> const std::vector<double>& {
      auto& p = Pxy[static_cast<std::size_t>(n) * O + m];
      if (p.empty()) {
         const auto& px = get_Px(n);
         p.assign(nz, 0.0);
         const double* h = &Hy[m * ny];
         for (std::size_t iz = 0; iz < nz; ++iz) {
            double acc = 0.0;
            for (std::size_t iy = 0; iy < ny; ++iy) acc += px[iz * ny + iy] * h[iy];
            p[iz] = acc;
         }
      }
      return p;
   };

   double P_accum = 0.0;
   out.coeffs.reserve(n_coeffs(max_order));
   for (int s = 0; s < max_order; ++s) {
      for (int l = 0; l <= s; ++l) {
         for (int m = 0; m <= s - l; ++m) {
            const int n = s - l - m;
            const auto& pxy = get_Pxy(n, m);
            const double* hz = &Hz[l * nz];
            double acc = 0.0;
            for (std::size_t iz = 0; iz < nz; ++iz) acc += pxy[iz] * hz[iz];
            const double c = acc * dv3;
            out.coeffs.push_back(c);
            P_accum += c * c;
         }
      }
      const double delta = std::sqrt(std::max(0.0, 1.0 - P_accum / P_total));
      out.deltas.push_back(delta);
      if (s > 2 && delta < tolerance) {
         out.order_used = s + 1;
         return out;
      }
   }
   out.order_used = max_order;
   return out;
}

std::vector<double> reconstruct(const Spectrum& spectrum, const Grid& grid, double vth,
                                const std::array<double, 3>& u) {
   if (!valid(grid)) return {};
   const std::size_t nx = grid.shape[0], ny = grid.shape[1], nz = grid.shape[2];
   std::vector<double> out(nx * ny * nz, 0.0);
   const int order = spectrum.order_used;
   if (order <= 0 || spectrum.coeffs.size() != n_coeffs(order) || !(vth > 0.0)) return out;

   const auto Hx = hermite_basis(grid_axis(grid, 0), order, vth, u[0]);
   const auto Hy = hermite_basis(grid_axis(grid, 1), order, vth, u[1]);
   const auto Hz = hermite_basis(grid_axis(grid, 2), order, vth, u[2]);
   const std::size_t O = order;

   // Dense h[l][m][n] from the canonical coefficient order
   std::vector<double> h(O * O * O, 0.0);
   std::size_t k = 0;
   for (int s = 0; s < order; ++s)
      for (int l = 0; l <= s; ++l)
         for (int m = 0; m <= s - l; ++m) h[(l * O + m) * O + (s - l - m)] = spectrum.coeffs[k++];

   // Pass 1: T1[l,m,x] = sum_n h[l,m,n] Hx[n,x]
   std::vector<double> T1(O * O * nx, 0.0);
   for (std::size_t lm = 0; lm < O * O; ++lm)
      for (std::size_t n = 0; n < O; ++n) {
         const double c = h[lm * O + n];
         if (c == 0.0) continue;
         for (std::size_t ix = 0; ix < nx; ++ix) T1[lm * nx + ix] += c * Hx[n * nx + ix];
      }
   // Pass 2: T2[l,y,x] = sum_m Hy[m,y] T1[l,m,x]
   std::vector<double> T2(O * ny * nx, 0.0);
   for (std::size_t l = 0; l < O; ++l)
      for (std::size_t m = 0; m < O; ++m)
         for (std::size_t iy = 0; iy < ny; ++iy) {
            const double w = Hy[m * ny + iy];
            for (std::size_t ix = 0; ix < nx; ++ix) T2[(l * ny + iy) * nx + ix] += w * T1[(l * O + m) * nx + ix];
         }
   // Pass 3: out[z,y,x] = sum_l Hz[l,z] T2[l,y,x]
   for (std::size_t l = 0; l < O; ++l)
      for (std::size_t iz = 0; iz < nz; ++iz) {
         const double w = Hz[l * nz + iz];
         double* dst = &out[iz * ny * nx];
         const double* src = &T2[l * ny * nx];
         for (std::size_t i = 0; i < ny * nx; ++i) dst[i] += w * src[i];
      }
   return out;
}

// ---- Cubic API: thin wrappers over the Grid API -------------------------------

std::array<double, 3> drift_velocity(const std::vector<double>& cube, double vlim, int vlen) {
   if (vlen <= 0) return {0.0, 0.0, 0.0};
   return drift_velocity(cube, cubic_grid(vlim, vlen));
}

double thermal_velocity(const std::vector<double>& cube, double vlim, int vlen, const std::array<double, 3>& u) {
   if (vlen <= 0) return 0.0;
   return thermal_velocity(cube, cubic_grid(vlim, vlen), u);
}

Spectrum adaptive_transform(const std::vector<double>& cube, double vlim, int vlen, double vth,
                            const std::array<double, 3>& u, int max_order, double tolerance) {
   if (vlen <= 0) return Spectrum{};
   return adaptive_transform(cube, cubic_grid(vlim, vlen), vth, u, max_order, tolerance);
}

std::vector<double> reconstruct(const Spectrum& spectrum, double vlim, int vlen, double vth,
                                const std::array<double, 3>& u) {
   if (vlen <= 0) return {};
   return reconstruct(spectrum, cubic_grid(vlim, vlen), vth, u);
}

} // namespace hermite
} // namespace ASTERIX
