#include "math.hpp"

// taken from pbrt
std::optional<Matrix4f> inverse(const Matrix4f &m) {
  int indxc[4], indxr[4];
  int ipiv[4] = {0};
  real minv[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      minv[i][j] = m.m[i][j];
  for (int i = 0; i < 4; i++) {
    int irow = 0, icol = 0;
    real big = 0.f;
    // Choose pivot
    for (int j = 0; j < 4; j++) {
      if (ipiv[j] != 1) {
        for (int k = 0; k < 4; k++) {
          if (ipiv[k] == 0) {
            if (std::abs(minv[j][k]) >= big) {
              big = std::abs(minv[j][k]);
              irow = j;
              icol = k;
            }
          } else if (ipiv[k] > 1)
            return {};  // singular
        }
      }
    }
    ++ipiv[icol];
    // Swap rows _irow_ and _icol_ for pivot
    if (irow != icol) {
      for (int k = 0; k < 4; ++k)
        std::swap(minv[irow][k], minv[icol][k]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (minv[icol][icol] == 0.f)
      return {};  // singular

    // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
    real pivinv = 1. / minv[icol][icol];
    minv[icol][icol] = 1.;
    for (int j = 0; j < 4; j++)
      minv[icol][j] *= pivinv;

    // Subtract this row from others to zero out their columns
    for (int j = 0; j < 4; j++) {
      if (j != icol) {
        real save = minv[j][icol];
        minv[j][icol] = 0;
        for (int k = 0; k < 4; k++)
          minv[j][k] = std::fma(-minv[icol][k], save, minv[j][k]);
      }
    }
  }
  // Swap columns to reflect permutation
  for (int j = 4 - 1; j >= 0; j--) {
    if (indxr[j] != indxc[j]) {
      for (int k = 0; k < 4; k++)
        std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
    }
  }
  return Matrix4f(minv);
}
