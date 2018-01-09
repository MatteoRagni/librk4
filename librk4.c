/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright 2018 - Matteo Ragni, Matteo Cocetti - University of Trento
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "librk4.h"

rk4_errno rk4(const rk4_opts *o, rk4_float *xp, const rk4_float t,
              const rk4_float *x, const rk4_float *u, const rk4_float **p, void *data) {

  ASSERT_PTR_NOT_NULL(o);
  ASSERT_PTR_NOT_NULL(xp);
  ASSERT_PTR_NOT_NULL(x);

  const rk4_float rk4_a[] = {0, 0.5, 0.5, 1};
  const rk4_float rk4_b[] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
  const rk4_float rk4_c[] = {0, 0.5, 0.5, 1};

  // Preparing the k matrix TODO: Static allocation of K
  rk4_float **k;
  k = (rk4_float **)malloc(RK4_ORDER * sizeof(rk4_float *));
  if (!k)
    return RK4_EMALLOC;
  for (size_t l = 0; l < RK4_ORDER; l++) {
    k[l] = (rk4_float *)calloc(o->f_size, sizeof(rk4_float));
    if (!(k[l])) {
      free(k);
      return RK4_EMALLOC;
    }
  }
  // Preparing support z vector
  rk4_float *z;
  z = (rk4_float *)calloc(o->f_size, sizeof(rk4_float));
  if (!z) {
    // Cleaning up if allocation fails
    for (size_t l = 0; l < RK4_ORDER; l++)
      free(k[l]);
    free(k);
    return RK4_EMALLOC;
  }

  // Evaluating k[0]
  o->f(k[0], t, x, u, p, data);
  // Evaluating other k[i]
  for (size_t l = 1; l < RK4_ORDER; l++) {
    for (size_t i = 0; i < o->f_size; i++)
      z[i] = x[i] + (o->h) * rk4_a[l - 1] * k[l - 1][i];
    for (size_t i = 0; i < o->f_size; i++) {
      o->f(k[l], t + rk4_c[l] * (o->h), z, u, p, data);
    }
  }
  // Evaluating output
  for (size_t i = 0; i < o->f_size; i++) {
    double df = 0;
    for (size_t l = 0; l < RK4_ORDER; l++)
      df += rk4_b[l] * k[l][i];
    xp[i] = x[i] + (o->h) * df;
  }

  // Clearing up k matrix
  for (size_t l = 0; l < RK4_ORDER; l++)
    free(k[l]);
  free(k);
  free(z);
  return RK4_SUCCESS;
}
