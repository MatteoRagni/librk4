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

/**
 * @brief Second order ODE with step input.
 *
 * The test evaluates the integral of a first order ODE with the following
 * formulation:
 *  \f{eqnarray*}{
 *    \dot{x}_1(t) &= x_2(t) \\
 *    \dot{x}_2(t) &= -k x_1(t) - c x_2(t) + u(t)
 *  \f}
 * and stores the result in a file. The integration is performed with
 * \f$h = 1\cdot10^{-3}\f$.
 */

#include <stdio.h>
#include <stdlib.h>
#include "../librk4.h"


/**
 * @brief ODE declaration.
 */
void ode(double *xdot, double t, const double *x, const double *u, const double **p) {
  xdot[0] = x[1];
  xdot[1] = -p[0][0] * x[0] - p[0][1] * x[1] + p[1][0] * u[0];
}

/**
 * @brief Step function.
 */
double step(double t, double st) { return t >= st ? 1.0 : 0.0; }

/**
 * @brief Integration options.
 */
rk4_opts options = { 1e-3, 2, ode };

char result[] = "test_ode_2.csv";

int main() {
  FILE *fp = fopen(result, "w");

  if (!fp)
    return 1;

  double t = 0;
  double u[] = { 0.0 };

  double **p = (double **)malloc(2 * sizeof(double *));
  p[0] = (double *)calloc(2, sizeof(double));
  p[1] = (double *)calloc(1, sizeof(double));
  p[0][0] = 1.0; /**< k */
  p[0][1] = 0.5; /**< c */
  p[1][0] = 1.0;

  double x[] = { 0.0, 0.0 };
  double xp[] = { 0.0, 0.0 };

  while (t <= 10.0) {

    fprintf(fp, "% 1.5f,% 3.5f,% 3.5f, % 3.5f\n", t, u[0], x[0], x[1]);

    t += options.h;
    u[0] = step(t, 1.0);
    rk4(&options, xp, t, x, u, (const double **)p);
    x[0] = xp[0];
    x[1] = xp[1];
  }

  fclose(fp);
  return 0;
}
