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
 * @brief Non linear, time variant, ODE with no input
 *
 * The test evaluates the integral of a first order ODE with the following
 * formulation:
 *  \f[
 *    \dot{x}(t) = (1 - 2 t) x^2
 *  \f]
 * and stores the result in a file. The ODE has dimension 3 because we
 * are passing 3 different initial conditions. There is no input and
 * no parameter.
 */

#include <stdio.h>
#include <stdlib.h>
#include "../librk4.h"


/**
 * @brief ODE declaration.
 */
void ode(double *xdot, double t, const double *x, const double *u, const double **p, void *q) {
  for (size_t i = 0; i < 3; i++)
    xdot[i] = (1 - 2 * t) * (x[i] * x[i]);
}

/**
 * @brief Integration options.
 *
 * The first option is the integration step. The second option is
 * the dimension of the ODE. The third option is the ode callback.
 */
rk4_opts options = { 1e-5, 3, ode };

char result[] = "test_ode_3.csv";

int main() {
  FILE *fp = fopen(result, "w");

  if (!fp)
    return 1;

  double t = 0;

  double x[] = { 1.0, 2.0, 3.0 };
  double xp[] = { 0.0, 0.0, 0.0 };

  while (t <= 3.0) {

    fprintf(fp, "% 1.5f,% 3.5f,% 3.5f, % 3.5f\n", t, x[0], x[1], x[2]);

    t += options.h;
    rk4(&options, xp, t, x, NULL, NULL, NULL);
    x[0] = xp[0];
    x[1] = xp[1];
    x[2] = xp[2];
  }

  fclose(fp);
  return 0;
}
