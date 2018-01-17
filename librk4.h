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
#ifndef LIBRK4_H_
#define LIBRK4_H_

#include <stddef.h>
#include <stdlib.h>

/**
 * @mainpage
 * @author Matteo Ragni, Matteo Cocetti
 * @date 4 Jan 2018
 *
 * The library implements a Runge-Kutta 4 scheme with the following tableau:
 * @code
 *      0  |  0   0   0   0
 *     1/2 | 1/2  0   0   0
 *     1/2 |  0  1/2  0   0
 *      1  |  0   0   1   0
 *     ----+----------------
 *         | 1/6 1/3 1/3 1/6
 * @endcode
 * For the ode in the form:
 *  \f[
 *    \dot{x}(t) = f(t, x(t), u(t), p), \quad
 *    x(t_0) = x_0
 *  \f]
 * The integration scheme evaluates:
 *  \f[
 *    x(t_k + h) = x(t_k) + h \sum_{i = 1}^{4} b_i k_i
 *  \f]
 * where:
 *  \f[
 *    k_i = f(t + c_i h, x(t_k) + h \sum_{j=1}^{i-1} a_{i,j} k_j, u(t_k), p)
 *  \f]
 *
 * The input \f$ u(t_k) \f$ is assumed constant between integration ministeps.
 * From an implementation point of view, \f$ p \f$ is a vector of vectors. This
 * is done for compatibility with MATLAB System Identification Toolbox model
 * file. Since we usually don't have an initialization and a termination
 * callback, everything must be handled inside the integrator step (malloc
 * and free included).
 */

#ifndef RK4_FLOAT_TYPE
#define RK4_FLOAT_TYPE                                                         \
  double /**< Pecision type override. Default is double.                       \
          */
#endif
typedef RK4_FLOAT_TYPE rk4_float; /**< Precision type declaration */

/**
 * @brief Ode callback function: implements \f$ f(t, x, u, p) \f$.
 *
 * The callback for the ODE model. The callback stores the output in
 * the first pointer (xdot). The function does not need to allocate nor free
 * the input vector, but may overflow if exceedes the dimension that is
 * declared in the rk4_opts structure.
 * When accessing an element of the parameter, please remember that the
 * implementation keeps in mind a MATLAB-like interface. If the parameters
 * (in MATLAB) are passed like:
 * @code
 * f(t, x, u, p1, [p2, p3])
 * @endcode
 * then, in the C code, p1, p2 and p3 may be accessed in p as:
 * @code
 * rk4_float p1, p2, p3;
 * p1 = p[0][0];
 * p2 = p[1][0];
 * p3 = p[1][1];
 * @endcode
 * If you need the time step **inside** the callback, you should directly define
 * it inside or pass it as a parameter (and remeber to update it in the rk4_opts
 * structure).
 * @param xdot otput vector for the ODE
 * @param t current time for the function
 * @param x state for the ODE evaluation
 * @param u external input for the ODE evaluation
 * @param p pointer to arrays of parameters
 * @param auxiliary data pointer to void for user data
 * @returns nothing
 * @warning The callback may produce a buffer overflow error if tries to
 *          write more than x_size elements in xdot.
 */
typedef void (*rk4_ode)(rk4_float *xdot, rk4_float t, const rk4_float *x,
                        const rk4_float *u, const rk4_float **p, void *data);

/**
 * @brief Return values for the integration step
 */
typedef enum rk4_errorcode {
  RK4_SUCCESS = 0, /**< Step seems good! */
  RK4_EMALLOC,     /**< During the step there was an allocation error */
  RK4_NULLPTR,     /**< Reveived a null pointer */
  RK4_GENERIC      /**< Unknown error generated */
} rk4_errorcode;

/**
 * @brief Assertion for pointer check. Can be disabled
 *
 * The next macro contains an assertion for checking for input
 * pointers. Can be disabled to speed up the code. If the assertion
 * fails, it makes rk4 function to return a RK4_NULLPTR. It may be disabled
 * after debugging phase, removing some if in the code.
 * @warning control (u) and parameter (p) pointer are not checked! This
 *          allows to pass NULL vector if the integrator step does not
 *          depend upon those two. Also user data is not checked.
 */
#define RK4_CHECK_NULL_PTRS
#ifdef RK4_CHECK_NULL_PTRS
#define ASSERT_PTR_NOT_NULL(ptrs)                                       \
  if (!(ptrs))                                                        \
    return RK4_NULLPTR;
#else
#define ASSERT_PTR_NOT_NULL(ptrs) ;
#endif

#define RK4_ORDER 4 /**< Order of the RK4 integrator (always 4) */

/**
 * @brief Integrator options.
 *
 * The structure contains some informations that are necessary for integration.
 * When used in MATLAB the structure may be initialized as a constant. The
 * structur contains the dimension of the vector field, the fundamental step,
 * and the vector field callback.
 */
typedef struct rk4_opts {
  rk4_float h;  /**< Fundamental time step. */
  size_t f_size; /**< Dimension for the vector field. */
  rk4_ode f;     /**< ODE vector field. This is a function pointer. */
} rk4_opts;

/**
 * @brief Performs an RK4 step.
 *
 * The function performs a step for the RK4 integration. The function receives
 * the vector field in the rk4_opts structure, and stores the result of the step
 * in the xp array.
 * @param o options structure for the integrator.
 * @param xp step result
 * @param t current time for the function
 * @param x state for the ODE evaluation
 * @param u external input for the ODE evaluation
 * @param p pointer to arrays of parameters
 * @param data space for user data, passed to ode callback
 * @return an exit error number
 */
rk4_errorcode rk4(const rk4_opts *o, rk4_float *xp, const rk4_float t,
              const rk4_float *x, const rk4_float *u, const rk4_float **p, void* data);

#endif /* LIBRK4_H_ */
