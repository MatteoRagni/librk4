# LibRK4

A simple step of Runge-Kutta 4 integration.

## Description

The library implements a Runge-Kutta 4 scheme with the following tableau:
```
      0  |  0   0   0   0
     1/2 | 1/2  0   0   0
     1/2 |  0  1/2  0   0
      1  |  0   0   1   0
     ----+----------------
         | 1/6 1/3 1/3 1/6
```
For the ode in the form:

<img src="https://latex.codecogs.com/gif.latex?\dot{x}(t)&space;=&space;f(t,&space;x(t),&space;u(t),&space;p),&space;\quad&space;x(t_0)&space;=&space;x_0" title="\dot{x}(t) = f(t, x(t), u(t), p), \quad x(t_0) = x_0" />

The integration scheme evaluates:

<img src="https://latex.codecogs.com/gif.latex?x(t_k&space;&plus;&space;h)&space;=&space;x(t_k)&space;&plus;&space;h&space;\sum_{i&space;=&space;1}^{4}&space;b_i&space;k_i" title="x(t_k + h) = x(t_k) + h \sum_{i = 1}^{4} b_i k_i" />

where:

<img src="https://latex.codecogs.com/gif.latex?k_i&space;=&space;f(t&space;&plus;&space;c_i&space;h,&space;x(t_k)&space;&plus;&space;h&space;\sum_{j=1}^{i-1}&space;a_{i,j}&space;k_j,&space;u(t_k),&space;p)" title="k_i = f(t + c_i h, x(t_k) + h \sum_{j=1}^{i-1} a_{i,j} k_j, u(t_k), p)" />

The input _u(t)_ is assumed constant between integration ministeps.
From an implementation point of view, _p_ is a vector of vectors. This
is done for compatibility with MATLAB System Identification Toolbox model
file. Since we usually don't have an initialization and a termination
callback, everything must be handled inside the integrator step (`malloc`
and `free` included).

## Examples

### Ode 1st order:

To represent an ODE in the form:

<img src="https://latex.codecogs.com/gif.latex?\dot{x}&space;=&space;-x&space;&plus;&space;u,&space;\quad&space;x(0)&space;=&space;0" title="\dot{x} = -x + u, \quad x(0) = 0" />

it is necessary to define the following vector field:

``` c
void ode(double *xdot,     /**< Output of the vector field */
         double t,         /**< Current time */
         const double *x,  /**< Current state vector */
         const double *u,  /**< Input vecotr */
         const double **p, /**< Vector of vectors of parameters */
         void *data)       /**< Space for user data */
{
  xdot[0] = -p[0][0] * x[0] + p[1][0] * u[0];
}
```
and the options for the integrator:
``` c
rk4_opts options = {
  1e-3, /**< Time steps */
  1,    /**< vector field dimensions */
  ode   /**< vector field callback (to evaluate it) */
};
```
a single integration step can be performed as follows:
``` c
rk4(&options, xp, t, x, u, p, NULL);
```
The last `NULL` is userdata placeholder (the callback does not use it).
The complete code is in [`test_ode_1`](test/test_ode_1.c). The result of the integration:

![Ode 1](.images/ode1.png?raw=true)

### Ode 2nd order:

To represent an ODE in the form:

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{rll}&space;\dot{x}_1&space;&=&space;x_2,&space;&&space;x_1(0)&space;=&space;0&space;\\&space;\dot{x}_2&space;&=&space;-&space;p_1\,&space;x_1&space;-&space;p_2\,&space;x_2&space;&plus;&space;p_3\,&space;u&space;,&space;&&space;x_2(0)&space;=&space;0&space;\\&space;\end{array}" title="\begin{array}{rll} \dot{x}_1 &= x_2, & x_1(0) = 0 \\ \dot{x}_2 &= - p_1\, x_1 - p_2\, x_2 + p_3\, u , & x_2(0) = 0 \\ \end{array}" />

it is necessary to define the following vector field:

``` c
void ode(double *xdot,     /**< Output of the vector field */
         double t,         /**< Current time */
         const double *x,  /**< Current state vector */
         const double *u,  /**< Input vecotr */
         const double **p, /**< Vector of vectors of parameters */
         void *data)       /**< Space for user data */
{
  xdot[0] = x[1];
  xdot[1] = -p[0][0] * x[0] - p[0][1] * x[1] + p[1][0] * u[0];
}
```
and the options for the integrator:
``` c
rk4_opts options = {
  1e-3, /**< Time steps */
  2,    /**< vector field dimensions */
  ode   /**< vector field callback (to evaluate it) */
};
```
a single integration step can be performed as follows:
``` c
rk4(&options, xp, t, x, u, p, NULL);
```
The complete code is in [`test_ode_2`](test/test_ode_2.c). The result of the integration:

![Ode 2](.images/ode2.png?raw=true)

### Non-linear ODE

To represent an ODE in the form:

<img src="https://latex.codecogs.com/gif.latex?\dot{x}&space;=&space;(1&space;-2\,&space;t)&space;\,x_1^2" title="\dot{x} = (1 -2\, t) \,x_1^2" />

it is necessary to define the following vector field (in this case we are evaluating the same ODE with 3 different initial conditions):

``` c
void ode(double *xdot,     /**< Output of the vector field */
         double t,         /**< Current time */
         const double *x,  /**< Current state vector */
         const double *u,  /**< Input vecotr */
         const double **p, /**< Vector of vectors of parameters */
         void *data)       /**< Space for user data */
{
  for (size_t i = 0; i < 3; i++)
    xdot[i] = (1 - 2 * t) * x[i] * x[i];
}
```
and the options for the integrator:
``` c
rk4_opts options = {
  1e-5, /**< Time steps */
  3,    /**< vector field dimensions */
  ode   /**< vector field callback (to evaluate it) */
};
```
a single integration step can be performed as follows (there is no input and no parameter):
``` c
rk4(&options, xp, t, x, NULL, NULL, NULL);
```
The complete code is in [`test_ode_3`](test/test_ode_3.c). The result of the integration:

![Ode 3](.images/ode3.png?raw=true)
