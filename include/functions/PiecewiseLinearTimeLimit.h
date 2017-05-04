/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef PIECEWISELINEARTIMELIMIT_H
#define PIECEWISELINEARTIMELIMIT_H

#include "Piecewise.h"

// Forward declarations
class PiecewiseLinearTimeLimit;

template<>
InputParameters validParams<PiecewiseLinearTimeLimit>();

/**
 * Base class for function objects.  Functions override value to supply a
 * value at a point.
 */
class PiecewiseLinearTimeLimit : public Piecewise
{
public:
  PiecewiseLinearTimeLimit(const InputParameters & parameters);
  virtual ~PiecewiseLinearTimeLimit();

  /**
   * Get the value of the function (based on time only)
   * \param t The time
   * \param pt The point in space (x,y,z) (unused)
   * \return The value of the function at the specified time
   */
  virtual Real value(Real t, const Point & pt);

  /**
   * Get the time derivative of the function (based on time only)
   * \param t The time
   * \param pt The point in space (x,y,z) (unused)
   * \return The time derivative of the function at the specified time
   */
  virtual Real timeDerivative(Real t, const Point & pt);

  virtual Real integral();

  virtual Real average();
  protected:
  Real _t_limit;
};

#endif
