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

#ifndef CONSTANTKERNEL_H
#define CONSTANTKERNEL_H

#include "Kernel.h"

//Forward Declarations
class ConstantKernel;
class Function;

template<>
InputParameters validParams<ConstantKernel>();

/**
 * This kernel implements a generic functional
 * body force term:
 * $ - c \cdof f \cdot \phi_i $
 *
 * The coefficient and function both have defaults
 * equal to 1.0.
 */
class ConstantKernel : public Kernel
{
public:
  ConstantKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  Real _t_limit;
  Real _val;
};

#endif
