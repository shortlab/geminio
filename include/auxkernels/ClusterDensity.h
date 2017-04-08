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

#ifndef CLUSTERDENSITY_H
#define CLUSTERDENSITY_H

#include "AuxKernel.h"

//Forward Declarations
class ClusterDensity;

/**
 * validParams returns the parameters that this Kernel accepts / needs
 * The actual body of the function MUST be in the .C file.
 */
template<>
InputParameters validParams<ClusterDensity>();

class ClusterDensity : public AuxKernel
{
public:

  ClusterDensity(const 
                   InputParameters & parameters);

protected:
  virtual Real computeValue();

  /**
   * This MooseArray will hold the reference we need to our
   * material property from the Material class
   */

  Real _scaling_factor;
  std::vector<unsigned int> _no_vars;
  std::vector<const VariableValue *> _val_vars;

};
#endif
