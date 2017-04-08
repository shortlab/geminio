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

#ifndef GSUMSIACLUSTERDENSITY_H
#define GSUMSIACLUSTERDENSITY_H

#include "AuxKernel.h"
#include "GGroup.h"

//Forward Declarations
class GSumSIAClusterDensity;

/**
 * validParams returns the parameters that this Kernel accepts / needs
 * The actual body of the function MUST be in the .C file.
 */
template<>
InputParameters validParams<GSumSIAClusterDensity>();

class GSumSIAClusterDensity : public AuxKernel
{
public:

  GSumSIAClusterDensity(const 
                   InputParameters & parameters);

protected:
  virtual Real computeValue();

  /**
   * This MooseArray will hold the reference we need to our
   * material property from the Material class
   */

  const GGroup & _gc;
  Real _scale_factor;
  int _lower_bound;
  int _upper_bound;
  std::vector<unsigned int> _no_vars;
  std::vector<const VariableValue *> _val_vars;

};
#endif
