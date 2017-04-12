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

#ifndef NODALCONSERVATIONCHECK_H
#define NODALCONSERVATIONCHECK_H

#include "GeneralPostprocessor.h"
#include "MooseVariable.h"

// Forward Declarations
class NodalConservationCheck;
class MooseMesh;

namespace libMesh
{
class Node;
}

template<>
InputParameters validParams<NodalConservationCheck>();

/**
 * Sums a nodal value across all processors and multiplies the result
 * by a scale factor.
 */
class NodalConservationCheck : public GeneralPostprocessor
{
public:
  NodalConservationCheck(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute() {}

  /**
   * This will return the degrees of freedom in the system.
   */
  virtual Real getValue();

protected:
  MooseMesh & _mesh;
  std::string _var_prefix;
  std::vector<int> _size_range;//cluster size range
  Node * _node_ptr;
  const Real _scale_factor;
};

#endif
