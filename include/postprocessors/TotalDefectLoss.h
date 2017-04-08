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

#ifndef TOTALDEFECTLOSS_H
#define TOTALDEFECTLOSS_H

#include "GeneralPostprocessor.h"
#include "MaterialConstants.h"

// Forward Declarations
class TotalDefectLoss;
class MooseMesh;

namespace libMesh
{
class Node;
}

template<>
InputParameters validParams<TotalDefectLoss>();

/**
 * Sums a nodal value across all processors and multiplies the result
 * by a scale factor.
 */
class TotalDefectLoss : public GeneralPostprocessor
{
public:
  TotalDefectLoss(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute() {}

  /**
   * This will return the degrees of freedom in the system.
   */
  virtual Real getValue();

protected:
  MooseMesh & _mesh;
  std::string _var_prefix_v;
  std::string _var_prefix_i;
  std::vector<int> _mobile_v;
  std::vector<int> _mobile_i;
  int _v_max;
  int _i_max;
  Real & _total_loss;//global variable to track defect loss to recombination and sinks
  Node * _node_ptr;
  const MaterialConstants * const _material;
  Real _T;
};

#endif
