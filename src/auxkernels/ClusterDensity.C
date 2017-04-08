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

//calculate void swelling from variable vaules with discrete cluster dynamics method 
#include "ClusterDensity.h"

template<>
InputParameters validParams<ClusterDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_vars","coupled vacancy type variables");
  params.addParam<Real>("scaling_factor",1.0,"scaling factor to overal cluster density, eg. atomic volume");
  return params;
}


ClusterDensity::ClusterDensity(const
                                   InputParameters & parameters)
  :AuxKernel(parameters),
  _scaling_factor(getParam<Real>("scaling_factor"))
{
  int ncoupled = coupledComponents("coupled_vars");
  _no_vars.resize(ncoupled);
  _val_vars.resize(ncoupled);

  for (int i=0; i < ncoupled; ++i)
  {
    _no_vars[i] = coupled("coupled_vars",i);
    _val_vars[i] = &coupledValue("coupled_vars",i);
  }
}

Real
ClusterDensity::computeValue()
{
  Real total = 0.0;//total defect cluster conentration
  int no_vars = _no_vars.size() ;
  for(int i=0;i<no_vars;i++){
      total += (*_val_vars[i])[_qp];
  }

  return total * _scaling_factor;
}

