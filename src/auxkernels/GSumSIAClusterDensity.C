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

//calculate total density of SIA cluster in range [lower_bounds,upper_bound] from variable vaules with grouping method
#include "GSumSIAClusterDensity.h"

template<>
InputParameters validParams<GSumSIAClusterDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_vars","coupled variables");
  params.addParam<int>("lower_bound",1,"starting size to count, inclusive");
  params.addParam<int>("upper_bound","ending size to count, inclusive");
  params.addParam<Real>("scale_factor", 1, "A scale factor to be applied to the variable");
  params.addRequiredParam<UserObjectName>("user_object","The name of user object providing interaction constants");
  return params;
}


GSumSIAClusterDensity::GSumSIAClusterDensity(const
                                   InputParameters & parameters)
  :AuxKernel(parameters),
  _gc(getUserObject<GGroup>("user_object")),
  _scale_factor(getParam<Real>("scale_factor")),
  _lower_bound(getParam<int>("lower_bound")),//[lower_bound,upper_bound],inclusive
  _upper_bound(isParamValid("upper_bound")?getParam<int>("upper_bound"):_gc.GroupScheme_i.back())
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
GSumSIAClusterDensity::computeValue()
{
  Real total_density = 0.0;//total cluster density in range [_lower_bound,_upper_bound]
  for(int i_size=_lower_bound;i_size<=_upper_bound;i_size++){
    int g = _gc.CurrentGroupI(i_size);
    total_density += (*_val_vars[2*(g-1)])[_qp]+(*_val_vars[2*(g-1)+1])[_qp]*(i_size-_gc.GroupScheme_i_avg[g-1]);
  }
  return total_density*_scale_factor;
}

