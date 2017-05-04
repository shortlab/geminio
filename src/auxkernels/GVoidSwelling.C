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

//calculate void swelling from variable vaules with grouping method
#include "GVoidSwelling.h"

template<>
InputParameters validParams<GVoidSwelling>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_v_vars","coupled vacancy type variables");
  params.addParam<int>("lower_bound",1,"starting size to count, inclusive");
  params.addParam<int>("upper_bound","ending size to count, inclusive");
  params.addParam<Real>("scale_factor",1.0,"Scaling factor applied to the variable");
  params.addRequiredParam<UserObjectName>("user_object","The name of user object providing interaction constants");
  return params;
}


GVoidSwelling::GVoidSwelling(const
                                   InputParameters & parameters)
  :AuxKernel(parameters),
  _gc(getUserObject<GGroup>("user_object")),
  _scale_factor(getParam<Real>("scale_factor")),
  _lower_bound(getParam<int>("lower_bound")),//[lower_bound,upper_bound],inclusive
  _upper_bound(isParamValid("upper_bound")?getParam<int>("upper_bound"):_gc.GroupScheme_v.back())
{
  int nvcoupled = coupledComponents("coupled_v_vars");
  _no_v_vars.resize(nvcoupled);
  _val_v_vars.resize(nvcoupled);

  for (int i=0; i < nvcoupled; ++i)
  {
    _no_v_vars[i] = coupled("coupled_v_vars",i);
    _val_v_vars[i] = &coupledValue("coupled_v_vars",i);
  }
}

Real
GVoidSwelling::computeValue()
{
  Real total_vacancy = 0.0;//total conentration in terms of point defects
  int low = _gc.CurrentGroupV(_lower_bound)-1;
  int high = _gc.CurrentGroupV(_upper_bound)-1;
  for(int i=low;i<high;i++){
  /*
    for(int j=_gc.GroupScheme_v[i]+1;j<=_gc.GroupScheme_v[i+1];j++){//(a,b]
      total_vacancy += ((*_val_v_vars[2*i])[_qp]+(*_val_v_vars[2*i+1])[_qp]*(j-_gc.GroupScheme_v_avg[i]))*j;
    }
  */
     Real a = (*_val_v_vars[2*i])[_qp];
     Real b = (*_val_v_vars[2*i+1])[_qp];
     Real avg = _gc.GroupScheme_v_avg[i];
     total_vacancy += ( (a-b*avg)*_gc.GroupScheme_v_del[i]*(_gc.GroupScheme_v[i]+1+_gc.GroupScheme_v[i+1])/2.0 + b * (_gc.GroupScheme_v[i+1]*(_gc.GroupScheme_v[i+1]+1.0)*(2.0*_gc.GroupScheme_v[i+1]+1) - _gc.GroupScheme_v[i]*(_gc.GroupScheme_v[i]+1.0)*(2.0*_gc.GroupScheme_v[i]+1))/6.0 );

  }

  return total_vacancy*_gc._atomic_vol;
}

