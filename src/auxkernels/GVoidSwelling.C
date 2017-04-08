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
  params.addRequiredParam<UserObjectName>("user_object","The name of user object providing interaction constants");
  return params;
}


GVoidSwelling::GVoidSwelling(const
                                   InputParameters & parameters)
  :AuxKernel(parameters),
  _gc(getUserObject<GGroup>("user_object"))
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
  Real total_vacancy = 0.0;//total vacancy conentration
  for(int i=0;i<_gc.GroupScheme_v.size()-1;i++){
    for(int j=_gc.GroupScheme_v[i]+1;j<=_gc.GroupScheme_v[i+1];j++){//(a,b]
      total_vacancy += ((*_val_v_vars[2*i])[_qp]+(*_val_v_vars[2*i+1])[_qp]*(j-_gc.GroupScheme_v_avg[i]))*j;
    }
  }

  return total_vacancy*_gc._atomic_vol;
}

