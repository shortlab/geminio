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

#include "SingleVariable.h"

template<>
InputParameters validParams<SingleVariable>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("secondVar","the second variable which will probably be used in kernel");
  params.addParam<Real>("coeff",1,"the coefficent before the product of these two variables");
  
  return params;
}

SingleVariable::SingleVariable(const
     InputParameters & parameters)
     :Kernel(parameters),
     _coeff(getParam<Real>("coeff"))
{
  int n = coupledComponents("secondVar");
  _vars.resize(n);
  _v_vals.resize(n);

  for (unsigned int i=0; i < _v_vals.size(); ++i)
  {
    _vars[i] = coupled("secondVar", i);
    _v_vals[i] = &coupledValue("secondVar", i);
  }
//  printf("coupled var num : %d\n",_v_vals.size());
/*
  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  std::vector<VariableName> vars_names = getParam<std::vector<VariableName> >("secondVar");
  printf("SingleVariable primary var: %s\n",cur_var_name.c_str());
  printf("coupled vars: "); 
  for (unsigned int i=0; i < vars_names.size(); ++i)
  {
  printf(" %s ",vars_names[i].c_str());
  }
  printf("\n");
 */ 
}

Real
SingleVariable::computeQpResidual()
{
  if (!_v_vals.size())
    return _coeff *_u[_qp]*_test[_i][_qp];
  return _coeff * (*_v_vals[0])[_qp]*_test[_i][_qp];
}

Real
SingleVariable::computeQpJacobian()
{
  if (!_v_vals.size())
    return _coeff * _phi[_j][_qp]*_test[_i][_qp];
  return 0.0;
}

Real 
SingleVariable::computeQpOffDiagJacobian(unsigned int jvar){
    if (_v_vals.size()==1 && jvar == _vars[0])
      return _coeff * _phi[_j][_qp]*_test[_i][_qp];
    else
      return 0.0;
}


