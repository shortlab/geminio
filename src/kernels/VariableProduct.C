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

#include "VariableProduct.h"

template<>
InputParameters validParams<VariableProduct>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("coupled_vars","a vector of coupled variables");
  params.addParam<Real>("coeff",1,"the coefficent before the product of these two variables");
  
  return params;
}

VariableProduct::VariableProduct(const
     InputParameters & parameters)
     :Kernel(parameters),
     _coeff(getParam<Real>("coeff"))
{
  int n = coupledComponents("coupled_vars");
  _vars.resize(n);
  _v_vals.resize(n);

  for (unsigned int i=0; i < _v_vals.size(); ++i)
  {
    _vars[i] = coupled("coupled_vars", i);
    _v_vals[i] = &coupledValue("coupled_vars", i);
  }

 // print var names, error-check
/*
  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  std::vector<VariableName> vars_names = getParam<std::vector<VariableName> >("coupled_vars");
  printf("VariableProduct primary var: %s\n",cur_var_name.c_str());
  printf("coupled vars: "); 
  for (unsigned int i=0; i < vars_names.size(); ++i)
  {
  printf(" %s ",vars_names[i].c_str());
  }
  printf("\n\n");
*/
}

Real
VariableProduct::computeQpResidual()
{
  Real res = 0.0;
  switch(_v_vals.size()){
  case 0:
  {
    res = _u[_qp]*_u[_qp];
    break;
  }
  case 1:
  {
    res = (*_v_vals[0])[_qp]*_u[_qp];
    break;
  }
  case 2:
  {
    res = (*_v_vals[0])[_qp]* (*_v_vals[1])[_qp];//second * third variable
    break;
  }
  default:
    mooseError("Error in VariableProduct kernel");
  }

  //printf("Coef: %f with res: %f\n",_coeff,res);
  return _coeff * res * _test[_i][_qp];
}

Real
VariableProduct::computeQpJacobian()
{
  Real jac=0.0;
  switch(_v_vals.size()){
    case 0:
    {
      jac = 2.0 * _phi[_j][_qp]*_u[_qp];
      break;
    }
    case 1:
    {
      jac = (*_v_vals[0])[_qp]*_phi[_j][_qp];
      break;
    }
    case 2:
    {
      jac = 0.0;//second * third variable
      break;
    }
    default:
      mooseError("Error in VariableProduct kernel");
  }

  //printf("Coef: %f with res: %f\n",_coeff,res);
  return _coeff * jac * _test[_i][_qp];
}

Real
VariableProduct::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int ncouples = _v_vals.size();
  Real off = 0.0;
  if(ncouples==1){
    if(jvar == _vars[0])
      off = _phi[_j][_qp] * _u[_qp];
  }
  else if(ncouples==2){
    if(jvar==_vars[0] && jvar==_vars[1])
      off = 2.0 * _phi[_j][_qp] * (*_v_vals[1])[_qp];
    else if(jvar==_vars[0])
      off = _phi[_j][_qp] * (*_v_vals[1])[_qp];
    else if(jvar==_vars[1])
      off = _phi[_j][_qp] * (*_v_vals[0])[_qp];
  }
  return _coeff * off * _test[_i][_qp];
}

