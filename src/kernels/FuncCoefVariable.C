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
#include "FuncCoefVariable.h"

template<>
InputParameters validParams<FuncCoefVariable>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<FunctionName>("coef", "0.5*x+0.5*y", "The function for conductivity");
  return params;
}

FuncCoefVariable::FuncCoefVariable(const InputParameters & parameters) :
    Kernel(parameters),
    _function(getFunction("coef"))
{
}

Real
FuncCoefVariable::computeQpResidual()
{
  Real k = _function.value(_t, _qp);
  return k*_test[_i][_qp]*_u[_qp];
}

Real
FuncCoefVariable::computeQpJacobian()
{
  Real k = _function.value(_t, _qp);
  return k*_test[_i][_qp]*_phi[_j][_qp];
}
