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

#include "UserObjectSingleVariable.h"

template<>
InputParameters validParams<UserObjectSingleVariable>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<UserObjectName>("user_object","the name of user object providing group constant");
  params.addRequiredParam<std::string>("call_function","map to function in the userobject");
  params.addCoupledVar("secondVar","the second variable which will probably be used in kernel");
  params.addParam<Real>("coeff",1,"the coefficent before the product of these two variables");
  
  return params;
}

UserObjectSingleVariable::UserObjectSingleVariable(const
     InputParameters & parameters)
     :Kernel(parameters),
     _gc(getUserObject<GroupConstant>("user_object")),
     _call_fun(getParam<std::string>("call_function")),
     _coeff(getParam<Real>("coeff"))
{
  int n = coupledComponents("secondVar");
  _vars.resize(n);
  _v_vals.resize(n);

  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  std::vector<VariableName> vars_names = getParam<std::vector<VariableName> >("secondVar");
  if(n==1)
  {
    groupNo = getGroupNumber(vars_names[0]);
    for (unsigned int i=0; i < _v_vals.size(); ++i)
    {
      _vars[i] = coupled("secondVar", i);
      _v_vals[i] = &coupledValue("secondVar", i);
    }
//  printf("coupled var num : %d\n",_v_vals.size());
  }
  else
    groupNo = getGroupNumber(cur_var_name);
  
  if(_call_fun.compare("emission") && _call_fun.compare("dislocation"))
    mooseError("Wrong call function name. Choices: emission, dislocation");
}

Real
UserObjectSingleVariable::computeQpResidual()
{
  Real gc;
  if(! _call_fun.compare("emission"))
      gc = _gc._emit(groupNo);
  else
      gc = _gc._disl(groupNo);

  if (!_v_vals.size())
    return _coeff * gc * _u[_qp]*_test[_i][_qp];
  return _coeff * gc * (*_v_vals[0])[_qp]*_test[_i][_qp];
}

Real
UserObjectSingleVariable::computeQpJacobian()
{
  Real gc;
  if(! _call_fun.compare("emission"))
      gc = _gc._emit(groupNo);
  else
      gc = _gc._disl(groupNo);
  if (!_v_vals.size())
    return _coeff * gc * _phi[_j][_qp]*_test[_i][_qp];
  return 0.0;
}

Real 
UserObjectSingleVariable::computeQpOffDiagJacobian(unsigned int jvar){
    Real gc;
    if(! _call_fun.compare("emission"))
        gc = _gc._emit(groupNo);
    else
        gc = _gc._disl(groupNo);

    if (_v_vals.size()==1 && jvar == _vars[0])
      return _coeff * gc * _phi[_j][_qp]*_test[_i][_qp];
    return 0.0;
}


int
UserObjectSingleVariable::getGroupNumber(std::string str)
{
  int len=str.length(),i=len;
  while(std::isdigit(str[i-1])) i--;
  int no = std::atoi((str.substr(i)).c_str());
  while(i>=0){
      i--;
      if(str[i]=='v'){no = no;break;}
      if(str[i]=='i'){no = -no;break;}
  }
  return no;
}
