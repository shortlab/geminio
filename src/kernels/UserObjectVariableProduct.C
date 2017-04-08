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

#include "UserObjectVariableProduct.h"

template<>
InputParameters validParams<UserObjectVariableProduct>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<UserObjectName>("user_object","the name of user object providing group constant");
  params.addCoupledVar("coupled_vars","the values of this variable will be used in kernel");
  params.addParam<Real>("coeff",1,"the coefficent before the product of these two variables");
  
  return params;
}

UserObjectVariableProduct::UserObjectVariableProduct(const
     InputParameters & parameters)
     :Kernel(parameters),
     _gc(getUserObject<GroupConstant>("user_object")),
     _coeff(getParam<Real>("coeff"))
{
  int ncoupled = coupledComponents("coupled_vars");

  int n = ncoupled;//only two variables, 
  _vars.resize(n);
  _v_vals.resize(n);

  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  std::vector<VariableName> vars_names = getParam<std::vector<VariableName> >("coupled_vars");
  switch(ncoupled){
  case 0:
  {
    groupa = getGroupNumber(cur_var_name.c_str()); 
    groupb = groupa;
    break;
  }
  case 1:
  {
    groupa = getGroupNumber(cur_var_name.c_str()); 
    groupb = getGroupNumber(vars_names[0].c_str());
    break;
  }
  case 2:
  {
    groupa = getGroupNumber(vars_names[0].c_str());
    groupb = getGroupNumber(vars_names[1].c_str());
    break;
  }
  default:
    mooseError("Wrong number of groups, should be 2");
  }
    
  for (unsigned int i=0; i < _v_vals.size(); ++i)
  {
    _vars[i] = coupled("coupled_vars", i);
    _v_vals[i] = &coupledValue("coupled_vars", i);
  }
}

Real
UserObjectVariableProduct::computeQpResidual()
{
  Real gc = _gc._absorb(groupa,groupb);
  Real res = 0.0;
//  printf("Absorption %d %d ---> %f\n",groupa,groupb,gc);
  switch(_v_vals.size()){
  case 0:
  {
     res =  _u[_qp]*_u[_qp];
     break;
  }
  case 1:
  {
    res = (*_v_vals[0])[_qp]*_u[_qp];
    break;
  }
  case 2:
  {
    res = (*_v_vals[1])[_qp]* (*_v_vals[0])[_qp];//second * third variable
    break;
  }
  default:
    mooseError("Wrong number of groups, should be 2");
  }
  //printf("Group: %d vs %d; coef: %f; gc: %f; res: %f\n",groupa,groupb,_coeff,gc,res);
  return _coeff * gc * res * _test[_i][_qp]; 
}

Real
UserObjectVariableProduct::computeQpJacobian()
{
  Real gc = _gc._absorb(groupa,groupb);
  switch(_v_vals.size()){
  case 0:
  {
    return _coeff * gc * 2.0 * _phi[_j][_qp] * _u[_qp] * _test[_i][_qp];
    break;
  }
  case 1:
  {
    return _coeff * gc * (*_v_vals[0])[_qp]*_phi[_j][_qp]*_test[_i][_qp];
    break;
  }
  case 2:
  {
    return 0.0;
    break;
  }
  default:
    mooseError("Wrong number of groups, should be 2");
  }
 
}

Real
UserObjectVariableProduct::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real gc = _gc._absorb(groupa,groupb);
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
  
  return _coeff * gc * off *  _test[_i][_qp] ; 

}


int
UserObjectVariableProduct::getGroupNumber(std::string str)
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
