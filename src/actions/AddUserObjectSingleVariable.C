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

#include "AddUserObjectSingleVariable.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "UserObjectSingleVariable.h"

#include <sstream>
#include <stdexcept>
#include <algorithm>
// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"
static unsigned int counter = 0;

template<>
InputParameters validParams<AddUserObjectSingleVariable>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<int> >("mobile_v_size", "A vector of mobile species");
  params.addRequiredParam<std::vector<int> >("mobile_i_size", "A vector of mobile species");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  params.addRequiredParam<std::string>("call_function", "map to function to call in the userobject");
  params.addParam<Real>("temperature",600.0,"Temperature");
  return params;
}


AddUserObjectSingleVariable::AddUserObjectSingleVariable(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddUserObjectSingleVariable::act()
{
  unsigned int number_v = getParam<unsigned int>("number_v");
  unsigned int number_i = getParam<unsigned int>("number_i");
  std::vector<int> v_size = getParam<std::vector<int> >("mobile_v_size");
  std::vector<int> i_size = getParam<std::vector<int> >("mobile_i_size");
  Real temp = getParam<Real>("temperature");
  std::string uo = getParam<std::string>("group_constant");
  std::string call_function = getParam<std::string>("call_function");

  for (unsigned int cur_num = 2; cur_num <= number_v; cur_num++)
  //vv emission
  {
    double emit_coef = 0.0;
    emit_coef = 1.0;//cur_num to vacancy,"+"
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("UserObjectSingleVariable");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<Real>("coeff") = emit_coef;//loss "+"
    params.set<UserObjectName>("user_object") = uo;
    params.set<std::string>("call_function") = call_function;
    _problem->addKernel("UserObjectSingleVariable", "SingleV_" + var_name_v+ "_" + Moose::stringify(counter), params);
    //printf("add UserObjectSingleVariable: %s (emit -), coef: %lf\n",var_name_v.c_str(),emit_coef);
    counter++;

    //point defects (+)
    std::string var_name_v2 = name() +"v"+ Moose::stringify(1);
    InputParameters params0 = _factory.getValidParams("UserObjectSingleVariable");
    params0.set<NonlinearVariableName>("variable") = var_name_v2;
    params0.set<std::vector<VariableName> > ("secondVar").push_back(var_name_v);
    params0.set<UserObjectName>("user_object") = uo;
    params0.set<std::string>("call_function") = call_function;
    emit_coef = (-1.0)*1.0;
    params0.set<Real>("coeff") = emit_coef; //emission coefficent from larger size to current size; gain should be negative in kernel
    _problem->addKernel("UserObjectSingleVariable", "SingleV_" + var_name_v2 + "_" +  Moose::stringify(counter), params0);
    //printf("add UserObjectSingleVariable: %s (%s emit +) , coef: %lf\n",var_name_v2.c_str(),var_name_v.c_str(),emit_coef);
    counter++;
    //cur_num-1 (+) (account for 2->1+1)
    std::string var_name_other = name() + "v" + Moose::stringify(cur_num-1);
    InputParameters params_other = _factory.getValidParams("UserObjectSingleVariable");
    params_other.set<NonlinearVariableName>("variable") = var_name_other;
    params_other.set<std::vector<VariableName> > ("secondVar").push_back(var_name_v);
    params_other.set<UserObjectName>("user_object") = uo;
    params_other.set<std::string>("call_function") = call_function;
    emit_coef = (-1.0)*1.0;
    params_other.set<Real>("coeff") = emit_coef;//emission coefficent from larger size to current size; gain should be negative in kernel
    _problem->addKernel("UserObjectSingleVariable", "SingleV_" + var_name_other + "_" +  Moose::stringify(counter), params_other);
    //printf("add UserObjectSingleVariable: %s (%s emit +) , coef: %lf\n",var_name_other.c_str(),var_name_v.c_str(),emit_coef);
    counter++;
  }

  for (unsigned int cur_num = 2; cur_num <= number_i; cur_num++)
  //ii emission
  {
    double emit_coef = 0.0;
    emit_coef = 1.0;//cur_num to vacancy,"+"
    std::string var_name_i = name() +"i"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("UserObjectSingleVariable");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<Real>("coeff") = emit_coef;//loss "+"
    params.set<UserObjectName>("user_object") = uo;
    params.set<std::string>("call_function") = call_function;
    _problem->addKernel("UserObjectSingleVariable", "SingleV_" + var_name_i+ "_" + Moose::stringify(counter), params);
    //printf("add UserObjectSingleVariable: %s (emit -), coef: %lf\n",var_name_i.c_str(),emit_coef);
    counter++;

    //point defects (+)
    std::string var_name_i2 = name() +"i"+ Moose::stringify(1);
    InputParameters params0 = _factory.getValidParams("UserObjectSingleVariable");
    params0.set<NonlinearVariableName>("variable") = var_name_i2;
    params0.set<std::vector<VariableName> > ("secondVar").push_back(var_name_i);
    params0.set<UserObjectName>("user_object") = uo;
    params0.set<std::string>("call_function") = call_function;
    emit_coef = -1.0*1.0;
    params0.set<Real>("coeff") = emit_coef;//emission coefficent from larger size to current size; gain should be negative in kernel
    _problem->addKernel("UserObjectSingleVariable", "SingleV_" + var_name_i2 + "_" +  Moose::stringify(counter), params0);
    //printf("add UserObjectSingleVariable: %s (%s emit +) , coef: %lf\n",var_name_i2.c_str(),var_name_i.c_str(),emit_coef);
    counter++;
    //cur_num-1 (+) (account for 2->1+1)
    std::string var_name_other = name() + "i" + Moose::stringify(cur_num-1);
    InputParameters params_other = _factory.getValidParams("UserObjectSingleVariable");
    params_other.set<NonlinearVariableName>("variable") = var_name_other;
    params_other.set<std::vector<VariableName> > ("secondVar").push_back(var_name_i);
    params_other.set<UserObjectName>("user_object") = uo;
    params_other.set<std::string>("call_function") = call_function;
    emit_coef = (-1.0)*1.0;
    params_other.set<Real>("coeff") = emit_coef;//emission coefficent from larger size to current size; gain should be negative in kernel
    _problem->addKernel("UserObjectSingleVariable", "SingleV_" + var_name_other + "_" +  Moose::stringify(counter), params_other);
    //printf("add UserObjectSingleVariable: %s (%s emit +) , coef: %lf\n",var_name_other.c_str(),var_name_i.c_str(),emit_coef);
    counter++;
  }
}
