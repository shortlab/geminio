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

#include "AddMobileDefects.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "MobileDefects.h"

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
static int counter = 0;

template<>
InputParameters validParams<AddMobileDefects>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addParam<std::vector<int> >("mobile_v_size", "A vector of mobile vacancy cluster");
  params.addParam<std::vector<int> >("mobile_i_size", "A vector of mobile intersitial cluster");
  params.addParam<int>("max_mobile_v", "maximum size of mobile vacancy cluster");
  params.addParam<int>("max_mobile_i", "maximum size of mobile intersitial cluster");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  params.addParam<Real>("dislocation",0.0,"dislocation density");
  params.addParam<std::vector<Real> >("disl_mobile_v", "A vector of dislocation bias for mobile species");
  params.addParam<std::vector<Real> >("disl_mobile_i", "A vector of dislocation bias for mobile species");
  return params;
}


AddMobileDefects::AddMobileDefects(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddMobileDefects::act()
{
  int number_v = getParam<int>("number_v");
  int number_i = getParam<int>("number_i");

  if ((isParamValid("mobile_v_size")^isParamValid("max_mobile_v") != 1) || (isParamValid("mobile_v_size")^isParamValid("max_mobile_v") != 1)){
    mooseError("Check the setting of mobile species, only one definition is allowed");
  }
  std::vector<int> v_size,i_size;
  if (isParamValid("mobile_v_size"))
    v_size = getParam<std::vector<int> >("mobile_v_size");
  else {
    int max_mobile_v = getParam<int>("max_mobile_v");
    for(int i=0;i<max_mobile_v;i++)
      v_size.push_back(i+1);
  }
  if (isParamValid("mobile_i_size"))
    i_size = getParam<std::vector<int> >("mobile_i_size");
  else {
    int max_mobile_i = getParam<int>("max_mobile_i");
    for(int i=0;i<max_mobile_i;i++)
      i_size.push_back(i+1);
  }

  Real disl = getParam<Real>("dislocation");
  std::string uo = getParam<std::string>("group_constant");
  std::vector<Real> disl_v = getParam<std::vector<Real> >("disl_mobile_v");
  std::vector<Real> disl_i = getParam<std::vector<Real> >("disl_mobile_i");


  int num_mobile_v = v_size.size();
  int num_mobile_i = i_size.size();
  int disl_v_size = disl_v.size();
  int disl_i_size = disl_i.size();

  std::vector<VariableName> coupled_v_vars;
  std::vector<VariableName> coupled_i_vars;
  std::string _prefix = name();
  std::string var_name;
  for (int i=1; i <= number_v; ++i){
    var_name = _prefix + "v" + Moose::stringify(i);
    coupled_v_vars.push_back(var_name);
  }
  for (int i=1; i <= number_i; ++i){
    var_name = _prefix + "i" + Moose::stringify(i);
    coupled_i_vars.push_back(var_name);
  }
//first add mobile v
  for(int cur_num=1; cur_num<=num_mobile_v; cur_num++){
    double coef = 1.0;// (disl_v_size>0)?disl_v[cur_num-1]:1.0;//dislocation absorption factor
    std::string var_name_v = name() +"v"+ Moose::stringify(v_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("MobileDefects");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<Real>("dislocation") = disl;
    params.set<Real>("dislocation_factor") = coef;
    params.set<UserObjectName>("user_object") = uo;
    params.set<int>("number_v") = number_v;
    params.set<int>("number_i") = number_i;
    params.set<std::vector<int> >("mobile_v_size") = v_size;
    params.set<std::vector<int> >("mobile_i_size") = i_size;
    _problem->addKernel("MobileDefects", "MobileDefects_" + var_name_v+ "_" + Moose::stringify(counter), params);
//    printf("add MobileDefects: %s \n",var_name_v.c_str());
    counter++;
  }
      
//Second add mobile i
  for(int cur_num=1; cur_num<=num_mobile_i; cur_num++){
    double coef = 1.0;//(disl_i_size>0)?disl_i[cur_num-1]:1.0;//dislocation absorption factor
    std::string var_name_i = name() +"i"+ Moose::stringify(i_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("MobileDefects");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<Real>("dislocation") = disl;
    params.set<Real>("dislocation_factor") = coef;
    params.set<UserObjectName>("user_object") = uo;
    params.set<int>("number_v") = number_v;
    params.set<int>("number_i") = number_i;
    params.set<std::vector<int> >("mobile_v_size") = v_size;
    params.set<std::vector<int> >("mobile_i_size") = i_size;
    _problem->addKernel("MobileDefects", "MobileDefects_" + var_name_i+ "_" + Moose::stringify(counter), params);
//    printf("add MobileDefects: %s \n",var_name_i.c_str());
    counter++;
  }
}
