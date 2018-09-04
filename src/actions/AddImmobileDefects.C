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

#include "AddImmobileDefects.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "ImmobileDefects.h"

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
InputParameters validParams<AddImmobileDefects>()
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
  return params;
}


AddImmobileDefects::AddImmobileDefects(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddImmobileDefects::act()
{
  int number_v = getParam<int>("number_v");
  int number_i = getParam<int>("number_i");
  if ((isParamValid("mobile_v_size")^isParamValid("max_mobile_v") != 1) || (isParamValid("mobile_v_size")^isParamValid("max_mobile_v") != 1)){
    mooseError("Check the setting of mobile species, only one definition is allowed");
  }
  std::vector<int> v_mobile,i_mobile;
  if (isParamValid("mobile_v_size"))
    v_mobile = getParam<std::vector<int> >("mobile_v_size");
  else {
    int max_mobile_v = getParam<int>("max_mobile_v");
    for(int i=0;i<max_mobile_v;i++)
      v_mobile.push_back(i+1);
  }
  if (isParamValid("mobile_i_size"))
    i_mobile = getParam<std::vector<int> >("mobile_i_size");
  else {
    int max_mobile_i = getParam<int>("max_mobile_i");
    for(int i=0;i<max_mobile_i;i++)
      i_mobile.push_back(i+1);
  }

  std::string uo = getParam<std::string>("group_constant");

  int v_size = v_mobile.size();
  int i_size = i_mobile.size();
  std::string _prefix = name();
  std::string var_name;

//first add immobile v
  for(int cur_size=v_size+1; cur_size<=number_v; cur_size++){
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_size);
    InputParameters params = _factory.getValidParams("ImmobileDefects");
    params.set<NonlinearVariableName>("variable") = var_name_v;

    std::vector<VariableName> coupled_v_vars;
    std::vector<VariableName> coupled_i_vars;
    for(int i=0;i<v_size;i++){
      var_name = _prefix + "v" + Moose::stringify(v_mobile[i]);
      coupled_v_vars.push_back(var_name);
    }
    for(int i=0;i<i_size;i++){
      var_name = _prefix + "i" + Moose::stringify(i_mobile[i]);
      coupled_i_vars.push_back(var_name);
    }
    for(int i=0;i<v_size;i++){//for vv reaction gain(+)
      var_name = _prefix + "v" + Moose::stringify(cur_size-v_mobile[i]);
      coupled_v_vars.push_back(var_name);
    }
    for(int i=0;i<std::min(i_size,number_v-cur_size);i++){//for vi reaction gain(+)
      var_name = _prefix + "v" + Moose::stringify(cur_size+i_mobile[i]);
      coupled_v_vars.push_back(var_name);
    } 
    if(i_size==0 && cur_size!=number_v){
      var_name = _prefix + "v" + Moose::stringify(cur_size+1);
      coupled_v_vars.push_back(var_name);
    } 

    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<UserObjectName>("user_object") = uo;
    params.set<int>("number_v") = number_v;
    params.set<int>("number_i") = number_i;
    params.set<std::vector<int> >("mobile_v_size") = v_mobile;
    params.set<std::vector<int> >("mobile_i_size") = i_mobile;
    _problem->addKernel("ImmobileDefects", "ImmobileDefects_" + var_name_v+ "_" + Moose::stringify(counter), params);
//    printf("add ImmobileDefects: %s \n",var_name_v.c_str());
    counter++;
  }
      
//Second add immobile i
  for(int cur_size=i_size+1; cur_size<=number_i; cur_size++){
    std::string var_name_i = name() +"i"+ Moose::stringify(cur_size);
    InputParameters params = _factory.getValidParams("ImmobileDefects");
    params.set<NonlinearVariableName>("variable") = var_name_i;


    std::vector<VariableName> coupled_v_vars;
    std::vector<VariableName> coupled_i_vars;
    for(int i=0;i<v_size;i++){
      var_name = _prefix + "v" + Moose::stringify(v_mobile[i]);
      coupled_v_vars.push_back(var_name);
    }
    for(int i=0;i<i_size;i++){
      var_name = _prefix + "i" + Moose::stringify(i_mobile[i]);
      coupled_i_vars.push_back(var_name);
    }
    for(int i=0;i<i_size;i++){
      var_name = _prefix + "i" + Moose::stringify(cur_size-i_mobile[i]);
      coupled_i_vars.push_back(var_name);
    }
    for(int i=0;i<std::min(v_size,number_i-cur_size);i++){
      var_name = _prefix + "i" + Moose::stringify(cur_size+v_mobile[i]);
      coupled_i_vars.push_back(var_name);
    } 
    if(v_size==0 && cur_size!=number_i){
      var_name = _prefix + "i" + Moose::stringify(cur_size+1);
      coupled_i_vars.push_back(var_name);
    } 

    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<UserObjectName>("user_object") = uo;
    params.set<int>("number_v") = number_v;
    params.set<int>("number_i") = number_i;
    params.set<std::vector<int> >("mobile_v_size") = v_mobile;
    params.set<std::vector<int> >("mobile_i_size") = i_mobile;
    _problem->addKernel("ImmobileDefects", "ImmobileDefects_" + var_name_i+ "_" + Moose::stringify(counter), params);
//    printf("add ImmobileDefects: %s \n",var_name_i.c_str());
    counter++;
  }
}
