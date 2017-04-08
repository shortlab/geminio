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

#include "AddGImmobile.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "GImmobileL0.h"
#include "GImmobileL1.h"

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
InputParameters validParams<AddGImmobile>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<int>("max_mobile_v", "maximum size of mobile vacancy cluster");
  params.addRequiredParam<int>("max_mobile_i", "maximum size of mobile intersitial cluster");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  return params;
}


AddGImmobile::AddGImmobile(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddGImmobile::act()
{
  int number_v = getParam<int>("number_v");
  int number_i = getParam<int>("number_i");
  int num_mobile_v = getParam<int>("max_mobile_v");
  int num_mobile_i = getParam<int>("max_mobile_i");
  
  std::string uo = getParam<std::string>("group_constant");
  std::string _prefix = name();
  std::string var_name;

//first add immobile v
  for(int cur_size=num_mobile_v+1; cur_size<=number_v; cur_size++){
    std::vector<VariableName> coupled_v_vars;
    std::vector<VariableName> coupled_i_vars;

    for(int i=1;i<=num_mobile_v;i++){
      var_name = _prefix + "0v" + Moose::stringify(i);
      coupled_v_vars.push_back(var_name);
      var_name = _prefix + "1v" + Moose::stringify(i);
      coupled_v_vars.push_back(var_name);
    }
    for(int i=1;i<=num_mobile_i;i++){
      var_name = _prefix + "0i" + Moose::stringify(i);
      coupled_i_vars.push_back(var_name);
      var_name = _prefix + "1i" + Moose::stringify(i);
      coupled_i_vars.push_back(var_name);
    }

    for(int i=num_mobile_v;i>=1;i--){//for vv reaction gain(+)
      var_name = _prefix + "0v" + Moose::stringify(cur_size-i);
      coupled_v_vars.push_back(var_name);
      var_name = _prefix + "1v" + Moose::stringify(cur_size-i);
      coupled_v_vars.push_back(var_name);
    }

    var_name = _prefix + "0v" + Moose::stringify(cur_size);
    coupled_v_vars.push_back(var_name);
    var_name = _prefix + "1v" + Moose::stringify(cur_size);
    coupled_v_vars.push_back(var_name);

    for(int i=1;i<=std::min(num_mobile_i,number_v-cur_size);i++){//for vi reaction gain(+)
      var_name = _prefix + "0v" + Moose::stringify(cur_size+i);
      coupled_v_vars.push_back(var_name);
      var_name = _prefix + "1v" + Moose::stringify(cur_size+i);
      coupled_v_vars.push_back(var_name);
    } 
    if(num_mobile_i==0 && cur_size!=number_v){
      var_name = _prefix + "0v" + Moose::stringify(cur_size+1);
      coupled_v_vars.push_back(var_name);
      var_name = _prefix + "1v" + Moose::stringify(cur_size+1);
      coupled_v_vars.push_back(var_name);
    } 

    var_name = name() +"0v"+ Moose::stringify(cur_size);
    InputParameters params = _factory.getValidParams("GImmobileL0");
    params.set<NonlinearVariableName>("variable") = var_name;
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<UserObjectName>("user_object") = uo;
    params.set<int>("number_v") = number_v;
    params.set<int>("number_i") = number_i;
    params.set<int>("max_mobile_v") = num_mobile_v;
    params.set<int>("max_mobile_i") = num_mobile_i;
    _problem->addKernel("GImmobileL0", "GImmobileL0_" + var_name+ "_" + Moose::stringify(counter), params);
    //printf("add GImmobileL0: %s \n",var_name.c_str());
    counter++;

    var_name = name() +"1v"+ Moose::stringify(cur_size);
    InputParameters params1 = _factory.getValidParams("GImmobileL1");
    params1.set<NonlinearVariableName>("variable") = var_name;
    params1.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params1.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params1.set<UserObjectName>("user_object") = uo;
    params1.set<int>("number_v") = number_v;
    params1.set<int>("number_i") = number_i;
    params1.set<int>("max_mobile_v") = num_mobile_v;
    params1.set<int>("max_mobile_i") = num_mobile_i;
    _problem->addKernel("GImmobileL1", "GImmobileL1_" + var_name+ "_" + Moose::stringify(counter), params1);
    //printf("add GImmobileL1: %s \n",var_name.c_str());
    counter++;

    //for(int i=0;i<coupled_v_vars.size();i++) printf("coupled v var: %s\n",coupled_v_vars[i].c_str());
    //for(int i=0;i<coupled_i_vars.size();i++) printf("coupled i var: %s\n",coupled_i_vars[i].c_str());
  }
      
//Second add immobile i
  for(int cur_size=num_mobile_i+1; cur_size<=number_i; cur_size++){

    std::vector<VariableName> coupled_v_vars;
    std::vector<VariableName> coupled_i_vars;
    for(int i=1;i<=num_mobile_v;i++){
      var_name = _prefix + "0v" + Moose::stringify(i);
      coupled_v_vars.push_back(var_name);
      var_name = _prefix + "1v" + Moose::stringify(i);
      coupled_v_vars.push_back(var_name);
    }
    for(int i=1;i<=num_mobile_i;i++){
      var_name = _prefix + "0i" + Moose::stringify(i);
      coupled_i_vars.push_back(var_name);
      var_name = _prefix + "1i" + Moose::stringify(i);
      coupled_i_vars.push_back(var_name);
    }
    for(int i=num_mobile_i;i>=1;i--){//for ii reaction gain(+)
      var_name = _prefix + "0i" + Moose::stringify(cur_size-i);
      coupled_i_vars.push_back(var_name);
      var_name = _prefix + "1i" + Moose::stringify(cur_size-i);
      coupled_i_vars.push_back(var_name);
    }

    var_name = _prefix + "0i" + Moose::stringify(cur_size);
    coupled_i_vars.push_back(var_name);
    var_name = _prefix + "1i" + Moose::stringify(cur_size);
    coupled_i_vars.push_back(var_name);

    for(int i=1;i<=std::min(num_mobile_v,number_i-cur_size);i++){//for iv reaction gain(+)
      var_name = _prefix + "0i" + Moose::stringify(cur_size+i);
      coupled_i_vars.push_back(var_name);
      var_name = _prefix + "1i" + Moose::stringify(cur_size+i);
      coupled_i_vars.push_back(var_name);
    } 
    if(num_mobile_v==0 && cur_size!=number_i){
      var_name = _prefix + "0i" + Moose::stringify(cur_size+1);
      coupled_i_vars.push_back(var_name);
      var_name = _prefix + "1i" + Moose::stringify(cur_size+1);
      coupled_i_vars.push_back(var_name);
    } 

    std::string var_name_i = name() +"0i"+ Moose::stringify(cur_size);
    InputParameters params = _factory.getValidParams("GImmobileL0");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<UserObjectName>("user_object") = uo;
    params.set<int>("number_v") = number_v;
    params.set<int>("number_i") = number_i;
    params.set<int>("max_mobile_v") = num_mobile_v;
    params.set<int>("max_mobile_i") = num_mobile_i;
    _problem->addKernel("GImmobileL0", "GImmobileL0_" + var_name_i+ "_" + Moose::stringify(counter), params);
    //printf("add GImmobileL0: %s \n",var_name_i.c_str());
    counter++;

    var_name_i = name() +"1i"+ Moose::stringify(cur_size);
    InputParameters params1 = _factory.getValidParams("GImmobileL1");
    params1.set<NonlinearVariableName>("variable") = var_name_i;
    params1.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params1.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params1.set<UserObjectName>("user_object") = uo;
    params1.set<int>("number_v") = number_v;
    params1.set<int>("number_i") = number_i;
    params1.set<int>("max_mobile_v") = num_mobile_v;
    params1.set<int>("max_mobile_i") = num_mobile_i;
    _problem->addKernel("GImmobileL1", "GImmobileL1_" + var_name_i+ "_" + Moose::stringify(counter), params1);
    //printf("add GImmobileL1: %s \n",var_name_i.c_str());
    counter++;

    //for(int i=0;i<coupled_v_vars.size();i++) printf("coupled v var: %s\n",coupled_v_vars[i].c_str());
    //for(int i=0;i<coupled_i_vars.size();i++) printf("coupled i var: %s\n",coupled_i_vars[i].c_str());

  }
}
