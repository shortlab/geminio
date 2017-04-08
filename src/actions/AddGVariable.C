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

#include "AddGVariable.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "MooseError.h"
#include <sstream>
#include <stdexcept>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

// class static initialization
const Real AddGVariable::_abs_zero_tol = 1e-12;

template<>
InputParameters validParams<AddGVariable>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();

  params.addRequiredParam<int>("number_v", "Total number of vacancy group to add");
  params.addRequiredParam<int>("number_i", "Total number of interstitial group to add");

  params.addParam<std::vector<int> >("IC_v_size","vacancy species number with initial concentration not ZERO");
  params.addParam<std::vector<int> >("IC_i_size","interstitial species number with initial concentration not ZERO");
  params.addParam<std::vector<Real> >("IC_v", "initial value for vacancy cluster correpsonding to IC_v_size");
  params.addParam<std::vector<Real> >("IC_i", "initial value for interstitial cluster corresponding to IC_i_size");

  params.addParam<std::string>("bc_type","neumann", "dirichlet or neumann, depending on w/t spatical dependence");

  params.addParam<Real>("boundary_value", 0.0, "Specifies the initial condition for this variable");
 // params.addParam<std::vector<SubdomainName> >("block", "The block id where this variable lives");
 // params.addParam<bool>("eigen", false, "True to make this variable an eigen variable");
  return params;
}


AddGVariable::AddGVariable(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddGVariable::act()
{
  int number_v = getParam<int>("number_v");
  int number_i = getParam<int>("number_i");
  std::vector<int> vv = getParam<std::vector<int> >("IC_v_size");
  std::vector<int> ii = getParam<std::vector<int> >("IC_i_size");
  std::vector<Real> initial_v = getParam<std::vector<Real> >("IC_v");
  std::vector<Real> initial_i = getParam<std::vector<Real> >("IC_i");
  if (vv.size() != initial_v.size() || ii.size() != initial_i.size())
    mooseError("IC_v_size and IC_v should have same length, so are IC_i_size and IC_i., groupsize = 1 ");
  
  std::string _bc_type = getParam<std::string>("bc_type");

  if (_current_task == "add_variable")
  {
    std::string var_name;

    for (int cur_num = 1; cur_num <= number_v; cur_num++)
    {
      var_name = name() +"0v" + Moose::stringify(cur_num);
      addVariable(var_name);
      var_name = name() +"1v" + Moose::stringify(cur_num);
      addVariable(var_name);
    }
    for (int cur_num = 1; cur_num <= number_i; cur_num++)
    {
      var_name = name() +"0i" + Moose::stringify(cur_num);
      addVariable(var_name);
      var_name = name() +"1i" + Moose::stringify(cur_num);
      addVariable(var_name);
    }
  }

  else if(_current_task == "add_bc")
  {
    Real bc_val = getParam<Real>("boundary_value");
    std::string bc_name;
    if(_bc_type == "dirichlet") bc_name = "DirichletBC";
    else if(_bc_type == "neumann") bc_name = "NeumannBC";
    else 
        mooseError("This bc name: ", bc_name, " does not exist");

    std::string var_name;
    for (int cur_num = 1; cur_num <= number_v; cur_num++)
    {
      var_name = name() +"0v" + Moose::stringify(cur_num);
      InputParameters params = _factory.getValidParams(bc_name);
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params.set<Real>("value") = bc_val;
      _problem->addBoundaryCondition(bc_name, var_name + "_left", params);
      params.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params.set<Real>("value") = bc_val;
      _problem->addBoundaryCondition(bc_name, var_name + "_right", params);

      var_name = name() +"1v" + Moose::stringify(cur_num);
      InputParameters params1 = _factory.getValidParams(bc_name);
      params1.set<NonlinearVariableName>("variable") = var_name;
      params1.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params1.set<Real>("value") = bc_val;
      _problem->addBoundaryCondition(bc_name, var_name + "_left", params1);
      params1.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params1.set<Real>("value") = bc_val;
      _problem->addBoundaryCondition(bc_name, var_name + "_right", params1);
    }

    for (int cur_num = 1; cur_num <= number_i; cur_num++)
    {
      var_name = name() +"0i" + Moose::stringify(cur_num);
      InputParameters params = _factory.getValidParams(bc_name);
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params.set<Real>("value") = bc_val;
      _problem->addBoundaryCondition(bc_name, var_name + "_left", params);
      params.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params.set<Real>("value") = bc_val;
      _problem->addBoundaryCondition(bc_name, var_name + "_right", params);

      var_name = name() +"1i" + Moose::stringify(cur_num);
      InputParameters params1 = _factory.getValidParams(bc_name);
      params1.set<NonlinearVariableName>("variable") = var_name;
      params1.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params1.set<Real>("value") = bc_val;
      _problem->addBoundaryCondition(bc_name, var_name + "_left", params1);
      params1.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params1.set<Real>("value") = bc_val;
      _problem->addBoundaryCondition(bc_name, var_name + "_right", params1);
    }
  }

  else if (_current_task == "add_ic")
  {
    std::string var_name;
    for (int cur_num = 1; cur_num <= number_v; cur_num++)
    {
      var_name = name()+ "0v" + Moose::stringify(cur_num);
      InputParameters params = _factory.getValidParams("ConstantIC");
      params.set<VariableName>("variable") = var_name;
      std::vector<int>::iterator it=find(vv.begin(),vv.end(),cur_num);
      params.set<Real>("value") = (it==vv.end()? 0.0: initial_v[it-vv.begin()]);
      _problem->addInitialCondition("ConstantIC", "ConstantIC_"+var_name, params);

      var_name = name() +"1v" + Moose::stringify(cur_num);
      InputParameters params1 = _factory.getValidParams("ConstantIC");
      params1.set<VariableName>("variable") = var_name;
      params1.set<Real>("value") = 0.0;
      _problem->addInitialCondition("ConstantIC", "ConstantIC_"+var_name, params1);
    }

    for (int cur_num = 1; cur_num <= number_i; cur_num++)
    {
      var_name = name()+ "0i" + Moose::stringify(cur_num);
      InputParameters params = _factory.getValidParams("ConstantIC");
      params.set<VariableName>("variable") = var_name;
      std::vector<int>::iterator it=find(ii.begin(),ii.end(),cur_num);
      params.set<Real>("value") = (it==ii.end()? 0.0: initial_i[it-ii.begin()]);
      _problem->addInitialCondition("ConstantIC", "ConstantIC_"+var_name, params);

      var_name = name() +"1i" + Moose::stringify(cur_num);
      InputParameters params1 = _factory.getValidParams("ConstantIC");
      params1.set<VariableName>("variable") = var_name;
      params1.set<Real>("value") = 0.0;
      _problem->addInitialCondition("ConstantIC", "ConstantIC_"+var_name, params1);
    }
  }

}
