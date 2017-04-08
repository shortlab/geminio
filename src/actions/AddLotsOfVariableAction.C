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

#include "AddLotsOfVariableAction.h"
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
const Real AddLotsOfVariableAction::_abs_zero_tol = 1e-12;

template<>
InputParameters validParams<AddLotsOfVariableAction>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addParam<MooseEnum>("family", families, "Specifies the family of FE shape functions to use for this variable");
  params.addParam<MooseEnum>("order", orders,  "Specifies the order of the FE shape function to use for this variable");
  params.addParam<Real>("initial_condition", 0.0, "Specifies the initial condition for this variable");
  params.addParam<std::string>("bc_type","neumann", "dirichlet or neumann, depending on w/t spatical dependence");
  params.addParam<Real>("boundary_value", 0.0, "Specifies the initial condition for this variable");
  params.addParam<Real>("scaling", 1.0, "Specifies a scaling factor to apply to this variable");
  params.addParam<bool>("use_constIC",false,"Specifies whether use sinlge value for all variables as initial condition");
 // params.addParam<std::vector<SubdomainName> >("block", "The block id where this variable lives");
 // params.addParam<bool>("eigen", false, "True to make this variable an eigen variable");
  return params;
}


AddLotsOfVariableAction::AddLotsOfVariableAction(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddLotsOfVariableAction::act()
{
  int number_v = getParam<int>("number_v");
  int number_i = getParam<int>("number_i");
  bool use_constIC = getParam<bool>("use_constIC");
  std::string _bc_type = getParam<std::string>("bc_type");

  if (_current_task == "add_variable")
  {
    for (int cur_num = 1; cur_num <= number_v; cur_num++)
    {
      std::string var_name_v = name()+ "v" + Moose::stringify(cur_num);
      addVariable(var_name_v);
    }
    for (int cur_num = 1; cur_num <= number_i; cur_num++)
    {
      std::string var_name_i = name() +"i" + Moose::stringify(cur_num);
      addVariable(var_name_i);
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
    for (int cur_num = 1; cur_num <= (number_v+number_i); cur_num++)
    {
      std::string var_name;
      if (cur_num <= number_v)
        var_name = name() +"v"+ Moose::stringify(cur_num);
      else
        var_name = name() +"i"+ Moose::stringify(cur_num-number_v);

      InputParameters params = _factory.getValidParams(bc_name);
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params.set<Real>("value") = bc_val;

      _problem->addBoundaryCondition(bc_name, var_name + "_left", params);

      params.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params.set<Real>("value") = bc_val;

      _problem->addBoundaryCondition(bc_name, var_name + "_right", params);
    }
  }

  else if ( use_constIC == true && _current_task == "add_ic")
  {
    Real initial = getParam<Real>("initial_condition");
    if (initial > _abs_zero_tol || initial < -_abs_zero_tol)
    {
        for (int cur_num = 1; cur_num <= (number_v+number_i); cur_num++)
        {
          std::string var_name;
          if (cur_num <= number_v)
            var_name = name() +"v"+ Moose::stringify(cur_num);
          else
            var_name = name() +"i"+ Moose::stringify(cur_num-number_v);
          if(!_scalar_var){//for non-scalar variable, scalasr variable should use ScalarConstantIC
              InputParameters params = _factory.getValidParams("ConstantIC");
              params.set<VariableName>("variable") = var_name;
              params.set<Real>("value") = initial;
              _problem->addBoundaryCondition("ConstantIC", var_name + "_ic", params);
          }
        }  
    }  
  }

}
