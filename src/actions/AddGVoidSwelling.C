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

#include "AddGVoidSwelling.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "Conversion.h"
#include "GVoidSwelling.h"
#include "AddVariableAction.h"

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

template<>
InputParameters validParams<AddGVoidSwelling>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<std::string>("aux_var","aux variable name to hold value");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  return params;
}


AddGVoidSwelling::AddGVoidSwelling(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddGVoidSwelling::act()
{
  int number_v = getParam<int>("number_v");

  std::string aux_var = getParam<std::string>("aux_var");
  std::string uo = getParam<std::string>("group_constant");

  std::vector<VariableName> coupled_v_vars;
  std::vector<VariableName> coupled_i_vars;

  std::string var_name;
  for (int cur_num = 1; cur_num <= number_v; cur_num++)
  {
    var_name = name() +"0v" + Moose::stringify(cur_num);
    coupled_v_vars.push_back(var_name);
    var_name = name() +"1v" + Moose::stringify(cur_num);
    coupled_v_vars.push_back(var_name);
  }

  InputParameters params = _factory.getValidParams("GVoidSwelling");
  params.set<AuxVariableName>("variable") = aux_var;
  params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
  params.set<UserObjectName>("user_object") = uo;
  _problem->addAuxKernel("GVoidSwelling", "GVoidSwelling_" + aux_var, params);
}
      
