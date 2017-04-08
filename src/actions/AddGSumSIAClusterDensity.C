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

#include "AddGSumSIAClusterDensity.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "Conversion.h"
#include "GSumSIAClusterDensity.h"
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
InputParameters validParams<AddGSumSIAClusterDensity>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::string>("aux_var","aux variable name to hold value");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  params.addParam<Real>("scale_factor",1.0,"scale factor used in the kernel");
  params.addParam<int>("lower_bound","starting size to count, inclusive");
  params.addParam<int>("upper_bound","ending size to count, inclusive");
  return params;
}


AddGSumSIAClusterDensity::AddGSumSIAClusterDensity(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddGSumSIAClusterDensity::act()
{
  int number_i = getParam<int>("number_i");
  Real scale_factor = getParam<Real>("scale_factor");

  std::string aux_var = getParam<std::string>("aux_var");
  std::string uo = getParam<std::string>("group_constant");

  std::vector<VariableName> coupled_i_vars;

  std::string var_name;
  for (int cur_num = 1; cur_num <= number_i; cur_num++)
  {
    var_name = name() +"0i" + Moose::stringify(cur_num);
    coupled_i_vars.push_back(var_name);
    var_name = name() +"1i" + Moose::stringify(cur_num);
    coupled_i_vars.push_back(var_name);
  }

  InputParameters params = _factory.getValidParams("GSumSIAClusterDensity");
  params.set<AuxVariableName>("variable") = aux_var;
  params.set<std::vector<VariableName> > ("coupled_vars") = coupled_i_vars;
  params.set<Real>("scale_factor") = scale_factor;
  if (isParamValid("lower_bound"))
    params.set<int>("lower_bound") = getParam<int>("lower_bound");
  if (isParamValid("upper_bound"))
    params.set<int>("upper_bound") = getParam<int>("upper_bound");
  params.set<UserObjectName>("user_object") = uo;
  _problem->addAuxKernel("GSumSIAClusterDensity", "GSumSIAClusterDensity_" + aux_var, params);
}
      
