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

#include "AddClusterDensity.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "Conversion.h"
#include "ClusterDensity.h"
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
InputParameters validParams<AddClusterDensity>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<std::string>("aux_var","aux variable name to hold value");
  params.addRequiredParam<std::string>("var_prefix","The prefix string of variables (in front of number)");
  params.addRequiredParam<std::vector<int> >("size_range","The number of vacancy variables to solve");
  params.addParam<Real>("scale_factor", 1.0, "A scale factor to be applied to the variable");
  return params;
}


AddClusterDensity::AddClusterDensity(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddClusterDensity::act()
{

  std::string aux_var = getParam<std::string>("aux_var");
  std::string  _var_prefix = getParam<std::string>("var_prefix");
  std::vector<int>  _size_range = getParam<std::vector<int> >("size_range");
  Real _scale_factor = getParam<Real>("scale_factor");
  std::vector<VariableName> coupled_vars;

  std::string var_name;
  for(int i=_size_range[0];i<=_size_range[1];i++){
      std::string var_name = _var_prefix + Moose::stringify(i);
      coupled_vars.push_back(var_name);
  }

  InputParameters params = _factory.getValidParams("ClusterDensity");
  params.set<AuxVariableName>("variable") = aux_var;
  params.set<std::vector<VariableName> > ("coupled_vars") = coupled_vars;
  params.set<Real>("scaling_factor") = _scale_factor;
  _problem->addAuxKernel("ClusterDensity", "ClusterDensity_" + aux_var, params);
}
      
