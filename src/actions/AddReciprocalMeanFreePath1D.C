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

#include "AddReciprocalMeanFreePath1D.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "Conversion.h"
#include "ReciprocalMeanFreePath1D.h"
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
InputParameters validParams<AddReciprocalMeanFreePath1D>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addParam<std::string>("aux_prefix","Rlambda","aux variable name prefix");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  params.addRequiredParam<int>("max_mobile_i", "maximum size of mobile intersitial cluster");
  return params;
}


AddReciprocalMeanFreePath1D::AddReciprocalMeanFreePath1D(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddReciprocalMeanFreePath1D::act()
{
  int number_v = getParam<int>("number_v");
  int number_i = getParam<int>("number_i");
  int max_mobile_i = getParam<int>("max_mobile_i");

  std::string aux_prefix = getParam<std::string>("aux_prefix");
  std::string uo = getParam<std::string>("group_constant");

  std::vector<VariableName> coupled_v_vars;
  std::vector<VariableName> coupled_i_vars;

  std::string var_name;
  for (int cur_num = 1; cur_num <= number_v; cur_num++)
  {
    var_name = name() +"0v" + Moose::stringify(cur_num);
    coupled_v_vars.push_back(var_name);
  }
  for (int cur_num = 1; cur_num <= number_i; cur_num++)
  {
    var_name = name() +"0i" + Moose::stringify(cur_num);
    coupled_i_vars.push_back(var_name);
  }


  for (int cur_num =1; cur_num <= max_mobile_i; cur_num++){
    InputParameters params = _factory.getValidParams("ReciprocalMeanFreePath1D");
    params.set<AuxVariableName>("variable") = aux_prefix + Moose::stringify(cur_num);
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<int>("mobile_SIA_size") = cur_num;
    params.set<UserObjectName>("user_object") = uo;
    _problem->addAuxKernel("ReciprocalMeanFreePath1D", "ReciprocalMeanFreePath1D_" + aux_prefix + Moose::stringify(cur_num), params);
  }

}
      
