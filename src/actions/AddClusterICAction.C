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

#include "AddClusterICAction.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddICAction.h"
#include "Conversion.h"
#include "ConstantIC.h"

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
static int counter=0;

template<>
InputParameters validParams<AddClusterICAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<unsigned int> >("IC_v_size","vacancy species number with initial concentration not ZERO");
  params.addRequiredParam<std::vector<unsigned int> >("IC_i_size","interstitial species number with initial concentration not ZERO");
  params.addRequiredParam<std::vector<Real> >("IC_v", "initial value for vacancy cluster correpsonding to IC_v_size");
  params.addRequiredParam<std::vector<Real> >("IC_i", "initial value for interstitial cluster corresponding to IC_i_size");
  return params;
}


AddClusterICAction::AddClusterICAction(const InputParameters & params) :
    Action(params)
{
}

void
AddClusterICAction::act()
{
  unsigned int number_v = getParam<unsigned int>("number_v");
  unsigned int number_i = getParam<unsigned int>("number_i");
  std::vector<unsigned int> vv = getParam<std::vector<unsigned int> >("IC_v_size");
  std::vector<unsigned int> ii = getParam<std::vector<unsigned int> >("IC_i_size");
  std::vector<Real> initial_v = getParam<std::vector<Real> >("IC_v");
  std::vector<Real> initial_i = getParam<std::vector<Real> >("IC_i");
  if (vv.size() != initial_v.size() || ii.size() != initial_i.size())
    mooseError("IC_v_size and IC_v should have same length, so are IC_i_size and IC_i.");
  
  unsigned int max_ic_v = (vv.size()>0)?(*std::max_element(vv.begin(),vv.end())):0; 
  unsigned int max_ic_i = (ii.size()>0)?(*std::max_element(ii.begin(),ii.end())):0; 
//initial values for v and i from IC_v_size and IC_i_size

  for (unsigned int cur_num = 1; cur_num <= max_ic_v; cur_num++)
  {
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = var_name_v;
    std::vector<unsigned int>::iterator it=find(vv.begin(),vv.end(),cur_num);
    params.set<Real>("value") = (it==vv.end()? 0.0: initial_v[it-vv.begin()]);
    //std::cout << "initial: " << (it==vv.end()? 0.0: initial_v[it-vv.begin()]) << std::endl;
    _problem->addInitialCondition("ConstantIC", "ConstantIC_"+var_name_v, params);
    counter++;
  }
  for (unsigned int cur_num = max_ic_v+1; cur_num <= number_v; cur_num++)
  {
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = var_name_v;
    params.set<Real>("value") = 0.0;
    _problem->addInitialCondition("ConstantIC", "ConstantIC_"+var_name_v, params);
    counter++;
  }


  for (unsigned int cur_num = 1; cur_num <= max_ic_i; cur_num++)
  {
    std::string var_name_i = name() +"i"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = var_name_i;
    std::vector<unsigned int>::iterator it=find(ii.begin(),ii.end(),cur_num);
    params.set<Real>("value") = (it==ii.end()? 0.0: initial_i[it-ii.begin()]);
    _problem->addInitialCondition("ConstantIC", "ConstantIC_"+var_name_i, params);
    counter++;
  }
  for (unsigned int cur_num = max_ic_i+1; cur_num <= number_i; cur_num++)
  {
    std::string var_name_i = name() +"i"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = var_name_i;
    params.set<Real>("value") = 0.0;
    _problem->addInitialCondition("ConstantIC", "ConstantIC_"+var_name_i, params);
    counter++;
  }
}
