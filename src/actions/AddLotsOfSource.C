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

#include "AddLotsOfSource.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "BodyForce.h"

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
static unsigned int counter = 0;

template<>
InputParameters validParams<AddLotsOfSource>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<Real> >("source_v_size", "A vector of distribution of creation along ion range");
  params.addRequiredParam<std::vector<Real> >("source_i_size", "A vector of distribution of creation along ion range");
  params.addParam<std::string>("func_pre_name", "the sub-block name of [functions]");
  params.addParam<bool>("custom_input",false,"true: use mannully input coefficients");
  params.addParam<std::vector<Real> >("source_i", "production of i clusters for source_i_size species");
  params.addParam<std::vector<Real> >("source_v", "production of v clusters for 1 source_v_size species"); //!!!if provided source_i or source_v (constant for each size), the function wouldn't be necessary. In other words, provide either function or source_i and source_v
  return params;
}


AddLotsOfSource::AddLotsOfSource(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddLotsOfSource::act()
{
  std::vector<Real> v_size = getParam<std::vector<Real> >("source_v_size");
  std::vector<Real> i_size = getParam<std::vector<Real> >("source_i_size");
  std::vector<Real> vv = getParam<std::vector<Real> >("source_v");
  std::vector<Real> ii = getParam<std::vector<Real> >("source_i");
  bool custom = getParam<bool>("custom_input");
  unsigned int _total_v = getParam<unsigned int>("number_v");
  unsigned int _total_i = getParam<unsigned int>("number_i");
  std::string func_pre_name = (isParamValid("func_pre_name") ? getParam<std::string>("func_pre_name") : "");

//ATTENTION: the emission of vacancy cluster emit an interstitial or interstitial cluster emit an vacancy is not considered

  for (unsigned int cur_num = 1; cur_num <= v_size.size(); cur_num++)
  {
    if(cur_num <= _total_v){
    std::string var_name_v = name() +"v"+ Moose::stringify(v_size[cur_num-1]);
    std::string fun_name_v = func_pre_name + "v" + Moose::stringify(v_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("BodyForce");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    if (!custom){
      params.set<FunctionName>("function") = fun_name_v;
      params.set<Real>("value") = 1.0;
    }
    else if (vv.size() != v_size.size())
      mooseError("Vacancy source number should be the same with provided variables");
    else params.set<Real>("value") = vv[cur_num-1];//Should be the production term of current size, gain should be negative in the kernel
    _problem->addKernel("BodyForce", "bodyforce_" +  var_name_v + Moose::stringify(counter), params);
    printf("add Source: %s\n",var_name_v.c_str());
    counter++;
    }
  }
  for (unsigned int cur_num = 1; cur_num <= i_size.size(); cur_num++)
  {
    if(cur_num <= _total_i){
    std::string var_name_i = name() +"i"+ Moose::stringify(i_size[cur_num-1]);
    std::string fun_name_i = func_pre_name + "i" + Moose::stringify(i_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("BodyForce");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    if (!custom){
      params.set<FunctionName>("function") = fun_name_i;
      params.set<Real>("value") = 1.0;
    }
    else if (ii.size() != i_size.size())
      mooseError("Intersitial source number should be the same with provided variable");
    else params.set<Real>("value") = ii[cur_num-1];// gain should be negative in the kernel
    _problem->addKernel("BodyForce", "bodyforce_"+ var_name_i+ Moose::stringify(counter), params);
    printf("add Source: %s\n",var_name_i.c_str());
    counter++;
    }
  }
}
