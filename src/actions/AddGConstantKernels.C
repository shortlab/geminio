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

#include "AddGConstantKernels.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "ConstantKernel.h"

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
InputParameters validParams<AddGConstantKernels>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addParam<std::vector<int> >("source_v_size", "A vector of distribution of creation along ion range");
  params.addParam<std::vector<int> >("source_i_size", "A vector of distribution of creation along ion range");
  params.addParam<std::vector<Real> >("source_v_value", "production of i clusters for source_i_size species");
  params.addParam<std::vector<Real> >("source_i_value", "production of v clusters for 1 source_v_size species");
  params.addParam<Real>("scaling_factor",1.0,"scaling factor to source rate");
  params.addParam<int>("number_single_v",0,"largest cluster size using group size of 1");
  params.addParam<int>("number_single_i",0,"largest cluster size using group size of 1");
  params.addParam<Real>("tlimit","set lifetime for the kernel");
  return params;
}


AddGConstantKernels::AddGConstantKernels(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddGConstantKernels::act()
{
  std::vector<int> v_size = getParam<std::vector<int> >("source_v_size");
  std::vector<int> i_size = getParam<std::vector<int> >("source_i_size");
  std::vector<Real> vv = getParam<std::vector<Real> >("source_v_value");
  std::vector<Real> ii = getParam<std::vector<Real> >("source_i_value");
  int max_single_v = getParam<int>("number_single_v");
  int max_single_i = getParam<int>("number_single_i");
  Real scaling_factor = getParam<Real>("scaling_factor");

//ATTENTION: the emission of vacancy cluster emit an interstitial or interstitial cluster emit an vacancy is not considered
  if((v_size.size() && *std::max_element(v_size.begin(),v_size.end()) > max_single_v) || (i_size.size() && *std::max_element(i_size.begin(),i_size.end()) > max_single_i))
    mooseError("Make sure number_single is larger than the largest source size");
   

  for (unsigned int cur_num = 1; cur_num <= v_size.size(); cur_num++)
  {
    std::string var_name_v = name() +"0v"+ Moose::stringify(v_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("ConstantKernel");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<Real>("value") = vv[cur_num-1]*scaling_factor;//Should be the production term of current size, gain should be negative in the kernel
    if (isParamValid("tlimit"))
      params.set<Real>("tlimit") = getParam<Real>("tlimit");
    _problem->addKernel("ConstantKernel", "ConstantKernel_" +  var_name_v + Moose::stringify(counter), params);
    counter++;
  }
  for (unsigned int cur_num = 1; cur_num <= i_size.size(); cur_num++)
  {
    std::string var_name_i = name() +"0i"+ Moose::stringify(i_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("ConstantKernel");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<Real>("value") = ii[cur_num-1]*scaling_factor;// gain should be negative in the kernel
    if (isParamValid("tlimit"))
      params.set<Real>("tlimit") = getParam<Real>("tlimit");
    _problem->addKernel("ConstantKernel", "ConstantKernel_"+ var_name_i+ Moose::stringify(counter), params);
    counter++;
  }
}
