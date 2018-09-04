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

#include "AddLotsOfCoeffDiffusion.h"
#include "MaterialParameters.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "CoeffDiffusion.h"

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
InputParameters validParams<AddLotsOfCoeffDiffusion>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addParam<std::vector<int> >("mobile_v_size", "A vector of mobile species sizes");
  params.addParam<std::vector<int> >("mobile_i_size", "A vector of mobile species sizes");
  params.addParam<int>("max_mobile_v", "maximum size of mobile vacancy cluster");
  params.addParam<int>("max_mobile_i", "maximum size of mobile intersitial cluster");
  params.addParam<bool>("custom_input",false,"mannually input coefficients");
  params.addParam<std::vector<Real> >("diff_i", "diffusivity of i clusters from 1 to n_i");
  params.addParam<std::vector<Real> >("diff_v", "diffusivity of v clusters from 1 to n_v");
  params.addParam<Real>("temperature",800,"system temperature [K]");
  return params;
}

AddLotsOfCoeffDiffusion::AddLotsOfCoeffDiffusion(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddLotsOfCoeffDiffusion::act()
{
  //unsigned int number_v = getParam<unsigned int>("number_v");
  //unsigned int number_i = getParam<unsigned int>("number_i");
  bool custom = getParam<bool>("custom_input");
  unsigned int _total_v = getParam<unsigned int>("number_v");
  unsigned int _total_i = getParam<unsigned int>("number_i");
  
  if ((isParamValid("mobile_v_size")^isParamValid("max_mobile_v") != 1) || (isParamValid("mobile_v_size")^isParamValid("max_mobile_v") != 1)){
    mooseError("Check the setting of mobile species, only one definition is allowed");
  }
  std::vector<int> v_size,i_size;
  if (isParamValid("mobile_v_size"))
    v_size = getParam<std::vector<int> >("mobile_v_size");
  else {
    int max_mobile_v = getParam<int>("max_mobile_v");
    for(int i=0;i<max_mobile_v;i++)
      v_size.push_back(i+1);
  }
  if (isParamValid("mobile_i_size"))
    i_size = getParam<std::vector<int> >("mobile_i_size");
  else {
    int max_mobile_i = getParam<int>("max_mobile_i");
    for(int i=0;i<max_mobile_i;i++)
      i_size.push_back(i+1);
  }

  std::vector<Real> vv,ii;
  Real temp = getParam<Real>("temperature");
  if (isParamValid("diff_v") && custom == true)
    vv = getParam<std::vector<Real> >("diff_v");
  else {
    for(int i=0;i<v_size.size();i++){
      vv.push_back(diff(v_size[i],"V",temp));}
  }
  if (isParamValid("diff_i") && custom == true)
    ii = getParam<std::vector<Real> >("diff_i");
  else {
    for(int i=0;i<i_size.size();i++){
      ii.push_back(diff(i_size[i],"I",temp));}
  }
  for (unsigned int cur_num = 1; cur_num <= v_size.size(); cur_num++)
  {
    if(cur_num <= _total_v){
      std::string var_name_v = name() +"v"+ Moose::stringify(v_size[cur_num-1]);
      InputParameters params0 = _factory.getValidParams("CoeffDiffusion");
      params0.set<NonlinearVariableName>("variable") = var_name_v;
      params0.set<Real>("diffusivity") =  vv[cur_num-1];
      _problem->addKernel("CoeffDiffusion", "diffusion_"+var_name_v+"_"+Moose::stringify(counter), params0);
      counter++;
    }
  }
  for (unsigned int cur_num = 1; cur_num <= i_size.size(); cur_num++)
  {
    if(cur_num <= _total_i){
      std::string var_name_i = name() +"i"+ Moose::stringify(i_size[cur_num-1]);
      InputParameters params0 = _factory.getValidParams("CoeffDiffusion");
      params0.set<NonlinearVariableName>("variable") = var_name_i;
      params0.set<Real>("diffusivity") =  vv[cur_num-1];
      _problem->addKernel("CoeffDiffusion", "diffusion_"+var_name_i+"_"+Moose::stringify(counter), params0);
      counter++;
    }
  }
}
