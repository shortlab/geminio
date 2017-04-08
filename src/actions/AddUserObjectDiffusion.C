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

#include "AddUserObjectDiffusion.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "UserObjectDiffusion.h"

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
static unsigned int counter = 0;

template<>
InputParameters validParams<AddUserObjectDiffusion>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<std::vector<int> >("mobile_v_size", "A vector of mobile species");
  params.addRequiredParam<std::vector<int> >("mobile_i_size", "A vector of mobile species");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  return params;
}


AddUserObjectDiffusion::AddUserObjectDiffusion(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddUserObjectDiffusion::act()
{
  std::vector<int> v_size = getParam<std::vector<int> >("mobile_v_size");
  std::vector<int> i_size = getParam<std::vector<int> >("mobile_i_size");
  std::string uo = getParam<std::string>("group_constant");
  std::vector<Real> vv,ii;
  int total_mobile_v = v_size.size();
  int total_mobile_i = i_size.size();
  for(int i=0;i<total_mobile_v;i++){
      vv.push_back(1.0);//account for sign; sink "+"
  }
  for(int i=0;i<total_mobile_i;i++){
      ii.push_back(1.0);//account for sign; sink "+"
  }

  for(int cur_num=1;cur_num<total_mobile_v;cur_num++){
    double coef = vv[cur_num-1];//cur_num to vacancy,"+"
    std::string var_name_v = name() +"v"+ Moose::stringify(v_size[cur_num-1]);//add kernel for mobile v
    InputParameters params = _factory.getValidParams("UserObjectDiffusion");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<Real>("coeff") = coef;//loss "+"
    params.set<UserObjectName>("user_object") = uo;
    _problem->addKernel("UserObjectDiffusion", "Diffusion_" + var_name_v+Moose::stringify(counter), params);
    printf("add UserObjectDiffusion disl: %s, coef: %lf\n",var_name_v.c_str(),coef);
    counter++;
  }

  for(int cur_num=1;cur_num<total_mobile_i;cur_num++){
    double coef = ii[cur_num-1];
    std::string var_name_i = name() +"i"+ Moose::stringify(i_size[cur_num-1]);//add kernel for mobile i
    InputParameters params = _factory.getValidParams("UserObjectDiffusion");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<Real>("coeff") = coef;//loss "+"
    params.set<UserObjectName>("user_object") = uo;
    _problem->addKernel("UserObjectDiffusion", "Diffusion_" + var_name_i+Moose::stringify(counter), params);
    printf("add UserObjectDiffusion disl: %s, coef: %lf\n",var_name_i.c_str(),coef);
    counter++;
  }
}
