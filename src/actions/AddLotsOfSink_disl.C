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

#include "AddLotsOfSink_disl.h"
#include "MaterialParameters.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
//#include "DislocationSinkRate.h"

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
static unsigned int counter=0;
template<>
InputParameters validParams<AddLotsOfSink_disl>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<int> >("mobile_v_size", "A vector of mobile species");
  params.addRequiredParam<std::vector<int> >("mobile_i_size", "A vector of mobile species");
  params.addParam<std::string>("dislocation","", "the name of dislocation variable");
  params.addParam<std::string>("const_dislocation","","the name of dislocation from materials");
  params.addParam<Real>("dislocation_density",-1, "dislocation density");
  params.addParam<bool>("custom_input",false,"true: use mannully input coefficients");
  params.addParam<std::vector<Real> >("diff_i", "diffusivity of i corresponding to the mobile_i_size list");
  params.addParam<std::vector<Real> >("diff_v", "diffusivity of v corresponding to the mobile_v_size list"); //if provided source_i or source_v, the function wouldn't be used.
  params.addParam<Real>("Rid",0.6,"capture radius of interstitials by dislocations[nm]");
  params.addParam<Real>("Rvd",0.5,"capture radius of interstitials by dislocations[nm]");
  params.addParam<Real>("i_disl_bias",1.1,"intersitial cluster efficiency factor");
  params.addParam<Real>("v_disl_bias",1.0,"vacancy cluster efficiency factor");
  params.addParam<Real>("temperature",800,"system temperature [K]");
  return params;
}


AddLotsOfSink_disl::AddLotsOfSink_disl(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddLotsOfSink_disl::act()
{
//selectively add mobile specise kernels
  std::vector<int> v_size = getParam<std::vector<int> >("mobile_v_size");
  std::vector<int> i_size = getParam<std::vector<int> >("mobile_i_size");
  Real rid = getParam<Real>("Rid");
  Real rvd = getParam<Real>("Rvd");
  Real i_disl_bias = getParam<Real>("i_disl_bias");
  Real v_disl_bias = getParam<Real>("v_disl_bias");
  int _total_v = getParam<int>("number_v");
  int _total_i = getParam<int>("number_i");
  Real temp = getParam<Real>("temperature");
  bool custom = getParam<bool>("custom_input");
  std::vector<Real> vv,ii;
  if (isParamValid("diff_v") && custom == true)
    vv = getParam<std::vector<Real> >("diff_v");
  else {
    for(unsigned int i=0;i<v_size.size();i++){
      vv.push_back(diff(v_size[i],"V",temp));}
  }
  if (isParamValid("diff_i") && custom == true)
    ii = getParam<std::vector<Real> >("diff_i");
  else {
    for(unsigned int i=0;i<i_size.size();i++){
      ii.push_back(diff(i_size[i],"I",temp));}
  }
  std::string varied_disl = getParam<std::string>("dislocation");
  std::string const_disl = getParam<std::string>("const_dislocation");
  Real rho_disl = getParam<Real>("dislocation_density");
  if(varied_disl.empty() && const_disl.empty() && rho_disl<0)
     mooseError("Error: dislocation density not specified" );

  for (int cur_num = 1; cur_num <= (int)v_size.size(); cur_num++)
  {
    if(cur_num <= _total_v){
    std::string var_name_v = name() +"v"+ Moose::stringify(v_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("DislocationSink");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<Real>("Diffusivity") = vv[cur_num-1];
    params.set<Real>("DislocationCoreSize") = rvd;//Rvd = 0.5 nm
    params.set<Real>("Coef") = v_disl_bias;
    if (!varied_disl.empty())
      params.set<std::vector<VariableName> > ("VariedDislocation").push_back(varied_disl);
      //params.set<NonlinearVariableName>("VariedDislocation") = varied_disl;
    else if(rho_disl>0){
      params.set<Real>("DislocationDensityValue") = rho_disl;
    }
    else{
      params.set<std::string>("DislocationDensity") = const_disl;
    }
    _problem->addKernel("DislocationSink", "Disl_" + var_name_v + Moose::stringify(counter), params);
    printf("add DislocationSink: %s\n",var_name_v.c_str());
    counter++;
    }
  }
  
  for (int cur_num = 1; cur_num <= (int)i_size.size(); cur_num++)
  {
    if(cur_num <= _total_i){
    std::string var_name_i = name() +"i"+ Moose::stringify(i_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("DislocationSink");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<Real>("Diffusivity") = ii[cur_num-1];
    params.set<Real>("DislocationCoreSize") = rid;//Rid = 0.5 nm
    params.set<Real>("Coef") = i_disl_bias;
    if (!varied_disl.empty())
      params.set<std::vector<VariableName> > ("VariedDislocation").push_back(varied_disl);
      //params.set<NonlinearVariableName>("VariedDislocation") = varied_disl;
    else if(rho_disl>0)
      params.set<Real>("DislocationDensityValue") = rho_disl;
    else
      params.set<std::string>("DislocationDensity") = const_disl;
      //params.set<NonlinearVariableName>("VariedDislocation") = varied_disl;
    _problem->addKernel("DislocationSink", "Disl_" + var_name_i + Moose::stringify(counter), params);
    printf("add DislocationSink: %s\n",var_name_i.c_str());
    counter++;
    }
  }
}
