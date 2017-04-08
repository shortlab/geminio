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

#include "AddLotsOfSingleVariable.h"
#include "MaterialParameters.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "SingleVariable.h"

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
InputParameters validParams<AddLotsOfSingleVariable>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<int> >("mobile_v_size", "A vector of mobile species");
  params.addRequiredParam<std::vector<int> >("mobile_i_size", "A vector of mobile species");
  params.addParam<bool>("custom_input",false,"true: use mannully input data");
  params.addParam<std::vector<Real> >("emit_vv", "emission coefficient from size > 1");
  params.addParam<std::vector<Real> >("emit_ii", "emission coefficient from size > 1");
  params.addParam<Real>("temperature",600.0,"Temperature");
  return params;
}


AddLotsOfSingleVariable::AddLotsOfSingleVariable(const InputParameters & params) :
    AddVariableAction(params)
{
}
//only emission of point defects of same type are considered
void
AddLotsOfSingleVariable::act()
{
  int number_v = getParam<int>("number_v");
  int number_i = getParam<int>("number_i");
  std::vector<int> v_size = getParam<std::vector<int> >("mobile_v_size");
  std::vector<int> i_size = getParam<std::vector<int> >("mobile_i_size");
  Real temp = getParam<Real>("temperature");
  bool custom = getParam<bool>("custom_input");
  std::vector<Real> vv,ii;
  int tagi = 0, tagj = 1;//used to mark mobility, 1: mobile; 0: immobile
  if(isParamValid("emit_vv") && custom == true)
    vv = getParam<std::vector<Real> >("emit_vv");
  else {
    for(int i=0;i<number_v;i++){
        if (std::find( v_size.begin(),v_size.end(),i+1) != v_size.end())//find mobile species
            tagi = 1;
        else tagi = 0;
        vv.push_back(emit(i+1,1,temp,"V","V",tagi,tagj));//emission of vacancy from i+1 cluster, a vector start with 0 (1->1)
   //     printf("size: %d, emission: %f\n",i+1,vv.back());//test
        //vv.push_back(1.0);//test
    }
  }
  if(isParamValid("emit_ii") && custom == true)
    ii = getParam<std::vector<Real> >("emit_ii");
  else {
    for(int i=0;i<number_i;i++){
        if (std::find(i_size.begin(),i_size.end(),i+1) != i_size.end())//find mobile species
            tagi = 1;
        else tagi = 0;
        ii.push_back(emit(i+1,1,temp,"I","I",tagi,tagj));//emission of vacancy from i+1 cluster, a vector start with 0 (1->1)
  //      printf("size: %d, emission: %f\n",-i-1,ii.back());//test
    }
  }


  for (int cur_num = 2; cur_num <= number_v; cur_num++)
  //vv emission
  {
    double emit_coef = 0.0;
    //self decrease (-)
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("SingleVariable");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    emit_coef = vv[cur_num-1];//cur_num to vacancy,"+"
    params.set<Real>("coeff") = emit_coef;//loss "+"
    _problem->addKernel("SingleVariable", "SingleV_" + var_name_v+ "_" +Moose::stringify(counter), params);
    //printf("add SingleVariable: %s (emit -), coef: %lf\n",var_name_v.c_str(),emit_coef);
    counter++;

    //point defects(+)
    std::string var_name_v2 = name() +"v"+ Moose::stringify(1);
    InputParameters params0 = _factory.getValidParams("SingleVariable");
    params0.set<NonlinearVariableName>("variable") = var_name_v2;
    params0.set<std::vector<VariableName> > ("secondVar").push_back(var_name_v);
    emit_coef = (-1.0)*vv[cur_num-1];
    params0.set<Real>("coeff") = emit_coef;//emission coefficent from larger size to current size; gain should be negative in kernel
    _problem->addKernel("SingleVariable", "SingleV_" + var_name_v2 + "_" + Moose::stringify(counter), params0);
    //printf("add SingleVariable: %s (%s emit +) , coef: %lf\n",var_name_v2.c_str(),var_name_v.c_str(),emit_coef);
    counter++;
    //cur_num-1 (+) (account for 2->1+1)
    std::string var_name_other = name() + "v" + Moose::stringify(cur_num-1);
    InputParameters params_other = _factory.getValidParams("SingleVariable");
    params_other.set<NonlinearVariableName>("variable") = var_name_other;
    params_other.set<std::vector<VariableName> > ("secondVar").push_back(var_name_v);
    emit_coef = (-1.0)*vv[cur_num-1];
    params_other.set<Real>("coeff") = emit_coef;//emission coefficent from larger size to current size; gain should be negative in kernel
    _problem->addKernel("SingleVariable", "SingleV_" + var_name_other + "_" + Moose::stringify(counter), params_other);
    //printf("add SingleVariable: %s (%s emit +) , coef: %lf\n",var_name_other.c_str(),var_name_v.c_str(),emit_coef);
    counter++;
  }

  for (int cur_num = 2; cur_num <= number_i; cur_num++)
  //ii emission
  {
    double emit_coef = 0.0;
    std::string var_name_i = name() +"i"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("SingleVariable");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    emit_coef = ii[cur_num-1];//cur_num to vacancy,"+"
    params.set<Real>("coeff") = emit_coef;//loss "+"
    _problem->addKernel("SingleVariable", "SingleV_" + var_name_i+ "_" +Moose::stringify(counter), params);
    //printf("add SingleVariable: %s (emit -), coef: %lf\n",var_name_i.c_str(),emit_coef);
    counter++;

    //point defects(+)
    std::string var_name_i2 = name() +"i"+ Moose::stringify(1);
    InputParameters params0 = _factory.getValidParams("SingleVariable");
    params0.set<NonlinearVariableName>("variable") = var_name_i2;
    params0.set<std::vector<VariableName> > ("secondVar").push_back(var_name_i);
    emit_coef = (-1.0)*ii[cur_num-1];
    params0.set<Real>("coeff") = emit_coef;//emission coefficent from larger size to current size; gain should be negative in kernel
    _problem->addKernel("SingleVariable", "SingleV_" + var_name_i2 + "_" + Moose::stringify(counter), params0);
    //printf("add SingleVariable: %s (%s emit +) , coef: %lf\n",var_name_i2.c_str(),var_name_i.c_str(),emit_coef);
    counter++;
    //cur_num-1 (+) (account for 2->1+1)
    std::string var_name_other = name() + "i" + Moose::stringify(cur_num-1);
    InputParameters params_other = _factory.getValidParams("SingleVariable");
    params_other.set<NonlinearVariableName>("variable") = var_name_other;
    params_other.set<std::vector<VariableName> > ("secondVar").push_back(var_name_i);
    emit_coef = (-1.0)*ii[cur_num-1];
    params_other.set<Real>("coeff") = emit_coef;//emission coefficent from larger size to current size; gain should be negative in kernel
    _problem->addKernel("SingleVariable", "SingleV_" + var_name_other + "_" + Moose::stringify(counter), params_other);
    //printf("add SingleVariable: %s (%s emit +) , coef: %lf\n",var_name_other.c_str(),var_name_i.c_str(),emit_coef);
    counter++;

    
  }
}
/*
void
AddLotsOfSingleVariable::act()
{
  unsigned int number_v = getParam<unsigned int>("number_v");
  unsigned int number_i = getParam<unsigned int>("number_i");
  std::vector<int> v_size = getParam<std::vector<int> >("mobile_v_size");
  std::vector<int> i_size = getParam<std::vector<int> >("mobile_i_size");
  Real temp = getParam<Real>("temperature");
  bool custom = getParam<bool>("custom_input");
  std::vector<Real> vv,ii;
  int tagi = 0, tagj = 0;//used to mark mobility, 1: mobile; 0: immobile
  if(isParamValid("emit_vv") && custom == true)
    vv = getParam<std::vector<Real> >("emit_vv");
  else {
    for(int i=0;i<number_v;i++){
        if (std::find( v_size.begin(),v_size.end(),i) != v_size.end())//only emit mobile species
            tagi = 1;
        else tagi = 0;
        for(int j=0;j<number_v;j++){
            if (std::find( v_size.begin(),v_size.end(),j) != v_size.end())//only emit mobile species
                tagj = 1;
            else tagj = 0;
            vv.push_back(emit(i+1,j+1,temp,"V","V",tagi,tagj));//emission of j+1 size from i+1 cluster, i+1 ought to be large than j+1, so vv essentially a low triangle matrix.
        }
    }
  }
  if(isParamValid("emit_ii") && custom == true)
    ii = getParam<std::vector<Real> >("emit_ii");
  else {
    for(int i=0;i<number_i;i++){
        if (std::find( i_size.begin(),i_size.end(),i+1) != i_size.end())//only emit mobile species
            tagi = 1;
        else tagi = 0;
        for(int j=0;j<number_i;j++){
            if (std::find( i_size.begin(),i_size.end(),j+1) != i_size.end())//only emit mobile species
                tagj = 1;
            else tagj = 0;
            double temp_coef = emit(i+1,j+1,temp,"I","I",tagi,tagj);
            ii.push_back(temp_coef);
        }
    }
  }
//ATTENTION: the emission of vacancy cluster emit an interstitial or interstitial cluster emit an vacancy is not considered
  for (unsigned int cur_num = 1; cur_num <= number_v; cur_num++)
  //vv emission
  {
    double emit_coef = 0;
    for (unsigned int j = 1; j < cur_num; j++){
      if (std::find( v_size.begin(),v_size.end(),j) != v_size.end())//only emit mobile species
        emit_coef += vv[(cur_num-1)*number_v+j-1];
    }
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("SingleVariable");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<Real>("coeff") = emit_coef;//the sum of all emission coefficent from current size to other sizes
    _problem->addKernel("SingleVariable", "SingleV_" + var_name_v+Moose::stringify(counter), params);
    printf("add SingleVariable: %s, coef: %lf\n",var_name_v.c_str(),emit_coef);
    counter++;
   //other sizes to current size if mobile
    if (std::find( v_size.begin(),v_size.end(),cur_num) != v_size.end())//mobile
    {
      for (unsigned int i = cur_num+1; i <=number_v; i++)
      {
          std::string var_name_v2 = name() +"v"+ Moose::stringify(i);
          InputParameters params0 = _factory.getValidParams("SingleVariable");
          params0.set<NonlinearVariableName>("variable") = var_name_v;
          params0.set<std::vector<VariableName> > ("secondVar").push_back(var_name_v2);
          if(2*cur_num == i)//double
            params0.set<Real>("coeff") = (-1.0)*2*vv[(i-1)*number_v+cur_num-1];//emission coefficent from larger size to current size; gain should be negative in kernel
          else
            params0.set<Real>("coeff") = (-1.0)*vv[(i-1)*number_v+cur_num-1];//emission coefficent from larger size to current size; gain should be negative in kernel
          _problem->addKernel("SingleVariable", "SingleV_" + var_name_v + Moose::stringify(counter), params0);
          printf("add SingleVariable: %s (%s) , coef: %lf\n",var_name_v.c_str(),var_name_v2.c_str(),(-1.0)*vv[(i-1)*number_v+cur_num-1]);
          counter++;
      }
    }
  }
  for (unsigned int cur_num = 1; cur_num <= number_i; cur_num++)
  //ii emission
  {
    double emit_coef = 0.0;
    for (unsigned int j = 1; j < cur_num; j++)
    {
      if (std::find( i_size.begin(),i_size.end(),j) != i_size.end()){
        emit_coef += ii[(cur_num-1)*number_i+j-1];
      }
    }
    std::string var_name_i = name() +"i"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("SingleVariable");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<Real>("coeff") = emit_coef;//Should be the sum of all emission coefficent from current size to other sizes
    _problem->addKernel("SingleVariable", "SingleV_" + var_name_i + Moose::stringify(counter), params);
    printf("add SingleVariable: %s, coef: %lf\n",var_name_i.c_str(),emit_coef);
    counter++;
    if (std::find( i_size.begin(),i_size.end(),cur_num) != i_size.end())//mobile
    {
      for (unsigned int i = cur_num+1; i <=number_i; i++)
      {
          std::string var_name_v2 = name() +"i"+ Moose::stringify(i);
          InputParameters params0 = _factory.getValidParams("SingleVariable");
          params0.set<NonlinearVariableName>("variable") = var_name_i;
          params0.set<std::vector<VariableName> > ("secondVar").push_back(var_name_v2);
          if(cur_num*2 == i)//double
            params0.set<Real>("coeff") = (-1.0)*2*ii[(i-1)*number_i+cur_num-1];//emission coefficent from larger sizes to current size; gain should be negative in kernel
          else
            params0.set<Real>("coeff") = (-1.0)*ii[(i-1)*number_i+cur_num-1];//emission coefficent from larger sizes to current size; gain should be negative in kernel
          _problem->addKernel("SingleVariable", "SingleV_" + var_name_i+ Moose::stringify(counter), params0);
          printf("add SingleVariable: %s (%s) , coef: %lf\n",var_name_i.c_str(),var_name_v2.c_str(),(-1.0)*ii[(i-1)*number_i+cur_num-1]);
          counter++;
      }
    }
  }
}
*/
