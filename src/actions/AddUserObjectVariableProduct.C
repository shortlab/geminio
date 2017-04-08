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

#include "AddUserObjectVariableProduct.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "UserObjectVariableProduct.h"

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
static unsigned int counter = 0;//define a counter for the name of kernels

template<>
InputParameters validParams<AddUserObjectVariableProduct>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<int> >("mobile_v_size", "A vector of mobile species sizes");
  params.addRequiredParam<std::vector<int> >("mobile_i_size", "A vector of mobile species sizes");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  params.addParam<Real>("temperature",600.0,"Temperature");
  return params;
}


AddUserObjectVariableProduct::AddUserObjectVariableProduct(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddUserObjectVariableProduct::act()
{
  int number_v = getParam<int>("number_v");
  int number_i = getParam<int>("number_i");
  std::vector<int> v_size = getParam<std::vector<int> >("mobile_v_size");
  std::vector<int> i_size = getParam<std::vector<int> >("mobile_i_size");
  Real temp = getParam<Real>("temperature");
  std::string uo = getParam<std::string>("group_constant");

  // First comes with the kernels for vacancy clusters
    for (int cur_num = 1; cur_num <= number_v; cur_num++)
    {
      std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
      for (int cur_num2 = 1; cur_num2 <= number_v-cur_num; cur_num2++)//ATTENTION: largest number_v-cur_num to ensure conservation (no outlier reaction)
      //vv reaction down(-)
      {
        Real  coef = 0.0;
        if(std::find(v_size.begin(),v_size.end(),cur_num2) != v_size.end() || std::find(v_size.begin(),v_size.end(),cur_num) != v_size.end())//confirm it as mobile
        {
          std::string var_name2;
          var_name2 = name() +"v"+ Moose::stringify(cur_num2);
          InputParameters params1 = _factory.getValidParams("UserObjectVariableProduct");
          params1.set<NonlinearVariableName>("variable") = var_name_v;
          params1.set<UserObjectName>("user_object") = uo;
          if (cur_num2 == cur_num){//if the same species, then no need to add second variable, the coefficient is multiplied by 2
            coef = 2*1.0;
            params1.set<Real>("coeff") = coef;
          }
          else{
            params1.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
            coef = 1.0;
            params1.set<Real>("coeff") = coef;
          }
          _problem->addKernel("UserObjectVariableProduct", "VarProd_" + var_name_v+ "_" + Moose::stringify(counter), params1);
          //printf("add UserObjectVariableProduct: %s (%s and %s), coef: %lf\n",var_name_v.c_str(),var_name_v.c_str(),var_name2.c_str(),coef);
          counter++;
        }
      }

      for (int cur_num2 = 1; cur_num2 <= number_i; cur_num2++)
      //vi reaction down(-)
      {
        Real coef = 0.0;
        if(std::find(i_size.begin(),i_size.end(),cur_num2) != i_size.end() || std::find(v_size.begin(),v_size.end(),cur_num) != v_size.end() )//confirm it as mobile
        {
          std::string var_name2;
          var_name2 = name() +"i"+ Moose::stringify(cur_num2);
          InputParameters params1 = _factory.getValidParams("UserObjectVariableProduct");
          params1.set<NonlinearVariableName>("variable") = var_name_v;
          params1.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
          coef =  1.0;
          params1.set<Real>("coeff") = coef;
          params1.set<UserObjectName>("user_object") = uo;
          _problem->addKernel("UserObjectVariableProduct", "VarProd_" + var_name_v + "_" + Moose::stringify(counter), params1);
          //printf("add UserObjectVariableProduct: %s and %s, coef: %lf\n",var_name_v.c_str(),var_name2.c_str(),coef);
          counter++;
        }
      }


    }

    int  pairs_v = (number_v+1)/2;

    for ( int cur_num = 1; cur_num <= number_v; cur_num++)
    {
      std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
      for (int j=1; j<= (int)(cur_num/2); j++)
      // vv reaction up(+)
      {
        Real coef = 0.0;
        if(std::find(v_size.begin(),v_size.end(),j) != v_size.end() || std::find(v_size.begin(),v_size.end(),cur_num-j) != v_size.end())//at least one is mobile
        { 
          std::string var_name2 = name() +"v"+ Moose::stringify(j);
          std::string var_name3 = name() +"v"+ Moose::stringify(cur_num-j);
          InputParameters params = _factory.getValidParams("UserObjectVariableProduct");
          params.set<NonlinearVariableName>("variable") = var_name_v;
          std::vector<VariableName> coupled_vars;
          coupled_vars.push_back(var_name2);
          coupled_vars.push_back(var_name3);
          params.set<std::vector<VariableName> > ("coupled_vars") = coupled_vars;
          coef = (-1.0)*1.0;
          params.set<Real>("coeff") = coef;//gain should be negative in the kernel;
          params.set<UserObjectName>("user_object") = uo;
          _problem->addKernel("UserObjectVariableProduct", "VarProd_" + var_name_v+ "_" + Moose::stringify(counter), params);
          //printf("add UserObjectVariableProduct: %s (%s and %s), coef: %lf\n",var_name_v.c_str(),var_name2.c_str(),var_name3.c_str(),coef);
          counter++;
        }
      }
      for (int j=cur_num+1;j<= number_v;j++)
      // vi reaction up(+)
      {
        Real coef = 0.0;
        if(std::find(v_size.begin(),v_size.end(),j) != v_size.end() || std::find(i_size.begin(),i_size.end(),j-cur_num) != i_size.end())//at least one is mobile
        { 
          if(j-cur_num <= number_i)//guarantee the maximum size of intersitital clusters satisfy
          {
            std::string var_name2 = name() +"v"+ Moose::stringify(j);
            std::string var_name3 = name() +"i"+ Moose::stringify(j-cur_num);
            InputParameters params = _factory.getValidParams("UserObjectVariableProduct");
            params.set<NonlinearVariableName>("variable") = var_name_v;

            std::vector<VariableName> coupled_vars;
            coupled_vars.push_back(var_name2);
            coupled_vars.push_back(var_name3);
            params.set<std::vector<VariableName> > ("coupled_vars") = coupled_vars;
            coef = (-1.0)*1.0;
            params.set<Real>("coeff") = coef;//gain should be negative in the kernel, factor 2 account for mutual position
            params.set<UserObjectName>("user_object") = uo;
            _problem->addKernel("UserObjectVariableProduct", "VarProd_" + var_name_v+ "_" + Moose::stringify(counter), params);
            //printf("add UserObjectVariableProduct: %s (%s and %s), coef: %lf\n",var_name_v.c_str(),var_name2.c_str(),var_name3.c_str(),coef);
            counter++;
          }
        }
      }
    }

  // Second comes with the kernels for interstitial clusters
    for (int cur_num = 1; cur_num <= number_i; cur_num++)
    {
      std::string var_name_i = name() +"i"+ Moose::stringify(cur_num);
      for (int cur_num2 = 1; cur_num2 <= number_i-cur_num; cur_num2++)//ATTENTION: largest number_v-cur_num to ensure conservation (no outlier reaction)
      // ii reaction down(-)
      {
        Real coef = 0.0;
        if(std::find(i_size.begin(),i_size.end(),cur_num2) != i_size.end() || std::find(i_size.begin(),i_size.end(),cur_num) != i_size.end() )//confirm it as mobile
        {
	    std::string var_name2;
	    var_name2 = name() +"i"+ Moose::stringify(cur_num2);
	    InputParameters params1 = _factory.getValidParams("UserObjectVariableProduct");
	    params1.set<NonlinearVariableName> ("variable") = var_name_i;
            params1.set<UserObjectName>("user_object") = uo;
	    if (cur_num == cur_num2){//same species, factor 2 to account for mutual position
              coef = 2*1.0;//ii[(cur_num-1)*number_i+cur_num2-1];
	      params1.set<Real>("coeff") = coef;
            }
            else{
              coef = 1.0;//ii[(cur_num-1)*number_i+cur_num2-1];
	      params1.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
	      params1.set<Real>("coeff") = coef;
            }
	    _problem->addKernel("UserObjectVariableProduct", "VarProd_" + var_name_i+ "_" + Moose::stringify(counter), params1);
            //printf("add UserObjectVariableProduct: %s (%s and %s), coef: %lf\n",var_name_i.c_str(),var_name_i.c_str(),var_name2.c_str(),coef);
            counter++;
        }
      }
      for (int cur_num2 = 1; cur_num2 <= number_v; cur_num2++)
      // iv reaction down(-)
      {
        Real coef = 0.0;
        if(std::find(v_size.begin(),v_size.end(),cur_num2) != v_size.end() || std::find(i_size.begin(),i_size.end(),cur_num) != i_size.end() )//confirm it as mobile
        {
          std::string var_name2;
          var_name2 = name() +"v"+ Moose::stringify(cur_num2);
          InputParameters params1 = _factory.getValidParams("UserObjectVariableProduct");
          params1.set<NonlinearVariableName> ("variable") = var_name_i;
          params1.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
          coef = 1.0;
          params1.set<Real>("coeff") = coef;
          params1.set<UserObjectName>("user_object") = uo;
          _problem->addKernel("UserObjectVariableProduct", "VarProd_" + var_name_i+ "_" + Moose::stringify(counter), params1);
          //printf("add UserObjectVariableProduct: %s (%s and %s), coef: %lf\n",var_name_i.c_str(),var_name_i.c_str(),var_name2.c_str(),coef);
          counter++;
        }
      }

    }

    int  pairs_i = (number_i+1)/2;
    for (int cur_num = 1; cur_num <= number_i; cur_num++)
    {
      std::string var_name_i = name() +"i"+ Moose::stringify(cur_num);
      //for (int j=1; j<= pairs_i; j++)
      //for (int j=1; j< cur_num; j++)
      for (int j=1; j<= (int)(cur_num/2); j++)
      //ii reaction up(+)
      {
        Real coef = 0.0;
        if(std::find(i_size.begin(),i_size.end(),j) != i_size.end() || std::find(i_size.begin(),i_size.end(),cur_num-j) != i_size.end())//at least one is mobile
        { 
            std::string var_name2 = name() +"i"+ Moose::stringify(j);
            std::string var_name3 = name() +"i"+ Moose::stringify(cur_num-j);
            InputParameters params = _factory.getValidParams("UserObjectVariableProduct");
            params.set<NonlinearVariableName>("variable") = var_name_i;

            std::vector<VariableName> coupled_vars;
            coupled_vars.push_back(var_name2);
            coupled_vars.push_back(var_name3);
            params.set<std::vector<VariableName> > ("coupled_vars") = coupled_vars;
            coef = (-1.0)*1.0;//ii[(j-1)*number_i+cur_num-j-1];
            params.set<Real>("coeff") = coef;//gain should be negative in the kernel; TODO: add 2 for test
            params.set<UserObjectName>("user_object") = uo;
            _problem->addKernel("UserObjectVariableProduct", "VarProd_" + var_name_i+ "_" + Moose::stringify(counter), params);
            //printf("add UserObjectVariableProduct: %s (%s and %s), coef: %lf\n",var_name_i.c_str(),var_name2.c_str(),var_name3.c_str(),coef);
            counter++;
        }
      }
      for (int j=cur_num+1;j<= number_i;j++)
      //iv reaction up(+)
      {
        Real coef = 0.0;
        if(std::find(i_size.begin(),i_size.end(),j) != i_size.end() || std::find(v_size.begin(),v_size.end(),j-cur_num) != v_size.end())//at least one is mobile
        { 
          if(j-cur_num <= number_v)//guarantee the maximum size of intersitital clusters satisfy
          {
            std::string var_name2 = name() +"i"+ Moose::stringify(j);
            std::string var_name3 = name() +"v"+ Moose::stringify(j-cur_num);
            InputParameters params = _factory.getValidParams("UserObjectVariableProduct");
            params.set<NonlinearVariableName>("variable") = var_name_i;

            std::vector<VariableName> coupled_vars;
            coupled_vars.push_back(var_name2);
            coupled_vars.push_back(var_name3);
            params.set<std::vector<VariableName> > ("coupled_vars") = coupled_vars;
            coef = (-1.0)*1.0;
            params.set<Real>("coeff") = coef;//gain should be negative in the kernel//*2
            params.set<UserObjectName>("user_object") = uo;
            _problem->addKernel("UserObjectVariableProduct", "VarProd_" + var_name_i+ "_" + Moose::stringify(counter), params);
            //printf("add UserObjectVariableProduct: %s (%s and %s), coef: %lf\n",var_name_i.c_str(),var_name2.c_str(),var_name3.c_str(),coef);//*2
            counter++;
          }
        }
      }
    }
}
