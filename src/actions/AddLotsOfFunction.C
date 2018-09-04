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

#include "AddLotsOfFunction.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddFunctionAction.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "PiecewiseLinear.h"
#include "Piecewise.h"
#include "PiecewiseLinearTimeLimit.h"

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

template<>
InputParameters validParams<AddLotsOfFunction>()
{
  //MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  //MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddFunctionAction>();
  params.addRequiredParam<std::vector<Real> >("source_v_size", "A vector of distribution of creation along ion range");
  params.addRequiredParam<std::vector<Real> >("source_i_size", "A vector of distribution of creation along ion range");
  params.addRequiredParam<std::string>("data_file", "File holding csv data for use with Piecewise");
  params.addRequiredParam<std::string>("format","Format of csv data file that is in either in columns or rows, but now columns are supported");
  params.addParam<Real>("tlimit",1.0,"timelimit for the function (lifetime) [s]");
  params.addParam<Real>("scale_factor",1.0,"scale factor for functions");
  return params;
}


AddLotsOfFunction::AddLotsOfFunction(const InputParameters &  params) :
    AddFunctionAction(params)
{
}

void
AddLotsOfFunction::act()
{
  std::vector<Real> vv = getParam<std::vector<Real> >("source_v_size");
  std::vector<Real> ii = getParam<std::vector<Real> >("source_i_size");
  std::string file_name = getParam<std::string>("data_file");
  std::string format = getParam<std::string>("format");
  Real t_limit = getParam<Real>("tlimit");
  Real _scale_factor = getParam<Real>("scale_factor");
  std::ifstream file(file_name.c_str());
  if (!file.good())
    mooseError("Error opening file '" + file_name );
  if (format.compare("columns") != 0)
    mooseError("Invalid option for format: "+format+" Use 'columns'.");

  std::string line;
  std::vector<Real> myvec;
  while (getline(file,line))
  {
    //Replace all commas with spaces
    while (size_t pos=line.find(','))
    {
      if (pos == line.npos)
        break;
      line.replace(pos,1,1,' ');
    }

    //Harvest floats separated by whitespace
    std::istringstream iss(line);
    Real f;
    while (iss>>f)
    {
      myvec.push_back(f);
    }
  }
  file.close();
  unsigned int cols = (vv.size()+ii.size()+1);
  if (myvec.size() % cols !=0)
  {
      mooseError("Check data format!(should be position + vacancy cluster sizes + intersitial cluster sizes with each in column)");
  } 
  unsigned int rows = myvec.size()/cols; //total rows;vv.size()+ii.size()+1 is the total column number; the first column is position along the ion range
  std::vector<Real> x;
  x.reserve(rows);
  for (unsigned int i=0; i<rows; ++i)
  {
      x.push_back(myvec[i*cols]);//position along the ion range
  }
  for (unsigned int cur_num = 1; cur_num <= vv.size(); cur_num++)
  {
    std::vector<Real> y;
    y.reserve(rows);
    for (unsigned int i=0; i<rows; ++i)
    {
      y.push_back(myvec[i*cols+cur_num]);//gain should be negative in the kernel
    }
    std::string fun_name_v = name() +"v"+ Moose::stringify(vv[cur_num-1]);
    InputParameters params = _factory.getValidParams("PiecewiseLinearTimeLimit");
    params.set<FunctionName>("function") = fun_name_v;
    params.set<int>("axis") = 0;//x axis, one dimension
    params.set<Real>("scale_factor") = _scale_factor;//x axis, one dimension
    params.set<std::vector<Real> >("x") = x;
    params.set<std::vector<Real> >("y") = y;
    params.set<Real>("tlimit") = t_limit;
    _problem->addFunction("PiecewiseLinearTimeLimit", fun_name_v, params);
  }
  for (unsigned int cur_num = 1; cur_num <= ii.size(); cur_num++)
  {
    std::vector<Real> y;
    y.reserve(rows);
    for (unsigned int i=0; i<rows; ++i)
    {
      y.push_back(myvec[i*cols+cur_num+vv.size()]);//gain should be negative in the kernel
    }    
    std::string fun_name_i = name() +"i"+ Moose::stringify(ii[cur_num-1]);
    InputParameters params = _factory.getValidParams("PiecewiseLinearTimeLimit");
    params.set<FunctionName>("function") = fun_name_i;
    params.set<int>("axis") = 0;//x axis, one dimension
    params.set<Real>("scale_factor") = _scale_factor;//x axis, one dimension
    params.set<std::vector<Real> >("x") = x;
    params.set<std::vector<Real> >("y") = y;
    params.set<Real>("tlimit") = t_limit;
    _problem->addFunction("PiecewiseLinearTimeLimit", fun_name_i, params);
  }
}
