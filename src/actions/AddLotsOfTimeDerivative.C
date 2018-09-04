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

#include "AddLotsOfTimeDerivative.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "TimeDerivative.h"

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
InputParameters validParams<AddLotsOfTimeDerivative>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  return params;
}

AddLotsOfTimeDerivative::AddLotsOfTimeDerivative(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddLotsOfTimeDerivative::act()
{
  unsigned int number_v = getParam<unsigned int>("number_v");
  unsigned int number_i = getParam<unsigned int>("number_i");

  for (unsigned int cur_num = 1; cur_num <= number_v; cur_num++)
  {
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("TimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    _problem->addKernel("TimeDerivative","dt_"+ var_name_v+Moose::stringify(counter), params);
    counter++;
  }
  for (unsigned int cur_num = 1; cur_num <= number_i; cur_num++)
  {
    std::string var_name_i = name() +"i"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("TimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    _problem->addKernel("TimeDerivative", "dt_"+ var_name_i+Moose::stringify(counter), params);
    counter++;
  }
}
