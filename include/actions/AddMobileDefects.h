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

#ifndef ADDMOBILEDEFECTS_H
#define ADDMOBILEDEFECTS_H

#include "AddVariableAction.h"

class AddMobileDefects;

template<>
InputParameters validParams<AddMobileDefects>();


class AddMobileDefects : public AddVariableAction
{
public:
  AddMobileDefects(const  InputParameters & parameters);

  virtual void act();
};

#endif // AddMobileDefects_H
