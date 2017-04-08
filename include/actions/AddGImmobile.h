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

#ifndef ADDGIMMOBILE_H
#define ADDGIMMOBILE_H

#include "AddVariableAction.h"

class AddGImmobile;

template<>
InputParameters validParams<AddGImmobile>();


class AddGImmobile : public AddVariableAction
{
public:
  AddGImmobile(const  InputParameters & parameters);

  virtual void act();
};

#endif // AddGImmobile_H
