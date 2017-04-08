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

#ifndef ADDUSEROBJECTDISLOCATIONSINK_H
#define ADDUSEROBJECTDISLOCATIONSINK_H

#include "AddVariableAction.h"

class AddUserObjectDislocationSink;

template<>
InputParameters validParams<AddUserObjectDislocationSink>();


class AddUserObjectDislocationSink : public AddVariableAction
{
public:
  AddUserObjectDislocationSink(const  InputParameters & parameters);

  virtual void act();
};

#endif // AddUserObjectDislocationSink_H
