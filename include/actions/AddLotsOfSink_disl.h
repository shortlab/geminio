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

#ifndef ADDLOTSOFSINK_DISL_H
#define ADDLOTSOFSINK_DISL_H

#include "AddVariableAction.h"

class AddLotsOfSink_disl;

template<>
InputParameters validParams<AddLotsOfSink_disl>();


class AddLotsOfSink_disl : public AddVariableAction
{
public:
  AddLotsOfSink_disl(const  InputParameters & parameters);

  virtual void act();
};

#endif // AddLotsOfSink_disl_H
