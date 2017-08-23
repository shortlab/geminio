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

#ifndef ADDRECIPROCALMEANFREEPATH1D_H
#define ADDRECIPROCALMEANFREEPATH1D_H

#include "AddVariableAction.h"

class AddReciprocalMeanFreePath1D;

template<>
InputParameters validParams<AddReciprocalMeanFreePath1D>();


class AddReciprocalMeanFreePath1D : public AddVariableAction
{
public:
  AddReciprocalMeanFreePath1D(const  InputParameters & parameters);

  virtual void act();
};

#endif // AddReciprocalMeanFreePath1D_H
