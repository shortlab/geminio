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

#ifndef USEROBJECTDIFFUSION_H
#define USEROBJECTDIFFUSION_H

#include "Diffusion.h"
#include "GroupConstant.h"
//Forward Declarations
class UserObjectDiffusion;


template<>
InputParameters validParams<UserObjectDiffusion>();

class UserObjectDiffusion : public Diffusion
{
public:
  
  UserObjectDiffusion(const 
                            InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  int getGroupNumber(std::string);

private:
  Real _coeff;
  const GroupConstant & _gc;
  int groupNo;
};
#endif 
