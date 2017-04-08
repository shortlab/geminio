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

#ifndef USEROBJECTSINGLEVARIABLE_H
#define USEROBJECTSINGLEVARIABLE_H

#include "Kernel.h"
#include "GroupConstant.h"
//Forward Declarations
class UserObjectSingleVariable;


template<>
InputParameters validParams<UserObjectSingleVariable>();

class UserObjectSingleVariable : public Kernel
{
public:
  
  UserObjectSingleVariable(const 
                            InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  int getGroupNumber(std::string);

private:

  const GroupConstant & _gc;
  std::string _call_fun;
  Real  _coeff;
  std::vector<unsigned int> _vars;
  std::vector<const VariableValue *> _v_vals;
  int groupNo;
};
#endif 
