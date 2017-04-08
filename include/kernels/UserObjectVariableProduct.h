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

#ifndef USEROBJECTVARIABLEPRODUCT_H
#define USEROBJECTVARIABLEPRODUCT_H

#include "Kernel.h"
#include "GroupConstant.h"

//Forward Declarations
class UserObjectVariableProduct;


template<>
InputParameters validParams<UserObjectVariableProduct>();

class UserObjectVariableProduct : public Kernel
{
public:
  
  UserObjectVariableProduct(const 
                            InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  int getGroupNumber(std::string);

private:
  const GroupConstant & _gc;
  Real  _coeff;
  std::vector<unsigned int> _vars;
  /// Coupled primary species concentrations.
  std::vector<const VariableValue *> _v_vals;
//  std::vector<std::string> _var_names;
  int groupa,groupb;
};
#endif 
