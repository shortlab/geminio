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

#ifndef IMMOBILEDEFECTS_H
#define IMMOBILEDEFECTS_H

#include "Kernel.h"
#include "GroupConstant.h"

//Forward Declarations
class ImmobileDefects;


template<>
InputParameters validParams<ImmobileDefects>();

class ImmobileDefects : public Kernel
{
public:
  
  ImmobileDefects(const 
                            InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  int getGroupNumber(std::string);

private:
  int _number_v;
  int _number_i;
  std::vector<int> _v_size; 
  std::vector<int> _i_size;
  const GroupConstant & _gc;
  std::vector<unsigned int> _no_v_vars;
  std::vector<const VariableValue *> _val_v_vars;
  std::vector<unsigned int> _no_i_vars;
  std::vector<const VariableValue *> _val_i_vars;
  int _cur_size;
};
#endif 
