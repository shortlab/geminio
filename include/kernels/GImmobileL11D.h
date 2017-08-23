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

#ifndef GIMMOBILEL11D_H
#define GIMMOBILEL11D_H

#include "Kernel.h"
#include "GGroup.h"

//Forward Declarations
class GImmobileL11D;


template<>
InputParameters validParams<GImmobileL11D>();

class GImmobileL11D : public Kernel
{
public:
  
  GImmobileL11D(const 
                            InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  int getGroupNumber(std::string);
  double getAbsorbCoef(int);

private:
  int _number_v;
  int _number_i;
  int _max_mobile_v; 
  int _max_mobile_i; 
  const GGroup & _gc;
  std::vector<unsigned int> _no_v_vars;
  std::vector<const VariableValue *> _val_v_vars;
  std::vector<unsigned int> _no_i_vars;
  std::vector<const VariableValue *> _val_i_vars;
  std::vector<const VariableValue *> _val_i_auxvars;
  int _cur_size;
};
#endif 
