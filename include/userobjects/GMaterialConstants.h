/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/

#ifndef GMATERIALCONSTANTS_H
#define GMATERIALCONSTANTS_H

#include "GeneralUserObject.h"

class GMaterialConstants : public GeneralUserObject
{
public:
  GMaterialConstants(const InputParameters & parameters);

  virtual ~GMaterialConstants(){}
  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  virtual Real absorb(int,int,std::string,std::string,double,int,int) const;
  virtual Real absorbVV(int,int,int,double) const;
  virtual Real absorbVI(int,int,int,double) const;
  virtual Real absorbII(int,int,int,double) const;
  virtual Real emit(int,int,double,std::string,std::string,int,int) const;
  virtual Real disl_ksq(int,std::string,double,int=1) const;//1 denotes mobile
  virtual Real diff(int,std::string,double) const;
  Real atomic_vol;

protected:
  Real _rho_d;
  Real _i_bias;
  Real _v_bias;
};

template<>
InputParameters validParams<GMaterialConstants>();

#endif
