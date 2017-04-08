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

#ifndef GROUPINGTEST_H
#define GROUPINGTEST_H

#include "GeneralUserObject.h"
#include "GMaterialConstants.h"

class GroupingTest : public GMaterialConstants
{
public:
  GroupingTest(const InputParameters & parameters);

  ~GroupingTest(){}
  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  Real absorb(int,int,std::string,std::string,double,int,int) const;
  Real emit(int,int,double,std::string,std::string,int,int) const;
  double energy(int,std::string,std::string) const;
  double D_prefactor(int,std::string) const;
  double diff(int, std::string,double) const;

private:
  Real _i_bias;
  Real _v_bias;
};

template<>
InputParameters validParams<GroupingTest>();

#endif
