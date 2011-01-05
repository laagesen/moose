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

#ifndef DIRICHLETPOSTPROCESSORBC_H
#define DIRICHLETPOSTPROCESSORBC_H

#include "BoundaryCondition.h"


//Forward Declarations
class DirichletPostprocessorBC;

template<>
InputParameters validParams<DirichletPostprocessorBC>();

/**
 * Implements a simple constant Dirichlet BC where u is a postprocessor value on the boundary.
 */
class DirichletPostprocessorBC : public BoundaryCondition
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  DirichletPostprocessorBC(const std::string & name, InputParameters parameters);
    
  virtual ~DirichletPostprocessorBC() {}

protected:
  virtual Real computeQpResidual();

private:
  /**
   * Value of u on the boundary (from postprocessor).
   */
  std::string _postprocessor_name;
  
  Real & _value;
};

#endif //DIRICHLETPOSTPROCESSORBC_H
