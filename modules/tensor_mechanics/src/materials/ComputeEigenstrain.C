/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeEigenstrain.h"

template<>
InputParameters validParams<ComputeEigenstrain>()
{
  InputParameters params = validParams<ComputeStressFreeStrainBase>();
  params.addClassDescription("Computes a constant Eigenstrain");
  params.addRequiredParam<std::vector<Real> >("eigen_base","Vector of values defining the constant base tensor for the Eigenstrain");
  params.addParam<MaterialPropertyName>("prefactor", 1.0, "Name of material defining the variable dependence");
  return params;
}

ComputeEigenstrain::ComputeEigenstrain(const InputParameters & parameters) :
    ComputeStressFreeStrainBase(parameters),
    _prefactor(getMaterialProperty<Real>("prefactor"))
{
  _eigen_base_tensor.fillFromInputVector(getParam<std::vector<Real> >("eigen_base"));
}

void
ComputeEigenstrain::computeQpStressFreeStrain()
{
  //Define Eigenstrain
  _stress_free_strain[_qp] = _eigen_base_tensor * _prefactor[_qp];
}

