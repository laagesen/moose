//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "NodeFaceConstraint.h"
#include "MooseEnum.h"

#include <vector>
#include <unordered_map>

class RANFSTieNode;
class PenetrationInfo;

namespace libMesh
{
template <typename>
class NumericVector;
}

template <>
InputParameters validParams<RANFSTieNode>();

class RANFSTieNode : public NodeFaceConstraint
{
public:
  RANFSTieNode(const InputParameters & parameters);

  bool shouldApply() override;
  void residualSetup() override;
  void jacobianSetup() override;
  void timestepSetup() override;
  void initialSetup() override;
  bool overwriteSlaveResidual() override;
  bool overwriteSlaveJacobian() override;
  void computeSlaveValue(NumericVector<Number> & solution) override;
  void computeJacobian() override;
  void computeOffDiagJacobian(unsigned int) override;

protected:
  virtual Real computeQpSlaveValue() override;

  virtual Real computeQpResidual(Moose::ConstraintType type) override;
  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;

  const MooseEnum _component;
  const unsigned int _mesh_dimension;
  NumericVector<Number> & _residual_copy;

  std::vector<unsigned int> _vars;
  std::vector<MooseVariable *> _var_objects;
  Real _lagrange_multiplier;
  const Node * _nearest_node;
  std::unordered_map<dof_id_type, Real> _node_to_lm;
  dof_id_type _master_index;

  dof_id_type _dof_number;
  dof_id_type _master_dof_number;
  std::vector<numeric_index_type> _master_cols;
  std::vector<Number> _master_values;
};
