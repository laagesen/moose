//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MECHANICALCONTACTLMTEST_H
#define MECHANICALCONTACTLMTEST_H

#include "MortarConstraint.h"

template <ComputeStage>
class MechanicalContactLMTest;

declareADValidParams(MechanicalContactLMTest);

template <ComputeStage compute_stage>
class MechanicalContactLMTest : public MortarConstraint<compute_stage>
{
public:
  MechanicalContactLMTest(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

  usingMortarConstraintMembers;
};

#endif /* MECHANICALCONTACTLMTEST_H */
