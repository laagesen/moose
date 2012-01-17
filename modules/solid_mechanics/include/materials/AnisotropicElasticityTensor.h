#ifndef ANISOTROPICELASTICITYTENSOR_H
#define ANISOTROPICELASTICITYTENSOR_H

#include "ElasticityTensor.h"
#include "libmesh.h"
#include "vector_value.h"
#include "dense_matrix.h"
#include "mesh.h"

/**
 * Defines an Isotropic Elasticity Tensor.
 *
 * The input is any two of the following:
 *
 * Youngs Modulus
 * Poissons Ration
 * First and Second Lame Coefficients (lambda and mu)
 * Bulk Modulus
 *
 * Internally this class uses the Lame Coefficients.
 * Everything is is transformed to these coefficients.
 *
 * Note that by default this tensor is _constant_... meaning
 * that it never changes after the first time it is computed.
 *
 * If you want to modify this behavior you can pass in
 * false to the constructor.
 */
class AnisotropicElasticityTensor : public ElasticityTensor
{
public:
  AnisotropicElasticityTensor();

  /**
   * Set the first euler angle
   */

  void setFirstEulerAngle(const Real a1);

  /**
   * Set the second euler angle
   */

  void setSecondEulerAngle(const Real a2);

  /**
   * Set the third euler angle
   */

  void setThirdEulerAngle(const Real a3);

  /**
   * Set the material constant c11
   */

  void setMaterialConstantc11(const Real c11);

  /**
   * Set the material constant c22
   */

  void setMaterialConstantc12(const Real c12);

  /**
   * Set the material constant c44
   */

  void setMaterialConstantc44(const Real c44);


protected:

  std::vector<Real> _euler_angle; // Stores Euler angeles

  Real _c11, _c12, _c44; // Material Constants

  /**
   * Fill in the matrix.
   */

  virtual void calculateEntries(unsigned int qp);

};

#endif //ANISOTROPICELASTICITYTENSOR_H
