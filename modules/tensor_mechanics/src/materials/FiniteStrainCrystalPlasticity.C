/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainCrystalPlasticity.h"
#include "petscblaslapack.h"

template<>
InputParameters validParams<FiniteStrainCrystalPlasticity>()
{
  InputParameters params = validParams<FiniteStrainMaterial>();
  params.addClassDescription("Crystal Plasticity base class: FCC system with power law flow rule implemented");
  params.addRequiredParam<int >("nss", "Number of slip systems");
  params.addParam<std::vector<Real> >("gprops", "Initial values of slip system resistances");
  params.addParam<std::vector<Real> >("hprops", "Hardening properties");
  params.addParam<std::vector<Real> >("flowprops", "Parameters used in slip rate equations");
  params.addRequiredParam<std::string>("slip_sys_file_name", "Name of the file containing the slip system");
  params.addParam<std::string>("slip_sys_res_prop_file_name", "", "Name of the file containing the initial values of slip system resistances");
  params.addParam<std::string>("slip_sys_flow_prop_file_name", "", "Name of the file containing the values of slip rate equation parameters");
  params.addParam<std::string>("slip_sys_hard_prop_file_name", "", "Name of the file containing the values of hardness evolution parameters");
  params.addParam<Real>("rtol", 1e-6, "Constitutive stress residue relative tolerance");
  params.addParam<Real>("abs_tol", 1e-6, "Constitutive stress residue absolute tolerance");
  params.addParam<Real>("gtol", 1e2, "Constitutive slip system resistance residual tolerance");
  params.addParam<Real>("slip_incr_tol", 2e-2, "Maximum allowable slip in an increment");
  params.addParam<unsigned int>("maxiter", 100 , "Maximum number of iterations for stress update");
  params.addParam<unsigned int>("maxitergss", 100 , "Maximum number of iterations for slip system resistance update");
  params.addParam<unsigned int>("num_slip_sys_flowrate_props", 2, "Number of flow rate properties for a slip system"); // Used for reading flow rate parameters
  params.addParam<UserObjectName>("read_prop_user_object","The ElementReadPropertyFile GeneralUserObject to read element specific property values from file");
  MooseEnum tan_mod_options("exact none","none");//Type of read
  params.addParam<MooseEnum>("tan_mod_type", tan_mod_options, "Type of tangent moduli for preconditioner: default elastic");
  MooseEnum intvar_read_options("slip_sys_file slip_sys_res_file none","none");
  params.addParam<MooseEnum>("intvar_read_type", intvar_read_options, "Read from options for initial value of internal variables: Default from .i file");
  params.addParam<unsigned int>("num_slip_sys_props", 0, "Number of slip system specific properties provided in the file containing slip system normals and directions");
  params.addParam<bool>("save_euler_angle", false , "Saves the Euler angles as Material Property if true");
  params.addParam<bool>("gen_random_stress_flag", false, "Flag to generate random stress to perform time cutback on constitutive failure");
  params.addParam<bool>("input_random_scaling_var", false, "Flag to input scaling variable: _Cijkl(0,0,0,0) when false");
  params.addParam<Real>("random_scaling_var", 1e9, "Random scaling variable: Large value can cause non-positive definiteness");
  params.addParam<unsigned int>("random_seed", 2000, "Random integer used to generate random stress when constitutive failure occurs");
  params.addParam<unsigned int>("maximum_substep_iteration", 1, "Maximum number of substep iteration");
  params.addParam<bool>("use_line_search", false, "Use line search in constitutive update");
  params.addParam<Real>("min_line_search_step_size", 0.01, "Minimum line search step size");

  return params;
}

FiniteStrainCrystalPlasticity::FiniteStrainCrystalPlasticity(const InputParameters & parameters) :
    FiniteStrainMaterial(parameters),
    _nss(getParam<int>("nss")),
    _gprops(getParam<std::vector<Real> >("gprops")),
    _hprops(getParam<std::vector<Real> >("hprops")),
    _flowprops(getParam<std::vector<Real> >("flowprops")),
    _slip_sys_file_name(getParam<std::string>("slip_sys_file_name")),
    _slip_sys_res_prop_file_name(getParam<std::string>("slip_sys_res_prop_file_name")),
    _slip_sys_flow_prop_file_name(getParam<std::string>("slip_sys_flow_prop_file_name")),
    _slip_sys_hard_prop_file_name(getParam<std::string>("slip_sys_hard_prop_file_name")),
    _rtol(getParam<Real>("rtol")),
    _abs_tol(getParam<Real>("abs_tol")),
    _gtol(getParam<Real>("gtol")),
    _slip_incr_tol(getParam<Real>("slip_incr_tol")),
    _maxiter(getParam<unsigned int>("maxiter")),
    _maxiterg(getParam<unsigned int>("maxitergss")),
    _num_slip_sys_flowrate_props(getParam<unsigned int>("num_slip_sys_flowrate_props")),
    _read_prop_user_object(isParamValid("read_prop_user_object") ? & getUserObject<ElementPropertyReadFile>("read_prop_user_object") : NULL),
    _tan_mod_type(getParam<MooseEnum>("tan_mod_type")),
    _intvar_read_type(getParam<MooseEnum>("intvar_read_type")),
    _num_slip_sys_props(getParam<unsigned int>("num_slip_sys_props")),
    _save_euler_angle(getParam<bool>("save_euler_angle")),
    _gen_rndm_stress_flag(getParam<bool>("gen_random_stress_flag")),
    _input_rndm_scale_var(getParam<bool>("input_random_scaling_var")),
    _rndm_scale_var(getParam<Real>("random_scaling_var")),
    _rndm_seed(getParam<unsigned int>("random_seed")),
    _max_substep_iter(getParam<unsigned int>("maximum_substep_iteration")),
    _use_line_search(getParam<bool>("use_line_search")),
    _min_lsrch_step(getParam<Real>("min_line_search_step_size")),
    _fp(declareProperty<RankTwoTensor>("fp")), // Plastic deformation gradient
    _fp_old(declarePropertyOld<RankTwoTensor>("fp")), // Plastic deformation gradient of previous increment
    _pk2(declareProperty<RankTwoTensor>("pk2")), // 2nd Piola Kirchoff Stress
    _pk2_old(declarePropertyOld<RankTwoTensor>("pk2")), // 2nd Piola Kirchoff Stress of previous increment
    _lag_e(declareProperty<RankTwoTensor>("lage")), // Lagrangian strain
    _lag_e_old(declarePropertyOld<RankTwoTensor>("lage")), // Lagrangian strain of previous increment
    _gss(declareProperty<std::vector<Real> >("gss")), // Slip system resistances
    _gss_old(declarePropertyOld<std::vector<Real> >("gss")), // Slip system resistances of previous increment
    _acc_slip(declareProperty<Real>("acc_slip")), // Accumulated slip
    _acc_slip_old(declarePropertyOld<Real>("acc_slip")), // Accumulated alip of previous increment
    _update_rot(declareProperty<RankTwoTensor>("update_rot")), // Rotation tensor considering material rotation and crystal orientation
    _update_rot_old(declarePropertyOld<RankTwoTensor>("update_rot")),
    _deformation_gradient_old(declarePropertyOld<RankTwoTensor>("deformation gradient"))
{
  if (_save_euler_angle)
  {
    _euler_ang = &declareProperty< std::vector<Real> >("euler_ang");
    _euler_ang_old = &declarePropertyOld< std::vector<Real> >("euler_ang");
  }

  if (!_input_rndm_scale_var)
    _rndm_scale_var = _Cijkl(0,0,0,0);

  _err_tol = false;

  _tau.resize(_nss);
  _slip_incr.resize(_nss);
  _dslipdtau.resize(_nss);

  _mo.resize(_nss*LIBMESH_DIM);
  _no.resize(_nss*LIBMESH_DIM);

  _s0.resize(_nss);

  if (_num_slip_sys_props > 0)
    _slip_sys_props.resize(_nss * _num_slip_sys_props);

  _pk2_tmp.zero();
  _pk2_tmp_old.zero();
  _gss_tmp.resize(_nss);
  _gss_tmp_old.resize(_nss);
  _delta_dfgrd.zero();

  _first_step_iter = false;
  _last_step_iter = false;
  //Initialize variables in the first iteration of substepping
  _first_substep = true;

  _read_from_slip_sys_file = false;
  if (_intvar_read_type == "slip_sys_file")
    _read_from_slip_sys_file = true;

  if (_read_from_slip_sys_file && ! ( _num_slip_sys_props > 0 ))
    mooseError("Crystal Plasticity Error: Specify number of internal variable's initial values to be read from slip system file");

  getSlipSystems();

  RankTwoTensor::initRandom( _rndm_seed );
}

void FiniteStrainCrystalPlasticity::initQpStatefulProperties()
{
  _stress[_qp].zero();

  _fp[_qp].zero();
  _fp[_qp].addIa(1.0);

  _pk2[_qp].zero();
  _acc_slip[_qp] = 0.0;
  _lag_e[_qp].zero();

  _update_rot[_qp].zero();
  _update_rot[_qp].addIa(1.0);

  if (_save_euler_angle)
  {
    (*_euler_ang)[_qp].resize(LIBMESH_DIM);
    (*_euler_ang_old)[_qp].resize(LIBMESH_DIM);
  }

  initSlipSysProps(); // Initializes slip system related properties
  initAdditionalProps();
}

void
FiniteStrainCrystalPlasticity::initSlipSysProps()
{
  switch (_intvar_read_type)
  {
  case 0:
    assignSlipSysRes();
    break;
  case 1:
    readFileInitSlipSysRes();
    break;
  default:
    getInitSlipSysRes();
  }

  if (_slip_sys_flow_prop_file_name.length()!=0)
    readFileFlowRateParams();
  else
    getFlowRateParams();

  if (_slip_sys_hard_prop_file_name.length()!=0)
    readFileHardnessParams();
  else
    getHardnessParams();
}

void
FiniteStrainCrystalPlasticity::assignSlipSysRes()
{
  _gss[_qp].resize(_nss);
  _gss_old[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
    _gss[_qp][i] = _gss_old[_qp][i] = _slip_sys_props[i];
}


// Read initial slip system resistances  from .txt file. See test.
void
FiniteStrainCrystalPlasticity::readFileInitSlipSysRes()
{
  _gss[_qp].resize(_nss);
  _gss_old[_qp].resize(_nss);

  MooseUtils::checkFileReadable(_slip_sys_res_prop_file_name);

  std::ifstream file;
  file.open(_slip_sys_res_prop_file_name.c_str());

  for (unsigned int i = 0; i < _nss; ++i)
    if (!(file >> _gss[_qp][i]))
      mooseError("Error FiniteStrainCrystalPlasticity: Premature end of slip_sys_res_prop file");

  file.close();
}

// Read initial slip system resistances  from .i file
void
FiniteStrainCrystalPlasticity::getInitSlipSysRes()
{
  if (_gprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading slip system resistance properties: Specify input in .i file or in slip_sys_res_prop_file or in slip_sys_file");

  _gss[_qp].resize(_nss, 0.0);
  _gss_old[_qp].resize(_nss, 0.0);

  unsigned int num_data_grp = 3; //Number of data per group e.g. start_slip_sys, end_slip_sys, value

  for (unsigned int i = 0; i < _gprops.size()/num_data_grp; ++i)
  {
    Real vs,ve;
    unsigned int is, ie;

    vs = _gprops[i * num_data_grp];
    ve = _gprops[i * num_data_grp + 1];

    if (vs <= 0 || ve <= 0)
      mooseError( "FiniteStrainCrystalPLasticity: Indices in gss property read must be positive integers: is = " << vs << " ie = " << ve );

    if (vs != floor(vs) || ve != floor(ve))
      mooseError("FiniteStrainCrystalPLasticity: Error in reading slip system resistances: Values specifying start and end number of slip system groups should be integer");

    is = static_cast<unsigned int>(vs);
    ie = static_cast<unsigned int>(ve);

    if (is > ie)
      mooseError("FiniteStrainCrystalPLasticity: Start index is = " << is << " should be greater than end index ie = " << ie << " in slip system resistance property read");

    for (unsigned int j = is; j <= ie; ++j)
      _gss[_qp][j-1] = _gprops[i * num_data_grp + 2];
  }

  for (unsigned int i = 0; i < _nss; ++i)
    if (_gss[_qp][i] <= 0.0)
      mooseError("FiniteStrainCrystalPLasticity: Value of resistance for slip system " << i + 1 << " non positive");
}

// Read flow rate parameters from .txt file. See test.
void
FiniteStrainCrystalPlasticity::readFileFlowRateParams()
{
  _a0.resize(_nss);
  _xm.resize(_nss);

  MooseUtils::checkFileReadable(_slip_sys_flow_prop_file_name);

  std::ifstream file;
  file.open(_slip_sys_flow_prop_file_name.c_str());

  std::vector<Real> vec;
  vec.resize(_num_slip_sys_flowrate_props);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    for (unsigned int j = 0; j < _num_slip_sys_flowrate_props; ++j)
      if (!(file >> vec[j]))
        mooseError("Error FiniteStrainCrystalPlasticity: Premature end of slip_sys_flow_rate_param file");

      _a0[i]=vec[0];
      _xm[i]=vec[1];
  }

  file.close();
}

// Read flow rate parameters from .i file
void
FiniteStrainCrystalPlasticity::getFlowRateParams()
{
  if (_flowprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading flow rate  properties: Specify input in .i file or a slip_sys_flow_prop_file_name");

  _a0.resize(_nss, 0.0);
  _xm.resize(_nss, 0.0);

  unsigned int num_data_grp = 2 + _num_slip_sys_flowrate_props; //Number of data per group e.g. start_slip_sys, end_slip_sys, value1, value2, ..

  for (unsigned int i = 0; i < _flowprops.size()/num_data_grp; ++i)
  {
    Real vs,ve;
    unsigned int is, ie;

    vs = _flowprops[i * num_data_grp];
    ve = _flowprops[i * num_data_grp + 1];

    if (vs <= 0 || ve <= 0)
      mooseError("FiniteStrainCrystalPLasticity: Indices in flow rate parameter read must be positive integers: is = " << vs << " ie = " << ve );

    if (vs != floor(vs) || ve != floor(ve))
      mooseError("FiniteStrainCrystalPLasticity: Error in reading flow props: Values specifying start and end number of slip system groups should be integer");

    is = static_cast<unsigned int>(vs);
    ie = static_cast<unsigned int>(ve);

    if (is > ie)
      mooseError("FiniteStrainCrystalPLasticity: Start index is = " << is << " should be greater than end index ie = " << ie << " in flow rate parameter read");

    for (unsigned int j = is; j <= ie; ++j)
    {
      _a0[j-1] = _flowprops[i * num_data_grp + 2];
      _xm[j-1] = _flowprops[i * num_data_grp + 3];
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (!(_a0[i] > 0.0 && _xm[i] > 0.0))
    {
      mooseWarning( "FiniteStrainCrystalPlasticity: Non-positive flow rate parameters " << _a0[i] << "," << _xm[i] );
      break;
    }
  }
}

// Read hardness parameters from .txt file
void
FiniteStrainCrystalPlasticity::readFileHardnessParams()
{
}

// Read hardness parameters from .i file
void
FiniteStrainCrystalPlasticity::getHardnessParams()
{
  if (_hprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading hardness properties: Specify input in .i file or a slip_sys_hard_prop_file_name");

  _r = _hprops[0];
  _h0 = _hprops[1];
  _tau_init = _hprops[2];
  _tau_sat = _hprops[3];
}

void
FiniteStrainCrystalPlasticity::getEulerAngles()
{
  if (_read_prop_user_object)
  {
    _Euler_angles(0) = _read_prop_user_object->getData(_current_elem, 0);
    _Euler_angles(1) = _read_prop_user_object->getData(_current_elem, 1);
    _Euler_angles(2) = _read_prop_user_object->getData(_current_elem, 2);
  }
}

// Calculate crystal rotation tensor from Euler Angles
void
FiniteStrainCrystalPlasticity::getEulerRotations()
{
  Real phi1, phi, phi2;
  Real cp, cp1, cp2, sp, sp1, sp2;
  RankTwoTensor RT;
  Real pi = libMesh::pi;

  phi1 = _Euler_angles(0) * (pi/180.0);
  phi =  _Euler_angles(1) * (pi/180.0);
  phi2 = _Euler_angles(2) * (pi/180.0);

  cp1 = std::cos(phi1);
  cp2 = std::cos(phi2);
  cp = std::cos(phi);

  sp1 = std::sin(phi1);
  sp2 = std::sin(phi2);
  sp = std::sin(phi);

  RT(0,0) = cp1 * cp2 - sp1 * sp2 * cp;
  RT(0,1) = sp1 * cp2 + cp1 * sp2 * cp;
  RT(0,2) = sp2 * sp;
  RT(1,0) = -cp1 * sp2 - sp1 * cp2 * cp;
  RT(1,1) = -sp1 * sp2 + cp1 * cp2 * cp;
  RT(1,2) = cp2 * sp;
  RT(2,0) = sp1 * sp;
  RT(2,1) = -cp1 * sp;
  RT(2,2) = cp;

  _crysrot = RT.transpose();
}

// Read slip systems from file
void
FiniteStrainCrystalPlasticity::getSlipSystems()
{
  Real vec[LIBMESH_DIM];
  std::ifstream fileslipsys;

  MooseUtils::checkFileReadable(_slip_sys_file_name);

  fileslipsys.open(_slip_sys_file_name.c_str());

  for (unsigned int i = 0; i < _nss; ++i)
  {
    // Read the slip normal
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      if ( ! (fileslipsys >> vec[j]) )
        mooseError("Crystal Plasticity Error: Premature end of file reading slip system file \n");

    // Normalize the vectors
    Real mag;
    mag = std::pow(vec[0], 2.0) + std::pow(vec[1], 2.0) + std::pow(vec[2], 2.0);
    mag = std::pow(mag, 0.5);

    for (unsigned j = 0; j < LIBMESH_DIM; ++j)
      _no[i*LIBMESH_DIM+j] = vec[j]/mag;

    // Read the slip direction
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      if ( ! (fileslipsys >> vec[j]) )
        mooseError("Crystal Plasticity Error: Premature end of file reading slip system file \n");

    // Normalize the vectors
    mag = std::pow(vec[0], 2.0) + std::pow(vec[1], 2.0) + std::pow(vec[2], 2.0);
    mag = std::pow(mag, 0.5);

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      _mo[i*LIBMESH_DIM+j] = vec[j] / mag;

    mag = 0.0;
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      mag += _mo[i*LIBMESH_DIM+j] * _no[i*LIBMESH_DIM+j];

    if (std::abs(mag) > 1e-8)
      mooseError("Crystal Plasicity Error: Slip direction and normal not orthonormal, System number = " << i << "\n");

    if (_read_from_slip_sys_file)
      for (unsigned int j = 0; j < _num_slip_sys_props; ++j)
        if (!(fileslipsys >> _slip_sys_props[i * _num_slip_sys_props + j]))
          mooseError("Crystal Plasticity Error: Premature end of file reading slip system file - check in slip system file read input options/values\n");
  }

  fileslipsys.close();
}

// Initialize addtional stateful material properties
void
FiniteStrainCrystalPlasticity::initAdditionalProps()
{
}

/**
 * Solves stress residual equation using NR.
 * Updates slip system resistances iteratively.
 */
void FiniteStrainCrystalPlasticity::computeQpStress()
{
  unsigned int substep_iter = 1;//Depth of substepping; Limited to maximum substep iteration
  unsigned int num_substep = 1;//Calculated from substep_iter as 2^substep_iter
  Real dt_original = _dt;//Stores original _dt; Reset at the end of solve
  _first_substep = true;//Initialize variables at substep_iter = 1

  if (_max_substep_iter > 1)
  {
    _dfgrd_tmp_old = _deformation_gradient_old[_qp];
    if (_dfgrd_tmp_old.det() == 0)
      _dfgrd_tmp_old.addIa(1.0);

    _delta_dfgrd = _deformation_gradient[_qp] - _dfgrd_tmp_old;
    _err_tol = true;//Indicator to continue substepping
  }

  //Substepping loop
  while (_err_tol && _max_substep_iter > 1)
  {
    _dt = dt_original/num_substep;

    for (unsigned int istep = 0; istep < num_substep; ++istep)
    {
      _first_step_iter = false;
      if (istep == 0)
        _first_step_iter = true;

      _last_step_iter = false;
      if (istep == num_substep - 1)
        _last_step_iter = true;

      _dfgrd_scale_factor = (static_cast<Real>(istep)+1)/num_substep;

      preSolveQp();
      solveQp();

      if (_err_tol)
      {
        substep_iter++;
        num_substep*=2;
        break;
      }
    }

    _first_substep = false;//Prevents reinitialization
    _dt = dt_original;//Resets dt

#ifdef DEBUG
    if (substep_iter > _max_substep_iter)
      mooseWarning("FiniteStrainCrystalPlasticity: Failure with substepping");
#endif

    if (!_err_tol || substep_iter > _max_substep_iter)
      postSolveQp();//Evaluate variables after successful solve or indicate failure
  }

  //No substepping
  if (_max_substep_iter == 1)
  {
    preSolveQp();
    solveQp();
    postSolveQp();
  }
}

void
FiniteStrainCrystalPlasticity::preSolveQp()
{
  //Initialize variable
  if (_first_substep)
  {
    _Jacobian_mult[_qp].zero();//Initializes jacobian for preconditioner
    getEulerAngles();
    getEulerRotations();

    calc_schmid_tensor();

    RealTensorValue rot;

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        rot(i,j) = _crysrot(i,j);

    _elasticity_tensor[_qp] = _Cijkl;
    _elasticity_tensor[_qp].rotate(rot);
  }

  if (_max_substep_iter == 1)
    _dfgrd_tmp = _deformation_gradient[_qp];//Without substepping
  else
    _dfgrd_tmp = _dfgrd_scale_factor * _delta_dfgrd + _dfgrd_tmp_old;

  _err_tol = false;
}

void
FiniteStrainCrystalPlasticity::solveQp()
{
  preSolveStatevar();
  solveStatevar();
  if (_err_tol)
    return;
  postSolveStatevar();
}

void
FiniteStrainCrystalPlasticity::postSolveQp()
{
  if (_err_tol)
  {
    _err_tol = false;
    if ( _gen_rndm_stress_flag )
      _stress[_qp] = RankTwoTensor::genRandomSymmTensor( _rndm_scale_var, 1.0 );
    else
      mooseError("FiniteStrainCrystalPlasticity: Constitutive failure");
  }
  else
  {
    _stress[_qp] = _fe * _pk2[_qp] * _fe.transpose()/_fe.det();

    _Jacobian_mult[_qp] += calcTangentModuli();//Calculate jacobian for preconditioner

    RankTwoTensor iden;
    iden.addIa(1.0);

    _lag_e[_qp] = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] - iden;
    _lag_e[_qp] = _lag_e[_qp] * 0.5;

    RankTwoTensor rot;
    rot = get_current_rotation(_deformation_gradient[_qp]); // Calculate material rotation
    _update_rot[_qp] = rot * _crysrot;

    if (_save_euler_angle)
      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
        (*_euler_ang)[_qp][i] = _Euler_angles(i);
  }
}

void
FiniteStrainCrystalPlasticity::preSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    for (unsigned i = 0; i < _nss; ++i)
      _gss_tmp[i] = _gss_old[_qp][i];
    _accslip_tmp_old = _acc_slip_old[_qp];
  }
  else
  {
    if (_first_step_iter)
    {
      for (unsigned i = 0; i < _nss; ++i)
        _gss_tmp[i] = _gss_tmp_old[i] = _gss_old[_qp][i];
      _accslip_tmp_old = _acc_slip_old[_qp];
    }
    else
      for (unsigned i = 0; i < _nss; ++i)
        _gss_tmp[i] = _gss_tmp_old[i];
  }
}

void
FiniteStrainCrystalPlasticity::solveStatevar()
{
  Real gmax, gdiff;
  unsigned int iterg;
  std::vector<Real> gss_prev(_nss);

  gmax = 1.1 * _gtol;
  iterg = 0;

  while (gmax > _gtol && iterg < _maxiterg) // Check for slip system resistance update tolerance
  {
    preSolveStress();
    solveStress();
    if (_err_tol)
      return;
    postSolveStress();

    for (unsigned i = 0; i < _nss; ++i)
      gss_prev[i] = _gss_tmp[i];

    update_slip_system_resistance(); // Update slip system resistance

    gmax = 0.0;
    for (unsigned i = 0; i < _nss; ++i)
    {
      gdiff = std::abs(gss_prev[i] - _gss_tmp[i]); // Calculate increment size

      if (gdiff > gmax)
        gmax = gdiff;
    }
    iterg++;
  }

  if (iterg == _maxiterg)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainCrystalPLasticity: Hardness Integration error gmax" << gmax << "\n");
#endif
    _err_tol = true;
  }
}

void
FiniteStrainCrystalPlasticity::postSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    for (unsigned i = 0; i < _nss; ++i)
      _gss[_qp][i] = _gss_tmp[i];
    _acc_slip[_qp] = _accslip_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      for (unsigned i = 0; i < _nss; ++i)
        _gss[_qp][i] = _gss_tmp[i];
      _acc_slip[_qp] = _accslip_tmp;
    }
    else
    {
      for (unsigned i = 0; i < _nss; ++i)
        _gss_tmp_old[i] = _gss_tmp[i];
      _accslip_tmp_old = _accslip_tmp;
    }
  }
}

void
FiniteStrainCrystalPlasticity::preSolveStress()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _pk2_tmp = _pk2_old[_qp];
    _fp_old_inv = _fp_old[_qp].inverse();
    _fp_inv = _fp_old_inv;
  }
  else
  {
    if (_first_step_iter)
    {
      _pk2_tmp = _pk2_tmp_old = _pk2_old[_qp];
      _fp_old_inv = _fp_old[_qp].inverse();
    }
    else
      _pk2_tmp = _pk2_tmp_old;

    _fp_inv = _fp_old_inv;
  }
}

void
FiniteStrainCrystalPlasticity::solveStress()
{
  unsigned int iter = 0;
  RankTwoTensor resid, dpk2;
  RankFourTensor jac;
  Real rnorm, rnorm0, rnorm_prev;

  calc_resid_jacob(resid, jac); // Calculate stress residual
  if (_err_tol)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainCrystalPLasticity: Slip increment exceeds tolerance - Element number " << _current_elem->id() << " Gauss point = " << _qp);
#endif
    return;
  }

  rnorm = resid.L2norm();
  rnorm0 = rnorm;

  while (rnorm > _rtol * rnorm0 && rnorm0 > _abs_tol && iter <  _maxiter) // Check for stress residual tolerance
  {
    dpk2 = - jac.invSymm() * resid; // Calculate stress increment
    _pk2_tmp = _pk2_tmp + dpk2; // Update stress

    calc_resid_jacob(resid, jac); // Calculate stress residual
    if (_err_tol)
    {
#ifdef DEBUG
      mooseWarning("FiniteStrainCrystalPLasticity: Slip increment exceeds tolerance - Element number " << _current_elem->id() << " Gauss point = " << _qp);
#endif
      return;
    }

    rnorm_prev = rnorm;
    rnorm = resid.L2norm();

    if (_use_line_search && rnorm > rnorm_prev && !line_search_update(rnorm_prev, dpk2))
    {
#ifdef DEBUG
      mooseWarning("FiniteStrainCrystalPLasticity: Failed with line search");
#endif
      _err_tol = true;
      return;
    }

    iter++;
  }

  if (iter >= _maxiter)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainCrystalPLasticity: Stress Integration error rmax = " << rnorm);
#endif
    _err_tol = true;
  }
}

void
FiniteStrainCrystalPlasticity::postSolveStress()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _fp[_qp] = _fp_inv.inverse();
    _pk2[_qp] = _pk2_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      _fp[_qp] = _fp_inv.inverse();
      _pk2[_qp] = _pk2_tmp;
    }
    else
    {
      _fp_old_inv = _fp_inv;
      _pk2_tmp_old = _pk2_tmp;
    }
  }
}

// Update slip system resistance. Overide to incorporate new slip system resistance laws
void
FiniteStrainCrystalPlasticity::update_slip_system_resistance()
{
  updateGss();
}

/**
 * Old function to update slip system resistances.
 * Kept to avoid code break at computeQpstress
 */
void
FiniteStrainCrystalPlasticity::updateGss()
{
  std::vector<Real> hb(_nss);
  Real qab;

  Real a = _hprops[4]; // Kalidindi

  _accslip_tmp = _accslip_tmp_old;
  for (unsigned int i=0; i < _nss; ++i)
    _accslip_tmp += std::abs(_slip_incr[i]);

  // Real val = std::cosh(_h0 * _accslip_tmp / (_tau_sat - _tau_init)); // Karthik
  // val = _h0 * std::pow(1.0/val,2.0); // Kalidindi

  for (unsigned int i = 0; i < _nss; ++i)
    // hb[i]=val;
    hb[i] = _h0 * std::pow(std::abs(1.0 - _gss_tmp[i]/_tau_sat),a) * copysign(1.0,1.0-_gss_tmp[i]/_tau_sat);

  for (unsigned int i=0; i < _nss; ++i)
  {
    if (_max_substep_iter == 1)//No substepping
      _gss_tmp[i] = _gss_old[_qp][i];
    else
      _gss_tmp[i] = _gss_tmp_old[i];

    for (unsigned int j = 0; j < _nss; ++j)
    {
      unsigned int iplane, jplane;
      iplane = i/3;
      jplane = j/3;

      if (iplane == jplane) // Kalidindi
        qab = 1.0;
      else
        qab = _r;

      _gss_tmp[i] += qab * hb[j] * std::abs(_slip_incr[j]);
    }
  }
}

// Calculates stress residual equation and jacobian
void
FiniteStrainCrystalPlasticity::calc_resid_jacob( RankTwoTensor & resid, RankFourTensor & jac)
{
  calcResidual( resid );
  if (_err_tol)
    return;
  calcJacobian( jac );
}

void
FiniteStrainCrystalPlasticity::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;

  iden.zero();
  iden.addIa(1.0);

  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    _tau[i] = ce_pk2.doubleContraction(_s0[i]);

  getSlipIncrements(); // Calculate dslip,dslipdtau

  if (_err_tol)
    return;

  eqv_slip_incr.zero();
  for (unsigned int i = 0; i < _nss; ++i)
    eqv_slip_incr += _s0[i] * _slip_incr[i];

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;

  pk2_new = _elasticity_tensor[_qp] * ee;

  resid = _pk2_tmp - pk2_new;
}

void
FiniteStrainCrystalPlasticity::calcJacobian( RankFourTensor &jac )
{
  RankFourTensor dfedfpinv, deedfe, dfpinvdpk2;

  std::vector<RankTwoTensor> dtaudpk2(_nss), dfpinvdslip(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    dtaudpk2[i] = _s0[i];
    dfpinvdslip[i] = - _fp_old_inv * _s0[i];
  }

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        dfedfpinv(i,j,k,j) = _dfgrd_tmp(i,k);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        deedfe(i,j,k,i) = deedfe(i,j,k,i) + _fe(k,j) * 0.5;
        deedfe(i,j,k,j) = deedfe(i,j,k,j) + _fe(k,i) * 0.5;
      }

  for (unsigned int i = 0; i < _nss; ++i)
    dfpinvdpk2 += (dfpinvdslip[i] * _dslipdtau[i]).outerProduct(dtaudpk2[i]);

  jac = RankFourTensor::IdentityFour() - (_elasticity_tensor[_qp] * deedfe * dfedfpinv * dfpinvdpk2);
}

// Calculate slip increment,dslipdtau. Override to modify.
void
FiniteStrainCrystalPlasticity::getSlipIncrements()
{
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_incr[i] = _a0[i] * std::pow(std::abs(_tau[i] / _gss_tmp[i]), 1.0 / _xm[i]) * copysign(1.0, _tau[i]) * _dt;
    if (std::abs(_slip_incr[i]) > _slip_incr_tol)
    {
      _err_tol = true;
#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded " << std::abs(_slip_incr[i]));
#endif
      return;
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
    _dslipdtau[i] = _a0[i] / _xm[i] * std::pow(std::abs(_tau[i] / _gss_tmp[i]), 1.0 / _xm[i] - 1.0) / _gss_tmp[i] * _dt;
}

// Calls getMatRot to perform RU factorization of a tensor.
RankTwoTensor
FiniteStrainCrystalPlasticity::get_current_rotation(const RankTwoTensor & a)
{
  return getMatRot(a);
}

// Performs RU factorization of a tensor
RankTwoTensor
FiniteStrainCrystalPlasticity::getMatRot(const RankTwoTensor & a)
{
  RankTwoTensor rot;
  RankTwoTensor c, diag, evec;
  PetscScalar cmat[LIBMESH_DIM][LIBMESH_DIM], work[10];
  PetscReal w[LIBMESH_DIM];
  PetscBLASInt nd = LIBMESH_DIM,
               lwork = 10,
               info;

  c = a.transpose() * a;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      cmat[i][j] = c(i,j);

  LAPACKsyev_("V", "U", &nd, &cmat[0][0], &nd, w, work, &lwork, &info);

  if (info != 0)
    mooseError("FiniteStrainCrystalPLasticity: DSYEV function call in getMatRot function failed");

  diag.zero();

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    diag(i,i) = std::pow(w[i], 0.5);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      evec(i,j) = cmat[i][j];

  rot = a * ((evec.transpose() * diag * evec).inverse());

  return rot;
}

// Calculates tangent moduli which is used for global solve
void
FiniteStrainCrystalPlasticity::computeQpElasticityTensor()
{

}

ElasticityTensorR4
FiniteStrainCrystalPlasticity::calcTangentModuli()
{
  ElasticityTensorR4 tan_mod;

  switch ( _tan_mod_type )
  {
  case 0:
    tan_mod = elastoPlasticTangentModuli();
    break;
  default:
    tan_mod = elasticTangentModuli();
  }

  return tan_mod;
}

void
FiniteStrainCrystalPlasticity::calc_schmid_tensor()
{
  std::vector<Real> mo(LIBMESH_DIM*_nss),no(LIBMESH_DIM*_nss);

  // Update slip direction and normal with crystal orientation
  for (unsigned int i = 0; i < _nss; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo[i*LIBMESH_DIM+j] = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo[i*LIBMESH_DIM+j] = mo[i*LIBMESH_DIM+j] + _crysrot(j,k) * _mo[i*LIBMESH_DIM+k];
    }

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no[i*LIBMESH_DIM+j] = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        no[i*LIBMESH_DIM+j] = no[i*LIBMESH_DIM+j] + _crysrot(j,k) * _no[i*LIBMESH_DIM+k];
    }
  }

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        _s0[i](j,k) = mo[i*LIBMESH_DIM+j] * no[i*LIBMESH_DIM+k];
}


ElasticityTensorR4
FiniteStrainCrystalPlasticity::elastoPlasticTangentModuli()
{
  ElasticityTensorR4 tan_mod;
  RankTwoTensor pk2fet, fepk2;
  RankFourTensor deedfe, dsigdpk2dfe;

  // Fill in the matrix stiffness material property

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        deedfe(i,j,k,i) = deedfe(i,j,k,i) + _fe(k,j) * 0.5;
        deedfe(i,j,k,j) = deedfe(i,j,k,j) + _fe(k,i) * 0.5;
      }

  dsigdpk2dfe = _fe.mixedProductIkJl(_fe) * _elasticity_tensor[_qp] * deedfe;

  pk2fet = _pk2_tmp * _fe.transpose();
  fepk2 = _fe * _pk2_tmp;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int l = 0; l < LIBMESH_DIM; ++l)
      {
        tan_mod(i,j,i,l) = tan_mod(i,j,i,l) + pk2fet(l,j);
        tan_mod(i,j,j,l) = tan_mod(i,j,j,l) + fepk2(i,l);
      }

  tan_mod += dsigdpk2dfe;

  Real je = _fe.det();
  if (je > 0.0)
    tan_mod /= je;

  return tan_mod;
}


ElasticityTensorR4
FiniteStrainCrystalPlasticity::elasticTangentModuli()
{
  return _elasticity_tensor[_qp];//update jacobian_mult
}

bool
FiniteStrainCrystalPlasticity::line_search_update(const Real rnorm_prev, const RankTwoTensor dpk2)
{
  Real rnorm;
  RankTwoTensor resid;
  Real step = 1.0;

  do
  {
    _pk2_tmp = _pk2_tmp - step * dpk2;
    step /= 2.0;
    _pk2_tmp = _pk2_tmp + step * dpk2;

    calcResidual(resid);
    rnorm = resid.L2norm();
  }
  while (rnorm > rnorm_prev && step > _min_lsrch_step);

  if (rnorm > rnorm_prev && step <= _min_lsrch_step)
    return false;

  return true;
}

