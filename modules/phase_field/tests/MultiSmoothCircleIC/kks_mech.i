#
# KKS coupled with elasticity
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 280
  ny = 280
  nz = 0
  xmin = 0
  xmax = 140
  ymin = 0
  ymax = 140
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]



[Variables]
  # order parameter
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]

  # solute concentration
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]

  # chemical potential
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]

  # solute phase concentration (matrix)
  [./cm]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
  # solute phase concentration (precipitate)
  [./cp]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./c_ic]
    type = MultiSmoothCircleIC
    variable = c
    numbub = 5
    bubspac = 2
    invalue = 0.235
    outvalue = 0.130
    int_width = 1
    radius_variation_type = normal
    radius_variation = 0.3
    numtries = 10000
    radius = 2
  [../]
  [./eta_ic]
    type = MultiSmoothCircleIC
    variable = eta
    numbub = 5
    bubspac = 2
    invalue = 1
    outvalue = 0
    int_width = 1
    radius_variation_type = normal
    radius_variation = 0.3
    numtries = 10000
    radius = 2
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  # Chemical free energy of the matrix
  [./fm]
    type = DerivativeParsedMaterial
    f_name = fm
    args = 'cm'
    function = '6.55*(cm-0.13)^2'
  [../]
  # Elastic energy of the matrix
  [./elastic_free_energy_m]
    type = ElasticEnergyMaterial
    base_name = matrix
    f_name = fe_m
    args = ' '
  [../]
# Total free energy of the matrix
  [./Total_energy_matrix]
    type = DerivativeSumMaterial
    f_name = f_total_matrix
    sum_materials = 'fm fe_m'
    args = 'cm'
  [../]

  # Free energy of the precipitate phase
  [./fp]
    type = DerivativeParsedMaterial
    f_name = fp
    args = 'cp'
    function = '6.55*(cp-0.235)^2'
  [../]

# Elastic energy of the precipitate
  [./elastic_free_energy_p]
    type = ElasticEnergyMaterial
    base_name = ppt
    f_name = fe_p
    args = ' '
  [../]

# Total free energy of the precipitate
  [./Total_energy_ppt]
    type = DerivativeSumMaterial
    f_name = f_total_ppt
    sum_materials = 'fp fe_p'
    args = 'cp'
  [../]



  # h(eta)
  [./h_eta]
    type = SwitchingFunctionMaterial
    h_order = HIGH
    eta = eta
  [../]

  # g(eta)
  [./g_eta]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'M   L   kappa    misfit'
    prop_values = '0.7 7.0 0.06382  0.0007'
  [../]

  #Mechanical properties
  [./Stiffness_matrix]
    type = ComputeElasticityTensor
    C_ijkl = '103.3 74.25 74.25 103.3 74.25 103.3 46.75 46.75 46.75'
    base_name = matrix
    fill_method = symmetric9
  [../]
  [./Stiffness_ppt]
    type = ComputeElasticityTensor
    C_ijkl = '100.7 71.45 71.45 100.7 71.45 100.7 50.10 50.10 50.10'
    base_name = ppt
    fill_method = symmetric9
  [../]
  [./stress_matrix]
    type = ComputeLinearElasticStress
    base_name = matrix
  [../]
  [./stress_ppt]
    type = ComputeLinearElasticStress
    base_name = ppt
  [../]
  [./strain_matrix]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = matrix
  [../]
  [./strain_ppt]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = ppt
  [../]
  [./eigen_strain]
    type = ComputeEigenstrain
    base_name = ppt
    eigen_base = '1 1 0 0 0 0'
    prefactor = misfit
  [../]
  [./global_stress]
    type = TwoPhaseStressMaterial
    base_A = matrix
    base_B = ppt
  [../]
  [./global_strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]

  # enforce c = (1-h(eta))*cm + h(eta)*cp
  [./PhaseConc]
    type = KKSPhaseConcentration
    ca       = cm
    variable = cp
    c        = c
    eta      = eta
  [../]

  # enforce pointwise equality of chemical potentials
  [./ChemPotVacancies]
    type = KKSPhaseChemicalPotential
    variable = cm
    cb       = cp
    fa_name  = f_total_matrix
    fb_name  = f_total_ppt
  [../]

  #
  # Cahn-Hilliard Equation
  #
  [./CHBulk]
    type = KKSSplitCHCRes
    variable = c
    ca       = cm
    cb       = cp
    fa_name  = f_total_matrix
    fb_name  = f_total_ppt
    w        = w
  [../]

  [./dcdt]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./ckernel]
    type = SplitCHWRes
    mob_name = M
    variable = w
  [../]

  #
  # Allen-Cahn Equation
  #
  [./ACBulkF]
    type = KKSACBulkF
    variable = eta
    fa_name  = f_total_matrix
    fb_name  = f_total_ppt
    w        = 0.1544
    args = 'cp cm'
  [../]
  [./ACBulkC]
    type = KKSACBulkC
    variable = eta
    ca       = cm
    cb       = cp
    fa_name  = f_total_matrix
    fb_name  = f_total_ppt
  [../]
  [./ACInterface]
    type = ACInterface
    variable = eta
    kappa_name = kappa
  [../]
  [./detadt]
    type = TimeDerivative
    variable = eta
  [../]
[]

[Postprocessors]
  [./sum_c]
    type =  ElementIntegralVariablePostprocessor
    variable = c
  [../]
  [./sum_eta]
    type =  ElementIntegralVariablePostprocessor
    variable = eta
  [../]
  [./n_particles]
    type = FeatureFloodCount
    variable = eta
    bubble_volume_file = bubble_volume
    threshold = 0.5
    #use_less_than_threshold_comparison = false
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       lu            nonzero'

  l_max_its = 30
  nl_max_its = 10
  num_steps = 2
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  dt = 0.1

  #[./TimeStepper]
  #  type = IterationAdaptiveDT
  #  dt = 0.1
  #  optimal_iterations = 6
  #  iteration_window = 2
  #[../]
[]

#
# Precondition using handcoded off-diagonal terms
#
[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]


[Outputs]
  [./exodus]
    type = Exodus
    #interval = 10
    sync_times = '2000 5000 10000 20000 30000'
  [../]
  checkpoint = true
  #interval = 20
  #print_perf_log = true
  #[./console]
  #  type = Console
  #  #output_on = 'FAILED NONLINEAR linear TIMESTEP_BEGIN TIMESTEP_END'
  #  all_variable_norms = true
  #  output_file = true
  #[../]
[]
