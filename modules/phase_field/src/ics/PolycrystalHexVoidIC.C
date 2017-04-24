/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PolycrystalHexVoidIC.h"
#include "MooseMesh.h"

InputParameters
PolycrystalHexVoidIC::actionParameters()
{
  InputParameters params = validParams<MultiSmoothCircleIC>();

  params.addRequiredParam<unsigned int>("op_num", "Number of order parameters");
  params.addRequiredParam<unsigned int>(
      "grain_num", "Number of grains being represented by the order parameters");
  params.addRequiredParam<unsigned int>("edge_bub",
                                        "Number of bubbles on grain edges (trijunctions)");

  params.addParam<unsigned int>("rand_seed", 12444, "The random seed");
  return params;
}

template <>
InputParameters
validParams<PolycrystalHexVoidIC>()
{
  InputParameters params = PolycrystalHexVoidIC::actionParameters();
  MooseEnum structure_options("grains voids");
  params.addRequiredParam<MooseEnum>("structure_type",
                                     structure_options,
                                     "Which structure type is being initialized, grains or voids");
  params.addParam<unsigned int>("op_index",
                                0,
                                "The index for the current "
                                "order parameter, not needed if "
                                "structure_type = voids");
  return params;
}

PolycrystalHexVoidIC::PolycrystalHexVoidIC(const InputParameters & parameters)
  : MultiSmoothCircleIC(parameters),
    _structure_type(getParam<MooseEnum>("structure_type")),
    _op_num(getParam<unsigned int>("op_num")),
    _grain_num(getParam<unsigned int>("grain_num")),
    _op_index(getParam<unsigned int>("op_index")),
    _edge_bub(getParam<unsigned int>("edge_bub")),
    _rand_seed(getParam<unsigned int>("rand_seed"))
{
  if (_invalue < _outvalue)
    mooseError("PolycrystalHexVoidIC requires that the voids be "
               "represented with invalue > outvalue");
  if (_numbub == 0)
    mooseError("PolycrystalHexVoidIC requires numbub > 0. If you want no voids to "
               "be "
               "represented, use invalue = outvalue. In general, you should use "
               "PolycrystalReducedIC to represent Hex grain structures without "
               "voids.");
  if (_grain_num != 4)
    mooseError("PolycrystalHexVoidIC currently only supports 4 grains. Check setting of grain_num "
               "parameter.");
}

void
PolycrystalHexVoidIC::initialSetup()
{
  if (_op_num <= _op_index)
    mooseError("op_index is too large in CircleGrainVoidIC");

  MooseRandom::seed(getParam<unsigned int>("rand_seed"));
  // Set up domain bounds with mesh tools
  for (unsigned int i = 0; i < LIBMESH_DIM; i++)
  {
    _bottom_left(i) = _mesh.getMinInDimension(i);
    _top_right(i) = _mesh.getMaxInDimension(i);
  }
  _range = _top_right - _bottom_left;

  // Create _centerpoints and _assigned_op vectors
  computeGrainCenters();

  // Call initial setup from MultiSmoothCircleIC to create _centers and _radii
  // for voids
  MultiSmoothCircleIC::initialSetup();
}

void
PolycrystalHexVoidIC::computeCircleRadii()
{
  _radii.resize(_numbub + _edge_bub);

  for (unsigned int i = 0; i < _numbub + _edge_bub; i++)
  {
    // Vary bubble radius
    switch (_radius_variation_type)
    {
      case 0: // Uniform distribution
        _radii[i] = _radius * (1.0 + (1.0 - 2.0 * _random.rand(_tid)) * _radius_variation);
        break;
      case 1: // Normal distribution
        _radii[i] = _random.randNormal(_tid, _radius, _radius_variation);
        break;
      case 2: // No variation
        _radii[i] = _radius;
    }

    _radii[i] = std::max(_radii[i], 0.0);
  }
}

void
PolycrystalHexVoidIC::computeCircleCenters()
{
  _centers.resize(_numbub + _edge_bub);

  // Point center_0, center_1;
  // center_0(0) = 0;
  // center_0(1) = 3.0/16.0 * _range(1);
  // center_0(2) = 0;
  // center_1(0) = 0.25 * _range(0);
  // center_1(1) = 5.0/16.0 * _range(1);
  // center_1(2) = 0;
  //
  // _centers[0] = center_0;
  // _centers[1] = center_1;

  // First place bubbles on grain edges (trijunctions in Hex geometry)
  for (unsigned int vp = 0; vp < _edge_bub; ++vp)
  {
    bool try_again;
    unsigned int num_tries = 0;

    do
    {
      try_again = false;
      num_tries++;

      if (num_tries > _max_num_tries)
        mooseError("Too many tries of assigning void centers on grain edges in "
                   "PolycrystalHexVoidIC");

      // Pick from 8 possible trijunctions in the 4-grain hexagonal arrangement
      const unsigned int trij_index =
          static_cast<unsigned int>(std::floor(MooseRandom::rand() * 8.0));

      Point rand_point;

      switch (trij_index)
      {
        case 0:
          rand_point(0) = _bottom_left(0);
          rand_point(1) = _bottom_left(1) + 3.0 / 16.0 * _range(1);
        case 1:
          rand_point(0) = 0;
          rand_point(1) = 3.0 / 16.0 * _range(1);
        case 2:
          rand_point(0) = 0;
          rand_point(1) = 3.0 / 16.0 * _range(1);
        case 3:
          rand_point(0) = 0;
          rand_point(1) = 3.0 / 16.0 * _range(1);
        case 4:
          rand_point(0) = 0;
          rand_point(1) = 3.0 / 16.0 * _range(1);
        case 5:
          rand_point(0) = 0;
          rand_point(1) = 3.0 / 16.0 * _range(1);
        case 6:
          rand_point(0) = 0;
          rand_point(1) = 3.0 / 16.0 * _range(1);
        case 7:
          rand_point(0) = 0;
          rand_point(1) = 3.0 / 16.0 * _range(1);

        default:
          mooseError("Invalid location for grain edge bubble in PolycrystalHexVoidIC");
          break;
      }

      rand_point(2) = 0;

    } while (try_again == true);
  }

  // Next place bubbles on grain boundaries
  for (unsigned int vp = _edge_bub; vp < _numbub + _edge_bub; ++vp)
  {
    bool try_again;
    unsigned int num_tries = 0;

    do
    {
      try_again = false;
      num_tries++;

      if (num_tries > _max_num_tries)
        mooseError("Too many tries of assigning void centers on grain boundaries in "
                   "PolycrystalHexVoidIC");

      Point rand_point;

      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
        rand_point(i) = _bottom_left(i) + _range(i) * MooseRandom::rand();

      _console << "random point" << rand_point << std::endl;
      // Allow the vectors to be sorted based on their distance from the
      // rand_point
      std::vector<PolycrystalHexVoidIC::DistancePoint> diff(_grain_num);

      for (unsigned int gr = 0; gr < _grain_num; ++gr)
      {
        diff[gr].d = _mesh.minPeriodicDistance(_var.number(), rand_point, _centerpoints[gr]);
        diff[gr].gr = gr;
      }

      std::sort(diff.begin(), diff.end(), _customLess);

      Point closest_point = _centerpoints[diff[0].gr];
      Point next_closest_point = _centerpoints[diff[1].gr];

      // If the closest_point or next_closest_point are across the periodic boundary
      // from rand_point, shift them to the equivalent position on the other
      // side of the domain
      if (rand_point(1) > _range(1) / 2 && closest_point(1) == 0.0)
        closest_point(1) += _range(1);
      if (rand_point(1) > _range(1) / 2 && next_closest_point(1) == 0.0)
        next_closest_point(1) += _range(1);
      if (rand_point(0) > _range(0) / 2 && closest_point(0) == 0.0)
        closest_point(0) += _range(0);
      if (rand_point(0) > _range(0) / 2 && next_closest_point(0) == 0.0)
        next_closest_point(0) += _range(0);

      _console << "closest point" << closest_point << std::endl;
      _console << "next closest point" << next_closest_point << std::endl;

      // Find Slope of Line in the plane orthogonal to the diff_centerpoint
      // vector
      Point diff_centerpoints =
          _mesh.minPeriodicVector(_var.number(), closest_point, next_closest_point);
      Point diff_rand_center = _mesh.minPeriodicVector(_var.number(), closest_point, rand_point);
      Point normal_vector = diff_centerpoints.cross(diff_rand_center);
      Point slope = normal_vector.cross(diff_centerpoints);

      _console << "diff centerpoints" << diff_centerpoints << std::endl;
      _console << "diff_rand_center" << diff_rand_center << std::endl;
      _console << "normal_vector" << normal_vector << std::endl;
      _console << "slope" << slope << std::endl;

      // Midpoint position vector between two center points
      Point midpoint = closest_point + (0.5 * diff_centerpoints);
      _console << "midpoint" << midpoint << std::endl;

      // Solve for the scalar multiplier solution on the line
      Real lambda = 0;
      Point mid_rand_vector = _mesh.minPeriodicVector(_var.number(), midpoint, rand_point);
      _console << "min_rand_vector" << mid_rand_vector << std::endl;

      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
        lambda += (mid_rand_vector(i) * slope(i)) /
                  (slope(0) * slope(0) + slope(1) * slope(1) + slope(2) * slope(2));

      // Assigning points to vector
      _centers[vp] = slope * lambda + midpoint;
      _console << "center" << _centers[vp] << std::endl;

      // Checking to see if points are in the domain ONLY WORKS FOR PERIODIC
      for (unsigned int i = 0; i < LIBMESH_DIM; i++)
        if ((_centers[vp](i) > _top_right(i)) || (_centers[vp](i) < _bottom_left(i)))
          try_again = true;

      for (unsigned int i = 0; i < vp; ++i)
      {
        Real dist = _mesh.minPeriodicDistance(_var.number(), _centers[vp], _centers[i]);

        if (dist < _bubspac)
          try_again = true;
      }

      // Two algorithms are available for screening bubbles falling in grain
      // interior. They produce
      // nearly identical results.
      // Here only one is listed. The other one is available upon request.

      // Use circle center for checking whether voids are at GBs
      if (try_again == false)
      {
        Real min_rij_1, min_rij_2, rij, rij_diff_tol;

        min_rij_1 = _range.norm();
        min_rij_2 = _range.norm();

        rij_diff_tol = 0.1 * _radius;

        for (unsigned int gr = 0; gr < _grain_num; ++gr)
        {
          rij = _mesh.minPeriodicDistance(_var.number(), _centers[vp], _centerpoints[gr]);

          if (rij < min_rij_1)
          {
            min_rij_2 = min_rij_1;
            min_rij_1 = rij;
          }
          else if (rij < min_rij_2)
            min_rij_2 = rij;
        }

        if (std::abs(min_rij_1 - min_rij_2) > rij_diff_tol)
          try_again = true;
      }

    } while (try_again == true);
  }
}

Real
PolycrystalHexVoidIC::value(const Point & p)
{
  Real value = 0.0;

  // Determine value for voids
  Real void_value = MultiSmoothCircleIC::value(p);

  // Determine value for grains
  Real grain_value = grainValueCalc(p);

  switch (_structure_type)
  {
    case 0:                 // assigning values for grains (order parameters)
      if (grain_value == 0) // Not in this grain
        value = grain_value;
      else                             // in this grain, but might be in a void
          if (void_value == _outvalue) // Not in a void
        value = grain_value;
      else if (void_value > _outvalue && void_value < _invalue) // On void interface
        value = 1.0 - (void_value - _outvalue) / (_invalue - _outvalue);
      else if (void_value == _invalue) // In a void, so op = 0
        value = 0.0;
      break;

    case 1: // assigning values for voids (concentration)
      value = void_value;
      break;
  }

  return value;
}

RealGradient
PolycrystalHexVoidIC::gradient(const Point & p)
{
  RealGradient gradient;
  RealGradient void_gradient = MultiSmoothCircleIC::gradient(p);

  // Order parameter assignment assumes zero gradient (sharp interface)
  switch (_structure_type)
  {
    case 1: // assigning gradient for voids
      gradient = void_gradient;
      break;
  }

  return gradient;
}

Real
PolycrystalHexVoidIC::grainValueCalc(const Point & p)
{
  Real val = 0.0;

  unsigned int min_index =
      PolycrystalICTools::assignPointToGrain(p, _centerpoints, _mesh, _var, _range.norm());

  // If the current order parameter index (_op_index) is equal to the min_index,
  // set the value to
  // 1.0
  if (_assigned_op[min_index] == _op_index)
    val = 1.0;

  if (val > 1.0)
    val = 1.0;

  if (val < 0.0)
    val = 0.0;

  return val;
}

void
PolycrystalHexVoidIC::computeGrainCenters()
{
  if (_op_num > _grain_num)
    mooseError("ERROR in PolycrystalHexVoidIC: Number of order parameters "
               "(op_num) can't be "
               "larger than the number of grains (grain_num)");

  // Initialize vectors
  _centerpoints.resize(_grain_num);
  _assigned_op.resize(_grain_num);

  // First grain
  _centerpoints[0](0) = _bottom_left(0) + _range(0) * 0.25;
  _centerpoints[0](1) = _bottom_left(1) + _range(1) * 0.0;
  _centerpoints[0](2) = _bottom_left(2) + _range(2) * 0.5;

  // Second grain
  _centerpoints[1](0) = _bottom_left(0) + _range(0) * 0.0;
  _centerpoints[1](1) = _bottom_left(1) + _range(1) * 0.5;
  _centerpoints[1](2) = _bottom_left(2) + _range(2) * 0.5;

  // Third grain
  _centerpoints[2](0) = _bottom_left(0) + _range(0) * 0.75;
  _centerpoints[2](1) = _bottom_left(1) + _range(1) * 0.0;
  _centerpoints[2](2) = _bottom_left(2) + _range(2) * 0.5;

  // Fourth grain
  _centerpoints[3](0) = _bottom_left(0) + _range(0) * 0.5;
  _centerpoints[3](1) = _bottom_left(1) + _range(1) * 0.5;
  _centerpoints[3](2) = _bottom_left(2) + _range(2) * 0.5;

  // Assign grains to specific order parameters in a way that maximizes the
  // distance
  _assigned_op = PolycrystalICTools::assignPointsToVariables(_centerpoints, _op_num, _mesh, _var);
}
