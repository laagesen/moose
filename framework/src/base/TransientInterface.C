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

#include "TransientInterface.h"
#include "FEProblem.h"
#include "TimePeriod.h"

template<>
InputParameters validParams<TransientInterface>()
{
  InputParameters params = emptyInputParameters();
  params.addParam<bool>("implicit", true, "Determines whether this object is calculated using an implicit or explicit form");

  params.addParamNamesToGroup("implicit", "Advanced");
  return params;
}


TransientInterface::TransientInterface(const InputParameters & parameters, const std::string & object_type) :
    _ti_feproblem(*parameters.get<FEProblem *>("_fe_problem")),
    _is_implicit(parameters.have_parameter<bool>("implicit") ? parameters.get<bool>("implicit") : true),
    _t(_is_implicit ? _ti_feproblem.time() : _ti_feproblem.timeOld()),
    _t_step(_ti_feproblem.timeStep()),
    _dt(_ti_feproblem.dt()),
    _dt_old(_ti_feproblem.dtOld()),
    _is_transient(_ti_feproblem.isTransient()),
    _object_type(object_type),
    _time_periods(_ti_feproblem.getTimePeriods()),
    _ti_name(MooseUtils::shortName(parameters.get<std::string>("name")))
{
}

TransientInterface::~TransientInterface()
{
}

bool
TransientInterface::isActive()
{
  // if we have zero or one time periods -> all objects are active all the time
  if (_time_periods.size() <= 1)
    return true;

  // look if _t lies in one of our time periods
  for (unsigned int i=1; i <= _time_periods.size(); ++i)  // Careful! We are purposely indexing one past the end of the array
    if ((i == _time_periods.size()) || // Are we in the last time period?
        ((_time_periods[i-1]->start() <= _t) && (_t < _time_periods[i]->start())))  // OR are we in one of the intermediate periods?
    {
      bool ret_value;
      const std::vector<std::string> & objects = _time_periods[i-1]->getObjectList(_object_type, ret_value);
      if (std::find(objects.begin(), objects.end(), _ti_name) != objects.end())
        return ret_value;
      else
        return !ret_value;
    }

  return false;
}
