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

// MOOSE includes
#include "DOFMapOutput.h"
#include "FEProblem.h"
#include "KernelBase.h"
#include "MooseApp.h"
#include "Moose.h"
#include "Conversion.h"

// libMesh includes
#include "libmesh/fe.h"

// compiler includes (for type demangling)
#include <cxxabi.h>

template<>
InputParameters validParams<DOFMapOutput>()
{
  // Get the parameters from the base class
  InputParameters params = validParams<BasicOutput<FileOutput> >();

  // Screen and file output toggles
  params.addParam<bool>("output_screen", false, "Output to the screen");
  params.addParam<bool>("output_file", true, "Output to the file");
  params.addParam<std::string>("system_name", "nl0", "System to output");

  // By default this only executes on the initial timestep
  params.set<MultiMooseEnum>("execute_on") = "initial";

  return params;
}

DOFMapOutput::DOFMapOutput(const InputParameters & parameters) :
    BasicOutput<FileOutput>(parameters),
    _write_file(getParam<bool>("output_file")),
    _write_screen(getParam<bool>("output_screen")),
    _system_name(getParam<std::string>("system_name")),
    _mesh(_problem_ptr->mesh())
{
  // Set output coloring
  if (Moose::_color_console)
  {
    char * term_env = getenv("TERM");
    if (term_env)
    {
      std::string term(term_env);
      if (term != "xterm-256color" && term != "xterm")
        Moose::_color_console = false;
    }
  }
}

std::string
DOFMapOutput::filename()
{
  if (_file_num > 0)
    return _file_base + "_" + Moose::stringify(_file_num) + ".json";
  else
    return _file_base + ".json";
}

std::string
DOFMapOutput::demangle(const std::string & name)
{
#if defined(LIBMESH_HAVE_GCC_ABI_DEMANGLE)
  return libMesh::demangle(name.c_str());
#else
  // at least remove leading digits
  std::string demangled(name);
  while (demangled.length() && demangled[0] >= '0' && demangled[0] <= '9')
    demangled.erase(0,1);

  return demangled;
#endif
}

void
DOFMapOutput::writeStreamToFile(bool /*append*/)
{
  // Create the stream
  std::ofstream output;

  // Open the file and write contents of file output stream and close the file
  output.open(filename().c_str(), std::ios::trunc);
  output << _file_output_stream.str();
  output.close();

  // Clear the file output stream
  _file_output_stream.str("");
  _file_num++;
}

template<typename T>
std::string
DOFMapOutput::join(const T & begin, const T & end, const char* const delim)
{
  std::ostringstream os;
  for (T it = begin; it != end; ++it)
    os << (it != begin ? delim : "") << *it;
  return os.str();
}

void
DOFMapOutput::output(const ExecFlagType & /*type*/)
{
  // Don't build this information if nothing is to be written
  if (!_write_screen && !_write_file)
    return;

  std::stringstream oss;

  // Get the DOF Map through the equation system
  const System & sys = _problem_ptr->es().get_system(_system_name); // TransientNonlinearImplicit
  const DofMap & dof_map = sys.get_dof_map();

  // fetch the KernelWarehouse through the NonlinearSystem
  NonlinearSystem & nl = _problem_ptr->getNonlinearSystem();
  const KernelWarehouse & kernels = nl.getKernelWarehouse(0);

  // get a set of all subdomains
  const std::set<SubdomainID> & subdomains = _mesh.meshSubdomains();

  bool first = true;
  oss << "{\"ndof\": " << sys.n_dofs() << ", \"demangled\": ";
#if defined(LIBMESH_HAVE_GCC_ABI_DEMANGLE)
  oss << "true";
#else
  oss << "false";
#endif
  oss << ", \"vars\": [";
  for (unsigned int vg = 0; vg < dof_map.n_variable_groups(); ++vg)
  {
    const VariableGroup &vg_description (dof_map.variable_group(vg));
    for (unsigned int vn = 0; vn < vg_description.n_variables(); ++vn)
    {
      unsigned int var = vg_description.number(vn);

      if (!first)
        oss << ", ";
      first = false;

      oss << "{\"name\": \"" << vg_description.name(vn) << "\", \"subdomains\": [";
      for (std::set<SubdomainID>::const_iterator sd = subdomains.begin(); sd != subdomains.end(); ++sd)
      {
        oss << (sd != subdomains.begin() ? ", " : "") << "{\"id\": " << *sd << ", \"kernels\": [";

        // build a list of all kernels in the current subdomain
        nl.updateActiveKernels(*sd, 0);

        // if this variable has active kernels output them
        if (kernels.hasActiveKernels(var))
        {
          const std::vector<KernelBase *> & active_kernels = kernels.activeVar(var);
          for (unsigned i = 0; i<active_kernels.size(); ++i)
            oss << (i>0 ? ", " : "") << "{\"name\": \""<< active_kernels[i]->name() << "\", \"type\": \"" << demangle(typeid(*active_kernels[i]).name()) << "\"}";
        }
        oss << "], \"dofs\": [";

        // get the list of unique DOFs for this variable
        std::set<dof_id_type> dofs;
        for (unsigned int i = 0; i < _mesh.nElem(); ++i)
          if (_mesh.elem(i)->subdomain_id() == *sd)
          {
            std::vector<dof_id_type> di;
            dof_map.dof_indices(_mesh.elem(i), di, var);
            dofs.insert(di.begin(), di.end());
          }
        oss << join(dofs.begin(), dofs.end(), ", ") << "]}";
      }
      oss << "]}";
    }
  }
  oss << "]}\n";

  // Write the message to file stream
  if (_write_file)
    _file_output_stream << oss.str() << std::endl;

  // Write message to the screen
  if (_write_screen)
    _console << oss.str() << std::flush;

  // Write the actual file
  if (_write_file)
    writeStreamToFile();
}
