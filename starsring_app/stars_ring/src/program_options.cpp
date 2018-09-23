#include <boost/program_options.hpp>
#include <iostream>

#include <stars_ring/program_options.hpp>

namespace {

void emit_help(std::ostream& s,
               const boost::program_options::options_description& desc) {
  s << "Program: stars-rings" << std::endl;
  s << desc << std::endl;
}
}  // namespace

ProgramOptions grep_program_options(int argc, char** argv) {
  ProgramOptions program_options;

  boost::program_options::options_description desc("Options");
  desc.add_options()
      // --help, -h:
      ("help,h", "Print help messages")
      // --n_cells,-c:
      ("n_cells,c",
       boost::program_options::value<unsigned>(&program_options.n_cells)
           ->default_value(4),
       "number of cells (half the number of sites).")
      // --n_max_stars,-s:
      ("n_max_stars,s",
       boost::program_options::value<unsigned>(&program_options.n_max_stars)
           ->default_value(2),
       "max number of excitations on one site.")
      // --perturbation_level,-p:
      ("perturbation_level,p",
       boost::program_options::value<unsigned>(
           &program_options.perturbation_level)
           ->default_value(3),
       "perturbation theory level used in basis generation.")
      // --multiplicity,-m:
      ("multiplicity,m",
       boost::program_options::value<unsigned>(&program_options.multiplicity)
           ->default_value(2),
       "spin multiplicity.")
      // --hamiltonian,-H:
      ("hamiltonian,H",
       boost::program_options::value<std::string>(
           &program_options.hamiltonian_type)
           ->default_value("spin"),
       "'spin' or 'osc' or 'jabcdz'.")
      // --print_basis:
      ("print_basis",
       boost::program_options::bool_switch(&program_options.print_basis)
           ->default_value(false),
       "print basis sets (localized basis, translation unique basis and "
       "ksubaspace basis).")
      // --omit_zero_star_space_calculations:
      ("omit_zero_star_space_calculations",
       boost::program_options::bool_switch(
           &program_options.omit_zero_star_space_calculations)
           ->default_value(false),
       "do not enter zero star subspace calculations block.")
      // --omit_one_star_space_calculations:
      ("omit_one_star_space_calculations",
       boost::program_options::bool_switch(
           &program_options.omit_one_star_space_calculations)
           ->default_value(false),
       "do not enter one star subspace calculations block.")
      // --n_extra_states_one_star_space:
      ("n_extra_states_one_star_space",
       boost::program_options::value<unsigned>(
           &program_options
                .n_extra_states_to_calculate_one_star_space_calculations)
           ->default_value(20),
       "Determine the number of calculate states.")
      // --hamiltonian_J,-J:
      ("hamiltonian_J,J",
       boost::program_options::value<double>(&program_options.hamiltonian_J)
           ->default_value(1.0),
       "Hamiltonian J (exchange, i.e. not spin) parameter value.")
      // --hamiltonian_Ez,-z:
      ("hamiltonian_Ez,z",
       boost::program_options::value<double>(&program_options.hamiltonian_Ez)
           ->default_value(0.0),
       "Hamiltonian Ez parameter value.")
      // --phi0,-0:
      ("phi0,0",
       boost::program_options::value<double>(&program_options.phi_0)
           ->default_value(0.0),
       "phi_0 in simple phis-establisher.")
      // --delta_phi,-d:
      ("delta_phi,d",
       boost::program_options::value<double>(&program_options.delta_phi)
           ->default_value(0.0),
       "delta_phi in simple phis-establisher.");
  boost::program_options::variables_map vm;
  try {
    boost::program_options::store(
        boost::program_options::command_line_parser(argc, argv)
            .options(desc)
            .run(),
        vm);  // may throw
    if (vm.count("help")) {
      emit_help(std::cout, desc);
      exit(0);
    }
    // sets auto variables (eq. class_specification_file_path),
    // throw is required variable is missing:
    boost::program_options::notify(vm);  // may throw
  } catch (boost::program_options::error& e) {
    std::cerr << "[GLOBAL ERROR] [PROGRAM OPTIONS ERROR]: " << e.what()
              << std::endl;
    std::cerr << std::endl;
    emit_help(std::cerr, desc);
    exit(1);
  }
  return program_options;
}