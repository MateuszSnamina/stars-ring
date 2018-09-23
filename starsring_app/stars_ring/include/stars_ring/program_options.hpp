#ifndef PROGRAM_OPTIONS_HPP
#define PROGRAM_OPTIONS_HPP

#include<string>

struct ProgramOptions {
    unsigned n_cells;
    unsigned n_max_stars;
    unsigned perturbation_level;
    unsigned multiplicity;
    std::string hamiltonian_type;
    bool print_basis;
    bool omit_zero_star_space_calculations;
    bool omit_one_star_space_calculations;
    unsigned n_extra_states_to_calculate_one_star_space_calculations;
    double hamiltonian_J;
    double hamiltonian_Ez;
    double phi_0;
    double delta_phi;
};

ProgramOptions grep_program_options(int argc, char** argv);

#endif
