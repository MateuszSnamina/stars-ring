#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include <armadillo>
#include <boost/optional/optional_io.hpp>

#include <maths_in_physic/oscilator_realm.hpp>
#include <maths_in_physic/spin_realm.hpp>
#include <stars_ring/program_options.hpp>
#include <stars_ring_analytical/analytical_formulas_box_af_oscilators.hpp>
#include <stars_ring_analytical/analytical_formulas_box_af_spins.hpp>
#include <stars_ring_analytical/analytical_formulas_box_jabcdzx.hpp>
#include <stars_ring_analytical/standard_calculator.hpp>
#include <stars_ring_basis/basis_print.hpp>
#include <stars_ring_basis/basis_toolbox.hpp>
#include <stars_ring_basis/basis_typedefs.hpp>
#include <stars_ring_basis/raw_state_coupled_elements_generator.hpp>
#include <stars_ring_basis/raw_state_operations.hpp>
#include <stars_ring_core/physical_system.hpp>
#include <stars_ring_numerical/hamiltonian_creator_jabcdzx.hpp>
#include <stars_ring_numerical/hamiltonian_creator_oscylator.hpp>
#include <stars_ring_numerical/hamiltonian_creator_spin.hpp>
#include <stars_ring_numerical/morphology_alalizer.hpp>
#include <stars_ring_numerical/phis_establisher.hpp>
#include <stars_ring_numerical/standard_calculator.hpp>
#include <utility_kit/print_utility.hpp>
#include <utility_kit/section_controller.hpp>

/*
std::complex<double> determine_coupling(
        const TranslationUniqueElement& lhs,
        const TranslationUniqueElement& rhs,
        double k) {
    std::complex<double> result = 0.0;
    assert(lhs.representative_state.size() == rhs.representative_state.size());
    unsigned size = lhs.representative_state.size();
    assert(size % 2 == 0);
    std::vector<unsigned> temp = rhs.representative_state;
    for (unsigned i = 0; i < size; i = i + 2) {
        std::rotate(temp.begin(), temp.begin() + 2, temp.end());
        result += std::exp(std::complex<double>(0, 1) * k * double(i)) *
determine_coupling(lhs.representative_state, temp);
    }
    double lhs_extra_norm_coef = 1 / std::sqrt(lhs.representative_state.size() /
2 / lhs.cycle); double rhs_extra_norm_coef = 1 /
std::sqrt(rhs.representative_state.size() / 2 / rhs.cycle); return result *
lhs_extra_norm_coef * rhs_extra_norm_coef;
}
 */

void print_basis_content(const stars_ring_basis::BasisBox& basis_box) {
  assert(basis_box.physical_system());
  std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system =
      basis_box.physical_system();
  std::cout << "[DATA   ] localized-basis" << std::endl;
  print_basis(basis_box.localized_basis());
  std::cout << "[DATA   ] translation-unique-basis" << std::endl;
  print_basis(basis_box.translation_unique_basis());
  for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
    std::cout << "[DATA   ] ksubspace-basis" << std::endl;
    print_basis(basis_box.ksubspace_basis(nk));
  }
}

void print_basis_sizes(const stars_ring_basis::BasisBox& basis_box) {
  assert(basis_box.physical_system());
  std::cout << utility::ValuePut<unsigned>("[INFO   ] localized-basis-set size",
                                           50,
                                           basis_box.localized_basis().size());
  std::cout << utility::ValuePut<unsigned>(
      "[INFO   ] translation-unique-basis size", 50,
      basis_box.translation_unique_basis().size());
  for (unsigned nk = 0; nk < basis_box.physical_system()->n_cells(); nk++)
    std::cout << utility::ValuePut<unsigned>(
        "[INFO   ] ksubspace-basis size (nk = " + std::to_string(nk) + ")", 50,
        basis_box.ksubspace_basis(nk).size());
  //    std::cout << utility::ValuePut<unsigned>(
  //            "[INFO   ] accumulated ksubspace-basis size", 50,
  //            std::accumulate(
  //            basis_box.ksubspace_basis().begin(),
  //            basis_box.ksubspace_basis().end(), 0u,
  //            [] (unsigned temp, const stars_ring_basis::KsubspaceBasis &
  //            k_subspace_basis) {
  //                return temp + k_subspace_basis.size();
  //            }));
}

void print_ground_state_basic_information(
    const stars_ring_numerical::StandardCalculator& calculator) {
  const double quantum_energy = calculator.energies(0);
  const double classical_energy = calculator.hamiltonian_sparse()(0, 0);
  const double correlation_energy = quantum_energy - classical_energy;
  std::cout << utility::ValuePut<double>(
      "[DATA   ] [NUMERICAL] ground state energy", 55, quantum_energy);
  std::cout << utility::ValuePut<double>(
      "[DATA   ] [NUMERICAL] ground state classical energy", 55,
      classical_energy);
  std::cout << utility::ValuePut<double>(
      "[DATA   ] [NUMERICAL] ground state correlation energy", 55,
      correlation_energy);
}

void print_subspace_indices(stars_ring_numerical::DegeneracySubspacesAnalizer&
                                degeneracy_subspace_analizer) {
  auto f = std::cout.flags();
  std::cout << std::fixed;
  std::cout << std::right;

  for (unsigned n_subspace = 0;
       n_subspace < degeneracy_subspace_analizer.n_subspaces(); n_subspace++) {
    std::cout << "[DATA   ] [NUMERICAL] ";
    std::cout << "n_subspace, [from, to), size: ";
    std::cout << std::setw(3) << n_subspace << ", ";
    std::cout << "[" << std::setw(3)
              << degeneracy_subspace_analizer.subspace_indices(n_subspace).first
              << ", ";
    std::cout
        << std::setw(3)
        << degeneracy_subspace_analizer.subspace_indices(n_subspace).second
        << ")"
        << ", ";
    std::cout
        << std::setw(3)
        << degeneracy_subspace_analizer.subspace_indices(n_subspace).second -
               degeneracy_subspace_analizer.subspace_indices(n_subspace).first
        << std::endl;
  }
  std::cout.flags(f);
}

void print_original_states_content(
    stars_ring_numerical::DegeneracySubspacesAnalizer&
        degeneracy_subspace_analizer,
    double energy_ground_state_zero_stars_subspace) {
  auto f = std::cout.flags();
  std::cout << std::fixed;
  std::cout << std::setprecision(10);
  std::cout << std::showpos;
  std::cout << std::right;
  for (unsigned n_state = 0;
       n_state < degeneracy_subspace_analizer.n_original_states(); n_state++) {
    std::cout << "[DATA   ] [NUMERICAL] ";
    std::cout << "n_state, energy, excitation_energy, (content...): ";
    std::cout << std::setw(5) << n_state << ", ";
    std::cout << std::setw(15)
              << degeneracy_subspace_analizer.original_state_energy(n_state)
              << ", ";
    std::cout << std::setw(15)
              << degeneracy_subspace_analizer.original_state_energy(n_state) -
                     energy_ground_state_zero_stars_subspace
              << ", ";
    std::cout << "(";
    for (double content :
         degeneracy_subspace_analizer.content_original_state(n_state)) {
      std::cout << std::setw(12) << content << ", ";
    }
    std::cout << ")";
    std::cout << std::endl;
  }
  std::cout.flags(f);
}

void print_subspace_content(stars_ring_numerical::DegeneracySubspacesAnalizer&
                                degeneracy_subspace_analizer,
                            double energy_ground_state_zero_stars_subspace) {
  auto f = std::cout.flags();
  std::cout << std::fixed;
  std::cout << std::setprecision(10);
  std::cout << std::showpos;
  std::cout << std::right;
  for (unsigned n_subspace = 0;
       n_subspace < degeneracy_subspace_analizer.n_subspaces(); n_subspace++) {
    std::cout << "[DATA   ] [NUMERICAL] ";
    std::cout << "n_subspace, energy, excitation_energy, (content...): ";
    std::cout << std::setw(3) << n_subspace << ", ";
    std::cout << std::setw(15)
              << degeneracy_subspace_analizer.subspace_energy(n_subspace)
              << ", ";
    std::cout << std::setw(15)
              << degeneracy_subspace_analizer.subspace_energy(n_subspace) -
                     energy_ground_state_zero_stars_subspace
              << ", ";
    std::cout << "(";
    for (double content :
         degeneracy_subspace_analizer.content_subspace(n_subspace)) {
      std::cout << std::setw(15) << content << ", ";
    }
    std::cout << ")";
    std::cout << std::endl;
  }
  std::cout.flags(f);
}

void print_subspace_filtered_content_version1(
    stars_ring_numerical::DegeneracySubspacesAnalizer&
        degeneracy_subspace_analizer,
    double energy_ground_state_zero_stars_subspace) {
  auto f = std::cout.flags();
  std::cout << std::fixed;
  std::cout << std::setprecision(10);
  std::cout << std::showpos;
  std::cout << std::right;
  for (unsigned n_subspace = 0;
       n_subspace < degeneracy_subspace_analizer.n_subspaces(); n_subspace++) {
    std::cout << "[DATA   ] [NUMERICAL] ";
    std::cout << "n_subspace, energy, excitation_energy, (content...): ";
    std::cout << std::setw(3) << n_subspace << ", ";
    std::cout << std::setw(15)
              << degeneracy_subspace_analizer.subspace_energy(n_subspace)
              << ", ";
    std::cout << std::setw(15)
              << degeneracy_subspace_analizer.subspace_energy(n_subspace) -
                     energy_ground_state_zero_stars_subspace
              << ", ";
    std::cout << "(";
    if (degeneracy_subspace_analizer.subspace_validity(n_subspace)) {
      for (unsigned nk :
           degeneracy_subspace_analizer.fileted_content_subspace(n_subspace))
        std::cout << std::setw(3) << nk << ", ";
    } else {
      std::cout << "not-valid-content";
    }
    std::cout << ")";
    std::cout << std::endl;
  }
  std::cout.flags(f);
}

void print_subspace_filtered_content_version2(
    stars_ring_numerical::DegeneracySubspacesAnalizer&
        degeneracy_subspace_analizer,
    double energy_ground_state_zero_stars_subspace) {
  auto f = std::cout.flags();
  std::cout << std::fixed;
  std::cout << std::setprecision(10);
  std::cout << std::showpos;
  std::cout << std::right;
  for (unsigned n_subspace = 0;
       n_subspace < degeneracy_subspace_analizer.n_subspaces(); n_subspace++) {
    if (degeneracy_subspace_analizer.subspace_validity(n_subspace)) {
      for (unsigned nk :
           degeneracy_subspace_analizer.fileted_content_subspace(n_subspace)) {
        std::cout << "[DATA   ] [NUMERICAL] ";
        std::cout << "n_subspace, nk, energy, excitation_energy: ";
        std::cout << std::setw(3) << n_subspace << ", ";
        std::cout << std::setw(3) << nk << ", ";
        std::cout << std::setw(15)
                  << degeneracy_subspace_analizer.subspace_energy(n_subspace)
                  << ", ";
        std::cout << std::setw(15)
                  << degeneracy_subspace_analizer.subspace_energy(n_subspace) -
                         energy_ground_state_zero_stars_subspace
                  << std::endl;
      }
    } else {
      for (unsigned idx =
               degeneracy_subspace_analizer.subspace_indices(n_subspace).first;
           idx <
           degeneracy_subspace_analizer.subspace_indices(n_subspace).second;
           ++idx) {
        std::cout << "[DATA   ] [NUMERICAL]";
        std::cout << "n_subspace, nk, energy, excitation_energy: ";
        std::cout << std::setw(3) << n_subspace << ", ";
        std::cout << std::setw(3) << "---"
                  << ", ";
        std::cout << std::setw(15)
                  << degeneracy_subspace_analizer.subspace_energy(n_subspace)
                  << ", ";
        std::cout << std::setw(15)
                  << degeneracy_subspace_analizer.subspace_energy(n_subspace) -
                         energy_ground_state_zero_stars_subspace
                  << std::endl;
      }
    }
  }
  std::cout.flags(f);
}

int main(int argc, char** argv) {
  utility::SectionController section_controller =
      utility::make_cyan_section_controller();

  // #######################################################################
  // ######################    Intro    ####################################
  // #######################################################################
  section_controller.new_section("Intro");
  std::cout << "[INFO   ] Program: stars-ring. Version: 1.1. Compilation date: "
            << __DATE__ << "." << std::endl;

  // #######################################################################
  // ######################     Program options    #########################
  // #######################################################################
  section_controller.new_section("Program options");

  const ProgramOptions program_options = grep_program_options(argc, argv);
  //    std::cout << utility::ValuePut<double>("[DATA   ] [PROGRAM_OPTIONS]
  //    n_cells", 50, program_options.n_cells); std::cout <<
  //    utility::ValuePut<double>("[DATA   ] [PROGRAM_OPTIONS] n_max_stars", 50,
  //    program_options.n_max_stars); std::cout <<
  //    utility::ValuePut<double>("[DATA   ] [PROGRAM_OPTIONS]
  //    perturbation_level", 50, program_options.perturbation_level); std::cout
  //    << utility::ValuePut<double>("[DATA   ] [PROGRAM_OPTIONS] multiplicity",
  //    50, program_options.multiplicity); std::cout <<
  //    utility::ValuePut<double>("[DATA   ] [PROGRAM_OPTIONS]
  //    hamiltonian_type", 50, program_options.hamiltonian_type);
  std::cout << "[DATA   ] [PROGRAM_OPTIONS] n_cells:            "
            << program_options.n_cells << std::endl;
  std::cout << "[DATA   ] [PROGRAM_OPTIONS] n_max_stars:        "
            << program_options.n_max_stars << std::endl;
  std::cout << "[DATA   ] [PROGRAM_OPTIONS] perturbation_level: "
            << program_options.perturbation_level << std::endl;
  std::cout << "[DATA   ] [PROGRAM_OPTIONS] multiplicity:       "
            << program_options.multiplicity << std::endl;
  std::cout << "[DATA   ] [PROGRAM_OPTIONS] hamiltonian_type:   "
            << program_options.hamiltonian_type << std::endl;

  // #######################################################################
  // ######################    Physical system    ##########################
  // #######################################################################
  section_controller.new_section("Physical system");

  std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system =
      stars_ring_core::PhysicalSystem::Builder()
          .n_cells(program_options.n_cells)
          .max_n_stars(program_options.n_max_stars)
          .build();

  std::cout << "[DATA   ] physical_system: " << *physical_system << std::endl;

  // #######################################################################
  // ######################    Interactions params    ######################
  // #######################################################################
  section_controller.new_section("Interactions parameters");

  std::shared_ptr<stars_ring_analytical::AnalyticalFormulasBox>
      analytical_formulas_box;
  std::shared_ptr<stars_ring_basis::CoupledRawStatesGenerator>
      coupled_raw_states_generator;
  std::shared_ptr<stars_ring_numerical::HamiltonianCreator> hamiltonian_creator;

  if (program_options.hamiltonian_type == "spin") {
    std::cout << "[INFO   ] Hamiltonian type -- spins" << std::endl;
    if (program_options.multiplicity < program_options.n_max_stars) {
      std::cerr << "[ERROR] The condition (multiplicity < n_max_stars) is "
                   "unphysical and will cause the program fatal error!"
                << std::endl;
      exit(102);
    }
    const maths_in_physic::SpinRealm spin_realm(program_options.multiplicity);
    // std::cout << "[DATA   ] spin_realm: " << std::endl
    //        << spin_realm << std::endl;
    analytical_formulas_box =
        std::make_shared<stars_ring_analytical::AnalyticalFormulasBoxAfSpins>(
            physical_system, program_options.multiplicity);
    coupled_raw_states_generator =
        std::make_shared<stars_ring_basis::CoupledRawStatesGenerator_AF>(
            program_options.n_max_stars);
    hamiltonian_creator =
        std::make_shared<stars_ring_numerical::SpinHamiltonianCreator>(
            spin_realm);
  } else if (program_options.hamiltonian_type == "osc") {
    std::cout << "[INFO   ] Hamiltonian type -- oscilators" << std::endl;
    const maths_in_physic::OscilatorRealm oscilator_realm(
        program_options.n_max_stars);
    // std::cout << "[DATA   ] oscilator_realm: " << std::endl
    //        << oscilator_realm << std::endl;
    analytical_formulas_box = std::make_shared<
        stars_ring_analytical::AnalyticalFormulasBoxAfOscylators>(
        physical_system);
    coupled_raw_states_generator =
        std::make_shared<stars_ring_basis::CoupledRawStatesGenerator_AF>(
            program_options.n_max_stars);
    hamiltonian_creator =
        std::make_shared<stars_ring_numerical::OscylatorHamiltonianCreator>(
            oscilator_realm);
  } else if (program_options.hamiltonian_type == "jabcdzx") {
    std::cout << "[INFO   ] Hamiltonian type -- jabcdzx" << std::endl;
    std::cout << "[DATA   ] [PROGRAM_OPTIONS] hamiltonian::J:     "
              << program_options.hamiltonian_J << std::endl;
    std::cout << "[DATA   ] [PROGRAM_OPTIONS] hamiltonian::Ez:    "
              << program_options.hamiltonian_Ez << std::endl;
    std::cout << "[DATA   ] [PROGRAM_OPTIONS] hamiltonian::Ex:    "
              << program_options.hamiltonian_Ex << std::endl;
    std::cout << "[DATA   ] [PROGRAM_OPTIONS] phis_establisher::phi0:       "
              << program_options.phi_0 << std::endl;
    std::cout << "[DATA   ] [PROGRAM_OPTIONS] phis_establisher::delta_phi:  "
              << program_options.delta_phi << std::endl;
    if (program_options.multiplicity < program_options.n_max_stars) {
      std::cerr << "[ERROR] The condition (multiplicity < n_max_stars) is "
                   "unphysical and will cause the program fatal error!"
                << std::endl;
      exit(102);
    }
    const maths_in_physic::SpinRealm spin_realm(program_options.multiplicity);
    // const double A = std::pow(spin_realm.S, 2);
    const double A = 0.0;  // debug!!
    const double B = +1.0;
    // const double C = 1.0 / 4.0;
    const double C = 0.0;  // debug!!
    const double D = +1.0;
    analytical_formulas_box =
        std::make_shared<stars_ring_analytical::AnalyticalFormulasBoxJABCDZX>(
            physical_system, program_options.multiplicity, A, B, C, D,
            program_options.hamiltonian_J, program_options.hamiltonian_Ez,
            program_options.hamiltonian_Ex, program_options.phi_0);
    coupled_raw_states_generator =
        std::make_shared<stars_ring_basis::CoupledRawStatesGenerator_AF>(
            program_options.n_max_stars);
    // std::shared_ptr<stars_ring_numerical::TrivialPhisEstablisher>
    // phis_establisher =
    //        std::make_shared<stars_ring_numerical::TrivialPhisEstablisher>(
    //        program_options.phi_0);
    std::shared_ptr<stars_ring_numerical::SimplePhisEstablisher>
        phis_establisher =
            std::make_shared<stars_ring_numerical::SimplePhisEstablisher>(
                program_options.phi_0, program_options.delta_phi);
    hamiltonian_creator =
        std::make_shared<stars_ring_numerical::JABCDZXHamiltonianCreator>(
            spin_realm, A, B, C, D, program_options.hamiltonian_J,
            program_options.hamiltonian_Ez, program_options.hamiltonian_Ex,
            phis_establisher);
  } else {
    std::cerr << "[ERROR] Invalid hamiltonian type program option."
              << std::endl;
    exit(101);
  }

  // #######################################################################
  // ######################    Theory    ###################################
  // #######################################################################
  {
    section_controller.new_section("ANALITYCAL-THEORY");
    // ---------------------------------------------------------------------
    if (program_options.hamiltonian_type == "jabcdzx") {
      auto analytical_formulas_box_jabcdzx = std::dynamic_pointer_cast<
          stars_ring_analytical::AnalyticalFormulasBoxJABCDZX>(
          analytical_formulas_box);
      std::cout << utility::ValuePut<double>(
          "[DATA   ] [ANALITYCAL-THEORY] [jabcdzx] theta", 70,
          analytical_formulas_box_jabcdzx->theta());
      std::cout << utility::ValuePut<double>(
          "[DATA   ] [ANALITYCAL-THEORY] [jabcdzx] mean_orbital_operator", 70,
          analytical_formulas_box_jabcdzx->mean_orbital_operator());
      std::cout << utility::ValuePut<double>(
          "[DATA   ] [ANALITYCAL-THEORY] [jabcdzx] J_spin", 70,
          analytical_formulas_box_jabcdzx->J_spin());

      //   std::cout << "[DATA   ] [ANALITYCAL-THEORY] theta_opt (nell) : "
      //             << analytical_formulas_box_jabcdzx->theta_opt_nell()
      //             << std::endl;
      //   std::cout << "[DATA   ] [ANALITYCAL-THEORY] theta_opt (corrected
      //   nell) : "
      //             <<
      //             analytical_formulas_box_jabcdzx->theta_opt_corrected_nell()
      //             << std::endl;
    }
    // ---------------------------------------------------------------------
    stars_ring_analytical::StandardCalculator standard_calculator(
        analytical_formulas_box);
    standard_calculator.calculate();
    std::cout << utility::ValuePut<double>(
        "[DATA   ] [ANALITYCAL-THEORY] ground state absolute classical energy",
        70, standard_calculator.ground_state_classical_energy());
    std::cout << utility::ValuePut<double>(
        "[DATA   ] [ANALITYCAL-THEORY] ground state correlation energy", 70,
        standard_calculator.ground_state_correlation_energy());
    std::cout << utility::ValuePut<double>(
        "[DATA   ] [ANALITYCAL-THEORY] ground state absolute energy", 70,
        standard_calculator.ground_state_energy());
    for (unsigned i = 0; i < standard_calculator.n_exc_states(); ++i)
      std::cout << utility::ValuePut<double>(
          "[DATA   ] [ANALITYCAL-THEORY] excitation energy (#" +
              std::to_string(i) + ")",
          70, standard_calculator.exc_state_relative_energy(i));
  }

  // #######################################################################
  // ################    Calculations in <n_A> - <n_B> = 0 subspace    #####
  // #######################################################################
  double energy_ground_state_zero_stars_subspace = arma::datum::nan;
  if (!program_options.omit_zero_star_space_calculations) {
    section_controller
        .new_section(
            "n_A - n_B = 0 subspace -- numerical exact diagonalization")
        .measure_time();
    // ---------------------------------------------------------------------
    const stars_ring_basis::RawState init_state =
        stars_ring_basis::operations::make_ground_state(
            physical_system->n_sites());
    stars_ring_numerical::StandardCalculator calculator(
        physical_system, {init_state}, coupled_raw_states_generator,
        program_options.perturbation_level, hamiltonian_creator);
    print_basis_sizes(calculator.basis_box());
    if (program_options.print_basis)
      print_basis_content(calculator.basis_box());
    // ---------------------------------------------------------------------
    const unsigned requested_states_to_calculate = 10;
    calculator.do_calculations_sparse(requested_states_to_calculate);
    if (program_options.print_hamiltonian) {
      std::cout << "[DATA   ] [NUMERICAL] hamiltonian: " << std::endl;
      std::cout << arma::mat(calculator.hamiltonian_sparse()) << std::endl;
    }
    print_ground_state_basic_information(calculator);
    // ---------------------------------------------------------------------
    stars_ring_numerical::MorphologyAlalizer morphology_analizator(
        calculator.basis_box());
    for (unsigned i = 0; i < calculator.n_calculated_states(); i++) {
      const double energy = calculator.energies(i);
      const arma::vec state = calculator.states(i);
      morphology_analizator.calculate_morphology(state);
      std::cout << "[DATA   ] [NUMERICAL] state nk, n_stars, n_stars_A, "
                   "n_stars_B, energy: "
                << std::setw(5) << morphology_analizator.main_nk() << ", "
                << std::setw(10) << morphology_analizator.average_n_stars()
                << ", " << std::setw(10)
                << morphology_analizator.average_n_stars_A() << ", "
                << std::setw(10) << morphology_analizator.average_n_stars_B()
                << ", " << std::setw(10) << energy << std::endl;
    }
    energy_ground_state_zero_stars_subspace = calculator.energies(0);
    // ---------------------------------------------------------------------
  }

  // #######################################################################
  // #####    Calculations in <n_A> - <n_B> = 1 subspace    ################
  // #######################################################################
  if (!program_options.omit_one_star_space_calculations) {
    section_controller
        .new_section(
            "n_A - n_B = 1 subspace -- numerical exact diagonalization")
        .measure_time();
    // ---------------------------------------------------------------------
    const stars_ring_basis::RawState init_state =
        stars_ring_basis::operations::make_one_star_state(
            physical_system->n_sites(), 0);
    stars_ring_numerical::StandardCalculator calculator(
        physical_system, {init_state}, coupled_raw_states_generator,
        program_options.perturbation_level, hamiltonian_creator);
    print_basis_sizes(calculator.basis_box());
    if (program_options.print_basis)
      print_basis_content(calculator.basis_box());
    // ---------------------------------------------------------------------
    const unsigned requested_states_to_calculate =
        physical_system->n_cells() +
        program_options.n_extra_states_to_calculate_one_star_space_calculations;
    calculator.do_calculations_sparse(requested_states_to_calculate);
    if (program_options.print_hamiltonian) {
      std::cout << "[DATA   ] [NUMERICAL] hamiltonian: " << std::endl;
      std::cout << arma::mat(calculator.hamiltonian_sparse()) << std::endl;
    }
    print_ground_state_basic_information(calculator);
    for (unsigned i = 0; i < calculator.n_calculated_states(); i++) {
      const double energy = calculator.energies(i);
      std::cout << "[DATA   ] [NUMERICAL] energy, excitation_energy: "
                << std::setw(10) << energy << ", " << std::setw(10)
                << energy - energy_ground_state_zero_stars_subspace
                << std::endl;
    }
    // ---------------------------------------------------------------------
    stars_ring_numerical::MorphologyAlalizer morphology_analizator(
        calculator.basis_box());
    for (unsigned i = 0; i < calculator.n_calculated_states(); i++) {
      const double energy = calculator.energies(i);
      const arma::vec state = calculator.states(i);
      morphology_analizator.calculate_morphology(state);
      std::cout << "[DATA   ] [NUMERICAL] state nk, n_stars, n_stars_A, "
                   "n_stars_B, energy, "
                   "excitation_energy: "
                << std::setw(5) << morphology_analizator.main_nk() << ", "
                << std::setw(10) << morphology_analizator.average_n_stars()
                << ", " << std::setw(10)
                << morphology_analizator.average_n_stars_A() << ", "
                << std::setw(10) << morphology_analizator.average_n_stars_B()
                << ", " << std::setw(10) << energy << ", " << std::setw(10)
                << energy - energy_ground_state_zero_stars_subspace
                << std::endl;
    }
    // ---------------------------------------------------------------------
    stars_ring_numerical::DegeneracySubspacesAnalizer
        degeneracy_subspace_analizer(calculator.energies(), calculator.states(),
                                     calculator.basis_box());
    degeneracy_subspace_analizer.determine_subspaces_indices();
    print_subspace_indices(degeneracy_subspace_analizer);
    degeneracy_subspace_analizer.determine_original_states_content();
    print_original_states_content(degeneracy_subspace_analizer,
                                  energy_ground_state_zero_stars_subspace);
    degeneracy_subspace_analizer.determine_subspaces_content();
    // print_subspace_content(degeneracy_subspace_analizer,
    //                       energy_ground_state_zero_stars_subspace);
    degeneracy_subspace_analizer.determine_subspaces_validity();
    degeneracy_subspace_analizer.determine_filtered_subspaces_content();
    // print_subspace_filtered_content_version1(
    //    degeneracy_subspace_analizer,
    //    energy_ground_state_zero_stars_subspace);
    print_subspace_filtered_content_version2(
        degeneracy_subspace_analizer, energy_ground_state_zero_stars_subspace);
    // ---------------------------------------------------------------------
  }
}
