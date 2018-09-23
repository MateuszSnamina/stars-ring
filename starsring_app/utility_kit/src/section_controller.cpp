#include"utility_kit/section_controller.hpp"

namespace utility {

    SectionController::SectionController(std::function<void(std::string) > print_function) :
    m_print_function(print_function),
    m_measure_time(false) {
    }

    SectionController& SectionController::new_section(std::string title) {
        if (m_measure_time) time_raport();
        m_title = title;
        m_print_function(title);
        m_measure_time = false;
        return *this;
    }

    void SectionController::measure_time() {
        m_measure_time = true;
        m_section_timer.tic();
    }

    void SectionController::time_raport() {
        double n_secs = m_section_timer.toc();
        std::cout << "[INFO   ] [TIME] section '" << m_title << "' took " << n_secs << " seconds." << std::endl;
        m_measure_time = false;
    }

    SectionController::~SectionController() {
        if (m_measure_time) time_raport();
    }

} // end of namespace utility
