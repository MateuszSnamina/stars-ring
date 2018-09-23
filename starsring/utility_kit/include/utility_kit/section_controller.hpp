#ifndef UTILITY_KIT_SECTION_CONTROLLER_HPP
#define UTILITY_KIT_SECTION_CONTROLLER_HPP

#include<armadillo>
#include<string>
#include<functional>
#include<iostream>
#include<iomanip>
#include"utility_kit/ansi_escape_code.hpp"

namespace utility {

    class SectionController {
    public:
        SectionController(std::function<void(std::string) > print_function);
        SectionController& new_section(std::string title);
        void measure_time();
        void time_raport();
        ~SectionController();
    private:
        const std::function<void(std::string) > m_print_function;
        std::string m_title;
        bool m_measure_time;
        arma::wall_clock m_section_timer;
    };

    inline SectionController make_cyan_section_controller() {
        auto print_function = [](std::string title) {
            std::cout << utility::BgColorCyan << std::left << std::setw(130) << "[SECTION]  " + title + ":"
                    << utility::ColorEnd << std::endl;
        };
        return SectionController(print_function);
    }

    inline SectionController make_green_section_controller() {
        auto print_function = [](std::string title) {
            std::cout << utility::BgColorGreen << std::left << std::setw(130) << "[SECTION]  " + title + ":"
                    << utility::ColorEnd << std::endl;
        };
        return SectionController(print_function);
    }

    inline SectionController make_magenta_section_controller() {
        auto print_function = [](std::string title) {
            std::cout << utility::BgColorMagenta << std::left << std::setw(130) << "[SECTION]  " + title + ":"
                    << utility::ColorEnd << std::endl;
        };
        return SectionController(print_function);
    }

}


#endif // UTILITY_KIT_SECTION_CONTROLLER_HPP
