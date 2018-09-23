#include<utility_kit/ansi_escape_code.hpp>

namespace utility {
    extern const std::string ColorGray("\033[30m");
    extern const std::string ColorRed("\033[31m");
    extern const std::string ColorGreen("\033[32m");
    extern const std::string ColorYellow("\033[33m");
    extern const std::string ColorBlue("\033[34m");
    extern const std::string ColorMagenta("\033[35m");
    extern const std::string ColorCyan("\033[36m");
    extern const std::string ColorWhite("\033[37m");

    extern const std::string BgColorGray("\033[40m");
    extern const std::string BgColorRed("\033[41m");
    extern const std::string BgColorGreen("\033[42m");
    extern const std::string BgColorYellow("\033[43m");
    extern const std::string BgColorBlue("\033[44m");
    extern const std::string BgColorMagenta("\033[45m");
    extern const std::string BgColorCyan("\033[46m");
    extern const std::string BgColorWhite("\033[47m");

    extern const std::string ColorEnd("\033[0m");
}
