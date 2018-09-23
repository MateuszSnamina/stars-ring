#include<utility_kit/check_memory.hpp>

extern "C" {
#include <malloc.h>
}
#include<iostream>

namespace utility {

    void check_memory() {
        struct mallinfo malloc_info = mallinfo();
        std::cout << "memory: " << malloc_info.uordblks / 1048576.0 << std::endl;
    }

}