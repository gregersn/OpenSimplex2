#include <iostream>
#include "OpenSimplex2s.hpp"

int main(int argc, char **argv) {
    OpenSimplex2S noise = OpenSimplex2S(1234);
    for(int i = 0; i < 100; i++) {
        std::cout << noise.noise2(i, 0) << "\n";
    }
}
