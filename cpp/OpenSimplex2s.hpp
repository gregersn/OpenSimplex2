/**
 * K.jpg's OpenSimplex 2, smooth variant ("SuperSimplex")
 *
 * - 2D is standard simplex, modified to support larger kernels.
 *   Implemented using a lookup table.
 * - 3D is "Re-oriented 8-point BCC noise" which constructs an
 *   isomorphic BCC lattice in a much different way than usual.
 * - 4D uses a naïve pregenerated lookup table, and averages out
 *   to the expected performance.
 *
 * Multiple versions of each function are provided. See the
 * documentation above each, for more info.
 */

#ifndef OPENSIMPLEX2S_H
#define OPENSIMPLEX2S_H

#include <vector>

class Grad2 {
    public:
    double dx, dy;
    Grad2() : dx(0), dy(0) {} ;
    Grad2(double dx, double dy) {
        this->dx = dx;
        this->dy = dy;
    }
};

class Grad3 {
    double dx, dy, dz;
    public:
    Grad3(double dx, double dy, double dz) {
        this->dx = dx;
        this->dy = dy;
        this->dz = dz;
    }
};

class Grad4 {
    double dx, dy, dz, dw;
    public:
    Grad4(double dx, double dy, double dz, double dw) {
        this->dx = dx;
        this->dy = dy;
        this->dz = dz;
        this->dw = dw;
    }
};

class LatticePoint2D {
    public:
    int xsv, ysv;
    double dx, dy;
    LatticePoint2D() : xsv(0), ysv(0) {}
    LatticePoint2D(int xsv, int ysv) {
        this->xsv = xsv;
        this->ysv = ysv;
        double ssv = (xsv + ysv) * -0.211324865405187;
        this->dx = -xsv - ssv;
        this->dy = -ysv - ssv;
    }
};

class OpenSimplex2S {
    static const int PSIZE = 2048;
    static const int PMASK = 2047;

    std::vector<short> perm;
    std::vector<Grad2> permGrad2;
    std::vector<Grad3> permGrad3;
    std::vector<Grad4> permGrad4;

    static const std::vector<Grad2> GRADIENTS_2D;
    //static const std::vector<Grad3> GRADIENTS_3D;
    //static const std::vector<Grad4> GRADIENTS_4D;

    static std::vector<LatticePoint2D> LOOKUP_2D;

    static int fastFloor(double x);

    static void init();

    double noise2_Base(double xs, double ys);

    public:
    OpenSimplex2S(long seed);
    double noise2(double x, double y);
    
};
#endif
