#include "OpenSimplex2s.hpp"

void OpenSimplex2S::init() {
    for(int i = 0; i < 8; i++) {
        int i1, j1, i2, j2;

        if ((i & 1) == 0) {
            if ((i & 2) == 0) { i1 = -1; j1 = 0; } else { i1 = 1; j1 = 0; }
            if ((i & 4) == 0) { i2 = 0; j2 = -1; } else { i2 = 0; j2 = 1; }
        } else {
            if ((i & 2) != 0) { i1 = 2; j1 = 1; } else { i1 = 0; j1 = 1; }
            if ((i & 4) != 0) { i2 = 1; j2 = 2; } else { i2 = 1; j2 = 0; }
        }

        LOOKUP_2D[i * 4 + 0] = LatticePoint2D(0, 0);
        LOOKUP_2D[i * 4 + 1] = LatticePoint2D(1, 1);
        LOOKUP_2D[i * 4 + 2] = LatticePoint2D(i1, j1);
        LOOKUP_2D[i * 4 + 3] = LatticePoint2D(i2, j2);
    }
}

OpenSimplex2S::OpenSimplex2S(long seed) : perm(PSIZE), permGrad2(PSIZE) {
    this->init();
    std::vector<short> source = std::vector<short>(PSIZE);
    for(short i = 0; i < PSIZE; i++) {
        source[i] = i;
    }

    for(int i = PSIZE - 1; i >= 0; i--) {
        seed = seed * 6364136223846793005L + 1442695040888963407L;
        int r = (int)((seed + 31) % (i + 1));
        if (r < 0)
            r += (i + 1);
        perm[i] = source[r];
        permGrad2[i] = GRADIENTS_2D[perm[i]];
        //permGrad3.push_back(GRADIENTS_3D[perm[i]]);
        //permGrad4.push_back(GRADIENTS_4D[perm[i]]);
        source[r] = source[i];
    }
}

/*
    * Utility
    */

int OpenSimplex2S::fastFloor(double x) {
    int xi = (int)x;
    return x < xi ? xi - 1 : xi;
}

/**
 * 2D SuperSimplex noise, standard lattice orientation.
 */
double OpenSimplex2S::noise2(double x, double y) {
		// Get points for A2* lattice
		double s = 0.366025403784439 * (x + y);
		double xs = x + s, ys = y + s;
		
		return noise2_Base(xs, ys);
}

/**
 * 2D SuperSimplex noise base.
 * Lookup table implementation inspired by DigitalShadow.
 */
double OpenSimplex2S::noise2_Base(double xs, double ys) {
    double value = 0;
    
    // Get base points and offsets
    int xsb = fastFloor(xs), ysb = fastFloor(ys);
    double xsi = xs - xsb, ysi = ys - ysb;
    
    // Index to point list
    int a = (int)(xsi + ysi);
    int index =
        (a << 2) |
        (int)(xsi - ysi / 2 + 1 - a / 2.0) << 3 |
        (int)(ysi - xsi / 2 + 1 - a / 2.0) << 4;
    
    double ssi = (xsi + ysi) * -0.211324865405187;
    double xi = xsi + ssi, yi = ysi + ssi;

    // Point contributions
    for (int i = 0; i < 4; i++) {
        LatticePoint2D c = LOOKUP_2D[index + i];

        double dx = xi + c.dx, dy = yi + c.dy;
        double attn = 2.0 / 3.0 - dx * dx - dy * dy;
        if (attn <= 0) continue;

        int pxm = (xsb + c.xsv) & PMASK, pym = (ysb + c.ysv) & PMASK;
        Grad2 grad = permGrad2[perm[pxm] ^ pym];
        double extrapolation = grad.dx * dx + grad.dy * dy;
        
        attn *= attn;
        value += attn * attn * extrapolation;
    }
    
    return value;
}

std::vector<LatticePoint2D> OpenSimplex2S::LOOKUP_2D = std::vector<LatticePoint2D>(8 * 4);



const std::vector<Grad2> OpenSimplex2S::GRADIENTS_2D({
			Grad2( 0.130526192220052,  0.99144486137381),
			Grad2( 0.38268343236509,   0.923879532511287),
			Grad2( 0.608761429008721,  0.793353340291235),
			Grad2( 0.793353340291235,  0.608761429008721),
			Grad2( 0.923879532511287,  0.38268343236509),
			Grad2( 0.99144486137381,   0.130526192220051),
			Grad2( 0.99144486137381,  -0.130526192220051),
			Grad2( 0.923879532511287, -0.38268343236509),
			Grad2( 0.793353340291235, -0.60876142900872),
			Grad2( 0.608761429008721, -0.793353340291235),
			Grad2( 0.38268343236509,  -0.923879532511287),
			Grad2( 0.130526192220052, -0.99144486137381),
			Grad2(-0.130526192220052, -0.99144486137381),
			Grad2(-0.38268343236509,  -0.923879532511287),
			Grad2(-0.608761429008721, -0.793353340291235),
			Grad2(-0.793353340291235, -0.608761429008721),
			Grad2(-0.923879532511287, -0.38268343236509),
			Grad2(-0.99144486137381,  -0.130526192220052),
			Grad2(-0.99144486137381,   0.130526192220051),
			Grad2(-0.923879532511287,  0.38268343236509),
			Grad2(-0.793353340291235,  0.608761429008721),
			Grad2(-0.608761429008721,  0.793353340291235),
			Grad2(-0.38268343236509,   0.923879532511287),
			Grad2(-0.130526192220052,  0.99144486137381)
});
