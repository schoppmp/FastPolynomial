#include <iostream>
#include <cmath>
#include <NTL/ZZ.h>


using namespace std;
using namespace NTL;

#include "iterative.h"
#include "recursive.h"

using namespace std;
using namespace NTL;
using namespace chrono;

//#define ITERATIVE

void usage_example_zp(){
    long degree = pow(2,20)-1;

    ZZ prime;
    GenPrime(prime, 128);
    ZZ_p::init(prime);

//  interpolation points:
    ZZ_p* X = new ZZ_p[degree+1];
    ZZ_p* Y = new ZZ_p[degree+1];
    for(unsigned int i=0;i<=degree; i++) {
        random(X[i]);
        random(Y[i]);
    }

    ZZ_pX P;
#ifdef ITERATIVE
    poly_interpolate_zp_iterative(degree, X, Y, P);
#else
    poly_interpolate_zp_recursive(degree, X, Y, P);
#endif

    // EVALUATE
    ZZ_p* X2 = new ZZ_p[degree+1];
    ZZ_p* Y2 = new ZZ_p[degree+1];
    for(unsigned int i=0;i<=degree; i++) {
        random(X[i]);
    }
#ifdef ITERATIVE
    poly_evaluate_zp_iterative(degree, P, X2, Y2);
#else
    poly_evaluate_zp_recursive(degree, P, X2, Y2);
#endif

    delete[] X;
    delete[] Y;
    delete[] X2;
    delete[] Y2;
}


int main() {

    usage_example_zp();

}
