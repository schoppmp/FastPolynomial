#include "recursive.h"
#include "utils.h"

#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/vec_GF2E.h>

#include <chrono>
#include <sys/resource.h>

using namespace std;
using namespace NTL;
using namespace chrono;

#define LEFT(X) (2*X+1)
#define RIGHT(X) (2*X+2)

namespace fastpoly {

/* A recursive function to build the tree of polynomials
 * (assumming a complete binary tree => size = 2*#leafs-1
 'tree_size' is the number of nodes (including leaves) in the tree =
 2*(degree+1)-1 = 2*degree+1 'root' is the index of the subtree in the array
 'tree'
 */
template<typename ElementType, typename PolynomialType>
void build_tree(PolynomialType *tree, const ElementType *points, unsigned int root,
                unsigned int tree_size) {
  // halting condition
  if (LEFT(root) >= tree_size) {
    unsigned int point_index = root - (tree_size - 1) / 2;
    // setting the polynomial to be x-m where m is points[point_index]
    ElementType negated;
    NTL::negate(negated, points[point_index]);
    SetCoeff(tree[root], 0, negated);
    SetCoeff(tree[root], 1, 1);
    return;
  }

  build_tree(tree, points, LEFT(root), tree_size);
  build_tree(tree, points, RIGHT(root), tree_size);
  tree[root] = tree[LEFT(root)] * tree[RIGHT(root)];
}

template<typename ElementType, typename PolynomialType>
void test_tree(PolynomialType &final_polynomial, ElementType *points, unsigned int npoints) {
  ElementType result;
  bool error = false;
  for (unsigned int i = 0; i < npoints; i++) {
    result = eval(final_polynomial, points[i]);
    if (0 != result) {
      cout << "FATAL ERROR: polynomials tree is incorrect!" << endl;
      error = true;
      break;
    }
  }
  if (!error)
    cout << "polynomials tree is correct." << endl;
}

/*
 * P - the polynomial to evaluate_zp_iterative
 * tree - the subproduct tree over the x points that we want to recursive_evaluate_zp root - the current subtree tree size - the size of a complete tree is 2*n-1 where n is the number of leafs results - the evaluation result over the x's (that are represented by the tree)
 */
template<typename ElementType, typename PolynomialType>
void recursive_evaluate_zp(const PolynomialType &P, PolynomialType *tree, unsigned int root,
                           unsigned int tree_size, ElementType *results) {
  // halting condition
  if (LEFT(root) >= tree_size) {
    PolynomialType R = P % tree[root];
    if (deg(R) > 0)
      cout << "ERROR: R should be constant...";
    unsigned int result_index = root - (tree_size - 1) / 2;
    results[result_index] = coeff(R, 0);
    return;
  }

  PolynomialType R = P % tree[root];
  recursive_evaluate_zp(R, tree, LEFT(root), tree_size, results);
  recursive_evaluate_zp(R, tree, RIGHT(root), tree_size, results);
}

template<typename ElementType, typename PolynomialType>
void test_evaluate_zp_recursive(PolynomialType &P, ElementType *points, ElementType *results,
                                unsigned int npoints) {
  bool error = false;
  for (unsigned int i = 0; i < npoints; i++) {
    ElementType y = eval(P, points[i]);
    if (y != results[i]) {
      cout << "y=" << y << " and results[i]=" << results[i] << endl;
      error = true;
    }
  }
  if (error)
    cout << "ERROR: evaluation results do not match real evaluation!" << endl;
  else
    cout << "All evaluation results computed correctly!" << endl;
}

template<typename ElementType, typename PolynomialType>
void multipoint_evaluate_zp_recursive(long degree, const PolynomialType &P,
                                      const ElementType *X, ElementType *Y) {
  // we want to recursive_evaluate_zp P on 'degree+1' values.
  PolynomialType *p_tree = new PolynomialType[degree * 2 + 1];
#ifndef NDEBUG
  steady_clock::time_point begin1 = steady_clock::now();
#endif
  build_tree(p_tree, X, 0, degree * 2 + 1);

#ifndef NDEBUG
  steady_clock::time_point end1 = steady_clock::now();
#endif
  //    test_tree_zp_iterative(p_tree[0], x, degree+1);

#ifndef NDEBUG
  steady_clock::time_point begin2 = steady_clock::now();
#endif
  recursive_evaluate_zp(P, p_tree, 0, degree * 2 + 1, Y);
#ifndef NDEBUG
  chrono::steady_clock::time_point end2 = steady_clock::now();
#endif

#ifndef NDEBUG
  cout << "Building tree: "
       << duration_cast<milliseconds>(end1 - begin1).count() << " ms" << endl;
  cout << "Evaluating points: "
       << duration_cast<milliseconds>(end2 - begin2).count() << " ms" << endl;
  cout << "Total: "
       << duration_cast<milliseconds>(end1 - begin1).count() +
              duration_cast<milliseconds>(end2 - begin2).count()
       << " ms" << endl;
#endif

  delete[] p_tree;
}

// void test_multipoint_eval_zp(ZZ prime, long degree)
//{
//    // init underlying prime field
//    ElementType::init(ZZ(prime));
//
//    // the given polynomial
//    PolynomialType P;
//    random(P, degree+1);
//    SetCoeff(P,degree,random_ElementType());
//
//    // evaluation points:
//    ElementType* x = new ElementType[degree+1];
//    ElementType* y = new ElementType[degree+1];
//
//    for(unsigned int i=0;i<=degree; i++) {
//        random(x[i]);
//    }
//
//    multipoint_evaluate_zp_iterative(P, x, y, degree);
//}

/*
 * expects an "empty" polynomial 'resultP'
 */
template<typename ElementType, typename PolynomialType>
void recursive_interpolate_zp(PolynomialType &resultP, unsigned int root, const ElementType *x,
                              const ElementType *y, ElementType *a, PolynomialType *M,
                              unsigned int tree_size) {
  // halting condition
  if (LEFT(root) >= tree_size) {
    unsigned int y_index = root - (tree_size - 1) / 2;
    ElementType inv_a;
    inv(inv_a, a[y_index]); // inv_a = 1/a
    SetCoeff(resultP, 0, y[y_index] * inv_a);
    return;
  }

  PolynomialType leftP, rightP;
  recursive_interpolate_zp(leftP, LEFT(root), x, y, a, M, tree_size);
  recursive_interpolate_zp(rightP, RIGHT(root), x, y, a, M, tree_size);

  resultP = leftP * M[RIGHT(root)] + rightP * M[LEFT(root)];
}

/*
 * We follow the algorithm and notation as in Moneck & Borodin '73
 */
template<typename ElementType, typename PolynomialType>
void interpolate_zp(long degree, const ElementType *X, const ElementType *Y, PolynomialType &resultP) {
#ifndef NDEBUG
  system_clock::time_point begin[4];
  system_clock::time_point end[4];
#endif

  // we first build the tree of the super moduli
  PolynomialType *M = new PolynomialType[degree * 2 + 1];
#ifndef NDEBUG
  begin[0] = system_clock::now();
#endif
  build_tree(M, X, 0, degree * 2 + 1);
#ifndef NDEBUG
  end[0] = system_clock::now();
#endif
  //    test_tree_zp_iterative(M[0], x, degree+1);

  // we construct a preconditioned global structure for the a_k for all 1<=k<=(degree+1)
  ElementType *a = new ElementType[degree + 1];
  PolynomialType d;
#ifndef NDEBUG
  begin[1] = system_clock::now();
#endif
  diff(d, M[0]);
#ifndef NDEBUG
  end[1] = system_clock::now();
#endif

  // recursive_evaluate_zp d(x) to obtain the results in the array a
#ifndef NDEBUG
  begin[2] = system_clock::now();
#endif
  recursive_evaluate_zp(d, M, 0, degree * 2 + 1, a);
#ifndef NDEBUG
  end[2] = system_clock::now();
#endif

  // now we can apply the recursive formula
#ifndef NDEBUG
  begin[3] = system_clock::now();
#endif
  recursive_interpolate_zp(resultP, 0, X, Y, a, M, degree * 2 + 1);
#ifndef NDEBUG
  end[3] = system_clock::now();
#endif

#ifndef NDEBUG
  cout << " -- Recursive --" << endl << endl;
  cout << "Building tree: "
       << duration_cast<milliseconds>(end[0] - begin[0]).count() << " ms"
       << endl;
  cout << "Differentiate: "
       << duration_cast<milliseconds>(end[1] - begin[1]).count() << " ms"
       << endl;
  cout << "Evaluate diff: "
       << duration_cast<milliseconds>(end[2] - begin[2]).count() << " ms"
       << endl;
  cout << "Interpolation: "
       << duration_cast<milliseconds>(end[3] - begin[3]).count() << " ms"
       << endl;
  cout << "Total: "
       << duration_cast<milliseconds>(end[0] - begin[0] + end[1] - begin[1] +
                                      end[2] - begin[2] + end[3] - begin[3])
              .count()
       << " ms" << endl;
#endif

  delete[] M;
  delete[] a;
}

template<typename ElementType, typename PolynomialType>
void test_interpolation_result_zp_recursive(long degree, ElementType *X, ElementType *Y,
                                            PolynomialType &P) {
  cout << "Testing result polynomial" << endl;
  ElementType res;
  for (long i = 0; i < degree + 1; i++) {
    eval(res, P, X[i]);
    if (res != Y[i]) {
      cout << "Error! x = " << X[i] << ", y = " << Y[i] << ", res = " << res
           << endl;
      return;
    }
  }
  cout << "Polynomial is interpolated correctly!" << endl;
}

template<typename ElementType, typename PolynomialType>
void poly_interpolate_zp_recursive(long degree, const ElementType *X, const ElementType *Y,
                                   PolynomialType &P) {

  interpolate_zp(degree, X, Y, P);
  // the next operation takes O(n^2) time, keep it commented out!
  //    test_interpolation_result_zp_recursive(degree, X, Y, P);
}

template<typename ElementType, typename PolynomialType>
void poly_evaluate_zp_recursive(long degree, const PolynomialType &P, const ElementType *X,
                                ElementType *Y) {
  multipoint_evaluate_zp_recursive(degree, P, X, Y);
  // the next operation takes O(n^2) time, keep it commented out!
  //    test_evaluate_zp_recursive(P,X,Y,degree+1);
}

// Explicit instantiations.

// Prime fields.
template
void poly_interpolate_zp_recursive(long degree, const NTL::ZZ_p *X,
                                   const ZZ_p *Y, ZZ_pX &P);

template
void poly_evaluate_zp_recursive(long degree, const ZZ_pX &P,
                                const ZZ_p *X, ZZ_p *Y);

// Binary extension fields.
template
void poly_interpolate_zp_recursive(long degree, const NTL::GF2E *X,
                                   const GF2E *Y, GF2EX &P);

template
void poly_evaluate_zp_recursive(long degree, const GF2EX &P,
                                const GF2E *X, GF2E *Y);

}