#ifndef FASTPOLY_RECURSIVE_H
#define FASTPOLY_RECURSIVE_H

namespace fastpoly {

template<typename ElementType, typename PolynomialType>
void poly_interpolate_zp_recursive(long degree, const ElementType *X,
                                   const ElementType *Y, PolynomialType &P);

template<typename ElementType, typename PolynomialType>
void poly_evaluate_zp_recursive(long degree, const PolynomialType &P,
                                const ElementType *X, ElementType *Y);

}

#endif //FASTPOLY_RECURSIVE_H
