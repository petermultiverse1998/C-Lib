//
// Created by peter on 12/27/2023.
//

#ifndef C_LIB_COMPLEX_H
#define C_LIB_COMPLEX_H

#define COMPLEX_EPSILON 1e-9f

typedef struct {
    float x, y;
} Complex;

struct ComplexControl {
    //////////////////////////BASIC//////////////////////////////////////////
    /**
     * Gives complex number
     * @param x real value
     * @param y imaginary value
     * @return complex number (x+iy)
     */
    Complex (*get)(float x, float y);

    /**
     * Gives complex number
     * @param x real value
     * @return complex number (r+i0)
     */
    Complex (*getFromReal)(float x);

    /**
     * Gives complex number
     * @param r magnitude value
     * @param theta angle
     * @return complex number (r<theta)
     */
    Complex (*getFromPolar)(float r, float theta);
    //////////////////////////COMPLEX OPERATIONS////////////////////////////
    /**
     * Gives conjugate
     * @param c complex number
     * @return conjugate of c
     */
    Complex (*conjugate)(Complex c);

    /**
     * Gives magnitude of complex number
     * @param c complex number
     * @return magnitude of c
     */
    float (*mag)(Complex c);

    /**
     * Gives angle of complex number
     * @param c complex number
     * @return angle of c
     */
    float (*angle)(Complex c);

    /**
     * Normalize complex number
     * @param c complex number
     * @return c/|c|
     */
    Complex (*normalise)(Complex c);
    //////////////////////////ARITHMETIC OPERATIONS//////////////////////////
    /**
     * Add two complex numbers
     * @param c1 first complex number
     * @param c2 second complex number
     * @return c1+c2
     */
    Complex (*add)(Complex c1, Complex c2);

    /**
     * Gives difference of two complex numbers
     * @param c1 first complex number
     * @param c2 second complex number
     * @return c1-c2
     */
    Complex (*sub)(Complex c1, Complex c2);

    /**
     * Multiply two complex numbers
     * @param c1 first complex number
     * @param c2 second complex number
     * @return c1xc2
     */
    Complex (*dot)(Complex c1, Complex c2);

    /**
     * Gives ratio of two complex numbers
     * @param c1 first complex number
     * @param c2 second complex number
     * @return c1/c2
     */
    Complex (*div)(Complex c1, Complex c2);
    ////////////////////////EXPONENTIAL,POWER AND LOGARITHMIC FUNCTIONS////////
    /**
     * Gives exponential of complex number
     * @param c complex number
     * @return exp(c)
     */
    Complex (*exp)(Complex c);

    /**
     * Gives natural log of complex number
     * @param c complex number
     * @return ln(c)
     */
    Complex (*ln)(Complex c);

    /**
     * Gives power of complex number raise to complex number
     * @param c1 first complex number
     * @param c2 second complex number
     * @return c1^c2
     */
    Complex (*pow)(Complex c1, Complex c2);
    //////////////////////////////TRIGONOMETRIC FUNCTIONS///////////////////////
    /**
     * Gives sine of complex number
     * @param c complex number
     * @return sin(c)
     */
    Complex (*sin)(Complex c);

    /**
     * Gives cosine of complex number
     * @param c complex number
     * @return cos(c)
     */
    Complex (*cos)(Complex c);

    /**
     * Gives tangent of complex number
     * @param c complex number
     * @return tan(c)
     */
    Complex (*tan)(Complex c);
    ///////////////////////////INVERSE TRIGONOMETRIC FUNCTIONS//////////////////
    /**
     * Gives inverse of sine of c
     * @param c complex number
     * @return asin(c)
     */
    Complex (*asin)(Complex c);

    /**
     * Gives inverse of cosine of c
     * @param c complex number
     * @return acos(c)
     */
    Complex (*acos)(Complex c);

    /**
     * Gives inverse of tangent of c
     * @param c complex number
     * @return atan(c)
     */
    Complex (*atan)(Complex c);

    /**
     * Gives inverse of tangent of c
     * @param y perpendicular complex number
     * @param x base complex number
     * @return atan2(y,x)
     */
    Complex (*atan2)(Complex y, Complex x);
    //////////////////////////////Hyperbolic FUNCTIONS//////////////////////////
    /**
     * Gives hyperbolic sine of complex number
     * @param c complex number
     * @return sinh(c)
     */
    Complex (*sinh)(Complex c);

    /**
     * Gives hyperbolic cosine of complex number
     * @param c complex number
     * @return cosh(c)
     */
    Complex (*cosh)(Complex c);

    /**
     * Gives hyperbolic tangent of complex number
     * @param c complex number
     * @return tanh(c)
     */
    Complex (*tanh)(Complex c);

    ///////////////////////////INVERSE HYPERBOLIC FUNCTIONS////////////////////
    /**
     * Gives inverse of hyperbolic sine of c
     * @param c complex number
     * @return asinh(c)
     */
    Complex (*asinh)(Complex c);

    /**
     * Gives inverse of hyperbolic cosine of c
     * @param c complex number
     * @return asinh(c)
     */
    Complex (*acosh)(Complex c);

    /**
     * Gives inverse of hyperbolic tangent of c
     * @param c complex number
     * @return atanh(c)
     */
    Complex (*atanh)(Complex c);

    //////////////////////////PRINTING//////////////////////////////////////////
    /**
     * Print complex number in form x+iy
     * @param c complex number
     */
    void (*print)(Complex c);

    /**
     * Print complex number in form r<theta
     * @param c complex number
     */
    void (*printA)(Complex c);
};

__attribute__((unused)) extern struct ComplexControl StaticComplex;


#endif //C_LIB_COMPLEX_H
