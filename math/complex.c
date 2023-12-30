//
// Created by peter on 12/27/2023.
//

#include <stdio.h>
#include "complex.h"
#include "math.h"


/////////////////////////BASIC//////////////////////////
/**
 * Gives complex number
 * @param x real value
 * @param y imaginary value
 * @return complex number (x+iy)
 */
static Complex get(float x, float y) {
    Complex complex = {.x=x, .y=y};
    return complex;
}

/**
 * Gives complex number
 * @param x real value
 * @return complex number (r+i0)
 */
static Complex getFromReal(float x) {
    Complex complex = {.x=x, .y=0};
    return complex;
}

/**
 * Gives complex number
 * @param r magnitude value
 * @param theta angle
 * @return complex number (r<theta)
 */
static Complex getFromPolar(float r, float theta) {
    return get(r * cosf(theta), r * sinf(theta));
}

//////////////////////////COMPLEX OPERATIONS////////////////////////////
/**
 * Gives conjugate
 * @param c complex number
 * @return conjugate of c
 */
static Complex conjugate(Complex c) {
    return get(c.x, -c.y);
}

/**
 * Gives magnitude of complex number
 * @param c complex number
 * @return magnitude of c
 */
static float mag(Complex c) {
    return sqrtf(c.x * c.x + c.y * c.y);
}

/**
 * Gives angle of complex number
 * @param c complex number
 * @return angle of c
 */
static float angle(Complex c) {
    return atan2f(c.y, c.x);
}

/**
 * Normalize complex number
 * @param c complex number
 * @return c/|c|
 */
static Complex normalise(Complex c) {
    float r = mag(c);
    return get(c.x / r, c.y / r);
}

//////////////////////////ARITHMETIC OPERATIONS//////////////////////////
/**
 * Add two complex numbers
 * @param c1 first complex number
 * @param c2 second complex number
 * @return c1+c2
 */
static Complex add(Complex c1, Complex c2) {
    return get(c1.x + c2.x, c1.y + c2.y);
}

/**
 * Gives difference of two complex numbers
 * @param c1 first complex number
 * @param c2 second complex number
 * @return c1-c2
 */
static Complex sub(Complex c1, Complex c2) {
    return get(c1.x - c2.x, c1.y - c2.y);
}

/**
 * Multiply two complex numbers
 * @param c1 first complex number
 * @param c2 second complex number
 * @return c1xc2
 */
static Complex dot(Complex c1, Complex c2) {
    return getFromPolar(mag(c1) * mag(c2), angle(c1) + angle(c2));
}

/**
 * Gives ratio of two complex numbers
 * @param c1 first complex number
 * @param c2 second complex number
 * @return c1/c2
 */
static Complex div(Complex c1, Complex c2) {
    return getFromPolar(mag(c1) / mag(c2), angle(c1) - angle(c2));
}

////////////////////////EXPONENTIAL,POWER AND LOGARITHMIC FUNCTIONS////////////////////
/**
 * Gives exponential of complex number
 * @param c complex number
 * @return exp(c)
 */
static Complex expC(Complex c) {
    return getFromPolar(expf(c.x), c.y);
}

/**
 * Gives natural log of complex number
 * @param c complex number
 * @return ln(c)
 */
static Complex lnC(Complex c) {
    return get(logf(mag(c)), angle(c));
}

/**
 * Gives power of complex number raise to complex number
 * @param c1 first complex number
 * @param c2 second complex number
 * @return c1^c2
 */
static Complex powC(Complex c1, Complex c2) {
    return expC(dot(c2, lnC(c1)));
}


//////////////////////////////TRIGONOMETRIC FUNCTIONS//////////////////////////
/**
 * Gives sine of complex number
 * @param c complex number
 * @return sin(c)
 */
static Complex sinC(Complex c) {
    return get(sinf(c.x) * coshf(c.y), cosf(c.x) * sinhf(c.y));
}

/**
 * Gives cosine of complex number
 * @param c complex number
 * @return cos(c)
 */
static Complex cosC(Complex c) {
    return get(cosf(c.x) * coshf(c.y), -sinf(c.x) * sinhf(c.y));
}

/**
 * Gives tangent of complex number
 * @param c complex number
 * @return tan(c)
 */
static Complex tanC(Complex c) {
    return div(sinC(c), cosC(c));
}
///////////////////////////INVERSE TRIGONOMETRIC FUNCTIONS////////////////////
/**
 * Gives inverse of sine of c
 * @param c complex number
 * @return asin(c)
 */
static Complex asinC(Complex c) {
    Complex iz = dot(get(0, 1), c);
    Complex sqrt_1_z2 = powC(sub(get(1, 0), powC(c, get(2, 0))), get(0.5f, 0));
    Complex x = add(iz, sqrt_1_z2);
    return dot(lnC(x), get(0, -1));
}

/**
 * Gives inverse of cos of c
 * @param c complex number
 * @return acos(c)
 */
static Complex acosC(Complex c) {
    Complex i_sqrt_1_z2 = dot(get(0, 1), powC(sub(get(1, 0), powC(c, get(2, 0))), get(0.5f, 0)));
    Complex x = add(c, i_sqrt_1_z2);
    return dot(lnC(x), get(0, -1));
}

/**
 * Gives inverse of tangent of c
 * @param c complex number
 * @return atan(c)
 */
static Complex atanC(Complex c) {
    Complex iz = dot(get(0, 1), c);
    Complex x = div(add(get(1, 0), iz), sub(get(1, 0), iz));
    return dot(lnC(x), get(0, -0.5f));
}

/**
 * Gives inverse of sine of c
 * @param c complex number
 * @return acos(c)
 */
static Complex atan2C(Complex y, Complex x) {
    Complex r = powC(add(powC(x, get(2, 0)), powC(y, get(2, 0))), get(0.5f, 0));
    return dot(get(2, 0), atanC(div(y, add(r, x))));
}

//////////////////////////////Hyperbolic FUNCTIONS//////////////////////////
/**
 * Gives hyperbolic sine of complex number
 * @param c complex number
 * @return sinh(c)
 */
static Complex sinhC(Complex c) {
    return get(sinhf(c.x) * cosf(c.y), coshf(c.x) * sinf(c.y));
}

/**
 * Gives hyperbolic cosine of complex number
 * @param c complex number
 * @return cosh(c)
 */
static Complex coshC(Complex c) {
    return get(coshf(c.x) * cosf(c.y), sinhf(c.x) * sinf(c.y));
}

/**
 * Gives hyperbolic tangent of complex number
 * @param c complex number
 * @return tanh(c)
 */
static Complex tanhC(Complex c) {
    return div(sinhC(c), coshC(c));
}
///////////////////////////INVERSE HYPERBOLIC FUNCTIONS////////////////////
/**
 * Gives inverse of hyperbolic sine of c
 * @param c complex number
 * @return asinh(c)
 */
static Complex asinhC(Complex c) {
    Complex sqrt_1_z2 = powC(add(powC(c, get(2, 0)),get(1, 0)), get(0.5f, 0));
    Complex x = add(c, sqrt_1_z2);
    return lnC(x);
}

/**
 * Gives inverse of hyperbolic cosine of c
 * @param c complex number
 * @return asinh(c)
 */
static Complex acoshC(Complex c) {
    Complex sqrt_1_z2 = powC(sub(powC(c, get(2, 0)),get(1, 0)), get(0.5f, 0));
    Complex x = add(c, sqrt_1_z2);
    return lnC(x);
}

/**
 * Gives inverse of hyperbolic tangent of c
 * @param c complex number
 * @return atanh(c)
 */
static Complex atanhC(Complex c) {
    Complex x = div(add(get(1, 0), c), sub(get(1, 0), c));
    return dot(lnC(x), get(0.5f, 0));
}


//////////////////////////PRINTING//////////////////////////////////////////
/**
 * Print complex number in form x+iy
 * @param c complex number
 */
static void print(Complex c) {
    if((c.y*c.y)<(COMPLEX_EPSILON))
        c.y = 0.0f;

    if (c.y >= 0.0f)
        printf("%.3f + i%.3f", c.x, c.y);
    else
        printf("%.3f - i%.3f", c.x, -c.y);
}

/**
 * Print complex number in form r<theta
 * @param c complex number
 */
static void printA(Complex c) {
    float r = mag(c);
    float t = angle(c);
    printf("%.3f <%.2f", r, t*180.0f/(float)M_PI);
}

__attribute__((unused)) struct ComplexControl StaticComplex = {
        .get = get,
        .getFromReal = getFromReal,
        .getFromPolar = getFromPolar,
        .conjugate = conjugate,
        .mag = mag,
        .angle = angle,
        .normalise = normalise,
        .add = add,
        .sub = sub,
        .dot = dot,
        .div =div,
        .exp = expC,
        .ln = lnC,
        .pow = powC,
        .sin = sinC,
        .cos = cosC,
        .tan = tanC,
        .asin = asinC,
        .acos = acosC,
        .atan = atanC,
        .atan2 = atan2C,
        .sinh = sinhC,
        .cosh = coshC,
        .tanh = tanhC,
        .asinh = asinhC,
        .acosh = acoshC,
        .atanh = atanhC,
        .print = print,
        .printA = printA
};