#include "FFTAux.h"

Complex multiply(Complex *a, Complex *b) {
    Complex c;
    c.real = (a->real * b->real) - (a->imag * b->imag);
    c.imag = (a->real * b->imag) + (b->real * a->imag);
    return c;
}

Complex add(Complex *a, Complex *b) {
    Complex c;
    c.real = a->real + b->real;
    c.imag = a->imag + b->imag;
    return c;
}


//***************************************************************************
//The function precomputes the euler values from x = 0 to N/2 of e^(-4*PI*x/N)
//to avoid recomputation of sines and cosines in the FFT
//***************************************************************************
void computeEulers(Complex *Euler) {
    int x = 0;
    float theta;
    float n = (4.0 * M_PI) / N;

    for (x = 0; x < (N >> 1); x++) {
        theta = x * n;
        Euler[x].real = cos(theta);
        Euler[x].imag = -sin(theta);
    }
}

void printHead() {
    printf("Total de valores processados: " STRONGT "%d\n\n"RESETT, N);
}


void printResult(Complex *Result) {

    printHead();

    int k;
    for (k = 0; k < PRINT_LIMIT; k++) {
        printf("X[%d]=\t  ( %.7f, %.7fi )\n", k, Result[k].real, Result[k].imag);
    }
}

void printResultSeparate(double *ResultR, double *ResultI) {

    printHead();

    int k;
    for (k = 0; k < PRINT_LIMIT; k++) {
        printf("X[%d]=\t  ( %.7f, %.7fi )\n", k, ResultR[k], ResultI[k]);
    }
}