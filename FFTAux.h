#ifndef FFTAUX_H
#define FFTAUX_H

#define N 65536/2
#define PRINT_LIMIT 16

#define STRONGT    "\033[1m\x1b[35m"
#define RESETT   "\x1b[0m"

#include <stdio.h>
#include<math.h>

typedef struct {
    double real;
    double imag;
} Complex;


Complex multiply(Complex *a, Complex *b);

Complex add(Complex *a, Complex *b);

void computeEulers(Complex *Euler);

void printResult(Complex *Result);

void printResultSeparate(double *ResultR, double *ResultI);

#endif