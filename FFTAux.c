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

//***************************************************************************
//This function define the 8 first elements of an input array as some values
//and the rest as 0
//***************************************************************************
void loadData(Complex *Input) {

    Input[0].real = 2.6;
    Input[1].real = 6.3;
    Input[2].real = 4.0;
    Input[3].real = 9.1;
    Input[4].real = 0.4;
    Input[5].real = 4.8;
    Input[6].real = 2.6;
    Input[7].real = 4.1;

    Input[0].imag = 3.6;
    Input[1].imag = 2.9;
    Input[2].imag = 5.6;
    Input[3].imag = 4.8;
    Input[4].imag = 3.3;
    Input[5].imag = 5.9;
    Input[6].imag = 5.0;
    Input[7].imag = 4.3;


    int n;
    for (n = 8; n < N; n++) {
        Input[n].real = Input[n].imag = 0.0;
    }
}


void printHead() {
    printf("Total de valores processados: " STRONGT "%d\n\n"RESETT, N);
}


void printResult(Complex *Result) {

    printHead();

    int k;
    for (k = 0; k < PRINT_LIMIT; k++) {
        printf("X[%d]=\t  ( %.7f, %.7f i )\n", k, Result[k].real, Result[k].imag);
    }
}

void printResultSeparate(double *ResultR, double *ResultI) {

    printHead();

    int k;
    for (k = 0; k < PRINT_LIMIT; k++) {
        printf("X[%d]=\t  ( %.7f, %.7f i )\n", k, ResultR[k], ResultI[k]);
    }
}