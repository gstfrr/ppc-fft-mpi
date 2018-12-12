//***************************************************************
// This is a serial implementation of the radix-2 FFT algorithm.
// This implementation first calculates all the euler values (e^(-itheta))
// and stores those in a complex array. Then the euler values are re-used in
// computing the FFT values from k = 0 to N-1. The main loop only goes from
// k = 1 to N/2 - 1 because the even and odd components for the first half
// can be reused to calculate the values for the last half.
//
// Compilation: gcc -o teste FFTS.c FFTAux.c
// Execution: ./teste
//
//*****************************************************************

#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include "FFTAux.h"


//Arrays for Input and Results
Complex Input[N];
Complex Result[N];

//Array for Euler values
Complex Euler[N / 2];

//Main
int main(int argc, char **argv) {
    //Redirect output to file
    //freopen("EnemKhalidSerial.txt", "w", stdout);

    //Input values
    loadData(Input);

    Complex even, odd;
    struct timespec now, tmstart;
    clock_gettime(CLOCK_REALTIME, &tmstart);

    //compute N/2 Euler values (cos(theta) - isin(theta))
    computeEulers(Euler);

    Complex twiddle, temp, euler;

    double theta, PI2_by_N = 2 * M_PI / N;

    int diff, idx;

    int n, k;
    for (k = 0; k < (N >> 1); k++) {
        even.real = even.imag = odd.real = odd.imag = 0.0;

        //Get difference between indices of Euler values for current k
        diff = (k - 1 + (N >> 1)) % (N >> 1);
//        printf("%d\t%d\n", k, diff);
        idx = 0; //start index is 0

        for (n = 0; n < (N >> 1); n++) {
            //get current euler component
            euler = Euler[idx];

            //multiply even input with euler component
            temp = multiply(&Input[n << 1], &euler);
            //add result to even
            even = add(&even, &temp);

            //multiply odd component with euler input
            temp = multiply(&Input[(n << 1) + 1], &euler);
            //add result to odd
            odd = add(&odd, &temp);

            //compute index for next euler component
            idx = (idx + diff + 1) % (N >> 1);
        }

        //Compute twiddle
        theta = k * PI2_by_N;
        twiddle.real = cos(theta);
        twiddle.imag = -sin(theta);
        //multiply twiddle with odd component
        temp = multiply(&odd, &twiddle);

        //Add even and odd to result
        Result[k] = add(&even, &temp);

        //Subtract odd from even to get k+N/2 component
        temp.real = -temp.real;
        temp.imag = -temp.imag;
        Result[k + (N >> 1)] = add(&even, &temp);
    }

    clock_gettime(CLOCK_REALTIME, &now);
    double seconds = (double) ((now.tv_sec + now.tv_nsec * 1e-9) - (double) (tmstart.tv_sec + tmstart.tv_nsec * 1e-9));

    printf("\nTempo sequencial:\t "STRONGT"%f"RESETT" segundos\n\n", seconds);

    //printResult(Result);

    return 0;

}

