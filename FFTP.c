//***************************************************************
// This is a parallel implementation of the radix-2 FFT algorithm.
// All processes first calculate (N/2)/comm_size of the euler values
// (e^(-theta)) and then distribute to every other process using
// MPI_Allgather. Then each process computes (N/2)/comm_size of
// the result and gathers their results to process 0 which ouputs it.
//
// Compilation:  mpicc -o teste FFTP.c FFTAux.c
// Execution: mpirun -np 2 teste
//
//*****************************************************************

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<mpi.h>
#include<stdlib.h>
#include "FFTAux.h"


//Main
int main(int argc, char **argv) {
    int size, rank;
    Complex Input[N];

    //Input values
    loadData(Input);
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int localN = (N / 2) / size;
    double EulerR[N / 2], EulerI[N / 2];
    double *ResultR = NULL, *ResultI = NULL;
    double tempEulerI[localN], tempEulerR[localN];
    double tempResultR[localN * 2], tempResultI[localN * 2];

    struct timespec now, tmstart;
    double mpist, mpiend;

    if (rank == 0) {
        ResultR = malloc(sizeof(double) * N);
        ResultI = malloc(sizeof(double) * N);

        //start timers
        clock_gettime(CLOCK_REALTIME, &tmstart);
        mpist = MPI_Wtime();
    }

    //Compute Euler values
    double theta, ang = 4.0 * M_PI / N;
    int x;
    for (x = 0; x < localN; x++) {
        theta = (localN * rank + x) * ang;
        tempEulerR[x] = cos(theta);
        tempEulerI[x] = -sin(theta);
    }

    //Distribute euler values to all processes
    MPI_Allgather(tempEulerR, localN, MPI_DOUBLE, EulerR, localN, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(tempEulerI, localN, MPI_DOUBLE, EulerI, localN, MPI_DOUBLE, MPI_COMM_WORLD);

    Complex even, odd;

    Complex twiddle, temp, euler, result;

    double PI2_by_N = 2 * M_PI / N;

    int diff, idx;

    int n, k;
    for (x = 0; x < localN; x++) {
        k = rank * localN + x;
        even.real = even.imag = odd.real = odd.imag = 0.0;

        //Get difference between indices of Euler values for current k
        diff = (k - 1 + (N >> 1)) % (N >> 1);
        idx = 0; //start index is 0

        for (n = 0; n < (N >> 1); n++) {
            //get current euler component
            euler.real = EulerR[idx];
            euler.imag = EulerI[idx];

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
        result = add(&even, &temp);
        tempResultR[x] = result.real;
        tempResultI[x] = result.imag;

        //Subtract odd from even to get k+N/2 component
        temp.real = -temp.real;
        temp.imag = -temp.imag;
        result = add(&even, &temp);
        tempResultR[x + localN] = result.real;
        tempResultI[x + localN] = result.imag;
    }

    //Gather results to process 0

    //First N/2 results
    MPI_Gather(tempResultR, localN, MPI_DOUBLE, ResultR, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(tempResultI, localN, MPI_DOUBLE, ResultI, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //Second N/2 results
    MPI_Gather(tempResultR + localN, localN, MPI_DOUBLE, ResultR + (N / 2), localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(tempResultI + localN, localN, MPI_DOUBLE, ResultI + (N / 2), localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        //Redirect output to file
//        freopen("EnemKhalidParallel.txt", "w", stdout);

        //End timers
        mpiend = MPI_Wtime();
        clock_gettime(CLOCK_REALTIME, &now);
        double seconds = (double) ((now.tv_sec + now.tv_nsec * 1e-9) -
                                   (double) (tmstart.tv_sec + tmstart.tv_nsec * 1e-9));

        printf("\n");
        printf("Running at: "STRONGT"%d"RESETT" proc:\n", size);
        printf("C time: "STRONGT"%f"RESETT" secs\n", seconds);
        printf("MPI time: "STRONGT"%f"RESETT" secs\n\n", mpiend - mpist);

        printResultSeparate(ResultR, ResultI);
    }
    MPI_Finalize();
}