//***************************************************************
// radix-2 FFT algorithm
// Todos os processos calculam (N/2)/np dos valores de Euler
// (e^(-theta))  e os distribui para todos os processos
// usando MPI_Allgather. Assim todos os processos computam sua parte
// e reunem para o processo 0 usando MPI_Gather. Os resultados imaginários
// são separados dos reais
//
// Compilação:  mpicc -o teste FFTP.c FFTAux.c
// Execução: mpirun -np 2 teste
//
//*****************************************************************

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<mpi.h>
#include<stdlib.h>
#include "FFTAux.h"


int main(int argc, char **argv) {
    int size, rank;
    Complex Input[N];

    //Entrada: array de tamanho N com valores
    // definidos apenas nas 8 entradas
    loadData(Input);

    //Iniciar MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    int localN = (N / 2) / size;
    double EulerR[N / 2], EulerI[N / 2];
    double *ResultR = NULL, *ResultI = NULL;
    double tempEulerI[localN], tempEulerR[localN];
    double tempResultR[localN * 2], tempResultI[localN * 2];


    // Variaveis para calcular o tempo
    struct timespec now, tmstart;
    double mpi_begin, mpi_end;

    if (rank == 0) {
        //Armazenar resultados apenas no processo root
        ResultR = malloc(sizeof(double) * N);
        ResultI = malloc(sizeof(double) * N);

        //Iniciar contagem do tempo
        clock_gettime(CLOCK_REALTIME, &tmstart);
        mpi_begin = MPI_Wtime();
    }

    //Calcula o valor de Euler (e^(-theta))
    double theta, ang = 4.0 * M_PI / N;
    int x;
    for (x = 0; x < localN; x++) {
        theta = (localN * rank + x) * ang;
        tempEulerR[x] = cos(theta);
        tempEulerI[x] = -sin(theta);
    }

    //Distribuir os valores para todos os processos
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
            // pegar o valor de Euler atual
            euler.real = EulerR[idx];
            euler.imag = EulerI[idx];

            //multiplicar os pares com o valor de Euler
            temp = multiply(&Input[n << 1], &euler);
            // somar resultado
            even = add(&even, &temp);

            //multiplicar os ímpares com o valor de Euler
            temp = multiply(&Input[(n << 1) + 1], &euler);
            //somar resultado
            odd = add(&odd, &temp);

            //compute index for next euler component
            idx = (idx + diff + 1) % (N >> 1);
        }

        //calcular o twiddle
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

    //Reunir os resultados para o processo 0

    //Primeiros N/2 resultados
    MPI_Gather(tempResultR, localN, MPI_DOUBLE, ResultR, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(tempResultI, localN, MPI_DOUBLE, ResultI, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //Segundos N/2 resultados
    MPI_Gather(tempResultR + localN, localN, MPI_DOUBLE, ResultR + (N / 2), localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(tempResultI + localN, localN, MPI_DOUBLE, ResultI + (N / 2), localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {


        // Calcular o tempo de processamento
        mpi_end = MPI_Wtime();
        clock_gettime(CLOCK_REALTIME, &now);
        double seconds = (double) ((now.tv_sec + now.tv_nsec * 1e-9) -
                                   (double) (tmstart.tv_sec + tmstart.tv_nsec * 1e-9));

        printf("\n");
        printf("Número de processos:\t "STRONGT"%d"RESETT" proc:\n", size);
        // Tempo total do programa
        printf("Tempo em paralelo:\t "STRONGT"%f"RESETT" secs\n", seconds);
        // Tempo do MPI
        printf("Tempo total do MPI: "STRONGT"%f"RESETT" secs\n\n", mpi_end - mpi_begin);

        printResultSeparate(ResultR, ResultI);
    }
    MPI_Finalize();
}