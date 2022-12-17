#include <iostream>
#include "mpi.h"


int main(int argc, char* argv[])
{
    int N = 2000;
    int my_rank;            /*My process rank*/
    int comm_sz;            /*Number of processes*/
    int local_N;
    int i, j, k;
    double start, finish;    /*timer*/
    int tem;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    //Количество строк, выделенных для каждой матрицы
    local_N = N / comm_sz;

    //Определить две матрицы
    int* A = new int[N * N];
    int* B = new int[N * N];

    //Матрица присваивается каждому процессу
    int* loc_mat_one = new int[local_N * N];


    //Матрица результатов в каждом процессе
    int* loc_res = new int[local_N * N];

    //Матрица результатов
    int* C = new int[N * N];

    if (my_rank == 0)
    {
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
                A[i * N + j] = (i + 2) * (j + 1);

        for (j = 0; j < N; j++)
            for (i = 0; i < N; i++)
                B[i * N + j] = (i + 1) + 2 * (j + 3);

        start = MPI_Wtime();
        // Распределение данных
        MPI_Scatter(A, local_N * N, MPI_INT, loc_mat_one, local_N * N, MPI_INT, 0, MPI_COMM_WORLD);

        //Передача данных
        MPI_Bcast(B, N * N, MPI_INT, 0, MPI_COMM_WORLD);

        //расчетlocalрезультат
        for (i = 0; i < local_N; i++)
            for (j = 0; j < N; j++) {
                tem = 0;
                for (k = 0; k < N; k++)
                    tem += loc_mat_one[i * N + k] * B[j * N + k];
                loc_res[i * N + j] = tem;
            }

        delete[] loc_mat_one;

        //Сбор результатов
        MPI_Gather(loc_res, local_N * N, MPI_INT, C, local_N * N, MPI_INT, 0, MPI_COMM_WORLD);

        //Обработка оставшейся строки (обработка ситуации, которая не может быть делимой)
        int rest = N % comm_sz;
        if (rest != 0)
            for (i = N - rest - 1; i < N; i++)
                for (j = 0; j < N; j++) {
                    tem = 0;
                    for (k = 0; k < N; k++)
                        tem += A[i * N + k] * B[j * N + k];
                    C[i * N + j] = tem;
                }
        finish = MPI_Wtime();

        delete[] A;
        delete[] B;
        delete[] loc_res;

        std::cout << "Elapsed time = " << finish - start << " seconds" << std::endl;
    }
    else {
        //Распределение данных
        MPI_Scatter(A, local_N * N, MPI_INT, loc_mat_one, local_N * N, MPI_INT, 0, MPI_COMM_WORLD);

        //Передача данных
        MPI_Bcast(B, N * N, MPI_INT, 0, MPI_COMM_WORLD);

        //расчетlocalрезультат
        for (i = 0; i < local_N; i++)
            for (j = 0; j < N; j++) {
                tem = 0;
                for (k = 0; k < N; k++)
                    tem += loc_mat_one[i * N + k] * B[j * N + k];
                loc_res[i * N + j] = tem;
            }

        delete[] loc_mat_one;
        delete[] B;

        //Сбор результатов
        MPI_Gather(loc_res, local_N * N, MPI_INT, C, local_N * N, MPI_INT, 0, MPI_COMM_WORLD);
        delete[] loc_res;
    }
    MPI_Finalize();
    return 0;
}