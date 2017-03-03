/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 * 
 *  MPI polynomial evaluation algorithm function implementations go here
 * 
 */

#include "mpi_evaluator.h"
#include "const.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>

void scatter(const int n, double* scatter_values, int &n_local, double* &local_values, int source_rank, const MPI_Comm comm){
    //Implementation
    int rank , size;
    MPI_Comm_size(comm , &size);
    MPI_Comm_rank(comm , &rank);
    if (rank == source_rank) {
        int m = int(n / size);
        int extra = n % size;
        int z = m + 1;
        for (int i = 0; i < size; i++) {
            if (i < extra) {
                if (i == rank) {
                    n_local = z;

                    local_values = (double *) malloc(z * sizeof(double));
                    if (local_values == NULL) printf("local_values is NULL on rank %d\n", rank);

                    for (int j = i * (m + 1), k = 0; k < z; j++, k++) { // j < (((i+1) * m) + 1) && 
                        local_values[k] = scatter_values[j];
                    }
                } else {
                    MPI_Send(&z, 1, MPI_INT, i, 0, comm);

                    double* local_values_i = (double *) malloc(z * sizeof(double));
                    if (local_values_i == NULL) printf("local_values_i is NULL on rank %d\n", rank);

                    for (int j = i * (m + 1), k = 0; k < z; j++, k++) { // j < ((i+1) * m + 1) && 
                        local_values_i[k] = scatter_values[j];
                    }
                    MPI_Send(local_values_i, z, MPI_DOUBLE, i, 1, comm);
                }
            } else {
                if (i == rank) {
                    n_local = m;

                    local_values = (double *) malloc(m * sizeof(double));
                    if (local_values == NULL) printf("local_values is NULL on rank %d\n", rank);

                    for (int j = (i * m) + extra, k = 0; k < m; j++, k++) {
                        local_values[k] = scatter_values[j];
                    }
                } else {
                    MPI_Send(&m, 1, MPI_INT, i, 0, comm);

                    double* local_values_i = (double *) malloc(m * sizeof(double));
                    if (local_values_i == NULL) printf("local_values_i is NULL on rank %d\n", rank);

                    for (int j = (i * m) + extra, k = 0; k < m; j++, k++) { // j < ((i+1) * m) && 
                        local_values_i[k] = scatter_values[j];
                    }
                    MPI_Send(local_values_i, m, MPI_DOUBLE, i, 1, comm);
                }
            }
        }
    } else {
        MPI_Status stat_n, stat_vals;
        MPI_Recv(&n_local, 1, MPI_INT, source_rank, 0, comm, &stat_n);
        local_values = (double *) malloc(n_local * sizeof(double));
        MPI_Recv(local_values, n_local, MPI_DOUBLE, source_rank, 1, comm, &stat_vals);
    }
    MPI_Barrier(comm);
}

double broadcast(double value, int source_rank, const MPI_Comm comm){
    //Implementation
    int rank , size, rank2;
    MPI_Comm_size (comm , &size );
    MPI_Comm_rank (comm , &rank );
    double val;
    if(rank == source_rank){
        val = value;
    } else{
        val = 0;
    }
    int d = ceil(log2(size));
    for(int i = 0; i < d; i++){
        rank2 = rank^(1<<(i));
        if((rank2 < size)){
            MPI_Send(&val,1,MPI_DOUBLE,rank2,111,comm);
            MPI_Status stat;
            MPI_Recv(&value,1,MPI_DOUBLE,rank2,111,comm,&stat);
            if(val == 0){
                val = value;
            }
        }
        MPI_Barrier(comm);
    }
    return val;
}

void parallel_prefix(const int n, const double* values, double* prefix_results, const int OP, const MPI_Comm comm){
    //Implementation
    int rank , size, rank2, total;
    MPI_Comm_size (comm , &size );
    MPI_Comm_rank (comm , &rank );

    prefix_results[0] = values[0];
    if (OP == PREFIX_OP_SUM) {
        for (int i = 1; i < n; i++) {
            prefix_results[i] = prefix_results[i - 1] + values[i];
        }
    } else if (OP == PREFIX_OP_PRODUCT) {
        for (int i = 1; i < n; i++) {
            prefix_results[i] = prefix_results[i - 1] * values[i];
        }
    }
    prefix_results[n] = prefix_results[n - 1];

    int steps = 0;
    int d = log2(size+1);
    //fix with??? int d = ceil(log2(size));
    if(d != log2(size+1)){
        d+=1;
    }
    for(int i = 0; i < d; i++){
        steps++;
        rank2 = rank^(1<<(i));
        total = prefix_results[n];
        if((rank2 < size)){
            MPI_Send(&total, 1, MPI_DOUBLE, rank2, 111, comm );
            MPI_Status stat;
            MPI_Recv(&total, 1, MPI_DOUBLE, rank2, MPI_ANY_TAG, comm, &stat);
            if (OP == PREFIX_OP_SUM) {
                if (rank > rank2) {
                    for (int i = 0; i <= n; i++) {
                        prefix_results[i] += total;
                    }
                } else {
                        prefix_results[n] += total;
                }
            } else if (OP == PREFIX_OP_PRODUCT) {
                if (rank > rank2) {
                    for (int i = 0; i <= n; i++) {
                        prefix_results[i] *= total;
                    }
                } else {
                        prefix_results[n] *= total;
                }
            }
        }
    }
}

double mpi_poly_evaluator(const double x, const int n, const double* constants, const MPI_Comm comm){
    //Implementation
    int rank;
    MPI_Comm_rank (comm , &rank );
    double answer = 0;
    double* pre_vec = (double *) malloc(n * sizeof(double));
    double* prefix_results = (double *) malloc((n + 1) * sizeof(double));
    if (pre_vec == NULL) { printf("pre_vec is null on %d", rank); }
    if (prefix_results == NULL) { printf("prefix_results is null on %d", rank); }

    if(rank == 0){
        pre_vec[0] = 1;
    } else {
        pre_vec[0] = x;
    }

    for (int i = 1; i < n; i++) {
        pre_vec[i] = x;
    }

    parallel_prefix(n, pre_vec, prefix_results, PREFIX_OP_PRODUCT, comm);

    for (int i = 0; i < n; i++) {
        pre_vec[i] = prefix_results[i] * constants[i];
    }
    parallel_prefix(n, pre_vec, prefix_results, PREFIX_OP_SUM, comm);
    answer = prefix_results[n];

    free(pre_vec);
    free(prefix_results);

    return answer;
}

