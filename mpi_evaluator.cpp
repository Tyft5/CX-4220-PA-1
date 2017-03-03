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
    //set the rank and size
    int rank , size;
    MPI_Comm_size(comm , &size);
    MPI_Comm_rank(comm , &rank);
    //if scattering from the source rank
    if (rank == source_rank) {
        //size of n/p
        int m = int(n / size);
        //determine extra elements if p isn't base 2
        int extra = n % size;
        int z = m + 1;
        //loop through each other rank
        for (int i = 0; i < size; i++) {
            //if a rank gets an extra element
            if (i < extra) {
                //if it's the source rank
                if (i == rank) {
                    //set the local value
                    n_local = z;
                    //allocate memory for the scattered array
                    local_values = (double *) malloc(z * sizeof(double));
                    if (local_values == NULL) printf("local_values is NULL on rank %d\n", rank);
                    //store the scattered values
                    for (int j = i * (m + 1), k = 0; k < z; j++, k++) { // j < (((i+1) * m) + 1) && 
                        local_values[k] = scatter_values[j];
                    }
                } else { //if not the source rank
                    //send the local value
                    MPI_Send(&z, 1, MPI_INT, i, 0, comm);
                    //allocate memory for a temporary array
                    double* local_values_i = (double *) malloc(z * sizeof(double));
                    if (local_values_i == NULL) printf("local_values_i is NULL on rank %d\n", rank);
                    //store the scatterred values in the local array
                    for (int j = i * (m + 1), k = 0; k < z; j++, k++) { // j < ((i+1) * m + 1) && 
                        local_values_i[k] = scatter_values[j];
                    }
                    //send the local array to the correct processor
                    MPI_Send(local_values_i, z, MPI_DOUBLE, i, 1, comm);
                }
            } else { //if a rank doesn't need an extra element
                //implement the same process as above but with m instead of z
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
    } else { //if not the source rank
        //recieve n_local and the scattered array
        MPI_Status stat_n, stat_vals;
        MPI_Recv(&n_local, 1, MPI_INT, source_rank, 0, comm, &stat_n);
        local_values = (double *) malloc(n_local * sizeof(double));
        MPI_Recv(local_values, n_local, MPI_DOUBLE, source_rank, 1, comm, &stat_vals);
    }
    MPI_Barrier(comm);
}

double broadcast(double value, int source_rank, const MPI_Comm comm){
    //Set size and rank
    int rank , size, rank2;
    MPI_Comm_size (comm , &size );
    MPI_Comm_rank (comm , &rank );
    //create a temporary value
    double val;
    //set to value if source, zero otherwise
    if(rank == source_rank){
        val = value;
    } else{
        val = 0;
    }
    //find the amound of steps
    int d = ceil(log2(size));
    //hypercubic communication
    for(int i = 0; i < d; i++){
        //use XOR to get the rank' to send/recieve from
        rank2 = rank^(1<<(i));
        //if rank' exists
        if((rank2 < size)){
            //send and recieve current value
            MPI_Send(&val,1,MPI_DOUBLE,rank2,111,comm);
            MPI_Status stat;
            MPI_Recv(&value,1,MPI_DOUBLE,rank2,111,comm,&stat);
            //if the recieved value isn't 0, update it
            if(val == 0){
                val = value;
            }
        }
        MPI_Barrier(comm);
    }
    //return the correct value
    return val;
}

void parallel_prefix(const int n, const double* values, double* prefix_results, const int OP, const MPI_Comm comm){
    //Set size and rank
    int rank , size, rank2, total;
    MPI_Comm_size (comm , &size );
    MPI_Comm_rank (comm , &rank );
    //set the initial value
    prefix_results[0] = values[0];
    //if the operator is addition get initial array
    if (OP == PREFIX_OP_SUM) {
        for (int i = 1; i < n; i++) {
            prefix_results[i] = prefix_results[i - 1] + values[i];
        }
    } else if (OP == PREFIX_OP_PRODUCT) { //get initial array for multiplication operator
        for (int i = 1; i < n; i++) {
            prefix_results[i] = prefix_results[i - 1] * values[i];
        }
    }
    //set the total value for each array
    prefix_results[n] = prefix_results[n - 1];

    //get the amount of steps
    int d = ceil(log2(size));
    //Hypercubic communication
    for(int i = 0; i < d; i++){
        //use XOR to get the rank' to communicate with
        rank2 = rank^(1<<(i));
        //get the total sum/product
        total = prefix_results[n];
        //if rank' is a valid processor
        if((rank2 < size)){
            //send and recieve the total from it
            MPI_Send(&total, 1, MPI_DOUBLE, rank2, 111, comm );
            MPI_Status stat;
            MPI_Recv(&total, 1, MPI_DOUBLE, rank2, MPI_ANY_TAG, comm, &stat);
            //addition operator
            if (OP == PREFIX_OP_SUM) {
                //if recieving from a smaller rank add recieved total to prefix and total values
                if (rank > rank2) {
                    for (int i = 0; i <= n; i++) {
                        prefix_results[i] += total;
                    }
                } else {
                    //if recieving from a larger rank add recieved total to total
                        prefix_results[n] += total;
                }
            } else if (OP == PREFIX_OP_PRODUCT) { //multiplication operator
                //if recieving from a smaller rank multiple recieved total to prefix and total values
                if (rank > rank2) {
                    for (int i = 0; i <= n; i++) {
                        prefix_results[i] *= total;
                    }
                } else {
                    //if recieving from a larger rank multiply recieved total with total
                    prefix_results[n] *= total;
                }
            }
        }
    }
}

double mpi_poly_evaluator(const double x, const int n, const double* constants, const MPI_Comm comm){
    //Set rank
    int rank;
    MPI_Comm_rank (comm , &rank );
    //set answer value
    double answer = 0;
    //alocate memory for prevec and prefix results arrays
    double* pre_vec = (double *) malloc(n * sizeof(double));
    double* prefix_results = (double *) malloc((n + 1) * sizeof(double));
    if (pre_vec == NULL) { printf("pre_vec is null on %d", rank); }
    if (prefix_results == NULL) { printf("prefix_results is null on %d", rank); }
    //if the rank is zero make the first element x^0
    if(rank == 0){
        pre_vec[0] = 1;
    } else { //otherwise make the first element x
        pre_vec[0] = x;
    }
    //make all other values x
    for (int i = 1; i < n; i++) {
        pre_vec[i] = x;
    }
    //run parallel prefix on the x values with multiplication as the operator
    parallel_prefix(n, pre_vec, prefix_results, PREFIX_OP_PRODUCT, comm);

    //set the new prevec array with the x prefix results times the a_i constants
    for (int i = 0; i < n; i++) {
        pre_vec[i] = prefix_results[i] * constants[i];
    }
    //run parallel prefix on the a_i x^i array
    parallel_prefix(n, pre_vec, prefix_results, PREFIX_OP_SUM, comm);
    //the final value
    answer = prefix_results[n];
    //free the allocated memory
    free(pre_vec);
    free(prefix_results);
    //return the answer
    return answer;
}

