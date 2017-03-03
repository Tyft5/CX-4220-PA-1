#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "const.h"

/*
 * TEST FILE- test each of the programs
 * 
 */

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
    printf("Rank %i was given the values %f", rank, local_values[0]);
    for (int i = 1; i < n_local; i++) {
    	printf(", %f", local_values[i]);
    }
    printf("\n");
}

double broadcast(double value, int source_rank, const MPI_Comm comm){
    //Implementation
	int rank , size, rank2;
	MPI_Comm_size (comm , &size );
	MPI_Comm_rank (comm , &rank );
    double val;
    // printf("rank %i  %f\n",rank,value);
    if(rank == source_rank){
        val = value;
    } else{
        val = 0;
    }
	int d = ceil(log2(size));
	for(int i = 0; i < d; i++){
        rank2 = rank^(1<<(i));
        if((rank2 < size)){
            //printf("rank: %i rank2: %i\n", rank,rank2);
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
    // int d = log2(size+1);
    // if(d != log2(size+1)){
    //     d+=1;
    // }
    int d = ceil(log2(size));
    for(int i = 0; i < d; i++){
        steps++;
        rank2 = rank^(1<<(i));
        total = prefix_results[n];
        if((rank2 < size)){
            //printf("rank: %i rank2: %i\n", rank,rank2);
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
        // if(rank == 5){
        //     printf("step %i total sum %f\n",steps, prefix_results[n]);
        // }
    }
    if (OP == PREFIX_OP_SUM) {
        for (int i = 1; i < n; i++) {
            prefix_results[i] = prefix_results[i - 1] + values[i];
        }
    } else if (OP == PREFIX_OP_PRODUCT) {
        for (int i = 1; i < n; i++) {
            prefix_results[i] = prefix_results[i - 1] * values[i];
        }
    }
    // printf("Rank %i has local sums %f %f\n",rank, prefix_results[0], prefix_results[1]);
}

double mpi_poly_evaluator(const double x, const int n, const double* constants, const MPI_Comm comm){
    //Implementation
    double answer = 0;
    double* pre_vec = (double *) malloc(n * sizeof(double));
    double* prefix_results = (double *) malloc((n + 1) * sizeof(double));
    int rank;
    MPI_Comm_rank (comm , &rank );
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

// int main(int argc, char *argv[]) {
//     MPI_Init(&argc, &argv);
//     int rank;
//     MPI_Comm comm = MPI_COMM_WORLD;
//     MPI_Comm_rank (comm , &rank );
//     int x;
//     int source = 3;
//     if(rank == source){
//         x = 5;
//     }
//     double val = broadcast(x,source,comm);
//     printf("Rank %i has value %f\n",rank,val);
//     MPI_Finalize();
//     return 0;
// }

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	const MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
	// double* loc_arr;
	int n_local, source = 2;
	const int n = 11;
	// double scatter_vals[8] = {0., 1., 2., 3., 4., 5., 6., 7.};
    double* scatter_vals = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        scatter_vals[i] = i;
    }
    double* local_values;
	scatter(n, scatter_vals, n_local, local_values, source, comm);
    // for (int i = 0; i < 2; i++) {
    //     printf("Rank %d: %f", rank, local_values[i]);
    // }
    free(scatter_vals);
    free(local_values);
    
    // double arr[2] = {2., 2.};
    // const double x = 2;
    //double* arr = (double *) malloc(8 * sizeof(double))
    // double* results = (double *) malloc(5 * sizeof(double));
    // parallel_prefix(2,arr,results,PREFIX_OP_SUM,comm);
    // double ans = mpi_poly_evaluator(x, 2, arr, comm);
    //printf("%f\n",ans);
    //broadcast(5,0,comm);
	MPI_Finalize();
    // free(scatter_vals);
    // free(local_values);
	return 0;
}

