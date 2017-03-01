#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

void scatter(const int n, double* scatter_values, int &n_local, double* &local_values, int source_rank, const MPI_Comm comm){
    //Implementation
    int rank , size;
    MPI_Comm_size(comm , &size);
    MPI_Comm_rank(comm , &rank);
    if (rank == source_rank) {
        int m = int(ceil(n / size));
        for (int i = 0; i < size; i++) {
            if (i == rank) {
                n_local = m;

                local_values = (double *) malloc(m * sizeof(double));
                // assert(local_values != NULL);
                if (local_values == NULL) printf("local_values is NULL on rank %d\n", rank);
                for (int j = i * m, k = 0; j < ((i+1) * m) && k < m; j++, k++) {
                    local_values[k] = scatter_values[j];
                }
            } else {
                MPI_Send(&m, 1, MPI_INT, i, 0, comm);

                double* local_values_i = (double *) malloc(n_local * sizeof(double));
                // assert(local_values_i != NULL);
                if (local_values_i == NULL) printf("local_values_i is NULL on rank %d\n", rank);
                for (int j = i * m, k = 0; j < ((i+1) * m) && k < m; j++, k++) {
                    local_values_i[k] = scatter_values[j];
                }
                MPI_Send(local_values_i, m, MPI_DOUBLE, i, 1, comm);
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
	int rank , size;
	MPI_Comm_size (comm , &size );
	MPI_Comm_rank (comm , &rank );

	int d = log2(size+1);
	for(int i = 0; i < d; i++){
		if(rank && 1<<(i)){
			if((rank^(1<<(i))) < size){
				MPI_Send(&value,1,MPI_DOUBLE,(rank^(1<<(i))),111,comm );
			}
		} else{
			MPI_Status stat;
			MPI_Recv(&value,1,MPI_DOUBLE,(rank^(1<<(i))),MPI_ANY_TAG,comm , &stat);
		}
		MPI_Barrier(comm);
	}
	printf("Rank %i has value %f\n",rank,value);
    return value;
}

// int main(int argc, char *argv[]) {
//     MPI_Init(&argc, &argv);
//     MPI_Comm comm = MPI_COMM_WORLD;
//     broadcast(5,0,comm);
//     MPI_Finalize();
//     return 0;
// }

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	const MPI_Comm comm = MPI_COMM_WORLD;
	double* loc_arr;
	int n_local, source = 0;
	const int n = 10;
	// double scatter_vals[10] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
    double* scatter_vals = (double *) malloc(10 * sizeof(double));
    for (int i = 0; i < 10; i++) {
        scatter_vals[i] = i;
    }
    double* local_values;
	scatter(n, scatter_vals, n_local, local_values, source, comm);
	MPI_Finalize();
    free(scatter_vals);
    free(local_values);
	return 0;
}

	// if(rank == source_rank){
	// 	for(int i = 0; i < size; i++){
	// 		if(i != source_rank){
	// 			MPI_Send(&value,1,MPI_INT,i,111,comm );
	// 		}
	// 	}
	// } else{
	// 	MPI_Status stat;
	// 	MPI_Recv(&value,1,MPI_INT,source_rank,MPI_ANY_TAG,comm , &stat);
		
	// }