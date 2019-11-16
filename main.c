#include<stdio.h>
#include<mpi.h>
#define MATRIX_ROWS 3
#define MATRIX_COLUMNS 3
#define VECTOR_ROWS 3
#define VECTOR_COLUMNS 1
#define MASTER_TO_SLAVE_TAG 1
#define SLAVE_TO_MASTER_TAG 4
void printArray();
int rank,size,i,j,k,low_bound,upper_bound,portion; 
float matrix[MATRIX_ROWS][MATRIX_COLUMNS],vector[VECTOR_ROWS][VECTOR_COLUMNS],result[MATRIX_ROWS][VECTOR_COLUMNS];
double start_time, end_time;
MPI_Status status;
MPI_Request request;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0) {
        printf("enter the values of the %dx%d matrix:\n",MATRIX_ROWS,MATRIX_COLUMNS);
        for(int i=0;i<MATRIX_ROWS;i++){
            for(int j=0;j<MATRIX_COLUMNS;j++){
                scanf("%f",&matrix[i][j]);
            }
        }
        printf("enter the values of the %dx%d vector:\n",VECTOR_ROWS,VECTOR_COLUMNS);
        for(int i=0;i<VECTOR_ROWS;i++){
            for(int j=0;j<VECTOR_COLUMNS;j++){
                scanf("%f",&vector[i][j]);
            }
        }
        start_time = MPI_Wtime();
        for (i = 1; i < size; i++) {
            portion = (MATRIX_ROWS / (size - 1));
            low_bound = (i - 1) * portion;
            if (((i + 1) == size) && ((MATRIX_ROWS % (size - 1)) != 0)) {
                upper_bound = MATRIX_ROWS;
            } else {
                upper_bound = low_bound + portion;
            }
            MPI_Isend(&low_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &request);
            MPI_Isend(&upper_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &request);
            MPI_Isend(&matrix[low_bound][0], (upper_bound - low_bound) * MATRIX_ROWS, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &request);
        }
        if(size==1){
        for (i = 0; i < MATRIX_ROWS; i++) {
            for (j = 0; j < VECTOR_COLUMNS; j++) {
                result[i][j]=0;
                for (k = 0; k < VECTOR_ROWS; k++) {
                    result[i][j] += (matrix[i][k] * vector[k][j]);
                }
            }
        }
    }
    }
    MPI_Bcast(&vector, VECTOR_ROWS*VECTOR_COLUMNS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank > 0) {
        printf("rank %d:",rank);
        MPI_Recv(&low_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&matrix[low_bound][0], (upper_bound - low_bound) * MATRIX_ROWS, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);
        for (i = low_bound; i < upper_bound; i++) {
            for (j = 0; j < VECTOR_COLUMNS; j++) {
                result[i][j]=0;
                for (k = 0; k < VECTOR_ROWS; k++) {
                    result[i][j] += (matrix[i][k] * vector[k][j]);
                }
            }
        }
        MPI_Isend(&low_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &request);
        MPI_Isend(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &request);
        MPI_Isend(&result[low_bound][0], (upper_bound - low_bound) * VECTOR_COLUMNS, MPI_DOUBLE, 0, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &request);
    }
    if (rank == 0) {
        for (i = 1; i < size; i++) {
            MPI_Recv(&low_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&result[low_bound][0], (upper_bound - low_bound) * VECTOR_COLUMNS, MPI_DOUBLE, i, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &status);
        }
        end_time = MPI_Wtime();
        printf("\nRunning Time = %f\n\n", end_time - start_time);
        for(i=0;i<MATRIX_ROWS;i++){
            printf("%f ",result[i][0]);
        }
        printf("\nsize = %d",size);
    }
    MPI_Finalize();
    return 0;
}
