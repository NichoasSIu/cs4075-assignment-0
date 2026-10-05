#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* Function prototypes */
int Read_input(int* data_count_p, double* a_p, double* b_p, int* bin_count_p, int my_rank, MPI_Comm comm);

void Find_distribution(int data_count, int comm_sz, int send_counts[], int displacements[]);

void Generate_measurements(double data[], int data_count, double a, double b);

void Count_local_bins( double local_data[], int local_data_count, int local_bin_counts[], int bin_count, double a, double b);

void Reduce_and_print(double data[], int data_count,int local_bin_counts[],int bin_count,double a,double b,int my_rank,MPI_Comm comm);


int main(void) {
    int my_rank;
    int comm_sz;
    int valid_input;

    int data_count;
    int bin_count;
    int local_data_count;

    double a;
    double b;

    double* data = NULL;
    double* local_data = NULL;

    int* send_counts = NULL;
    int* displacements = NULL;
    int* local_bin_counts = NULL;

    MPI_Comm comm;

    // Start MPI 
    MPI_Init(NULL, NULL);

    comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &comm_sz);

    /* Read and broadcast the input */
    valid_input = Read_input(&data_count, &a, &b, &bin_count, my_rank, comm);

    if (!valid_input) {
        if (my_rank == 0) {
            printf("Invalid input.\n");
        }

        MPI_Finalize();
        return 1;
    }

    /*
     * Every process calculates how many measurements it receives.
     * Some processes may receive one additional measurement.
     */
    local_data_count = data_count / comm_sz + (my_rank < data_count % comm_sz);

    /*
     * Allocate the current process's local arrays.
     * Allocate at least one position even if local_data_count is zero.
     */
    local_data = malloc((local_data_count > 0 ? local_data_count : 1) * sizeof(double));

    local_bin_counts = calloc(bin_count, sizeof(int));

    /*
     * Only process 0 needs the complete data and distribution arrays.
     */
    if (my_rank == 0) {
        data = malloc(data_count * sizeof(double));
        send_counts = malloc(comm_sz * sizeof(int));
        displacements = malloc(comm_sz * sizeof(int));

        Find_distribution(data_count, comm_sz, send_counts, displacements);

        Generate_measurements(data, data_count, a, b);
    }

    /*
     * Distribute the measurements from process 0.
     */
    MPI_Scatterv(data, send_counts, displacements, MPI_DOUBLE, local_data, local_data_count, MPI_DOUBLE, 0, comm);

    /*
     * Each process counts its local measurements.
     */
    Count_local_bins(local_data, local_data_count, local_bin_counts, bin_count, a, b);

    /*
     * Combine the local counts and print the result.
     */
    Reduce_and_print(data, data_count, local_bin_counts, bin_count, a, b, my_rank, comm);

    //Free memory
    free(local_bin_counts);
    free(local_data);

    /* Free memory used only by process 0 */
    if (my_rank == 0) {
        free(displacements);
        free(send_counts);
        free(data);
    }

    MPI_Finalize();
    return 0;
}

/*
 * Process 0 reads the input.
 * MPI_Bcast sends the input to every process.
 */
int Read_input(int* data_count_p, double* a_p, double* b_p, int* bin_count_p, int my_rank, MPI_Comm comm){

    int valid_input = 1;

    *data_count_p = 0;
    *a_p = 0.0;
    *b_p = 0.0;
    *bin_count_p = 0;

    if (my_rank == 0) {
        printf("Enter number of measurements: ");

        if (scanf("%d", data_count_p) != 1) {
            valid_input = 0;
        }

        printf("Enter lower bound a: ");

        if (scanf("%lf", a_p) != 1) {
            valid_input = 0;
        }

        printf("Enter upper bound b: ");

        if (scanf("%lf", b_p) != 1) {
            valid_input = 0;
        }

        printf("Enter number of bins: ");

        if (scanf("%d", bin_count_p) != 1) {
            valid_input = 0;
        }

        if (
            *data_count_p <= 0 ||
            *bin_count_p <= 0 ||
            *b_p <= *a_p
        ) {
            valid_input = 0;
        }
    }

    /* Broadcast the input status */
    MPI_Bcast(&valid_input, 1, MPI_INT, 0, comm);

    /* Broadcast the input values */
    MPI_Bcast(data_count_p, 1, MPI_INT, 0, comm);

    MPI_Bcast(a_p, 1, MPI_DOUBLE, 0, comm);

    MPI_Bcast(b_p, 1, MPI_DOUBLE, 0, comm);

    MPI_Bcast(bin_count_p, 1, MPI_INT, 0, comm);

    return valid_input;
}


/*-------------------------------------------------------------------*/
/*
 * Calculates how many measurements each process receives
 * and where each process's section begins.
 */
void Find_distribution(int data_count, int comm_sz, int send_counts[], int displacements[]){

    int i;
    int base_count;
    int remainder;

    base_count = data_count / comm_sz;
    remainder = data_count % comm_sz;

    for (i = 0; i < comm_sz; i++) {
        send_counts[i] =
            base_count + (i < remainder);

        if (i == 0) {
            displacements[i] = 0;
        } else {
            displacements[i] =
                displacements[i - 1] +
                send_counts[i - 1];
        }
    }
}


/*-------------------------------------------------------------------*/
/*
 * Process 0 generates random measurements from a up to b.
 */
void Generate_measurements(double data[], int data_count, double a, double b){

    int i;

    /*
     * A fixed seed makes repeated tests produce the same sequence
     * on the same system.
     */
    srand(1);

    for (i = 0; i < data_count; i++) {
        data[i] =
            a +
            (b - a) *
            ((double) rand() / ((double) RAND_MAX + 1.0));
    }
}


/*-------------------------------------------------------------------*/
/*
 * Each process places its local measurements into histogram bins.
 */
void Count_local_bins( double local_data[], int local_data_count, int local_bin_counts[], int bin_count, double a, double b){

    double bin_width;
    int bin;
    int i;

    bin_width = (b - a) / bin_count;

    for (i = 0; i < local_data_count; i++) {
        /*
         * Calculate the bin containing this measurement.
         */
        bin = (int) ((local_data[i] - a) / bin_width);

        /*
         * If a value equals b, put it in the final bin.
         */
        if (bin == bin_count) {
            bin = bin_count - 1;
        }

        if (bin >= 0 && bin < bin_count) {
            local_bin_counts[bin]++;
        }
    }
}


/*-------------------------------------------------------------------*/
/*
 * MPI_Reduce adds corresponding bin counts from every process.
 * Process 0 prints the measurements, ranges, and counts.
 */
void Reduce_and_print(double data[], int data_count,int local_bin_counts[],int bin_count,double a,double b,int my_rank,MPI_Comm comm){

    int* global_bin_counts = NULL;
    double bin_width;
    double bin_low;
    double bin_high;
    int i;

    if (my_rank == 0) {
        global_bin_counts = malloc(
            bin_count * sizeof(int)
        );
    }

    MPI_Reduce(local_bin_counts, global_bin_counts, bin_count, MPI_INT, MPI_SUM, 0, comm);

    if (my_rank == 0) {
        printf("\nMeasurements:\n");

        for (i = 0; i < data_count; i++) {
            printf("%.2f ", data[i]);
        }

        printf("\n\nHistogram:\n");

        bin_width = (b - a) / bin_count;

        for (i = 0; i < bin_count; i++) {
            bin_low = a + i * bin_width;
            bin_high = bin_low + bin_width;

            if (i == bin_count - 1) {
                printf(
                    "Bin %d [%.2f, %.2f]: %d\n",
                    i,
                    bin_low,
                    bin_high,
                    global_bin_counts[i]
                );
            } else {
                printf(
                    "Bin %d [%.2f, %.2f): %d\n",
                    i,
                    bin_low,
                    bin_high,
                    global_bin_counts[i]
                );
            }
        }

        free(global_bin_counts);
    }
}