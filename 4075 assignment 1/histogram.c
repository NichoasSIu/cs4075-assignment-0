#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

//have process 0 read input data. and distribuite among process, also it must print out program
//create bins ranging from the given input, n = amount of bins, b-a/n will be the width



//function prototype
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

    //read input
    valid_input = Read_input(&data_count, &a, &b, &bin_count, my_rank, comm);

    //sanity check, program ends when input is invalid
    if (!valid_input) {
        if (my_rank == 0) {
            printf("Invalid input.\n");
        }

        MPI_Finalize();
        return 1;
    }

    //calculates how much shouyld it take per
    local_data_count = data_count / comm_sz + (my_rank < data_count % comm_sz);

    //allocate
    local_data = malloc((local_data_count > 0 ? local_data_count : 1) * sizeof(double));

    local_bin_counts = calloc(bin_count, sizeof(int));

    //process 0 does heavy lifinting thefore only malloc is needed for it
    if (my_rank == 0) {
        data = malloc(data_count * sizeof(double));
        send_counts = malloc(comm_sz * sizeof(int));
        displacements = malloc(comm_sz * sizeof(int));

        Find_distribution(data_count, comm_sz, send_counts, displacements);

        Generate_measurements(data, data_count, a, b);
    }

    //distribute measurements from process 0
    MPI_Scatterv(data, send_counts, displacements, MPI_DOUBLE, local_data, local_data_count, MPI_DOUBLE, 0, comm);

    //count local bins
    Count_local_bins(local_data, local_data_count, local_bin_counts, bin_count, a, b);

    //combine 
    Reduce_and_print(data, data_count, local_bin_counts, bin_count, a, b, my_rank, comm);

    //Free memory
    free(local_bin_counts);
    free(local_data);

    if (my_rank == 0) {
        free(displacements);
        free(send_counts);
        free(data);
    }

    MPI_Finalize();
    return 0;
}

// process 0 reads input and also distribuites it
int Read_input(int* data_count_p, double* a_p, double* b_p, int* bin_count_p, int my_rank, MPI_Comm comm){
    //saninty check when user inputs nothing it will return 0 else return 1 and end program
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

        if (*data_count_p <= 0 || *bin_count_p <= 0 || *b_p <= *a_p) {
            valid_input = 0;
        }
    }

    //ditribute proper measurements to every process coresposnding to its name upper, lower, etc
    MPI_Bcast(&valid_input, 1, MPI_INT, 0, comm);

    MPI_Bcast(data_count_p, 1, MPI_INT, 0, comm);

    MPI_Bcast(a_p, 1, MPI_DOUBLE, 0, comm);

    MPI_Bcast(b_p, 1, MPI_DOUBLE, 0, comm);

    MPI_Bcast(bin_count_p, 1, MPI_INT, 0, comm);

    return valid_input;
}


// Calculates how many measurements each process receives and where each process's section begins.
void Find_distribution(int data_count, int comm_sz, int send_counts[], int displacements[]){

    int i;
    int base_count;
    int remainder;
    
    //calcs how much per section, ex: data cout of 22 and comm_sz of 4, therefore it takes 5 increments
    base_count = data_count / comm_sz;
    remainder = data_count % comm_sz;

    for (i = 0; i < comm_sz; i++) {
        send_counts[i] = base_count + (i < remainder);

        //find where each process should start if index 0 its starts process 0
        //esle would be calc with previous index plus send_counts to find its index
        if (i == 0) {
            displacements[i] = 0;
        } else {
            displacements[i] = displacements[i - 1] + send_counts[i - 1];
        }
    }
}


//randomly generate numbers bewteen a-b
void Generate_measurements(double data[], int data_count, double a, double b){

    int i;

    srand(1);

    for (i = 0; i < data_count; i++) {
        data[i] = a + (b - a) * ((double) rand() / ((double) RAND_MAX + 1.0));
    }
}

//place process into the correct bin

void Count_local_bins( double local_data[], int local_data_count, int local_bin_counts[], int bin_count, double a, double b){

    double bin_width;
    int bin;
    int i;

    //calcs the width between interval
    bin_width = (b - a) / bin_count;

    for (i = 0; i < local_data_count; i++) {
        //sub lower bound and divide by width this gives the range away from value a
        //example: value 5.7 and the ranges are 2(0-2)(2-4)... this would be placed in (4-6) because its 5.7
        bin = (int) ((local_data[i] - a) / bin_width);

        //because index start 0 we cannot find last bin ex: bin 5 with ranges 8-10 therefore we explictly state it
        if (bin == bin_count) {
            bin = bin_count - 1;
        }

        //prevents it from accessing memory outside of bin counts
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

    //process 0 must print it out
    if (my_rank == 0) {
        global_bin_counts = malloc(bin_count * sizeof(int));
    }

    //combine all counts
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
                printf("Bin %d [%.2f, %.2f]: %d\n", i, bin_low, bin_high, global_bin_counts[i]);
            } else {
                printf("Bin %d [%.2f, %.2f): %d\n", i, bin_low, bin_high, global_bin_counts[i]);
            }
        }

        free(global_bin_counts);
    }
}