#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <fstream>
#include <string>
#include <mpi.h>


unsigned overlapMC(const double x2, const double R1, const double R2, size_t n, int rank, int procs)
{
	unsigned pts_inside = 0;

	std::default_random_engine g(42 + rank);  // random generator with seed 42
	std::uniform_real_distribution<double> u; // uniform distribution in [0, 1]

	// TODO_b: split the amount of work as equally as possible for each process.
	size_t n_local_size = n/procs;
	size_t n_start = 0;
	size_t n_end   = n;

	n_start = n_local_size*rank;
	if (rank+1==procs) {
		n_end = n;
	}
	else n_end = n_local_size*(rank+1);


	double box_width = 2*(R1+R2)-(R1-(x2-R2));
	double box_height = std::max(R1, R2);

	for (size_t i = n_start; i < n_end; ++i)
	{
		// TODO_a: implement the MC integration part here!
		std::pair<double, double> point (box_width*u(g)-R1, 2*box_height*u(g)-box_height);
		bool in_one = std::sqrt(std::pow(point.first, 2)+ std::pow(point.second, 2)) <= R1;
		bool in_two = std::sqrt(std::pow(12-point.first, 2)+ std::pow(point.second, 2)) <= R2;
		
		if(in_one && in_two) pts_inside++;
	}

	return  pts_inside;
}




int main(int argc, char *argv[])
{
	// TODO_b: Start-up the MPI environment and determine this process' rank ID as
    // well as the total number of processes (=ranks) involved in the
    // communicator

    int rank, procs;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);


	const double R1 = 5.0;		// Radius of first circle
	const double R2 = 10.0;		// Radius of second circle
	const double x2 = 12.0;		// x2 coordinate of second circle center

	// TODO_a: calculate the rectangle area for which you uniformly sample x & y
	const double area_rectangle = (2*(R1+R2)-(R1-(x2-R2)))*2*std::max(R1, R2);

	size_t n = 1e9 + 1;// default number of MC samples

	double ref = 17.0097776; // reference solution

	double t0 = MPI_Wtime(); 

	unsigned local_sum = overlapMC(x2, R1, R2, n, rank, procs);

	unsigned global_sum = (procs == -1) ? local_sum : 0;

	// TODO_b: Sum up and average all the local_sums to the master ranks global_sum using MPI_Reduce


	double area = global_sum / double(n) * area_rectangle;

	double t1 = MPI_Wtime();

	if(rank == 0){
		double error = std::abs(area - ref);
		if(error > 1e-2){
			printf("ERROR: you should get pi, but got: %.20f\n", area);
		}
		else{
			printf("result:  %.20f\nref: %.20e\ntime: %.20f\n", area, ref, t1 - t0);

			std::string file_name = "out/";
			file_name += std::to_string(procs);
			file_name += ".txt";
			//output time in a file
			std::ofstream myfile;
			myfile.open (file_name);
			myfile << procs << " " << t1 - t0 << "\n";
			myfile.close();
		}
	}

	// TODO_b: Shutdown the MPI environment
	MPI_Finalize();

	return 0;
}
