#include "head.H"
#include <mpi.h>

int main()
{
	MPI_Init(NULL, NULL);
	const int    N			= 200 ;
	const double Velocity_of_Upwall	= 1.0;
	const double edgeLength		= 1.0 ;
	const double Re			= 400.;
	const double dt			= 1.0e-3;

	double begin_time, end_time;
	begin_time = MPI_Wtime();
	
	struct Problem pb;
	
	initializeProblem(pb, Velocity_of_Upwall, edgeLength, Re, N, dt);
	 
	Solver(pb);

	end_time = MPI_Wtime();
	
	outputSolution(pb);

	freeMemories(pb);
	
	MPI_Finalize();
	printf("Total time = %lf\n", end_time - begin_time);

	return 0;
}
