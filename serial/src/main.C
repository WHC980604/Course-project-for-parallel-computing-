#include "head.H"
#include "omp.h"

int main()
{
	const int    N			= 200 ;
	const double Velocity_of_Upwall	= 1.0;
	const double edgeLength		= 1.0 ;
	const double Re			= 400.;
	const double dt			= 1.0e-3;

	double begin_time, end_time;
	begin_time = omp_get_wtime();
	
	struct Problem pb;
	
	initializeProblem(pb, Velocity_of_Upwall, edgeLength, Re, N, dt);
	 
	Solver(pb);

	end_time = omp_get_wtime();
	printf("Total time = %lf\n", end_time - begin_time);
	
	outputSolution(pb);

	freeMemories(pb);
	
	return 0;
}
