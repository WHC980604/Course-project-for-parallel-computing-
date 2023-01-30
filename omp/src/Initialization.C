#include "head.H"

void initializeProblem
(
	struct	Problem& pb,
	double	Velocity_of_Upwall,
	double	edgeLength,
	double	Re,
	int	N,
	double	dt
)
{
	pb.Velocity_of_Upwall = Velocity_of_Upwall;
	pb.L = edgeLength;
	pb.Re = Re;
	pb.N = N;
	pb.dt = dt;
	pb.dh = edgeLength / (N - 1.0);
	
	pb.omega	= (double*) calloc( N * N, sizeof(double) );
	pb.psi		= (double*) calloc( N * N, sizeof(double) );
	pb.R		= (double*) calloc( N * N, sizeof(double) );
	pb.S		= (double*) calloc( N * N, sizeof(double) );
	pb.u		= (double*) calloc( N * N, sizeof(double) );
	pb.v		= (double*) calloc( N * N, sizeof(double) );
}
