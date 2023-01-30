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
	pb.globalsize = N * N;

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Blocking //
	const double section	= static_cast<double>(pb.N) / static_cast<double>(size);
	pb.upRowID		= floor(section * rank);
	pb.lowRowID		= floor(section * (rank + 1)) - 1;
	//pb.localsize	= (pb.lowRowID - upRowID + 1) * pb.N;
	pb.left			= pb.upRowID * N;
	pb.right		= (pb.lowRowID + 1) * N - 1;
	pb.localsize	= pb.right - pb.left + 1;
	
	pb.local_omega	= (double*) calloc( pb.localsize + N * 2, sizeof(double) );
	pb.local_psi	= (double*) calloc( pb.localsize + N * 2, sizeof(double) );
	pb.local_R		= (double*) calloc( pb.localsize + N * 2, sizeof(double) );
	pb.local_S		= (double*) calloc( pb.localsize + N * 2, sizeof(double) );
	pb.local_u		= (double*) calloc( pb.localsize + N * 2, sizeof(double) );
	pb.local_v		= (double*) calloc( pb.localsize + N * 2, sizeof(double) );

	pb.global_omega	= (double*) calloc( pb.globalsize, sizeof(double) );
	pb.global_psi	= (double*) calloc( pb.globalsize, sizeof(double) );
	pb.global_R		= (double*) calloc( pb.globalsize, sizeof(double) );
	pb.global_S		= (double*) calloc( pb.globalsize, sizeof(double) );
	pb.global_u		= (double*) calloc( pb.globalsize, sizeof(double) );
	pb.global_v		= (double*) calloc( pb.globalsize, sizeof(double) );

}
