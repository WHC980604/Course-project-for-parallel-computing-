#include "head.H"

void UpdateBC
(
	struct Problem& pb
)
{
	int N = pb.N;
	for(int k = N; k < pb.localsize + N; k++)
	{
		int globalIndex = k + pb.left - N;
		if(globalIndex > 0 && globalIndex < N - 1)
		{
			pb.local_omega[k] = 3.0 * pb.local_psi[k + N] / (pb.dh * pb.dh) - 0.5 * pb.local_omega[k + N];
		}
		if(globalIndex > (N * N - N) && globalIndex < (N * N - 1))
		{
			pb.local_omega[k] = 3.0 * pb.local_psi[k - N] / (pb.dh * pb.dh) - 0.5 * pb.local_omega[k - N];
		}
		if(globalIndex % N == 0)
		{
			pb.local_omega[k] = 3.0 * pb.local_psi[k + 1] / (pb.dh * pb.dh) - 0.5 * pb.local_omega[k + 1];
		}
		if(globalIndex % N == (N - 1))
		{
			pb.local_omega[k] = 3.0 * (pb.local_psi[k - 1] + pb.dh) / (pb.dh * pb.dh) - 0.5 * pb.local_omega[k - 1];
		}
	}
}

void ComputeVelocity
(
	struct Problem& pb
)
{
	int N = pb.N;
	for(int k = N; k < pb.localsize + N; k++)
	{
		int globalIndex = k + pb.left - N;
		if((globalIndex >= 0 && k <= N - 1) || (globalIndex >= (N * N - N) && globalIndex <= (N * N - 1)) || globalIndex % pb.N == 0)
		{
			pb.local_u[k] = 0.0;
			pb.local_v[k] = 0.0;
		}
		if(globalIndex % N == (N - 1))
		{
			pb.local_u[k] = 1.0;
			pb.local_v[k] = 0.0;
		}
		else
		{
			pb.local_u[k] = 0.5 * (pb.local_psi[k + 1] - pb.local_psi[k - 1]) / pb.dh;
			pb.local_v[k] = 0.5 * (pb.local_psi[k - N] - pb.local_psi[k + N]) / pb.dh;
		}
	}
}


double ComputeResidual
(
	int length,
	double* vec
)
{
	double sum = 0;
	for( int i = 0 ; i < length ; i++ )
	{
	      sum = sum + vec[i] * vec[i];
	}
	return sum;
}

void Solver
(
	struct Problem& pb
)
{
	const double residualTolerance 	= 1.e-6 ;
	const int    maxIterations     	= 1000000 ;
	const int N			= pb.N;
	const double dt			= pb.dt; 
	const double dh			= pb.dh;
	const double Re			= pb.Re;
	int upRowID			= pb.upRowID;
	int lowRowID			= pb.lowRowID;	
	int left				= pb.left;
	int right				= pb.right;
	int localsize		= pb.localsize;
	int globalsize		= pb.globalsize;
	
	double bP = -4.0 / (pb.dh * pb.dh);
	double errR, errS, local_errR, local_errS;

	double timeRS = 0.0;
	double timeCom = 0.0;
	double timeBC = 0.0;
	double timeVelocity = 0.0;
	double timeResidual = 0.0;
	double start, end;

	int rank, size;
	int up, low;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;

	int* recvCounts		= (int*) calloc(size, sizeof(int));
	int* displacement	= (int*) calloc(size, sizeof(int));
	// get the receive counts for later gatherV //
	MPI_Allgather(&pb.localsize, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);
	// get the displacement for later gatherV //
	MPI_Allgather(&pb.left, 1, MPI_INT, displacement, 1, MPI_INT, MPI_COMM_WORLD);

	// Set process index of up and low
	if( rank == 0 )
	{
		up = MPI_PROC_NULL;
		low = rank + 1;
	}
	else if( rank == size - 1 )
	{
		up = rank - 1;
		low = MPI_PROC_NULL;
	}
	else
	{
		up = rank - 1;
		low = rank + 1;
	}

	for(int nIters = 0; nIters <= maxIterations; nIters++)
	{
		start = MPI_Wtime();
		// up to low comunicate
		MPI_Sendrecv(	&pb.local_omega[localsize],	N,	MPI_DOUBLE, 	low, 	0, 
				&pb.local_omega[0],		N,	MPI_DOUBLE,	up,	0,
				MPI_COMM_WORLD, 		&status					);
		MPI_Sendrecv(	&pb.local_psi[localsize],	N,	MPI_DOUBLE, 	low, 	1, 
				&pb.local_psi[0],		N,	MPI_DOUBLE,	up,	1,
				MPI_COMM_WORLD, 		&status					);
		MPI_Sendrecv(	&pb.local_u[localsize],		N,	MPI_DOUBLE, 	low, 	2, 
				&pb.local_u[0],			N,	MPI_DOUBLE,	up,	2,
				MPI_COMM_WORLD, 		&status					);
		MPI_Sendrecv(	&pb.local_v[localsize],		N,	MPI_DOUBLE, 	low, 	3, 
				&pb.local_v[0],			N,	MPI_DOUBLE,	up,	3,
				MPI_COMM_WORLD, 		&status					);
		// low to up comunicate
		MPI_Sendrecv(	&pb.local_omega[N],		N,	MPI_DOUBLE, 	up, 	0, 
				&pb.local_omega[localsize + N],	N,	MPI_DOUBLE,	low,	0,
				MPI_COMM_WORLD, 		&status					);
		MPI_Sendrecv(	&pb.local_psi[N],		N,	MPI_DOUBLE, 	up, 	1, 
				&pb.local_psi[localsize + N],	N,	MPI_DOUBLE,	low,	1,
				MPI_COMM_WORLD, 		&status					);
		MPI_Sendrecv(	&pb.local_u[N],			N,	MPI_DOUBLE, 	up, 	2, 
				&pb.local_u[localsize + N],	N,	MPI_DOUBLE,	low,	2,
				MPI_COMM_WORLD, 		&status					);
		MPI_Sendrecv(	&pb.local_v[N],			N,	MPI_DOUBLE, 	up, 	3, 
				&pb.local_v[localsize + N],	N,	MPI_DOUBLE,	low,	3,
				MPI_COMM_WORLD, 		&status					);
		end = MPI_Wtime();
		timeCom += (end - start);

		//if(rank == 0){printf("rank = %d\n", rank);}		

		//Inner grid iteration //
		//compute local R, S and omega, psi
		start = MPI_Wtime();
		for(int k = N; k < localsize + N; k++)
		{
			int globalIndex = k + pb.left - N;
			if( (globalIndex > N && globalIndex < (N * N - N - 1)) && (globalIndex % N != 0) && (globalIndex % N != (N - 1)))
			{
				// compute local R[i] and omega[i] //
				pb.local_R[k] 	= ((pb.local_omega[k + N] - 2. * pb.local_omega[k] + pb.local_omega[k - N]) / ((dh) * (dh)) 
						+ (pb.local_omega[k + 1] - 2. * pb.local_omega[k] + pb.local_omega[k - 1]) / ((dh) * (dh))) / Re 
						- (pb.local_omega[k + N] * pb.local_u[k + N] - pb.local_omega[k - N] * pb.local_u[k - N]) / (2. * dh) 
						- (pb.local_omega[k + 1] * pb.local_v[k + 1] - pb.local_omega[k - 1] * pb.local_v[k - 1]) / (2. * dh);
				pb.local_omega[k] += dt * pb.local_R[k];

				// compute local S[i] and psi[i] //
				pb.local_S[k] 	= pb.local_omega[k] 
						- (pb.local_psi[k + N] - 2. * pb.local_psi[k] + pb.local_psi[k - N]) / (dh * dh) 
						- (pb.local_psi[k + 1] - 2. * pb.local_psi[k] + pb.local_psi[k - 1]) / (dh * dh);
				pb.local_psi[k] += (0.8 / bP) * pb.local_S[k];
			}
		}
		end = MPI_Wtime();
		timeRS += (end - start);

		//Update the Vortex boundary condition //
		start = MPI_Wtime();
		UpdateBC(pb);
		end = MPI_Wtime();
		timeBC += (end - start);

		//Calculate the velocity field //
		start = MPI_Wtime();
		ComputeVelocity(pb);
		end = MPI_Wtime();
		timeVelocity += (end - start);

		//Calculate residual
		start = MPI_Wtime();
		local_errR = ComputeResidual(localsize, pb.local_R);
		local_errS = ComputeResidual(localsize, pb.local_S);
		MPI_Allreduce(&local_errR, &errR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&local_errS, &errS, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		errR = sqrt(errR);
		errS = sqrt(errS);
		end = MPI_Wtime();
		timeResidual += (end - start);
		
		if(rank == 0)
		{
			printf( "Total interations = %d, R residual = %.10e, S residual = %.10e \n", nIters, errR, errS ) ;
		}

		if(errR < residualTolerance && errS < residualTolerance && nIters > 2) break;
	}

	if( rank == 0 )
	{
		printf("timeRS = %lf, timeBC = %lf, timeVelocity = %lf, timeResidual = %lf, timeCom = %lf\n", timeRS, timeBC, timeVelocity, timeResidual, timeCom);
	}
	MPI_Allgatherv(&pb.local_omega[N], localsize, MPI_DOUBLE, pb.global_omega, recvCounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgatherv(&pb.local_psi[N], localsize, MPI_DOUBLE, pb.global_psi, recvCounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgatherv(&pb.local_u[N], localsize, MPI_DOUBLE, pb.global_u, recvCounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgatherv(&pb.local_v[N], localsize, MPI_DOUBLE, pb.global_v, recvCounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);
}
