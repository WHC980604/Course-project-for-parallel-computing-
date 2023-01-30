#include "head.H"

void UpdateBC
(
	struct Problem& pb
)
{
	int N = pb.N;
	#pragma omp parallel for collapse(2)
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			int k = i * N + j;
			if( i == 0 )
			{
				pb.omega[k] = 3.0 * pb.psi[k + N] / (pb.dh * pb.dh) - 0.5 * pb.omega[k + N];
			}
			if( i == N - 1 )
			{
				pb.omega[k] = 3.0 * pb.psi[k - N] / (pb.dh * pb.dh) - 0.5 * pb.omega[k - N];
			}
			if( j == 0 )
			{
				pb.omega[k] = 3.0 * pb.psi[k + 1] / (pb.dh * pb.dh) - 0.5 * pb.omega[k + 1];
			}
			if( j == N - 1 )
			{
				pb.omega[k] = 3.0 * (pb.psi[k - 1] + pb.dh) / (pb.dh * pb.dh) - 0.5 * pb.omega[k - 1];
			}
		}	
	}
}

void ComputeVelocity
(
	struct Problem& pb
)
{
	int N = pb.N;
	#pragma omp parallel for collapse(2)
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			int k = i * N + j;
			if( i == 0 || i == pb.N - 1 )
			{
				pb.u[k] = 0.0;
				pb.v[k] = 0.0;
			}
			if( j == 0 )
			{
				pb.u[k] = 0.0;
				pb.v[k] = 0.0;
			}
			if( j == pb.N - 1 )
			{
				pb.u[k] = 1.0;
				pb.v[k] = 0.0;
			}
			else
			{
				pb.u[k] = 0.5 * (pb.psi[k + 1] - pb.psi[k - 1]) / pb.dh;
				pb.v[k] = 0.5 * (pb.psi[k - N] - pb.psi[k + N]) / pb.dh;
			}
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
	#pragma omp parallel for reduction(+:sum)
	for( int i = 0 ; i < length ; i++ )
	{
	      sum = sum + vec[i] * vec[i];
	}
	return sqrt(sum);
}

void Solver
(
	struct Problem& pb
)
{
	const double residualTolerance = 1.e-6 ;
	const int    maxIterations     = 1000000 ;
	const int N			= pb.N;
	const double dt		= pb.dt; 
	const double dh		= pb.dh;
	const double Re		= pb.Re;

	double timeRS = 0.0;
	double timeBC = 0.0;
	double timeVelocity = 0.0;
	double timeResidual = 0.0;
	double start, end;
	
	double bP = -4.0 / (pb.dh * pb.dh); 
	double errR = 1.0;
	double errS = 1.0;
	
	for(int nIters = 0; nIters <= maxIterations; nIters++)
	{
		//Inner grid iteration
		start = omp_get_wtime();
		#pragma omp parallel for collapse(2) 
		for(int i = 1; i < N - 1; i++)
		{
			for(int j = 1; j < N - 1; j++)
			{
				int k = i * N + j;
				//calculate the vortex
				pb.R[k] = ((pb.omega[k + N] - 2. * pb.omega[k] + pb.omega[k - N]) / ((dh) * (dh)) + (pb.omega[k + 1] - 2. * pb.omega[k] + pb.omega[k - 1]) / ((dh) * (dh))) / Re - (pb.omega[k + N] * pb.u[k + N] - pb.omega[k - N] * pb.u[k - N]) / (2. * dh) - (pb.omega[k + 1] * pb.v[k + 1] - pb.omega[k - 1] * pb.v[k - 1]) / (2. * dh);
				pb.omega[k] += dt * pb.R[k];
				
				//calculate the Stream Function
				pb.S[k] = pb.omega[k] - (pb.psi[k + N] - 2. * pb.psi[k] + pb.psi[k - N]) / (dh * dh) - (pb.psi[k + 1] - 2. * pb.psi[k] + pb.psi[k - 1]) / (dh * dh);
				pb.psi[k] += (0.8 / bP) * pb.S[k];
			}
		}
		end = omp_get_wtime();
		timeRS += (end - start);
		
		//Update the Vortex boundary condition
		start = omp_get_wtime();
		UpdateBC(pb);
		end = omp_get_wtime();
		timeBC += (end - start);
		
		//Calculate the velocity field
		start = omp_get_wtime();
		ComputeVelocity(pb);
		end = omp_get_wtime();
		timeVelocity += (end - start);
		
		//Calculate residual
		start = omp_get_wtime();
		errR = ComputeResidual(N * N, pb.R);
		errS = ComputeResidual(N * N, pb.S);
		end = omp_get_wtime();
		timeResidual += (end - start);

		printf( "Total interations = %d, R residual = %.10e, S residual = %.10e \n", nIters, errR, errS ) ;
		printf("timeRS = %lf, timeBC = %lf, timeVelocity = %lf, timeResidual = %lf\n", timeRS, timeBC, timeVelocity, timeResidual);

		if(errR < residualTolerance && errS < residualTolerance && nIters > 2) break;
	}
}
