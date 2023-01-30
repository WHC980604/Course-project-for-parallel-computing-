#ifndef HEAD_H
#define HEAD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// define structs
struct Problem
{
	int N;
	double Re;
	double L;
	double	Velocity_of_Upwall;

	int upRowID;
	int lowRowID;	
	int left;
	int right;
	int localsize;
	int globalsize;

	double* local_omega;
	double* local_psi;
	double* local_u;
	double* local_v;
	double* local_R;
	double* local_S;

	double* global_omega;
	double* global_psi;
	double* global_u;
	double* global_v;
	double* global_R;
	double* global_S;

	double dh;
	double dt;
};


//initialize the problem
void initializeProblem
(
	struct	Problem& problem,
	double	Velocity_of_Upwall,
	double	edgeLength,
	double	Re,
	int	N,
	double	dt
);

void Solver
(
	struct Problem& pb
);

void UpdateBC
(
	struct Problem& pb
);

void ComputeVelocity
(
	struct Problem& pb
);

double ComputeResidual
(
	int length,
	double* vec
);

void outputSolution
(
    struct Problem& pb
);

inline void freeMemories
(
    struct Problem& pb
)
{
    free(pb.global_omega);
    free(pb.global_psi);
    free(pb.global_u);
    free(pb.global_v);
    free(pb.global_R);
    free(pb.global_S);
    free(pb.local_omega);
    free(pb.local_psi);
    free(pb.local_u);
    free(pb.local_v);
    free(pb.local_R);
    free(pb.local_S);
}

#endif
