#ifndef HEAD_H
#define HEAD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// define structs
struct Problem
{
	int N;
	double Re;
	double L;
	double	Velocity_of_Upwall;
	double* omega;
	double* psi;
	double* u;
	double* v;
	double* R;
	double* S;
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
    free(pb.omega);
    free(pb.psi);
    free(pb.u);
    free(pb.v);
    free(pb.R);
    free(pb.S);
}

#endif
