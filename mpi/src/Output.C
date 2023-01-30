#include "head.H"

void outputSolution
(
    struct Problem&	pb
)
{
	//Output PLT file
	double w_out, PSI_out, u_out, v_out;
	FILE* op;
	op = fopen("output.plt", "w");
	fprintf(op, "variables=x,y,u,v,omega,psi\n zone i=%d,j=%d,f=point\n", pb.N, pb.N);
	for (int i = 0; i < pb.N; i++)
	{
		for (int j = 0; j < pb.N; j++)
		{
			int k = i * pb.N + j;
			w_out = pb.global_omega[k];
			PSI_out = pb.global_psi[k];
			u_out = pb.global_u[k];
			v_out = pb.global_v[k];
			fprintf(op, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", i * pb.dh, j * pb.dh, u_out, v_out, w_out, PSI_out);
		}
	}
	fclose(op);
}
