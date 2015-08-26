__device__ void mc_c(vertex* in, vertex* out, float T, int seed)
{
	float E0, E1, p;

	curandState s;
	curand_init(seed, 0, 0, &s);

	out->cc = in->cc;
	for (int d = 0; d < in->nn; d++)
		out->c[d] = in->c[d];

	for (int d = 0; d < in->nn; d++)
	{
		E0 = dev_wc(T)*dev_Ei_c(out) + dev_wx(T)*dev_Ei_x(out, in->x, in->y, in->z);

		E1 = dev_wc(T)*dev_Ei_c(out, d) + dev_wx(T)*dev_Ei_x(out, in->x, in->y, in->z);

		p = fminf(1., expf((E0 - E1) / T));

		if (curand_uniform(&s) < p)
		{
			out->c[d] = 1 - in->c[d];
			out->cc += -1 + 2 * in->c[d];

			int j = 0;
			for (; j < in->n[d]->nn; j++)
			if (in == in->n[d]->n[j])
			{
				out->n[d]->c[j] = 1 - out->n[d]->c[j];
				break;
			}
			out->n[d]->cc += -1 + 2 * out->c[d];
		}
	}
}

__device__ void mc_x(vertex* in, vertex* out, float T, int seed)
{
	float x0, y0, z0, x, y, z;
	float E0, E1, p;

	curandState s;
	curand_init(seed, 0, 0, &s);

	x0 = in->x;
	y0 = in->y;
	z0 = in->z;

	for (int i = 0; i < dev_nx; i++)
	{
		E0 = dev_wx(T)*dev_Ei_x(in, x0, y0, z0);

		x = x0 - (x0 - in->pos_x + dev_pixdim_0 * 0.5) * dev_delta_x + curand_uniform(&s)* dev_delta_x;
		y = y0 - (y0 - in->pos_y + dev_pixdim_1 * 0.5) * dev_delta_x + curand_uniform(&s)* dev_delta_x;
		z = z0 - (z0 - in->pos_z + dev_pixdim_2 * 0.5) * dev_delta_x + curand_uniform(&s)* dev_delta_x;

		E1 = dev_wx(T)*dev_Ei_x(in, x, y, z);

		p = fminf(1., expf((E0 - E1) / T));

		if (curand_uniform(&s) < p)
		{
			x0 = x;
			y0 = y;
			z0 = z;

			out->x = x;
			out->y = y;
			out->z = z;
		}
	}
}

__global__ void mc(vertex* in, vertex* out,float T,int scount)
{
	int id = blockDim.x*blockIdx.x + threadIdx.x;

	if (id < scount)
	{
		mc_c(&in[id], &out[id], T, id);
		mc_x(&in[id], &out[id], T, id);
	}
}