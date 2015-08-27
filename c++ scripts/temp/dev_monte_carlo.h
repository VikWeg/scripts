__device__ void mc_c(
		float* in_x, float* in_y, float* in_z, int* in_cc, int* in_c,

		float* out_x, float* out_y, float* out_z, int* out_cc, int* out_c,

		float* pos_x, float* pos_y, float* pos_z,
		float* T0, float* T1, float* T2, float* T3, float* T4, float* T5,
		float* Emin, float* Emax, float* delta_E,
		int* sig,
		int* nc,
		int* n_id, int* n,

		float T, int id
	)
	{
		float E0, E1;
		float p;

		curandState s;
		curand_init(id, 0, 0, &s);

		out_cc[id] = in_cc[id];
		for (int d = n_id[id]; d < n_id[id] + nc[id]; d++)
			out_c[d] = in_c[d];

		for (int d = n_id[id]; d < n_id[id] + nc[id]; d++)
		{
			E0 =	dev_wc(T)*dev_Ei_c(	
											-1,
											in_x, in_y, in_z, in_cc, in_c,

											out_x, out_y, out_z, out_cc, out_c,

											pos_x, pos_y, pos_z,
											T0, T1, T2, T3, T4, T5,
											Emin, Emax, delta_E,
											sig,
											nc,
											n_id, n,

											T, id
										)
				+	dev_wx(T)*dev_Ei_x(	
											in_x[id], in_y[id], in_z[id],

											in_x, in_y, in_z,
											in_cc, in_c,

											pos_x, pos_y, pos_z,
											T0, T1, T2, T3, T4, T5,
											Emin, Emax, delta_E,
											sig,
											nc,
											n_id, n,

											T, id
										);

			E1 =	dev_wc(T)*dev_Ei_c(	
											d,
											in_x, in_y, in_z, in_cc, in_c,

											out_x, out_y, out_z, out_cc, out_c,

											pos_x, pos_y, pos_z,
											T0, T1, T2, T3, T4, T5,
											Emin, Emax, delta_E,
											sig,
											nc,
											n_id, n,

											T, id
										) 
				+	dev_wx(T)*dev_Ei_x(
											in_x[id], in_y[id], in_z[id],

											in_x, in_y, in_z,
											in_cc, in_c,

											pos_x, pos_y, pos_z,
											T0, T1, T2, T3, T4, T5,
											Emin, Emax, delta_E,
											sig,
											nc,
											n_id, n,

											T, id
										);

			p = fminf(1., expf((E0 - E1) / T));

			if (curand_uniform(&s) < p)
			{
				out_c[d] = 1 - in_c[d];
				out_cc[id] += -1 + 2 * in_c[d];

				for (int j = n_id[n[d]]; j < n_id[n[d]] + nc[n[d]]; j++)
				if (id == n[j])
				{			
					out_c[j] = 1 - out_c[j];
					break;
				}
				out_cc[n[d]] += -1 + 2 * out_c[d];
			}
		}
	}

__device__ void mc_x(
		float* in_x, float* in_y, float* in_z, int* in_cc, int* in_c,

		float* out_x, float* out_y, float* out_z, int* out_cc, int* out_c,

		float* pos_x, float* pos_y, float* pos_z,
		float* T0, float* T1, float* T2, float* T3, float* T4, float* T5,
		float* Emin, float* Emax, float* delta_E,
		int* sig,
		int* nc,
		int* n_id, int* n,

		float T, int id
	)
	{
		float x0, y0, z0, x, y, z;
		float E0, E1, p;

		curandState s;
		curand_init(id, 0, 0, &s); //id used as seed;

		x0 = in_x[id];
		y0 = in_y[id];
		z0 = in_z[id];

		for (int i = 0; i < *dev_nx; i++)
		{
			E0 = dev_wx(T)*	dev_Ei_x
							(
								x0, y0, z0,
								
								in_x, in_y, in_z,
								in_cc, in_c,

								pos_x, pos_y, pos_z,
								T0, T1, T2, T3, T4, T5,
								Emin, Emax, delta_E,
								sig,
								nc,
								n_id, n,

								T, id
							);

			x = x0 - (x0 - pos_x[id] + dev_pixdim[0] * 0.5) * *dev_delta_x + curand_uniform(&s)* *dev_delta_x;
			y = y0 - (y0 - pos_y[id] + dev_pixdim[1] * 0.5) * *dev_delta_x + curand_uniform(&s)* *dev_delta_x;
			z = z0 - (z0 - pos_z[id] + dev_pixdim[2] * 0.5) * *dev_delta_x + curand_uniform(&s)* *dev_delta_x;

			E1 = dev_wx(T)*	dev_Ei_x(
								x, y, z,

								in_x, in_y, in_z,
								in_cc, in_c,

								pos_x, pos_y, pos_z,
								T0, T1, T2, T3, T4, T5,
								Emin, Emax, delta_E,
								sig,
								nc,
								n_id, n,

								T, id
							);

			p = fminf(1., expf((E0 - E1) / T));

			if (curand_uniform(&s) < p)
			{
				x0 = x;
				y0 = y;
				z0 = z;

				out_x[id] = x;
				out_y[id] = y;
				out_z[id] = z;
			}
		}
	}

__global__ void mc(	
		float* in_x, float* in_y, float* in_z, int* in_cc, int* in_c,

		float* out_x, float* out_y, float* out_z, int* out_cc, int* out_c,

		float* pos_x, float* pos_y, float* pos_z,
		float* T0, float* T1, float* T2, float* T3, float* T4, float* T5,
		float* Emin, float* Emax, float* delta_E,
		int* sig,
		int* nc,
		int* n_id, int* n,

		float T, int scount
		)
		{
			int id = blockDim.x*blockIdx.x + threadIdx.x;

			while (id < scount)
			{
				mc_c(
					in_x, in_y, in_z, in_cc, in_c,

					out_x, out_y, out_z, out_cc, out_c,

					pos_x, pos_y, pos_z,
					T0, T1, T2, T3, T4, T5,
					Emin, Emax, delta_E,
					sig,
					nc,
					n_id, n,
			
					T, id
				);

				__syncthreads();

				mc_x(
					in_x, in_y, in_z, in_cc, in_c,

					out_x, out_y, out_z, out_cc, out_c,

					pos_x, pos_y, pos_z,
					T0, T1, T2, T3, T4, T5,
					Emin, Emax, delta_E,
					sig,
					nc,
					n_id, n,

					T, id);

				__syncthreads();

				id += blockDim.x * gridDim.x;
			}
		}