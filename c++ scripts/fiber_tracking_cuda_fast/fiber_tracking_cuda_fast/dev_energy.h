__device__ float dev_Edata
(
	int nj,

	float* in_x, float* in_y, float* in_z,
	int* in_cc, int* in_c,

	float* pos_x, float* pos_y, float* pos_z,
	float* T0, float* T1, float* T2, float* T3, float* T4, float* T5,
	float* Emin, float* Emax, float* delta_E,
	int* sig,
	int* nc,
	int* n_id, int* n,

	float T, int id
)
{
	float E1, E2;
	float xij[3];

	xij[0] = in_x[id] - in_x[nj];
	xij[1] = in_y[id] - in_y[nj];
	xij[2] = in_z[id] - in_z[nj];

	float norm = 1 / (xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2]);

	E1 = (    T0[id] * xij[0] * xij[0]
		+ 2 * T1[id] * xij[0] * xij[1]
		+ 2 * T2[id] * xij[0] * xij[2]
		+	  T3[id] * xij[1] * xij[1]
		+ 2 * T4[id] * xij[1] * xij[2]
		+	  T5[id] * xij[2] * xij[2]) * norm;

	E2 = (
		+     T0[nj] * xij[0] * xij[0]
		+ 2 * T1[nj] * xij[0] * xij[1]
		+ 2 * T2[nj] * xij[0] * xij[2]
		+     T3[nj] * xij[1] * xij[1]
		+ 2 * T4[nj] * xij[1] * xij[2]
		+     T5[nj] * xij[2] * xij[2]) * norm;

	return  (E1 - Emin[id]) / (Emax[id] - Emin[id]) + (E2 - Emin[nj]) / (Emax[nj] - Emin[nj]);
}

__device__ float dev_Eint
(
	int nj, int nk,

	float* in_x, float* in_y, float* in_z,
	int* in_cc, int* in_c,

	float* pos_x, float* pos_y, float* pos_z,
	float* T0, float* T1, float* T2, float* T3, float* T4, float* T5,
	float* Emin, float* Emax, float* delta_E,
	int* sig,
	int* nc,
	int* n_id, int* n,

	float T, int id
)
{
	float xij[3];
	xij[0] = in_x[nj] - in_x[id];
	xij[1] = in_y[nj] - in_y[id];
	xij[2] = in_z[nj] - in_z[id];

	float xik[3];
	xik[0] = in_x[nk] - in_x[id];
	xik[1] = in_y[nk] - in_y[id];
	xik[2] = in_z[nk] - in_z[id];

	float norm_ij = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
	float norm_ik = xik[0] * xik[0] + xik[1] * xik[1] + xik[2] * xik[2];
	float norm = sqrtf(norm_ij*norm_ik);
	float dot = xij[0] * xik[0] + xij[1] * xik[1] + xij[2] * xik[2];

	float cos = dot / (norm + 0.0000001);

	return (1 + cos) *0.5;
}

__device__ float dev_Ei_x
(
	float x, float y, float z,
	
	float* in_x, float* in_y, float* in_z,
	int* in_cc, int* in_c,

	float* pos_x, float* pos_y, float* pos_z,
	float* T0, float* T1, float* T2, float* T3, float* T4, float* T5,
	float* Emin, float* Emax, float* delta_E,
	int* sig,
	int* nc,
	int* n_id, int* n,

	float T, int id
)
{
	int count = 0;
	float E = 0;

	for (int j = n_id[id]; j < n_id[id] + nc[id]; j++)
	{
		E += in_c[j] * dev_Edata
			(
				n[j],

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
		count += in_c[j];
	}

	for (int j = n_id[id];	j < n_id[id] + nc[id]; j++)
	for (int k = j + 1;		k < n_id[id] + nc[id]; k++)
	{
		E += in_c[j] * in_c[k] * dev_Eint
			(
				n[j], n[k],

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

		count += in_c[j] * in_c[k];
	}

	for (int j = n_id[id]; j < n_id[id] + nc[id]; j++)
	for (int k = n_id[n[j]]; k < n_id[n[j]] + nc[n[j]]; k++)
	if (n[k] != id)
	{
		E += in_c[j] * in_c[k] * dev_Eint
			(
				n[j], n[k],

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

		count += in_c[j] * in_c[k];
	}

	return E / (count + 0.0000001);
}

__device__ float dev_Ei_c
(
	int d,
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
	float E = 0;
	if (d > 0)
	{
		int cc = in_cc[id] - 1 + 2 * in_c[d];

		E += (1 - sig[id])*fabsf(cc - 2) / (nc[id] - 2 + 0.0000001);
		E += sig[id]*      fabsf(cc - 1) / (nc[id] - 1 + 0.0000001);

		int cj;
		for (int j = n_id[id]; j < n_id[id] + nc[id]; j++)
		{
			if (j == d)
			{
				cj = in_cc[n[j]] + 1 - 2 * in_c[d];
				E += (1 - sig[n[j]])*	fabsf(cj - 2) / (nc[n[j]] - 2 + 0.0000001);
				E +=      sig[n[j]] *	fabsf(cj - 1) / (nc[n[j]] - 1 + 0.0000001);
			}
			else
			{
				cj = in_cc[n[j]];
				E += (1 - sig[n[j]])*	fabsf(cj - 2) / (nc[n[j]] - 2 + 0.0000001);
				E +=      sig[n[j]] *	fabsf(cj - 1) / (nc[n[j]] - 1 + 0.0000001);
			}
		}

		return E / (nc[id] + 1);
	}
	else
	{
		E += (1 - sig[id])*	fabsf(in_cc[id] - 2)	/ (nc[id] - 2 + 0.0000001);
		E += sig[id] *		fabsf(in_cc[id] - 1)	/ (nc[id] - 1 + 0.0000001);

		int cj;
		for (int j = n_id[id]; j < n_id[id] + nc[id]; j++)
		{
			cj = in_cc[n[j]];

			E += (1 - sig[n[j]])* fabsf(cj - 2) / (nc[n[j]] - 2 + 0.0000001);
			E +=	  sig[n[j]] * fabsf(cj - 1) / (nc[n[j]] - 1 + 0.0000001);
		}

		return E / (nc[id] + 1);
	}
}