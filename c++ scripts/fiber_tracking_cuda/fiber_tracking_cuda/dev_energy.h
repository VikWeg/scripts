__device__ float dev_Edata(vertex* vi, vertex* vj, float x, float y, float z)
{
	float E1, E2;
	float xij[3];

	xij[0] = x - vj->x;
	xij[1] = y - vj->y;
	xij[2] = z - vj->z;

	float norm = 1 / (xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2]);

	E1 = (    (*vi).T0 * xij[0] * xij[0]
		+ 2 * (*vi).T1 * xij[0] * xij[1]
		+ 2 * (*vi).T2 * xij[0] * xij[2]
		+	  (*vi).T3 * xij[1] * xij[1]
		+ 2 * (*vi).T4 * xij[1] * xij[2]
		+	  (*vi).T5 * xij[2] * xij[2]) * norm;

	E2 = (
		+     (*vj).T0 * xij[0] * xij[0]
		+ 2 * (*vj).T1 * xij[0] * xij[1]
		+ 2 * (*vj).T2 * xij[0] * xij[2]
		+     (*vj).T3 * xij[1] * xij[1]
		+ 2 * (*vj).T4 * xij[1] * xij[2]
		+     (*vj).T5 * xij[2] * xij[2]) * norm;

	return  0.5*((E1 - vi->Emin) / (vi->Emax - E1 + 0.0000001) + (E2 - vj->Emin) / (vj->Emax - E2 + 0.0000001));
}

__device__ float dev_Eint(vertex* vj, vertex* vi, vertex* vk)
{
	float xij[3];
	xij[0] = vj->x - vi->x;
	xij[1] = vj->y - vi->y;
	xij[2] = vj->z - vi->z;

	float xik[3];
	xik[0] = vk->x - vi->x;
	xik[1] = vk->y - vi->y;
	xik[2] = vk->z - vi->z;

	float norm_ij = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
	float norm_ik = xik[0] * xik[0] + xik[1] * xik[1] + xik[2] * xik[2];
	float norm = sqrtf(norm_ij*norm_ik);
	float dot = xij[0] * xik[0] + xij[1] * xik[1] + xij[2] * xik[2];

	float cos = dot / norm;

	return (1 + cos) / (1.0001 - cos);
}

__device__ float dev_Eint(vertex* vj, vertex* vi, vertex* vk, float x, float y, float z)
{
	float xij[3];
	xij[0] = vj->x - x;
	xij[1] = vj->y - y;
	xij[2] = vj->z - z;

	float xik[3];
	xik[0] = vk->x - x;
	xik[1] = vk->y - y;
	xik[2] = vk->z - z;

	float norm_ij = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
	float norm_ik = xik[0] * xik[0] + xik[1] * xik[1] + xik[2] * xik[2];
	float norm = sqrtf(norm_ij*norm_ik);
	float dot = xij[0] * xik[0] + xij[1] * xik[1] + xij[2] * xik[2];

	float cos = dot / norm;

	return dev_wint(cos);
}

__device__ float dev_Ei_x(vertex* vi, float x, float y, float z)
{
	int count = 0;
	float E = 0;

	for (int j = 0; j < vi->nn; j++)
	{
		E += vi->c[j] * dev_Edata(vi, vi->n[j], x, y, z);
		count += vi->c[j];
	}

	for (int j = 0; j < vi->nn; j++)
	for (int k = j + 1; k < vi->nn; k++)
	{
		E += vi->c[j] * vi->c[k] * dev_Eint(vi->n[j], vi, vi->n[k], x, y, z);
		count += vi->c[j] * vi->c[k];
	}

	for (int j = 0; j < vi->nn; j++)
	for (int k = 0; k < vi->n[j]->nn; k++)
	if (vi->n[j]->n[k] != vi)
	{
		E += vi->c[j] * vi->n[j]->c[k] * dev_Eint(vi, vi->n[j], vi->n[j]->n[k], x, y, z);
		count += vi->c[j] * vi->n[j]->c[k];
	}

	return E / (count + 0.0000001);
}

__device__ float dev_Ei_c(vertex* vi)
{
	float E = 0;

	E += (1 - vi->sig)*fabsf(vi->cc - 2) / (vi->nn - 2 + 0.0000001);
	E += (vi->sig)*    fabsf(vi->cc - 1) / (vi->nn - 1 + 0.0000001);

	int cj;
	for (int j = 0; j < vi->nn; j++)
	{
		cj = vi->n[j]->cc;

		E += (1 - vi->n[j]->sig)*fabsf(cj - 2) / (float)(vi->n[j]->nn - 2 + 0.0000001);
		E += (vi->n[j]->sig)*    fabsf(cj - 1) / (float)(vi->n[j]->nn - 1 + 0.0000001);
	}

	return E / (vi->nn + 1);
}

__device__ float dev_Ei_c(vertex* vi, int d)
{
	float E = 0;

	int cc = vi->cc - 1 + 2 * vi->c[d];

	E += (1 - vi->sig)*fabsf(cc - 2) / (vi->nn - 2 + 0.0000001);
	E += (vi->sig)*    fabsf(cc - 1) / (vi->nn - 1 + 0.0000001);

	int cj;
	for (int j = 0; j < vi->nn; j++)
	{
		if (j == d)
		{
			cj = vi->n[j]->cc + 1 - 2 * vi->c[d];
			E += (1 - vi->n[j]->sig)*fabsf(cj - 2) / (float)(vi->n[j]->nn - 2 + 0.0000001);
			E += (vi->n[j]->sig)*    fabsf(cj - 1) / (float)(vi->n[j]->nn - 1 + 0.0000001);
		}
		else
		{
			cj = vi->n[j]->cc;
			E += (1 - vi->n[j]->sig)*fabsf(cj - 2) / (float)(vi->n[j]->nn - 2 + 0.0000001);
			E += (vi->n[j]->sig)*    fabsf(cj - 1) / (float)(vi->n[j]->nn - 1 + 0.0000001);
		}
	}

	return E / (vi->nn + 1);
}