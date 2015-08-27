float Edata(int s, int nj)
{
	float E1, E2;
	float xij[3];

	xij[0] = x[s] - x[nj];
	xij[1] = y[s] - y[nj];
	xij[2] = z[s] - z[nj];

	float norm = 1 / (xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2]);

	E1 = (    T0[s] * xij[0] * xij[0]
		+ 2 * T1[s] * xij[0] * xij[1]
		+ 2 * T2[s] * xij[0] * xij[2]
		+	  T3[s] * xij[1] * xij[1]
		+ 2 * T4[s] * xij[1] * xij[2]
		+     T5[s] * xij[2] * xij[2]) * norm;

	E2 = (
		+     T0[nj] * xij[0] * xij[0]
		+ 2 * T1[nj] * xij[0] * xij[1]
		+ 2 * T2[nj] * xij[0] * xij[2]
		+     T3[nj] * xij[1] * xij[1]
		+ 2 * T4[nj] * xij[1] * xij[2]
		+     T5[nj] * xij[2] * xij[2]) * norm;

	return 0.5*((E1 - Emin[s]) / (Emax[s] - E1 + 0.0000001) + (E2 - Emin[nj]) / (Emax[nj] - E2 + 0.0000001));
}

float Eint(int nj, int s, int nk)
{
	float xij[3];
	xij[0] = x[nj] - x[s];
	xij[1] = y[nj] - y[s];
	xij[2] = z[nj] - z[s];

	float xik[3];
	xik[0] = x[nk] - x[s];
	xik[1] = y[nk] - y[s];
	xik[2] = z[nk] - z[s];

	float norm_ij = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
	float norm_ik = xik[0] * xik[0] + xik[1] * xik[1] + xik[2] * xik[2];
	float norm = sqrtf(norm_ij*norm_ik);
	float dot = xij[0] * xik[0] + xij[1] * xik[1] + xij[2] * xik[2];

	float cos = dot / norm;

	return (1 + cos) / (1.0001 - cos);
}

float Ei_C(int s)
{
	float E = 0;

	E += (1 - sig[s])*fabsf(cc[s] - 2) / (nc[s] - 2 + 0.0000001);
	E += sig[s] * fabsf(cc[s] - 1) / (nc[s] - 1 + 0.0000001);

	return E;
}

float Ei_D(int s)
{
	int count = 0;
	float E = 0;

	for (int j = n_id[s]; j < n_id[s] + nc[s]; j++)
	{
		E += c[j] * Edata(s, n[j]);
		count += c[j];
	}

	return E / (count + 0.0000001);
}

float Ei_I(int s)
{
	int count = 0;
	float E = 0;

	for (int j = n_id[s]; j < n_id[s] + nc[s]; j++)
	for (int k = j + 1; k < n_id[s] + nc[s]; k++)
	{
		E += c[j] * c[k] * Eint(n[j], s, n[k]);
		count += c[j] * c[k];
	}

	return E / (count + 0.0000001);
}