
float Edata(vertex* vi, vertex* vj)
{
	float E1,E2;
	float xij [3];
	
	for (int n = 0; n < 3; n++)
		xij[n] = (*vi).x[n] - (*vj).x[n];

	float norm = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];

	E1 = (		(*vi).T[0] * xij[0] * xij[0]
		+ 2 *	(*vi).T[1] * xij[0] * xij[1]
		+ 2 *	(*vi).T[2] * xij[0] * xij[2]
		+		(*vi).T[3] * xij[1] * xij[1]
		+ 2 *	(*vi).T[4] * xij[1] * xij[2]
		+		(*vi).T[5] * xij[2] * xij[2]) / norm;

	E2 = (
		+		(*vj).T[0] * xij[0] * xij[0]
		+ 2 *	(*vj).T[1] * xij[0] * xij[1]
		+ 2 *	(*vj).T[2] * xij[0] * xij[2]
		+		(*vj).T[3] * xij[1] * xij[1]
		+ 2 *	(*vj).T[4] * xij[1] * xij[2]
		+		(*vj).T[5] * xij[2] * xij[2]) / norm;

	return (E1 - vi->Emin) / (vi->Emax - vi->Emin) + (E2 - vj->Emin) / (vj->Emax - vj->Emin);
}

float Eint(vertex* vj, vertex* vi, vertex* vk)
{
	float xij[3];
	for (int n = 0; n < 3; n++)
		xij[n] = (*vj).x[n] - (*vi).x[n];

	float xik[3];
	for (int n = 0; n < 3; n++)
		xik[n] = (*vk).x[n] - (*vi).x[n];

	float norm_ij = sqrtf(xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2]);
	float norm_ik = sqrtf(xik[0] * xik[0] + xik[1] * xik[1] + xik[2] * xik[2]);

	//return  logf((cosmax + 1)/(cosmax - (xij[0] * xik[0] + xij[1] * xik[1] + xij[2] * xik[2]) / norm_ij / norm_ik));
	return ((xij[0] * xik[0] + xij[1] * xik[1] + xij[2] * xik[2]) / norm_ij / norm_ik + 1)/2;
}

float Ei_x(vertex* vi)
{
	int count = 0;
	float E = 0;

	for (int j = 0; j < vi->nn; j++)
	{
		E += vi->c[j] * Edata(vi, vi->n[j]);
		count += vi->c[j];
	}

	for (int j = 0; j < vi->nn; j++)
		for (int k = j + 1; k < vi->nn; k++)
		{
			E += vi->c[j] * vi->c[k] * Eint(vi->n[j], vi, vi->n[k]);
			count += vi->c[j] * vi->c[k];
		}

	for (int j = 0; j < vi->nn; j++)
		for (int k = 0; k < vi->n[j]->nn; k++)
			if (vi->n[j]->n[k] != vi)
			{
				E += vi->c[j] * vi->n[j]->c[k] * Eint(vi, vi->n[j], vi->n[j]->n[k]);
				count += vi->c[j] * vi->n[j]->c[k];
			}

	return E/(count+0.00001);
}

float Ei_c(vertex* vi)
{

	float E = 0;

	E += (1 - vi->sig)*fabsf(vi->cc - 2);
	E += (vi->sig)*    fabsf(vi->cc - 1);

	E /= vi->nn;

	int cj;
	for (int j = 0; j < vi->nn; j++)
	{
		cj = vi->n[j]->cc;

		E += (1 - vi->n[j]->sig)*fabsf(cj - 2) / (float)(vi->n[j]->nn);
		E += (vi->n[j]->sig)*    fabsf(cj - 1) / (float)(vi->n[j]->nn);
	}

	return E/(vi->nn + 1);
}

float Ei_C(vertex* vi)
{
	float E = 0;

	E += (1 - vi->sig)*fabsf(vi->cc - 2);
	E += (vi->sig)*    fabsf(vi->cc - 1);

	return E / vi->nn;
}

float Ei_D(vertex* vi)
{
	int count = 0;
	float E = 0;

	for (int j = 0; j < vi->nn; j++)
	{
		E += vi->c[j] * Edata(vi, vi->n[j]);
		count += vi->c[j];
	}

	return 0.5 *E /(count + 0.00001);
}

float Ei_I(vertex* vi)
{
	int count=0;
	float E = 0;

	for (int j = 0; j < vi->nn; j++)
		for (int k = j + 1; k < vi->nn; k++)
		{
			E += vi->c[j] * vi->c[k] * Eint(vi->n[j], vi, vi->n[k]);
			count += vi->c[j] * vi->c[k];
		}

	return E / (count+0.00001);
}