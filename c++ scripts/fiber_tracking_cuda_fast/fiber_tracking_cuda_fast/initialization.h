int offset(int i, int j, int k)
{
	return cube_size[2] * cube_size[1] * i + cube_size[2] * j + k;
}

int vertex_offset(int i, int j, int k)
{
	int count = 0;

	for (int ii = 0; ii < i; ii++)
	for (int jj = 0; jj < cube_size[1]; jj++)
	for (int kk = 0; kk < cube_size[2]; kk++)
		count += snum[offset(ii, jj, kk)];

	for (int ii = i, int jj = 0; jj < j; jj++)
	for (int kk = 0; kk < cube_size[2]; kk++)
		count += snum[offset(ii, jj, kk)];

	for (int ii = i, int jj = j, int kk = 0; kk < k; kk++)
		count += snum[offset(ii, jj, kk)];

	return count;
}

void initialize()
{	
	dim[1] = hdr.dim[1];
	dim[2] = hdr.dim[2];
	dim[3] = hdr.dim[3];

	host_params = new parameters;

	host_params->delta_x = delta_x;
	host_params->nx = nx;
	host_params->pix_dim_x = hdr.pixdim[1];
	host_params->pix_dim_y = hdr.pixdim[2];
	host_params->pix_dim_z = hdr.pixdim[3];

	long coo;
	float ten[6] = { 0, 0, 0, 0, 0, 0 };

	w_vox_num = 0;
	wmask = new int[cube_size[0] * cube_size[1] * cube_size[2]];

	for (int i = 0; i < cube_size[0]; i++)
		for (int j = 0; j < cube_size[1]; j++)
			for (int k = 0; k < cube_size[2]; k++)
			{
				coo = (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
					+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
					+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

				if (mask[coo]>0)
				{
					wmask[offset(i,j,k)] = 1;
					w_vox_num++;
				}
				else
					wmask[offset(i,j,k)] = 0;
			}

	surf_vox_num = 0;
	surf_mask = new int[cube_size[0] * cube_size[1] * cube_size[2]];

	for (int i = 0; i < cube_size[0]; i++)
		for (int j = 0; j < cube_size[1]; j++)
			for (int k = 0; k < cube_size[2]; k++)
			{
				int wn = 0;

				for (int ii = fmax(0, i - 1); ii <= fmin(cube_size[0] - 1, i + 1); ii++)
				for (int jj = fmax(0, j - 1); jj <= fmin(cube_size[1] - 1, j + 1); jj++)
				for (int kk = fmax(0, k - 1); kk <= fmin(cube_size[2] - 1, k + 1); kk++)
				if ((ii != i || jj != j || kk != k) && wmask[offset(ii,jj,kk)] == 1)
					wn++;

				if (wmask[offset(i,j,k)] == 0 && wn > 0)
				{
					surf_mask[offset(i,j,k)] = 1;
					surf_vox_num++;
				}
				else
					surf_mask[offset(i,j,k)] = 0;
			}
	
	//******** snum, scount ********
	scount = 0;
	snum = new int[cube_size[0] * cube_size[1] * cube_size[2]];
	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
		if ((wmask[offset(i, j, k)] == 1 || surf_mask[offset(i, j, k)] == 1))
		{
			coo = (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
				+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
				+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

			if (L1data[coo] / L3data[coo] < cutoff && surf_mask[offset(i, j, k)] == 0)
			{
				scount += 3;
				snum[offset(i, j, k)] = 3;
			}
			else if ((L2data[coo] - L3data[coo]) / (L1data[coo] - L3data[coo]) > 0.5 && surf_mask[offset(i, j, k)] == 0)
			{
				scount += 2;
				snum[offset(i, j, k)] = 2;
			}
			else
			{
				scount += 1;
				snum[offset(i, j, k)] = 1;
			}
		}
		else
			snum[offset(i, j, k)] = 0;

	host_params->scount = scount;

	//======== spin variables ========

	T0 = new float[scount];
	T1 = new float[scount];
	T2 = new float[scount];
	T3 = new float[scount];
	T4 = new float[scount];
	T5 = new float[scount];

	Emin = new float[scount];
	Emax = new float[scount];
	delta_E = new float[scount];

	pos_x = new float[scount];
	pos_y = new float[scount];
	pos_z = new float[scount];

	sig = new int[scount];

	x = new float[scount];
	y = new float[scount];
	z = new float[scount];

	cc = new int[scount];

	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	for (int s = 0; s < snum[offset(i, j, k)]; s++)
	{
		coo = (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
			+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
			+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

		for (int t = 0; t < 6; t++)
		{
			long coo_t = t*dim[1] * dim[2] * dim[3]
				+ (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
				+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
				+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

			ten[t] = data[coo_t];
		}
		float norm = 1. / (ten[1] * ten[1] + ten[2] * ten[2] + ten[4] * ten[4] - ten[0] * ten[3] - ten[5] * (ten[0] + ten[3]));

		int id = vertex_offset(i, j, k) + s;

		//******** T ********
		T0[id] = (ten[4] * ten[4] - ten[3] * ten[5])* norm;
		T1[id] = (ten[1] * ten[5] - ten[2] * ten[4])* norm;
		T2[id] = (ten[2] * ten[3] - ten[1] * ten[4])* norm;
		T3[id] = (ten[2] * ten[2] - ten[0] * ten[5])* norm;
		T4[id] = (ten[0] * ten[4] - ten[1] * ten[2])* norm;
		T5[id] = (ten[1] * ten[1] - ten[0] * ten[3])* norm;

		//******** Emin & Emax ********
		float norm2 = ten[2] * ten[2] * ten[3] + ten[1] * ten[1] * ten[5] - 2 * ten[1] * ten[2] * ten[4] + ten[0] * (ten[4] * ten[4] - ten[3] * ten[5]);

		Emin[id] = norm*norm2 / L1data[coo];
		Emax[id] = norm*norm2 / L3data[coo];
		delta_E[id] = 1. / (Emax[id] - Emin[id] + 0.000001);

		//******** pos ********
		pos_x[id] = i*hdr.pixdim[1];
		pos_y[id] = j*hdr.pixdim[2];
		pos_z[id] = k*hdr.pixdim[3];

		//******** sig ********
		if (surf_mask[offset(i, j, k)] == 1 || i == 0 || i == cube_size[0] - 1 || j == 0 || j == cube_size[1] - 1 || k == 0 || k == cube_size[2] - 1)
			sig[id] = 1;
		else
			sig[id] = 0;

		//******** x, y, z ********
		x[id] = i*hdr.pixdim[1];
		y[id] = j*hdr.pixdim[2];
		z[id] = k*hdr.pixdim[3];
	}

	//******** nc ********
	nc = new int[scount];

	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	for (int s = 0; s < snum[offset(i, j, k)]; s++)
	{
		int id = vertex_offset(i, j, k) + s;

		nc[id] = 0;
		for (int ii = fmax(0, i - 1); ii < fmin(cube_size[0], i + 2); ii++)
		for (int jj = fmax(0, j - 1); jj < fmin(cube_size[1], j + 2); jj++)
		for (int kk = fmax(0, k - 1); kk < fmin(cube_size[2], k + 2); kk++)
		if ((ii != i || jj != j || kk != k) && (wmask[offset(ii, jj, kk)] == 1 || surf_mask[offset(ii, jj, kk)] == 1))
			nc[id] += snum[offset(ii, jj, kk)];
	}

	free(data);
	free(L1data);
	free(L2data);
	free(L3data);

	//******** n_id ********
	n_id = new int[scount];
	for (int s = 0; s < scount; s++)
	{
		n_id[s] = 0;
		for (int ss = 0; ss < s; ss++)
			n_id[s] += nc[ss];
	}

	//******** n_num ********
	n_num = n_id[scount - 1] + nc[scount - 1];

	//******** n ********
	n = new int[n_num];

	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	for (int s = 0; s < snum[offset(i, j, k)]; s++)
	{
		int id = vertex_offset(i, j, k) + s;

		int m = 0;
		for (int ii = fmax(0, i - 1); ii < fmin(cube_size[0], i + 2); ii++)
		for (int jj = fmax(0, j - 1); jj < fmin(cube_size[1], j + 2); jj++)
		for (int kk = fmax(0, k - 1); kk < fmin(cube_size[2], k + 2); kk++)
		if (ii != i || jj != j || kk != k)
		for (int ss = 0; ss < snum[offset(ii, jj, kk)]; ss++)
		{
			n[n_id[id] + m] = vertex_offset(ii,jj,kk) + ss;
			m++;
		}
	}

	//******** cc ********
	for (int s = 0; s < scount; s++)
		cc[s] = 0;

	//******** c ********
	c = new int[n_num];
	for (int s = 0; s < n_num; s++)
		c[s] = 0;
}