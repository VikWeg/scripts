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

	cudaMalloc(&dev_pixdim, 3 * sizeof(float));
	cudaMemcpy(dev_pixdim, &hdr.pixdim[1], 3 * sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&dev_nx, sizeof(int));
	cudaMemset(dev_nx, nx, sizeof(int));

	cudaMalloc(&dev_delta_x, sizeof(float));
	cudaMemset(dev_delta_x, delta_x, sizeof(float));

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
	
	//******** box_num ********
	box_num = surf_vox_num + w_vox_num;

	//******** snum, scount ********
	scount = 0;
	snum = new int[box_num];
	int m = 0;
	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	if (wmask[offset(i, j, k)] == 1 || surf_mask[offset(i, j, k)] == 1)
	{
		coo = (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
			+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
			+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

		if (L1data[coo] / L3data[coo] < cutoff && surf_mask[offset(i, j, k)] == 0)
			snum[m] = 3;
		else if ((L2data[coo] - L3data[coo]) / (L1data[coo] - L3data[coo]) > 0.5 && surf_mask[offset(i, j, k)] == 0)
			snum[m] = 2;
		else
			snum[m] = 1;

		m++;
	}

	//======== spin variables ========

	id = new int [box_num];

	T0 = new float[box_num];
	T1 = new float[box_num];
	T2 = new float[box_num];
	T3 = new float[box_num];
	T4 = new float[box_num];
	T5 = new float[box_num];

	Emin = new float[box_num];
	Emax = new float[box_num];
	delta_E = new float[box_num];

	pos_x = new float[box_num];
	pos_y = new float[box_num];
	pos_z = new float[box_num];

	sig = new int[box_num];

	x = new float[scount];
	y = new float[scount];
	z = new float[scount];

	nc = new int[box_num];

	m = 0;
	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
		if ((wmask[offset(i, j, k)] == 1 || surf_mask[offset(i, j, k)] == 1))
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

			//******** id ********
			id[m] = vertex_offset(i, j, k);

			//******** T ********
			T0[m] = (ten[4] * ten[4] - ten[3] * ten[5])* norm;
			T1[m] = (ten[1] * ten[5] - ten[2] * ten[4])* norm;
			T2[m] = (ten[2] * ten[3] - ten[1] * ten[4])* norm;
			T3[m] = (ten[2] * ten[2] - ten[0] * ten[5])* norm;
			T4[m] = (ten[0] * ten[4] - ten[1] * ten[2])* norm;
			T5[m] = (ten[1] * ten[1] - ten[0] * ten[3])* norm;

			//******** Emin & Emax ********
			float norm2 = ten[2] * ten[2] * ten[3] + ten[1] * ten[1] * ten[5] - 2 * ten[1] * ten[2] * ten[4] + ten[0] * (ten[4] * ten[4] - ten[3] * ten[5]);

			Emin[m] = norm*norm2 / L1data[coo];
			Emax[m] = norm*norm2 / L3data[coo];
			delta_E[m] = 1. / (Emax[m] - Emin[m] + 0.000001);

			//******** pos ********
			pos_x[m] = i*hdr.pixdim[1];
			pos_y[m] = j*hdr.pixdim[2];
			pos_z[m] = k*hdr.pixdim[3];

			//******** sig ********
			if (surf_mask[offset(i, j, k)] == 1 || i == 0 || i == cube_size[0] - 1 || j == 0 || j == cube_size[1] - 1 || k == 0 || k == cube_size[2] - 1)
				sig[m] = 1;
			else
				sig[m] = 0;

			//******** x, y, z ********
			for (int s = 0; s < snum[m]; s++)
			{
				x[id[m] + s] = i*hdr.pixdim[1];
				y[id[m] + s] = j*hdr.pixdim[2];
				z[id[m] + s] = k*hdr.pixdim[3];
			}

			//******** nc ********
			nc[m] = 0;
			for (int mm = 0; mm < box_num; mm++)
			if (
				   abs(pos_x[mm] - pos_x[m]) <= 1
				&& abs(pos_y[mm] - pos_y[m]) <= 1
				&& abs(pos_z[mm] - pos_z[m]) <= 1
				&& (pos_x[mm] != pos_x[m] || pos_y[mm] != pos_y[m] || pos_z[mm] != pos_z[m])
				)
				nc[m] += snum[mm];

			m++;
		}

	free(data);
	free(L1data);
	free(L2data);
	free(L3data);

	//******** n_id ********
	n_id = new int[box_num];
	for (int i = 0; i < box_num; i++)
	{
		n_id[i] = 0;
		for (int j = 0; j < i; j++)
			n_id[i] += nc[j];
	}

	//******** n ********
	n = new int[ n_id[box_num-1] + nc[box_num-1] ];
	for (int m = 0; m < box_num; m++)
		for (int mm = 0; mm < box_num; mm++)
		if (
			      abs(pos_x[mm] - pos_x[m]) <= 1
			   && abs(pos_y[mm] - pos_y[m]) <= 1
			   && abs(pos_z[mm] - pos_z[m]) <= 1
			   && (pos_x[mm] != pos_x[m] || pos_y[mm] != pos_y[m] || pos_z[mm] != pos_z[m])
			)
			for (int s = 0; s < snum[mm]; s++)
				n[n_id[m] + s] = id[mm] + s;

	//******** c ********
	for (int m = 0; m < box_num; m++)
		for (int n = 0; n < nc[m]; n++)
			c[n_id[m] + n] = 0;
}