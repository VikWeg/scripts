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
	cudaMemset(&dev_pixdim_0, hdr.pixdim[1], sizeof(float));
	cudaMemset(&dev_pixdim_1, hdr.pixdim[2], sizeof(float));
	cudaMemset(&dev_pixdim_2, hdr.pixdim[3], sizeof(float));

	cudaMemset(&dev_nx, nx, sizeof(int));
	cudaMemset(&dev_delta_x, delta_x, sizeof(float));

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
	
	scount = 0;
	snum = new int[cube_size[0] * cube_size[1] * cube_size[2]];
	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	{
		coo = (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
			+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
			+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

		if ((wmask[offset(i, j, k)] == 1 || surf_mask[offset(i, j, k)] == 1))
		{
			//******** scount + snum ********
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
	}

	ensemble = new vertex[scount];
	dev_in_ensemble = new vertex[scount];
	dev_out_ensemble = new vertex[scount];

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

					for (int s = 0; s < snum[offset(i, j, k)]; s++)
					{
						//******** T ********
						ensemble[vertex_offset(i, j, k) + s].T0 = (ten[4] * ten[4] - ten[3] * ten[5])* norm;
						ensemble[vertex_offset(i, j, k) + s].T1 = (ten[1] * ten[5] - ten[2] * ten[4])* norm;
						ensemble[vertex_offset(i, j, k) + s].T2 = (ten[2] * ten[3] - ten[1] * ten[4])* norm;
						ensemble[vertex_offset(i, j, k) + s].T3 = (ten[2] * ten[2] - ten[0] * ten[5])* norm;
						ensemble[vertex_offset(i, j, k) + s].T4 = (ten[0] * ten[4] - ten[1] * ten[2])* norm;
						ensemble[vertex_offset(i, j, k) + s].T5 = (ten[1] * ten[1] - ten[0] * ten[3])* norm;

						dev_in_ensemble[vertex_offset(i, j, k) + s].T0 = (ten[4] * ten[4] - ten[3] * ten[5])* norm;
						dev_in_ensemble[vertex_offset(i, j, k) + s].T1 = (ten[1] * ten[5] - ten[2] * ten[4])* norm;
						dev_in_ensemble[vertex_offset(i, j, k) + s].T2 = (ten[2] * ten[3] - ten[1] * ten[4])* norm;
						dev_in_ensemble[vertex_offset(i, j, k) + s].T3 = (ten[2] * ten[2] - ten[0] * ten[5])* norm;
						dev_in_ensemble[vertex_offset(i, j, k) + s].T4 = (ten[0] * ten[4] - ten[1] * ten[2])* norm;
						dev_in_ensemble[vertex_offset(i, j, k) + s].T5 = (ten[1] * ten[1] - ten[0] * ten[3])* norm;

						dev_out_ensemble[vertex_offset(i, j, k) + s].T0 = (ten[4] * ten[4] - ten[3] * ten[5])* norm;
						dev_out_ensemble[vertex_offset(i, j, k) + s].T1 = (ten[1] * ten[5] - ten[2] * ten[4])* norm;
						dev_out_ensemble[vertex_offset(i, j, k) + s].T2 = (ten[2] * ten[3] - ten[1] * ten[4])* norm;
						dev_out_ensemble[vertex_offset(i, j, k) + s].T3 = (ten[2] * ten[2] - ten[0] * ten[5])* norm;
						dev_out_ensemble[vertex_offset(i, j, k) + s].T4 = (ten[0] * ten[4] - ten[1] * ten[2])* norm;
						dev_out_ensemble[vertex_offset(i, j, k) + s].T5 = (ten[1] * ten[1] - ten[0] * ten[3])* norm;

						//******** Emin & Emax ********
						float norm2 = ten[2] * ten[2] * ten[3] + ten[1] * ten[1] * ten[5] - 2 * ten[1] * ten[2] * ten[4] + ten[0] * (ten[4] * ten[4] - ten[3] * ten[5]);

						ensemble[vertex_offset(i, j, k) + s].Emin = norm*norm2 / L1data[coo];
						ensemble[vertex_offset(i, j, k) + s].Emax = norm*norm2 / L3data[coo];
						ensemble[vertex_offset(i, j, k) + s].delta_E = 1. / (ensemble[vertex_offset(i, j, k) + s].Emax - ensemble[vertex_offset(i, j, k) + s].Emin + 0.000001);

						dev_in_ensemble[vertex_offset(i, j, k) + s].Emin = norm*norm2 / L1data[coo];
						dev_in_ensemble[vertex_offset(i, j, k) + s].Emax = norm*norm2 / L3data[coo];
						dev_in_ensemble[vertex_offset(i, j, k) + s].delta_E = 1. / (ensemble[vertex_offset(i, j, k) + s].Emax - ensemble[vertex_offset(i, j, k) + s].Emin + 0.000001);

						dev_out_ensemble[vertex_offset(i, j, k) + s].Emin = norm*norm2 / L1data[coo];
						dev_out_ensemble[vertex_offset(i, j, k) + s].Emax = norm*norm2 / L3data[coo];
						dev_out_ensemble[vertex_offset(i, j, k) + s].delta_E = 1. / (ensemble[vertex_offset(i, j, k) + s].Emax - ensemble[vertex_offset(i, j, k) + s].Emin + 0.000001);

						//******** pos ********
						ensemble[vertex_offset(i, j, k) + s].x = i*hdr.pixdim[1];
						ensemble[vertex_offset(i, j, k) + s].y = j*hdr.pixdim[2];
						ensemble[vertex_offset(i, j, k) + s].z = k*hdr.pixdim[3];
						ensemble[vertex_offset(i, j, k) + s].pos_x = i*hdr.pixdim[1];
						ensemble[vertex_offset(i, j, k) + s].pos_y = j*hdr.pixdim[2];
						ensemble[vertex_offset(i, j, k) + s].pos_z = k*hdr.pixdim[3];

						dev_in_ensemble[vertex_offset(i, j, k) + s].x = i*hdr.pixdim[1];
						dev_in_ensemble[vertex_offset(i, j, k) + s].y = j*hdr.pixdim[2];
						dev_in_ensemble[vertex_offset(i, j, k) + s].z = k*hdr.pixdim[3];
						dev_in_ensemble[vertex_offset(i, j, k) + s].pos_x = i*hdr.pixdim[1];
						dev_in_ensemble[vertex_offset(i, j, k) + s].pos_y = j*hdr.pixdim[2];
						dev_in_ensemble[vertex_offset(i, j, k) + s].pos_z = k*hdr.pixdim[3];

						dev_out_ensemble[vertex_offset(i, j, k) + s].x = i*hdr.pixdim[1];
						dev_out_ensemble[vertex_offset(i, j, k) + s].y = j*hdr.pixdim[2];
						dev_out_ensemble[vertex_offset(i, j, k) + s].z = k*hdr.pixdim[3];
						dev_out_ensemble[vertex_offset(i, j, k) + s].pos_x = i*hdr.pixdim[1];
						dev_out_ensemble[vertex_offset(i, j, k) + s].pos_y = j*hdr.pixdim[2];
						dev_out_ensemble[vertex_offset(i, j, k) + s].pos_z = k*hdr.pixdim[3];

						//******** sig ********
						if (surf_mask[offset(i, j, k)] == 1 || i == 0 || i == cube_size[0] - 1 || j == 0 || j == cube_size[1] - 1 || k == 0 || k == cube_size[2] - 1)
						{
							ensemble[vertex_offset(i, j, k) + s].sig = 1;
							dev_in_ensemble[vertex_offset(i, j, k) + s].sig = 1;
							dev_out_ensemble[vertex_offset(i, j, k) + s].sig = 1;
						}
						else
						{
							ensemble[vertex_offset(i, j, k) + s].sig = 0;
							dev_in_ensemble[vertex_offset(i, j, k) + s].sig = 0;
							dev_out_ensemble[vertex_offset(i, j, k) + s].sig = 0;
						}
					}
			}

	free(data);
	free(L1data);
	free(L2data);
	free(L3data);

	int nn;
	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
		if (wmask[offset(i, j, k)] == 1 || surf_mask[offset(i, j, k)] == 1)
		for (int s = 0; s < snum[offset(i, j, k)]; s++)
		{
			//******** nn ********
			nn = 0;
			for (int ii = fmax(0, i - 1); ii < fmin(cube_size[0], i + 2); ii++)
			for (int jj = fmax(0, j - 1); jj < fmin(cube_size[1], j + 2); jj++)
			for (int kk = fmax(0, k - 1); kk < fmin(cube_size[2], k + 2); kk++)
			if ((ii != i || jj != j || kk != k) && (wmask[offset(ii, jj, kk)] == 1 || surf_mask[offset(ii, jj, kk)] == 1))
				nn += snum[offset(ii, jj, kk)];

			ensemble[vertex_offset(i, j, k) + s].nn = nn;
			dev_in_ensemble[vertex_offset(i, j, k) + s].nn = nn;
			dev_out_ensemble[vertex_offset(i, j, k) + s].nn = nn;

			ensemble[vertex_offset(i, j, k) + s].n = new vertex* [nn];

			int m = 0;
			for (int ii = fmax(0, i - 1); ii < fmin(cube_size[0], i + 2); ii++)
			for (int jj = fmax(0, j - 1); jj < fmin(cube_size[1], j + 2); jj++)
			for (int kk = fmax(0, k - 1); kk < fmin(cube_size[2], k + 2); kk++)
			if ((ii != i || jj != j || kk != k) && (wmask[offset(ii, jj, kk)] == 1 || surf_mask[offset(ii, jj, kk)] == 1))
				for (int ss = 0; ss < snum[offset(ii, jj, kk)]; ss++)
				{
					ensemble[vertex_offset(i, j, k) + s].n[m] = &ensemble[vertex_offset(ii, jj, kk) + ss];
					m++;
				}

			//******** c ********
			ensemble[vertex_offset(i, j, k) + s].cc = 0;
			ensemble[vertex_offset(i, j, k) + s].c = new int[nn];

			dev_in_ensemble[vertex_offset(i, j, k) + s].c = new int[nn];
			dev_in_ensemble[vertex_offset(i, j, k) + s].cc = 0;

			dev_out_ensemble[vertex_offset(i, j, k) + s].c = new int[nn];
			dev_out_ensemble[vertex_offset(i, j, k) + s].cc = 0;

			for (int n = 0; n < nn; n++)
			{
				ensemble[vertex_offset(i, j, k) + s].c[n] = 0;
				dev_in_ensemble[vertex_offset(i, j, k) + s].c[n] = 0;
				dev_out_ensemble[vertex_offset(i, j, k) + s].c[n] = 0;
			}
		}

		cudaMalloc((void**)&dev_in_ensemble, scount*sizeof(vertex));
		cudaMalloc((void**)&dev_out_ensemble, scount*sizeof(vertex));

		cudaMemcpy(dev_in_ensemble, ensemble, scount*sizeof(vertex), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_out_ensemble, ensemble, scount*sizeof(vertex), cudaMemcpyHostToDevice);

		dev_in_c = new int*[scount];
		dev_out_c = new int*[scount];

		//dev_in_n = new vertex**[scount];
		//dev_out_n = new vertex**[scount];

		for (int i = 0; i < cube_size[0]; i++)
		for (int j = 0; j < cube_size[1]; j++)
		for (int k = 0; k < cube_size[2]; k++)
		if (wmask[offset(i, j, k)] == 1 || surf_mask[offset(i, j, k)] == 1)
		for (int s = 0; s < snum[offset(i, j, k)]; s++)
		{
			int id = vertex_offset(i, j, k) + s;
			int nn = ensemble[id].nn;
			//int* c = ensemble[id].c;

			cudaMemset(&(dev_in_ensemble[id].c), 0, nn*sizeof(int));

			/*
			dev_in_c[id] = new int[nn];
			dev_out_c[id] = new int[nn];

			cudaMalloc(&(dev_in_c[id]), nn * sizeof(int));
			cudaMalloc(&(dev_out_c[id]), nn * sizeof(int));

			cudaMemcpy(dev_in_c[id], c, nn*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_out_c[id], c, nn*sizeof(int), cudaMemcpyHostToDevice);

			cudaMemcpy(&(dev_in_ensemble[id].c), &(dev_in_c[id]), sizeof(int*), cudaMemcpyHostToDevice);
			cudaMemcpy(&(dev_out_ensemble[id].c), &(dev_out_c[id]), sizeof(int*), cudaMemcpyHostToDevice);

			//cudaMemcpy(&(dev_in_ensemble->c), c, nn*sizeof(int), cudaMemcpyHostToDevice);
			//cudaMemcpy(&(dev_out_ensemble->c), c, nn*sizeof(int), cudaMemcpyHostToDevice);

			
			cudaMemcpy(dev_in_c[id], c, nn*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_out_c[id], c, nn*sizeof(int), cudaMemcpyHostToDevice);

			//=======================================================================================//

			vertex** n = ensemble[id].n;

			dev_in_n[id] = new vertex*[nn];
			dev_out_n[id] = new vertex*[nn];

			cudaMalloc(&dev_in_n[id], nn * sizeof(vertex*));
			cudaMalloc(&dev_out_n[id], nn * sizeof(vertex*));

			cudaMemcpy(dev_in_n[id], n, nn*sizeof(vertex*), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_out_n[id], n, nn*sizeof(vertex*), cudaMemcpyHostToDevice);

			cudaMemcpy(&dev_in_ensemble->n, &dev_in_n[id], sizeof(vertex*), cudaMemcpyHostToDevice);
			cudaMemcpy(&dev_out_ensemble->n, &dev_out_n[id], sizeof(vertex*), cudaMemcpyHostToDevice);
			*/
		}
}