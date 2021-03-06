int VoxNum(int i, int j, int k)
{
	return cube_size[1] * cube_size[0] * k + cube_size[0] * j + i;
}

int VoxId(int i, int j, int k)
{
	if (snum[VoxNum(i, j, k)] == 0) return -1;

	int count = 0;

	for (int kk = 0; kk < k; kk++)
	for (int jj = 0; jj < cube_size[1]; jj++)
	for (int ii = 0; ii < cube_size[0]; ii++)
		count += snum[VoxNum(ii, jj, kk)];

	for (int jj = 0; jj < j; jj++)
	for (int ii = 0; ii < cube_size[0]; ii++)
		count += snum[VoxNum(ii, jj, k)];

	for (int ii = 0; ii < i; ii++)
		count += snum[VoxNum(ii, j, k)];

	return count;
}

void init_ensemble()
{
	int coo;
	int dim[4] = { 0, hdr.dim[1], hdr.dim[2], hdr.dim[3] };
	float TenDat[6] = { 0, 0, 0, 0, 0, 0 };
	
	//******** WMask ********
	w_vox_num = 0;
	wmask = new int**[cube_size[0]];
	for (int i = 0; i < cube_size[0]; i++)
	{
		wmask[i] = new int*[cube_size[1]];

		for (int j = 0; j < cube_size[1]; j++)
		{
			wmask[i][j] = new int[cube_size[2]];

			for (int k = 0; k < cube_size[2]; k++)
			{
				coo = (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
					+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
					+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

				if (mask[coo]>0)
				{
					wmask[i][j][k] = 1;
					w_vox_num++;
				}
				else
					wmask[i][j][k] = 0;
			}
		}
	}

	//******** SurfMask ********
	surf_vox_num = 0;
	surf_mask = new int**[cube_size[0]];
	int wn;
	for (int i = 0; i < cube_size[0]; i++)
	{
		surf_mask[i] = new int*[cube_size[1]];

		for (int j = 0; j < cube_size[1]; j++)
		{
			surf_mask[i][j] = new int[cube_size[2]];

			for (int k = 0; k < cube_size[2]; k++)
			{
				wn = 0;

				for (int ii = fmax(0, i - 1); ii <= fmin(cube_size[0] - 1, i + 1); ii++)
				for (int jj = fmax(0, j - 1); jj <= fmin(cube_size[1] - 1, j + 1); jj++)
				for (int kk = fmax(0, k - 1); kk <= fmin(cube_size[2] - 1, k + 1); kk++)
				if ((ii != i || jj != j || kk != k) && wmask[ii][jj][kk] == 1)
					wn++;

				if (wmask[i][j][k] == 0 && wn > 0)
				{
					surf_mask[i][j][k] = 1;
					surf_vox_num++;
				}
				else
					surf_mask[i][j][k] = 0;
			}
		}
	}

	//******** Snum + Scount ********
	snum = new int[cube_size[0] * cube_size[1] * cube_size[2]];
	scount = 0;
	for (int i = 0; i < cube_size[0]; i++)
		for (int j = 0; j < cube_size[1]; j++)
			for (int k = 0; k < cube_size[2]; k++)
				if (wmask[i][j][k] == 1 || surf_mask[i][j][k] == 1)
				{
					coo = (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
						+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
						+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

					if (L1data[coo] / L3data[coo] < cutoff && surf_mask[i][j][k] == 0)
					{
						snum[VoxNum(i,j,k)] = 3; scount += 3;
					}
					else if ((L2data[coo] - L3data[coo]) / (L1data[coo] - L3data[coo]) > 0.5 && surf_mask[i][j][k] == 0)
					{
						snum[VoxNum(i, j, k)] = 2; scount += 2;
					}
					else
					{
						snum[VoxNum(i, j, k)] = 1; scount += 1;
					}
				}
				else
					snum[VoxNum(i, j, k)] = 0;

	//******** x, y, z, VoxIds, c, sig, ten ********

	VoxIds = new int[cube_size[0] * cube_size[1] * cube_size[2]];
	ten = new float[8 * cube_size[0] * cube_size[1] * cube_size[2]];
	sig = new int[cube_size[0] * cube_size[1] * cube_size[2]];
	x = new float[scount];
	y = new float[scount];
	z = new float[scount];
	c = new unsigned long long[scount];

	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	{
		coo = (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
			+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
			+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

		int voxnum = VoxNum(i, j, k);
		int voxid = VoxId(i, j, k);

		if (wmask[i][j][k] == 0 && surf_mask[i][j][k] == 0) {VoxIds[voxnum] = -1; sig[voxnum] = 0; continue; }
		else												 VoxIds[voxnum] = voxid;

		for (int t = 0; t < 6; t++)
		{
			long coo_t = t*dim[1] * dim[2] * dim[3]
				+ (dim[3] - (vox_origin[3] + k - cube_size[2] / 2)) * dim[2] * dim[1]
				+ (dim[2] - (vox_origin[2] + j - cube_size[1] / 2)) * dim[1]
				+ (dim[1] - (vox_origin[1] + i - cube_size[0] / 2));

			TenDat[t] = data[coo_t];
		}
		
		float Li3 = 1. / L1data[coo];
		float Li2 = 1. / L2data[coo];
		float Li1 = 1. / L3data[coo];

		float norm = (Li1 + Li2 + Li3)*(  TenDat[2] * TenDat[2] * TenDat[3]
									- 2 * TenDat[1] * TenDat[2] * TenDat[4]
										+ TenDat[0] * TenDat[4] * TenDat[4]
										+ TenDat[1] * TenDat[1] * TenDat[5]
										- TenDat[0] * TenDat[3] * TenDat[5]);

		ten[8 * voxnum + 0] = (TenDat[4] * TenDat[4] - TenDat[3] * TenDat[5]) / norm;
		ten[8 * voxnum + 1] = (TenDat[1] * TenDat[5] - TenDat[2] * TenDat[4]) / norm;
		ten[8 * voxnum + 2] = (TenDat[2] * TenDat[3] - TenDat[1] * TenDat[4]) / norm;
		ten[8 * voxnum + 3] = (TenDat[2] * TenDat[2] - TenDat[0] * TenDat[5]) / norm;
		ten[8 * voxnum + 4] = (TenDat[0] * TenDat[4] - TenDat[1] * TenDat[2]) / norm;
		ten[8 * voxnum + 5] = (TenDat[1] * TenDat[1] - TenDat[0] * TenDat[3]) / norm;

		ten[8 * voxnum + 6] = Li3 / (Li1 + Li2 + Li3);
		ten[8 * voxnum + 7] = Li1 / (Li1 + Li2 + Li3);

		if (surf_mask[i][j][k] == 1 || i == 0 || i == cube_size[0] - 1 || j == 0 || j == cube_size[1] - 1 || k == 0 || k == cube_size[2] - 1)
			sig[voxnum] = 1;
		else
			sig[voxnum] = 0;

		for (int s = 0; s < snum[voxnum]; s++)
		{
			x[voxid + s] = i*hdr.pixdim[1];
			y[voxid + s] = j*hdr.pixdim[2];
			z[voxid + s] = k*hdr.pixdim[3];

			c[voxid + s] = (unsigned long long)0;
		}
	}

	free(data);
	free(L1data);
	free(L2data);
	free(L3data);
}