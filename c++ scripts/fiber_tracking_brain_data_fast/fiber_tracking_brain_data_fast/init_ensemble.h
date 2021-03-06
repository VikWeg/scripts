void init_ensemble()
{
	long coo;
	int dim[4] = { 0, hdr.dim[1], hdr.dim[2], hdr.dim[3] };
	float ten[6] = { 0, 0, 0, 0, 0, 0 };
	
	w_vox_num = 0;
	wmask = new int**[size[0]];

	for (int i = 0; i < size[0]; i++)
	{
		wmask[i] = new int*[size[1]];

		for (int j = 0; j < size[1]; j++)
		{
			wmask[i][j] = new int [size[2]];

			for (int k = 0; k < size[2]; k++)
			{
				coo = (dim[3] - (center[3] + k - size[2] / 2)) * dim[2] * dim[1]
					+ (dim[2] - (center[2] + j - size[1] / 2)) * dim[1]
					+ (dim[1] - (center[1] + i - size[0] / 2));

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

	surf_vox_num = 0;
	surf_mask = new int**[size[0]];
	int wn;

	for (int i = 0; i < size[0]; i++)
	{
		surf_mask[i] = new int*[size[1]];

		for (int j = 0; j < size[1]; j++)
		{
			surf_mask[i][j] = new int[size[2]];

			for (int k = 0; k < size[2]; k++)
			{
				wn = 0;

				for (int ii = fmax(0, i - 1); ii <= fmin(size[0]-1, i + 1); ii++)
				for (int jj = fmax(0, j - 1); jj <= fmin(size[1]-1, j + 1); jj++)
				for (int kk = fmax(0, k - 1); kk <= fmin(size[2]-1, k + 1); kk++)
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
	
	ensemble = new vertex*** [size[0]];
	snum = new int**[size[0]];

	for (int i = 0; i < size[0]; i++)
	{
		ensemble[i] = new vertex**[size[1]];
		snum[i] = new int* [size[1]];

		for (int j = 0; j < size[1]; j++)
		{
			ensemble[i][j] = new vertex*[size[2]];
			snum[i][j] = new int[size[2]];

			for (int k = 0; k < size[2]; k++)
			{
				coo = (dim[3] - (center[3] + k - size[2] / 2)) * dim[2] * dim[1]
					+ (dim[2] - (center[2] + j - size[1] / 2)) * dim[1]
					+ (dim[1] - (center[1] + i - size[0] / 2));

				if ((wmask[i][j][k] == 1 || surf_mask[i][j][k] == 1))// && L1data[coo] > 0 && L2data[coo] > 0 && L3data[coo] > 0 && fabs(L1data[coo] / L3data[coo]) < 10)
				{
					//******** Snum ********
					if (L1data[coo] / L3data[coo] < cutoff && surf_mask[i][j][k] == 0)
						snum[i][j][k] = 3;
					else if ((L2data[coo] - L3data[coo]) / (L1data[coo] - L3data[coo]) > 0.5 && surf_mask[i][j][k] == 0)
						snum[i][j][k] = 2;
					else
						snum[i][j][k] = 1;

					ensemble[i][j][k] = new vertex[snum[i][j][k]];

					//******** Tensor ********
					for (int t = 0; t < 6; t++)
					{
						long coo_t = t*dim[1] * dim[2] * dim[3]
							+ (dim[3] - (center[3] + k - size[2] / 2)) * dim[2] * dim[1]
							+ (dim[2] - (center[2] + j - size[1] / 2)) * dim[1]
							+ (dim[1] - (center[1] + i - size[0] / 2));

						ten[t] = data[coo_t];
					}
					float norm = 1. / (ten[1] * ten[1] + ten[2] * ten[2] + ten[4] * ten[4] - ten[0] * ten[3] - ten[5] * (ten[0] + ten[3]));

					for (int s = 0; s < snum[i][j][k]; s++)
					{
						ensemble[i][j][k][s].T0 = (ten[4] * ten[4] - ten[3] * ten[5])* norm;
						ensemble[i][j][k][s].T1 = (ten[1] * ten[5] - ten[2] * ten[4])* norm;
						ensemble[i][j][k][s].T2 = (ten[2] * ten[3] - ten[1] * ten[4])* norm;
						ensemble[i][j][k][s].T3 = (ten[2] * ten[2] - ten[0] * ten[5])* norm;
						ensemble[i][j][k][s].T4 = (ten[0] * ten[4] - ten[1] * ten[2])* norm;
						ensemble[i][j][k][s].T5 = (ten[1] * ten[1] - ten[0] * ten[3])* norm;

						//******** Emin & Emax ********
						float norm2 = ten[2] * ten[2] * ten[3] + ten[1] * ten[1] * ten[5] - 2 * ten[1] * ten[2] * ten[4] + ten[0] * (ten[4] * ten[4] - ten[3] * ten[5]);

						ensemble[i][j][k][s].Emin = norm*norm2 / L1data[coo];
						ensemble[i][j][k][s].Emax = norm*norm2 / L3data[coo];
						ensemble[i][j][k][s].delta_E = 1. / (ensemble[i][j][k][s].Emax - ensemble[i][j][k][s].Emin + 0.000001);

						//******** Position ********
						ensemble[i][j][k][s].x = i;
						ensemble[i][j][k][s].y = j;
						ensemble[i][j][k][s].z = k;

						ensemble[i][j][k][s].pos_x = i;
						ensemble[i][j][k][s].pos_y = j;
						ensemble[i][j][k][s].pos_z = k;

						//******** Surface ********
						if (surf_mask[i][j][k] == 1)
							ensemble[i][j][k][s].sig = 1;
						else
							ensemble[i][j][k][s].sig = 0;
					}
				}
				else
				snum[i][j][k] = 0;
			}
		}
	}
	free(data);
	free(L1data);
	free(L2data);
	free(L3data);

	int nn;
	for (int i = 0; i < size[0]; i++)
		for (int j = 0; j < size[1]; j++)
			for (int k = 0; k < size[2]; k++)
			if (wmask[i][j][k] == 1 || surf_mask[i][j][k] == 1)
			{
				for (int s = 0; s < snum[i][j][k]; s++)
				{
					//******** Neighbours ********
					nn = 0;
					for (int ii = fmax(0, i - 1); ii < fmin(size[0], i + 2); ii++)
					for (int jj = fmax(0, j - 1); jj < fmin(size[1], j + 2); jj++)
					for (int kk = fmax(0, k - 1); kk < fmin(size[2], k + 2); kk++)
					for (int ss = 0; ss < snum[ii][jj][kk]; ss++)
					if ((ii != i || jj != j || kk != k) && (wmask[ii][jj][kk] == 1 || surf_mask[ii][jj][kk] == 1))
						nn++;

					ensemble[i][j][k][s].nn = nn;
					ensemble[i][j][k][s].n = new vertex*[nn];

					int m = 0;
					for (int ii = fmax(0, i - 1); ii < fmin(size[0], i + 2); ii++)
					for (int jj = fmax(0, j - 1); jj < fmin(size[1], j + 2); jj++)
					for (int kk = fmax(0, k - 1); kk < fmin(size[2], k + 2); kk++)
					for (int ss = 0; ss < snum[ii][jj][kk]; ss++)
					if ((ii != i || jj != j || kk != k) && (wmask[ii][jj][kk] == 1 || surf_mask[ii][jj][kk] == 1))
					{
						ensemble[i][j][k][s].n[m] = &ensemble[ii][jj][kk][ss];
						m++;
					}

					//******** Connected ********
					ensemble[i][j][k][s].cc = 0;
					ensemble[i][j][k][s].c = new int[nn];
					for (int n = 0; n < nn; n++)
						ensemble[i][j][k][s].c[n] = 0;
				}
			}
}