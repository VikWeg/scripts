void init_ensemble()
{
	ensemble = new vertex***;
	snum = new int**;

	int coo;
	int center[4] = { 0, 48, 64, 64 };
	int dim[4] = { 0, hdr.dim[1], hdr.dim[2], hdr.dim[3] };
	float ten[6] = { 0, 0, 0, 0, 0, 0 };

	for (int i = 0; i < size; i++)
	{
		ensemble[i] = new vertex**;
		snum[i] = new int*;

		for (int j = 0; j < size; j++)
		{
			ensemble[i][j] = new vertex*;
			snum[i][j] = new int;

			for (int k = 0; k < size; k++)
			{
				//******** Tensor ********
				for (int t = 0; t < 6; t++)
					{
						coo = t*dim[1] * dim[2] * dim[3]
							+ (dim[3] - (center[3] + k - size / 2)) * dim[2] * dim[1]
							+ (dim[2] - (center[2] + j - size / 2)) * dim[1]
							+ (dim[1] - (center[1] + i - size / 2));

						ten[t] = data[coo];
					}

				float norm = 1. / (ten[1] * ten[1] + ten[2] * ten[2] + ten[4] * ten[4] - ten[0] * ten[3] - ten[5] * (ten[0] + ten[3]));

				coo = (dim[3] - (center[3] + k - size / 2)) * dim[2] * dim[1]
					+ (dim[2] - (center[2] + j - size / 2)) * dim[1]
					+ (dim[1] - (center[1] + i - size / 2));

				if (L1data[coo] / L3data[coo] < cutoff)
					snum[i][j][k] = 3;
				else if ((L2data[coo]  - L3data[coo] ) / (L1data[coo]  - L3data[coo] ) > 0.5)
					snum[i][j][k] = 2;
				else
					snum[i][j][k] = 1;

				ensemble[i][j][k] = new vertex[snum[i][j][k]];

				for (int s = 0; s < snum[i][j][k]; s++)
				{
					ensemble[i][j][k][s].T = new float[6];
					ensemble[i][j][k][s].T[0] = (ten[4] * ten[4] - ten[3] * ten[5])* norm;
					ensemble[i][j][k][s].T[1] = (ten[1] * ten[5] - ten[2] * ten[4])* norm;
					ensemble[i][j][k][s].T[2] = (ten[2] * ten[3] - ten[1] * ten[4])* norm;
					ensemble[i][j][k][s].T[3] = (ten[2] * ten[2] - ten[0] * ten[5])* norm;
					ensemble[i][j][k][s].T[4] = (ten[0] * ten[4] - ten[1] * ten[2])* norm;
					ensemble[i][j][k][s].T[5] = (ten[1] * ten[1] - ten[0] * ten[3])* norm;

					//******** Emin & Emax ********
					float norm2 = ten[2] * ten[2] * ten[3] + ten[1] * ten[1] * ten[5] - 2 * ten[1] * ten[2] * ten[4] + ten[0] * (ten[4] * ten[4] - ten[3] * ten[5]);

					ensemble[i][j][k][s].Emin = 1. / L1data[coo] * norm * norm2;
					ensemble[i][j][k][s].Emax = 1. / L3data[coo] * norm * norm2;

					//******** Position ********
					ensemble[i][j][k][s].x = new float[3];
					ensemble[i][j][k][s].x[0] = i;
					ensemble[i][j][k][s].x[1] = j;
					ensemble[i][j][k][s].x[2] = k;

					ensemble[i][j][k][s].pos = new int[3];
					ensemble[i][j][k][s].pos[0] = i;
					ensemble[i][j][k][s].pos[1] = j;
					ensemble[i][j][k][s].pos[2] = k;

					//******** Surface ********
					if (i == 0 || j == 0 || k == 0 || i == size - 1 || j == size - 1 || k == size - 1)
						ensemble[i][j][k][s].sig = 1;
					else
						ensemble[i][j][k][s].sig = 0;
				}
			}
		}
	}
	free(data);
	free(L1data);
	free(L2data);
	free(L3data);

	int nn;
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			for (int k = 0; k < size; k++)
				for (int s = 0; s < snum[i][j][k]; s++)
				{
					//******** Neighbours ********
					nn = 0;
					for (int ii = fmax(0, i - 1); ii < fmin(size, i + 2); ii++)
					for (int jj = fmax(0, j - 1); jj < fmin(size, j + 2); jj++)
					for (int kk = fmax(0, k - 1); kk < fmin(size, k + 2); kk++)
					for (int ss = 0; ss < snum[ii][jj][kk]; ss++)
					if (ii != i || jj != j || kk != k)
						nn++;

					ensemble[i][j][k][s].nn = nn;
					ensemble[i][j][k][s].n = new vertex*[nn];

					int m = 0;
					for (int ii = fmax(0, i - 1); ii < fmin(size, i + 2); ii++)
						for (int jj = fmax(0, j - 1); jj < fmin(size, j + 2); jj++)
							for (int kk = fmax(0, k - 1); kk < fmin(size, k + 2); kk++)
								for (int ss = 0; ss < snum[ii][jj][kk]; ss++)
									if (ii != i || jj != j || kk != k)
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