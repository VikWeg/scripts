void init_test_ensemble()
{
	ensemble = new vertex***;

	float ten[6] = { 0, 0, 0, 0, 0, 0 };

	for (int i = 0; i < 10; i++)
	{
		ensemble[i] = new vertex**;
		snum[i] = new int*;

		for (int j = 0; j < 10; j++)
		{
			ensemble[i][j] = new vertex*;
			snum[i][j] = new int;

			for (int k = 0; k < 1; k++)
			{
				//******** Tensor ********
				if (j<3 && i<7 && i>2 || j>6 && i<7 && i>2)
				{
					ten[0] = 1;
					ten[1] = 0;
					ten[2] = 0;
					ten[3] = 5;
					ten[4] = 0;
					ten[5] = 1;

					snum[i][j][k] = 1;
				}
				else if (j>2 && j < 7 && i<3 || j>2 && j<7 && i>6)
				{
					ten[0] = 5;
					ten[1] = 0;
					ten[2] = 0;
					ten[3] = 1;
					ten[4] = 0;
					ten[5] = 1;

					snum[i][j][k] = 1;
				}
				else if (i>2 && i<7 && j>2 && j < 7)
				{
					ten[0] = 1.5;
					ten[1] = 0;
					ten[2] = 0;
					ten[3] = 1.5;
					ten[4] = 0;
					ten[5] = 1;

					snum[i][j][k] = 2;
				}
				else
					snum[i][j][k] = 1;

				float norm = 1. / (ten[1] * ten[1] + ten[2] * ten[2] + ten[4] * ten[4] - ten[0] * ten[3] - ten[5] * (ten[0] + ten[3]));

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

					if (i<3 && j<7 && j>2 || i>6 && j<7 && j>2)
					{
						ensemble[i][j][k][s].Emin = 1. / 5 * norm * norm2;
						ensemble[i][j][k][s].Emax = 1. / 1 * norm * norm2;
					}
					else if (i > 2 && i < 7 && j<3 || i>2 && i<7 && j>6)
					{
						ensemble[i][j][k][s].Emin = 1. / 5 * norm * norm2;
						ensemble[i][j][k][s].Emax = 1. / 1 * norm * norm2;
					}
					else if (i>2 && i<7 && j>2 && j < 7)
					{
						ensemble[i][j][k][s].Emin = 1. / 1.5 * norm * norm2;
						ensemble[i][j][k][s].Emax = 1. / 1.0 * norm * norm2;
					}

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
					if (   (i == 0 && j<7 && j>2)
						|| (i == 9 && j<7 && j>2)
						|| (j == 0 && i<7 && i>2)
						|| (j == 9 && i<7 && i>2))
						ensemble[i][j][k][s].sig = 1;
					else if( (i<3 && (j<3 || j>6))
						  || (i>6 && (j<3 || j>6)) )
						ensemble[i][j][k][s].sig = -1;
					else
						ensemble[i][j][k][s].sig = 0;
				}
			}
		}
	}

	int nn;
	for (int i = 0; i < 10; i++)
	for (int j = 0; j < 10; j++)
	for (int k = 0; k < 1; k++)
	{
		for (int s = 0; s < snum[i][j][k]; s++)
		{
			//******** Neighbours ********
			nn = 0;
			for (int ii = fmax(0, i - 1); ii < fmin(10, i + 2); ii++)
			for (int jj = fmax(0, j - 1); jj < fmin(10, j + 2); jj++)
			for (int kk = 0; kk < 1; kk++)
			for (int ss = 0; ss < snum[ii][jj][kk]; ss++)
			if ((ii != i || jj != j) && ensemble[ii][jj][kk][ss].sig >= 0)
				nn++;

			ensemble[i][j][k][s].nn = nn;
			ensemble[i][j][k][s].n = new vertex*[nn];

			int m = 0;
			for (int ii = fmax(0, i - 1); ii < fmin(10, i + 2); ii++)
			for (int jj = fmax(0, j - 1); jj < fmin(10, j + 2); jj++)
			for (int kk = 0; kk < 1; kk++)
			for (int ss = 0; ss < snum[ii][jj][kk]; ss++)
			if ((ii != i || jj != j) && ensemble[ii][jj][kk][ss].sig >= 0)
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