struct EigVals
{
	float L1;
	float L3;
};

EigVals calcEigVals(int VoxNum, float* ten)
{
	EigVals EigVals;

	float a11 = ten[8 * VoxNum];
	float a12 = ten[8 * VoxNum+1];
	float a13 = ten[8 * VoxNum+2];
	float a22 = ten[8 * VoxNum+3];
	float a23 = ten[8 * VoxNum+4];
	float a33 = ten[8 * VoxNum+5];

	float p1 = a12*a12 + a13*a13 + a23*a23;
	float q = (a11 + a22 + a33)*0.33333333f;
	float p2 = (a11 - q)*(a11 - q) + (a22 - q)*(a22 - q) + (a33 - q)*(a33 - q) + 2 * p1;
	float p = sqrtf(p2*0.16666666f);

	float r = 0.5*(1 / (p*p*p)) * (2 * a12*a13*a23 + (a11-q)*( (a22 - q)*(a33 - q) - a23*a23 ) + a13*a13*(q - a22) + a12*a12*(q - a33) );

	//Check for r missing

	float phi = acosf(r)*0.33333333f;

	EigVals.L1 = q + 2 * p*cosf(phi);
	EigVals.L3 = q + 2 * p*cosf(phi + 4.71238898f);

	return EigVals;
}


__host__ __device__ float Edata(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds,
	int VoxNum, int SpinId,
	int NeighborVoxNum, int NeighborSpinId, Constant Const)
{
	float E1,E2;
	float xij [3];
	
	xij[0] = x[SpinId] - x[NeighborSpinId];
	xij[1] = y[SpinId] - y[NeighborSpinId];
	xij[2] = z[SpinId] - z[NeighborSpinId];

	float norm = 1/(xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2]);

	E1 = (		ten[8 * VoxNum    ] * xij[0] * xij[0]
		  + 2 * ten[8 * VoxNum + 1] * xij[0] * xij[1]
		  + 2 * ten[8 * VoxNum + 2] * xij[0] * xij[2]
		  +		ten[8 * VoxNum + 3] * xij[1] * xij[1]
		  + 2 * ten[8 * VoxNum + 4] * xij[1] * xij[2]
		  +		ten[8 * VoxNum + 5] * xij[2] * xij[2]) * norm;

	E2 = (		ten[8 * NeighborVoxNum    ] * xij[0] * xij[0]
		  + 2 * ten[8 * NeighborVoxNum + 1] * xij[0] * xij[1]
		  + 2 * ten[8 * NeighborVoxNum + 2] * xij[0] * xij[2]
		  +		ten[8 * NeighborVoxNum + 3] * xij[1] * xij[1]
		  + 2 * ten[8 * NeighborVoxNum + 4] * xij[1] * xij[2]
		  +		ten[8 * NeighborVoxNum + 5] * xij[2] * xij[2]) * norm;

	//Calculate Emin, Emax each time from tensor (Problem with eigenvalues though...)
	//see 3x3 matrices in https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
	// http://arxiv.org/pdf/physics/0610206v3.pdf

	float E1min = ten[8 * VoxNum + 6];
	float E1max = ten[8 * VoxNum + 7];
	float E2min = ten[8 * NeighborVoxNum + 6];
	float E2max = ten[8 * NeighborVoxNum + 7];

	return	((E1 > E1min) ? (E1 - E1min) : 0) / (((E1 < E1max) ? (E1max - E1) : 0) + Const.eps)
		+ ((E2 > E2min) ? (E2 - E2min) : 0) / (((E2 < E2max) ? (E2max - E2) : 0) + Const.eps);
}

__host__ __device__ float Eint(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds,
			int VoxNum, int SpinId,
			int NeighborVoxNum, int NeighborSpinId,
			int OtherNeighborVoxNum, int OtherNeighborSpinId,
			Constant Const)
{
	float xij[3];
	xij[0] = x[NeighborSpinId] - x[SpinId];
	xij[1] = y[NeighborSpinId] - y[SpinId];
	xij[2] = z[NeighborSpinId] - z[SpinId];

	float xik[3];
	xik[0] = x[OtherNeighborSpinId] - x[SpinId];
	xik[1] = y[OtherNeighborSpinId] - y[SpinId];
	xik[2] = z[OtherNeighborSpinId] - z[SpinId];

	float norm_ij = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
	float norm_ik = xik[0] * xik[0] + xik[1] * xik[1] + xik[2] * xik[2];
	float norm = sqrtf(norm_ij*norm_ik);
	float dot = xij[0] * xik[0] + xij[1] * xik[1] + xij[2] * xik[2];

	float cos = dot / norm;

	return wint(cos,Const);
}

__device__ float Ei_x(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int VoxNum, int SpinId, Constant Const)
{
	float E = 0;
	int count = 0;

	//int nc = GetNeighborSpinCount(VoxNum, VoxIds, Const);

	for (int i = 0; i < 64; i++)
	if (getBit(i, c[SpinId]))
	{
		Next Neighbor = GetNextSpin(VoxIds, VoxNum, i, Const);

		//Add Edata
		E += Edata(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Neighbor.VoxNum, Neighbor.SpinId, Const);
		count++;

		//Add Eint_1
		for (int j = i + 1; j < 64; j++)
		if (getBit(j, c[SpinId]))
		{
			Next OtherNeighbor = GetNextSpin(VoxIds, VoxNum, j, Const);

			E += Eint(x, y, z, c, ten, sig, VoxIds,
				VoxNum, SpinId,
				Neighbor.VoxNum, Neighbor.SpinId,
				OtherNeighbor.VoxNum, OtherNeighbor.SpinId, Const);
			count++;
		}

		//Add Eint_2
		//int nnc = GetNeighborSpinCount(Neighbor.VoxNum, VoxIds, Const);
		for (int j = 0; j < 64; j++)
		if (getBit(j, c[Neighbor.SpinId]))
		{
			Next NextNeighbor = GetNextSpin(VoxIds, Neighbor.VoxNum, j, Const);

			if (NextNeighbor.SpinId != SpinId)
			{
				E += Eint(x, y, z, c, ten, sig, VoxIds,
					Neighbor.VoxNum, Neighbor.SpinId,
					VoxNum, SpinId,
					NextNeighbor.VoxNum, NextNeighbor.SpinId, Const);
				count++;
			}
		}
	}

	return E / (count + Const.eps);
}

__device__ float Ei_c(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int VoxNum, int SpinId, int NeighborVoxNum, int NeighborSpinId, Constant Const)
{
	float E = 0;

	int ConnectionCount = sumBits(c[SpinId]);
	int NearSpinsCount = GetNeighborSpinCount(VoxNum, VoxIds,Const);

	if (sig[VoxNum])
		E += fabsf(ConnectionCount - 1) / (NearSpinsCount - 1 + Const.eps);
	else
		E += fabsf(ConnectionCount - 2) / (NearSpinsCount - 2 + Const.eps);

	int NeighborConnectionCount = sumBits(c[NeighborSpinId]);
	int NeighbourNearSpinsCount = GetNeighborSpinCount(NeighborVoxNum, VoxIds, Const);

	if (sig[NeighborVoxNum])
		E += fabsf(NeighborConnectionCount - 1) / (NeighbourNearSpinsCount - 1 + Const.eps);
	else
		E += fabsf(NeighborConnectionCount - 2) / (NeighbourNearSpinsCount - 2 + Const.eps);

	return E / NearSpinsCount;
}


__host__ float Ei_C(int VoxNum, Constant Const)
{
	int VoxId = VoxIds[VoxNum];

	if (VoxId < 0) return 0;

	float E = 0;
	int	SpinsInVoxel;

	if (VoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
	{
		int next = VoxNum + 1;
		while (VoxIds[next] < 0) next++;
		SpinsInVoxel = VoxIds[next] - VoxId;
	}
	else SpinsInVoxel = scount - VoxId;

	int NearSpinsCount = GetNeighborSpinCount(VoxNum, VoxIds,Const);

	for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
	{
		int ConnectionCount = sumBits(c[SpinId]);

		E += (1 - sig[VoxNum])	*fabsf(ConnectionCount - 2.0f)/(NearSpinsCount - 2.0f + eps);
		E += sig[VoxNum] * fabsf(ConnectionCount - 1.0f)/(NearSpinsCount - 1.0f + eps);
	}

	return E;
}

__host__ float Ei_D(int VoxNum, Constant Const)
{
	int VoxId = VoxIds[VoxNum];
	int	SpinsInVoxel;

	if (VoxId < 0) return 0;

	if (VoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
	{
		int next = VoxNum + 1;
		while (VoxIds[next] < 0) next++;
		SpinsInVoxel = VoxIds[next] - VoxId;
	}
	else SpinsInVoxel = scount - VoxId;

	float E = 0;

	int nc = GetNeighborSpinCount(VoxNum, VoxIds, Const);

	for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
	{
		float Ei = 0;
		int count = 0;

		for (int i = 0; i < nc; i++)
		if (getBit(i, c[SpinId]))
		{
			Next Neighbor = GetNextSpin(VoxIds, VoxNum, i, Const);
			Ei += Edata(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Neighbor.VoxNum, Neighbor.SpinId, Const);
			count++;
		}

		E += Ei / (count + eps);
	}

	return E;
}

__host__ float Ei_I(int VoxNum, Constant Const)
{
	int VoxId = VoxIds[VoxNum];
	int	SpinsInVoxel;

	if (VoxId < 0) return 0;

	if (VoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
	{
		int next = VoxNum + 1;
		while (VoxIds[next] < 0) next++;
		SpinsInVoxel = VoxIds[next] - VoxId;
	}
	else SpinsInVoxel = scount - VoxId;

	float E = 0;

	int nc = GetNeighborSpinCount(VoxNum, VoxIds, Const);

	//Loop over spins in VoxNum

	for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
	{
		float Ei = 0;
		int count = 0;

		for (int i = 0; i < nc; i++)
		{
			if (getBit(i, c[SpinId]))
			{
				Next Neighbor = GetNextSpin(VoxIds, VoxNum, i, Const);

				for (int j = i + 1; j < nc; j++)
				if (getBit(j, c[SpinId]))
				{
					Next OtherNeighbor = GetNextSpin(VoxIds, VoxNum, j, Const);

					Ei += Eint(x, y, z, c, ten, sig, VoxIds,
						VoxNum, SpinId,
						Neighbor.VoxNum, Neighbor.SpinId,
						OtherNeighbor.VoxNum, OtherNeighbor.SpinId, Const);
					count++;
				}

			}
		}

		E += Ei / (count + eps);
	}

	return E;
}
