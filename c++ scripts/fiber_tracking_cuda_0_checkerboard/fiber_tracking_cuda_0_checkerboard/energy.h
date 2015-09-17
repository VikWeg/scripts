float Edata(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds,
	int VoxNum, int SpinId,
	int NeighborVoxNum, int NeighborSpinId)
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

	float E1min = ten[8 * VoxNum + 6];
	float E1max = ten[8 * VoxNum + 7];
	float E2min = ten[8 * NeighborVoxNum + 6];
	float E2max = ten[8 * NeighborVoxNum + 7];

	return	((E1 > E1min) ? (E1 - E1min) : 0) / (((E1 < E1max) ? (E1max - E1) : 0) + eps)
		+	((E2 > E2min) ? (E2 - E2min) : 0) / (((E2 < E2max) ? (E2max - E2) : 0) + eps);
}

float Eint(	float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds,
			int VoxNum, int SpinId,
			int NeighborVoxNum, int NeighborSpinId,
			int OtherNeighborVoxNum, int OtherNeighborSpinId)
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

	float cos = dot / (norm+eps);

	return wint(cos);
}

float Ei_x(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int VoxNum, int SpinId)
{
	float E = 0;
	int count = 0;

	int nc = GetNeighborSpinCount(VoxNum, VoxIds);

	for (int i = 0; i < nc; i++)
	if (getBit(i, c[SpinId]))
	{
		Next Neighbor = GetNextSpin(VoxIds, VoxNum, i);

		//Add Edata
		E += Edata(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Neighbor.VoxNum, Neighbor.SpinId);
		count++;

		//Add Eint_1
		for (int j = i + 1; j < nc; j++)
		if (getBit(j, c[SpinId]))
		{
			Next OtherNeighbor = GetNextSpin(VoxIds, VoxNum, j);

			E += Eint(x, y, z, c, ten, sig, VoxIds,
				VoxNum, SpinId,
				Neighbor.VoxNum, Neighbor.SpinId,
				OtherNeighbor.VoxNum, OtherNeighbor.SpinId);
			count++;
		}

		//Add Eint_2
		int nnc = GetNeighborSpinCount(Neighbor.VoxNum, VoxIds);
		for (int j = 0; j < nnc; j++)
		if (getBit(j, c[Neighbor.SpinId]))
		{
			Next NextNeighbor = GetNextSpin(VoxIds, Neighbor.VoxNum, j);

			if (NextNeighbor.SpinId != SpinId)
			{
				E += Eint(x, y, z, c, ten, sig, VoxIds,
					Neighbor.VoxNum, Neighbor.SpinId,
					VoxNum, SpinId,
					NextNeighbor.VoxNum, NextNeighbor.SpinId);
				count++;
			}
		}
	}

	return E /(count + eps);
}

float Ei_c(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds,int VoxNum, int SpinId, int NeighborVoxNum, int NeighborSpinId)
{
	float E = 0;

	int ConnectionCount = sumBits(c[SpinId]);
	int NearSpinsCount = GetNeighborSpinCount(VoxNum, VoxIds);

	if (sig[VoxNum])
	E += fabsf(ConnectionCount - 1) / (NearSpinsCount - 1 + eps);
	else
	E += fabsf(ConnectionCount - 2) / (NearSpinsCount - 2 + eps);

	int NeighborConnectionCount = sumBits(c[NeighborSpinId]);
	int NeighbourNearSpinsCount = GetNeighborSpinCount(NeighborVoxNum, VoxIds);

	if (sig[NeighborVoxNum])
	E += fabsf(NeighborConnectionCount - 1) / (NeighbourNearSpinsCount - 1 + eps);
	else
	E += fabsf(NeighborConnectionCount - 2) / (NeighbourNearSpinsCount - 2 + eps);

	return E / NearSpinsCount;
}


float Ei_C(int VoxNum)
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

	int NearSpinsCount = GetNeighborSpinCount(VoxNum, VoxIds);

	for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
	{
		int ConnectionCount = sumBits(c[SpinId]);

		E += (1 - sig[VoxNum])	*fabsf(ConnectionCount - 2.0f)/(NearSpinsCount - 2.0f + eps);
		E += sig[VoxNum] * fabsf(ConnectionCount - 1.0f)/(NearSpinsCount - 1.0f + eps);
	}

	return E;
}

float Ei_D(int VoxNum)
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

	int nc = GetNeighborSpinCount(VoxNum, VoxIds);

	for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
	{
		float Ei = 0;
		int count = 0;

		for (int i = 0; i < nc; i++)
		if (getBit(i, c[SpinId]))
		{
			Next Neighbor = GetNextSpin(VoxIds, VoxNum, i);
			Ei += Edata(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Neighbor.VoxNum, Neighbor.SpinId);
			count++;
		}

		E += Ei / (count + eps);
	}

	return E;
}

float Ei_I(int VoxNum)
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

	int nc = GetNeighborSpinCount(VoxNum, VoxIds);

	//Loop over spins in VoxNum

	for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
	{
		float Ei = 0;
		int count = 0;

		for (int i = 0; i < nc; i++)
		{
			if (getBit(i, c[SpinId]))
			{
				Next Neighbor = GetNextSpin(VoxIds, VoxNum, i);

				for (int j = i + 1; j < nc; j++)
				if (getBit(j, c[SpinId]))
				{
					Next OtherNeighbor = GetNextSpin(VoxIds, VoxNum, j);

					Ei += Eint(x, y, z, c, ten, sig, VoxIds,
						VoxNum, SpinId,
						Neighbor.VoxNum, Neighbor.SpinId,
						OtherNeighbor.VoxNum, OtherNeighbor.SpinId);
					count++;
				}

			}
		}

		E += Ei / (count + eps);
	}

	return E;
}
