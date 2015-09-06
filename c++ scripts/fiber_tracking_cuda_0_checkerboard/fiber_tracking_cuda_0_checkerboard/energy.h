float Edata(float* x, float* y, float* z, unsigned long long* c, float* ten, char* sig, int* VoxIds,
	int VoxNum, int VoxId, int SpinId,
	int NeighborVoxNum, int NeighborVoxId, int NeighborSpinId)
{
	float E1,E2;
	float xij [3];
	
	xij[0] = x[SpinId] - x[NeighborSpinId];
	xij[1] = y[SpinId] - y[NeighborSpinId];
	xij[2] = z[SpinId] - z[NeighborSpinId];

	float norm = 1/(xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2]);

	E1 = (		ten[VoxNum] * xij[0] * xij[0]
		+ 2 * ten[VoxNum+1] * xij[0] * xij[1]
		+ 2 * ten[VoxNum+2] * xij[0] * xij[2]
		+	  ten[VoxNum+3] * xij[1] * xij[1]
		+ 2 * ten[VoxNum+4] * xij[1] * xij[2]
		+	  ten[VoxNum+5] * xij[2] * xij[2]) * norm;

	E2 = (	  ten[NeighborVoxNum]     * xij[0] * xij[0]
		+ 2 * ten[NeighborVoxNum + 1] * xij[0] * xij[1]
		+ 2 * ten[NeighborVoxNum + 2] * xij[0] * xij[2]
		+     ten[NeighborVoxNum + 3] * xij[1] * xij[1]
		+ 2 * ten[NeighborVoxNum + 4] * xij[1] * xij[2]
		+     ten[NeighborVoxNum + 5] * xij[2] * xij[2]) * norm;

	//Calculate Emin, Emax each time from tensor (Problem with eigenvalues though...)
	//see 3x3 matrices in https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices

	return 0.5*(E1 + E2);
}

float Eint(	float* x, float* y, float* z, unsigned long long* c, float* ten, char* sig, int* VoxIds,
			int VoxNum, int VoxId, int SpinId,
			int NeighborVoxNum, int NeighborVoxId, int NeighborSpinId,
			int OtherNeighborVoxNum, int OtherNeighborVoxId, int OtherNeighborSpinId)
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

	return wint(cos);
}

float Ei_x(float* x, float* y, float* z, unsigned long long* c, float* ten, char* sig, int* VoxIds, int VoxNum, int VoxId, int SpinId)
{
	float E = 0;

	int NonEmptyNeighborVoxsCount = GetNonEmptyNeighborVoxsCount(VoxNum, VoxIds);

	int NeighborSpinsCovered = 0;

	//Loop over neighbor voxels
	for (int NeighborNumber = 0; NeighborNumber < NonEmptyNeighborVoxsCount; NeighborNumber++)
	{
		int NeighborVoxNum = GetVoxNumFromNeighborNum(VoxNum, NeighborNumber, VoxIds);
		int NeighborVoxId = VoxIds[NeighborVoxNum];
		int	SpinsInNeighborVoxel;

		if (NeighborVoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
		{
			int next = NeighborVoxNum + 1;
			while (VoxIds[next] < 0) next++;
			SpinsInNeighborVoxel = VoxIds[next] - NeighborVoxId;
		}
		else SpinsInNeighborVoxel = scount - NeighborVoxId; //global scount

		//Loop over spins in current neighbor voxel
		for (int NeighborSpinNumber = 0; NeighborSpinNumber < SpinsInNeighborVoxel; NeighborSpinNumber++)
		{
			int NeighborSpinId = NeighborVoxId + NeighborSpinNumber;

			//Add data and first+second order interaction energy
			if (getBit(NeighborSpinsCovered + NeighborSpinNumber, c[SpinId]))
			{
				E += Edata(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxId, SpinId, NeighborVoxNum, NeighborVoxId, NeighborSpinId);

				//Add first order smoothness energy
				int OtherNeighborSpinsCovered = 0;
				for (int OtherNeighborNumber = NeighborNumber + 1; OtherNeighborNumber < NonEmptyNeighborVoxsCount; OtherNeighborNumber++)
				{
					int OtherNeighborVoxNum = GetVoxNumFromNeighborNum(VoxNum, OtherNeighborNumber, VoxIds);
					int OtherNeighborVoxId = VoxIds[OtherNeighborVoxNum];

					int	SpinsInOtherNeighborVoxel;
					if (OtherNeighborVoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
					{
						int next = OtherNeighborVoxNum + 1;
						while (VoxIds[next] < 0) next++;
						SpinsInOtherNeighborVoxel = VoxIds[next] - OtherNeighborVoxId;
					}
					else SpinsInOtherNeighborVoxel = scount - OtherNeighborVoxId; //global scount

					//Loop over spins in current other neighbor voxel
					for (int OtherNeighborSpinNumber = 0; OtherNeighborSpinNumber < SpinsInOtherNeighborVoxel; OtherNeighborSpinNumber++)
					{
						int OtherNeighborSpinId = OtherNeighborVoxId + OtherNeighborSpinNumber;

						if (getBit(OtherNeighborSpinsCovered + OtherNeighborSpinNumber, c[SpinId]))
							E += Eint(	x, y, z, c, ten, sig, VoxIds,
										VoxNum, VoxId, SpinId,
										NeighborVoxNum, NeighborVoxId, NeighborSpinId,
										OtherNeighborVoxNum, OtherNeighborVoxId, OtherNeighborSpinId);
					}
					OtherNeighborSpinsCovered += SpinsInOtherNeighborVoxel;
				}

				//Add second order smoothness energy
				int NextNeighborSpinsCovered = 0;
				int VoxNeighborNumber = GetNeighborNumFromVoxNum(NeighborVoxNum, VoxNum, VoxIds);
				int NeighborNonEmptyNeighborVoxsCount = GetNonEmptyNeighborVoxsCount(NeighborVoxNum, VoxIds);

				for (int NextNeighborNumber = 0; NextNeighborNumber < NeighborNonEmptyNeighborVoxsCount; NextNeighborNumber++)
				{
					int NextNeighborVoxNum = GetVoxNumFromNeighborNum(NeighborVoxNum, NextNeighborNumber, VoxIds);
					int NextNeighborVoxId = VoxIds[NextNeighborVoxNum];
					int	SpinsInNextNeighborVoxel;

					if (NextNeighborVoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
					{
						int next = NextNeighborVoxNum + 1;
						while (VoxIds[next] < 0) next++;
						SpinsInNextNeighborVoxel = VoxIds[next] - NextNeighborVoxId;
					}
					else SpinsInNextNeighborVoxel = scount - NextNeighborVoxId;

					//Loop over spins in current next neighbor voxel
					if (NextNeighborNumber != VoxNeighborNumber)
					for (int NextNeighborSpinNumber = 0; NextNeighborSpinNumber < SpinsInNextNeighborVoxel; NextNeighborSpinNumber++)
					{
						int NextNeighborSpinId = NextNeighborVoxId + NextNeighborSpinNumber;

						if (getBit(NextNeighborSpinsCovered + NextNeighborSpinNumber, c[NeighborSpinId]))
							E += Eint(x, y, z, c, ten, sig, VoxIds,
										NeighborVoxNum, NeighborVoxId, NeighborSpinId,
										VoxNum, VoxId, SpinId,
										NextNeighborVoxNum, NextNeighborVoxId, NextNeighborSpinId);
					}
					NextNeighborSpinsCovered += SpinsInNextNeighborVoxel;
				}
			}
		}
		NeighborSpinsCovered += SpinsInNeighborVoxel;
	}

	return E;
}

float Ei_c(float* x, float* y, float* z, unsigned long long* c, float* ten, char* sig, int* VoxIds,int VoxNum, int VoxId, int SpinId, int NeighborVoxNum, int NeighborSpinId)
{
	int c0 = cube_size[0];
	int c1 = cube_size[1];

	float E = 0;

	int ConnectionCount = sumBits(c[SpinId]);
	int NearSpinsCount = GetNeighborSpinCount(VoxNum, VoxIds);

	E += (1 - sig[VoxNum])	*fabsf(ConnectionCount - 2) / (NearSpinsCount - 2 + 0.0000001);
	E += sig[VoxNum] * fabsf(ConnectionCount - 1) / (NearSpinsCount - 1 + 0.0000001);

	int NeighborConnectionCount = sumBits(c[NeighborSpinId]);
	int NeighbourNearSpinsCount = GetNeighborSpinCount(NeighborVoxNum, VoxIds);

	E += (1 - sig[NeighborVoxNum])	*fabsf(NeighborConnectionCount - 2) / (NeighbourNearSpinsCount - 2 + 0.0000001);
	E += sig[NeighborVoxNum] * fabsf(NeighborConnectionCount - 1) / (NeighbourNearSpinsCount - 1 + 0.0000001);

	return E / (NearSpinsCount + 1);
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

		E += (1 - sig[VoxNum])	*fabsf(ConnectionCount - 2) / (NearSpinsCount - 2 + 0.0000001);
		E += sig[VoxNum] * fabsf(ConnectionCount - 1) / (NearSpinsCount - 1 + 0.0000001);
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

	int NonEmptyNeighborVoxsCount = GetNonEmptyNeighborVoxsCount(VoxNum, VoxIds);

	int NeighborSpinsCovered = 0;

	//Loop over neighbor voxels
	for (int NeighborNumber = 0; NeighborNumber < NonEmptyNeighborVoxsCount; NeighborNumber++)
	{
		int NeighborVoxNum = GetVoxNumFromNeighborNum(VoxNum, NeighborNumber, VoxIds);
		int NeighborVoxId = VoxIds[NeighborVoxNum];

		int	SpinsInNeighborVoxel;
		if (NeighborVoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
		{
			int next = NeighborVoxNum + 1;
			while (VoxIds[next] < 0) next++;
			SpinsInNeighborVoxel = VoxIds[next] - NeighborVoxId;
		}
		else SpinsInNeighborVoxel = scount - NeighborVoxId; //global scount

		//Loop over spins in current neighbor voxel
		for (int NeighborSpinNumber = 0; NeighborSpinNumber < SpinsInNeighborVoxel; NeighborSpinNumber++)
		{
			int NeighborSpinId = NeighborVoxId + NeighborSpinNumber;

			//Add data energy for each spin in VoxNum
			for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
			if (getBit(NeighborSpinsCovered + NeighborSpinNumber, c[SpinId]))
				E += Edata(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxId, SpinId, NeighborVoxNum, NeighborVoxId, NeighborSpinId);
			
		}
		NeighborSpinsCovered += SpinsInNeighborVoxel;
	}

	return E;
}

float Ei_I(int VoxNum)
{
	int VoxId = VoxIds[VoxNum];

	if (VoxId < 0) return 0;

	int	SpinsInVoxel;
	if (VoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
	{
		int next = VoxNum + 1;
		while (VoxIds[next] < 0) next++;
		SpinsInVoxel = VoxIds[next] - VoxId;
	}
	else SpinsInVoxel = scount - VoxId;

	float E = 0;

	int NonEmptyNeighborVoxsCount = GetNonEmptyNeighborVoxsCount(VoxNum, VoxIds);

	int NeighborSpinsCovered = 0;

	//Loop over neighbor voxels
	for (int NeighborNumber = 0; NeighborNumber < NonEmptyNeighborVoxsCount; NeighborNumber++)
	{
		int NeighborVoxNum = GetVoxNumFromNeighborNum(VoxNum, NeighborNumber, VoxIds);
		int NeighborVoxId = VoxIds[NeighborVoxNum];

		int	SpinsInNeighborVoxel;
		if (NeighborVoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
		{
			int next = NeighborVoxNum + 1;
			while (VoxIds[next] < 0) next++;
			SpinsInNeighborVoxel = VoxIds[next] - NeighborVoxId;
		}
		else SpinsInNeighborVoxel = scount - NeighborVoxId; //global scount

		//Loop over spins in current neighbor voxel
		for (int NeighborSpinNumber = 0; NeighborSpinNumber < SpinsInNeighborVoxel; NeighborSpinNumber++)
		{
			int NeighborSpinId = NeighborVoxId + NeighborSpinNumber;

			//Add first order interaction energy for each spin in VoxNum
			for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
			if (getBit(NeighborSpinsCovered + NeighborSpinNumber, c[SpinId]))
			{
				int OtherNeighborSpinsCovered = 0;
				for (int OtherNeighborNumber = NeighborNumber + 1; OtherNeighborNumber < NonEmptyNeighborVoxsCount; OtherNeighborNumber++)
				{
					int OtherNeighborVoxNum = GetVoxNumFromNeighborNum(VoxNum, OtherNeighborNumber, VoxIds);
					int OtherNeighborVoxId = VoxIds[OtherNeighborVoxNum];

					int	SpinsInOtherNeighborVoxel;
					if (OtherNeighborVoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
					{
						int next = OtherNeighborVoxNum + 1;
						while (VoxIds[next] < 0) next++;
						SpinsInOtherNeighborVoxel = VoxIds[next] - OtherNeighborVoxId;
					}
					else SpinsInOtherNeighborVoxel = scount - OtherNeighborVoxId; //global scount

					//Loop over spins in current other neighbor voxel
					for (int OtherNeighborSpinNumber = 0; OtherNeighborSpinNumber < SpinsInOtherNeighborVoxel; OtherNeighborSpinNumber++)
					{
						int OtherNeighborSpinId = OtherNeighborVoxId + OtherNeighborSpinNumber;

						if (getBit(OtherNeighborSpinsCovered + OtherNeighborSpinNumber, c[SpinId]))
							E += Eint(x, y, z, c, ten, sig, VoxIds,
							VoxNum, VoxId, SpinId,
							NeighborVoxNum, NeighborVoxId, NeighborSpinId,
							OtherNeighborVoxNum, OtherNeighborVoxId, OtherNeighborSpinId);
					}
					OtherNeighborSpinsCovered += SpinsInOtherNeighborVoxel;
				}
			}
		}
		NeighborSpinsCovered += SpinsInNeighborVoxel;
	}

	return E;
}
