void mc_c(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int VoxNum, int VoxId, int SpinId)
{
	float E0, E1, p;
	
	int nc = GetNeighborSpinCount(VoxNum, VoxIds);

	for (int i = 0; i < nc; i++)
	{
		Next Neighbor = GetNextSpin(VoxIds, VoxNum, i);

		E0 = wc(T)*Ei_c(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Neighbor.VoxNum, Neighbor.SpinId)
			+ wx(T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId);

		c[SpinId] = changeBit(i, c[SpinId]);

		int VoxNeighborNum = GetNeighborNumFromVoxNum(Neighbor.VoxNum, VoxNum, VoxIds);
		int SpinNumber = SpinId - VoxId;
		int ConnectivityPosition = GetConnectivityOffset(Neighbor.VoxNum, VoxNeighborNum, VoxIds) + SpinNumber;

		c[Neighbor.SpinId] = changeBit(ConnectivityPosition, c[Neighbor.SpinId]);

		E1 = wc(T)*Ei_c(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Neighbor.VoxNum, Neighbor.SpinId)
			+ wx(T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId);

		p = fminf(1., expf((E0 - E1) / T));
		std::bernoulli_distribution acceptQ(p);

		if (acceptQ(generate))
			;
		else
		{
			c[SpinId] = changeBit(i, c[SpinId]);
			c[Neighbor.SpinId] = changeBit(ConnectivityPosition, c[Neighbor.SpinId]);
		}
	}
}

void mc_x(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int VoxNum, int SpinId)
{

	float x0, y0, z0;
	float E0, E1, p;

	float pos_x = round(x[SpinId] / hdr.pixdim[1]) * hdr.pixdim[1];
	float pos_y = round(y[SpinId] / hdr.pixdim[2]) * hdr.pixdim[2];
	float pos_z = round(z[SpinId] / hdr.pixdim[3]) * hdr.pixdim[3];

		for (int i = 0; i < nx; i++)
		{
			E0 = wx(T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId);

			x0 = x[SpinId];
			y0 = y[SpinId];
			z0 = z[SpinId];

			std::uniform_real_distribution<float> u_x(x0 - (x0 - pos_x + hdr.pixdim[1] * 0.5f) * delta_x, x0 + (pos_x + hdr.pixdim[1] * 0.5f - x0) * delta_x);
			x[SpinId] = u_x(generate);
			std::uniform_real_distribution<float> u_y(y0 - (y0 - pos_y + hdr.pixdim[2] * 0.5f) * delta_x, y0 + (pos_y + hdr.pixdim[2] * 0.5f - y0) * delta_x);
			y[SpinId] = u_y(generate);
			std::uniform_real_distribution<float> u_z(z0 - (z0 - pos_z + hdr.pixdim[3] * 0.5f) * delta_x, z0 + (pos_z + hdr.pixdim[3] * 0.5f - z0) * delta_x);
			z[SpinId] = u_z(generate);

			E1 = wx(T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId);

			p = fminf(1., expf((E0 - E1) / T));
			std::bernoulli_distribution acceptQ(p);

			if (acceptQ(generate))
				;
			else
			{
				x[SpinId] = x0;
				y[SpinId] = y0;
				z[SpinId] = z0;
			}
		}
}

void mc(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int LatticeId)
{
	//for (int i = 0; i < cube_size[0]; i++)
	//for (int j = 0; j < cube_size[1]; j++)
	//for (int k = 0; k < cube_size[2]; k++)

	//int MaxThreadId = getMaxThreadId(LatticeId);

	//for (int threadId = 0; threadId < MaxThreadId; threadId++)
	//{
		//int VoxNum = threadIdToVoxNum(threadId, LatticeId);
		//int VoxNum = threadId;
		//int VoxNum = u_k(generate)*cube_size[0] * cube_size[1] + u_j(generate)*cube_size[0] + u_i(generate);
		//int VoxNum = k*cube_size[0] * cube_size[1] + j*cube_size[0] + i;
	for (int VoxNum = 0; VoxNum < cube_size[0] * cube_size[1] * cube_size[2]; VoxNum++)
	if (VoxIds[VoxNum] >= 0)
	{
		int SpinsInVoxel;

		if (VoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
		{
			int next = VoxNum + 1;
			while (VoxIds[next] < 0) next++;
			SpinsInVoxel = VoxIds[next] - VoxIds[VoxNum];
		}
		else
			SpinsInVoxel = scount - VoxIds[VoxNum]; //scount is global

		for (int SpinId = VoxIds[VoxNum]; SpinId < VoxIds[VoxNum] + SpinsInVoxel; SpinId++)
		{
			mc_c(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxIds[VoxNum], SpinId);
			mc_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId);
		}
	}
}