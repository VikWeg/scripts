void mc_c(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* VoxIds, int VoxNum, int VoxId, int SpinId)
{
	float E0, E1, p;

	int c0 = cube_size[0];
	int c1 = cube_size[1];
	
	int NonEmptyNeighborVoxsCount = GetNonEmptyNeighborVoxsCount(VoxNum,VoxIds);

	int NeighborSpinsCovered = 0;
	for (int NeighborNumber = 0; NeighborNumber < NonEmptyNeighborVoxsCount; NeighborNumber++)
	{
		int NeighborVoxNum = GetVoxNumFromNeighborNum(VoxNum,NeighborNumber,VoxIds);
		int NeighborVoxId = VoxIds[NeighborVoxNum];

		int next = NeighborVoxNum + 1;
		while (VoxIds[next] < 0) next++;
		int	SpinsInNeighborVoxel = VoxIds[next] - NeighborVoxId;

		for (int NeighborSpinNumber = 0; NeighborSpinNumber < SpinsInNeighborVoxel; NeighborSpinNumber++)
			{
			int NeighborSpinId = NeighborVoxId + NeighborSpinNumber;

			E0 =		wc(T)*Ei_c(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxId, SpinId, NeighborVoxNum, NeighborSpinId)
				+		wx(T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxId, SpinId);

				//=== Change connection between spin and neighbour ===

				c[SpinId] = changeBit(NeighborSpinsCovered + NeighborSpinNumber, c[SpinId]);
				
				int VoxNeighborNum = GetNeighborNumFromVoxNum(NeighborVoxNum,VoxNum,VoxIds);
				int SpinNumber = SpinId - VoxId;
				int ConnectivityPosition = GetConnectivityOffset(NeighborVoxNum, VoxNeighborNum,VoxIds) + SpinNumber;

				c[NeighborSpinId] = changeBit(ConnectivityPosition, c[NeighborSpinId]);

				//=====================================================

				E1 =	wc(T)*Ei_c(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxId, SpinId, NeighborVoxNum, NeighborSpinId)
					+	wx(T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxId, SpinId);

				p = fminf(1., expf((E0 - E1) / T));
				std::bernoulli_distribution acceptQ(p);

				if (acceptQ(generate))
					;
				else
				{
					c[SpinId] = changeBit(NeighborSpinsCovered + NeighborSpinNumber, c[SpinId]);
					c[NeighborSpinId] = changeBit(ConnectivityPosition, c[NeighborSpinId]);
				}
			}
			NeighborSpinsCovered += SpinsInNeighborVoxel;
		}
}

void mc_x(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* VoxIds, int VoxNum, int VoxId, int SpinId)
{

	float x0, y0, z0, x1, y1, z1;
	float E0, E1, p;

		for (int i = 0; i < nx; i++)
		{
			E0 = wx(T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxId, SpinId);

			x0 = x[SpinId];
			y0 = y[SpinId];
			z0 = z[SpinId];

			int pos_x = round(x0 / hdr.pixdim[1]);
			int	pos_y = round(y0 / hdr.pixdim[2]);
			int	pos_z = round(z0 / hdr.pixdim[3]);

			std::uniform_real_distribution<float> u_x(x0 - (x0 - pos_x + hdr.pixdim[1] * 0.5) * delta_x, x0 + (pos_x + hdr.pixdim[1] * 0.5 - x0) * delta_x);
			x[SpinId] = u_x(generate);
			std::uniform_real_distribution<float> u_y(y0 - (y0 - pos_y + hdr.pixdim[2] * 0.5) * delta_x, y0 + (pos_y + hdr.pixdim[2] * 0.5 - y0) * delta_x);
			y[SpinId] = u_y(generate);
			std::uniform_real_distribution<float> u_z(z0 - (z0 - pos_z + hdr.pixdim[3] * 0.5) * delta_x, z0 + (pos_z + hdr.pixdim[3] * 0.5 - z0) * delta_x);
			z[SpinId] = u_z(generate);

			E1 = wx(T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxId, SpinId);

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

void mc(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* VoxIds, int lattice_id)
{
	for (int threadId = 0; threadId <= cube_size[0] * cube_size[1] * cube_size[2]; threadId++)
	{
		int VoxNum = threadIdToVoxNum(threadId, lattice_id);

		if (VoxIds[VoxNum] >= 0)
		{
			int next = VoxNum + 1;
			while (VoxIds[next] < 0) next++;

			for (char SpinId = VoxIds[VoxNum]; SpinId < VoxIds[next]; SpinId++)
			{
				mc_c(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxIds[VoxNum], SpinId);
				mc_x(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxIds[VoxNum], SpinId);
			}
		}
	}
}