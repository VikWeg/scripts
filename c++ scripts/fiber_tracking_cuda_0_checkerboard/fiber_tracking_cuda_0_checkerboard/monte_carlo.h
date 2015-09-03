void mc_c(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* vox_ids, int spin_id, int vox_id, int vox_num)
{
	float E0, E1, p;

	int c0 = cube_size[0];
	int c1 = cube_size[1];
	
	int NonEmptyNeighborVoxsCount = GetNonEmptyNeighborVoxsCount(vox_num,vox_ids);

	int NeighborSpinsCovered = 0;
	for (int NeighborNumber = 0; NeighborNumber < NonEmptyNeighborVoxsCount; NeighborNumber++)
	{
		int NeighborVoxNum = GetVoxNumFromNeighborNum(vox_num,NeighborNumber,vox_ids);
		int NeighborVoxId = vox_ids[NeighborVoxNum];

		int next = NeighborVoxNum + 1;
		while (vox_ids[next] < 0) next++;
		int	SpinsInNeighborVoxel = vox_ids[next] - NeighborVoxId;

		for (int NeighborSpinNumber = 0; NeighborSpinNumber < SpinsInNeighborVoxel; NeighborSpinNumber++)
			{
			int NeighborSpinId = NeighborVoxId + NeighborSpinNumber;

				E0 =	wc(T)*Ei_c(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num)
					+	wx(T)*Ei_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num);

				//=== Change connection between spin and neighbour ===

				c[spin_id] = changeBit(NeighborSpinsCovered + NeighborSpinNumber, c[spin_id]);
				
				int VoxNeighborNum = GetNeighborNumFromVoxNum(NeighborVoxNum,vox_num,vox_ids);
				int SpinNumber = spin_id - vox_id;
				int ConnectivityPosition = GetConnectivityOffset(NeighborVoxNum, VoxNeighborNum,vox_ids) + SpinNumber;

				c[NeighborSpinId] = changeBit(ConnectivityPosition, c[NeighborSpinId]);

				//OCTREE???

				//=====================================================

				E1 =	wc(T)*Ei_c(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num)
					+	wx(T)*Ei_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num);

				p = fminf(1., expf((E0 - E1) / T));
				std::bernoulli_distribution acceptQ(p);

				if (acceptQ(generate))
					;
				else
				{
					c[spin_id] = changeBit(NeighborSpinsCovered + NeighborSpinNumber, c[spin_id]);
					c[NeighborSpinId] = changeBit(ConnectivityPosition, c[NeighborSpinId]);
				}
			}
			NeighborSpinsCovered += SpinsInNeighborVoxel;
		}
}

void mc_x(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* vox_ids, int spin_id, int vox_id, int vox_num)
{

	float x0, y0, z0, x1, y1, z1;
	float E0, E1, p;

		for (int i = 0; i < nx; i++)
		{
			E0 = wx(T)*Ei_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num);

			x0 = x[vox_id];
			y0 = y[vox_id];
			z0 = z[vox_id];

			int pos_x = round(x0 / hdr.pixdim[1]);
			int	pos_y = round(y0 / hdr.pixdim[2]);
			int	pos_z = round(z0 / hdr.pixdim[3]);

			std::uniform_real_distribution<float> u_x(x0 - (x0 - pos_x + hdr.pixdim[1] * 0.5) * delta_x, x0 + (pos_x + hdr.pixdim[1] * 0.5 - x0) * delta_x);
			x[vox_id] = u_x(generate);
			std::uniform_real_distribution<float> u_y(y0 - (y0 - pos_y + hdr.pixdim[2] * 0.5) * delta_x, y0 + (pos_y + hdr.pixdim[2] * 0.5 - y0) * delta_x);
			y[vox_id] = u_y(generate);
			std::uniform_real_distribution<float> u_z(z0 - (z0 - pos_z + hdr.pixdim[3] * 0.5) * delta_x, z0 + (pos_z + hdr.pixdim[3] * 0.5 - z0) * delta_x);
			z[vox_id] = u_z(generate);

			E1 = wx(T)*Ei_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num);

			p = fminf(1., expf((E0 - E1) / T));
			std::bernoulli_distribution acceptQ(p);

			if (acceptQ(generate))
				;
			else
			{
				x[vox_id] = x0;
				y[vox_id] = y0;
				z[vox_id] = z0;
			}
		}
}

void mc(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* vox_ids, int lattice_id)
{
	for (int threadId = 0; threadId <= cube_size[0] * cube_size[1] * cube_size[2]; threadId++)
	{
		int vox_num = threadIdToVoxNum(threadId, lattice_id);

		if (vox_id[vox_num] >= 0)
		{
			int next = vox_num + 1;
			while (vox_id[next] < 0) next++;

			for (char spin_id = vox_id[vox_num]; spin_id < vox_id[next]; spin_id++)
			{
				mc_c(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id[vox_num], vox_num);
				mc_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id[vox_num], vox_num);
			}
		}
	}
}