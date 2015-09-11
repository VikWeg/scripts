__device__ void mc_c(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int VoxNum, int VoxId, int SpinId, Constant Const)
{
	curandState s;
	curand_init(SpinId, 0, 0, &s);

	float E0=0, E1=0, p;
	
	//int nc = GetNeighborSpinCount(VoxNum, VoxIds, Const);

	for (int i = 0; i < 64; i++)
	{
		Next Neighbor = GetNextSpin(VoxIds, VoxNum, i, Const);

		E0 = wc(Const.T)*Ei_c(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Neighbor.VoxNum, Neighbor.SpinId, Const)
			+ wx(Const.T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Const);

		c[SpinId] = changeBit(i, c[SpinId]);

		int VoxNeighborNum = GetNeighborNumFromVoxNum(Neighbor.VoxNum, VoxNum, VoxIds, Const);
		int SpinNumber = SpinId - VoxId;
		int ConnectivityPosition = GetConnectivityOffset(Neighbor.VoxNum, VoxNeighborNum, VoxIds, Const) + SpinNumber;

		c[Neighbor.SpinId] = changeBit(ConnectivityPosition, c[Neighbor.SpinId]);

		E1 = wc(Const.T)*Ei_c(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Neighbor.VoxNum, Neighbor.SpinId, Const)
			+ wx(Const.T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Const);

		p = fminf(1., expf((E0 - E1) / Const.T));

		if (curand_uniform(&s) < p)
			;
		else
		{
			c[SpinId] = changeBit(i, c[SpinId]);
			c[Neighbor.SpinId] = changeBit(ConnectivityPosition, c[Neighbor.SpinId]);
		}
	}
}

__device__ void mc_x(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int VoxNum, int SpinId, Constant Const)
{
	curandState s;
	curand_init(SpinId, 0, 0, &s);

	float x0, y0, z0;
	float E0, E1, p;

	float pos_x = round(x[SpinId] / Const.pixdim0) * Const.pixdim0;
	float pos_y = round(y[SpinId] / Const.pixdim1) * Const.pixdim1;
	float pos_z = round(z[SpinId] / Const.pixdim2) * Const.pixdim2;

	for (int i = 0; i < Const.nx; i++)
		{
			E0 = wx(Const.T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Const);

			x0 = x[SpinId];
			y0 = y[SpinId];
			z0 = z[SpinId];

			x[SpinId] = x0 - (x0 - pos_x + Const.pixdim0 * 0.5) * Const.delta + curand_uniform(&s)* Const.delta;
			y[SpinId] = y0 - (y0 - pos_y + Const.pixdim0 * 0.5) * Const.delta + curand_uniform(&s)* Const.delta;
			z[SpinId] = z0 - (z0 - pos_z + Const.pixdim0 * 0.5) * Const.delta + curand_uniform(&s)* Const.delta;

			E1 = wx(Const.T)*Ei_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Const);

			p = fminf(1., expf((E0 - E1) / Const.T));

			if (curand_uniform(&s) < p)
				;
			else
			{
				x[SpinId] = x0;
				y[SpinId] = y0;
				z[SpinId] = z0;
			}
		}
}

__global__ void mc(float* x, float* y, float* z, unsigned long long* c, float* ten, int* sig, int* VoxIds, int LatticeId, Constant Const)
{
	int MaxThreadId = getMaxThreadId(LatticeId,Const);
	int threadId = threadIdx.x + blockIdx.x * blockDim.x;

	while (threadId < MaxThreadId)
	{
		int VoxNum = threadIdToVoxNum(threadId, LatticeId, Const);

		if (VoxIds[VoxNum] >= 0)
		{
			int SpinsInVoxel;

			if (VoxNum < Const.VoxCount - 1)
			{
				int next = VoxNum + 1;
				while (VoxIds[next] < 0) next++;
				SpinsInVoxel = VoxIds[next] - VoxIds[VoxNum];
			}
			else
				SpinsInVoxel = Const.SpinCount - VoxIds[VoxNum];

			for (int SpinId = VoxIds[VoxNum]; SpinId < VoxIds[VoxNum] + SpinsInVoxel; SpinId++)
			{
				mc_c(x, y, z, c, ten, sig, VoxIds, VoxNum, VoxIds[VoxNum], SpinId, Const);
				mc_x(x, y, z, c, ten, sig, VoxIds, VoxNum, SpinId, Const);
			}
		}

		threadId += gridDim.x*blockDim.x;
	}
}