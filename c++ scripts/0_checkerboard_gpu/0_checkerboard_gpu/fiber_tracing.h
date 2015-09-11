void get_fiber(int* VoxIds, int PrevVoxNum, int PrevSpinId, int CurrVoxNum, int CurrSpinId, int n, Constant Const)
{
	Next Next = { -1, -1 };

	int nc = GetNeighborSpinCount(CurrVoxNum, VoxIds,Const);

	for (int i = 0; i < nc; i++)
	if (getBit(i, c[CurrSpinId]))
	{
		Next = GetNextSpin(VoxIds, CurrVoxNum, i, Const);
		if (Next.VoxNum != PrevVoxNum) break;
		else Next = { -1, -1 };
	};

	fwrite((char*)&(x[CurrSpinId]), 1, 4, fiber_file);
	fwrite((char*)&(y[CurrSpinId]), 1, 4, fiber_file);
	fwrite((char*)&(z[CurrSpinId]), 1, 4, fiber_file);

	if (Next.VoxNum != -1 && n < 50)
	{
		if (sig[Next.VoxNum] != 1)
			get_fiber(VoxIds, CurrVoxNum, CurrSpinId, Next.VoxNum, Next.SpinId, n + 1, Const);
		else
		{
			fwrite((char*)&(x[Next.SpinId]), 1, 4, fiber_file);
			fwrite((char*)&(y[Next.SpinId]), 1, 4, fiber_file);
			fwrite((char*)&(z[Next.SpinId]), 1, 4, fiber_file);
		}
	}
}

void get_fiber(int* VoxIds, int VoxNum, int SpinId,Constant Const)
{
	Next Next = { -1, -1 };

	int nc = GetNeighborSpinCount(VoxNum, VoxIds, Const);

	for (int i = 0; i < nc; i++)
	if (getBit(i, c[SpinId])) { Next = GetNextSpin(VoxIds, VoxNum, i, Const); break; };

	if (Next.VoxNum != -1 && sig[Next.VoxNum] != 1)
	{
		fwrite((char*)&(x[SpinId]), 1, 4, fiber_file);
		fwrite((char*)&(y[SpinId]), 1, 4, fiber_file);
		fwrite((char*)&(z[SpinId]), 1, 4, fiber_file);

		get_fiber(VoxIds, VoxNum, SpinId, Next.VoxNum, Next.SpinId, 1, Const);
	}
	else if (Next.VoxNum != -1 && sig[Next.VoxNum] == 1)
	{
		fwrite((char*)&(x[SpinId]), 1, 4, fiber_file);
		fwrite((char*)&(y[SpinId]), 1, 4, fiber_file);
		fwrite((char*)&(z[SpinId]), 1, 4, fiber_file);

		fwrite((char*)&(x[Next.SpinId]), 1, 4, fiber_file);
		fwrite((char*)&(y[Next.SpinId]), 1, 4, fiber_file);
		fwrite((char*)&(z[Next.SpinId]), 1, 4, fiber_file);
	}
}

int get_fiber_length(int* VoxIds, int PrevVoxNum, int PrevSpinId, int CurrVoxNum, int CurrSpinId, int n,Constant Const)
{
	Next Next = { -1, -1 };

	int nc = GetNeighborSpinCount(CurrVoxNum, VoxIds, Const);

	for (int i = 0; i < nc; i++)
		if (getBit(i, c[CurrSpinId]))
		{ 
			Next = GetNextSpin(VoxIds, CurrVoxNum, i, Const);
			if (Next.VoxNum != PrevVoxNum) break;
			else Next = { -1, -1 };
		};

	if (Next.VoxNum != -1 && n < 50)
	{
		if (sig[Next.VoxNum] != 1)
			return get_fiber_length(VoxIds, CurrVoxNum, CurrSpinId, Next.VoxNum, Next.SpinId, n + 1, Const);
		else
			return n + 2;
	}
	else
		return n + 1;
}

int get_fiber_length(int* VoxIds, int VoxNum, int SpinId,Constant Const)
{
	Next Next = { -1, -1 };

	int nc = GetNeighborSpinCount(VoxNum, VoxIds, Const);

	for (int i = 0; i < nc; i++)
	if (getBit(i, c[SpinId])) { Next = GetNextSpin(VoxIds, VoxNum, i, Const); break; };

	if (Next.VoxNum != -1 && sig[Next.VoxNum] != 1)
		return get_fiber_length(VoxIds, VoxNum, SpinId, Next.VoxNum, Next.SpinId, 1, Const);
	else if (Next.VoxNum != -1 && sig[Next.VoxNum] == 1)
		return 2;
	else
		return 1;
}

// GetNextSpin ?
// NextVOxNum != PrevVoxNUm ?
//getBit ?