// **** Map Thread to Vox ****
int threadIdToVoxNum(int threadId, int latticeId)
{
	int sliceOff = latticeId / 9;
	int colOff = latticeId % 3;
	int rowOff = (latticeId / 3) % 3;

	int voxsPerRow = (cube_size[0] - colOff - 1) / 3 + 1;
	int rowsPerSlice = (cube_size[1] - colOff - 1) / 3 + 1;

	int voxsPerSlice = voxsPerRow*rowsPerSlice;

	int voxsToFill = threadId + 1;

	int voxsLeft = voxsToFill % voxsPerSlice;

	int lastRowLeft = voxsLeft % voxsPerRow;

	return	cube_size[0] * cube_size[1] * sliceOff

			+cube_size[0] * cube_size[1] * 3 * (voxsToFill / voxsPerSlice)

			+(voxsLeft == 0) ? 1 : 0 *
									(
										rowOff*cube_size[0]
										+(rowsPerSlice - 1) * 3 * cube_size[0]
										+colOff
										+(voxsPerRow - 1) * 3
										-3 * cube_size[0] * cube_size[1]
									)

			+(voxsLeft != 0) ? 1 : 0 *
									(
										rowOff*cube_size[0]
										+3 * cube_size[0] * (voxsLeft / voxsPerRow)
										-(lastRowLeft == 0) ? 1 : 0 *
																	(
																		2 * cube_size[0]
																		+(cube_size[0] - colOff - 1) % 3
																		+ 1
																	)

										+(lastRowLeft != 0) ? 1 : 0 *
																	(
																		colOff
																		+3 * (lastRowLeft - 1)
																	)
									);

}

// **** Get i-th bit ****

int getBit(int i, unsigned long long c)
{
	unsigned long long mask = 1 << i;

	return ((c&mask) == 0) ? 0 : 1;
}

// **** Change i-th bit ****

int changeBit(int i, unsigned long long c)
{
	unsigned long long mask = 1 << i;

	return c ^ mask;

}
// **** Sum all bits **** 'Hamming Weight', 'popcount' or 'sideways addition' http://bisqwit.iki.fi/source/misc/bitcounting/

int sumBits(unsigned long long c)
{
	const unsigned TEST_BITS = sizeof(unsigned long long)* sizeof(char)* 8;

	unsigned long long m1 = (~(unsigned long long)0) / 3;
	unsigned long long m2 = (~(unsigned long long)0) / 5;
	unsigned long long m4 = (~(unsigned long long)0) / 17;
	unsigned long long h01 = (~(unsigned long long)0) / 255;

	c -= (c >> 1) & m1;
	c = (c & m2) + ((c >> 2) & m2);
	c = (c + (c >> 4)) & m4;

	return (c * h01) >> (TEST_BITS - 8);
}

// Find first set bit in c https://en.wikipedia.org/wiki/Find_first_set

int FindFirstSet(unsigned long long c)
{
	if (c == 0) return 0;
	else return sumBits(c ^ (~(-c)));
}

int GetVoxNumFromNeighborNum(int VoxNum, int NeighborNum, int* VoxIds)
{//Returns the voxel number of the NeighborNum-th neighbor of the voxel VoxNum

	int c0 = cube_size[0];
	int c1 = cube_size[1];
	int c2 = cube_size[2];

	int x = VoxNum % c0;
	int y = (VoxNum / c0) % c1;
	int z = VoxNum / (c0*c1);

	int x0 = x - 1;
	int y0 = y - 1;
	int z0 = z - 1;

	int i, j, k;
	int CurrentNeighborNum = 0;
	for (k = 0; k < 3;k++)
		for (j = 0; j < 3; j++)
			for (i = 0; i < 3; i++)
				if (
						x0 + i >= 0 && x0 + i < c0 && y0 + j >= 0 && y0 + j < c1 && z0 + k >= 0 && z0 + k < c2
					&&	(x0 + i != x || y0 + j != y || z0 + k != z)
					&&	VoxIds[(z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i)] >= 0
					)
				if (CurrentNeighborNum == NeighborNum) return (z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i);
				else CurrentNeighborNum++;

	return -1;
}

int GetNonEmptyNeighborVoxsCount(int BaseVoxNum, int* VoxIds)
{//Returns the number of voxels next to the voxel BaseVoxNum, that have voxId > -1

	int NonEmptyNeighborVoxsCount = 0;

	int c0 = cube_size[0];
	int c1 = cube_size[1];
	int c2 = cube_size[2];

	int x = BaseVoxNum % c0;
	int y = (BaseVoxNum / c0) % c1;
	int z = BaseVoxNum / (c0*c1);

	int x0 = x - 1;
	int y0 = y - 1;
	int z0 = z - 1;

	int i, j, k;
	for (k = 0; k < 3; k++)
	for (j = 0; j < 3; j++)
	for (i = 0; i < 3; i++)
	if (
		x0 + i >= 0 && x0 + i < c0 && y0 + j >= 0 && y0 + j < c1 && z0 + k >= 0 && z0 + k < c2
		&& (x0 + i != x || y0 + j != y || z0 + k != z)
		&& VoxIds[(z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i)] >= 0
		)
		NonEmptyNeighborVoxsCount++;

	return NonEmptyNeighborVoxsCount;
}

int GetNeighborNumFromVoxNum(int BaseVoxNum, int VoxNum, int* VoxIds)
{//Returns the neighbor number of the voxel VoxNum, relative to the voxel BaseVoxNum. Assumes that BaseVox and Vox are both valid voxels.

	int c0 = cube_size[0];
	int c1 = cube_size[1];
	int c2 = cube_size[2];

	int x = BaseVoxNum % c0;
	int y = (BaseVoxNum / c0) % c1;
	int z = BaseVoxNum / (c0*c1);

	int x0 = x - 1;
	int y0 = y - 1;
	int z0 = z - 1;

	int i, j, k;
	int NeighborNum = 0;
	for (k = 0; k < 3; k++)
		for (j = 0; j < 3; j++)
			for (i = 0; i < 3; i++)
				if (
						x0 + i >= 0 && x0 + i < c0 && y0 + j >= 0 && y0 + j < c1 && z0 + k >= 0 && z0 + k < c2
						&& (x0 + i != x || y0 + j != y || z0 + k != z)
						&& VoxIds[(z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i)] >= 0
					)
				if ((z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i) == VoxNum) return NeighborNum;
				else NeighborNum++;

	return -1;
}

int GetConnectivityOffset(int VoxNum, int NeighborNum, int* VoxIds)
{//Returns the offset of the NeighborNum-th voxel in the connectivity of voxel VoxNum

	int c0 = cube_size[0];
	int c1 = cube_size[1];
	int c2 = cube_size[2];

	int x = VoxNum % c0;
	int y = (VoxNum / c0) % c1;
	int z = VoxNum / (c0*c1);

	int x0 = x - 1;
	int y0 = y - 1;
	int z0 = z - 1;

	int i, j, k;
	int CurrentNeighborNum = 0;
	int Offset = 0;
	for (k = 0; k < 3; k++)
		for (j = 0; j < 3; j++)
			for (i = 0; i < 3; i++)
			if (
				x0 + i >= 0 && x0 + i < c0 && y0 + j >= 0 && y0 + j < c1 && z0 + k >= 0 && z0 + k < c2
				&& (x0 + i != x || y0 + j != y || z0 + k != z)
				&& VoxIds[(z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i)] >= 0
				)
				{
					if (CurrentNeighborNum == NeighborNum) return Offset;

					int vox = (z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i);

					if (vox < c0*c1*c2 - 1)
					{
						int next = vox + 1;
						while (VoxIds[next] < 0) next++;
						Offset += VoxIds[next] - VoxIds[vox];
					}
					else Offset += scount - VoxIds[vox];

					CurrentNeighborNum++;
				}

	return -1;
}

int GetNeighborSpinCount(int VoxNum, int* VoxIds)
{//Returns number of spins around VoxNum
	int c0 = cube_size[0];
	int c1 = cube_size[1];
	int c2 = cube_size[2];

	int x = VoxNum % c0;
	int y = (VoxNum / c0) % c1;
	int z = VoxNum / (c0*c1);

	int x0 = x - 1;
	int y0 = y - 1;
	int z0 = z - 1;

	int i, j, k;
	int count = 0;
	for (k = 0; k < 3; k++)
	for (j = 0; j < 3; j++)
	for (i = 0; i < 3; i++)
	if (
		x0 + i >= 0 && x0 + i < c0 && y0 + j >= 0 && y0 + j < c1 && z0 + k >= 0 && z0 + k < c2
		&& (x0 + i != x || y0 + j != y || z0 + k != z)
		&& VoxIds[(z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i)] >= 0
		)
		{
			int vox = (z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i);

			if (vox < c0*c1*c2 - 1)
			{
				int next = vox + 1;
				while (VoxIds[next] < 0) next++;
				count += VoxIds[next] - VoxIds[vox];
			}
			else count += scount - VoxIds[vox];
		}

	return count;
}

Next GetNextSpin(int* VoxIds, int VoxNum, int NextSpinOffset)
{//Returns the voxnum and spin id of the spin with spin offset NextSpinOffset relative to VoxNum

	Next Next = { -1, -1 };
	int SpinsInVoxel;

	int c0 = cube_size[0];
	int c1 = cube_size[1];
	int c2 = cube_size[2];

	int x = VoxNum % c0;
	int y = (VoxNum / c0) % c1;
	int z = VoxNum / (c0*c1);

	int x0 = x - 1;
	int y0 = y - 1;
	int z0 = z - 1;

	int i, j, k;
	int Offset = 0;
	for (k = 0; k < 3; k++)
	for (j = 0; j < 3; j++)
	for (i = 0; i < 3; i++)
	if (
		x0 + i >= 0 && x0 + i < c0 && y0 + j >= 0 && y0 + j < c1 && z0 + k >= 0 && z0 + k < c2
		&& (x0 + i != x || y0 + j != y || z0 + k != z)
		&& VoxIds[(z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i)] >= 0
		)
		{
			int voxnum = (z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i);

			if (voxnum < c0*c1*c1 - 1)
			{
				int next = voxnum + 1;
				while (VoxIds[next] < 0) next++;
				SpinsInVoxel = VoxIds[next] - VoxIds[voxnum];
			}
			else SpinsInVoxel = scount - VoxIds[voxnum];

			for (int SpinNumber = 0; SpinNumber < SpinsInVoxel; SpinNumber++)
			{
				if (Offset == NextSpinOffset)
				{
					Next.VoxNum = voxnum;
					Next.SpinId = VoxIds[voxnum] + SpinNumber;
					return Next;
				}

				Offset++;
			}

		}

	return Next;
}