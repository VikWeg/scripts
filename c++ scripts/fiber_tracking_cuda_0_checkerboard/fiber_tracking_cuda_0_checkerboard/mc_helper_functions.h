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

int getBit(int i, long long int c)
{
	long long int mask = 1 << i;

	return (c&mask == 0) ? 0 : 1;
}

// **** Change i-th bit ****

int changeBit(int i, long long int c)
{
	long long int mask = 1 << i;

	return c ^ mask;

}
// **** Sum all bits **** 'Hamming Weight', 'popcount' or 'sideways addition' http://bisqwit.iki.fi/source/misc/bitcounting/

int sumBits(long long int c)
{
	const unsigned TEST_BITS = sizeof(long long int)* sizeof(char)* 8;

	long long int m1 = (~(long long int)0) / 3;
	long long int m2 = (~(long long int)0) / 5;
	long long int m4 = (~(long long int)0) / 17;
	long long int h01 = (~(long long int)0) / 255;

	c -= (c >> 1) & m1;
	c = (c & m2) + ((c >> 2) & m2);
	c = (c + (c >> 4)) & m4;

	return (c * h01) >> (TEST_BITS - 8);
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
				if (CurrentNeighborNum = NeighborNum) return (z0 + k)*c0*c1 + (y0 + j)*c0 + (x0 + i);
				else CurrentNeighborNum++;

	return -1;
}

int GetNonEmptyNeighborVoxsCount(int BaseVoxNum, int* VoxIds) //TO BE REPLACE
{//Returns the number of voxels next to the voxel BaseVoxNum, that have voxId != -1

	int NonEmptyNeighborVoxsCount = 0;

	for (int NeighborNum = 0; NeighborNum < 26; NeighborNum++)
		NonEmptyNeighborVoxsCount += VoxIds[GetVoxNumFromNeighborNum(BaseVoxNum, NeighborNum, VoxIds)] >= 0 ? 1 : 0;

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



}