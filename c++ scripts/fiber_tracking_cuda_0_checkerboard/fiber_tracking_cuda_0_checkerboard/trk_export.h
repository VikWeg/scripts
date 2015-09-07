void export_fibers(int* VoxIds, int n)
{
	chdir("fibers");
	fiber_file = fopen(("fibers_" + std::to_string(n) + ".trk").c_str(),"w+b");

	write_header(fiber_file);

	int len, fiber_num = 0;
	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	{
		int VoxNum = k*cube_size[0] * cube_size[1] + j*cube_size[0] + i;

		if (sig[VoxNum])
		{
			int VoxId = VoxIds[VoxNum];

			int	SpinsInVoxel;
			if (VoxNum < cube_size[0] * cube_size[1] * cube_size[2] - 1)
			{
				int next = VoxNum + 1;
				while (VoxIds[next] < 0) next++;
				SpinsInVoxel = VoxIds[next] - VoxId;
			}
			else SpinsInVoxel = scount - VoxId;

			for (int SpinId = VoxId; SpinId < VoxId + SpinsInVoxel; SpinId++)
			{
				len = get_fiber_length(VoxIds,VoxNum, SpinId);

				if (len>1)
				{
					fwrite((char*)&len, 1, 4, fiber_file);
					get_fiber(VoxIds,VoxNum, SpinId);
					fiber_num++;
				}
			}
		}
	}

	fseek(fiber_file, 1000 - 12, SEEK_SET);
	fwrite((char*)&fiber_num, 1, 4, fiber_file);

	fclose(fiber_file);
}