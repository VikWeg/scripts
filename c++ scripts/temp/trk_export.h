void export_fibers(int n)
{
	chdir("fibers");
	fiber_file = fopen(("fibers_" + std::to_string(n) + ".trk").c_str(),"w+b");

	trackvis_header trk_hdr;
	write_header(&trk_hdr, fiber_file);

	int len, fiber_num = 0;
	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	for (int s = 0; s < snum[offset(i, j, k)]; s++)
	if (ensemble[vertex_offset(i, j, k) + s].sig)
	{
		len = get_fiber_length(&ensemble[vertex_offset(i, j, k) + s]);

		if (len>1)
		{
			fwrite((char*)&len, 1, 4, fiber_file);
			get_fiber(&ensemble[vertex_offset(i, j, k) + s]);
			fiber_num++;
		}
	}

	fseek(fiber_file, 1000 - 12, SEEK_SET);
	fwrite((char*)&fiber_num, 1, 4, fiber_file);

	fclose(fiber_file);
}