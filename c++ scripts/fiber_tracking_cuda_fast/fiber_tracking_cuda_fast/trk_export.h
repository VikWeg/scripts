void export_fibers(int n)
{
	chdir("fibers");
	fiber_file = fopen(("fibers_" + std::to_string(n) + ".trk").c_str(),"w+b");

	trackvis_header trk_hdr;
	write_header(&trk_hdr, fiber_file);

	int len, fiber_num = 0;
	for (int s = 0; s < scount; s++)
	if (sig[s])
	{
		len = get_fiber_length(s);

		if (len>1)
		{
			fwrite((char*)&len, 1, 4, fiber_file);
			get_fiber(s);
			fiber_num++;
		}
	}

	fseek(fiber_file, 1000 - 12, SEEK_SET);
	fwrite((char*)&fiber_num, 1, 4, fiber_file);

	fclose(fiber_file);
}