void annealing()
{
	import_and_init_data();
	create_log();
	init_Efile();
	stats_out();

	if (TestFiberTrace == 1)
		FiberTraceTest();
	else
	{

		RDTSC(start_all);
		T = Ti;
		for (; T > Tf; T *= etha)
		{
			std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";

			/*========*/ RDTSC(start); /*========*/
			for (int n = 1; n <= S; n++)
			{
				//for (int lattice_id = 0; lattice_id < 27; lattice_id++)
				mc(x, y, z, c, ten, sig, VoxIds, 0);

				std::cout << "% done: " << std::setprecision(2) << std::setw(6) << std::left << std::fixed << (100. * n) / S << "\r";
			}
			/*========*/ RDTSC(stop); /*========*/

			export_fibers(VoxIds,tstep);
			write_E_file(tstep, T);
			plot_E();
			time_out();
		}
		RDTSC(stop_all);

		write_par_file(T / etha);
		LogFile.close();
		cudaDeviceReset();
	}
}