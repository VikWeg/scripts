void annealing()
{
	import_and_init_data();
	create_log();
	init_Efile();
	stats_out();

		RDTSC(start_all);
		for (T = Ti; T > Tf; T *= etha)
		{
			Const.T = T;
			std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";

			/*========*/ RDTSC(start); /*========*/
			for (int n = 1; n <= S; n++)
			{
				for (int lattice_id = 0; lattice_id < 27; lattice_id++)
					mc << <10, 32 >> >(dev_x, dev_y, dev_z, dev_c, dev_ten, dev_sig, dev_VoxIds, lattice_id, Const);

				std::cout << "% done: " << std::setprecision(2) << std::setw(6) << std::left << std::fixed << (100. * n) / S << "\r";
			}
			/*========*/ RDTSC(stop); /*========*/

			cudaDeviceSynchronize();
			CopyDeviceToHost();
			export_fibers(VoxIds,tstep,Const);
			write_E_file(tstep, T);
			plot_E();
			time_out();
		}
		RDTSC(stop_all);

		write_par_file(T / etha);
		LogFile.close();
		cudaDeviceReset();
}