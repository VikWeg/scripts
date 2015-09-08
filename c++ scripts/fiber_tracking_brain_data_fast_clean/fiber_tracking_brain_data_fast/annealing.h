void fiber_tracking_cube()
{
	import_and_init_data();
	create_log();
	init_Efile();
	stats_out();

	RDTSC(start_all);
	float T = Ti;
	for (; T > Tf; T *= etha)
	{
		std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";

		/*========*/ RDTSC(start); /*========*/
		for (int n = 1; n <= S; n++)
		{
			for (int i = 0; i < cube_size[0]; i++)
			for (int j = 0; j < cube_size[1]; j++)
			for (int k = 0; k < cube_size[2]; k++)
			for (int s = 0; s < snum[i][j][k]; s++)
			{
				mc_c(&ensemble[i][j][k][s], T);
				mc_x(&ensemble[i][j][k][s], T);
			}
			
			std::cout << "% done: " << std::setprecision(2) << std::setw(6) << std::left << std::fixed << (100. * n) / S << "\r";
		}
		/*========*/ RDTSC(stop); /*========*/

		export_fibers(tstep);
		write_E_file(tstep,T);
		plot_E();
		time_out();
	}
	RDTSC(stop_all);

	write_par_file(T/etha);
}