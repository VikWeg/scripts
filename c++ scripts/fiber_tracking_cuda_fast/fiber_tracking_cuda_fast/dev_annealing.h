void annealing()
{
	import_and_init_data();
	create_log();
	init_Efile();
	stats_out();

	RDTSC(start_all);
	float T = Ti;
	int toggle = 0;
	for (; T > Tf; T *= etha)
	{
		std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";

		/*========*/ RDTSC(start); /*========*/

		for (int i = 1; i <= sweeps; i++)
		{
			if (toggle == 0)
				mc << <scount / 32, 64 >> > (dev_in_, dev_out_ensemble, T, scount);
			else
				mc << <scount / 32, 64 >> > (dev_out_ensemble, dev_in_ensemble, T, scount);

			toggle = 1 - toggle;
		}

		/*========*/ RDTSC(stop); /*========*/

		export_fibers(tstep);
		write_E_file(tstep, T);
		plot_E();
		time_out();
	}
	RDTSC(stop_all);

	write_par_file(T / etha);
}