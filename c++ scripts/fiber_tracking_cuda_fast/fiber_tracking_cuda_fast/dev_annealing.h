void annealing()
{
	import_and_init_data();
	create_log();
	init_Efile();
	stats_out();

	RDTSC(start_all);
	float T = Ti;
	int toggle = 1;
	for (; T > Tf; T *= etha)
	{
		std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";

		/*========*/ RDTSC(start); /*========*/

		for (int i = 1; i <= sweeps; i++)
		if (toggle % 2)
		{
			mc << <scount/32, 64 >> > (dev_in_ensemble, dev_out_ensemble, T,scount);
			toggle++;
		}
		else
		{
			mc << <scount/32, 64>> > (dev_out_ensemble, dev_in_ensemble, T,scount);
			toggle++;
		}

		/*========*/ RDTSC(stop); /*========*/

		cpyDev2Host(toggle);

		export_fibers(tstep);
		write_E_file(tstep, T);
		plot_E();
		time_out();
	}
	RDTSC(stop_all);

	write_par_file(T / etha);
}