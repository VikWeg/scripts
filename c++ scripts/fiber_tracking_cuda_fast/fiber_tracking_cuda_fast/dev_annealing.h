void annealing()
{
	import_and_init_data();
	create_log();
	init_Efile();
	stats_out();

	cudaEvent_t dev_start, dev_stop;
	cudaEventCreate(&dev_start);
	cudaEventCreate(&dev_stop);

	RDTSC(start_all);

	float T = Ti;
	for (; T > Tf; T *= etha)
	{
		std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";

		/*========*/ RDTSC(start); /*========*/

		for (int i = 1; i <= sweeps; i++)
		{
			cudaEventRecord(dev_start, 0);
			if (toggle == 0)
				mc << <scount / 32, 64 >> >
				(
					dev_in_x, dev_in_y, dev_in_z, dev_in_cc, dev_in_c,

					dev_out_x, dev_out_y, dev_out_z, dev_out_cc, dev_out_c,

					dev_pos_x, dev_pos_y, dev_pos_z,
					dev_T0, dev_T1, dev_T2, dev_T3, dev_T4, dev_T5,
					dev_Emin, dev_Emax, dev_delta_E,
					dev_sig,
					dev_nc,
					dev_n_id, dev_n,

					T, scount
				);
			else
				mc << <scount / 32, 64 >> >
				(
					dev_out_x, dev_out_y, dev_out_z, dev_out_cc, dev_out_c,

					dev_in_x, dev_in_y, dev_in_z, dev_in_cc, dev_in_c,

					dev_pos_x, dev_pos_y, dev_pos_z,
					dev_T0, dev_T1, dev_T2, dev_T3, dev_T4, dev_T5,
					dev_Emin, dev_Emax, dev_delta_E,
					dev_sig,
					dev_nc,
					dev_n_id, dev_n,

					T, scount
				);
			cudaEventRecord(dev_stop, 0);
			cudaEventSynchronize(dev_stop);

			toggle = 1 - toggle;
		}

		/*========*/ RDTSC(stop); /*========*/

		cpyEnsembleDevToHost();
		export_fibers(tstep);
		write_E_file(tstep, T);
		plot_E();
		time_out();
	}
	RDTSC(stop_all);

	write_par_file(T / etha);
}