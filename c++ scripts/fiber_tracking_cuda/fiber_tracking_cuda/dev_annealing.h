void annealing()
{
	import_and_init_data();
	create_log();
	init_Efile();
	stats_out();

	dim3 grid(scount);

	RDTSC(start_all);
	float T = Ti;
	for (; T > Tf; T *= etha)
	{
		std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";

		/*========*/ RDTSC(start); /*========*/

		for (int i = 1; i <= sweeps; i++)
		if (i%2)
			mc <<<grid, 1>>> (dev_out_ensemble, dev_in_ensemble,T);
		else
			mc <<<grid, 1 >> > (dev_in_ensemble, dev_out_ensemble,T);

		/*========*/ RDTSC(stop); /*========*/

		//if (sweeps%2)
		//cudaMemcpy(ensemble, dev_in_ensemble, scount*sizeof(vertex), cudaMemcpyDeviceToHost);
		//else
		//cudaMemcpy(ensemble, dev_out_ensemble, scount*sizeof(vertex), cudaMemcpyDeviceToHost);

		export_fibers(tstep);
		write_E_file(tstep, T);
		plot_E();
		time_out();
	}
	RDTSC(stop_all);

	write_par_file(T / etha);
}