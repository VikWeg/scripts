void fiber_tracking_cube()
{
	import_and_init_data();
	create_log();
	init_Efile();
	stats_out();

	cudaMalloc((void**)&dev_in_ensemble, scount*sizeof(vertex));
	cudaMalloc((void**)&dev_out_ensemble, scount*sizeof(vertex));

	cudaMemcpy(dev_in_ensemble, ensemble, scount*sizeof(vertex), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_ensemble, ensemble, scount*sizeof(vertex), cudaMemcpyHostToDevice);

	dim3 grid(scount);

	RDTSC(start_all);
	float T = Ti;
	for (; T > Tf; T *= etha)
	{
		std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";

		/*========*/ RDTSC(start); /*========*/

		mc <<<grid, 1>>> (dev_in_ensemble, dev_out_ensemble);

		/*========*/ RDTSC(stop); /*========*/

		cudaMemcpy(ensemble, dev_out_ensemble, scount*sizeof(vertex), cudaMemcpyDeviceToHost);

		export_fibers(tstep);
		write_E_file(tstep, T);
		plot_E();
		time_out();
	}
	RDTSC(stop_all);

	write_par_file(T / etha);
}