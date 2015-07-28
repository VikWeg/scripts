void fiber_tracking_cube()
{
	import_data();
	std::cout << "Data Import done\n";
	
	init_ensemble();
	std::cout << "Ensemble Initialization done\n";

	int i=0, j=0, k=0, s, scount=0;

	for (int i = 0; i < size[0]; i++)
	for (int j = 0; j < size[1]; j++)
	for (int k = 0; k < size[2]; k++)
		if (wmask[i][j][k] == 1 || surf_mask[i][j][k] == 1)
		for (int s = 0; s < snum[i][j][k]; s++)
			scount++;
	
	float EiC;
	float EiD;
	float EiI;
	std::ofstream Efile;
	Efile.open("C:/ETH/Neuro/GlobalTracking/energy_.txt");
	Efile << "{";

	std::cout << "Annealing\n\n";
	int tsteps_tot = ceilf(logf(Tf / Ti) / logf(etha));
	int tstep = 1;
	for (float T = Ti; T > Tf; T *= etha)
	{
		std::cout << std::setprecision(2) << "Tstep = " << tstep << "/" << tsteps_tot << "\n";
		RDTSC(start);
		for (int n = 0; n < S * scount; n++)
		{
			i = u_i(generate);
			j = u_j(generate);
			k = u_k(generate);
			while (wmask[i][j][k] == 0 && surf_mask[i][j][k] == 0)
			{
				i = u_i(generate);
				j = u_j(generate);
				k = u_k(generate);
			}

			for (int s = 0; s < snum[i][j][k]; s++)
			{
				mc_c(&ensemble[i][j][k][s], T);
				mc_x(&ensemble[i][j][k][s], T);
			}

			std::cout << "% done: " << std::setprecision(2) << std::setw(6) << std::left << std::fixed << (100. * n) / (S*scount) << "\r";
		}
		RDTSC(stop);

		float time = (float)(stop.int64 - start.int64) / CPU_freq / 60;
		std::cout << std::setprecision(2) << "\nElapsed time = " << std::fixed << time << " min\n";
		std::cout << "Expected time left: " << (tsteps_tot - tstep)*(int)time / 60 << " h : " << fmod((tsteps_tot - tstep)*time , 60) << " min\n";
		tstep++,

		EiC = 0;
		EiD = 0;
		EiI = 0;
		for (int i = 0; i < size[0]; i++)
		for (int j = 0; j < size[1]; j++)
		for (int k = 0; k < size[2]; k++)
		for (int s = 0; s < snum[i][j][k]; s++)
		if (wmask[i][j][k] == 1 || surf_mask[i][j][k] == 1)
		{
			EiC += Ei_C(&ensemble[i][j][k][s]);
			EiD += Ei_D(&ensemble[i][j][k][s]);
			EiI += Ei_I(&ensemble[i][j][k][s]);
		}
		Efile << std::fixed << EiC / scount << "," << EiD / scount << "," << EiI / scount << ",";

		std::cout << std::setprecision(3) << "Econstr = " << EiC / scount << " Edata = " << EiD / scount << " Eint = " << EiI / scount << "\n\n";
	}

	long pos = Efile.tellp();
	Efile.seekp(pos - 1);
	Efile << "}";
	Efile.close();

	std::cout << "\nStarting Fiber Export\n";
	export_fibers();
	std::cout << "Finished Fiber Export\n";
}