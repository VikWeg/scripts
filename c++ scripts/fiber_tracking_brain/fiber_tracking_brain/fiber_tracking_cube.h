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

	std::cout << "Annealing\n";
	for (float T = Ti; T > Tf; T *= etha)
	{
		RDTSC(start);
		for (int n = 0; n < S * scount; n++)
		{
			i = 0; j = 0; k = 0;
			while (wmask[i][j][k] == 0 && surf_mask[i][j][k] == 0)
			{
				i = u_i(generate);
				j = u_j(generate);
				k = u_k(generate);
			}

			std::uniform_int_distribution<int> u_s(0, snum[i][j][k]-1);
			s = u_s(generate);

			mc_c(&ensemble[i][j][k][s], T);
			mc_x(&ensemble[i][j][k][s], T);
		}
		RDTSC(stop);

		std::cout << std::setprecision(2) << "T = " << T << " Number of cycles = " << (stop.int64 - start.int64)/(S*scount) << "\xd";


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
	}
	Efile << "0}";
	Efile.close();

	std::cout << "\nStarting Fiber Export\n";
	export_fibers();
	std::cout << "Finished Fiber Export\n";
}