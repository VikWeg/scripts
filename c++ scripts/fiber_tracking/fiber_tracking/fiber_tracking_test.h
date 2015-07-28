void fiber_tracking_test()
{
	std::cout << "Ensemble Initialization\n";
	init_test_ensemble();

	int i=0, j=0, k, s;

	float EiC;
	float EiD;
	float EiI;

	std::uniform_int_distribution<int> u_i_t(0, 9);
	std::uniform_int_distribution<int> u_j_t(0, 9);

	std::ofstream Efile;
	Efile.open("C:/ETH/Neuro/GlobalTracking/energy_cross_.txt");
	Efile << "{";

	std::cout << "Annealing\n";
	for (float T = Ti; T > Tf; T *= etha)
	{
		for (int n = 0; n < S * 2 * 64; n++)
		{
			i = 0; j = 0;
			while ((i<3 && (j<3 || j>6)) || (i>6 && (j<3 || j>6)))
			{
				i = u_i_t(generate);
				j = u_j_t(generate);
			}

			k = 0;

			std::uniform_int_distribution<int> u_s_t(0, snum[i][j][k]-1);
			s = u_s_t(generate);

			mc_c(&ensemble[i][j][k][s], T);
			mc_x(&ensemble[i][j][k][s], T);
		}
		std::cout << std::setprecision(2) << "T = " << T << "\xd";

		EiC = 0;
		EiD = 0;
		EiI = 0;
		for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
		for (int k = 0; k < 1;  k++)
		for (int s = 0; s < snum[i][j][k];  s++)
		{
			if (ensemble[i][j][k][s].sig >= 0)
			{
				EiC += Ei_C(&ensemble[i][j][k][s]);
				EiD += Ei_D(&ensemble[i][j][k][s]);
				EiI += Ei_I(&ensemble[i][j][k][s]);
			}
		}
		Efile << std::fixed << EiC / (2 * 64) << "," << EiD / (2 * 64) << "," << EiI / (2 * 64) << ",";
	}
	Efile << "0}";
	Efile.close();

	std::cout << "\nStarting Fiber Export\n";
	export_test_fibers();
	std::cout << "Finished Fiber Export\n";
}