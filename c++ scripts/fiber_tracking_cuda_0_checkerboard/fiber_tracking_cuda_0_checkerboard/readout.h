float EiC, EiD, EiI;
int tstep = 1;

void init_Efile()
{
	Efile.open("energy.dat");
	Efile << "Tstep\tTemp\tEc\tEd\tEi\tEt\n";
	Efile.close();
}

void calc_E(float T)
{
	EiC = 0; EiD = 0; EiI = 0;
	for (int i = 0; i < cube_size[0]; i++)
	for (int j = 0; j < cube_size[1]; j++)
	for (int k = 0; k < cube_size[2]; k++)
	{
		int VoxNum = k*cube_size[0] * cube_size[1] + j*cube_size[0] + i;
		EiC += Ei_C(VoxNum);
		EiD += Ei_D(VoxNum);
		EiI += Ei_I(VoxNum);
	}

	//EiC *= wc(T);
	//EiD *= wx(T);
	//EiI *= wx(T);
}

void stats_out()
{
	std::cout << "Annealing\n\n";
	std::cout << "WM Voxels: \t" << w_vox_num;
	std::cout << "\nSurf Voxels: \t" << surf_vox_num;
	std::cout << "\nVertices: \t" << scount << "\n\n";
}

void time_out()
{
	float time = (float)(stop.int64 - start.int64) / CPU_freq / 60;
	std::cout << std::setprecision(2) << "\nElapsed time = " << std::fixed << time << " min\n";
	std::cout << "Expected time left: " << (tsteps_tot - tstep)*(int)time / 60 << " h : " << fmod((tsteps_tot - tstep)*time, 60) << " min\n";
	tstep++;
}

void import_and_init_data()
{
	import_data();
	std::cout << "Data Import done\n";
	init_ensemble();
	std::cout << "Ensemble Initialization done\n";
}

void write_E_file(int t, float T)
{
	chdir("..\\");
	calc_E(T);
	Efile.open("energy.dat", std::ofstream::app);

	Efile	<< std::fixed << t << "\t" << T << "\t" << EiC / scount << "\t"
			<< 0.5*EiD / scount << "\t" << EiI / scount << "\t"
			<< (EiC + 0.5*EiD + EiI)/scount << "\n";

	Efile.close();
	std::cout << std::setprecision(3) << "Econstr = " << EiC / scount << " Edata = " << 0.5*EiD / scount << " Eint = " << EiI / scount << "\n\n";
}

void write_par_file(float T)
{
	calc_E(T);

	par_file.open("parameters.txt", std::ofstream::app);

	par_file << "Constraint Energy: \t" << EiC / scount << "\n";
	par_file << "Data Energy: \t\t" << 0.5*EiD / scount << "\n";
	par_file << "Interaction Energy: \t" << EiI / scount << "\n";
	par_file << "Total Energy: \t\t" << (EiC + 0.5*EiD + EiI) / scount << "\n\n";

	float time_all = (float)(stop_all.int64 - start_all.int64) / CPU_freq / 60;
	par_file << "Elapsed Time: \t\t" << std::setprecision(2) << (int)time_all / 60 << " h : " << fmod(time_all, 60) << " min\n\n";
	par_file << "Comments:";
	par_file.close();
}


void plot_E()
{
	Gnuplot g;
	g.plot_E();
}
