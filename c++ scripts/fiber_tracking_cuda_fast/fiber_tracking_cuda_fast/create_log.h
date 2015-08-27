void create_log()
{
	time_t rawtime;
	struct tm * timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, ("C:\\ETH\\Neuro\\GlobalTracking\\results\\%Y_%m_%d_%H_%M" + comment).c_str(), timeinfo);

	mkdir(buffer);
	chdir(buffer);

	mkdir("fibers");

	par_file.open("parameters.txt");

	par_file << "Subject: " << "\t\t" << subject.c_str() << "\n";
	par_file << "Cube center: " << "\t\t(" << vox_origin[1] << "," << vox_origin[2] << "," << vox_origin[3] << ")\n";
	par_file << "Cube size: " << "\t\t" << cube_size[0] << "x" << cube_size[1] << "x" << cube_size[2] << "\n\n";

	par_file << "Anisotropy cutoff: " << "\t" << cutoff << "\n\n";

	par_file << "Initial Temperature: " << "\t" << Ti << "\n";
	par_file << "Stop Temperature: " << "\t" << Tf << "\n";
	par_file << "Cooling rate: " << "\t\t" << etha << "\n";
	par_file << "Number T steps: " << "\t" << tsteps_tot << "\n\n";

	par_file << "Number sweeps: " << "\t\t" << sweeps << "\n";
	par_file << "Number x-proposals: " << "\t" << nx << "\n";
	par_file << "Range x-proposals: " << "\t" << delta_x << "\n\n";

	par_file << "Constraint weight: " << "\t" << wc_str << "\n";
	par_file << "Position weight: " << "\t" << wx_str << "\n\n";

	par_file << "Data function: " << "\t\t" << wdata_str << "\n";
	par_file << "Angle function: " << "\t" << wint_str << "\n\n";

	par_file << "Number WM voxels: " << "\t" << w_vox_num << "\n";
	par_file << "Number surface voxels: " << "\t" << surf_vox_num << "\n";
	par_file << "Number vertices: " << "\t" << scount << "\n\n";

	par_file.close();
}