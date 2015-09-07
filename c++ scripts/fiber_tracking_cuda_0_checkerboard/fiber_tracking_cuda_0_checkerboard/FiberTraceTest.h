void FiberTraceTest()
{
	cube_size[0] = 2;
	cube_size[1] = 2;
	cube_size[2] = 2;

	scount = 8;

	delete[] VoxIds;
	delete[] c;
	delete[] sig;
	delete[] x;
	delete[] y;
	delete[] z;

	VoxIds = new int[8] {0,1,2,3,4,5,6,7};
	c = new unsigned long long[8] {1, 5, 36, 6, 16, 80, 68, 96};

	sig = new int[8] {1, 0, 0, 0, 0, 0, 0, 0};

	x = new float[8] { 0, 1, 0, 1, 0, 1, 0, 1 };
	y = new float[8] { 0, 0, 1, 1, 0, 0, 1, 1 };
	z = new float[8] { 0, 0, 0, 0, 1, 1, 1, 1 };

	export_fibers(VoxIds,0);
}