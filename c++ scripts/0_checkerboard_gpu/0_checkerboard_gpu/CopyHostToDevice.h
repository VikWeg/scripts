void cpyHostToDevice()
{
	cudaMalloc((void**)&dev_x, scount*sizeof(float));
	cudaMalloc((void**)&dev_y, scount*sizeof(float));
	cudaMalloc((void**)&dev_z, scount*sizeof(float));
	cudaMalloc((void**)&dev_c, scount*sizeof(unsigned long long));

	cudaMalloc((void**)&dev_sig, VoxCount*sizeof(int));
	cudaMalloc((void**)&dev_ten, VoxCount*sizeof(float));
	cudaMalloc((void**)&dev_VoxIds, VoxCount*sizeof(int));

	cudaMemcpy(dev_x, x, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y, y, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_z, z, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_c, c, scount*sizeof(unsigned long long), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_sig, sig, VoxCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_ten, ten, VoxCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_VoxIds, VoxIds, VoxCount*sizeof(int), cudaMemcpyHostToDevice);

	//=======================
}