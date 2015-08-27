void dev_init()
{
	// =========  Malloc ===========

	//read & write
	cudaMalloc(&dev_in_x,  scount*sizeof(float));
	cudaMalloc(&dev_out_x, scount*sizeof(float));

	cudaMalloc(&dev_in_y, scount*sizeof(float));
	cudaMalloc(&dev_out_y, scount*sizeof(float));

	cudaMalloc(&dev_in_z, scount*sizeof(float));
	cudaMalloc(&dev_out_z, scount*sizeof(float));

	cudaMalloc(&dev_in_cc, scount*sizeof(int));
	cudaMalloc(&dev_out_cc, scount*sizeof(int));

	cudaMalloc(&dev_in_c, n_num*sizeof(int));
	cudaMalloc(&dev_out_c, n_num*sizeof(int));


	// read-only
	cudaMalloc(&dev_pos_x, scount*sizeof(float));
	cudaMalloc(&dev_pos_y, scount*sizeof(float));
	cudaMalloc(&dev_pos_z, scount*sizeof(float));

	cudaMalloc(&dev_T0, scount*sizeof(float));
	cudaMalloc(&dev_T1, scount*sizeof(float));
	cudaMalloc(&dev_T2, scount*sizeof(float));
	cudaMalloc(&dev_T3, scount*sizeof(float));
	cudaMalloc(&dev_T4, scount*sizeof(float));
	cudaMalloc(&dev_T5, scount*sizeof(float));

	cudaMalloc(&dev_Emin, scount*sizeof(float));
	cudaMalloc(&dev_Emax, scount*sizeof(float));
	cudaMalloc(&dev_delta_E, scount*sizeof(float));

	cudaMalloc(&dev_nc, scount*sizeof(int));
	cudaMalloc(&dev_sig, scount*sizeof(int));

	cudaMalloc(&dev_n_id, scount*sizeof(int));
	cudaMalloc(&dev_n, n_num*sizeof(int));

	// =======   MemCpy ==============

	//read & write
	cudaMemcpy(dev_in_x, x, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_x, x, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_in_y, y, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_y, y, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_in_z, z, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_z, z, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_in_cc, cc, scount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_cc, cc, scount*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_in_c, c, n_num*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_c, c, n_num*sizeof(int), cudaMemcpyHostToDevice);

	// read-only
	cudaMemcpy(dev_pos_x, pos_x, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pos_y, pos_y, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pos_y, pos_z, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&dev_T0, scount*sizeof(float));
	cudaMalloc(&dev_T1, scount*sizeof(float));
	cudaMalloc(&dev_T2, scount*sizeof(float));
	cudaMalloc(&dev_T3, scount*sizeof(float));
	cudaMalloc(&dev_T4, scount*sizeof(float));
	cudaMalloc(&dev_T5, scount*sizeof(float));

	cudaMemcpy(dev_T0, T0, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T1, T1, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T2, T2, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T3, T3, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T4, T4, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T5, T5, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_Emin, Emin, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Emax, Emin, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_delta_E, Emin, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_nc, nc, scount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sig, sig, scount*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_n_id, n_id, scount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_n, n, n_num*sizeof(int), cudaMemcpyHostToDevice);
}