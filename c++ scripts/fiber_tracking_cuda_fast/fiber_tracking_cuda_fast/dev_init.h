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

	cudaMalloc(&dev_in_cc, box_num*sizeof(int));
	cudaMalloc(&dev_out_cc, box_num*sizeof(int));

	cudaMalloc(&dev_in_c, n_num*sizeof(int));
	cudaMalloc(&dev_out_c, n_num*sizeof(int));


	// read-only
	cudaMalloc(&dev_id, box_num*sizeof(int));

	cudaMalloc(&dev_pos_x, box_num*sizeof(float));
	cudaMalloc(&dev_pos_y, box_num*sizeof(float));
	cudaMalloc(&dev_pos_z, box_num*sizeof(float));

	cudaMalloc(&dev_T0, box_num*sizeof(float));
	cudaMalloc(&dev_T1, box_num*sizeof(float));
	cudaMalloc(&dev_T2, box_num*sizeof(float));
	cudaMalloc(&dev_T3, box_num*sizeof(float));
	cudaMalloc(&dev_T4, box_num*sizeof(float));
	cudaMalloc(&dev_T5, box_num*sizeof(float));

	cudaMalloc(&dev_Emin, box_num*sizeof(float));
	cudaMalloc(&dev_Emax, box_num*sizeof(float));
	cudaMalloc(&dev_delta_E, box_num*sizeof(float));

	cudaMalloc(&dev_nc, box_num*sizeof(int));
	cudaMalloc(&dev_sig, box_num*sizeof(int));

	cudaMalloc(&dev_n_id, box_num*sizeof(int));
	cudaMalloc(&dev_n, n_num*sizeof(int));

	// =======   MemCpy ==============

	//read & write
	cudaMemcpy(dev_in_x, x, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_x, x, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_in_y, y, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_y, y, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_in_z, z, scount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_z, z, scount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_in_cc, cc, box_num*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_cc, cc, box_num*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_in_c, c, n_num*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_out_c, c, n_num*sizeof(int), cudaMemcpyHostToDevice);

	// read-only

	cudaMemcpy(dev_id, id, box_num*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_pos_x, pos_x, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pos_y, pos_y, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pos_y, pos_z, box_num*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&dev_T0, box_num*sizeof(float));
	cudaMalloc(&dev_T1, box_num*sizeof(float));
	cudaMalloc(&dev_T2, box_num*sizeof(float));
	cudaMalloc(&dev_T3, box_num*sizeof(float));
	cudaMalloc(&dev_T4, box_num*sizeof(float));
	cudaMalloc(&dev_T5, box_num*sizeof(float));

	cudaMemcpy(dev_T0, T0, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T1, T1, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T2, T2, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T3, T3, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T4, T4, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_T5, T5, box_num*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_Emin, Emin, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Emax, Emin, box_num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_delta_E, Emin, box_num*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_nc, nc, box_num*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sig, sig, box_num*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_n_id, n_id, box_num*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_n, n, n_num*sizeof(int), cudaMemcpyHostToDevice);
}