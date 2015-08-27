void cpyEnsembleDevToHost()
{
	if (toggle == 0)
	{
		cudaMemcpy(x, dev_in_x, scount*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(y, dev_in_y, scount*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(z, dev_in_z, scount*sizeof(float), cudaMemcpyHostToDevice);

		cudaMemcpy(cc, dev_in_cc, box_num*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(c, dev_in_c, n_num*sizeof(int), cudaMemcpyHostToDevice);
	}
	else
	{
		cudaMemcpy(x, dev_out_x, scount*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(y, dev_out_y, scount*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(z, dev_out_z, scount*sizeof(float), cudaMemcpyHostToDevice);

		cudaMemcpy(cc, dev_out_cc, box_num*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(c, dev_out_c, n_num*sizeof(int), cudaMemcpyHostToDevice);
	}
}