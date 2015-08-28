void cpyEnsembleDevToHost()
{
	if (toggle == 0)
	{
		cudaMemcpy(x, dev_in_x, scount*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(y, dev_in_y, scount*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(z, dev_in_z, scount*sizeof(float), cudaMemcpyDeviceToHost);

		cudaMemcpy(cc, dev_in_cc, scount*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(c, dev_in_c, n_num*sizeof(int), cudaMemcpyDeviceToHost);
	}
	else
	{
		cudaMemcpy(x, dev_out_x, scount*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(y, dev_out_y, scount*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(z, dev_out_z, scount*sizeof(float), cudaMemcpyDeviceToHost);

		cudaMemcpy(cc, dev_out_cc, scount*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(c, dev_out_c, n_num*sizeof(int), cudaMemcpyDeviceToHost);
	}
}