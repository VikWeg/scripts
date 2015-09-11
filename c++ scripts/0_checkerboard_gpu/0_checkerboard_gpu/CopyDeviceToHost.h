void CopyDeviceToHost()
{
	cudaMemcpy(x, dev_x, scount*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(y, dev_y, scount*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(z, dev_z, scount*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(c, dev_c, scount*sizeof(unsigned long long), cudaMemcpyDeviceToHost);

}