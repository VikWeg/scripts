void cpyDev2Host(int toggle)
{
	if (toggle % 2)
	{


		for (int i = 0; i < scount; i++)
		{
			int nn = ensemble[i].nn;

			cudaMemcpy(ensemble[i].c, temp_in_ensemble[i].c, nn*sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(&ensemble[i].x, &temp_in_ensemble[i].x, sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(&ensemble[i].y, &temp_in_ensemble[i].y, sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(&ensemble[i].z, &temp_in_ensemble[i].z, sizeof(float), cudaMemcpyDeviceToHost);
		}
	}
	else
	{
		for (int i = 0; i < scount; i++)
		{
			int nn = ensemble[i].nn;
			cudaMemcpy(ensemble[i].c, temp_out_ensemble[i].c, nn*sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(&ensemble[i].x, &dev_out_ensemble[i].x, sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(&ensemble[i].y, &dev_out_ensemble[i].y, sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(&ensemble[i].z, &dev_out_ensemble[i].z, sizeof(float), cudaMemcpyDeviceToHost);
		}
	}
}