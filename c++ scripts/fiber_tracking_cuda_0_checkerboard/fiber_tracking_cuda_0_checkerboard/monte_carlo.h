void mc_c(vertex* vi, float T)
{
	int d;
	float E0, E1, p;

		for (int d = 0; d < vi->nn; d++)
		{
			E0 = wc(T)*Ei_c(vi) + wx(T)*Ei_x(vi);

			vi->c[d] = 1 - vi->c[d];
			vi->cc += -1 + 2 * vi->c[d];

			int j = 0;
			for (; j < vi->n[d]->nn; j++)
				if (vi == vi->n[d]->n[j])
				{
					vi->n[d]->c[j] = 1 - vi->n[d]->c[j];
					break;
				}
			vi->n[d]->cc += -1 + 2 * vi->c[d];

			E1 = wc(T)*Ei_c(vi) + wx(T)*Ei_x(vi);

			p = fminf(1., expf((E0 - E1) / T));
			std::bernoulli_distribution acceptQ(p);

			if (acceptQ(generate))
				;
			else
			{
				vi->cc -= -1 + 2 * vi->c[d];
				vi->n[d]->cc -= -1 + 2 * vi->c[d];

				vi->c[d] = 1 - vi->c[d];
				vi->n[d]->c[j] = 1 - vi->n[d]->c[j];
			}
		}
}

void mc_x(vertex* vi, float T)
{
	float x0, y0, z0, x1, y1, z1;
	float E0, E1, p;

		for (int i = 0; i < nx; i++)
		{
			E0 = wx(T)*Ei_x(vi);

			x0 = vi->x;
			y0 = vi->y;
			z0 = vi->z;

			std::uniform_real_distribution<float> u_x(x0 - (x0 - vi->pos_x + hdr.pixdim[1] * 0.5) * delta_x, x0 + (vi->pos_x + hdr.pixdim[1] * 0.5 - x0) * delta_x);
			vi->x = u_x(generate);
			std::uniform_real_distribution<float> u_y(y0 - (y0 - vi->pos_y + hdr.pixdim[2] * 0.5) * delta_x, y0 + (vi->pos_y + hdr.pixdim[2] * 0.5 - y0) * delta_x);
			vi->y = u_y(generate);
			std::uniform_real_distribution<float> u_z(z0 - (z0 - vi->pos_z + hdr.pixdim[3] * 0.5) * delta_x, z0 + (vi->pos_z + hdr.pixdim[3] * 0.5 - z0) * delta_x);
			vi->z = u_z(generate);

			E1 = wx(T)*Ei_x(vi);

			p = fminf(1., expf((E0 - E1) / T));
			std::bernoulli_distribution acceptQ(p);

			if (acceptQ(generate))
				;
			else
			{
				vi->x = x0;
				vi->y = y0;
				vi->z = z0;
			}
		}
}

void mc(float* x, float* y, float* z, double* c, float* ten, char* sig, int* vox_id, int lattice_id)
{
	for (int vox = 0; vox <= cube_size[0] * cube_size[1] * cube_size[2]; vox++)


	if (vox_id[vox] > 0)
	{
		int next = vox + 1;
		while (vox_id[next] < 0) next++;

		for (char spin_id = vox_id[vox]; spin_id < vox_id[next]; spin_id++)
		{
			mc_c(x, y, z, c, ten, sig, spin_id, vox_id);
			mc_x(x, y, z, c, ten, sig, spin_id, vox_id);
		}
	}
}