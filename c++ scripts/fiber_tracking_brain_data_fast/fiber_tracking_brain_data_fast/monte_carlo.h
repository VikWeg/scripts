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

			std::uniform_real_distribution<float> u_x(x0 - (x0 - vi->pos_x + 0.5) * delta_x, x0 + (vi->pos_x + 0.5 - x0) * delta_x);
			vi->x = u_x(generate);
			std::uniform_real_distribution<float> u_y(y0 - (y0 - vi->pos_y + 0.5) * delta_x, y0 + (vi->pos_y + 0.5 - y0) * delta_x);
			vi->y = u_y(generate);
			std::uniform_real_distribution<float> u_z(z0 - (z0 - vi->pos_z + 0.5) * delta_x, z0 + (vi->pos_z + 0.5 - z0) * delta_x);
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