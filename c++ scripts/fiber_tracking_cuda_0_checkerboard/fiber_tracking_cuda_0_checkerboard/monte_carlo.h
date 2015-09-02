void mc_c(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* vox_ids, int spin_id, int vox_id, int vox_num)
{
	float E0, E1, p;

	int c0 = cube_size[0];
	int c1 = cube_size[1];

	int nvoxnums[26] = {vox_num - c0*c1 - c0 - 1,	vox_num - c0*c1 - c0,	vox_num - c0*c1 - c0 + 1,
						vox_num - c0*c1 - 1,		vox_num - c0*c1,		vox_num - c0*c1 + 1,
						vox_num - c0*c1 + c0 - 1,	vox_num - c0*c1 + c0,	vox_num - c0*c1 + c0 + 1,

						vox_num - c0 - 1,			vox_num - c0,			vox_num - c0 + 1,
						vox_num - 1,										vox_num + 1,
						vox_num + c0 - 1,			vox_num + c0,			vox_num + c0 + 1,

						vox_num + c0*c1 - c0 - 1,	vox_num + c0*c1 - c0,	vox_num + c0*c1 - c0 + 1,
						vox_num + c0*c1 - 1,		vox_num + c0*c1,		vox_num + c0*c1 + 1,
						vox_num + c0*c1 + c0 - 1,	vox_num + c0*c1 + c0,	vox_num + c0*c1 + c0 + 1};

	/* Neighbor numbering:

		0, 1, 2,
		3, 4, 5
		6, 7, 8

		9,10, 11
		12,  ,13
		14,15,16

		17,18,19
		20,21,22
		23,24,25
	
	*/
	
	int spin_off = 0;
	for (int n = 0; n < 26; n++)
	{
		int n_vox_id = vox_ids[nvoxnums[n]];

		int nsc = 0;
		if (n_vox_id > 0)
		{
			int next = nvoxnums[n] + 1;
			while (vox_ids[next] < 0) next++;

			nsc = vox_ids[next] - n_vox_id;

			for (int ns = 0; ns < nsc; ns++)
			{
				int n_spin_id = n_vox_id + ns;

				E0 =	wc(T)*Ei_c(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num)
					+	wx(T)*Ei_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num);

				//=== Change connection between spin and neighbour ===

				c[spin_id] = changeBit(spin_off + ns, c[spin_id]);
				
				int diff = nvoxnums[n] - vox_num;
				int n_nsnum = 13 - diff - (diff < 0) ? 1 : 0;

				int s_spos = spin_id - vox_id;

				int n_n0num;
				int n_n0voxid = -1;
				int m = 0;
				while (n_n0voxid < 0)
				{
					n_n0num =	nvoxnums[n]
								- c0*c1*((m<9) ? 1 : 0)
								+ c0*c1*((m>16) ? 1 : 0)
								- c0*((m + m>12?1:0)%9 < 3 ? 1 : 0)
								+ c0*((m + m>12?1:0)%9 > 5 ? 1 : 0)
								+	  (m + m>12?1:0)%3 > 1 ? 1 : 0
								-	  (m + m>12?1:0)%3 < 1 ? 1 : 0
								;
					n_n0num = (n_n0num>0) ? n_n0num : 0;
					n_n0voxid = vox_ids[n_n0num];
					m++;
				}

				c[n_spin_id] = changeBit(spin_id - n_n0voxid, c[n_spin_id]);  //bullshit spin_id - n_n0voxid... mind the gaps...

				//OCTREE???

				//=====================================================

				E1 =	wc(T)*Ei_c(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num)
					+	wx(T)*Ei_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num);

				p = fminf(1., expf((E0 - E1) / T));
				std::bernoulli_distribution acceptQ(p);

				if (acceptQ(generate))
					;
				else
				{
					c[spin_id] =	changeBit(spin_off + ns, c[spin_id]);
					c[n_spin_id] =	changeBit(spin_id - n_n0voxid, c[n_spin_id]);
				}
			}
			spin_off += nsc;
		}
	}
}

void mc_x(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* vox_ids, int spin_id, int vox_id, int vox_num)
{

	float x0, y0, z0, x1, y1, z1;
	float E0, E1, p;

		for (int i = 0; i < nx; i++)
		{
			E0 = wx(T)*Ei_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num);

			x0 = x[vox_id];
			y0 = y[vox_id];
			z0 = z[vox_id];

			int pos_x = round(x0 / hdr.pixdim[1]);
			int	pos_y = round(y0 / hdr.pixdim[2]);
			int	pos_z = round(z0 / hdr.pixdim[3]);

			std::uniform_real_distribution<float> u_x(x0 - (x0 - pos_x + hdr.pixdim[1] * 0.5) * delta_x, x0 + (pos_x + hdr.pixdim[1] * 0.5 - x0) * delta_x);
			x[vox_id] = u_x(generate);
			std::uniform_real_distribution<float> u_y(y0 - (y0 - pos_y + hdr.pixdim[2] * 0.5) * delta_x, y0 + (pos_y + hdr.pixdim[2] * 0.5 - y0) * delta_x);
			y[vox_id] = u_y(generate);
			std::uniform_real_distribution<float> u_z(z0 - (z0 - pos_z + hdr.pixdim[3] * 0.5) * delta_x, z0 + (pos_z + hdr.pixdim[3] * 0.5 - z0) * delta_x);
			z[vox_id] = u_z(generate);

			E1 = wx(T)*Ei_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id, vox_num);

			p = fminf(1., expf((E0 - E1) / T));
			std::bernoulli_distribution acceptQ(p);

			if (acceptQ(generate))
				;
			else
			{
				x[vox_id] = x0;
				y[vox_id] = y0;
				z[vox_id] = z0;
			}
		}
}

void mc(float* x, float* y, float* z, long long int* c, float* ten, char* sig, int* vox_ids, int lattice_id)
{
	for (int threadId = 0; threadId <= cube_size[0] * cube_size[1] * cube_size[2]; threadId++)
	{
		int vox_num = threadIdToVoxNum(threadId, lattice_id);

		if (vox_id[vox_num] > 0)
		{
			int next = vox_num + 1;
			while (vox_id[next] < 0) next++;

			for (char spin_id = vox_id[vox_num]; spin_id < vox_id[next]; spin_id++)
			{
				mc_c(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id[vox_num], vox_num);
				mc_x(x, y, z, c, ten, sig, vox_ids, spin_id, vox_id[vox_num], vox_num);
			}
		}
	}
}