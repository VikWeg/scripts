void get_fiber(int prev, int curr, int nn)
{
	int next = -1;
	for (int i = n_id[curr]; i < n_id[curr] + nc[curr]; i++)
	if (c[i] == 1 && n[i] != prev)
	{
		next = n[i];
		break;
	}

	fwrite((char*)&(x[curr]), 1, 4, fiber_file);
	fwrite((char*)&(y[curr]), 1, 4, fiber_file);
	fwrite((char*)&(z[curr]), 1, 4, fiber_file);

	if (next != -1 && nn < 50)
	{
		if (sig[next] != 1)
			get_fiber(curr, next, nn + 1);
		else
		{
			fwrite((char*)&(x[next]), 1, 4, fiber_file);
			fwrite((char*)&(y[next]), 1, 4, fiber_file);
			fwrite((char*)&(z[next]), 1, 4, fiber_file);
		}
	}
}

void get_fiber(int s)
{
	int next = -1;
	for (int i = n_id[s]; i < n_id[s] + nc[s]; i++)
	if (c[i] == 1)
	{
		next = n[i];
		break;
	}

	if (next != -1 && sig[next] != 1)
	{
		fwrite((char*)&(x[s]), 1, 4, fiber_file);
		fwrite((char*)&(y[s]), 1, 4, fiber_file);
		fwrite((char*)&(z[s]), 1, 4, fiber_file);

		get_fiber(s, next, 1);
	}
	else if (next != -1 && sig[next] == 1)
	{
		fwrite((char*)&(x[s]), 1, 4, fiber_file);
		fwrite((char*)&(y[s]), 1, 4, fiber_file);
		fwrite((char*)&(z[s]), 1, 4, fiber_file);

		fwrite((char*)&(x[next]), 1, 4, fiber_file);
		fwrite((char*)&(y[next]), 1, 4, fiber_file);
		fwrite((char*)&(z[next]), 1, 4, fiber_file);
	}
}

int get_fiber_length(int prev, int curr, int nn)
{
	int next = -1;
	for (int i = n_id[curr]; i < n_id[curr] + nc[curr]; i++)
	if (c[i] == 1 && n[i] != prev)
	{
		next = n[i];
		break;
	}

	if (next != -1 && nn < 50)
	{
		if (sig[next] != 1)
			return get_fiber_length(curr, next, nn + 1);
		else
			return nn + 2;
	}
	else
		return nn + 1;
}

int get_fiber_length(int s)
{
	int next = -1;
	for (int i = n_id[s]; i < n_id[s] + nc[s]; i++)
	if (c[i] == 1)
	{
		next = n[i];
		break;
	}

	if (next != -1 && sig[next] != 1)
		return get_fiber_length(s, next, 1);
	else if (next != -1 && sig[next] == 1)
		return 2;
	else
		return 1;
}