void get_fiber(vertex* prev, vertex* curr, int n)
{
	int next = -1;
	for (int i = 0; i < (*curr).nn; i++)
	if (curr->c[i] == 1 && curr->n[i] != prev)
	{
		next = i;
		break;
	}

	fwrite((char*)&(curr->x), 1, 4, fiber_file);
	fwrite((char*)&(curr->y), 1, 4, fiber_file);
	fwrite((char*)&(curr->z), 1, 4, fiber_file);

	if (next != -1 && n < 50)
	{
		if (curr->n[next]->sig != 1)
			get_fiber(curr, curr->n[next], n + 1);
		else
		{
			fwrite((char*)&(curr->n[next]->x), 1, 4, fiber_file);
			fwrite((char*)&(curr->n[next]->y), 1, 4, fiber_file);
			fwrite((char*)&(curr->n[next]->z), 1, 4, fiber_file);
		}
	}
}

void get_fiber(vertex* v)
{
	int next = -1;
	for (int i = 0; i < v->nn; i++)
	if (v->c[i] == 1)
	{
		next = i;
		break;
	}

	if (next != -1 && v->n[next]->sig != 1)
	{
		fwrite((char*)&(v->x), 1, 4, fiber_file);
		fwrite((char*)&(v->y), 1, 4, fiber_file);
		fwrite((char*)&(v->z), 1, 4, fiber_file);

		get_fiber(v, v->n[next], 1);
	}
	else if (next != -1 && v->n[next]->sig == 1)
	{
		fwrite((char*)&(v->x), 1, 4, fiber_file);
		fwrite((char*)&(v->y), 1, 4, fiber_file);
		fwrite((char*)&(v->z), 1, 4, fiber_file);

		fwrite((char*)&(v->n[next]->x), 1, 4, fiber_file);
		fwrite((char*)&(v->n[next]->y), 1, 4, fiber_file);
		fwrite((char*)&(v->n[next]->z), 1, 4, fiber_file);
	}
}

int get_fiber_length(vertex* prev, vertex* curr, int n)
{
	int next = -1;
	for (int i = 0; i < curr->nn; i++)
	if (curr->c[i] == 1 && curr->n[i] != prev)
	{
		next = i;
		break;
	}

	if (next != -1 && n < 50)
	{
		if (curr->n[next]->sig != 1)
			return get_fiber_length(curr, curr->n[next], n + 1);
		else
			return n + 2;
	}
	else
		return n + 1;
}

int get_fiber_length(vertex* v)
{
	int next = -1;
	for (int i = 0; i < v->nn; i++)
	if (v->c[i] == 1)
	{
		next = i;
		break;
	}

	if (next != -1 && v->n[next]->sig != 1)
		return get_fiber_length(v, v->n[next], 1);
	else if (next != -1 && v->n[next]->sig == 1)
		return 2;
	else
		return 1;
}