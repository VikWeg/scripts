void export_fibers_1(vertex* prev, vertex* curr, int n)
{
	int next = -1;
	for (int i = 0; i < (*curr).nn; i++)
	if ((*curr).c[i] == 1 && (*curr).n[i] != prev)
	{
		next = i;
		break;
	}

	fiber_file_1 << ",{";
	fiber_file_1 << curr->x + 1 << ",";
	fiber_file_1 << curr->y + 1 << ",";
	fiber_file_1 << curr->z + 1;
	fiber_file_1 << "}";

	if (next != -1 && n < 50)
		if( (*(*curr).n[next]).sig != 1)
			export_fibers_1(curr, (*curr).n[next], n + 1);
		else
		{
			fiber_file_1 << ",{";
			fiber_file_1 << (*(*curr).n[next]).x + 1 << ",";
			fiber_file_1 << (*(*curr).n[next]).y + 1 << ",";
			fiber_file_1 << (*(*curr).n[next]).z + 1;
			fiber_file_1 << "}";
		}
}

void export_fibers_1(vertex* v)
{
	int next = -1;
	for (int i = 0; i < v->nn; i++)
		if (v->c[i] == 1)
		{
			next = i;
			break;
		}

	if (next != -1)
	{
		fiber_file_1 << "{";
		fiber_file_1 << v->x + 1 << ",";
		fiber_file_1 << v->y + 1 << ",";
		fiber_file_1 << v->z + 1;
		fiber_file_1 << "}";
		export_fibers_1(v, (*v).n[next], 0);
	}
}

void export_fibers_1(int n)
{
	chdir("fibers");
	fiber_file_1.open("fibers_" + std::to_string(n) + ".txt");

	fiber_file_1 << "{";

	for (int i = 0; i < size[0]; i++)
	for (int j = 0; j < size[1]; j++)
	for (int k = 0; k < size[2]; k++)
	for (int s = 0; s < snum[offset(i, j, k)]; s++)
		if (ensemble[offset(i, j, k) + s].sig)
		{
			fiber_file_1 << "{";
			export_fibers_1(&ensemble[offset(i, j, k) + s]);
			fiber_file_1 << "},";
		}

		long pos = fiber_file_1.tellp();
		fiber_file_1.seekp(pos - 1);
		fiber_file_1 << "}";
		fiber_file_1.close();
}