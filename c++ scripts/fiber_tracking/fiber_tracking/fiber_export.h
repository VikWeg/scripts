void export_fibers(vertex* prev, vertex* curr,int n)
{
	int next = -1;
	for (int i = 0; i < (*curr).nn; i++)
	if ((*curr).c[i] == 1 && (*curr).n[i] != prev)
	{
		next = i;
		break;
	}

	fiber_file << ",{";
	fiber_file << (*curr).x[0] + 1 << ",";
	fiber_file << (*curr).x[1] + 1 << ",";
	fiber_file << (*curr).x[2] + 1;
	fiber_file << "}";

	if (next != -1 && n < 50)
		if( (*(*curr).n[next]).sig != 1)
			export_fibers(curr, (*curr).n[next],n+1);
		else
		{
			fiber_file << ",{";
			fiber_file << (*(*curr).n[next]).x[0] + 1 << ",";
			fiber_file << (*(*curr).n[next]).x[1] + 1 << ",";
			fiber_file << (*(*curr).n[next]).x[2] + 1;
			fiber_file << "}";
		}
}

void export_fibers(vertex* v)
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
		fiber_file << "{";
		fiber_file << v->x[0] + 1 << ",";
		fiber_file << v->x[1] + 1 << ",";
		fiber_file << v->x[2] + 1;
		fiber_file << "}";
		export_fibers(v, (*v).n[next],0);
	}
}

void export_fibers()
{
	fiber_file.open("C:/ETH/Neuro/GlobalTracking/fibers_.txt");

	fiber_file << "{";

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			for (int k = 0; k < size; k++)
				for (int s = 0; s < snum[i][j][k]; s++)
					if (ensemble[i][j][k][s].sig == 1)
					{
						fiber_file << "{";
						export_fibers(&ensemble[i][j][k][s]);

						if (i == size - 1 && j == size - 1 && k == size - 1 && s == snum[size-1][size-1][size-1]-1)
							fiber_file << "}";
						else
							fiber_file << "},";
						std::cout << "i = " << i << " j = " << j << " k = " << k << " s = " << s << "\xd";
					}
	
	std::cout << "\nend of export function loop\n";

	fiber_file << "}";
	fiber_file.close();

	std::cout << "end of export function\n";
}