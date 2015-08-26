void import_data()
{
	FILE* header_file = fopen(("C:/ETH/Neuro/GlobalTracking/subjects/" + subject + "/tensor.nii").c_str(), "r");
	if (header_file == NULL)
	{
		fprintf(stderr, "\nError opening header file\n");
		exit(1);
	}

	long ret = fread(&hdr, MIN_HEADER_SIZE, 1, header_file);
	if (ret != 1)
	{
		fprintf(stderr, "\nError reading header file %s\n", header_file);
		exit(1);
	}

	fclose(header_file);

	//=============================================================================

	//******** Open data file ********
	FILE* data_file = fopen(("C:/ETH/Neuro/GlobalTracking/subjects/" + subject + "/tensor.nii").c_str(), "rb");
	if (data_file == NULL)
	{
		fprintf(stderr, "\nError opening data file\n");
	}
	//******** Set read head ********
	ret = fseek(data_file, (long)(hdr.vox_offset), SEEK_SET);
	if (ret != 0)
	{
		fprintf(stderr, "\nError doing fseek() to data file\n");
	}
	//******** Allocate data memory ********
	data = (float *)malloc(sizeof(float)* hdr.dim[1] * hdr.dim[2] * hdr.dim[3] * hdr.dim[4]);
	if (data == NULL)
	{
		fprintf(stderr, "\nError allocating data buffer\n");
	}
	//******** Read data from file ********
	ret = fread(data, sizeof(float), hdr.dim[1] * hdr.dim[2] * hdr.dim[3] * hdr.dim[4], data_file);
	if (ret != hdr.dim[1] * hdr.dim[2] * hdr.dim[3] * hdr.dim[4])
	{
		fprintf(stderr, "\nError reading data\n");
	}

	fclose(data_file);

	//=============================================================================

	//******** Open data file ********
	FILE* L1_file = fopen(("C:/ETH/Neuro/GlobalTracking/subjects/" + subject + "/L1.nii").c_str(), "rb");
	if (L1_file == NULL)
	{
		fprintf(stderr, "\nError opening data file\n");
	}
	//******** Set read head ********
	ret = fseek(L1_file, (long)(hdr.vox_offset), SEEK_SET);
	if (ret != 0)
	{
		fprintf(stderr, "\nError doing fseek() to data file\n");
	}
	//******** Allocate data memory ********
	L1data = (float *)malloc(sizeof(float)* hdr.dim[1] * hdr.dim[2] * hdr.dim[3]);
	if (L1data == NULL)
	{
		fprintf(stderr, "\nError allocating data buffer\n");
	}
	//******** Read data from file ********
	ret = fread(L1data, sizeof(float), hdr.dim[1] * hdr.dim[2] * hdr.dim[3], L1_file);
	if (ret != hdr.dim[1] * hdr.dim[2] * hdr.dim[3])
	{
		fprintf(stderr, "\nError reading data\n");
	}

	fclose(L1_file);

	//=============================================================================

	//******** Open data file ********
	FILE* L3_file = fopen(("C:/ETH/Neuro/GlobalTracking/subjects/" + subject + "/L3.nii").c_str(), "rb");
	if (L3_file == NULL)
	{
		fprintf(stderr, "\nError opening data file\n");
	}
	//******** Set read head ********
	ret = fseek(L3_file, (long)(hdr.vox_offset), SEEK_SET);
	if (ret != 0)
	{
		fprintf(stderr, "\nError doing fseek() to data file\n");
	}
	//******** Allocate data memory ********
	L3data = (float *)malloc(sizeof(float)* hdr.dim[1] * hdr.dim[2] * hdr.dim[3]);
	if (L3data == NULL)
	{
		fprintf(stderr, "\nError allocating data buffer\n");
	}
	//******** Read data from file ********
	ret = fread(L3data, sizeof(float), hdr.dim[1] * hdr.dim[2] * hdr.dim[3], L3_file);
	if (ret != hdr.dim[1] * hdr.dim[2] * hdr.dim[3])
	{
		fprintf(stderr, "\nError reading data\n");
	}

	fclose(L3_file);

	//=============================================================================

	//******** Open data file ********
	FILE* L2_file = fopen(("C:/ETH/Neuro/GlobalTracking/subjects/" + subject + "/L2.nii").c_str(), "rb");
	if (L2_file == NULL)
	{
		fprintf(stderr, "\nError opening data file\n");
	}
	//******** Set read head ********
	ret = fseek(L2_file, (long)(hdr.vox_offset), SEEK_SET);
	if (ret != 0)
	{
		fprintf(stderr, "\nError doing fseek() to data file\n");
	}
	//******** Allocate data memory ********
	L2data = (float *)malloc(sizeof(float)* hdr.dim[1] * hdr.dim[2] * hdr.dim[3]);
	if (L2data == NULL)
	{
		fprintf(stderr, "\nError allocating data buffer\n");
	}
	//******** Read data from file ********
	ret = fread(L2data, sizeof(float), hdr.dim[1] * hdr.dim[2] * hdr.dim[3], L2_file);
	if (ret != hdr.dim[1] * hdr.dim[2] * hdr.dim[3])
	{
		fprintf(stderr, "\nError reading data\n");
	}

	fclose(L2_file);

	//=============================================================================

	//******** Open data file ********
	FILE* mask_file = fopen(("C:/ETH/Neuro/GlobalTracking/subjects/" + subject + "/wm_mask.nii").c_str(), "rb");
	if (mask_file == NULL)
	{
		fprintf(stderr, "\nError opening data file\n");
	}
	//******** Set read head ********
	ret = fseek(mask_file, (long)(hdr.vox_offset), SEEK_SET);
	if (ret != 0)
	{
		fprintf(stderr, "\nError doing fseek() to data file\n");
	}
	//******** Allocate data memory ********
	mask = (float *)malloc(sizeof(float)* hdr.dim[1] * hdr.dim[2] * hdr.dim[3]);
	if (mask == NULL)
	{
		fprintf(stderr, "\nError allocating data buffer\n");
	}
	//******** Read data from file ********
	ret = fread(mask, sizeof(float), hdr.dim[1] * hdr.dim[2] * hdr.dim[3], mask_file);
	if (ret != hdr.dim[1] * hdr.dim[2] * hdr.dim[3])
	{
		fprintf(stderr, "\nError reading data\n");
	}

	fclose(mask_file);
}