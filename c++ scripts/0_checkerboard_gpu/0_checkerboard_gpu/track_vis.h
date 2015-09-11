struct trackvis_header
{
	char 		id_string[6];
	short int 	dim[3];
	float 		voxel_size[3];
	float 		origin[3];
	short int 	n_scalars;
	char 		scalar_name[10][20];
	short int 	n_properties;
	char 		property_name[10][20];
	char 		reserved[508];
	char 		voxel_order[4];
	char 		pad2[4];
	float 		image_orientation_patient[6];
	char 		pad1[2];
	unsigned char 	invert_x;
	unsigned char 	invert_y;
	unsigned char 	invert_z;
	unsigned char 	swap_xy;
	unsigned char 	swap_yz;
	unsigned char 	swap_zx;
	int 	n_count;
	int 	version;
	int 	hdr_size;
};

void write_header(FILE* fiber_file)
{
	trackvis_header trk_hdr;

	sprintf(trk_hdr.id_string, "TRACK");

	trk_hdr.dim[0] = cube_size[0];
	trk_hdr.dim[1] = cube_size[1];
	trk_hdr.dim[2] = cube_size[2];

	trk_hdr.voxel_size[0] = hdr.pixdim[1];
	trk_hdr.voxel_size[1] = hdr.pixdim[2];
	trk_hdr.voxel_size[2] = hdr.pixdim[3];

	trk_hdr.origin[0] = 0;
	trk_hdr.origin[1] = 0;
	trk_hdr.origin[2] = 0;

	trk_hdr.n_scalars = 0;
	trk_hdr.n_properties = 0;
	sprintf(trk_hdr.voxel_order, "RAS");
	sprintf(trk_hdr.pad2, "LPS");
	trk_hdr.image_orientation_patient[0] = 1.0;
	trk_hdr.image_orientation_patient[1] = 0.0;
	trk_hdr.image_orientation_patient[2] = 0.0;
	trk_hdr.image_orientation_patient[3] = 0.0;
	trk_hdr.image_orientation_patient[4] = 1.0;
	trk_hdr.image_orientation_patient[5] = 0.0;
	trk_hdr.pad1[0] = 0;
	trk_hdr.pad1[1] = 0;
	trk_hdr.invert_x = 0;
	trk_hdr.invert_y = 0;
	trk_hdr.invert_z = 0;
	trk_hdr.swap_xy = 0;
	trk_hdr.swap_yz = 0;
	trk_hdr.swap_zx = 0;
	trk_hdr.n_count = 0;
	trk_hdr.version = 1;
	trk_hdr.hdr_size = 1000;

	fwrite((char*)&trk_hdr, 1, 1000, fiber_file);
}