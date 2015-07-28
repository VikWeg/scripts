#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nifti1.h"

typedef float MY_DATATYPE;

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

/**************  MAIN  ****************************/

/*
main(int argc,char *argv[]) 
{

	if(argc != 3)
		printf( "\nUsage: %s <header file> <data file>\n",argv[0]);
	else
	{
		char *hdr_file = argv[1];
		char *data_file = argv[2];

		read_nifti_file(hdr_file, data_file);
	}
}
*/

/************ read_nifti_file ****************/

int read_nifti_file(char *hdr_file, char *data_file)
{
	nifti_1_header hdr;
	FILE *fp;
	int ret,i;
	double total;
	MY_DATATYPE *data=NULL;


	/********** open and read header */
	fp = fopen(hdr_file,"r");
	if (fp == NULL)
	{
		fprintf(stderr, "\nError opening header file %s\n",hdr_file);
		exit(1);
	}
	ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
	if (ret != 1)
	{
		fprintf(stderr, "\nError reading header file %s\n",hdr_file);
		exit(1);
	}
	fclose(fp);


	/********** print a little header information 
	fprintf(stderr, "\n%s header information:",hdr_file);
	fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
	fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
	fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
	fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
	fprintf(stderr, "\n");
	*/

	/********** open the datafile, jump to data offset */
	fp = fopen(data_file,"r");
	if (fp == NULL)
	{
		fprintf(stderr, "\nError opening data file %s\n",data_file);
		exit(1);
	}

	ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);
	if (ret != 0)
	{
		fprintf(stderr, "\nError doing fseek() to %ld in data file %s\n",(long)(hdr.vox_offset), data_file);
		exit(1);
	}


	/********** allocate buffer and read first 3D volume from data file */
	data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]);

	if (data == NULL)
	{
		fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
		exit(1);
	}

	ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4], fp);

	if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4])
	{
		fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",data_file,ret);
		exit(1);
	}
	fclose(fp);


/********** scale the data buffer  
if (hdr.scl_slope != 0) {
        for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
                data[i] = (data[i] * hdr.scl_slope) + hdr.scl_inter;
}
*/

/********** print mean of data */
/*total = 0;
for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
        total += data[i];
total /= (hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
fprintf(stderr, "\nMean of volume 1 in %s is %.3f\n",data_file,total);
*/
int x = 17;
int y = 57;
int z = 32;
int t = 8;

int coo = (t-1)*128*128*72 + (72-z)*128*128 + (128-y)*128 + (128-x);

//for(i = 0; i < hdr.dim[1]; i++)
printf("\ndata[%d] = %f\n",coo,data[coo]);

return(0);
}

