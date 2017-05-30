#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>


void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
	uchar *y_old = new uchar[xSize * ySize]();
	uchar *y_new = new uchar[newXSize * newYSize]();

	char *v_old = new char[xSize * ySize / 4]();
	char *u_old = new char[xSize * ySize / 4]();
	char *v_new = new char[newXSize * newYSize / 4]();
	char *u_new = new char[newXSize * newYSize / 4]();

	/*
		Convert to YUV from RGB
	*/
	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);
	
	/*
		Calculate scale factors so i can apply formula
	*/
	const double scale_horizontal = (double)newXSize / xSize;
	const double scale_vertical = (double)newYSize / ySize;

	/*
		Walking through imput image and applying formula
		NewImage(p,q) = Image([(p - 1) / SCALE + 1], [(q - 1) / SCALE + 1])
		for Y component
	*/
	for (int i = 0; i < newYSize; i++) {
		
		for (int j = 0; j < newXSize; j++) {
			
			int ii = (i - 1) / scale_vertical + 1;
			int jj = (j - 1) / scale_horizontal + 1;
			y_new[i * newXSize + j] = y_old[ii * xSize + jj];
		
		}
	}

	/*
		Same thing for U and V components
	*/
	for (int i = 0; i < newYSize / 2; i++) {

		for (int j = 0; j < newXSize / 2; j++) {

			int ii = (i - 1) / scale_vertical + 1;
			int jj = (j - 1) / scale_horizontal + 1;
			
			u_new[i * newXSize / 2 + j] = u_old[ii * xSize / 2 + jj];
			v_new[i * newXSize / 2 + j] = v_old[ii * xSize / 2 + jj];
		}
	}

	/*
		Go back to RGB color space
	*/
	YUV420toRGB(y_new, u_new, v_new, newXSize, newYSize, output);

	/*
		Free allocated resources
	*/
	delete[] y_old;
	delete[] y_new;
	delete[] u_old;
	delete[] u_new;
	delete[] v_old;
	delete[] v_new;
}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
	uchar *y_old = new uchar[xSize * ySize]();
	uchar *y_new = new uchar[newXSize * newYSize]();

	char *v_old = new char[xSize * ySize / 4]();
	char *u_old = new char[xSize * ySize / 4]();
	char *v_new = new char[newXSize * newYSize / 4]();
	char *u_new = new char[newXSize * newYSize / 4]();

	/*
		Convert to YUV from RGB
	*/
	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);

	/*
		Calculate scale factors so i can apply formula
	*/
	const double scale_horizontal = (double)newXSize / xSize;
	const double scale_vertical = (double)newYSize / ySize;

	for (int i = 0; i < newYSize; i++) {

		for (int j = 0; j < newXSize; j++) {
			
			int a = i / scale_vertical - floor(i / scale_vertical);
			int b = j / scale_horizontal - floor(j / scale_horizontal);

			int ii = i / scale_vertical;
			int jj = j / scale_horizontal;

			y_new[i * newXSize + j] =
				(1 - a) * (1 - b) * y_old[ii * xSize + jj] +
				(1 - a) * b * y_old[(ii + 1) * xSize + jj] +
				a * (1 - b) * y_old[ii * xSize + (jj + 1)] +
				a * b * y_old[(ii + 1) * xSize + (jj + 1)];

		}
	}

	for (int i = 0; i < newYSize / 2; i++) {

		for (int j = 0; j < newXSize / 2; j++) {

			int a = i / scale_vertical - floor(i / scale_vertical);
			int b = j / scale_horizontal - floor(j / scale_horizontal);

			int ii = i / scale_vertical;
			int jj = j / scale_horizontal;

			u_new[i * newXSize / 2 + j] =
				(1 - a) * (1 - b) * u_old[ii * xSize / 2 + jj] +
				(1 - a) * b * u_old[(ii + 1) * xSize / 2 + jj] +
				a * (1 - b) * u_old[ii * xSize / 2 + (jj + 1)] +
				a * b * u_old[(ii + 1) * xSize / 2 + (jj + 1)];

			v_new[i * newXSize / 2 + j] =
				(1 - a) * (1 - b) * v_old[ii * xSize / 2 + jj] +
				(1 - a) * b * v_old[(ii + 1) * xSize / 2 + jj] +
				a * (1 - b) * v_old[ii * xSize / 2 + (jj + 1)] +
				a * b * v_old[(ii + 1) * xSize / 2 + (jj + 1)];
		}
	}

	/*
		Go back to RGB color space
	*/
	YUV420toRGB(y_new, u_new, v_new, newXSize, newYSize, output);

	/*
		Free allocated resources
	*/
	delete[] y_old;
	delete[] y_new;
	delete[] u_old;
	delete[] u_new;
	delete[] v_old;
	delete[] v_new;
}

void bicubicInterpolate(uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}