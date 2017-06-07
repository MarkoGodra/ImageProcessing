#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>


void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
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

			int ii;
			int jj;

			ii = (i - 1) / scale_vertical;
			jj = (j - 1) / scale_horizontal;

			if (ii < ySize - 1)
				ii = (i - 1) / scale_vertical + 1;

			if (jj < xSize - 1)
				jj = (j - 1) / scale_horizontal + 1;
			
			y_new[i * newXSize + j] = y_old[ii * xSize + jj];
		
		}
	}

	/*
		Same thing for U and V components
	*/
	for (int i = 0; i < newYSize / 2; i++) {

		for (int j = 0; j < newXSize / 2; j++) {
			
			int ii;
			int jj;

			ii = (i - 1) / scale_vertical;
			jj = (j - 1) / scale_horizontal;
			
			if (ii < ySize / 2 - 1)
				ii = (i - 1) / scale_vertical + 1;

			if (jj < xSize / 2 - 1)
				jj = (j - 1) / scale_horizontal + 1;
			
			
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
			
			double a = i / scale_vertical - floor(i / scale_vertical);
			double b = j / scale_horizontal - floor(j / scale_horizontal);

			int ii = i / scale_vertical;
			int jj = j / scale_horizontal;

			int iii = ii;
			int jjj = jj;

			/*
				Bottom tearing protection
			*/
			if (ii < ySize - 1)
				iii = ii + 1;
			if (jj < xSize - 1)
				jjj = jj + 1;

			y_new[i * newXSize + j] =
				(1 - a) * (1 - b) * y_old[ii * xSize + jj] +
				(1 - a) * b * y_old[ii * xSize + jjj] +
				a * (1 - b) * y_old[iii * xSize + jj] +
				a * b * y_old[iii * xSize + jjj];

		}
	}

	for (int i = 0; i < newYSize / 2; i++) {

		for (int j = 0; j < newXSize / 2; j++) {

			double a = i / scale_vertical - floor(i / scale_vertical);
			double b = j / scale_horizontal - floor(j / scale_horizontal);

			int ii = i / scale_vertical;
			int jj = j / scale_horizontal;

			int iii = ii;
			int jjj = jj;

			/*
			Bottom tearing protection
			*/
			if (ii < ySize / 2 - 1)
				iii = ii + 1;
			if (jj < xSize / 2 - 1)
				jjj = jj + 1;

			u_new[i * newXSize / 2 + j] =
				(1 - a) * (1 - b) * u_old[ii * xSize / 2 + jj] +
				(1 - a) * b * u_old[ii * xSize / 2 + jjj] +
				a * (1 - b) * u_old[iii * xSize / 2 + jj] +
				a * b * u_old[iii * xSize / 2 + jjj];

			v_new[i * newXSize / 2 + j] =
				(1 - a) * (1 - b) * v_old[ii * xSize / 2 + jj] +
				(1 - a) * b * v_old[ii * xSize / 2 + jjj] +
				a * (1 - b) * v_old[iii * xSize / 2 + jj] +
				a * b * v_old[iii * xSize / 2 + jjj];
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
	Processing
	*/
	for (int i = 0; i < newYSize; i++) {
		for (int j = 0; j < newXSize; j++) {

			double d_vertical = i / scale_vertical - floor(i / scale_vertical);
			double d_horizontal = j / scale_horizontal - floor(j / scale_horizontal);

			int ii = i / scale_vertical;
			int jj = j / scale_horizontal;
			
			uchar temp_array[4];
			uchar temp_array_2[4];

			int l;
			for (int k = ii - 1, l = 0; k < ii + 3; k++, l++) {
				temp_array[0] = y_old[k * xSize + jj - 1];
				temp_array[1] = y_old[k * xSize + jj];
				temp_array[2] = y_old[k * xSize + jj + 1];
				temp_array[3] = y_old[k * xSize + jj + 2];

				temp_array_2[l] = cubicInterpolation(temp_array, d_horizontal);
			}

			y_new[i * newXSize + j] = cubicInterpolation(temp_array_2, d_vertical);
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

double weight(double d) {
	if (abs(d) < 1.0)
		return (3.0 / 2) * abs(d) * abs(d) * abs(d) - (5.0 / 2) * abs(d) * abs(d) + 1;
	else if (abs(d) >= 1.0 && abs(d) < 2.0)
		return (-1.0 / 2) * abs(d) * abs(d) * abs(d) + 5.0 / 2 * abs(d) * abs(d) - 4.0 * abs(d) + 2;
	else
		return 0;
}

uchar cubicInterpolation(uchar points[4], double d) {

	double *w = new double[4];
	w[0] = weight(d + 1);
	w[1] = weight(d);
	w[2] = weight(1 - d);
	w[3] = weight(2 - d);

	uchar point = 0;
	for (int i = 0; i < 4; i++) {
		point += w[i] * points[i];
	}

	delete[] w;

	return point;

}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}