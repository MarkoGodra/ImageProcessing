
#ifndef IMAGEINTERPOLATION_H_
#define IMAGEINTERPOLATION_H_

#include <QString>
#include <QVector>
#include <QImage>
#include "ImageFilter.h"

void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize);

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize);

void bicubicInterpolate(uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize);

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle);

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle);

uchar cubicInterpolation(uchar points[4], double d);

char cubicInterpolation_char(char points[4], double d);

double weight(double d);

#endif // IMAGEINTERPOLATION_H_
