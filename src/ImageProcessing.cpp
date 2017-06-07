
#include "ImageProcessing.h"
#include "ImageInterpolation.h"
#include <math.h>
#include <QDebug>

int upto(int x, int y) {
	return (x + y - 1) & ~(y - 1);
}

void imageProcessingFun(const QString& progName, QImage* const outImgs, const QImage* const inImgs, const QVector<double>& params) 
{
	int X_SIZE = inImgs->width();
	int Y_SIZE = inImgs->height();

	int X_SIZE_NEW, Y_SIZE_NEW;
	
	if(progName == "Sample and hold") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */

		X_SIZE_NEW = upto(X_SIZE * params[1], 4);
		Y_SIZE_NEW = upto(X_SIZE * params[0], 4);

		new (outImgs) QImage(X_SIZE_NEW, Y_SIZE_NEW, inImgs->format());
		
		sampleAndHold(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), X_SIZE_NEW, Y_SIZE_NEW);
	}
	else if (progName == "Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */

		X_SIZE_NEW = upto(X_SIZE * params[1], 4);
		Y_SIZE_NEW = upto(X_SIZE * params[0], 4);

		new (outImgs) QImage(X_SIZE_NEW, Y_SIZE_NEW, inImgs->format());

		bilinearInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), X_SIZE_NEW, Y_SIZE_NEW);
	}
	else if (progName == "Bicubic")
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */
		
		X_SIZE_NEW = upto(X_SIZE * params[1], 4);
		Y_SIZE_NEW = upto(X_SIZE * params[0], 4);

		new (outImgs) QImage(X_SIZE_NEW, Y_SIZE_NEW, inImgs->format());

 		bicubicInterpolate((uchar*)inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), X_SIZE_NEW, Y_SIZE_NEW);
	}
	else if(progName == "Rotation") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Rotation angle in degrees is params[0]*/
		/* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */

		new (outImgs) QImage(X_SIZE, Y_SIZE, inImgs->format());

		imageRotate((uchar*)inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), X_SIZE / 2, Y_SIZE / 2, params[0]);
	
	}
	else if (progName == "Rotation Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Rotation angle in degrees is params[0]*/
		/* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */

		new (outImgs) QImage(X_SIZE, Y_SIZE, inImgs->format());

		imageRotateBilinear((uchar*)inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), X_SIZE / 2, Y_SIZE / 2, (double)params[0]);
	}

}

