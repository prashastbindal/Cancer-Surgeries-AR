#ifndef FIDUCIAL_H
#define FIDUCIAL_H
#include <jni.h>
#include<math.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <vector>
#include <fstream>
#include <iostream>




using namespace std;
using namespace cv;

#ifdef __cplusplus
extern "C"
#endif
{


/*
struct attributesFiducial
{
	int *pattern;
};*/



	//short ProjectimagePoints[1024*576];
	//jbyte temp[230400];


class Fiducial
{

	private:
	Point2f fiducial[4];




	public:
	bool getPoints( Mat &src, double *outputfiducial);
	bool find_glyphs(const Mat &img, Point2f *glyph_center) ;

};

#ifdef __cplusplus
}
#endif





#endif
