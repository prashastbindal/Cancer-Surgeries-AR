#ifndef PROJECT3D_H
#define PROJECT3D_H
#include <jni.h>
#include<math.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <vector>
#include <iostream>
#include"Camera.h"
#include"Fiducial.h"



using namespace std;
using namespace cv;

#ifdef __cplusplus
extern "C"
#endif
{


	class Project3D
	{
		public:
			bool Project3Dto2D(Camera cam, vector<Point3f> kinect3dCoord,  vector< Point2f>  & ProjectimagePoints);


	};




#ifdef __cplusplus
}
#endif





#endif
