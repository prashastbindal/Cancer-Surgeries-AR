#ifndef CAMERA_H
#define CAMERA_H
#include <jni.h>
#include<math.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <vector>
#include <iostream>
#include"Fiducial.h"




using namespace std;
using namespace cv;

#ifdef __cplusplus
extern "C"
#endif
{


	class Camera
	{

		public:
			Mat rvec;
			Mat tvec;
			vector<float> distorVector;
			Mat cameraMatrix;
			bool estimateCameraParameters(double Fiducial[]);
			//bool calculateProjectPoints(double Fiducial[], vector< Point2f>  & ProjectimagePoints, vector<Point3f> kinect3dCoord);


	};




#ifdef __cplusplus
}
#endif





#endif
