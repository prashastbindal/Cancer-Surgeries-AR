#include"Camera.h"


using namespace std;
using namespace cv;




bool Camera::estimateCameraParameters(double Fiducial[])
{
	Point3f objectPoints[4] = { Point3i(00.00, 0.00, 0.0), Point3i( 0.0,45.00, 0.0), Point3i(25.0,0.0,  0.0), Point3i(25.0,45.0,  00.0) }; //this numbers are in centimeter
	Mat objectPointMat = Mat(4, 1, CV_32FC3, objectPoints);

	Point2f imagePoint[4] = { Point(Fiducial[0], Fiducial[1]), Point(Fiducial[2], Fiducial[3]), Point(Fiducial[4], Fiducial[5]), Point(Fiducial[6], Fiducial[7]) };
	Mat imagePointMat = Mat(4, 1, CV_32FC2, imagePoint);

	//float cameraParameter[3][3] = { { 1003.40559, 0, 511.5 }, { 0, 511.5, 287.5}, { 0, 0, 1 } };
	float cameraParameter[3][3] = { { 988, 0, 537.5 }, { 0, 986, 328}, { 0, 0, 1 } };
	//float cameraParameter[3][3] = { { 481.995, 0, 399.5 }, { 0, 481.9952, 239.5 }, { 0, 0, 1 } };

	cameraMatrix = Mat(3, 3, CV_32FC1, cameraParameter);

	// (&distC[0], &distC[5]);
	/*distorVector.push_back(0.031733);
	distorVector.push_back(0.1256);
	distorVector.push_back(0);
	distorVector.push_back(0);
	distorVector.push_back(-0.597547193);*/

	distorVector.clear();

	distorVector.push_back(0.1677);
	distorVector.push_back(-2.109);
	distorVector.push_back(0.02306);
	distorVector.push_back(0.006997);
	distorVector.push_back(10.365);






	solvePnP(objectPointMat, imagePointMat, cameraMatrix, distorVector, rvec, tvec);







	/////////// project 3D points to 2D iamge/////////////////////////
	if (rvec.empty() && tvec.empty())
	{
		return false;
	}


	return true;


}
