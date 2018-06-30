#pragma once
#include <Windows.h>
#include <Kinect.h>
#include <opencv2/opencv.hpp>
#include <NuiKinectFusionApi.h>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/core/core.hpp>
//#include <WINSOCK2.H>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iostream>

using namespace std;
using namespace cv;

extern IplImage* img1;
extern IplImage* img2;
extern IplImage* DispImage;
extern char* name;

extern IplImage *updateGui(IplImage *DispImage, IplImage *img1, IplImage *img2);

class Kincet3D
{
public:
	Kincet3D(void);
	~Kincet3D(void);



	/////////////////////////////////////add by me///////////////
	IColorFrameSource* pColorSource;
	IFrameDescription* pColorDescription;
	Mat colorBufferMat;
	Mat colorMat;
	unsigned int colorBufferSize;
	bool FindFiducials(const Mat& mat);
	bool find_glyphs(const Mat &img, Point2f *glyph_center);
	int Kincet3D::KinectRun();
	int Kincet3D::KinectInitial();
	int KinectCoordinateCalibrate();
	void Kincet3D::thetransform();
	bool Kincet3D::kinect_image( Mat& image);
	Mat rgb_dst;
	Mat gray;
	Mat gs_src;
	Mat bw_src;
	Mat edges;
	Mat src;
	Point2f fiducial[8];
	bool NewCalculate;
	Mat depthBufferMat;
	UINT16*                     m_pDepthRawPixelBuffer;
	Point3f first[8];
	bool isFinddepth;
	bool isSaved;
	bool isSave;
	Mat DataTransform;
	/////////////////////////////////


	Mat cloudPoints;
	bool saveOriginalData;
	bool isSaveCloudPoints;
	int minRow, maxRow, minCol, maxCol;
	vector <Point3f> Target3D;
	vector< char> data;
	bool dataReady, dataReserved, askData, getWrong;
	Mat depthMat;

	double glassFiducial[8];
	vector<Point> flouce_points_kinect;
	bool is_flouce_points_kinect_using;
	bool is_flouce_points_kinect_ready;
	int trackPos1;
	int trackPos2;


	

	

};

