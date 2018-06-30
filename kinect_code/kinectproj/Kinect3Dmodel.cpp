#include <Windows.h>
#include <Kinect.h>
#include <opencv2/opencv.hpp>
#include <NuiKinectFusionApi.h>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/core/core.hpp>
#include <fstream>
//#include <WINSOCK2.H>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include"projectPoints.h"
#include "Kincet3D.h"
#include "CameraControl.h"
#pragma comment(lib,"ws2_32.lib")


using namespace std;
using namespace cv;

//template<class Interface>
//inline void SafeRelease( Interface *& pInterfaceToRelease )
//{
//	if( pInterfaceToRelease != NULL ){
//		pInterfaceToRelease->Release();
//		pInterfaceToRelease = NULL;
//	}
//}


 IplImage* img1;
 IplImage* img2;
 IplImage* DispImage;
char* name = "Demo Window";

IplImage *updateGui(IplImage *DispImage, IplImage *img1, IplImage *img2);


bool calculateProjectPoints(double inputdata[], vector<Point2f> & ProjectimagePoints);
bool FindFloucePoints(const Mat& inputmat, vector<Point>&  flouce_poins);
void RegisterImages();
bool FindFiducials(const Mat& inputmat, Point2f* fiducial);
bool find_glyphs(const Mat &img, Point2f *glyph_center);
void GetFlouce();
void TestWithFlouce();
void mouseHandler(int event, int x, int y, int flags, void *param);
bool GetBoxByMouse(Mat& inputimage);
Kincet3D mKinect;

IplImage*updateGui(IplImage *DispImage, IplImage *img1, IplImage *img2)
{
	
	
	int size;
	int i;
	int m, n;
	int x, y;

	// w - Maximum number of images in a row 
	// h - Maximum number of images in a column 
	int w, h;

	// scale - How much we have to resize the image
	float scale;
	int max;
	w = 2; h = 1;
	size = 600;




	m = 80;
	n = 160;
	// Get the Pointer to the IplImage
	

	// Find the width and height of the image
	x = img1->width;
	y = img1->height;

	// Find whether height or width is greater in order to resize the image
	max = (x > y) ? x : y;

	// Find the scaling factor to resize the image
	scale = (float)((float)max / size);
	// Used to Align the images




	// Set the image ROI to display the current image
	cvSetImageROI(DispImage, cvRect(m, n, (int)(x / scale), (int)(y / scale)));

	// Resize the input image and copy the it to the Single Big Image
	cvResize(img1, DispImage);
	// Reset the ROI in order to display the next image
	cvResetImageROI(DispImage);




	
	m += 40 + size;

	// Find the width and height of the image
	x = img2->width;
	y = img2->height;

	// Find whether height or width is greater in order to resize the image
	max = (x > y) ? x : y;

	// Find the scaling factor to resize the image
	scale = (float)((float)max / size);
	// Used to Align the images




	// Set the image ROI to display the current image
	cvSetImageROI(DispImage, cvRect(m, n, (int)(x / scale), (int)(y / scale)));

	// Resize the input image and copy the it to the Single Big Image
	cvResize(img2, DispImage);
	// Reset the ROI in order to display the next image
	cvResetImageROI(DispImage);		
	return DispImage;
}

int colorInt = 0;
int colorInt2 = 0;



void switch_callback(int position){
	colorInt = position;
}

void switch_callback2(int position){
	colorInt2 = position;
}

void SetIdentityMatrix(Matrix4 &mat)
{
	mat.M11 = 1; mat.M12 = 0; mat.M13 = 0; mat.M14 = 0;
	mat.M21 = 0; mat.M22 = 1; mat.M23 = 0; mat.M24 = 0;
	mat.M31 = 0; mat.M32 = 0; mat.M33 = 1; mat.M34 = 0;
	mat.M41 = 0; mat.M42 = 0; mat.M43 = 0; mat.M44 = 1;
}
void rigidTransform();
int sendData();
void test();
void	calculate();
int KinectRun();
void affine3D(cv::Mat& DataTransform);
void thetransform();

int mLoadData(string fileName, cv::Mat& matData, int matRows, int matCols, int matChns);
bool* dataReady;

DWORD  WINAPI  AnswerThreadSendModel(LPVOID  lparam);
DWORD  WINAPI  KinectThread(LPVOID  lparam)
{
	//	mKinect.KinectCoordinateCalibrate();
	mKinect.KinectRun();
	return 0;
};
char  sendbuf[3473408];
vector<short> sendPoints;
BYTE imageDataMatrix[576][1024];
BYTE tempDataMatrix[576][1024];

/////for get points using mouse
bool is_need_manual_camera_points = true;
bool is_get_fiducials = false;
Vector<Point> ficucial_centers;

///// for the threshould of the camera with flouce 
int light_threshold = 150;




int g_switch_value = 0;
int g_switch_value2 = 0;




int main()
{
	//mKinect.KinectCoordinateCalibrate();
	int radius = 30;
	int thickness = 2;
	int connectivity = 8;

	img1 = cvLoadImage("abc.jpg");
	img2 = cvLoadImage("abc.jpg");
	IplImage* src1 = cvLoadImage("abc.jpg");
	IplImage *img;
	
	int size;
	int i;
	int m, n;
	int x, y;

	// w - Maximum number of images in a row 
	// h - Maximum number of images in a column 
	int w, h;

	// scale - How much we have to resize the image
	float scale;
	int max;
	w = 2; h = 1;
	size = 600;

	DispImage = cvCreateImage(cvSize(200 + size*w, 120 + size*h), 8, 3);
	DispImage = updateGui(DispImage, img1, img2);


	// Create a new window, and show the Single Big Image
	cvNamedWindow(name, 1);




	// Create trackbar
	cvCreateTrackbar("Threshold", name, &g_switch_value, 255, switch_callback);
	cvCreateTrackbar("Calibrate", name, &g_switch_value, 1, switch_callback);
	cvCreateTrackbar("Control", name, &g_switch_value2, 2, switch_callback2); 
	
	Mat d(DispImage);
	cv::putText(d, "Control:", Point(20, 20), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0));
	cv::putText(d,"1 for Adjusting Threshold", Point(20, 40), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0.));
	cv::putText(d, "2 for Transfer of Points to Glass ", Point(20, 60), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0.));
	
	cv::putText(d, "Captured Image", Point(160, 140), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0.));
	cv::putText(d, "Flourescent Points", Point(800, 140), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0.));
	//Mat newImage(mKinect.DispImage);
	// Loop to update the circle color
	//while (1) {
		imshow(name, d);
		cvWaitKey(5000);
		
		
		cvSetTrackbarPos("Threshold", name, 50);
		light_threshold = colorInt;
		
		
			
		if (colorInt2==1) 
		{
			
			TestWithFlouce();
			
			

			////	KinectRun();
			////	Point2f fiducial[8];
			////int test;
			////	affine3D();
			////	sendData();
			//test();
				//HANDLE handle = CreateThread(NULL, 0, KinectThread, NULL, 0, NULL);
			//if (handle == NULL)
			//{
			//	printf("CreatThread  KinectThread  failed.\n");
			//}
			//else
			//{
			//	printf("CreateThread KinectThread OK.\n");
			//}

			////	KinectRun();
			////	Point2f fiducial[8];
			////int test;
			////	affine3D();
			////	sendData();
			
			//break;
		}
		



	cvReleaseImage(&src1);

	
	//cvDestroyWindow(mKinect.name);
	cvReleaseImage(&DispImage);
	return 0;
}






DWORD  WINAPI  AnswerThreadSendModel(LPVOID  lparam)
{



	SOCKET  ClientSocket = (SOCKET)(LPVOID)lparam;
	double receviedFiducial[8];
	vector<float> send3dPoints;

	int  bytesRecv;


	char  recvbuf[100];
	char  sendfirst[10];
	//int dataPoint=0;
	//char str_android[200];
	//int i=0;
	while (1)
	{
		bytesRecv = SOCKET_ERROR;
		// first time get the inform from glass
		//Receiving Data
		bytesRecv = recv(ClientSocket, recvbuf, 100, 0);
		if (bytesRecv == -1)
		{
			break;
		}
		else
		{
			// the first data from glass should its 4 fiducials. 
			if (bytesRecv == 64)
			{
				//using the flouce camera to get target area
				//TestWithFlouce();
				//GetFlouce();
				////
				printf("Fidcials received\n");
				memcpy(&mKinect.glassFiducial[0], recvbuf, bytesRecv);
				mKinect.askData = true;//after this, another thread will to get 3D surface data
				mKinect.dataReady = false;

				//		printf("%s\n",recvbuf);
				vector< Point2f> ProjectPoints;
				//bool getData = calculateProjectPoints(mKinect.glassFiducial, ProjectPoints);// in here, the function will calculate the glass camera position based on the glass fiducials, then wait the 3D surface data ready, then project them to 2d image with glass position.
				//bool getData = true;
				//sendPoints.push_back(30);
				//sendPoints.push_back(40);
				//sendPoints.push_back(50);
				//sendPoints.push_back(60);
				bool getData = true;

				if (mKinect.getWrong)
					getData = false;

				while (1)
					if (mKinect.dataReady && !mKinect.getWrong)
						break;


				if (getData)
				{
					send3dPoints.clear();
					for (int k = 0; k < mKinect.Target3D.size(); k++)
					{
						send3dPoints.push_back(mKinect.Target3D.at(k).x);
						send3dPoints.push_back(mKinect.Target3D.at(k).y);
						send3dPoints.push_back(mKinect.Target3D.at(k).z);
					}


					//long bb = sendPoints.size() * 2;// get the total length in bytes
					long bb = send3dPoints.size() * sizeof(float);
					memcpy(sendfirst, &bb, 4);
					bytesRecv = send(ClientSocket, sendfirst, 4, 0);// send the length of the total data to glass, here we are checking whether the google glass is ready to accept data because the size of the data to be sent after this is too much





					if (bytesRecv == SOCKET_ERROR)
					{
						printf("send error");
						break;
					}
					//while(bytesRecv==SOCKET_ERROR)  
					//{  //Receiving Data





					bytesRecv = recv(ClientSocket, recvbuf, 100, 0);// wait the glass answer, the answer length is 2 means not ready, other 4 means ready.



					//printf("bytes received=%d\n", bytesRecv);
					if (bytesRecv == 2)
						continue;
					if (bytesRecv == 4)
					{

						//} 




						memcpy(sendbuf, &send3dPoints[0], size_t(bb));


						cout << send3dPoints.at(0) << " and " << send3dPoints.at(1) << " and " << send3dPoints.at(2) << "\n";
						cout << send3dPoints.size() << "\n";

						//	short testPoints[640 * 360];
						//	memcpy(&testPoints[0], sendbuf, size_t(bb));



						bytesRecv = send(ClientSocket, sendbuf, bb, 0);// send the project 2d image data( only need to be changed points)




						//printf("send data success\n");
						//			bytesRecv=send(ClientSocket,sendbuf,100,0); 
						if (bytesRecv == SOCKET_ERROR)
						{
							printf("send error");
							break;
						}
					}


					// printf("%s\n",sendbuf);

				}
				else
				{
					bytesRecv = send(ClientSocket, "wrong", 5, 0);
					bytesRecv = recv(ClientSocket, recvbuf, 100, 0);
				}
				if (bytesRecv == SOCKET_ERROR)
				{
					printf("send error");
					break;
				}
			}
		}
		//char sendBuf[50];  
		//       sendBuf[0]='w';//elcome';
		//		memcpy(sendbuf,&inData[0],size_t(dataLenth));
		//        send(ClientSocket,sendbuf,dataLenth+1,0); 

	}
	printf("thread error and end");
	return  0;


}
void test()
{
	//creat socket
	WORD myVersionRequest;
	WSADATA wsaData;
	myVersionRequest = MAKEWORD(1, 1);
	int err;
	err = WSAStartup(myVersionRequest, &wsaData);
	if (!err)
	{
		printf("Socket opened\n");
	}
	else
	{
		//
		printf("Socket unopened!");
		return;
	}
	SOCKET serSocket = socket(AF_INET, SOCK_STREAM, 0);//creat socket
	if (INVALID_SET_FILE_POINTER == serSocket)
		printf("Socket unopened!");

	//bind parameters
	SOCKADDR_IN addr;
	addr.sin_family = AF_INET;
	addr.sin_addr.S_un.S_addr = htonl(INADDR_ANY);//ip
	addr.sin_port = htons(6000);//port

	int bindid;
	::bind(serSocket, (SOCKADDR*)&addr, sizeof(SOCKADDR));
	//	printf("bind erro");
	err = listen(serSocket, 2);//

	//////////////////////////////////////////////////////////////////////////
	//beging listenning
	//////////////////////////////////////////////////////////////////////////
	SOCKADDR_IN clientsocket;
	int len = sizeof(SOCKADDR);
	while (1)
	{
		SOCKET clientConn = accept(serSocket, (SOCKADDR*)&clientsocket, &len);//
		DWORD  dwThreadId;
		HANDLE  hThread;

		hThread = CreateThread(NULL, NULL, AnswerThreadSendModel,
			(LPVOID)clientConn, 0, &dwThreadId);
		if (hThread == NULL)
		{
			printf("CreatThread  AnswerThread()  failed.\n");
		}
		else
		{
			printf("CreateThread  OK.\n");
		}
	}
}



bool calculateProjectPoints(double Fiducial[], vector< Point2f>  & ProjectimagePoints)
{
	Point3f objectPoints[4] = { Point3i(00.00, 0.00, 0.0), Point3i(45.00, 0.0, 0.0), Point3i(0.0, 25.0, 0.0), Point3i(45.0, 25.0, 00.0) }; //this numbers are in centimeter
	Mat objectPointMat = Mat(4, 1, CV_32FC3, objectPoints);
	//vector<Point2f> imagePointMat;
	//vector<Point3f> objectPointMat;
	//objectPointMat.push_back(Point3f(0.0f, 0.0f, 0.0f));
	//objectPointMat.push_back(Point3f(45.0f, 0.0f, 0.0f));
	//objectPointMat.push_back(Point3f(0.0f, 25.0f, 0.0f));
	//objectPointMat.push_back(Point3f(45.0f, 25.0f, 0.0f));
	Point2f imagePoint[4] = { Point(Fiducial[0], Fiducial[1]), Point(Fiducial[2], Fiducial[3]), Point(Fiducial[4], Fiducial[5]), Point(Fiducial[6], Fiducial[7]) };
	Mat imagePointMat = Mat(4, 1, CV_32FC2, imagePoint);
	/*imagePointMat.push_back(Point2f((float)Fiducial[0], (float)Fiducial[1]));
	imagePointMat.push_back(Point2f((float)Fiducial[2], (float)Fiducial[3]));
	imagePointMat.push_back(Point2f((float)Fiducial[4], (float)Fiducial[5]));
	imagePointMat.push_back(Point2f((float)Fiducial[6], (float)Fiducial[7]));*/
	//float cameraParameter[3][3] = { { 481.995, 0, 399.5 }, { 0, 481.9952, 239.5 }, { 0, 0, 1 } }; // this is in pixel coordinates assuming pixel size is 5micron and focal length is 2.4 mm
	//float cameraParameter[3][3] = { { 481.995, 0, 512 }, { 0, 481.9952, 288 }, { 0, 0, 1 } }; // this is in pixel coordinates assuming pixel size is 5micron and focal length is 2.4 mm
	//float cameraParameter[3][3] = { { 481.995, 0, 0 }, { 0, 481.9952, 0 }, { 399.5, 239.5, 1 } };
	float cameraParameter[3][3] = { { 1003.40559, 0, 511.5 }, { 0, 511.5, 287.5 }, { 0, 0, 1 } };

	Mat cameraMatrix = Mat(3, 3, CV_32FC1, cameraParameter);
	/*
	for (int j = 0; j < 3; j++)
	for (int i = 0; i < 3;i++)
	printf("%f	", cameraMatrix.at<float>(i,j));

	printf("\n\n\n");*/

	//float distC[5] = { 0.031733, 0.1256, 0.0, 0.0, -0.597547193 };
	vector<float> distorVector;// (&distC[0], &distC[5]);
	distorVector.push_back(0.031733);
	distorVector.push_back(0.1256);
	distorVector.push_back(0);
	distorVector.push_back(0);
	distorVector.push_back(-0.597547193);

	/*for (int i = 0; i <distorVector.size();i++)
	printf("%f	", distorVector.at(i));

	printf("\n\n\n");*/

	//Mat cameraMatrix = Mat(3, 3, CV_32FC1, Scalar::all(0));
	//cameraMatrix.at<float>(0, 0) = 481.995f;
	//cameraMatrix.at<float>(0, 2) = 399.5f;
	//cameraMatrix.at<float>(1, 1) = 481.9952f;
	//cameraMatrix.at<float>(1, 2) = 239.5f;
	//cameraMatrix.at<float>(2, 2) = 1.0f;


	//vector<float> distorVector;
	//distorVector.push_back(0.031733);
	//distorVector.push_back(0.1256);
	//distorVector.push_back(0.0);
	//distorVector.push_back(0.0);
	//distorVector.push_back(-0.597547193);
	Mat rvec;
	Mat tvec;

	/*FileStorage fs1("data.xml", FileStorage::WRITE);
	fs1 << "objectPointMat" << objectPointMat;
	fs1 << "imagePointMat" << imagePointMat;
	fs1 << "cameraMatrix" << cameraMatrix;
	fs1 << "distorVector" << distorVector;
	fs1.release();*/
	solvePnP(objectPointMat, imagePointMat, cameraMatrix, distorVector, rvec, tvec);

	/*for (int i = 0; i< 3; i++)
	for (int j = 0; j < 1; j++)
	printf("%d	", rvec.at<float>(0, 0));

	printf("\n\n\n");*/
	Mat tvec2, rvec2;
	/*cout << "tvec is:\n";
	for (int i = 0; i < 1; i++)
	{

	tvec.convertTo(tvec2, CV_32F);
	for (int j = 0; j < tvec2.rows; j++)
	{
	printf("%f	", tvec2.at<float>(j, i));
	}
	}

	printf("\n\n\n");

	cout << "rvec is:\n";
	for (int i = 0; i < 1; i++)
	{

	rvec.convertTo(rvec2, CV_32F);
	for (int j = 0; j < rvec2.rows; j++)
	{
	printf("%f	", rvec2.at<float>(j, i));
	}
	}

	printf("\n\n\n");

	cout << "camera matrix is:\n";
	for (int i = 0; i < cameraMatrix.rows; i++)
	{

	//rvec.convertTo(rvec2, CV_32F);
	for (int j = 0; j < cameraMatrix.cols; j++)
	{
	printf("%f	", cameraMatrix.at<float>(j, i));
	}
	printf("\n");
	}

	printf("\n\n\n");


	cout << "distortion vector is:\n";
	for (int j = 0; j < distorVector.size(); j++)
	{
	printf("%f	", distorVector.at(j));
	}


	printf("\n\n\n");*/



	/////////// project 3D points to 2D iamge/////////////////////////		
	if (rvec.empty() && tvec.empty())
	{
		return false;
	}

	while (1)
	{
		if (mKinect.getWrong)
			return false;
		if (mKinect.dataReady && !mKinect.getWrong)
			break;

	}


	vector<double> projpts;
	//cout << "3d point" << mKinect.Target3D.at(0)<<"\n";
	projectPoints(mKinect.Target3D, rvec, tvec, cameraMatrix, distorVector, ProjectimagePoints);
	//getPoints(Fiducial, projpts,mKinect.Target3D);

	if (ProjectimagePoints.size() == 0)
		return false;




	memset(tempDataMatrix, 0, sizeof(tempDataMatrix));
	ofstream myfile;
	//myfile.open("print.txt");
	//cout << "size=" << projpts.size();
	cout << "project points size=" << projpts.size() << "\n";
	for (int i = 0; i < ProjectimagePoints.size(); i = i + 5)
	{
		if (ProjectimagePoints[i].x > 1020 || ProjectimagePoints[i].y > 570 || ProjectimagePoints[i].x < 0 || ProjectimagePoints[i].y < 0)
			continue;

		//myfile << ProjectimagePoints[i].x << " , " << ProjectimagePoints[i].y << "\n";
		int a = ProjectimagePoints[i].x;
		int b = ProjectimagePoints[i].y;
		tempDataMatrix[(short)(ProjectimagePoints[i].y + 0.5)][(short)(ProjectimagePoints[i].x + 0.5)] = 1;
		//tempDataMatrix[(short)(projpts.at(2*i+1) + 0.5)][(short)(projpts.at(2*i) + 0.5)] = 1;
		// myfile<< "projected points\n";
		//myfile << ProjectimagePoints[i].x << "," << ProjectimagePoints[i].y << "\n";

		//printf("%d,%d\n", ProjectimagePoints[i].x, ProjectimagePoints[i].y);
	}
	//ProjectimagePoints.clear();
	//tempDataMatrix[30][40] = 1;
	//tempDataMatrix[50][60] = 1;

	//myfile.close();
	sendPoints.clear();
	for (short i = 0; i < 576; i++)// [640][360]
		for (short j = 0; j < 1024; j++)
		{
			if (imageDataMatrix[i][j] != tempDataMatrix[i][j])
				//if (tempDataMatrix[i][j])
			{
				sendPoints.push_back(i);
				sendPoints.push_back(j);
				//printf("pushing i %d\n", i);
				//printf("pushing j %d\n", j);
				if (imageDataMatrix[i][j] == 0)
					imageDataMatrix[i][j] = 1;
				else
					imageDataMatrix[i][j] = 0;
			}
		}
	//printf("exiting calculate\n");
	return true;
}

void RegisterImages()
{
	Mat homography_transform;
	vector<Point2f> kinect_points_vector;
	vector<Point2f> camear_points_vector;
	////get kinect color images' fiducials points
	Mat kinect_image;
	Point2f kinect_points[4];
	bool is_find_kinect_points = false;
	while (!is_find_kinect_points)
	{
		mKinect.kinect_image(kinect_image);
		if (FindFiducials(kinect_image, kinect_points))
			is_find_kinect_points = true;
		else
			is_find_kinect_points = false;

		for (int i = 0; i < 4; i++)
		{
			circle(kinect_image, kinect_points[i], 4, Scalar(1, 255, 9));
		}
		imshow("kinect image", kinect_image);
		waitKey(30);
	}

	///get the manta flouce image fiducials' points 
	CameraControl mCameraControl;
	mCameraControl.CameraInitialize();
	mCameraControl.CamerasOpening();
	mCameraControl.ImagePrepar();
	vector <vector <IplImage*>> mImages;
	bool is_get_manta_points = false;
	Point2f manta_points[4];
	while (!is_get_manta_points)
	{
		mImages.clear();
		//	mCameraControl.CameraGetFrames(&mImages, 1);
		//ShowImages(mCameraControl.CameraGetFrames(&mImages, 1));
		mImages = mCameraControl.CameraGetFrames(&mImages, 1);
		IplImage* image1;
		image1 = mImages[0][0];
		Mat image = Mat(image1);

		imshow("captured image", image);
		cvWaitKey(30);
		//	resize(image,image,Size(image.size().width / 2, image.size().height / 2));
		if (is_need_manual_camera_points)
		{
			while (!GetBoxByMouse(image));
			is_get_manta_points=true;
			for (int i = 0; i < 4; i++)
			{
				circle(image, manta_points[i], 4, Scalar(1, 255, 9));
				camear_points_vector.push_back(ficucial_centers[i]);
			}
			imshow("captured image", image);
			cvWaitKey(10);
		}
		else if (FindFiducials(image, manta_points))
		{
			is_get_manta_points = true;
			for (int i = 0; i < 4; i++)
			{
				circle(image, manta_points[i], 4, Scalar(1, 255, 9));
				camear_points_vector.push_back(manta_points[i]);
			}
			imshow("captured image", image);
			cvWaitKey(10);
		}
		else
			is_get_manta_points = false;

	}
	mCameraControl.CameraStopFrameStream();


	/// i manually measured the coordinate of four fiducials

	//	camear_points_vector.push_back(Point(378, 299));
	//	camear_points_vector.push_back(Point(1140, 298));
	//	camear_points_vector.push_back(Point(386, 719));
	//	camear_points_vector.push_back(Point(1141, 712));

	for (int i = 0; i < 4; i++)
	{
		kinect_points_vector.push_back(kinect_points[i]);

	}

	homography_transform = findHomography(camear_points_vector, kinect_points_vector);

	FileStorage fs1("homography_transform.xml", FileStorage::WRITE);
	fs1 << "homography_transform" << homography_transform;
	fs1.release();
}

bool FindFloucePoints(const Mat& inputmat, vector<Point>&  flouce_poins, vector<unsigned char>&  magnitute)
{
	Mat gray_image;
	vector<Point> points;
	if (inputmat.channels() == 3)
	{
		cvtColor(inputmat, gray_image, COLOR_BGR2GRAY);
	}
	else
		inputmat.copyTo(gray_image);

	//imshow("gray debugh image", gray_image);
	for (int i = 0; i < gray_image.cols; i++)
	{
		for (int j = 0; j < gray_image.rows; j++)
		{

			if (gray_image.at<unsigned char>(j, i)>light_threshold) {
				points.push_back(Point(i, j));
				magnitute.push_back(gray_image.at<unsigned char>(j, i));

			}
		}
	}
	flouce_poins = points;
	return true;
}

bool FindFiducials(const Mat& inputmat, Point2f* fiducial)
{
	bool first_time = true;
	bool NewCalculate = true;
	int widowsize = 30;

	double outputfiducial[16];

	Mat src = inputmat;
	if (src.empty())
		return false;
	Mat gray;
	if (src.channels() == 1)
		src.copyTo(gray);
	else
		cvtColor(src, gray, COLOR_BGR2GRAY);

	//================================================first time========================
	if (NewCalculate == true)
	{
		if (find_glyphs(gray, fiducial))
		{

			return true;
		}
		else
			return false;

	}
}

bool find_glyphs(const Mat &img, Point2f *glyph_center)
{
	int pattern1[] = { 0, 1, 0, 1, 1, 1, 1, 0, 1,
		1, 1, 0, 0, 1, 1, 1, 1, 0,
		1, 0, 1, 1, 1, 1, 0, 1, 0,
		0, 1, 1, 1, 1, 0, 0, 1, 1 };

	int pattern2[] = { 0, 0, 1, 1, 0, 1, 0, 1, 0,
		0, 1, 0, 1, 0, 0, 0, 1, 1,
		0, 1, 0, 1, 0, 1, 1, 0, 0,
		1, 1, 0, 0, 0, 1, 0, 1, 0 };

	int pattern3[] = { 0, 0, 1, 0, 1, 0, 1, 0, 1,
		1, 0, 0, 0, 1, 0, 1, 0, 1,
		1, 0, 1, 0, 1, 0, 1, 0, 0,
		1, 0, 1, 0, 1, 0, 0, 0, 1 };

	int pattern4[] = { 0, 1, 0, 0, 1, 1, 0, 0, 1,
		0, 0, 0, 0, 1, 1, 1, 1, 0,
		1, 0, 0, 1, 1, 0, 0, 1, 0,
		0, 1, 1, 1, 1, 0, 0, 0, 0 };

	//int pattern5[] = {0, 0, 0, 0, 1, 0, 0, 0, 0,
	//	0, 0, 0, 0, 1, 0, 0, 0, 0,
	//	0, 0, 0, 0, 1, 0, 0, 0, 0,
	//	0, 0, 0, 0, 1, 0, 0, 0, 0};

	//int pattern6[] = {0, 1, 0, 0, 1, 0, 0, 1, 0,
	//	0, 0, 0, 1, 1, 1, 0, 0, 0,
	//	0, 1, 0, 0, 1, 0, 0, 1, 0,
	//	0, 0, 0, 1, 1, 1, 0, 0, 0};

	//int pattern7[] = {0, 0, 0, 0, 0, 0, 1, 1, 1,
	//	0, 0, 1, 0, 0, 1, 0, 0, 1,
	//	1, 1, 1, 0, 0, 0, 0, 0, 0,
	//	1, 0, 0, 1, 0, 0, 1, 0, 0};

	//int pattern8[] = {1, 0, 0, 1, 0, 0, 1, 1, 1,
	//	0, 0, 1, 0, 0, 1, 1, 1, 1,
	//	1, 1, 1, 0, 0, 1, 0, 0, 1,
	//	1, 1, 1, 1, 0, 0, 1, 0, 0};

	const int *pattern[] = { pattern1, pattern2, pattern3, pattern4 };//,pattern5};//, pattern6, pattern7, pattern8};

	const int APPROX_POLY_EPSILON = 3;
	const int MIN_CONTOUR_AREA = 50;
	const int CELL_NUM_ROW = 3;
	const int GLYPH_SIZE = 30;
	const int CELL_NUM = CELL_NUM_ROW * CELL_NUM_ROW;
	const int CELL_SIZE = GLYPH_SIZE / CELL_NUM_ROW;
	const int MAX_DELAY_FRAME = 3;
	Mat gs_src;
	gs_src = img;

	Point2f cw_dst_points[4];
	cw_dst_points[0] = Point2f(-CELL_SIZE, -CELL_SIZE);
	cw_dst_points[1] = Point2f(GLYPH_SIZE + CELL_SIZE, -CELL_SIZE);
	cw_dst_points[2] = Point2f(GLYPH_SIZE + CELL_SIZE, GLYPH_SIZE + CELL_SIZE);
	cw_dst_points[3] = Point2f(-CELL_SIZE, GLYPH_SIZE + CELL_SIZE);

	Point2f ccw_dst_points[4];
	ccw_dst_points[0] = Point2f(-CELL_SIZE, -CELL_SIZE);
	ccw_dst_points[1] = Point2f(-CELL_SIZE, GLYPH_SIZE + CELL_SIZE);
	ccw_dst_points[2] = Point2f(GLYPH_SIZE + CELL_SIZE, GLYPH_SIZE + CELL_SIZE);
	ccw_dst_points[3] = Point2f(GLYPH_SIZE + CELL_SIZE, -CELL_SIZE);

	Mat bw_src;
	double high_threshold = threshold(gs_src, bw_src, 0, 255, THRESH_BINARY | THRESH_OTSU);
	Mat edges;
	Canny(gs_src, edges, 0.5*high_threshold, high_threshold);
	vector<vector<Point> > contours;
	findContours(edges, contours, RETR_LIST, CHAIN_APPROX_SIMPLE);

	vector<vector<Point> > vertices;
	vertices.clear();
	for (int i = 0; i < contours.size(); i++)
	{


		vector<Point> cnt = contours.at(i);
		vector<Point> v;
		approxPolyDP(cnt, v, APPROX_POLY_EPSILON, true);
		if (v.size() == 4 && contourArea(v) > MIN_CONTOUR_AREA && isContourConvex(v))
		{
			vertices.push_back(v);
			//	polylines(src,contours[i], true, cv::Scalar(155,255,55));			
		}
	}
	//	imshow("plylines",src);
	//	waitKey(20);

	if (vertices.size() < 4)
		return false;
	Mat m;
	Size dsize(GLYPH_SIZE, GLYPH_SIZE);
	Mat glyph(dsize, CV_8UC1);
	int p[CELL_NUM + 1];
	int found[4] = { 0, 0, 0, 0 };//, 0, 0, 0, 0};
	Point2f fv[4];

	for (int g = 0; g < vertices.size(); g++)
	{
		for (int i = 0; i < 4; i++)
		{
			Point2f f = Point(vertices.at(g).at(i));
			fv[i] = f;
		}

		if (contourArea(vertices.at(g), true) < 0)
			m = getPerspectiveTransform(fv, cw_dst_points);
		else
			m = getPerspectiveTransform(fv, ccw_dst_points);

		warpPerspective(gs_src, glyph, m, dsize);
		threshold(glyph, glyph, 0, 255, THRESH_BINARY | THRESH_OTSU);
		imshow("glyph", glyph);
		waitKey(20);
		for (int i = 0; i < CELL_NUM; i++)
		{
			int sum = 0;
			int row = (i / 3) * CELL_SIZE;
			int col = (i % 3) * CELL_SIZE;
			for (int r = row; r < row + CELL_SIZE; r++)
				for (int c = col; c<col + CELL_SIZE; c++)
					sum += (int)glyph.at<uchar>(r, c);
			if (sum / 255 > CELL_SIZE*CELL_SIZE / 2)
				p[i] = 1;
			else
				p[i] = 0;
		}

		for (int j = 4, i = 0; i < 4 && j == 4; i++)
		{
			for (j = 0; j < 4; j++)
			{
				int k;
				for (k = 0; k < CELL_NUM; k++)
				{
					if (pattern[i][j*CELL_NUM + k] != p[k])
						break;
				}
				if (k == CELL_NUM)
				{
					Point2f center;
					center.x = (fv[0].x + fv[1].x + fv[2].x + fv[3].x)*0.25;
					center.y = (fv[0].y + fv[1].y + fv[2].y + fv[3].y)*0.25;

					if (!found[i])
					{
						glyph_center[i] = center;
						found[i]++;
						break;
					}
				}
			}
		}
	}

	if (found[0] && found[1] && found[2] && found[3])// && found[4] && found[5] && found[6] && found[7])
	{
		return true;
	}
	else
		return false;

}

void GetFlouce()
{
	Mat DataTransform;
	FileStorage fsmy2("homography_transform.xml", FileStorage::READ);
	fsmy2["homography_transform"] >> DataTransform;
	fsmy2.release();
	double a1 = DataTransform.at<double>(0, 0);
	double a2 = DataTransform.at<double>(0, 1);
	double a3 = DataTransform.at<double>(0, 2);

	double b1 = DataTransform.at<double>(1, 0);
	double b2 = DataTransform.at<double>(1, 1);
	double b3 = DataTransform.at<double>(1, 2);

	double c1 = DataTransform.at<double>(2, 0);
	double c2 = DataTransform.at<double>(2, 1);
	double c3 = DataTransform.at<double>(2, 2);

	CameraControl mCameraControl;
	mCameraControl.CameraInitialize();
	mCameraControl.CamerasOpening();
	mCameraControl.ImagePrepar();
	vector <vector <IplImage*>> mImages;
	for (;;)

	{
		light_threshold = colorInt;
		mImages.clear();
		//	mCameraControl.CameraGetFrames(&mImages, 1);
		//ShowImages(mCameraControl.CameraGetFrames(&mImages, 1));
		mImages = mCameraControl.CameraGetFrames(&mImages, 1);
		IplImage* image1;
		image1 = mImages[0][0];
		Mat image = Mat(image1);
		//resize(image,image,Size(image.size().width / 2, image.size().height / 2));
		vector<Point> flouce_poins_image;
		vector<unsigned char> flouce_magnitute;




		//Mat test_image=imread("2.tif",1);
		//imshow("sfsdf",test_image);
		//waitKey(20);
		//Mat test_image(image);
		//FindFloucePoints(test_image, flouce_poins_image); // use threshold value to find the fluorescence part. the threshold value in the function is 10;

		//for (int i = 0; i<flouce_poins_image.size(); i++)
		//{
		//circle(test_image, flouce_poins_image[i], 1, Scalar(0, 255, 0));
		//}
		//imshow("image", test_image);
		//waitKey(20);

		//	imshow("image", image);// original image from camera prash
		/*if (image.type() == CV_8UC1)
			cout << "type of image=yes\n";

			else
			cout << "type of image=no\n";*/
		Mat outImage;
		Mat temp; 
		cvtColor(image,outImage , CV_GRAY2BGR );
		//merge(temp, 3, outImage);

		cvReleaseImage(&img1);
		img1 = cvCloneImage(&(IplImage)outImage);

		//DispImage = updateGui(DispImage, img1, img2);
		/*cvShowImage(name, DispImage);
		cvWaitKey(15);*/


		//waitKey(20);

		//	imwrite("original image.jpg", image);
		FindFloucePoints(image, flouce_poins_image, flouce_magnitute);// use threshold value to find the fluorescence part. the threshold value in the function is 10;

		//	Mat image_show;
		Mat image_show(outImage);
		
		//image_show = Mat(image.cols, image.rows, );
		//	image.copyTo(image_show);
		//	image_show = Scalar(125, 135, 124);
		for (int i = 0; i < flouce_poins_image.size(); i++)
		{
			float magnitude = (flouce_magnitute[i] / 255.0 * 50) + 205;
			circle(image_show, flouce_poins_image[i], 1, Scalar(0, magnitude, 0));
			//circle(image_show, flouce_poins_image[i], 1, Scalar(0, 255, 0));
		}
		//		imshow("image_show", image_show);
		//		waitKey(20);
		///////////transform camera points to kinect points
		while (mKinect.is_flouce_points_kinect_using);
		mKinect.is_flouce_points_kinect_ready = false;
		mKinect.flouce_points_kinect.clear();
		for (int i = 0; i < flouce_poins_image.size(); i++)
		{
			float x = flouce_poins_image[i].x;
			float y = flouce_poins_image[i].y;
			Point F;
			F.x = (x * a1 + y * a2 + a3) / (c1*x + c2*y + c3);

			F.y = (x * b1 + y * b2 + b3) / (c1*x + c2*y + c3);

			mKinect.flouce_points_kinect.push_back(F);
		}

		mKinect.is_flouce_points_kinect_ready = true;

		/////////
		//cvShowImage(mKinect.name, mKinect.DispImage);
		//imshow("captured image", image_show);
		//cvWaitKey(20);
		//	imwrite("captured image.jpg", image_show);
	//	Mat outImage2;
		//Mat temp2[] = { image_show, image_show, image_show };
		//merge(temp2, 3, outImage2);

		
		//cvtColor(image_show,outImage2, CV_GRAY2BGR);
		

		cvReleaseImage(&img2);
		img2 = cvCloneImage(&(IplImage)image_show);

		DispImage = updateGui(DispImage, img1, img2);
		//cvShowImage(name, DispImage);
		//cvShowImage("new", DispImage);
		Mat tempIm(DispImage);
		imshow(name, tempIm);
		//imshow("changllin", tempIm);
		cvWaitKey(5);


		if (colorInt2 == 2)
		{
			HANDLE handle = CreateThread(NULL, 0, KinectThread, NULL, 0, NULL);
			if (handle == NULL)
			{
				printf("CreatThread  KinectThread  failed.\n");
			}
			else
			{
				printf("CreateThread KinectThread OK.\n");
			}
			test();
			break;
		}

	}
	
	mCameraControl.CameraStopFrameStream();
	
	//cvShowImage(mKinect.name, mKinect.DispImage);
	//	return false;
}

void TestWithFlouce()
{

	

	GetFlouce();

	//mKinect.askData = true;
	////cvShowImage(mKinect.name, mKinect.DispImage);
	//
	//HANDLE handle = CreateThread(NULL, 0, KinectThread, NULL, 0, NULL);

	//if (handle == NULL)
	//{
	//	printf("CreatThread  KinectThread  failed.\n");
	//}
	//else
	//{
	//	printf("CreateThread KinectThread OK.\n");
	//}

	//while (1)
	//{
	//	//cvShowImage(mKinect.name, mKinect.DispImage);
	//}
}

bool GetBoxByMouse(Mat& inputimage)
{

	Mat temporary_image;
	inputimage.copyTo(temporary_image);
	cv::namedWindow("CameraImage", 1);
	cvSetMouseCallback("CameraImage", mouseHandler, NULL);
	while (!is_get_fiducials)
	{
		for (int i = 0; i < ficucial_centers.size(); i++)
		{
			circle(temporary_image, ficucial_centers[i], 4, Scalar(0, 212, 111), 1);
		}
		imshow("CameraImage", temporary_image);
		cv::waitKey(10);
	}
	cvSetMouseCallback("FirstWindow", NULL, NULL);
	return true;

}

void mouseHandler(int event, int x, int y, int flags, void *param){
	switch (event){
	case CV_EVENT_LBUTTONDOWN:
	{
		ficucial_centers.push_back(Point(x, y));
		if (ficucial_centers.size() == 4)
			is_get_fiducials = true;
		break;
	}
	}
}