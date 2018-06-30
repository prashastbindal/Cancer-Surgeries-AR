#include"Project3D.h"


using namespace std;
using namespace cv;

bool Project3D:: Project3Dto2D(Camera cam, vector<Point3f> kinect3dCoord,  vector< Point2f>  & ProjectimagePoints)
{

	vector<Point2f> tempProjectimagePoints;
	projectPoints(kinect3dCoord, cam.rvec, cam.tvec, cam.cameraMatrix, cam.distorVector, tempProjectimagePoints);


	if (tempProjectimagePoints.size() == 0)
		return false;




	//memset(tempDataMatrix, 0, sizeof(tempDataMatrix));

	ProjectimagePoints.clear();
	for (int i = 0; i<tempProjectimagePoints.size(); i=i+5)
	{
		if (!(tempProjectimagePoints[i].x>1020 || tempProjectimagePoints[i].y > 570 || tempProjectimagePoints[i].x<0|| tempProjectimagePoints[i].y<0))
			ProjectimagePoints.push_back(tempProjectimagePoints.at(i));


		//tempDataMatrix[(short)(ProjectimagePoints[i].y + 0.5)][(short)(ProjectimagePoints[i].x + 0.5)] = 1;

	}


	//myfile.close();
/*	sendPoints.clear();
	for (short i = 0; i<576; i++)// [640][360]
		for (short j = 0; j<1024; j++)
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
		}*/
	//printf("exiting calculate\n");
	return true;

}
