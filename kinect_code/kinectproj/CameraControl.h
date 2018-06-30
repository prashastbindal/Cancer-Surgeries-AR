#pragma once
#include <PvApi.h>
#include "windows.h"		// for Sleep()
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/core/core.hpp>
#include <ppl.h>			// Parallel processing
#include <vector>

using std::vector;
	typedef struct 
	{
    unsigned long   UID;
    tPvHandle       Handle;
    tPvFrame        Frame;
    tPvUint32       Counter;
    char            Filename[20];

	}tCamera;
class CameraControl
{


public:
	CameraControl(void);
	~CameraControl(void);

private:
	int showingratio ;
	static const unsigned long nCams=1;
	int IsPacketOK ;
	unsigned long nFrames ;
	DWORD time_start, time_end;
	DWORD time_start_saving, time_end_saving;
	double elapedtime_saving;
	double fps;
	tPvErr errCode;
	tCamera* Camera;
	unsigned long	numCameras;

	vector <vector <IplImage*>> origImg;
	vector <IplImage*> outImag[nCams];
	
	IplImage** img_ori;
	IplImage** img_sml_BW;
	IplImage *img_shw_big;

	unsigned long StreamBytesPerSecondPerCamera;		// StreamBytesPerSecond
	unsigned long ExposureValue;
public:
	void CameraControl::CameraSetParemeter();
	void CameraControl::CameraInitialize();
	bool CameraControl::CamerasOpening();
	void CameraControl::ImagePrepar();
	void CameraControl::ImageSaving();
	vector <vector <IplImage*>>  CameraControl::CameraGetFrames(vector <vector <IplImage*>> * mvector, int times);
	void CameraControl::CameraStarFrameStream();
	void CameraControl::CameraStopFrameStream();
};

