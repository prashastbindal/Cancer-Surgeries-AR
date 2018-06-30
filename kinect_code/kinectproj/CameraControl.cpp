#include "CameraControl.h"


CameraControl::CameraControl(void)
{
	showingratio = 2;
	IsPacketOK = 1;
	nFrames = 100000000;
	fps = 1.0;
	Camera = new tCamera [nCams];

}


CameraControl::~CameraControl(void)
{
}

void CameraControl::CameraSetParemeter()
	{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// *I can change this value according to network condition*
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// StreamBytesPerSecond (MegaBytes per second)
//	double NetworkStreamMegaBytesPerSecond = 100;		//	double NetworkStreamMegaBytesPerSecond = 60;
	double NetworkStreamMegaBytesPerSecond =  90;		//	double NetworkStreamMegaBytesPerSecond = 60;
	
	// Camera exposure time [sec] (Auto: ExporeTimeInSecond = -1; manual-> auto; needs to restart cameras)
//	double ExposureTimeInSecond = 0.02;	
	double ExposureTimeInSecond = 0.01;	

	// Calculating Stream bytes per seconds
//	StreamBytesPerSecondPerCamera = (unsigned long)(NetworkStreamMegaBytesPerSecond*1000000.0/(double)nCams);
	StreamBytesPerSecondPerCamera = (unsigned long)(NetworkStreamMegaBytesPerSecond*1024*1024/(double)nCams);

	// Calculating ExposureValue
	ExposureValue =  (unsigned long)(ExposureTimeInSecond*1000000);
	}


void CameraControl::CameraInitialize()
{
	CameraSetParemeter();
	for(int i=0;i<nCams;i++)
	{
		memset(&Camera[i], 0, sizeof(tCamera));
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialize the PvAPI
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	PvInitialize();   

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for cameras
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	printf("Waiting for a camera");
    
	while(PvCameraCount()<nCams)
    {
        printf(".");
		Sleep(100);
    }
    printf("\n\n");
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Getting cameras
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	tPvCameraInfoEx*	list;
	list = new tPvCameraInfoEx [nCams];
	tPvUint32    connected;
	numCameras = PvCameraListEx(list, nCams, &connected, sizeof(tPvCameraInfoEx));
	printf("Number of detected camera: %d\n\n", numCameras);
	for(int i=0;i<numCameras;i++)
		printf("%s [ID %u]\n", list[i].SerialNumber, list[i].UniqueId);
	printf("\n");
	int selcount = 0;
	if(numCameras>0)
    {
		for(int i=0;i<numCameras;i++)
		{
			Camera[selcount].UID = list[selcount].UniqueId;
			printf("Selected camera%d: %s\n",selcount, list[selcount].SerialNumber);
			selcount++;
		}
		printf("\n");
    }
	delete[] list;
}

bool CameraControl::CamerasOpening()
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Opening cameras
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(int i=0;i<numCameras;i++)
	{
		printf("Opening camera%d: ", i);
		if ((errCode = PvCameraOpen(Camera[i].UID,ePvAccessMaster,&(Camera[i].Handle))) != ePvErrSuccess)
		{
			if (errCode == ePvErrAccessDenied)
			{
				printf("PvCameraOpen returned ePvErrAccessDenied:\nCamera already open as Master, or camera wasn't properly closed and still waiting to HeartbeatTimeout.");
				return false;
			}
			else
			{
				printf("PvCameraOpen err: %u\n", errCode);
				return false;
			}
		}
		else
		{
			printf("success\n");
		}
	}
	return true;
}

void CameraControl::ImageSaving()
{
	// Calculate frame buffer size & allocate image buffer
	unsigned long FrameSize = 0;
	for(int i=0;i<numCameras;i++)
	{
		PvAttrUint32Get(Camera[i].Handle,"TotalBytesPerFrame",&FrameSize);
		Camera[i].Frame.ImageBuffer = new char[FrameSize];
		Camera[i].Frame.ImageBufferSize = FrameSize;	
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Image acquisition
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	for(int i=0;i<numCameras;i++)
	{
		PvCaptureAdjustPacketSize(Camera[i].Handle, 8228);
		PvCaptureStart(Camera[i].Handle);
		PvCaptureQueueFrame(Camera[i].Handle, &(Camera[i].Frame), NULL);
		PvAttrUint32Set(Camera[i].Handle, "StreamBytesPerSecond", StreamBytesPerSecondPerCamera);	
		PvAttrEnumSet(Camera[i].Handle, "PixelFormat", "Mono8");	
	//	PvAttrEnumSet(Camera[i].Handle, "ExposureMode", "Auto");	
		PvAttrEnumSet(Camera[i].Handle, "ExposureMode", "Manual");	
		PvAttrUint32Set(Camera[i].Handle, "ExposureValue", ExposureValue);
		PvAttrEnumSet(Camera[i].Handle, "GainMode", "Auto");	
	//	PvAttrEnumSet(Camera[i].Handle, "GainMode", "AutoOnce");
		PvAttrEnumSet(Camera[i].Handle, "FrameStartTriggerMode", "Software");
		PvAttrEnumSet(Camera[i].Handle, "AcquisitionMode", "Continuous");
		PvCommandRun(Camera[i].Handle,"AcquisitionStart");
	}
	
	unsigned long nSavedFrames = 1000;
	
	IplImage** img_ori;
	img_ori = new IplImage* [numCameras];

	IplImage** img_sml_BW;
	img_sml_BW = new IplImage* [numCameras];

	for(int i=0;i<numCameras;i++)
	{
		img_ori[i]		= cvCreateImage(cvSize(Camera[i].Frame.Width, Camera[i].Frame.Height), 8, 1);
		img_sml_BW[i]	= cvCreateImage(cvSize(img_ori[i]->width/showingratio, img_ori[i]->height/showingratio), 8, 1);
	}
	
	// Images for showing	
	IplImage *img_shw_big = cvCreateImage(cvSize(img_sml_BW[0]->width*2 + 10, img_sml_BW[0]->height), IPL_DEPTH_8U, 1); 

	// cvRect for copying images to a big image to show
	CvRect* rtShwBig = new CvRect[numCameras]; 	
	for(int i=0;i<numCameras;i++)
		rtShwBig[i] = cvRect(0, 0, img_sml_BW[0]->width,  img_sml_BW[0]->height);
	 
	// Setting values for cvRect
	if(numCameras>0)
	{
		rtShwBig[0].x = 0;							
		rtShwBig[0].y = 0;				
	}
	if(numCameras>1)
	{
		rtShwBig[1].x = img_sml_BW[0]->width + 10;		
		rtShwBig[1].y = 0;	
	}
	if(numCameras>2)
	{
		rtShwBig[2].x = 0;							
		rtShwBig[2].y = img_sml_BW[0]->height + 10;	
	}
	if(numCameras>3)
	{
		rtShwBig[3].x = img_sml_BW[0]->width + 10;		
		rtShwBig[3].y = img_sml_BW[0]->height + 10;	
	}

	int Sizeof_SavedImageFileName = 50;

	char** SavedImageFileName = new char* [numCameras];
	for(int i=0;i<numCameras;i++)
		SavedImageFileName[i] = new char [Sizeof_SavedImageFileName];
}

void CameraControl::ImagePrepar()
{
	// Calculate frame buffer size & allocate image buffer
	unsigned long FrameSize = 0;
	for(int i=0;i<numCameras;i++)
	{
		PvAttrUint32Get(Camera[i].Handle,"TotalBytesPerFrame",&FrameSize);
		Camera[i].Frame.ImageBuffer = new char[FrameSize];
		Camera[i].Frame.ImageBufferSize = FrameSize;	
	}
	
	for(int i=0;i<numCameras;i++)
	{
		PvCaptureAdjustPacketSize(Camera[i].Handle, 8228);
		PvCaptureStart(Camera[i].Handle);
		PvCaptureQueueFrame(Camera[i].Handle, &(Camera[i].Frame), NULL);
		PvAttrUint32Set(Camera[i].Handle, "StreamBytesPerSecond", StreamBytesPerSecondPerCamera);	
		PvAttrEnumSet(Camera[i].Handle, "PixelFormat", "Mono8");	
	//	PvAttrEnumSet(Camera[i].Handle, "ExposureMode", "Auto");	
		PvAttrEnumSet(Camera[i].Handle, "ExposureMode", "Manual");	
		PvAttrUint32Set(Camera[i].Handle, "ExposureValue", ExposureValue);
		PvAttrEnumSet(Camera[i].Handle, "GainMode", "Auto");	
	//	PvAttrEnumSet(Camera[i].Handle, "GainMode", "AutoOnce");
		PvAttrEnumSet(Camera[i].Handle, "FrameStartTriggerMode", "Software");
		PvAttrEnumSet(Camera[i].Handle, "AcquisitionMode", "Continuous");
		PvCommandRun(Camera[i].Handle,"AcquisitionStart");
	}

	img_ori = new IplImage* [numCameras];

	
	img_sml_BW = new IplImage* [numCameras];
	
	for(int i=0;i<numCameras;i++)
	{
		img_ori[i]		= cvCreateImage(cvSize(Camera[i].Frame.Width, Camera[i].Frame.Height), 8, 1);
		img_sml_BW[i]	= cvCreateImage(cvSize(img_ori[i]->width/showingratio, img_ori[i]->height/showingratio), 8, 1);
	}
		// Images for showing	
	img_shw_big = cvCreateImage(cvSize(img_sml_BW[0]->width*2 + 10, img_sml_BW[0]->height), IPL_DEPTH_8U, 1); 

	// cvRect for copying images to a big image to show
	CvRect* rtShwBig = new CvRect[numCameras]; 	
	for(int i=0;i<numCameras;i++)
		rtShwBig[i] = cvRect(0, 0, img_sml_BW[0]->width,  img_sml_BW[0]->height);
	 
	// Setting values for cvRect
	if(numCameras>0)
	{
		rtShwBig[0].x = 0;							
		rtShwBig[0].y = 0;				
	}
	if(numCameras>1)
	{
		rtShwBig[1].x = img_sml_BW[0]->width + 10;		
		rtShwBig[1].y = 0;	
	}
	if(numCameras>2)
	{
		rtShwBig[2].x = 0;							
		rtShwBig[2].y = img_sml_BW[0]->height + 10;	
	}
	if(numCameras>3)
	{
		rtShwBig[3].x = img_sml_BW[0]->width + 10;		
		rtShwBig[3].y = img_sml_BW[0]->height + 10;	
	}






}



vector <vector <IplImage*>>  CameraControl::CameraGetFrames(vector <vector <IplImage*>> * mvector, int times)
{
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Image acquisition
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	origImg.clear();
	outImag[0].clear();
//	outImag[1].clear();
	int IsCaptureTrue = 0;
	int IsFrameCaptureTrue = 0;
	int FrameCount = 0;
	int IsCamPacketOK[4];
	for(int time=0;time<times;time++)
	{
		// Software trigger
		for(int i=0;i<numCameras;i++)
			if(PvCommandRun(Camera[i].Handle,"FrameStartTriggerSoftware")!=ePvErrSuccess)
				printf("\nFrameStartTriggerSoftware: Cam%d", i);

		// Wait for camera to return to host
		for(int i=0;i<numCameras;i++)
			if(PvCaptureWaitForFrameDone(Camera[i].Handle, &(Camera[i].Frame), PVINFINITE)!=ePvErrSuccess)
				printf("\nPvCaptureWaitForFrameDone: Cam%d", i);
			
		//check returned Frame.Status
		IsPacketOK = 1;
		for(int i=0;i<numCameras;i++)
		{			
			if(Camera[i].Frame.Status!=ePvErrSuccess)
				IsCamPacketOK[i] = 0;			
			else
				IsCamPacketOK[i] = 1;

			IsPacketOK *= IsCamPacketOK[i];
		}
		// Print out Packet status
		printf("\nFrame-%06d, Packet chk: ", FrameCount);
		for(int i=0;i<numCameras;i++)
			printf("%s", IsCamPacketOK[i]?"o":"x");
		printf("=%s,", IsPacketOK?"GOOD":"FAIL");
		
		// Copy images to memory

	//	for(int i=0;i<numCameras;i++)

		Concurrency::parallel_for((unsigned long)0, numCameras, [&](unsigned long i)
		{
			memcpy(img_ori[i]->imageDataOrigin, (char*)Camera[i].Frame.ImageBuffer, Camera[i].Frame.ImageBufferSize);

//			cvResize(img_ori[i], img_sml_BW[i], CV_INTER_NN);

			// Cameras are mounted upside down
		//	cvFlip(img_sml_RGB[i], img_sml_RGB[i], -1);
		});
		for (int i=0;i<numCameras;i++)
		{
			outImag[i].push_back(img_ori[i]);
		}
		for(int i=0;i<numCameras;i++)
		{
			if(PvCaptureQueueFrame(Camera[i].Handle, &(Camera[i].Frame), NULL)!=ePvErrSuccess)
			{
				printf("PvCaptureQueueFrame: Cam%d", i);
			}
		}
		
		// Showing additional info
		if(!IsPacketOK)
			printf(" **PACKET LOSS**");

	}
	for (int i=0;i<numCameras;i++)
	{
		origImg.push_back(outImag[i]);
//		* mvector.push_back(outImag[i]);
	}

	///////////////////////////////////
	// end the camera
	////////////////////////////////////////////
	//for(int i=0;i<numCameras;i++)
	//{
	//	PvCommandRun(Camera[i].Handle, "AcquisitionStop");
	//	PvCaptureQueueClear(Camera[i].Handle);
	//	PvCaptureEnd(Camera[i].Handle);
	//	PvCameraClose(Camera[i].Handle);
	//	delete [] (char*)Camera[i].Frame.ImageBuffer;
	//}

	//PvUnInitialize();

	//cvDestroyAllWindows();




	//delete[] Camera;
	return origImg;

}

void CameraControl::CameraStopFrameStream()
{
	for(int i=0;i<numCameras;i++)
	{
		PvCommandRun(Camera[i].Handle, "AcquisitionStop");
		PvCaptureQueueClear(Camera[i].Handle);
		PvCaptureEnd(Camera[i].Handle);
		PvCameraClose(Camera[i].Handle);
		delete [] (char*)Camera[i].Frame.ImageBuffer;
	}

	PvUnInitialize();

	//cvDestroyAllWindows();
	for(int i=0;i<numCameras;i++)
	{
		cvReleaseImage(&img_ori[i]);
		cvReleaseImage(&img_sml_BW[i]);
	}
	delete[] img_ori;
	delete[] img_sml_BW;
	delete[] Camera;
}