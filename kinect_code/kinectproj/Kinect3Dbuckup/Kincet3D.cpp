#include "Kincet3D.h"
using namespace std;
using namespace cv;

template<class Interface>
inline void SafeRelease(Interface *& pInterfaceToRelease)
{
	if (pInterfaceToRelease != NULL){
		pInterfaceToRelease->Release();
		pInterfaceToRelease = NULL;
	}
}
Kincet3D::Kincet3D(void)
{

	dataReady = false;
	dataReserved = false;
	isFinddepth = false;
	isSaved = false;
	isSave = false;
	askData = false;
	saveOriginalData = false;
	isSaveCloudPoints = false;
	getWrong = false;
	minRow = 140;//80;//512 * 424//down
	maxRow = 170;//250;//up
	minCol = 250;//170;//
	maxCol = 300;//400;
}


Kincet3D::~Kincet3D(void)
{
}

int Kincet3D::KinectCoordinateCalibrate()
{
	cv::setUseOptimized(true);

	// Sensor
	IKinectSensor* pSensor;
	HRESULT hResult = S_OK;
	hResult = GetDefaultKinectSensor(&pSensor);
	if (FAILED(hResult)){
		std::cerr << "Error : GetDefaultKinectSensor" << std::endl;
		return -1;
	}

	hResult = pSensor->Open();
	if (FAILED(hResult)){
		std::cerr << "Error : IKinectSensor::Open()" << std::endl;
		return -1;
	}

	IColorFrameSource* pColorSource;
	hResult = pSensor->get_ColorFrameSource(&pColorSource);
	if (FAILED(hResult)){
		std::cerr << "Error : IKinectSensor::get_ColorFrameSource()" << std::endl;
		return -1;
	}

	IDepthFrameSource* pDepthSource;
	hResult = pSensor->get_DepthFrameSource(&pDepthSource);
	if (FAILED(hResult)){
		std::cerr << "Error : IKinectSensor::get_DepthFrameSource()" << std::endl;
		return -1;
	}

	// Reader
	IColorFrameReader* pColorReader;
	hResult = pColorSource->OpenReader(&pColorReader);
	if (FAILED(hResult)){
		std::cerr << "Error : IColorFrameSource::OpenReader()" << std::endl;
		return -1;
	}

	IDepthFrameReader* pDepthReader;
	hResult = pDepthSource->OpenReader(&pDepthReader);
	if (FAILED(hResult)){
		std::cerr << "Error : IDepthFrameSource::OpenReader()" << std::endl;
		return -1;
	}

	// Description
	IFrameDescription* pColorDescription;
	hResult = pColorSource->get_FrameDescription(&pColorDescription);
	if (FAILED(hResult)){
		std::cerr << "Error : IColorFrameSource::get_FrameDescription()" << std::endl;
		return -1;
	}

	int colorWidth = 0;
	int colorHeight = 0;
	pColorDescription->get_Width(&colorWidth); // 1920
	pColorDescription->get_Height(&colorHeight); // 1080
	unsigned int colorBufferSize = colorWidth * colorHeight * 4 * sizeof(unsigned char);

	cv::Mat colorBufferMat(colorHeight, colorWidth, CV_8UC4);
	cv::Mat colorMat(colorHeight / 2, colorWidth / 2, CV_8UC4);
	cv::namedWindow("Color");

	IFrameDescription* pDepthDescription;
	hResult = pDepthSource->get_FrameDescription(&pDepthDescription);
	if (FAILED(hResult)){
		std::cerr << "Error : IDepthFrameSource::get_FrameDescription()" << std::endl;
		return -1;
	}

	int depthWidth = 0;
	int depthHeight = 0;
	pDepthDescription->get_Width(&depthWidth); // 512
	pDepthDescription->get_Height(&depthHeight); // 424
	unsigned int depthBufferSize = depthWidth * depthHeight * sizeof(unsigned short);

	cv::Mat depthBufferMat(depthHeight, depthWidth, CV_16SC1);
	cv::Mat depthMat(depthHeight, depthWidth, CV_8UC1);
	cv::namedWindow("Depth");

	// Coordinate Mapper
	ICoordinateMapper* pCoordinateMapper;
	hResult = pSensor->get_CoordinateMapper(&pCoordinateMapper);
	if (FAILED(hResult)){
		std::cerr << "Error : IKinectSensor::get_CoordinateMapper()" << std::endl;
		return -1;
	}

	cv::Mat coordinateMapperMat(depthHeight, depthWidth, CV_8UC4);
	string name = "jjjkk";
	cv::namedWindow(name);

	unsigned short minDepth, maxDepth;
	pDepthSource->get_DepthMinReliableDistance(&minDepth);
	pDepthSource->get_DepthMaxReliableDistance(&maxDepth);

	while (1)
	{
		// Depth Frame
		IDepthFrame* pDepthFrame = nullptr;
		hResult = pDepthReader->AcquireLatestFrame(&pDepthFrame);
		if (SUCCEEDED(hResult))
		{
			hResult = pDepthFrame->AccessUnderlyingBuffer(&depthBufferSize, reinterpret_cast<UINT16**>(&depthBufferMat.data));
			if (SUCCEEDED(hResult))
			{
				depthBufferMat.convertTo(depthMat, CV_8U, -255.0f / 4500.0f, 255.0f);
			}
		}

		// Color Frame
		IColorFrame* pColorFrame = nullptr;
		hResult = pColorReader->AcquireLatestFrame(&pColorFrame);
		if (SUCCEEDED(hResult))
		{
			hResult = pColorFrame->CopyConvertedFrameDataToArray(colorBufferSize, reinterpret_cast<BYTE*>(colorBufferMat.data), ColorImageFormat::ColorImageFormat_Bgra);
			if (SUCCEEDED(hResult))
			{
				cv::resize(colorBufferMat, colorMat, cv::Size(), 0.5, 0.5);
				colorMat.copyTo(src);
				if(FindFiducials(colorMat))
				{

					for (int i = 0; i<8; i++)
					{

						circle(colorMat, fiducial[i], 10, cv::Scalar(250.));
						fiducial[i].x = fiducial[i].x * 2;						//because the fiducials were found in the image which was of half size. To find the actual position we have to double it
						fiducial[i].y = fiducial[i].y * 2;
						fiducial[i].x = static_cast<int>(std::floor(fiducial[i].x + 0.5));
						fiducial[i].y = static_cast<int>(std::floor(fiducial[i].y + 0.5));


					}
					int foundNumber = 0;
					if (!isFinddepth)
					{

						vector<CameraSpacePoint> referencePoint(colorBufferMat.rows*colorBufferMat.cols);
						hResult = pCoordinateMapper->MapColorFrameToCameraSpace(512 * 424, reinterpret_cast<UINT16*>(depthBufferMat.data), colorBufferMat.rows*colorBufferMat.cols, &referencePoint[0]);
						if (SUCCEEDED(hResult))
						{

							for (int i = 0; i<8; i++)
							{
								Point3f F = Point3f(0, 0, 0);
								first[i] = F;
								int number = 0;

								for (int m = -10; m <= 10; m++)
									for (int n = -10; n <= 10; n++)
									{
										int aa = (fiducial[i].y + m)*colorBufferMat.cols + fiducial[i].x + n;
										double a1 = referencePoint[aa].Z;
										if (referencePoint[aa].Z>0 && referencePoint[aa].Z<10)
										{
											F.x = referencePoint[aa].X + F.x;
											F.y = referencePoint[aa].Y + F.y;
											F.z = referencePoint[aa].Z + F.z;
											number++;
										}
									}
								if (number >= 10)
								{
									F.x = F.x / number * 100;								//taking the average of x coordinates in 10 neighbourhood of the pixel
									F.y = F.y / number * 100;
									F.z = F.z / number * 100;
									first[i] = F;
									foundNumber++;

								}

							}


						}
					}
					if (foundNumber = 8)
					{
						isFinddepth = true;
					}
				}

			}

		}



		///////////////////////////////////
		// Mapping (Depth to Color)


		SafeRelease(pColorFrame);
		SafeRelease(pDepthFrame);

		cv::imshow("Color", colorMat);
		cv::imshow("Depth", depthMat);

		if (cv::waitKey(10) == VK_SPACE)
		{
			if (isFinddepth)

				break;
		}

	}

	/////////////////
	std::vector<cv::Point3f> mfirst, msecond;
	std::vector<uchar> inliers;
	cv::Mat aff(3, 4, CV_32F);
	for (int i = 0; i<8; i++)
	{
		mfirst.push_back(first[i]);
		
	}


	/*msecond.push_back(Point3f(00.00, 0.00, 0.0));
	msecond.push_back(Point3f(45.00, 0.0, 0.0));
	msecond.push_back(Point3f(0.0, 25.0, 0.0));
	msecond.push_back(Point3f(45.0, 25.0, 00.0));
	msecond.push_back(Point3f(10.0, 25.0, 6.0));
	msecond.push_back(Point3f(0.0, 12.5, 6.0));
	msecond.push_back(Point3f(25.0, 0.0, 6.0));
	msecond.push_back(Point3f(35.0, 45.0, 6.0));*/

	/*msecond.push_back(Point3f(00.00, 0.00, 0.0));
	msecond.push_back(Point3f(45.00, 0.0, 0.0));
	msecond.push_back(Point3f(0.0, 25.0, 0.0));
	msecond.push_back(Point3f(45.0, 25.0, 00.0));
	msecond.push_back(Point3f(45.0, 12.5, 6.7));
	msecond.push_back(Point3f(0.0, 12.5, 6.7));
	msecond.push_back(Point3f(22.5, 25.0, 6.0));
	msecond.push_back(Point3f(22.5, 0.0, 6.0));*/


	msecond.push_back(Point3f(00.00, 0.00, 0.0));
	msecond.push_back(Point3f(0.0, 45.00, 0.0));
	msecond.push_back(Point3f(25.0, 0.0, 0.0));
	msecond.push_back(Point3f(25.0, 45.0, 00.0));
	msecond.push_back(Point3f(12.5, 45.0, 6.7));
	msecond.push_back(Point3f(12.5, 0.0, 6.7));
	msecond.push_back(Point3f(25.0, 22.5, 6.0));
	msecond.push_back(Point3f(0.0, 22.5, 6.0));

	


	int cdc = aff.type();
	int ret = cv::estimateAffine3D(mfirst, msecond, aff, inliers);
	std::cout << mfirst << std::endl;
	std::cout << msecond << std::endl;
	std::cout << aff << std::endl;
	int cc = aff.type();
	DataTransform = aff;
	FileStorage fs1("estimateAffine3D.xml", FileStorage::WRITE);
	fs1 << "DataTransform" << DataTransform;
	fs1 << "inliers" << inliers;
	fs1.release();


	//////////////////////




	SafeRelease(pColorSource);
	SafeRelease(pDepthSource);
	SafeRelease(pColorReader);
	SafeRelease(pDepthReader);
	SafeRelease(pColorDescription);
	SafeRelease(pDepthDescription);
	SafeRelease(pCoordinateMapper);
	if (pSensor){
		pSensor->Close();
	}
	SafeRelease(pSensor);
	cv::destroyAllWindows();
	return 0;
}
bool Kincet3D::FindFiducials(const Mat& inputmat)
{
	bool first_time = true;
	NewCalculate = true;
	int widowsize = 30;

	double outputfiducial[16];

	Mat src = inputmat;
	if (src.empty())
		return false;
	if (first_time)
	{
		gs_src = Mat(src.size(), CV_8UC1);
		bw_src = Mat(src.size(), CV_8UC1);
		edges = Mat(src.size(), CV_8UC1);
		first_time = false;
	}
	cvtColor(src, gray, COLOR_BGR2GRAY);

	//=======================tracking========================
	if (NewCalculate == false)
	{
		//			widowsize=0;
		for (int m = 0; m < 8; m++)
		{
			Point2f potentianltest[8];
			double x = fiducial[m].x;
			double y = fiducial[m].y;
			double cols = gray.cols;
			double rows = gray.rows;
			if (x - widowsize>0 && y + widowsize<cols && y - widowsize>0 && y + widowsize<rows)
			{
				Rect roi(x - widowsize, y - widowsize, widowsize * 2, widowsize * 2);
				Mat subMatrix(gray, roi);
				Mat subMatrixshow(gray, roi);
				//					cv::imshow("submatrix", subMatrixshow);
				//					cv::waitKey(10);
				if (find_glyphs(subMatrix, potentianltest))
				{
					fiducial[m].x = potentianltest[m].x + x - widowsize;
					fiducial[m].y = potentianltest[m].y + y - widowsize;
				}
				else
				{
					NewCalculate = true;
					break;
				}
			}
		}

		for (int i = 0; i<8; i++)
		{
			outputfiducial[i * 2] = (double)fiducial[i].x;
			outputfiducial[i * 2 + 1] = (double)fiducial[i].y;
		}
		return true;
	}
	//================================================first time========================
	if (NewCalculate == true)
	{
		if (find_glyphs(gray, fiducial))
		{
			NewCalculate = false;
			for (int i = 0; i<8; i++)
			{
				outputfiducial[i * 2] = (double)fiducial[i].x;
				outputfiducial[i * 2 + 1] = (double)fiducial[i].y;
			}
			return true;
		}
		else
			return false;

	}
}

bool Kincet3D::find_glyphs(const Mat &img, Point2f *glyph_center)
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

	int pattern5[] = { 0, 1, 0, 1, 0, 1, 1, 1, 1,
		0, 1, 1, 1, 0, 1, 0, 1, 1,
		1, 1, 1, 1, 0, 1, 0, 1, 0,
		1, 1, 0, 1, 0, 1, 1, 1, 0 };

	int pattern6[] = { 0, 1, 0, 1, 1, 0, 0, 1, 0,
		0, 0, 0, 1, 1, 1, 0, 1, 0,
		0, 1, 0, 0, 1, 1, 0, 1, 0,
		0, 1, 0, 1, 1, 1, 0, 0, 0 };



	int pattern7[] = { 1, 0, 1, 1, 0, 0, 1, 1, 1,
		1, 0, 1, 0, 0, 1, 1, 1, 1,
		1, 1, 1, 0, 0, 1, 1, 0, 1,
		1, 1, 1, 1, 0, 0, 1, 0, 1 };

	int pattern8[] = { 1, 1, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 0, 1, 0, 1, 1, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 1, 1,
		0, 0, 1, 1, 0, 1, 0, 0, 0 };


	const int *pattern[] = { pattern1, pattern2, pattern3, pattern4, pattern5, pattern6, pattern7, pattern8 };

	const int APPROX_POLY_EPSILON = 3;
	const int MIN_CONTOUR_AREA = 50;
	const int CELL_NUM_ROW = 3;
	const int GLYPH_SIZE = 30;
	const int CELL_NUM = CELL_NUM_ROW * CELL_NUM_ROW;
	const int CELL_SIZE = GLYPH_SIZE / CELL_NUM_ROW;
	const int MAX_DELAY_FRAME = 3;

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

	double high_threshold = threshold(gs_src, bw_src, 0, 255, THRESH_BINARY | THRESH_OTSU);
	Canny(gs_src, edges, 0.5*high_threshold, high_threshold);
	vector<vector<Point> > contours;
	findContours(edges, contours, RETR_LIST, CHAIN_APPROX_SIMPLE);

	vector<vector<Point> > vertices;
	vertices.clear();
	for (int i = 0; i<contours.size(); i++)
	{


		vector<Point> cnt = contours.at(i);
		vector<Point> v;
		approxPolyDP(cnt, v, APPROX_POLY_EPSILON, true);
		if (v.size() == 4 && contourArea(v)>MIN_CONTOUR_AREA && isContourConvex(v))
		{
			vertices.push_back(v);
			polylines(src, contours[i], true, cv::Scalar(155, 255, 55));
		}
	}
	imshow("plylines", src);
	waitKey(20);

	bool submatrix = false;

	if (!vertices.empty() && NewCalculate == false)
	{
		submatrix = true;
	}


	else if (vertices.size() < 8)
		return false;
	Mat m;
	Size dsize(GLYPH_SIZE, GLYPH_SIZE);
	Mat glyph(dsize, CV_8UC1);
	int p[CELL_NUM + 1];
	int found[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	Point2f fv[4];

	for (int g = 0; g<vertices.size(); g++)
	{
		for (int i = 0; i<4; i++)
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
		for (int i = 0; i<CELL_NUM; i++)
		{
			int sum = 0;
			int row = (i / 3) * CELL_SIZE;
			int col = (i % 3) * CELL_SIZE;
			for (int r = row; r<row + CELL_SIZE; r++)
				for (int c = col; c<col + CELL_SIZE; c++)
					sum += (int)glyph.at<uchar>(r, c);
			if (sum / 255 > CELL_SIZE*CELL_SIZE / 2)
				p[i] = 1;
			else
				p[i] = 0;
		}

		for (int j = 4, i = 0; i<8 && j == 4; i++)
		{
			for (j = 0; j<4; j++)
			{
				int k;
				for (k = 0; k<CELL_NUM; k++)
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
	cout << "found glyphs" << "\n";
	if (submatrix && (found[0] + found[1] + found[2] + found[3] + found[4] + found[5] + found[6] + found[7]))
	{
		return true;
	}

	else if (found[0] && found[1] && found[2] && found[3] && found[4] && found[5] && found[6] && found[7])
	{
		return true;
	}
	else
		return false;

}

int Kincet3D::KinectRun()

{

	cv::setUseOptimized(true);

	// Sensor
	IKinectSensor* pSensor;
	HRESULT hResult = S_OK;
	hResult = GetDefaultKinectSensor(&pSensor);
	if (FAILED(hResult)){
		std::cerr << "Error : GetDefaultKinectSensor" << std::endl;
		return -1;
	}


	hResult = pSensor->Open();
	if (FAILED(hResult)){
		std::cerr << "Error : IKinectSensor::Open()" << std::endl;
		return -1;
	}


	IDepthFrameSource* pDepthSource;
	hResult = pSensor->get_DepthFrameSource(&pDepthSource);
	if (FAILED(hResult)){
		std::cerr << "Error : IKinectSensor::get_DepthFrameSource()" << std::endl;
		return -1;
	}

	// Reader
	IDepthFrameReader* pDepthReader;
	hResult = pDepthSource->OpenReader(&pDepthReader);
	if (FAILED(hResult)){
		std::cerr << "Error : IDepthFrameSource::OpenReader()" << std::endl;
		return -1;
	}

	//	/////////////////////////////////////////////////////////////this is the test of infrared image, we don't need it///////
	//
	//	ILongExposureInfraredFrameSource* pInfraredSource = NULL;
	//	hResult = pSensor->get_LongExposureInfraredFrameSource( &pInfraredSource );
	//	if( FAILED( hResult ) ){
	//		std::cerr << "Error : IKinectSensor::pInfraredSource()" << std::endl;
	//		return -1;
	//	}
	//
	//	// Reader
	//	ILongExposureInfraredFrameReader* pInfraredReader;
	//	hResult = pInfraredSource->OpenReader( &pInfraredReader );
	//	if( FAILED( hResult ) ){
	//		std::cerr << "Error : IDepthFrameSource::OpenReader()" << std::endl;
	//		return -1;
	//	}
	//
	//	IFrameDescription* pInfraredDescription;
	//	hResult = pInfraredSource->get_FrameDescription( &pInfraredDescription );
	//	if( FAILED( hResult ) ){
	//		std::cerr << "Error : IDepthFrameSource::get_FrameDescription()" << std::endl;
	//		return -1;
	//	}
	//
	//	int mInfraredWidth = 0;
	//	int mInfraredHeight = 0;
	//	pInfraredDescription->get_Width( &mInfraredWidth ); // 512
	//	pInfraredDescription->get_Height( &mInfraredHeight ); // 424
	//	
	//	unsigned int mInfraredBufferSize;// = mInfraredWidth * mInfraredHeight * sizeof( unsigned short );
	////	pInfraredDescription->get_BytesPerPixel(&mInfraredBufferSize);
	//	cv::Mat mInfraredBufferMat(mInfraredHeight,mInfraredWidth, CV_16SC1 );
	//	Mat mInfraredMat;
	//	mInfraredMat=Mat(mInfraredHeight,mInfraredWidth, CV_8UC1 );
	//	cv::namedWindow( "mInfrared" );
	//	UINT16 *infraData=NULL;
	//	UINT number;
	//	byte infraMatData[512 * 424];
	//	
	//	UINT16 *pBuffer = NULL;
	//	while (1)
	//	{
	//		ILongExposureInfraredFrame* pInfraredFrame = NULL;
	//
	//		HRESULT hr = pInfraredReader->AcquireLatestFrame(&pInfraredFrame);
	//		if( SUCCEEDED( hr ) )
	//		{
	//			pInfraredFrame->CopyFrameDataToArray(512 * 424, infraData);
	//			pInfraredFrame->AccessUnderlyingBuffer(&mInfraredBufferSize, reinterpret_cast<UINT16**>( &mInfraredBufferMat.data ) );
	//			pInfraredFrame->AccessUnderlyingBuffer(&mInfraredBufferSize, &pBuffer );
	////			pInfraredFrame->AccessUnderlyingBuffer(&mInfraredBufferSize, &infraData);
	//	
	//
	//			for (int i = 0; i < 512 * 424; ++i)
	//			{
	//				// Get infrared value
	//				UINT16 ir = *pBuffer;
	//				if (ir>0)
	//					int aa=33;
	//
	//				// Bitshift
	//				byte intensity = (byte)(ir >> 8);
	//				if (intensity!=0)
	//					int aa=33;
	//				// Assign infrared intensity
	//				infraMatData[i] = intensity;
	//				pBuffer++;
	//
	//			}
	//
	//	//		memcpy( mInfraredMat.data ,infraMatData,mInfraredBufferSize );
	//			{
	//		//		mInfraredBufferMat.convertTo( mInfraredMat, CV_8U, -255.0f / 4500.0f, 255.0f );	
	//			}
	//		}
	//		cv::imshow( "mInfrared", mInfraredMat );
	//
	//		if( cv::waitKey( 1 ) == VK_ESCAPE )
	//		{
	//			break;
	//		}
	//		SafeRelease( pInfraredFrame );
	//	}
	//
	//	SafeRelease( pInfraredSource );
	//	SafeRelease( pInfraredDescription );
	//	SafeRelease( pInfraredReader );
	/////////////////////////////////////////////////////////////////////

	// Description

	IFrameDescription* pDepthDescription;
	hResult = pDepthSource->get_FrameDescription(&pDepthDescription);
	if (FAILED(hResult)){
		std::cerr << "Error : IDepthFrameSource::get_FrameDescription()" << std::endl;
		return -1;
	}

	int depthWidth = 0;
	int depthHeight = 0;
	pDepthDescription->get_Width(&depthWidth); // 512
	pDepthDescription->get_Height(&depthHeight); // 424
	unsigned int depthBufferSize = depthWidth * depthHeight * sizeof(unsigned short);

	cv::Mat depthBufferMat(depthHeight, depthWidth, CV_16SC1);
	depthMat = Mat(depthHeight, depthWidth, CV_8UC1);
	cv::namedWindow("Depth");

	// Coordinate Mapper
	ICoordinateMapper* pCoordinateMapper;
	hResult = pSensor->get_CoordinateMapper(&pCoordinateMapper);
	if (FAILED(hResult)){
		std::cerr << "Error : IKinectSensor::get_CoordinateMapper()" << std::endl;
		return -1;
	}

	unsigned short minDepth, maxDepth;
	pDepthSource->get_DepthMinReliableDistance(&minDepth);
	pDepthSource->get_DepthMaxReliableDistance(&maxDepth);

	FileStorage fsmy2("estimateAffine3D.xml", FileStorage::READ);
	fsmy2["DataTransform"] >> DataTransform;
	fsmy2.release();

	double a1 = DataTransform.at<double>(0, 0);
	double a2 = DataTransform.at<double>(0, 1);
	double a3 = DataTransform.at<double>(0, 2);
	double a4 = DataTransform.at<double>(0, 3);

	double b1 = DataTransform.at<double>(1, 0);
	double b2 = DataTransform.at<double>(1, 1);
	double b3 = DataTransform.at<double>(1, 2);
	double b4 = DataTransform.at<double>(1, 3);

	double c1 = DataTransform.at<double>(2, 0);
	double c2 = DataTransform.at<double>(2, 1);
	double c3 = DataTransform.at<double>(2, 2);
	double c4 = DataTransform.at<double>(2, 3);


	vector<CameraSpacePoint> referencePoint(512 * 424);

	DWORD time_start, time_end;
	Point3f F = Point3f(0, 0, 0);
	Target3D.reserve(600 * 500);

	int times = 0;
	while (1)
	{
		//printf("kinect thread\n"); 
		if (times>4)
		{
			askData = false;
			getWrong = true;
			times = 0;
		}

		if (times == 0)
		{
			time_start = GetTickCount();

		}
		// Depth Frame
		IDepthFrame* pDepthFrame = nullptr;
		hResult = pDepthReader->AcquireLatestFrame(&pDepthFrame);
		if (SUCCEEDED(hResult))
		{
			hResult = pDepthFrame->AccessUnderlyingBuffer(&depthBufferSize, reinterpret_cast<UINT16**>(&depthBufferMat.data));
			if (SUCCEEDED(hResult))
			{
				//printf("1----\n");
				depthBufferMat.convertTo(depthMat, CV_8U, -255.0f / 4500.0f, 255.0f);

				if (askData)
				{
					//printf("asked data in kinect thread\n");
					printf("2----\n");
					hResult = pCoordinateMapper->MapDepthFrameToCameraSpace(512 * 424, reinterpret_cast<UINT16*>(depthBufferMat.data), 512 * 424, &referencePoint[0]);
					if (SUCCEEDED(hResult))
					{
						printf("3----\n");
						Target3D.clear(); 
						ofstream myfile;
						myfile.open("print2.txt");
						
							

						
						
						for (int i = minRow; i<maxRow; i++)
						{
							for (int j = minCol; j<maxCol; j++)
							{
								int aa = i * 512 + j;
								float x = referencePoint[aa].X;
								float y = referencePoint[aa].Y;
								float z = referencePoint[aa].Z;
								if (z>0 && z<10)
								{
									//cout << x << " " << y << " " << z << " " << a1 << " " << a2 << " " << a3 << " " << a4 << " " << b1 << " " << b2 << " " << b3 << " " << b4 << " " << c1 << " " << c2 << " " << c3 << " " << c4;

									F.x = x * 100 * a1 + y * 100 * a2 + z * 100 * a3 + 1 * a4;
									int xx = int(F.x);
									F.y = x * 100 * b1 + y * 100 * b2 + z * 100 * b3 + 1 * b4;
									int yy = int(F.y);
									F.z = x * 100 * c1 + y * 100 * c2 + z * 100 * c3 + 1 * c4;
									int zz = int(F.z);
									if (F.z < 2.3)
										continue;
									Target3D.push_back(F);
									myfile << F.x << "," << F.y <<"," <<F.z<<"\n";
									//cout <<"previous size=" <<Target3D.size()<<"\n"; 

								}
							}
						}
						myfile.close();

						getWrong = false;
						dataReady = true;
						printf("data is ready\n");
						askData = false;
						times = 0;
					}
					else times++;
				}
				time_end = GetTickCount();
				double fps = 1000.0 / (double)(time_end - time_start);
				cv::putText(depthMat, to_string((float)fps), Point(20, 20), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(150.));

				cv::imshow("Depth", depthMat);

				if (cv::waitKey(1) == VK_ESCAPE)
				{
					break;
				}

			}
			else
			{
				if (askData)
					times++;
			}
		}
		else
		{
			if (askData)
				times++;
		}
		SafeRelease(pDepthFrame);
	}





	SafeRelease(pDepthSource);
	SafeRelease(pDepthReader);
	SafeRelease(pDepthDescription);
	SafeRelease(pCoordinateMapper);


	if (pSensor)
	{
		pSensor->Close();
	}
	SafeRelease(pSensor);
	cv::destroyAllWindows();
	return 0;
}

int Kincet3D::KinectInitial()
{
	//cv::setUseOptimized( true );

	//// Sensor
	//IKinectSensor* pSensor;
	//HRESULT hResult = S_OK;
	//hResult = GetDefaultKinectSensor( &pSensor );
	//if( FAILED( hResult ) ){
	//	std::cerr << "Error : GetDefaultKinectSensor" << std::endl;
	//	return -1;
	//}

	//hResult = pSensor->Open();
	//if( FAILED( hResult ) ){
	//	std::cerr << "Error : IKinectSensor::Open()" << std::endl;
	//	return -1;
	//}


	//IDepthFrameSource* pDepthSource;
	//hResult = pSensor->get_DepthFrameSource( &pDepthSource );
	//if( FAILED( hResult ) ){
	//	std::cerr << "Error : IKinectSensor::get_DepthFrameSource()" << std::endl;
	//	return -1;
	//}

	//// Reader
	//IDepthFrameReader* pDepthReader;
	//hResult = pDepthSource->OpenReader( &pDepthReader );
	//if( FAILED( hResult ) ){
	//	std::cerr << "Error : IDepthFrameSource::OpenReader()" << std::endl;
	//	return -1;
	//}

	//// Description

	//IFrameDescription* pDepthDescription;
	//hResult = pDepthSource->get_FrameDescription( &pDepthDescription );
	//if( FAILED( hResult ) ){
	//	std::cerr << "Error : IDepthFrameSource::get_FrameDescription()" << std::endl;
	//	return -1;
	//}

	//int depthWidth = 0;
	//int depthHeight = 0;
	//pDepthDescription->get_Width( &depthWidth ); // 512
	//pDepthDescription->get_Height( &depthHeight ); // 424
	//unsigned int depthBufferSize = depthWidth * depthHeight * sizeof( unsigned short );

	//cv::Mat depthBufferMat( depthHeight, depthWidth, CV_16SC1 );
	//cv::Mat depthMat( depthHeight, depthWidth, CV_8UC1 );
	//cv::namedWindow( "Depth" );

	//// Coordinate Mapper
	//ICoordinateMapper* pCoordinateMapper;
	//hResult = pSensor->get_CoordinateMapper( &pCoordinateMapper );
	//if( FAILED( hResult ) ){
	//	std::cerr << "Error : IKinectSensor::get_CoordinateMapper()" << std::endl;
	//	return -1;
	//}

	//unsigned short minDepth, maxDepth;
	//pDepthSource->get_DepthMinReliableDistance( &minDepth );
	//pDepthSource->get_DepthMaxReliableDistance( &maxDepth );

	//FileStorage fsmy2("estimateAffine3D.xml", FileStorage::READ);
	//fsmy2 ["DataTransform"]>>DataTransform;
	//fsmy2.release();

	//double a1=DataTransform.at<double>(0,0);
	//double a2=DataTransform.at<double>(0,1);
	//double a3=DataTransform.at<double>(0,2);
	//double a4=DataTransform.at<double>(0,3);

	//double b1=DataTransform.at<double>(1,0);
	//double b2=DataTransform.at<double>(1,1);
	//double b3=DataTransform.at<double>(1,2);
	//double b4=DataTransform.at<double>(1,3);

	//double c1=DataTransform.at<double>(2,0);
	//double c2=DataTransform.at<double>(2,1);
	//double c3=DataTransform.at<double>(2,2);
	//double c4=DataTransform.at<double>(2,3);


	//vector<CameraSpacePoint> referencePoint(512 * 424);

	//DWORD time_start,time_end;
	//Point3f F=Point3f(0,0,0);
	//Target3D.reserve(600 * 500);
	//while( 1 )
	//{

	//	Target3D.clear();
	//	time_start = GetTickCount();
	//	// Depth Frame
	//	IDepthFrame* pDepthFrame = nullptr;
	//	hResult = pDepthReader->AcquireLatestFrame( &pDepthFrame );
	//	if( SUCCEEDED( hResult ) )
	//	{
	//		hResult = pDepthFrame->AccessUnderlyingBuffer( &depthBufferSize, reinterpret_cast<UINT16**>( &depthBufferMat.data ) );
	//		if( SUCCEEDED( hResult ) )
	//		{
	//			depthBufferMat.convertTo( depthMat, CV_8U, -255.0f / 4500.0f, 255.0f );
	//			hResult=pCoordinateMapper->MapDepthFrameToCameraSpace(512 * 424 , reinterpret_cast<UINT16*>( depthBufferMat.data ),512 * 424,&referencePoint[0]);
	//			if (SUCCEEDED(hResult))
	//			{
	//				if ( saveOriginalData)
	//				{
	//					ofstream fout;
	//					fout.open("OriginalCloudPoints.txt");
	//					for (int i=0;i<512 * 424;i++)
	//					{						
	//						if (referencePoint[i].Z>0 && referencePoint[i].Z<10)
	//						{
	//							fout << referencePoint[i].X*100<<' '<< referencePoint[i].Y*100<<' '<<referencePoint[i].Z*100<<endl; 
	//						}



	//					}
	//					fout << flush; 
	//					fout.close();

	//				}
	//				////////////////////////////////////

	//				dataReady=false;
	//				if(1)
	//				{


	//					for (int i=minRow;i<maxRow;i++)
	//					{
	//						for (int j=minCol;j<maxCol;j++)

	//						{						

	//							int aa=i*512+j;
	//							double x=referencePoint[aa].X;
	//							double y=referencePoint[aa].Y;
	//							double z=referencePoint[aa].Z;
	//							if (z>0 && z<10)
	//							{
	//								F.x=x*100*a1+y*100*a2+z*100*a3+4*a4;
	//								F.y=x*100*b1+y*100*b2+z*100*b3+4*b4;
	//								F.z=x*100*c1+y*100*c2+z*100*c3+4*c4;
	//								Target3D.push_back(F);
	//							}
	//						}
	//					}
	//					dataReady=true;

	//					/////////////////////////////////////////

	//					//Mat Target3DH;
	//					//cv::convertPointsToHomogeneous(Target3D,Target3DH);

	//					////				int cc=Target3DH.channels();
	//					////				int ccc=Target3DH.cols;
	//					////				int cccc=Target3DH.rows;

	//					//Mat matData=Target3DH.reshape(1).clone();

	//					////	std::cout << Target3DH << std::endl ;

	//					////			int cols=DataTransform.type();
	//					////			int rows=matData.type();
	//					//Mat T;
	//					//DataTransform.convertTo(T,5);
	//					////			int tcols=T.type();
	//					//cloudPoints=T*matData.t();
	//					//cloudPoints=cloudPoints.t();
	//					////std::cout << Model << std::endl ;
	//				}


	//				if (isSaveCloudPoints)
	//				{
	//					ofstream fout;
	//					fout.open("cloudPoints.txt");
	//					for (int i=0;i<Target3D.size();i++)
	//					{
	//						fout << Target3D[i].x<<' '<< Target3D[i].y<<' '<<Target3D[i].z<<endl; 
	//					}
	//					fout << flush; 
	//					fout.close();
	//				}

	//			}
	//			time_end=GetTickCount();
	//			double fps = 1000.0/(double)(time_end - time_start);	
	//			putText(depthMat,to_string((float) fps),Point(20,20),FONT_HERSHEY_SIMPLEX,0.5,Scalar( 150. ));

	//			cv::imshow( "Depth", depthMat );



	//			if( cv::waitKey( 1 ) == VK_ESCAPE )
	//			{
	//				break;
	//			}			


	//		}


	//		///////////////////////////////////
	//	}



	//	SafeRelease( pDepthFrame );


	//}


	//SafeRelease( pDepthSource );
	//SafeRelease( pDepthReader );
	//SafeRelease( pDepthDescription );
	//SafeRelease( pCoordinateMapper );
	//if( pSensor ){
	//	pSensor->Close();
	//}
	//SafeRelease( pSensor );
	//cv::destroyAllWindows();
	return 0;

}

void Kincet3D::thetransform()
{

}