
#include"Fiducial.h"

using namespace std;
using namespace cv;

#ifdef __cplusplus
extern "C"
#endif
{


bool first_time = true;
extern bool NewCalculate= true;




//const int MAX_DELAY_FRAME     = 3;

//Mat rgb_dst;

	//Mat gs_src;
	//Mat bw_src;
	//Mat edges;
	//Mat src;

	//jdouble outputfiducial[8];


Mat gray;
bool Fiducial::getPoints(Mat &src,double *outputfiducial)
//int main(int argc, char** argv)
{
//	memcopy(&targetModel2[0],receved2,1000);

	int widowsize=30;
//	string filename = "rectange.mp4";
//    VideoCapture capture(filename);

//	cv::namedWindow("output");
//	cv::namedWindow("submatrix");
//    if( !capture.isOpened() )
//        throw "Error when reading steam_avi";
//	for(; ; )
    {
 //       capture >> src;
        if(src.empty())
            return false;
		if (first_time)
		{
			//gs_src = Mat(src.size(), CV_8UC1);
			//bw_src = Mat(src.size(), CV_8UC1);
			//edges  = Mat(src.size(), CV_8UC1);
			first_time = false;
		}
		cvtColor(src, gray, COLOR_BGR2GRAY);

		//=======================tracking========================
		if (NewCalculate==false )
	    {
//			widowsize=0;
           for (int m = 0; m < 4; m++)
           {
			   Point2f potentianltest[4];
        	   double x= fiducial[m].x;
        	   double y= fiducial[m].y;
        	   double cols= gray.cols;
        	   double rows= gray.rows;
        	   if (x-widowsize>0 & y+widowsize<cols & y-widowsize>0 & y+widowsize<rows )
        	   {
					Rect roi( x-widowsize, y-widowsize, widowsize*2, widowsize*2 );
					Mat subMatrix(gray,roi);
					Mat subMatrixshow(gray,roi);
//					cv::imshow("submatrix", subMatrixshow);
//					cv::waitKey(10);
					if(find_glyphs(subMatrix, potentianltest))
					{
						fiducial[m].x=potentianltest[m].x+x-widowsize;
						fiducial[m].y=potentianltest[m].y+y-widowsize;
					}
					else
					{
						NewCalculate=true;
						break;
					}
			    }
		    }

           for(int i=0;i<4;i++)
           {
        	   outputfiducial[i*2  ] = (double)fiducial[i].x;
        	   outputfiducial[i*2+1] = (double)fiducial[i].y;
           }
           return true;


//			for (int i=0;i<4;i++)
//			{
//				circle(src, Point(fiducial[i]), (i+1)*5,  cv::Scalar(255,255,0));
//			}
//		cv::imshow("output", src);
//		waitKey(5);
		 }
	//================================================first time========================
		if (NewCalculate==true )
		{
			if(find_glyphs(gray, fiducial))
			{
				NewCalculate=false;
		           for(int i=0;i<4;i++)
		           {
		        	   outputfiducial[i*2  ] = (double)fiducial[i].x;
		        	   outputfiducial[i*2+1] = (double)fiducial[i].y;
		           }
		           return true;

//				for (int i=0;i<4;i++)
//				{
//					circle(src, fiducial[i], (i+1)*5,  cv::Scalar(255,255,0));
//				}

			}
//		cv::imshow("output", src);
//		waitKey(5);

		}
	}

}

bool Fiducial::find_glyphs(const Mat &img, Point2f *glyph_center) {



	const int APPROX_POLY_EPSILON = 3;
	const int MIN_CONTOUR_AREA    = 50;
	const int CELL_NUM_ROW        = 3;
	const int GLYPH_SIZE          = 30;
	const int CELL_NUM            = CELL_NUM_ROW * CELL_NUM_ROW;
	const int CELL_SIZE           = GLYPH_SIZE / CELL_NUM_ROW;





	int pattern1[] = {/*0, 1, 0, 1, 1, 1, 1, 0, 1,*/378,
						 /* 1, 1, 0, 0, 1, 1, 1, 1, 0,*/243,
						 /* 1, 0, 1, 1, 1, 1, 0, 1, 0,*/189,
						 /* 0, 1, 1, 1, 1, 0, 0, 1, 1};*/414};

		int pattern2[] = {/*0, 0, 1, 1, 0, 1, 0, 1, 0,*/172,
						  /*0, 1, 0, 1, 0, 0, 0, 1, 1,*/394,
						  /*0, 1, 0, 1, 0, 1,	1, 0, 0,*/106,
						  /*1, 1, 0, 0, 0, 1, 0, 1, 0};*/163};

		int pattern3[] = {/*0, 0, 1, 0, 1, 0, 1, 0, 1,*/340,
						  /*1, 0, 0, 0, 1, 0, 1, 0, 1,*/337,
						  /*1, 0, 1, 0, 1, 0, 1, 0, 0,*/85,
						  /*1, 0, 1, 0, 1, 0, 0, 0, 1};*/277};

		int pattern4[] = {/*0, 1, 0, 0, 1, 1, 0, 0, 1,*/306,
						  /*0, 0, 0, 0, 1, 1, 1, 1, 0,*/240,
						  /*1, 0, 0, 1, 1, 0, 0, 1, 0,*/153,
						  /*0, 1, 1, 1, 1, 0, 0, 0, 0}*/30};

	const int *pattern[] = {pattern1, pattern2, pattern3, pattern4};



	Mat gs_src (img);
	Mat bw_src;
	Mat edges;

	Point2f cw_dst_points[4];
	cw_dst_points[0] = Point2f(-CELL_SIZE,			 -CELL_SIZE);
	cw_dst_points[1] = Point2f(GLYPH_SIZE + CELL_SIZE, -CELL_SIZE);
	cw_dst_points[2] = Point2f(GLYPH_SIZE + CELL_SIZE, GLYPH_SIZE + CELL_SIZE);
	cw_dst_points[3] = Point2f(-CELL_SIZE,             GLYPH_SIZE + CELL_SIZE);

	Point2f ccw_dst_points[4];
	ccw_dst_points[0] = Point2f(-CELL_SIZE,             -CELL_SIZE);
	ccw_dst_points[1] = Point2f(-CELL_SIZE,             GLYPH_SIZE + CELL_SIZE);
	ccw_dst_points[2] = Point2f(GLYPH_SIZE + CELL_SIZE, GLYPH_SIZE + CELL_SIZE);
	ccw_dst_points[3] = Point2f(GLYPH_SIZE + CELL_SIZE, -CELL_SIZE);

	double high_threshold = threshold(gs_src, bw_src, 0, 255, THRESH_BINARY | THRESH_OTSU);
	Canny(gs_src, edges, 0.5*high_threshold, high_threshold);
	vector<vector<Point> > contours;
	findContours(edges, contours, RETR_LIST, CHAIN_APPROX_SIMPLE);

	vector<vector<Point> > vertices;
	vertices.clear();
	for(int i=0;i<contours.size();i++)
	{


		vector<Point> cnt = contours.at(i);
		vector<Point> v;
		approxPolyDP(cnt, v, APPROX_POLY_EPSILON, true);
		if(v.size()==4 && contourArea(v)>MIN_CONTOUR_AREA && isContourConvex(v))
		{
			vertices.push_back(v);
///			polylines(src,contours[i], true, cv::Scalar(255,255,0));
		}

	}

//	cv::imshow("submatrix", src);
//	cv::waitKey(10);

	bool submatrix=false;

	if (!vertices.empty() && NewCalculate==false)
	{
		submatrix=true;
	}


	else if (vertices.size() < 4)
		return false;
	Mat m;
	Size dsize(GLYPH_SIZE, GLYPH_SIZE);
	Mat glyph(dsize, CV_8UC1);
	int p[CELL_NUM+1];
	int found[4] = {0, 0, 0, 0};
	Point2f fv[4];

	for(int g=0;g<vertices.size();g++)
	{
		for (int i=0; i<4; i++)
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

		int candidatePattern=0;
		for (int i=0; i<CELL_NUM; i++)
		{
			int sum = 0;
			int row = (i/3) * CELL_SIZE;
			int col = (i%3) * CELL_SIZE;
			for (int r=row; r<row+CELL_SIZE; r++)
				for (int c=col; c<col+CELL_SIZE; c++)
					sum += (int)glyph.at<uchar>(c, r);
			if (sum / 255 > CELL_SIZE*CELL_SIZE/2)
				p[i] = 1;
			else
				p[i] = 0;
			candidatePattern+=p[i]*pow(2,i);

		}

		int checkEqual;
		for (int j=4, i=0; i<4 && j==4; i++)
		{
			for (j=0; j<4; j++)
			{

				checkEqual=pattern[i][j] ^ candidatePattern;
				if (checkEqual==0)
				{
					Point2f center;
					center.x = (fv[0].x+fv[1].x+fv[2].x+fv[3].x)*0.25;
					center.y = (fv[0].y+fv[1].y+fv[2].y+fv[3].y)*0.25;

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
	if (submatrix && (found[0] + found[1] + found[2] +found[3]))
	{
		return true;
	}

	else if (found[0] && found[1] && found[2] && found[3])
	{
		return true;
	}
	else
		return false;
}
/*
JNIEXPORT void JNICALL Java_com_pcvlab_glasslocalization_MainsActivity_setSrcPoints(JNIEnv*, jobject, jint jx, jint jy, jint ji)
{
	int i = int(ji);
	int x = int(jx);
	int y = int(jy);
//	src_points[i] = Point2f(x, y);
}
*/


#ifdef __cplusplus
}
#endif
