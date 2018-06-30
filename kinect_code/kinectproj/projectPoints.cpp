
#include"projectPoints.h"

using namespace std;
using namespace cv;

#ifdef __cplusplus
extern "C"
#endif
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

	const int *pattern[] = { pattern1, pattern2, pattern3, pattern4 };

	const int APPROX_POLY_EPSILON = 3;
	const int MIN_CONTOUR_AREA = 50;
	const int CELL_NUM_ROW = 3;
	const int GLYPH_SIZE = 30;
	const int CELL_NUM = CELL_NUM_ROW * CELL_NUM_ROW;
	const int CELL_SIZE = GLYPH_SIZE / CELL_NUM_ROW;
	const int MAX_DELAY_FRAME = 3;

	Mat rgb_dst;
	Mat gs_src;
	Mat bw_src;
	Mat edges;
	bool first_time = true;

	bool find_glyphs(const Mat &src, Point2f *glyph_center);
	//JNIEXPORT void JNICALL Java_com_pcvlab_glasslocalization_MainActivity_setSrcPoints(JNIEnv*, jobject, jint jx, jint jy, jint ji);

	struct InteriorOrientationParameters
	{
		double f;
		Point2f PP;
		double k1, k2, k3, p1, p2, p3;

		InteriorOrientationParameters(){};
		InteriorOrientationParameters(double f, double x, double y,
			double k1, double k2, double k3,
			double p1, double p2, double p3)
		{
			this->f = f;
			this->PP.x = x;
			this->PP.y = y;
			this->k1 = k1;
			this->k2 = k2;
			this->k3 = k3;
			this->p1 = p1;
			this->p2 = p2;
			this->p3 = p3;
		}
		/*
		void Set_IOP(double focallength, double xo, double yo,
		double k_1, double k_2, double k_3,
		double p_1, double p_2, double p_3)
		{
		f		= focallength;
		PP.x	= xo;
		PP.y	= yo;
		k1		= k_1;
		k2		= k_2;
		k3		= k_3;
		p1		= p_1;
		p2		= p_2;
		p3		= p_3;
		}
		*/
	};

	struct OrientationParameters
	{
		Point3f Point;
		double o;
		double p;
		double k;
		/*
		void Set_EOP(double X, double Y, double Z, double omega, double phi, double kappa)
		{
		Point.x = X;
		Point.y = Y;
		Point.z = Z;
		o = omega;
		p = phi;
		k = kappa;
		}
		*/
	};

	class CPhoto
	{
	public:
		CPhoto(){}
		virtual ~CPhoto(){};

	public:
		InteriorOrientationParameters	IOP;
		OrientationParameters			EOP;
		double	PixelSize;
		int		nRows, nCols;
		double 	half_nRows, half_nCols;

	protected:
		// Elements of the rotation matrix
		double m11, m12, m13;
		double m21, m22, m23;
		double m31, m32, m33;
		Mat R;

	public:
		void Set_PixelSize(double pixelsize)
		{
			PixelSize = pixelsize;
		}
		void Set_ImageSize(int nrows, int ncols)
		{
			nRows = nrows;
			nCols = ncols;
			half_nRows = nRows / 2;
			half_nCols = nCols / 2;
		}
		void Set_IOP(InteriorOrientationParameters iop)
		{
			IOP = iop;
		}
		void Set_IOP(double focallength, double xo, double yo,
			double k_1, double k_2, double k_3,
			double p_1, double p_2, double p_3)
		{
			IOP.f = focallength;
			IOP.PP.x = xo;
			IOP.PP.y = yo;
			IOP.k1 = k_1;
			IOP.k2 = k_2;
			IOP.k3 = k_3;
			IOP.p1 = p_1;
			IOP.p2 = p_2;
			IOP.p3 = p_3;
		}
		void Set_EOP(OrientationParameters eop)
		{
			EOP = eop;
		}
		void Set_EOP(double x, double y, double z, double omega, double phi, double kappa)
		{
			EOP.Point.x = x;
			EOP.Point.y = y;
			EOP.Point.z = z;
			EOP.o = omega;
			EOP.p = phi;
			EOP.k = kappa;
		}

		InteriorOrientationParameters	Get_IOP(){ return IOP; }
		OrientationParameters 			Get_EOP(){ return EOP; }

		Point2f Image2Photo(Point2f imagecoord)
		{
			Point2f PhotoCoord;

			PhotoCoord.x = ((double)imagecoord.x - half_nCols)*PixelSize;
			PhotoCoord.y = ((double)half_nRows - imagecoord.y)*PixelSize;

			return PhotoCoord;
		}
		Point2f Photo2Image(Point2f photocoord)
		{
			Point2f ImageCoord;

			ImageCoord.x = photocoord.x / PixelSize + half_nCols;
			ImageCoord.y = half_nRows - photocoord.y / PixelSize;

			return ImageCoord;
		}
		Point2f Photo2CalibratedPhoto(Point2f photocoord)
		{
			Point2f CPhotoCoord;

			double drx, dry;
			double ddx, ddy;
			double x_bar, y_bar;
			double xy_bar;
			double x_bar2, y_bar2;
			double r2, r4, r6;

			x_bar = photocoord.x - IOP.PP.x;
			y_bar = photocoord.y - IOP.PP.y;
			xy_bar = x_bar*y_bar;
			x_bar2 = x_bar*x_bar;
			y_bar2 = y_bar*y_bar;
			r2 = x_bar2 + y_bar2;
			r4 = r2*r2;
			r6 = r4*r2;

			// Radial lens distortion
			drx = x_bar*(IOP.k1*r2 + IOP.k2*r4 + IOP.k3*r6);
			dry = y_bar*(IOP.k1*r2 + IOP.k2*r4 + IOP.k3*r6);

			// Decentering lens distortion
			ddx = IOP.p1*(r2 + 2 * x_bar2) + 2 * IOP.p2*xy_bar;
			ddy = 2 * IOP.p1*xy_bar + IOP.p2*(r2 + 2 * y_bar2);

			CPhotoCoord.x = photocoord.x + (drx + ddx);
			CPhotoCoord.y = photocoord.y + (dry + ddy);

			return CPhotoCoord;
		}
		Point2f CalibratedPhoto2Photo(Point2f cphotocoord)
		{
			CvMat* A = cvCreateMat(2, 2, CV_64FC1);
			CvMat* AT = cvCreateMat(2, 2, CV_64FC1);
			CvMat* N = cvCreateMat(2, 2, CV_64FC1);
			CvMat* Ninv = cvCreateMat(2, 2, CV_64FC1);
			CvMat* csi_hat = cvCreateMat(2, 1, CV_64FC1);
			CvMat* y = cvCreateMat(2, 1, CV_64FC1);
			CvMat* c = cvCreateMat(2, 1, CV_64FC1);
			CvMat* e_tilda = cvCreateMat(2, 1, CV_64FC1);
			CvMat* e_tilda_T = cvCreateMat(1, 2, CV_64FC1);
			CvMat* eTe = cvCreateMat(1, 1, CV_64FC1);

			double xo, yo;	// initial values
			double x_bar, y_bar;
			double xy_bar;
			double x_bar2, y_bar2;
			double r2, r4, r6;
			double K1, K2;
			double drx, dry;
			double ddx, ddy;

			double a00, a01, a11, y0, y1;

			int iteration_stop = 0;
			int iteration_count = 0;

			double sigma_o_2_hat_change;
			double sigma_o_2_hat;
			double sigma_o_2_hat_old = 9999999999;

			xo = cphotocoord.x;
			yo = cphotocoord.y;

			while (iteration_stop == 0)
			{
				x_bar = xo - IOP.PP.x;
				y_bar = yo - IOP.PP.y;
				xy_bar = x_bar*y_bar;
				x_bar2 = x_bar*x_bar;
				y_bar2 = y_bar*y_bar;
				r2 = x_bar2 + y_bar2;
				r4 = r2*r2;
				r6 = r4*r2;

				K1 = IOP.k1*r2 + IOP.k2*r4 + IOP.k3*r6;
				K2 = IOP.k1 + 2 * IOP.k2*r2 + 3 * IOP.k3*r4;

				// Radial lens distortion
				drx = x_bar*K1;
				dry = y_bar*K1;

				// Decentering lens distortion
				ddx = IOP.p1*(r2 + 2 * x_bar2) + 2 * IOP.p2*xy_bar;
				ddy = 2 * IOP.p1*xy_bar + IOP.p2*(r2 + 2 * y_bar2);

				a00 = 1 + K1 + 2 * x_bar2*K2 + 6 * IOP.p1*x_bar + 2 * IOP.p2*y_bar;
				a01 = 2 * xy_bar*K2 + 2 * IOP.p1*y_bar + 2 * IOP.p2*x_bar;
				a11 = 1 + K1 + 2 * y_bar2*K2 + 2 * IOP.p1*x_bar + 6 * IOP.p2*y_bar;

				y0 = cphotocoord.x - (xo + drx + ddx);
				y1 = cphotocoord.y - (yo + dry + ddy);

				// Matrix
				cvmSet(A, 0, 0, a00);
				cvmSet(A, 0, 1, a01);
				cvmSet(A, 1, 0, a01);
				cvmSet(A, 1, 1, a11);

				cvmSet(y, 0, 0, y0);
				cvmSet(y, 1, 0, y1);

				// LESS
				cvTranspose(A, AT);
				cvMatMul(AT, A, N);
				cvMatMul(AT, y, c);
				cvInvert(N, Ninv);
				cvMatMul(Ninv, c, csi_hat);
				cvGEMM(A, csi_hat, 1.0, y, -1.0, e_tilda);
				cvTranspose(e_tilda, e_tilda_T);
				cvMatMul(e_tilda_T, e_tilda, eTe);

				sigma_o_2_hat = cvmGet(eTe, 0, 0) / 2;

				// Update
				xo += cvmGet(csi_hat, 0, 0);
				yo += cvmGet(csi_hat, 1, 0);

				iteration_count++;

				//Iteration Termination Conditions///////////////////////////////
				//Wolf, Ghilani, Adjustment Computations, pp. 243~245

				//1. Maximum iteration
				if (iteration_count == 10)
					iteration_stop = -1;

				//2. Maximum correction


				//3. Minitoring the Reference Variance
				sigma_o_2_hat_change = fabs((sigma_o_2_hat - sigma_o_2_hat_old) / sigma_o_2_hat_old);

				if (sigma_o_2_hat_change<0.01)
					iteration_stop = 1;

				sigma_o_2_hat_old = sigma_o_2_hat;
			}

			Point2f PhotoCoord;

			PhotoCoord.x = xo;
			PhotoCoord.y = yo;

			return PhotoCoord;
		}
		Point2f Image2CalibratedPhoto(Point2f imagecoord)
		{
			Point2f PhotoCoord;
			PhotoCoord = Photo2CalibratedPhoto(Image2Photo(imagecoord));
			return PhotoCoord;
		}
		Point2f CalibratedPhoto2Image(Point2f cphotocoord)
		{
			Point2f ImageCoord;
			ImageCoord = Photo2Image(CalibratedPhoto2Photo(cphotocoord));
			return ImageCoord;
		}
	};

	class CPhotogrammetry
	{
	public:
		CPhotogrammetry(){};
		virtual ~CPhotogrammetry(){};

	public:
		CPhoto Photo;

	protected:
		int n, m;
		int iteration_stop;
		int iteration_count;
		int iteration_tolerence;

		double sigma_o_2_hat_change_tolerence;
		double sigma_o_2_hat;
		double sigma_o_2_hat_change;
		double sigma_o_2_hat_old;

		double sin_o, sin_p, sin_k;
		double cos_o, cos_p, cos_k;

		// Elements of the rotation matrix
		double m11, m12, m13;
		double m21, m22, m23;
		double m31, m32, m33;

		// Parameters for linearization of the collinearity equation
		double q, r, s;
		double b11, b12, b13, b14, b15, b16, J;					//for x
		double b21, b22, b23, b24, b25, b26, K;					//for y

	protected:
		void Get_Elements_of_RotationMatrix(double o, double p, double k)
		{
			//For calculation efficiency
			sin_o = sin(o);
			sin_p = sin(p);
			sin_k = sin(k);
			cos_o = cos(o);
			cos_p = cos(p);
			cos_k = cos(k);

			//Calculation of elements of rotation matrix
			m11 = cos_p*cos_k;
			m12 = sin_o*sin_p*cos_k + cos_o*sin_k;
			m13 = -cos_o*sin_p*cos_k + sin_o*sin_k;

			m21 = -cos_p*sin_k;
			m22 = -sin_o*sin_p*sin_k + cos_o*cos_k;
			m23 = cos_o*sin_p*sin_k + sin_o*cos_k;

			m31 = sin_p;
			m32 = -sin_o*cos_p;
			m33 = cos_o*cos_p;
		}

		void Calculate_Parameters(InteriorOrientationParameters iop,
			OrientationParameters eop,
			Point2f photopoint,
			Point3f gcp)
		{
			Get_Elements_of_RotationMatrix(eop.o, eop.p, eop.k);

			Point3f Delta;
			double temp1, temp2, temp3, temp4, temp5;
			double fq;
			double fqq;

			Delta.x = gcp.x - eop.Point.x;
			Delta.y = gcp.y - eop.Point.y;
			Delta.z = gcp.z - eop.Point.z;

			q = m31*Delta.x + m32*Delta.y + m33*Delta.z;
			r = m11*Delta.x + m12*Delta.y + m13*Delta.z;
			s = m21*Delta.x + m22*Delta.y + m23*Delta.z;

			fq = iop.f / q;
			fqq = iop.f / (q*q);


			//Calculation derivatives//////////////////////////////////////////////////

			//For calculation efficiency
			temp1 = -m33*Delta.y + m32*Delta.z;
			temp2 = cos_p*Delta.x + sin_o*sin_p*Delta.y - cos_o*sin_p*Delta.z;
			temp3 = sin_o*cos_p*Delta.y;
			temp4 = cos_o*cos_p*Delta.z;
			temp5 = sin_p*Delta.x;

			//For x
			b11 = fqq*(r*(temp1)-q*(-m13*Delta.y + m12*Delta.z));				// dF/d(omega)
			b12 = fqq*(r*(temp2)-q*(-cos_k*temp5 + cos_k*temp3 - cos_k*temp4));  // dF/d(phi)
			b13 = -fq*(m21*Delta.x + m22*Delta.y + m23*Delta.z);					// dF/d(kappa)
			b14 = fqq*(r*m31 - q*m11);												// dF/d(XL)
			b15 = fqq*(r*m32 - q*m12);												// dF/d(YL)
			b16 = fqq*(r*m33 - q*m13);												// dF/d(ZL)

			J = photopoint.x - iop.PP.x + fq*r;

			//For y
			b21 = fqq*(s*(temp1)-q*(-m23*Delta.y + m22*Delta.z));				// dG/d(omega)
			b22 = fqq*(s*(temp2)-q*(sin_k*temp5 - sin_k*temp3 + sin_k*temp4));	// dG/d(phi)
			b23 = fq*(m11*Delta.x + m12*Delta.y + m13*Delta.z);					// dG/d(kappa)
			b24 = fqq*(s*m31 - q*m21);												// dG/d(XL)
			b25 = fqq*(s*m32 - q*m22);												// dG/d(YL)
			b26 = fqq*(s*m33 - q*m23);												// dG/d(ZL)

			K = photopoint.y - iop.PP.y + fq*s;
		}

	public:
		// Not Tested yet
		Point2f CollinearityEquations(Point3f point, InteriorOrientationParameters iop, OrientationParameters eop)
		{
			Point2f PhotoPoint;
			double r, s, q;

			Get_Elements_of_RotationMatrix(eop.o, eop.p, eop.k);

			r = m11*(point.x - eop.Point.x) + m12*(point.y - eop.Point.y) + m13*(point.z - eop.Point.z);
			s = m21*(point.x - eop.Point.x) + m22*(point.y - eop.Point.y) + m23*(point.z - eop.Point.z);
			q = m31*(point.x - eop.Point.x) + m32*(point.y - eop.Point.y) + m33*(point.z - eop.Point.z);

			PhotoPoint.x = iop.PP.x - iop.f*r / q;
			PhotoPoint.y = iop.PP.y - iop.f*s / q;

			return PhotoPoint;
		}
		int SpaceResection(CPhoto& photo, int ngcps, Point2f* photopoint, Point3f* gcp, double& sigma, double precision[], double residual[])
		{
			int i, j;
			int trow;

			//Number of unoknowns, knowns/////////////////////////////////

			//Number of Observations
			//(x,y)*number_of_points*number_of_photos
			n = 2 * ngcps;

			//Number of Unknown parameters
			m = 6;

			//For the iteration termination conditions////////////////////

			int iteration_stop = 0;
			int iteration_count = 0;

			double angle_tolerence = 0.000001;
			double position_tolerence = 0.000001;

			double change_sigma_o_2;
			double new_sigma_o_2;
			double old_sigma_o_2 = 9999999999;

			//Save Initial values////////////////////////////////////////
			OrientationParameters EOP;
			OrientationParameters EOP_initial;
			EOP_initial = photo.Get_EOP();

			//Matrix settings/////////////////////////////////////////////
			CvMat* A = cvCreateMat(n, m, CV_64FC1);
			CvMat* AT = cvCreateMat(m, n, CV_64FC1);
			CvMat* y = cvCreateMat(n, 1, CV_64FC1);
			CvMat* N = cvCreateMat(m, m, CV_64FC1);
			CvMat* Ninv = cvCreateMat(m, m, CV_64FC1);
			CvMat* c = cvCreateMat(m, 1, CV_64FC1);
			CvMat* csi_hat = cvCreateMat(m, 1, CV_64FC1);
			CvMat* e_tilda = cvCreateMat(n, 1, CV_64FC1);
			CvMat* e_tilda_T = cvCreateMat(1, n, CV_64FC1);
			CvMat* eTe = cvCreateMat(1, 1, CV_64FC1);

			for (i = 0; i<n; i++)
				for (j = 0; j<m; j++)
					cvmSet(A, i, j, 0.0);

			while (iteration_stop == 0)
			{
				for (i = 0; i<ngcps; i++)
				{
					Calculate_Parameters(photo.Get_IOP(), photo.Get_EOP(), photopoint[i], gcp[i]);
					trow = i * 2;

					//A matrix/////////////////////////////////////////////
					cvmSet(A, trow, 0, b11);	cvmSet(A, trow + 1, 0, b21);
					cvmSet(A, trow, 1, b12);	cvmSet(A, trow + 1, 1, b22);
					cvmSet(A, trow, 2, b13);	cvmSet(A, trow + 1, 2, b23);
					cvmSet(A, trow, 3, -b14);	cvmSet(A, trow + 1, 3, -b24);
					cvmSet(A, trow, 4, -b15);	cvmSet(A, trow + 1, 4, -b25);
					cvmSet(A, trow, 5, -b16);	cvmSet(A, trow + 1, 5, -b26);

					//y matrix////////////////////////////////////////////////////////////
					cvmSet(y, trow, 0, J);
					cvmSet(y, trow + 1, 0, K);

				}

				cvTranspose(A, AT);
				cvMatMul(AT, A, N);
				cvMatMul(AT, y, c);
				cvInvert(N, Ninv);
				cvMatMul(Ninv, c, csi_hat);
				cvGEMM(A, csi_hat, 1.0, y, -1.0, e_tilda);
				cvTranspose(e_tilda, e_tilda_T);
				cvMatMul(e_tilda_T, e_tilda, eTe);


				//Update////////////////////////////////////////////////////////////////////

				EOP = photo.Get_EOP();

				EOP.o += cvmGet(csi_hat, 0, 0);
				EOP.p += cvmGet(csi_hat, 1, 0);
				EOP.k += cvmGet(csi_hat, 2, 0);
				EOP.Point.x += cvmGet(csi_hat, 3, 0);
				EOP.Point.y += cvmGet(csi_hat, 4, 0);
				EOP.Point.z += cvmGet(csi_hat, 5, 0);

				photo.Set_EOP(EOP);

				iteration_count++;


				//Iteration Termination Conditions///////////////////////////////
				//Wolf, Ghilani, Adjustment Computations, pp. 243~245

				//1. Maximum iteration
				if (iteration_count == 10)
					iteration_stop = -1;

				//2. Maximum correction


				//3. Minitoring the Reference Variance
				new_sigma_o_2 = cvmGet(eTe, 0, 0) / (n - m);
				change_sigma_o_2 = fabs((new_sigma_o_2 - old_sigma_o_2) / old_sigma_o_2);

				if (change_sigma_o_2<0.01)
					iteration_stop = 1;

				old_sigma_o_2 = new_sigma_o_2;

				/*for (i = 0; i<ngcps; i++)
				{
					residual[i * 2] = cvmGet(eTe, i * 2, 0);
					residual[i * 2 + 1] = cvmGet(eTe, i * 2 + 1, 0);
				}*/
				EOP = photo.Get_EOP();
			}

			cvReleaseMat(&A);
			cvReleaseMat(&y);
			cvReleaseMat(&AT);
			cvReleaseMat(&N);
			cvReleaseMat(&Ninv);
			cvReleaseMat(&c);
			cvReleaseMat(&csi_hat);
			cvReleaseMat(&e_tilda);
			cvReleaseMat(&e_tilda_T);
			cvReleaseMat(&eTe);

			sigma = new_sigma_o_2;

			return iteration_stop;
		}
	};
	

	double dist_2P(Point2f A, Point2f B)
	{
		return sqrt((A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y));
	}

	double dist3D(Point3f A, Point3f B)
	{
		return sqrt((A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y) + (A.z - B.z)*(A.z - B.z));
	}

	void getPoints(double* fiducial, vector<double> &projpts, vector<Point3f> kinectpts)
	{
		

		
			

			CPhotogrammetry SR;
			//	SR.Photo.Set_IOP(2.0585e+000, 1.1326e-002, 3.1120e-001,-1.9197e-002, 1.8684e-002,-3.3957e-003, -4.2030e-004,  5.8225e-004, 0.0);
			SR.Photo.Set_IOP(2.0585e+000, 1.1326e-002, 3.1120e-001, 0, 0, 0, 0, 0, 0.0);

			SR.Photo.Set_PixelSize(0.005);
			SR.Photo.Set_ImageSize(576, 1024);

			double XL, YL, ZL;
			double omega, phi, kappa;


			int ngcps = 4;
			Point2f photopoint[4];
			Point2f imagepoint[4];
			Point3f gcp[4];
			double sigma;
			double precision[6];
			double residual[8];

			//	gcp[0].x = 1000;	gcp[0].y = 2000;	gcp[0].z = 3000;
			//	gcp[1].x = 1153;	gcp[1].y = 2000;	gcp[1].z = 3000;
			//	gcp[2].x = 1153;	gcp[2].y = 2133;	gcp[2].z = 3000;
			//	gcp[3].x = 1000;	gcp[3].y = 2133;	gcp[3].z = 3000;

			//	gcp[0].x = 1000;	gcp[0].y = 2000;	gcp[0].z = 3000;
			//	gcp[1].x = 1300;	gcp[1].y = 2000;	gcp[1].z = 3000;
			//	gcp[2].x = 1300;	gcp[2].y = 2200;	gcp[2].z = 3000;
			//	gcp[3].x = 1000;	gcp[3].y = 2200;	gcp[3].z = 3000;

			gcp[0].x = 0;	gcp[0].y = 0;	gcp[0].z = 0;
			gcp[1].x = 0;	gcp[1].y = 440;	gcp[1].z = 0;
			gcp[2].x = 230;	gcp[2].y = 0;	gcp[2].z =0;
			gcp[3].x = 230;	gcp[3].y = 440;	gcp[3].z = 0;

			imagepoint[0].x = fiducial[0];	imagepoint[0].y = fiducial[1];
			imagepoint[1].x = fiducial[2];	imagepoint[1].y = fiducial[3];
			imagepoint[2].x = fiducial[4];	imagepoint[2].y = fiducial[5];
			imagepoint[3].x = fiducial[6];	imagepoint[3].y = fiducial[7];

			for (int i = 0; i<4; i++)
			{
				photopoint[i] = SR.Photo.Image2CalibratedPhoto(imagepoint[i]);
			}

			omega = 0;
			phi = 0;
			kappa = -atan2(photopoint[1].y - photopoint[0].y, photopoint[1].x - photopoint[0].x);
			//	XL = 1070;
			//	YL = 2060;
			//	ZL = 3500;
			XL = 1150;
			YL = 1700;
			ZL = 3500;

			if (dist_2P(imagepoint[0], imagepoint[1])>dist_2P(imagepoint[2], imagepoint[3]))
			{
				omega = 20 * CV_PI / 180.0;
				YL -= 200;
			}
			else if (dist_2P(imagepoint[0], imagepoint[1])<dist_2P(imagepoint[2], imagepoint[3]))
			{
				omega = -20 * CV_PI / 180.0;
				YL += 200;
			}

			if (dist_2P(imagepoint[1], imagepoint[2])>dist_2P(imagepoint[3], imagepoint[0]))
			{
				phi = 20 * CV_PI / 180.0;
				XL += 200;
			}
			else if (dist_2P(imagepoint[0], imagepoint[1])<dist_2P(imagepoint[2], imagepoint[3]))
			{
				phi = -20 * CV_PI / 180.0;
				XL -= 200;
			}

			//	SR.Photo.Set_EOP(1070, 2060, 3500, -62.900,-0.740,-64.300);

			SR.Photo.Set_EOP(XL, YL, ZL, omega, phi, kappa);
			SR.SpaceResection(SR.Photo, ngcps, photopoint, gcp, sigma, precision, residual);


			double eop[6];

			eop[0] = SR.Photo.EOP.Point.x;
			eop[1] = SR.Photo.EOP.Point.y;
			eop[2] = SR.Photo.EOP.Point.z;
			eop[3] = SR.Photo.EOP.o;
			eop[4] = SR.Photo.EOP.p;
			eop[5] = SR.Photo.EOP.k;

			Point2f targetphoto[3];
			Point2f targetimage[3];
			Point3f target_ground[3];
			double target_r[3];
			double target_r_pixel[3];
			double targetcoord[6];
			double targetradious[3];

			//target_ground[0].x = 1090;		target_ground[0].y = 2080;		target_ground[0].z = 3020;
			//target_ground[1].x = 1071;		target_ground[1].y = 1969;		target_ground[1].z = 3001;
			//target_ground[2].x = 1215;		target_ground[2].y = 2163;		target_ground[2].z = 3008.5;

			target_ground[0].x = 110;		target_ground[0].y = 110;		target_ground[0].z = 20;
			target_ground[1].x = 200;		target_ground[1].y =90 ;		target_ground[1].z = 30;
			target_ground[2].x = 150;		target_ground[2].y = 300;		target_ground[2].z = 10;

			target_r[0] = 7.0;
			target_r[1] = 7.0;
			target_r[2] = 10;

			int i = 0;
			Point2f tempvar;
			for (int k = 0; k<kinectpts.size(); k++)
			{
				//targetphoto[i] = SR.CollinearityEquations(target_ground[i], SR.Photo.Get_IOP(), SR.Photo.Get_EOP());
				
				tempvar = SR.CollinearityEquations(kinectpts.at(k), SR.Photo.Get_IOP(), SR.Photo.Get_EOP());
				//targetimage[i] = SR.Photo.CalibratedPhoto2Image(targetphoto[i]);
				tempvar = SR.Photo.CalibratedPhoto2Image(tempvar);
				//targetcoord[i * 2] = tempvar.x;
				//targetcoord[i * 2 + 1] = tempvar.y;
				projpts.push_back(tempvar.y);
				projpts.push_back(tempvar.x);
				//double dist = dist3D(SR.Photo.Get_EOP().Point, target_ground[i]);
				//target_r_pixel[i] = SR.Photo.Get_IOP().f*target_r[i] / (dist*0.005);
				//targetradious[i] = target_r_pixel[i];
			}

			//	jintArray result = env->NewIntArray(8);
			//	jdoubleArray result = env->NewDoubleArray(8);
			/*jdoubleArray result = env->NewDoubleArray(23);
			if (result == NULL)
			{
				return NULL;
			}
			//	env->SetIntArrayRegion(result, 0, 8, fiducial);
			env->SetDoubleArrayRegion(result, 0, 8, fiducial);
			env->SetDoubleArrayRegion(result, 8, 6, eop);
			env->SetDoubleArrayRegion(result, 14, 6, targetcoord);
			env->SetDoubleArrayRegion(result, 20, 3, targetradious);
			*/
			//return result;
		
			//return NULL;
	}

	
	
#ifdef __cplusplus
}
#endif
