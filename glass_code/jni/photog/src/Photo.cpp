// Photo.cpp: implementation of the CPhoto class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Photo.h"
#include "Matrix.h"
#include <math.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPhoto::CPhoto()
{
	IOP.k1 =	0;
	IOP.k2 =	0;
	IOP.k3 =	0;
	IOP.p1 =	0;
	IOP.p2 =	0;
	IOP.p3 =	0;
	IOP.PP.x =	0;
	IOP.PP.y =	0;

	EOP.o = NO_VALUE;
	EOP.p = NO_VALUE;
	EOP.k = NO_VALUE;
	EOP.Point.x = NO_VALUE;
	EOP.Point.y = NO_VALUE;
	EOP.Point.z = NO_VALUE;
}

CPhoto::~CPhoto()
{
}

void CPhoto::Set_ImagePoint(ImageCoordinates <double> imagepoint)
{
	ImagePoint = imagepoint;
}

void CPhoto::Set_ImagePoint(double row, double col)
{
	ImagePoint.Row = row;
	ImagePoint.Col = col;
}

void CPhoto::Set_IOP(InteriorOrientationParameters iop)
{
	IOP = iop;
}

void CPhoto::Set_IOP(double  f, double  x, double  y, 
			 		 double k1, double k2, double k3,
					 double p1, double p2, double p3)
{
	IOP.f = f;
	IOP.PP.x = x;
	IOP.PP.y = y;
	IOP.k1 = k1;
	IOP.k2 = k2;
	IOP.k3 = k3;
	IOP.p1 = p1;
	IOP.p2 = p2;
	IOP.p3 = p3;
}

void CPhoto::Set_EOP(OrientationParameters eop)
{
	EOP = eop;

	Get_Elements_of_RotationMatrix();

	R.Resize(3,3);
	R(0,0) = m11;
	R(0,1) = m12;
	R(0,2) = m13;

	R(1,0) = m21;
	R(1,1) = m22;
	R(1,2) = m23;
	
	R(2,0) = m31;
	R(2,1) = m32;
	R(2,2) = m33;
}

void CPhoto::Set_EOP(double x, double y, double z, double o, double p, double k)
{
	EOP.Point.x = x;
	EOP.Point.y = y;
	EOP.Point.z = z;
	EOP.o = o;
	EOP.p = p;
	EOP.k = k;

	Get_Elements_of_RotationMatrix();

	R.Resize(3,3);
	R(0,0) = m11;
	R(0,1) = m12;
	R(0,2) = m13;

	R(1,0) = m21;
	R(1,1) = m22;
	R(1,2) = m23;
	
	R(2,0) = m31;
	R(2,1) = m32;
	R(2,2) = m33;
}

void CPhoto::Set_PhotoPoint(Point2D <double> photopoint)
{
	PhotoPoint = photopoint;
}

void CPhoto::Set_PhotoPoint(double x, double y)
{
	PhotoPoint.x = x;
	PhotoPoint.y = y;
}

int CPhoto::Image2Photo()
{
	PhotoPoint = Image2Photo(ImagePoint); 	
	
	return 1;
}

Point2D <double> CPhoto::Image2Photo(const ImageCoordinates <double>& image)
{
	Point2D <double> Photo;

	if (image.Row!=NO_VALUE && image.Col!=NO_VALUE)
	{
		Photo.x = (image.Col - half_nCols)*PixelSize;
		Photo.y = (half_nRows - image.Row)*PixelSize;
	}
	else
	{
		Photo.x = NO_VALUE;
		Photo.y = NO_VALUE;
	}
	
	return Photo;
}

Point2D <double> CPhoto::Image2Photo(const double& row, const double& col)
{
	Point2D <double> Photo;

	if (row!=NO_VALUE && col!=NO_VALUE)
	{
		Photo.x = (col - half_nCols)*PixelSize;
		Photo.y = (half_nRows - row)*PixelSize;
	}
	else
	{
		Photo.x = NO_VALUE;
		Photo.y = NO_VALUE;
	}
	
	return Photo;
}

int CPhoto::Photo2Image()
{
	ImagePoint = Photo2Image(PhotoPoint);		

	return 1;
}

ImageCoordinates <double> CPhoto::Photo2Image(const Point2D <double>& photo)
{
	ImageCoordinates <double> Image;

	if (photo.x!=NO_VALUE && photo.y!=NO_VALUE)
	{
		Image.Col = photo.x/PixelSize + half_nCols;
		Image.Row = half_nRows - photo.y/PixelSize;
	}
	else
	{
		Image.Col = NO_VALUE;
		Image.Row = NO_VALUE;
	}

	return Image;
}


ImageCoordinates <double> CPhoto::Photo2Image(const double &x, const double& y)
{
	Point2D <double> photo;
	photo.x = x;
	photo.y = y;

	return Photo2Image(photo);
}

Point2D <double> CPhoto::Photo2CalibratedPhoto(Point2D <double> photo)
{
	Point2D <double> CPhoto;

	if (photo.x!=NO_VALUE && photo.y!=NO_VALUE)
	{
		double drx, dry;
		double ddx, ddy;
		double x_bar,  y_bar;
		double xy_bar;
		double x_bar2, y_bar2;
		double r2, r4, r6;

		x_bar	= photo.x - IOP.PP.x;
		y_bar	= photo.y - IOP.PP.y;
		xy_bar  = x_bar*y_bar;
		x_bar2  = x_bar*x_bar;
		y_bar2  = y_bar*y_bar;
		r2		= x_bar2 + y_bar2;
		r4		= r2*r2;
		r6		= r4*r2;

		// Radial lens distortion
		drx  = x_bar*(IOP.k1*r2 + IOP.k2*r4 + IOP.k3*r6);
		dry  = y_bar*(IOP.k1*r2 + IOP.k2*r4 + IOP.k3*r6);

		// Decentering lens distortion
		ddx = IOP.p1*(r2 + 2*x_bar2) + 2*IOP.p2*xy_bar;
		ddy = 2*IOP.p1*xy_bar + IOP.p2*(r2 + 2*y_bar2);

		CPhoto.x = photo.x + (drx + ddx);
	   	CPhoto.y = photo.y + (dry + ddy);
	}
	else
	{
		CPhoto.x = NO_VALUE;
		CPhoto.y = NO_VALUE;
	}

	return CPhoto;
}

int CPhoto::CalibratedPhoto2Image()
{
	ImagePoint = Photo2Image(CalibratedPhoto2Photo(PhotoPoint));

 	return 1;
}

int CPhoto::CalibratedPhoto2Photo()
{
	PhotoPoint = CalibratedPhoto2Photo(PhotoPoint);
	
	return 1;
}

int CPhoto::Photo2CalibratedPhoto()
{	
	PhotoPoint = Photo2CalibratedPhoto(PhotoPoint);

	return 1;
}

int CPhoto::Image2CalibratedPhoto()
{	
	PhotoPoint = Image2CalibratedPhoto(ImagePoint);

	return 1;
}

Point2D <double> CPhoto::Image2CalibratedPhoto(double row, double col)
{
	ImageCoordinates <double> Image;
	Image.Row = row;
	Image.Col = col;

	return Image2CalibratedPhoto(Image);
}

Point2D <double> CPhoto::Image2CalibratedPhoto(ImageCoordinates <double> image)
{
	Point2D <double> Photo;

	if (image.Row!=NO_VALUE && image.Col!=NO_VALUE)
	{
		Photo = Photo2CalibratedPhoto(Image2Photo(image));
	}
	else
	{
		Photo.x = NO_VALUE;
		Photo.y = NO_VALUE;
	}

	return Photo;
}

ImageCoordinates <double> CPhoto::CalibratedPhoto2Image(Point2D <double> cphoto)
{
	ImageCoordinates <double> Image;

//	bool control = (cphoto.x!=NO_VALUE && cphoto.y!=NO_VALUE);

	if (cphoto.x!=NO_VALUE && cphoto.y!=NO_VALUE)
	{
		Image = Photo2Image(CalibratedPhoto2Photo(cphoto));
	}
	else
	{
		Image.Row = NO_VALUE;
		Image.Col = NO_VALUE;
	}

	if (cphoto.x!=NO_VALUE && cphoto.y!=NO_VALUE)
	{
		Image = Photo2Image(CalibratedPhoto2Photo(cphoto));
	}
	else
	{
		Image.Row = NO_VALUE;
		Image.Col = NO_VALUE;
	}

	return Image;
}

Point2D <double> CPhoto::CalibratedPhoto2Photo(Point2D <double> cphoto)
{
	Point2D <double> Photo;

	if (cphoto.x!=NO_VALUE && cphoto.y!=NO_VALUE)
	{
		Matrix <double> A(2,2);
		Matrix <double> csi_hat(2,1);
		Matrix <double> y(2,1);
		Matrix <double> e_tilda;

		double xo, yo;	// initial values
		double x_bar,  y_bar;
		double xy_bar;
		double x_bar2, y_bar2;
		double r2, r4, r6;
		double K1,K2;
		double drx, dry;
		double ddx, ddy;

		xo = cphoto.x;
		yo = cphoto.y;


		//For the iteration termination conditions////////////////////

		int iteration_stop  = 0;
		int iteration_count = 0;

		double sigma_o_2_hat_change;
		double sigma_o_2_hat;
		double sigma_o_2_hat_old	= 9999999999;
		
		while(iteration_stop==0)
		{
			x_bar	= xo - IOP.PP.x;
			y_bar	= yo - IOP.PP.y;
			xy_bar  = x_bar*y_bar;
			x_bar2  = x_bar*x_bar;
			y_bar2  = y_bar*y_bar;
			r2		= x_bar2 + y_bar2;
			r4		= r2*r2;
			r6		= r4*r2;
			
			K1 = IOP.k1*r2 + IOP.k2*r4 + IOP.k3*r6;
			K2 = IOP.k1 + 2*IOP.k2*r2 + 3*IOP.k3*r4;

			// Radial lens distortion
			drx  = x_bar*K1;
			dry  = y_bar*K1;

			// Decentering lens distortion
			ddx = IOP.p1*(r2 + 2*x_bar2) + 2*IOP.p2*xy_bar;
			ddy = 2*IOP.p1*xy_bar + IOP.p2*(r2 + 2*y_bar2);

			// Matrix
			A(0,0) = 1 + K1 + 2*x_bar2*K2 + 6*IOP.p1*x_bar + 2*IOP.p2*y_bar;
			A(0,1) = 2*xy_bar*K2 + 2*IOP.p1*y_bar + 2*IOP.p2*x_bar;
			A(1,0) = A(0,1);
			A(1,1) = 1 + K1 + 2*y_bar2*K2 + 2*IOP.p1*x_bar + 6*IOP.p2*y_bar;

			y(0,0) = cphoto.x - (xo + drx + ddx);
			y(1,0) = cphoto.y - (yo + dry + ddy);

			// LESS
			csi_hat			= (A.Transpose()*A).Inverse()*A.Transpose()*y;
			e_tilda			= y - A*csi_hat;
			sigma_o_2_hat	= (e_tilda.Transpose()*e_tilda)(0,0)/2;

			// Update				
			xo += csi_hat(0,0);
			yo += csi_hat(1,0);

			iteration_count++;

			//Iteration Termination Conditions///////////////////////////////
			//Wolf, Ghilani, Adjustment Computations, pp. 243~245

			//1. Maximum iteration
			if(iteration_count==10)
				iteration_stop = -1;

			//2. Maximum correction


			//3. Minitoring the Reference Variance
			sigma_o_2_hat_change = fabs((sigma_o_2_hat - sigma_o_2_hat_old)/sigma_o_2_hat_old);

			if(sigma_o_2_hat_change<0.01)
				iteration_stop = 1;

			sigma_o_2_hat_old = sigma_o_2_hat;
		}

		Photo.x = xo;
		Photo.y = yo;
	}
	else
	{
		Photo.x = NO_VALUE;
		Photo.y = NO_VALUE;
	}

	return Photo;
}

void CPhoto::Get_Elements_of_RotationMatrix()
{
	//For calculation efficiency
	double sin_o = sin(EOP.o);
	double sin_p = sin(EOP.p);
	double sin_k = sin(EOP.k);
	double cos_o = cos(EOP.o);
	double cos_p = cos(EOP.p);
	double cos_k = cos(EOP.k);

	//Calculation of elements of rotation matrix
	m11 =  cos_p*cos_k;
	m12 =  sin_o*sin_p*cos_k + cos_o*sin_k;
	m13 = -cos_o*sin_p*cos_k + sin_o*sin_k;

	m21 = -cos_p*sin_k;
	m22 = -sin_o*sin_p*sin_k + cos_o*cos_k;
	m23 =  cos_o*sin_p*sin_k + sin_o*cos_k;

	m31 =  sin_p;
	m32 = -sin_o*cos_p;
	m33 =  cos_o*cos_p;
}

Matrix <double> CPhoto::Get_RotationMatrix()
{
	return R;
}

Matrix <double> CPhoto::GetCameraMatrix()
{
	Matrix <double> P(3,4);
	/*
	P(0,0) = -IOP.f*0.001*m11; 
    P(0,1) = -IOP.f*0.001*m12;
    P(0,2) = -IOP.f*0.001*m13;
	P(0,3) =  IOP.f*0.001*(m11*EOP.Point.x + m12*EOP.Point.y + m13*EOP.Point.z);
    
    P(1,0) = -IOP.f*0.001*m21;
    P(1,1) = -IOP.f*0.001*m22;
    P(1,2) = -IOP.f*0.001*m23;
    P(1,3) =  IOP.f*0.001*(m21*EOP.Point.x + m22*EOP.Point.y + m23*EOP.Point.z);
*/
	P(0,0) = -IOP.f*m11; 
    P(0,1) = -IOP.f*m12;
    P(0,2) = -IOP.f*m13;
	P(0,3) =  IOP.f*(m11*EOP.Point.x + m12*EOP.Point.y + m13*EOP.Point.z);
    
    P(1,0) = -IOP.f*m21;
    P(1,1) = -IOP.f*m22;
    P(1,2) = -IOP.f*m23;
    P(1,3) =  IOP.f*(m21*EOP.Point.x + m22*EOP.Point.y + m23*EOP.Point.z);
  
    P(2,0) =   m31;
    P(2,1) =   m32;
    P(2,2) =   m33;
    P(2,3) =  -(m31*EOP.Point.x + m32*EOP.Point.y + m33*EOP.Point.z);

	return P;
}

 Matrix <double> CPhoto::GetCameraCenterMatrix()
 {
	Matrix <double> C(4,1);

	C(0,0) = EOP.Point.x;
	C(1,0) = EOP.Point.y;
	C(2,0) = EOP.Point.z;
	C(3,0) = 1.0;

	return C;
 }

Matrix <double> CPhoto::GetPseudoInverseCameraMatrix()
{
	Matrix <double> P = GetCameraMatrix();	
	Matrix <double> Ppinv = P.Transpose()*(P*P.Transpose()).Inverse();
	
	return Ppinv;
}