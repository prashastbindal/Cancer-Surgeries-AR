#if !defined(AFX_PHOTO_H__363B3453454353434AE_923E_4730_A86A_6714832811C9__INCLUDED_)
#define AFX_PHOTO_H__363B3453454353434AE_923E_4730_A86A_6714832811C9__INCLUDED_

#define M_PI		3.14159265358979323846
#define M_PI_2		1.57079632679489661923
#define NUMBER_ID   5
#define NO_VALUE  -99999
#define BIG_VALUE  99999

#include "Matrix.h"

template<class T>
class ImageCoordinates
{
public:
	T Row;
	T Col;

	ImageCoordinates() {Row=0; Col=0;}		
	ImageCoordinates(T row, T col) { Row=row; Col=col; }
	ImageCoordinates(const ImageCoordinates &other) { Row=other.Row; Col=other.Col; }	
	virtual ~ImageCoordinates() {}

	void SetPoint(T row, T col) { Row=row; Col=col; }
};

template<class T>
class Point2D  
{
public:
	T x;
	T y;

	Point2D() {	x=0; y=0; }
	Point2D(T xx,T yy) { x=xx; y=yy; }
	Point2D(const Point2D &other) { x=other.x; y=other.y; }	
	virtual ~Point2D() {}

	void SetPoint(T xx,T yy) { x=xx; y=yy; }
};

template<class T>
class Point3D  
{
public:
	T x;
	T y;
	T z;

	Point3D() {	x=0; y=0; z=0; }
	Point3D(T xx,T yy, T zz) { x=xx; y=yy; z=zz; }
	Point3D(const Point3D &other) { x=other.x; y=other.y; z=other.z; }	
	virtual ~Point3D() {}

	void SetPoint(T xx,T yy, T zz) { x=xx; y=yy; z=zz;}
	void SetPoint(Point3D <T> point3d) {x = point3d.x; y = point3d.y; z = point3d.z;} 

	Point3D<T> operator = (const Point3D<T>& TempPoint)
	{
		x = TempPoint.x;
		y = TempPoint.y;
		z = TempPoint.z;
		return *this;
	}

	Point3D<T> operator - (const Point3D<T>& TempPoint)
	{
		Point3D <T> returnPoint;
		returnPoint.x = x - TempPoint.x;
		returnPoint.y = y - TempPoint.y;
		returnPoint.z = z - TempPoint.z;
		return returnPoint;
	}
};

template<class T>
class Line2D  
{
	//ax + by + c = 0
public:
	T a;
	T b;
	T c;

	Line2D() {	a=0; b=0; c=0; }
	Line2D(T aa,T bb, T cc) { a=aa; b=bb; c=cc; }
	Line2D(const Line2D &other) { a=other.a; b=other.b; c=other.c; }	
	virtual ~Line2D() {}

	void SetParam(T aa,T bb, T cc) { x=aa; y=bb; z=cc;}
};

struct OrientationParameters
{
	Point3D <double> Point;
	double o;
	double p;
	double k;

	void Set_EOP(double x, double y, double z, double omega, double phi, double kappa)
	{
		Point.x = x;
		Point.y = y;
		Point.z = z;
		o = omega;
		p = phi;
		k = kappa;
	}

	void Set_EOP(OrientationParameters eop)
	{
		Point.x = eop.Point.x;
		Point.y = eop.Point.y;
		Point.z = eop.Point.z;
		o = eop.o;
		p = eop.p;
		k = eop.k;
	}
};

struct InteriorOrientationParameters
{
	double f;
	Point2D <double> PP;
	double p1;
	double p2;
	double p3;
	double k1;
	double k2;
	double k3;	

	InteriorOrientationParameters(){};
	InteriorOrientationParameters(double f, double x, double y, 
								  double k1, double k2, double k3,
								  double p1, double p2, double p3)
	{
		this->f		= f;
		this->PP.x	= x;
		this->PP.y	= y;
		this->k1	= k1;
		this->k2	= k2;
		this->k3	= k3;
		this->p1	= p1;
		this->p2	= p2;
		this->p3	= p3;
	}
	
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
};

class CPhoto
{
public:
	CPhoto();
	virtual ~CPhoto();

protected:
	InteriorOrientationParameters	IOP;
	OrientationParameters			EOP;
	double	PixelSize;				//Single pixel size
	int		nRows, nCols;			//Row, column
	double  half_nRows, half_nCols;	// for calculation efficiency

protected:
	// Elements of the rotation matrix
	double m11,m12,m13;
	double m21,m22,m23;
	double m31,m32,m33;
	Matrix <double> R;

protected:
	void Get_Elements_of_RotationMatrix();

public:	
	InteriorOrientationParameters	Get_IOP(){return IOP;}
	OrientationParameters			Get_EOP(){return EOP;}
	void	Set_IOP(InteriorOrientationParameters iop);
	void	Set_IOP(double  f, double  x, double  y, 
					double k1, double k2, double k3,
					double p1, double p2, double p3);
	void	Set_EOP(OrientationParameters eop);
	void	Set_EOP(double x, double y, double z,
					double o, double p, double k);

	double	Get_PixelSize(){return PixelSize;}
	int		Get_nRows(){return nRows;}
	int		Get_nCols(){return nCols;}
	void	Set_PixelSize(double pixel_size){PixelSize = pixel_size;}
	void	Set_nRows(int nrows){nRows = nrows;	half_nRows  = (double)nRows/2.0;}
	void	Set_nCols(int ncols){nCols = ncols; half_nCols  = (double)nCols/2.0;}

	// For coordinates observations
	Point2D <double>	PhotoPoint;
	void Set_PhotoPoint(Point2D <double> photopoint);
	void Set_PhotoPoint(double x, double y);

	ImageCoordinates <double> ImagePoint;
	void	Set_ImagePoint(ImageCoordinates <double> imagepoint);
	void	Set_ImagePoint(double row, double col);

	// Coordinate transformation 
	int Image2Photo();				// Image2Photo
	int Photo2Image();				// Photo2Image
	int Photo2CalibratedPhoto();	// Photo2CalibratedPhoto
	int CalibratedPhoto2Photo();	// CalibratedPhoto2Photo
	int Image2CalibratedPhoto();	// Image2CalibratedPhoto
	int CalibratedPhoto2Image();	// CalibratedPhoto2Image

	Point2D <double> Image2Photo(const ImageCoordinates <double>& image);
	Point2D <double> Image2Photo(const double& row, const double& col);
	ImageCoordinates <double> Photo2Image(const Point2D <double>& photo);
	ImageCoordinates <double> Photo2Image(const double &x, const double& y);

	Point2D <double> Photo2CalibratedPhoto(Point2D <double> photo);
	Point2D <double> CalibratedPhoto2Photo(Point2D <double> cphoto);
	
	Point2D <double> Image2CalibratedPhoto(double row, double col);
	Point2D <double> Image2CalibratedPhoto(ImageCoordinates <double> image);
	ImageCoordinates <double> CalibratedPhoto2Image(Point2D <double> cphoto);

	Matrix <double> Get_RotationMatrix();

	// Computer Vision
	// Camera matrix which is not considering lens distortion parameters
	Matrix <double> GetCameraCenterMatrix();
	Matrix <double> GetCameraMatrix();
	Matrix <double> GetPseudoInverseCameraMatrix();
};

#endif // !defined(AFX_PHOTO_H__363B3453454353434AE_923E_4730_A86A_6714832811C9__INCLUDED_)

