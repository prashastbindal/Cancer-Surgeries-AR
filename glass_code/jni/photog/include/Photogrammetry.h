// Photogrammetry.h: interface for the CPhotogrammetryTrmble class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PHOTOGRAMMETRY_H__363B34AE_923E_4730_A86A_6714832811C9__INCLUDED_)
#define AFX_PHOTOGRAMMETRY_H__363B34AE_923E_4730_A86A_6714832811C9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Photo.h"
#include "Matrix.h"
//#include "cv.h"

struct IntersectionParam4
{
	int CheckBA;
	int nCams;
	int Count;
	Point3D <double> ObjCoordEst;
	double STD[3];
	double RMSE;
	double Residual[8];
	ImageCoordinates <double> ImgCoordObs[4];
	ImageCoordinates <double> ImgCoordEst[4];
	int tempInt0[4];
	int tempInt1[4];
	int IsTrackingSuccessful[4];
};

class IntersectionParam
{
public:
	IntersectionParam()
	{
		Residual = NULL;
		ImgCoordObs = NULL;
		ImgCoordEst = NULL;
		tempInt0 = NULL;
		tempInt1 = NULL;
	};
	IntersectionParam(int ncams)
	{
		Count = 0;
		nCams = ncams;
		Residual = new double [nCams*2];
		ImgCoordObs = new ImageCoordinates <double> [nCams];
		ImgCoordEst = new ImageCoordinates <double> [nCams];
		tempInt0 = new int [nCams];
		tempInt1 = new int [nCams];
	};
	virtual ~IntersectionParam()
	{
		if(Residual!=NULL)
			delete[] Residual;
		if(ImgCoordObs!=NULL)
			delete[] ImgCoordObs;
		if(ImgCoordEst!=NULL)
			delete[] ImgCoordEst;
		if(tempInt0!=NULL)
			delete[] tempInt0;
		if(tempInt1!=NULL)
			delete[] tempInt1;
	};

	void Set_nCams(int ncams)
	{
		Count = 0;
		nCams = ncams;
		Residual = new double [nCams*2];
		ImgCoordObs = new ImageCoordinates <double> [nCams];
		ImgCoordEst = new ImageCoordinates <double> [nCams];
		tempInt0 = new int [nCams];
		tempInt1 = new int [nCams];
	};

	int CheckBA;
	int nCams;
	int Count;
	Point3D <double> ObjCoordEst;
	double STD[3];
	double RMSE;
	double* Residual;
	ImageCoordinates <double>* ImgCoordObs;
	ImageCoordinates <double>* ImgCoordEst;
	int* tempInt0;
	int* tempInt1;	
};

class CPhotogrammetry  
{
public:
	CPhotogrammetry();
	CPhotogrammetry(int nphotos);
	virtual ~CPhotogrammetry();

protected:
	void Construction();
	void Initialize_Photo();
	void Initialize_FundamentalMatrix();

protected:
	void Get_Elements_of_RotationMatrix(double o, double p, double k);
	void Calculate_Parameters(InteriorOrientationParameters iop, OrientationParameters eop, Point2D <double> photopoint, Point3D <double> gcp);
	void Calculate_Parameters(int photoID, Point2D <double> photopoint, Point3D <double> gcp);

protected:
	int nPhotos;
	int n, m;

	// For iteration termination conditions
	int iteration_stop;
	int iteration_count;
	int iteration_tolerence;

	double angle_tolerence;
	double position_tolerence;
	double sigma_o_2_hat_change_tolerence;

	double sigma_o_2_hat;
	double sigma_o_2_hat_change;
	double sigma_o_2_hat_old;
	
	//sin & cos
	double sin_o, sin_p, sin_k;
	double cos_o, cos_p, cos_k;

	// Elements of the rotation matrix
	double m11,m12,m13;
	double m21,m22,m23;
	double m31,m32,m33;

	// Parameters for linearization of the collinearity equation
	double q, r, s;
	double b11,  b12,  b13,  b14,  b15,  b16,  J;					//for x
	double b21,  b22,  b23,  b24,  b25,  b26,  K;					//for y

	//Matrix formation
	Matrix <double> A;
	Matrix <double> AT;
	Matrix <double> y;
	Matrix <double> csi_hat;
	Matrix <double> c;
	Matrix <double> e_tilda;
	Matrix <double> N;
	Matrix <double> Ninv;
	Matrix <double> P;
	Matrix <double> Precision;

	// Computer vision (Fundamental matrix)
	Matrix <double>** FMatrix;
	Matrix <double>** Epipole;

public:	
	CPhoto* Photo;
	void Set_Photo(int nphotos);	
	int Get_nPhotos();
	void Get_3DPoints(IntersectionParam& IP);
	void Get_3DPoints(IntersectionParam4& IP);
	int Get_3DPoints(Point3D <double>& initial, double residual[], double& rmse);
	int Get_3DPoints(Point3D <double>& initial, double residual[], double& rmse, double std[]);
	ImageCoordinates <double> Get_ImageCoordinates(CPhoto photo, Point3D <double> point);
//	CvPoint2D Get_ImageCoordinates(CPhoto photo, double x, double y, double z);
	Point2D <double> CollinearityEquations( Point3D <double> point, InteriorOrientationParameters iop, OrientationParameters eop); 
	
	// Computer vision (Fundamental matrix)
	void EstimateFundamentalMatrix();
	Matrix <double> Get_FundamentalMatrix(int fromCAM, int toCAM);
	Matrix <double> Get_Epipole(int fromCAM, int toCAM);

	int Space_Resection(CPhoto& photo, int ngcps, Point2D <double>* photopoint, Point3D <double>* gcp, double& sigma, double precision[], double residual[]);
};

#endif // !defined(AFX_PHOTOGRAMMETRY_H__363B34AE_923E_4730_A86A_6714832811C9__INCLUDED_)
