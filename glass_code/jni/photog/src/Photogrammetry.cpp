// Photogrammetry.cpp: implementation of the CPhotogrammetry class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Photogrammetry.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPhotogrammetry::CPhotogrammetry()
{
	Construction();
}

CPhotogrammetry::CPhotogrammetry(int nphotos)
{
	Construction();
	Set_Photo(nphotos);
}

void CPhotogrammetry::Construction()
{
	Photo = NULL;

	// For iteration control
	iteration_stop		= 0;
	iteration_count		= 0;
	sigma_o_2_hat_old	= 9999999999;
	
	// For iteration termination conditions
	iteration_tolerence				= 10;
	angle_tolerence					= 0.000001;
	position_tolerence				= 0.0001;
	sigma_o_2_hat_change_tolerence	= 0.01;	
	sigma_o_2_hat_change_tolerence	= 0.001;	

	// Computer vision (Fundamental matrix)
	FMatrix = NULL;
	Epipole = NULL;
}

CPhotogrammetry::~CPhotogrammetry()
{
	Initialize_Photo();
}

void CPhotogrammetry::Initialize_Photo()
{
	if(Photo!=NULL)
		delete[] Photo;
	Photo = NULL;
}

void CPhotogrammetry::Initialize_FundamentalMatrix()
{
	int i;

	if(FMatrix!=NULL)
	{
		for(i=0;i<nPhotos;i++)
		{
			delete[] FMatrix[i];
		}
		delete[] FMatrix;
	}
	if(Epipole!=NULL)
	{
		for(i=0;i<nPhotos;i++)
		{
			delete[] Epipole[i];
		}
		delete[] Epipole;
	}
}

void CPhotogrammetry::Set_Photo(int nphotos)
{
	Initialize_Photo();

	if(nphotos>0)
	{
		nPhotos = nphotos;
		Photo = new CPhoto [nphotos];
	}
}

int CPhotogrammetry::Get_nPhotos()
{
	return nPhotos;
}

Point2D <double> CPhotogrammetry::CollinearityEquations(Point3D <double> point, 
														InteriorOrientationParameters iop, 
														OrientationParameters eop)
{
	Point2D <double> PhotoPoint;
	double r, s, q;

	Get_Elements_of_RotationMatrix(eop.o, eop.p, eop.k);
	
	r = m11*(point.x - eop.Point.x) + m12*(point.y - eop.Point.y) + m13*(point.z - eop.Point.z);
	s = m21*(point.x - eop.Point.x) + m22*(point.y - eop.Point.y) + m23*(point.z - eop.Point.z);
	q = m31*(point.x - eop.Point.x) + m32*(point.y - eop.Point.y) + m33*(point.z - eop.Point.z);

	PhotoPoint.x = iop.PP.x - iop.f*r/q;
	PhotoPoint.y = iop.PP.y - iop.f*s/q;

	return PhotoPoint;
}

ImageCoordinates <double> CPhotogrammetry::Get_ImageCoordinates(CPhoto photo, Point3D <double> point)
{
	ImageCoordinates <double> imagepoint;
	Point2D <double> photopoint;

	photopoint = CollinearityEquations(point, photo.Get_IOP(), photo.Get_EOP());
	imagepoint = photo.CalibratedPhoto2Image(photopoint);

	return imagepoint;
}

void CPhotogrammetry::Get_Elements_of_RotationMatrix(double o, double p, double k)
{
	//For calculation efficiency
	sin_o = sin(o);
	sin_p = sin(p);
	sin_k = sin(k);
	cos_o = cos(o);
	cos_p = cos(p);
	cos_k = cos(k);

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

void CPhotogrammetry::Calculate_Parameters(	int photoID,
											Point2D <double> photopoint,
											Point3D <double> gcp)
{
	//Calculation of Elements of Rotation matrix////////////////////////////

	Matrix <double> R = Photo[photoID].Get_RotationMatrix();

	m11  = R(0,0);	m12 = R(0,1);	m13 = R(0,2);
	m21  = R(1,0);	m22 = R(1,1);	m23 = R(1,2);
	m31  = R(2,0);	m32 = R(2,1);	m33 = R(2,2);

	//For calculation efficiency////////////////////////////////////////////

	Point3D <double> Delta;
	double temp1, temp2, temp3, temp4, temp5;
	double fq;
	double fqq;	

	Delta.x = gcp.x - Photo[photoID].Get_EOP().Point.x;
	Delta.y = gcp.y - Photo[photoID].Get_EOP().Point.y;
	Delta.z = gcp.z - Photo[photoID].Get_EOP().Point.z;

	q = m31*Delta.x + m32*Delta.y + m33*Delta.z;
	r = m11*Delta.x + m12*Delta.y + m13*Delta.z;
	s = m21*Delta.x + m22*Delta.y + m23*Delta.z;

	fq	= Photo[photoID].Get_IOP().f/q;
	fqq	= Photo[photoID].Get_IOP().f/(q*q);


	//Calculation derivatives//////////////////////////////////////////////////

	//For calculation efficiency
	temp1 = -m33*Delta.y + m32*Delta.z;
	temp2 = cos_p*Delta.x + sin_o*sin_p*Delta.y - cos_o*sin_p*Delta.z;
	temp3 = sin_o*cos_p*Delta.y;
	temp4 = cos_o*cos_p*Delta.z;
	temp5 = sin_p*Delta.x;

	//For x
	b11 =  fqq*(r*(temp1) - q*(-m13*Delta.y + m12*Delta.z));				// dF/d(omega)
	b12 =  fqq*(r*(temp2) - q*(-cos_k*temp5 + cos_k*temp3 - cos_k*temp4));  // dF/d(phi)
	b13 =  -fq*(m21*Delta.x + m22*Delta.y + m23*Delta.z);					// dF/d(kappa)
	b14 =  fqq*(r*m31 - q*m11);												// dF/d(XL)
	b15 =  fqq*(r*m32 - q*m12);												// dF/d(YL)
	b16 =  fqq*(r*m33 - q*m13);												// dF/d(ZL)
	
	J	=  photopoint.x - Photo[photoID].Get_IOP().PP.x + fq*r;

	//For y
	b21 =  fqq*(s*(temp1) - q*(-m23*Delta.y + m22*Delta.z));				// dG/d(omega)
	b22 =  fqq*(s*(temp2) - q*(sin_k*temp5 - sin_k*temp3 + sin_k*temp4));	// dG/d(phi)
	b23 =   fq*(m11*Delta.x + m12*Delta.y + m13*Delta.z);					// dG/d(kappa)
	b24 =  fqq*(s*m31 - q*m21);												// dG/d(XL)
	b25 =  fqq*(s*m32 - q*m22);												// dG/d(YL)
	b26 =  fqq*(s*m33 - q*m23);												// dG/d(ZL)

	K	=  photopoint.y - Photo[photoID].Get_IOP().PP.y + fq*s;
}

void CPhotogrammetry::Calculate_Parameters(	InteriorOrientationParameters iop,
											OrientationParameters eop,
											Point2D <double> photopoint,
											Point3D <double> gcp)
{
	//Calculation of Elements of Rotation matrix////////////////////////////

	Get_Elements_of_RotationMatrix(eop.o, eop.p, eop.k);


	//For calculation efficiency////////////////////////////////////////////

	Point3D <double> Delta;
	double temp1, temp2, temp3, temp4, temp5;
	double fq;
	double fqq;	

	Delta.x = gcp.x - eop.Point.x;
	Delta.y = gcp.y - eop.Point.y;
	Delta.z = gcp.z - eop.Point.z;

	q = m31*Delta.x + m32*Delta.y + m33*Delta.z;
	r = m11*Delta.x + m12*Delta.y + m13*Delta.z;
	s = m21*Delta.x + m22*Delta.y + m23*Delta.z;

	fq	= iop.f/q;
	fqq	= iop.f/(q*q);


	//Calculation derivatives//////////////////////////////////////////////////

	//For calculation efficiency
	temp1 = -m33*Delta.y + m32*Delta.z;
	temp2 = cos_p*Delta.x + sin_o*sin_p*Delta.y - cos_o*sin_p*Delta.z;
	temp3 = sin_o*cos_p*Delta.y;
	temp4 = cos_o*cos_p*Delta.z;
	temp5 = sin_p*Delta.x;

	//For x
	b11 =  fqq*(r*(temp1) - q*(-m13*Delta.y + m12*Delta.z));				// dF/d(omega)
	b12 =  fqq*(r*(temp2) - q*(-cos_k*temp5 + cos_k*temp3 - cos_k*temp4));  // dF/d(phi)
	b13 =  -fq*(m21*Delta.x + m22*Delta.y + m23*Delta.z);					// dF/d(kappa)
	b14 =  fqq*(r*m31 - q*m11);												// dF/d(XL)
	b15 =  fqq*(r*m32 - q*m12);												// dF/d(YL)
	b16 =  fqq*(r*m33 - q*m13);												// dF/d(ZL)
	
	J	=  photopoint.x - iop.PP.x + fq*r;

	//For y
	b21 =  fqq*(s*(temp1) - q*(-m23*Delta.y + m22*Delta.z));				// dG/d(omega)
	b22 =  fqq*(s*(temp2) - q*(sin_k*temp5 - sin_k*temp3 + sin_k*temp4));	// dG/d(phi)
	b23 =   fq*(m11*Delta.x + m12*Delta.y + m13*Delta.z);					// dG/d(kappa)
	b24 =  fqq*(s*m31 - q*m21);												// dG/d(XL)
	b25 =  fqq*(s*m32 - q*m22);												// dG/d(YL)
	b26 =  fqq*(s*m33 - q*m23);												// dG/d(ZL)

	K	=  photopoint.y - iop.PP.y + fq*s;
}

int CPhotogrammetry::Space_Resection(CPhoto& photo, int ngcps, Point2D <double>* photopoint, Point3D <double>* gcp, double& sigma, double precision[], double residual[])
{
	int i, j;
	int trow;

	//Number of unoknowns, knowns/////////////////////////////////

	//Number of Observations
	//(x,y)*number_of_points*number_of_photos
	n = 2*ngcps;
	
	//Number of Unknown parameters
	m = 6;

	//For the iteration termination conditions////////////////////

	int iteration_stop  = 0;
	int iteration_count = 0;

	double angle_tolerence		= 0.000001;
	double position_tolerence	= 0.000001;

	double change_sigma_o_2;
	double new_sigma_o_2;
	double old_sigma_o_2		= 9999999999;			

	
	//Save Initial values////////////////////////////////////////
	
	OrientationParameters  EOP;
	OrientationParameters EOP_initial;
	EOP_initial = photo.Get_EOP();


	//Matix settings/////////////////////////////////////////////

	//Matrix resize
	A.Resize(n,m);
	y.Resize(n,1);

	//Set to Zero
	for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			A(i,j)=0;

	while(iteration_stop==0)
	{
		//For the First (Reference) photo/////////////////////////
		//ROP[0]
		
		for(i=0;i<ngcps;i++)
		{
			Calculate_Parameters(photo.Get_IOP(), photo.Get_EOP(), photopoint[i], gcp[i]);
			trow = i*2;

			
			//A matrix/////////////////////////////////////////////

			//for ROP
			A(trow, 0) =  b11;			A(trow+1, 0) =  b21;
			A(trow, 1) =  b12;			A(trow+1, 1) =  b22;
			A(trow, 2) =  b13;			A(trow+1, 2) =  b23;
			A(trow, 3) = -b14;			A(trow+1, 3) = -b24;
			A(trow, 4) = -b15;			A(trow+1, 4) = -b25;
			A(trow, 5) = -b16;			A(trow+1, 5) = -b26;

			//y matrix////////////////////////////////////////////////////////////
			y(trow  , 0) = J;
			y(trow+1, 0) = K;
		}

		AT = A.Transpose();
		N = AT*A;
		Ninv = N.Inverse();
		csi_hat = Ninv*AT*y;
		e_tilda = A*csi_hat-y;	


		//Update////////////////////////////////////////////////////////////////////

		EOP = photo.Get_EOP();

		EOP.o		+= csi_hat(0, 0);
		EOP.p		+= csi_hat(1, 0);
		EOP.k		+= csi_hat(2, 0);
		EOP.Point.x	+= csi_hat(3, 0);
		EOP.Point.y	+= csi_hat(4, 0);
		EOP.Point.z	+= csi_hat(5, 0);

		photo.Set_EOP(EOP);

		iteration_count++;


		//Iteration Termination Conditions///////////////////////////////
		//Wolf, Ghilani, Adjustment Computations, pp. 243~245

		//1. Maximum iteration
		if(iteration_count==10)
			iteration_stop = -1;

		//2. Maximum correction


		//3. Minitoring the Reference Variance
		new_sigma_o_2 = (e_tilda.Transpose()*e_tilda)(0,0)/(n-m);
		change_sigma_o_2 = fabs((new_sigma_o_2-old_sigma_o_2)/old_sigma_o_2);

		if(change_sigma_o_2<0.01)
			iteration_stop = 1;

		old_sigma_o_2 = new_sigma_o_2;
		Precision = Ninv*new_sigma_o_2;

		precision[0] = Precision(0,0);
		precision[1] = Precision(1,1);
		precision[2] = Precision(2,2);
		precision[3] = Precision(3,3);
		precision[4] = Precision(4,4);
		precision[5] = Precision(5,5);

		for(i=0;i<ngcps;i++)
		{
			residual[i*2  ] = e_tilda(i*2,  0);
			residual[i*2+1] = e_tilda(i*2+1,0);
		}
		
		EOP = photo.Get_EOP();
	}
	sigma = new_sigma_o_2;
	return iteration_stop;
}


int CPhotogrammetry::Get_3DPoints(Point3D <double>& initial, double residual[], double& rmse)
{
	int i;
	int count = 0;
	int count_row = 0;
	
	for(i=0;i<nPhotos;i++)
	{
//		Photo[i].Image2CalibratedPhoto();

		if(Photo[i].PhotoPoint.x != NO_VALUE)
		{			
			count++;
		}
	}

	//Number of Unknown parameters
	m = 3;

	//Number of Observations
	n = 2*count;

	if(n<4)
		return -100;

	//For the iteration termination conditions
	iteration_stop		= 0;
	iteration_count		= 0;
	sigma_o_2_hat_old	= 9999999999;			

	//Matrix resize
	A.Resize(n,m);
	y.Resize(n,1);

	while(iteration_stop==0)
	{
		count_row = 0;

		for(i=0;i<nPhotos;i++)
		{
			if(Photo[i].PhotoPoint.x != NO_VALUE)
			{
			//	Calculate_Parameters(Photo[i].Get_IOP(), Photo[i].Get_EOP(), Photo[i].PhotoPoint, initial);
				Calculate_Parameters(i, Photo[i].PhotoPoint, initial);

				A(count_row*2,  0) = b14;	A(count_row*2,  1) = b15;	A(count_row*2,  2) = b16;	
				A(count_row*2+1,0) = b24;	A(count_row*2+1,1) = b25;	A(count_row*2+1,2) = b26;

				y(count_row*2,  0) = J;
				y(count_row*2+1,0) = K;

				count_row++;
			}
		}
	
		AT = A.Transpose();
		N = AT*A;
		Ninv = N.Inverse();
		csi_hat = Ninv*AT*y;
		e_tilda = y - A*csi_hat;	

		sigma_o_2_hat = (e_tilda.Transpose()*e_tilda)(0,0)/(n - m);
//		Precision = Ninv*sigma_o_2_hat;

		initial.x		+= csi_hat(0, 0);
		initial.y		+= csi_hat(1, 0);
		initial.z		+= csi_hat(2, 0);		

		iteration_count++;

		if(iteration_count==10)
			iteration_stop = -1;

		sigma_o_2_hat_change = fabs((sigma_o_2_hat - sigma_o_2_hat_old)/sigma_o_2_hat_old);

		if(sigma_o_2_hat_change<sigma_o_2_hat_change_tolerence)
			iteration_stop = 1;

		sigma_o_2_hat_old = sigma_o_2_hat;
	}

	// Residual in [Pixel]
	if(residual!=NULL)
	{
		count = 0;
		for(i=0;i<nPhotos;i++)
		{
			if(Photo[i].PhotoPoint.x != NO_VALUE)
			{
			//	residual[i*2  ] = sqrt(pow(e_tilda(count*2  , 0),2)*pow(e_tilda(count*2+1, 0),2))/Photo[0].Get_PixelSize();
			//	residual[i*2+1] = residual[i*2  ];
			//	count++;

				residual[i] = sqrt(pow(e_tilda(count*2  , 0),2)+pow(e_tilda(count*2+1, 0),2))/Photo[0].Get_PixelSize();
				count++;
			}
			else
			{
			//	residual[i*2  ] = -1;
			//	residual[i*2+1] = -1; 
				residual[i] = NO_VALUE;			
			}
		}
	}

	rmse = sqrt(sigma_o_2_hat)/Photo[0].Get_PixelSize();

	/*
	count = 0;
	for(i=0;i<nPhotos;i++)
	{
		if(Photo[i].PhotoPoint.x != NO_VALUE)
		{
			residual[i*2  ] = e_tilda(count*2  , 0);
			residual[i*2+1] = e_tilda(count*2+1, 0);
			count++;
		}
		else
		{
			residual[i*2  ] = -1;
			residual[i*2+1] = -1; 
		}
	}
	*/

	return iteration_stop;
}

void CPhotogrammetry::Get_3DPoints(IntersectionParam& IP)
{
	IP.CheckBA = Get_3DPoints(IP.ObjCoordEst, IP.Residual, IP.RMSE, IP.STD);
}

void CPhotogrammetry::Get_3DPoints(IntersectionParam4& IP)
{
	IP.CheckBA = Get_3DPoints(IP.ObjCoordEst, IP.Residual, IP.RMSE, IP.STD);
}

int CPhotogrammetry::Get_3DPoints(Point3D <double>& initial, double residual[], double& rmse, double std[])
{
	int i;
	int count = 0;
	int count_row = 0;
	
	for(i=0;i<nPhotos;i++)
	{
//		Photo[i].Image2CalibratedPhoto();

		if(Photo[i].PhotoPoint.x != NO_VALUE)
		{			
			count++;
		}
	}

	//Number of Unknown parameters
	m = 3;

	//Number of Observations
	n = 2*count;

	if(n<4)
		return -100;

	//For the iteration termination conditions
	iteration_stop		= 0;
	iteration_count		= 0;
	sigma_o_2_hat_old	= 9999999999;			

	//Matrix resize
	A.Resize(n,m);
	y.Resize(n,1);

	while(iteration_stop==0)
	{
		count_row = 0;

		for(i=0;i<nPhotos;i++)
		{
			if(Photo[i].PhotoPoint.x != NO_VALUE)
			{
			//	Calculate_Parameters(Photo[i].Get_IOP(), Photo[i].Get_EOP(), Photo[i].PhotoPoint, initial);
				Calculate_Parameters(i, Photo[i].PhotoPoint, initial);

				A(count_row*2,  0) = b14;	A(count_row*2,  1) = b15;	A(count_row*2,  2) = b16;	
				A(count_row*2+1,0) = b24;	A(count_row*2+1,1) = b25;	A(count_row*2+1,2) = b26;

				y(count_row*2,  0) = J;
				y(count_row*2+1,0) = K;

				count_row++;
			}
		}
	
		AT = A.Transpose();
		N = AT*A;
		Ninv = N.Inverse();
		csi_hat = Ninv*AT*y;
		e_tilda = y - A*csi_hat;	

		sigma_o_2_hat = (e_tilda.Transpose()*e_tilda)(0,0)/(n - m);
//		Precision = Ninv*sigma_o_2_hat;

		initial.x		+= csi_hat(0, 0);
		initial.y		+= csi_hat(1, 0);
		initial.z		+= csi_hat(2, 0);		

		iteration_count++;

		if(iteration_count==10)
			iteration_stop = -1;

		sigma_o_2_hat_change = fabs((sigma_o_2_hat - sigma_o_2_hat_old)/sigma_o_2_hat_old);

		if(sigma_o_2_hat_change<sigma_o_2_hat_change_tolerence)
			iteration_stop = 1;

		sigma_o_2_hat_old = sigma_o_2_hat;
	}

	// Residual in [Pixel]
	if(residual!=NULL)
	{		
		count = 0;

		for(i=0;i<nPhotos;i++)
		{
			if(Photo[i].PhotoPoint.x != NO_VALUE)
			{
				residual[i*2  ] = e_tilda(count*2  , 0)/Photo[i].Get_PixelSize();
				residual[i*2+1] = e_tilda(count*2+1, 0)/Photo[i].Get_PixelSize();

				count++;
			}
			else
			{
				residual[i*2  ] = NO_VALUE;
				residual[i*2+1] = NO_VALUE; 	
			}
		}
	}

	rmse = sqrt(sigma_o_2_hat)/Photo[0].Get_PixelSize();

	if(std!=NULL)
	{
		for(i=0;i<m;i++)
		{
			std[i] = sqrt(sigma_o_2_hat*Ninv(i, i));
		} 
	}

	return iteration_stop;
}

void CPhotogrammetry::EstimateFundamentalMatrix()
{
	if(nPhotos>1)
	{
		int i, j;

		// Mem set
		Initialize_FundamentalMatrix();

		FMatrix = new Matrix <double>* [nPhotos];
		Epipole = new Matrix <double>* [nPhotos];

		for(i=0;i<nPhotos;i++)
		{
			FMatrix[i] = new Matrix <double> [nPhotos];
			Epipole[i] = new Matrix <double> [nPhotos];
		}

		// Camera, camera center, pseudo inverse of camera matrix
		Matrix <double>* P		= new Matrix <double> [nPhotos];
		Matrix <double>* C		= new Matrix <double> [nPhotos];
		Matrix <double>* Ppinv	= new Matrix <double> [nPhotos];

		for(i=0;i<nPhotos;i++)
		{
			P[i]	 = Photo[i].GetCameraMatrix();
			C[i]	 = Photo[i].GetCameraCenterMatrix();
			Ppinv[i] = Photo[i].GetPseudoInverseCameraMatrix();
		}

		// Calculating Fundamental matrices
		for(i=0;i<nPhotos;i++)
		{
			for(j=0;j<nPhotos;j++)
			{
				if(i!=j)
				{
					// From ith to jth camera
					// e of ith cam onto jth cam
					Epipole[i][j] = P[j]*C[i];
					FMatrix[i][j] = (Epipole[i][j].Skew())*P[j]*Ppinv[i];
				}
			}
		}

		delete[] P;
		delete[] C;
		delete[] Ppinv;
	}
}

Matrix <double> CPhotogrammetry::Get_FundamentalMatrix(int fromCAM, int toCAM)
{
	return FMatrix[fromCAM][toCAM];
}

Matrix <double> CPhotogrammetry::Get_Epipole(int fromCAM, int toCAM)
{
	return Epipole[fromCAM][toCAM];
}