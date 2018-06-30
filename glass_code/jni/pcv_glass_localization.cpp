#include "pcv_glass_localization.h"
#include"Fiducial.h"
#include"Camera.h"
#include"Project3D.h"

using namespace std;
using namespace cv;

#ifdef __cplusplus
extern "C"
#endif
{


	jbyte temp[1024*576];
	float kinect3dBuf[12000000];
	vector<Point3f> kinect3dCoord;
	Fiducial fid;
	Camera cam;
	Project3D proj;
	//double *fiducial;
//
	jdouble outfiducial[8];
JNIEXPORT jdoubleArray JNICALL Java_com_pcvlab_glasslocalization_MainActivity_getPointsWrapper(JNIEnv* env, jobject, jlong addrsrc)
//int main(int argc, char** argv)
{
//	memcopy(&targetModel2[0],receved2,1000);


	Mat &src = *(Mat*)addrsrc;



	if(!fid.getPoints(src,outfiducial))
		return NULL;


	   jdoubleArray result = env->NewDoubleArray(8);
	 if(result==NULL)
	 {
	   return NULL;
	 }

	 env->SetDoubleArrayRegion(result, 0, 8, outfiducial);
	 return result;


}


/*
JNIEXPORT void JNICALL Java_com_pcvlab_glasslocalization_MainActivity_setSrcPoints(JNIEnv*, jobject, jint jx, jint jy, jint ji)
{
	int i = int(ji);
	int x = int(jx);
	int y = int(jy);
//	src_points[i] = Point2f(x, y);
}
*/




JNIEXPORT jdouble JNICALL Java_com_pcvlab_glasslocalization_MainActivity_dataTransforme(JNIEnv* env, jobject, jlong addrsrc,jbyteArray receved2,  jint receNumber2, jbyteArray templ)
{
//	result.clear();
	Mat &src = *(Mat*)addrsrc;
	int length=receNumber2;


	char* data=(char*)env->GetByteArrayElements(receved2,0);
	//memcpy(&ProjectimagePoints[0],data,length);
	memcpy(&kinect3dBuf[0],data,length);
	//ofstream myfile;
	//myfile.open("print.txt");
	kinect3dCoord.clear();
	Point3f tempCoord;
	for(int k=0;k<length/sizeof(float);k++)
	{

		if(k%3==0)
			tempCoord.x=kinect3dBuf[k];



		if(k%3==1)
					tempCoord.y=kinect3dBuf[k];

		if(k%3==2)
		{
					tempCoord.z=kinect3dBuf[k];
					kinect3dCoord.push_back(tempCoord);

		}

		//myfile<<kinect3dBuf[k]<<"\n";
	}







	vector<Point2f> ProjectimagePoints;


	//cam.calculateProjectPoints(outfiducial, ProjectimagePoints, kinect3dCoord);
	cam.estimateCameraParameters(outfiducial);

	proj.Project3Dto2D(cam, kinect3dCoord, ProjectimagePoints);

	//receNumber2=ProjectimagePoints.size();

	//char* data2=(char*)env->GetByteArrayElements(templ,0);
	//memcpy(&temp[0],data2,2304000);
	//memcpy(&temp[0],data2,1024*576);



//	circle(src, Point(mx,my),  1, Scalar(0,0,255),1);
	for (int i=0;i<ProjectimagePoints.size();i++)
		{


			short y=(short)ProjectimagePoints[i].y;
			short x=(short)ProjectimagePoints[i].x;
		//	cout<<x,y;
		//	circle(src, Point(x,y),  1, Scalar(0,0,255),1);
			//if (temp[y*1024+x]==1)
			//	temp[y*1024+x]=0;
			//else

			circle(src, Point(x,y),  2, Scalar(55,255,80),2);
			//	temp[y*1024+x]=1;
		}
	/*env->SetByteArrayRegion(templ, 0, 1024*576,&temp[0]);



					for (short i=0; i<576;i++)
						for(short j=0;j<1024;j++)
						{
							if(temp[i*1024+j]==1)
								circle(src, Point(j,i),  5, Scalar(0,0,255),1);
						}*/

	return 0;


}
JNIEXPORT jint JNICALL Java_com_pcvlab_glasslocalization_MainActivity_dataTransforme2Long(JNIEnv* env, jobject, jbyteArray receved)
{
	int result2;
	char* data=(char*)env->GetByteArrayElements(receved,0);
	memcpy(&result2,data,4);
			           if(result2==NULL)
			           {
			        	   return NULL;
			           }

			           return result2;

}

JNIEXPORT void JNICALL Java_com_pcvlab_glasslocalization_MainActivity_doubleArray2ByteArray(JNIEnv* env, jobject, jdoubleArray senddoubleArray,jbyteArray datain)
{
	double* data=(double*)env->GetDoubleArrayElements(senddoubleArray,0);

	jbyte jreturnByteArray[64];// =env->NewByteArray(64);

//	jbyte *jby =env->GetByteArrayElements(jreturnByteArray, 0);

	memcpy(&jreturnByteArray[0], data, 64);

	env->SetByteArrayRegion(datain, 0,64, &jreturnByteArray[0]);



}

#ifdef __cplusplus
}
#endif
