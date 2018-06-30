#ifndef PCV_GLASS
#define PCV_GLASS
#include <jni.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <vector>
#include <iostream>
#include <fstream>



using namespace std;
using namespace cv;

#ifdef __cplusplus
extern "C"
#endif
{

JNIEXPORT jdoubleArray JNICALL Java_com_pcvlab_glasslocalization_MainActivity_getPointsWrapper(JNIEnv* env, jobject, jlong addrsrc);
JNIEXPORT jdouble JNICALL Java_com_pcvlab_glasslocalization_MainActivity_dataTransforme(JNIEnv* env, jobject, jlong addrsrc,jbyteArray receved2,  jint receNumber2, jbyteArray templ);


#ifdef __cplusplus
}
#endif
#endif
