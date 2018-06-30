#ifndef PROJECT_POINTS_H
#define PROJECT_POINTS_H
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <vector>

#if __cplusplus
extern "C" {
#endif

using namespace std;
using namespace cv;


void getPoints(double* fiducial, vector<double> &projpts, vector<Point3f> kinectpts);

#if __cplusplus
}
#endif

#endif