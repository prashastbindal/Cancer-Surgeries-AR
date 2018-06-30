LOCAL_PATH := $(call my-dir)


#include $(CLEAR_VARS)

#LOCAL_MODULE    := photogrammetry
#LOCAL_C_INCLUDES := \
#        $(LOCAL_PATH)/photog/include
#LOCAL_CFLAGS := $(LOCAL_C_INCLUDES:%=-I%)
#LOCAL_LDLIBS := -L$(SYSROOT)/usr/lib -ldl

#LOCAL_SRC_FILES := \
#        src\Photo.cpp \
#        src\Photogrammetry.cpp
#include $(BUILD_STATIC_LIBRARY)


include $(CLEAR_VARS)

# Path to OpenCV Android SDK
#include ~/Downloads/OpenCV-2.4.8-android-sdk/sdk/native/jni/OpenCV.mk
include C:/NVPACK/OpenCV-2.4.8.2-Tegra-sdk/sdk/native/jni/OpenCV.mk

# For custom OpenCV modules
LOCAL_MODULE    := pcv_glass_localization
LOCAL_SRC_FILES := pcv_glass_localization.cpp Fiducial.cpp Camera.cpp Project3D.cpp
LOCAL_LDLIBS +=  -llog -ldl

include $(BUILD_SHARED_LIBRARY)
