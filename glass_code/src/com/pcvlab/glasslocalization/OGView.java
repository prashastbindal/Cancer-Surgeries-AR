package com.pcvlab.glasslocalization;

import java.io.File;
import java.io.FileOutputStream;
import java.util.List;

import org.opencv.android.JavaCameraView;

import android.content.Context;
import android.hardware.Camera;
import android.hardware.Camera.PictureCallback;
import android.hardware.Camera.Size;
import android.util.AttributeSet;
import android.util.Log;

public class OGView extends JavaCameraView implements PictureCallback{
	private String mPictureFileName;

	public OGView(Context context, AttributeSet attrs){
		super(context, attrs);
	}

	public List<String> getEffectList(){
		return mCamera.getParameters().getSupportedColorEffects();
	}

	public boolean isEffectSupported(){
		return (mCamera.getParameters().getColorEffect() != null);
	}

	public String getEffect(){
		return mCamera.getParameters().getColorEffect();
	}

	public void setEffect(String effect){
		Camera.Parameters params = mCamera.getParameters();
		params.setColorEffect(effect);
		mCamera.setParameters(params);
	}

	public List<Size> getResolutionList(){
		return mCamera.getParameters().getSupportedPreviewSizes();
	}

	public void setResolution(Size resolution){
		disconnectCamera();
		mMaxHeight = resolution.height;
		mMaxWidth = resolution.width;
		connectCamera(getWidth(), getHeight());
	}

	public Size getResolution(){
		return mCamera.getParameters().getPreviewSize();
	}

	public void takePicture(final String fileName){
		this.mPictureFileName = fileName;

		mCamera.setPreviewCallback(null);
		mCamera.takePicture(null, null, this);
	}
	
	@Override
	public void onPictureTaken(byte[] data, Camera camera){	
		mCamera.startPreview();
		mCamera.setPreviewCallback(this);

		try{
			FileOutputStream fos = new FileOutputStream(mPictureFileName);

			fos.write(data);
			fos.close();
		} 
		catch (java.io.IOException e){
			Log.e("PictureDemo", "Exception in photoCallback", e);
		}
	}

	protected boolean initializeCamera(int width, int height){
		super.initializeCamera(width, height);

		cameraFix();

		return true;
	}

	private void cameraFix(){
		if (mCamera != null){
			Camera.Parameters params = mCamera.getParameters();
			params.setFocusMode(Camera.Parameters.FOCUS_MODE_INFINITY);
			//float fl=params.getFocalLength();
			
			//try{
			
			//FileOutputStream fos = new FileOutputStream("focal.txt");

			//fos.write((int)fl);
			//fos.close();
			//} 
			//catch (java.io.IOException e){
			//	Log.e("PictureDemo", "Exception in photoCallback", e);
			//}
			
			
			//System.out.println("Focal length="+fl);
			params.setPreviewFpsRange(30000, 30000);
			//List<Camera.Size> previewSizes = params.getSupportedPreviewSizes();
			//Camera.Size previewSize =previewSizes.get(5);
			
			
			// Google Glass Only!!!
			//params.setPreviewSize(640,360);
			params.setPreviewSize(1024,576);
			
			mCamera.setParameters(params);		
		}
	}
}
