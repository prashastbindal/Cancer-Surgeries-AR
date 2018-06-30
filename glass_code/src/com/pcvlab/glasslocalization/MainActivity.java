package com.pcvlab.glasslocalization;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.InetAddress;
import java.net.Socket;
import java.util.Date;

import org.opencv.android.BaseLoaderCallback;
import org.opencv.android.CameraBridgeViewBase.CvCameraViewFrame;
import org.opencv.android.CameraBridgeViewBase.CvCameraViewListener2;
import org.opencv.android.LoaderCallbackInterface;
import org.opencv.android.OpenCVLoader;
import org.opencv.core.Core;
import org.opencv.core.Mat;
import org.opencv.core.Point;
import org.opencv.core.Point3;
import org.opencv.core.Scalar;
import org.opencv.core.Size;
import org.opencv.highgui.Highgui;
import org.opencv.imgproc.Imgproc;

import android.app.Activity;
import android.os.Bundle;
import android.os.Environment;
import android.util.Log;
import android.view.KeyEvent;
import android.view.SurfaceView;
import android.view.WindowManager;


public class MainActivity extends Activity implements CvCameraViewListener2
{		
	/** OpenCV View **/
	private OGView 		 		 mOpenCvCameraView;
	
	private Mat scr;
	private Mat dst;

	double[] Fiducial=null;
	byte [] sendFiducial=new byte[64];
	Boolean IsgetModel=false;

	byte[] receivedBuf=new byte[12000000];
	byte[] receivedData=new byte[12000000];
	byte[] temp=new byte[1024*576];
	int dataLenth=0;	
	boolean firstTime=true;
	int receNumber;
	int totalReceived=0;
	
	Socket socket=null;
	InetAddress ia;
	BufferedReader in = null;
	DataOutputStream out=null;
	int num=0;
		
	
	
	public native double[] getPointsWrapper(long addrsrc);
	public native double dataTransforme(long addrsrc,byte[] receved2, int receNumber2, byte[] temp2);
	public native int dataTransforme2Long(byte[] receved3);
	public native void doubleArray2ByteArray(double[] senddoubleArray,byte[] inputdata);
	/** Callback for loading OpenCV **/
	private BaseLoaderCallback mLoaderCallback = new BaseLoaderCallback(this){
		@Override
		public void onManagerConnected(int status){
			switch (status){
			case LoaderCallbackInterface.SUCCESS:
				/** Native JNI Library **/
				System.loadLibrary("pcv_glass_localization");
				mOpenCvCameraView.enableView();
				break;

			default:
				super.onManagerConnected(status);
				break;
			}
		}
	};
	
	@Override
	public void onCreate(Bundle savedInstanceState){
		super.onCreate(savedInstanceState);
		if(!OpenCVLoader.initAsync(OpenCVLoader.OPENCV_VERSION_2_4_8, this, mLoaderCallback)){	
			Log.e("opencv", "Cannot connect to OpenCV Manager");
		}

		this.getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
		this.setContentView(R.layout.surface_view);

		mOpenCvCameraView = (OGView) findViewById(R.id.pcvlab_surface_view);
		mOpenCvCameraView.setVisibility(SurfaceView.VISIBLE);
		mOpenCvCameraView.setMaxFrameSize(1280, 640);
		mOpenCvCameraView.enableFpsMeter();
		mOpenCvCameraView.setCvCameraViewListener(this);
		
		
				
	}

	@Override
	protected void onStop() {
		// TODO Auto-generated method stub
		super.onStop();
		try
		{
			socket.close();
			out.close();
			in.close();
		}catch(IOException e)
		{
			e.printStackTrace();
		}
	};	

	@Override
	public void onPause(){
		super.onPause();

		if (mOpenCvCameraView != null){
			mOpenCvCameraView.disableView();
		}
	}
	
	@Override
	public void onResume(){
		super.onResume();
		OpenCVLoader.initAsync(OpenCVLoader.OPENCV_VERSION_2_4_8, this, mLoaderCallback);
	}

	@Override
	public void onDestroy(){
		super.onDestroy();
		android.os.Process.killProcess(android.os.Process.myPid());
		System.exit(0);
		if (mOpenCvCameraView != null){
			mOpenCvCameraView.disableView();
		}
	}
	
	public void onCameraViewStarted(int width, int height){

	}

	public void onCameraViewStopped(){ 
	}

	public Mat onCameraFrame(CvCameraViewFrame inputFrame)
	{
		//File path = new File(Environment.getExternalStorageDirectory() + "/Images2/");
		 // path.mkdirs();
		//  File file = new File(path, "image"+num+".png");

		//String filename = file.toString();
		dst = inputFrame.rgba();
		//Mat colr;
		//Imgproc.cvtColor(dst, dst, Imgproc.COLOR_BayerBG2BGR);
		//Boolean bool=Highgui.imwrite(filename, dst);
		
		num++;
		if (dst.empty())
			return dst;
		//System.out.println("writing the image");
		//System.out.println(bool);
		System.out.println("find fiducial begin at : " + System.currentTimeMillis() );
		Fiducial = getPointsWrapper(dst.getNativeObjAddr());	
		//System.out.println("get points returns");
		System.out.println("find fiducial end at : " + System.currentTimeMillis());
/*		if (Fiducial==null)
		{
			double[] Ficucial2= new double[8];
			Ficucial2[0]=22.22;
			Ficucial2[1]=22.22;
			Ficucial2[2]=322.22;
			Ficucial2[3]=22.22;
			Ficucial2[4]=322.22;
			Ficucial2[5]=222.22;
			Ficucial2[6]=22.22;
			Ficucial2[7]=122.22;
			Fiducial=Ficucial2;
		}
*/	
		if(Fiducial!=null )
		{
			System.out.println("sending fiducial");
			totalReceived=0;
			receNumber=0;
			
			if (firstTime)
			{
				try 
				{
					ia = InetAddress.getByName("10.0.1.2");
					socket = new Socket(ia, 6000);
				       System.out.println(socket.getSendBufferSize());
				       System.out.println(socket.getReceiveBufferSize());
				     
				       socket.setSendBufferSize(1024*32);
				       socket.setReceiveBufferSize(1024*32);
				       
				       System.out.println(socket.getSendBufferSize());
				       System.out.println(socket.getReceiveBufferSize()); 
				}
				catch (IOException e) 
				{					
					e.printStackTrace();
					System.out.println("Client:Connection error: " + e.toString());
					System.exit(-1);
				}
				try
				{
					out = new DataOutputStream(socket.getOutputStream());
					in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
				} 
				catch (IOException e) 
				{
					e.printStackTrace();
					System.out.println("Client:Connection error: " + e.toString());
					System.exit(-1);
				}
				
				firstTime=false;
					
			
				
			//if(!IsgetModel) //we don't need it this time
			

				System.out.println("communicate begin at : " + System.currentTimeMillis());
				
				doubleArray2ByteArray(Fiducial,sendFiducial);
	
				try 
				{	
					
					out.write(sendFiducial);
	
				} 
				catch (IOException e1) 
				{
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
	
				System.out.println("Client:Connectted!");
				
				try
				{
					receNumber=socket.getInputStream().read(receivedBuf);
//					System.out.println("str_received");		    	
				}

				catch (IOException e1) 
				{
					e1.printStackTrace();
				}	

				if (receNumber!=4)
				{
					try {
						out.writeUTF("wrong");

					} catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}

				}
				if (receNumber==4)
				{
								
					//Log.v("check123","check123");
					//System.out.println("check123");
				dataLenth=dataTransforme2Long(receivedBuf );

				if (dataLenth==0)
				{
					System.out.println("communicate end at : " + System.currentTimeMillis());
				try {
						out.writeUTF("O");
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				if (dataLenth>0)
				{
				try {
					out.writeUTF("OK");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				while (totalReceived<dataLenth)
				{
					try
					{
						receNumber=socket.getInputStream().read(receivedBuf);

						System.arraycopy(receivedBuf,0,receivedData,totalReceived,receNumber);
						totalReceived=totalReceived+receNumber;					    	
					}
					// while(!str_received.equals("End"));

					catch (IOException e1) 
					{
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}	
					Log.i("app","receiving the project points");
				}
				
				System.out.println("communicate end and data process begain at : " + System.currentTimeMillis());
				System.out.println("datalength="+dataLenth);
				double testvar=dataTransforme(dst.getNativeObjAddr(),receivedData,dataLenth,temp );
				System.out.println("datalength="+dataLenth);
				System.out.println("testvar="+testvar);
				/*for (int i=0;i<dataLenth/4;i++)
				{
					int x=receivedBuf[i*4+1]*255+receivedBuf[i*4];
					System.out.println("x="+x);
					
					int y=receivedBuf[i*4+3]*255+receivedBuf[i*4+2];
					System.out.println("y="+y);
					if (temp[x*1024+y]==1)
						temp[x*1024+y]=0;
					else
						temp[x*1024+y]=1;
					
				}*/
				
			}
		
		
				
				//System.out.println("data process end and draw circle at : " + System.currentTimeMillis());
			/*	for (short i=0; i<1024;i++)
					for(short j=0;j<576;j++)
					{
						if(temp[i*576+j]==1)
							Core.circle(dst, new Point(i,j),  5, new Scalar(0,0,255),1);			//marking the points of the tumor
					}*/
			
	
				}
				//System.out.println("draw circle end at : " + System.currentTimeMillis());
			return dst;	
			
			}
			double testvar=dataTransforme(dst.getNativeObjAddr(),receivedData,dataLenth,temp );
			//Boolean bool=Highgui.imwrite(filename, dst);
			
			//System.out.println(bool);
			//num++;
		
		}
	
		return dst;
		
	}
	@Override
	public boolean onKeyUp(int keyCode, KeyEvent event){
		// TODO Auto-generated method stub
		
		return super.onKeyUp(keyCode, event);
	}
}

