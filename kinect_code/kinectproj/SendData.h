#pragma once
#include <stdio.h>
#include <Windows.h>
class SendData
{
public:
	SendData();
	~SendData(void);

private:
	char sendBuf[10000];
	char recvBuf[100];
	char tempData[4000000];

public:
	int sendRate;

public:
	bool deliver(void * inputData, long dataLenth, LPVOID  lparam);


};

