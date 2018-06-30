#include "SendData.h"
#include <stdio.h>
#include <iostream>

SendData::SendData(void)
{
}


SendData::~SendData(void)
{
}
/*
bool SendData::deliver(void * inputData, long dataLenth, LPVOID  lparam)
{
	if (dataLenth > 4000000)
		return false;
	memcpy(tempData, inputData, dataLenth);
	SOCKET mSocket = (SOCKET)(LPVOID)lparam;

	int sendLenth = 0;
	int copyLenth = 0;
	long totalSend = 0;
	char* sendFirst;
	sprintf(sendFirst, "%ld", dataLenth);
	int bytesRecv = send(mSocket, sendFirst, 4, 0);
	if (bytesRecv == -1)
		return false;

	while (1)
	{
		if (dataLenth - totalSend<sendRate)
			copyLenth = dataLenth - totalSend;
		else
			copyLenth = sendRate;

		if (copyLenth>0)
			memcpy(sendBuf, tempData + totalSend, copyLenth);

		if (copyLenth<0)
		{
			bytesRecv = recv(mSocket, recvBuf, 100, 0);
			if (bytesRecv == -1 || bytesRecv != 2)
				return false;
			else
				return true;
		}
		bytesRecv = recv(mSocket, recvBuf, 100, 0);
		if (bytesRecv != 2 || bytesRecv == -1)
			return false;

		bytesRecv = send(mSocket, sendBuf, copyLenth, 0);
		if (bytesRecv == -1)
			return false;
		totalSend = totalSend + bytesRecv;

	}

}*/