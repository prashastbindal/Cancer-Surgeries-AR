// Preparing TCP/IP


#include <WINSOCK2.H>
#include <stdio.h>
#pragma comment(lib,"ws2_32.lib")

#ifdef	ANDROID_TCPIP

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Android TCP/IP Communication
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void main()
{
	WSADATA wsa;
	SOCKET soc_server, soc_client;
	struct sockaddr_in add_server, add_client;
	int c;
	int iResult;
	int iSendResult;
	int port = 85;

	// Initializing WS2_32.dll
	printf("\nInitialising Winsock...");
	iResult = WSAStartup(MAKEWORD(2, 2), &wsa);
	if (iResult != 0)
	{
		printf("WSAStartup failed: %d\n", iResult);
		return 1;
	}
	printf("Initialised.\n");

	// Create a socket
	if ((soc_server = socket(AF_INET, SOCK_STREAM, 0)) == INVALID_SOCKET)
		printf("Could not create socket : %d", WSAGetLastError());
	printf("Socket created.\n");

	// Prepare the sockaddr_in structure
	add_server.sin_family = AF_INET;
	add_server.sin_addr.s_addr = INADDR_ANY;
	add_server.sin_port = htons(port);

	// Bind
	if (bind(soc_server, (struct sockaddr *)&add_server, sizeof(add_server)) == SOCKET_ERROR)
		printf("Bind failed with error code : %d", WSAGetLastError());
	printf("Bind done\n");

	//Listen to incoming connections
	listen(soc_server, 3);

	//Accept and incoming connection
	printf("Waiting for incoming connections...\n");

	c = sizeof(struct sockaddr_in);
	soc_client = accept(soc_server, (struct sockaddr *)&add_client, &c);
	if (soc_client == INVALID_SOCKET)
	{
		printf("accept failed with error code : %d\n", WSAGetLastError());
	}
	printf("Connection accepted\n");




	// Sending data in loop

	char str_android[200];
	int t = 0;
	if ((int)t>0)
	{
		sprintf_s(str_android, "%d_%10.3lf_%10.3lf_%10.3lf_%10.3lf_%10.3lf_%10.3lf\n", 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		iSendResult = send(soc_client, str_android, strlen(str_android), 0);
	}
	else
	{
		sprintf_s(str_android, "%d_%10.3lf_%10.3lf_%10.3lf_%10.3lf_%10.3lf_%10.3lf\n", 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		iSendResult = send(soc_client, str_android, strlen(str_android), 0);
	}


	// Closing communication

	closesocket(soc_client);
	closesocket(soc_server);
	puts("Socket closed");
	WSACleanup();
}

#endif