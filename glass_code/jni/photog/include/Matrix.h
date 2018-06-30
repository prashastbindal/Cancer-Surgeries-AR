#ifndef mymatrix_aefwenkfwejriqjnlsdnklfasdlkdsalkfdnqipwjepodsknaksmNkclanwioerhqip3wenfkds
#define mymatrix_aefwenkfwejriqjnlsdnklfasdlkdsalkfdnqipwjepodsknaksmNkclanwioerhqip3wenfkds


#if !defined PI 
#define PI 3.141592653589793
#endif

//////////////////////////////////////////////////////////////////////
//								
//	Matrix Class
//	(Matrix.h)			
//							
//	inha univ.				
//	spacematics lab.			
//	Lee, Young-Jin			
//
//	2002-01-15
//	2008-10-27	
//  2009-05-18 Resize(int row,int col, T value), eye(int row)		
//  2009-05-23 zeros(int row, int col);
//  2012-03-23 MFC removed
//							
//////////////////////////////////////////////////////////////////////

template<class T>
class Matrix  
{
private:	
	int Row;
	int Col;
	T **Element;
	
public:
	// Constructor/Destructor
	Matrix()	
	{
		int i;

		Row = 1;
		Col = 1;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];
	}

	Matrix(int row,int col)
	{	
		int i;

		Row = row;
		Col = col;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];
	}
	
	Matrix(const Matrix<T>& other)
	{
		int i;

		Row=other.Row;
		Col=other.Col;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];

		for(i=0;i<Row;i++)
			for(int j=0;j<Col;j++)
				Element[i][j]=other.Element[i][j];
	}

	Matrix(int row)
	{	
		int i,j;

		Row = row;
		Col = row;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];	

		for(i=0;i<Row;i++)
			for (j=0;j<Col;j++)	
				if(i==j)
					Element[i][j]=1;
				else
					Element[i][j]=0;
	}//for Identity Matrix

	virtual ~Matrix()
	{
		if(Element)
		{
			int i;

			for(i=0;i<Row;i++)
				delete [] Element[i];	
			delete [] Element;
		}
	}	
	
	// Get/Set Element
	T& operator()(int row,int col){return Element[row][col];};

	// Function
	Matrix<T> Transpose(void)
	{
		Matrix<T> ReturnMatrix(Col,Row);

		int i,j;

		for(i=0;i<Row;i++)
			for (j=0;j<Col;j++)	
				ReturnMatrix.Element[j][i]=Element[i][j];	
			
		return ReturnMatrix;
	}

	Matrix<T> Inverse(void)
	{
		if(Row!=Col)
			ShowErrorMsg(1);

		Matrix<T> ReturnMatrix(*this);

		int i,j,k;

		for(k=0;k<Row ;k++)				
		{
			for(j=0;j<Row ;j++)
				if(j!=k) 
					ReturnMatrix.Element[k][j]=ReturnMatrix.Element[k][j]/ReturnMatrix.Element[k][k];
			ReturnMatrix.Element [k][k]= 1.0/ReturnMatrix.Element [k][k];
			
			for(i=0;i<Row;i++)
				if(i!=k)
				{
					for(j=0;j<Row;j++)
						if(j!=k)
							ReturnMatrix.Element[i][j]=ReturnMatrix.Element[i][j]-ReturnMatrix.Element[i][k]*ReturnMatrix.Element[k][j];
					ReturnMatrix.Element [i][k]=-ReturnMatrix.Element [i][k]*ReturnMatrix.Element [k][k];
				}
		}
		
		return ReturnMatrix;
	}

	void Resize(int row,int col)
	{
		int i;

		if(Element)
		{
			for(i=0;i<Row;i++)
				delete [] Element[i];	
			delete [] Element;
		}

		Row = row;
		Col = col;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];
	}

	// 2009-05-17
	void Resize(int row,int col, T value)
	{
		int i,j;

		if(Element)
		{
			for(i=0;i<Row;i++)
				delete [] Element[i];	
			delete [] Element;
		}

		Row = row;
		Col = col;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];

		for(i=0;i<Row;i++)
			for (j=0;j<Col;j++)	
				Element[i][j]=value;
	}

	// 2009-05-17
	void eye(int row)
	{
		int i,j;

		if(Element)
		{
			for(i=0;i<Row;i++)
				delete [] Element[i];	
			delete [] Element;
		}

		Row = row;
		Col = row;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];

		for(i=0;i<Row;i++)
			for (j=0;j<Col;j++)	
				if(i==j)
					Element[i][j]=1.0;
				else
					Element[i][j]=0.0;
	}

	// 2009-05-23
	void zeros(int row, int col)
	{
		int i,j;

		if(Element)
		{
			for(i=0;i<Row;i++)
				delete [] Element[i];	
			delete [] Element;
		}

		Row = row;
		Col = col;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];

		for(i=0;i<Row;i++)
			for (j=0;j<Col;j++)	
				Element[i][j]=0.0;
	}
	
	T Sum()
	{
		T sum = 0;

		int i,j;

		for(i=0;i<Row;i++)
		{
			for(j=0;j<Col;j++)
			{
				sum += Element[i][j];
			}
		}
		return sum;
	}

	T Ave()
	{
		sum = Sum();
		return sum/(Row*Col);
	}

	// Operator
	Matrix<T> operator + (const Matrix<T>& TempMatrix)
	{
		if(Row!=TempMatrix.Row||Col!=TempMatrix.Col)
			ShowErrorMsg(2);

		int i,j;

		Matrix<T> ReturnMatrix(Row,Col);

		for(i=0;i<Row;i++)
			for(j=0;j<Col;j++)
				ReturnMatrix.Element[i][j]=Element[i][j]+TempMatrix.Element[i][j];

		return ReturnMatrix;
	}

	Matrix<T> operator - (const Matrix<T>& TempMatrix)
	{
		if(Row!=TempMatrix.Row||Col!=TempMatrix.Col)
			ShowErrorMsg(3);

		int i,j;

		Matrix<T> ReturnMatrix(Row,Col);

		for(i=0;i<Row;i++)
			for(j=0;j<Col;j++)
				ReturnMatrix.Element[i][j]=Element[i][j]-TempMatrix.Element[i][j];

		return ReturnMatrix;
	}

	Matrix<T> operator * (const Matrix<T>& TempMatrix)
	{
		if(Col!=TempMatrix.Row)
			ShowErrorMsg(4);

		int i,j,k;

		Matrix<T> ReturnMatrix(Row,TempMatrix.Col);

		for(i=0;i<ReturnMatrix.Row;i++)
			for(j=0;j<ReturnMatrix.Col;j++)
			{
				ReturnMatrix.Element[i][j]=0;
				for(k=0;k<Col;k++)
					ReturnMatrix.Element[i][j]+=Element[i][k]*TempMatrix.Element [k][j];
			}

		return ReturnMatrix;
	}

	Matrix<T> operator * (double Scalar)
	{
		int i,j;

		Matrix<T> ReturnMatrix(*this);

		for(i=0;i<Row;i++)
			for(j=0;j<Col;j++)
				ReturnMatrix.Element[i][j]=Element[i][j]*Scalar;
			
		return ReturnMatrix;
	}
	
	Matrix<T> operator = (const Matrix<T>& TempMatrix)
	{
		int i,j;

		if(Element)
		{
			for(i=0;i<Row;i++)
				delete [] Element[i];	
			delete [] Element;
		}

		Row = TempMatrix.Row;
		Col = TempMatrix.Col;

		Element = new T *[Row];
		for(i=0;i<Row;i++)
			Element[i] = new T [Col];

		for(i=0;i<Row;i++)
			for(j=0;j<Col;j++)
				Element[i][j]=TempMatrix.Element[i][j];

		return *this;
	}	

	// Error
	void ShowErrorMsg(int ErrorCode)
	{
		// Removed
	}	

	// Skew symmtric matrix
	// 10/16/2013
	Matrix<T> Skew(void)
	{
		Matrix<T> SSM(3,3);

		SSM(0,0) =  0;
		SSM(0,1) = -Element[2][0];
		SSM(0,2) =  Element[1][0];
		
		SSM(1,0) =  Element[2][0];
		SSM(1,1) =  0;
		SSM(1,2) = -Element[0][0];
		
		SSM(2,0) = -Element[1][0];
		SSM(2,1) =  Element[0][0];
		SSM(2,2) =  0;
		
		return SSM;
	}
};
#endif