#define  Dim_C      10000                        //ά��
#define  PI         3.141592653589793238
#define  Problem    37                         // Total of problems ��������
#define  TotalRun   3                          // The times of executions per algorithm ÿ���㷨��ִ�д���
#define  Maxgen     3*pow(10,6)         // Maximal evaluation number  �����������������������
//#define  Maxgen     30000  
#define  Ro_N       100                        // Rotional matrix's order     ��ת����Ľ�

/********************************Network's structural parameters����ṹ����*****************************/
//#define  Total     5                        // Layer 1 T
#define  Total     6                        // Layer 1
int      evaluation = 0;                        // Evaluation number  ��������
//double   Current_Best_value;                  // Best value at the current iteration  ��ǰ���������ֵ
double   Learning_rate = 1;                     // Membrance votage  Ĥ��ѹ
double   w1[3][3] = { { -1,-1, -1 }, { -1, -9, -1 }, { -1, -1, -1 } };//convolution weights ���Ȩ��
double   Diversity;                        //������
double  sum_Lon1, sum_Ron1, sum_Uon1, sum_Don1, sum_Loff1, sum_Roff1, sum_Uoff1, sum_Doff1,sum_L, sum_R;
/*****************************Problems' parameter settings�����������******************/
int      problem;                              // Problem number ������
int      dim = 0;                                // Problem's dimension  �����ά��
int      iter = 0;                               // Iterative number   ��������
double   a[Dim_C], b[Dim_C];                    // Variants interval bound    ���������
double   w[30];                                // Decomposable sub-objective weight   �ɷֽ���Ŀ��Ȩ��
double   z_tans[Dim_C];                        // Rotional transformation    ��ת�任
int      Division_no;              //���      // Variable segment nimber     �ɱ��������
int      size1[201];               // Variable size for each subsegment   ÿ���ӶεĿɱ��С
int      Y_left_number[100], Y_right_number[100]; // Size of each cross subsegment  ÿ�������ӶεĴ�С
//double   threshold;	
//double   Average_value=0;
/********************************Network's structural parameters*****************************/

struct   indiv
{
	double var[Dim_C];
	double Value;
	double Yvar[Dim_C];
} x[Total + 1][Total + 1], x1[Total + 1][Total + 1], x2[Total + 1][Total + 1], x3[Total + 1][Total + 1], x4[Total + 1][Total + 1], 
Last_DFly[Total + 1][Total + 1], Memory_value, x_1[Total + 1][Total + 1];

struct  Indiv
{
	double dir[Total + 2][Total + 2];
} m2h_on, m2v_on, m2h_off, m2v_off, Delay_error_on, Delay_error_off, Excite, X_cart_on, X_cart_off,
X_oo_Lon, X_oo_Loff, X_oo_Ron, X_oo_Roff, X_oo_Uon, X_oo_Uoff, X_oo_Don, X_oo_Doff, Cart_D_off, Cart_D_on,
Delay_on, Delay_off, OO_D_on, OO_D_off, Error, Error_on, Error_off;

struct indi
{
	double var[Ro_N][Ro_N];
} Result, Mturn, R[21];
//Mת��
/************************Rotational transformation statement��ת�任�׶�*******************/
void    Mrot(int, int, int);
void    MulMatrix(double X[][Ro_N], double Y[][Ro_N], int);   //����˷�(MatrixMul)������
void    Create_Matrix(int, int n);                           //�������� 
void    Unit_Matrix(int);                                   //��λ����

/****************************************Network's modules ����ģ��*********************************/
void    Initialization();                 // Create initial population
void    melanization();                   // Stochastic disturbance	  ����Ŷ��׶�
void    Retina_Layer();                   // Retina
void    Lamina_Layer();                   // Lamina
void    Medulla_Layer();                  // Medulla
void    Update_DFV(int);                  // Update

/****************************************Problems' modules ����ģ��*********************************/
void    Problem_Parameter(int);   //��������
double  fun(indiv X, int);       //����Ӧ��Ϊ�����ļ���
double  Derivatefun(indiv, int, int, double, double, int, int, int); //���������ļ���

//���Ժ���
double  Schwefel_function(double x[], int);
double  Rastrigin_function(double x[], int);
double  Ackley_function(double x[], int);
double  Elliptic_function(double x[], int); //��Բ����
double  Sphere_function(double x[], int);
double  Rosenbrock_function(double x[], int);

// Transformation function and its derivate funtion     �任�������䵼������
double  Tosz(double);
double  Tosz_derivate(double);
double  Tasy(double, double, int, int);
double  Tasy_derivate(double, double, int, int);
double  sign(double);                                            //���ź���
double  Overlapedsum(indiv, int, int);                           //�ص�����
void    weight_value_combination(int, int);                      //

void    strand(unsigned int seed);

// Performance statistic variables    ����ͳ�Ʊ���
double Iter_value[Problem][200000] = { 0, 0 };    //������ֵ��
double Objective_best[Problem][TotalRun];        //Ŀ�����ֵ
double Max_iter[Problem][TotalRun];             //����ֵ��������������
double Problem_Min_Iter[Problem];              //�������С������
double Problem_Run_Min_Value[Problem][TotalRun];  //�������е���Сֵ��
void   Output_Result();                            //������
double Z[Total + 1][Total + 1], Z_on[Total + 1][Total + 1], Z_off[Total + 1][Total + 1], L[Total + 1][Total + 1];

int Rand(int i) { return rand() % i; }

double sigmod(double a)
{
	return(double(1 / (1 + pow(2.718, -a))));
	//return(sin(a));
}

void IndexSort(double  *p, int length, int * ind_diff)
{
	for (int m = 0; m < length; m++)
	{
		ind_diff[m] = m;
	}
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length - i - 1; j++)
		{
			if (p[j] > p[j + 1])
			{
				double temp = p[j];
				p[j] = p[j + 1];
				p[j + 1] = temp;
				int ind_temp = ind_diff[j];
				ind_diff[j] = ind_diff[j + 1];
				ind_diff[j + 1] = ind_temp;
			}
		}
	}
}

void weight_value_combination(int Division_no, int problem)
{
	int i = 0, j = 0;
	for (j = 0; j < Division_no; j++)
	{
		double sum = 0;
		for (i = 0; i < 12; i++)
			sum += (rand() % 1001) / 1000.0;
		sum = 1 * (sum - 6);
		w[j] = pow(10, 3 * sum);
	}
}

