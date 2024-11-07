#define  Dim_C      10000                        //维度
#define  PI         3.141592653589793238
#define  Problem    37                         // Total of problems 问题总数
#define  TotalRun   3                          // The times of executions per algorithm 每个算法的执行次数
#define  Maxgen     3*pow(10,6)         // Maximal evaluation number  最大评价数（最大进化次数）
//#define  Maxgen     30000  
#define  Ro_N       100                        // Rotional matrix's order     旋转矩阵的阶

/********************************Network's structural parameters网络结构参数*****************************/
//#define  Total     5                        // Layer 1 T
#define  Total     6                        // Layer 1
int      evaluation = 0;                        // Evaluation number  进化代数
//double   Current_Best_value;                  // Best value at the current iteration  当前迭代的最佳值
double   Learning_rate = 1;                     // Membrance votage  膜电压
double   w1[3][3] = { { -1,-1, -1 }, { -1, -9, -1 }, { -1, -1, -1 } };//convolution weights 卷积权重
double   Diversity;                        //多样性
double  sum_Lon1, sum_Ron1, sum_Uon1, sum_Don1, sum_Loff1, sum_Roff1, sum_Uoff1, sum_Doff1,sum_L, sum_R;
/*****************************Problems' parameter settings问题参数设置******************/
int      problem;                              // Problem number 问题编号
int      dim = 0;                                // Problem's dimension  问题的维度
int      iter = 0;                               // Iterative number   迭代次数
double   a[Dim_C], b[Dim_C];                    // Variants interval bound    变量区间界
double   w[30];                                // Decomposable sub-objective weight   可分解子目标权重
double   z_tans[Dim_C];                        // Rotional transformation    旋转变换
int      Division_no;              //编号      // Variable segment nimber     可变段敏捷器
int      size1[201];               // Variable size for each subsegment   每个子段的可变大小
int      Y_left_number[100], Y_right_number[100]; // Size of each cross subsegment  每个交叉子段的大小
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
//M转弯
/************************Rotational transformation statement旋转变换阶段*******************/
void    Mrot(int, int, int);
void    MulMatrix(double X[][Ro_N], double Y[][Ro_N], int);   //矩阵乘法(MatrixMul)？？？
void    Create_Matrix(int, int n);                           //创建矩阵 
void    Unit_Matrix(int);                                   //单位矩阵

/****************************************Network's modules 网络模块*********************************/
void    Initialization();                 // Create initial population
void    melanization();                   // Stochastic disturbance	  随机扰动阶段
void    Retina_Layer();                   // Retina
void    Lamina_Layer();                   // Lamina
void    Medulla_Layer();                  // Medulla
void    Update_DFV(int);                  // Update

/****************************************Problems' modules 问题模块*********************************/
void    Problem_Parameter(int);   //参数问题
double  fun(indiv X, int);       //以下应该为函数的计算
double  Derivatefun(indiv, int, int, double, double, int, int, int); //衍生函数的计算

//测试函数
double  Schwefel_function(double x[], int);
double  Rastrigin_function(double x[], int);
double  Ackley_function(double x[], int);
double  Elliptic_function(double x[], int); //椭圆函数
double  Sphere_function(double x[], int);
double  Rosenbrock_function(double x[], int);

// Transformation function and its derivate funtion     变换函数及其导数函数
double  Tosz(double);
double  Tosz_derivate(double);
double  Tasy(double, double, int, int);
double  Tasy_derivate(double, double, int, int);
double  sign(double);                                            //符号函数
double  Overlapedsum(indiv, int, int);                           //重叠函数
void    weight_value_combination(int, int);                      //

void    strand(unsigned int seed);

// Performance statistic variables    性能统计变量
double Iter_value[Problem][200000] = { 0, 0 };    //迭代优值？
double Objective_best[Problem][TotalRun];        //目标最佳值
double Max_iter[Problem][TotalRun];             //最优值（最大迭代次数）
double Problem_Min_Iter[Problem];              //问题的最小迭代？
double Problem_Run_Min_Value[Problem][TotalRun];  //问题运行的最小值？
void   Output_Result();                            //输出结果
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

