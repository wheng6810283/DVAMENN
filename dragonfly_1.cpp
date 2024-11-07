#include<iostream>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include <fstream>
#include<string>
#include"dragonfly_1.h"
#include<algorithm>
#include<cstdio>
#include<ctime>
#include <vector>
using namespace std;

/*����������*/
void main()
{
	srand((unsigned)time(NULL));
	char fileName[500];
	ofstream fp;
	//ofstream mycout("F:\\ʵ����\\2013f35.txt");
	for (problem = 20; problem < 21; problem++)//F16,                                      //������11������35
	{
		Problem_Parameter(problem);//���������Problem��

		//for (int div = 0; div < Division_no; div++)
		//	weight_value_combination(Division_no, problem);                             //Ȩ��ֵ��ϣ�
		//evaluation = 0;
		sprintf_s(fileName, "F:\\ʵ����\\F%d.txt", problem);
		for (int run = 1; run < 2; run++)//ÿ����������д�����һ����25��
		{
			fp.open(fileName, ios::app);

			cout << endl << "Problem=" << problem << endl;
			//cout <<endl   //ʹ�������л��еĲ���
			srand((unsigned)time(NULL));
			evaluation = 0;
			Initialization();
			cout << endl << "Problem=" << problem << ",Run=" << run << endl;
			//cout << "Problem=" << problem << ",Run=" << run << endl;
			iter = 0;                                                                //  �����Ĺ��� 

			while (evaluation < Maxgen)
			{
				if (/*evaluation == Total * Total ||*/ evaluation % 1000 == 0 /*|| evaluation == Maxgen*/)
				{

					fp << Memory_value.Value << " ";
					if (Memory_value.Value == 0) {
						break;
					}
				}
				/*melanization();*/
				Retina_Layer();
				Lamina_Layer();
				Medulla_Layer();
				Update_DFV(problem);
				Iter_value[problem][iter] += Memory_value.Value;                              //����ֵ+----����ֵ
				iter++;
				//evaluation++;
			}
			fp << endl;
			fp.close();
			Problem_Run_Min_Value[problem][run] = Memory_value.Value;                         //�������е���Сֵ
			Max_iter[problem][run] = iter;	/**/                                              //��������Ϊ�����ĵ�����

		}
		//mycout.close();
	}
	/*Output_Result();*/
	system("pause");
}


/*��ʼ��*/
void Initialization()
{
	/*1.����ṹ�ĳ�ʼ��*/
	int i = 0, j = 0, k = 0, l = 0;
	for (i = 0; i <= Total + 1; i++)     //�����߼��㸳��ֵ
	{
		Z[i][0] = 0, Z[i][Total + 1] = 0;
		L[i][0] = 0, L[i][Total + 1] = 0;
		X_oo_Lon.dir[i][0] = 0, X_oo_Lon.dir[i][Total + 1] = 0;
		X_oo_Ron.dir[i][0] = 0, X_oo_Ron.dir[i][Total + 1] = 0;
		X_oo_Uon.dir[i][0] = 0, X_oo_Uon.dir[i][Total + 1] = 0;
		X_oo_Don.dir[i][0] = 0, X_oo_Don.dir[i][Total + 1] = 0;
		Cart_D_on.dir[i][0] = 0, Cart_D_on.dir[i][Total + 1] = 0;

		X_oo_Loff.dir[i][0] = 0, X_oo_Loff.dir[i][Total + 1] = 0;
		X_oo_Roff.dir[i][0] = 0, X_oo_Roff.dir[i][Total + 1] = 0;
		X_oo_Uoff.dir[i][0] = 0, X_oo_Uoff.dir[i][Total + 1] = 0;
		X_oo_Doff.dir[i][0] = 0, X_oo_Doff.dir[i][Total + 1] = 0;
		Cart_D_off.dir[i][0] = 0, Cart_D_off.dir[i][Total + 1] = 0;

	}
	for (j = 0; j <= Total + 1; j++) //�����߼��㸳��ֵ
	{
		Z[0][j] = 0, Z[Total + 1][j] = 0;
		L[0][j] = 0, L[Total + 1][j] = 0;
		X_oo_Lon.dir[0][j] = 0, X_oo_Lon.dir[Total + 1][j] = 0;
		X_oo_Ron.dir[0][j] = 0, X_oo_Ron.dir[Total + 1][j] = 0;
		X_oo_Uon.dir[0][j] = 0, X_oo_Uon.dir[Total + 1][j] = 0;
		X_oo_Don.dir[0][j] = 0, X_oo_Don.dir[Total + 1][j] = 0;
		Cart_D_on.dir[0][j] = 0, Cart_D_on.dir[Total + 1][j] = 0;


		X_oo_Loff.dir[0][j] = 0, X_oo_Loff.dir[Total + 1][j] = 0;
		X_oo_Roff.dir[0][j] = 0, X_oo_Roff.dir[Total + 1][j] = 0;
		X_oo_Uoff.dir[0][j] = 0, X_oo_Uoff.dir[Total + 1][j] = 0;
		X_oo_Doff.dir[0][j] = 0, X_oo_Doff.dir[Total + 1][j] = 0;
		Cart_D_off.dir[0][j] = 0, Cart_D_off.dir[Total + 1][j] = 0;

	}

	for (i = 0; i <= Total; i++) {
		for (j = 0; j <= Total; j++) {
			X_oo_Lon.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;          //������Ϊ-1��1֮������������8*8�Ĵ�С����ֵ
			X_oo_Loff.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;

			X_oo_Ron.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;
			X_oo_Roff.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;

			X_oo_Uon.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;
			X_oo_Uoff.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;

			X_oo_Don.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;
			X_oo_Doff.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;
		}
	}
	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++)
		{
			Cart_D_on.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;

			Cart_D_off.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;

		}
	}

	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++) {
			Delay_error_on.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;
			Delay_error_off.dir[i][j] = -1 + 2 * (rand() % 1000) / 1000.0;
		}
	}

	/*2.�����������������ĳ�ʼ��*/
	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++)
		{
			for (k = 1; k < dim; k++)
			{
				x[i][j].var[k] = a[k] + (b[k] - a[k])*(rand() % 1001) / 1000.0;
				//  cout<<x[i][j].var[k]<<",";
			}
			x[i][j].Value = fun(x[i][j], problem);                              //���㺯��ֵ

																  //x[i][j]��x[i][j].Value�����ֵ�Ƿ���� ����һ�£���
			x_1[i][j] = x[i][j];                                                //�����x1�ƺ�û�õ�
			Last_DFly[i][j] = x[i][j];                                          //��ǰ֡��ֵ����һ֡ 
		}
	}

	Memory_value = x[1][1];                                        //��¼��ʼ�ļ���ֵ������ļ���ֵ�Ǵ�[1][1]��ʼ�� ѭ���Ǵ�1��1��ʼ�� ϸ�����⣨��Ϊ�������������˳�ֵ��0��1-8����ʵ�ʵ�ֵ��

	/* for(k=0;k<dim;k++)
	cout<<x[0][0].var[k]<<",";*/

	cout << endl << endl;                                          //���в���
}

void Retina_Layer()
{
	int i = 0, j = 0, k = 0, l = 0, u = 0, v = 0;
	double weight = 0;
	double sum_on = 0, sum_off = 0, Sum = 0;
	double Max = -1000, min = pow(10, 30), max = -pow(10, 8);
	double mean = 0;
	int P = 0;
	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++)
		{
			Z[i][j] = Last_DFly[i][j].Value - x[i][j].Value;                //����ƫ������֡��֣����˴���Last_Flyδ��������xΪ���������õ���

			if (min > x[i][j].Value)
			{
				min = x[i][j].Value;                                                              //��¼����ֵ����Сֵ������¼�����
				u = i, v = j;
			}
			if (Max < fabs(Z[i][j]))                                                                //��¼��ƫ������ֵ
				Max = fabs(Z[i][j]);
		}
		for (k = 1; k < dim; k++) {
			if (Memory_value.Value > x[u][v].Value)                                  //�����ֵ�ȼ�¼ֵС�����������¼��ֵ(��¼��Сֵ��λ��)
			{
				Memory_value = x[u][v];
				/*else
					a[k] + (b[k] - a[k])*(rand() % 1001) / 1000.0;
				Memory_value = x[u][v];*/
			}
		}
		for (i = 1; i <= Total; i++) {
			for (j = 1; j <= Total; j++) {
				for (k = 0; k < dim; k++) {
					if ((x[u][v].var[k] < a[k]) || (x[u][v].var[k] > b[k]))
					{
						Memory_value.var[k] = a[k] + (b[k] - a[k])*(rand() % 1001) / 1000.0;
					}
				}
			}
		}
	}

	if (evaluation % 1000 == 0) {
		cout << "Evaluation=" << evaluation << ",Memory_value=" << Memory_value.Value << endl;
	}


}

void Lamina_Layer()
{
	int i = 0, j = 0, k = 1, l = 0, u = 0, v = 0;
	double temp = 0.0, weight1 = 0.0;
	double Max = -1000, min = pow(10, 30), max = -pow(10, 8);
	for (i = 1; i <= Total; i++)
		for (j = 1; j <= Total; j++)
		{
			for (u = 0; u < 3; u++)
				for (v = 0; v < 3; v++)
				{
					weight1 += w1[u][v];
					temp += w1[u][v] * Z[i + u - 1][j + v - 1] / (Max + 0.1);
					L[i][j] = temp / weight1;
					/*cout << "L[i][j]=" << L[i][j] << endl;*/
				}
			/*L[i][j] = temp;*/
			/*cout << "L[i][j]=" << L[i][j] << endl;*/
			if (L[i][j] < 0)
			{
				X_cart_off.dir[i][j] = abs(L[i][j]);
				/*cout << "L2[i][j]=" << X_cart_off.dir[i][j] << endl;*/
			}
			else if (L[i][j] > 0)
			{
				X_cart_on.dir[i][j] = L[i][j];
				/*cout << "L1[i][j]=" << X_cart_off.dir[i][j] << endl;*/
			}
		}

	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++) {
			Cart_D_on.dir[i][j] = X_cart_on.dir[i][j];           //��ǰ��ֵ��Ϊ��һ֡�����
			Cart_D_off.dir[i][j] = X_cart_off.dir[i][j];
			/*cout << "L2[i][j]=" << Cart_D_off.dir[i][j] << endl;*/
		}
	}

	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++)
		{
			double sum_LEon = 0, sum_LIon = 0, sum_LEoff = 0, sum_LIoff = 0;
			for (k = 0; k < 2; k++)
			{
				sum_LEon += X_cart_on.dir[i][j] * Cart_D_on.dir[i + k][j];
				sum_LIon += X_cart_on.dir[i + k][j] * Cart_D_on.dir[i][j];            //��һ���У�������������ʱ֮�ͣ���Ӧ�����ĸ����֣�������һ�����֣�
				sum_LEoff += X_cart_off.dir[i][j] * Cart_D_off.dir[i + k][j];
				sum_LIoff += X_cart_off.dir[i + k][j] * Cart_D_off.dir[i][j];
			}
			X_oo_Lon.dir[i][j] = sum_LEon - 0.4*sum_LIon;                            //sum/8Ϊ8��С��ƽ���ļ������ֵ
			X_oo_Loff.dir[i][j] = sum_LEoff - 0.4*sum_LIoff;
			/*cout << "L1[i][j]=" << X_oo_Lon.dir[i][j] << endl;
			cout << "L2[i][j]=" << X_oo_Loff.dir[i][j] << endl;*/

		}
	}
	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++)
		{
			double sum_REon = 0, sum_RIon = 0, sum_REoff = 0, sum_RIoff = 0;
			for (k = 0; k < 2; k++)
			{
				sum_REon += X_cart_on.dir[i][j] * Cart_D_on.dir[i - k][j];
				sum_RIon += X_cart_on.dir[i - k][j] * Cart_D_on.dir[i][j];            //��һ���У�������������ʱ֮�ͣ���Ӧ�����ĸ����֣�������һ�����֣�
				sum_REoff += X_cart_off.dir[i][j] * Cart_D_off.dir[i - k][j];
				sum_RIoff += X_cart_off.dir[i - k][j] * Cart_D_off.dir[i][j];
			}
			X_oo_Ron.dir[i][j] = sum_REon - 0.4*sum_RIon;                            //sum/8Ϊ8��С��ƽ���ļ������ֵ
			X_oo_Roff.dir[i][j] = sum_REoff - 0.4*sum_RIoff;
		}
	}

	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++)
		{
			double sum_UEon = 0, sum_UIon = 0, sum_UEoff = 0, sum_UIoff = 0;
			for (k = 0; k < 2; k++)
			{
				sum_UEon += X_cart_on.dir[i][j] * Cart_D_on.dir[i][j + k];
				sum_UIon += X_cart_on.dir[i][j + k] * Cart_D_on.dir[i][j];            //��һ���У�������������ʱ֮�ͣ���Ӧ�����ĸ����֣�������һ�����֣�
				sum_UEoff += X_cart_off.dir[i][j] * Cart_D_off.dir[i][j + k];
				sum_UIoff += X_cart_off.dir[i][j + k] * Cart_D_off.dir[i][j];
			}
			X_oo_Uon.dir[i][j] = sum_UEon - 0.4*sum_UIon;                            //sum/8Ϊ8��С��ƽ���ļ������ֵ
			X_oo_Uoff.dir[i][j] = sum_UEoff - 0.4*sum_UIoff;
		}
	}
	for (i = 1; i <= Total; i++) {
		for (j = 1; j <= Total; j++)
		{
			double sum_DEon = 0, sum_DIon = 0, sum_DEoff = 0, sum_DIoff = 0;
			for (k = 0; k < 2; k++)
			{
				sum_DEon += X_cart_on.dir[i][j] * Cart_D_on.dir[i][j - k];
				sum_DIon += X_cart_on.dir[i][j - k] * Cart_D_on.dir[i][j];            //��һ���У�������������ʱ֮�ͣ���Ӧ�����ĸ����֣�������һ�����֣�
				sum_DEoff += X_cart_off.dir[i][j] * Cart_D_off.dir[i][j - k];
				sum_DIoff += X_cart_off.dir[i][j - k] * Cart_D_off.dir[i][j];
			}
			X_oo_Don.dir[i][j] = sum_DEon - 0.4*sum_DIon;                            //sum/8Ϊ8��С��ƽ���ļ������ֵ
			X_oo_Doff.dir[i][j] = sum_DEoff - 0.4*sum_DIoff;
		}
	}

}

void Medulla_Layer()
{
	int i = 0, j = 0;
	double sum_Lon = 0, sum_Loff = 0;
	double sum_Ron = 0, sum_Roff = 0;
	double sum_Don = 0, sum_Doff = 0;
	double sum_Uon = 0, sum_Uoff = 0;
	for (i = 1; i <= Total; i++)
	{

		for (j = 1; j <= Total; j++)
		{
			sum_Lon += X_oo_Lon.dir[i][j];                                //�����Ƿ���Ҫ���ݸ���һ֡���� ����Ҳ���������forѭ������
			sum_Loff += X_oo_Loff.dir[i][j];
		}

	}
	sum_Lon1 = 2 / (1 + sigmod(sum_Lon)) - 1;
	sum_Loff1 = 2 / (1 + sigmod(sum_Loff)) - 1;
	/*cout<< "sum_Lon="<< sum_Lon1[i] <<endl;*/
	for (i = 1; i <= Total; i++)
	{

		for (j = 1; j <= Total; j++)
		{
			sum_Ron += X_oo_Ron.dir[i][j];                                //�����Ƿ���Ҫ���ݸ���һ֡���� ����Ҳ���������forѭ������
			sum_Roff += X_oo_Roff.dir[i][j];
		}

	}
	sum_Ron1 = 2 / (1 + sigmod(sum_Ron)) - 1;
	sum_Roff1 = 2 / (1 + sigmod(sum_Roff)) - 1;

	for (i = 1; i <= Total; i++)
	{

		for (j = 1; j <= Total; j++)
		{
			sum_Uon += X_oo_Uon.dir[i][j];                                //�����Ƿ���Ҫ���ݸ���һ֡���� ����Ҳ���������forѭ������
			sum_Uoff += X_oo_Uoff.dir[i][j];
		}
	}
	sum_Uon1 = 2 / (1 + sigmod(sum_Uon)) - 1;
	sum_Uoff1 = 2 / (1 + sigmod(sum_Uoff)) - 1;
	for (i = 1; i <= Total; i++)
	{

		for (j = 1; j <= Total; j++)
		{
			sum_Don += X_oo_Don.dir[i][j];                                //�����Ƿ���Ҫ���ݸ���һ֡���� ����Ҳ���������forѭ������
			sum_Doff += X_oo_Doff.dir[i][j];
		}

	}//ÿһ���ˮƽ�������ֱ����֮�͵��ܺ� //��Ϊ�����ѧϰ��//˫�����м�����������Ҳ�����ò�ͬ������Ч����
	/*cout << "1213Learning=" << Learning_rate << endl;*/
	sum_Don1 = 2 / (1 + sigmod(sum_Don)) - 1;
	sum_Doff1 = 2 / (1 + sigmod(sum_Doff)) - 1;

	sum_L = (sum_Lon1 + sum_Ron1 + sum_Uon1 + sum_Don1) / 4;
	sum_R = (sum_Loff1+ sum_Roff1+ sum_Uoff1+ sum_Doff1) / 4;

}


void Update_DFV(int problem)                                                  //���¹�ʽ��������������Ȼ���ڽ��и���
{
	int i = 0, j = 0, k = 0;
	int d[64] = { 0 };
	for (i = 0; i < 64; i++)
	{
		d[i] = i;
	}

	for (i = 1; i <= Total; i++)
		for (j = 1; j <= Total; j++)
		{
			Last_DFly[i][j] = x[i][j];                                              //��ǰֵ����ֵ����һ��ֵ
		}

	for (i = 1; i <= Total; i++)
		for (j = 1; j <= Total; j++)
		{
			srand(unsigned(time(0)));
			random_shuffle(d, d + 64, Rand);//������Ԫ�ش��ң���ÿ�ζ���ͬһ�ִ���˳��
			int e = 0.0, e1 = 0.0, e2 = 0.0;
			e = d[0] / 8; e1 = d[1] / 8; e2 = d[2] / 8;
			int f = 0, f1 = 0, f2 = 0;
			f = d[0] % 8; f1 = d[1] % 8; f2 = d[2] % 8;
			double r = -1 + 2 * (rand() % 1001) / 1000.0;                            //���������-1��1֮��������
			for (int k = 0; k < dim; k++)
			{
				double s = 1 * (rand() % 1001) / 1000.0;                                   //������Ϊ0��1֮��������                   

				double u = 1 * (rand() % 1001) / 1000.0; 
				
				
				if(sum_L> sum_R)
				{
					//x1[i][j].var[k] += x[i][j].var[k] + r * sum_Lon1[i] * (Memory_value.var[k] - x[e][f].var[k]) + r * sum_Loff1[i] * (x[e1][f1].var[k] - x[e2][f2].var[k]);
					x1[i][j].var[k] += Memory_value.var[k] + r * (sum_L / sum_R) * (x[i][j].var[k] - sum_L * Memory_value.var[k]);
					x2[i][j].var[k] += x[i][j].var[k] + r * (sum_R / sum_L) * (x[i][j].var[k] - sum_R * Memory_value.var[k]);
				}

				
				
				if (sum_L < sum_R)
				{
					//x1[i][j].var[k] += x[i][j].var[k] + r * sum_Lon1[i] * (Memory_value.var[k] - x[e][f].var[k]) + r * sum_Loff1[i] * (x[e1][f1].var[k] - x[e2][f2].var[k]);
					x1[i][j].var[k] += Memory_value.var[k] + r * (sum_R / sum_L) * (x[i][j].var[k] - sum_L * Memory_value.var[k]);
					x2[i][j].var[k] += x[i][j].var[k] + r * (sum_L / sum_R) * (x[i][j].var[k] - sum_R * Memory_value.var[k]);
				}
				
				
				
				if (sum_L = sum_R )
				{
					//x1[i][j].var[k] += x[i][j].var[k] + r * sum_Lon1[i] * (Memory_value.var[k] - x[e][f].var[k]) + r * sum_Loff1[i] * (x[e1][f1].var[k] - x[e2][f2].var[k]);
					x1[i][j].var[k] += (x1[i][j].var[k] + x2[i][j].var[k]) / 2 + r * (sum_L / (sum_L + sum_R)) * (x1[i][j].var[k] - Memory_value.var[k]) + r * (sum_L / (sum_L + sum_R)) * (x2[i][j].var[k] - Memory_value.var[k]);
				}


				
				x1[i][j].var[k] += x[i][j].var[k] + r * sum_Lon1 * (Memory_value.var[k] - x[e][f].var[k]) + r * sum_Loff1 * (x[e1][f1].var[k] - x[e2][f2].var[k]);
				x2[i][j].var[k] += x[i][j].var[k] + r * sum_Ron1 * (Memory_value.var[k] - x[e][f].var[k]) + r * sum_Roff1 * (x[e1][f1].var[k] - x[e2][f2].var[k]);
				x3[i][j].var[k] += x[i][j].var[k] + r * sum_Uon1 * (Memory_value.var[k] - x[e][f].var[k]) + r * sum_Uoff1 * (x[e1][f1].var[k] - x[e2][f2].var[k]);
				x4[i][j].var[k] += x[i][j].var[k] + r * sum_Don1 * (Memory_value.var[k] - x[e][f].var[k]) + r * sum_Doff1 * (x[e1][f1].var[k] - x[e2][f2].var[k]);

				if ((x1[i][j].var[k] < a[k]) || (x1[i][j].var[k] > b[k]))                     //���ж��Ƿ�Խ�磬ûԽ�����ԭ����ֵ��Խ��ʱ�������ٸ÷�Χ�����������
				{
					x1[i][j].var[k] = a[k] + (b[k] - a[k]) * (rand() % 1001) / 1000.0;
				}
				if ((x2[i][j].var[k] < a[k]) || (x2[i][j].var[k] > b[k]))
				{
					x2[i][j].var[k] = a[k] + (b[k] - a[k]) * (rand() % 1001) / 1000.0;
				}
				if ((x3[i][j].var[k] < a[k]) || (x3[i][j].var[k] > b[k]))
				{
					x3[i][j].var[k] = a[k] + (b[k] - a[k]) * (rand() % 1001) / 1000.0;
				}
				if ((x4[i][j].var[k] < a[k]) || (x4[i][j].var[k] > b[k]))
				{
					x4[i][j].var[k] = a[k] + (b[k] - a[k]) * (rand() % 1001) / 1000.0;
				}
				/*cout << "x1 =" << x1[i][j].var[k] << endl;
				cout << "x2 = " << x2[i][j].var[k] << endl;
				cout << "x3 = " << x3[i][j].var[k] << endl;
				cout <<"x4 = " << x4[i][j].var[k] << endl;
				cout << "x = " << x[i][j].var[k] << endl;*/
			}
			x1[i][j].Value = fun(x1[i][j], problem);
			x2[i][j].Value = fun(x2[i][j], problem);
			x3[i][j].Value = fun(x3[i][j], problem);
			x4[i][j].Value = fun(x4[i][j], problem);
			/*cout << "x1.Value =" << x1[i][j].Value << endl;
			cout << "x2.Value = " << x2[i][j].Value << endl;
			cout << "x3.Value = " << x3[i][j].Value << endl;
			cout << "x4.Value = " << x4[i][j].Value << endl;*/

			double arr[] = { x1[i][j].Value, x2[i][j].Value, x3[i][j].Value, x4[i][j].Value };
			int ind[4] = { 0 };
			IndexSort(arr, 4, ind);
			if (ind[0] = 0)
				x[i][j] = x1[i][j];
			else if (ind[0] = 1)
				x[i][j] = x2[i][j];
			else if (ind[0] = 2)
				x[i][j] = x3[i][j];
			else
				x[i][j] = x4[i][j];
			/*x[i][j] = x1[i][j];  */                                       //�������һ��ģ�飬����ȫ�����Ž�
			x[i][j].Value = fun(x[i][j], problem);                                         //������Ӧ��ֵ  
			/*cout << "x[i][j] = " << x[i][j].var[k] << endl;*/
			evaluation++;
		}//1
}


//����Ϊ��������ļܹ�����





//�����ļ���
double fun(indiv X, int problem)
{
	int i = 0, j = 0, k = 0, p = 0, sum0 = 0;
	double sum = 0, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, Sum = 0, z0[1001];

	if (problem == 0)
	{
		for (j = 0; j < dim; j++)
			sum += pow(X.var[j] - 1, 2);

		/*for (j = 0; j<dim; j++)
		sum += pow(10, 6 * j / (dim - 1.0))*pow(Tosz(X.var[j]), 2);        */        //ToszΪ����������ĺ��������
	}
	//
	if (problem == 1)
	{
		for (j = 0; j < dim; j++)
		{
			double z = pow(10, 0.5*(j*1.0) / (dim))*Tasy(Tosz(X.var[j]), 0.2, j, dim);  //����������ĺ���
			sum += pow(z, 2) - 10 * cos(2 * PI*z) + 10;
		}
	}

	if (problem == 2)
	{
		for (j = 0; j < dim; j++)
		{
			double z = pow(10, 0.5*(j*1.0) / (dim - 1))*Tasy(Tosz(X.var[j]), 0.2, j - 0, dim);
			sum1 += pow(z, 2);
			sum2 += cos(2 * PI*z);
		}
		sum = -20 * exp(-0.2*sqrt(sum1 / dim)) - exp(sum2 / dim) + 20 + 2.718281828459045;
	}

	if (problem == 3)
	{
		sum = 0, sum0 = 0;
		for (i = 0; i < 7; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			p = 0;
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z0[p] = Tosz(z_tans[j]);
				p++;
			}
			double x = Elliptic_function(z0, size1[i]);
			sum += w[i] * x;
			sum0 += size1[i];
		}
		p = 0;
		for (j = sum0; j < sum0 + size1[7]; j++)
		{
			z0[p] = Tosz(X.var[j]);
			p++;
		}
		sum += Elliptic_function(z0, size1[7]);
	}


	if (problem == 4)
	{
		int i = 0;
		for (i = 0; i < 7; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			p = 0;

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z0[p] = pow(10, 0.5*((j - sum0)*1.0) / (size1[i] - 1))*Tasy(Tosz(z_tans[j]), 0.2, j - sum0, size1[i]);
				p++;
			}
			sum0 += size1[i];
			sum += w[i] * Rastrigin_function(z0, size1[i]);
		}
		p = 0;
		for (j = sum0; j < sum0 + size1[7]; j++)
		{
			z0[p] = pow(10, 0.5*((j - sum0)*1.0) / (size1[7] - 1))*Tasy(Tosz(X.var[j]), 0.2, j - sum0, size1[7]);
			p++;
		}
		sum += Rastrigin_function(z0, size1[7]);
	}
	//
	if (problem == 5)
	{
		sum = 0;
		for (i = 0; i < 7; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			p = 0;
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z0[p] = pow(10, 0.5*((j - sum0)*1.0) / (size1[i] - 1))*Tasy(Tosz(z_tans[j]), 0.2, j - sum0, size1[i]);
				p++;
			}
			sum += w[i] * Ackley_function(z0, size1[i]);
			sum0 += size1[i];
		}
		p = 0;
		for (j = sum0; j < sum0 + size1[7]; j++)
		{
			z0[p] = pow(10, 0.5*((j - sum0)*1.0) / (size1[7] - 1))*Tasy(Tosz(X.var[j]), 0.2, j - sum0, size1[7]);//No Rotation
			p++;
		}
		sum += Ackley_function(z0, size1[7]);
	}
	//
	if (problem == 6)
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no - 1; i++)
		{

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			p = 0;
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z0[p] = Tasy(Tosz(z_tans[j]), 0.2, j - sum0, size1[i]);
				p++;
			}
			sum0 += size1[i];
			sum += w[i] * Schwefel_function(z0, size1[i]);

		}
		p = 0;
		for (j = sum0; j < sum0 + size1[Division_no - 1]; j++)
		{
			double Sum = 0;
			double z = Tasy(Tosz(X.var[j]), 0.2, j - sum0, size1[Division_no - 1]);
			Sum += pow(z, 2);
		}
		sum += Sum;
	}
	//
	if (problem == 7)
	{
		int j = 0, sum0 = 0;
		sum = 0;
		for (i = 0; i < Division_no; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			p = 0;
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z0[p] = Tosz(z_tans[j]);
				p++;
			}
			sum += w[i] * Elliptic_function(z0, size1[i]);
			sum0 += size1[i];
		}
	}
	//
	if (problem == 8)
	{
		for (i = 0; i < Division_no; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			p = 0;
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z0[p] = pow(10, 0.5*((j - sum0)*1.0) / (size1[i] - 1))*Tasy(Tosz(z_tans[j]), 0.2, j - sum0, size1[i]);
				p++;
			}
			sum += w[i] * Rastrigin_function(z0, size1[i]);
			sum0 += size1[i];
		}
	}
	//
	if (problem == 9)
	{
		int sum0 = 0;
		for (i = 0; i < Division_no; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}

			p = 0;
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z0[p] = pow(10, 0.5*((j - sum0)*1.0) / (size1[i] - 1))*Tasy(Tosz(z_tans[j]), 0.2, j - sum0, size1[i]);
				p++;
			}
			sum += w[i] * Ackley_function(z0, size1[i]);
			sum0 += size1[i];
		}
	}

	//
	if (problem == 10)
	{
		for (i = 0; i < Division_no; i++)
		{

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			p = 0;
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z0[p] = pow(10, 0.5*((j - sum0)*1.0) / (size1[i] - 1))*Tasy(Tosz(z_tans[j]), 0.2, j - sum0, size1[i]);
				p++;
			}
			sum += w[i] * Schwefel_function(z0, size1[i]);
			sum0 += size1[i];
		}
	}

	//
	if (problem == 11)
	{
		sum = 0;
		for (int j = 0; j < size1[0] - 1; j++)
		{
			double x10 = X.var[j], x20 = X.var[j + 1];
			sum += 100 * (x10*x10 - x20)*(x10*x10 - x20) + (x10 - 1)*(x10 - 1);
		}
	}

	if (problem == 12)   //��ֵ��̫��̫��ʱ����Ҫ  �����ԵĴ�������������
	{
		sum = 0, sum0 = 0;
		for (i = 0; i < Division_no; i++)
		{
			for (int s = Y_left_number[i]; s < Y_right_number[i]; s++)
			{
				Sum = 0;
				for (int t = Y_left_number[i]; t < Y_right_number[i]; t++)
					Sum += R[i].var[s - Y_left_number[i]][t - Y_left_number[i]] * X.var[t];
				z_tans[s] = Sum;
			}

			p = 0;
			for (j = Y_left_number[i]; j < Y_right_number[i]; j++)
			{
				z0[p] = Tasy(Tosz(z_tans[j]), 0.2, j - Y_left_number[i], Y_right_number[i] - Y_left_number[i]);
				p++;
			}
			sum += w[i] * Schwefel_function(z0, Y_right_number[i] - Y_left_number[i]);
		}
	}

	if (problem == 13)
	{
		sum = 0;
		for (i = 0; i < Division_no; i++)
		{
			for (int s = Y_left_number[i]; s < Y_right_number[i]; s++)
			{
				Sum = 0;
				for (int t = Y_left_number[i]; t < Y_right_number[i]; t++)
					Sum += R[i].var[s - Y_left_number[i]][t - Y_left_number[i]] * X.var[t];
				z_tans[s] = Sum;
			}

			p = 0;
			double Z = Y_right_number[i] - Y_left_number[i] - 1;
			for (j = Y_left_number[i]; j < Y_right_number[i]; j++)
			{
				double X0 = pow(10, 0.5*(j - Y_left_number[i]) / Z);
				z0[p] = X0 * Tasy(Tosz(z_tans[j]), 0.2, j - Y_left_number[i], Y_right_number[i] - Y_left_number[i]);
				p++;
			}
			sum += w[i] * Schwefel_function(z0, Y_right_number[i] - Y_left_number[i]);
		}
	}
	if (problem == 14)
	{
		for (i = 0; i < dim; i++)
		{
			z0[i] = Tasy(Tosz(X.var[i]), 0.2, i, dim);
		}

		sum = Schwefel_function(z0, dim);
	}

	/**************************************************************************************************/
	/*���ԵĲ�λ������*/


	if (problem == 15)///F1
	{
		for (j = 0; j < dim; j++)
			sum += pow(10, 6 * j / (dim - 1.0))*pow(X.var[j], 2);
	}
	//
	if (problem == 16)//F2
	{
		for (j = 0; j < dim; j++)
		{
			sum += pow(X.var[j], 2) - 10 * cos(2 * PI*X.var[j]) + 10;
		}
	}

	if (problem == 17)//F3
	{
		for (j = 0; j < dim; j++)
		{
			z_tans[j] = X.var[j];
		}
		sum = Ackley_function(z_tans, dim);                                  //Ackley_function �����ں��涨��  3/5
	}

	if (problem == 18)//F4
	{
		sum = 0, sum0 = 0;
		for (i = 0; i < Division_no - 1; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
		}
		double x = Elliptic_function(z_tans, size1[0]);

		for (j = size1[0]; j < dim; j++)
		{
			z_tans[j] = X.var[j];
		}
		double y = Elliptic_function(z_tans, size1[Division_no - 1]);
		sum = x * pow(10, 6) + y;
	}
	//
	if (problem == 19)//F5
	{
		sum = 0, sum0 = 0;
		for (i = 0; i < Division_no - 1; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
		}
		double x = Rastrigin_function(z_tans, size1[0]);

		for (j = size1[0]; j < dim; j++)
		{
			z_tans[j] = X.var[j];
		}
		double y = Rastrigin_function(z_tans, size1[Division_no - 1]);
		sum = x * pow(10, 6) + y;
	}
	//
	if (problem == 20)//F6
	{
		sum = 0, sum0 = 0;
		for (i = 0; i < 1; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
		}
		double x = Ackley_function(z_tans, size1[0]);

		for (j = size1[0]; j < dim; j++)
		{
			z_tans[j] = X.var[j];
		}
		double y = Ackley_function(z_tans, size1[Division_no - 1]);
		sum = x * pow(10, 6) + y;
	}

	//
	if (problem == 21)//F7
	{
		sum = 0, sum0 = 0;
		for (i = 0; i < Division_no - 1; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
		}
		double x = Schwefel_function(z_tans, size1[0]);

		for (j = size1[0]; j < dim; j++)
		{
			z_tans[j] = X.var[j];
		}
		double y = Sphere_function(z_tans, size1[Division_no - 1]);
		sum = x * pow(10, 6) + y;
	}

	//
	if (problem == 22)//F8
	{
		sum = 0, sum0 = 0;
		for (j = 0; j < size1[0]; j++)
		{
			z_tans[j] = X.var[j];
		}
		double x = Rosenbrock_function(z_tans, size1[0]);
		for (j = size1[0]; j < dim; j++)
		{
			z_tans[j] = X.var[j];
		}
		double y = Sphere_function(z_tans, size1[Division_no - 1]);
		sum = x * pow(10, 6) + y;
	}

	//
	if (problem == 23)//F9
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no - 1; i++)
		{

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			sum0 += size1[i];
			sum += Elliptic_function(z_tans, size1[i]);
		}

		for (j = 0; j < dim / 2; j++)
		{
			z_tans[j] = X.var[j + dim / 2];
		}
		sum += Sphere_function(z_tans, dim / 2);
	}
	//
	if (problem == 24)//F10
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no - 1; i++)
		{

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			sum += Rastrigin_function(z_tans, size1[i]);
			sum0 += size1[i];
		}

		for (j = 0; j < dim / 2; j++)
		{
			z_tans[j] = X.var[j + dim / 2];
		}
		sum += Rastrigin_function(z_tans, dim / 2);
	}

	//
	if (problem == 25)//F11
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no - 1; i++)
		{

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			sum0 += size1[i];
			sum += Ackley_function(z_tans, size1[i]);
		}

		for (j = 0; j < dim / 2; j++)
		{
			z_tans[j] = X.var[j + dim / 2];
		}
		sum += Ackley_function(z_tans, dim / 2);
	}

	//
	if (problem == 26)//F12
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no - 1; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z_tans[j] = X.var[j];
			}
			sum += Schwefel_function(z_tans, size1[i]);
			sum0 += size1[i];
		}

		for (j = 0; j < dim / 2; j++)
		{
			z_tans[j] = X.var[j + dim / 2];
		}
		sum += Sphere_function(z_tans, dim / 2);
	}

	//
	if (problem == 27)//F13
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no - 1; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z_tans[j] = X.var[j];
			}
			sum += Rosenbrock_function(z_tans, size1[i]);
			sum0 += size1[i];
		}

		for (j = 0; j < dim / 2; j++)
		{
			z_tans[j] = X.var[j + dim / 2];
		}
		sum += Sphere_function(z_tans, dim / 2);
	}

	//
	if (problem == 28)//F14
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			sum += Elliptic_function(z_tans, size1[i]);
			sum0 += size1[i];
		}
	}

	//
	if (problem == 29)//F15
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no; i++)
		{

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			sum += Rastrigin_function(z_tans, size1[i]);
			sum0 += size1[i];
			// cout<<sum<<",";
		}
	}

	//
	if (problem == 30)//F16
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no; i++)
		{

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				Sum = 0;
				for (int k = sum0; k < sum0 + size1[i]; k++)
					Sum += R[i].var[j - sum0][k - sum0] * X.var[k];
				z_tans[j] = Sum;
			}
			sum += Ackley_function(z_tans, size1[i]);
			sum0 += size1[i];
		}
	}

	//
	if (problem == 31)//F17
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no; i++)
		{

			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z_tans[j] = X.var[j];
			}
			sum += Schwefel_function(z_tans, size1[i]);
			sum0 += size1[i];
		}
	}

	//
	if (problem == 32)//F18
	{
		sum0 = 0, sum = 0;
		for (i = 0; i < Division_no; i++)
		{
			for (j = sum0; j < sum0 + size1[i]; j++)
			{
				z_tans[j] = X.var[j];
			}
			sum += Rosenbrock_function(z_tans, size1[i]);
			sum0 += size1[i];
		}
	}

	if (problem == 33)//F19
	{
		for (i = 0; i < dim; i++)
		{
			z_tans[i] = X.var[i];
		}
		sum = Schwefel_function(z_tans, dim);
	}

	if (problem == 34)//F20
	{
		for (i = 0; i < dim; i++)
		{
			z_tans[i] = X.var[i];
		}
		sum = Rosenbrock_function(z_tans, dim);
	}

	return sum;
}



/*����Ϊ�����ÿ�����㺯�������*/
double Rosenbrock_function(double x[], int Dim)
{
	double sum = 0;
	for (int j = 0; j < Dim - 1; j++)
		sum += 100 * (x[j] * x[j] - x[j + 1])*(x[j] * x[j] - x[j + 1]) + (x[j] - 1)*(x[j] - 1);
	return sum;
}

double Elliptic_function(double x[], int Dim)
{
	int i = 0;
	double sum = 0;
	for (i = 0; i < Dim; i++)
	{
		double Q1 = pow(10, (6 * i / (Dim - 1.0)));
		double Q2 = pow(x[i], 2);
		sum += Q1 * Q2;
	}
	return sum;
}
double Sphere_function(double x[], int size)
{
	double F = 0;
	for (int i = 0; i < size; i++)
		F += pow(x[i], 2);
	return F;
}
double Ackley_function(double x[], int Dim)
{
	int i = 0;
	double sum = 0, sum1 = 0, sum2 = 0;
	for (int j = 0; j < Dim; j++)
	{
		sum1 += pow(x[j], 2);
		sum2 += cos(2 * PI*x[j]);
	}
	sum1 = exp(-0.2*sqrt(sum1 / Dim));
	sum2 = exp(sum2 / Dim);
	sum = -20 * sum1 - sum2 + 20 + 2.718281828459045;
	return sum;
}

double Rastrigin_function(double x[], int Dim)
{
	int i = 0;
	double sum = 0;
	for (i = 0; i < Dim; i++)
		sum += pow(x[i], 2) - 10 * cos(2 * PI*x[i]) + 10;
	return sum;
}

double Schwefel_function(double x[], int Dim)
{
	int i = 0, j = 0;
	double sum = 0, sum1 = 0, sum2 = 0;
	for (j = 0; j < Dim; j++)
	{
		sum1 = 0;
		for (int r = 0; r <= j; r++)
			sum1 += x[r];
		sum += pow(sum1, 2);
	}
	return sum;
}

double sign(double F)
{
	double x = 0;
	if (F > 0)
		x = 1;
	if (F == 0)
		x = 0;
	if (F < 0)
		x = -1;
	return x;
}

double Tosz(double x)
{
	double c1 = 0, c2 = 0, W = 0;
	double F = 0, Y = 0;
	if (x > 0)
	{
		c1 = 10, c2 = 7.9;
	}
	else
	{
		c1 = 5.5, c2 = 3.1;
	}

	if (!(x == 0))
	{
		F = log(fabs(x));
	}
	else
	{
		F = 0;
	};
	W = sign(x)*exp(F + 0.049*(sin(c1*F) + sin(c2*F)));
	return W;
}

double Tasy(double x, double beta, int var_number, int Module_size)
{
	double c1 = 0, c2 = 0, W = 0;
	double F = 0, Y = 0;
	if (x > 0)
	{
		F = pow(x, 1 + beta * var_number*sqrt(x) / (Module_size - 1));
	}
	else
	{
		F = x;
	}
	return F;
}

void Problem_Parameter(int problem)                    //��������
{
	int k = 0;
	int left, right;
	for (int j = 0; j < 20; j++)
	{
		double sum = 0;
		for (int i = 0; i < 12; i++)
		{
			sum += (rand() % 1001) / 1000.0;
		}
		w[j] = pow(10, 3 * (sum / 12.0 - 0.5));
	}

	//ÿһ����������½�  ����ά�ȵ���������
	if (problem == 0)
	{
		dim = 1000;
		Division_no = 1;
		size1[0] = dim;
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
	if (problem == 1)
	{
		dim = 1000;
		Division_no = 1;
		size1[0] = dim;
		for (k = 0; k < dim; k++)
		{
			a[k] = -5, b[k] = 5;
		}
	}
	if (problem == 2)
	{
		dim = 1000;
		Division_no = 1;
		size1[0] = dim;
		for (k = 0; k < dim; k++)
		{
			a[k] = -32, b[k] = 32;
		}
	}
	if (problem == 3)
	{
		dim = 1000;
		Division_no = 8;
		size1[0] = 50, size1[1] = 25, size1[2] = 25, size1[3] = 100, size1[4] = 50, size1[5] = 25, size1[6] = 25, size1[7] = 700;
		int sum0 = 0;
		for (int i = 0; i < 7; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			Create_Matrix(i, right - left);
			sum0 += size1[i];
		}
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
	if (problem == 4)
	{
		dim = 1000;
		Division_no = 8;
		size1[0] = 50, size1[1] = 25, size1[2] = 25, size1[3] = 100, size1[4] = 50, size1[5] = 25, size1[6] = 25, size1[7] = 700;
		int sum0 = 0;
		for (int i = 0; i < 7; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -5, b[k] = 5;
		}
	}
	if (problem == 5)
	{
		dim = 1000;
		Division_no = 8;
		size1[0] = 50, size1[1] = 25, size1[2] = 25, size1[3] = 100, size1[4] = 50, size1[5] = 25, size1[6] = 25, size1[7] = 700;
		int sum0 = 0;
		for (int i = 0; i < 7; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			Create_Matrix(i, right - left);
			sum0 += size1[i];
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -32, b[k] = 32;
		}
	}
	if (problem == 6)
	{
		dim = 1000;
		Division_no = 8;
		size1[0] = 50, size1[1] = 25, size1[2] = 25, size1[3] = 100, size1[4] = 50, size1[5] = 25, size1[6] = 25, size1[7] = 700;
		int sum0 = 0;
		for (int i = 0; i < 7; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
	if (problem == 7)
	{
		dim = 1000;
		Division_no = 20;
		size1[0] = 50, size1[1] = 50, size1[2] = 25, size1[3] = 25, size1[4] = 100, size1[5] = 100, size1[6] = 25;
		size1[7] = 25, size1[8] = 50, size1[9] = 25, size1[10] = 100, size1[11] = 25, size1[12] = 100, size1[13] = 50;
		size1[14] = 25, size1[15] = 25, size1[16] = 25, size1[17] = 100, size1[18] = 50, size1[19] = 25;
		int sum0 = 0;
		for (int i = 0; i < 20; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
	if (problem == 8)
	{
		dim = 1000;
		Division_no = 20;
		size1[0] = 50, size1[1] = 50, size1[2] = 25, size1[3] = 25, size1[4] = 100, size1[5] = 100, size1[6] = 25;
		size1[7] = 25, size1[8] = 50, size1[9] = 25, size1[10] = 100, size1[11] = 25, size1[12] = 100, size1[13] = 50;
		size1[14] = 25, size1[15] = 25, size1[16] = 25, size1[17] = 100, size1[18] = 50, size1[19] = 25;
		int sum0 = 0;
		for (int i = 0; i < 20; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -5, b[k] = 5;
		}
	}
	if (problem == 9)
	{
		dim = 1000;
		Division_no = 20;
		size1[0] = 50, size1[1] = 50, size1[2] = 25, size1[3] = 25, size1[4] = 100, size1[5] = 100, size1[6] = 25;
		size1[7] = 25, size1[8] = 50, size1[9] = 25, size1[10] = 100, size1[11] = 25, size1[12] = 100, size1[13] = 50;
		size1[14] = 25, size1[15] = 25, size1[16] = 25, size1[17] = 100, size1[18] = 50, size1[19] = 25;
		int sum0 = 0;
		for (int i = 0; i < 20; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		for (k = 0; k < dim; k++)
		{
			a[k] = -32, b[k] = 32;
		}
	}
	if (problem == 10)
	{
		dim = 1000;
		Division_no = 20;
		size1[0] = 50, size1[1] = 50, size1[2] = 25, size1[3] = 25, size1[4] = 100, size1[5] = 100, size1[6] = 25;
		size1[7] = 25, size1[8] = 50, size1[9] = 25, size1[10] = 100, size1[11] = 25, size1[12] = 100, size1[13] = 50;
		size1[14] = 25, size1[15] = 25, size1[16] = 25, size1[17] = 100, size1[18] = 50, size1[19] = 25;
		int sum0 = 0;
		for (int i = 0; i < 20; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
	if (problem == 11)
	{
		dim = 1000;
		Division_no = 1;
		size1[0] = dim;
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
	if (problem == 12)
	{
		dim = 905;
		Division_no = 20;
		size1[0] = 50, size1[1] = 50, size1[2] = 25, size1[3] = 25, size1[4] = 100, size1[5] = 100, size1[6] = 25;
		size1[7] = 25, size1[8] = 50, size1[9] = 25, size1[10] = 100, size1[11] = 25, size1[12] = 100, size1[13] = 50;
		size1[14] = 25, size1[15] = 25, size1[16] = 25, size1[17] = 100, size1[18] = 50, size1[19] = 25;
		int i = 0, C[20], sum0 = 0;
		for (i = 0; i < 20; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		sum0 = 0;
		for (i = 0; i < Division_no; i++)
		{
			sum0 += size1[i];
			C[i] = sum0;
		}
		Y_left_number[0] = 0;//Y_left_number[0]=1
		Y_right_number[0] = 50;//��һ�������1~50
		for (i = 1; i < Division_no; i++)
		{
			Y_left_number[i] = C[i - 1] - i * 5;//Y_left_number[i]=C[i-1]-i*5;
			Y_right_number[i] = C[i] - i * 5;
		}
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
	if (problem == 13)
	{
		dim = 905;
		Division_no = 20;
		size1[0] = 50, size1[1] = 50, size1[2] = 25, size1[3] = 25, size1[4] = 100, size1[5] = 100, size1[6] = 25;
		size1[7] = 25, size1[8] = 50, size1[9] = 25, size1[10] = 100, size1[11] = 25, size1[12] = 100, size1[13] = 50;
		size1[14] = 25, size1[15] = 25, size1[16] = 25, size1[17] = 100, size1[18] = 50, size1[19] = 25;
		int i = 0, C[20], sum0 = 0;
		for (i = 0; i < 20; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}

		sum0 = 0;
		for (i = 0; i < Division_no; i++)
		{
			sum0 += size1[i];
			C[i] = sum0;
		}
		Y_left_number[0] = 0;//Y_left_number[0]=1
		Y_right_number[0] = 50;//��һ�������1~50
		for (i = 1; i < Division_no; i++)
		{
			Y_left_number[i] = C[i - 1] - i * 5;//Y_left_number[i]=C[i-1]-i*5;
			Y_right_number[i] = C[i] - i * 5;
		}
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 14)
	{
		dim = 1000;
		Division_no = 1;
		size1[0] = dim;
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	/**************************************************************************************/

	if (problem == 15)//F1
	{
		dim = 10000;
		Division_no = 1;
		size1[0] = dim;
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 16)//F2
	{
		dim = 10000;
		Division_no = 1;
		size1[0] = dim;
		for (k = 0; k < dim; k++)
		{
			a[k] = -5, b[k] = 5;
		}
	}

	if (problem == 17)//F3
	{
		dim = 10000;
		Division_no = 1;
		size1[0] = dim;
		for (k = 0; k < dim; k++)
		{
			a[k] = -32, b[k] = 32;
		}
	}


	if (problem == 18)//F4
	{
		dim = 10000;
		Division_no = 2;
		size1[0] = 50, size1[1] = 950;
		int sum0 = 0;
		for (int i = 0; i < Division_no - 1; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 19)//F5
	{
		dim = 10000;
		Division_no = 2;
		size1[0] = 50, size1[1] = 950;
		int sum0 = 0;
		for (int i = 0; i < Division_no - 1; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		for (k = 0; k < dim; k++)
		{
			a[k] = -5, b[k] = 5;
		}
	}

	if (problem == 20)//F6
	{
		dim = 10000;
		Division_no = 2;
		size1[0] = 50, size1[1] = 950;
		int sum0 = 0;
		for (int i = 0; i < Division_no - 1; i++)
		{
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		for (k = 0; k < dim; k++)
		{
			a[k] = -32, b[k] = 32;
		}
	}

	if (problem == 21)//F7
	{
		dim = 10000;
		Division_no = 2;
		size1[0] = 50, size1[1] = 950;
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 22)//F8
	{
		dim = 10000;
		Division_no = 2;
		size1[0] = 50, size1[1] = 950;
		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 23)//F9
	{
		dim = 10000;
		Division_no = 11;
		int sum0 = 0;
		for (int i = 0; i < Division_no - 1; i++)
		{
			size1[i] = 50;
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		size1[Division_no - 1] = 500;

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 24)//F10
	{
		dim = 10000;
		Division_no = 11;
		int sum0 = 0;
		for (int i = 0; i < Division_no - 1; i++)
		{
			size1[i] = 50;
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		size1[Division_no - 1] = 500;

		for (k = 0; k < dim; k++)
		{
			a[k] = -5, b[k] = 5;
		}
	}

	if (problem == 25)//F11
	{
		dim = 10000;
		Division_no = 11;
		int sum0 = 0;
		for (int i = 0; i < Division_no - 1; i++)
		{
			size1[i] = 50;
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		size1[Division_no - 1] = 500;

		for (k = 0; k < dim; k++)
		{
			a[k] = -32, b[k] = 32;
		}
	}

	if (problem == 26)//F12
	{
		dim = 10000;
		Division_no = 11;
		int sum0 = 0;
		for (int i = 0; i < Division_no - 1; i++)
		{
			size1[i] = 50;
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}
		size1[Division_no - 1] = 500;

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 27)//F13
	{
		dim = 10000;
		Division_no = 11;
		int sum0 = 0;
		for (int i = 0; i < Division_no - 1; i++)
		{
			size1[i] = 50;
		}
		size1[Division_no - 1] = 500;

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 28)//F14
	{
		dim = 10000;
		Division_no = 20;
		int sum0 = 0;
		for (int i = 0; i < Division_no; i++)
		{
			size1[i] = 50;
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 29)//F15
	{
		dim = 10000;
		Division_no = 20;
		int sum0 = 0;
		for (int i = 0; i < Division_no; i++)
		{
			size1[i] = 50;
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		} /**/

		for (k = 0; k < dim; k++)
		{
			a[k] = -5, b[k] = 5;
		}
	}

	if (problem == 30)//F16
	{
		dim = 10000;
		Division_no = 20;
		int sum0 = 0;
		for (int i = 0; i < Division_no; i++)
		{
			size1[i] = 50;
			if (i == 0)
			{
				left = 0;
				right = size1[0];
			}
			else
			{
				left = sum0;
				right = sum0 + size1[i];
			}
			sum0 += size1[i];
			Create_Matrix(i, right - left);
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -32, b[k] = 32;
		}
	}

	if (problem == 31)//F17
	{
		dim = 10000;
		Division_no = 20;
		int sum0 = 0;
		for (int i = 0; i < Division_no; i++)
		{
			size1[i] = 50;
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
	if (problem == 32)//F18
	{
		dim = 10000;
		Division_no = 20;
		int sum0 = 0;
		for (int i = 0; i < Division_no; i++)
		{
			size1[i] = 50;
		}

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 33)//F19
	{
		dim = 10000;
		Division_no = 1;
		int sum0 = 0;
		size1[0] = 1000;

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}

	if (problem == 34)//F20
	{
		dim = 10000;
		Division_no = 1;
		int sum0 = 0;
		size1[0] = 1000;

		for (k = 0; k < dim; k++)
		{
			a[k] = -100, b[k] = 100;
		}
	}
}


void Mrot(int i, int j, int Rot_dim)
{
	int k, l;
	for (k = 0; k < Rot_dim; k++)
	{
		for (l = 0; l < Rot_dim; l++)
			if (k == l)
			{
				Mturn.var[k][l] = 1;
			}
			else
			{
				Mturn.var[k][l] = 0;
			}
	}
	double alpha = (rand() % 1001) / 1000.0;
	alpha = (alpha - 0.5)*PI*0.5;
	Mturn.var[i][i] = cos(alpha);
	Mturn.var[j][j] = cos(alpha);
	Mturn.var[i][j] = (-1)*sin(alpha);
	Mturn.var[j][i] = sin(alpha);
}

void MulMatrix(indi A, indi B, int Rot_dim)
{
	double sum = 0;
	for (int i = 0; i < Rot_dim; i++)
		for (int j = 0; j < Rot_dim; j++)
		{
			sum = 0;
			for (int k = 0; k < Rot_dim; k++)
				sum += A.var[i][k] * B.var[k][j];
			Result.var[i][j] = sum;
		}
}

void Create_Matrix(int Module_no, int Rot_dim)
{
	int i, j;
	Unit_Matrix(Rot_dim);
	for (i = 1; i < Rot_dim; i++)
	{
		Mrot(0, i, Rot_dim);
		MulMatrix(Result, Mturn, Rot_dim);
	}
	for (i = 1; i < Rot_dim - 1; i++)
	{
		Mrot(i, Rot_dim - 1, Rot_dim);
		MulMatrix(Result, Mturn, Rot_dim);
	}

	for (i = 0; i < Rot_dim; i++)
		for (j = 0; j < Rot_dim; j++)
			R[Module_no].var[i][j] = Result.var[i][j];
}

void Unit_Matrix(int Rot_dim)
{
	int i, j;
	for (i = 0; i < Rot_dim; i++)
	{
		for (j = 0; j < Rot_dim; j++)
			if (i == j)
			{
				Result.var[i][j] = 1;
			}
			else
			{
				Result.var[i][j] = 0;
			}
	}
}


void Output_Result()      //�������׶�
{
	int i = 0, run = 0;
	fstream Zhang_Result_Table, Zhang_Result_Table0;
	Zhang_Result_Table.open("F:\\ʵ����\\4.txt", ios::out);
	Zhang_Result_Table0.open("F:\\ʵ����\\0.txt", ios::out);
	for (problem = 0; problem < Problem; problem++)
	{
		double min_iter = pow(10, 9);
		for (run = 0; run < TotalRun; run++)
		{
			if (min_iter > Max_iter[problem][run])
				min_iter = Max_iter[problem][run];
		}

		Problem_Min_Iter[problem] = min_iter;
		cout << "min_iter=" << min_iter << endl;
		for (i = 0; i < Problem_Min_Iter[problem]; i++)
			Iter_value[problem][i] /= TotalRun;
	}
	Zhang_Result_Table << "All the results acquired with total of 25" << endl << endl;

	for (problem = 0; problem < Problem; problem++)
	{
		Zhang_Result_Table << "Problem=" << problem << endl;
		Zhang_Result_Table << "Minimization[" << problem << "]=[";
		for (run = 0; run < TotalRun; run++)
			Zhang_Result_Table << Problem_Run_Min_Value[problem][run] << ",";
		Zhang_Result_Table << "]" << endl;

		Zhang_Result_Table0 << "Average_curve[" << problem << "]=[";
		for (i = 0; i < Problem_Min_Iter[problem]; i++)
			Zhang_Result_Table0 << Iter_value[problem][i] << ",";
		Zhang_Result_Table0 << "];" << endl << endl;
	}
}

//void Output_Result()
//{
//	int i = 0, run = 0;
//	fstream LIU_Result_Table, LIU_Result_Table0;
//	LIU_Result_Table.open("F:\\ʵ����\\4.txt", ios::out);
//	LIU_Result_Table0.open("F:\\ʵ����\\0.txt", ios::out);
//	for (problem = 0; problem < Problem; problem++)
//	{
//		double min_iter = pow(10, 9);
//		for (run = 0; run < TotalRun; run++)
//		{
//			if (min_iter > Max_iter[problem][run])
//				min_iter = Max_iter[problem][run];
//		}
//
//		Problem_Min_Iter[problem] = min_iter;
//		//cout << "min_iter=" << min_iter << endl;
//		for (i = 0; i < Problem_Min_Iter[problem]; i++)
//			Iter_value[problem][i] /= TotalRun;
//	}
//	LIU_Result_Table << "All the results acquired with total of 25" << endl << endl;
//
//	for (problem = 0; problem < Problem; problem++)
//	{
//		LIU_Result_Table << "Problem=" << problem << endl;
//		LIU_Result_Table << "Minimization[" << problem << "]=[";
//		for (run = 0; run < TotalRun; run++)
//			LIU_Result_Table << fixed << setprecision(5) << Problem_Run_Min_Value[problem][run] << ",";
//		LIU_Result_Table << "]" << endl;
//
//		LIU_Result_Table0 << "Problem_Run_constraint_Value[" << problem << "]=[";
//		for (run = 0; run < TotalRun; run++)
//			LIU_Result_Table0 << fixed << setprecision(3) << Problem_Run_constraint_Value[problem][run] << ",";
//		LIU_Result_Table0 << "];" << endl << endl;
//	}
//}