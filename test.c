#define _CRT_SECURE_NO_WARNINGS
#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//给定常量
extern double froller_radius = 0.01;       //前滚轮半径
extern double fheight = 0.01;              //前滚轮高度
extern double stroke = 0.02;               //平均冲程
extern double aver_stressradius = 0.045;   //平均压力半径
extern double diameter = 0.02;             //气缸直径
extern double piston_mess = 0.2;           //活塞质量
extern double p_q = 1e+6;                  //凸轮箱中压力
extern double u1 = 0.3, u2 = 0.3;          //凸轮与圆柱材料泊松比u1，u2
extern double E1 = 2e+11, E2 = 2e+11;      //凸轮与滚轮材料弹性模量E1，E2 

int main() {

	//根据赫兹接触理论可得凸轮前工作曲面所承受的接触应力为 :
	//double contact_stress = sqrt((Fj() / fheight * abs(Kng())) / (PI * ( (1 - pow(u1, 2)) / E1 + (1 - pow(u2, 2)) / E2)));
	double angle[100];
	int i, j;
	for (i = 100, j = 0; i > 0; i--, j++) {      //角线的变量
		angle[j] = PI / (float)i;
	}
	//定义角线的坐标，并将相应坐标存在列表中
	double *p_x = (double*)malloc(sizeof(double) * 100);
	double *p_y = (double*)malloc(sizeof(double) * 100);
	double *p_z = (double*)malloc(sizeof(double) * 100);
	double x[100], y[100], z[100];
	p_x = &x, p_y = &y, p_z = &z;
	
	struct Vector R[100];
	
	//角线方程的函数声明
	void fwork_angle_function_X(double stressradius, double angle, double* X);
	void fwork_angle_function_Y(double stressradius, double angle, double* Y);
	void fwork_angle_function_Z(double stressradius, double angle, double* Z);

	for (i = 0; i < 100; i++) {
		fwork_angle_function_X(aver_stressradius, angle[i], p_x+i);
		fwork_angle_function_Y(aver_stressradius, angle[i], p_y+i);
		fwork_angle_function_Z(aver_stressradius, angle[i], p_z+i);
		//printf("%d	%e,%e,%e	,%lf",i, *(p_x+i), *(p_y+i), *(p_z+i),angle[i]);
		//printf("%e,%e,%e\n", x[i], y[i], z[i]);
		//放入向量R中
		R[i].position[0] = x[i];
		R[i].position[1] = y[i];
		R[i].position[2] = z[i];
		
	}

	double x_der(double);
	printf("%lf\n", x_der(aver_stressradius,angle[0]));
	struct Vector v_OuterProduct(struct Vector R1, struct Vector R2);
	struct Vector o = v_OuterProduct(R[0], R[1]);
	/*printf("%e,%e,%e\n", o.position[0], o.position[1], o.position[2]);
	printf("%e,%e,%e\n", R[0].position[0], R[0].position[1], R[0].position[2]);
	printf("%e,%e,%e\n", R[1].position[0], R[1].position[1], R[1].position[2]);*/
	return 0;
}



