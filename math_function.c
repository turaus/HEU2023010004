#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constant.h"
extern double froller_radius;       //前滚轮半径
extern double fheight;              //前滚轮高度
extern double stroke;               //平均冲程
extern double aver_stressradius;    //平均压力半径
extern double diameter;             //气缸直径
extern double piston_mess;          //活塞质量
extern double p_q;                 //凸轮箱中压力
extern double u1, u2;          //凸轮与圆柱材料泊松比u1，u2
extern double E1, E2;          //凸轮与滚轮材料弹性模量E1，E2 

//角线的曲线方程
//传入参数：转过角度，x，y，z的地址
void fwork_angle_function_X(double stressradius,double angle,double *X) {
	
	*X = -stroke * pow(sin(angle), 2) - froller_radius * stressradius / sqrt(pow(stressradius, 2) + pow(stroke * sin(2 * angle), 2));
	
}
void fwork_angle_function_Y(double stressradius,double angle,double *Y) {

	*Y = stressradius * cos(angle) + froller_radius * sin(angle) * stroke * sin(2 * angle) / sqrt(pow(stressradius, 2) + pow(stroke * sin(2 * angle), 2));

}
void fwork_angle_function_Z(double stressradius,double angle, double* Z) {

	*Z = stressradius * sin(angle) - froller_radius * cos(angle) * stroke * sin(2 * angle) / sqrt(pow(stressradius, 2) + pow(stroke * sin(2 * angle), 2));

}

//(法向）接触应力模型
//传入参数：缸内工作气体压强，转过角度，内外轴相对转速
double n_force(double gas_pressure, double angle, double n0) {
	double z_f = (PI / 4) * diameter * (gas_pressure - p_q) - 2 * piston_mess * stroke * pow(2 * PI * n0 / 60, 2) * cos(2 * angle);
	double n_f = z_f * sqrt(pow(aver_stressradius, 2) + pow(stroke * sin(2 * angle), 2)) / aver_stressradius;
	return n_f;
}

//向量求导数（参数方程）
//根据定义所求出的近似值
double move = 1e-2;
//下列三个为对角度（angle）求导
double x_der(double stressradius,double angle) {
	double der_x;     //x坐标导数值
	//坐标起始、终止点
	double* x1 = (double*)malloc(sizeof(double));
	double* x2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //变化量
	fwork_angle_function_X(stressradius,angle,x1);
	fwork_angle_function_X(stressradius, angle + move, x2);
	der_x = (*x2 - *x1) / move;
	return der_x;
}

double y_der(double stressradius,double angle) {
	double der_y;     //y坐标导数值
	//坐标起始、终止点
	double* y1 = (double*)malloc(sizeof(double));
	double* y2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //变化量
	fwork_angle_function_X(stressradius, angle, y1);
	fwork_angle_function_X(stressradius, angle + move, y2);
	der_y = (*y2 - *y1) / move;
	return der_y;
}

double z_der(double stressradius,double angle) {
	double der_z;     //z坐标导数值
	//坐标起始、终止点
	double* z1 = (double*)malloc(sizeof(double));
	double* z2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //变化量
	fwork_angle_function_X(stressradius, angle, z1);
	fwork_angle_function_X(stressradius, angle + move, z2);
	der_z = (*z2 - *z1) / move;
	return der_z;
}

//对压力半径求导
double x_der_r(double stressradius, double angle) {
	double der_x;     //x坐标导数值
	//坐标起始、终止点
	double* x1 = (double*)malloc(sizeof(double));
	double* x2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //变化量
	fwork_angle_function_X(stressradius, angle, x1);
	fwork_angle_function_X(stressradius + move, angle , x2);
	der_x = (*x2 - *x1) / move;
	return der_x;
}

double y_der_r(double stressradius, double angle) {
	double der_y;     //y坐标导数值
	//坐标起始、终止点
	double* y1 = (double*)malloc(sizeof(double));
	double* y2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //变化量
	fwork_angle_function_X(stressradius, angle, y1);
	fwork_angle_function_X(stressradius + move, angle, y2);
	der_y = (*y2 - *y1) / move;
	return der_y;
}

double z_der_r(double stressradius, double angle) {
	double der_z;     //z坐标导数值
	//坐标起始、终止点
	double* z1 = (double*)malloc(sizeof(double));
	double* z2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //变化量
	fwork_angle_function_X(stressradius, angle, z1);
	fwork_angle_function_X(stressradius + move, angle, z2);
	der_z = (*z2 - *z1) / move;
	return der_z;
}

//求解工作曲面接触点法曲率
//向量的模长
double v_length(struct Vector R) {
	return sqrt(pow(R.position[0], 2) + pow(R.position[1], 2) + pow(R.position[2], 2));
}
//向量的外积,得到向量
struct Vector v_OuterProduct(struct Vector R1,struct Vector R2) {
	struct Vector R3;
	R3.position[0] = R1.position[1] * R2.position[2] - R2.position[1] * R1.position[2];
	R3.position[1] = R1.position[2] * R2.position[0] - R2.position[2] * R1.position[0];
	R3.position[2] = R1.position[0] * R2.position[1] - R2.position[0] * R1.position[1];
	return R3;
}

double v_InnerProduct(struct Vector R1, struct Vector R2) {
	return R1.position[0] * R2.position[0] + R1.position[1] * R2.position[1] + R1.position[2] * R2.position[2];
}
//角线在P点的诱导法曲率
double ng_K(double stressradius,double angle) {
	//K部分
	struct Vector single_der, twice_der;
	//向量的一阶导数
	single_der.position[0] = x_der(aver_stressradius,angle);
	single_der.position[1] = y_der(aver_stressradius, angle);
	single_der.position[2] = z_der(aver_stressradius, angle);
	//向量的二阶导数
	twice_der.position[0] = (x_der(aver_stressradius, angle + move) - x_der(aver_stressradius, angle)) / move;
	twice_der.position[1] = (y_der(aver_stressradius, angle + move) - y_der(aver_stressradius, angle)) / move;
	twice_der.position[2] = (z_der(aver_stressradius, angle + move) - z_der(aver_stressradius, angle)) / move;
	double devident = v_length(v_OuterProduct(single_der, twice_der));
	double devidor = pow(v_length(single_der), 3);
	//β矢量的系数
	double temp = devident / devidor;
	
	//β矢量
	struct Vector B;
	double a = v_length(single_der) / v_length(v_OuterProduct(single_der, twice_der));
	double b = v_InnerProduct(single_der, twice_der) / (v_length(single_der) * v_length(v_OuterProduct(single_der, twice_der)));
	for (int i = 0; i < 3; i++) {
		B.position[i] = a * twice_der.position[i] - b * single_der.position[i];
	}

	struct Vector K;
	for (int j = 0; j < 3; j++) {
		K.position[j] = temp * B.position[j];
	}
	//K部分结束
	//法矢量n部分
	struct Vector n, n0, n1, n2;
	//对压力半径导
	n1.position[0] = x_der_r(stressradius, angle);
	n1.position[1] = y_der_r(stressradius, angle);
	n1.position[2] = z_der_r(stressradius, angle);
	//对角度导
	n2.position[0] = x_der(stressradius, angle);
	n2.position[1] = y_der(stressradius, angle);
	n2.position[2] = z_der(stressradius, angle);

	n0 = v_OuterProduct(n1, n2);
	double t = v_length(n0);
	for (int i = 1; i < 3; i++) {
		n.position[i] = n0.position[i] / t;
	}

	double n_K;
	n_K = v_InnerProduct(K, n);

	//诱导法曲率
	return n_K - 1 / froller_radius;
}

//根据赫兹接触理论可得凸轮前工作曲面所承受的接触应力为 :
//double contact_stress = sqrt((Fj() / fheight * abs(Kng())) / (PI * ( (1 - pow(u1, 2)) / E1 + (1 - pow(u2, 2)) / E2)));
double Contact_Stress(double stressradius, double angle, double gas_pressure, double n0) {
	double a = n_force(gas_pressure, angle, n0);
	double b = (1 - pow(u1, 2)) / E1 + (1 - pow(u2, 2)) / E2;
	double contact_stress = sqrt(a / fheight * abs(ng_K(stressradius, angle)) / (PI * b));
}