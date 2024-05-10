#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constant.h"
extern double froller_radius;       //ǰ���ְ뾶
extern double fheight;              //ǰ���ָ߶�
extern double stroke;               //ƽ�����
extern double aver_stressradius;    //ƽ��ѹ���뾶
extern double diameter;             //����ֱ��
extern double piston_mess;          //��������
extern double p_q;                 //͹������ѹ��
extern double u1, u2;          //͹����Բ�����ϲ��ɱ�u1��u2
extern double E1, E2;          //͹������ֲ��ϵ���ģ��E1��E2 

//���ߵ����߷���
//���������ת���Ƕȣ�x��y��z�ĵ�ַ
void fwork_angle_function_X(double stressradius,double angle,double *X) {
	
	*X = -stroke * pow(sin(angle), 2) - froller_radius * stressradius / sqrt(pow(stressradius, 2) + pow(stroke * sin(2 * angle), 2));
	
}
void fwork_angle_function_Y(double stressradius,double angle,double *Y) {

	*Y = stressradius * cos(angle) + froller_radius * sin(angle) * stroke * sin(2 * angle) / sqrt(pow(stressradius, 2) + pow(stroke * sin(2 * angle), 2));

}
void fwork_angle_function_Z(double stressradius,double angle, double* Z) {

	*Z = stressradius * sin(angle) - froller_radius * cos(angle) * stroke * sin(2 * angle) / sqrt(pow(stressradius, 2) + pow(stroke * sin(2 * angle), 2));

}

//(���򣩽Ӵ�Ӧ��ģ��
//������������ڹ�������ѹǿ��ת���Ƕȣ����������ת��
double n_force(double gas_pressure, double angle, double n0) {
	double z_f = (PI / 4) * diameter * (gas_pressure - p_q) - 2 * piston_mess * stroke * pow(2 * PI * n0 / 60, 2) * cos(2 * angle);
	double n_f = z_f * sqrt(pow(aver_stressradius, 2) + pow(stroke * sin(2 * angle), 2)) / aver_stressradius;
	return n_f;
}

//�����������������̣�
//���ݶ���������Ľ���ֵ
double move = 1e-2;
//��������Ϊ�ԽǶȣ�angle����
double x_der(double stressradius,double angle) {
	double der_x;     //x���굼��ֵ
	//������ʼ����ֹ��
	double* x1 = (double*)malloc(sizeof(double));
	double* x2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //�仯��
	fwork_angle_function_X(stressradius,angle,x1);
	fwork_angle_function_X(stressradius, angle + move, x2);
	der_x = (*x2 - *x1) / move;
	return der_x;
}

double y_der(double stressradius,double angle) {
	double der_y;     //y���굼��ֵ
	//������ʼ����ֹ��
	double* y1 = (double*)malloc(sizeof(double));
	double* y2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //�仯��
	fwork_angle_function_X(stressradius, angle, y1);
	fwork_angle_function_X(stressradius, angle + move, y2);
	der_y = (*y2 - *y1) / move;
	return der_y;
}

double z_der(double stressradius,double angle) {
	double der_z;     //z���굼��ֵ
	//������ʼ����ֹ��
	double* z1 = (double*)malloc(sizeof(double));
	double* z2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //�仯��
	fwork_angle_function_X(stressradius, angle, z1);
	fwork_angle_function_X(stressradius, angle + move, z2);
	der_z = (*z2 - *z1) / move;
	return der_z;
}

//��ѹ���뾶��
double x_der_r(double stressradius, double angle) {
	double der_x;     //x���굼��ֵ
	//������ʼ����ֹ��
	double* x1 = (double*)malloc(sizeof(double));
	double* x2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //�仯��
	fwork_angle_function_X(stressradius, angle, x1);
	fwork_angle_function_X(stressradius + move, angle , x2);
	der_x = (*x2 - *x1) / move;
	return der_x;
}

double y_der_r(double stressradius, double angle) {
	double der_y;     //y���굼��ֵ
	//������ʼ����ֹ��
	double* y1 = (double*)malloc(sizeof(double));
	double* y2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //�仯��
	fwork_angle_function_X(stressradius, angle, y1);
	fwork_angle_function_X(stressradius + move, angle, y2);
	der_y = (*y2 - *y1) / move;
	return der_y;
}

double z_der_r(double stressradius, double angle) {
	double der_z;     //z���굼��ֵ
	//������ʼ����ֹ��
	double* z1 = (double*)malloc(sizeof(double));
	double* z2 = (double*)malloc(sizeof(double));
	double move = 1e-2;    //�仯��
	fwork_angle_function_X(stressradius, angle, z1);
	fwork_angle_function_X(stressradius + move, angle, z2);
	der_z = (*z2 - *z1) / move;
	return der_z;
}

//��⹤������Ӵ��㷨����
//������ģ��
double v_length(struct Vector R) {
	return sqrt(pow(R.position[0], 2) + pow(R.position[1], 2) + pow(R.position[2], 2));
}
//���������,�õ�����
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
//������P����յ�������
double ng_K(double stressradius,double angle) {
	//K����
	struct Vector single_der, twice_der;
	//������һ�׵���
	single_der.position[0] = x_der(aver_stressradius,angle);
	single_der.position[1] = y_der(aver_stressradius, angle);
	single_der.position[2] = z_der(aver_stressradius, angle);
	//�����Ķ��׵���
	twice_der.position[0] = (x_der(aver_stressradius, angle + move) - x_der(aver_stressradius, angle)) / move;
	twice_der.position[1] = (y_der(aver_stressradius, angle + move) - y_der(aver_stressradius, angle)) / move;
	twice_der.position[2] = (z_der(aver_stressradius, angle + move) - z_der(aver_stressradius, angle)) / move;
	double devident = v_length(v_OuterProduct(single_der, twice_der));
	double devidor = pow(v_length(single_der), 3);
	//��ʸ����ϵ��
	double temp = devident / devidor;
	
	//��ʸ��
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
	//K���ֽ���
	//��ʸ��n����
	struct Vector n, n0, n1, n2;
	//��ѹ���뾶��
	n1.position[0] = x_der_r(stressradius, angle);
	n1.position[1] = y_der_r(stressradius, angle);
	n1.position[2] = z_der_r(stressradius, angle);
	//�ԽǶȵ�
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

	//�յ�������
	return n_K - 1 / froller_radius;
}

//���ݺ��ȽӴ����ۿɵ�͹��ǰ�������������ܵĽӴ�Ӧ��Ϊ :
//double contact_stress = sqrt((Fj() / fheight * abs(Kng())) / (PI * ( (1 - pow(u1, 2)) / E1 + (1 - pow(u2, 2)) / E2)));
double Contact_Stress(double stressradius, double angle, double gas_pressure, double n0) {
	double a = n_force(gas_pressure, angle, n0);
	double b = (1 - pow(u1, 2)) / E1 + (1 - pow(u2, 2)) / E2;
	double contact_stress = sqrt(a / fheight * abs(ng_K(stressradius, angle)) / (PI * b));
}