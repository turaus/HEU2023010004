#define _CRT_SECURE_NO_WARNINGS
#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//��������
extern double froller_radius = 0.01;       //ǰ���ְ뾶
extern double fheight = 0.01;              //ǰ���ָ߶�
extern double stroke = 0.02;               //ƽ�����
extern double aver_stressradius = 0.045;   //ƽ��ѹ���뾶
extern double diameter = 0.02;             //����ֱ��
extern double piston_mess = 0.2;           //��������
extern double p_q = 1e+6;                  //͹������ѹ��
extern double u1 = 0.3, u2 = 0.3;          //͹����Բ�����ϲ��ɱ�u1��u2
extern double E1 = 2e+11, E2 = 2e+11;      //͹������ֲ��ϵ���ģ��E1��E2 

int main() {

	//���ݺ��ȽӴ����ۿɵ�͹��ǰ�������������ܵĽӴ�Ӧ��Ϊ :
	//double contact_stress = sqrt((Fj() / fheight * abs(Kng())) / (PI * ( (1 - pow(u1, 2)) / E1 + (1 - pow(u2, 2)) / E2)));
	double angle[100];
	int i, j;
	for (i = 100, j = 0; i > 0; i--, j++) {      //���ߵı���
		angle[j] = PI / (float)i;
	}
	//������ߵ����꣬������Ӧ��������б���
	double *p_x = (double*)malloc(sizeof(double) * 100);
	double *p_y = (double*)malloc(sizeof(double) * 100);
	double *p_z = (double*)malloc(sizeof(double) * 100);
	double x[100], y[100], z[100];
	p_x = &x, p_y = &y, p_z = &z;
	
	struct Vector R[100];
	
	//���߷��̵ĺ�������
	void fwork_angle_function_X(double stressradius, double angle, double* X);
	void fwork_angle_function_Y(double stressradius, double angle, double* Y);
	void fwork_angle_function_Z(double stressradius, double angle, double* Z);

	for (i = 0; i < 100; i++) {
		fwork_angle_function_X(aver_stressradius, angle[i], p_x+i);
		fwork_angle_function_Y(aver_stressradius, angle[i], p_y+i);
		fwork_angle_function_Z(aver_stressradius, angle[i], p_z+i);
		//printf("%d	%e,%e,%e	,%lf",i, *(p_x+i), *(p_y+i), *(p_z+i),angle[i]);
		//printf("%e,%e,%e\n", x[i], y[i], z[i]);
		//��������R��
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



