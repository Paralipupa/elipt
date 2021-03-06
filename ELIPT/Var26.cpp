#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;

double** point = NULL;
double** gamma_betta = NULL;
int** finit = NULL;
int* finit_in_area = NULL;
int** kraev = NULL;
int* point_in_area = NULL;

int* ig = NULL;
int* jg = NULL;
double* di = NULL;
double* ggl = NULL;
double* ggu = NULL;
double* F = NULL;
double* x = NULL;

double* z, * r, * p, * t, * r1, * l, * l1, * f;

int num_points = 0;//количество точек
int num_finit_elements = 0;//количество конечных элементов
int num_areas = 0;
int num_kraev = 0;

double func(double x, double y, int i)
{
	if (i == 0)
		return -20;
	return 0;
}

double func_kraev1(double* x, int k)
{
	switch (k) {
	case 0:  return x[1] * x[1];
	default: return 0;
	}
}

double func_kraev2(double* x, int k)
{
	switch (k)
	{
	case 0: return 20;
	case 1: return 0;
		//case 2: return (2.);
	default: return 0;
	}
}

double func_kraev3(double* x, int k)
{
	if (k == 0)
		return (20 * x[1] - 27);
	return 0;
}

double resh(double x, double y, int k)
{
	switch (k)
	{
	case 0: return (y * y);
	case 1: return 20.0 * y - 19.0;
	default: return 0;
	}
}

int input()
{
	ifstream _file("coords.txt");

	if (!_file.is_open())
	{
		_file.close();
		return 1;
	}

	_file >> num_points;

	if (num_points <= 0)
	{
		_file.close();
		return 2;
	}

	point = new double* [num_points];
	for (int i = 0; i < num_points; i++) {
		point[i] = new double[2];
		for (int j = 0; j < 2; j++)
			_file >> point[i][j];
	}
	_file.close();

	ifstream _file1("finit_elements.txt");
	if (!_file1.is_open()) {
		_file1.close();
		return 1;
	}
	_file1 >> num_finit_elements;
	if (num_finit_elements <= 0) {
		_file.close();
		return 2;
	}
	finit = new int* [num_finit_elements];
	for (int i = 0; i < num_finit_elements; i++) {
		finit[i] = new int[3];
		for (int j = 0; j < 3; j++)
			_file1 >> finit[i][j];
	}
	_file1.close();

	ifstream _file2("areas.txt");
	if (!_file2.is_open()) {
		_file2.close();
		return 1;
	}
	_file2 >> num_areas;
	if (num_areas <= 0) {
		_file.close();
		return 2;
	}
	gamma_betta = new double* [num_areas];
	finit_in_area = new int[num_finit_elements];
	for (int i = 0; i < num_areas; i++) {
		gamma_betta[i] = new double[2];
		for (int j = 0; j < 2; j++)
			_file2 >> gamma_betta[i][j];
	}
	for (int i = 0; i < num_finit_elements; i++)
		_file2 >> finit_in_area[i];
	_file2.close();

	ifstream _file3("board.txt");
	if (!_file3.is_open()) {
		_file3.close();
		return 1;
	}
	_file3 >> num_kraev;
	kraev = new int* [num_kraev];
	for (int i = 0; i < num_kraev; i++) {
		kraev[i] = new int[4];
		for (int j = 0; j < 5; j++)
			_file3 >> kraev[i][j];
	}
	_file3.close();

	ifstream _file4("point_in_area.txt");
	if (!_file4.is_open()) {
		_file4.close();
		return 1;
	}
	point_in_area = new int[num_points];
	for (int i = 0; i < num_points; i++)
		_file4 >> point_in_area[i];
	_file4.close();
	return 0;
}
void sort(int* mas, int k)
{
	int l = 0;
	for (int i = k - 1; i >= 0; i--)
		for (int j = 0; j <= i; j++)
			if (mas[j] > mas[j + 1]) {
				l = mas[j];
				mas[j] = mas[j + 1];
				mas[j + 1] = l;
			}
}

void portret()
{
	int i = 0;
	int j = 0;
	int k = 0;
	int kk = 0;
	int key = 0;
	int position = 0;//позиция в массиве jg, в которую надо добавлять
	int* mas = new int[num_points];
	struct List {
		int num;
		List* next;
	};

	List* list = new List[num_points];
	List* p = NULL;

	ig = new int[num_points + 1];

	for (i = 0; i < num_points; i++)
		list[i].next = NULL;

	//составление "массива", содержащего номер точки и смежные с ней 
	for (i = 0; i < num_finit_elements; i++) {
		for (j = 0; j < 3; j++) {
			key = 0;
			k = finit[i][(j == 2)];//0 0 1
			kk = finit[i][((j != 0) + 1)];// 1 2 2
			if (k < kk) {
				k += kk;
				kk = k - kk;
				k -= kk;
			}
			p = &list[k];
			while (p->next) {
				if (p->next->num == kk) {
					key = 1;
					break;
				}
				p = p->next;
			}
			if (!key) {
				p->next = new List;
				p->next->num = kk;
				p->next->next = NULL;
			}
		}
	}

	//составление массива ig
	ig[0] = 0;
	for (i = 0; i < num_points; i++) {
		k = 0;
		p = &list[i];
		while (p = p->next)
			k++;
		ig[i + 1] = ig[i] + k;
	}

	jg = new int[ig[i] - 1];
	//составление массива jg
	for (i = 0; i < num_points; i++) {
		k = 0;
		key = 0;
		p = &list[i];
		while (p = p->next) {
			mas[k] = p->num;
			k++;
			key = 1;
		}
		if (key) {
			sort(mas, --k);//сортировка
			int ii = 0;//добавляет в jg
			int jj = 0;
			for (ii = position, jj = 0; ii <= k + position; ii++, jj++)
				jg[ii] = mas[jj];

			position += k + 1;
		}
	}
}

void tochnoe()
{
	for (int i = 0; i < num_points; i++)
		x[i] = resh(point[i][0], point[i][1], point_in_area[i]);
}

double lambda(int k)
{
	if (!k) return 10;
	return 1;
}

double mes_G(double* x, double* y)
{
	return (sqrt((y[0] - x[0]) * (y[0] - x[0]) + (y[1] - x[1]) * (y[1] - x[1])));
}

void M_matrix(double* p1, double* p2, double* p3, double gamma, double** M_matr, double* local_F, int num_of_area)
{
	int i = 0;
	int j = 0;
	double det = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
	double mnoz = fabs(det) / 24;
	double* f = new double[3];
	double mnoz2 = mnoz * gamma;

	f[0] = mnoz * func(p1[0], p1[1], num_of_area);
	f[1] = mnoz * func(p2[0], p2[1], num_of_area);
	f[2] = mnoz * func(p3[0], p3[1], num_of_area);

	local_F[0] = 2 * f[0] + f[1] + f[2];
	local_F[1] = f[0] + 2 * f[1] + f[2];
	local_F[2] = f[0] + f[1] + 2 * f[2];

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			if (i == j)
				M_matr[i][j] = 2 * mnoz2;
			else
				M_matr[i][j] = mnoz2;

}

void G_matrix(double* p1, double* p2, double* p3, double** G_matr, int k)
{
	int i = 0;
	int j = 0;
	double ck = 0;
	double det = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
	double* p12 = new double[2];
	double* p23 = new double[2];
	double* p31 = new double[2];
	for (i = 0; i < 2; i++) {
		p12[i] = (p1[i] + p2[i]) / 2;
		p23[i] = (p2[i] + p3[i]) / 2;
		p31[i] = (p3[i] + p1[i]) / 2;
	}

	ck = lambda(k) * fabs(det) / 2;
	G_matr[0][0] = ck * ((p2[1] - p3[1]) * (p2[1] - p3[1]) / (det * det) + (p3[0] - p2[0]) * (p3[0] - p2[0]) / (det * det));
	G_matr[0][1] = ck * ((p2[1] - p3[1]) * (p3[1] - p1[1]) / (det * det) + (p3[0] - p2[0]) * (p1[0] - p3[0]) / (det * det));
	G_matr[0][2] = ck * ((p2[1] - p3[1]) * (p1[1] - p2[1]) / (det * det) + (p3[0] - p2[0]) * (p2[0] - p1[0]) / (det * det));
	G_matr[1][0] = ck * ((p3[1] - p1[1]) * (p2[1] - p3[1]) / (det * det) + (p1[0] - p3[0]) * (p3[0] - p2[0]) / (det * det));
	G_matr[1][1] = ck * ((p3[1] - p1[1]) * (p3[1] - p1[1]) / (det * det) + (p1[0] - p3[0]) * (p1[0] - p3[0]) / (det * det));
	G_matr[1][2] = ck * ((p3[1] - p1[1]) * (p1[1] - p2[1]) / (det * det) + (p1[0] - p3[0]) * (p2[0] - p1[0]) / (det * det));
	G_matr[2][0] = ck * ((p2[1] - p3[1]) * (p1[1] - p2[1]) / (det * det) + (p3[0] - p2[0]) * (p2[0] - p1[0]) / (det * det));
	G_matr[2][1] = ck * ((p3[1] - p1[1]) * (p1[1] - p2[1]) / (det * det) + (p1[0] - p3[0]) * (p2[0] - p1[0]) / (det * det));
	G_matr[2][2] = ck * ((p1[1] - p2[1]) * (p1[1] - p2[1]) / (det * det) + (p2[0] - p1[0]) * (p2[0] - p1[0]) / (det * det));
}

void local_matrix(int num_of_finit_element, double** local_matr, double* local_F)
{
	int k = finit[num_of_finit_element][0];
	int l = finit[num_of_finit_element][1];
	int m = finit[num_of_finit_element][2];
	double** M_matr = new double* [3];
	double** G_matr = new double* [3];

	for (int i = 0; i < 3; i++) {
		M_matr[i] = new double[3];
		G_matr[i] = new double[3];
	}


	M_matrix(point[k], point[l], point[m],
		gamma_betta[finit_in_area[num_of_finit_element]][0],
		M_matr, local_F,
		finit_in_area[num_of_finit_element]);

	G_matrix(point[k], point[l], point[m],
		G_matr,
		finit_in_area[num_of_finit_element]);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			local_matr[i][j] = M_matr[i][j] + G_matr[i][j];
}

int uchet_kraev(int* current_kraev, double** a, double* b, double betta)
{
	double ck = 0;
	if (current_kraev[3] == 1)
		return 0;
	else
		if (current_kraev[3] == 2)
		{
			ck = mes_G(point[current_kraev[1]], point[current_kraev[2]]) / 6.0;
			b[0] = ck * (2 * func_kraev2(point[current_kraev[1]], current_kraev[4]) + func_kraev2(point[current_kraev[2]], current_kraev[4]));
			b[1] = ck * (func_kraev2(point[current_kraev[1]], current_kraev[4]) + 2 * func_kraev2(point[current_kraev[2]], current_kraev[4]));
			return 1;
		}
		else
		{
			ck = betta * mes_G(point[current_kraev[1]], point[current_kraev[2]]) / 6;
			a[0][0] = 2 * ck;
			a[0][1] = ck;
			a[1][0] = ck;
			a[1][1] = 2 * ck;
			b[0] = ck * (2 * func_kraev3(point[current_kraev[1]], current_kraev[4]) + func_kraev3(point[current_kraev[2]], current_kraev[4]));
			b[1] = ck * (func_kraev3(point[current_kraev[1]], current_kraev[4]) + 2 * func_kraev3(point[current_kraev[2]], current_kraev[4]));
			return 2;
		}
}

void pervoe_kraevoe(int* current_kraev)
{
	int kol = 0, m = 0;
	int lbeg;
	int lend;
	di[current_kraev[1]] = 1;
	di[current_kraev[2]] = 1;

	F[current_kraev[1]] = func_kraev1(point[current_kraev[1]], current_kraev[4]);
	F[current_kraev[2]] = func_kraev1(point[current_kraev[2]], current_kraev[4]);

	kol = ig[current_kraev[1] + 1] - ig[current_kraev[1]];

	for (int i = 0; i < kol; i++)
		ggl[ig[current_kraev[1]] + i] = 0;

	kol = ig[current_kraev[2] + 1] - ig[current_kraev[2]];
	for (int i = 0; i < kol; i++)
		ggl[ig[current_kraev[2]] + i] = 0;

	for (int i = current_kraev[1] + 1; i < num_points; i++) {
		lbeg = ig[i];
		lend = ig[i + 1];
		for (int p = lbeg; p < lend; p++)
			if (jg[p] == current_kraev[1]) {
				ggu[p] = 0;
				continue;
			}
	}
	for (int i = current_kraev[2] + 1; i < num_points; i++)
	{
		lbeg = ig[i];
		lend = ig[i + 1];
		for (int p = lbeg; p < lend; p++)
			if (jg[p] == current_kraev[2]) {
				ggu[p] = 0;
				continue;
			}
	}
}

void global_matrix()
{
	int i = 0;
	int j = 0;
	int k = 0;
	int p = 0;
	int ibeg = 0;
	int iend = 0;
	int h = 0;
	int key = 0;
	int kol = 0;
	int* L = new int[3];
	int* L2 = new int[2];
	int* K = new int[num_points / 2];

	double* local_F = new double[3];
	double** local_matr = new double* [3];
	double* b = new double[2];//вектор для краевых
	double** a = new double* [2];//матрица для краевых

	for (i = 0; i < 2; i++) {
		a[i] = new double[2];
		b[i] = 0;
		for (j = 0; j < 2; j++)
			a[i][j] = 0;
	}

	for (i = 0; i < 3; i++)
		local_matr[i] = new double[3];

	for (k = 0; k < num_finit_elements; k++) {
		local_matrix(k, local_matr, local_F);
		memcpy(L, finit[k], 3 * sizeof(double));
		//локальная в глобальную//
		for (i = 0; i < 3; i++) {
			//h=0;
			ibeg = L[i];
			for (j = i + 1; j < 3; j++) {
				iend = L[j];
				if (ibeg < iend) {
					h = ig[iend];
					while (jg[h++] - ibeg);
					h--;
					ggl[h] += local_matr[i][j];
					ggu[h] += local_matr[j][i];
				}
				else {
					h = ig[ibeg];
					while (jg[h++] - iend);
					h--;
					ggl[h] += local_matr[i][j];
					ggu[h] += local_matr[j][i];
				}
			}
			di[ibeg] += local_matr[i][i];
		}
		//правая часть//
		for (i = 0; i < 3; i++)
			F[L[i]] += local_F[i];
	}

	for (k = 0; k < num_kraev; k++) {
		key = uchet_kraev(kraev[k], a, b, gamma_betta[kraev[k][0]][1]);
		if (!key) {
			K[p] = k;
			p++;
			continue;
		}
		L2[0] = kraev[k][1];
		L2[1] = kraev[k][2];
		if (key == 2) {
			for (i = 0; i < 2; i++) {
				ibeg = L2[i];
				for (j = i + 1; j < 2; j++) {
					iend = L2[j];
					if (ibeg < iend) {
						h = ig[iend];
						while (jg[h++] - ibeg);
						h--;
						ggl[h] += a[i][j];
						ggu[h] += a[j][i];
					}
					else {
						h = ig[ibeg];
						while (jg[h++] - iend);
						h--;
						ggl[h] += a[i][j];
						ggu[h] += a[j][i];
					}
				}
				di[ibeg] += a[i][i];
			}
		}
		for (i = 0; i < 2; i++)
			F[L2[i]] += b[i];

		for (i = 0; i < 2; i++) {
			b[i] = 0;
			for (j = 0; j < 2; j++)
				a[i][j] = 0;
		}
	}
	for (i = 0; i < p; i++)
		pervoe_kraevoe(kraev[K[i]]);
}

double norma(double* w)//НОРМА ВЕКТОРА
{
	double s = 0;

	for (int i = 0; i < num_points; i++)
		s += w[i] * w[i];

	return sqrt(s);
}

double calc(int i, int j, double* gl, double* gu, int kl)//по формуле это сумма которую вычитаем
{
	double s = 0;
	int k, J = jg[kl], p;
	for (k = j; k > 0; k--)
		for (p = ig[J]; p < ig[J + 1]; p++)
			if (jg[p] == jg[kl - k])
				s += gl[kl - k] * gu[p];
	return s;
}

double calcD(int j, double* gl, double* gu, int kl)//аналогично только для диагонали
{
	double s = 0;
	for (int k = kl - j; k < kl; k++)
		s += gl[k] * gu[k];
	return s;
}

void lulu(double* gl, double* gu, double* gd)
{
	int i, j, kol, kl = 0, ku = 0;
	for (i = 0; i < num_points; i++) {
		kol = ig[i + 1] - ig[i];
		for (j = 0; j < kol; j++, kl++)
			gl[kl] = (ggl[kl] - calc(i, j, gl, gu, kl)) / gd[jg[kl]];

		for (j = 0; j < kol; j++, ku++)
			gu[ku] = (ggu[ku] - calc(i, j, gu, gl, ku)) / gd[jg[ku]];

		gd[i] = sqrt(di[i] - calcD(j, gu, gl, kl));
	}
}

void Mult_A_Vect(double* xn)//перемножение исходной матрицы А на вектор
{
	long i, j, st;
	for (i = 0; i < num_points; i++) {
		f[i] = di[i] * xn[i];
		for (j = ig[i]; j < ig[i + 1]; j++) {
			st = jg[j];
			f[i] += ggl[j] * xn[st];
			f[st] += ggu[j] * xn[i];
		}
	}
}//на выходе вектор f который является глобальным

double sk_pr(double* a, double* b)//скалярное произведение векторов.
{
	double s = 0;
	for (int i = 0; i < num_points; i++)
		s += a[i] * b[i];
	return s;
}

void LOC()
{
	double nvzk, alfa, beta, skp, eps = 9.999999682655226e-030;
	int i;

	double lastnvzk;

	Mult_A_Vect(x);

	for (i = 0; i < num_points; i++)
		z[i] = r[i] = F[i] - f[i];

	Mult_A_Vect(z);

	for (i = 0; i < num_points; i++)
		p[i] = f[i];

	nvzk = sqrt(sk_pr(r, r)) / sqrt(sk_pr(F, F));

	for (int k = 1; k<10000 && nvzk > eps; k++)
	{
		lastnvzk = nvzk;
		skp = sk_pr(p, p);
		alfa = sk_pr(p, r) / skp;
		for (i = 0; i < num_points; i++) {
			x[i] += alfa * z[i];
			r[i] -= alfa * p[i];
		}

		Mult_A_Vect(r);
		beta = -sk_pr(p, f) / skp;
		for (i = 0; i < num_points; i++) {
			z[i] = r[i] + beta * z[i];
			p[i] = f[i] + beta * p[i];
		}
		nvzk = sqrt(sk_pr(r, r)) / sqrt(sk_pr(F, F));
	}
}

int main()
{
	int i = 0;
	int key = 0;

	if (input())
		cout << "Input error!" << endl;

	double* tr;

	x = new double[num_points];
	di = new double[num_points];
	F = new double[num_points];
	z = new double[num_points];
	r = new double[num_points];
	p = new double[num_points];
	l = new double[num_points];
	l1 = new double[num_points];
	f = new double[num_points];
	tr = new double[num_points];

	portret();

	for (i = 0; i < num_points; i++)
		di[i] = F[i] = x[i] = 0;

	ggu = new double[ig[num_points] - 1];
	ggl = new double[ig[num_points] - 1];

	for (i = 0; i < ig[num_points]; i++)
		ggu[i] = ggl[i] = 0;

	global_matrix();

	LOC();
	for (int i = 0; i < num_points; i++)
		tr[i] = x[i];

	tochnoe();

	ofstream file11("1.txt");
	for (i = 0; i < num_points; i++)
		file11 << setprecision(20) << x[i] << "	" << tr[i] << endl;
	file11.close();

	return 0;
}

