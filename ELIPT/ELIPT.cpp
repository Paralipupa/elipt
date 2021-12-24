#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;
class ELIPT
{
public:
	double det, mes_G, b0, b1, a0, a1, ** dot, ** params, ** M, ** G, * F, * di, * diag, *
		ggl, * ggu, * x, * z, * r, * p, * buff, * buff_2, * p12, * p23, * p31, * local_F;
	int** elem, ** kr, * fa, * pa, * ig, * jg, * m;
	int n, n_f, n_a, n_k;
	void input(ifstream& f, ifstream& f1, ifstream& f2, ifstream& f3, ifstream& f4,
		ifstream& f5);
	void memory();
	void portrait();
	void sort(int k);
	void global_m();
	int uchet_kr(int* kr, double b);
	void kraevoe1(int* kr);
	void M_local(double* v, double* v1, double* v2, double g, int num);
	void G_local(double* v, double* v1, double* v2, int k);
	double kraevye(double* x, int k, int flag);
	double function(double x, double y, int i);
	double L(double* v, double* v1, double* v2);
	double lambda(double* p1, double* p2, double* p3, int k);
	void deleteall();
	double scal(double* w, double* t);
	void multiplication_matrix_vector(double* f, double* xn);
	void vectorMinus(double* z, double* x, double* y);
	void vectorPlus(double* z, double* x, double* y);
	void vectorDiv(double* z, double* x, double* y);
	void vectorMult(double* z, double x, double* y);
	void vectorMultv(double* z, double* x, double* y);
	void LOS();
	void output(ofstream& f);
};
void ELIPT::input(ifstream& f, ifstream& f1, ifstream& f2, ifstream& f3, ifstream&
	f4, ifstream& f5)
{
	f >> n;
	f >> n_f;
	dot = new double* [n];
	for (int i = 0; i < n; i++)
	{
		dot[i] = new double[2];
		for (int j = 0; j < 2; j++)
			f1 >> dot[i][j];
	}
	elem = new int* [n_f];
	for (int i = 0; i < n_f; i++)
	{
		elem[i] = new int[3];
		for (int j = 0; j < 3; j++)
			f2 >> elem[i][j];
	}
	f3 >> n_a;
	params = new double* [n_a];
	fa = new int[n_f];
	for (int i = 0; i < n_a; i++)
	{
		params[i] = new double[2];
		for (int j = 0; j < 2; j++)
			f3 >> params[i][j];
	}
	for (int i = 0; i < n_f; i++)
		f3 >> fa[i];
	f4 >> n_k;
	kr = new int* [n_k];
	for (int i = 0; i < n_k; i++)
	{
		kr[i] = new int[4];
		for (int j = 0; j < 5; j++)
			f4 >> kr[i][j];
	}
	pa = new int[n];
	for (int i = 0; i < n; i++)
		f5 >> pa[i];
	memory();
}
void ELIPT::memory()
{
	x = new double[n];
	di = new double[n];
	diag = new double[n];
	z = new double[n];
	r = new double[n];
	p = new double[n];
	buff = new double[n];
	buff_2 = new double[n];
	F = new double[n];
	ig = new int[n + 1];
	M = new double* [3];
	G = new double* [3];
	p12 = new double[2];
	p23 = new double[2];
	p31 = new double[2];
	local_F = new double[3];
	m = new int[n];
	for (int i = 0; i < n; i++)
		di[i] = F[i] = x[i] = 0;
	for (int i = 0; i < 3; i++)
	{
		M[i] = new double[3];
		G[i] = new double[3];
	}
}
void ELIPT::portrait()
{
	int s = 0, st = 0, f = 0, position = 0;
	struct List
	{
		int N;
		List* next;
	};
	List* list = new List[n];
	List* a = NULL;
	for (int i = 0; i < n; i++)
		list[i].next = NULL;
	for (int i = 0; i < n_f; i++)
		for (int j = 0; j < 3; j++)
		{
			f = 0;
			s = elem[i][(j == 2)];
			st = elem[i][((j != 0) + 1)];
			if (s < st)
			{
				s += st;
				st = s - st;
				s -= st;
			}
			a = &list[s];
			while (a->next)
			{
				if (a->next->N == st)
				{
					f = 1;
					break;
				}
				a = a->next;
			}
			if (!f)
			{
				a->next = new List;
				a->next->N = st;
				a->next->next = NULL;
			}
		}
	ig[0] = 0;
	for (int i = 0; i < n; i++)
	{
		s = 0;
		a = &list[i];
		while (a = a->next)
			s++;
		ig[i + 1] = ig[i] + s;
	}
	jg = new int[ig[n]];
	for (int i = 0; i < n; i++)
	{
		s = 0;
		f = 0;
		a = &list[i];
		while (a = a->next)
		{
			m[s] = a->N;
			s++;
			f = 1;
		}
		if (f)
		{
			sort(--s);
			for (int i = position, j = 0; i <= s + position; i++, j++)
				jg[i] = m[j];
				position += s + 1;
		}
	}
	ggu = new double[ig[n]];
	ggl = new double[ig[n]];
	for (int i = 0; i < ig[n]; i++)
		ggu[i] = ggl[i] = 0;
}
void ELIPT::sort(int k)
{
	for (int i = k - 1; i >= 0; i--)
		for (int j = 0; j <= i; j++)
			if (m[j] > m[j + 1])
			{
				int l = m[j];
				m[j] = m[j + 1];
				m[j + 1] = l;
			}
}
void ELIPT::global_m()
{
	int s = 0, beg = 0, end = 0, h = 0, f = 0;
	for (int k = 0; k < n_f; k++)
	{
		G_local(dot[elem[k][0]], dot[elem[k][1]], dot[elem[k][2]], fa[k]);
		M_local(dot[elem[k][0]], dot[elem[k][1]], dot[elem[k][2]], params[fa[k]][0],
			fa[k]);
		for (int i = 0; i < 3; i++)
		{
			beg = elem[k][i];
			for (int j = i + 1; j < 3; j++)
			{
				end = elem[k][j];
				if (beg < end)
				{
					h = ig[end];
					while (jg[h++] - beg);
					h--;
					ggl[h] += M[i][j] + G[i][j];
					ggu[h] += M[j][i] + G[j][i];
				}
				else
				{
					h = ig[beg];
					while (jg[h++] - end);
					h--;
					ggl[h] += M[i][j] + G[i][j];
					ggu[h] += M[j][i] + G[j][i];
				}
			}
			di[beg] += M[i][i] + G[i][i];
		}
		for (int i = 0; i < 3; i++)
			F[elem[k][i]] += local_F[i];
	}
	for (int k = 0; k < n_k; k++)
	{
		elem[n_f - 1][k] = 0;
		f = uchet_kr(kr[k], params[kr[k][0]][1]);
		if (!f)
		{
			elem[n_f - 1][s] = k;
			s++;
			continue;
		}
		if (f == 2)
			for (int i = 0; i < 2; i++)
			{
				beg = kr[k][i + 1];
				for (int j = i + 1; j < 2; j++)
				{
					end = kr[k][2];
					if (beg < end)
					{
						h = ig[end];
						while (jg[h++] - beg);
						h--;
						ggl[h] += a1;
						ggu[h] += a1;
					}
					else
					{
						h = ig[beg];
						while (jg[h++] - end);
						h--;
						ggl[h] += a1;
						ggu[h] += a1;
					}
				}
				di[beg] += a0;
			}
		F[kr[k][1]] += b0;
		F[kr[k][2]] += b1;
	}
	for (int i = 0; i < s; i++)
		kraevoe1(kr[elem[n_f - 1][i]]);
}
int ELIPT::uchet_kr(int* kr, double b)
{
	mes_G = sqrt((dot[kr[2]][0] - dot[kr[1]][0]) * (dot[kr[2]][0] - dot[kr[1]][0]) +
		(dot[kr[2]][1] - dot[kr[1]][1]) * (dot[kr[2]][1] - dot[kr[1]][1]));
	if (kr[3] == 1)
		return 0;
	else
		if (kr[3] == 2)
		{
			double coeff = mes_G / 6.;
			b0 = coeff * (2 * kraevye(dot[kr[1]], kr[4], 2) + kraevye(dot[kr[2]], kr[4],
				2));
			b1 = coeff * (kraevye(dot[kr[1]], kr[4], 2) + 2 * kraevye(dot[kr[2]], kr[4],
				2));
			return 1;
		}
		else
		{
			double coeff = b * mes_G / 6.;
			a0 = 2 * coeff;
			a1 = coeff;
			b0 = coeff * (2 * kraevye(dot[kr[1]], kr[4], 3) + kraevye(dot[kr[2]], kr[4],
				3));
			b1 = coeff * (kraevye(dot[kr[1]], kr[4], 3) + 2 * kraevye(dot[kr[2]], kr[4],
				3));
			return 2;
		}
}
void ELIPT::kraevoe1(int* kr)
{
	for (int i = 1; i <= 2; i++)
	{
		di[kr[i]] = 1;
		F[kr[i]] = kraevye(dot[kr[i]], kr[4], 1);
	}
	for (int i = 0; i < ig[kr[1] + 1] - ig[kr[1]]; i++)
		ggl[ig[kr[1]] + i] = 0;
	for (int i = 0; i < ig[kr[2] + 1] - ig[kr[2]]; i++)
		ggl[ig[kr[2]] + i] = 0;
	for (int i = kr[1] + 1; i < n; i++)
		for (int p = ig[i]; p < ig[i + 1]; p++)
			if (jg[p] == kr[1])
			{
				ggu[p] = 0;
				continue;
			}
	for (int i = kr[2] + 1; i < n; i++)
		for (int p = ig[i]; p < ig[i + 1]; p++)
			if (jg[p] == kr[2])
			{
				ggu[p] = 0;
				continue;
			}
}
double ELIPT::L(double* v, double* v1, double* v2)
{
	return ((v[0] * v1[1] - v1[0] * v[1] + (v[1] - v1[1]) * v2[0] + (v1[0] - v[0]) * v2[1]) /
		det);
}
double ELIPT::lambda(double* v, double* v1, double* v2, int k)
{
	double l = (v[0] * v[0]) * L(v1, v2, v) * (2 * L(v1, v2, v) - 1) +
		+(p12[0] * p12[0]) * L(v2, v, p12) * (2 * L(v2, v, p12) - 1) +
		+(v1[0] * v1[0]) * L(v, v1, v1) * (2 * L(v, v1, v1) - 1) +
		+(p23[0] * p23[0]) * 4 * L(v1, v2, p23) * L(v2, v, p23) +
		+(v2[0] * v2[0]) * 4 * L(v2, v, v2) * L(v, v1, v2) +
		+(p31[0] * p31[0]) * 4 * L(v1, v2, p31) * L(v, v1, p31);
	//return l; //если лямбда представлена как х*х
	return 10; //если лямбда число
}
void ELIPT::M_local(double* v, double* v1, double* v2, double g, int num)
{
	double mnoz = fabs(det) / 24.;
	double c1 = mnoz * function(v[0], v[1], num);
	double c2 = mnoz * function(v1[0], v1[1], num);
	double c3 = mnoz * function(v2[0], v2[1], num);
	local_F[0] = 2 * c1 + c2 + c3;
	local_F[1] = c1 + 2 * c2 + c3;
	local_F[2] = c1 + c2 + 2 * c3;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			if (i == j)
				M[i][j] = 2 * mnoz * g;
			else
				M[i][j] = mnoz * g;
}
void ELIPT::G_local(double* v, double* v1, double* v2, int k)
{
	for (int i = 0; i < 2; i++)
	{
		p12[i] = (v[i] + v1[i]) / 2;
		p23[i] = (v1[i] + v2[i]) / 2;
		p31[i] = (v2[i] + v[i]) / 2;
	}
	double a = v1[1] - v2[1];
	double b = v2[0] - v1[0];
	double c = v2[1] - v[1];
	double d = v[0] - v2[0];
	double e = v[1] - v1[1];
	double f = v1[0] - v[0];
	det = f * c - (v2[0] - v[0]) * (v1[1] - v[1]);
	double coeff = lambda(v, v1, v2, k) / (6 * fabs(det));
	G[0][0] = coeff * (a * a + b * b);
	G[0][1] = coeff * (a * c + b * d);
	G[0][2] = coeff * (a * e + b * f);
	G[1][0] = coeff * (c * a + d * b);
	G[1][1] = coeff * (c * c + d * d);
	G[1][2] = coeff * (c * e + d * f);
	G[2][0] = coeff * (a * e + b * f);
	G[2][1] = coeff * (c * e + d * f);
	G[2][2] = coeff * (e * e + f * f);
}
double ELIPT::kraevye(double* x, int k, int flag)
{
	switch (flag)
	{
	case 1: //первые краевые
		switch (k)
		{
		case 0: return 0;
		case 1: return sin(x[0]);
		case 2: return 0;
		case 3: return sin(x[1]);
		}
	case 2: //вторые краевые
		switch (k)
		{
		case 0: return x[0] / 400.;
		case 1: return 0;
		default: return 0;
		case 3: //трерьи краевые
			if (!k)
				return 81 / x[1] / x[1];
			return 0;
		}
	}
	return 0;
}
double ELIPT::function(double x, double y, int i)
{
	if (!i)
		return x * x - 2;
	return 0;
}
void ELIPT::deleteall()
{
	delete[] x;
	delete[] di;
	delete[] diag;
	delete[] z;
	delete[] r;
	delete[] p;
	delete[] buff;
	delete[] buff_2;
	delete[] F;
	delete[] ig;
	delete[] M;
	delete[] G;
	delete[] p12;
	delete[] p23;
	delete[] p31;
	delete[] local_F;
	delete[] m;
}
//////решатель/////////////////////////////////////////////////////////////////////////////////////
double ELIPT::scal(double* w, double* t)
{
	double s = 0;
	for (int i = 0; i < n; i++)
		s += w[i] * t[i];
	return s;
}
void ELIPT::multiplication_matrix_vector(double* f, double* xn)
{
	for (int i = 0; i < n; i++)
	{
		f[i] = di[i] * xn[i];
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			f[i] += ggl[j] * xn[jg[j]];
			f[jg[j]] += ggu[j] * xn[i];
		}
	}
}
void ELIPT::vectorMinus(double* z, double* x, double* y)
{
	for (int i = 0; i < n; i++)
		z[i] = x[i] - y[i];
}
void ELIPT::vectorPlus(double* z, double* x, double* y)
{
	for (int i = 0; i < n; i++)
		z[i] = x[i] + y[i];
}
void ELIPT::vectorDiv(double* z, double* x, double* y)
{
	for (int i = 0; i < n; i++)
		z[i] = x[i] / y[i];
}
void ELIPT::vectorMult(double* z, double x, double* y)
{
	for (int i = 0; i < n; i++)
		z[i] = x * y[i];
}
void ELIPT::vectorMultv(double* z, double* x, double* y)
{
	for (int i = 0; i < n; i++)
		z[i] = x[i] * y[i];
}
void ELIPT::LOS() //ЛОС c диаг.предобусл.
{
	double eps = 1e-40;
	for (int i = 0; i < n; ++i)
		diag[i] = sqrt(di[i]);
	multiplication_matrix_vector(buff, x);
	vectorMinus(r, F, buff);
	vectorDiv(r, r, diag);
	vectorDiv(z, r, diag);
	multiplication_matrix_vector(p, z);
	vectorDiv(p, p, diag);
	double normf = sqrt(scal(F, F));
	double nev = 1000;
	for (int i = 0; nev > eps && i <= 1000; ++i)
	{
		double length = scal(p, p);
		double alpha = scal(p, r) / length;
		vectorMult(buff, alpha, z);
		vectorPlus(x, x, buff);
		vectorMult(buff, alpha, p);
		vectorMinus(r, r, buff);
		vectorDiv(r, r, diag);
		multiplication_matrix_vector(buff, r);
		vectorDiv(buff, buff, diag);
		double beta = -scal(p, buff) / length;
		vectorMult(z, beta, z);
		vectorPlus(z, r, z);
		vectorMult(p, beta, p);
		vectorPlus(p, buff, p);
		vectorMultv(r, r, diag);
		nev = sqrt(scal(r, r)) / normf;
	}
}
void ELIPT::output(ofstream& f)
{
	for (int i = 0; i < n; i++)
		f << setprecision(20) << x[i] << endl;
}
void main()
{
	ELIPT A;
	ifstream f("info.txt");
	ifstream f1("xy.txt");
	ifstream f2("elem.txt");
	ifstream f3("mat.txt");
	ifstream f4("kraevye.txt");
	ifstream f5("in_area.txt");
	A.input(f, f1, f2, f3, f4, f5);
	f.close();
	f1.close();
	f2.close();
	f3.close();
	f4.close();
	A.portrait();
	A.global_m();
	A.LOS();
	ofstream fout("1.txt");
	A.output(fout);
	fout.close();
	A.deleteall();
}
