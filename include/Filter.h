#ifndef FILTER_H_
#define FILTER_H_
#include "tool.h"
inline void movingAverage2D(double** grad, int nx, int ny, int wx, int wy, int times)
{
	//expand
	double** grad_ex;
	int bond = 30;
	int nx1 = nx + 2 * bond;
	int ny1 = ny + 2 * bond;
	grad_ex = alloc_double_2d(nx1, ny1);
	exmodel(grad, grad_ex, nx, ny, bond);
	for (int it = 0; it < times; it++)
	{
		for (int i = 0; i < nx1; ++i)
		{
			for (int j = 0; j < ny1; ++j)
			{
				double sum = 0.0;
				int count = 0;
				for (int ii = max(0, i - wx / 2); ii <= min(nx1 - 1, i + wx / 2); ++ii)
				{
					for (int jj = max(0, j - wy / 2); jj <= min(ny1 - 1, j + wy / 2); ++jj)
					{
						sum += grad_ex[ii][jj];
						count++;
					}
				}
				grad_ex[i][j] = sum / count;
			}
		}
	}
	for (int ix = bond; ix < nx + bond; ix++)
		for (int iy = bond; iy < ny + bond; iy++)
		{
			grad[ix - bond][iy - bond] = grad_ex[ix][iy];
		}
	free_2d(grad_ex, nx1);
}
inline void movingAverage2D_vec(vector<vector<double>>&grad1, int nx, int ny, int wx, int wy, int times)
{
	//expand
	double** grad;
	grad = alloc_double_2d(nx, ny);
	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{
			grad[i][j] = grad1[i][j];
		}
	}
	double** grad_ex;
	int bond = 30;
	int nx1 = nx + 2 * bond;
	int ny1 = ny + 2 * bond;
	grad_ex = alloc_double_2d(nx1, ny1);
	exmodel(grad, grad_ex, nx, ny, bond);
	for (int it = 0; it < times; it++)
	{
		for (int i = 0; i < nx1; ++i)
		{
			for (int j = 0; j < ny1; ++j)
			{
				double sum = 0.0;
				int count = 0;
				for (int ii = max(0, i - wx / 2); ii <= min(nx1 - 1, i + wx / 2); ++ii)
				{
					for (int jj = max(0, j - wy / 2); jj <= min(ny1 - 1, j + wy / 2); ++jj)
					{
						sum += grad_ex[ii][jj];
						count++;
					}
				}
				grad_ex[i][j] = sum / count;
			}
		}
	}
	for (int ix = bond; ix < nx + bond; ix++)
		for (int iy = bond; iy < ny + bond; iy++)
		{
			grad[ix - bond][iy - bond] = grad_ex[ix][iy];
		}
	free_2d(grad_ex, nx1);
}
inline void movingAverage1D(double* grad, int nx, int ny, int wx, int wy, int times)
{
	//expand
	double** tmp;
	tmp = alloc_double_2d(nx, ny);
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
		{
			tmp[i][j] = grad[j + i * ny];
		}
	//
	double** grad_ex;
	int bond = 30;
	int nx1 = nx + 2 * bond;
	int ny1 = ny + 2 * bond;
	grad_ex = alloc_double_2d(nx1, ny1);
	exmodel(tmp, grad_ex, nx, ny, bond);
	for (int it = 0; it < times; it++)
	{
		for (int i = 0; i < nx1; ++i)
		{
			for (int j = 0; j < ny1; ++j)
			{
				double sum = 0.0;
				int count = 0;
				for (int ii = max(0, i - wx / 2); ii <= min(nx1 - 1, i + wx / 2); ++ii)
				{
					for (int jj = max(0, j - wy / 2); jj <= min(ny1 - 1, j + wy / 2); ++jj)
					{
						sum += grad_ex[ii][jj];
						count++;
					}
				}
				grad_ex[i][j] = sum / count;
			}
		}
	}
	for (int ix = bond; ix < nx + bond; ix++)
		for (int iy = bond; iy < ny + bond; iy++)
		{
			tmp[ix - bond][iy - bond] = grad_ex[ix][iy];
		}
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
		{
			grad[j + i * ny] = tmp[i][j];
		}
	free_2d(grad_ex, nx1);
	free_2d(tmp, nx);
}
inline void Gauss2D(double** mod, int nx, int nz, int h_size, double sigma)
{
	int bond, ix, iz, i, j;
	bond = (h_size - 1) / 2;
	double** exmod;
	exmod = alloc_double_2d(nx + 2 * bond, nz + 2 * bond);
	exmodel(mod, exmod, nx, nz, bond);
	double** gau, tmp;
	gau = alloc_double_2d(h_size, h_size);
	MakeGauss(gau, h_size, sigma);
	for (ix = bond; ix < nx + bond; ix++)
		for (iz = bond; iz < nz + bond; iz++)
		{
			tmp = 0;
			for (i = 0; i < h_size; i++)
				for (j = 0; j < h_size; j++)
					tmp += gau[i][j] * exmod[ix - bond + i][iz - bond + j];
			mod[ix - bond][iz - bond] = tmp;
		}
	free_2d(exmod, nx + 2 * bond);
	free_2d(gau, h_size);
}
inline void smoothgauss2dtopo(double** v, int nz, int nx, double dz, double freq, vector<int>& itopo, int maxn, double fracx, double fracz, double vref, int nptt)
{
	double* beta1, * beta2, ** vf;
	double wl, taux, tauz, taux2, tauz2, h, xl2, xl1, d, betatot;
	int ix, iz, il2, il1, k1, k2, j1, j2, jj1, jj2;

	beta1 = new double[maxn];
	beta2 = new double[maxn];
	vf = new double* [nx];
	for (ix = 0; ix < nx; ix++)vf[ix] = new double[nz];

	// Correlation lengths of Gaussian filter are defined as fraction of the wavelength 
	// as defined by the inverted frequency component and a reference velocity

	wl = vref / freq;
	taux = fracx * wl;
	tauz = fracz * wl;
	h = dz;
	if (taux == 0 && tauz == 0)
		for (ix = 0; ix < nx; ix++)
			for (iz = 0; iz < nz; iz++)
				vf[ix][iz] = v[ix][iz];
	else
	{
		for (ix = 0; ix < nx; ix++)
			for (iz = 0; iz < itopo[ix] + nptt; iz++)
				vf[ix][iz] = v[ix][iz];
		xl2 = 3. * taux;
		xl1 = 3. * tauz;
		il2 = int(xl2 / h) + 1;
		il1 = int(xl1 / h + 1);
		if (2 * il1 + 1 > maxn)
		{
			cout << "*******increase maxn in smoothgauss*******" << endl;
			exit(0);
		}
		if (2 * il2 + 1 > maxn)
		{
			cout << "*******increase maxn in smoothgauss*******" << endl;
			exit(0);
		}
		tauz2 = tauz * tauz;
		taux2 = taux * taux;
		k2 = 0;
		for (ix = -il2; ix <= il2; ix++)
		{
			d = ix * h;
			beta2[k2] = exp(-d * d / taux2);
			k2++;
		}
		k1 = 0;
		for (iz = -il1; iz <= il1; iz++)
		{
			d = iz * h;
			beta1[k1] = exp(-d * d / tauz2);
			k1++;
		}

		for (ix = 0; ix < nx; ix++)
			for (iz = itopo[ix] + nptt; iz < nz; iz++)
			{
				vf[ix][iz] = 0;
				betatot = 0;
				jj2 = -1;
				for (j2 = -il2; j2 <= il2; j2++)
				{
					k2 = j2 + ix;
					jj2++;
					if (k2 < 0 || k2 >= nx)continue;
					jj1 = -1;
					for (j1 = -il1; j1 <= il1; j1++)
					{
						k1 = j1 + iz;
						jj1++;
						if (k1 < itopo[ix] + nptt || k1 >= nz)continue;
						vf[ix][iz] += beta1[jj1] * beta2[jj2] * v[k2][k1];
						betatot += beta1[jj1] * beta2[jj2];
					}
				}
				vf[ix][iz] /= betatot;
			}
	}

	for (ix = 0; ix < nx; ix++)
		for (iz = 0; iz < nz; iz++)
			v[ix][iz] = vf[ix][iz];

	delete[] beta1;
	delete[] beta2;
	for (ix = 0; ix < nx; ix++)delete[] vf[ix];
	delete[] vf;
}

inline void horizontal_smooth_fast(double** grad, int nx, int ny, int smx, const vector<int>& topo_ynum, double** output) {
	int half = smx / 2;
	for (int iy = 0; iy < ny; ++iy) {
		// 前缀和数组
		vector<double> sum_pre(nx + 1, 0.0);
		vector<int> cnt_pre(nx + 1, 0);
		for (int ix = 0; ix < nx; ++ix) {
			double val = (iy >= topo_ynum[ix]) ? grad[ix][iy] : 0.0;
			int cnt = (iy >= topo_ynum[ix]) ? 1 : 0;
			sum_pre[ix + 1] = sum_pre[ix] + val;
			cnt_pre[ix + 1] = cnt_pre[ix] + cnt;
		}
		// 对每个地下网格计算窗口均值
		for (int ix = 0; ix < nx; ++ix) {
			if (iy < topo_ynum[ix]) {
				output[ix][iy] = 0.0;
				continue;
			}
			int L = max(0, ix - half);
			int R = min(nx - 1, ix + half);
			double sum = sum_pre[R + 1] - sum_pre[L];
			int cnt = cnt_pre[R + 1] - cnt_pre[L];
			if (cnt > 0) output[ix][iy] = sum / cnt;
			else output[ix][iy] = grad[ix][iy];
		}
	}
}

inline void vertical_smooth_fast(double** grad, int nx, int ny, int smy, const vector<int>& topo_ynum, double** output) {
	int half = smy / 2;
	for (int ix = 0; ix < nx; ++ix) {
		// 前缀和数组
		vector<double> sum_pre(ny + 1, 0.0);
		vector<int> cnt_pre(ny + 1, 0);
		for (int iy = 0; iy < ny; ++iy) {
			double val = (iy >= topo_ynum[ix]) ? grad[ix][iy] : 0.0;
			int cnt = (iy >= topo_ynum[ix]) ? 1 : 0;
			sum_pre[iy + 1] = sum_pre[iy] + val;
			cnt_pre[iy + 1] = cnt_pre[iy] + cnt;
		}
		for (int iy = 0; iy < ny; ++iy) {
			if (iy < topo_ynum[ix]) {
				output[ix][iy] = 0.0;
				continue;
			}
			int L = max(0, iy - half);
			int R = min(ny - 1, iy + half);
			double sum = sum_pre[R + 1] - sum_pre[L];
			int cnt = cnt_pre[R + 1] - cnt_pre[L];
			if (cnt > 0) output[ix][iy] = sum / cnt;
			else output[ix][iy] = grad[ix][iy];
		}
	}
}

inline void movingAverage2D_topo(double** grad, int nx, int ny, int smx, int smy, int smtimes, const vector<int>& topo_ynum) {
	double** tmp = alloc_double_2d(nx, ny);
	for (int iter = 0; iter < smtimes; ++iter) {
		// 水平平滑
		horizontal_smooth_fast(grad, nx, ny, smx, topo_ynum, tmp);
		// 垂直平滑
		vertical_smooth_fast(tmp, nx, ny, smy, topo_ynum, grad);
	}
	free_2d(tmp, nx);
}
#endif