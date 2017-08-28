#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

const double c = 2.9979246e8; //! 光速
const double eps0 = 8.8541878e-12; //! 真空の誘電率
const double mu0 = 1.256637e-6; //! 真空の透磁率
const double z0 = 376.73031; //! 真空のインピーダンス

const std::string filename = "../gnuplot_src/fdtd_data.dat";

const int nz0 = 1000; //! 解析領域の分割数
const int nstep = 2000; //! 計算のステップ数
const int lpml = 8; //! PMLの総数
const int order = 4; //! PMLの次数
const double rmax = -120; //! PMLの要求精度dB
const double copml = -1.5280063e-4; // ln(10)/(40Z0)
const int nz = nz0 + 2 * lpml; //! 全分割数

const double d = 0.1; //! 誘電体厚さ [m]
const double epsr = 3.0; //! 誘電体の比誘電率
const int nd = 50; //! 誘電体の分割数
const int kd = nz / 2; //! 誘電体の位置
const double dz = d / nd; //! セルサイズ
const double dt = 0.99999*dz / c; //! 時間ステップ [s]
const double zp = (lpml + 100.0)*dz; //! 初期パルスの初期位置
const double a = 20.0*dz; //! 初期パルスの幅

std::vector<double> ex(nz + 1, 0.0), hy(nz + 1, 0.0); //! 電界, 磁界を記憶する配列
std::vector<double> ae(nz + 1, 0.0), be(nz + 1, 0.0), am(nz + 1, 0.0), bm(nz + 1, 0.0); //! 電界磁界更新式の係数
double z, t; //! 位置, 時刻

std::vector<double> epsd(nz + 1, 0.0), sgmd(nz + 1, 0.0), mud(nz + 1, 0.0), msgmd(nz + 1, 0.0);
double eps, sgm, mu, msgm; //! 領域内の誘電率, 導電率, 透磁率, 導磁率
double eps_l,eps_r,smax0e,sgme,sgmm;

void setup(void);
void pmlcoef(void);
double pulse(double z, double tm);
void e_cal(void);
void h_cal(void);

int main() {

	std::ofstream outputfile(filename);

	setup();
	pmlcoef();

	t = dt;

	for (int n = 1; n <= nstep; n++) {
		e_cal();
		t += 0.5*dt;
		h_cal();
		t += 0.5*dt;

		if (n % 10 == 0) {
			for (int k = 0; k <= nz; k++) {
				if (k % 10 == 0) {
					z = k*dz;
					outputfile << z << " " << ex[k] << std::endl;
				}
			}
			outputfile << std::endl << std::endl;
			std::cout << n << "/" << nstep << std::endl;
		}
		
	}

	outputfile.close();

	return 0;
}

void setup(void) {

	for (int k = 0; k <= nz; k++) {
		if ((k < kd) || (k >= kd + nd)) {
			epsd[k] = 1.0;
			sgmd[k] = 0.0;
			mud[k] = 1.0;
			msgmd[k] = 0.0;
		}
		else {
			epsd[k] = epsr;
			sgmd[k] = 0.0;
			mud[k] = 1.0;
			msgmd[k] = 0.0;
		}
	}

	for (int k = 1; k <= nz; k++) {
		eps = 0.5*(epsd[k] + epsd[k - 1])*eps0;
		sgm = 0.5*(sgmd[k] + sgmd[k - 1]);
		mu = mud[k] * mu0;
		msgm = msgmd[k];

		ae[k] = (1.0 - (sgm*dt / (2.0*eps))) / (1.0 + (sgm*dt / (2.0*eps)));
		be[k] = (dt/eps) / (1.0 + (sgm*dt / (2.0*eps))) / dz;

		am[k] = (1.0 - (msgm*dt / (2.0*mu))) / (1.0 + (msgm*dt / (2.0*mu)));
		bm[k] = (dt/mu) / (1.0 + (msgm*dt / (2.0*mu))) / dz;
	}

	for (int k = 0; k <= nz; k++) {
		z = k * dz;
		ex[k] = pulse(z, 0.0);
	}

	for (int k = 0; k <= nz; k++) {
		z = (k+0.5) * dz;
		hy[k] = pulse(z, 0.5*dt)/z0;
	}
}

void pmlcoef(void) {
	eps_l = epsd[lpml + 1] *eps0;
	eps_r = epsd[nz - lpml - 1] * eps0;

	smax0e = copml * rmax * (order + 1) / (lpml*dz);

	for (int k = 0; k <= lpml - 1; k++) {
		sgme = pow(((double)(lpml - k) / (double)(lpml)), order)*smax0e*epsd[lpml + 1];
		sgmm = pow((((double)(lpml - k)-0.5) / (double)(lpml)), order)*smax0e*epsd[lpml + 1];

		ae[k] = (1.0 - ((sgme*dt) / (2.0*eps_l))) / (1.0 + ((sgme*dt) / (2.0*eps_l)));
		be[k] = dt / eps_l / (1.0 + ((sgme*dt) / (2.0*eps_l))) / dz;

		am[k] = (1.0 - ((sgmm*dt) / (2.0*eps_l))) / (1.0 + ((sgmm*dt) / (2.0*eps_l)));
		bm[k] = dt / (mud[lpml + 1] * mu0) / (1.0 + ((sgmm*dt) / (2.0*eps_l))) / dz;
	}

	for (int k = nz-lpml; k <= nz-1; k++) {
		sgme = pow(((double)(k-nz+lpml) / (double)(lpml)), order)*smax0e*epsd[nz - lpml - 1];
		sgmm = pow((((double)(k-nz+lpml) + 0.5) / (double)(lpml)), order)*smax0e*epsd[nz - lpml - 1];

		ae[k] = (1.0 - ((sgme*dt) / (2.0*eps_r))) / (1.0 + ((sgme*dt) / (2.0*eps_r)));
		be[k] = (dt / eps_r) / (1.0 + ((sgme*dt) / (2.0*eps_r))) / dz;

		am[k] = (1.0 - ((sgmm*dt) / (2.0*eps_r))) / (1.0 + ((sgmm*dt) / (2.0*eps_r)));
		bm[k] = (dt / (mud[nz - lpml - 1] * mu0)) / (1.0 + ((sgmm*dt) / (2.0*eps_r))) / dz;
	}
}

double pulse(double z, double tm) {
	return exp(-1.0*pow((z-zp-c*tm)/a,2));
}

void e_cal(void) {
	for (int k = 1; k <= nz - 1; k++)
		ex[k] = ae[k] * ex[k] - be[k] * (hy[k] - hy[k - 1]);
}

void h_cal(void) {
	for (int k = 0; k <= nz - 1; k++)
		hy[k] = am[k] * hy[k] - bm[k] * (ex[k+1] - ex[k]);
}