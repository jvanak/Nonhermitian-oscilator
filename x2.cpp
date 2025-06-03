#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>

#define _USE_MATH_DEFINES

using namespace std;
using ldouble = long double;
using lcomplex = complex<long double>;

ldouble theta = 0.05L;
ldouble q = 1.L;

ldouble om = 0.057L;
ldouble fr = 0.06L;
ldouble sigma = 0.00001L;
// ldouble sigma = 1.L;
ldouble c = 137.03599L;
ldouble mu = fr / c;

int n_0 = 20000000;
int n = 300000;
ldouble a_0 = 10.L;
ldouble a_1 = 0.00015L;
ldouble a_2 = 0.00015L;
int n_plot = 4000;
ldouble a_plot = 10000.L;
ldouble sw = 0.0001L;

ldouble t = -a_plot;
ldouble q2 = q * q;

/*ldouble sgn_re(lcomplex x) {
    if (real(x) > 0.0L)
        return 1.0L;
    if (real(x) < 0.0L)
        return -1.0L;
    return real(x);
}
ldouble sgn_im(lcomplex x) {
    if (imag(x) > 0.0L)
        return 1.0L;
    if (imag(x) < 0.0L)
        return -1.0L;
    return imag(x);
}*/
lcomplex integrate_0(lcomplex f(ldouble x)) {
    ldouble step = 2.L * a_0 / n_0;
    lcomplex result = 0.L;

    for (int l = 0; l <= n_0; l++) {
        result += f(l * step) * step;
    }
    return 2.L * result;
}
lcomplex integrate_1(lcomplex f(ldouble x)) {
    ldouble zn = 1.L;
    ldouble step = 2.L * a_1 / n;
    lcomplex result = 0.L;

    for (int l = 0; l <= n; l++) {
        /*  if ((f(-a_1 + mu + l * step) == 0.L || sgn_im(f(-a_1 + mu + l * step)) != sgn_im(f(-a_1 + mu + (l - 1.L) * step)) || sgn_re(f(-a_1 + mu + l * step)) != sgn_re(f(-a_1 + mu + (l - 1.L) * step))) && (abs(real(f(-a_1 + mu + l * step)) - real(f(-a_1 + mu + (l - 1.L) * step))) > sw || abs(imag(f(-a_1 + mu + l * step)) - imag(f(-a_1 + mu + (l - 1.L) * step))) > sw)) {
              zn = -zn;
          }*/
        if (abs(f(-a_1 + mu + l * step) - f(-a_1 + mu + (l - 1.L) * step)) > sw && abs(f(-a_1 + mu + l * step) - f(-a_1 + mu + (l - 1.L) * step)) > abs(f(-a_1 + mu + l * step) + f(-a_1 + mu + (l - 1.L) * step))) {
            zn = -zn;
        }
        result += zn * f(-a_2 + mu + l * step) * step;
    }
    return result;
}
lcomplex integrate_2(lcomplex f(ldouble x)) {
    ldouble zn = 1.L;
    ldouble step = 2.L * a_2 / n;
    lcomplex result = 0.L;

    for (int l = 0; l <= n; l++) {
        if (abs(f(-a_2 + mu + l * step) - f(-a_2 + mu + (l - 1.L) * step)) > sw && abs(f(-a_2 + mu + l * step) - f(-a_2 + mu + (l - 1.L) * step)) > abs(f(-a_2 + mu + l * step) + f(-a_2 + mu + (l - 1.L) * step))) {
            zn = -zn;
        }
        result += zn * f(-a_2 + mu + l * step) * step;
    }
    return result;
}
lcomplex z_0(ldouble k) {
    return k / (q2 * q2 * exp(-2.il * theta) * c * c * k * k + 4.L * M_PI * M_PI * c * c * (om * om - exp(-2.il * theta) * c * c * k * k) * (om * om - exp(-2.il * theta) * c * c * k * k));
}
lcomplex z_1(ldouble k) {
    return exp(1.il * exp(-1.il * theta) * c * fabs(k) * t - 3.il * theta / 2.L - (k * exp(-1.il * theta) - mu) * (k * exp(-1.il * theta) - mu) / (2.L * sigma * sigma)) * sqrt(fabs(k) / (exp(-2.il * theta) * c * c * k * k * q2 * q2 / 4.L + M_PI * M_PI * c * c * (om * om - exp(-2.il * theta) * c * c * k * k) * (om * om - exp(-2.il * theta) * c * c * k * k)));
}
lcomplex z_2(ldouble k) {
    return exp(-1.il * exp(-1.il * theta) * c * fabs(k) * t - 3.il * theta / 2.L - (k * exp(-1.il * theta) - mu) * (k * exp(-1.il * theta) - mu) / (2.L * sigma * sigma)) * sqrt(fabs(k) / (exp(-2.il * theta) * c * c * k * k * q2 * q2 / 4.L + M_PI * M_PI * c * c * (om * om - exp(-2.il * theta) * c * c * k * k) * (om * om - exp(-2.il * theta) * c * c * k * k)));
}
int main() {
    /*
        t = -472.L;
        ldouble a = a_1;
        ldouble kstep_plot = 2.L * a / n;
        lcomplex int1;
        ldouble znam = -1.L;

        ofstream kout;
        kout.open("int1.dat");
        ofstream kout2;
        kout2.open("int1_re.dat");
        ofstream kout3;
        kout3.open("int1_im.dat");

        for (int j = 0; j <= n; j++) {
            ldouble k = mu + j * kstep_plot - a;

            if (abs(z_1(k) - z_1(k - kstep_plot)) > sw && abs(z_1(k) - z_1(k - kstep_plot)) > abs(z_1(k) + z_1(k - kstep_plot))) {
                znam = -znam;
            }
            int1 = znam * z_1(k);

            kout << k << ", " << setprecision(20) << fabs(int1) << endl;
            kout2 << k << ", " << setprecision(20) << real(int1) << endl;
            kout3 << k << ", " << setprecision(20) << imag(int1) << endl;
        }
        kout.close();
        kout2.close();
        kout3.close();

        cout << "int plot done" << endl;
    */
    lcomplex Omega = sqrt(om * om - q2 * q2 / (16.L * M_PI * M_PI * c * c)) - 1il * q2 / (4.L * M_PI * c);
    ldouble K = q2 * c * c * c / (2.L * sqrt(M_PI) * sigma);
    lcomplex x2_0 = Omega / (2.L * om * om - 1.il * q2 * Omega / (2.L * M_PI * c)) + q2 * c * c * c * exp(-2.il * theta) * integrate_0(z_0);

    cout << "arg(Omega) = " << arg(Omega) << endl;
    cout << "mu = " << mu << endl;
    cout << "K = " << K << endl;
    cout << "x^2_0 = " << x2_0 << endl;

    ldouble step_plot = 2.L * a_plot / n_plot;
    lcomplex x2;
    lcomplex int1_full;
    lcomplex int2_full;

    ofstream out;
    out.open("x2_t.dat");
    ofstream out2;
    out2.open("x2_t_re.dat");
    ofstream out3;
    out3.open("x2_t_im.dat");
    ofstream out4;
    out4.open("int1_full_re.dat");
    ofstream out5;
    out5.open("int1_full_im.dat");
    ofstream out6;
    out6.open("int2_full_re.dat");
    ofstream out7;
    out7.open("int2_full_im.dat");

    for (int j = 0; j <= n_plot; j++) {
        t = j * step_plot - a_plot;

        int1_full = integrate_1(z_1);
        int2_full = integrate_2(z_2);

        x2 = x2_0 + K * int1_full * int2_full;

        if (real(x2) < 0.L) {
            x2 = -x2;
        }

        out << t << ", " << setprecision(20) << fabs(x2) << endl;
        out2 << t << ", " << setprecision(20) << real(x2) << endl;
        out3 << t << ", " << setprecision(20) << imag(x2) << endl;
        out4 << t << ", " << setprecision(20) << real(int1_full) << endl;
        out5 << t << ", " << setprecision(20) << imag(int1_full) << endl;
        out6 << t << ", " << setprecision(20) << real(int2_full) << endl;
        out7 << t << ", " << setprecision(20) << imag(int2_full) << endl;
        cout << j << endl;

        cout << t << ", int1 = " << int1_full << endl;
        cout << t << ", int2 = " << int2_full << endl;
        cout << t << ", x^2 = " << x2 << endl;
    }
    out.close();
    out2.close();
    out3.close();
    out4.close();
    out5.close();
    out6.close();
    out7.close();

    return 0;
}