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
ldouble c = 137.03599L;
ldouble mu = fr / c;

int n_0 = 20000000;
int n = 300000;
ldouble a_0 = 0.5L;
ldouble a = 0.00015L;
int n_plot = 5000;
ldouble a_plot = 10000.L;
ldouble sw = 0.0001L;

ldouble t = -a_plot;
ldouble q2 = q * q;

lcomplex integrate_0(lcomplex f(ldouble x)) {
    ldouble step = 2.L * a_0 / n_0;
    lcomplex result = 0.;

    for (int l = 0; l <= n_0; l++) {
        result += f(l * step) * step;
    }
    return 2.L * result;
}
lcomplex integrate(lcomplex f(ldouble x)) {
    ldouble zn = 1.L;
    ldouble step = 2.L * a / n;
    lcomplex result = 0.L;

    for (int l = 0; l <= n; l++) {
        if (abs(f(-a + mu + l * step) - f(-a + mu + (l - 1.L) * step)) > sw && abs(f(-a + mu + l * step) - f(-a + mu + (l - 1.L) * step)) > abs(f(-a + mu + l * step) + f(-a + mu + (l - 1.L) * step))) {
            zn = -zn;
        }
        result += zn * f(-a + mu + l * step) * step;
    }
    return result;
}
lcomplex z_0(ldouble k) {
    return k / (q2 * q2 * exp(-2.il * theta) * c * c * k * k + 4.L * M_PI * M_PI * c * c * (om * om - exp(-2.il * theta) * c * c * k * k) * (om * om - exp(-2.il * theta) * c * c * k * k));
}
lcomplex z_a(ldouble k) {
    return exp(-3.il * theta / 2.L - (k * exp(-1.il * theta) - mu) * (k * exp(-1.il * theta) - mu) / (2.L * sigma * sigma)) * sin(exp(-1.il * theta) * c * abs(k) * t) * sqrt(abs(k) / (q2 * q2 * exp(-2.il * theta) * c * c * k * k + 4.L * M_PI * M_PI * c * c * (om * om - exp(-2.il * theta) * c * c * k * k) * (om * om - exp(-2.il * theta) * c * c * k * k)));
}
int main() {

    ldouble kstep_plot = 2.L * a / n;
    t = 10000.L;
    lcomplex int1;
    ldouble znam = 1.L;

    ofstream kout;
    kout.open("int1_a.dat");
    ofstream kout2;
    kout2.open("int1_a_re.dat");
    ofstream kout3;
    kout3.open("int1_a_im.dat");

    for (int j = 0; j <= n; j++) {
        ldouble k = mu + j * kstep_plot - a;

        if (abs(real(z_a(k)) - real(z_a(k - kstep_plot))) > sw && abs(z_a(k) - z_a(k - kstep_plot)) > abs(z_a(k) + z_a(k - kstep_plot))) {
            znam = -znam;
        }
        int1 = znam * z_a(k);

        kout << k << ", " << setprecision(20) << fabs(int1) << endl;
        kout2 << k << ", " << setprecision(20) << real(int1) << endl;
        kout3 << k << ", " << setprecision(20) << imag(int1) << endl;
    }
    kout.close();
    kout2.close();
    kout3.close();

    cout << "int plot done" << endl;

    lcomplex Omega = sqrt(om * om - q2 * q2 / (16.L * M_PI * M_PI * c * c)) - 1il * q2 / (4.L * M_PI * c);
    lcomplex x2_0 = Omega / (2.L * om * om - 1.il * q2 * Omega / (2.L * M_PI * c)) + q2 * c * c * c * exp(-2.il * theta) * integrate_0(z_0);
    lcomplex K_a = -2.L * q * c * sqrt(c) / (sqrt(sqrt(M_PI)) * sqrt(sigma));
    ldouble step_plot = 2.L * a_plot / n_plot;
    lcomplex x_a;
    lcomplex x2_a;

    cout << "arg(Omega) = " << arg(Omega) << endl;
    cout << "mu = " << mu << endl;
    cout << "K_a = " << K_a << endl;
    cout << "x^2_0 = " << x2_0 << endl;

    ofstream out;
    out.open("x_a_t.dat");
    ofstream out2;
    out2.open("x_a_t_re.dat");
    ofstream out3;
    out3.open("x_a_t_im.dat");
    ofstream out4;
    out4.open("x2_a_t.dat");
    ofstream out5;
    out5.open("x2_a_t_re.dat");
    ofstream out6;
    out6.open("x2_a_t_im.dat");

    for (int j = 0; j <= n_plot; j++) {
        t = j * step_plot - a_plot;

        x_a = K_a * integrate(z_a);
        x2_a = x2_0 + x_a * x_a;

        out << t << ", " << setprecision(20) << fabs(x_a) << endl;
        out2 << t << ", " << setprecision(20) << real(x_a) << endl;
        out3 << t << ", " << setprecision(20) << imag(x_a) << endl;
        out4 << t << ", " << setprecision(20) << fabs(x2_a) << endl;
        out5 << t << ", " << setprecision(20) << real(x2_a) << endl;
        out6 << t << ", " << setprecision(20) << imag(x2_a) << endl;
        cout << j << endl;

        cout << t << ", x_a = " << x_a << endl;
        cout << t << ", x^2_a = " << x2_a << endl;
    }
    out.close();
    out2.close();
    out3.close();
    out4.close();
    out5.close();
    out6.close();

    return 0;
}