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
ldouble c = 137.03599L;
ldouble mu = 0.L;

int n = 2100000001;
ldouble a_x = 14.L;
ldouble a_p = 1.L;
ldouble a_0 = 4.L;
ldouble a_N = 4.L;
ldouble kp = 5.L;

ldouble q2 = q * q;
lcomplex Omega = sqrt(om * om - q2 * q2 / (16.L * M_PI * M_PI * c * c)) - 1il * q2 / (4.L * M_PI * c);

ldouble integrate_0(ldouble f(ldouble x)) {
    ldouble step = a_0 / n;
    ldouble result = 0.L;

    for (int l = 1; l <= n; l++) {
        result += f(l * step) * step;
    }
    return 2.L * result;
}
lcomplex integrate_x(lcomplex f(ldouble x)) {
    ldouble step = a_x / n;
    lcomplex result = 0.L;

    for (int l = 0; l <= n; l++) {
        result += f(l * step) * step;
    }
    return 2.L * result;
}
lcomplex integrate_p(lcomplex f(ldouble x)) {
    ldouble step = a_p / n;
    lcomplex result = 0.L;

    for (int l = 1; l <= n; l++) {
        result += f(l * step) * step;
    }
    return 2.L * result;
}
lcomplex integrate_N(lcomplex f(ldouble x)) {
    ldouble step = 2.L * a_N / n;
    lcomplex result = 0.L;

    for (int l = 1; l <= n; l++) {
        kp = l * step;
        result += f(kp) * step;

        cout << l << endl;
    }
    return 2.L * result;
}
ldouble omk(ldouble k) {
    return sqrt(c * c * k * k + mu * mu);
}
ldouble z_0(ldouble k) {
    return 1.L / (omk(k) * (omk(kp) + omk(k)) * (omk(kp) + omk(k)));
}
lcomplex z_N(ldouble k) {
    return Omega * Omega * Omega / (om * om * omk(k) * (Omega + exp(-1.il * theta) * omk(k)) * (Omega + exp(-1.il * theta) * omk(k)) * (2.L - 1.il * q2 * Omega / (2.L * M_PI * c * om * om))) + q2 * exp(-2.il * theta) * omk(k) * omk(k) * omk(k) * integrate_0(z_0) / (q2 * q2 * exp(-2.il * theta) * omk(k) * omk(k) / (c * c) + 4.L * M_PI * M_PI * (om * om - exp(-2.il * theta) * omk(k) * omk(k)) * (om * om - exp(-2.il * theta) * omk(k) * omk(k)));
}
lcomplex z_x(ldouble k) {
    return omk(k) / (q2 * q2 * exp(-2.il * theta) * omk(k) * omk(k) / (c * c) + 4.L * M_PI * M_PI * (om * om - exp(-2.il * theta) * omk(k) * omk(k)) * (om * om - exp(-2.il * theta) * omk(k) * omk(k)));
}
lcomplex z_p(ldouble k) {
    return 1.L / (omk(k) * (q2 * q2 * exp(-2.il * theta) * omk(k) * omk(k) / (c * c) + 4.L * M_PI * M_PI * (om * om - exp(-2.il * theta) * omk(k) * omk(k)) * (om * om - exp(-2.il * theta) * omk(k) * omk(k))));
}
int main() {
    /*   int n_plot = 10000;
       ldouble a_plot = 0.01;
       ldouble step_plot = a_plot / n_plot;

       ofstream out;
       out.open("z_N.dat");

       for (int j = 0; j <= n_plot; j++) {
           ldouble k = j * step_plot;

           // cout << k << ", " << setprecision(20) << abs(z_N(j * step_plot)) << endl;
           out << k << ", " << setprecision(20) << abs(z_p(j * step_plot)) << endl;
       }
       out.close();
       cout << "done" << endl;*/

    lcomplex N_st = q2 * z_N(kp) / (4.L * M_PI * M_PI);
    cout << "<aa^+> = " << N_st << endl;

    lcomplex x2 = Omega / (2.L * om * om - 1.il * q2 * Omega / (2.L * M_PI * c)) + q2 * exp(-2.il * theta) * integrate_x(z_x);
    cout << "<x^2> = " << x2 << endl;

    lcomplex p2 = om * om / (Omega * (2.L - 1.il * q2 * Omega / (2.L * M_PI * c * om * om))) + q2 * om * om * om * om * integrate_p(z_p);
    cout << "<p^2> = " << p2 << endl;

    return 0;
}