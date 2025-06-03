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
ldouble mu = 10.L;

int n = 21000001;
ldouble a_0 = 1.L;
ldouble a_1 = 1300.L;
ldouble kp = 0.L;

ldouble q2 = q * q;

ldouble omk(ldouble k) {
    return sqrt(c * c * k * k + mu * mu);
}
lcomplex integrate_0(lcomplex f(ldouble x)) {
    ldouble step = 2.L * a_0 / n;
    lcomplex result = 0.L;

    for (int l = 0; l <= n; l++) {
        result += f(l * step) * step;
    }
    return 2.L * result;
}
lcomplex integrate_1(lcomplex f(ldouble x)) {
    ldouble step = 2.L * a_1 / n;
    lcomplex result = 0.L;

    for (int l = 0; l <= n; l++) {
        if (abs(kp - l * step) >= step) {
            result += f(l * step) * step;
        } else if (kp < l * step) {
            result += f(kp + step) * step;
        } else {
            result += f(kp - step) * step;
        }
    }
    return 2.L * result;
}
lcomplex z_0(ldouble k) {
    return 1.L / (q2 * q2 * exp(-2.il * theta) * omk(k) * omk(k) / 2.L + 2.L * M_PI * M_PI * c * c * (om * om - exp(-2.il * theta) * omk(k) * omk(k)) * (om * om - exp(-2.il * theta) * omk(k) * omk(k)));
}
lcomplex z_1(ldouble k) {
    return omk(k) * omk(k) * (1.L / (kp - k) + 1.L / (kp + k)) / (q2 * q2 * exp(-2.il * theta) * omk(k) * omk(k) / (c * c) + 4.L * M_PI * M_PI * (om * om - exp(-2.il * theta) * omk(k) * omk(k)) * (om * om - exp(-2.il * theta) * omk(k) * omk(k)));
}
int main() {
    lcomplex Omega = sqrt(om * om - q2 * q2 / (16.L * M_PI * M_PI * c * c)) - 1il * q2 / (4.L * M_PI * c);

    lcomplex K_0 = 1.L / (1.L - 1.il * q2 * Omega / (4.L * M_PI * c * om * om)) + exp(-1.il * theta) * c * c * om * om * q2 * integrate_0(z_0);
    cout << "(1, 0) = " << setprecision(10) << K_0 << endl;

    lcomplex K_1 = Omega * Omega * omk(kp) / ((Omega * Omega - exp(-2.il * theta) * omk(kp) * omk(kp)) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + omk(kp) * (exp(-2.il * theta) * omk(kp) * omk(kp) - om * om) / (q2 * q2 * exp(-2.il * theta) * omk(kp) * omk(kp) / (c * c * 4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * omk(kp) * omk(kp)) * (om * om - exp(-2.il * theta) * omk(kp) * omk(kp))) + q2 * exp(-1.il * theta) * integrate_1(z_1) / c;
    cout << endl
         << "k = 0:    (0, 0) = " << setprecision(10) << K_1 << endl;

    kp = a_1 / 2000.L;
    K_1 = Omega * Omega * omk(kp) / ((Omega * Omega - exp(-2.il * theta) * omk(kp) * omk(kp)) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + omk(kp) * (exp(-2.il * theta) * omk(kp) * omk(kp) - om * om) / (q2 * q2 * exp(-2.il * theta) * omk(kp) * omk(kp) / (c * c * 4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * omk(kp) * omk(kp)) * (om * om - exp(-2.il * theta) * omk(kp) * omk(kp))) + q2 * exp(-1.il * theta) * integrate_1(z_1) / c;
    cout << "k = 1:    (0, 0) = " << scientific << setprecision(10) << K_1 << endl;

    kp = a_1 / 2.L;
    K_1 = Omega * Omega * omk(kp) / ((Omega * Omega - exp(-2.il * theta) * omk(kp) * omk(kp)) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + omk(kp) * (exp(-2.il * theta) * omk(kp) * omk(kp) - om * om) / (q2 * q2 * exp(-2.il * theta) * omk(kp) * omk(kp) / (c * c * 4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * omk(kp) * omk(kp)) * (om * om - exp(-2.il * theta) * omk(kp) * omk(kp))) + q2 * exp(-1.il * theta) * integrate_1(z_1) / c;
    cout << "k = 500:  (0, 0) = " << scientific << setprecision(10) << K_1 << endl;

    kp = a_1;
    K_1 = Omega * Omega * omk(kp) / ((Omega * Omega - exp(-2.il * theta) * omk(kp) * omk(kp)) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + omk(kp) * (exp(-2.il * theta) * omk(kp) * omk(kp) - om * om) / (q2 * q2 * exp(-2.il * theta) * omk(kp) * omk(kp) / (c * c * 4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * omk(kp) * omk(kp)) * (om * om - exp(-2.il * theta) * omk(kp) * omk(kp))) + q2 * exp(-1.il * theta) * integrate_1(z_1) / c;
    cout << "k = 1000: (0, 0) = " << scientific << setprecision(10) << K_1 << endl;

    return 0;
}