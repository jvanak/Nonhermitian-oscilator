#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>

using namespace std;

#define _USE_MATH_DEFINES

int m = 100;
int n = 1000000;
int n_plot = 10000;
double a = 1000.;
double a_plot = 3.;
double a_qplot = 0.8;

double theta = 0.2;
complex<double> Omega = 1.;
complex<double> Omega_q2 = 1.;
complex<double> Omega_q4 = 1.;

// double kp;
double x;
complex<double> Omega2 = 1.;
double step = 2. * a / n;
double Q = 1.;
double q2;
double c = 137.03599;

complex<double> integrate(complex<double> f(double x)) {
    complex<double> result = 0.;
    for (int j = 0; j < n; j++) {
        result += f(-a + (j + 0.5) * step) * step;
        // cout << j + 1 << endl;
    }
    return result;
}
/*complex<double> integrate_2(complex<double> f(double x))
{
    complex<double> result = 0.;
    for (int l = 0; l <= n; l++)
    {
        kp = -a + l * step;
        result += f(kp) * step;

        cout << l + 1 << endl;
    }
    return result;
}*/
complex<double> f(double k) {
    return 1. / (Omega2 - k * k * exp(-2.i * theta));
}
complex<double> g(double k) {
    return (1. / ((Omega - exp(-1.i * theta) * fabs(k)) * (Omega - exp(-1.i * theta) * fabs(k))) - 1. / ((Omega + exp(-1.i * theta) * fabs(k)) * (Omega + exp(-1.i * theta) * fabs(k)))) / fabs(k);
}
/*complex<double> h(double k)
{
    return 1. / (kp * kp - k * k);
}*/
/*complex<double> z(double k) {
    return (1. + exp(-1.i * theta) * fabs(k)) / (exp(2.i * theta) * q2 * q2 / (2. * k * k) + 2. * M_PI * M_PI * (1. - exp(-2.i * theta) * k * k) * (1. - exp(-2.i * theta) * k * k));
    //  return (1. + exp(-1.i * theta) * fabs(k)) / (exp(2.i * theta) * q2 * q2 / (2. * k * k) + 2. * M_PI * M_PI * (1. - exp(-2.i * theta) * k * k + exp(1.i * theta) * q2 * M_PI * 1.i / (2. * M_PI * M_PI * fabs(k))) * (1. - exp(-2.i * theta) * k * k + exp(1.i * theta) * q2 * M_PI * 1.i / (2. * M_PI * M_PI * fabs(k))));
}*/
void solve() {
    for (int j = 0; j < m; j++) {
        Omega2 = 1. + q2 * exp(-1.i * theta) * integrate(f) / (2. * M_PI * M_PI);
    }
}
/*void solve2(double Q2)
{
    for (int j = 0; j < m; j++)
    {
        Omega2 = 1. + Q2 * exp(-1.i * theta) * integrate(f) / (8. * M_PI * M_PI);
    }
}*/
complex<double> z(double k) {
    return (1. + exp(-1.i * theta) * c * fabs(k)) / (q2 * q2 * exp(-2.i * theta) * c * c * k * k / 2. + 2. * M_PI * M_PI * c * c * (1. - exp(-2.i * theta) * c * c * k * k) * (1. - exp(-2.i * theta) * c * c * k * k));
}
complex<double> zp(double k) {
    return (1. + exp(-1.i * theta) * c * fabs(k)) / (c * fabs(k) * (q2 * q2 * exp(-2.i * theta) * c * c * k * k / 2. + 2. * M_PI * M_PI * c * c * (1. - exp(-2.i * theta) * c * c * k * k) * (1. - exp(-2.i * theta) * c * c * k * k)));
} /*
 complex<double> zN_0(double k) {
     return 1 / ((c * fabs(k) * (c * fabs(kp) + c * fabs(k)) * (c * fabs(kp) + c * fabs(k))));
 }
 complex<double> zN(double k) {
     return 1. / ((c * kp * (Omega + exp(-1.i * theta) * c * fabs(kp)) * (Omega + exp(-1.i * theta) * c * fabs(kp))) * (2. * Omega - 1.i * q2 * Omega * Omega / (2. * M_PI * c))) + q2 * c * c * exp(2.i * theta) * integrate(zN_0) / (exp(-2.i * theta) * q2 * q2 * c * c * c * fabs(k) * fabs(k) * fabs(k) + 4. * M_PI * M_PI * c * c * c * k * (1. - exp(-2.i * theta * c * c * k * k)) * (1. - exp(-2.i * theta * c * c * k * k)));
 }*/
int main() {
    double step_plot = 2. * a_plot / n_plot;
    double step_qplot = a_qplot / n_plot;
    /* ofstream out;
      out.open("Omega_re_0.05.dat");
      ofstream out2;
      out2.open("Omega_im_0.05.dat");
      ofstream out3;
      out3.open("error_0.05.dat");

      for (int j = 0; j < n_plot; j++)
      {
          Q = j * step_plot;

          solve(Q);
          Omega = sqrt(Omega2);
          double error = sqrt(norm(Omega2 - 1. - q2 * exp(-i * theta) * integrate(f) / (8. * M_PI * M_PI)));

          out << Q << ", " << setprecision(20) << real(Omega) << endl;
          out2 << Q << ", " << setprecision(20) << imag(Omega) << endl;
          out3 << Q << ", " << error << endl;

          cout << theta << " " << j + 1 << endl;
          // cout << "theta = " << theta << " : Q = " << Q << " : Omega = " << Omega << " : error = " << error << endl;
      }
      out.close();
      out2.close();
      out3.close();*/
    /*
        ofstream out_z;
        out_z.open("Omega_re_0.1_z.dat");
        ofstream out2_z;
        out2_z.open("Omega_im_0.1_z.dat");
        ofstream out3_z;
        out3_z.open("error_0.1_z.dat");

        theta = 0.1L;
        Omega2 = 1.;

        for (int j = 0; j < n_plot; j++)
        {
            double Q2 = j * step_plot;

            solve2(Q2);
            double error2 = norm(Omega2 - 1. - Q2 * exp(-i * theta) * integrate(f) / (8. * M_PI * M_PI));

            out_z << Q2 << ", " << setprecision(20) << real(Omega2) << endl;
            out2_z << Q2 << ", " << setprecision(20) << imag(Omega2) << endl;
            out3_z << Q2 << ", " << error2 << endl;

            cout << theta << " " << j + 1 << endl;
            // cout << "theta = " << theta << " : Q = " << Q << " : Omega = " << Omega << " : error = " << error << endl;
        }
        out_z.close();
        out2_z.close();
    */

    //   ofstream out_0;
    //  out_0.open("Omega.dat");
    //   ofstream out_0_b;
    //   out_0_b.open("Omega_q2.dat");
    //    ofstream out_0_c;
    //   out_0_c.open("Omega_q4.dat");
    ofstream out_x;
    out_x.open("x2.dat");
    ofstream out_p;
    out_p.open("p2.dat");
    ofstream out_x_q2;
    out_x_q2.open("x2_q2.dat");
    //  ofstream out_p_q2;
    // out_x_q2.open("p2_q2.dat");

    // ofstream out2_0;
    //  out2_0.open("Omega_im_old.dat");
    // ofstream out3_0;
    // out3_0.open("Omega_error_old.dat");
    double x2;
    complex<double> p2;
    q2 = Q * Q;
    Omega = (sqrt(16. * M_PI * M_PI * c * c - q2 * q2) - 1.i * q2) / (4. * M_PI * c);

    p2 = (-1. + (1. + 1. / Omega) / (1. - 1.i * q2 * Omega / (4. * M_PI * c)) + q2 * c * c * integrate(zp)) / 2.;
    cout << p2 << endl;

    for (int j = 0; j <= n_plot; j++) {
        Q = j * step_qplot;
        q2 = Q * Q;
        /* if (16. * M_PI * M_PI * c * c > q2 * q2) {
              Omega = (sqrt(16. * M_PI * M_PI * c * c - q2 * q2) - 1.i * q2) / (4. * M_PI * c);
          } else {
              Omega = (1.i * sqrt(q2 * q2 - 16. * M_PI * M_PI * c * c) - 1.i * q2) / (4. * M_PI * c);
          }*/
        // Omega_q2 = 1. - 1.i * q2 / (4. * M_PI * c);
        // Omega_q4 = 1. - 1.i * q2 / (4. * M_PI * c) - q2 * q2 / (32. * M_PI * M_PI * c * c);

        Omega = (sqrt(16. * M_PI * M_PI * c * c - q2 * q2) - 1.i * q2) / (4. * M_PI * c);
        //  x2 = real(0.5 * (-1. + (1. + Omega) / (1. - 1.i * q2 * Omega / (4 * M_PI * c)) + q2 * c * c * exp(-1.i * theta) * integrate(z)));
        //  p2 = real(0.5 * (-1. + (1. + 1. / Omega) / (1. - 1.i * q2 * Omega / (4. * M_PI * c)) + q2 * c * c * integrate(zp)));

        x2 = real(0.5 * (-1. + (1. + Omega) / (1. - 1.i * q2 * Omega / (4 * M_PI * c)) + q2 * c * c * exp(-1.i * theta) * integrate(z)));

        //   solve();
        //  Omega = sqrt(Omega2);
        //  double error = sqrt(norm(Omega * Omega - 1. - q2 * exp(-1.i * theta) * integrate(f) / (2. * M_PI * M_PI)));
        //  double vartheta = atan(q2 / sqrt(16. * M_PI * M_PI * c * c - q2 * q2));
        //   double x2 = real(0.5 * ((1. + exp(-1.i * vartheta)) / (1. - 1.i * q2 * exp(-1.i * vartheta) / (4. * M_PI * c)) + (1.i * q2 - 8 * c * atan(q2 / (sqrt(16. * M_PI * M_PI * c * c - q2 * q2)))) / (sqrt(16. * M_PI * M_PI * c * c - q2 * q2)) + 0.25 * (1.i * q2 * (1. + exp(-1.i * vartheta))) / (M_PI * c * (exp(-3.i * vartheta) - exp(-1.i * vartheta)))));
        //    out_0 << setprecision(25) << real(Omega) << ", " << setprecision(25) << imag(Omega) << endl;
        //  out_0_b << setprecision(25) << real(Omega_q2) << ", " << setprecision(25) << imag(Omega_q2) << endl;
        //   out_0_c << setprecision(25) << real(Omega_q4) << ", " << setprecision(25) << imag(Omega_q4) << endl;
        out_x << setprecision(25) << Q << ", " << x2 << endl;
        //    out_p << setprecision(25) << Q << ", " << p2 << endl;
        out_x_q2 << setprecision(25) << q2 << x2 << endl;
        //  out_p_q2 << setprecision(25) << q2 << p2 << endl;
        // out2_0 << setprecision(25) << imag(Omega) << endl;
        // out3_0 << Q << ", " << error << endl;
        cout << j << endl;
    }
    out_x_q2.close();
    //  out_p_q2.close();
    //  out_0.close();
    // out_0_b.close();
    //  out_0_c.close();
    out_x.close();
    out_p.close();
    //   out2_0.close();
    // out3_0.close();
    /*
        ofstream out_a;
        out_a.open("rho_x_0.1.dat");
        ofstream out2_a;
        out2_a.open("rho_x_0.dat");
        // ofstream out3_a;
        // out3_a.open("error_0.1.dat");

        // theta = 0.1;
        // Omega2 = 1.;
        //  double error;
        Q = 0.1;
        q2 = Q * Q;
        Omega = 1.000002190744742724959337 - 0.0007955191066086031957682589i;

        complex<double> x2 = 0.5 * (-1. + (1. + Omega) / (1. + q2 * integrate(g) / (8. * M_PI * M_PI * Omega)) + exp(-1.i * theta) * q2 * integrate(z));
        cout << "<x^2> (q = 0.1) = " << x2 << endl;
        double x2_re = real(x2);

        for (int j = 0; j <= n_plot; j++)
        {
            x = j * step_plot - a_plot;

            // solve(Q);
            // Omega = sqrt(Omega2);
            //  error = sqrt(norm(Omega2 - 1. - q2 * exp(-1.i * theta) * integrate(f) / (8. * M_PI * M_PI)));

            // out_a << Q << ", " << setprecision(20) << real(Omega) << endl;
            // out2_a << Q << ", " << setprecision(20) << imag(Omega) << endl;
            // out3_a << Q << ", " << error << endl;

            double rho_x = 1. / sqrt(2. * M_PI * x2_re) * exp(-0.5 * x * x / x2_re);
            double rho_x_0 = 1. / sqrt(M_PI) * exp(-x * x);

            out_a << x << ", " << rho_x << endl;
            out2_a << x << ", " << rho_x_0 << endl;
            // cout << "rho(x = " << x << ") = " << rho_x << endl;
        }
        out_a.close();
        out2_a.close();
    /*
        Q = 0.05;
        q2 = Q * Q;
        Omega = 1.000000369660443144681494 - 0.0001988806191521279168391012i;

        x2 = 0.5 * (-1. + (1. + Omega) / (1. + q2 * integrate(g) / (8. * M_PI * M_PI * Omega)) + exp(-1.i * theta) * q2 * integrate(z));
        cout << "<x^2> (q = 0.05) = " << x2 << endl;
        x2_re = real(x2);

        ofstream out_b;
        out_b.open("rho_x_0.05.dat");

        for (int j = 0; j <= n_plot; j++)
        {
            x = j * step_plot - a_plot;

            // solve(Q);
            // Omega = sqrt(Omega2);
            //  error = sqrt(norm(Omega2 - 1. - q2 * exp(-1.i * theta) * integrate(f) / (8. * M_PI * M_PI)));

            // out_a << Q << ", " << setprecision(20) << real(Omega) << endl;
            // out2_a << Q << ", " << setprecision(20) << imag(Omega) << endl;
            // out3_a << Q << ", " << error << endl;

            double rho_x = 1. / sqrt(2. * M_PI * x2_re) * exp(-0.5 * x * x / x2_re);

            out_b << x << ", " << rho_x << endl;
            // cout << "rho(x = " << x << ") = " << rho_x << endl;
        }
        out_b.close();
    */
    //  out3_a.close();
    /*
        ofstream out_b;
        out_b.open("Omega_re_0.15.dat");
        ofstream out2_b;
        out2_b.open("Omega_im_0.15.dat");
        ofstream out3_b;
        out3_b.open("error_0.15.dat");

        theta = 0.15L;
        Omega2 = 1.;

        for (int j = 0; j < n_plot; j++)
        {
            double Q = j * step_plot;

            solve(Q);
            complex<double> Omega = sqrt(Omega2);
            double error = sqrt(norm(Omega2 - 1. - q2 * exp(-i * theta) * integrate(f) / (8. * M_PI * M_PI)));

            out_b << Q << ", " << setprecision(20) << real(Omega) << endl;
            out2_b << Q << ", " << setprecision(20) << imag(Omega) << endl;
            out3_b << Q << ", " << error << endl;

            cout << theta << " " << j + 1 << endl;
            // cout << "theta = " << theta << " : Q = " << Q << " : Omega = " << Omega << " : error = " << error << endl;
        }
        out_b.close();
        out2_b.close();
        out3_b.close();

        ofstream out_c;
        out_c.open("Omega_re_0.2.dat");
        ofstream out2_c;
        out2_c.open("Omega_im_0.2.dat");
        ofstream out3_c;
        out3_c.open("error_0.2.dat");

        theta = 0.2L;
        Omega2 = 1.;

        for (int j = 0; j < n_plot; j++)
        {
            double Q = j * step_plot;

            solve(Q);
            complex<double> Omega = sqrt(Omega2);
            double error = sqrt(norm(Omega2 - 1. - q2 * exp(-i * theta) * integrate(f) / (8. * M_PI * M_PI)));

            out_c << Q << ", " << setprecision(20) << real(Omega) << endl;
            out2_c << Q << ", " << setprecision(20) << imag(Omega) << endl;
            out3_c << Q << ", " << error << endl;

            cout << theta << " " << j + 1 << endl;
            // cout << "theta = " << theta << " : Q = " << Q << " : Omega = " << Omega << " : error = " << error << endl;
        }
        out_c.close();
        out2_c.close();
        out3_c.close();

        ofstream out_d;
        out_d.open("Omega_re_0.3.dat");
        ofstream out2_d;
        out2_d.open("Omega_im_0.3.dat");
        ofstream out3_d;
        out3_d.open("error_0.3.dat");

        theta = 0.3L;
        Omega2 = 1.;

        for (int j = 0; j < n_plot; j++)
        {
            double Q = j * step_plot;

            solve(Q);
            complex<double> Omega = sqrt(Omega2);
            double error = sqrt(norm(Omega2 - 1. - q2 * exp(-i * theta) * integrate(f) / (8. * M_PI * M_PI)));

            out_d << Q << ", " << setprecision(20) << real(Omega) << endl;
            out2_d << Q << ", " << setprecision(20) << imag(Omega) << endl;
            out3_d << Q << ", " << error << endl;

            cout << theta << " " << j + 1 << endl;
            // cout << "theta = " << theta << " : Q = " << Q << " : Omega = " << Omega << " : error = " << error << endl;
        }
        out_d.close();
        out2_d.close();
        out3_d.close();

        ofstream out_e;
        out_e.open("Omega_re_0.4.dat");
        ofstream out2_e;
        out2_e.open("Omega_im_0.4.dat");
        ofstream out3_e;
        out3_e.open("error._0.4.dat");

        theta = 0.4L;
        Omega2 = 1.;

        for (int j = 0; j < n_plot; j++)
        {
            double Q = j * step_plot;

            solve(Q);
            complex<double> Omega = sqrt(Omega2);
            double error = sqrt(norm(Omega2 - 1. - q2 * exp(-i * theta) * integrate(f) / (8. * M_PI * M_PI)));

            out_e << Q << ", " << setprecision(20) << real(Omega) << endl;
            out2_e << Q << ", " << setprecision(20) << imag(Omega) << endl;
            out3_e << Q << ", " << error << endl;

            cout << theta << " " << j + 1 << endl;
            // cout << "theta = " << theta << " : Q = " << Q << " : Omega = " << Omega << " : error = " << error << endl;
        }
        out_e.close();
        out2_e.close();
        out3_e.close();

        return 0;*/
}