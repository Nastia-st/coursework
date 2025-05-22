#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define T 1.0
#define S 1e-4
#define l0 0.01
#define k0 0.005
#define m 0.01
#define i0 0.1
#define n 0.01
#define w 1000
#define dt 0.01

double i_ob(double t) {
    if (t >= 0 && t <= T/3.0) return i0 * (1 + k0);
    else if (t > T/3.0 && t <= 2*T/3.0) return i0 * (1 + m);
    else return i0 * (1 + n);
}

double l(double t) {
    if (t >= 0 && t <= T/4.0) return l0 * (1 + t / T);
    else if (t > T/4.0 && t <= 3*T/4.0) return l0 * (1 + 2*t / T);
    else return l0 * (1 + 3*t / T);
}

double h(double t) {
    if (t >= 0 && t <= T/2.0) return k0 * (1 + t);
    else return k0 * (1 + t / T);
}

double I(double t) {
    if (t >= 0 && t <= T/3.0)
        return i0 * (1 + h(t));
    else if (t > T/3.0 && t <= 2*T/3.0)
        return i0 * (1 + m) * (1 - h(t));
    else
        return i0 * (1 + n) + h(t) * t * t;
}

double F(double t) {
    double mu0 = 4 * M_PI * 1e-7;
    return 1e6 * (mu0 * pow(w, 2) * S) / (8 * l(t)) * I(t);
}

int main() {
    system("chcp 65001");
    FILE *out = fopen("results.txt", "w");
    if (!out) {
        perror("Не вдалося відкрити файл");
        return 1;
    }

    fprintf(out, "t\tI(t)\tl(t)\th(t)\tF(t)\n");
    for (double t = 0; t <= T; t += dt) {
        double i_t = I(t);
        double l_t = l(t);
        double h_t = h(t);
        double f_t = F(t);
        fprintf(out, "%.2f\t%.6f\t%.6f\t%.6f\t%.6f\n", t, i_t, l_t, h_t, f_t);
    }

    fclose(out);
    printf("Результати збережено у файл 'results.txt'\n");
    return 0;
}
