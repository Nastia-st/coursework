#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592653589793
#define T 100
#define DT 5

typedef struct {
    double Sc;
    double w;
    double i0;
    double k;
    double m;
    double n;
    double l0;
} MagnetParams;

double current(double t, MagnetParams p) {
    double t_frac = t / T;
    if (t_frac >= 0 && t_frac < 1.0/8.0)
        return p.i0 * (1 + p.k * t);
    else if (t_frac >= 1.0/8.0 && t_frac < 3.0/8.0)
        return p.i0 * (1 + p.k);
    else if (t_frac >= 3.0/8.0 && t_frac < 5.0/8.0)
        return p.i0 * (1 + p.k - p.k * (t - 3*T/8.0));
    else if (t_frac >= 5.0/8.0 && t_frac < 7.0/8.0)
        return p.i0 * (1 - p.k);
    else
        return p.i0 * (1 - p.k + p.k * (t - 7*T/8.0));
}
double air_gap(double t, MagnetParams p) {
    if (t <= T / 2.0)
        return p.l0 * (1 + p.m * exp(-t / T));
    else
        return p.l0 * (1 + p.n * exp(-(T - t) / T));
}
double force(double I, double w, double Sc, double lb) {
    return (1e-5 * 0.4 * I * I * w * w * Sc * Sc) / (8 * PI * PI * lb);
}
void simulate(MagnetParams p, const char* filename) {
    FILE *out = fopen(filename, "w");
    if (!out) {
        perror("Помилка відкриття файлу");
        exit(1);
    }
    fprintf(out, "t\tI(t)\tl(t)\tF(t)\n");
    for (double t = 0; t <= T; t += DT) {
        double I_t = current(t, p);
        double l_t = air_gap(t, p);
        double F_t = force(I_t, p.w, p.Sc, l_t);
        fprintf(out, "%.2f\t%.6f\t%.6f\t%.6f\n", t, I_t, l_t, F_t);
    }
    fclose(out);
    printf("Результати збережено у файл: %s\n", filename);
}
int main() {
    system("chcp 65001");
    MagnetParams variant2 = {
        .Sc = 0.0212,
        .w = 2000,
        .i0 = 0.2,
        .k = 0.005,
        .m = 0.01,
        .n = 0.01,
        .l0 = 0.01
    };

    simulate(variant2, "output_variant_2.txt");
    return 0;
}
