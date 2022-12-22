#include <stdio.h>
int main(void) {
	double a[3][3] = { {4,12,-10},{12,37,-43},{-10,-43,98} };

    //All intermediate matrix should have the same number of row&column as A
	double c[sizeof(a) / sizeof(a[0])][sizeof(a) / sizeof(a[0])] = { 0 };
    double d[sizeof(a) / sizeof(a[0])][sizeof(a) / sizeof(a[0])] = { 0 };
    double l[sizeof(a) / sizeof(a[0])][sizeof(a) / sizeof(a[0])] = { 0 };


    //initialize L, set the diagonal entry in the diagonal to 1
    for (int i = 0;i < sizeof(a) / sizeof(a[0]);i++) {
        l[i][i] = 1;
    }
    for (int j = 0;j < sizeof(a) / sizeof(a[0]);j++) {
        double sum = 0.0;
        for (int s = 0;s < j;s++) {
            sum += (d[s][s] * l[j][s]*l[j][s]);
        }
        c[j][j] = a[j][j] - sum;
        d[j][j] = c[j][j];
        for (int i = j + 1;i < sizeof(a) / sizeof(a[0]);i++) {
            sum = 0.0;
            for (int s = 0;s < j;s++) {
                sum += (d[s][s] * l[i][s] * l[j][s]);
            }
            c[i][j] = a[i][j] - sum;
            l[i][j] = c[i][j] / d[j][j];
        }
    }
    for (int i = 0;i < sizeof(l) / sizeof(l[0]);i++) {
        for (int j = 0;j < sizeof(l) / sizeof(l[0]);j++) {
        printf("%f ", l[i][j]);
        }
        printf("\n");
	}
    printf("\n");
    for (int i = 0;i < sizeof(d) / sizeof(d[0]);i++) {
        for (int j = 0;j < sizeof(d) / sizeof(d[0]);j++) {
            printf("%f ", d[i][j]);
        }
        printf("\n");
    }
}