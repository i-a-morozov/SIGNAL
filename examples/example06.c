// gfortran -c -cpp -std=f2018 -Wall -pedantic -O3 -ffast-math -march=native -Wno-unused-function signal.f90
// ar rcs libsignal.a signal.o
// gcc -o example06 -std=c99 -Wall -pedantic -O3 -ffast-math -march=native -L. -I. example06.c -lsignal  -llapack -lblas -lm -lfftw3 -lgsl -lgslcblas -lgfortran

// EXAMPLE-06: C INTERFACE

#include <stdio.h>
#include <math.h>

#include "signal.h"

// SIGNAL FUNCTIONS WITH C INTERFACE

// void    convert_real_(int*, double*, double*) ;
// void    convert_complex_(int*, double*, double*, double*) ;
// void    fft_radix_two_(int*, int*, double*) ;
// void    fft_radix_eight_(int*, int*, double*) ;
// void    fft_external_(int*, int*, double*) ;
// void    ffrft_(int*, double*, double*) ;
// void    window_(int*, int*, double*) ;
// int     peak_(int*, double*, int*) ;
// double  frequency_(int*, int*, int*, double*, double*, double*) ;
// void    decomposition_(int*, int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
// void    frequency_list_(int*, int*, int*, double*, double*, double*, int*, double*) ;
// void    amplitude_list_(int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
// void    fit_(int*, double*, int*, double*, double*, double*, double*, double*) ;

int main(){

    // define test sample parameters (five harmonics)
    double f1 = 1.0*0.123456789, a1=1.0E+0, b1 = 1.0E-1 ;
    double f2 = 2.0*0.123456789, a2=1.0E-1, b2 = 5.0E-2 ;
    double f3 = 3.0*0.123456789, a3=1.0E-3, b3 = 0.0E+0 ;
    double f4 = 4.0*0.123456789, a4=0.0E+0, b4 = 1.0E-4 ;
    double f5 = 5.0*0.123456789, a5=1.0E-6, b5 = 1.0E-9 ;

    // set test sample length, real and imaginary parts
    double pi = 2.0*acos(0.0) ;
    int length = 2048 ;
    double sample_r[length] ;
    double sample_i[length] ;
    for(int i=1;i<=length;i++){
        sample_r[i-1]  = a1*cos(2.0*pi*f1*i) + b1*sin(2.0*pi*f1*i) +
                         a2*cos(2.0*pi*f2*i) + b2*sin(2.0*pi*f2*i) +
                         a3*cos(2.0*pi*f3*i) + b3*sin(2.0*pi*f3*i) +
                         a4*cos(2.0*pi*f4*i) + b4*sin(2.0*pi*f4*i) +
                         a5*cos(2.0*pi*f5*i) + b5*sin(2.0*pi*f5*i) ;
        sample_i[i-1]  = b1*cos(2.0*pi*f1*i) - a1*sin(2.0*pi*f1*i) +
                         b2*cos(2.0*pi*f2*i) - a2*sin(2.0*pi*f2*i) +
                         b3*cos(2.0*pi*f3*i) - a3*sin(2.0*pi*f3*i) +
                         b4*cos(2.0*pi*f4*i) - a4*sin(2.0*pi*f4*i) +
                         b5*cos(2.0*pi*f5*i) - a5*sin(2.0*pi*f5*i) ;
    }

    // print input sample parameters
    printf("input \n") ;
    printf("%.15f %.15f %.15f\n", f1, a1, b1) ;
    printf("%.15f %.15f %.15f\n", f2, a2, b2) ;
    printf("%.15f %.15f %.15f\n", f3, a3, b3) ;
    printf("%.15f %.15f %.15f\n", f4, a4, b4) ;
    printf("%.15f %.15f %.15f\n", f5, a5, b5) ;
    printf("\n") ;

    // set complex flag and format test sample
    int flag = 1 ;
    double sample[2*length] ;
    if (flag == 0) {
        convert_real_(&length, sample_r, sample) ;
    } else {
        convert_complex_(&length, sample_r, sample_i, sample) ;
    }

    // set window
    int order = 4 ;
    double window[length] ;
    window_(&length, &order, window) ;
    double total = 0.0 ;
    for (int i = 0; i < length; i += 1) {
        total = total + window[i] ;
    }

    int peak ;
    double frequency ;

    // frequency_ (bin)
    peak = 0 ;
    frequency = frequency_(&flag, &peak, &length, &total, window, sample) ;
    printf("frequency_\n") ;
    printf("%.15f\n", frequency) ;
    printf("\n") ;

    // frequency_ (peak)
    printf("frequency_\n") ;
    for (int i = 1; i <= 5 ; i++)
    {
        peak = i ;
        frequency = frequency_(&flag, &peak, &length, &total, window, sample) ;
        printf("%.15f\n", frequency) ;
    }
    printf("\n") ;

    int method ;
    int loop = 5 ;
    double fre_amp[loop], cos_amp[loop], sin_amp[loop] ;

    // decomposition_ (subtract)
    method = 0 ;
    decomposition_(&flag, &method, &length, &total, window, sample, &loop, fre_amp, cos_amp, sin_amp) ;
    printf("decomposition_\n") ;
    for(int i = 0; i < loop; i++){
        printf("%.15f %.15f %.15f\n", fre_amp[i], cos_amp[i], sin_amp[i]) ;
        fre_amp[i] = 0.0 ;
        cos_amp[i] = 0.0 ;
        sin_amp[i] = 0.0 ;
    }
    printf("\n") ;

    // decomposition_ (peak)
    method = 1 ;
    decomposition_(&flag, &method, &length, &total, window, sample, &loop, fre_amp, cos_amp, sin_amp) ;
    printf("decomposition_\n") ;
    for(int i = 0; i < loop; i++){
        printf("%.15f %.15f %.15f\n", fre_amp[i], cos_amp[i], sin_amp[i]) ;
        fre_amp[i] = 0.0 ;
        cos_amp[i] = 0.0 ;
        sin_amp[i] = 0.0 ;
    }
    printf("\n") ;

    // frequency_list_ (peak)
    frequency_list_(&flag, &method, &length, &total, window, sample, &loop, fre_amp) ;
    printf("frequency_list_\n") ;
    for(int i = 0; i < loop; i++){
        printf("%.15f\n", fre_amp[i]) ;
    }
    printf("\n") ;

    // amplitude_list_
    amplitude_list_(&flag, &length, &total, window, sample, &loop, fre_amp, cos_amp, sin_amp) ;
    printf("amplitude_list_\n") ;
    for(int i = 0; i < loop; i++){
        printf("%.15f %.15f %.15f\n", fre_amp[i], cos_amp[i], sin_amp[i]) ;
    }
    printf("\n") ;

    // fit_ (least squares)
    printf("fit_ \n") ;
    for (int i = 1; i <= loop ; i++)
    {
        peak = i ;
        frequency = frequency_(&flag, &peak, &length, &total, window, sample) ;
        fre_amp[i-1] = frequency ;
    }
    double mean ;
    double error ;
    fit_(&length, sample_r, &loop, fre_amp, &mean, cos_amp, sin_amp, &error) ;
    printf("%.15f\n", error) ;
    printf("%.15f\n", mean) ;
    for (int i = 0; i < loop ; i++)
    {
        printf("%.15f %.15f %.15f\n", fre_amp[i], cos_amp[i], sin_amp[i]) ;
    }

    return 0 ;
}