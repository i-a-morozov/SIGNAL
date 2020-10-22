// EXAMPLE-06: C INTERFACE

#include <stdio.h>
#include <math.h>

#include "signal.h"

// SIGNAL FUNCTIONS WITH EXPLICIT C INTERFACE

// void    generate_signal_(int* flag, int* length, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
// void    ffrft_(int* length, double* argument, double* sequence) ;
// void    fft_external_(int* length, int* direction, double* sequence) ;
// void    fft_radix_two_(int* length, int* direction, double* sequence) ;
// void    fft_radix_eight_(int* length, int* direction, double* sequence) ;
// void    compute_table_(int* length, int* pad) ;
// void    destroy_table_() ;
// void    convert_real_(int* length, double* r_part, double* sequence) ;
// void    convert_complex_(int* length, double* r_part, double* i_part, double* sequence) ;
// int     round_up_(int* number) ;
// void    pad_(int* linput, int* loutput, double* input, double* output) ;
// void    remove_mean_(int* length, double* input, double* output) ;
// void    remove_window_mean_(int* length, double* total, double* window, double* input, double* output) ;
// void    apply_window_(int* length, double* window, double* input, double* output) ;
// void    filter_(int* length, double* sequence, int* limit, double* svd_list) ;
// int     peak_(int* length, double* sequence, int* id) ;
// void    window_cos_(int* length, int* order, double* window) ;
// void    window_cos_generic_(int* length, double* order, double* window) ;
// void    window_kaiser_(int* length, double* order, double* window) ;
// double  frequency_initial_(double* range_min, double* range_max, int* peak, int* length, int* pad, double* sequence) ;
// double  frequency_initial__(double* range_min, double* range_max, int* peak, int* length, int* pad, double* sequence) ;
// double  frequency_refine_(int* method, int* length, double* sequence, double* initial) ;
// double  frequency_refine__(int* method, int* length, double* sequence, double* initial) ;
// double  binary_amplitude_(int* flag, int* length, double* total, double* window, double* sequence, double* initial) ;
// double  golden_amplitude_(int* flag, int* length, double* total, double* window, double* sequence, double* initial) ;
// double  frequency_(int* flag, double* range_min, double* range_max, int* peak, int* method, int* length, int* pad, double* total, double* window, double* sequence) ;
// double  frequency__(int* flag, double* range_min, double* range_max, int* peak, int* method, int* length, int* pad, double* total, double* window, double* sequence) ;
// void    amplitude_(int* flag, int* length, double* total, double* window, double* sequence, double* frequency, double* cos_amp, double* sin_amp, double* amp) ;
// void    decomposition_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
// void    decomposition__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
// void    frequency_list_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency) ;
// void    frequency_list__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency) ;
// void    amplitude_list_(int* flag, int* length, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
// void    frequency_correction_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
// void    frequency_correction__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
// void    fit_(int* length, double* sequence, int* loop, double* frequency, double* mean, double* cos_amp, double* sin_amp, double* error) ;
// void    fit_parabola_(int* length, double* x, double* y, double* a, double* b, double* c, double* maximum) ;

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

    // set frequency range
    double range_min = 0.00 ;
    double range_max = 0.99 ;

    // set window
    int order = 4 ;
    double window[length] ;
    window_cos_(&length, &order, window) ;
    double total = 0.0 ;
    for (int i = 0; i < length; i += 1) {
        total = total + window[i] ;
    }

    int peak ;
    int method = 2 ;
    double frequency ;

    // frequency_ (bin)
    peak = 0 ;
    frequency = frequency_(&flag, &range_min, &range_max, &peak, &method, &length, &length, &total, window, sample) ;
    printf("frequency_\n") ;
    printf("%.15f\n", frequency) ;
    printf("\n") ;

    // frequency_ (peak)
    printf("frequency_\n") ;
    for (int i = 1; i <= 5 ; i++)
    {
        peak = i ;
        frequency = frequency_(&flag, &range_min, &range_max, &peak, &method, &length, &length, &total, window, sample) ;
        printf("%.15f\n", frequency) ;
    }
    printf("\n") ;

    int mode ;
    int loop = 5 ;
    double fre_amp[loop], cos_amp[loop], sin_amp[loop] ;

    // decomposition_ (subtract)
    mode = 0 ;
    decomposition_(&flag, &range_min, &range_max, &method, &mode, &length, &length, &total, window, sample, &loop, fre_amp, cos_amp, sin_amp) ;
    printf("decomposition_\n") ;
    for(int i = 0; i < loop; i++){
        printf("%.15f %.15f %.15f\n", fre_amp[i], cos_amp[i], sin_amp[i]) ;
        fre_amp[i] = 0.0 ;
        cos_amp[i] = 0.0 ;
        sin_amp[i] = 0.0 ;
    }
    printf("\n") ;

    // frequency_list_ (subtract)
    frequency_list_(&flag, &range_min, &range_max, &method, &mode, &length, &length, &total, window, sample, &loop, fre_amp) ;
    printf("frequency_list_\n") ;
    for(int i = 0; i < loop; i++){
        printf("%.15f\n", fre_amp[i]) ;
    }
    printf("\n") ;

    // decomposition_ (peak)
    mode = 1 ;
    decomposition_(&flag, &range_min, &range_max, &method, &mode, &length, &length, &total, window, sample, &loop, fre_amp, cos_amp, sin_amp) ;
    printf("decomposition_\n") ;
    for(int i = 0; i < loop; i++){
        printf("%.15f %.15f %.15f\n", fre_amp[i], cos_amp[i], sin_amp[i]) ;
        fre_amp[i] = 0.0 ;
        cos_amp[i] = 0.0 ;
        sin_amp[i] = 0.0 ;
    }
    printf("\n") ;

    // frequency_list_ (peak)
    frequency_list_(&flag, &range_min, &range_max, &method, &mode, &length, &length, &total, window, sample, &loop, fre_amp) ;
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
        cos_amp[i] = 0.0 ;
        sin_amp[i] = 0.0 ;
    }
    printf("\n") ;

    // fit_ (least squares)
    printf("fit_ \n") ;
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