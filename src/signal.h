#ifndef _SIGNAL__H__
#define _SIGNAL__H__
void    generate_signal_(int*, int*, double*, int*, double*, double*, double*, double*) ;
void    ffrft_(int*, double*, double*) ;
void    fft_external_(int*, int*, double*) ;
void    fft_radix_two_(int*, int*, double*) ;
void    fft_radix_eight_(int*, int*, double*) ;
void    compute_table_(int*) ;
void    destroy_table_() ;
void    convert_real_(int*, double*, double*) ;
void    convert_complex_(int*, double*, double*, double*) ;
int     round_up_(int*) ;
void    pad_(int*, int*, double*, double*) ;
void    remove_mean_(int*, double*, double*) ;
void    remove_window_mean_(int*, double*, double*, double*, double*) ;
void    apply_window_(int*, double*, double*, double*) ;
void    filter_(int*, double*, int*) ;
int     peak_(int*, double*, int*) ;
void    window_cos_(int*, int*, double*) ;
void    window_cos_generic_(int*, double*, double*) ;
void    window_kaiser_(int*, double*, double*) ;
double  frequency_(int*, int*, int*, int*, double*) ;
double  frequency__(int*, int*, int*, int*, double*) ;
void    decomposition_(int*, int*, int*, int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
void    frequency_list_(int*, int*, int*, int*, int*, double*, double*, double*, int*, double*) ;
void    amplitude_(int*, int*, double*, double*, double*, double*, double*, double*, double*) ;
void    amplitude_list_(int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
void    decomposition__(int*, int*, int*, int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
void    frequency_list__(int*, int*, int*, int*, int*, double*, double*, double*, int*, double*) ;
void    frequency_correction_(int*, int*, int*, int*, int*, double*, double*, int*, double*, double*, double*) ;
void    fit_(int*, double*, int*, double*, double*, double*, double*, double*) ;
#endif