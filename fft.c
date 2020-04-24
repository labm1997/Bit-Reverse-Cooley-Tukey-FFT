#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#define N 1024
#define log2N 10
#define hlength 2047
#define PI 3.14159265358979323846264338327950288
#define n_bench_loops 100000

typedef struct {
    double real;
    double imag;
} complex;

#define cmr(x,y) (((x).real * (y).real - (x).imag * (y).imag))
#define cmi(x,y) (((x).real * (y).imag + (x).imag * (y).real))
#define cmul(x,y,u) (u).real = cmr((x),(y)); (u).imag = cmi((x),(y));
#define csum(x,y,u) (u).real = (x).real + (y).real; (u).imag = (x).imag + (y).imag;
#define csub(x,y,u) (u).real = (x).real - (y).real; (u).imag = (x).imag - (y).imag;

void init_samples(complex *samples, double sampleWindow){
    for(uint16_t i=0; i<N ; i++){
        double t = ((double)i)*sampleWindow/N;
        samples[i].real = cos(2*PI*t*1)+2*cos(2*PI*t*10)+3*cos(2*PI*t*20)+3*cos(2*PI*t*30)+2*cos(2*PI*t*40)+cos(2*PI*t*50);
        samples[i].imag = 0;
        // samples[i].real = (i%64 > 32) ? 1 : -1;
        // samples[i].imag = 0;

    }
}

void init_halves(complex *halves){
    for(uint16_t p=1,k=1 ; p<hlength ; p=(p<<1)+1,k++){
        uint16_t s = 1 << k; // Tamanho de cada sub-array
        for(uint16_t i=0 ; i<s ; i++){
            halves[p+i].real = cos(-2*PI*(double)i/(double)s);
            halves[p+i].imag = sin(-2*PI*(double)i/(double)s);
        }
    }
}

uint16_t reverse_index(uint16_t i){
    return  ((i & (1 << 0)) << 9) | 
            ((i & (1 << 1)) << 7) |
            ((i & (1 << 2)) << 5) |
            ((i & (1 << 3)) << 3) |
            ((i & (1 << 4)) << 1) |
            ((i & (1 << 5)) >> 1) |
            ((i & (1 << 6)) >> 3) |
            ((i & (1 << 7)) >> 5) |
            ((i & (1 << 8)) >> 7) |
            ((i & (1 << 9)) >> 9);
}

void bit_reverse_copy(complex *input, complex *output){
    for(uint16_t i=0 ; i<N ; i++){
        output[reverse_index(i)] = input[i];
    }
}

void fft(complex *input, complex *output, complex *halves){
    bit_reverse_copy(input, output);
    for(uint16_t k=1 ; k<=log2N ; k++){
        uint16_t s = 1 << k; // Tamanho de cada sub-array
        uint16_t halfs = (s>>1);

        for(uint16_t j=0 ; j<N ; j+=s){ // Percorre os sub-arrays, j é o índicie inicial no superarray
            for(uint16_t i=0 ; i<halfs ; i++){ // Percorre metade dos elementos
                complex t = output[i+j];
                complex u;

                // Multiplicação complexa
                cmul(output[i+j+halfs], halves[s-1 + i], u);

                // Soma complexa
                csum(t, u, output[i+j]);

                // Subtração complexa
                csub(t, u, output[i+j+halfs]);
            }
        }
    }
}

void print_fft(complex *fft, complex *samples, double sampleWindow){
    printf("{");
    printf("\"data\": [");
    for(uint16_t i=0 ; i<N/2 ; i++){
        double magpos = sqrt(pow(fft[i].real, 2)+pow(fft[i].imag,2));
        double magneg = sqrt(pow(fft[N-i].real, 2)+pow(fft[N-i].imag,2));
        double mag = (magpos+magneg)/1024;

        printf("[%lf,%lf]", (double)i/sampleWindow, mag);
        if(i<N/2-1) printf(",");
    }
    printf("]");
    printf(",\"samples\": [");
    for(uint16_t i=0 ; i<N ; i++){
        printf("[%lf,%lf]", (double)i*sampleWindow/N, samples[i].real);
        if(i<N-1) printf(",");
    }
    printf("]");
    printf("}\n");
}

int main(){
    complex samples[N];
    complex halves[hlength];
    complex samples_fft[N];
    init_samples(samples, 10);
    init_halves(halves);

    clock_t begin = clock();
    for(uint32_t i=0 ; i<n_bench_loops ; i++){
        fft(samples, samples_fft, halves);
    }
    clock_t end = clock();
    
    printf("%lf us/fft\n", (double)(end - begin) / CLOCKS_PER_SEC / n_bench_loops * 1000000);
    //print_fft(samples_fft, samples, 10);
    return 0;
}