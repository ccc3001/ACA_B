#include "new_aca.h"
#include "aca_b.h"
#include "fullmatrix.h"
#include <stdio.h>
#include <time.h>
#include "rkmatrix.h"

int main(){
    printf("Testing the functionality of :\n"
            "- ACA \n"
            "- BACA \n");
    
    clock_t start_time; 
    clock_t end_time;
    double elapsed; 
    double *rk_;
    double *elapsed_;
    rk_ =allocate_doubles(11);
    elapsed_ =allocate_doubles(11);
    for (int i = 1; i < 11; i++)
    {
        int n= i * 1000 ;
        printf("\nGenerating  Matrix of size %dx%d\n",n,n);
        pfullmatrix A = new_random_fullmatrix(n,n,12);
        
        printf("Calculating ACA Matrix\n");
        start_time =clock(); 
        prkmatrix RK_ACA = b_aca_rkmatrix_new(0.001,2,A);
        end_time =clock(); 
        elapsed= (double)(end_time - start_time) / CLOCKS_PER_SEC;
        printf("time:%.4g\n",elapsed);
        printf("rank:%d\n",RK_ACA->kt);
        printf("time/rank:%.4g",elapsed/RK_ACA->kt);
        elased_[i-1]=elapsed;
        rk_[i-1]=RK_ACA->kt;
        
    }
    printf("rk:[");
    for (int i = 0; i < n; i++) {
        printf("%d,", arr[i]);
    }
     
    return 0;
}
