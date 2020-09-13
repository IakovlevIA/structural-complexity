#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iostream>

using namespace std;

typedef struct point_t {
    float sx;
    float sy;
    float sz;
} spin;

void complexity_func(spin* A1, spin* A2, int n, int lambda, float* over, int l);

int main()
{
	using namespace std;
	
    std::ofstream outMatrix;
    std::ofstream out;
    std::ifstream inp;
    std::ifstream inp_fi1;
    std::ifstream inp_fi2;
    out.open("complexity_128.dat");

	char buf1[100], buf2[100];
	float overA[10];
    int n, exam, starting_point;
    int t0=time(0);
    float complexityA;
    int n_bl, num_conf, conf_count;
    int block_count = 0;
    int lambda, lambda0, compl_steps, compl_start;
    div_t div_count1, div_count2, start;
    
    float *compMatrix;

    inp.open("inp.dat");
    inp >> n >> compl_steps >> compl_start >> lambda0 >> num_conf;
    inp.close();
    
    spin *A1 = NULL;
    spin *A2 = NULL;

    size_t lsize = n * n * sizeof(spin);

    // Host memory allocation
    A1 = (spin*) malloc(lsize);
    A2 = (spin*) malloc(lsize);

    
    //--------------------conf loop--------------------------
    for (int c_i = 0; c_i < num_conf; c_i++){
            
        sprintf(buf1,"conf%d", c_i);
        inp_fi1.open(buf1);
        
        for (int i = 0; i < n * n; i++){
            inp_fi1 >> A1[i].sz;
        }
            
        inp_fi1.close();
        
        for (int i = 0; i < compl_steps; i++){
            lambda = pow(lambda0, i+1);
            if (i % 2 == 0){
                complexity_func(A1,  A2,  n, lambda, overA, i);
            }
            else{
                complexity_func(A2,  A1,  n, lambda, overA, i);
            }
        }
    
        complexityA = 0;
        for (int comp = compl_start; comp < compl_steps; comp++){
            complexityA += fabs(overA[comp]);
        }
        out << 4.5 - c_i * 4.0/150 << ' ' << complexityA << endl;

    }//end of conf loop
    out.close();
}


void complexity_func(spin* A1, spin* A2, int n, int lambda, float* over, int l){
    
    float overlap2 = 0, overlap1 = 0, overlap3 = 0, ren_r, o1 = 0, o2 = 0, o3 = 0;
    spin renorm;
    int n_bl = n/lambda;
    
    for (int n1_counti = 0; n1_counti < n_bl; n1_counti++){
        for (int n1_countj = 0; n1_countj < n_bl; n1_countj++){
            renorm.sz = 0;
            for (int count_i=0; count_i < lambda; count_i++){
                for (int count_j=0; count_j < lambda; count_j++){
                    renorm.sz += A1[(n1_counti * lambda + count_i) * n + n1_countj * lambda + count_j].sz;
                }
            }
            ren_r = lambda * lambda;
            for (int count_i=0; count_i < lambda; count_i++){
                for (int count_j=0; count_j < lambda; count_j++){
                    A2[(n1_counti * lambda + count_i) * n + n1_countj * lambda + count_j].sz = renorm.sz/ren_r;
                }
            }
        }
    }
    
    for (int i = 0; i < n * n; i++){
        overlap1 += A1[i].sz * A1[i].sz;
        overlap2 += A1[i].sz * A2[i].sz;
        overlap3 += A2[i].sz * A2[i].sz;
    }
    
    o1 = overlap1/n/n;
    o2 = overlap2/n/n;
    o3 = overlap3/n/n;
    over[l]= o2 - 0.5 * (o1 +  o3);
}

