#include "tfhe/tfhe.h"
#include <iostream>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe/polynomials.h"
#include "tfhe/lwesamples.h"
#include "tfhe/lweparams.h"

using namespace std;
int main(){
    LweParams *params = new_LweParams(512, 0.2, 0.5); 
    int32_t n = params->n;
    LweKey *key = new_LweKey(params);
    LweSample *cipher = new_LweSample(params);
    Torus32 mu = dtot32(0.5);
    // noise tolerance
    double alpha = 0.0625;
    Torus32 phi;
    double message;
    int32_t Msize = 2;
    lweKeyGen(key);
    lweSymEncrypt(cipher, mu, alpha, key);
    cout << "a = [";
    for (int32_t i = 0; i < n - 1; ++i) cout << t32tod(cipher->a[i]) << ", ";
    cout << t32tod(cipher->a[n - 1]) << "]" << endl;
    cout << "b = " << t32tod(cipher->b) << endl;
}
