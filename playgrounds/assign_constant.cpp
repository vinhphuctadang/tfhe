#include "tfhe/tfhe.h"
#include <iostream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe/tfhe.h"
#include "tfhe/polynomials.h"
#include "tfhe/lwesamples.h"
#include "tfhe/lwekey.h"
#include "tfhe/lweparams.h"
#include "tfhe/tlwe.h"
#include "tfhe/tgsw.h"
using namespace std;

int main(){
     // generate param and secret keys
    TFheGateBootstrappingParameterSet* parameneterSet = new_default_gate_bootstrapping_parameters(128);
    const LweParams* in_out_param = parameneterSet->in_out_params;

    // generate secret key based on paramenter set
    TFheGateBootstrappingSecretKeySet* secretKeySet = new_random_gate_bootstrapping_secret_keyset(parameneterSet);

    int bit_count = 16;
    LweSample* a = new_LweSample_array(bit_count, in_out_param);
    bootsCONSTANT(a, 1, &secretKeySet->cloud);
    bootsCONSTANT(a+1, 0, &secretKeySet->cloud);
    bootsCONSTANT(a+2, 1, &secretKeySet->cloud);

    for(int i = 0; i<bit_count; ++i) {
        cout << bootsSymDecrypt(a+i, secretKeySet);
        // expected '1010000000000000'
    }
}