#include "tfhe/tfhe.h"
#include <iostream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe/polynomials.h"
#include "tfhe/lwesamples.h"
#include "tfhe/lweparams.h"

using namespace std;

void print_sample(LweSample sample){
    cout << sample.
}
int main(){
    // generate param and secret keys
    TFheGateBootstrappingParameterSet* parameneterSet = new_default_gate_bootstrapping_parameters(16);
    const LweParams* in_out_param = parameneterSet->in_out_params;

    // generate secret key based on paramenter set
    TFheGateBootstrappingSecretKeySet* secretKeySet = new_random_gate_bootstrapping_secret_keyset(parameneterSet);

    // encrypt the boolean '1'
    int bit = 0;
    LweSample* sample = new_LweSample_array(1, in_out_param); // 1 bit
    bootsSymEncrypt(sample, bit, secretKeySet);
    
    LweSample* notSample = new_LweSample_array(1, in_out_param); // 1 bit
    bootsNOT(notSample, sample, &secretKeySet->cloud);
    
    // decrypt
    int decrypted_bit = bootsSymDecrypt(notSample, secretKeySet);
    cout << "Decrypted bit:" << decrypted_bit << endl;
    cout << "Expected bit:" << !bit << endl;

    delete_LweSample(sample);
    delete_LweSample(notSample);
}
