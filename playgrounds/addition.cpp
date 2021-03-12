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

void full_adder(LweSample *sum, 
    const LweSample *x, const LweSample *y, 
    const int32_t nb_bits, const LweParams *in_out_params, LweSample *carry, const TFheGateBootstrappingCloudKeySet* cloud_key) {
    LweSample *temp = new_LweSample_array(3, in_out_params);
    for (int32_t i = 0; i < nb_bits; ++i) {
        //sumi = xi XOR yi XOR carry(i-1) 
        bootsXOR(temp, x + i, y + i, cloud_key); // temp = xi XOR yi
        bootsXOR(sum + i, temp, carry, cloud_key);

        // carry = (xi AND yi) XOR (carry(i-1) AND (xi XOR yi))
        bootsAND(temp + 1, x + i, y + i, cloud_key); // temp1 = xi AND yi
        bootsAND(temp + 2, carry, temp, cloud_key); // temp2 = carry AND temp
        bootsXOR(carry + 1, temp + 1, temp + 2, cloud_key);
        bootsCOPY(carry, carry + 1, cloud_key);
    }
    bootsCOPY(sum + nb_bits, carry, cloud_key);

    delete_LweSample_array(3, temp);
    delete_LweSample_array(2, carry);
}

int main(){
    // generate param and secret keys
    TFheGateBootstrappingParameterSet* parameneterSet = new_default_gate_bootstrapping_parameters(16);
    const LweParams* in_out_param = parameneterSet->in_out_params;

    // generate secret key based on paramenter set
    TFheGateBootstrappingSecretKeySet* secretKeySet = new_random_gate_bootstrapping_secret_keyset(parameneterSet);

    // encrypt the boolean '1'
    LweSample* x = new_LweSample_array(1, in_out_param); // 1 bit
    bootsSymEncrypt(x, 1, secretKeySet);
    
    // encrypt the boolean '0'
    LweSample* y = new_LweSample_array(1, in_out_param); // 1 bit
    bootsSymEncrypt(y, 1, secretKeySet);

    // add placeholder for result
    LweSample* z = new_LweSample_array(2, in_out_param); // 2 bit
    LweSample* carry = new_LweSample_array(2, in_out_param); // 1 bit
    bootsSymEncrypt(carry, 0, secretKeySet);
    
    full_adder(z, x, y, 1, in_out_param, carry, &secretKeySet->cloud);

    for(int i = 0; i<2; ++i) {
        int decrypted_bit = bootsSymDecrypt(z+i, secretKeySet);
        cout << decrypted_bit;
        // expected to see 01, which will be 2 in DEC
    }

    delete_gate_bootstrapping_parameters(parameneterSet);
    delete_gate_bootstrapping_secret_keyset(secretKeySet);
}