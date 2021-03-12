#include "tfhe/tfhe.h"
#include <iostream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include "tfhe/tfhe.h"
#include "tfhe/polynomials.h"
#include "tfhe/lwesamples.h"
#include "tfhe/lwekey.h"
#include "tfhe/lweparams.h"
#include "tfhe/tlwe.h"
#include "tfhe/tgsw.h"

using namespace std;
static const Torus32 MU = modSwitchToTorus32(1, 8);

time_t now(){
    struct timeval tv;
    gettimeofday(&tv, 0);
    return (tv.tv_sec * 1000 + tv.tv_usec/1000);
}

void nonbootsXOR(LweSample *result, const LweSample *ca, const LweSample *cb, const TFheGateBootstrappingCloudKeySet *bk) {
    // static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;

    // LweSample *temp_result = new_LweSample(in_out_params);
    //compute: (0,1/4) + 2*(ca + cb)
    static const Torus32 XorConst = modSwitchToTorus32(1, 4);
    lweNoiselessTrivial(result, XorConst, in_out_params);
    lweAddMulTo(result, 2, ca, in_out_params);
    lweAddMulTo(result, 2, cb, in_out_params);
}

void nonbootsAND(LweSample *result, const LweSample *ca, const LweSample *cb, const TFheGateBootstrappingCloudKeySet *bk) {
    
    const LweParams *in_out_params = bk->params->in_out_params;
    //compute: (0,-1/8) + ca + cb
    static const Torus32 AndConst = modSwitchToTorus32(-1, 8);
    lweNoiselessTrivial(result, AndConst, in_out_params);
    lweAddTo(result, ca, in_out_params);
    lweAddTo(result, cb, in_out_params);
    // delete_LweSample(temp_result);
}

void nonbootsOR(LweSample *result, const LweSample *ca, const LweSample *cb, const TFheGateBootstrappingCloudKeySet *bk) {
    
    const LweParams *in_out_params = bk->params->in_out_params;
    LweSample *temp_result = new_LweSample(in_out_params);
    //compute: (0,1/8) + ca + cb
    static const Torus32 OrConst = modSwitchToTorus32(1, 8);
    lweNoiselessTrivial(result, OrConst, in_out_params);
    lweAddTo(result, ca, in_out_params);
    lweAddTo(result, cb, in_out_params);
    // tfhe_bootstrap_woKS_FFT(result, bk->bkFFT, MU, temp_result);
    // delete_LweSample(temp_result);
}


void nonbootsCONSTANT(LweSample *result, int32_t value, const TFheGateBootstrappingCloudKeySet *bk) {
    // const LweParams *extracted_params = &bk->params->tgsw_params->tlwe_params->extracted_lweparams;
    const LweParams *in_out_params = bk->params->in_out_params;
    // static const Torus32 MU = modSwitchToTorus32(1, 8);
    lweNoiselessTrivial(result, value ? MU : -MU, in_out_params);
}

void full_adder_optimized(LweSample *sum, 
    const LweSample *x, const LweSample *y, 
    const int32_t nb_bits, const LweParams *in_out_params, const TFheGateBootstrappingCloudKeySet* cloud_key) {
    const LweParams *extracted_params = &cloud_key->params->tgsw_params->tlwe_params->extracted_lweparams;

    LweSample *carry = new_LweSample_array(2, in_out_params); // 1 bit
    LweSample *temp = new_LweSample_array(4, in_out_params);
    LweSample *exported_carry = new_LweSample(in_out_params);
    LweSample *u = new_LweSample(extracted_params);
    LweSample *tmp = new_LweSample(extracted_params);

    nonbootsCONSTANT(carry, 0, cloud_key);

    for (int32_t i = 0; i < nb_bits; ++i) {

        time_t marker = now();
        //sumi = xi XOR yi XOR carry(i-1) 
        nonbootsXOR(temp, x + i, y + i, cloud_key); // temp = xi XOR yi
        nonbootsXOR(temp + 3, temp, carry, cloud_key);
        tfhe_bootstrap_woKS_FFT(u, cloud_key->bkFFT, MU, temp + 3);

        // carry = (xi AND yi) XOR (carry(i-1) AND (xi XOR yi))
        nonbootsAND(temp + 1, x + i, y + i, cloud_key); // temp1 = xi AND yi
        nonbootsAND(temp + 2, carry, temp, cloud_key); // temp2 = carry AND temp
        nonbootsXOR(carry + 1, temp + 1, temp + 2, cloud_key);
        static const Torus32 AdderConst = modSwitchToTorus32(1, 8);
        lweNoiselessTrivial(tmp, AdderConst, extracted_params);
        lweAddTo(tmp, u, extracted_params);
        // Key switching
        lweKeySwitch(sum + i, cloud_key->bkFFT->ks, tmp);
        // copy carry flag
        lweCopy(carry, carry + 1, in_out_params);
        cout << "1 bit calculation cosume: " << now() - marker << endl;
    }
    lweKeySwitch(exported_carry, cloud_key->bkFFT->ks, carry);
    bootsCOPY(sum + nb_bits, carry, cloud_key);
    delete_LweSample(exported_carry);
    delete_LweSample_array(3, temp);
    delete_LweSample_array(2, carry);
}

void full_adder(LweSample *sum, 
    const LweSample *x, const LweSample *y, 
    const int32_t nb_bits, const LweParams *in_out_params, const TFheGateBootstrappingCloudKeySet* cloud_key) {
    LweSample *carry = new_LweSample_array(2, in_out_params);
    bootsCONSTANT(carry, 0, cloud_key);
    LweSample *temp = new_LweSample_array(3, in_out_params);
    for (int32_t i = 0; i < nb_bits; ++i) {

        time_t marker = now();
        //sumi = xi XOR yi XOR carry(i-1) 
        bootsXOR(temp, x + i, y + i, cloud_key); // temp = xi XOR yi
        bootsXOR(sum + i, temp, carry, cloud_key);

        // carry = (xi AND yi) XOR (carry(i-1) AND (xi XOR yi))
        bootsAND(temp + 1, x + i, y + i, cloud_key); // temp1 = xi AND yi
        bootsAND(temp + 2, carry, temp, cloud_key); // temp2 = carry AND temp
        bootsXOR(carry + 1, temp + 1, temp + 2, cloud_key);
        bootsCOPY(carry, carry + 1, cloud_key);

        cout << "1 bit calculation cosume: " << now() - marker << endl;
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

    cout << "Running: " << endl;

    cout << "Optimize with non-keyswitch:" << endl;
    // add placeholder for result
    LweSample* z = new_LweSample_array(2, in_out_param); // 2 bit
    
    full_adder_optimized(z, x, y, 1, in_out_param, &secretKeySet->cloud);
    for(int i = 0; i<2; ++i) {
        int decrypted_bit = bootsSymDecrypt(z+i, secretKeySet);
        cout << decrypted_bit;
    }
    cout << endl << endl;
    cout << "Keyswitch:" << endl;
    full_adder(z, x, y, 1, in_out_param, &secretKeySet->cloud);
    for(int i = 0; i<2; ++i) {
        int decrypted_bit = bootsSymDecrypt(z+i, secretKeySet);
        cout << decrypted_bit;
    }

    cout << endl;
    delete_gate_bootstrapping_parameters(parameneterSet);
    delete_gate_bootstrapping_secret_keyset(secretKeySet);
}