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
    LweSample *temp_result = new_LweSample(in_out_params);
    //compute: (0,1/4) + 2*(ca + cb)
    static const Torus32 XorConst = modSwitchToTorus32(1, 4);
    lweNoiselessTrivial(temp_result, XorConst, in_out_params);
    lweAddMulTo(temp_result, 2, ca, in_out_params);
    lweAddMulTo(temp_result, 2, cb, in_out_params);
    lweCopy(result, temp_result, in_out_params);
}

void nonbootsXOR3(LweSample *result, const LweSample *ca, const LweSample *cb, const LweSample *cc, const TFheGateBootstrappingCloudKeySet *bk) {
    // static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;
    const LweParams *extracted_params = &bk->params->tgsw_params->tlwe_params->extracted_lweparams;

    LweSample *temp_result = new_LweSample(in_out_params);
    LweSample 
        *u1 = new_LweSample(extracted_params),
        *u2 = new_LweSample(extracted_params),
        *u = new_LweSample(extracted_params);

    //compute: (0,1/4) + 2*(ca + cb)
    static const Torus32 XorConst = modSwitchToTorus32(1, 4);
    lweNoiselessTrivial(temp_result, XorConst, in_out_params);
    lweAddMulTo(temp_result, 2, ca, in_out_params);
    lweAddMulTo(temp_result, 2, cb, in_out_params);
    tfhe_bootstrap_woKS_FFT(u1, bk->bkFFT, MU, temp_result);

    lweNoiselessTrivial(temp_result, XorConst, in_out_params);
    lweAddMulTo(temp_result, -2, cc, in_out_params);
    tfhe_bootstrap_woKS_FFT(u2, bk->bkFFT, MU, temp_result);
    
    lweNoiselessTrivial(u, XorConst, extracted_params);
    lweAddMulTo(u, 2, u1, extracted_params);
    lweAddMulTo(u, 2, u2, extracted_params);
    // sum bit
    lweKeySwitch(result, bk->bkFFT->ks, u);
}

void nonbootsAND(LweSample *result, const LweSample *ca, const LweSample *cb, const TFheGateBootstrappingCloudKeySet *bk) {
    const LweParams *in_out_params = bk->params->in_out_params;
    //compute: (0,-1/8) + ca + cb
    static const Torus32 AndConst = modSwitchToTorus32(-1, 8);
    LweSample *temp_result = new_LweSample(in_out_params);
    lweNoiselessTrivial(temp_result, AndConst, in_out_params);
    lweAddTo(temp_result, ca, in_out_params);
    lweAddTo(temp_result, cb, in_out_params);
    // delete_LweSample(temp_result);
    lweCopy(result, temp_result, in_out_params);
}

void nonbootsOR(LweSample *result, const LweSample *ca, const LweSample *cb, const TFheGateBootstrappingCloudKeySet *bk) {
    
    const LweParams *in_out_params = bk->params->in_out_params;
    LweSample *temp_result = new_LweSample(in_out_params);
    //compute: (0,1/8) + ca + cb
    static const Torus32 OrConst = modSwitchToTorus32(1, 8);
    lweNoiselessTrivial(result, OrConst, in_out_params);
    lweAddTo(temp_result, ca, in_out_params);
    lweAddTo(temp_result, cb, in_out_params);
    // tfhe_bootstrap_woKS_FFT(result, bk->bkFFT, MU, temp_result);
    // delete_LweSample(temp_result);
    lweCopy(result, temp_result, in_out_params);
}


void nonbootsCONSTANT(LweSample *result, int32_t value, const TFheGateBootstrappingCloudKeySet *bk) {
    // const LweParams *extracted_params = &bk->params->tgsw_params->tlwe_params->extracted_lweparams;
    const LweParams *in_out_params = bk->params->in_out_params;
    // static const Torus32 MU = modSwitchToTorus32(1, 8);
    lweNoiselessTrivial(result, value ? MU : -MU, in_out_params);
}

void full_adder_gate(
    LweSample *sum, 
    LweSample *carry, // in-out
    const LweSample *x, const LweSample *y, 
    const TFheGateBootstrappingCloudKeySet* cloud_key) {
    
    static const Torus32 AdderConst = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    const LweParams *extracted_params = &cloud_key->params->tgsw_params->tlwe_params->extracted_lweparams;
    LweSample *tmp = new_LweSample(extracted_params);
    LweSample *temp = new_LweSample_array(3, in_out_params);

    LweSample *temp_result = new_LweSample(in_out_params);
    LweSample 
        *u1 = new_LweSample(extracted_params),
        *u2 = new_LweSample(extracted_params),
        *u3 = new_LweSample(extracted_params),
        *u = new_LweSample(extracted_params);

    // compute: XOR(XOR(x,y), carry)
    static const Torus32 XorConst = modSwitchToTorus32(1, 4);
    // u1 = x xor y 
    lweNoiselessTrivial(temp_result, XorConst, in_out_params);
    lweAddMulTo(temp_result, 2, x, in_out_params);
    lweAddMulTo(temp_result, 2, y, in_out_params);
    tfhe_bootstrap_woKS_FFT(u1, cloud_key->bkFFT, MU, temp_result);

    // u2 = carry
    lweNoiselessTrivial(temp_result, XorConst, in_out_params);
    lweAddMulTo(temp_result, 2, carry, in_out_params);
    tfhe_bootstrap_woKS_FFT(u2, cloud_key->bkFFT, MU, temp_result);
    
    // u = u1 xor u2
    lweNoiselessTrivial(u, XorConst, extracted_params);
    lweAddMulTo(u, 2, u1, extracted_params);
    lweSubMulTo(u, 2, u2, extracted_params);
    // produce sum bit
    lweKeySwitch(sum, cloud_key->bkFFT->ks, u);

    // ---
    // compute carry = (x and y) or ((x xor y) and carry)
    // u1 = x xor y
    lweNoiselessTrivial(temp_result, XorConst, in_out_params);
    lweAddMulTo(temp_result, 2, x, in_out_params);
    lweAddMulTo(temp_result, 2, y, in_out_params);
    tfhe_bootstrap_woKS_FFT(u1, cloud_key->bkFFT, MU, temp_result);

    // u3 = x and y
    static const Torus32 AndConst = modSwitchToTorus32(-1, 8);
    lweNoiselessTrivial(temp_result, AndConst, in_out_params);
    lweAddTo(temp_result, x, in_out_params);
    lweAddTo(temp_result, y, in_out_params);
    tfhe_bootstrap_woKS_FFT(u3, cloud_key->bkFFT, MU, temp_result);
    
    // // u2 = carry
    lweNoiselessTrivial(temp_result, AndConst, in_out_params);
    lweAddTo(temp_result, carry, in_out_params);
    tfhe_bootstrap_woKS_FFT(u2, cloud_key->bkFFT, MU, temp_result);
    
    // // u = u1 and u2
    lweNoiselessTrivial(u, AndConst, extracted_params);
    lweAddTo(u, u1, extracted_params);
    lweAddTo(u, u2, extracted_params);
    lweKeySwitch(carry, cloud_key->bkFFT->ks, u);
    tfhe_bootstrap_woKS_FFT(u2, cloud_key->bkFFT, MU, carry);

    static const Torus32 OrConst = modSwitchToTorus32(1, 8);
    lweNoiselessTrivial(u, OrConst, extracted_params);
    lweAddTo(u, u3, extracted_params);
    lweAddTo(u, u2, extracted_params);
    lweKeySwitch(carry, cloud_key->bkFFT->ks, u);

    delete_LweSample(u);
    delete_LweSample(u1);
    delete_LweSample(u2);
    delete_LweSample(u3);
    delete_LweSample(tmp);
    delete_LweSample_array(3, temp);
}

// void full_adder_optimized(LweSample *sum, 
//     const LweSample *x, const LweSample *y, 
//     const int32_t nb_bits, const LweParams *in_out_params, const TFheGateBootstrappingCloudKeySet* cloud_key) {
//     const LweParams *extracted_params = &cloud_key->params->tgsw_params->tlwe_params->extracted_lweparams;

//     LweSample *carry = new_LweSample_array(2, in_out_params); // 1 bit
//     LweSample *temp = new_LweSample_array(4, in_out_params);
//     LweSample *exported_carry = new_LweSample(in_out_params);
//     LweSample *u = new_LweSample(extracted_params);
//     LweSample *tmp = new_LweSample(extracted_params);
//     nonbootsCONSTANT(carry, 0, cloud_key);
//     static const Torus32 AdderConst = modSwitchToTorus32(1, 8);

//     for (int32_t i = 0; i < nb_bits; ++i) {

//         time_t marker = now();
//         //sumi = xi XOR yi XOR carry(i-1) 
//         // nonbootsXOR(temp, x + i, y + i, cloud_key); // temp = xi XOR yi
//         // // nonbootsXOR(temp + 3, temp, carry, cloud_key);
//         // tfhe_bootstrap_woKS_FFT(u, cloud_key->bkFFT, MU, temp);

//         // // carry = (xi AND yi) XOR (carry(i-1) AND (xi XOR yi))
//         // nonbootsAND(temp + 1, x + i, y + i, cloud_key); // temp1 = xi AND yi
//         // nonbootsAND(temp + 2, carry, temp, cloud_key); // temp2 = carry AND temp
//         // nonbootsXOR(carry + 1, temp + 1, temp + 2, cloud_key);
        
//         // lweNoiselessTrivial(tmp, AdderConst, extracted_params);
//         // lweAddTo(tmp, u, extracted_params);
//         // // Key switching
//         // lweKeySwitch(sum + i, cloud_key->bkFFT->ks, tmp);
//         // // copy carry flag
//         // lweCopy(carry, carry + 1, in_out_params);
//         cout << "1 bit calculation cosume: " << now() - marker << endl;
//     }
//     lweKeySwitch(exported_carry, cloud_key->bkFFT->ks, carry);
//     bootsCOPY(sum + nb_bits, carry, cloud_key);
//     delete_LweSample(exported_carry);
//     delete_LweSample_array(3, temp);
//     delete_LweSample_array(2, carry);
// }

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

    for(int i = 0; i<2; ++i) 
        for(int j = 0; j<2; ++j) 
            for(int k = 0; k<2; ++k) {

                        
                // encrypt the boolean '1'
                LweSample* x = new_LweSample(in_out_param); // 1 bit
                bootsSymEncrypt(x, i, secretKeySet);
                
                // encrypt the boolean '0'
                LweSample* y = new_LweSample(in_out_param); // 1 bit
                bootsSymEncrypt(y, j, secretKeySet);

                cout << "Running: " << endl;
                cout << "Optimize with non-keyswitch:" << endl;

                // add placeholder for result
                LweSample 
                    *z = new_LweSample(in_out_param), 
                    *carry = new_LweSample(in_out_param),
                    *result = new_LweSample(in_out_param);

                bootsCONSTANT(carry, k, &secretKeySet->cloud);
                // full_adder_gate(z, carry, x, y, &secretKeySet->cloud);
                full_adder_gate(result, carry, x, y, &secretKeySet->cloud);
                int decrypted_bit = bootsSymDecrypt(result, secretKeySet);
                int carry_bit = bootsSymDecrypt(carry, secretKeySet);
                cout << "Sum bit: " << i << " " << j << " " << k << " == " << decrypted_bit << ", expected: " << (i^j^k) << ", actual carry: " << carry_bit << " " << endl;
            }

    // full_adder_optimized(z, x, y, 1, in_out_param, &secretKeySet->cloud);
    // for(int i = 0; i<2; ++i) {
    //     int decrypted_bit = bootsSymDecrypt(z+i, secretKeySet);
    //     cout << decrypted_bit;
    // }
    // cout << endl << endl;
    // cout << "Keyswitch:" << endl;
    // full_adder(z, x, y, 1, in_out_param, &secretKeySet->cloud);
    // for(int i = 0; i<2; ++i) {
    //     int decrypted_bit = bootsSymDecrypt(z+i, secretKeySet);
    //     cout << decrypted_bit;
    // }

    delete_gate_bootstrapping_parameters(parameneterSet);
    delete_gate_bootstrapping_secret_keyset(secretKeySet);
}