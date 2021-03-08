#include "tfhe/tfhe.h"
#include <cstring>
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

// define 32 bit integer number manipulation
struct FheUInt32 {
    LweSample *sample;
    const LweParams *params;
    FheUInt32(const LweParams *params){
        this->params = params; // where param is CONSTANTLY READ-ONLY
        this->sample = new_LweSample_array(32, params); // 32 bit
    }

    FheUInt32(const FheUInt32 &x) {
        // deep clone 
        // cloning this type for future change
        // const int32_t n;
        // const double alpha_min;//le plus petit bruit tq sur
        // const double alpha_max;
        this->params = x.params; // where param is CONSTANTLY READ-ONLY
        for(int i = 0; i<32; ++i) {
            memcpy(&this->sample[i], x.sample + i, sizeof(Torus32) * x.params->n);
        }
    }
    ~FheUInt32() {
        delete_LweSample_array(32, this->sample);
    }
};

FheUInt32 fromUInt32(unsigned int x, TFheGateBootstrappingSecretKeySet* secretKeySet){
    FheUInt32 result(secretKeySet->params->in_out_params);
    for(int i = 0; i<32; ++i) {
        bootsSymEncrypt(result.sample+i, x & 1, secretKeySet);
        x >>= 1;
    }
    return result;
}

unsigned int toUInt32(FheUInt32 &x, TFheGateBootstrappingSecretKeySet* secretKeySet) {
    unsigned int result = 0;
    for(int i = 31; i>=0; --i) {
        int u = bootsSymDecrypt(x.sample+i, secretKeySet);
        result = 2*result + u;
    }
    return result;
}

FheUInt32 add(
    FheUInt32 &a, FheUInt32 &b, // avoid copying of object 
    const TFheGateBootstrappingCloudKeySet *cloud_key, const FheUInt32 &car){
    const int nb_bits = 32;
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    FheUInt32 result(in_out_params);
    LweSample *carry = car.sample; // new_LweSample_array(2, in_out_params);
    LweSample *temp = new_LweSample_array(3, in_out_params);

    LweSample *sum = result.sample;
    LweSample *x = a.sample, *y = b.sample;

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
    // bootsCOPY(sum + nb_bits, carry, cloud_key);

    delete_LweSample_array(3, temp);
    return result;
}

FheUInt32 sub(
    FheUInt32 &a, FheUInt32 &b, // avoid copying of object 
    const TFheGateBootstrappingCloudKeySet *cloud_key, const FheUInt32 &car){
    const int nb_bits = 32;
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    FheUInt32 result(in_out_params);
    LweSample *carry = car.sample; // new_LweSample_array(2, in_out_params);
    LweSample *temp = new_LweSample_array(3, in_out_params);

    LweSample *sum = result.sample;
    LweSample *x = a.sample, *y = b.sample;

    for (int32_t i = 0; i < nb_bits; ++i) {
        //sumi = xi XOR yi XOR carry
        bootsXOR(temp, x + i, y + i, cloud_key); // temp = xi XOR yi
        bootsXOR(sum + i, temp, carry, cloud_key);
        // carry = (xi AND yi) XOR (carry(i-1) AND (xi XOR yi))
        // = (!A & carry) | (!A & B) | (B & carry)
        bootsNOT(temp, x + i, cloud_key); // !A
        bootsAND(temp + 1, temp, carry, cloud_key); // !A & carry
        bootsAND(temp + 2, temp, y + i, cloud_key); // !A & B
        bootsOR(temp, temp + 1, temp + 2, cloud_key); // (!A & carry) | (!A & B)
        bootsAND(temp + 1, y + i, carry, cloud_key); // (B & carry)
        bootsOR(carry, temp, temp + 1, cloud_key); // (!A & carry) | (!A & B) | (B & carry)
    }
    delete_LweSample_array(3, temp);
    return result;
}

// void full_adder_circuit(LweSample* a, LweSample *b, LweSample *carry, LweSample *sum, const TFheGateBootstrappingCloudKeySet *cloud_key) {
//     const LweParams *in_out_params = cloud_key->params->in_out_params;
//     const LweSample *temp = new_LweSample_array(3, in_out_params);
//     // a xor b => temp, 
//     //            temp xor carry => sum
//     //            temp and carry => first_carry
//     // a and b => second_carry
//     // carry = first_carry or second_carry
//     bootsXOR(temp, a, b, cloud_key);
//     bootsXOR(sum, temp, carry, cloud_key);
//     bootsAND(temp+1, temp, carry, cloud_key);
//     bootsAND(carry, a, b, cloud_key);
//     bootsAND(carry, carry, temp+1, cloud_key);

//     delete_LweSample_array(2, temp);
// }
FheUInt32 mul(
    FheUInt32 &a, FheUInt32 &b, // avoid copying of object 
    const TFheGateBootstrappingCloudKeySet *cloud_key, const FheUInt32 &car){
    const int nb_bits = 32;
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    FheUInt32 result(in_out_params);
    LweSample *carry = car.sample; // new_LweSample_array(2, in_out_params);
    LweSample *temp = new_LweSample_array(3, in_out_params);
    LweSample *sum = result.sample;
    LweSample *x = a.sample, *y = b.sample;

    // use binary shifting stragegy
    for (int32_t j = 0; j < nb_bits; ++j) { // iterate over B
        bootsCONSTANT(carry, 0, cloud_key);
        for (int32_t i = 0; i < nb_bits; ++i) { // iterate over A
            LweSample *Xi = x+i;
            LweSample *Yj = y+j;
            if (j+i > 31) break;
            // x[i] * y[j]
            bootsAND(temp, Xi, Yj, cloud_key);
            LweSample *c = sum+(j+i);
            // Add to sum[j+i] and Carry
            bootsXOR(temp, Xi, Yj, cloud_key); // temp = a xor b
            bootsXOR(c, temp, carry, cloud_key); // c = temp xor carry 
            bootsAND(temp, temp, carry, cloud_key); // temp = (a xor b) and carry
            bootsAND(carry, Xi, Yj, cloud_key); // carry = a and b
            bootsAND(carry, carry, temp, cloud_key); // carry = (a and b) and ((a xor b) and carry)
        }
        // dont need to consider "redudant" bit
    }
    delete_LweSample_array(3, temp);
    return result;
}


int main(){

    TFheGateBootstrappingParameterSet *parameneterSet = new_default_gate_bootstrapping_parameters(16);
    const LweParams* in_out_param = parameneterSet->in_out_params;
    TFheGateBootstrappingSecretKeySet *keySet = new_random_gate_bootstrapping_secret_keyset(parameneterSet);
    // cout << in_out_param->n;
    // client decrypt data with its own key
    FheUInt32 
            x = fromUInt32(37, keySet), 
            y = fromUInt32(20, keySet);
    FheUInt32 car = fromUInt32(0, keySet);

    cout << "Start mul 2 encrypted number:" << endl;
    auto marked = time(0);
    // computation takes place on cloud over encrypted data
    FheUInt32 result = mul(x, y, &keySet->cloud, car);
    // client get result and decrypt
    cout << toUInt32(x, keySet) << " * " << toUInt32(y, keySet) << " = " << toUInt32(result, keySet) << endl;

    cout << "Time ellapsed:" << time(0)-marked << " second(s)";
}