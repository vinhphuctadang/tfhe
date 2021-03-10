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

    FheUInt32(const FheUInt32 &x) { // only allow copying from the same params
        // deep clone 
        this->params = x.params; // where param is CONSTANTLY READ-ONLY
        for(int i = 0; i<32; ++i) {
            lweCopy(this->sample + i, x.sample + i, this->params);
        }
    }
    ~FheUInt32() {
        // cout << "DEALLOCATE PERFORMED TO " << this->sample << endl;
        delete_LweSample_array(32, this->sample);
    }

    FheUInt32& operator=(const FheUInt32 &x)=delete;
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
    const TFheGateBootstrappingCloudKeySet *cloud_key){
    const int nb_bits = 32;
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    FheUInt32 result(in_out_params);
    LweSample *carry = new_LweSample(in_out_params);
    bootsCONSTANT(carry, 0, cloud_key);

    LweSample *temp = new_LweSample_array(2, in_out_params);
    LweSample *sum = result.sample;
    LweSample *x = a.sample, *y = b.sample;

    for (int32_t i = 0; i < nb_bits; ++i) {
        bootsXOR(temp, x + i, y + i, cloud_key); // temp = xi XOR yi
        bootsAND(temp + 1, x + i, y + i, cloud_key); // temp1 = xi AND yi
        bootsXOR(sum + i, temp, carry, cloud_key); // sum = temp XOR carry = (xi XOR yi) XOR carry
        bootsAND(carry, temp, carry, cloud_key); // carry = (xi XOR yi) and carry
        bootsOR(carry, carry, temp + 1, cloud_key); // carry = (((xi XOR yi) and carry) or (xi AND yi))
    }

    delete_LweSample_array(3, temp);
    return result;
}

FheUInt32 sub(
    FheUInt32 &a, FheUInt32 &b, // avoid copying of object 
    const TFheGateBootstrappingCloudKeySet *cloud_key){
    const int nb_bits = 32;
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    FheUInt32 result(in_out_params);

    LweSample *carry = new_LweSample(in_out_params);
    bootsCONSTANT(carry, 0, cloud_key);

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
        bootsAND(temp + 1, temp, carry, cloud_key);     // !A & carry
        bootsAND(temp + 2, temp, y + i, cloud_key);    // !A & B
        bootsOR(temp, temp + 1, temp + 2, cloud_key); // (!A & carry) | (!A & B)
        bootsAND(temp + 1, y + i, carry, cloud_key); // (B & carry)
        bootsOR(carry, temp, temp + 1, cloud_key);  // (!A & carry) | (!A & B) | (B & carry)
    }
    delete_LweSample_array(3, temp);
    return result;
}


FheUInt32 mul(
    FheUInt32 &a, FheUInt32 &b, // avoid copying of object 
    const TFheGateBootstrappingCloudKeySet *cloud_key){
    const int nb_bits = 8;
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    FheUInt32 result(in_out_params);
    LweSample *carry = new_LweSample(in_out_params);

    LweSample *temp = new_LweSample_array(3, in_out_params);
    LweSample *sum = result.sample;
    LweSample *x = a.sample, *y = b.sample;
    for (int32_t j = 0; j < nb_bits; ++j) { // iterate over B
        bootsCONSTANT(sum + j, 0, cloud_key);
    }
    // use binary shifting stragegy
    for (int32_t j = 0; j < nb_bits; ++j) { // iterate over B
        bootsCONSTANT(carry, 0, cloud_key);
        for (int32_t i = 0; i < nb_bits; ++i) { // iterate over A
            LweSample *Xi = x+i;
            LweSample *Yj = y+j;
            if (j+i >= nb_bits) break;
            // x[i] * y[j]
            bootsAND(temp, Xi, Yj, cloud_key);
            LweSample *c = sum+(j+i);
            // sum[j+i] = sum[j+i] + temp + carry
            bootsXOR(temp + 1, c, temp, cloud_key);
            bootsAND(temp + 2, c, temp, cloud_key);
            bootsXOR(c, temp + 1, carry, cloud_key);
            bootsAND(carry, temp + 1, carry, cloud_key); 
            bootsOR(carry, carry, temp + 2, cloud_key);
        }
        // dont need to consider "redudant" bit
    }
    delete_LweSample_array(3, temp);
    delete_LweSample(carry);
    return result;
}

FheUInt32 sub_out_carry(
    FheUInt32 &a, FheUInt32 &b, // avoid copying of object 
    const TFheGateBootstrappingCloudKeySet *cloud_key, LweSample* carry){ 
    const int nb_bits = 32;
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    FheUInt32 result(in_out_params);

    // LweSample *carry = new_LweSample(in_out_params);
    bootsCONSTANT(carry, 0, cloud_key);
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
        bootsAND(temp + 1, temp, carry, cloud_key);     // !A & carry
        bootsAND(temp + 2, temp, y + i, cloud_key);    // !A & B
        bootsOR(temp, temp + 1, temp + 2, cloud_key); // (!A & carry) | (!A & B)
        bootsAND(temp + 1, y + i, carry, cloud_key); // (B & carry)
        bootsOR(carry, temp, temp + 1, cloud_key);  // (!A & carry) | (!A & B) | (B & carry)
    }
    delete_LweSample_array(3, temp);

    return result;
}

void shift_left_abit(FheUInt32 &a, const TFheGateBootstrappingCloudKeySet *cloud_key){
    for(int i = 30; i>=0; --i) {
        // copy memory only
        bootsCOPY(a.sample + (i + 1), a.sample + i, cloud_key);
    }
    bootsCONSTANT(a.sample, 0, cloud_key);
}

FheUInt32 div(
    FheUInt32 &a, FheUInt32 &b, // avoid copying of object 
    const TFheGateBootstrappingCloudKeySet *cloud_key){
    const int nb_bits = 32;
    const LweParams *in_out_params = cloud_key->params->in_out_params;

    FheUInt32 result(in_out_params), remainder(in_out_params);
    LweSample *carry = new_LweSample(in_out_params), *not_carry = new_LweSample(in_out_params);

    // x = yq + r
    LweSample *q = result.sample;
    LweSample *r = remainder.sample;
    LweSample *x = a.sample, *y = b.sample;

    for (int32_t j = 0; j < nb_bits; ++j) { 
        bootsCONSTANT(r + j, 0, cloud_key);
    }

    // use binary shifting stragegy
    for (int32_t i = nb_bits-1; i >= 0; --i) { 
        shift_left_abit(remainder, cloud_key);
        bootsCOPY(remainder.sample, x + i, cloud_key);
        FheUInt32 temp = sub_out_carry(remainder, b, cloud_key, carry);
        shift_left_abit(result, cloud_key);
        bootsNOT(not_carry, carry, cloud_key);
        bootsCOPY(result.sample, not_carry, cloud_key);
        // reassign if carry is 0
        for(int32_t i = nb_bits-1; i >= 0; --i) {
            bootsMUX(remainder.sample + i, carry, remainder.sample + i, temp.sample + i, cloud_key);
        }
    }
    delete_LweSample(carry);
    delete_LweSample(not_carry);
    return result;
}

FheUInt32 mod(
    FheUInt32 &a, FheUInt32 &b, // avoid copying of object 
    const TFheGateBootstrappingCloudKeySet *cloud_key){
    const int nb_bits = 32;
    const LweParams *in_out_params = cloud_key->params->in_out_params;

    FheUInt32 result(in_out_params), remainder(in_out_params);
    LweSample *carry = new_LweSample(in_out_params), *not_carry = new_LweSample(in_out_params);

    // x = yq + r
    LweSample *q = result.sample;
    LweSample *r = remainder.sample;
    LweSample *x = a.sample, *y = b.sample;

    for (int32_t j = 0; j < nb_bits; ++j) { 
        bootsCONSTANT(r + j, 0, cloud_key);
    }

    // use binary shifting stragegy
    for (int32_t i = nb_bits-1; i >= 0; --i) { 
        shift_left_abit(remainder, cloud_key);
        bootsCOPY(remainder.sample, x + i, cloud_key);
        FheUInt32 temp = sub_out_carry(remainder, b, cloud_key, carry);
        shift_left_abit(result, cloud_key);
        bootsNOT(not_carry, carry, cloud_key);
        bootsCOPY(result.sample, not_carry, cloud_key);
        // reassign if carry is 0
        for(int32_t i = nb_bits-1; i >= 0; --i) {
            bootsMUX(remainder.sample + i, carry, remainder.sample + i, temp.sample + i, cloud_key);
        }
    }
    delete_LweSample(carry);
    delete_LweSample(not_carry);
    return remainder;
}

int main(){
    TFheGateBootstrappingParameterSet *parameneterSet = new_default_gate_bootstrapping_parameters(16);
    const LweParams* in_out_param = parameneterSet->in_out_params;
    TFheGateBootstrappingSecretKeySet *keySet = new_random_gate_bootstrapping_secret_keyset(parameneterSet);
    // cout << in_out_param->n;
    // client decrypt data with its own key
    FheUInt32 
            x = fromUInt32(7, keySet), 
            y = fromUInt32(3, keySet);
    // cout << "Start mul 2 encrypted number:" << endl;
    auto marked = time(0);
    // FheUInt32 x(in_out_param), y(in_out_param);
    cout << "Start division" << endl;
    auto z = mod(x, y, &keySet->cloud);
    cout << toUInt32(z, keySet) << endl;
    cout << "Time ellapsed:" << time(0)-marked << " second(s)";
}