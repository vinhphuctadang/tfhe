/*
 * fhe interpreter
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include "tfhe/tfhe.h"
#include "tfhe/polynomials.h"
#include "tfhe/lwesamples.h"
#include "tfhe/lwekey.h"
#include "tfhe/lweparams.h"
#include "tfhe/tlwe.h"
#include "tfhe/tgsw.h"
using namespace std;

const int NB_BITS = 8;
// global cloud key
const TFheGateBootstrappingCloudKeySet *cloud_key;
// global carry flag
const LweSample* carry;

// define n bit integer number manipulation
struct FheUInt {
    LweSample *sample;
    const LweParams *params;
    FheUInt(const LweParams *params){
        this->params = params; // param is READ-ONLY
        this->sample = new_LweSample_array(NB_BITS, params); // NB_BITS bit
    }

    FheUInt(const FheUInt &x) { // only allow copying from the same params
        // deep clone 
        this->params = x.params; // param is READ-ONLY
        for(int i = 0; i<NB_BITS; ++i) {
            lweCopy(this->sample + i, x.sample + i, this->params);
        }
    }
    ~FheUInt() {
        delete_LweSample_array(NB_BITS, this->sample);
    }

    FheUInt& operator=(const FheUInt &x) {
        // deep clone 
        this->params = x.params; // param is READ-ONLY
        for(int i = 0; i<NB_BITS; ++i) {
            lweCopy(this->sample + i, x.sample + i, this->params);
        }
        return *this;
    }
};

void add(FheUInt &a, FheUInt &b){
    const int nb_bits = NB_BITS;
    const LweParams *in_out_params = cloud_key->params->in_out_params;
    // FheUInt result(in_out_params);
    LweSample *carry = new_LweSample(in_out_params);
    bootsCONSTANT(carry, 0, cloud_key);

    LweSample *temp = new_LweSample_array(2, in_out_params);
    LweSample *sum = a.sample;
    LweSample *x = a.sample, *y = b.sample;

    for (int32_t i = 0; i < nb_bits; ++i) {
        bootsXOR(temp, x + i, y + i, cloud_key); // temp = xi XOR yi
        bootsAND(temp + 1, x + i, y + i, cloud_key); // temp1 = xi AND yi
        bootsXOR(x + i, temp, carry, cloud_key); // sum = temp XOR carry = (xi XOR yi) XOR carry
        bootsAND(carry, temp, carry, cloud_key); // carry = (xi XOR yi) and carry
        bootsOR(carry, carry, temp + 1, cloud_key); // carry = (((xi XOR yi) and carry) or (xi AND yi))
    }

    delete_LweSample_array(3, temp);
}

void mov(FheUInt &a, int x){
    for(int i = 0; i<NB_BITS; ++i) {
        bootsCONSTANT(a.sample + i, (x >> i) & 1, cloud_key);
    }
}

void mov(FheUInt &a, FheUInt &b){
    a = b;
}

int main(int n, char* argv[]){
    if (n == 1) {
        return cout << "No code to run" << endl, -1;
    }
    const char* scriptName = argv[0];
    ifstream file(scriptName);

    // File reader
    string line;
    while (getline(file, line)) {
        cout << line << endl;
    }
    // preprocessor

    // executor
    file.close();
}