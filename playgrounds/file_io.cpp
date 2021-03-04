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
    TFheGateBootstrappingParameterSet* parameneterSet = new_default_gate_bootstrapping_parameters(100);
    const LweParams* in_out_param = parameneterSet->in_out_params;

    // generate secret key based on paramenter set
    TFheGateBootstrappingSecretKeySet* secretKeySet = new_random_gate_bootstrapping_secret_keyset(parameneterSet);

    int bit_count = 4;
    LweSample* a = new_LweSample_array(bit_count, in_out_param);
    // bootsCONSTANT(a, 1, &secretKeySet->cloud);
    // bootsCONSTANT(a+1, 0, &secretKeySet->cloud);
    // bootsCONSTANT(a+2, 0, &secretKeySet->cloud);
    // bootsCONSTANT(a+3, 1, &secretKeySet->cloud);
    int plain = 1101;
    cout << "Plain:";
    for(int i = 0; i<bit_count; ++i) {
        bootsSymEncrypt(a+i, plain % 10, secretKeySet);
        cout << plain % 10;
        plain /= 10;
    }
    cout << endl;
    // output cipher
    FILE* f = fopen("cipher.txt", "w");
    for(int i = 0; i<bit_count; ++i) {
        export_gate_bootstrapping_ciphertext_toFile(f, a+i, parameneterSet);
    }
    fclose(f);

    // output secret
    f = fopen("secret.txt", "w");
    export_tfheGateBootstrappingSecretKeySet_toFile(f, secretKeySet);
    fclose(f);

    LweSample *b = new_LweSample_array(bit_count, in_out_param); // b = a, but read from file
    f = fopen("cipher.txt", "r");
    for(int i = 0; i<bit_count; ++i) {
        import_gate_bootstrapping_ciphertext_fromFile(f, b+i, parameneterSet);
    }
    fclose(f);

    cout << "Decrypted received text:";
    for(int i = 0; i<bit_count; ++i) {
        int bit = bootsSymDecrypt(b+i, secretKeySet);
        cout << bit;
    }
}