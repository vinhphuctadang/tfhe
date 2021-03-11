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

int main(){
    TFheGateBootstrappingParameterSet *parameneterSet = new_default_gate_bootstrapping_parameters(128);
    const LweParams* in_out_param = parameneterSet->in_out_params;
    TFheGateBootstrappingSecretKeySet *keyset = new_random_gate_bootstrapping_secret_keyset(parameneterSet);
    
    // export secret keys
    // export cloud key
    FILE *cloud_key_file = fopen("cloudkey.fhe", "w");
    export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key_file, &keyset->cloud);
    fclose(cloud_key_file);

    // export secret key;
    FILE *secret_key_file = fopen("secretkey.fhe", "w");
    export_tfheGateBootstrappingSecretKeySet_toFile(secret_key_file, keyset);
    fclose(secret_key_file);
}