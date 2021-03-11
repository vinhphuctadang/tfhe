/*
 * fhe interpreter
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <sstream>
#include "tfhe/tfhe.h"
#include "tfhe/polynomials.h"
#include "tfhe/lwesamples.h"
#include "tfhe/lwekey.h"
#include "tfhe/lweparams.h"
#include "tfhe/tlwe.h"
#include "tfhe/tgsw.h"
using namespace std;

const char* CLOUD_KEY_FILE = "cloudkey.fhe";
const int NB_BITS = 8;

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

// global cloud key
TFheGateBootstrappingCloudKeySet *cloud_key;
// global carry flag
LweSample *carry;
// variables
unordered_map<string, FheUInt*> variables;


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

void out(FheUInt &a) {
    for(int i = 0; i<NB_BITS; ++i) {
        export_gate_bootstrapping_ciphertext_toFile(stdout, a.sample+i, cloud_key->params);
    }
}

bool is_digit(const string& s) {
    return !s.empty() && all_of(s.begin(), s.end(), ::isdigit);
}

int to_int(const string& s) {
    stringstream ss;
    ss << s;
    int result;
    ss >> result;
    return result;
}

int main(int n, char* argv[]){
    if (n == 1) {
        return cout << "No code to run" << endl, -1;
    }
    const char* scriptName = argv[1];
    ifstream file(scriptName);

    // import cloud key
    FILE *cloud_key_file = fopen("cloudkey.fhe", "r");
    cloud_key = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key_file);
    fclose(cloud_key_file);

    // File reader
    string line;
    while (getline(file, line)) {
        //cout << line << endl;
        // simple syntax analyse
        
        vector<string> instruction_detail;
        // split string
        int pivot = 0;
        while (pivot < line.length()) {
            int p;
            p = line.find(" ", pivot);
            if (p != string::npos) {

                string operand = line.substr(pivot, p-pivot);

                if (operand[operand.length()-1] == ',') operand.pop_back(); // remove trailing comma
                instruction_detail.push_back(operand);
                pivot = p + 1;
                continue;
            }
    
            // not found
            instruction_detail.push_back(line.substr(pivot, line.length()-pivot));
            break;
        }

        // print out
        // for(auto detail : instruction_detail) {
        //     cout << detail << "|";
        // }

        if (instruction_detail[0] == "mov"){
            string a_name = instruction_detail[1];
            string b = instruction_detail[1];
            
            // mov mem, const
            if (is_digit(b)){
                if (variables.find(a_name) == variables.end()) {
                    variables[a_name] = new FheUInt(cloud_key->params->in_out_params);
                }
                mov(*variables[a_name], to_int(b));
            }
            // mov mem, mem
            else {
                if (variables.find(a_name) == variables.end()) {
                    variables[a_name] = new FheUInt(cloud_key->params->in_out_params);
                }
                mov(*variables[a_name], *variables[b]);
            }
        }
        else if (instruction_detail[0] == "add") {
            string a_name = instruction_detail[1];
            string b = instruction_detail[1];            
            if (variables.find(a_name) == variables.end()) {
                variables[a_name] = new FheUInt(cloud_key->params->in_out_params);
            }
            add(*variables[a_name], *variables[b]);
        }
        else if (instruction_detail[0] == "out") {
            string a_name = instruction_detail[1];
            out(*variables[a_name]);
        }

        // cout << endl;
    }
    file.close();
    // delete_gate_bootstrapping_secret_keyset(cloud_key);
    // TODO: Destruct arrays
}