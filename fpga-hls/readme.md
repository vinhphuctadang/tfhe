# C synthesis problem:

- pointer to pointer is not supported:
Description: 
    Pointer to a struct which has a pointer field inside
    Example:
        StructX* x;
        struct StructX {
            int *a,
            int *b
        }

Bottleneck part:
``
    tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
```

Paramenter types:

```
result: LweSample
    Torus32* a; //-- the n coefs of the mask
    Torus32 b;  //
   	double current_variance; //-- average noise of the sample

bk: TFheGateBootstrappingCloudKeySet:
    const TFheGateBootstrappingParameterSet *const params;
    const LweBootstrappingKey *const bk;
    const LweBootstrappingKeyFFT *const bkFFT;

    LweBootstrappingKeyFFT:
        const LweParams* in_out_params; 
        const TGswParams* bk_params; 
        const TLweParams* accum_params; 
        const LweParams* extract_params;
        const TGswSampleFFT* bkFFT; 
        const LweKeySwitchKey* ks;

        LweParams:
            LweParams:
            // no inner pointer field
            const int32_t n;
            const double alpha_min;
            const double alpha_max;

        TGswParams:
            const int32_t l; 
            const int32_t Bgbit;
            const int32_t Bg;// (must be a power of 2)
            const int32_t halfBg; 
            const uint32_t maskMod; 
            const TLweParams *tlwe_params; // Params of each row
            const int32_t kpl; 
            Torus32 *h; (Size ?)
            uint32_t offset; // offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))

        TLweParams:
            const int32_t N; /// a power of 2: degree of the polynomials
            const int32_t k; // number of polynomials in the mask
            const double alpha_min; // minimal noise s.t. the sample is secure
            const double alpha_max; // maximal noise s.t. we can decrypt
            const LweParams extracted_lweparams; // lwe params if one extracts

        TGswSampleFFT:
            TLweSampleFFT *all_samples;
            TLweSampleFFT **sample;
            //double current_variance;
            const int32_t k;
            const int32_t l;

        LweKeySwitchKey:
            int32_t n; 
            int32_t t; 
            int32_t basebit; // == ceil(log_2(base))
            int32_t base; // decomposition base: power of 2 
            const LweParams* out_params;
            LweSample* ks0_raw; 
            LweSample** ks1_raw;// de taille nl  pointe vers un tableau ks0_raw dont les cases sont espaceés de base positions
            LweSample*** ks; (3 level pointer, represents 3D array)

MU: Torus32:
temp_result: LweSample

```


80-bit gate:
```
    static const int32_t N = 1024;
    static const int32_t k = 1;
    static const int32_t n = 500;
    > static const int32_t bk_l = 2;
    > static const int32_t bk_Bgbit = 10;
    static const int32_t ks_basebit = 2;
    static const int32_t ks_length = 8;
    static const double ks_stdev = 2.44e-5; //standard deviation
    static const double bk_stdev = 7.18e-9; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space


    LweParams *params_in = new_LweParams(n, ks_stdev, max_stdev);

    --> (easy)
        n(n),
		alpha_min(alpha_min),
		alpha_max(alpha_max) {}

    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);

    --> N(N),
        k(k),
        alpha_min(alpha_min),
        alpha_max(alpha_max),
        extracted_lweparams(N * k, alpha_min, alpha_max) {} // LweParams constructor (not TLwe)

    TGswParams *params_bk = new_TGswParams(bk_l, bk_Bgbit, params_accum);

TGswParams(int32_t l, int32_t Bgbit, const TLweParams *tlwe_params) :
        l(l),
        Bgbit(Bgbit),
        Bg(1 << Bgbit),
        halfBg(Bg / 2),
        maskMod(Bg - 1),
        tlwe_params(tlwe_params),
        kpl(int32_t((tlwe_params->k + 1) * l)) {

    h = new Torus32[l];
    for (int32_t i = 0; i < l; ++i) {
        int32_t kk = (32 - (i + 1) * Bgbit);
        h[i] = 1 << kk; // 1/(Bg^(i+1)) as a Torus32
    }

    // offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))
    uint32_t temp1 = 0;
    for (int32_t i = 0; i < l; ++i) {
        uint32_t temp0 = 1 << (32 - (i + 1) * Bgbit);
        temp1 += temp0;
    }
    offset = temp1 * halfBg;

    }

    ks_t(ks_t),
    ks_basebit(ks_basebit),
    in_out_params(in_out_params),
    tgsw_params(tgsw_params)
```


new_TGswParams(bk_l, bk_Bgbit, params_accum)


static TFheGateBootstrappingParameterSet *default_80bit_gate_bootstrapping_parameters() {
    // These are the historic parameter set provided in 2016,
    // that were analyzed against lattice enumeration attacks
    // Currently (in 2020), the security of these parameters is estimated to lambda = **80-bit security**
    // (w.r.t bkz-sieve model, + hybrid attack model)
    static const int32_t N = 1024;
    static const int32_t k = 1;
    static const int32_t n = 500;
    static const int32_t bk_l = 2;
    static const int32_t bk_Bgbit = 10;
    static const int32_t ks_basebit = 2;
    static const int32_t ks_length = 8;
    static const double ks_stdev = 2.44e-5; //standard deviation
    static const double bk_stdev = 7.18e-9; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space

    LweParams *params_in = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}