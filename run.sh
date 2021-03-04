if [[ $OSTYPE == "darwin"*  ]]; then 
    EXT="dylib"
else
    EXT="so"
fi
# run the script including the tfhe lib
if g++ ${1} /usr/local/lib/*.${EXT} -o tmp; then
    ./tmp
    rm tmp
fi
