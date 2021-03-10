if [[ $OSTYPE == "darwin"*  ]]; then 
    EXT="dylib"
else
    EXT="so"
fi
# run the script including the tfhe lib
if g++ --std=c++11 interpreter.cpp /usr/local/lib/*.${EXT} -o interpreter; then
    ./interpreter $@
fi
