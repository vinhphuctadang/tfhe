# run the script including the tfhe lib
if g++ ${1} /usr/local/lib/*.so -o tmp; then
    ./tmp
    rm tmp
fi
