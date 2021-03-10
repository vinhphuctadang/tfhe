#include <iostream>

using namespace std;

#include <ctime>
#include <sys/time.h>

time_t now(){
    struct timeval _now;
    gettimeofday(&_now, 0);
    return _now.tv_sec * 1000 + _now.tv_usec / 1000;
}

int main(){
    time_t marked = now();
    for(int i = 0; i<1E7; ++i) {
        
    }
    cout << now() - marked;
}