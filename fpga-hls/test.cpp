#include <stdio.h>

typedef struct {
	const int a[5];
} Param;

Param assign(Param param) {
	int b[5] = { param.a };
	
	return param;
}

int main(){
	printf("Hello world");
	Param x;
	assign(x);
}
