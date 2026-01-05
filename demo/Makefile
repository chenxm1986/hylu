all:
	gcc demo.c   -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu    -o demo   -lpthread -lm -ldl
	gcc democ.c  -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_c  -o democ  -lpthread -lm -ldl
	gcc demol.c  -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_l  -o demol  -lpthread -lm -ldl
	gcc democl.c -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_cl -o democl -lpthread -lm -ldl

clean:
	rm demo democ demol democl