all:
	gcc ./unsym_r/unsym_r.c     -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu       -o ./unsym_r/unsym_r     -lpthread -lm -ldl
	gcc ./unsym_r/unsym_r_l.c   -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_l     -o ./unsym_r/unsym_r_l   -lpthread -lm -ldl
	gcc ./unsym_c/unsym_c.c     -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_c     -o ./unsym_c/unsym_c     -lpthread -lm -ldl
	gcc ./unsym_c/unsym_c_l.c   -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_cl    -o ./unsym_c/unsym_c_l   -lpthread -lm -ldl
	gcc ./spd_r/spd_r.c         -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu       -o ./spd_r/spd_r         -lpthread -lm -ldl
	gcc ./spd_r/spd_r_l.c       -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_l     -o ./spd_r/spd_r_l       -lpthread -lm -ldl
	gcc ./spd_c/spd_c.c         -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_c     -o ./spd_c/spd_c         -lpthread -lm -ldl
	gcc ./spd_c/spd_c_l.c       -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_cl    -o ./spd_c/spd_c_l       -lpthread -lm -ldl
	gcc ./sym_r/sym_r.c         -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu       -o ./sym_r/sym_r         -lpthread -lm -ldl
	gcc ./sym_r/sym_r_l.c       -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_l     -o ./sym_r/sym_r_l       -lpthread -lm -ldl
	gcc ./sym_c/sym_c.c         -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_c     -o ./sym_c/sym_c         -lpthread -lm -ldl
	gcc ./sym_c/sym_c_l.c       -O3 -std=gnu11 -I ../include -L ../centos7_x64_mkl -lhylu_cl    -o ./sym_c/sym_c_l       -lpthread -lm -ldl

clean:
	rm ./unsym_r/unsym_r ./unsym_r/unsym_r_l ./unsym_c/unsym_c ./unsym_c/unsym_c_l ./spd_r/spd_r ./spd_r/spd_r_l ./spd_c/spd_c ./spd_c/spd_c_l ./sym_r/sym_r ./sym_r/sym_r_l ./sym_c/sym_c ./sym_c/sym_c_l