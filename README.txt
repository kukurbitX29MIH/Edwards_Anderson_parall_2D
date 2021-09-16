This programs calculates with parallelization ground state energy of 2D rectangle lattice for Ising/Edwards Anderson model.
It can calculate ground state spin configuration as well.
Distribution of coupling constant is bimodal (J_bond_i= +1 or -1).
But program can be modified easily to treat floating point constant.
The type of boundary condition is FBC (free boundary conditions).
But one can set out boundary spins around the lattice and bonds. 
Program returns as a result ground state energy, up row ground state spin configuration.
Result is given in out file.
Knowing up row ground state spin configuration one can start new calculation with less height lattice with given up boundary spin frame corresponding obtained on the previous calculation up row gs spin state.
It gives a possibility to get ground state spin configuration.

The program must get 3 or 5 arguments of command line.
Argument format:   1) *argc[1]: Ly    2) *argc[2]: Lx    3) *argc[3]: J_bond_plus_perc    4) *argc[4]:  bool__1_data_from_command_line_args__0_from_file    5) *argc[5]:  chstr_filename\n"

There are two options:

  1)  if (bool__1_data_from_command_line_args__0_from_file not equals '1') or (amount of argument == 3):         
  
      In this case J_bonds data will be generated automatically with given (1) Ly, (2) Lx and percent of positive bonds (3) J_bond_plus_perc.
      Outer spin frame will not be (full fbc)           
          
  2)  if (bool__1_data_from_command_line_args__0_from_file     equals '1') or (amount of argument == 5):         
       In this case all aprameters will be taken from file in 5th argument of command line.
       Ly, Lx, J_bond_plus_perc,  outer frame of spins and bonds are taken from file.            
 
 
Format input file has view:

#  Ly
6
#  Lx
6
#  J_plus_perc
50
#  J_min_perc
50
#  down edge spin state
-1 +1 -1 +1 -1 -1
#  up edge spin state
-1 -1 -1 -1 -1 -1
#  left edge spin state
-1 -1 +1 -1 -1 -1
#  right edge spin state
-1 -1 -1 -1 -1 -1
#  down edge spin existing
-1 +1 -1 +1 -1 +1
#  up edge spin existing
+1 -1 -1 -1 -1 -1
#  left edge spin existing
-1 -1 +1 +1 +1 -1
#  right edge spin existing
+1 -1 -1 -1 -1 -1
#  down edge J_ver
-1 -1 -1 -1 +1 -1
#  up edge J_ver
-1 +1 +1 +1 +1 +1
#  left edge J_hor
-1 -1 -1 -1 -1 +1
#  right edge J_hor
-1 +1 -1 -1 -1 +1
#  J table
#  s1_index_y	s1_index_x	s2_index_y	s2_index_x	bond_mean
0		0		1		0		-1
0		1		1		1		+1
0		2		1		2		+1
0		3		1		3		+1
0		4		1		4		+1
0		5		1		5		-1
1		0		2		0		-1
1		1		2		1		-1
1		2		2		2		-1
1		3		2		3		-1
1		4		2		4		-1
1		5		2		5		-1
2		0		3		0		-1
2		1		3		1		+1
2		2		3		2		+1
2		3		3		3		+1
2		4		3		4		+1
2		5		3		5		+1
3		0		4		0		-1
3		1		4		1		-1
3		2		4		2		+1
3		3		4		3		-1
3		4		4		4		+1
3		5		4		5		+1
4		0		5		0		+1
4		1		5		1		+1
4		2		5		2		+1
4		3		5		3		-1
4		4		5		4		+1
4		5		5		5		+1
0		0		0		1		+1
0		1		0		2		-1
0		2		0		3		+1
0		3		0		4		-1
0		4		0		5		+1
1		0		1		1		+1
1		1		1		2		-1
1		2		1		3		-1
1		3		1		4		+1
1		4		1		5		-1
2		0		2		1		+1
2		1		2		2		-1
2		2		2		3		+1
2		3		2		4		+1
2		4		2		5		-1
3		0		3		1		-1
3		1		3		2		-1
3		2		3		3		-1
3		3		3		4		+1
3		4		3		5		+1
4		0		4		1		+1
4		1		4		2		-1
4		2		4		3		+1
4		3		4		4		-1
4		4		4		5		-1
5		0		5		1		-1
5		1		5		2		+1
5		2		5		3		+1
5		3		5		4		+1
5		4		5		5		+1
5		0		6		0		-1
5		1		6		1		+1
5		2		6		2		+1
5		3		6		3		+1
5		4		6		4		+1
5		5		6		5		+1

This data corresponds to

matrix struct:
-------------------------------------------------------------------
|                                                                 |
|            d                                                    |
|                                                                 |
|            -                                                    |
|                                                                 |
|            S   -   S   +   S   +   S   +   S   +   S            |
|                                                                 |
|            +       +       +       -       +       +            |
|                                                                 |
|    d   -   S   +   S   -   S   +   S   -   S   -   S            |
|                                                                 |
|            -       -       +       -       +       +            |
|                                                                 |
|    d   -   S   -   S   -   S   -   S   +   S   +   S            |
|                                                                 |
|            -       +       +       +       +       +            |
|                                                                 |
|    u   -   S   +   S   -   S   +   S   +   S   -   S            |
|                                                                 |
|            -       -       -       -       -       -            |
|                                                                 |
|            S   +   S   -   S   -   S   +   S   -   S            |
|                                                                 |
|            -       +       +       +       +       -            |
|                                                                 |
|            S   +   S   -   S   +   S   -   S   +   S   -   d    |
|                                                                 |
|                    -               -               -            |
|                                                                 |
|                    u               u               d            |
|                                                                 |
-------------------------------------------------------------------

Thus one can understand easily meanings of parameters of input file. We can set existing of outer frame spins and bonds. The name of parameters have clear means.
Data is read on the next line after key phrase excepting J-table.
In J-table bonds are set as value +1 or -1 between spins with coordinate (y1, x1) and (y2, x2).
It gives undertanding of bond structure. Spin's row coordinate numerations starts with 0.
Numeration way of outer frame spins is given below
//
//            up (spins, existing/notexisting)
//        s1  s2  s3
//        |   |   |
// t s3 - *   *   * - s3 t
// f                      h
// e s2 - *   *   * - s2 g (spins, existing/notexisting)
// l                      i
//   s1 - *   *   * - s1 r
//        |   |   |
//        s1  s2  s3 
//          down (spins, existing/notexisting)



Program generates at the end of its work the file with same structure as input. So one can take out file and use it as input in the next launch.
Examples of the input files are given in the main directory.




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  V  V  V
!!
!!
For program compilation one needs to have OS linux and installed libraries:    mpfr    gmp    lgmpxx
Also one needs installed and tuned Open MPI environment.
Examples

how to compile:   mpic++ code.cpp   -o prg.exe -lmpfr -lgmp -lgmpxx

how to launch for 4 threads lattice 6x6, J-bond is generated automatically, without outer spins:    mpirun -n 4 prg.exe 6 6 50 
  or
 ow to launch for 4 threads lattice 6x6, J-bond is generated automatically (50% of positive bonds), without outer spins:    mpirun -n 4 prg.exe 6 6 50 0 Ising_EA_2D_gbc_bond_pl50_min50_0006x0006.txt
  or
how to launch for 4 threads,  all data will be taken from file:    mpirun -n 4 prg.exe 6 6 50 1 aIsing_EA_2D_gbc_bond_pl50_min50_0006x0006.txt
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ^  ^  ^





