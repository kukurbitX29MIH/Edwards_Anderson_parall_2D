// Last modified 16.09.2021 PadalkoMA
// how to compile:   mpic++ code.cpp   -o prg.exe -lmpfr -lgmp -lgmpxx
// how to launch:    mpirun -n 4 prg.exe 6 6 50 1 aIsing_EA_2D_gbc_bond_pl50_min50_0008x0008.txt
//  or
// how to launch:    mpirun -n 4 prg.exe 6 6 50 0 aIsing_EA_2D_gbc_bond_pl50_min50_0008x0008.txt
//  or
// how to launch:    mpirun -n 4 prg.exe 6 6 50 1 


#include <iostream>
#include <cstdlib>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <random>
#include <string.h>
#include <gmp.h>
#include <gmpxx.h>
#include <mpi.h>
//#include <mpfr.h>
#include "anyfunctions.cpp"
#include "rdk_io_formatting.cpp"
#include "rdk_io_formatting_ext.cpp"
#include "rdk_io_lattice_models.cpp"



#define __MACR_GS_THREAD_NUM_TO_COMMENT 0
#define __MACR_PRINT_PROGRESS 1
#define __MACR_TO_EXTRACT_DATA_FROM_FILE 1      
#define __MACR_TO_DO_RECORD_INTO_FILE 1




#define __MACR_MATR_TYPE_1__MPZ
//#define __MACR_CHECK_ACCESS_MEM_VIOLATION
#define __MACR_LATTICE_LY 8
#define __MACR_LATTICE_LX 8

#define __MACR_MPZ_MAX_CHSTR_SIZE 30
#define __MACR_MPZ_BIT_SIZE 100

#define __MACR_SPIN_CONFIG_AR_AMOUNT_MAX_PER_ROW_STATE 256


#define __MACR_GS_ANOUNT_MAX_L_09 725336
#define __MACR_GS_ANOUNT_MAX_L_10 24103424
#define __MACR_GS_ANOUNT_MAX_L_11 134627770
#define __MACR_GS_ANOUNT_MAX_L_12 134627770  //
#define __MACR_GS_ANOUNT_MAX_L_13 134627770  //
#define __MACR_GS_ANOUNT_MAX_L_14 134627770  //
#define __MACR_GS_ANOUNT_MAX_L_15 134627770  //
#define __MACR_GS_ANOUNT_MAX_L_16 134627770  //
#define __MACR_GS_ANOUNT_MAX_L_17 134627770  //





// 3 x 3 takes __MACR_MPZ_BIT_SIZE=9    bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=3   using memory for storing of matrices is about:    0 Mb +   0 kb + 702 b      time T=        0 ms   ratio= 
// 4 x 4 takes __MACR_MPZ_BIT_SIZE=16   bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=5   using memory for storing of matrices is about:    0 Mb +   7 kb + 832 b      time T=        2 ms   ratio= 4
// 5 x 5 takes __MACR_MPZ_BIT_SIZE=25   bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=8   using memory for storing of matrices is about:    0 Mb +  56 kb +  56 b      time T=        8 ms   ratio= 5.5
// 6 x 6 takes __MACR_MPZ_BIT_SIZE=36   bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=11  using memory for storing of matrices is about:    0 Mb + 343 kb + 128 b      time T=       44 ms   ratio= 4.11
// 7 x 7 takes __MACR_MPZ_BIT_SIZE=49   bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=15  using memory for storing of matrices is about:    1 Mb + 668 kb +  32 b      time T=      181 ms   ratio= 4.21
// 8 x 8 takes __MACR_MPZ_BIT_SIZE=64   bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=20  using memory for storing of matrices is about:    7 Mb + 516 kb +   0 b      time T=      763 ms   ratio= 3.65

// 9 x 9 takes __MACR_MPZ_BIT_SIZE=81   bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=25  using memory for storing of matrices is about:   30 Mb + 110 kb + 640 b      time T=    2 782 ms   ratio= 3.83
// 10x10 takes __MACR_MPZ_BIT_SIZE=100  bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=30  using memory for storing of matrices is about:  114 Mb + 914 kb +   0 b      time T=   10 643 ms   ratio= 2,90
// 11x11 takes __MACR_MPZ_BIT_SIZE=121  bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=37  using memory for storing of matrices is about:  404 Mb + 789 kb + 512 b      time T=   38 952 ms   ratio=
// 12x12 takes __MACR_MPZ_BIT_SIZE=144  bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=44  using memory for storing of matrices is about: 1378 Mb + 848 kb +   0 b      time T=  133 513 ms   ratio=
// 13x13 takes __MACR_MPZ_BIT_SIZE=169  bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=51  using memory for storing of matrices is about: 5152 Mb + 280 kb +   0 b      time T=   ms   ratio=  
// 14x14 takes __MACR_MPZ_BIT_SIZE=196  bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=60
// 15x15 takes __MACR_MPZ_BIT_SIZE=225  bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=68
// 16x16 takes __MACR_MPZ_BIT_SIZE=256  bits on 1 element of hist-matrix  __MACR_MPZ_MAX_CHSTR_SIZE=78



//                    ___    ___________________________
//                   \    ________\_______   _______   /
//        ________    ______________________________/ /
//       \  ____________\_________   ____________/                                              
//            _______________   _______          
//     _____________________\___  ___________________\             
//                        
//                      





//  
//   +--^----------,--------,-----,--------^-,
//   | |||||||||   `--------'     |          O
//   `+---------------------------^----------|
//     `\_,---------,---------,--------------'
//       / XXXXXX /'|       /'
//      / XXXXXX /  `\    /'
//     / XXXXXX /`-------'
//    / XXXXXX /
//   / XXXXXX /
//  (________(                
//   `------' 

using namespace std;
using namespace _NMSP_RDK_;

MPI_Status Status;                                 // ====================  MPI
int mpisize, mpirank;                              // ====================  MPI





// *************************************************************** (1) in/out/repres
// ---------------------------------------------------------------------------------

//  display aerror message
void display_on_screen_error_msg(char * chstr_error_text); 

inline void print_progress(unsigned long long int progress, double & mustachieve_progress_perc, unsigned long long int max_progress, ofstream &out);

inline void print_progress(unsigned long long int progress, unsigned long long int max_progress, ofstream &out);



//  write chstr text to file
//--// int write_chstr_to_file(char * text, const char * const filename, bool B1New0Add=true); 

//  concatenation
//--// void add_int_to_string(string * string_ptr, int num);

//  concatenation
//--// void add_uint_to_string(string * string_ptr, unsigned int num);

//  concatenation
//--// void add_llint_to_string(string * string_ptr, long long int num);

//  concatenation
//--// void add_ullint_to_string(string * string_ptr, unsigned long long int num);

//  concatenation
//--// void add_double_to_string(string * string_ptr, double num);

//  calculate amount of bit1 in num
inline int amount_of_bit1(unsigned long long int num);

//  calculate amount of pairs with dif means of bit.   pair is only 2 bit-neighbours within of len.   last and first bit is connected in a cyclic way.
inline int amount_of_difdir_pairs_pbc(unsigned long long int row_state, int len);

//  calculate amount of pairs with dif means of bit.   pair is only 2 bit-neighbours within of len.   last and first bit is connected in a cyclic way.
inline int amount_of_difdir_pairs_fbc(unsigned long long int row_state, int len);

//    calculate amount of pairs with dif means of bit in opposite located spins in 2 rows.   pair is only 2 bit-neighbours within of len. 
inline int amount_of_difdir_pairs_betw_2_rows(unsigned long long int row1_state1, unsigned long long int row2_state2, int len);



// J treatment procedures
//
bool generate_2D_lattice_J_v1(bool * J_ver,  bool * J_hor,  int Ly,  int Lx, int plus_count, int negat_count, bool pbc_1_fbs_0=true);

bool generate_2D_lattice_J_v3(bool * J_ver,  bool * J_hor,  int Ly,  int Lx, int plus_count, int negat_count);

bool invert_J(bool * J_ver,  bool * J_hor,  int Ly,  int Lx, bool pbc_1_fbs_0=true);

inline void inver_J(bool *J_ver,  bool *J_hor,  int Ly,  int Lx);

inline void read_ullint_J_ver_line(unsigned long long int & J_ver_line,  bool * J_ver_array,  int ver_line_num, int amount_of_bond_read, int Lx);

inline void read_ullint_J_hor_line(unsigned long long int & J_hor_line,  bool * J_hor_array,  int ver_line_num, int amount_of_bond_read, int Lx);

inline int amount_of_bit1_logistic_xor_3_nums(unsigned long long int J_line,  unsigned long long int num1,  unsigned long long int num2, int len);

inline int amount_of_bit1_logistic_xor_1_num_with_shift_pbc(unsigned long long int J_line,  unsigned long long int spin_line_state, int len);

inline int amount_of_bit1_logistic_xor_1_num_with_shift_fbc(unsigned long long int J_line,  unsigned long long int spin_line_state, int len);


//  neg pair if e.g.:    s1=d  s2=d  J=-1   (fer)           ---=-                      pos pair if e.g.:    s1=d  s2=u  J=-1   (fer)           -+-=+
//
void calc_pospair_and_upspin_EA_2D_lattice_PBC(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount=0,  int * up_amount=0);

void calc_pospair_and_upspin_EA_2D_lattice_FBC(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount=0,  int * up_amount=0);

void calc_pospair_and_upspin_EA_2D_lattice_PBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount=0,  int * up_amount=0);

void calc_pospair_and_upspin_EA_2D_lattice_FBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount=0,  int * up_amount=0);

void calc_pospair_and_upspin_EA_2D_lattice_with_edges_without_down_edge_v3(bool * spin_array,  bool * bool_edge_spin_exist_d_u_l_r_ar, bool * J_ver, bool * J_hor, int Ly, int Lx, int * pospair_amount=0, int * up_amount=0, bool calc_edge_down=false);

void calc_pospair_and_upspin_EA_2D_lattice_cylind_PBC_FBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount=0,  int * up_amount=0);

void calc_pospair_and_upspin_EA_2D_lattice_cylind_hor_PBC_ver_FBC_v3(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount=0,  int * up_amount=0);

inline int amount_of_positiv_E_par_revers_order_edge_v3(bool *J_line,  bool * bool_edge_spin_exist_d_u_l_r_ar, unsigned long long int spin_line_state, int len);

inline int positi_E_of_1_par(bool spin_1,  bool spin_2, bool J);

//  calc will be started with 0-index
inline int amount_of_positiv_E_par_edge_v3(bool *J_line,  bool * bool_edge_spin_exist_d_u_l_r_ar, unsigned long long int spin_line_state, int len);

//  calc will be started with 0-index,         J_line[0] <--> bool_edge_spin_exist_d_u_l_r_ar[2*0] <--> spin_line_state[len-i-1 bit]
inline int amount_of_positiv_E_par_revers_order_edge_v3(bool *J_line,  bool * bool_edge_spin_exist_d_u_l_r_ar, unsigned long long int spin_line_state, int len);

inline int amount_of_positiv_E_par_bool_edge_v3(bool *J_line,  bool * bool_edge_spin_exist_d_u_l_r_ar,  bool * bool_spin_array, int len);




//  -----------------------
//  
//           a a
//         a * * a
//         a * * a
//           a a   
//            
bool calc_pospair_and_upspin_amount_EA_2D_lat_interact_with_bound(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  unsigned long long int *boundstate,  int * energy=0,  int * magnetization=0); 

void calc_energy_and_magnetiz_EA_2D_lat_PBC(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * energy=0,  int * magnetization=0);

void calc_energy_and_magnetiz_EA_2D_lat_FBC(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * energy=0,  int * magnetization=0);

void calc_energy_and_magnetiz_EA_2D_lat_PBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * energy=0,  int * magnetization=0);

void calc_energy_and_magnetiz_EA_2D_lat_FBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * energy=0,  int * magnetization=0);

int calc_magnetiz_lat(bool * spin_array,  int Ly,  int Lx);




//  s[1]=   -1,  0(not spin),  +1
//
//               J[1]
//               s[1]
//   s[2] J[2]   up1dn  s[0] J[0]
//               s[3]
//               J[3]
inline int pospair_amount_1_cell_EA(bool up1dn0 , bool *J, char *s);      //  !  !  !




// ============= (( (( < gmp

//  stat variant for chstr_num,      convert mpz_t to chstr var:  necessary len will be set with symb ch_empt_symb       desired_len - len we want 
//            if (desired_len >  0) {it will be:  filling with ch_empt_symb of rest pos;}
//            if (desired_len == 0) {it will be: desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}
//            if (desired_len <  0) {it will be: desired_len=real_len;}  
int mpz_t__to__chstr__stat(mpz_t * mpz_t_num,   char * const chstr_num,  int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);

//  dyn variant for chstr_num,      convert mpz_t to chstr var:  necessary len will be set with symb ch_empt_symb       desired_len - len we want 
//            if (desired_len >  0) {it will be:  filling with ch_empt_symb of rest pos;}
//            if (desired_len == 0) {it will be: desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}
//            if (desired_len <  0) {it will be: desired_len=real_len;}  
int mpz_t__to__chstr__dyn(mpz_t * mpz_t_num,   char * & chstr_num,  int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);

//  convert mpz_t to screen:  necessary len will be set with symb ch_empt_symb    (1) desired_len - len we want, if (desired_len <= 0) {it will be: desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}  
//            if (desired_len >  0) {it will be:  filling with ch_empt_symb of rest pos;}
//            if (desired_len == 0) {it will be: desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}
//            if (desired_len <  0) {it will be: desired_len=real_len;}  
void mpz_t__to__screen(mpz_t * mpz_t_num,   int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);



// *************************************************************** (2) additional
// ------------------------------------------------------------------------------



unsigned long long int ullint_2_pow_n(int n);
inline unsigned int bit_inverse_int(unsigned int num,  int len);
inline unsigned long long int bit_inverse_ullint(unsigned long long int num,  int len);
inline void write_bits_uint(unsigned int & num_dest, int start_pos_dest,   unsigned int num_src, int start_pos_src,   int len);
inline void write_bits_ullint(unsigned long long int & num_dest, int start_pos_dest,   unsigned long long int num_src, int start_pos_src,   int len);

void fill_random_array_with_0_or_1(bool * const array,  int N_loc);

inline int cycl_correction(int mean, int period); 
void shift_J_bonds(bool *J_ver,  bool *J_hor,  int shift_dy,  int shift_dx,  int Ly,  int Lx);


struct matr_t1
{

matr_t1();
void free();
~matr_t1() {free();}  
bool init(int size_y_loc, int size_x_loc, unsigned int seting_mean);

void update_min_max_for_adding_elem(int  index_y,  int index_x);
void set_val_for_all_el(unsigned int value);
void set_val_for_one_el(int  index_y,  int index_x,  unsigned int value);
void assign_val_from_one_el_to_other(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src);
void set_zero_val_for_all_el();
void add_val_from_one_el_to_other(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src); 
void mul_matr_by_num(unsigned int mult);
void swap(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src);
void swap_with_add(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src);

bool shift_elements(int dy, int dx);
inline bool shift_elements_spec_fast_case(int dy);
inline bool shift_elements_spec_fast_case_M_div_2(int dy);
inline bool shift_elements_spec_fast_case_M_div_2_v20(int dy);
inline bool shift_elements_spec_without_M_v22(int dy);
unsigned long long int get_element_as_ullint_if_possible(int  index_y,  int index_x);
int get_element_as_chstr_stat(int  index_y,  int index_x,  char * const chstr_matr_element,     int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);      //  //  ? ? ?
int get_element_as_chstr_dyn(int  index_y,  int index_x,  char * & chstr_matr_element,     int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);      //  //  ? ? ?
void display_on_screen_one_element(int  index_y,  int index_x,      int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);
void display_on_screen_not_0_area_size_param();
void display_on_screen_whole_matr(char * interval_betw_elem=0,  int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);
int  save_to_file_whole_matr(const char * const chstr_filename,  char * interval_betw_elem=0,  int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);


#ifdef __MACR_MATR_TYPE_1__MPZ
mpz_t * element;
#endif

int size_y,  size_x, size_N;
int min_not_0_x_index, max_not_0_x_index,  not_0_len_x;
int min_not_0_y_index, max_not_0_y_index,  not_0_len_y;

}; 


//  Sum    elements matr x1 (len_y * len_x)   and    elements matr x2 (len_y * len_x).    Put result matr to result.
void mpz_t_matr_sum(mpz_t *result,  mpz_t *x1,  mpz_t *x2,  int len_y,  int len_x);

//  Dif    elements matr x1 (len_y * len_x)   and    elements matr x2 (len_y * len_x).    Put result matr to result.
void mpz_t_matr_dif(mpz_t *result,  mpz_t *x1,  mpz_t *x2,  int len_y,  int len_x);

//  Shift    elements of matr (len_y * len_x) on shift_y and shift_x.    
void mpz_t_matr_shift(mpz_t *matr,  int len_y,  int len_x,  int shift_y,  int shift_x);


//  Add elements of matr_add_src to matr_dest.
void matr1_add(matr_t1 *matr_dest,  matr_t1 *matr_add_src);

// copy matrix
inline void matr_t1_copy(matr_t1 *matr_dest,  matr_t1 *matr_src);

// copy matrix with shift
inline void matr_t1_copy_with_shift(matr_t1 *matr_dest,  matr_t1 *matr_src, int shf_y, int shf_x);

//  Add elements of matr_add_src to matr_dest.
void matr1_add_dif_sizes(matr_t1 *matr_dest,  matr_t1 *matr_add_src);

//  add given area of elements of matr_src with preliminary shift  to  matr_dest
inline void matr1_area_is_added_to_matr2(matr_t1 *matr_dest, int matr_dest_start_index,    matr_t1 *matr_src, int matr_src_start_index,    int matr_src_subarea_len_y, int matr_src_subarea_len_x);

//  assign given area of elements of matr_src with preliminary shift  to  matr_dest
inline void matr1_area_is_assign_to_matr2(matr_t1 * matr_dest, int matr_dest_start_index,    matr_t1 *matr_src, int matr_src_start_index,    int matr_src_subarea_len_y, int matr_src_subarea_len_x);





//  ! ! ! not-0 elements are not modified of matr ! ! !
//                         
//                   
// 
//          F F * * *          * * * * *
//          F F * * *          * * * * *
//    matr1 F @ A A *          * * * * *  matr2
//    src   * A A A *          * ^ * * *  dest
//          * | * * *          * | * * *
//            |                  |
//            |__________________|
//                        add
//
//   A - coping area of matr1-src
//   F - forbidden  for coping area of matr1-src
//   ! ! ! not-0 elements are not modified of matr ! ! !
//
inline void matr1_area_spec_add_to_matr2(matr_t1 * matr_dest,  int matr_dest_start_pos,                matr_t1 *matr_src,  int matr_src_subarea_y0,  int matr_src_subarea_x0,  int matr_src_subarea_len_y,  int matr_src_subarea_len_x,  int matr_src_forbidsubarea_y0,  int matr_src_forbidsubarea_x0,  int matr_src_forbidsubarea_len_y,  int matr_src_forbidsubarea_len_x);

//
//                         add
//                   _________________
//                  |                |
//          * * * * V          * * * | *
//          * * * * *          * * * B *
//    matr1 * * * * *          * * * * *  matr2
//          * A * * *          * ^ * * *
//          * | * * *          * | * * *
//            | elem1            |
//            |__________________|
//                        add
//
//  we have 4 points: 2 src,  2 dest (where add to)
//  ! ! ! key problem is that any points can coincide ! ! !
//  it's the main reason why this proc is needed
//  ! ! ! not-0 elements are not modified of matr ! ! !
//
inline void two_matrices_mutual_add_with_1_elem(matr_t1 *matr1, int elem1_index_from, int elem1_index_to,    matr_t1 *matr2, int elem2_index_from, int elem2_index_to);


//  swap with adding 2 the same quad chunks between 2 matrices 
//
//  * * * * *     * * * * *
//  * B B * *     * * * * *
//  * B B * *     * * A A *
//  * * * * *     * * A A *
//  * * * * *     * * * * *
//    ^ | ____add_____^ |
//    |_______add_______|
//
//  add assign crd are:      int y0_to_matr2, int x0_to_matr2                 int y0_to_matr1, int x0_to_matr1
//  src crd of chunk:        int y0_from_matr1, int x0_from_matr1             int y0_from_matr2, int x0_from_matr2
//  subarea sizes are the same
//  ! ! ! base 0-point is left down point ! ! !
//  ! ! ! matr1_index_from=matr2_index_to ! ! !
//  ! ! ! not-0 elements are not modified of matr ! ! !
//
//  or other scheme
//  matr1                    matr2
//  * * * * * *              * * * * * *
//        |____________________^   |
//        ^________________________|
//        int y0_matr2_to 
//        int x0_matr2_to
//
inline void matr1_area_spec_mutual_add_swap(matr_t1 * matr_1, int y0_matr1_to, int x0_matr1_to,     matr_t1 * matr_2, int y0_matr2_from, int x0_matr2_from, int y0_matr2_to, int x0_matr2_to,  int subarea_len_y, int subarea_len_x);
inline void matr1_area_spec_mutual_add_swap_v22(matr_t1 * matr_1, int y0_matr1_to, int x0_matr1_to,     matr_t1 * matr_2, int y0_matr2_from, int x0_matr2_from, int y0_matr2_to, int x0_matr2_to,  int subarea_len_y, int subarea_len_x);


//
//   ========= matr_cross_add ============
//
//  using procedures
//
//  level 1
//  -------
//  1) bool matr_t1 :: shift_elements_spec_fast_case(int dy, int dx)
//  2) void matr1_area_spec_add_to_matr2(matr_t1 * matr_dest,  int matr_dest_start_pos,          matr_t1 *matr_src,  int matr_src_subarea_y0,  int matr_src_subarea_x0,  int matr_src_subarea_len_y,  int matr_src_subarea_len_x,  int matr_src_forbidsubarea_y0,  int matr_src_forbidsubarea_x0,  int matr_src_forbidsubarea_len_y,  int matr_src_forbidsubarea_len_x)
//  3) matr1_area_spec_mutual_add_swap(matr_t1 * matr_1, int y0_matr1_to, int x0_matr1_to,     matr_t1 * matr_2, int y0_matr2_from, int x0_matr2_from, int y0_matr2_to, int x0_matr2_to,  int subarea_len_y, int subarea_len_x)
//
//  level 2
//  -------
//  1) void two_matrices_mutual_add_with_1_elem(matr_t1 *matr1, int elem1_index_from, int elem1_index_to,    matr_t1 *matr2, int elem2_index_from, int elem2_index_to)
//
//
//  common scheme
//   _____________________________________
//
//  | was 
//  |    matr1                    matr2
//  |    @                            @
//  |    * *                       x  x 
//  |    *    *_shift1          x     x_shift2
//  |    *       *           x        x
//  |    *          *     x           x
//  |    *             @              x
//  |    *          x     *           x
//  |    * _add_ x            * _add_ x 
//  |    *    x                   *   x
//  |    * x                         *x   
//  |    @                            @ 
//  |    matr1                    matr2
//  | became
//  |
//  V t
// 
//    shift up and right only (not left and not down)
//

//    
//   notation system
//
//     AA - matr1 aubarea    BBB - matr2 subarea 
//     AA   (sm1)            BBB   (sm2)
//
//   * A A * *             * * B B B
//   * A A * *             * * B B B
//   * * * * *             * * * * *
//   * * * * *             * * * * *
//   * * * * *             * * * * *
//   matr 1                matr 2
//
//   ism1_inm2:    image of matr 1 subarea in matr 2
//   ism2_inm1:    image of matr 2 subarea in matr 1
//
//   ism2acrsm1_inm1:    image of matr 1 subarea in matr 2
//   ism1acrsm2_inm2:    image of acrossing (subarea ofmatr 1 and subarea ofmatr 2   in matr 1)  in matr 2

//       ________________________
//      V                        |
//   * * * * *             * * * * *
//   * * * * *             * * f f *
//   * a a * *             * t @ f *
//   * a a * *             * t t * *
//   * * * * *             * * * * *
//      |____________________^
//
//
void matr_cross_add_v1(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y,  int shf1_x,   int shf2_y,  int shf2_x);      // ! ! ! ======= very important and difficult proc  ======= ! ! !
void matr_cross_add_v2(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y);
void matr_cross_add_v3(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y);
void matr_cross_add_v4(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y,  int shf2_y);
bool matr_cross_add_v5(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y,  int shf2_y);
bool matr_cross_add_v6(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y,  int shf2_y);
//
//  optimal enumeration one line of bits, returns multiplier
//  ! ! ! == very important proc impacted on optimization== ! ! !
//  if (start_neo_count_and_do_multipl == true) 
//  {
//  multiply having matrix on already storing multiplier (it was returned by this proc on previous teps)
//  we must store-rewrite returned multiplier of mean of proc
//  we must zero matrix  
//  }
// 
//  profit in performance:  size_line:   7: 
//
int optim_enumer_one_line_pbc(unsigned long long int & state, int len, bool & start_neo_count_and_do_multipl);
//  matr_dest and matr_src have different sizes ! !
void matr_symmetr_reflect(matr_t1 * matr_dest,  matr_t1 * matr_src,   bool is_refl_last_col=true);
//  matr have dif sizes
void matr_half_to_whole_fill(matr_t1 * matr_dest,  matr_t1 * matr_src);
void matr_yaxis_invers_symmetr_selfadd(matr_t1 * matr); 




// ----------------------------------------------------
// in out procedures
// ----------------------------------------------------
//
//
//     [ O ]
//       \ \      p
//        \ \  \o/
//         \ \--'---_
//         /\ \   / ~~\_
//   ./---/__|=/_/------//~~~\
//  /___________________/O   O \
//  (===(\_________(===(Oo o o O)        
//   \~~~\____/     \---\Oo__o--
//     ~~~~~~~       ~~~~~~~~~~
//
//
//
//  format of out data
//
//  #  dim
//  2
//  #  lattice type
//  rect
//  #  bound cond
//  fbc
//  #  Ly
//  6
//  #  Lx
//  6
//  #  amount of threads 
//  4
//  #  cur thread number 
//  -1
//  #  last completed row main state
//  1636111
//  #  lattice struct
//  |    |    |    |
//  * -- * -- * -- * -- 
//  |    |    |    |
//  * -- * -- * -- * --
//  |    |    |    |
//  * -- * -- * -- * --
//  |    |    |    |
//  * -- * -- * -- * --
//  #  min energy
//  -125
//  #  J+ percent
//  20
//  #  J- percent
//  80
//  #  J table 
//  #  s1_index_y	s1_index_x	s2_index_y	s2_index_x	bond_mean
//  5		        1		0		1		-1
//  5		        1		5		2		+1
//  5		        2		0		2		+1
//  #  J-matr
//    +   +   +   +   +   +   +   +
//  + * - * - * + * - * - * + * - * +
//    -   +   -   +   -   -   +   -
//  + * - * + * - * + * - * - * - * +
//    -   -   +   +   +   -   +   +
//  + * + * - * - * + * + * + * + * +
//    -   +   -   +   +   +   +   -
//  + * + * + * - * - * - * + * + * +
//    -   +   -   +   +   -   -   +
//  + * + * + * - * - * + * + * + * +
//    -   +   +   -   +   +   +   -
//  + * - * + * + * - * - * - * + * +
//    +   -   +   -   +   -   +   +
//  + * + * - * - * - * - * + * + * +
//    -   +   -   -   -   -   +   +
//  + * - * + * - * - * - * - * + * +
//    +   +   +   +   +   +   +   +
//  #  spin configuration for
//  ?????
//  #  spin configuration:
//  # s_y   s_x     spin
//  0       0       -1
//  0       1       -1
//  1       0       +1
//  1       1       -1
//  #  hist E M g 
//  #  ================================================================
//  -13   -19   971274912840123402394
//  -35   -62   43895723845244
//  -23   -65   19284078234982347823
//  -18   -90   87395782398752982347823
//
//
//
//  #  subhist for boundary state
//  0
//  #  subhist start 
//  #  ================================================================
//  -2   -9    092397856199010
//  -13  -73   3284
//  #  subhist end
//  #  subhist for boundary state
//  1
//  #  subhist start 
//  #  ================================================================
//  -3   -15    52311234
//  -24  -88    1122
//  # subhist end 
//  #  subhist for boundary state
//  2
//  #  subhist start 
//  #  ================================================================
//  -89   -75    121
//  -12   -67    8564543
//  # subhist end 
//
//  writing and extracting data
//  save J matr to file:       izing_EA_2D_rect_fbc_06x06__J                  +postfix    + .txt
//  save hist matr to file:    izing_EA_2D_rect_fbc_06x06___thrd_06_of_16__hist   +postfix    + .txt
//

int  matr_t1_to_chstr_hist(matr_t1 * matr_hist,  char *& chstr_text,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment=true,  int Ly=0,  int Lx=0);
int  matr_t1_to_file(matr_t1 * matr_t1_var,  char * chstr_filename,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment=true,  int Ly=0,  int Lx=0,  bool newfile=true);
void matr_t1_to_screen(matr_t1 * matr_hist,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment=true,  int Ly=0,  int Lx=0);
int  matr_t1_to_short_display(matr_t1 * matr_hist,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment=true,  int Ly=0,  int Lx=0);
bool extract_2D_rect_lattice_hist_from_file(char * filename, int Ly, int Lx, matr_t1 *matr, int E_min, int E_max, int M_min, int M_max, int dE_dcell, int dM_dcell);
bool read_J_bonds_from_file(char * chstr_filename_with_J_bonds,  bool * & J_ver,  bool * & J_hor,  int & Ly,  int & Lx, bool is_1pbc_0fbc=false);







void debug_display_all_matr(matr_t1 * matr_hist, int amount_of_matr,  int int_start_pos)
{

cout<<"===== debug analysis start ====="<<endl; 

  for (int i=0;  i < amount_of_matr;  i++)
  {
  cout<<i<<") state:"<<endl;
  VDisplayFixLenULLIntToBin((unsigned long long int) 3, 3, false, true); cout<<endl;
  VDisplayFixLenULLIntToBin((unsigned long long int) i, 3, false, true); cout<<endl;
  VDisplayFixLenULLIntToBin((unsigned long long int) 0, 3, false, true); cout<<endl;
  cout<<endl; 
  matr_hist[i].display_on_screen_whole_matr();
  }

cout<<"===== debug analysis end   ====="<<endl;

}








struct st_spin_config_int_ar
{

  st_spin_config_int_ar() {bit_row_state=0;  max_amount_of_spin_config=0;  current_cycle_out=0;  wished_cycle_out=0;   E_pospair=0;  current_g=0;}
  void set(unsigned long long int max_amount_of_spin_config_arg) 
  {free(); max_amount_of_spin_config=max_amount_of_spin_config_arg; bit_row_state=new unsigned int[max_amount_of_spin_config*__MACR_LATTICE_LY]; current_cycle_out=0;  wished_cycle_out=0;  E_pospair=0;}
  void set_states_to_zero() {std::memset(bit_row_state, max_amount_of_spin_config*__MACR_LATTICE_LY*sizeof(unsigned int), (char) 0);    E_pospair=0;}
  void free() {  if (bit_row_state) {delete[] bit_row_state;  bit_row_state=0;}  was_overfilled=false;  E=0;    E_pospair=0;  current_g=0;  max_amount_of_spin_config=0;    current_cycle_out=0;  wished_cycle_out=0;}
  ~st_spin_config_int_ar() {free();}

  unsigned int *bit_row_state;
  unsigned long long int max_amount_of_spin_config;

  int E, E_pospair;
  unsigned long long int current_g;
  bool was_overfilled;
  unsigned long long int current_cycle_out;
  unsigned long long int wished_cycle_out;
  
  
  //  #define __MACR_LATTICE_LY 13
  //  #define __MACR_LATTICE_LX 13
 
  void set_default() {was_overfilled=false;  E=0;  current_g=0;    current_cycle_out=0;  wished_cycle_out=0;}
  inline void set_row0_and_row1_state_to_end_for_1st_matr(unsigned short int spin_config_row0_state, unsigned short int spin_config_row1_state, int E_)  
                                                          {bit_row_state[0]=(unsigned int) spin_config_row0_state;  bit_row_state[1]=(unsigned int) spin_config_row1_state;  E=E_;  current_g=1;  was_overfilled=false;}
  inline void set_row0_state_to_end_for_1st_matr(unsigned short int spin_config_row0_state, int E_)  
                                                          {bit_row_state[0]=(unsigned int) spin_config_row0_state;  E=E_;  current_g=1;  was_overfilled=false;}                                                          
  inline void sort_by_M( );  

                                                            
  inline void set_rowi(int row_num,  unsigned int spin_config_row_state)  {  for (unsigned long long int i=0;  i < current_g;  i++) {bit_row_state[i*__MACR_LATTICE_LY+row_num]=spin_config_row_state;}  }
  inline void set_rowi(int row_num,  unsigned int spin_config_row_state, unsigned long long int start_config_index, unsigned long long int amount_adding_configs)  
              {  for (unsigned long long int i=start_config_index;  ((i < (start_config_index+amount_adding_configs)) && (i < current_g));  i++) {bit_row_state[i*__MACR_LATTICE_LY+row_num]=spin_config_row_state;}  }    
  inline void set_rowi(int row_num,  unsigned int spin_config_row_state, unsigned long long int config_num)  
              {  if ((config_num <  current_g) && (row_num < __MACR_LATTICE_LY)) {bit_row_state[config_num*__MACR_LATTICE_LY+row_num]=spin_config_row_state;}  }                  
  inline unsigned long long int replace_row_state_ar_with_nouv(const st_spin_config_int_ar & o_spin_config_int_ar_arg);
  inline unsigned long long int add_row_state_nouv(const st_spin_config_int_ar & o_spin_config_int_ar_arg);
  inline unsigned long long int add_row_state_nouv_rand(const st_spin_config_int_ar & o_spin_config_int_ar_arg);  
  //int E_

  bool convert_bit_row_state_to_bool_ar(unsigned long long int bit_row_state_number, bool * const spin_array, int *Ly_arg=0,  int *Lx_arg=0);
  void write_to_screen(unsigned long long int bit_row_state_number,  char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, bool is_pbc1_fbc0=true);
  int write_to_file(char * chstr_filename, unsigned long long int bit_row_state_number,  bool is_new,  char ch_spin_up,  char ch_spin_dn='d', bool add_info=true, bool is_pbc1_fbc0=true);
  void write_to_screen_all_config(char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, bool is_pbc1_fbc0=true);
  int write_to_file_all_config(char * chstr_filename, bool is_new,  char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, bool is_pbc1_fbc0=true);
  //  here  J_ver and J_hor can be 0
  void write_to_screen_with_J_bonds_one_config(unsigned long long int config_num, bool * J_ver, bool * J_hor, char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, char ch_ver_frame=' ', char ch_hor_frame='=', bool is_pbc1_fbc0=true, bool is_do_frame=true, char bond_bkg_symb=' ', bool is_check_E=true, int *Ly_arg=0,  int *Lx_arg=0);
  void write_to_screen_with_J_bonds_all_config(bool * J_ver, bool * J_hor, char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, char ch_ver_frame=' ', char ch_hor_frame='=', bool is_pbc1_fbc0=true);
  int write_to_file_with_J_bonds_all_config(char * chstr_filename, bool is_new, bool * J_ver, bool * J_hor, char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, char ch_ver_frame=' ', char ch_hor_frame='=', bool is_pbc1_fbc0=true);
  bool write_to_screen_with_J_bonds_all_config_with_E_check(bool * J_ver, bool * J_hor, char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, char ch_ver_frame=' ', char ch_hor_frame='=', bool is_pbc1_fbc0=true);
  bool check_energy(bool * J_ver, bool * J_hor, short int E_must_be, bool is_pbc1_fbc0, unsigned long long int start_state=0, unsigned long long int amount_of_states=(~0));
  bool check_energy_edge_v3(bool *J_ver, bool * J_hor, int  amount_of_null_edge_bonds,  bool * bool_edge_spin_exist_d_u_l_r_ar,  short int E_must_be, unsigned long long int start_state=0, unsigned long long int amount_of_states=(~((unsigned long long int) 0)), int *Ly_arg=0,  int *Lx_arg=0);
  
  //void write_to_screen(bool * J_ver, bool * J_hor, unsigned long long int bit_row_state_number,  char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, bool is_pbc1_fbc0=true);
  
  
  //void write_to_screen_2_state_cmp(unsigned long long int bit_row_state_number,  st_spin_config_int_ar & o_spin_config_int_ar_inst_2,  unsigned long long int bit_row_state_number_inst_2);
  void write_to_screen_2_state_cmp(unsigned long long int bit_row_state_number,  st_spin_config_int_ar & o_spin_config_int_ar_inst_2,  unsigned long long int bit_row_state_number_inst_2);  
  bool write_to_screen_2_state_cmp_2(char *chstr_filename_J_bond,  unsigned long long int bit_row_state_number,  st_spin_config_int_ar & o_spin_config_int_ar_inst_2,  unsigned long long int bit_row_state_number_inst_2,   int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x,  bool is_flip_bond_1ver_0hor, int flip_bond_index_y, int flip_bond_index_x); 


  void write_to_screen_with_J_bonds_config_edge_v3(bool * J_ver, bool * J_hor, int config_num, int amount_of_null_edge_bonds, bool *bool_edge_spin_exist_d_u_l_r_ar=0, bool check_E=false, char ch_spin_up='u', char ch_spin_dn='d', bool add_info=true, char ch_ver_frame=' ', char ch_hor_frame='=',  int *Ly_arg=0,  int *Lx_arg=0);
  
  bool write_to_screen_2_state_edge_cmp_2_v3( st_spin_config_int_ar * o_spin_config_int_ar_inst_to_cmp,  unsigned long long int number_of_gs_1,  unsigned long long int number_of_gs_2_to_cmp,  bool * J_ver,  bool * J_hor,  bool * array_spin_exist_d_u_l_r,   bool is_flip_bond_1ver_0hor, int flip_bond_index_y, int flip_bond_index_x, bool is_draw_domain_wall=true); 
  
};

 


void check_histogr_matr_sum_all_elem(matr_t1 * matr_hist)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

mpz_t sum;
mpz_init2(sum, __MACR_MPZ_BIT_SIZE*2);
mpz_set_ui(sum, (unsigned long long int) 0);  //  cout<<"g: "<<element[i]<<" :g"<<endl;

  for (int y_index=matr_hist->min_not_0_y_index;  y_index <= matr_hist->max_not_0_y_index;  y_index++)
  {
    for (int x_index=matr_hist->min_not_0_x_index;  x_index <= matr_hist->max_not_0_x_index;  x_index++)
    {
      if (mpz_sgn(matr_hist->element[y_index*matr_hist->size_x+x_index]) != 0)
      {mpz_add(sum,  sum,  matr_hist->element[y_index*matr_hist->size_x+x_index]);}
    }
  }

cout<<"cheking sumof hist is: "<<sum<<endl;

mpz_clear(sum);

#endif

} 



void check_pow_of(int mypow)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

mpz_t result;
mpz_init2(result, __MACR_MPZ_BIT_SIZE);
mpz_set_ui(result, (int) 1);  //  cout<<"g: "<<element[i]<<" :g"<<endl;

  if (mypow > 1)
  {
  mpz_t two;
  mpz_init2(two, __MACR_MPZ_BIT_SIZE);
  mpz_set_ui(two, (unsigned int) 2);  //  cout<<"g: "<<element[i]<<" :g"<<endl;
  mpz_set_ui(result, (unsigned int) 2);
  
    for (unsigned int pow_loc=1;  pow_loc < mypow;  pow_loc++)
    {mpz_mul(result,  result,  two);}   

  mpz_clear(two);
  }

cout<<"total number of all states of lattice must be: "<<result<<endl;

mpz_clear(result);

#endif

}




//  02.08.2021
//  frame bond calc proj v03
//
//

// (1)
//
//void shift_J_bonds(bool *J_ver,  bool *J_hor,  int shift_dy,  int shift_dx,  int Ly,  int Lx)

// (2)
//
//  extracting as spins, so crds are pointed as spin indexes
//  example:  we point crd of 4 spins but all 12 bonds are extracted  
//
//      |      |
//    - * -  - * -
//      |      |
//    - * -  - * -
//      |      |

bool extract_J_bond_submatr_from_ver(bool * J_bond_ver_src, int Ly_src, int Lx_src, bool *& J_bond_ver_dest, int spin_start_pos_y_src,  int spin_start_pos_x_src,   int spin_dy,  int spin_dx);

// (3)
//
bool extract_J_bond_submatr_from_hor(bool * J_bond_hor_src, int Ly_src, int Lx_src, bool *& J_bond_hor_dest, int spin_start_pos_y_src,  int spin_start_pos_x_src,   int spin_dy,  int spin_dx);

// (4)
//

//  
//  spin_dy=3   int spin_dx=3
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

bool extract_environment_4(bool *spin_array_src, int Ly, int Lx,  int down_y, int left_x, int dy, int dx,  bool *& bool_edge_spin_exist_d_u_l_r_ar, bool pbc_1_fbc_0); 


// (5)
//
bool rand_generate_edges_spin_val(bool * bool_edge_spin_exist_d_u_l_r_ar, int dy, int dx);


// (6)
//
//
//   example:  dy=4,  dx=3          
//   down: spin_array_src[0]  spin_array_src[1]  spin_array_src[2]          up: spin_array_src[3]  spin_array_src[4]  spin_array_src[5]            left: spin_array_src[6]  spin_array_src[7]            right: spin_array_src[8]  spin_array_src[9]  
//
//
//   -------------------------
//   | *         *         * |
//   -------------------------
//   -----               -----
//   |   |               |   |
//   | * |               | * |
//   |   |               |   |
//   |   |               |   |
//   | * |               | * |
//   |   |               |   |
//   -----               ----- 
//   -------------------------
//   | *         *         * |
//   -------------------------

bool extract_frame_spin_value(bool *spin_array_src, int Ly, int Lx,  bool *& spin_frame_value_4, int down_y, int left_x, int dy, int dx, bool pbc_1_fbc_0);

bool paste_spin_sublat_into_spin_lat(bool * spin_array_dest,  int Ly_dest,  int Lx_dest,  bool * spin_array_src, int dest_down_y, int dest_left_x, int dy, int dx, bool pbc_1_fbc_0);

bool calc_amount_of_null_edge_bonds_edge_v3(int Ly,  int Lx,  bool * bool_edge_spin_exist_d_u_l_r_ar);   

bool print_struct_with_J_edges_v3(bool *J_ver, bool *J_hor, int Ly,  int Lx,  bool * bool_edge_spin_exist_d_u_l_r_ar);

int match_test_for_0_radius_struct(int N);      //  !  !  !

bool print_bool_edge_spin_exist_d_u_l_r_ar_for_debug(int Ly,  int Lx,  bool * bool_edge_spin_exist_d_u_l_r_ar);      //  !  !  !






unsigned long long int ullint_mylog_base_2(unsigned long long int num) {  if (num < 2) {return 0;}   unsigned long long int answer=1;  for (unsigned long long int i=0;  i < 1000; i++) {answer*=2;  if (num <= answer) {return (i+1);}  }  return 0;}


bool calc_min_E_partially_into_mem_chunk__for_given_first_down_row_and_col(short int *min_E_array,  int num_thrd,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  unsigned long long int amount_of_states_per_thrd,  int Ly,  int Lx,  bool * J_ver,  bool * J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, int amount_of_null_edge_bonds=0);

bool calc_min_E_partially_into_mem_chunk__for_given_row_and_col(short int *min_E_array,  unsigned long long int amount_of_states_per_thrd_div_2,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  int amount_of_steps_into_thread, int row,  int col,  int Ly,  int Lx,  bool *& J_ver,  bool *& J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, short int *result_in_thrd_min_E=0,  unsigned long long int *up_row_gs_state_in_thrd=0, int amount_of_null_edge_bonds=0);

bool swap_mem(short int *min_E_array,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  unsigned long long int mem_per_thrd_in_unit_sint,  int thrd_num_1, int thrd_num_2);


bool swap_mem_paral(short int *min_E_array,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  unsigned long long int mem_per_thrd_in_unit_sint,  int thrd_num_1, int thrd_num_2,  int row,  int col,  int thread_id_code_shift=0);























//  given_amount_of_gs_match_states_to_write:  amount of gs we want with matching frame,   if (wi  will get less) {we change given_amount_of_gs_match_states_to_write  on number of found gs configs}
//
//  4 refs only:    1) int &up_choiced_state,  2) int &up_choiced_state,  3) unsigned long long int & amount_of_gs_frame_matching=0,   4) unsigned long long int & amount_of_gs_match_states_to_write=1


//  numeration notation: 
//  r0col3 r0col2 r0col1 r0col0
//  r1col3 r1col2 r1col1 r1col0  
//  r2col3 r2col2 r2col1 r2col0
//
//  bit3   bit2   bit1   bit0
//
//  
//  row_state_range_out=3  row_state_range_in_div_2=4
//
//  *  *  *  * | *  *  *  *      *  *  *  * | *  *  *  *      *  *  *  * | *  *  *  *              :  * means one state of row
//  |__|__|__|___|  |  |  |
//     |__|__|______|  |  |
//        |__|_________|  |
//           |____________|
//          

//


bool izing_2D_exact_gbc_EAnd_v55_minE_sr_parallel_export_variant(int Ly,  int Lx,  bool *&J_ver,  bool *&J_hor,  bool * & bool_edge_spin_exist_d_u_l_r_ar, int posit_bond_perc_arg,  int negat_bond_perc_arg, bool is_to_do_record_in_file=true,  bool is_print_comments=true, bool do_enum_check=true, short int * result_min_E_arg=0,  bool *bool_one_up_row_gs_state=0,  bool *check_is_ok_arg=0,  bool print_progress_arg=true)
{ 

//  if (bool_edge_spin_exist_d_u_l_r_ar == 0) {return false;}   
//  if ((J_ver == 0) || (J_hor == 0)) {return false;} 


bool is_write_J_config_to_file=false;   
int thrd_amount=mpisize; 
int  result_min_E=((int) (~((unsigned int) 0)));
short int  result_min_E_pospair=256*128-1;
unsigned long long int ullint_one_up_row_gs_state=0;   

unsigned long long int state_size_in_array=sizeof(unsigned long long int)/sizeof(result_min_E_pospair);
    
   
   

  if ((mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT) && (Ly > 1) && (Lx > 1))
  {            
    if (is_print_comments == true) 
    {
    cout<<endl;
    cout<<endl<<"===============  ---------------------------    izing_2D_exact_gbc_EAnd_v54_minE_sr_mult_thrd_parallel   is launched    ------------------------   ============"<<endl;
    cout<<"===============  --------------------------          ================================================================   ------------------------   ============"<<endl<<endl;  
    }
  }





// ====================================================================================================  // ====================  MPI
// ============================  J-matr transmit from 0 to other threads ==============================  // ====================  MPI
//
//

//  transfer:  Ly
//  transfer:  Lx
//  transfer:  posit_bond_perc_arg
//  transfer:  negat_bond_perc_arg
//  transfer:  bool_edge_spin_exist_d_u_l_r_ar
//  transfer:  J_ver
//  transfer:  J_hor


  {  
  //  transfer:  Ly
  //  transfer:  Lx
  //  transfer:  posit_bond_perc_arg
  //  transfer:  negat_bond_perc_arg

    if (mpisize > 1)
    {
    char *chstr_array=0;
    int chstr_size=sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int);
    chstr_array=new char[chstr_size+1];    for (int i=0;  i < (chstr_size);  i++) {chstr_array[i]=0;}    chstr_array[0]='\0';
    
      if (mpirank > 0) 
      {
      MPI_Recv(chstr_array, chstr_size, MPI_CHAR, 0, 3, MPI_COMM_WORLD, &Status);      //  !  !  !
      std::memcpy(& Ly,                  & chstr_array[0            ], sizeof(int));  std::memcpy(& Lx,                  & chstr_array[1*sizeof(int)], sizeof(int));  
      std::memcpy(& posit_bond_perc_arg, & chstr_array[2*sizeof(int)], sizeof(int));  std::memcpy(& negat_bond_perc_arg, & chstr_array[3*sizeof(int)], sizeof(int)); 
      }      //  ! ! !
      else
      {
      std::memcpy(& chstr_array[0            ], & Ly,                  sizeof(int));  std::memcpy(& chstr_array[1*sizeof(int)], & Lx,                  sizeof(int));  
      std::memcpy(& chstr_array[2*sizeof(int)], & posit_bond_perc_arg, sizeof(int));  std::memcpy(& chstr_array[3*sizeof(int)], & negat_bond_perc_arg, sizeof(int)); 
      chstr_array[4*sizeof(int)]='\0';
      
        for (int i=1; i < mpisize; i++)      //  thrd 0  send array to every other thrd
        {MPI_Send(chstr_array,  chstr_size, MPI_CHAR, i, 3, MPI_COMM_WORLD);}      //  !  !  !
    
      }  

      if (chstr_array) {delete[] chstr_array;  chstr_array=0;}
    }
    
  }


  if (Lx < 2) {return false;}
  if (Ly < 2) {return false;}    


  {  
  //  transfer:  bool_edge_spin_exist_d_u_l_r_ar

    if (mpisize > 1)
    {
    char *chstr_array=0;
    int chstr_size=4*Ly+4*Lx;
    chstr_array=new char[chstr_size+1];    for (int i=0;  i < (chstr_size+1);  i++) {chstr_array[i]='0';}    chstr_array[0]='\0';
    

      if (mpirank > 0) 
      {
      MPI_Recv(chstr_array, chstr_size, MPI_CHAR, 0, 4, MPI_COMM_WORLD, &Status);      //  !  !  !
      bool_edge_spin_exist_d_u_l_r_ar=new bool[4*Ly+4*Lx];  for (int i=0;  i < 4*Ly+4*Lx;  i++)  {bool_edge_spin_exist_d_u_l_r_ar[i]=false;}
      
        for (int i=0;       i < chstr_size;  i++)  {  if (chstr_array[i]          == '1') {bool_edge_spin_exist_d_u_l_r_ar[i]=true;} else {bool_edge_spin_exist_d_u_l_r_ar[i]=false;}  } 
      }      //  ! ! !
      else
      {
        for (int i=0;       i < chstr_size;  i++)  {  if (bool_edge_spin_exist_d_u_l_r_ar[i] == true) {chstr_array[i]='1';         } else {chstr_array[i]='0';         }  }
      chstr_array[chstr_size]='\0';

        for (int i=1; i < mpisize; i++)      //  thrd 0  send array to every other thrd
        {MPI_Send(chstr_array,  chstr_size, MPI_CHAR, i, 4, MPI_COMM_WORLD);}      //  !  !  !
    
      }  

      if (chstr_array) {delete[] chstr_array;  chstr_array=0;}
    }
    
  }
  

  {  
  //  transfer:  J_ver
  //  transfer:  J_hor

    if (mpisize > 1)
    {
    char *chstr_array=0;
    int chstr_size=Ly*Lx+Lx + Ly*Lx+Ly;  chstr_array=new char[chstr_size+1];  chstr_array[0]='\0';

      if (mpirank > 0) 
      {
      MPI_Recv(chstr_array, chstr_size, MPI_CHAR, 0, 5, MPI_COMM_WORLD, &Status);      //  !  !  !
            
      J_ver=new bool[Lx*Ly+Lx];    for (int i=0;  i < Ly*Lx+Lx;  i++)  {J_ver[i]=false;}
      J_hor=new bool[Lx*Ly+Ly];    for (int i=0;  i < Ly*Lx+Ly;  i++)  {J_hor[i]=false;}   

  
        for (int i=0;       i < Ly*Lx+Lx;  i++)  {  if (chstr_array[i]          == '1') {J_ver[i]=true;} else {J_ver[i]=false;}  }
        for (int i=0;       i < Ly*Lx+Ly;  i++)  {  if (chstr_array[Ly*Lx+Lx+i] == '1') {J_hor[i]=true;} else {J_hor[i]=false;}  }
      }      //  ! ! !
      else
      {
        for (int i=0;       i < Ly*Lx+Lx;  i++)  {  if (J_ver[i] == true) {chstr_array[i]='1';         } else {chstr_array[i]='0';         }  }
        for (int i=0;       i < Ly*Lx+Ly;  i++)  {  if (J_hor[i] == true) {chstr_array[Ly*Lx+Lx+i]='1';} else {chstr_array[Ly*Lx+Lx+i]='0';}  }
      chstr_array[chstr_size]='\0';

        for (int i=1; i < mpisize; i++)      //  thrd 0  send array to every other thrd
        {MPI_Send(chstr_array,  chstr_size, MPI_CHAR, i, 5, MPI_COMM_WORLD);}      //  !  !  !
    
      }  

      if (chstr_array) {delete[] chstr_array;  chstr_array=0;}
    }
    
  }
//  
// ====================================================================================================  // ====================  MPI
// ====================================================================================================  // ====================  MPI
  







//  ========  *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *
//  prepare
//  ========
//  vars and arrays prepare 
//
ofstream out_progress("progres2D.txt"); 
double mustachieve_progress_perc=0;
unsigned long long int max_row_state=ullint_2_pow_n(Lx), max_row_state_div_2=max_row_state/2;







//  ========  *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *     *
//  prepare
//  ========
//
//  ---------------------------------- v51 additive ------------- V
//
  if (bool_edge_spin_exist_d_u_l_r_ar == 0)
  {
  bool_edge_spin_exist_d_u_l_r_ar=new bool[Ly+Lx+Ly+Lx];
  
    for (int y=0;  y < Ly;  y++)
    {
    bool_edge_spin_exist_d_u_l_r_ar[Lx*4+2*y]=((bool) (rand() % 2));  bool_edge_spin_exist_d_u_l_r_ar[Lx*4+2*y+1]=true;    bool_edge_spin_exist_d_u_l_r_ar[Lx*4+Ly*2+2*y]=((bool) (rand() % 2));   bool_edge_spin_exist_d_u_l_r_ar[Lx*4+Ly*2+2*y+1]=true;
    
    }

    for (int x=0;  x < Lx;  x++)
    {
    bool_edge_spin_exist_d_u_l_r_ar[2*x]=((bool) (rand() % 2));  bool_edge_spin_exist_d_u_l_r_ar[2*x+1]=true;    bool_edge_spin_exist_d_u_l_r_ar[Lx*2+2*x]=true;   bool_edge_spin_exist_d_u_l_r_ar[Lx*2+2*x+1]=((bool) (rand() % 2));
    }          
  }


int amount_of_null_edge_bonds=0;      //  calc here

  {
    for (int y=0;  y < Ly;  y++)
    {
      if (bool_edge_spin_exist_d_u_l_r_ar[Lx*4+2*y+1] == false)  {amount_of_null_edge_bonds++;}    if (bool_edge_spin_exist_d_u_l_r_ar[Lx*4+Ly*2+2*y+1] == false)  {amount_of_null_edge_bonds++;}  
    }

    for (int x=0;  x < Lx;  x++)
    {
      if (bool_edge_spin_exist_d_u_l_r_ar[2*x+1] == false)  {amount_of_null_edge_bonds++;}         if (bool_edge_spin_exist_d_u_l_r_ar[Lx*2+2*x+1] == false)  {amount_of_null_edge_bonds++;}  
    }          

  }
//
//  ---------------------------------- v51 additive ------------- ^





//==DEBUG==//  cout<<"mpirank="<<mpirank<<" amount_of_null_edge_bonds="<<amount_of_null_edge_bonds<<endl;  //==DEBUG==//









//  prepare
//  ========
//  vars and arrays for stage 1 prepare 
//
unsigned long long int amount_of_states_per_thrd      =max_row_state/thrd_amount;
unsigned long long int amount_of_states_per_thrd_div_2=amount_of_states_per_thrd/2;

unsigned long long int mem_per_thrd_in_unit_sint=3*(state_size_in_array+amount_of_states_per_thrd/2);                    //  the first 4 short int (8 byte) stores unsigned long long int state_code for every third
unsigned long long int mem_part_third_in_thrd_in_unit_sint=(mem_per_thrd_in_unit_sint)/3;
unsigned long long int memory_for_all_thrd_emul_in_unit_sint=mem_per_thrd_in_unit_sint*thrd_amount;

int amount_of_steps_into_thread=ullint_mylog_base_2(amount_of_states_per_thrd);
int amount_of_steps_with_mem_thread_swap=Lx-amount_of_steps_into_thread;







  if (mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT)
  {
    //if (is_print_comments == true)
    {
    cout<<"boundary conditions: "<<"given edges (d, u, l, r)"<<endl;
    cout<<"Ly="<<Ly<<", Lx="<<Lx<<endl;
    cout<<"posit_bond_perc="<<posit_bond_perc_arg<<", negat_bond_perc="<<negat_bond_perc_arg<<endl;
    cout<<"amount of threads = "<<thrd_amount<<endl;
    cout<<"amount_of_states_per_thrd = "<<amount_of_states_per_thrd<<endl;
    
    cout<<"using memory per one thread for storing energy array: ";

    long long int mem_size_bit=mem_per_thrd_in_unit_sint*sizeof(result_min_E_pospair);
    long long int mem_size_Gb=mem_size_bit/8/1024/1024/1024; 
    mem_size_bit-=mem_size_Gb*1024*1024*1024*8;
    long long int mem_size_Mb=mem_size_bit/8/1024/1024; 
    mem_size_bit-=mem_size_Mb*1024*1024*8;
    long long int mem_size_kb=mem_size_bit/8/1024;
    mem_size_bit-=mem_size_kb*1024*8;
    long long int mem_size_b =mem_size_bit/8;
    cout<<mem_size_Gb<<" Gb + "<<mem_size_Mb<<" Mb + "<<mem_size_kb<<" kb + "<<mem_size_b<<" b"<<endl;

    cout<<"using memory for all for storing of matrices is about: ";
    mem_size_bit=mem_per_thrd_in_unit_sint*sizeof(result_min_E_pospair)*((long long int) mpisize);
    mem_size_Gb=mem_size_bit/8/1024/1024/1024; 
    mem_size_bit-=mem_size_Gb*1024*1024*1024*8;
    mem_size_Mb=mem_size_bit/8/1024/1024; 
    mem_size_bit-=mem_size_Mb*1024*1024*8;
    mem_size_kb=mem_size_bit/8/1024;
    mem_size_bit-=mem_size_kb*1024*8;
    mem_size_b =mem_size_bit/8;
    cout<<mem_size_Gb<<" Gb + "<<mem_size_Mb<<" Mb + "<<mem_size_kb<<" kb + "<<mem_size_b<<" b"<<endl;
    
    cout<<endl<<"matrix struct:"<<endl;
    print_struct_with_J_edges_v3(J_ver, J_hor,  Ly,  Lx,  bool_edge_spin_exist_d_u_l_r_ar);
    cout<<endl;   
    }
  } 









short int *min_E_array=0;
min_E_array=new short int[mem_per_thrd_in_unit_sint]; 

  for (unsigned long long int row_state_in_row=0;  row_state_in_row < mem_per_thrd_in_unit_sint;  row_state_in_row++)
  {min_E_array[row_state_in_row]=0;}


//  filling with state code in begin of the third part
//
  //<>//for (int thrd_num=0;  thrd_num < thrd_amount;  thrd_num++)
  {
  unsigned long long int start_state=((unsigned long long int) mpirank)*amount_of_states_per_thrd+amount_of_states_per_thrd_div_2*0;
  std::memcpy(& min_E_array[0*mem_per_thrd_in_unit_sint+mem_part_third_in_thrd_in_unit_sint*0], & start_state, sizeof(unsigned long long int));

  start_state=((unsigned long long int) mpirank)*amount_of_states_per_thrd+amount_of_states_per_thrd_div_2*1; 
  std::memcpy(& min_E_array[0*mem_per_thrd_in_unit_sint+mem_part_third_in_thrd_in_unit_sint*1], & start_state, sizeof(unsigned long long int)); 
  }



//==DEBUG==//  cout<<"DEBUG:  ";       //==DEBUG==//  
//==DEBUG==//    for (int k=0;  k < memory_for_all_thrd_emul_in_unit_sint;  k++) {cout<<"min_E_array["<<k<<"]="<<min_E_array[k]<<"  ";}       //==DEBUG==//
//==DEBUG==//  cout<<endl;       //==DEBUG==//
  

//  bool calc_min_E_partially_into_mem_chunk__for_given_first_down_row_and_col(short int *min_E_array,  int num_thrd,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  unsigned long long int amount_of_states_per_thrd,  int Ly,  int Lx,  bool * J_ver,  bool * J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, int amount_of_null_edge_bonds=0);

//  bool calc_min_E_partially_into_mem_chunk__for_given_row_and_col(short int *min_E_array,  unsigned long long int amount_of_states_per_thrd_div_2,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  int amount_of_steps_into_thread, int row,  int col,  int Ly,  int Lx,  bool *& J_ver,  bool *& J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, short int *result_in_thrd_min_E, int amount_of_null_edge_bonds=0);

//  bool swap_mem(short int *min_E_array,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  unsigned long long int mem_per_thrd_in_unit_sint,  int thrd_num_1, int thrd_num_2);








//  a  a  a  a  a    a  a  a  a  a =========================================  V  V
//  ------------------------------------------------------
//  a  a  a  a  a    a  a  a  a  a =========================================  V  V
//  ------------------------------------------------------
//  a  a  a  a  a    a  a  a  a  a =========================================  V  V
//  ------------------------------------------------------
//  a  a  a  a  a    a  a  a  a  a =========================================  V  V
//
//

short int result_min_E_pospair_per_thrd=result_min_E_pospair;
unsigned long long int ullint_one_up_row_gs_state_per_thread=0;  
int interval_amount_of_thrds=0, interval_amount_of_thrds_div_2=0,  amount_of_bunchs_of_connected_thrds=0;


//  calc row 0 (first down row):
//
//==DEBUG==//  cout<<"DEBUG:  "<<"0-row"<<endl;       //==DEBUG==//

  //<>//for (int thrd_num=0; thrd_num < thrd_amount;  thrd_num++)
  {
  int thrd_num=0;
  calc_min_E_partially_into_mem_chunk__for_given_first_down_row_and_col(& min_E_array[mem_per_thrd_in_unit_sint*thrd_num],  thrd_num, mem_part_third_in_thrd_in_unit_sint,  amount_of_states_per_thrd,  Ly,  Lx,  J_ver,  J_hor,  bool_edge_spin_exist_d_u_l_r_ar, amount_of_null_edge_bonds);
  }      //  !  !  !        !  !  !   

  for (int row=1;  row < Ly;  row++) 
  {
  //==DEBUG==// cout<<endl<<"DEBUG:  "<<"row="<<row<<endl;       //==DEBUG==//  
  
  //  !  !  !  into thread:  V
  //
    for (int col=0;  col < amount_of_steps_into_thread;  col++)      //  ==========================================  1) A A  A A  A A  A A  A A  A A  A A  A A  A A  A A  A A  A A
    {
    //==DEBUG==//  cout<<"DEBUG:  "<<"col="<<col<<endl;       //==DEBUG==//
           
      //<>//for (int thrd_num=0; thrd_num < thrd_amount;  thrd_num++)
      {
      int thrd_num=0;
      calc_min_E_partially_into_mem_chunk__for_given_row_and_col(& min_E_array[thrd_num*mem_per_thrd_in_unit_sint],  amount_of_states_per_thrd_div_2,  mem_part_third_in_thrd_in_unit_sint,  amount_of_steps_into_thread,  row,  col,  Ly,  Lx,  J_ver,  J_hor,  bool_edge_spin_exist_d_u_l_r_ar,  & result_min_E_pospair_per_thrd, & ullint_one_up_row_gs_state_per_thread, amount_of_null_edge_bonds);      //  !  !  !        !  !  !    
            
      //<>//  if ( ((col+1 == Lx) && (row+1 == Ly)) == true) {  if (result_min_E_pospair_per_thrd <= result_min_E_pospair) {result_min_E_pospair=result_min_E_pospair_per_thrd;  ullint_one_up_row_gs_state=ullint_one_up_row_gs_state_per_thread;}  }
      }    
            
      if ((mpirank == 0) && (print_progress_arg == true)) {print_progress(row*Lx+col, mustachieve_progress_perc, Ly*Lx-1, out_progress);}      
    }
   
   

     


  //  !  !  !  out thread:  V
  //
  interval_amount_of_thrds=2; interval_amount_of_thrds_div_2=1;  amount_of_bunchs_of_connected_thrds=thrd_amount/2;      if (amount_of_bunchs_of_connected_thrds < 1) {amount_of_bunchs_of_connected_thrds=1;}


    for (int col=amount_of_steps_into_thread;  col < Lx;  col++)      //  ==========================================  2) A A  A A  A A  A A  A A  A A  A A  A A  A A  A A  A A  A A
    {
    //==DEBUG==//  cout<<"DEBUG:  "<<"col="<<col<<endl;       //==DEBUG==//
    
    
    // mem swap     ======  v  v  v  !  !  !
    //
      for (int thrd_chunk_num=0;  thrd_chunk_num < amount_of_bunchs_of_connected_thrds;  thrd_chunk_num++)
      {
        for (int thrd_local_num=0;  thrd_local_num < interval_amount_of_thrds_div_2;  thrd_local_num++)
        {
        int thrd_num_1=thrd_chunk_num*interval_amount_of_thrds+thrd_local_num,  thrd_num_2=thrd_num_1+interval_amount_of_thrds_div_2;
          if ((mpirank == thrd_num_1) || (mpirank == thrd_num_2))
          {swap_mem_paral(min_E_array,  mem_part_third_in_thrd_in_unit_sint,  mem_per_thrd_in_unit_sint,  thrd_num_1, thrd_num_2, row, col, 0);}      //  !  !  !        !  !  !        
        }      
      }   //--//cout<<"debug22: "<<"  "<<"mpirank="<<mpirank<<"  "<<"min_E_array[0]="<<min_E_array[0]<<"  "<<"min_E_array[5]="<<min_E_array[5]<<"  "<<"min_E_array[4]="<<min_E_array[4]<<"  "<<"min_E_array[9]="<<min_E_array[9]<<endl;
    // mem swap     ======  ^  ^  ^  !  !  !  

    
      //<>//for (int thrd_num=0; thrd_num < thrd_amount;  thrd_num++)
      {
      int thrd_num=0;
      calc_min_E_partially_into_mem_chunk__for_given_row_and_col(& min_E_array[thrd_num*mem_per_thrd_in_unit_sint],  amount_of_states_per_thrd_div_2,  mem_part_third_in_thrd_in_unit_sint,  amount_of_steps_into_thread,  row,  col,  Ly,  Lx,  J_ver,  J_hor,  bool_edge_spin_exist_d_u_l_r_ar,  & result_min_E_pospair_per_thrd, & ullint_one_up_row_gs_state_per_thread, amount_of_null_edge_bonds);      //  !  !  !        !  !  !    
            
        //<>//if ( ((col+1 == Lx) && (row+1 == Ly)) == true) 
        //<>//{  if (result_min_E_pospair_per_thrd <= result_min_E_pospair) {result_min_E_pospair=result_min_E_pospair_per_thrd;    ullint_one_up_row_gs_state=ullint_one_up_row_gs_state_per_thread;}  }    

      }


      if (col < Lx)  {interval_amount_of_thrds*=2;  interval_amount_of_thrds_div_2=interval_amount_of_thrds/2;  amount_of_bunchs_of_connected_thrds/=2;}
      if ((mpirank == 0) && (print_progress_arg == true)) {print_progress(row*Lx+col, mustachieve_progress_perc, Ly*Lx-1, out_progress);}  
      
    }      //  for (int col=amount_of_steps_into_thread;  col < Lx;  col++)
         
     
  // mem zuruck (restore order of bytes)     ======  v  v  v  !  !  !
  //
    {
    interval_amount_of_thrds=thrd_amount;  interval_amount_of_thrds_div_2=interval_amount_of_thrds/2;  amount_of_bunchs_of_connected_thrds=1;      if (amount_of_bunchs_of_connected_thrds < 1) {amount_of_bunchs_of_connected_thrds=1;}
    
      for (int col=Lx-1;  col >= amount_of_steps_into_thread;  col--)
      {       
        for (int thrd_chunk_num=0;  thrd_chunk_num < amount_of_bunchs_of_connected_thrds;  thrd_chunk_num++)
        {
          for (int thrd_local_num=0;  thrd_local_num < interval_amount_of_thrds_div_2;  thrd_local_num++)
          {
          int thrd_num_1=thrd_chunk_num*interval_amount_of_thrds+thrd_local_num,  thrd_num_2=thrd_num_1+interval_amount_of_thrds_div_2;
          //--//cout<<"debug2: "<<"  "<<"mpirank="<<mpirank<<"  "<<"thrd_num_1="<<thrd_num_1<<"  "<<"thrd_num_2="<<thrd_num_2<<endl;

            if ((mpirank == thrd_num_1) || (mpirank == thrd_num_2))
            {swap_mem_paral(min_E_array,  mem_part_third_in_thrd_in_unit_sint,  mem_per_thrd_in_unit_sint,  thrd_num_1, thrd_num_2, row, col, 1);}      //  !  !  !        !  !  !
          }      
        }
        
        if (col < Lx)  {interval_amount_of_thrds/=2;  interval_amount_of_thrds_div_2/=2;  amount_of_bunchs_of_connected_thrds*=2;}
        
      }
      
    }
  // mem zuruck (restore order of bytes)     ======  v  v  v  !  !  !
  //     
    
  }      //  for (int row=1;  row < Ly;  row++)






// ====================================================================================================  // ====================  MPI
// ============================  J-matr transmit from 0 to other threads ==============================  // ====================  MPI
//
//
  {  
  char *chstr_array=0;
  int chstr_size=sizeof(result_min_E_pospair_per_thrd)+sizeof(ullint_one_up_row_gs_state_per_thread);
  chstr_array=new char[chstr_size+1];    for (int i=0;  i < (chstr_size);  i++) {chstr_array[i]=0;}    chstr_array[0]='\0';  
  
    if (mpisize > 1)
    {
      if (mpirank > 0) 
      {
      std::memcpy(& chstr_array[0                                       ], & result_min_E_pospair_per_thrd,          sizeof(result_min_E_pospair_per_thrd        ) );
      std::memcpy(& chstr_array[sizeof(result_min_E_pospair_per_thrd)   ], & ullint_one_up_row_gs_state_per_thread,  sizeof(ullint_one_up_row_gs_state_per_thread) );
      
      MPI_Send(chstr_array,  chstr_size, MPI_CHAR, 0, 1, MPI_COMM_WORLD);       //  ! ! !
      }      //  ! ! !
      else
      {
      result_min_E_pospair=result_min_E_pospair_per_thrd;
      ullint_one_up_row_gs_state=ullint_one_up_row_gs_state_per_thread;
     
        for (int i=1; i < mpisize; i++)
        {
        chstr_array[0]='\0';
        MPI_Recv(chstr_array, chstr_size, MPI_CHAR, i, 1, MPI_COMM_WORLD, &Status);
        
        std::memcpy(& result_min_E_pospair_per_thrd,           & chstr_array[0                                       ], sizeof(result_min_E_pospair_per_thrd        ) );
        std::memcpy(& ullint_one_up_row_gs_state_per_thread,   & chstr_array[sizeof(result_min_E_pospair_per_thrd)   ], sizeof(ullint_one_up_row_gs_state_per_thread) );
        
          if (result_min_E_pospair_per_thrd  <=  result_min_E_pospair ) {result_min_E_pospair=result_min_E_pospair_per_thrd;    ullint_one_up_row_gs_state=ullint_one_up_row_gs_state_per_thread;}
        
        }     
            
      }  
    }
    else
    {
    result_min_E_pospair=result_min_E_pospair_per_thrd;
    ullint_one_up_row_gs_state=ullint_one_up_row_gs_state_per_thread;
    }

    if (chstr_array) {delete[] chstr_array;  chstr_array=0;}

  }
//  
// ====================================================================================================  // ====================  MPI
// ====================================================================================================  // ====================  MPI








  if (bool_one_up_row_gs_state)   
  {
  unsigned long long int shift_one=1;
    for (int i=0;  i < Lx;  i++)
    {  if ((ullint_one_up_row_gs_state & (1 << i)) != 0) {bool_one_up_row_gs_state[i]=true;} else {bool_one_up_row_gs_state[i]=false;}  }  
  }
  
  
result_min_E=(-2*Ly*Lx+Ly+Lx-2*Ly-2*Lx+amount_of_null_edge_bonds)+2*result_min_E_pospair;      //  ---------------------------------- v53 modif ------------- <
  
  if (result_min_E_arg) {*result_min_E_arg=result_min_E;}

//
//  a  a  a  a  a    a  a  a  a  a =========================================  ^  ^
//  ------------------------------------------------------
//  a  a  a  a  a    a  a  a  a  a =========================================  ^  ^
//  ------------------------------------------------------
//  a  a  a  a  a    a  a  a  a  a =========================================  ^  ^
//  ------------------------------------------------------
//  a  a  a  a  a    a  a  a  a  a =========================================  ^  ^











//  comment
//  ========
//  out results of stage 1 
//
  if (mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT)
  {    
    if (is_print_comments == true)
    { 
    cout<<endl<<"program is completed   for lattice "<<Ly<<"x"<<Lx<<endl;
    cout<<"result_min_E="<<result_min_E<<"   result_min_E_pospair="<<result_min_E_pospair<<"  amount_of_null_edge_bonds="<<amount_of_null_edge_bonds<<"      "<<endl;
      
  
      if (bool_one_up_row_gs_state)
      {
      cout<<"  up_row_gs_state:"<<endl;  
    
        for (int i=0;  i < Lx;  i++)
        {  if (bool_one_up_row_gs_state[i] == true) {cout<<"+1 ";} else {cout<<"-1 ";}  }
      }  
    cout<<endl<<endl;
    }
  }    





//  recording to file info and result
//  ========
//
//
  if (mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT)
  {

    if (is_to_do_record_in_file == true)
    {
    //char *chstr_filename_with_J_bonds_do_record=0;
    bool is_write_lat_struct=true;
    bool pbc1fbc0=false;
    bool print_J_list=true;
    bool *one_up_row_gs_state_ptr=0;
  
      if (bool_one_up_row_gs_state) {one_up_row_gs_state_ptr=bool_one_up_row_gs_state;}  //if (bool_one_up_row_gs_state) {cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<endl;}
  
    //-- prototype --//
    //  bool print_to_file_2D_rect_gbc_lattice_info_and_ext_edge_and_up_gs_row_v54(int Ly, int Lx, bool * J_ver,  bool * J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, int *min_energ, bool *bool_up_row_gs_ar,  char * dest_filename_arg, int J_plus_perc, int J_min_perc,  bool is_write_lat_struct,  bool pbc1fbc0,  bool print_J_list)  
   
    print_to_file_2D_rect_gbc_lattice_info_and_ext_edge_and_up_gs_row_v54(Ly, Lx, J_ver,  J_hor,  bool_edge_spin_exist_d_u_l_r_ar, & result_min_E,  bool_one_up_row_gs_state,  0, posit_bond_perc_arg,  negat_bond_perc_arg,  is_write_lat_struct,  pbc1fbc0,  print_J_list, true);
    one_up_row_gs_state_ptr=0;
    }

  }














//  checking
//  ========
//  to display some results for checking and self-control
//
//   8 7 6
//   5 4 3
//   2 1 0
//
  if (mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT)
  {
  
    if ((do_enum_check == true) && (Ly*Lx <= 30) && (is_print_comments == true) && (1 == 1)) 
    {
    cout<<endl<<endl<<"========== exact enumertaion total check =========="<<endl;
    unsigned long long int lack_spins_shift_one_in_lattice=1;
    unsigned long long int lack_spins_in_lattice_mask=0;  

      //for (int i=0;  i < Lx*Ly;  i++)  
      for (int i=0;  i < Ly;  i++)
      {  
      lack_spins_shift_one_in_lattice=(((unsigned long long int) 1) << ((i+1)*Lx-1));
        for (int j=0;  j < Lx;  j++)  
        {  
        //  if (vacan_array[i*Lx+j] == true) {lack_spins_in_lattice_mask+=lack_spins_shift_one_in_lattice;}
        //==DEBUG==//   cout<<"here yes "<<i<<" "<<j<<" "<<lack_spins_shift_one_in_lattice<<endl;} else {cout<<"here no "<<i<<" "<<j<<endl;}   
        lack_spins_shift_one_in_lattice=(lack_spins_shift_one_in_lattice >> 1); 
        }
      }  

    matr_t1 matr_check;
    matr_check.init(2*Ly*Lx-Ly-Lx+2*Ly+2*Lx+1, Ly*Lx+1, 0);  
    //matr_check.init(2*Ly*Lx+1, Ly*Lx+1, 0); 
    //  ---------------------------------- v51 additive ------------- V
    //
    bool *out_spin_array=0;    out_spin_array=new bool[Ly*2+Lx*2];
    bool *out_J_array   =0;    out_J_array   =new bool[Ly*2+Lx*2];  
      for (int i=0;  i < Lx;  i++) {out_J_array[i]=J_ver[i];                out_J_array[i+Lx]     =J_ver[Ly*Lx+i];}    
      for (int i=0;  i < Ly;  i++) {out_J_array[Lx*2+i]=J_hor[i*(Lx+1)];    out_J_array[Lx*2+Ly+i]=J_hor[i*(Lx+1)+Lx];}
    //     
    //  ---------------------------------- v51 additive ------------- ^

  
    //==DEBUG==//cout<<"debug: matr_check.size_y="<<matr_check.size_y<<" matr_check.size_x="<<matr_check.size_x<<" lack_spins_in_lattice_mask="<<lack_spins_in_lattice_mask<<endl;  //==DEBUG==//
    int spin_main=0, spin_neigh=0, bond=0;
    int index_neigh=0, index_neigh_=0;
    unsigned long long int my_one=1;
    int dE=0;

    int E_pos_pair_amount=0, M_pos_amount=0;
    unsigned long long int max_spin_code=ullint_2_pow_n(Ly*Lx)-1;

      for (unsigned long long int spin_code=0;  spin_code <= max_spin_code;  spin_code++)
      {
      //==DEBUG==//cout<<"DEBUG:  spin_code="<<spin_code<<endl;    //==DEBUG==//
        if ((lack_spins_in_lattice_mask & spin_code) == 0)
        {
        E_pos_pair_amount=0; M_pos_amount=0;

          for (int y=0;  y < Ly;  y++)
          {
            for (int x_=0;  x_ < Lx;  x_++)
            {
            int x=Lx-x_-1;
            //==DEBUG==//cout<<"DEBUG:  y="<<y<<"  x_="<<x_<<" my_one ="<<my_one<<endl;    //==DEBUG==// 

            {spin_main=-1;    if (((my_one << Lx*y+x_) & spin_code) != 0) {spin_main=+1;}    }

              //  right neigh
              if (x < Lx-1) 
              {  
              index_neigh=Lx*y+x+1;  index_neigh_=Lx*y+x_-1;
                if (x >= Lx-1) {index_neigh=Lx*y;  index_neigh_=Lx*y+Lx-1;}
                {spin_neigh=-1;    if (((my_one << index_neigh_) & spin_code) != 0) {spin_neigh=+1;}    }
                if (J_hor[(Lx+1)*y+x+1] == true) {bond=1;} else {bond=-1;}
                dE=bond*spin_main*spin_neigh;  if (dE > 0) {E_pos_pair_amount++;}    // cout<<"debug:  y="<<"   x="<<x<<"  energy add ver!"<<endl;    //==DEBUG==//
                //==DEBUG==//cout<<"  DEBUG:  "<<"  spin_main="<<spin_main<<" neigh hor spin="<<spin_neigh<<" bond="<<bond<<"  index_neigh_="<<index_neigh_<<endl;    //==DEBUG==//
              }
            
              //  up neigh
              if (y < Ly-1) 
              {          
              index_neigh=Lx*y+x+Lx;  index_neigh_=Lx*y+x_+Lx;
                if (y >= Ly-1) {index_neigh=x;  index_neigh_=x_;}
                {spin_neigh=-1;    if (((my_one << index_neigh_) & spin_code) != 0) {spin_neigh=+1;}    }
                if (J_ver[Lx*(y+1)+x] == true) {bond=1;} else {bond=-1;}  
              dE=bond*spin_main*spin_neigh;  if (dE > 0) {E_pos_pair_amount++;}  // cout<<"debug:  y="<<"   x="<<x<<"  energy add hor!"<<endl;    //==DEBUG==//
              //==DEBUG==//cout<<"  DEBUG:  "<<"  spin_main="<<spin_main<<" neigh ver spin="<<spin_neigh<<" bond="<<bond<<"  index_neigh_="<<index_neigh_<<endl;   //==DEBUG==//
              }
                      
            //  magnetiz
              if (spin_main == 1) {M_pos_amount++;}
            }
          }      //  for (int y=0;  y < Ly;  y++)


        // edges
        //  ---------------------------------- v51 additive ------------- V
          {
          //   if (((my_one << Lx*y+x_) & spin_code) != 0) {spin_main=+1;}
                                  
          unsigned long long int part_of_spin_lat_code=0, spin_code_temp=0;
          //  for (int edge_i=0;  edge_i < Lx;  edge_i++)
          //  {spin_main=-1;  if (((my_one << Lx*y+x_) & spin_code) != 0) {spin_main=+1;}  }
          //spin_code_temp=spin_code;
          //spin_code_temp=(spin_code_temp << (64-(Lx)));  spin_code_temp=(spin_code_temp >> (64-(Lx)));      
              
            for (int edge_i=0;  edge_i < Lx;  edge_i++)
            {my_one=1;  if (((my_one << (Lx-edge_i-1)) & spin_code)     != 0) {out_spin_array[edge_i]=true;} else {out_spin_array[edge_i]=false;}  }
              
            for (int edge_i=0;  edge_i < Lx;  edge_i++)
            {my_one=1;  if (((my_one << (Ly*Lx-edge_i-1)) & spin_code)  != 0) {out_spin_array[Lx+edge_i]=true;} else {out_spin_array[Lx+edge_i]=false;}  }
              
            for (int edge_i=0;  edge_i < Ly;  edge_i++)
            {my_one=1;  if (((my_one << (Lx*(edge_i+1)-1)) & spin_code) != 0) {out_spin_array[2*Lx+edge_i]=true;} else {out_spin_array[2*Lx+edge_i]=false;}  }
              
            for (int edge_i=0;  edge_i < Ly;  edge_i++)
            {my_one=1;  if (((my_one << (Lx*edge_i)) & spin_code)       != 0) {out_spin_array[2*Lx+Ly+edge_i]=true;} else {out_spin_array[2*Lx+Ly+edge_i]=false;}  }

          //--//  int dEE=amount_of_positiv_E_par_bool_edge_v3(out_J_array,  bool_edge_spin_exist_d_u_l_r_ar,  out_spin_array, 2*Lx+2*Ly);
          //--//  cout<<"dEE="<<dEE<<endl;
          E_pos_pair_amount+=amount_of_positiv_E_par_bool_edge_v3(out_J_array,  bool_edge_spin_exist_d_u_l_r_ar,  out_spin_array, 2*Lx+2*Ly);      //   !  !  !
                
          //part_of_spin_lat_code+=spin_code_temp;  
          //E_pos_pair_amount+=spin_code_temp;
          }
        //  ---------------------------------- v51 additive ------------- ^



        //==DEBUG==//cout<<"debug:  E_pos_pair_amount="<<E_pos_pair_amount<<"  M_pos_amount="<<M_pos_amount<<" E_pos_pair_amount*matr_check.size_x+M_pos_amount="<<(E_pos_pair_amount*matr_check.size_x+M_pos_amount)<<endl;
        mpz_add_ui(matr_check.element[E_pos_pair_amount*matr_check.size_x+M_pos_amount],  matr_check.element[E_pos_pair_amount*matr_check.size_x+M_pos_amount],  (unsigned int) 1);
        }      //  if ((lack_spins_in_lattice_mask & spin_code) == 0)


      }      //  for (unsigned long long int spin_code=0;  spin_code <= max_spin_code;  spin_code++)

    //matr_check.min_not_0_y_index=0;  matr_check.max_not_0_y_index=matr_check.size_y-1;  matr_check.not_0_len_y=matr_check.size_y;
    //matr_check.min_not_0_x_index=0;  matr_check.max_not_0_x_index=matr_check.size_x-1;  matr_check.not_0_len_x=matr_check.size_x;
    //cout<<"testing main matr:"<<endl;
    //matr_result_hist.display_on_screen_whole_matr((char *) "  ", 5);      //==DEBUG==//
    cout<<"right checked matr:"<<endl;
    matr_check.display_on_screen_whole_matr((char *) "  ", 5);            //==DEBUG==//
    cout<<"lowest row in matr points E="<<(-(2*Ly*Lx-Ly-Lx+2*Ly+2*Lx-amount_of_null_edge_bonds))<<",    dif E between rows equals 2"<<endl;  
    matr_check.free();
  
      if (out_spin_array) {delete[] out_spin_array;  out_spin_array=0;}
      if (out_J_array   ) {delete[] out_J_array;     out_J_array=0;}
        
    } 
  }






//  cleaning
//  ========
//
//

  if (min_E_array   ) {delete[] min_E_array;     min_E_array=0;  }



  if (mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT)
  {
    if (is_print_comments == true) 
    {    
    cout<<endl;
    cout<<endl<<"===============  ---------------------------    izing_2D_exact_gbc_EAnd_v54_minE_sr_mult_thrd_parallel   is launched    ------------------------   ============"<<endl;
    cout<<"===============  --------------------------          ==================================================================   ------------------------   ============"<<endl<<endl;  
    }
  }



return true;

}








bool swap_mem_paral(short int *min_E_array,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  unsigned long long int mem_per_thrd_in_unit_sint,  int thrd_num_1, int thrd_num_2,  int row,  int col, int thread_id_code_shift)
{

  if (min_E_array == 0) {return false;}


//std::memcpy(& min_E_array[thrd_num_1*mem_per_thrd_in_unit_sint+mem_part_third_in_thrd_in_unit_sint*2],  & min_E_array[thrd_num_1*mem_per_thrd_in_unit_sint+mem_part_third_in_thrd_in_unit_sint*1],  sizeof(short int)*mem_part_third_in_thrd_in_unit_sint);
//std::memcpy(& min_E_array[thrd_num_1*mem_per_thrd_in_unit_sint+mem_part_third_in_thrd_in_unit_sint*1],  & min_E_array[thrd_num_2*mem_per_thrd_in_unit_sint+mem_part_third_in_thrd_in_unit_sint*0],  sizeof(short int)*mem_part_third_in_thrd_in_unit_sint);
//std::memcpy(& min_E_array[thrd_num_2*mem_per_thrd_in_unit_sint+mem_part_third_in_thrd_in_unit_sint*0],  & min_E_array[thrd_num_1*mem_per_thrd_in_unit_sint+mem_part_third_in_thrd_in_unit_sint*2],  sizeof(short int)*mem_part_third_in_thrd_in_unit_sint);


  if (thrd_num_1 == mpirank) 
  {std::memcpy(& min_E_array[mem_part_third_in_thrd_in_unit_sint*2],  & min_E_array[mem_part_third_in_thrd_in_unit_sint*1],  sizeof(short int)*mem_part_third_in_thrd_in_unit_sint);}

  {  
  //  transfer:  min_E_array

    if (mpisize > 1)
    {

      if (mpirank == thrd_num_1)  
      {MPI_Recv(& min_E_array[mem_part_third_in_thrd_in_unit_sint*1],  mem_part_third_in_thrd_in_unit_sint, MPI_SHORT, thrd_num_2, 16+4*row*50*mpisize+thread_id_code_shift*2*50*mpisize+2*col*mpisize+2*thrd_num_1+0, MPI_COMM_WORLD, &Status);}  
      if (mpirank == thrd_num_2)  
      {MPI_Send(& min_E_array[mem_part_third_in_thrd_in_unit_sint*0],  mem_part_third_in_thrd_in_unit_sint, MPI_SHORT, thrd_num_1, 16+4*row*50*mpisize+thread_id_code_shift*2*50*mpisize+2*col*mpisize+2*thrd_num_1+0, MPI_COMM_WORLD);}  

      if (mpirank == thrd_num_2)  
      {MPI_Recv(& min_E_array[mem_part_third_in_thrd_in_unit_sint*0],  mem_part_third_in_thrd_in_unit_sint, MPI_SHORT, thrd_num_1, 16+4*row*50*mpisize+thread_id_code_shift*2*50*mpisize+2*col*mpisize+2*thrd_num_1+1, MPI_COMM_WORLD, &Status);}  
      if (mpirank == thrd_num_1)  
      {MPI_Send(& min_E_array[mem_part_third_in_thrd_in_unit_sint*2],  mem_part_third_in_thrd_in_unit_sint, MPI_SHORT, thrd_num_2, 16+4*row*50*mpisize+thread_id_code_shift*2*50*mpisize+2*col*mpisize+2*thrd_num_1+1, MPI_COMM_WORLD);}  

    }
    
  }
  
  


return true;

}














bool calc_min_E_partially_into_mem_chunk__for_given_first_down_row_and_col(short int *min_E_array,  int num_thrd,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  unsigned long long int amount_of_states_per_thrd,  int Ly,  int Lx,  bool * J_ver,  bool * J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, int amount_of_null_edge_bonds)
{

  if (bool_edge_spin_exist_d_u_l_r_ar == 0) {return false;}
  if (J_ver == 0) {return false;}
  if (J_hor == 0) {return false;}

unsigned long long int state_size_in_array=sizeof(unsigned long long int)/sizeof(min_E_array[0]);

unsigned long long int amount_of_states_per_thrd_div_2=amount_of_states_per_thrd/2;

//  starting main calc
//  ========
// 
//



//  calc 2, 3, 4, ..., last row         lattice energy (all lattice without 0 and 1 row).
//
int shE_0=0, shM_0=0;

  {
  bool * spin_array_loc=0;
  spin_array_loc=new bool[Lx*(Ly-1)];
    for (int i=0;  i < Lx*(Ly-1);  i++)
    {spin_array_loc[i]=false;}
  bool *J_ver_loc=0,  *J_hor_loc=0;
  J_ver_loc=J_ver+Lx*(sizeof(bool));  J_hor_loc=J_hor+(Lx+1)*(sizeof(bool));
  //  ---------------------------------- v51 modif ------------- V
  //
  bool *bool_edge_spin_exist_d_u_l_r_ar_loc=0;                                         //  temp edge array without, without the most down objects in left and right edges
  bool_edge_spin_exist_d_u_l_r_ar_loc=new bool[Lx*4+Ly*4-4];
    for (int i=0;  i < Lx*4;  i++)
    {bool_edge_spin_exist_d_u_l_r_ar_loc[i]=bool_edge_spin_exist_d_u_l_r_ar[i];}
    for (int i=0;  i < (Ly*2-2);  i++)
    {bool_edge_spin_exist_d_u_l_r_ar_loc[Lx*4+i]=bool_edge_spin_exist_d_u_l_r_ar[Lx*4+i+2];  bool_edge_spin_exist_d_u_l_r_ar_loc[Lx*4+Ly*2-2+i]=bool_edge_spin_exist_d_u_l_r_ar[Lx*4+Ly*2+i+2];}
  
  calc_pospair_and_upspin_EA_2D_lattice_with_edges_without_down_edge_v3(spin_array_loc,  bool_edge_spin_exist_d_u_l_r_ar_loc,  J_ver_loc,  J_hor_loc,  Ly-1,  Lx,  & shE_0);
  //==DEBUG==//  cout<<" shE_0="<<shE_0<<" shM_0="<<shM_0<<"  amount_of_null_edge_bonds="<<amount_of_null_edge_bonds<<endl;      //==DEBUG==//
  
    if (spin_array_loc) {delete[] spin_array_loc;  spin_array_loc=0;}
    if (bool_edge_spin_exist_d_u_l_r_ar_loc) {delete[] bool_edge_spin_exist_d_u_l_r_ar_loc;  bool_edge_spin_exist_d_u_l_r_ar_loc=0;}
  //
  //  ---------------------------------- v51 modif ------------- ^  
  } 


unsigned long long int max_row_state_local=num_thrd*amount_of_states_per_thrd+amount_of_states_per_thrd;
unsigned long long int row_0_state=0;      std::memcpy(& row_0_state,  & min_E_array[0], sizeof(unsigned long long int));
unsigned long long int row_0_state_index_in_array=state_size_in_array;



  for (unsigned long long int row_0_state_num=0;  row_0_state_num < amount_of_states_per_thrd;  row_0_state_num++,  row_0_state++,  row_0_state_index_in_array++)   
  {
  
    if (row_0_state_num == amount_of_states_per_thrd_div_2) 
    {//min_E_array[mem_part_third_in_thrd_in_unit_sint]=1;       
    row_0_state_index_in_array+=state_size_in_array;   
    std::memcpy(& row_0_state,  & min_E_array[mem_part_third_in_thrd_in_unit_sint], sizeof(unsigned long long int)); 
    //==DEBUG==//  cout<<"  row_0_state="<<row_0_state<<"  amount_of_states_per_thrd_div_2="<<amount_of_states_per_thrd_div_2<<" row_0_state_index_in_array="<<row_0_state_index_in_array<<endl;  //==DEBUG==//
    //==DEBUG==//  cout<<"  mem_part_third_in_thrd_in_unit_sint="<<mem_part_third_in_thrd_in_unit_sint<<endl;  //==DEBUG==//
    }
  
  
  //for (unsigned long long int row_0_state=num_thrd*amount_of_states_per_thrd;  row_0_state < max_row_state_local;  row_0_state++)   
  //{
  //==DEBUG==//  cout<<"DEBUG: start value: row_0_state="<<row_0_state<<"  ";  VDisplayULLIntToBin(row_0_state, Lx, true);  cout<<endl;  //==DEBUG==//  

  //  int posit_pair_number=0;
  //  int up_number=0;  
  int shE=0, shM=0;
  char var1=0, var2=0, var3=0;
  unsigned long long int cur_spin=1;   
 
    if ((row_0_state & cur_spin) != 0) {shM++;}

  //  we start with right spin     *--*--*
  //                                     ^
    for (int cell_num=0;  cell_num < Lx-1;  cell_num++)
    {
    var1=0; var2=0; var3=0;

    //  hor bond (left of main)
      if (J_hor[Lx-cell_num-1] == true) {var1=1;}
      if ((row_0_state & cur_spin) != 0) {var2=1;}    // spin 1
    cur_spin=(cur_spin << 1);
      if ((row_0_state & cur_spin) != 0) {var3=1; shM++;}    // spin 2
      if ((var1 ^ var2 ^ var3) != 0) {shE++;}      //==DEBUG==//  cout<<"shE="<<shE<<endl;      //==DEBUG==//

    //  (1)------------------------------- v51 additive ------------- V
    //  right down hor edge (corner)
      if ((cell_num == 0) && (bool_edge_spin_exist_d_u_l_r_ar[Lx*4+Ly*2+1] == true))     // interact  first col bottom row with abov spins
      {  
      char var6=0, var7=0;
        if (J_hor[Lx] == true) {var6=1;} else {var6=0;}
      //  var1                //  already defined - current left spin
        if (bool_edge_spin_exist_d_u_l_r_ar[Lx*4+Ly*2] == true) {var7=1;} else {var7=0;}
        if ((var2 ^ var6 ^ var7) != 0) {shE++;}      //==DEBUG==//  cout<<"shE="<<shE<<endl;      //==DEBUG==//
      }
    //  ---------------------------------- v51 additive ------------- ^  


    //  ver bond (above main)
      if (Ly > 1)   // interact  bottom row with abov spins excluding first col
      {
        if (J_ver[Lx+Lx-cell_num-1] == true) {var1=1;} else {var1=0;}
      int var4=0;    //  spin above in row 1
        if ((var1 ^ var2 ^ var4) != 0) {shE++;}      //==DEBUG==//  cout<<"shE="<<shE<<endl;      //==DEBUG==//
      }
      
      
    //  (2)------------------------------- v51 additive ------------- V  
    //  down edge,   ver bond
      if (bool_edge_spin_exist_d_u_l_r_ar[2*(Lx-cell_num-1)+1] == true)
      {
      char var5=0, var6=0;
        if (bool_edge_spin_exist_d_u_l_r_ar[2*(Lx-cell_num-1)] == true) {var5=1;}
        if (J_ver[(Lx-cell_num-1)] == true) {var6=1;}
        if ((var2 ^ var5 ^ var6) != 0) {shE++;}      //==DEBUG==//  cout<<"shE="<<shE<<endl;      //==DEBUG==//
      }
    //  ---------------------------------- v51 additive ------------- ^  
      

      if ((Ly > 1) && (cell_num == Lx-2))     // interact  first col bottom row with abov spins
      {  
        if (J_ver[Lx+Lx-cell_num-2] == true) {var1=1;} else {var1=0;}
      int var4=0;    //  spin above  (the first left spin in row 0)
      //  var3                //  already defined - current left spin
        if ((var1 ^ var3 ^ var4) != 0) {shE++;}      //==DEBUG==//  cout<<"shE="<<shE<<endl;      //==DEBUG==//
      }


    //  (3)--------------------------------- v51 additive ------------- V
    //  left down ver edge (corner)
      if ((cell_num == Lx-2) && (bool_edge_spin_exist_d_u_l_r_ar[1] == true))     // interact  first col bottom row with abov spins
      {
        if (J_ver[0] == true) {var1=1;} else {var1=0;}
      //  var3                //  already defined - current left spin
        if (bool_edge_spin_exist_d_u_l_r_ar[0] == true) {var2=1;} else {var2=0;}
        if ((var1 ^ var2 ^ var3) != 0) {shE++;}      //==DEBUG==//  cout<<"shE="<<shE<<endl;      //==DEBUG==//
      }
    //  ---------------------------------- v51 additive ------------- ^  


    //  (4)--------------------------------- v51 additive ------------- V
    //  left down hor edge (corner)
      if ((cell_num == Lx-2) && (bool_edge_spin_exist_d_u_l_r_ar[Lx*4+1] == true))     // interact  first col bottom row with abov spins
      {  
        if (J_hor[0] == true) {var1=1;} else {var1=0;}
      //  var3                //  already defined - current left spin
        if (bool_edge_spin_exist_d_u_l_r_ar[Lx*4] == true) {var2=1;} else {var2=0;}
        if ((var1 ^ var2 ^ var3) != 0) {shE++;}      //==DEBUG==//  cout<<"shE="<<shE<<endl;      //==DEBUG==//
      }
    //  ---------------------------------- v51 additive ------------- ^  


    }
    

  min_E_array[row_0_state_index_in_array]=shE_0+shE;   //==DEBUG==//==DEBUG==//if ((row_0_state) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;}            //==DEBUG==//==DEBUG==//
  //==DEBUG==//  matr_hist[row_0_state].display_on_screen_whole_matr();  //==DEBUG==//  
  //==DEBUG==//  cout<<"row_0_state="<<row_0_state<<"  shE_0="<<(shE_0)<<"  "<<"  shE="<<(shE)<<"  "<<"min_E_array[row_0_state_index_in_array]="<<min_E_array[row_0_state_index_in_array]<<endl;    //==DEBUG==//    
  
  }      //  for (unsigned long long int row_0_state_num=0;  row_0_state_num < amount_of_states_per_thrd;  row_0_state_num++,  row_0_state++)    

//==DEBUG==//  cout<<endl;    //==DEBUG==//  




return true;

}










bool calc_min_E_partially_into_mem_chunk__for_given_row_and_col(short int *min_E_array,  unsigned long long int amount_of_states_per_thrd_div_2,  unsigned long long int mem_part_third_in_thrd_in_unit_sint,  int amount_of_steps_into_thread, int row,  int col,  int Ly,  int Lx,  bool *& J_ver,  bool *& J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, short int *result_in_thrd_min_E,  unsigned long long int *up_row_gs_state_in_thrd, int amount_of_null_edge_bonds)
{


unsigned long long int state_size_in_array=sizeof(unsigned long long int)/sizeof(min_E_array[0]);

short int result_min_E_pospair=0;
unsigned long long int start_state_1=0,  start_state_2=0;

  //for (int row=1;  row < Ly;  row++) 
  {
  //==DEBUG==//  cout<<"DEBUG:  ==============================================   row="<<row<<endl;      //==DEBUG==//
   
  
  unsigned long long int shifting_one=1;
  bool J[4];
  char s[4];
  char cur_spin=0;
  int shE_matr1=0, shE_matr2=0;


  //  unsigned long long int start_state_1=0,  start_state_2=0;      //  already done above
  unsigned long long int local_number_of_state_1=0,  local_number_of_state_2=0;  
  unsigned long long int cur_state1_index_in_array=state_size_in_array, cur_state2_index_in_array=state_size_in_array;  
  int row_state_range_out=amount_of_states_per_thrd_div_2,  row_state_range_in=2,  row_state_range_in_div_2=1;
  std::memcpy(& start_state_1,  & min_E_array[0],                                   sizeof(unsigned long long int));
  std::memcpy(& start_state_2,  & min_E_array[mem_part_third_in_thrd_in_unit_sint], sizeof(unsigned long long int));
  
    for (int i=0;  (i < col) && (i < amount_of_steps_into_thread-1);  i++)                                                        //  !  !  !       !   IMPORTANT STEP,   DANGER  !       !  !  !  
    {row_state_range_out/=2;  row_state_range_in_div_2*=2;  row_state_range_in*=2;}                                               //  !  !  !       !   IMPORTANT STEP,   DANGER  !       !  !  !  

  //==DEBUG==// cout<<"DEBUG: "<<" row_state_range_out="<<row_state_range_out<<" row_state_range_in_div_2="<<row_state_range_in_div_2<<"  row_state_range_in="<<row_state_range_in<<"  start_state_1="<<start_state_1<<"  start_state_2="<<start_state_2<<endl;
    
  unsigned long long int cur_state1=0, cur_state2=0;      // down  and  up  states

    //for (int col=0;  col < Lx;  col++)
    { 
    //==DEBUG==//  cout<<"DEBUG:  ======================     col="<<col<<endl;      //==DEBUG==//  
    shifting_one=1;
    shifting_one=(shifting_one << col);

      for (unsigned long long int index_out=0;  index_out < row_state_range_out;  index_out++)
      {
      //==DEBUG==//  cout<<"DEBUG:     index_out="<<index_out<<endl;      //==DEBUG==//  
              
        for (unsigned long long int index_in=0;  index_in < row_state_range_in_div_2;  index_in++)
        {
        //==DEBUG==//  cout<<"DEBUG:     index_in="<<index_in<<endl;      //==DEBUG==//
        
        local_number_of_state_1=index_out*row_state_range_in+index_in;        local_number_of_state_2=local_number_of_state_1+row_state_range_in_div_2;
        
          if (local_number_of_state_1 >= amount_of_states_per_thrd_div_2) 
          {local_number_of_state_1-=amount_of_states_per_thrd_div_2;  cur_state1=start_state_2+local_number_of_state_1;  cur_state1_index_in_array=mem_part_third_in_thrd_in_unit_sint+4+local_number_of_state_1;} 
          else {cur_state1=start_state_1+local_number_of_state_1;  cur_state1_index_in_array=state_size_in_array+local_number_of_state_1;}
          
          if (local_number_of_state_2 >= amount_of_states_per_thrd_div_2) 
          {local_number_of_state_2-=amount_of_states_per_thrd_div_2;  cur_state2=start_state_2+local_number_of_state_2;  cur_state2_index_in_array=mem_part_third_in_thrd_in_unit_sint+4+local_number_of_state_2;} 
          else {cur_state2=start_state_1+local_number_of_state_2;  cur_state2_index_in_array=state_size_in_array+local_number_of_state_2;}

        //==DEBUG==//  cout<<"DEBUG: "<<"  local_number_of_state_1="<<local_number_of_state_1<<" cur_state1="<<cur_state1<<" cur_state2="<<cur_state2<<"    "<<endl;
        //==DEBUG==//  cout<<"DEBUG: "<<"  local_number_of_state_2="<<local_number_of_state_2<<" cur_state1_index_in_array="<<cur_state1_index_in_array<<"  cur_state2_index_in_array="<<cur_state2_index_in_array<<endl;

          
        
        //==DEBUG==//  cout<<"index_in="<<index_in<<endl;      //==DEBUG==//
        //cur_state1=start_state_1+index_out*row_state_range_in+index_in;   cur_state2=cur_state1+row_state_range_in_div_2;
        //==DEBUG==//  cout<<"DEBUG:     cur_state1="<<cur_state1<<"  cur_state2="<<cur_state2<<endl;      //==DEBUG==//  

        // --------\    
        //          -------------v        
        // neigbours analising
        // input  data:  cur_state1, cur_state2,  shifting_one, col, rowJ_ver[..]       
        // output data:  s[..], J[..]
        
        // left neigh
        J[2]=J_hor[(Lx+1)*row+Lx-col-1  ];
          if (col < Lx-1) 
          {
          //shifting_one=(shifting_one << 1);  if ((cur_state1 & shifting_one) != 0) {s[2]=1;} else {s[2]=-1;}  shifting_one=(shifting_one >> 1);
          s[2]=-1;
          }
          else {  if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+row*2+1] == true) {if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+row*2] == true) {s[2]=1;} else {s[2]=-1;}} else {s[2]=0;}  }      //  edge left            //  -------- v51 modif ----------- < 

        // right neigh
        J[0]=J_hor[(Lx+1)*row+Lx-col];
          if (col > 0) 
          {
          shifting_one=(shifting_one >> 1);  if ((cur_state1 & shifting_one) != 0) {s[0]=1;} else {s[0]=-1;}  shifting_one=(shifting_one << 1);
          }
          else {  if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*Ly+row*2+1] == true) {if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*Ly+row*2] == true) {s[0]=1;} else {s[0]=-1;}} else {s[0]=0;}  }      //  edge right     //  -------- v51 modif ----------- <

        // down
        J[3]=J_ver[Lx*row+Lx-col-1 ];

        //  up
        //  ---------------------------------- v51 modif ------------- V
        //        
        J[1]=J_ver[Lx*row+Lx+Lx-col-1 ];
          if (row+1 < Ly) { s[1]=-1;} else {  if (bool_edge_spin_exist_d_u_l_r_ar[2*Lx+(Lx-col-1)*2+1] == true) {if (bool_edge_spin_exist_d_u_l_r_ar[2*Lx+(Lx-col-1)*2] == true) {s[1]=1;} else {s[1]=-1;}} else {s[1]=0;}  }      //  edge up
        //
        //  ---------------------------------- v51 modif ------------- ^
          

        s[3]=+1;
        //==DEBUG==//    if (s[0] == -1) {cout<<"s[0]=dn ";} else {  if (s[0] == +1) {cout<<"s[0]=up ";} else {cout<<"s[0]=0 ";}  }
        //==DEBUG==//    if (s[1] == -1) {cout<<"s[1]=dn ";} else {  if (s[1] == +1) {cout<<"s[1]=up ";} else {cout<<"s[1]=0 ";}  }
        //==DEBUG==//    if (s[2] == -1) {cout<<"s[2]=dn ";} else {  if (s[2] == +1) {cout<<"s[2]=up ";} else {cout<<"s[2]=0 ";}  }
        //==DEBUG==//    if (s[3] == -1) {cout<<"s[3]=dn ";} else {  if (s[3] == +1) {cout<<"s[3]=up ";} else {cout<<"s[3]=0 ";}  }
        //==DEBUG==//  cout<<"      ";
        //==DEBUG==//    if (J[0] == false) {cout<<"J[0]=-1 ";} else {cout<<"J[0]=+1 ";}
        //==DEBUG==//    if (J[1] == false) {cout<<"J[1]=-1 ";} else {cout<<"J[1]=+1 ";}
        //==DEBUG==//    if (J[2] == false) {cout<<"J[2]=-1 ";} else {cout<<"J[2]=+1 ";}
        //==DEBUG==//    if (J[3] == false) {cout<<"J[3]=-1 ";} else {cout<<"J[3]=+1 ";}
        //==DEBUG==//  cout<<"      "; 
        //shE_matr2 =pospair_amount_1_cell_EA(true  , J, s);
        //cout<<"aa"<<J[0]<<J[1]<<J[2]<<J[3]<<"aa";
        int E11=pospair_amount_1_cell_EA(false  , J, s);
        //==DEBUG==//  cout<<" ex 1: E11="<<E11<<"  ";  //==DEBUG==//  
        //==DEBUG==//shE_matr2-=pospair_amount_1_cell_EA(false , J, s);  //==DEBUG==//
        int E12=pospair_amount_1_cell_EA(true , J, s);
        //==DEBUG==//  cout<<"   ex 2: E12="<<E12<<"  ";  //==DEBUG==//
        shE_matr2=E12-E11;
        //==DEBUG==//  cout<<endl;
        //matr_gather_hist.display_on_screen_whole_matr((char *) "  ", 2);
        //==DEBUG==//  if ((col >= 0) && (index_in >= 0)) {cout<<"DEBUG PAIR: until:"<<endl;  matr_hist[row_state_1].display_on_screen_whole_matr();  matr_hist[row_state_2].display_on_screen_whole_matr();} 
        // matr_cross_add_v3(& matr_hist[row_state_1],  & matr_hist[row_state_2],  shE);     //  !!! ======== v3 M_div_2 line ======== !!!
        //==DEBUG==//  if ((col >= 0) && (index_in >= 0)) {cout<<"DEBUG PAIR: after:"<<endl;  matr_hist[row_state_1].display_on_screen_whole_matr();  matr_hist[row_state_2].display_on_screen_whole_matr();}

        s[3]=-1;
        int E21=pospair_amount_1_cell_EA(false , J, s);
        //==DEBUG==//cout<<" ex 1: E21="<<E21<<"  ";  //==DEBUG==//  
        int E22=pospair_amount_1_cell_EA(true  , J, s);
        //==DEBUG==//  cout<<"   ex 2: E22="<<E22<<"  ";  //==DEBUG==//  
        shE_matr1=E22-E21;
        //==DEBUG==//  cout<<endl;
        // shE_matr1 =pospair_amount_1_cell_EA(true  , J, s);
        // shE_matr1-=pospair_amount_1_cell_EA(false , J, s);

        //==DEBUG==//  cout<<"DEBUG:     shE_matr1="<<shE_matr1<<"  shE_matr2="<<shE_matr2<<endl;      //==DEBUG==// 
        

        //
        // neigbours analising
        //          -------------^
        // --------/        
          //==DEBUG==//  if ((col >= 0) && (index_in >= 0)) 
          //==DEBUG==//  {cout<<"DEBUG PAIR: until:"<<endl;  matr_hist[cur_state1].display_on_screen_whole_matr((char *) "  ", 2);  matr_hist[cur_state2].display_on_screen_whole_matr((char *) "  ", 2);} 
          //matr_cross_add_v4(& matr_hist[cur_state1],  & matr_hist[cur_state2],  shE_matr1,  shE_matr2);     //  !!! ======== v3 M_div_2 line ======== !!!
      

        //cur_state1         shift_index_old   shift_index_neo
        short int inst_1_copy_min_E=min_E_array[cur_state1_index_in_array];    
        
          if (min_E_array[cur_state2_index_in_array] < min_E_array[cur_state1_index_in_array]) 
          {
          min_E_array[cur_state1_index_in_array]=min_E_array[cur_state2_index_in_array]; 
            //==DEBUG==//==DEBUG==//if ((cur_state2) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;} 
            //==DEBUG==//==DEBUG==//if ((cur_state1) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;}             
          } //  processing inst 1
          else 
          {  
            if (min_E_array[cur_state2_index_in_array] == min_E_array[cur_state1_index_in_array]) 
            {
            min_E_array[cur_state1_index_in_array]=min_E_array[cur_state2_index_in_array];    
              //==DEBUG==//==DEBUG==//if ((cur_state1) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;}                  
              //==DEBUG==//==DEBUG==//if ((cur_state1) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;}            //==DEBUG==//==DEBUG==//   
              //==DEBUG==//==DEBUG==//if ((cur_state2) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;}            //==DEBUG==//==DEBUG==//   
              
              
                                                                                                   //  insurance for out of ullint,   it can be commented
            }
          } 
        short int inst_2_min_E_variant1=inst_1_copy_min_E+shE_matr1;  // short int inst_2_min_E_g_variant1=inst_1_copy_min_E_g;
        short int inst_2_min_E_variant2=min_E_array[cur_state2_index_in_array]+shE_matr2;  // short int inst_2_min_E_g_variant2=min_E_g_array[cur_state2]+shE_matr2;

          if (inst_2_min_E_variant1 < inst_2_min_E_variant2) {min_E_array[cur_state2_index_in_array]=inst_2_min_E_variant1; } 
          else 
          {
            
            if (inst_2_min_E_variant1 == inst_2_min_E_variant2) 
            {                                                                                          
              //==DEBUG==//==DEBUG==//if ((cur_state2) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;}            //==DEBUG==//==DEBUG==//        
              //==DEBUG==//==DEBUG==//if ((cur_state2) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;}            //==DEBUG==//==DEBUG==//                      
                                                                              //  insurance for out of ullint,   it can be commented
            }            
          min_E_array[cur_state2_index_in_array]=inst_2_min_E_variant2;
            //==DEBUG==//==DEBUG==//if ((cur_state2) >= max_row_state) {cout<<"ERROR!!"<<endl;  char ch;  cin>>ch;}            //==DEBUG==//==DEBUG==//
          }     //)))))))))))))))))))))))))))//  if (cur_state2 == 3) {cout<<"min_E_array[3]="<<min_E_array[3]<<endl;}
                                
          
        }

      }


    //==DEBUG==//  cout<<endl;   //==DEBUG==//
    //==DEBUG==//    for (unsigned long long int state_temp=0;  state_temp < max_row_state ;  state_temp++)  //==DEBUG==//
    //==DEBUG==//    {cout<<"row="<<row<<"  "<<"min_E_array[shift+"<<state_temp<<"]="<<min_E_array[state_temp]<<endl;}  //==DEBUG==//
    //==DEBUG==//  cout<<endl;  //==DEBUG==//
      
    //row_state_range_out/=2;  row_state_range_in_div_2*=2;  row_state_range_in*=2;
    //
    }      //  for (int col=0;  col < Lx;  col++)


  //  ----------------
  //
  //)________________)____________________________________)
  //
  
 
  //  ==================================================================
  //  search for minimum for the last row
  //
    if ((row+1 == Ly) && (col == Lx-1) && (result_in_thrd_min_E))
    {
    //==DEBUG==// cout<<"debug55: "<<"  "<<"mpirank="<<mpirank<<"  "<<"min_E_array[0]="<<min_E_array[0]<<"  "<<"min_E_array[5]="<<min_E_array[5]<<"  "<<"min_E_array[4]="<<min_E_array[4]<<"  "<<"min_E_array[9]="<<min_E_array[9]<<"  "<<endl; //==DEBUG==//
    
    unsigned long long int up_row_gs_state_in_thrd_candidate=start_state_1;  
    unsigned long long int index_in_array=state_size_in_array;
    result_min_E_pospair=min_E_array[index_in_array];    //if (mpirank == 0) {cout<<"VVVVVVVVV  result_min_E_pospair1="<<result_min_E_pospair<<endl;}
    index_in_array++;
  
      for (unsigned long long int row_state_num=1;  row_state_num < amount_of_states_per_thrd_div_2*2;  row_state_num++,  index_in_array++)
      {
      up_row_gs_state_in_thrd_candidate++;
        if (row_state_num == amount_of_states_per_thrd_div_2)  {index_in_array+=state_size_in_array;  up_row_gs_state_in_thrd_candidate=start_state_2;}
              
                            
        if (min_E_array[index_in_array] < result_min_E_pospair) 
        {  //  if (mpirank == 0) {cout<<"VVVVVVVVV  result_min_E_pospair2="<<result_min_E_pospair<<endl;}
        result_min_E_pospair=min_E_array[index_in_array];
          if (up_row_gs_state_in_thrd) {*up_row_gs_state_in_thrd=up_row_gs_state_in_thrd_candidate;}            
        }
        
      }      //  for (unsigned long long int row_state=1;  row_state < max_row_state;  row_state++)  

      if (result_in_thrd_min_E) {*result_in_thrd_min_E=result_min_E_pospair;} 
    }      //  if (row+1 == Ly)


    //if (is_print_comments == true) {print_progress(row, mustachieve_progress_perc, Ly-1, out_progress);}  //cout<<" min_E_array[0]="<<min_E_array[3]<<endl;

  }      //  for (int row=1;  row < Ly;  row++)


  


return true;

}

























//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *


//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *


//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *


//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *















int main(int argn,  char **argc)
{

// if (mpirank == 0) {cout<<"argn="<<argn<<"   ";  for (int i=0;  i < argn;  i++) {cout<<"argc["<<i<<"]="<<argc[i]<<"     ";}   cout<<endl;}



//MPI_Status Status;                               // ====================  MPI
//int mpisize, mpirank;                            // ====================  MPI
MPI_Init(& argn, & argc);                          // ====================  MPI           //  ************ Start MPI ***********
MPI_Comm_size(MPI_COMM_WORLD, &mpisize);           // ====================  MPI
MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);           // ====================  MPI


//  check mpisize parameter
//
  if (mpisize > 512) {cout<<"Error! Amount of threads must less or equal 512!"<<endl;  MPI_Finalize(); return 0;}
  
  if (mpisize > 2)
  {
  unsigned long long int my_power=1; 
    for (unsigned long long int i=0;  i < 512;  i++)
    {my_power=my_power*2;  if (mpisize == my_power) {break;}   if (mpisize < my_power) {cout<<"Error! Amount of threads must be power of 2!"<<endl;  MPI_Finalize(); return 0;}  }
  }





//mpf_set_default_prec(100);

  if (mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT)
  {
  cout<<endl<<"=============================================================== start =============================================================== "<<endl;
  cout<<endl; 
  }
  
clock_t start=0, end=0;  
start = clock();
srand( (unsigned) time(NULL) ); 






//int Ly=__MACR_LATTICE_LY, Lx=__MACR_LATTICE_LX;
int Ly=0, Lx=0;
int J_plus_perc= 50, J_min_perc= 50;
            
bool *J_ver=0; 
bool *J_hor=0;

bool * bool_edge_spin_exist_d_u_l_r_ar=0;
bool *  bool_one_up_row_gs_state_arg=0;


//  ____________________________________________________v
//
//
bool is_print_comments=false; 
bool print_progress=true;       if (__MACR_PRINT_PROGRESS != 1) {print_progress=false;}
bool do_enum_check=false;    
short int  result_min_E_arg=0; 
unsigned long long int  one_up_row_gs_state_arg=0;
bool check_is_ok_arg=true; 
bool is_to_do_record_in_file=true;    if (__MACR_TO_DO_RECORD_INTO_FILE != 1) {is_to_do_record_in_file=false;}
bool to_do_extract_initial_data_from_file=true;      if (__MACR_TO_EXTRACT_DATA_FROM_FILE != 1) {to_do_extract_initial_data_from_file=false;}
bool is_to_record_gs_up_row=true;        
bool to_do_extract_up_row_gs=true; 

char filename_to_extract_data[127];  filename_to_extract_data[0]='\0';

//
//  ____________________________________________________^




// ----------------------------//   v
//


bool is_error_in_command_line_arguments=false;

  if (mpirank == 0) 
  {  


    if (is_error_in_command_line_arguments  == false)
    {      
      if ((argn != 4) && (argn != 6)) 
      {is_error_in_command_line_arguments=true;  cout<<"Error in command line arguments!  At least 3 first arguments must be."<<endl;}      
    }


    if (is_error_in_command_line_arguments  == false)
    {      
      if ((argn > 4) && (strlen(argc[4]) < 1)) 
      {is_error_in_command_line_arguments=true;  cout<<"Error in command line arguments! 4th argument must be equaled 1 or 0  if  you want to extract data from file."<<endl;  }      
    }

  
    if ((is_error_in_command_line_arguments  == false))
    {      
      if (((argn > 5) && (strlen(argc[5]) < 1)) || (argn == 5))
      {
      is_error_in_command_line_arguments=true;
      cout<<"Error in command line arguments! Last 5th argument must be not null variable with name of file if you want extract data from file!"<<endl;  
      }      
    }


    if ((is_error_in_command_line_arguments  == false) && (argn == 4))
    {     
      if ((strlen(argc[1]) < 1) || (strlen(argc[2]) < 1) || (strlen(argc[3]) < 1))
      {
      is_error_in_command_line_arguments=true;
      cout<<"Error in command line arguments! One of three arguments is empty char string!"<<endl; 
      }    
    }
  
    if ((is_error_in_command_line_arguments  == false) && (argn == 6))
    {     
      if ((strlen(argc[4]) < 1) || (strlen(argc[5]) < 1))
      {
      is_error_in_command_line_arguments=true;
      cout<<"Error in command line arguments! One of two last arguments is empty char string!"<<endl;  
      }    
    } 


    if ((is_error_in_command_line_arguments  == false) && (argn == 6))  
    {
      if (argc[4][0] == '1') {to_do_extract_initial_data_from_file=true;} else {to_do_extract_initial_data_from_file=false;}
      
      for (int i=0;  i < strlen(argc[5]);  i++)
      {filename_to_extract_data[i]=argc[5][i];}
       
    filename_to_extract_data[strlen(argc[5])]='\0';   
    }    
   
   
    if ((is_error_in_command_line_arguments  == false) && ((argn == 4) || (argn == 6)))
    {    
      if (argn == 4) {to_do_extract_initial_data_from_file=false;}
      
      if ((argn == 4) || (argc[4][0] != '1'))
      {
      Ly=atoi(argc[1]);  Lx=atoi(argc[2]);  J_plus_perc=atoi(argc[3]);    J_min_perc=100-J_plus_perc;
  
        if (Ly < 2) {is_error_in_command_line_arguments=true;  cout<<"Error in command line arguments! Ly must be greater then 1!"<<endl;}
        if ((Lx < 2) || (Lx > 100)) {is_error_in_command_line_arguments=true;  cout<<"Error in command line arguments! Lx must be:   1 < Lx < 100    !"<<endl;}
        if ((J_plus_perc < 0) || (J_plus_perc > 100)) {is_error_in_command_line_arguments=true;  cout<<"Error in command line arguments!  Fraction of positive bond bonds (J_plus_perc) must be  must be:   0 <= J_plus_perc <= 100    !"<<endl;}
      }      
    }   
    


  char chstr_addit_1_to_error_comments[]="  argument format:   1) *argc[1]: Ly    2) *argc[2]: Lx    3) *argc[3]: J_bond_plus_perc    4) *argc[4]:  bool__1_data_from_command_line_args__0_from_file    5) *argc[5]:  chstr_filename\n";
  char chstr_addit_2_to_error_comments[]="  if (bool__1_data_from_command_line_args__0_from_file == '1') {all data will be taken from file, 3 first args are ignored}\n";
  char chstr_addit_3_to_error_comments[]="  if (bool__1_data_from_command_line_args__0_from_file != '1') {all data will be generated automatically, 5 arg is ignored}\n";
  char chstr_addit_4_to_error_comments[]="  it's allowed to set only 3 first arguments as well\n";
  char chstr_addit_5_to_error_comments[]="  Example 1 of launching using 4 threads on linux with generating J bonds (+50%, -50%) for 10x10 lattice: \n";
  char chstr_addit_6_to_error_comments[]="  mpirun -n 4  prg.exe  10 10 50\n";
  char chstr_addit_7_to_error_comments[]="  Example 2 of launching using 4 threads on linux with reading data from file with name 'myfile.txt',  3 first args are ingnored: \n";
  char chstr_addit_8_to_error_comments[]="  mpirun -n 4  prg.exe  10 10 50  1  myfile.txt\n";


    if (is_error_in_command_line_arguments  == true) 
    {
    cout<<chstr_addit_1_to_error_comments;    cout<<chstr_addit_2_to_error_comments;    cout<<chstr_addit_3_to_error_comments;    cout<<chstr_addit_4_to_error_comments;
    cout<<chstr_addit_5_to_error_comments;    cout<<chstr_addit_6_to_error_comments;    cout<<chstr_addit_7_to_error_comments;    cout<<chstr_addit_8_to_error_comments;
    cout<<endl<<endl;
    }

  }
  
//
// ----------------------------//   ^









//    ----    ----    ----    ----    generation in data  ----    ----    ----    ---- 
//     V  V  V
//

  if ((mpirank == 0) && (is_error_in_command_line_arguments  == false) && (to_do_extract_initial_data_from_file == false))
  {
  
    if (to_do_extract_initial_data_from_file == false)
    {
    J_ver=new bool[Lx*Ly+Lx];
    J_hor=new bool[Lx*Ly+Ly];   


      for (int i=0;  i < Ly*Lx+Lx;  i++)
      {J_ver[i]=false;}

      for (int i=0;  i < Ly*Lx+Ly;  i++)
      {J_hor[i]=false;}                      
    

    bool_edge_spin_exist_d_u_l_r_ar=new bool[4*Ly+4*Lx];

      {
      bool ok=generate_2D_lattice_J_v3(J_ver,  J_hor,  Ly,  Lx, J_plus_perc, J_min_perc);
        if (ok == false) {return false;}        
      }    



      for (int i=0;  i < (2*Ly+2*Lx);  i++)
      {
      bool_edge_spin_exist_d_u_l_r_ar[2*i]=false;    
        //if ((rand() % 2) == 0) 
        //{bool_edge_spin_exist_d_u_l_r_ar[2*i]=true;} 
        //else 
        //{bool_edge_spin_exist_d_u_l_r_ar[2*i]=false;}  
      }    

      for (int i=0;  i < (2*Ly+2*Lx);  i++)
      {
      //bool_edge_spin_exist_d_u_l_r_ar[2*i+1]=true;    
        //if ((rand() % 2) == 0) 
        //{bool_edge_spin_exist_d_u_l_r_ar[2*i+1]=true;} 
        //else 
        {bool_edge_spin_exist_d_u_l_r_ar[2*i+1]=false;}  
      } 

  
    rand_generate_edges_spin_val(bool_edge_spin_exist_d_u_l_r_ar, Ly, Lx);
    }
    
  }   

//
//     ^  ^  ^
//    ----    ----    ----    ----    generation in data  ----    ----    ----    ---- 


//print_bool_edge_spin_exist_d_u_l_r_ar_for_debug(Ly,  Lx, bool_edge_spin_exist_d_u_l_r_ar);

  
   




//  ..............................  v  v  v  v   .....................................  
//  ..............................  v  v  v  v   .....................................  
//  extract data from file
//
bool is_lat_info_extracted_ok=false;

  if (mpirank == 0)
  {

    if ((to_do_extract_initial_data_from_file == true) && (is_error_in_command_line_arguments  == false))
    {
    int min_energ=0;
    //char filename_to_extract_data[]="Ising_EA_2D_gbc_bond_pl50_min50_0020x0020.txt";
    bool is_up_gs_found=false;
    bool is_edges_extracted=false;

    // --  -- prototype --  -- //
    // 
    //  bool extract_lattice_gbc_data_v54(char * filename, int &Ly, int &Lx, bool * & J_ver,  bool * &J_hor,  int  & J_plus_perc, int & J_min_perc, bool & is_edges_extracted, bool & is_up_gs_found, bool *& bool_edge_spin_exist_d_u_l_r_ar,  bool * bool_up_row_gs_ar,  int *min_energ, bool pbc1fbc0,  bool is_print_comments,  bool is_to_extract_up_row_gs);
    //
  
      if (bool_edge_spin_exist_d_u_l_r_ar) {delete[] bool_edge_spin_exist_d_u_l_r_ar;  bool_edge_spin_exist_d_u_l_r_ar=0;}    
      if (J_ver) {delete[] J_ver;  J_ver=0;}
      if (J_hor) {delete[] J_hor;  J_hor=0;}    
  
    is_lat_info_extracted_ok=extract_lattice_gbc_data_v54(filename_to_extract_data, Ly, Lx,  J_ver, J_hor,  J_plus_perc, J_min_perc,  is_edges_extracted, is_up_gs_found, bool_edge_spin_exist_d_u_l_r_ar,  bool_one_up_row_gs_state_arg,  & min_energ, false, is_print_comments, to_do_extract_up_row_gs);      //  !  !  !
  
  
      if (is_to_record_gs_up_row == true) 
      {
      cout<<"is_up_gs_found="<<((int) is_up_gs_found)<<"  "<<"bool_one_up_row_gs_state_arg: "<<bool_one_up_row_gs_state_arg<<endl;
        if (bool_one_up_row_gs_state_arg == 0) {bool_one_up_row_gs_state_arg=new bool[Lx]; }
      
        if ((bool_one_up_row_gs_state_arg) && (is_up_gs_found == false))
        {
          for (int i=0;  i < Lx;  i++)
          {bool_one_up_row_gs_state_arg[i]=false;}
        }
      
      }
  
  
      if (is_edges_extracted == false) {is_lat_info_extracted_ok=false;}
      
      if (is_print_comments == true)
      {  
      cout<<"lattice info is extracted successfully: ";    
        if (is_lat_info_extracted_ok == true) {cout<<"yes"<<endl;} else {cout<<"no"<<endl;}
      }


      if (is_lat_info_extracted_ok == false) {Ly=0;  Lx=0;}
      
    }      //  if ((to_do_extract_initial_data_from_file == true) && (is_error_in_command_line_arguments  == false))
  
  //  cout<<endl<<"print bolls ar: "<<endl;  print_bool_edge_spin_exist_d_u_l_r_ar_for_debug(Ly,  Lx, bool_edge_spin_exist_d_u_l_r_ar); cout<<end<<endl;      //==DEBUG==//


  } 

//
//
//  ..............................   ^  ^  ^  ^   .....................................  
//  ..............................   ^  ^  ^  ^   ..................................... 





//____II____II____II____II____II____II____II____II____II  
//
//  *  *  *    *  *  *  !  !  !    *  *  *    *  *  *      
//  -------------------  prototype   -------------------
// 


//  bool izing_2D_exact_gbc_EAnd_v54_minE_sr_mult_thrd_emul(int Ly,  int Lx,  bool *J_ver,  bool *J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, int posit_bond_perc_arg,  int negat_bond_perc_arg, bool is_to_do_record_in_file=true,  bool is_print_comments=true, bool do_enum_check=true, short int * result_min_E_arg=0,  bool *bool_one_up_row_gs_state=0,  bool *check_is_ok_arg=0, bool  print_progress=true)

bool min_E_multthread_emul_done_ok=false;

  if  ((mpirank == 0) && (bool_one_up_row_gs_state_arg == 0) && (Lx > 0))  {bool_one_up_row_gs_state_arg=new bool[Lx];}

min_E_multthread_emul_done_ok=izing_2D_exact_gbc_EAnd_v55_minE_sr_parallel_export_variant(Ly,  Lx,  J_ver,  J_hor,  bool_edge_spin_exist_d_u_l_r_ar, J_plus_perc, J_min_perc,  is_to_do_record_in_file,  is_print_comments, do_enum_check, & result_min_E_arg, bool_one_up_row_gs_state_arg, & check_is_ok_arg,  print_progress);


  if (mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT)
  {  
  cout<<"............................................................"<<endl;
    if (min_E_multthread_emul_done_ok == true)
    {
    cout<<"calculation is completed successfully"<<endl; 
    cout<<"min_E="<<result_min_E_arg<<endl<<endl;
    }
    else  
    {
    cout<<"calculation is completed with error"<<endl;
      if ((is_error_in_command_line_arguments  == false) && (to_do_extract_initial_data_from_file == true)) {cout<<"May be error occured because of reading wrong parameters from file!"<<endl;}
    }
  cout<<"............................................................"<<endl;  
  }





           
  if (bool_one_up_row_gs_state_arg) {delete[] bool_one_up_row_gs_state_arg;  bool_one_up_row_gs_state_arg=0;}
  if (bool_edge_spin_exist_d_u_l_r_ar) {delete[] bool_edge_spin_exist_d_u_l_r_ar;  bool_edge_spin_exist_d_u_l_r_ar=0;}    
  if (J_ver) {delete[] J_ver;  J_ver=0;}
  if (J_hor) {delete[] J_hor;  J_hor=0;}    


end = clock();


  if (mpirank == __MACR_GS_THREAD_NUM_TO_COMMENT)
  {
  printf("work took:  ");   print_time((int) ((end-start)*1E3/CLOCKS_PER_SEC));   cout<<endl;

  cout<<endl<<"=============================================================== end =============================================================== "<<endl<<endl;
  }




MPI_Finalize();           //  ************ End MPI ***********



return 0;

}







    //--DEBUG--//  cout<<"DEBUG:  "<<"spin1_crd_betw_flip_bond_y="<<spin1_crd_betw_flip_bond_y<<"  spin1_crd_betw_flip_bond_x="<<spin1_crd_betw_flip_bond_x<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<"spin2_crd_betw_flip_bond_y="<<spin2_crd_betw_flip_bond_y<<"  spin2_crd_betw_flip_bond_x="<<spin2_crd_betw_flip_bond_x<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<"flip_bond_2d_index="<<flip_bond_2d_index<<endl;
    //--DEBUG--//  cout<<"DEBUG:  "<<endl;
    //--DEBUG--//  cout<<"DEBUG:  "<<"frame_wished_dy="<<frame_wished_dy<<"  frame_wished_dx="<<frame_wished_dx<<"    ";      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<"frame_shift="<<frame_shift<<"    "<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<endl;
    //--DEBUG--//  cout<<"DEBUG:  "<<"frame_wished_y_dn="<<frame_wished_y_dn<<"  frame_wished_y_up="<<frame_wished_y_up<<"  frame_wished_x_lt="<<frame_wished_x_lt<<"  frame_wished_x_rt="<<frame_wished_x_rt<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<"inner_frame_area_y_dn="<<inner_frame_area_y_dn<<"  inner_frame_area_y_up="<<inner_frame_area_y_up<<"  inner_frame_area_x_lt="<<inner_frame_area_x_lt<<"  inner_frame_area_x_rt="<< inner_frame_area_x_rt<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<endl;
    //--DEBUG--//  cout<<"DEBUG:  "<<"frame_wished_Ly="<<frame_wished_Ly<<"  frame_wished_Lx="<<frame_wished_Lx<<"  frames_border_coincide_yes="<<frames_border_coincide_yes<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<endl;
    //--DEBUG--//  cout<<"DEBUG:  "<<"down_row__x_start="<<down_row__x_start<<"  down_row__x_end="<<down_row__x_end<<"  down_row__x_len="<<down_row__x_len<<"  down_row__y="<<down_row__y<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<"up_row__x_start="<<up_row__x_start<<"  up_row__x_end="<<up_row__x_end<<"  up_row__x_len="<<up_row__x_len<<"  up_row__y="<<up_row__y<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<endl;
    //--DEBUG--//  cout<<"DEBUG:  "<<"left_col__y_start="<<left_col__y_start<<"  left_col__y_end="<<left_col__y_end<<"  left_col__y_len="<<left_col__y_len<<"  left_col__x="<<left_col__x<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<"right_col__y_start="<<right_col__y_start<<"  right_col__y_end="<<right_col__y_end<<"  right_col__y_len="<<right_col__y_len<<"  right_col__x="<<right_col__x<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<endl;
    //--DEBUG--//  cout<<"DEBUG:  "<<"flip_bond_index_y="<<flip_bond_index_y<<"  flip_bond_index_x="<<flip_bond_index_x<<endl;      //--DEBUG--//
    //--DEBUG--//  cout<<"DEBUG:  "<<endl;

  



//  ==============================================================================
//
//  
//                           .-----------------TTTT_-----_______
//                         /''''''''''(______O] ----------____  \______/]_
//      __...---'"""\_ --''   Q                               ___________@
//  |'''                   ._   _______________=---------"""""""
//  |                ..--''|   l L |_l   |
//  |          ..--''      .  /-___j '   '
//  |    ..--''           /  ,       '   '
//  |--''                /           `    \
//                       L__'         \    -
//                                     -    '-.
//                                      '.    /
//                                        '-./
//
//  ==============================================================================




//  display aerror message
void display_on_screen_error_msg(char * chstr_error_text)
{

cout<<endl<<endl;
cout<<"==== ERROR ====  ==== ERROR ====  ==== ERROR ====  V ==== ERROR ==== V  ==== ERROR ====  ==== ERROR ====  ==== ERROR ===="<<endl;
cout<<"-------------------------  ";
cout<<chstr_error_text;
cout<<"  -------------------------"<<endl;
cout<<"==== ERROR ====  ==== ERROR ====  ==== ERROR ====  ^ ==== ERROR ==== ^  ==== ERROR ====  ==== ERROR ====  ==== ERROR ===="<<endl;
cout<<endl;

} 


// ==============================================================================
// ==============================================================================


inline void print_progress(unsigned long long int progress, double & mustachieve_progress_perc, unsigned long long int max_progress, ofstream &out)
{

progress=((double) progress)/((double) max_progress)*100.0;

  if (progress >= mustachieve_progress_perc) 
  {  
  cout<<"progress = "<<progress<<" % "<<endl;  out<<"progress = "<<progress<<" % "<<endl;  mustachieve_progress_perc+=1.0;  
  }

}



// ==============================================================================
// ==============================================================================



inline void print_progress(unsigned long long int progress, unsigned long long int max_progress, ofstream &out)
{

progress=((double) progress)/((double) max_progress)*100.0;  
cout<<"progress = "<<progress<<" % "<<endl;  out<<"progress = "<<progress<<" % "<<endl;  

}



// ==============================================================================
// ==============================================================================



//       write char array to file
int write_chstr_to_file(char * text, const char * const filename, bool B1New0Add)
{

  if (filename == 0) {return 0;}
  if (text == 0) {return 0;}

int filenamelen=strlen(filename);
int textlen=strlen(text);

  if (filenamelen == 0) {return 0;}
  if (textlen == 0) {return 0;}

FILE *file=0;

  if (B1New0Add == true) {file=fopen(filename, "wt");} else {file=fopen(filename, "at");}
  if (file == 0) {return 0;}

int filesize=fwrite(text, textlen, 1, file);
fclose(file);

return filesize;

}









inline int amount_of_bit1(unsigned long long int num)
{

int posit_pair_count=0;

  while (num)
  {
  num &= num-1;
  ++posit_pair_count;
  }

return posit_pair_count;

}



inline int amount_of_difdir_pairs_pbc(unsigned long long int row_state, int len)
{

//  we will compare:   row_state  and  row_state_copy (cyclic shifted)
//                                                                  //   0000 0011 0111  len=6
unsigned long long int row_state_copy=row_state << 1;               //   0000 0110 1110
row_state_copy &=~(((unsigned long long int) 1) << len);            //   0000 0100 0000;   1111 1011 1111;    0000 0010 1110
row_state_copy+=(row_state >> (len-1));                             //   0000 0010 1110 +  0000 0000 0001  =  0000 0010 1111
//                                                                  //   was:     0000 0011 0111
//                                                                  //   becaome: 0000 0010 1111 cycle shift completed
row_state_copy=(row_state_copy ^ row_state);                        //   now row_state_copy contains 1 bits     where bit1 means difdir pair 

return amount_of_bit1(row_state_copy);

}





inline int amount_of_difdir_pairs_fbc(unsigned long long int row_state, int len)
{

//  we will compare:   row_state  and  row_state_copy (cyclic shifted)
//                                                                  //   0000 0011 0111  len=6
unsigned long long int row_state_copy=row_state;                    //   0000 0110 1110
row_state_copy &=~(((unsigned long long int) 1) << (len-1));        //   0000 0100 0000;   1111 1011 1111;    0000 0010 1110
row_state=(row_state >> 1);                                         //   0000 0010 1110 +  0000 0000 0001  =  0000 0010 1111
//                                                                  //   was:     0000 0011 0111
//                                                                  //   becaome: 0000 0010 1111 cycle shift completed
row_state_copy=(row_state_copy ^ row_state);                        //   now row_state_copy contains 1 bits     where bit1 means difdir pair 

return amount_of_bit1(row_state_copy);
return 0;

}




//    calculate amount of pairs with dif means of bit in opposite located spins in 2 rows.   pair is only 2 bit-neighbours within of len. 
inline int amount_of_difdir_pairs_betw_2_rows(unsigned long long int row1_state1, unsigned long long int row2_state2, int len)
{

//  we will compare:   row_state  and  row_state_copy (cyclic shifted)
//                                                                  //   0000 0011 0111  len=6
unsigned long long int row1_state1_copy=row1_state1 << 1;           //   0000 0110 1110
row1_state1_copy &=~(((unsigned long long int) 1) << len);          //   0000 0100 0000;   1111 1011 1111;    0000 0010 1110
row1_state1_copy+=(row1_state1 >> (len-1));                         //   0000 0010 1110 +  0000 0000 0001  =  0000 0010 1111
//                                                                  //   was:     0000 0011 0111
//                                                                  //   becaome: 0000 0010 1111 cycle shift completed
row1_state1_copy=(row1_state1 ^ row2_state2);                         //   now row_state_copy contains 1 bits     where bit1 means difdir pair 

return amount_of_bit1(row1_state1_copy);

}




bool generate_2D_lattice_J_v1(bool * J_ver,  bool * J_hor,  int Ly,  int Lx, int plus_percent, int negat_percent, bool pbc_1_fbs_0)
{

  if ((Ly < 2) || (Lx < 2)) {return false;}
  if (J_ver == 0) {return false;}
  if (J_hor == 0) {return false;}
  if (plus_percent+negat_percent != 100) {return false;}

  if ((plus_percent == 0) && (negat_percent == 0)) {return false;}



int bond_amount=0;
int ver_bond_amount=0;
int hor_bond_amount=0;

  if (pbc_1_fbs_0 == true) 
  {
    if ((Ly != 2) && (Lx != 2)) 
    {ver_bond_amount=Ly*Lx;  hor_bond_amount=Ly*Lx;} 
    else 
    {  if ((Ly == 2) && (Lx == 2)) {ver_bond_amount=2;  hor_bond_amount=2;} else {  if (Ly == 2) {ver_bond_amount=Lx; hor_bond_amount=2*Lx;}  if (Lx == 2) {hor_bond_amount=Ly; ver_bond_amount=2*Ly;}  }   }
  bond_amount=ver_bond_amount+hor_bond_amount;  
  }
  else {ver_bond_amount=(Ly-1)*Lx;  hor_bond_amount=(Lx-1)*Ly;  bond_amount=ver_bond_amount+hor_bond_amount;}


double dbl_plus_bond_count_ratio_must_be=((double) plus_percent)/((double) (plus_percent+negat_percent));
double dbl_plus_bond_count_must_be=(((double) bond_amount)*dbl_plus_bond_count_ratio_must_be);
int    int_plus_bond_count_must_be=lround(dbl_plus_bond_count_must_be);

//  if ((dbl_plus_bond_count_must_be-((double) int_plus_bond_count_must_be)) > 0.5) {int_plus_bond_count_must_be++;}     //   <----- what is it ???!!!
  if (int_plus_bond_count_must_be > bond_amount) {int_plus_bond_count_must_be=bond_amount;}  



bool * J_temp=0;
J_temp=new bool[bond_amount];

  for (int i=0;  i < bond_amount;  i++)
  {J_temp[i]=false;}



int number_bonds_within__out_ver_var=(Ly+1)*Lx,  number_bonds_within__out_hor_var=Ly*(Lx+1);      // they can be repeated,  it's number the same for fbc and pbc

  for (int i=0;  i < number_bonds_within__out_ver_var;  i++)
  {J_ver[i]=false;}

  for (int i=0;  i < number_bonds_within__out_hor_var;  i++)
  {J_hor[i]=false;}



bool exit_of_random_cycle=false;

  if ((negat_percent == 0) || (plus_percent == 0))
  {
    
    if (negat_percent == 0)
    {
      for (int i=0;  i < bond_amount;  i++)
      {J_temp[i]=true;}
    }
    
    if (plus_percent == 0)
    {
      for (int i=0;  i < bond_amount;  i++)
      {J_temp[i]=false;}
    }

  exit_of_random_cycle=true;
  }



// --------- main cycle
//
  {
    if (exit_of_random_cycle == false)
    {
      for (int i=0;  i < int_plus_bond_count_must_be;  i++)
      {
      int rand_index=(rand() % (bond_amount-i));
      int counter=0;

        for (int j=0;  j < bond_amount;  j++)
        {   
          if (J_temp[j] == false) {  if (counter == rand_index) {J_temp[j]=true;  break;} else {counter++;}  }
        }

      }
    }
  }

//
// --------- main cycle





//  carring random values into arrays
//  
  if (pbc_1_fbs_0 == true)   //  pbc
  {

    if (Ly == 2)
    {
      for (int i=0;  i < Lx;  i++)      //  filling main rect area
      {J_ver[i]=J_temp[i];  J_ver[i+Lx]=J_temp[i];  J_ver[i+Lx+Lx]=J_temp[i];}
    }
    else
    {
      for (int i=0;  i < Ly*Lx;  i++)      //  filling main rect area
      {J_ver[i+Lx]=J_temp[i];}

      for (int i=0;  i < Lx;  i++)         //  filling down last row (copy of up last)
      {J_ver[i]=J_ver[Ly*Lx+i];}
    }

  //  *** hor
  //
    if (Lx == 2)
    {
      for (int i=0;  i < Ly;  i++)
      {
      J_hor[i*(Lx+1)+0]=J_temp[ver_bond_amount+i];  J_hor[i*(Lx+1)+1]=J_temp[ver_bond_amount+i];  J_hor[i*(Lx+1)+2]=J_temp[ver_bond_amount+i];
      }
    }
    else
    {
      for (int i=0;  i < Ly;  i++)
      {
        for (int j=0;  j < Lx;  j++) 
        {
        J_hor[i*(Lx+1)+j+1]=J_temp[ver_bond_amount+i*Lx+j];      //  filling main rect area
        //==DEBUG==//   cout<<" A "<<bond_amount<<"  "<<(i*(Lx+1)+j+1)<<"   "<<(ver_bond_amount+i*Lx+j)<<endl;  //==DEBUG==//
        }
      }

      for (int i=0;  i < Ly;  i++)
      {
      J_hor[i*(Lx+1)]=J_hor[i*(Lx+1)+Lx];          //  filling left column (copy of right)
      }
    }

  }
  else  //  fbc
  {

  //  *** ver
    for (int i=0;  i < (Ly-1)*Lx;  i++)      //  filling main rect area
    {J_ver[i+Lx]=J_temp[i];}

  //  *** hor
    for (int i=0;  i < Ly;  i++)
    {
      for (int j=0;  j < Lx-1;  j++) 
      {
      J_hor[i*(Lx+1)+j+1]=J_temp[ver_bond_amount+i*(Lx-1)+j];      //  filling main rect area
      }
    }

  }


  if (J_temp) {delete[] J_temp;  J_temp=0;}


return true;

}







bool generate_2D_lattice_J_v2(bool * J_ver,  bool * J_hor,  int Ly,  int Lx, int plus_percent, int negat_percent, bool pbc_1_fbs_0)
{


  if (J_ver == 0) {return false;}
  if (J_hor == 0) {return false;}

  if ((plus_percent == 0) && (negat_percent == 0)) {return false;}


bool * J_temp=0;

int bond_amount=Ly*Lx*2;
int ver_bond_amount=Ly*Lx;
int hor_bond_amount=Ly*Lx;

  if (pbc_1_fbs_0 == false) {bond_amount=Ly*Lx*2-Lx-Ly;  ver_bond_amount=(Ly-1)*Lx;  hor_bond_amount=(Lx-1)*Ly; }



J_temp=new bool[bond_amount];

  for (int i=0;  i < bond_amount;  i++)
  {J_temp[i]=false;}


double dbl_plus_bond_amount_ratio_must_be=((double) plus_percent)/((double) (plus_percent+negat_percent));
double dbl_plus_bond_amount_must_be=((double) (plus_percent*bond_amount))/((double) (plus_percent+negat_percent));
int    int_plus_bond_amount_must_be=(int)     (((double) (plus_percent*bond_amount))/((double) (plus_percent+negat_percent)));

  if ((dbl_plus_bond_amount_must_be-((double) int_plus_bond_amount_must_be)) > 0.5) {int_plus_bond_amount_must_be++;}
  if (int_plus_bond_amount_must_be > bond_amount) {int_plus_bond_amount_must_be=bond_amount;}  

  if (pbc_1_fbs_0 == true)
  {
    for (int i=0;  i < bond_amount/2+Lx; i++)
    {J_ver[i]=false;}
    for (int i=0;  i < bond_amount/2+Ly; i++)
    {J_hor[i]=false;}

  }
  else
  {
    for (int i=0;  i < bond_amount/2; i++)
    {J_ver[i]=false;  J_hor[i]=false;}
  }


int true_counter=0;
bool exit=false;    


  if ((negat_percent == 0) || (plus_percent == 0))
  {
    
    if (negat_percent == 0)
    {
      for (int i=0;  i < bond_amount;  i++)
      {J_temp[i]= true;}
    }
    
    if (plus_percent == 0)
    {
      for (int i=0;  i < bond_amount;  i++)
      {J_temp[i]=false;}
    }

  exit=true;
  }


  if (true_counter >= int_plus_bond_amount_must_be) {exit=true;}  

  while (exit == false)
  {
  int rand_index=(rand() % (bond_amount));      //  choose random index of bond
  
    if (J_temp[rand_index] == false) {J_temp[rand_index]=true;  true_counter++;}  
    if (true_counter >= int_plus_bond_amount_must_be) {exit=true;}
  
  //cout<<"ggggggggggg   "<<int_plus_bond_amount_must_be<<"      hhhhhhh"<<endl;
  }
  

  if (pbc_1_fbs_0 == true)
  {

    for (int i=0;  i < Ly*Lx;  i++)
    {J_ver[i+Lx]=J_temp[i];}

    for (int i=0;  i < Lx;  i++)
    {J_ver[i]=J_ver[Ly*Lx+i];}
  

    for (int i=0;  i < Ly;  i++)
    {
      for (int j=0;  j < Lx;  j++)
      {
      J_hor[i*Lx+j+i+1]=J_temp[Ly*Lx+i*Lx+j];
      }
    }

    for (int i=0;  i < Ly;  i++)
    {
    J_hor[i*(Lx+1)]=J_hor[(i+1)*(Lx+1)-1];
    }

  }
  else
  {

    for (int i=0;  i < ver_bond_amount;  i++)
    {J_ver[i]=J_temp[i];}

    for (int i=0;  i < hor_bond_amount;  i++)
    {J_hor[i]=J_temp[i+ver_bond_amount];}  

  }


  if (J_temp) {delete[] J_temp;  J_temp=0;}


return true;

}






bool generate_2D_lattice_J_v3(bool * J_ver,  bool * J_hor,  int Ly,  int Lx, int plus_percent, int negat_percent)
{

  if ((Ly < 2) || (Lx < 2)) {return false;}
  if (J_ver == 0) {return false;}
  if (J_hor == 0) {return false;}
  if (plus_percent+negat_percent != 100) {return false;}

  if ((plus_percent == 0) && (negat_percent == 0)) {return false;}



int bond_amount=Ly*(Lx+1)+Lx*(Ly+1);
int ver_bond_amount=Lx*(Ly+1);
int hor_bond_amount=Ly*(Lx+1);


double dbl_plus_bond_count_ratio_must_be=((double) plus_percent)/((double) (plus_percent+negat_percent));
double dbl_plus_bond_count_must_be=(((double) bond_amount)*dbl_plus_bond_count_ratio_must_be);
int    int_plus_bond_count_must_be=lround(dbl_plus_bond_count_must_be);

//  if ((dbl_plus_bond_count_must_be-((double) int_plus_bond_count_must_be)) > 0.5) {int_plus_bond_count_must_be++;}     //   <----- what is it ???!!!
  if (int_plus_bond_count_must_be > bond_amount) {int_plus_bond_count_must_be=bond_amount;}  



bool * J_temp=0;
J_temp=new bool[bond_amount];

  for (int i=0;  i < bond_amount;  i++)
  {J_temp[i]=false;}



int number_bonds_within__out_ver_var=(Ly+1)*Lx,  number_bonds_within__out_hor_var=Ly*(Lx+1);      // they can be repeated,  it's number the same for fbc and pbc

  for (int i=0;  i < number_bonds_within__out_ver_var;  i++)
  {J_ver[i]=false;}

  for (int i=0;  i < number_bonds_within__out_hor_var;  i++)
  {J_hor[i]=false;}



bool exit_of_random_cycle=false;

  if ((negat_percent == 0) || (plus_percent == 0))
  {
    
    if (negat_percent == 0)
    {
      for (int i=0;  i < bond_amount;  i++)
      {J_temp[i]=true;}
    }
    
    if (plus_percent == 0)
    {
      for (int i=0;  i < bond_amount;  i++)
      {J_temp[i]=false;}
    }

  exit_of_random_cycle=true;
  }



// --------- main cycle
//
  {
    if (exit_of_random_cycle == false)
    {
      for (int i=0;  i < int_plus_bond_count_must_be;  i++)
      {
      int rand_index=(rand() % (bond_amount-i));
      int counter=0;

        for (int j=0;  j < bond_amount;  j++)
        {   
          if (J_temp[j] == false) {  if (counter == rand_index) {J_temp[j]=true;  break;} else {counter++;}  }
        }

      }
    }
  }

//
// --------- main cycle





unsigned int counter=0;

  for (int i=0;  i < Ly+1;  i++)
  {
    for (int j=0;  j < Lx;  j++) 
    {
    J_ver[i*Lx+j]=J_temp[counter];      //  filling main rect area
    counter++;
    }
  }


  for (int i=0;  i < Ly;  i++)
  {
    for (int j=0;  j < Lx+1;  j++) 
    {
    J_hor[i*(Lx+1)+j]=J_temp[counter];      //  filling main rect area
    counter++;
    }
  }




  if (J_temp) {delete[] J_temp;  J_temp=0;}


return true;

}






bool invert_J(bool * J_ver,  bool * J_hor,  int Ly,  int Lx, bool pbc_1_fbs_0)
{

  if ((J_ver == 0) || (J_hor == 0)) {return false;}
  if ((Ly < 1) ||  (Lx < 1)) {return false;}

  if (pbc_1_fbs_0 == true)
  {
  
    for (int i=0; i < Lx*(Ly+1); i++)
    {
      if (J_ver[i] == false) {J_ver[i]=true;} else {J_ver[i]=false;}
    }
  
    for (int i=0; i < Ly*(Lx+1); i++)
    {
      if (J_hor[i] == false) {J_hor[i]=true;} else {J_hor[i]=false;}
    }

  }
  else
  {

    for (int i=0; i < Lx*(Ly+1); i++)
    {
      if (J_ver[i] == false) {J_ver[i]=true;} else {J_ver[i]=false;}
    }
  
    for (int i=0; i < Ly*(Lx+1); i++)
    {
      if (J_hor[i] == false) {J_hor[i]=true;} else {J_hor[i]=false;}
    }

  }

return true;

}




inline void inver_J(bool *J_ver,  bool *J_hor,  int Ly,  int Lx)
{

  for (int index_y=0;  index_y < Ly+1;  index_y++)
  {
    for (int index_x=0;  index_x < Lx;  index_x++)
    {
    J_ver[index_y*Lx+index_x]=!J_ver[index_y*Lx+index_x];
    }
  }

  for (int index_y=0;  index_y < Ly;  index_y++)
  {
    for (int index_x=0;  index_x < Lx+1;  index_x++)
    {
    J_hor[index_y*(Lx+1)+index_x]=!J_hor[index_y*(Lx+1)+index_x];
    }
  }

}





inline void read_ullint_J_ver_line(unsigned long long int & J_ver_line,  bool * J_ver_array,  int ver_line_num, int amount_of_bond_read, int Lx)
{

unsigned long long int moving_1_bit1=1;
J_ver_line=0;

  for (int index_y=0;  index_y < amount_of_bond_read;  index_y++)
  {
    if (J_ver_array[(Lx+1)*index_y+ver_line_num] == true) {J_ver_line+=moving_1_bit1;}

  moving_1_bit1=moving_1_bit1 << 1;
  }

}




inline void read_ullint_J_hor_line(unsigned long long int & J_hor_line,  bool * J_hor_array,  int hor_line_num, int amount_of_bond_read, int Lx)
{

unsigned long long int moving_1_bit1=1;
J_hor_line=0;

  for (int index_x=0;  index_x < amount_of_bond_read;  index_x++)
  {
    if (J_hor_array[(Lx+1)*hor_line_num+index_x] == true) {J_hor_line+=moving_1_bit1;}

  moving_1_bit1=moving_1_bit1 << 1;
  }

}



//  xor
//  0 ^ 0 ^ 0 = 0      1 ^ 1 ^ 1 = 1
//  0 ^ 0 ^ 1 = 1      0 ^ 1 ^ 0 = 1     1 ^ 0 ^ 0 = 1
//  0 ^ 1 ^ 1 = 0      1 ^ 0 ^ 1 = 0     1 ^ 1 ^ 0 = 0
//

inline int amount_of_bit1_logistic_xor_3_nums(unsigned long long int num1,  unsigned long long int num2,  unsigned long long int num3, int len)
{

unsigned long long int mult_result=(num1 ^ num2 ^ num3);      //  !  !  !

int bit_1_count=0;
                  
  while (mult_result)
  {
  mult_result &= mult_result-1;
  ++bit_1_count;
  }

return bit_1_count;

}




inline int amount_of_bit1_logistic_xor_1_num_with_shift_pbc(unsigned long long int J_line,  unsigned long long int spin_line_state, int len)
{

unsigned long long int spin_line_state_copy=spin_line_state << 1;               //   0000 0110 1110
spin_line_state_copy &=~(((unsigned long long int) 1) << len);            //   0000 0100 0000;   1111 1011 1111;    0000 0010 1110
spin_line_state_copy+=(spin_line_state >> (len-1));        

unsigned long long int mult_result=(spin_line_state ^ spin_line_state_copy ^ J_line);      //  !  !  !

int bit_1_count=0;
                  
  while (mult_result)
  {
  mult_result &= mult_result-1;
  ++bit_1_count;
  }

return bit_1_count;

}




inline int amount_of_bit1_logistic_xor_1_num_with_shift_fbc(unsigned long long int J_line,  unsigned long long int spin_line_state, int len)
{

unsigned long long int spin_line_state_copy=spin_line_state;                    //   0000 0110 1110
spin_line_state_copy &=~(((unsigned long long int) 1) << (len-1));        //   0000 0100 0000;   1111 1011 1111;    0000 0010 1110
spin_line_state=(spin_line_state >> 1);             

unsigned long long int mult_result=(spin_line_state ^ spin_line_state_copy ^ J_line);      //  !  !  !

int bit_1_count=0;
                  
  while (mult_result)
  {
  mult_result &= mult_result-1;
  ++bit_1_count;
  }

return bit_1_count;

}




void calc_pospair_and_upspin_EA_2D_lattice_PBC(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount,  int * up_amount)
{

int energy_loc=0,  magnetization_loc=0;

  for (int index_y=0;  index_y < Ly; index_y++)
  {

    for (int index_x=0;  index_x < Lx; index_x++)
    {
    int main_index=index_y*Lx+index_x;
    int next_ver_index=main_index+Lx;    // with respecting of PBC
      if (index_y >= Ly-1) {next_ver_index=index_x;}
    int next_hor_index=main_index+1;   // with respecting of PBC
      if (index_x >= Lx-1) {next_hor_index=index_y*Lx;}

    // energy: ver bond
    bool pair_energy=false;
      if (spin_array[main_index] != spin_array[next_ver_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_ver[index_y*Lx+index_x+Lx]);
      if (pair_energy == true) {energy_loc++;}      //  ! ! !
    

    // energy: hor bond
    pair_energy=false;
      if (spin_array[main_index] != spin_array[next_hor_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_hor[index_y*(Lx+1)+index_x+1]);
        if (pair_energy == true) {energy_loc++;}      //  ! ! !

    // magnetization:
      if (spin_array[main_index] == true)  {magnetization_loc++;}
    
    }

  }

  if (pospair_amount ) {*pospair_amount=energy_loc;       }
  if (up_amount      ) {*up_amount     =magnetization_loc;}

}




void calc_pospair_and_upspin_EA_2D_lattice_FBC(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount,  int * up_amount)
{

int energy_loc=0,  magnetization_loc=0;

  for (int index_y=0;  index_y < Ly-1; index_y++)
  {

    for (int index_x=0;  index_x < Lx-1; index_x++)
    {
    int main_index=index_y*Lx+index_x;
    int next_ver_index=main_index+Lx;    // with respecting of PBC
      if (index_y >= Ly-1) {next_ver_index=index_x;}
    int next_hor_index=main_index+1;   // with respecting of PBC
      if (index_x >= Lx-1) {next_hor_index=index_y*Lx;}

    // energy: ver bond
    bool pair_energy=false;
      if (spin_array[main_index] != spin_array[next_ver_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_ver[index_y*(Lx-1)+index_x]);
      if (pair_energy == true) {energy_loc++;}      //  ! ! !

    // energy: hor bond
    pair_energy=false;
      if (spin_array[main_index] != spin_array[next_hor_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_hor[index_y*(Lx-1)+index_x]);
        if (pair_energy == true) {energy_loc++;}      //  ! ! !

    // magnetization:
      if (spin_array[main_index] == true)  {magnetization_loc++;}
    
    }

  }


  if (pospair_amount ) {*pospair_amount=energy_loc;       }
  if (up_amount      ) {*up_amount     =magnetization_loc;}

}







void calc_pospair_and_upspin_EA_2D_lattice_PBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount,  int * up_amount)
{

int energy_loc=0,  magnetization_loc=0;

  for (int index_y=0;  index_y < Ly; index_y++)
  {

    for (int index_x=0;  index_x < Lx; index_x++)
    {
    int main_index=index_y*Lx+index_x;
    int next_ver_index=main_index+Lx;    // with respecting of PBC
      if (index_y >= Ly-1) {next_ver_index=index_x;}
    int next_hor_index=main_index+1;   // with respecting of PBC
      if (index_x >= Lx-1) {next_hor_index=index_y*Lx;}

    // energy: ver bond
    bool pair_energy=false;
      if (spin_array[main_index] != spin_array[next_ver_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_ver[index_y*Lx+index_x+Lx]);
      if (pair_energy == true) {energy_loc++;}      //  ! ! !
    

    // energy: hor bond
    pair_energy=false;
      if (spin_array[main_index] != spin_array[next_hor_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_hor[index_y*(Lx+1)+index_x+1]);
        if (pair_energy == true) {energy_loc++;}      //  ! ! !

    // magnetization:
      if (spin_array[main_index] == true)  {magnetization_loc++;}
    
    }

  }

  if (pospair_amount ) {*pospair_amount=energy_loc;       }
  if (up_amount      ) {*up_amount     =magnetization_loc;}

}







void calc_pospair_and_upspin_EA_2D_lattice_FBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount,  int * up_amount)
{
//  neg pair if e.g.:    s1=d  s2=d  J=-1   (fer)           ---=-                      pos pair if e.g.:    s1=d  s2=u  J=-1   (fer)           -+-=+

int energy_loc=0,  magnetization_loc=0;

  for (int index_y=0;  index_y < Ly; index_y++)
  {

    for (int index_x=0;  index_x < Lx; index_x++)
    {
    int main_index=index_y*Lx+index_x;
    int next_ver_index=main_index+Lx;    // with respecting of PBC
      if (index_y >= Ly-1) {next_ver_index=index_x;}
    int next_hor_index=main_index+1;   // with respecting of PBC
      if (index_x >= Lx-1) {next_hor_index=index_y*Lx;}

    // energy: ver bond
    bool pair_energy=false;
      if (index_y < Ly-1)
      {
        if (spin_array[main_index] != spin_array[next_ver_index]) {pair_energy=true;}
      pair_energy=(pair_energy != J_ver[index_y*Lx+index_x+Lx]);
        if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }

    // energy: hor bond
    pair_energy=false;
      if (index_x < Lx-1)
      {
        if (spin_array[main_index] != spin_array[next_hor_index]) {pair_energy=true;}
      pair_energy=(pair_energy != J_hor[index_y*(Lx+1)+index_x+1]);
          if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }

    // magnetization:
      if (spin_array[main_index] == true)  {magnetization_loc++;}
    
    }

  }


  if (pospair_amount ) {*pospair_amount=energy_loc;       }
  if (up_amount      ) {*up_amount     =magnetization_loc;}

}




void calc_pospair_and_upspin_EA_2D_lattice_with_edges_without_down_edge_v3(bool * spin_array,  bool * bool_edge_spin_exist_d_u_l_r_ar,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount,  int * up_amount,  bool calc_edge_down)
{
//  neg pair if e.g.:    s1=d  s2=d  J=-1   (fer)           ---=-                      pos pair if e.g.:    s1=d  s2=u  J=-1   (fer)           -+-=+

int energy_loc=0,  magnetization_loc=0;

  for (int index_y=0;  index_y < Ly; index_y++)
  {

    for (int index_x=0;  index_x < Lx; index_x++)
    {
    int main_index=index_y*Lx+index_x;
    int next_ver_index=main_index+Lx;    // with respecting of PBC
      if (index_y >= Ly-1) {next_ver_index=index_x;}
    int next_hor_index=main_index+1;   // with respecting of PBC
      if (index_x >= Lx-1) {next_hor_index=index_y*Lx;}

    // energy: ver bond
    bool pair_energy=false;
      if (index_y < Ly-1)
      {
        if (spin_array[main_index] != spin_array[next_ver_index]) {pair_energy=true;}
      pair_energy=(pair_energy != J_ver[index_y*Lx+index_x+Lx]);
        if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }

    // energy: hor bond
    pair_energy=false;
      if (index_x < Lx-1)
      {
        if (spin_array[main_index] != spin_array[next_hor_index]) {pair_energy=true;}
      pair_energy=(pair_energy != J_hor[index_y*(Lx+1)+index_x+1]);
          if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }

    // magnetization:
      if (spin_array[main_index] == true)  {magnetization_loc++;}
    
    }

  }



//  edges
  {
  //  ===
  //
  bool pair_energy=false;
  
    for (int index_x=0;  index_x < Lx; index_x++)
    {
    
    //  down edge
      if (calc_edge_down == true)
      {
      pair_energy=false;
        if (bool_edge_spin_exist_d_u_l_r_ar[2*index_x+1] == true) 
        {
          if (bool_edge_spin_exist_d_u_l_r_ar[2*index_x] != spin_array[index_x]) {pair_energy=true;}
        pair_energy=(pair_energy != J_ver[index_x]);   
          if (pair_energy == true) {energy_loc++;}      //  ! ! !
        }
      }
      
    //  up edge
    pair_energy=false;
      if (bool_edge_spin_exist_d_u_l_r_ar[Lx*2+2*index_x+1] == true) 
      {
        if (bool_edge_spin_exist_d_u_l_r_ar[Lx*2+2*index_x] != spin_array[(Ly-1)*Lx+index_x]) {pair_energy=true;}
      pair_energy=(pair_energy != J_ver[Ly*Lx+index_x]);   
        if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }
      
    }      //  for (int index_x=0;  index_x < Lx; index_x++)


  //  ===
  // 
  pair_energy=false;
    
    for (int index_y=0;  index_y < Ly; index_y++)
    {
    
    //  left edge
    pair_energy=false;
      if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*index_y+1] == true) 
      {
        if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*index_y] != spin_array[index_y*Lx]) {pair_energy=true;}
      pair_energy=(pair_energy != J_hor[index_y*(Lx+1)]);   
        if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }
      
    //  up edge
    pair_energy=false;
      if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*Ly+2*index_y+1] == true) 
      {
        if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*Ly+2*index_y] != spin_array[index_y*Lx+Lx-1]) {pair_energy=true;}
      pair_energy=(pair_energy != J_hor[index_y*(Lx+1)+Lx]);   
        if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }    
    }      //cout<<"DEBUG:  "<<" energy_loc="<<energy_loc<<endl;  //  for (int index_x=0;  index_x < Lx; index_x++)

  }


  if (pospair_amount ) {*pospair_amount=energy_loc;       }
  if (up_amount      ) {*up_amount     =magnetization_loc;}

}










void calc_pospair_and_upspin_EA_2D_lattice_cylind_PBC_FBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount,  int * up_amount)
{
//  neg pair if e.g.:    s1=d  s2=d  J=-1   (fer)           ---=-                      pos pair if e.g.:    s1=d  s2=u  J=-1   (fer)           -+-=+

int energy_loc=0,  magnetization_loc=0;

  for (int index_y=0;  index_y < Ly; index_y++)
  {

    for (int index_x=0;  index_x < Lx; index_x++)
    {
    int main_index=index_y*Lx+index_x;
    int next_ver_index=main_index+Lx;    // with respecting of PBC
      if (index_y >= Ly-1) {next_ver_index=index_x;}
    int next_hor_index=main_index+1;   // with respecting of PBC
      if (index_x >= Lx-1) {next_hor_index=index_y*Lx;}

    // energy: ver bond
    bool pair_energy=false;
      if (spin_array[main_index] != spin_array[next_ver_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_ver[index_y*Lx+index_x+Lx]);
      if (index_y < Ly-1)       //  cylindric condition, we don't consider the top with bottom
      {
        if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }

    // energy: hor bond
    pair_energy=false;
      if (spin_array[main_index] != spin_array[next_hor_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_hor[index_y*(Lx+1)+index_x+1]);
        if (pair_energy == true) {energy_loc++;}      //  ! ! !

    // magnetization:
      if (spin_array[main_index] == true)  {magnetization_loc++;}
    
    }

  }

  if (pospair_amount ) {*pospair_amount=energy_loc;       }
  if (up_amount      ) {*up_amount     =magnetization_loc;}

}





void calc_pospair_and_upspin_EA_2D_lattice_cylind_hor_PBC_ver_FBC_v3(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount,  int * up_amount)
{
//  neg pair if e.g.:    s1=d  s2=d  J=-1   (fer)           ---=-                      pos pair if e.g.:    s1=d  s2=u  J=-1   (fer)           -+-=+

int energy_loc=0,  magnetization_loc=0;

  for (int index_y=0;  index_y < Ly; index_y++)
  {

    for (int index_x=0;  index_x < Lx; index_x++)
    {
    int main_index=index_y*Lx+index_x;
    int next_ver_index=main_index+Lx;    // with respecting of PBC
      if (index_y >= Ly-1) {next_ver_index=index_x;}
    int next_hor_index=main_index+1;   // with respecting of PBC
      if (index_x >= Lx-1) {next_hor_index=index_y*Lx;}

    // energy: ver bond
    bool pair_energy=false;
      if (spin_array[main_index] != spin_array[next_ver_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_ver[index_y*Lx+index_x+Lx]);
      if (index_y < Ly-1)       //  cylindric condition, we don't consider the top with bottom
      {
        if (pair_energy == true) {energy_loc++;}      //  ! ! !
      }

    // energy: hor bond
    pair_energy=false;
      if (spin_array[main_index] != spin_array[next_hor_index]) {pair_energy=true;}
    pair_energy=(pair_energy != J_hor[index_y*(Lx+1)+index_x+1]);
        if (pair_energy == true) {energy_loc++;}      //  ! ! !

    // magnetization:
      if (spin_array[main_index] == true)  {magnetization_loc++;}
    
    }

  }

  if (pospair_amount ) {*pospair_amount=energy_loc;       }
  if (up_amount      ) {*up_amount     =magnetization_loc;}

}







inline int positi_E_of_1_par(bool spin_1,  bool spin_2, bool J)
{

bool pair_energy=false;

  if (spin_1 != spin_2) {pair_energy=true;}
pair_energy=(pair_energy != J);  
  if (pair_energy == true) {return 1;}

return 0;

}





inline int amount_of_positiv_E_par_edge_v3(bool *J_line,  bool * bool_edge_spin_exist_d_u_l_r_ar, unsigned long long int spin_line_state, int len)
{

int energy_loc=0;
unsigned long long int shifting_one=1;
bool spin_of__spin_line_state=false;
bool pair_energy=false;

  for (int i=0;  i < len;  i++)
  {
  pair_energy=false;
  
    if (bool_edge_spin_exist_d_u_l_r_ar[2*i+1] == true)
    {
    spin_of__spin_line_state=(((shifting_one << (i)) & spin_line_state) != 0);
  
      if (spin_of__spin_line_state != bool_edge_spin_exist_d_u_l_r_ar[2*i]) {pair_energy=true;}
    pair_energy=(pair_energy != J_line[i]);  
      if (pair_energy == true) {energy_loc++;}
    }   
    
  }

return energy_loc;

}





inline int amount_of_positiv_E_par_revers_order_edge_v3(bool *J_line,  bool * bool_edge_spin_exist_d_u_l_r_ar, unsigned long long int spin_line_state, int len)
{

int energy_loc=0;
unsigned long long int shifting_one=1;
bool spin_of__spin_line_state=false;
bool pair_energy=false;

  for (int i=0;  i < len;  i++)
  {
  pair_energy=false;
  
    if (bool_edge_spin_exist_d_u_l_r_ar[2*i+1] == true)
    {
    spin_of__spin_line_state=(((shifting_one << (len-i-1)) & spin_line_state) != 0);
  
      if (spin_of__spin_line_state != bool_edge_spin_exist_d_u_l_r_ar[2*i]) {pair_energy=true;}
    pair_energy=(pair_energy != J_line[i]);  
      if (pair_energy == true) {energy_loc++;}
    }   
    
  }

return energy_loc;

}




inline int amount_of_positiv_E_par_bool_edge_v3(bool *J_line,  bool * bool_edge_spin_exist_d_u_l_r_ar,  bool * bool_spin_array, int len)
{

int energy_loc=0;
bool pair_energy=false;

  for (int i=0;  i < len;  i++)
  {
  pair_energy=false;
  
    if (bool_edge_spin_exist_d_u_l_r_ar[2*i+1] == true)
    {  
      if (bool_spin_array[i] != bool_edge_spin_exist_d_u_l_r_ar[2*i]) {pair_energy=true;}
    pair_energy=(pair_energy != J_line[i]);
      if (pair_energy == true) {energy_loc++;}  
    }   
    
  }

return energy_loc;

}



void calc_energy_and_magnetiz_EA_2D_lat_PBC(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * energy,  int * magnetization)
{

int total_amount_of_bond=2*Ly*Lx;
int spin_amount=Ly*Lx;
int pospair_amount=0,  up_amount=0;
calc_pospair_and_upspin_EA_2D_lattice_PBC(spin_array,  J_ver,  J_hor,  Ly,  Lx,  & pospair_amount,  & up_amount);

  if (energy       ) {*energy=pospair_amount*2-total_amount_of_bond;}
  if (magnetization) {*magnetization=up_amount*2-spin_amount;}

}



void calc_energy_and_magnetiz_EA_2D_lat_FBC(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * energy,  int * magnetization)
{

int total_amount_of_bond=2*Ly*Lx-Lx-Ly;
int spin_amount=Ly*Lx;
int pospair_amount=0,  up_amount=0;
calc_pospair_and_upspin_EA_2D_lattice_FBC(spin_array,  J_ver,  J_hor,  Ly,  Lx,  & pospair_amount,  & up_amount);

  if (energy       ) {*energy=pospair_amount*2-total_amount_of_bond;}
  if (magnetization) {*magnetization=up_amount*2-spin_amount;}

}




void calc_energy_and_magnetiz_EA_2D_lat_PBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * energy,  int * magnetization)
{

int total_amount_of_bond=2*Ly*Lx;
int spin_amount=Ly*Lx;
int pospair_amount=0,  up_amount=0;
calc_pospair_and_upspin_EA_2D_lattice_PBC_v2(spin_array,  J_ver,  J_hor,  Ly,  Lx,  & pospair_amount,  & up_amount);

  if (energy       ) {*energy=pospair_amount*2-total_amount_of_bond;}
  if (magnetization) {*magnetization=up_amount*2-spin_amount;}

}



int calc_magnetiz_lat(bool * spin_array,  int Ly,  int Lx)
{

  if (spin_array == 0) {return 0;}

int magnetization=0;

  for (int i=0;  i < Ly;  i++)
  {
    for (int j=0;  j < Lx;  j++)
    {  
      if (spin_array[i*Lx+j] == true) {magnetization++;} else {magnetization--;}
    }
  }

return magnetization;

}




void calc_energy_and_magnetiz_EA_2D_lat_FBC_v2(bool * spin_array,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * energy,  int * magnetization)
{

int total_amount_of_bond=2*Ly*Lx-Lx-Ly;
int spin_amount=Ly*Lx;
int pospair_amount=0,  up_amount=0;
calc_pospair_and_upspin_EA_2D_lattice_FBC_v2(spin_array,  J_ver,  J_hor,  Ly,  Lx,  & pospair_amount,  & up_amount);

  if (energy       ) {*energy=pospair_amount*2-total_amount_of_bond;}
  if (magnetization) {*magnetization=up_amount*2-spin_amount;}

}




inline int pospair_amount_1_cell_EA(bool up1dn0 , bool *J, char *s)
{

//          s[1]
//   s[2]   up1dn  s[0]
//          s[3]
//DEBUG==DEBUG//  cout<<"aa"<<J[0]<<J[1]<<J[2]<<J[3]<<"aa";    //DEBUG==DEBUG//


int energy=0;
//  char energy_1_pair=0;
int var1=0, var2=0, var3=0;

  if (up1dn0 == true) {var1=1;}


  if (s[0] != 0) 
  {
    if (J[0] == true) {var2=1;} else {var2=0;}
    if (s[0] > 0) {var3=1;} else {var3=0;}
    if ((var1 ^ var2 ^ var3) != 0)  {energy++;}
  //DEBUG==DEBUG//  cout<<" 0) "<<var1<<" "<<var2<<" "<<var3<<" en="<<energy<"  ";  //DEBUG==DEBUG//  
  }

  if (s[1] != 0) 
  {
    if (J[1] == true) {var2=1;} else {var2=0;}
    if (s[1] > 0) {var3=1;} else {var3=0;}
    if ((var1 ^ var2 ^ var3) != 0) {energy++;}
  //DEBUG==DEBUG//  cout<<" 1) "<<var1<<" "<<var2<<" "<<var3<<" en="<<energy<"  ";  //DEBUG==DEBUG//  
  }

  if (s[2] != 0) 
  {
    if (J[2] == true) {var2=1;} else {var2=0;}
    if (s[2] > 0) {var3=1;} else {var3=0;}
    if ((var1 ^ var2 ^ var3) != 0) {energy++;}
  //DEBUG==DEBUG//  cout<<" 2) "<<var1<<" "<<var2<<" "<<var3<<" en="<<energy<"  ";  //DEBUG==DEBUG//  
  }

  if (s[3] != 0) 
  {
    if (J[3] == true) {var2=1;} else {var2=0;}
    if (s[3] > 0) {var3=1;} else {var3=0;}
    if ((var1 ^ var2 ^ var3) != 0) {energy++;}
  //DEBUG==DEBUG//  cout<<" 3) "<<var1<<" "<<var2<<" "<<var3<<" en="<<energy<"  ";  //DEBUG==DEBUG//  
  }

return ((int) energy);

}





// ============= (( (( < gmp

//  stat variant for chstr_num,      convert mpz_t to chstr var:  necessary len will be set with symb ch_empt_symb       desired_len - len we want
//            if (desired_len >  0) {it will be:  filling with ch_empt_symb of rest pos;} 
//            if (desired_len == 0) {it will be: desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}
//            if (desired_len <  0) {it will be: desired_len=real_len;} 
// 
int mpz_t__to__chstr__stat(mpz_t * mpz_t_num,   char * const chstr_num,  int desired_len,  char ch_empt_symb,  bool fill_empt_symb_left1_right0)
{ 

bool is_real_len=false;

  if ((chstr_num == 0) || (mpz_t_num == 0)) {return 0;}

int not_nul_num_len=mpz_sizeinbase(*mpz_t_num, 10);

  if (desired_len <  0) {desired_len=not_nul_num_len;    is_real_len=true;}
  if (desired_len == 0) {desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}
  if (desired_len < not_nul_num_len) {desired_len=not_nul_num_len;}


char *chstr_temp=0;
chstr_temp=new char[desired_len+1+64];
chstr_temp[0]='\0';
mpz_get_str(chstr_temp, 10,  *mpz_t_num); 
//chstr_temp[desired_len]='\0';


not_nul_num_len=strlen(chstr_temp);  
  if (is_real_len == true) {desired_len=not_nul_num_len;}


int  empt_symb_len=desired_len-not_nul_num_len;
  if (empt_symb_len < 0) {empt_symb_len=0;}
int empt_symb_start_pos=0,  empt_symb_end_pos=empt_symb_start_pos+empt_symb_len-1;

int  num_symb_len=not_nul_num_len;
  if (num_symb_len < 0) {num_symb_len=0;}
int  num_symb_start_pos=empt_symb_end_pos+1; 
  if (empt_symb_len < 1) {num_symb_start_pos=0;}
int  num_symb_end_pos=num_symb_start_pos+num_symb_len-1; 
  

  if (fill_empt_symb_left1_right0 == false) 
  {
  num_symb_len=not_nul_num_len;  num_symb_start_pos=0;  num_symb_end_pos=num_symb_start_pos+num_symb_len-1;
  empt_symb_len=desired_len-not_nul_num_len;
    if (empt_symb_len < 0) {empt_symb_len=0;}  
  empt_symb_start_pos=num_symb_end_pos+1;  empt_symb_end_pos=empt_symb_start_pos+empt_symb_len-1;
  }


chstr_num[0]='\0';

  for (int pos=empt_symb_start_pos, counter=0;  counter < empt_symb_len;  counter++, pos++)
  {chstr_num[pos]=ch_empt_symb;}

  for (int pos=num_symb_start_pos, counter=0;  pos <= num_symb_end_pos;  counter++, pos++)
  {chstr_num[pos]=chstr_temp[pos-num_symb_start_pos];}

chstr_num[desired_len]='\0';

  if (chstr_temp) {delete[] chstr_temp;  chstr_temp=0;}


return strlen(chstr_num);

}




//  dyn variant for chstr_num,      convert mpz_t to chstr var:  necessary len will be set with symb ch_empt_symb       desired_len - len we want 
//            if (desired_len >  0) {it will be:  filling with ch_empt_symb of rest pos;}
//            if (desired_len == 0) {it will be: desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}
//            if (desired_len <  0) {it will be: desired_len=real_len;}  
//
int mpz_t__to__chstr__dyn(mpz_t * mpz_t_num,   char * & chstr_num,  int desired_len,  char ch_empt_symb,  bool fill_empt_symb_left1_right0)
{   

  if ((chstr_num) || (mpz_t_num == 0)) {return 0;}

int not_nul_num_len=mpz_sizeinbase(*mpz_t_num, 10);


  if (desired_len <  0) {desired_len=not_nul_num_len;}
  if (desired_len == 0) {desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}
  if (desired_len < not_nul_num_len) {desired_len=not_nul_num_len;}   
  

chstr_num=new char[desired_len+1];
chstr_num[0]='\0';
int len=mpz_t__to__chstr__stat(mpz_t_num,   chstr_num,   desired_len,  ch_empt_symb,   fill_empt_symb_left1_right0);

return len;

}




//  convert mpz_t to screen:  necessary len will be set with symb ch_empt_symb    (1) desired_len - len we want, if (desired_len <= 0) {it will be: desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}  
//            if (desired_len >  0) {it will be:  filling with ch_empt_symb of rest pos;}
//            if (desired_len == 0) {it will be: desired_len=__MACR_MPZ_MAX_CHSTR_SIZE;}
//            if (desired_len <  0) {it will be: desired_len=real_len;}  
//
void mpz_t__to__screen(mpz_t * mpz_t_num,   int desired_len,  char ch_empt_symb,  bool fill_empt_symb_left1_right0)
{   

char * chstr_num=0;

int len=mpz_t__to__chstr__dyn(mpz_t_num,   chstr_num,  desired_len,  ch_empt_symb,  fill_empt_symb_left1_right0);

  if (chstr_num) 
  {
    if (len > 0) {cout<<chstr_num;}

  delete[] chstr_num;  chstr_num=0;
  }

}




// ------------------------------------------------------------------------------------------------------------




unsigned long long int ullint_2_pow_n(int n)
{
    
return (  ((unsigned long long int) 1) << (n) );

}



unsigned int bit_inverse(unsigned int num,  int len)
{

bool left_bit_is_1=false, right_bit_is_1=false;
unsigned int bit_1_left=1, bit_1_right=1;
bit_1_left=(bit_1_left << (len-1));

int len_div_2=len/2;

  for (int i=0;  i < len_div_2;  i++)
  {
    if ((bit_1_left  & num) != 0) {left_bit_is_1 =true;} else {left_bit_is_1 =false;}
    if ((bit_1_right & num) != 0) {right_bit_is_1=true;} else {right_bit_is_1=false;}

    if (left_bit_is_1 != right_bit_is_1) 
    {
      if (right_bit_is_1 == true) {num=~((~num)+bit_1_right);} else {num=num+bit_1_right;}
      if (left_bit_is_1  == true) {num=~((~num)+bit_1_left); } else {num=num+bit_1_left; }
    }

  bit_1_right=(bit_1_right << 1);
  bit_1_left =(bit_1_left  >> 1); 
  }

return num;

}



unsigned long long int bit_inverse_ullint(unsigned long long int num,  int len)
{

bool left_bit_is_1=false, right_bit_is_1=false;
unsigned long long int bit_1_left=1, bit_1_right=1;
bit_1_left=(bit_1_left << (len-1));

int len_div_2=len/2;

  for (int i=0;  i < len_div_2;  i++)
  {
    if ((bit_1_left  & num) != 0) {left_bit_is_1 =true;} else {left_bit_is_1 =false;}
    if ((bit_1_right & num) != 0) {right_bit_is_1=true;} else {right_bit_is_1=false;}

    if (left_bit_is_1 != right_bit_is_1) 
    {
      if (right_bit_is_1 == true) {num=~((~num)+bit_1_right);} else {num=num+bit_1_right;}
      if (left_bit_is_1  == true) {num=~((~num)+bit_1_left); } else {num=num+bit_1_left; }
    }

  bit_1_right=(bit_1_right << 1);
  bit_1_left =(bit_1_left  >> 1); 
  }

return num;

}




inline void write_bits_uint(unsigned int & num_dest, int start_pos_dest,   unsigned int num_src, int start_pos_src,   int len)
{

unsigned int bit_1_src=1;
bit_1_src=bit_1_src << start_pos_src;
unsigned int bit_1_dest=1;
bit_1_dest=bit_1_dest << start_pos_dest;
bool src_bit_is_1=false, dest_bit_is_1=false;

  for (int i=0;  i < len;  i++)
  {
    if ((bit_1_src  & num_src ) != 0) { src_bit_is_1=true;} else { src_bit_is_1=false;}
    if ((bit_1_dest & num_dest) != 0) {dest_bit_is_1=true;} else {dest_bit_is_1=false;}

    if (src_bit_is_1 != dest_bit_is_1) 
    {
      if (src_bit_is_1  == true) {num_dest=num_dest+bit_1_dest;} else {num_dest=~((~num_dest)+bit_1_dest);}
    }

  bit_1_src =(bit_1_src  >> 1);
  bit_1_dest=(bit_1_dest >> 1);
  }

}




inline void write_bits_ullint(unsigned long long int & num_dest, int start_pos_dest,   unsigned long long int num_src, int start_pos_src,   int len)
{

unsigned long long int bit_1_src=1;
bit_1_src=bit_1_src << start_pos_src;
unsigned long long int bit_1_dest=1;
bit_1_dest=bit_1_dest << start_pos_dest;
bool src_bit_is_1=false, dest_bit_is_1=false;

  for (int i=0;  i < len;  i++)
  {
    if ((bit_1_src  & num_src ) != 0) { src_bit_is_1=true;} else { src_bit_is_1=false;}
    if ((bit_1_dest & num_dest) != 0) {dest_bit_is_1=true;} else {dest_bit_is_1=false;}

    if (src_bit_is_1 != dest_bit_is_1) 
    {
      if (src_bit_is_1  == true) {num_dest=num_dest+bit_1_dest;} else {num_dest=~((~num_dest)+bit_1_dest);}
    }

  bit_1_src =(bit_1_src  >> 1);
  bit_1_dest=(bit_1_dest >> 1);
  }

}




void fill_random_array_with_0_or_1(bool * const array,  int N_loc)
{

  if (array != 0) 
  {

    for (int i=0;  i < N_loc;  i++)
    {
    if (rand() % 2 == 0) {array[i]=false;} else {array[i]=true;} 
    }
  }  

}





inline int cycl_correction(int mean, int period)
{

  if (mean >= period) {return (mean % period);}

  if (mean < 0) {return (period-(abs(mean) % period));}

return mean;

}



void shift_J_bonds(bool *J_ver,  bool *J_hor,  int shift_dy,  int shift_dx,  int Ly,  int Lx)
{

int nouv_y=0, nouv_x=0;


bool *J_ver_temp=0;
J_ver_temp=new bool[(Ly+1)*Lx];

  for (int i=0;  i < Ly;  i++)
  {
  nouv_y=cycl_correction(i+shift_dy, Ly);
    for (int j=0;  j < Lx;  j++)
    {nouv_x=cycl_correction(j+shift_dx, Lx);  J_ver_temp[nouv_y*Lx+nouv_x]=J_ver[i*Lx+j];} 
  }

  for (int i=0;  i < Ly;  i++)
  {
    for (int j=0;  j < Lx;  j++)
    {J_ver[i*Lx+j]=J_ver_temp[i*Lx+j];} 
  }
  
  for (int j=0;  j < Lx;  j++)
  {J_ver[Ly*Lx+j]=J_ver[j];} 
  
  if (J_ver_temp) {delete[] J_ver_temp;  J_ver_temp=0;}  


bool *J_hor_temp=0;
J_hor_temp=new bool[Ly*(Lx+1)];

  for (int i=0;  i < Ly;  i++)
  {
  nouv_y=cycl_correction(i+shift_dy, Ly);
    for (int j=0;  j < Lx;  j++)
    {nouv_x=cycl_correction(j+shift_dx, Lx);  J_hor_temp[nouv_y*(Lx+1)+nouv_x]=J_hor[i*(Lx+1)+j];} 
  }

  for (int i=0;  i < Ly;  i++)
  {
    for (int j=0;  j < Lx;  j++)
    {J_hor[i*(Lx+1)+j]=J_hor_temp[i*(Lx+1)+j];} 
  }
  
  for (int j=0;  j < Ly;  j++)
  {J_hor[j*(Lx+1)+Lx]=J_hor[j*(Lx+1)];} 
    
  if (J_hor_temp) {delete[] J_hor_temp;  J_hor_temp=0;}  

}





//======================================================================================
//                             ______
//          |\_______________ (_____\\______________
//  HH======#H###############H#######################
//          ' ~""""""""""""""`##(_))#H\"""""Y########
//                            ))    \#H\       `"Y###
//                            "      }#H)
//
//=======================================================================================













/*

struct matr_t1
{

matr_t1();
void free();
~matr_t1() {free();}  
bool init(int size_y_loc, int size_x_loc, unsigned int seting_mean);

void update_min_max_for_adding_elem(int  index_y,  int index_x);
void set_val_for_all_el(unsigned int value);
void set_val_for_one_el(int  index_y,  int index_x,  unsigned int value);
void assign_val_from_one_el_to_other(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src);
void add_val_from_one_el_to_other(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src);
void swap(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src);
void swap_with_add(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src);


bool shift_elements(int dy, int dx);
inline bool shift_elements_spec_fast_case(int dy, int dx);
unsigned long long int get_element_as_ullint_if_possible(int  index_y,  int index_x);
int get_element_as_chstr_stat(int  index_y,  int index_x,  char * const chstr_matr_element,     int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);      //  //  ? ? ?
int get_element_as_chstr_dyn(int  index_y,  int index_x,  char * & chstr_matr_element,     int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);      //  //  ? ? ?
void display_on_screen_one_element(int  index_y,  int index_x,      int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);
void display_on_screen_whole_matr(char * interval_betw_elem=0,  int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);
int  save_to_file_whole_matr(const char * const chstr_filename,  char * interval_betw_elem=0,  int desired_len=0,  char ch_empt_symb=' ',  bool fill_empt_symb_left1_right0=true);


#ifdef __MACR_MATR_TYPE_1__MPZ
mpz_t * element;
#endif

int size_y,  size_x, size_N;
int min_not_0_x_index, max_not_0_x_index,  not_0_len_x;
int min_not_0_y_index, max_not_0_y_index,  not_0_len_y;

};
*/




matr_t1 :: matr_t1()
{

#ifdef __MACR_MATR_TYPE_1__MPZ
element=0;
#endif

size_y=0;  size_x=0;  size_N=0;
min_not_0_x_index=0;  max_not_0_x_index=0;  not_0_len_x=0;
min_not_0_y_index=0;  max_not_0_y_index=0;  not_0_len_y=0;

}




void matr_t1 :: free()
{


#ifdef __MACR_MATR_TYPE_1__MPZ

  if (element)
  {

    for (int i=0;  i < size_N;  i++)
    {  
    mpz_clear(element[i]);
    }

  delete[] element;  element=0;
  }

#endif



size_y=0;  size_x=0;  size_N=0;
min_not_0_x_index=0;  max_not_0_x_index=0;  not_0_len_x=0;
min_not_0_y_index=0;  max_not_0_y_index=0;  not_0_len_y=0;

}




bool matr_t1 :: init(int size_y_loc, int size_x_loc, unsigned int seting_mean)
{   

free(); 

size_y=size_y_loc;    size_x=size_x_loc;
size_N=size_y*size_x;  

  if (size_N == 0) {return false;}


#ifdef __MACR_MATR_TYPE_1__MPZ

element=new mpz_t[size_N];

  if (seting_mean != 0)  {min_not_0_y_index=0;  max_not_0_y_index=size_y-1;    min_not_0_x_index=0;  max_not_0_x_index=size_x-1;           not_0_len_y=size_y;  not_0_len_x=size_x;}

  for (int i=0;  i < size_N;  i++)
  {  
  mpz_init2(element[i], __MACR_MPZ_BIT_SIZE);
  mpz_set_ui(element[i], seting_mean);  //  cout<<"g: "<<element[i]<<" :g"<<endl;
  }

#endif


return true;

}




void matr_t1 :: update_min_max_for_adding_elem(int  index_y,  int index_x)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

  if (not_0_len_y == 0)
  {min_not_0_y_index=index_y;  max_not_0_y_index=index_y;  not_0_len_y=1;}
  else
  {
    if (index_y < min_not_0_y_index)
    {min_not_0_y_index=index_y;}
    else
    {
      if (index_y > max_not_0_y_index) {max_not_0_y_index=index_y;}
    }

  not_0_len_y=max_not_0_y_index-min_not_0_y_index+1;
  }

  if (not_0_len_x == 0)
  {min_not_0_x_index=index_x;  max_not_0_x_index=index_x;  not_0_len_x=1;}
  else
  {
    if (index_x < min_not_0_x_index)
    {min_not_0_x_index=index_x;}
    else
    {
      if (index_x > max_not_0_x_index) {max_not_0_x_index=index_x;}
    }

  not_0_len_x=max_not_0_x_index-min_not_0_x_index+1;
  }

#endif

}




void matr_t1 :: set_val_for_all_el(unsigned int value)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

  if (element)
  {
   
    for (int i=0;  i < size_N;  i++)
    {mpz_set_ui(element[i],  value);}

    if (value != 0) {min_not_0_y_index=0;  max_not_0_y_index=size_y-1;    min_not_0_x_index=0;  max_not_0_x_index=size_x-1;           not_0_len_y=size_y;  not_0_len_x=size_x;}
    else {min_not_0_y_index=0;  max_not_0_y_index=0;    min_not_0_x_index=0;  max_not_0_x_index=0;           not_0_len_y=0;  not_0_len_x=0;}
    //else {min_not_0_y_index=0;  max_not_0_y_index=size_y-1;    min_not_0_x_index=0;  max_not_0_x_index=size_x-1;           not_0_len_y=0;  not_0_len_x=0;}

  }

#endif

}




void matr_t1 :: set_zero_val_for_all_el()
{
unsigned int zero=0;

  if ((not_0_len_y != 0) && (not_0_len_x != 0))
  {
    for (int y=min_not_0_y_index;  y <= max_not_0_y_index;  y++)
    {
      for (int x=min_not_0_x_index;  x <= max_not_0_x_index;  x++)
      {
      int index=size_x*y+x;  
        if (mpz_sgn(element[index]) != 0) {mpz_set_ui(element[index], zero);}   
      }
    }

  min_not_0_y_index=0;  max_not_0_y_index=0;    min_not_0_x_index=0;  max_not_0_x_index=0;           not_0_len_y=0;  not_0_len_x=0;
  }

}




void matr_t1 :: set_val_for_one_el(int  index_y,  int index_x,  unsigned int value)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

  if (element)
  {
  mpz_set_ui(element[size_x*index_y+index_x],  value);
  
    if (value != 0) {update_min_max_for_adding_elem(index_y,  index_x);}

  }

#endif

}




void matr_t1 :: assign_val_from_one_el_to_other(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

  if (element)
  {
  mpz_set(element[size_x*index_y_dest+index_x_dest], element[size_x*index_y_src+index_x_src]);
  
    if (mpz_sgn(element[size_x*index_y_dest+index_x_dest]) != 0) {update_min_max_for_adding_elem(index_y_dest,  index_x_dest);}

  }

#endif

}



void matr_t1 :: add_val_from_one_el_to_other(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

  if (element)
  {
  mpz_add(element[size_x*index_y_dest+index_x_dest], element[size_x*index_y_dest+index_x_dest], element[size_x*index_y_src+index_x_src]);

    if (mpz_sgn(element[size_x*index_y_src+index_x_src]) != 0) {update_min_max_for_adding_elem(index_y_dest,  index_x_dest);}

  }

#endif

}



void matr_t1 :: mul_matr_by_num(unsigned int mult) 
{

#ifdef __MACR_MATR_TYPE_1__MPZ    

  if (element)
  {

    for (int y=min_not_0_y_index;  y <= max_not_0_y_index;  y++)
    {
      for (int x=min_not_0_x_index;  x <= max_not_0_x_index;  x++)
      {
      int index=size_x*y+x; 

        if (mpz_sgn(element[index]) != 0) {mpz_mul_ui(element[index],  element[index],  mult);}  

      }
    }

  }

#endif

}



void matr_t1 :: swap(int index_y_dest,  int index_x_dest,  int  index_y_src,  int index_x_src)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

  if (element)
  {
  mpz_swap(element[size_x*index_y_dest+index_x_dest], element[size_x*index_y_src+index_x_src]);
  }

#endif

}




bool matr_t1 :: shift_elements(int dy, int dx)
{

  if ((dy == 0) && (dx == 0)) {return false;}
  if ((not_0_len_y == 0) || (not_0_len_x == 0)) {return false;}  

#ifdef __MACR_MATR_TYPE_1__MPZ

int  y_start_src=0,  x_start_src=0;
int  y_pos_src=0,  x_pos_src=0;
int  y_pos_dest=0,  x_pos_dest=0;
int  dy_on_1=0,  dx_on_1=0; 

  if (dy < 0) {y_start_src=min_not_0_y_index;  dy_on_1=+1;} else {y_start_src=max_not_0_y_index;  dy_on_1=-1;}
  if (dx < 0) {x_start_src=min_not_0_x_index;  dx_on_1=+1;} else {x_start_src=max_not_0_x_index;  dx_on_1=-1;} 
  //for (int h=0; h < 36; h++) {mpz_set_ui(element[h], 0);}

int copy_size_y=max_not_0_y_index-min_not_0_y_index+1;
int copy_size_x=max_not_0_x_index-min_not_0_x_index+1;
  
y_pos_src=y_start_src;

  for (int y=0;  y < copy_size_y;  y++)
  {
  x_pos_src=x_start_src;
   
    for (int x=0;  x < copy_size_x;  x++)
    {
    y_pos_dest=y_pos_src+dy;  x_pos_dest=x_pos_src+dx;

      if (mpz_sgn(element[y_pos_src*size_x+x_pos_src]) != 0)
      {

        if ((y_pos_dest >= 0) && (y_pos_dest < size_y) && (x_pos_dest >= 0) && (x_pos_dest < size_x)) 
        {
        mpz_set(element[y_pos_dest*size_x+x_pos_dest], element[y_pos_src*size_x+x_pos_src]);
        }    
        #ifdef __MACR_CHECK_ACCESS_MEM_VIOLATION
        else
        { 
        //  error 
        display_on_screen_error_msg((char *) "in proc    matr_t1 :: shift_elements:  mem_access_viol");
        cout<<"error details start:"<<endl;   
        cout<<"dy="<<dy<<"  dx="<<dx<<endl;
        cout<<"min_not_0_y_index="<<min_not_0_y_index<<"  max_not_0_y_index="<<max_not_0_y_index<<"  min_not_0_x_index="<<min_not_0_x_index<<"  max_not_0_x_index="<<max_not_0_x_index<<endl;
        cout<<"y_pos_dest="<<y_pos_dest<<"  x_pos_dest="<<x_pos_dest<<"   y_pos_src="<<y_pos_src<<"  x_pos_src="<<x_pos_src<<endl;
        cout<<"error details end:"<<endl;  
        return false; 
        }
        #endif      
     
      mpz_set_ui(element[y_pos_src*size_x+x_pos_src], 0);
      }

    x_pos_src+=dx_on_1;
    }

  y_pos_src+=dy_on_1;
  }




min_not_0_y_index+=dy;  max_not_0_y_index+=dy;
min_not_0_x_index+=dx;  max_not_0_x_index+=dx;

  if ((min_not_0_y_index > size_y-1) || (max_not_0_y_index < 0)) {min_not_0_y_index=0;  max_not_0_y_index=0;  not_0_len_y=0;}
  if (min_not_0_y_index < 0) {min_not_0_y_index=0;  not_0_len_y=max_not_0_y_index-min_not_0_y_index+1;}
  if (max_not_0_y_index > size_y-1) {max_not_0_y_index=size_y-1;  not_0_len_y=max_not_0_y_index-min_not_0_y_index+1;}

  if ((min_not_0_x_index > size_x-1) || (max_not_0_x_index < 0)) {min_not_0_x_index=0;  max_not_0_x_index=0;  not_0_len_x=0;}
  if (min_not_0_x_index < 0) {min_not_0_x_index=0;  not_0_len_x=max_not_0_x_index-min_not_0_x_index+1;}
  if (max_not_0_x_index > size_x-1) {max_not_0_x_index=size_x-1;  not_0_len_x=max_not_0_x_index-min_not_0_x_index+1;}


#endif


return true;

}




//  in this case    dy  and  dx should be positive

/*bool matr_t1 :: shift_elements_spec_fast_case(int dy, int dx)
{  
  
#ifdef __MACR_MATR_TYPE_1__MPZ

int index_src=max_not_0_y_index*size_x+max_not_0_x_index;
  
  for (int area_index_y=0;  area_index_y < not_0_len_y;  area_index_y++)
  {
   
    for (int area_index_x=0;  area_index_x < not_0_len_x;  area_index_x++)
    {
      if (element[index_src] != 0)
      {
      mpz_set(element[index_src+dy*size_x+dx], element[index_src]);
      mpz_set_ui(element[index_src], (unsigned int) 0);
      }
    index_src--;
    }
  
  index_src-=size_x-not_0_len_x;
  }


min_not_0_y_index+=dy;  max_not_0_y_index+=dy;
min_not_0_x_index+=dx;  max_not_0_x_index+=dx;

#endif


return true;

}*/



//  shift on dy,  dx=1

inline bool matr_t1 :: shift_elements_spec_fast_case(int dy)
{  
  
#ifdef __MACR_MATR_TYPE_1__MPZ

int index_src=max_not_0_y_index*size_x+max_not_0_x_index;
int src_dif=-size_x+not_0_len_x;

  if (dy < 0) {index_src=min_not_0_y_index*size_x+max_not_0_x_index;  src_dif=size_x+not_0_len_x;} 
  
  for (int area_index_y=0;  area_index_y < not_0_len_y;  area_index_y++)
  {
   
    for (int area_index_x=0;  area_index_x < not_0_len_x;  area_index_x++)
    {
      if (mpz_sgn(element[index_src]) != 0)
      {
      //cout<<"dest_index="<<(index_src+dy*size_x+1)<<"  index_src="<<index_src<<endl;
      mpz_set(element[index_src+dy*size_x+1], element[index_src]);
      mpz_set_ui(element[index_src], (unsigned int) 0);
      }
    index_src--;
    }
  
  index_src+=src_dif;
  }


min_not_0_y_index+=dy;  max_not_0_y_index+=dy;
min_not_0_x_index++;    max_not_0_x_index++;

#endif


return true;

}




inline bool matr_t1 :: shift_elements_spec_fast_case_M_div_2(int dy)
{  
  
#ifdef __MACR_MATR_TYPE_1__MPZ

int index_src=min_not_0_y_index*size_x+max_not_0_x_index;

  if (max_not_0_x_index+1 >= size_x) 
  {
    for (int y_num=0;  y_num < not_0_len_y;  y_num++)
    {
    mpz_set_ui(element[index_src], (unsigned int) 0);
    index_src+=size_x;
    }

  
  max_not_0_x_index--;  not_0_len_x--;

    if (not_0_len_x <= 0) {min_not_0_y_index=0; max_not_0_y_index=0; not_0_len_y=0;    min_not_0_x_index=0; max_not_0_x_index=0; not_0_len_x=0;}

  }



  if ((not_0_len_y > 0) && (not_0_len_x > 0))
  {
  index_src=max_not_0_y_index*size_x+max_not_0_x_index;
  int src_dif=-size_x+not_0_len_x;

    if (dy < 0) {index_src=min_not_0_y_index*size_x+max_not_0_x_index;  src_dif=size_x+not_0_len_x;} 
  
    for (int area_index_y=0;  area_index_y < not_0_len_y;  area_index_y++)
    {
   
      for (int area_index_x=0;  area_index_x < not_0_len_x;  area_index_x++)
      {
        if (mpz_sgn(element[index_src]) != 0)
        {
        //cout<<"dest_index="<<(index_src+dy*size_x+1)<<"  index_src="<<index_src<<endl;
          if ((index_src+dy*size_x+1 >= 0) && (index_src+dy*size_x+1 < size_y*size_x)) {mpz_set(element[index_src+dy*size_x+1], element[index_src]);}
        mpz_set_ui(element[index_src], (unsigned int) 0);
        }
      index_src--;
      }
  
    index_src+=src_dif;
    }


  min_not_0_y_index+=dy;  max_not_0_y_index+=dy;
  min_not_0_x_index++;    max_not_0_x_index++;
  }

#endif


return true;

}




inline bool matr_t1 :: shift_elements_spec_fast_case_M_div_2_v20(int dy)
{  
  
#ifdef __MACR_MATR_TYPE_1__MPZ

int index_src=min_not_0_y_index*size_x+max_not_0_x_index;

  if (max_not_0_x_index+1 >= size_x) 
  {
    for (int y_num=0;  y_num < not_0_len_y;  y_num++)
    {
    mpz_set_ui(element[index_src], (unsigned int) 0);
    index_src+=size_x;
    }

  
  max_not_0_x_index--;  not_0_len_x--;

    if (not_0_len_x <= 0) {min_not_0_y_index=0; max_not_0_y_index=0; not_0_len_y=0;    min_not_0_x_index=0; max_not_0_x_index=0; not_0_len_x=0;}

  }


  if ((not_0_len_y > 0) && (not_0_len_x > 0))
  {

    if (dy >= 0)
    {
    int size=size_y*size_x;
    int index_dest=0;

      for (int area_index_y=0;  area_index_y < not_0_len_y;  area_index_y++)
      {
      index_src = (max_not_0_y_index-area_index_y   )*size_x+max_not_0_x_index;
      index_dest= (max_not_0_y_index-area_index_y+dy)*size_x+max_not_0_x_index+1;

        for (int area_index_x=0;  area_index_x < not_0_len_x;  area_index_x++)
        {
          if (mpz_sgn(element[index_src]) != 0)
          {
          //cout<<"dest_index="<<(index_src+dy*size_x+1)<<"  index_src="<<index_src<<endl;
            if ((index_dest >= 0) && (index_dest < size)) {mpz_set(element[index_dest], element[index_src]);}
          mpz_set_ui(element[index_src], (unsigned int) 0);
          }
        index_src--;
        index_dest--;
        }
  
      }

    }
    else
    {      //  if (dy < 0)
    int size=size_y*size_x;
    int index_dest=0;

      for (int area_index_y=0;  area_index_y < not_0_len_y;  area_index_y++)
      {
      index_src = (min_not_0_y_index+area_index_y   )*size_x+max_not_0_x_index;
      index_dest= (min_not_0_y_index+area_index_y+dy)*size_x+max_not_0_x_index+1;

        for (int area_index_x=0;  area_index_x < not_0_len_x;  area_index_x++)
        {
          if (mpz_sgn(element[index_src]) != 0)
          {
          //cout<<"dest_index="<<(index_src+dy*size_x+1)<<"  index_src="<<index_src<<endl;
            if ((index_dest >= 0) && (index_dest < size)) {mpz_set(element[index_dest], element[index_src]);}
          mpz_set_ui(element[index_src], (unsigned int) 0);
          }
        index_src--;
        index_dest--;
        }
  
      }


    }

  min_not_0_y_index+=dy;  max_not_0_y_index+=dy;
  min_not_0_x_index++;    max_not_0_x_index++;

    if (min_not_0_y_index < 0) {min_not_0_y_index=0;}  if (max_not_0_y_index >= size_y) {max_not_0_y_index=size_y-1;}  not_0_len_y=max_not_0_y_index-min_not_0_y_index+1; 
    if (min_not_0_x_index < 0) {min_not_0_x_index=0;}  if (max_not_0_x_index >= size_x) {max_not_0_x_index=size_x-1;}  not_0_len_x=max_not_0_x_index-min_not_0_x_index+1; 

  }


#endif


return true;

}




inline bool matr_t1 :: shift_elements_spec_without_M_v22(int dy)
{  
  
int index_src=0,  index_dest=0;
  if (dy == 0) {return true;}


  if ((not_0_len_y > 0) && (not_0_len_x > 0))
  {

    if (dy >= 0)
    {
    int size=size_y;
    int index_dest=0;

      for (int area_index_y=0;  area_index_y < not_0_len_y;  area_index_y++)
      {
      index_src = (max_not_0_y_index-area_index_y   );
      index_dest= (max_not_0_y_index-area_index_y+dy);

      int area_index_x=0;
        {
          if (mpz_sgn(element[index_src]) != 0)
          {
          //cout<<"dest_index="<<(index_src+dy*size_x+1)<<"  index_src="<<index_src<<endl;
            if ((index_dest >= 0) && (index_dest < size)) {mpz_set(element[index_dest], element[index_src]);}
          mpz_set_ui(element[index_src], (unsigned int) 0);
          }
        }
  
      }

    }
    else
    {      //  if (dy < 0)
    int size=size_y;
    int index_dest=0;

      for (int area_index_y=0;  area_index_y < not_0_len_y;  area_index_y++)
      {
      index_src = (min_not_0_y_index+area_index_y   );
      index_dest= (min_not_0_y_index+area_index_y+dy);

      int area_index_x=0;
        {
          if (mpz_sgn(element[index_src]) != 0)
          {
          //cout<<"dest_index="<<(index_src+dy*size_x+1)<<"  index_src="<<index_src<<endl;
            if ((index_dest >= 0) && (index_dest < size)) {mpz_set(element[index_dest], element[index_src]);}
          mpz_set_ui(element[index_src], (unsigned int) 0);
          }
        }
  
      }


    }

  min_not_0_y_index+=dy;  max_not_0_y_index+=dy;

    if (min_not_0_y_index < 0) {min_not_0_y_index=0;}  if (max_not_0_y_index >= size_y) {max_not_0_y_index=size_y-1;}  not_0_len_y=max_not_0_y_index-min_not_0_y_index+1; 
  min_not_0_x_index=0;  max_not_0_x_index=0;  not_0_len_x=1; 

  }



return true;

}




unsigned long long int matr_t1 :: get_element_as_ullint_if_possible(int index_y,  int index_x)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

mpz_t  mpz_t_max_ullint;
mpz_init2(mpz_t_max_ullint, __MACR_MPZ_BIT_SIZE); 
mpz_set_ui(mpz_t_max_ullint, ULLONG_MAX);  

  if (element[index_y*size_x+index_x] <= mpz_t_max_ullint)
  {
  char chstr_mpz_t[__MACR_MPZ_MAX_CHSTR_SIZE+1];
  chstr_mpz_t[0]='\0';
  mpz_get_str(chstr_mpz_t, 10, element[index_y*size_x+index_x]);
  chstr_mpz_t[__MACR_MPZ_MAX_CHSTR_SIZE]='\0';
  return atoi(chstr_mpz_t); 
  }
  else
  {
  return 7777777;
  }

#endif

return 0;

}




int matr_t1 :: get_element_as_chstr_stat(int  index_y,  int index_x,  char * const chstr_matr_element,     int desired_len,  char ch_empt_symb,  bool fill_empt_symb_left1_right0)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

return mpz_t__to__chstr__stat(  & element[index_y*size_x+index_x],   chstr_matr_element,  desired_len,  ch_empt_symb,  fill_empt_symb_left1_right0);

#endif

return 0;

}




int matr_t1 :: get_element_as_chstr_dyn(int  index_y,  int index_x,  char * & chstr_matr_element,     int desired_len,  char ch_empt_symb,  bool fill_empt_symb_left1_right0)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

return mpz_t__to__chstr__dyn(  & element[index_y*size_x+index_x],   chstr_matr_element,  desired_len,  ch_empt_symb,  fill_empt_symb_left1_right0);

#endif

return 0;

}




void matr_t1 :: display_on_screen_one_element(int  index_y,  int index_x,      int desired_len,  char ch_empt_symb,  bool fill_empt_symb_left1_right0)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

mpz_t__to__screen(  & element[index_y*size_x+index_x],   desired_len,  ch_empt_symb,  fill_empt_symb_left1_right0);

#endif

}



void matr_t1 :: display_on_screen_not_0_area_size_param()
{

cout<<"min_not_0_y_index="<<min_not_0_y_index<<"  max_not_0_y_index="<<max_not_0_y_index<<"  not_0_len_y="<<not_0_len_y<<"  size_y="<<size_y<<endl;
cout<<"min_not_0_x_index="<<min_not_0_x_index<<"  max_not_0_x_index="<<max_not_0_x_index<<"  not_0_len_x="<<not_0_len_x<<"  size_x="<<size_x<<endl;

}



void matr_t1 :: display_on_screen_whole_matr(char *interval_betw_elem, int desired_len,  char ch_empt_symb,  bool fill_empt_symb_left1_right0)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

char * chstr_interval_betw_elem_loc=0;

  if (interval_betw_elem) 
  {chstr_interval_betw_elem_loc=new char[strlen(interval_betw_elem)+1];  chstr_interval_betw_elem_loc[0]='\0';  strcpy(chstr_interval_betw_elem_loc, interval_betw_elem);}
  else
  {chstr_interval_betw_elem_loc=new char[2];  chstr_interval_betw_elem_loc[0]=' ';  chstr_interval_betw_elem_loc[1]='\0';}


cout<<"******************************************* matr *******************************************"<<endl;

cout<<"size_y="<<size_y<<"    "<<"size_x="<<size_x<<endl;

  for (int index_y_inver=0;  index_y_inver < size_y;  index_y_inver++)
  {
  int index_y=size_y-index_y_inver-1;

    for (int index_x=0;  index_x < size_x;  index_x++)
    {
    mpz_t__to__screen(  & element[index_y*size_x+index_x],   desired_len,  ch_empt_symb,  fill_empt_symb_left1_right0);
    cout<<chstr_interval_betw_elem_loc;
    }

  cout<<endl;
  }

cout<<"********************************************************************************************"<<endl;


  if (chstr_interval_betw_elem_loc) {delete[] chstr_interval_betw_elem_loc;  chstr_interval_betw_elem_loc=0;}

#endif

}




int matr_t1 :: save_to_file_whole_matr(const char * const chstr_filename,  char * interval_betw_elem,  int desired_len,  char ch_empt_symb,  bool fill_empt_symb_left1_right0)
{ 

#ifdef __MACR_MATR_TYPE_1__MPZ

  if (chstr_filename == 0) {return 0;}  

char * chstr_interval_betw_elem_loc=0;  

  if (interval_betw_elem) 
  {chstr_interval_betw_elem_loc=new char[strlen(interval_betw_elem)+1];  chstr_interval_betw_elem_loc[0]='\0';  strcpy(chstr_interval_betw_elem_loc, interval_betw_elem);}
  else
  {chstr_interval_betw_elem_loc=new char[2];  chstr_interval_betw_elem_loc[0]=' ';  chstr_interval_betw_elem_loc[1]='\0';}  


char * chstr_text=0;
string string_var;
char chstr_num[__MACR_MPZ_MAX_CHSTR_SIZE+1+1000];   chstr_num[0]='\0';


string_var.append((char *)  "******************************************* matr *******************************************\n");


  for (int index_y_inver=0;  index_y_inver < size_y;  index_y_inver++)
  {
  int index_y=size_y-index_y_inver-1;

    for (int index_x=0;  index_x < size_x;  index_x++)
    {
    get_element_as_chstr_stat(index_y,  index_x,  chstr_num,    0,  ch_empt_symb,  fill_empt_symb_left1_right0);
    string_var.append(chstr_num); 
    string_var.append(chstr_interval_betw_elem_loc); 
    }

  string_var.append((char *)  "\n");
  }  

string_var.append((char *)  "********************************************************************************************\n");


chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';



int len=write_chstr_to_file(chstr_text,  chstr_filename,  true);
        

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

  if (chstr_interval_betw_elem_loc) {delete[] chstr_interval_betw_elem_loc;  chstr_interval_betw_elem_loc=0;}

return len;

#endif


return 0;

}





//  Sum    elements matr x1 (len_y * len_x)   and    elements matr x2 (len_y * len_x).    Put result matr to result.
void mpz_t_matr_sum(mpz_t *result,  mpz_t *x1,  mpz_t *x2,  int len_y,  int len_x)
{

  for (int index_y=0;  index_y < len_y;  index_y++)
  {
    for (int index_x=0;  index_x < len_x;  index_x++)
    {
    mpz_add(result[len_x*index_y+index_x],  x1[len_x*index_y+index_x],  x2[len_x*index_y+index_x]);
    }
  }

}


//  Dif    elements matr x1 (len_y * len_x)   and    elements matr x2 (len_y * len_x).    Put result matr to result.
void mpz_t_matr_dif(mpz_t *result,  mpz_t *x1,  mpz_t *x2,  int len_y,  int len_x)
{

  for (int index_y=0;  index_y < len_y;  index_y++)
  {
    for (int index_x=0;  index_x < len_x;  index_x++)
    {
    mpz_sub(result[len_x*index_y+index_x],  x1[len_x*index_y+index_x],  x2[len_x*index_y+index_x]);
    }
  }

}




//  Shift    elements of matr (len_y * len_x) on shift_y and shift_x.    
void mpz_t_matr_shift(mpz_t *matr,  int len_y,  int len_x,  int shift_y,  int shift_x)
{

  if ((shift_y != 0) || (shift_x != 0))
  {
  int  y_start_src=0,  x_start_src=0;
  int  y_pos_src=0,  x_pos_src=0;
  int  y_pos_dest=0,  x_pos_dest=0;
  int  dy_on_1=0,  dx_on_1=0;

    if (shift_y < 0) {y_start_src=0;  dy_on_1=+1;} else {y_start_src=len_y-1;  dy_on_1=-1;}
    if (shift_x < 0) {x_start_src=0;  dx_on_1=+1;} else {x_start_src=len_x-1;  dx_on_1=-1;} 

  
  y_pos_src=y_start_src;

    for (int y=0;  y < len_y;  y++)
    {
    x_pos_src=x_start_src;
   
      for (int x=0;  x < len_x;  x++)
      {
      y_pos_dest=y_pos_src+shift_y;  x_pos_dest=x_pos_src+shift_x;

        if ((y_pos_dest >= 0) && (y_pos_dest < len_y) && (x_pos_dest >= 0) && (x_pos_dest < len_x)) 
        {mpz_set(matr[y_pos_dest*len_x+x_pos_dest], matr[y_pos_src*len_x+x_pos_src]);}
     
      mpz_set_ui(matr[y_pos_src*len_x+x_pos_src], 0);
      x_pos_src+=dx_on_1;
      }

    y_pos_src+=dy_on_1;
    }
  }

}






//  Add elements of matr_add_src into matr_dest.
//
inline void matr1_add(matr_t1 *matr_dest,  matr_t1 *matr_add_src)
{


  for (int y=matr_add_src->min_not_0_y_index;  y <= matr_add_src->max_not_0_y_index;  y++)
  {

    for (int x=matr_add_src->min_not_0_x_index;  x <= matr_add_src->max_not_0_x_index;  x++)
    {
    int index=matr_add_src->size_x*y+x;  
      if (mpz_sgn(matr_add_src->element[index]) != 0) {mpz_add(matr_dest->element[index], matr_dest->element[index], matr_add_src->element[index]);}
        
    }

  }

  if ((matr_dest->not_0_len_x == 0) || (matr_dest->not_0_len_y == 0))
  {
  matr_dest->min_not_0_y_index=matr_add_src->min_not_0_y_index;  matr_dest->max_not_0_y_index=matr_add_src->max_not_0_y_index;  matr_dest->not_0_len_y=matr_add_src->not_0_len_y;
  matr_dest->min_not_0_x_index=matr_add_src->min_not_0_x_index;  matr_dest->max_not_0_x_index=matr_add_src->max_not_0_x_index;  matr_dest->not_0_len_x=matr_add_src->not_0_len_x;
  }
  else
  {

    if (matr_dest->min_not_0_y_index > matr_add_src->min_not_0_y_index) {matr_dest->min_not_0_y_index=matr_add_src->min_not_0_y_index;}
    if (matr_dest->max_not_0_y_index < matr_add_src->max_not_0_y_index) {matr_dest->max_not_0_y_index=matr_add_src->max_not_0_y_index;}

    if (matr_dest->min_not_0_x_index > matr_add_src->min_not_0_x_index) {matr_dest->min_not_0_x_index=matr_add_src->min_not_0_x_index;}
    if (matr_dest->max_not_0_x_index < matr_add_src->max_not_0_x_index) {matr_dest->max_not_0_x_index=matr_add_src->max_not_0_x_index;}
  
  matr_dest->not_0_len_y=matr_dest->max_not_0_y_index-matr_dest->min_not_0_y_index+1;  if (matr_dest->not_0_len_y < 0) {matr_dest->not_0_len_y=0;}
  matr_dest->not_0_len_x=matr_dest->max_not_0_x_index-matr_dest->min_not_0_x_index+1;  if (matr_dest->not_0_len_x < 0) {matr_dest->not_0_len_x=0;}
  }

}






//  matr copy
//
inline void matr_t1_copy(matr_t1 *matr_dest,  matr_t1 *matr_src)
{

  if ((matr_src->not_0_len_x != 0) && (matr_src->not_0_len_y != 0))
  {
    for (int y=matr_src->min_not_0_y_index;  y <= matr_src->max_not_0_y_index;  y++)
    {
      for (int x=matr_src->min_not_0_x_index;  x <= matr_src->max_not_0_x_index;  x++)
      {
      int index=matr_src->size_x*y+x;  
        if (mpz_sgn(matr_src->element[index]) != 0) {mpz_set(matr_dest->element[index], matr_src->element[index]);}   
      }
    }
  }


  if ((matr_dest->not_0_len_x == 0) || (matr_dest->not_0_len_y == 0))
  {
  matr_dest->min_not_0_y_index=matr_src->min_not_0_y_index;  matr_dest->max_not_0_y_index=matr_src->max_not_0_y_index;  matr_dest->not_0_len_y=matr_src->not_0_len_y;
  matr_dest->min_not_0_x_index=matr_src->min_not_0_x_index;  matr_dest->max_not_0_x_index=matr_src->max_not_0_x_index;  matr_dest->not_0_len_x=matr_src->not_0_len_x;
  }
  else
  {

    if (matr_dest->min_not_0_y_index > matr_src->min_not_0_y_index) {matr_dest->min_not_0_y_index=matr_src->min_not_0_y_index;}
    if (matr_dest->max_not_0_y_index < matr_src->max_not_0_y_index) {matr_dest->max_not_0_y_index=matr_src->max_not_0_y_index;}

    if (matr_dest->min_not_0_x_index > matr_src->min_not_0_x_index) {matr_dest->min_not_0_x_index=matr_src->min_not_0_x_index;}
    if (matr_dest->max_not_0_x_index < matr_src->max_not_0_x_index) {matr_dest->max_not_0_x_index=matr_src->max_not_0_x_index;}
  
  matr_dest->not_0_len_y=matr_dest->max_not_0_y_index-matr_dest->min_not_0_y_index+1;  if (matr_dest->not_0_len_y < 0) {matr_dest->not_0_len_y=0;}
  matr_dest->not_0_len_x=matr_dest->max_not_0_x_index-matr_dest->min_not_0_x_index+1;  if (matr_dest->not_0_len_x < 0) {matr_dest->not_0_len_x=0;}
  }


}






//  matr copy with shift
//
inline void matr_t1_copy_with_shift(matr_t1 *matr_dest,  matr_t1 *matr_src, int shf_y, int shf_x)
{

  if ((matr_src->not_0_len_x != 0) && (matr_src->not_0_len_y != 0))
  {
  int dest_y=0, dest_x=0;
  int index_src=0, index_dest=0;

    for (int y=matr_src->min_not_0_y_index;  y <= matr_src->max_not_0_y_index;  y++)
    {
    dest_y=y+shf_y;
      for (int x=matr_src->min_not_0_x_index;  x <= matr_src->max_not_0_x_index;  x++)
      {
      index_src=matr_src->size_x*y+x;
      dest_x=x+shf_x;
        if ((0 <= dest_y) && (dest_y < matr_dest->size_y) && (0 <= dest_x) && (dest_x < matr_dest->size_x))
        {
        index_dest=matr_dest->size_x*dest_y+dest_x;
          if (mpz_sgn(matr_src->element[index_src]) != 0) {mpz_set(matr_dest->element[index_dest], matr_src->element[index_src]);}
        }   
      }
    }



  matr_dest->min_not_0_y_index=matr_src->min_not_0_y_index+shf_y;  matr_dest->max_not_0_y_index=matr_src->max_not_0_y_index+shf_y;
  matr_dest->min_not_0_x_index=matr_src->min_not_0_x_index+shf_x;  matr_dest->max_not_0_x_index=matr_src->max_not_0_x_index+shf_x;  

    if ((matr_dest->min_not_0_y_index >= matr_dest->size_y) || (matr_dest->min_not_0_x_index >= matr_dest->size_x) || (matr_dest->max_not_0_y_index < 0) || (matr_dest->max_not_0_x_index < 0)) 
    {
    matr_dest->min_not_0_y_index=0;  matr_dest->max_not_0_y_index=0;  matr_dest->not_0_len_y=0;
    matr_dest->min_not_0_x_index=0;  matr_dest->max_not_0_x_index=0;  matr_dest->not_0_len_x=0;
    }
    else
    {
      if (matr_dest->max_not_0_y_index >=  matr_dest->size_y) {matr_dest->max_not_0_y_index=matr_dest->size_y-1;}
      if (matr_dest->max_not_0_x_index >=  matr_dest->size_x) {matr_dest->max_not_0_x_index=matr_dest->size_x-1;}
      if (matr_dest->min_not_0_y_index <   0) {matr_dest->min_not_0_y_index=0;}
      if (matr_dest->min_not_0_x_index <   0) {matr_dest->min_not_0_x_index=0;}
 
    matr_dest->not_0_len_y=matr_dest->max_not_0_y_index-matr_dest->min_not_0_y_index+1;
    matr_dest->not_0_len_x=matr_dest->max_not_0_x_index-matr_dest->min_not_0_x_index+1;
    }

  }


}




//  Add elements of matr_add_src into matr_dest.
//
inline void matr1_add_dif_sizes(matr_t1 *matr_dest,  matr_t1 *matr_add_src)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

  for (int y=matr_add_src->min_not_0_y_index;  y <= matr_add_src->max_not_0_y_index;  y++)
  {

    for (int x=matr_add_src->min_not_0_x_index;  x <= matr_add_src->max_not_0_x_index;  x++)
    {
    int index_src =matr_add_src->size_x*y+x;
    int index_dest=matr_dest->size_x*y+x;

      if (mpz_sgn(matr_add_src->element[index_src]) != 0) {mpz_add(matr_dest->element[index_dest], matr_dest->element[index_dest], matr_add_src->element[index_src]);}
        
    }

  }

  if ((matr_dest->not_0_len_x == 0) || (matr_dest->not_0_len_y == 0))
  {
  matr_dest->min_not_0_y_index=matr_add_src->min_not_0_y_index;  matr_dest->max_not_0_y_index=matr_add_src->max_not_0_y_index;  matr_dest->not_0_len_y=matr_add_src->not_0_len_y;
  matr_dest->min_not_0_x_index=matr_add_src->min_not_0_x_index;  matr_dest->max_not_0_x_index=matr_add_src->max_not_0_x_index;  matr_dest->not_0_len_x=matr_add_src->not_0_len_x;
  }
  else
  {

    if (matr_dest->min_not_0_y_index > matr_add_src->min_not_0_y_index) {matr_dest->min_not_0_y_index=matr_add_src->min_not_0_y_index;}
    if (matr_dest->max_not_0_y_index < matr_add_src->max_not_0_y_index) {matr_dest->max_not_0_y_index=matr_add_src->max_not_0_y_index;}

    if (matr_dest->min_not_0_x_index > matr_add_src->min_not_0_x_index) {matr_dest->min_not_0_x_index=matr_add_src->min_not_0_x_index;}
    if (matr_dest->max_not_0_x_index < matr_add_src->max_not_0_x_index) {matr_dest->max_not_0_x_index=matr_add_src->max_not_0_x_index;}
  
  matr_dest->not_0_len_y=matr_dest->max_not_0_y_index-matr_dest->min_not_0_y_index+1;
  matr_dest->not_0_len_x=matr_dest->max_not_0_x_index-matr_dest->min_not_0_x_index+1;
  }

#endif



}






inline void matr1_area_is_assign_to_matr2(matr_t1 * matr_dest, int matr_dest_start_index,    matr_t1 *matr_src, int matr_src_start_index,    int matr_src_subarea_len_y, int matr_src_subarea_len_x)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

//int index_dest=matr_dest_start_index;
//int index_src=matr_src_start_index;
//
//  for (int area_index_y=0;  area_index_y=0;  area_index_y++)
//  {
//
//    for (int area_index_x=0;  area_index_x=0;  area_index_x++)
//    {
//    mpz_set(matr_dest->element[index_dest], matr_src->element[index_src]);
//    index_dest++;  index_src++;
//    }

//  index_dest+=matr_dest->size_x-matr_src_subarea_len_x;
//  index_src +=matr_src->size_x-matr_src_subarea_len_x;
//  }


  for (int area_index_y=0;  area_index_y=0;  area_index_y++)
  {

    for (int area_index_x=0;  area_index_x=0;  area_index_x++)
    {
    mpz_set(matr_dest->element[matr_dest_start_index], matr_src->element[matr_src_start_index]);
    matr_dest_start_index++;  matr_src_start_index++;
    }

  matr_dest_start_index+=matr_dest->size_x-matr_src_subarea_len_x;
  matr_src_start_index +=matr_src->size_x-matr_src_subarea_len_x;
  }

#endif

}




inline void matr1_area_is_added_to_matr2(matr_t1 * matr_dest, int matr_dest_start_index,    matr_t1 *matr_src, int matr_src_start_index,    int matr_src_subarea_len_y, int matr_src_subarea_len_x)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

//int index_dest=matr_dest_start_index;
//int index_src=matr_src_start_index;
//
//  for (int area_index_y=0;  area_index_y=0;  area_index_y++)
//  {
//
//    for (int area_index_x=0;  area_index_x=0;  area_index_x++)
//    {
//    mpz_set(matr_dest->element[index_dest], matr_src->element[index_src]);
//    index_dest++;  index_src++;
//    }

//  index_dest+=matr_dest->size_x-matr_src_subarea_len_x;
//  index_src +=matr_src->size_x-matr_src_subarea_len_x;
//  }


  for (int area_index_y=0;  area_index_y < matr_src_subarea_len_y;  area_index_y++)
  {

    for (int area_index_x=0;  area_index_x < matr_src_subarea_len_x;  area_index_x++)
    {
    //mpz_set(matr_dest->element[matr_dest_start_index], matr_src->element[matr_src_start_index]);
    mpz_add(matr_dest->element[matr_dest_start_index], matr_dest->element[matr_dest_start_index], matr_src->element[matr_src_start_index]);
    matr_dest_start_index++;  matr_src_start_index++;
    }

  matr_dest_start_index+=matr_dest->size_x-matr_src_subarea_len_x;
  matr_src_start_index +=matr_src->size_x-matr_src_subarea_len_x;
  }

#endif

}





//  ! ! ! not-0 elements are not modified of matr ! ! !
//                         
//                   
// 
//          F F * * *          * * * * *
//          F F * * *          * * * * *
//    matr1 F @ A A *          * * * * *  matr2
//    src   * A A A *          * ^ * * *  dest
//          * | * * *          * | * * *
//            |                  |
//            |__________________|
//                        add
//
//   A - coping area of matr1-src
//   F - forbidden  for coping area of matr1-src
//   ! ! ! not-0 elements are not modified of matr ! ! !
//
inline void matr1_area_spec_add_to_matr2(matr_t1 * matr_dest,  int matr_dest_start_pos,          matr_t1 *matr_src,  int matr_src_subarea_y0,  int matr_src_subarea_x0,  int matr_src_subarea_len_y,  int matr_src_subarea_len_x,  int matr_src_forbidsubarea_y0,  int matr_src_forbidsubarea_x0,  int matr_src_forbidsubarea_len_y,  int matr_src_forbidsubarea_len_x)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

int end_pos_src_y=matr_src_subarea_y0+matr_src_subarea_len_y-1,  end_pos_src_x=matr_src_subarea_x0+matr_src_subarea_len_x-1;

int pos_src_forbidarea_max_y=matr_src_forbidsubarea_y0+matr_src_forbidsubarea_len_y-1,   pos_src_forbidarea_max_x=matr_src_forbidsubarea_x0+matr_src_forbidsubarea_len_x-1;

  if (matr_src_forbidsubarea_len_y < 1) {matr_src_forbidsubarea_y0=matr_src_subarea_y0+matr_src_subarea_len_y;  pos_src_forbidarea_max_y=matr_src_subarea_y0-1;}
  if (matr_src_forbidsubarea_len_x < 1) {matr_src_forbidsubarea_x0=matr_src_subarea_x0+matr_src_subarea_len_x;  pos_src_forbidarea_max_x=matr_src_subarea_x0-1;}

int pos_dest_and_src_dif=matr_dest_start_pos-matr_src_subarea_y0*matr_src->size_x-matr_src_subarea_x0;     


  for (int pos_src_y=matr_src_subarea_y0;  pos_src_y <= end_pos_src_y;  pos_src_y++)
  {

    for (int pos_src_x=matr_src_subarea_x0;  pos_src_x <= end_pos_src_x;  pos_src_x++)
    {  
    int pos_src=pos_src_y*matr_src->size_x+pos_src_x;

      if (mpz_sgn(matr_src->element[pos_src]) != 0)
      {
        if ((matr_src_forbidsubarea_len_x < 1) || (pos_src_x < matr_src_forbidsubarea_x0) || (pos_src_forbidarea_max_x < pos_src_x) || (pos_src_y < matr_src_forbidsubarea_y0) || (pos_src_forbidarea_max_y < pos_src_y)) 
        {    
        int pos_dest=pos_src+pos_dest_and_src_dif;

        #ifdef __MACR_CHECK_ACCESS_MEM_VIOLATION
          if ((pos_dest >= matr_dest->size_y*matr_dest->size_x) || (pos_src >= matr_src->size_y*matr_src->size_x))
          { 
          //  error 
          display_on_screen_error_msg((char *) "in proc    matr1_area_spec_add_to_matr2:  mem_access_viol");
          cout<<"error details start:"<<endl;   
          cout<<"matr_dest->size_y="<<matr_dest->size_y<<"  matr_dest->size_x="<<matr_dest->size_x<<"  MATR1_SIZE=matr_1->size_y*matr_1->size_x="<<(matr_dest->size_y*matr_dest->size_x)<<endl;
          cout<<"matr_src->size_y="<<matr_src->size_y<<"  matr_src->size_x="<<matr_src->size_x<<"  MATR1_SIZE=matr_src->size_y*matr_src->size_x="<<(matr_src->size_y*matr_src->size_x)<<endl;
          cout<<"pos_dest="<<pos_dest<<"  pos_src="<<pos_src<<endl;
          cout<<"error details end:"<<endl;  
          //return false; 
          }
          else
        #endif
          {
          mpz_add(matr_dest->element[pos_dest],  matr_dest->element[pos_dest],  matr_src->element[pos_src]);
          }  
          
        }
      }

    }


  }

  /*for (int pos_src_y=matr_src_subarea_y0;  pos_src_y <= end_pos_src_y;  pos_src_y++)
  {

    if ((pos_src_y < matr_src_forbidsubarea_y0) || (pos_src_forbidarea_max_y < pos_src_y)) 
    {

      for (int pos_src_x=matr_src_subarea_x0;  pos_src_x <= end_pos_src_x;  pos_src_x++)
      {  

        if ((pos_src_x < matr_src_forbidsubarea_x0) || (pos_src_forbidarea_max_x < pos_src_x)) 
        {  
        int pos_src=pos_src_y*matr_src->size_x+pos_src_x;
          
          if (matr_src->element[pos_src] != 0)
          {
          int pos_dest=pos_src+pos_dest_and_src_dif;
          mpz_add(matr_dest->element[pos_dest],  matr_dest->element[pos_dest],  matr_src->element[pos_src]);          
          }

        }

      }

    }


  }*/



#endif

}





//
//                         add
//                   _________________
//                  |                |
//          * * * * V          * * * | *
//          * * * * *          * * * B *
//    matr1 * * * * *          * * * * *  matr2
//          * A * * *          * ^ * * *
//          * | * * *          * | * * *
//            | elem1            |
//            |__________________|
//                        add
//
//  we have 4 points: 2 src,  2 dest (where add to)
//  ! ! ! key problem is that any points can coincide ! ! !
//  it's the main reason why this proc is needed
//  ! ! ! not-0 elements are not modified of matr ! ! !
//
inline void two_matrices_mutual_add_with_1_elem(matr_t1 *matr1, int elem1_index_from, int elem1_index_to,    matr_t1 *matr2, int elem2_index_from, int elem2_index_to)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

//int  index_1_from=matr1->size_x*elem1_y_from+elem1_x_from,  index_1_to=matr1->size_x*elem1_y_to+elem1_x_to;
//int  index_2_from=matr2->size_x*elem2_y_from+elem2_x_from,  index_2_to=matr2->size_x*elem2_y_to+elem2_x_to;

  if (elem1_index_to != elem2_index_from)
  {   
  mpz_add(matr2->element[elem1_index_to],  matr2->element[elem1_index_to],  matr1->element[elem1_index_from]);
  mpz_add(matr1->element[elem2_index_to],  matr1->element[elem2_index_to],  matr2->element[elem2_index_from]);
  }
  else
  {
  mpz_t temp;
  mpz_init2(temp, __MACR_MPZ_BIT_SIZE);
  mpz_set(temp, matr1->element[elem1_index_from]);
  mpz_add(matr1->element[elem2_index_to],  matr1->element[elem2_index_to],  matr2->element[elem2_index_from]);
  mpz_add(matr2->element[elem1_index_to],  matr2->element[elem1_index_to],  temp);
  mpz_clear(temp);
  }

  //if ((elem1_index_to != elem2_index_from) && (elem2_index_to != elem1_index_from))
  //{ 
  //mpz_add(matr1->element[elem2_index_to],  matr1->element[elem2_index_to],  matr2->element[elem2_index_from]);
  //mpz_add(matr2->element[elem1_index_to],  matr2->element[elem1_index_to],  matr1->element[elem1_index_from]);
  //}
  //else
  //{
  //mpz_t temp;
  //mpz_init2(temp, __MACR_MPZ_BIT_SIZE);
  //mpz_set(temp, matr1->element[elem1_index_from]);
  //mpz_add(matr1->element[elem2_index_to],  matr1->element[elem2_index_to],  matr2->element[elem2_index_from]);
  //mpz_add(matr2->element[elem1_index_to],  matr2->element[elem1_index_to],  temp);
  //mpz_clear(temp);
  //}

#endif

}



//  swap with adding 2 the same quad chunks between 2 matrices 
//
//  using procedures:
//  1) void two_matrices_mutual_add_with_1_elem(matr_t1 *matr1, int elem1_index_from, int elem1_index_to,    matr_t1 *matr2, int elem2_index_from, int elem2_index_to)
//
//  * * * * *     * * * * *
//  * B B * *     * * * * *
//  * B B * *     * * A A *
//  * * * * *     * * A A *
//  * * * * *     * * * * *
//    ^ | ____add_____^ |
//    |_______add_______|
//
//  add assign crd are:      int y0_to_matr2, int x0_to_matr2                 int y0_to_matr1, int x0_to_matr1
//  src crd of chunk:        int y0_from_matr1, int x0_from_matr1             int y0_from_matr2, int x0_from_matr2
//  subarea sizes are the same
//  ! ! ! base 0-point is left down point ! ! !
//  ! ! ! matr1_index_from=matr2_index_to ! ! !
//  ! ! ! not-0 elements are not modified of matr ! ! !
//
//  or other scheme
//  matr1                    matr2
//  * * * * * *              * * * * * *
//        |____________________^   |
//        ^________________________|
//        int y0_matr2_to 
//        int x0_matr2_to
//
inline void matr1_area_spec_mutual_add_swap(matr_t1 * matr_1, int y0_matr1_to, int x0_matr1_to,     matr_t1 * matr_2, int y0_matr2_from, int x0_matr2_from, int y0_matr2_to, int x0_matr2_to,  int subarea_len_y, int subarea_len_x) 
{

#ifdef __MACR_MATR_TYPE_1__MPZ

int                     matr1_index_to;      //  matr1_index_from=matr2_index_to
int  matr2_index_from,  matr2_index_to;

matr2_index_from=y0_matr2_from*matr_2->size_x+x0_matr2_from;
matr2_index_to=y0_matr2_to*matr_1->size_x+x0_matr2_to;
matr1_index_to=y0_matr1_to*matr_2->size_x+x0_matr1_to;
 
int  index_dy,  index_dx;

  if (y0_matr1_to < y0_matr2_from) {index_dy=+1;} 
  else {index_dy=-1;  matr2_index_from+=(subarea_len_y-1)*matr_2->size_x;  matr1_index_to+=(subarea_len_y-1)*matr_2->size_x;  matr2_index_to+=(subarea_len_y-1)*matr_1->size_x;}

  if (x0_matr1_to < x0_matr2_from) {index_dx=+1;} 
  else {index_dx=-1;  matr2_index_from+=subarea_len_x-1;  matr1_index_to+=subarea_len_x-1;  matr2_index_to+=subarea_len_x-1;}

int  index_new_line_d=-subarea_len_x*index_dx+matr_1->size_x*index_dy; 
//==debug==//  cout<<"DEBUG: index_new_line_d="<<index_new_line_d<<endl;

  for (int num_y=0;  num_y < subarea_len_y;  num_y++)
  { 

    for (int num_x=0;  num_x < subarea_len_x;  num_x++)
    {

      if ((mpz_sgn(matr_1->element[matr2_index_to]) != 0) || (mpz_sgn(matr_2->element[matr2_index_from]) != 0))
      {

      #ifdef __MACR_CHECK_ACCESS_MEM_VIOLATION
        if ((matr2_index_to >= matr_1->size_y*matr_1->size_x) || (matr1_index_to >= matr_2->size_y*matr_2->size_x) || (matr2_index_from >= matr_2->size_y*matr_2->size_x))
        { 
        //  error 
        display_on_screen_error_msg((char *) "in proc    matr1_area_spec_mutual_add_swap:  mem_access_viol");
        cout<<"error details start:"<<endl;   
        cout<<"matr_1->size_y="<<matr_1->size_y<<"  matr_1->size_x="<<matr_1->size_x<<"  MATR1_SIZE=matr_1->size_y*matr_1->size_x="<<(matr_1->size_y*matr_1->size_x)<<endl;
        cout<<"matr_2->size_y="<<matr_2->size_y<<"  matr_2->size_x="<<matr_2->size_x<<"  MATR1_SIZE=matr_2->size_y*matr_1->size_x="<<(matr_1->size_y*matr_1->size_x)<<endl;
        cout<<"matr2_index_to="<<matr2_index_to<<"  matr1_index_to="<<matr1_index_to<<"  matr2_index_from="<<matr2_index_from<<endl;
        cout<<"error details end:"<<endl;  
        //return false; 
        }
        else
      #endif
        {
        two_matrices_mutual_add_with_1_elem(matr_1, matr2_index_to, matr1_index_to,    matr_2, matr2_index_from, matr2_index_to);
        }
   
      //==debug==//  cout<<"DEBUG: "<<matr2_index_to<<"  "<<matr1_index_to<<"  "<<matr2_index_from<<"  "<<matr2_index_to<<endl;
      }

    matr2_index_to+=index_dx;  matr2_index_from+=index_dx;  matr1_index_to+=index_dx;
    }

  matr2_index_to+=index_new_line_d;  matr2_index_from+=index_new_line_d;  matr1_index_to+=index_new_line_d;
  } 

#endif 

}



inline void matr1_area_spec_mutual_add_swap_v22(matr_t1 * matr_1, int y0_matr1_to, int x0_matr1_to,     matr_t1 * matr_2, int y0_matr2_from, int x0_matr2_from, int y0_matr2_to, int x0_matr2_to,  int subarea_len_y, int subarea_len_x) 
{

#ifdef __MACR_MATR_TYPE_1__MPZ

int                     matr1_index_to;      //  matr1_index_from=matr2_index_to
int  matr2_index_from,  matr2_index_to;

matr2_index_from=y0_matr2_from*matr_2->size_x+x0_matr2_from;
matr2_index_to=y0_matr2_to*matr_1->size_x+x0_matr2_to;
matr1_index_to=y0_matr1_to*matr_2->size_x+x0_matr1_to;
 
int  index_dy,  index_dx;

  if (y0_matr1_to < y0_matr2_from) {index_dy=+1;} 
  else {index_dy=-1;  matr2_index_from+=(subarea_len_y-1)*matr_2->size_x;  matr1_index_to+=(subarea_len_y-1)*matr_2->size_x;  matr2_index_to+=(subarea_len_y-1)*matr_1->size_x;}

//  if (x0_matr1_to < x0_matr2_from) {index_dx=+1;} 
//  else {index_dx=-1;  matr2_index_from+=subarea_len_x-1;  matr1_index_to+=subarea_len_x-1;  matr2_index_to+=subarea_len_x-1;}

int  index_new_line_d=-subarea_len_x*index_dx+matr_1->size_x*index_dy; 
//==debug==//  cout<<"DEBUG: index_new_line_d="<<index_new_line_d<<endl;
index_dx=0;

  for (int num_y=0;  num_y < subarea_len_y;  num_y++)
  { 

    {

      if ((mpz_sgn(matr_1->element[matr2_index_to]) != 0) || (mpz_sgn(matr_2->element[matr2_index_from]) != 0))
      {
        {
        two_matrices_mutual_add_with_1_elem(matr_1, matr2_index_to, matr1_index_to,    matr_2, matr2_index_from, matr2_index_to);
        }
   
      //==debug==//  cout<<"DEBUG: "<<matr2_index_to<<"  "<<matr1_index_to<<"  "<<matr2_index_from<<"  "<<matr2_index_to<<endl;
      }

    //matr2_index_to+=index_dx;  matr2_index_from+=index_dx;  matr1_index_to+=index_dx;
    }

  matr2_index_to+=index_new_line_d;  matr2_index_from+=index_new_line_d;  matr1_index_to+=index_new_line_d;
  } 

#endif 

}



//
//   ========= matr_cross_add ============
//
//  using procedures
//
//  level 1
//  -------
//  1) bool matr_t1 :: shift_elements_spec_fast_case(int dy, int dx)
//  2) void matr1_area_spec_add_to_matr2(matr_t1 * matr_dest,  int matr_dest_start_pos,          matr_t1 *matr_src,  int matr_src_subarea_y0,  int matr_src_subarea_x0,  int matr_src_subarea_len_y,  int matr_src_subarea_len_x,  int matr_src_forbidsubarea_y0,  int matr_src_forbidsubarea_x0,  int matr_src_forbidsubarea_len_y,  int matr_src_forbidsubarea_len_x)
//  3) matr1_area_spec_mutual_add_swap(matr_t1 * matr_1, int y0_matr1_to, int x0_matr1_to,     matr_t1 * matr_2, int y0_matr2_from, int x0_matr2_from, int y0_matr2_to, int x0_matr2_to,  int subarea_len_y, int subarea_len_x)
//
//  level 2
//  -------
//  1) void two_matrices_mutual_add_with_1_elem(matr_t1 *matr1, int elem1_index_from, int elem1_index_to,    matr_t1 *matr2, int elem2_index_from, int elem2_index_to)
//
//
//  common scheme
//   _____________________________________
//
//  | was 
//  |    matr1                    matr2
//  |    @                            @
//  |    * *                       x  x 
//  |    *    *_shift1          x     x_shift2
//  |    *       *           x        x
//  |    *          *     x           x
//  |    *             @              x
//  |    *          x     *           x
//  |    * _add_ x            * _add_ x 
//  |    *    x                   *   x
//  |    * x                         *x   
//  |    @                            @ 
//  |    matr1                    matr2
//  | became
//  |
//  V t
// 
//    shift up and right only (not left and not down)
//

//    
//   notation system
//
//     AA - matr1 aubarea    BBB - matr2 subarea 
//     AA   (sm1)            BBB   (sm2)
//
//   * A A * *             * * B B B
//   * A A * *             * * B B B
//   * * * * *             * * * * *
//   * * * * *             * * * * *
//   * * * * *             * * * * *
//   matr 1                matr 2
//
//   ism1_inm2:    image of matr 1 subarea in matr 2
//   ism2_inm1:    image of matr 2 subarea in matr 1
//
//   ism2acrsm1_inm1:    image of matr 1 subarea in matr 2
//   ism1acrsm2_inm2:    image of acrossing (subarea ofmatr 1 and subarea ofmatr 2   in matr 1)  in matr 2

//       ________________________
//      V                        |
//   * * * * *             * * * * *
//   * * * * *             * * f f *
//   * a a * *             * t @ f *
//   * a a * *             * t t * *
//   * * * * *             * * * * *
//      |____________________^
//
//  ! ! ! not-0 elements are not modified of matr ! ! !
//
void matr_cross_add_v1(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y,  int shf1_x,   int shf2_y,  int shf2_x)      // ! ! ! ======= very importand and difficult proc  ======= ! ! !
{

#ifdef __MACR_MATR_TYPE_1__MPZ

//  _(1)_

//matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x);
//==v1==//matr2->shift_elements(shf2_y,  shf2_x);
matr2->shift_elements_spec_fast_case(shf2_y);

//  2
//  define main parameters

int  ism2_inm1_y0=matr2->min_not_0_y_index-shf2_y,  ism2_inm1_y_end=ism2_inm1_y0+matr2->not_0_len_y-1,   ism2_inm1_size_y=ism2_inm1_y_end-ism2_inm1_y0+1;   //matr2->not_0_len_y;
int  ism2_inm1_x0=matr2->min_not_0_x_index-shf2_x,  ism2_inm1_x_end=ism2_inm1_x0+matr2->not_0_len_x-1,   ism2_inm1_size_x=ism2_inm1_x_end-ism2_inm1_x0+1;   //matr2->not_0_len_x;  
// debug // cout<<"debug:  ism2_inm1_y0="<<ism2_inm1_y0<<",  ism2_inm1_x0="<<ism2_inm1_x0<<",  ism2_inm1_y_end="<<ism2_inm1_y_end<<",  ism2_inm1_x_end="<<ism2_inm1_x_end<<endl;

int  ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y,  ism1_inm2_y_end=ism1_inm2_y0+matr1->not_0_len_y-1,   ism1_inm2_size_y=ism1_inm2_y_end-ism1_inm2_y0+1;      //matr1->not_0_len_y;
int  ism1_inm2_x0=matr1->min_not_0_x_index+shf1_x,  ism1_inm2_x_end=ism1_inm2_x0+matr1->not_0_len_x-1,   ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;      //matr1->not_0_len_x;
// debug // cout<<"debug:  ism1_inm2_y0="<<ism1_inm2_y0<<",  ism1_inm2_x0="<<ism1_inm2_x0<<",  ism1_inm2_y_end="<<ism1_inm2_y_end<<",  ism1_inm2_x_end="<<ism1_inm2_x_end<<endl;

//int  ism2acrsm1_inm1_y0=ism2_inm1_y0,      ism2acrsm1_inm1_y_end=ism1_inm2_y_end,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
//int  ism2acrsm1_inm1_x0=ism2_inm1_x0,      ism2acrsm1_inm1_x_end=ism1_inm2_x_end,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_y0=0,      ism2acrsm1_inm1_y_end=0,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_x0=0,      ism2acrsm1_inm1_x_end=0,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !

//int  ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//int  ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

int  ism2acrsm1_puttom2_y0=0,  ism2acrsm1_puttom2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_puttom2_x0=0,  ism2acrsm1_puttom2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity

int  ism2acrsm1_fromm2_y0=0,  ism2acrsm1_fromm2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_fromm2_x0=0,  ism2acrsm1_fromm2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity



//  across section handling
//
//  case (1)  *  case (2)  *  case (3)  *  case (4)  *  case (5)  *   here
//            *        |   *            *            *            *
//            *        |   *            *            *            *
//            *            *        |   *            *            *  matr1->max_not_0_y_index
//      |     *      |     *      | |   *      |     *      |     *  |  
//      |     *      |     *      |     *      |     *      | |   *  |
//      |     *      |     *      |     *      |     *      | |   *  |                           ism2acrsm1_inm1_y_end  ***************** ism2acrsm1_inm1_y_end (cross area)
//      |     *      |     *      |     *      | |   *      |     *  |                           |                      ***************** ism2acrsm1_inm1_y0    (cross area)
//            *            *            *        |   *            *  matr1->min_not_0_y_index    |
//        |   *            *            *            *            *                              ism2acrsm1_inm1_y0 
//        |   *            *            *            *            *                           
//is_across_ar*is_across_ar*is_across_ar*is_across_ar*is_across_ar*
// =false     * =false     * =true      * =true      * =true      *

bool is_across_area=true;

  if ((is_across_area == true) && (ism2_inm1_y0          <= matr1->max_not_0_y_index)) {ism2acrsm1_inm1_y0=ism2_inm1_y0;       } else {is_across_area=false;}
  if ((is_across_area == true) && (ism2_inm1_y_end       >= matr1->min_not_0_y_index)) {ism2acrsm1_inm1_y_end=ism2_inm1_y_end; } else {is_across_area=false;}
  if ((is_across_area == true) && (ism2_inm1_x0          <= matr1->max_not_0_x_index)) {ism2acrsm1_inm1_x0=ism2_inm1_x0;       } else {is_across_area=false;}
  if ((is_across_area == true) && (ism2_inm1_x_end       >= matr1->min_not_0_x_index)) {ism2acrsm1_inm1_x_end=ism2_inm1_x_end; } else {is_across_area=false;}

  if (is_across_area == true)
  {
  ism2acrsm1_inm1_y0=matr1->min_not_0_y_index;  ism2acrsm1_inm1_y_end=matr1->max_not_0_y_index;  ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1;
  ism2acrsm1_inm1_x0=matr1->min_not_0_x_index;  ism2acrsm1_inm1_x_end=matr1->max_not_0_x_index;  ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;

  ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_inm1_y_end+shf1_y;
  ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+shf1_x;    ism2acrsm1_puttom2_x_end=ism2acrsm1_inm1_x_end+shf1_x;

  ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
  ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
  }


//==//  //  it' the same but in more details
//==//  bool is_across_area=true;

//==//    if ((is_across_area == true) && (ism2_inm1_y0          <= matr1->max_not_0_y_index)) {ism2acrsm1_inm1_y0=ism2_inm1_y0;} else {is_across_area=false;}
//==//    if ((is_across_area == true) && (ism2acrsm1_inm1_y0    <  matr1->min_not_0_y_index)) {ism2acrsm1_inm1_y0=atr1->min_not_0_y_index;}
//==//    if ((is_across_area == true) && (ism2_inm1_y_end       >= matr1->min_not_0_y_index)) {ism2acrsm1_inm1_y_end=ism2_inm1_y_end;} else {is_across_area=false;}
//==//    if ((is_across_area == true) && (ism2acrsm1_inm1_y_end >  matr1->max_not_0_y_index)) {ism2acrsm1_inm1_y_end=matr1->max_not_0_y_index;}

//==//    if (is_across_area == true) {ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1;} else {ism2acrsm1_inm1_size_y=0; ism2acrsm1_inm1_y0=ism2acrsm1_inm1_y_end+1;}

//==//    if ((is_across_area == true) && (ism2_inm1_x0          <= matr1->max_not_0_x_index)) {ism2acrsm1_inm1_x0=ism2_inm1_x0;} else {is_across_area=false;}
//==//    if ((is_across_area == true) && (ism2acrsm1_inm1_x0    <  matr1->min_not_0_x_index)) {ism2acrsm1_inm1_x0=atr1->min_not_0_x_index;}
//==//    if ((is_across_area == true) && (ism2_inm1_x_end       >= matr1->min_not_0_x_index)) {ism2acrsm1_inm1_x_end=ism2_inm1_x_end;} else {is_across_area=false;}
//==//    if ((is_across_area == true) && (ism2acrsm1_inm1_x_end >  matr1->max_not_0_x_index)) {ism2acrsm1_inm1_x_end=matr1->max_not_0_x_index;}

//==//    if (is_across_area == true) {ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;} else {ism2acrsm1_inm1_size_x=0; ism2acrsm1_inm1_x0=ism2acrsm1_inm1_x_end+1;}

//==//    if (is_across_area == true)
//==//    {
//==//    ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_inm1_y_end+shf1_y;
//==//    ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+shf1_x;    ism2acrsm1_puttom2_x_end=ism2acrsm1_inm1_x_end+shf1_x;

//==//    ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
//==//    ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
//==//    }



/*  if (ism2_inm1_y0 < matr1->min_not_0_y_index   ) {ism2acrsm1_inm1_y0   =matr1->min_not_0_y_index;}       //  else {ism2acrsm1_inm1_y0=ism2_inm1_y0;      }    //  already done
  if (matr1->max_not_0_y_index < ism2_inm1_y_end) {ism2acrsm1_inm1_y_end=matr1->max_not_0_y_index;}       //  else {ism2acrsm1_inm1_y_end=ism1_inm2_y_end;}    //  already done

  if (ism2_inm1_x0 < matr1->min_not_0_x_index   ) {ism2acrsm1_inm1_x0   =matr1->min_not_0_x_index;}       //  else {ism2acrsm1_inm1_x0=ism2_inm1_x0;      }    //  already done
  if (matr1->max_not_0_x_index < ism2_inm1_x_end) {ism2acrsm1_inm1_x_end=matr1->max_not_0_x_index;}       //  else {ism2acrsm1_inm1_x_end=ism1_inm2_x_end;}    //  already done

ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1; cout<<"  kkk "<<ism2acrsm1_inm1_y_end<<"  "<<ism2acrsm1_inm1_y0<<"  kkk"<<endl;
ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;

  if ((ism2acrsm1_inm1_y_end < matr1->min_not_0_y_index) || (matr1->max_not_0_y_index < ism2acrsm1_inm1_y0)) {ism2acrsm1_inm1_size_y=0;}
  if ((ism2acrsm1_inm1_x_end < matr1->min_not_0_x_index) || (matr1->max_not_0_x_index < ism2acrsm1_inm1_x0)) {ism2acrsm1_inm1_size_x=0;}

  if ((ism2acrsm1_inm1_size_y > 0) && (ism2acrsm1_inm1_size_x > 0))
  {
  ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_inm1_y_end+shf1_y;
  ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+shf1_x;    ism2acrsm1_puttom2_x_end=ism2acrsm1_inm1_x_end+shf1_x;

  ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
  ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
  }
  else
  {
  ism2acrsm1_inm1_size_y=0;  ism2acrsm1_inm1_size_x=0;
  }*/


//  _(3)_   
matr1_area_spec_add_to_matr2(matr1,  ism2_inm1_y0*matr1->size_x+ism2_inm1_x0,          matr2,  matr2->min_not_0_y_index,  matr2->min_not_0_x_index,  ism2_inm1_size_y,  ism2_inm1_size_x,    ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);


//  _(4)_
//==debug==//  cout<<"y0_matr1_to="<<ism2acrsm1_puttom2_y0<<", x0_matr1_to="<<ism2acrsm1_puttom2_x0<<endl;
//==debug==//  cout<<"y0_matr2_from="<<ism2acrsm1_fromm2_y0<<", x0_matr2_from="<<ism2acrsm1_fromm2_x0<<endl;
//==debug==//  cout<<"y0_matr2_to="<<ism2acrsm1_inm1_y0<<", x0_matr2_to="<<ism2acrsm1_inm1_x0<<endl;
//==debug==//  cout<<"subarea_len_y="<<ism2acrsm1_inm1_size_y<<"  "<<"subarea_len_x="<<ism2acrsm1_inm1_size_x<<"  "<<endl;
matr1_area_spec_mutual_add_swap(matr1, ism2acrsm1_puttom2_y0, ism2acrsm1_puttom2_x0,     matr2, ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,                   ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,  ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);


//  _(5)_
matr1_area_spec_add_to_matr2(matr2,  ism1_inm2_y0*matr2->size_x+ism1_inm2_x0,          matr1,  matr1->min_not_0_y_index,  matr1->min_not_0_x_index,  ism1_inm2_size_y,  ism1_inm2_size_x,    ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);

//  6:    update not-0-area indicators 

//==//matr2->min_not_0_y_index+=shf2_y;  matr2->max_not_0_y_index+=shf2_y;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)
//==//matr2->min_not_0_x_index+=shf2_x;  matr2->max_not_0_x_index+=shf2_x;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)

  if (ism2_inm1_y0                < matr1->min_not_0_y_index ) {matr1->min_not_0_y_index=ism2_inm1_y0;   }
  if (matr1->max_not_0_y_index    < ism2_inm1_y_end          ) {matr1->max_not_0_y_index=ism2_inm1_y_end;}
  if (ism2_inm1_x0                < matr1->min_not_0_x_index ) {matr1->min_not_0_x_index=ism2_inm1_x0;   }
  if (matr1->max_not_0_x_index    < ism2_inm1_x_end          ) {matr1->max_not_0_x_index=ism2_inm1_x_end;}

  if (ism1_inm2_y0                < matr2->min_not_0_y_index ) {matr2->min_not_0_y_index=ism1_inm2_y0;   }
  if (matr2->max_not_0_y_index    < ism1_inm2_y_end          ) {matr2->max_not_0_y_index=ism1_inm2_y_end;}
  if (ism1_inm2_x0                < matr2->min_not_0_x_index ) {matr2->min_not_0_x_index=ism1_inm2_x0;   }
  if (matr2->max_not_0_x_index    < ism1_inm2_x_end          ) {matr2->max_not_0_x_index=ism1_inm2_x_end;}

matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1;
matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1;
matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;   


#ifdef __MACR_CHECK_ACCESS_MEM_VIOLATION

  if ((matr1->min_not_0_y_index < 0) || (matr1->min_not_0_y_index >= matr1->size_y) || (matr1->min_not_0_x_index < 0) || (matr1->min_not_0_x_index >= matr1->size_x) ||
      (matr1->max_not_0_y_index < 0) || (matr1->max_not_0_y_index >= matr1->size_y) || (matr1->max_not_0_x_index < 0) || (matr1->max_not_0_x_index >= matr1->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr1"<<endl;  
  cout<<"matr1:"<<endl;
  matr1->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr1->min_not_0_y_index=0;  matr1->max_not_0_y_index=matr1->size_y-1;  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1; 
  matr1->min_not_0_x_index=0;  matr1->max_not_0_x_index=matr1->size_x-1;  matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
  //return false; 
  }

  if ((matr2->min_not_0_y_index < 0) || (matr2->min_not_0_y_index >= matr2->size_y) || (matr2->min_not_0_x_index < 0) || (matr2->min_not_0_x_index >= matr2->size_x) ||
      (matr2->max_not_0_y_index < 0) || (matr2->max_not_0_y_index >= matr2->size_y) || (matr2->max_not_0_x_index < 0) || (matr2->max_not_0_x_index >= matr2->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr2"<<endl;  
  cout<<"matr2:"<<endl;
  matr2->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr2->min_not_0_y_index=0;  matr2->max_not_0_y_index=matr2->size_y-1;  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1; 
  matr2->min_not_0_x_index=0;  matr2->max_not_0_x_index=matr2->size_x-1;  matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;
  //return false; 
  }

#endif


#endif


//==//ism2_inm1_y0=matr1->min_not_0_y_index-shf2_y;  ism2_inm1_y_end=matr2->min_not_0_y_index-shf2_y;   ism2_inm1_size_y=matr2->not_0_len_y;
//==//ism2_inm1_x0=matr1->min_not_0_x_index-shf2_x;  ism2_inm1_x_end=matr2->min_not_0_x_index-shf2_x;   ism2_inm1_size_x=matr2->not_0_len_x;

//==//ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y;  ism1_inm2_y_end=matr1->max_not_0_y_index+shf1_y;   ism1_inm2_size_y=matr1->not_0_len_y;
//==//ism1_inm2_x0=matr1->min_not_0_x_index+shf1_x;  ism1_inm2_x_end=matr1->max_not_0_x_index+shf1_x;   ism1_inm2_size_x=matr1->not_0_len_x;


//==//ism2acrsm1_inm1_y0=;  ism2acrsm1_inm1_y_end=;       ism2acrsm1_inm1_size_y=;
//==//ism2acrsm1_inm1_x0=;  ism2acrsm1_inm1_x_end=;       ism2acrsm1_inm1_size_x=;

//==//ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//==//ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

//==//ism2acrsm1_inm2_y0=;  ism2acrsm1_inm2_y_end=;      // ism2acrsm1_inm2_size_y=;
//==//ism2acrsm1_inm2_x0=;  ism2acrsm1_inm2_x_end=;      // ism2acrsm1_inm2_size_x=;



}



//  --------------------------------------------------------------------------------------------------------------------------------------
//  ======================================================================================================================================
//  --------------------------------------------------------------------------------------------------------------------------------------


void matr_cross_add_v2(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y)      // ! ! ! ======= very importand and difficult proc  ======= ! ! !
{

//int shf1_y,  int shf1_x,   int shf2_y,  int shf2_x

#ifdef __MACR_MATR_TYPE_1__MPZ

//  _(1)_

//matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x);
//==v1==//matr2->shift_elements(shf2_y,  shf2_x);
matr2->shift_elements_spec_fast_case(shf1_y-1);

//  2
//  define main parameters

int  ism2_inm1_y0=matr2->min_not_0_y_index-(shf1_y-1),  ism2_inm1_y_end=ism2_inm1_y0+matr2->not_0_len_y-1,   ism2_inm1_size_y=ism2_inm1_y_end-ism2_inm1_y0+1;   //matr2->not_0_len_y;
int  ism2_inm1_x0=matr2->min_not_0_x_index-1,           ism2_inm1_x_end=ism2_inm1_x0+matr2->not_0_len_x-1,   ism2_inm1_size_x=ism2_inm1_x_end-ism2_inm1_x0+1;   //matr2->not_0_len_x;  
// debug // cout<<"debug:  ism2_inm1_y0="<<ism2_inm1_y0<<",  ism2_inm1_x0="<<ism2_inm1_x0<<",  ism2_inm1_y_end="<<ism2_inm1_y_end<<",  ism2_inm1_x_end="<<ism2_inm1_x_end<<endl;

int  ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y,  ism1_inm2_y_end=ism1_inm2_y0+matr1->not_0_len_y-1,   ism1_inm2_size_y=ism1_inm2_y_end-ism1_inm2_y0+1;      //matr1->not_0_len_y;
int  ism1_inm2_x0=matr1->min_not_0_x_index+1,       ism1_inm2_x_end=ism1_inm2_x0+matr1->not_0_len_x-1,   ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;      //matr1->not_0_len_x;
// debug // cout<<"debug:  ism1_inm2_y0="<<ism1_inm2_y0<<",  ism1_inm2_x0="<<ism1_inm2_x0<<",  ism1_inm2_y_end="<<ism1_inm2_y_end<<",  ism1_inm2_x_end="<<ism1_inm2_x_end<<endl;

//int  ism2acrsm1_inm1_y0=ism2_inm1_y0,      ism2acrsm1_inm1_y_end=ism1_inm2_y_end,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
//int  ism2acrsm1_inm1_x0=ism2_inm1_x0,      ism2acrsm1_inm1_x_end=ism1_inm2_x_end,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_y0=0,      ism2acrsm1_inm1_y_end=0,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_x0=0,      ism2acrsm1_inm1_x_end=0,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !

//int  ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//int  ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

int  ism2acrsm1_puttom2_y0=0,  ism2acrsm1_puttom2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_puttom2_x0=0,  ism2acrsm1_puttom2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity

int  ism2acrsm1_fromm2_y0=0,  ism2acrsm1_fromm2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_fromm2_x0=0,  ism2acrsm1_fromm2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity



//  across section handling
//
//  case (1)  *  case (2)  *  case (3)  *  case (4)  *  case (5)  *   here
//            *        |   *            *            *            *
//            *        |   *            *            *            *
//            *            *        |   *            *            *  matr1->max_not_0_y_index
//      |     *      |     *      | |   *      |     *      |     *  |  
//      |     *      |     *      |     *      |     *      | |   *  |
//      |     *      |     *      |     *      |     *      | |   *  |                           ism2acrsm1_inm1_y_end  ***************** ism2acrsm1_inm1_y_end (cross area)
//      |     *      |     *      |     *      | |   *      |     *  |                           |                      ***************** ism2acrsm1_inm1_y0    (cross area)
//            *            *            *        |   *            *  matr1->min_not_0_y_index    |
//        |   *            *            *            *            *                              ism2acrsm1_inm1_y0 
//        |   *            *            *            *            *                           
//is_across_ar*is_across_ar*is_across_ar*is_across_ar*is_across_ar*
// =false     * =false     * =true      * =true      * =true      *

bool is_across_area=true;

  if ((is_across_area == true) && (ism2_inm1_y0          <= matr1->max_not_0_y_index)) {ism2acrsm1_inm1_y0=ism2_inm1_y0;       } else {is_across_area=false;}
  if ((is_across_area == true) && (ism2_inm1_y_end       >= matr1->min_not_0_y_index)) {ism2acrsm1_inm1_y_end=ism2_inm1_y_end; } else {is_across_area=false;}
  if ((is_across_area == true) && (ism2_inm1_x0          <= matr1->max_not_0_x_index)) {ism2acrsm1_inm1_x0=ism2_inm1_x0;       } else {is_across_area=false;}
  if ((is_across_area == true) && (ism2_inm1_x_end       >= matr1->min_not_0_x_index)) {ism2acrsm1_inm1_x_end=ism2_inm1_x_end; } else {is_across_area=false;}

  if (is_across_area == true)
  {
  ism2acrsm1_inm1_y0=matr1->min_not_0_y_index;  ism2acrsm1_inm1_y_end=matr1->max_not_0_y_index;  ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1;
  ism2acrsm1_inm1_x0=matr1->min_not_0_x_index;  ism2acrsm1_inm1_x_end=matr1->max_not_0_x_index;  ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;

  ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_inm1_y_end+shf1_y;
  ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+1;         ism2acrsm1_puttom2_x_end=ism2acrsm1_inm1_x_end+1;

  ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
  ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
  }


//==//  //  it' the same but in more details
//==//  bool is_across_area=true;

//==//    if ((is_across_area == true) && (ism2_inm1_y0          <= matr1->max_not_0_y_index)) {ism2acrsm1_inm1_y0=ism2_inm1_y0;} else {is_across_area=false;}
//==//    if ((is_across_area == true) && (ism2acrsm1_inm1_y0    <  matr1->min_not_0_y_index)) {ism2acrsm1_inm1_y0=atr1->min_not_0_y_index;}
//==//    if ((is_across_area == true) && (ism2_inm1_y_end       >= matr1->min_not_0_y_index)) {ism2acrsm1_inm1_y_end=ism2_inm1_y_end;} else {is_across_area=false;}
//==//    if ((is_across_area == true) && (ism2acrsm1_inm1_y_end >  matr1->max_not_0_y_index)) {ism2acrsm1_inm1_y_end=matr1->max_not_0_y_index;}

//==//    if (is_across_area == true) {ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1;} else {ism2acrsm1_inm1_size_y=0; ism2acrsm1_inm1_y0=ism2acrsm1_inm1_y_end+1;}

//==//    if ((is_across_area == true) && (ism2_inm1_x0          <= matr1->max_not_0_x_index)) {ism2acrsm1_inm1_x0=ism2_inm1_x0;} else {is_across_area=false;}
//==//    if ((is_across_area == true) && (ism2acrsm1_inm1_x0    <  matr1->min_not_0_x_index)) {ism2acrsm1_inm1_x0=atr1->min_not_0_x_index;}
//==//    if ((is_across_area == true) && (ism2_inm1_x_end       >= matr1->min_not_0_x_index)) {ism2acrsm1_inm1_x_end=ism2_inm1_x_end;} else {is_across_area=false;}
//==//    if ((is_across_area == true) && (ism2acrsm1_inm1_x_end >  matr1->max_not_0_x_index)) {ism2acrsm1_inm1_x_end=matr1->max_not_0_x_index;}

//==//    if (is_across_area == true) {ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;} else {ism2acrsm1_inm1_size_x=0; ism2acrsm1_inm1_x0=ism2acrsm1_inm1_x_end+1;}

//==//    if (is_across_area == true)
//==//    {
//==//    ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_inm1_y_end+shf1_y;
//==//    ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+1;         ism2acrsm1_puttom2_x_end=ism2acrsm1_inm1_x_end+1;

//==//    ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
//==//    ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
//==//    }



/*  if (ism2_inm1_y0 < matr1->min_not_0_y_index   ) {ism2acrsm1_inm1_y0   =matr1->min_not_0_y_index;}       //  else {ism2acrsm1_inm1_y0=ism2_inm1_y0;      }    //  already done
  if (matr1->max_not_0_y_index < ism2_inm1_y_end) {ism2acrsm1_inm1_y_end=matr1->max_not_0_y_index;}       //  else {ism2acrsm1_inm1_y_end=ism1_inm2_y_end;}    //  already done

  if (ism2_inm1_x0 < matr1->min_not_0_x_index   ) {ism2acrsm1_inm1_x0   =matr1->min_not_0_x_index;}       //  else {ism2acrsm1_inm1_x0=ism2_inm1_x0;      }    //  already done
  if (matr1->max_not_0_x_index < ism2_inm1_x_end) {ism2acrsm1_inm1_x_end=matr1->max_not_0_x_index;}       //  else {ism2acrsm1_inm1_x_end=ism1_inm2_x_end;}    //  already done

ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1; cout<<"  kkk "<<ism2acrsm1_inm1_y_end<<"  "<<ism2acrsm1_inm1_y0<<"  kkk"<<endl;
ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;

  if ((ism2acrsm1_inm1_y_end < matr1->min_not_0_y_index) || (matr1->max_not_0_y_index < ism2acrsm1_inm1_y0)) {ism2acrsm1_inm1_size_y=0;}
  if ((ism2acrsm1_inm1_x_end < matr1->min_not_0_x_index) || (matr1->max_not_0_x_index < ism2acrsm1_inm1_x0)) {ism2acrsm1_inm1_size_x=0;}

  if ((ism2acrsm1_inm1_size_y > 0) && (ism2acrsm1_inm1_size_x > 0))
  {
  ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_inm1_y_end+shf1_y;
  ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+1;         ism2acrsm1_puttom2_x_end=ism2acrsm1_inm1_x_end+1;

  ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
  ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
  }
  else
  {
  ism2acrsm1_inm1_size_y=0;  ism2acrsm1_inm1_size_x=0;
  }*/


//  _(3)_   
matr1_area_spec_add_to_matr2(matr1,  ism2_inm1_y0*matr1->size_x+ism2_inm1_x0,          matr2,  matr2->min_not_0_y_index,  matr2->min_not_0_x_index,  ism2_inm1_size_y,  ism2_inm1_size_x,    ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);


//  _(4)_
//==debug==//  cout<<"y0_matr1_to="<<ism2acrsm1_puttom2_y0<<", x0_matr1_to="<<ism2acrsm1_puttom2_x0<<endl;
//==debug==//  cout<<"y0_matr2_from="<<ism2acrsm1_fromm2_y0<<", x0_matr2_from="<<ism2acrsm1_fromm2_x0<<endl;
//==debug==//  cout<<"y0_matr2_to="<<ism2acrsm1_inm1_y0<<", x0_matr2_to="<<ism2acrsm1_inm1_x0<<endl;
//==debug==//  cout<<"subarea_len_y="<<ism2acrsm1_inm1_size_y<<"  "<<"subarea_len_x="<<ism2acrsm1_inm1_size_x<<"  "<<endl;
matr1_area_spec_mutual_add_swap(matr1, ism2acrsm1_puttom2_y0, ism2acrsm1_puttom2_x0,     matr2, ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,                   ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,  ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);


//  _(5)_
matr1_area_spec_add_to_matr2(matr2,  ism1_inm2_y0*matr2->size_x+ism1_inm2_x0,          matr1,  matr1->min_not_0_y_index,  matr1->min_not_0_x_index,  ism1_inm2_size_y,  ism1_inm2_size_x,    ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);

//  6:    update not-0-area indicators 

//==//matr2->min_not_0_y_index+=shf1_y-1;  matr2->max_not_0_y_index+=shf1_y;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf1_y-1,  1)
//==//matr2->min_not_0_x_index+=1;         matr2->max_not_0_x_index+=1;        //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf1_y-1,  1)

  if (ism2_inm1_y0                < matr1->min_not_0_y_index ) {matr1->min_not_0_y_index=ism2_inm1_y0;   }
  if (matr1->max_not_0_y_index    < ism2_inm1_y_end          ) {matr1->max_not_0_y_index=ism2_inm1_y_end;}
  if (ism2_inm1_x0                < matr1->min_not_0_x_index ) {matr1->min_not_0_x_index=ism2_inm1_x0;   }
  if (matr1->max_not_0_x_index    < ism2_inm1_x_end          ) {matr1->max_not_0_x_index=ism2_inm1_x_end;}

  if (ism1_inm2_y0                < matr2->min_not_0_y_index ) {matr2->min_not_0_y_index=ism1_inm2_y0;   }
  if (matr2->max_not_0_y_index    < ism1_inm2_y_end          ) {matr2->max_not_0_y_index=ism1_inm2_y_end;}
  if (ism1_inm2_x0                < matr2->min_not_0_x_index ) {matr2->min_not_0_x_index=ism1_inm2_x0;   }
  if (matr2->max_not_0_x_index    < ism1_inm2_x_end          ) {matr2->max_not_0_x_index=ism1_inm2_x_end;}

matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1;
matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1;
matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;   


#ifdef __MACR_CHECK_ACCESS_MEM_VIOLATION

  if ((matr1->min_not_0_y_index < 0) || (matr1->min_not_0_y_index >= matr1->size_y) || (matr1->min_not_0_x_index < 0) || (matr1->min_not_0_x_index >= matr1->size_x) ||
      (matr1->max_not_0_y_index < 0) || (matr1->max_not_0_y_index >= matr1->size_y) || (matr1->max_not_0_x_index < 0) || (matr1->max_not_0_x_index >= matr1->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr1"<<endl;  
  cout<<"matr1:"<<endl;
  matr1->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr1->min_not_0_y_index=0;  matr1->max_not_0_y_index=matr1->size_y-1;  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1; 
  matr1->min_not_0_x_index=0;  matr1->max_not_0_x_index=matr1->size_x-1;  matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
  //return false; 
  }

  if ((matr2->min_not_0_y_index < 0) || (matr2->min_not_0_y_index >= matr2->size_y) || (matr2->min_not_0_x_index < 0) || (matr2->min_not_0_x_index >= matr2->size_x) ||
      (matr2->max_not_0_y_index < 0) || (matr2->max_not_0_y_index >= matr2->size_y) || (matr2->max_not_0_x_index < 0) || (matr2->max_not_0_x_index >= matr2->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr2"<<endl;  
  cout<<"matr2:"<<endl;
  matr2->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr2->min_not_0_y_index=0;  matr2->max_not_0_y_index=matr2->size_y-1;  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1; 
  matr2->min_not_0_x_index=0;  matr2->max_not_0_x_index=matr2->size_x-1;  matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;
  //return false; 
  }

#endif


#endif


//==//ism2_inm1_y0=matr1->min_not_0_y_index-(shf1_y-1);  ism2_inm1_y_end=matr2->min_not_0_y_index-(shf1_y-1);   ism2_inm1_size_y=matr2->not_0_len_y;
//==//ism2_inm1_x0=matr1->min_not_0_x_index-1;           ism2_inm1_x_end=matr2->min_not_0_x_index-1;            ism2_inm1_size_x=matr2->not_0_len_x;

//==//ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y;  ism1_inm2_y_end=matr1->max_not_0_y_index+shf1_y;   ism1_inm2_size_y=matr1->not_0_len_y;
//==//ism1_inm2_x0=matr1->min_not_0_x_index+1;       ism1_inm2_x_end=matr1->max_not_0_x_index+1;        ism1_inm2_size_x=matr1->not_0_len_x;


//==//ism2acrsm1_inm1_y0=;  ism2acrsm1_inm1_y_end=;       ism2acrsm1_inm1_size_y=;
//==//ism2acrsm1_inm1_x0=;  ism2acrsm1_inm1_x_end=;       ism2acrsm1_inm1_size_x=;

//==//ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//==//ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

//==//ism2acrsm1_inm2_y0=;  ism2acrsm1_inm2_y_end=;      // ism2acrsm1_inm2_size_y=;
//==//ism2acrsm1_inm2_x0=;  ism2acrsm1_inm2_x_end=;      // ism2acrsm1_inm2_size_x=;



}




//  --------------------------------------------------------------------------------------------------------------------------------------
//  ======================================================================================================================================
//  --------------------------------------------------------------------------------------------------------------------------------------



void matr_cross_add_v3(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y)      // ! ! ! ======= very importand and difficult proc  ======= ! ! !
{
//==DEBUG==//  cout<<"DEBUG:   matr_cross_add_v3"<<endl;  //==DEBUG==//  

#ifdef __MACR_MATR_TYPE_1__MPZ

int   shf1_x=1,  shf2_y=shf1_y-2,  shf2_x=1;

//  _(3)_  
//matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x);
//==v1==//matr2->shift_elements(shf2_y,  shf2_x);
//==DEBUG==//  cout<<"DEBUG:   "<<"  shf2_y="<<shf2_y<<endl;  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  WAS:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  WAS:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
  if ((matr2->not_0_len_y != 0) && (matr2->not_0_len_x != 0)) {matr2->shift_elements_spec_fast_case_M_div_2(shf2_y);}      //  !!! ======== v3 M_div_2 line modif ======== !!!
//==DEBUG==//  cout<<"DEBUG  INTERMED:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  INTERMED:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  

     

//  1
//  define main parameters

//==DEBUG==//  cout<<"debug:  matr2->not_0_len_x="<<matr2->not_0_len_x<<endl;  //==DEBUG==//   
int  ism2_inm1_y0=matr2->min_not_0_y_index-shf2_y,  ism2_inm1_y_end=ism2_inm1_y0+matr2->not_0_len_y-1,   ism2_inm1_size_y=ism2_inm1_y_end-ism2_inm1_y0+1;   //matr2->not_0_len_y;
int  ism2_inm1_x0=matr2->min_not_0_x_index-shf2_x,  ism2_inm1_x_end=ism2_inm1_x0+matr2->not_0_len_x-1,   ism2_inm1_size_x=ism2_inm1_x_end-ism2_inm1_x0+1;   //matr2->not_0_len_x; 
  if ((matr2->not_0_len_y == 0) || (matr2->not_0_len_x == 0)) 
  {
  ism2_inm1_y0=0;  ism2_inm1_y_end=0;   ism2_inm1_size_y=0;   //matr2->not_0_len_y;
  ism2_inm1_x0=0;  ism2_inm1_x_end=0;   ism2_inm1_size_x=0;   //matr2->not_0_len_x; 
  }  
//==DEBUG==//  cout<<"debug:  ism2_inm1_y0="<<ism2_inm1_y0<<",  ism2_inm1_x0="<<ism2_inm1_x0<<",  ism2_inm1_y_end="<<ism2_inm1_y_end<<",  ism2_inm1_x_end="<<ism2_inm1_x_end<<endl;  //==DEBUG==//  


int  ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y,  ism1_inm2_y_end=ism1_inm2_y0+matr1->not_0_len_y-1,   ism1_inm2_size_y=ism1_inm2_y_end-ism1_inm2_y0+1;      //matr1->not_0_len_y;
int  ism1_inm2_x0=matr1->min_not_0_x_index+shf1_x,  ism1_inm2_x_end=ism1_inm2_x0+matr1->not_0_len_x-1,   ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;      //matr1->not_0_len_x;
  if ((matr1->not_0_len_y == 0) || (matr1->not_0_len_x == 0)) 
  {
  ism1_inm2_y0=0;  ism1_inm2_y_end=0;   ism1_inm2_size_y=0;   //matr2->not_0_len_y;
  ism1_inm2_x0=0;  ism1_inm2_x_end=0;   ism1_inm2_size_x=0;   //matr2->not_0_len_x; 
  } 

  if (ism1_inm2_x_end >= matr2->size_x) {ism1_inm2_x_end=matr2->size_x-1; ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;}       //  !!! ======== v3 M_div_2 line ======== !!!
  if (ism1_inm2_size_x < 1) {ism1_inm2_size_x=0;  ism1_inm2_x0=0;  ism1_inm2_x_end=0;}                                            //  !!! ======== v3 M_div_2 line ======== !!!

//==DEBUG==//   cout<<"debug:  ism1_inm2_y0="<<ism1_inm2_y0<<",  ism1_inm2_x0="<<ism1_inm2_x0<<",  ism1_inm2_y_end="<<ism1_inm2_y_end<<",  ism1_inm2_x_end="<<ism1_inm2_x_end<<endl;  //==DEBUG==//  

//int  ism2acrsm1_inm1_y0=ism2_inm1_y0,      ism2acrsm1_inm1_y_end=ism1_inm2_y_end,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
//int  ism2acrsm1_inm1_x0=ism2_inm1_x0,      ism2acrsm1_inm1_x_end=ism1_inm2_x_end,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_y0=0,      ism2acrsm1_inm1_y_end=0,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_x0=0,      ism2acrsm1_inm1_x_end=0,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !

//int  ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//int  ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

int  ism2acrsm1_puttom2_y0=0,  ism2acrsm1_puttom2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_puttom2_x0=0,  ism2acrsm1_puttom2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity

int  ism2acrsm1_fromm2_y0=0,  ism2acrsm1_fromm2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_fromm2_x0=0,  ism2acrsm1_fromm2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity



//  across section handling
//
//  case (1)  *  case (2)  *  case (3)  *  case (4)  *  case (5)  *   here
//            *        |   *            *            *            *
//            *        |   *            *            *            *
//            *            *        |   *            *            *  matr1->max_not_0_y_index
//      |     *      |     *      | |   *      |     *      |     *  |  
//      |     *      |     *      |     *      |     *      | |   *  |
//      |     *      |     *      |     *      |     *      | |   *  |                           ism2acrsm1_inm1_y_end  ***************** ism2acrsm1_inm1_y_end (cross area)
//      |     *      |     *      |     *      | |   *      |     *  |                           |                      ***************** ism2acrsm1_inm1_y0    (cross area)
//            *            *            *        |   *            *  matr1->min_not_0_y_index    |
//        |   *            *            *            *            *                              ism2acrsm1_inm1_y0 
//        |   *            *            *            *            *                           
//is_across_ar*is_across_ar*is_across_ar*is_across_ar*is_across_ar*
// =false     * =false     * =true      * =true      * =true      *

bool is_across_area=true;

  if ((ism2_inm1_size_y > 0) && (ism2_inm1_size_x > 0))        //  !!! ======== v3 M_div_2 line ======== !!!
  {
    if ((is_across_area == true) && (ism2_inm1_y0          <= matr1->max_not_0_y_index)) {ism2acrsm1_inm1_y0=ism2_inm1_y0;       } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_y_end       >= matr1->min_not_0_y_index)) {ism2acrsm1_inm1_y_end=ism2_inm1_y_end; } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_x0          <= matr1->max_not_0_x_index)) {ism2acrsm1_inm1_x0=ism2_inm1_x0;       } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_x_end       >= matr1->min_not_0_x_index)) {ism2acrsm1_inm1_x_end=ism2_inm1_x_end; } else {is_across_area=false;}

    if (is_across_area == true)
    {
      if (ism2acrsm1_inm1_y0    < matr1->min_not_0_y_index) {ism2acrsm1_inm1_y0   =matr1->min_not_0_y_index;}
      if (ism2acrsm1_inm1_y_end > matr1->max_not_0_y_index) {ism2acrsm1_inm1_y_end=matr1->max_not_0_y_index;}  
    ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1;

      if (ism2acrsm1_inm1_x0    < matr1->min_not_0_x_index) {ism2acrsm1_inm1_x0   =matr1->min_not_0_x_index;}
      if (ism2acrsm1_inm1_x_end > matr1->max_not_0_x_index) {ism2acrsm1_inm1_x_end=matr1->max_not_0_x_index;}  
    ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;

    ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_puttom2_y0+ism2acrsm1_inm1_size_y-1;
    ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+shf1_x;    ism2acrsm1_puttom2_x_end=ism2acrsm1_puttom2_x0+ism2acrsm1_inm1_size_x-1;

    ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
    ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
    }
  }
  else
  {is_across_area=false;}



//  _(2)_   
  if ((ism2_inm1_size_y != 0) && (ism2_inm1_size_x != 0))
  {
  matr1_area_spec_add_to_matr2(matr1,  ism2_inm1_y0*matr1->size_x+ism2_inm1_x0,          matr2,  matr2->min_not_0_y_index,  matr2->min_not_0_x_index,  ism2_inm1_size_y,  ism2_inm1_size_x,    ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  }

//  _(4)_
//==DEBUG==//    cout<<"y0_matr1_to="<<ism2acrsm1_puttom2_y0<<", x0_matr1_to="<<ism2acrsm1_puttom2_x0<<endl;
//==DEBUG==//    cout<<"y0_matr2_from="<<ism2acrsm1_fromm2_y0<<", x0_matr2_from="<<ism2acrsm1_fromm2_x0<<endl;
//==DEBUG==//    cout<<"y0_matr2_to="<<ism2acrsm1_inm1_y0<<", x0_matr2_to="<<ism2acrsm1_inm1_x0<<endl;
//==DEBUG==//    cout<<"subarea_len_y="<<ism2acrsm1_inm1_size_y<<"  "<<"subarea_len_x="<<ism2acrsm1_inm1_size_x<<"  "<<endl;
  if ((ism2acrsm1_inm1_size_y !=0) && (ism2acrsm1_inm1_size_x !=0))
  { 
  matr1_area_spec_mutual_add_swap(matr1, ism2acrsm1_puttom2_y0, ism2acrsm1_puttom2_x0,     matr2, ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,                   ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,  ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  } 

//  _(5)_
  if (( ism1_inm2_size_y != 0) && (ism1_inm2_size_x != 0))
  {
  matr1_area_spec_add_to_matr2(matr2,  ism1_inm2_y0*matr2->size_x+ism1_inm2_x0,          matr1,  matr1->min_not_0_y_index,  matr1->min_not_0_x_index,  ism1_inm2_size_y,  ism1_inm2_size_x,    ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  }

//  6:    update not-0-area indicators 

//==//matr2->min_not_0_y_index+=shf2_y;  matr2->max_not_0_y_index+=shf2_y;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)
//==//matr2->min_not_0_x_index+=shf2_x;  matr2->max_not_0_x_index+=shf2_x;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)


  if ((ism2_inm1_size_y != 0) && (ism2_inm1_size_x != 0))
  {
    if ((matr1->not_0_len_y != 0) && (matr1->not_0_len_x != 0))
    {
      if (ism2_inm1_y0                < matr1->min_not_0_y_index ) {matr1->min_not_0_y_index=ism2_inm1_y0;   }
      if (matr1->max_not_0_y_index    < ism2_inm1_y_end          ) {matr1->max_not_0_y_index=ism2_inm1_y_end;}
      if (ism2_inm1_x0                < matr1->min_not_0_x_index ) {matr1->min_not_0_x_index=ism2_inm1_x0;   }
      if (matr1->max_not_0_x_index    < ism2_inm1_x_end          ) {matr1->max_not_0_x_index=ism2_inm1_x_end;}
    }
    else {matr1->min_not_0_y_index=ism2_inm1_y0;  matr1->max_not_0_y_index=ism2_inm1_y_end;     matr1->min_not_0_x_index=ism2_inm1_x0;  matr1->max_not_0_x_index=ism2_inm1_x_end;}

  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1;
  matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
  }


  if ((ism1_inm2_size_y != 0) && (ism1_inm2_size_x != 0))
  {
    if ((matr2->not_0_len_y != 0) && (matr2->not_0_len_x != 0))
    {
      if (ism1_inm2_y0                < matr2->min_not_0_y_index ) {matr2->min_not_0_y_index=ism1_inm2_y0;   }
      if (matr2->max_not_0_y_index    < ism1_inm2_y_end          ) {matr2->max_not_0_y_index=ism1_inm2_y_end;}
      if (ism1_inm2_x0                < matr2->min_not_0_x_index ) {matr2->min_not_0_x_index=ism1_inm2_x0;   }
      if (matr2->max_not_0_x_index    < ism1_inm2_x_end          ) {matr2->max_not_0_x_index=ism1_inm2_x_end;}
    }
    else
    {matr2->min_not_0_y_index=ism1_inm2_y0;  matr2->max_not_0_y_index=ism1_inm2_y_end;    matr2->min_not_0_x_index=ism1_inm2_x0;  matr2->max_not_0_x_index=ism1_inm2_x_end;}

  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1;
  matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;   
  }


//==DEBUG==//  cout<<"DEBUG  became:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  became:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  


#ifdef __MACR_CHECK_ACCESS_MEM_VIOLATION

  if ((matr1->min_not_0_y_index < 0) || (matr1->min_not_0_y_index >= matr1->size_y) || (matr1->min_not_0_x_index < 0) || (matr1->min_not_0_x_index >= matr1->size_x) ||
      (matr1->max_not_0_y_index < 0) || (matr1->max_not_0_y_index >= matr1->size_y) || (matr1->max_not_0_x_index < 0) || (matr1->max_not_0_x_index >= matr1->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr1"<<endl;  
  cout<<"matr1:"<<endl;
  matr1->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr1->min_not_0_y_index=0;  matr1->max_not_0_y_index=matr1->size_y-1;  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1; 
  matr1->min_not_0_x_index=0;  matr1->max_not_0_x_index=matr1->size_x-1;  matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
  //return false; 
  }

  if ((matr2->min_not_0_y_index < 0) || (matr2->min_not_0_y_index >= matr2->size_y) || (matr2->min_not_0_x_index < 0) || (matr2->min_not_0_x_index >= matr2->size_x) ||
      (matr2->max_not_0_y_index < 0) || (matr2->max_not_0_y_index >= matr2->size_y) || (matr2->max_not_0_x_index < 0) || (matr2->max_not_0_x_index >= matr2->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr2"<<endl;  
  cout<<"matr2:"<<endl;
  matr2->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr2->min_not_0_y_index=0;  matr2->max_not_0_y_index=matr2->size_y-1;  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1; 
  matr2->min_not_0_x_index=0;  matr2->max_not_0_x_index=matr2->size_x-1;  matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;
  //return false; 
  }

#endif


#endif


//==//ism2_inm1_y0=matr1->min_not_0_y_index-shf2_y;  ism2_inm1_y_end=matr2->min_not_0_y_index-shf2_y;   ism2_inm1_size_y=matr2->not_0_len_y;
//==//ism2_inm1_x0=matr1->min_not_0_x_index-shf2_x;  ism2_inm1_x_end=matr2->min_not_0_x_index-shf2_x;   ism2_inm1_size_x=matr2->not_0_len_x;

//==//ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y;  ism1_inm2_y_end=matr1->max_not_0_y_index+shf1_y;   ism1_inm2_size_y=matr1->not_0_len_y;
//==//ism1_inm2_x0=matr1->min_not_0_x_index+shf1_x;  ism1_inm2_x_end=matr1->max_not_0_x_index+shf1_x;   ism1_inm2_size_x=matr1->not_0_len_x;


//==//ism2acrsm1_inm1_y0=;  ism2acrsm1_inm1_y_end=;       ism2acrsm1_inm1_size_y=;
//==//ism2acrsm1_inm1_x0=;  ism2acrsm1_inm1_x_end=;       ism2acrsm1_inm1_size_x=;

//==//ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//==//ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

//==//ism2acrsm1_inm2_y0=;  ism2acrsm1_inm2_y_end=;      // ism2acrsm1_inm2_size_y=;
//==//ism2acrsm1_inm2_x0=;  ism2acrsm1_inm2_x_end=;      // ism2acrsm1_inm2_size_x=;


}




//  --------------------------------------------------------------------------------------------------------------------------------------
//  ======================================================================================================================================
//  --------------------------------------------------------------------------------------------------------------------------------------



void matr_cross_add_v4(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y,  int shf2_y)      // ! ! ! ======= very importand and difficult proc  ======= ! ! !
{
//==DEBUG==//  cout<<"DEBUG:   matr_cross_add_v3"<<endl;  //==DEBUG==//  

#ifdef __MACR_MATR_TYPE_1__MPZ

int  shf1_x=1, shf2_x=1;

//  _(3)_  
//matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x);
//==v1==//matr2->shift_elements(shf2_y,  shf2_x);
//==DEBUG==//  cout<<"DEBUG:   "<<"  shf2_y="<<shf2_y<<endl;  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  WAS:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  WAS:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
  if ((matr2->not_0_len_y != 0) && (matr2->not_0_len_x != 0)) {matr2->shift_elements_spec_fast_case_M_div_2(shf2_y);}      //  !!! ======== v3 M_div_2 line modif ======== !!!
//==DEBUG==//  cout<<"DEBUG  INTERMED:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  INTERMED:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  

     

//  1
//  define main parameters

//==DEBUG==//  cout<<"debug:  matr2->not_0_len_x="<<matr2->not_0_len_x<<endl;  //==DEBUG==//   
int  ism2_inm1_y0=matr2->min_not_0_y_index-shf2_y,  ism2_inm1_y_end=ism2_inm1_y0+matr2->not_0_len_y-1,   ism2_inm1_size_y=ism2_inm1_y_end-ism2_inm1_y0+1;   //matr2->not_0_len_y;
int  ism2_inm1_x0=matr2->min_not_0_x_index-shf2_x,  ism2_inm1_x_end=ism2_inm1_x0+matr2->not_0_len_x-1,   ism2_inm1_size_x=ism2_inm1_x_end-ism2_inm1_x0+1;   //matr2->not_0_len_x; 
  if ((matr2->not_0_len_y == 0) || (matr2->not_0_len_x == 0)) 
  {
  ism2_inm1_y0=0;  ism2_inm1_y_end=0;   ism2_inm1_size_y=0;   //matr2->not_0_len_y;
  ism2_inm1_x0=0;  ism2_inm1_x_end=0;   ism2_inm1_size_x=0;   //matr2->not_0_len_x; 
  }  
//==DEBUG==//  cout<<"debug:  ism2_inm1_y0="<<ism2_inm1_y0<<",  ism2_inm1_x0="<<ism2_inm1_x0<<",  ism2_inm1_y_end="<<ism2_inm1_y_end<<",  ism2_inm1_x_end="<<ism2_inm1_x_end<<endl;  //==DEBUG==//  


int  ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y,  ism1_inm2_y_end=ism1_inm2_y0+matr1->not_0_len_y-1,   ism1_inm2_size_y=ism1_inm2_y_end-ism1_inm2_y0+1;      //matr1->not_0_len_y;
int  ism1_inm2_x0=matr1->min_not_0_x_index+shf1_x,  ism1_inm2_x_end=ism1_inm2_x0+matr1->not_0_len_x-1,   ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;      //matr1->not_0_len_x;
  if ((matr1->not_0_len_y == 0) || (matr1->not_0_len_x == 0)) 
  {
  ism1_inm2_y0=0;  ism1_inm2_y_end=0;   ism1_inm2_size_y=0;   //matr2->not_0_len_y;
  ism1_inm2_x0=0;  ism1_inm2_x_end=0;   ism1_inm2_size_x=0;   //matr2->not_0_len_x; 
  } 

  if (ism1_inm2_x_end >= matr2->size_x) {ism1_inm2_x_end=matr2->size_x-1; ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;}       //  !!! ======== v3 M_div_2 line ======== !!!
  if (ism1_inm2_size_x < 1) {ism1_inm2_size_x=0;  ism1_inm2_x0=0;  ism1_inm2_x_end=0;}                                            //  !!! ======== v3 M_div_2 line ======== !!!

//==DEBUG==//   cout<<"debug:  ism1_inm2_y0="<<ism1_inm2_y0<<",  ism1_inm2_x0="<<ism1_inm2_x0<<",  ism1_inm2_y_end="<<ism1_inm2_y_end<<",  ism1_inm2_x_end="<<ism1_inm2_x_end<<endl;  //==DEBUG==//  

//int  ism2acrsm1_inm1_y0=ism2_inm1_y0,      ism2acrsm1_inm1_y_end=ism1_inm2_y_end,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
//int  ism2acrsm1_inm1_x0=ism2_inm1_x0,      ism2acrsm1_inm1_x_end=ism1_inm2_x_end,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_y0=0,      ism2acrsm1_inm1_y_end=0,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_x0=0,      ism2acrsm1_inm1_x_end=0,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !

//int  ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//int  ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

int  ism2acrsm1_puttom2_y0=0,  ism2acrsm1_puttom2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_puttom2_x0=0,  ism2acrsm1_puttom2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity

int  ism2acrsm1_fromm2_y0=0,  ism2acrsm1_fromm2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_fromm2_x0=0,  ism2acrsm1_fromm2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity



//  across section handling
//
//  case (1)  *  case (2)  *  case (3)  *  case (4)  *  case (5)  *   here
//            *        |   *            *            *            *
//            *        |   *            *            *            *
//            *            *        |   *            *            *  matr1->max_not_0_y_index
//      |     *      |     *      | |   *      |     *      |     *  |  
//      |     *      |     *      |     *      |     *      | |   *  |
//      |     *      |     *      |     *      |     *      | |   *  |                           ism2acrsm1_inm1_y_end  ***************** ism2acrsm1_inm1_y_end (cross area)
//      |     *      |     *      |     *      | |   *      |     *  |                           |                      ***************** ism2acrsm1_inm1_y0    (cross area)
//            *            *            *        |   *            *  matr1->min_not_0_y_index    |
//        |   *            *            *            *            *                              ism2acrsm1_inm1_y0 
//        |   *            *            *            *            *                           
//is_across_ar*is_across_ar*is_across_ar*is_across_ar*is_across_ar*
// =false     * =false     * =true      * =true      * =true      *

bool is_across_area=true;

  if ((ism2_inm1_size_y > 0) && (ism2_inm1_size_x > 0))        //  !!! ======== v3 M_div_2 line ======== !!!
  {
    if ((is_across_area == true) && (ism2_inm1_y0          <= matr1->max_not_0_y_index)) {ism2acrsm1_inm1_y0=ism2_inm1_y0;       } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_y_end       >= matr1->min_not_0_y_index)) {ism2acrsm1_inm1_y_end=ism2_inm1_y_end; } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_x0          <= matr1->max_not_0_x_index)) {ism2acrsm1_inm1_x0=ism2_inm1_x0;       } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_x_end       >= matr1->min_not_0_x_index)) {ism2acrsm1_inm1_x_end=ism2_inm1_x_end; } else {is_across_area=false;}

    if (is_across_area == true)
    {
      if (ism2acrsm1_inm1_y0    < matr1->min_not_0_y_index) {ism2acrsm1_inm1_y0   =matr1->min_not_0_y_index;}
      if (ism2acrsm1_inm1_y_end > matr1->max_not_0_y_index) {ism2acrsm1_inm1_y_end=matr1->max_not_0_y_index;}  
    ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1;

      if (ism2acrsm1_inm1_x0    < matr1->min_not_0_x_index) {ism2acrsm1_inm1_x0   =matr1->min_not_0_x_index;}
      if (ism2acrsm1_inm1_x_end > matr1->max_not_0_x_index) {ism2acrsm1_inm1_x_end=matr1->max_not_0_x_index;}  
    ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;

    ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_puttom2_y0+ism2acrsm1_inm1_size_y-1;
    ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+shf1_x;    ism2acrsm1_puttom2_x_end=ism2acrsm1_puttom2_x0+ism2acrsm1_inm1_size_x-1;

    ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
    ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
    }
  }
  else
  {is_across_area=false;}



//  _(2)_   
  if ((ism2_inm1_size_y != 0) && (ism2_inm1_size_x != 0))
  {
  matr1_area_spec_add_to_matr2(matr1,  ism2_inm1_y0*matr1->size_x+ism2_inm1_x0,          matr2,  matr2->min_not_0_y_index,  matr2->min_not_0_x_index,  ism2_inm1_size_y,  ism2_inm1_size_x,    ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  }

//  _(4)_
//==DEBUG==//    cout<<"y0_matr1_to="<<ism2acrsm1_puttom2_y0<<", x0_matr1_to="<<ism2acrsm1_puttom2_x0<<endl;
//==DEBUG==//    cout<<"y0_matr2_from="<<ism2acrsm1_fromm2_y0<<", x0_matr2_from="<<ism2acrsm1_fromm2_x0<<endl;
//==DEBUG==//    cout<<"y0_matr2_to="<<ism2acrsm1_inm1_y0<<", x0_matr2_to="<<ism2acrsm1_inm1_x0<<endl;
//==DEBUG==//    cout<<"subarea_len_y="<<ism2acrsm1_inm1_size_y<<"  "<<"subarea_len_x="<<ism2acrsm1_inm1_size_x<<"  "<<endl;
  if ((ism2acrsm1_inm1_size_y !=0) && (ism2acrsm1_inm1_size_x !=0))
  { 
  matr1_area_spec_mutual_add_swap(matr1, ism2acrsm1_puttom2_y0, ism2acrsm1_puttom2_x0,     matr2, ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,                   ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,  ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  } 

//  _(5)_
  if (( ism1_inm2_size_y != 0) && (ism1_inm2_size_x != 0))
  {
  matr1_area_spec_add_to_matr2(matr2,  ism1_inm2_y0*matr2->size_x+ism1_inm2_x0,          matr1,  matr1->min_not_0_y_index,  matr1->min_not_0_x_index,  ism1_inm2_size_y,  ism1_inm2_size_x,    ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  }

//  6:    update not-0-area indicators 

//==//matr2->min_not_0_y_index+=shf2_y;  matr2->max_not_0_y_index+=shf2_y;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)
//==//matr2->min_not_0_x_index+=shf2_x;  matr2->max_not_0_x_index+=shf2_x;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)


  if ((ism2_inm1_size_y != 0) && (ism2_inm1_size_x != 0))
  {
    if ((matr1->not_0_len_y != 0) && (matr1->not_0_len_x != 0))
    {
      if (ism2_inm1_y0                < matr1->min_not_0_y_index ) {matr1->min_not_0_y_index=ism2_inm1_y0;   }
      if (matr1->max_not_0_y_index    < ism2_inm1_y_end          ) {matr1->max_not_0_y_index=ism2_inm1_y_end;}
      if (ism2_inm1_x0                < matr1->min_not_0_x_index ) {matr1->min_not_0_x_index=ism2_inm1_x0;   }
      if (matr1->max_not_0_x_index    < ism2_inm1_x_end          ) {matr1->max_not_0_x_index=ism2_inm1_x_end;}
    }
    else {matr1->min_not_0_y_index=ism2_inm1_y0;  matr1->max_not_0_y_index=ism2_inm1_y_end;     matr1->min_not_0_x_index=ism2_inm1_x0;  matr1->max_not_0_x_index=ism2_inm1_x_end;}

  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1;
  matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
  }


  if ((ism1_inm2_size_y != 0) && (ism1_inm2_size_x != 0))
  {
    if ((matr2->not_0_len_y != 0) && (matr2->not_0_len_x != 0))
    {
      if (ism1_inm2_y0                < matr2->min_not_0_y_index ) {matr2->min_not_0_y_index=ism1_inm2_y0;   }
      if (matr2->max_not_0_y_index    < ism1_inm2_y_end          ) {matr2->max_not_0_y_index=ism1_inm2_y_end;}
      if (ism1_inm2_x0                < matr2->min_not_0_x_index ) {matr2->min_not_0_x_index=ism1_inm2_x0;   }
      if (matr2->max_not_0_x_index    < ism1_inm2_x_end          ) {matr2->max_not_0_x_index=ism1_inm2_x_end;}
    }
    else
    {matr2->min_not_0_y_index=ism1_inm2_y0;  matr2->max_not_0_y_index=ism1_inm2_y_end;    matr2->min_not_0_x_index=ism1_inm2_x0;  matr2->max_not_0_x_index=ism1_inm2_x_end;}

  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1;
  matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;   
  }


//==DEBUG==//  cout<<"DEBUG  became:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  became:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  


#ifdef __MACR_CHECK_ACCESS_MEM_VIOLATION

  if ((matr1->min_not_0_y_index < 0) || (matr1->min_not_0_y_index >= matr1->size_y) || (matr1->min_not_0_x_index < 0) || (matr1->min_not_0_x_index >= matr1->size_x) ||
      (matr1->max_not_0_y_index < 0) || (matr1->max_not_0_y_index >= matr1->size_y) || (matr1->max_not_0_x_index < 0) || (matr1->max_not_0_x_index >= matr1->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr1"<<endl;  
  cout<<"matr1:"<<endl;
  matr1->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr1->min_not_0_y_index=0;  matr1->max_not_0_y_index=matr1->size_y-1;  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1; 
  matr1->min_not_0_x_index=0;  matr1->max_not_0_x_index=matr1->size_x-1;  matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
  //return false; 
  }

  if ((matr2->min_not_0_y_index < 0) || (matr2->min_not_0_y_index >= matr2->size_y) || (matr2->min_not_0_x_index < 0) || (matr2->min_not_0_x_index >= matr2->size_x) ||
      (matr2->max_not_0_y_index < 0) || (matr2->max_not_0_y_index >= matr2->size_y) || (matr2->max_not_0_x_index < 0) || (matr2->max_not_0_x_index >= matr2->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr2"<<endl;  
  cout<<"matr2:"<<endl;
  matr2->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr2->min_not_0_y_index=0;  matr2->max_not_0_y_index=matr2->size_y-1;  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1; 
  matr2->min_not_0_x_index=0;  matr2->max_not_0_x_index=matr2->size_x-1;  matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;
  //return false; 
  }

#endif


#endif


//==//ism2_inm1_y0=matr1->min_not_0_y_index-shf2_y;  ism2_inm1_y_end=matr2->min_not_0_y_index-shf2_y;   ism2_inm1_size_y=matr2->not_0_len_y;
//==//ism2_inm1_x0=matr1->min_not_0_x_index-shf2_x;  ism2_inm1_x_end=matr2->min_not_0_x_index-shf2_x;   ism2_inm1_size_x=matr2->not_0_len_x;

//==//ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y;  ism1_inm2_y_end=matr1->max_not_0_y_index+shf1_y;   ism1_inm2_size_y=matr1->not_0_len_y;
//==//ism1_inm2_x0=matr1->min_not_0_x_index+shf1_x;  ism1_inm2_x_end=matr1->max_not_0_x_index+shf1_x;   ism1_inm2_size_x=matr1->not_0_len_x;


//==//ism2acrsm1_inm1_y0=;  ism2acrsm1_inm1_y_end=;       ism2acrsm1_inm1_size_y=;
//==//ism2acrsm1_inm1_x0=;  ism2acrsm1_inm1_x_end=;       ism2acrsm1_inm1_size_x=;

//==//ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//==//ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

//==//ism2acrsm1_inm2_y0=;  ism2acrsm1_inm2_y_end=;      // ism2acrsm1_inm2_size_y=;
//==//ism2acrsm1_inm2_x0=;  ism2acrsm1_inm2_x_end=;      // ism2acrsm1_inm2_size_x=;


}




//  --------------------------------------------------------------------------------------------------------------------------------------
//  ======================================================================================================================================
//  --------------------------------------------------------------------------------------------------------------------------------------



bool matr_cross_add_v5(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y,  int shf2_y)      // ! ! ! ======= very importand and difficult proc  ======= ! ! !
{
//==DEBUG==//  cout<<"DEBUG:   matr_cross_add_v3"<<endl;  //==DEBUG==//  

#ifdef __MACR_MATR_TYPE_1__MPZ

int  shf1_x=1, shf2_x=1;

//  _(3)_  
//matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x);
//==v1==//matr2->shift_elements(shf2_y,  shf2_x);
//==DEBUG==//  cout<<"DEBUG:   "<<"  shf2_y="<<shf2_y<<endl;  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  WAS:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  WAS:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
  if ((matr2->not_0_len_y != 0) && (matr2->not_0_len_x != 0)) {matr2->shift_elements_spec_fast_case_M_div_2_v20(shf2_y);}      //  !!! ======== v3 M_div_2 line modif ======== !!!
//==DEBUG==//  cout<<"DEBUG  INTERMED:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  INTERMED:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  

     

//  1
//  define main parameters

//==DEBUG==//  cout<<"debug:  matr2->min_not_0_y_index="<<matr2->min_not_0_y_index<<" matr2->max_not_0_y_index="<<matr2->max_not_0_y_index<<" matr2->not_0_len_y="<<matr2->not_0_len_y<<endl;  //==DEBUG==//
//==DEBUG==//  cout<<"debug:  matr2->min_not_0_x_index="<<matr2->min_not_0_x_index<<" matr2->max_not_0_x_index="<<matr2->max_not_0_x_index<<" matr2->not_0_len_x="<<matr2->not_0_len_x<<endl;  //==DEBUG==//
//cout<<"debug:  matr2->min_not_0_y_index" matr2->not_0_len_x="<<matr2->not_0_len_x<<endl;  //==DEBUG==//   
int  ism2_inm1_y0=matr2->min_not_0_y_index-shf2_y,  ism2_inm1_y_end=ism2_inm1_y0+matr2->not_0_len_y-1,   ism2_inm1_size_y=ism2_inm1_y_end-ism2_inm1_y0+1;   //matr2->not_0_len_y;
int  ism2_inm1_x0=matr2->min_not_0_x_index-shf2_x,  ism2_inm1_x_end=ism2_inm1_x0+matr2->not_0_len_x-1,   ism2_inm1_size_x=ism2_inm1_x_end-ism2_inm1_x0+1;   //matr2->not_0_len_x; 
  if ((matr2->not_0_len_y == 0) || (matr2->not_0_len_x == 0) || (ism2_inm1_y_end < 0) || (ism2_inm1_x_end < 0) || (ism2_inm1_y0 >= matr2->size_y) || (ism2_inm1_x0 >= matr2->size_x)) 
  {
  ism2_inm1_y0=0;  ism2_inm1_y_end=0;   ism2_inm1_size_y=0;   //matr2->not_0_len_y;
  ism2_inm1_x0=0;  ism2_inm1_x_end=0;   ism2_inm1_size_x=0;   //matr2->not_0_len_x; 
  }  
  else
  {
    if (ism2_inm1_y0 < 0) {ism2_inm1_y0=0;}    if (ism2_inm1_y0 >= matr2->size_y) {ism2_inm1_y0=matr2->size_y-1;}  ism2_inm1_size_y=ism2_inm1_y_end-ism2_inm1_y0+1;
    if (ism2_inm1_x0 < 0) {ism2_inm1_x0=0;}    if (ism2_inm1_x0 >= matr2->size_x) {ism2_inm1_x0=matr2->size_x-1;}  ism2_inm1_size_x=ism2_inm1_x_end-ism2_inm1_x0+1;
  }
//==DEBUG==//  cout<<"debug:  ism2_inm1_y0="<<ism2_inm1_y0<<",  ism2_inm1_x0="<<ism2_inm1_x0<<",  ism2_inm1_y_end="<<ism2_inm1_y_end<<",  ism2_inm1_x_end="<<ism2_inm1_x_end<<endl;  //==DEBUG==//  



int  ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y,  ism1_inm2_y_end=ism1_inm2_y0+matr1->not_0_len_y-1,   ism1_inm2_size_y=ism1_inm2_y_end-ism1_inm2_y0+1;      //matr1->not_0_len_y;
int  ism1_inm2_x0=matr1->min_not_0_x_index+shf1_x,  ism1_inm2_x_end=ism1_inm2_x0+matr1->not_0_len_x-1,   ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;      //matr1->not_0_len_x;
  if ((matr1->not_0_len_y == 0) || (matr1->not_0_len_x == 0) || (ism1_inm2_y_end < 0) || (ism1_inm2_x_end < 0) || (ism1_inm2_y0 >= matr2->size_y) || (ism1_inm2_x0 >= matr1->size_x))
  {
  ism1_inm2_y0=0;  ism1_inm2_y_end=0;   ism1_inm2_size_y=0;   //matr2->not_0_len_y;
  ism1_inm2_x0=0;  ism1_inm2_x_end=0;   ism1_inm2_size_x=0;   //matr2->not_0_len_x; 
  } 
  else
  {
    if (ism1_inm2_y0 < 0) {ism1_inm2_y0=0;}    if (ism1_inm2_y0 >= matr1->size_y) {ism1_inm2_y0=matr1->size_y-1;}  ism1_inm2_size_y=ism1_inm2_y_end-ism1_inm2_y0+1;
    if (ism1_inm2_x0 < 0) {ism1_inm2_x0=0;}    if (ism1_inm2_x0 >= matr1->size_x) {ism1_inm2_x0=matr1->size_x-1;}  ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;
  }

  if (ism1_inm2_x_end >= matr2->size_x) {ism1_inm2_x_end=matr2->size_x-1; ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;}       //  !!! ======== v3 M_div_2 line ======== !!!
  if (ism1_inm2_size_x < 1) {ism1_inm2_size_x=0;  ism1_inm2_x0=0;  ism1_inm2_x_end=0;}                                            //  !!! ======== v3 M_div_2 line ======== !!!

int m1_area_will_be_car_to_m2_y0=ism1_inm2_y0-shf1_y,  m1_area_will_be_car_to_m2_y_end=ism1_inm2_y_end-shf1_y;      //  taking into account getting out of matrix 1 within m2 because of negative shf_y
int m1_area_will_be_car_to_m2_x0=ism1_inm2_x0-1,       m1_area_will_be_car_to_m2_x_end=ism1_inm2_x_end-1;           //  taking into account getting out of matrix 1 within m2 because of negative shf_y
//==DEBUG==//   cout<<"debug:  m1_area_will_be_car_to_m2_y0="<<m1_area_will_be_car_to_m2_y0<<" m1_area_will_be_car_to_m2_y_end="<<m1_area_will_be_car_to_m2_y_end<<endl;  //==DEBUG==//  
//==DEBUG==//   cout<<"debug:  m1_area_will_be_car_to_m2_x0="<<m1_area_will_be_car_to_m2_x0<<" m1_area_will_be_car_to_m2_x_end="<<m1_area_will_be_car_to_m2_x_end<<endl;  //==DEBUG==//  



//==DEBUG==//  cout<<"debug:  ism1_inm2_y0="<<ism1_inm2_y0<<",  ism1_inm2_x0="<<ism1_inm2_x0<<",  ism1_inm2_y_end="<<ism1_inm2_y_end<<",  ism1_inm2_x_end="<<ism1_inm2_x_end<<endl;  //==DEBUG==//  

//int  ism2acrsm1_inm1_y0=ism2_inm1_y0,      ism2acrsm1_inm1_y_end=ism1_inm2_y_end,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
//int  ism2acrsm1_inm1_x0=ism2_inm1_x0,      ism2acrsm1_inm1_x_end=ism1_inm2_x_end,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_y0=0,      ism2acrsm1_inm1_y_end=0,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_x0=0,      ism2acrsm1_inm1_x_end=0,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !

//int  ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//int  ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

int  ism2acrsm1_puttom2_y0=0,  ism2acrsm1_puttom2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_puttom2_x0=0,  ism2acrsm1_puttom2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity

int  ism2acrsm1_fromm2_y0=0,  ism2acrsm1_fromm2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_fromm2_x0=0,  ism2acrsm1_fromm2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity



//  across section handling
//
//  case (1)  *  case (2)  *  case (3)  *  case (4)  *  case (5)  *   here
//            *        |   *            *            *            *
//            *        |   *            *            *            *
//            *            *        |   *            *            *  matr1->max_not_0_y_index
//      |     *      |     *      | |   *      |     *      |     *  |  
//      |     *      |     *      |     *      |     *      | |   *  |
//      |     *      |     *      |     *      |     *      | |   *  |                           ism2acrsm1_inm1_y_end  ***************** ism2acrsm1_inm1_y_end (cross area)
//      |     *      |     *      |     *      | |   *      |     *  |                           |                      ***************** ism2acrsm1_inm1_y0    (cross area)
//            *            *            *        |   *            *  matr1->min_not_0_y_index    |
//        |   *            *            *            *            *                              ism2acrsm1_inm1_y0 
//        |   *            *            *            *            *                           
//is_across_ar*is_across_ar*is_across_ar*is_across_ar*is_across_ar*
// =false     * =false     * =true      * =true      * =true      *

bool is_across_area=true;

  if ((ism2_inm1_size_y > 0) && (ism2_inm1_size_x > 0))        //  !!! ======== v3 M_div_2 line ======== !!!
  {
    if ((is_across_area == true) && (ism2_inm1_y0          <= m1_area_will_be_car_to_m2_y_end)) {ism2acrsm1_inm1_y0   =ism2_inm1_y0;    } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_y_end       >= m1_area_will_be_car_to_m2_y0   )) {ism2acrsm1_inm1_y_end=ism2_inm1_y_end; } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_x0          <= m1_area_will_be_car_to_m2_x_end)) {ism2acrsm1_inm1_x0   =ism2_inm1_x0;    } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_x_end       >= m1_area_will_be_car_to_m2_x0   )) {ism2acrsm1_inm1_x_end=ism2_inm1_x_end; } else {is_across_area=false;}

    if (is_across_area == true)      //  further limitation
    {
      if (ism2acrsm1_inm1_y0    < m1_area_will_be_car_to_m2_y0   ) {ism2acrsm1_inm1_y0   =m1_area_will_be_car_to_m2_y0;}
      if (ism2acrsm1_inm1_y_end > m1_area_will_be_car_to_m2_y_end) {ism2acrsm1_inm1_y_end=m1_area_will_be_car_to_m2_y_end;}  
    ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1;

      if (ism2acrsm1_inm1_x0    < m1_area_will_be_car_to_m2_x0   ) {ism2acrsm1_inm1_x0   =m1_area_will_be_car_to_m2_x0;}
      if (ism2acrsm1_inm1_x_end > m1_area_will_be_car_to_m2_x_end) {ism2acrsm1_inm1_x_end=m1_area_will_be_car_to_m2_x_end;}
    ism2acrsm1_inm1_size_x=ism2acrsm1_inm1_x_end-ism2acrsm1_inm1_x0+1;

    ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_inm1_y_end+shf1_y;
    ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0+shf1_x;    ism2acrsm1_puttom2_x_end=ism2acrsm1_inm1_x_end+shf1_x;

    //ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
    //ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
    ism2acrsm1_fromm2_y0=ism2acrsm1_inm1_y0+shf2_y;    ism2acrsm1_fromm2_y_end=ism2acrsm1_inm1_y_end+shf2_y;
    ism2acrsm1_fromm2_x0=ism2acrsm1_inm1_x0+shf2_x;    ism2acrsm1_fromm2_x_end=ism2acrsm1_inm1_x_end+shf2_x;
    }
  }
  else
  {is_across_area=false;}


//  _(2)_   
  if ((ism2_inm1_size_y != 0) && (ism2_inm1_size_x != 0))
  {
  matr1_area_spec_add_to_matr2(matr1,  ism2_inm1_y0*matr1->size_x+ism2_inm1_x0,          matr2,  matr2->min_not_0_y_index,  matr2->min_not_0_x_index,  ism2_inm1_size_y,  ism2_inm1_size_x,    ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  }

//  _(4)_
  //==DEBUG==//  cout<<"y0_matr1_to="<<ism2acrsm1_puttom2_y0<<", x0_matr1_to="<<ism2acrsm1_puttom2_x0<<endl;                  //==DEBUG==//
  //==DEBUG==//  cout<<"y0_matr2_from="<<ism2acrsm1_fromm2_y0<<", x0_matr2_from="<<ism2acrsm1_fromm2_x0<<endl;                //==DEBUG==//
  //==DEBUG==//  cout<<"y0_matr2_to="<<ism2acrsm1_inm1_y0<<", x0_matr2_to="<<ism2acrsm1_inm1_x0<<endl;                        //==DEBUG==//
  //==DEBUG==//  cout<<"subarea_len_y="<<ism2acrsm1_inm1_size_y<<"  "<<"subarea_len_x="<<ism2acrsm1_inm1_size_x<<"  "<<endl;  //==DEBUG==//
  if ((ism2acrsm1_inm1_size_y !=0) && (ism2acrsm1_inm1_size_x !=0))
  { 
  matr1_area_spec_mutual_add_swap(matr1, ism2acrsm1_puttom2_y0, ism2acrsm1_puttom2_x0,     matr2, ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,                   ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,  ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  } 

//  _(5)_
  if (( ism1_inm2_size_y != 0) && (ism1_inm2_size_x != 0))
  {
  matr1_area_spec_add_to_matr2(matr2,  ism1_inm2_y0*matr2->size_x+ism1_inm2_x0,          matr1,  m1_area_will_be_car_to_m2_y0,  m1_area_will_be_car_to_m2_x0,  ism1_inm2_size_y,  ism1_inm2_size_x,    ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  }

//  6:    update not-0-area indicators 

//==//matr2->min_not_0_y_index+=shf2_y;  matr2->max_not_0_y_index+=shf2_y;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)
//==//matr2->min_not_0_x_index+=shf2_x;  matr2->max_not_0_x_index+=shf2_x;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)


  if ((ism2_inm1_size_y != 0) && (ism2_inm1_size_x != 0))
  {
    if ((matr1->not_0_len_y != 0) && (matr1->not_0_len_x != 0))
    {
      if (ism2_inm1_y0                < matr1->min_not_0_y_index ) {matr1->min_not_0_y_index=ism2_inm1_y0;   }
      if (matr1->max_not_0_y_index    < ism2_inm1_y_end          ) {matr1->max_not_0_y_index=ism2_inm1_y_end;}
      if (ism2_inm1_x0                < matr1->min_not_0_x_index ) {matr1->min_not_0_x_index=ism2_inm1_x0;   }
      if (matr1->max_not_0_x_index    < ism2_inm1_x_end          ) {matr1->max_not_0_x_index=ism2_inm1_x_end;}
    }
    else {matr1->min_not_0_y_index=ism2_inm1_y0;  matr1->max_not_0_y_index=ism2_inm1_y_end;     matr1->min_not_0_x_index=ism2_inm1_x0;  matr1->max_not_0_x_index=ism2_inm1_x_end;}

  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1;
  matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
  }


  if ((ism1_inm2_size_y != 0) && (ism1_inm2_size_x != 0))
  {
    if ((matr2->not_0_len_y != 0) && (matr2->not_0_len_x != 0))
    {
      if (ism1_inm2_y0                < matr2->min_not_0_y_index ) {matr2->min_not_0_y_index=ism1_inm2_y0;   }
      if (matr2->max_not_0_y_index    < ism1_inm2_y_end          ) {matr2->max_not_0_y_index=ism1_inm2_y_end;}
      if (ism1_inm2_x0                < matr2->min_not_0_x_index ) {matr2->min_not_0_x_index=ism1_inm2_x0;   }
      if (matr2->max_not_0_x_index    < ism1_inm2_x_end          ) {matr2->max_not_0_x_index=ism1_inm2_x_end;}
    }
    else
    {matr2->min_not_0_y_index=ism1_inm2_y0;  matr2->max_not_0_y_index=ism1_inm2_y_end;    matr2->min_not_0_x_index=ism1_inm2_x0;  matr2->max_not_0_x_index=ism1_inm2_x_end;}

  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1;
  matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;   
  }


//==DEBUG==//  cout<<"DEBUG  became:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  became:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  


#ifdef __MACR_CHECK_ACCESS_MEM_VIOLATION

  if ((matr1->min_not_0_y_index < 0) || (matr1->min_not_0_y_index >= matr1->size_y) || (matr1->min_not_0_x_index < 0) || (matr1->min_not_0_x_index >= matr1->size_x) ||
      (matr1->max_not_0_y_index < 0) || (matr1->max_not_0_y_index >= matr1->size_y) || (matr1->max_not_0_x_index < 0) || (matr1->max_not_0_x_index >= matr1->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr1"<<endl;  
  cout<<"matr1:"<<endl;
  matr1->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr1->min_not_0_y_index=0;  matr1->max_not_0_y_index=matr1->size_y-1;  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1; 
  matr1->min_not_0_x_index=0;  matr1->max_not_0_x_index=matr1->size_x-1;  matr1->not_0_len_x=matr1->max_not_0_x_index-matr1->min_not_0_x_index+1;
  //return false; 
  }

  if ((matr2->min_not_0_y_index < 0) || (matr2->min_not_0_y_index >= matr2->size_y) || (matr2->min_not_0_x_index < 0) || (matr2->min_not_0_x_index >= matr2->size_x) ||
      (matr2->max_not_0_y_index < 0) || (matr2->max_not_0_y_index >= matr2->size_y) || (matr2->max_not_0_x_index < 0) || (matr2->max_not_0_x_index >= matr2->size_x)    )
  { 
  //  error 
  display_on_screen_error_msg((char *) "in proc    matr_cross_add_v1:  mem_access_viol");
  cout<<"error details start:"<<endl; 
  cout<<"not 0-area is out of accessible area of matrix matr2"<<endl;  
  cout<<"matr2:"<<endl;
  matr2->display_on_screen_not_0_area_size_param();
  cout<<"error details end"<<endl; 
  matr2->min_not_0_y_index=0;  matr2->max_not_0_y_index=matr2->size_y-1;  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1; 
  matr2->min_not_0_x_index=0;  matr2->max_not_0_x_index=matr2->size_x-1;  matr2->not_0_len_x=matr2->max_not_0_x_index-matr2->min_not_0_x_index+1;
  //return false; 
  }

#endif


#endif


//==//ism2_inm1_y0=matr1->min_not_0_y_index-shf2_y;  ism2_inm1_y_end=matr2->min_not_0_y_index-shf2_y;   ism2_inm1_size_y=matr2->not_0_len_y;
//==//ism2_inm1_x0=matr1->min_not_0_x_index-shf2_x;  ism2_inm1_x_end=matr2->min_not_0_x_index-shf2_x;   ism2_inm1_size_x=matr2->not_0_len_x;

//==//ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y;  ism1_inm2_y_end=matr1->max_not_0_y_index+shf1_y;   ism1_inm2_size_y=matr1->not_0_len_y;
//==//ism1_inm2_x0=matr1->min_not_0_x_index+shf1_x;  ism1_inm2_x_end=matr1->max_not_0_x_index+shf1_x;   ism1_inm2_size_x=matr1->not_0_len_x;


//==//ism2acrsm1_inm1_y0=;  ism2acrsm1_inm1_y_end=;       ism2acrsm1_inm1_size_y=;
//==//ism2acrsm1_inm1_x0=;  ism2acrsm1_inm1_x_end=;       ism2acrsm1_inm1_size_x=;

//==//ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//==//ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

//==//ism2acrsm1_inm2_y0=;  ism2acrsm1_inm2_y_end=;      // ism2acrsm1_inm2_size_y=;
//==//ism2acrsm1_inm2_x0=;  ism2acrsm1_inm2_x_end=;      // ism2acrsm1_inm2_size_x=;

return true;

}




//  --------------------------------------------------------------------------------------------------------------------------------------
//  ======================================================================================================================================
//  --------------------------------------------------------------------------------------------------------------------------------------



bool matr_cross_add_v6(matr_t1 *matr1,  matr_t1 *matr2,  int shf1_y,  int shf2_y)      // ! ! ! ======= very importand and difficult proc  ======= ! ! !
{
//==DEBUG==//  cout<<"DEBUG:   matr_cross_add_v3"<<endl;  //==DEBUG==//  


int  shf1_x=0, shf2_x=0;

//  _(3)_  
//matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x);
//==v1==//matr2->shift_elements(shf2_y,  shf2_x);
//==DEBUG==//  cout<<"DEBUG:   "<<"  shf2_y="<<shf2_y<<endl;  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  WAS:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  WAS:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
  if ((matr2->not_0_len_y != 0) && (matr2->not_0_len_x != 0)) {matr2->shift_elements_spec_without_M_v22(shf2_y);}      //  !!! ======== v3 M_div_2 line modif ======== !!!
//==DEBUG==//  cout<<"DEBUG  INTERMED:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  INTERMED:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  

     

//  1
//  define main parameters

//==DEBUG==//  cout<<"debug:  matr2->min_not_0_y_index="<<matr2->min_not_0_y_index<<" matr2->max_not_0_y_index="<<matr2->max_not_0_y_index<<" matr2->not_0_len_y="<<matr2->not_0_len_y<<endl;  //==DEBUG==//
//==DEBUG==//  cout<<"debug:  matr2->min_not_0_x_index="<<matr2->min_not_0_x_index<<" matr2->max_not_0_x_index="<<matr2->max_not_0_x_index<<" matr2->not_0_len_x="<<matr2->not_0_len_x<<endl;  //==DEBUG==//
//cout<<"debug:  matr2->min_not_0_y_index" matr2->not_0_len_x="<<matr2->not_0_len_x<<endl;  //==DEBUG==//   
int  ism2_inm1_y0=matr2->min_not_0_y_index-shf2_y,  ism2_inm1_y_end=ism2_inm1_y0+matr2->not_0_len_y-1,   ism2_inm1_size_y=ism2_inm1_y_end-ism2_inm1_y0+1;   //matr2->not_0_len_y;
int  ism2_inm1_x0=0,  ism2_inm1_x_end=0,   ism2_inm1_size_x=1;   //matr2->not_0_len_x; 
  if ((matr2->not_0_len_y == 0) || (matr2->not_0_len_x == 0) || (ism2_inm1_y_end < 0) || (ism2_inm1_x_end < 0) || (ism2_inm1_y0 >= matr2->size_y) || (ism2_inm1_x0 >= matr2->size_x)) 
  {
  ism2_inm1_y0=0;  ism2_inm1_y_end=0;   ism2_inm1_size_y=0;   //matr2->not_0_len_y;
  ism2_inm1_x0=0;  ism2_inm1_x_end=0;   ism2_inm1_size_x=0;   //matr2->not_0_len_x; 
  }  
  else
  {
    if (ism2_inm1_y0 < 0) {ism2_inm1_y0=0;}    if (ism2_inm1_y0 >= matr2->size_y) {ism2_inm1_y0=matr2->size_y-1;}  ism2_inm1_size_y=ism2_inm1_y_end-ism2_inm1_y0+1;
  }
ism2_inm1_x0=0;  ism2_inm1_x_end=0;  ism2_inm1_size_x=1;
//==DEBUG==//  cout<<"debug:  ism2_inm1_y0="<<ism2_inm1_y0<<",  ism2_inm1_x0="<<ism2_inm1_x0<<",  ism2_inm1_y_end="<<ism2_inm1_y_end<<",  ism2_inm1_x_end="<<ism2_inm1_x_end<<endl;  //==DEBUG==//  



int  ism1_inm2_y0=matr1->min_not_0_y_index+shf1_y,  ism1_inm2_y_end=ism1_inm2_y0+matr1->not_0_len_y-1,   ism1_inm2_size_y=ism1_inm2_y_end-ism1_inm2_y0+1;      //matr1->not_0_len_y;
int  ism1_inm2_x0=0,  ism1_inm2_x_end=0,   ism1_inm2_size_x=1;      //matr1->not_0_len_x;
  if ((matr1->not_0_len_y == 0) || (matr1->not_0_len_x == 0) || (ism1_inm2_y_end < 0) || (ism1_inm2_x_end < 0) || (ism1_inm2_y0 >= matr2->size_y) || (ism1_inm2_x0 >= matr1->size_x))
  {
  ism1_inm2_y0=0;  ism1_inm2_y_end=0;   ism1_inm2_size_y=0;   //matr2->not_0_len_y;
  ism1_inm2_x0=0;  ism1_inm2_x_end=0;   ism1_inm2_size_x=0;   //matr2->not_0_len_x; 
  } 
  else
  {
    if (ism1_inm2_y0 < 0) {ism1_inm2_y0=0;}    if (ism1_inm2_y0 >= matr1->size_y) {ism1_inm2_y0=matr1->size_y-1;}  ism1_inm2_size_y=ism1_inm2_y_end-ism1_inm2_y0+1;
    if (ism1_inm2_x0 < 0) {ism1_inm2_x0=0;}    if (ism1_inm2_x0 >= matr1->size_x) {ism1_inm2_x0=matr1->size_x-1;}  ism1_inm2_size_x=ism1_inm2_x_end-ism1_inm2_x0+1;
  }

  if (ism1_inm2_size_x < 1) {ism1_inm2_size_x=0;  ism1_inm2_x0=0;  ism1_inm2_x_end=0;}   
ism1_inm2_x0=0;  ism1_inm2_x_end=0;  ism1_inm2_size_x=1;                                         //  !!! ======== v3 M_div_2 line ======== !!!

int m1_area_will_be_car_to_m2_y0=ism1_inm2_y0-shf1_y,  m1_area_will_be_car_to_m2_y_end=ism1_inm2_y_end-shf1_y;      //  taking into account getting out of matrix 1 within m2 because of negative shf_y
int m1_area_will_be_car_to_m2_x0=0,                    m1_area_will_be_car_to_m2_x_end=0;           //  taking into account getting out of matrix 1 within m2 because of negative shf_y
//==DEBUG==//   cout<<"debug:  m1_area_will_be_car_to_m2_y0="<<m1_area_will_be_car_to_m2_y0<<" m1_area_will_be_car_to_m2_y_end="<<m1_area_will_be_car_to_m2_y_end<<endl;  //==DEBUG==//  
//==DEBUG==//   cout<<"debug:  m1_area_will_be_car_to_m2_x0="<<m1_area_will_be_car_to_m2_x0<<" m1_area_will_be_car_to_m2_x_end="<<m1_area_will_be_car_to_m2_x_end<<endl;  //==DEBUG==//  



//==DEBUG==//  cout<<"debug:  ism1_inm2_y0="<<ism1_inm2_y0<<",  ism1_inm2_x0="<<ism1_inm2_x0<<",  ism1_inm2_y_end="<<ism1_inm2_y_end<<",  ism1_inm2_x_end="<<ism1_inm2_x_end<<endl;  //==DEBUG==//  

//int  ism2acrsm1_inm1_y0=ism2_inm1_y0,      ism2acrsm1_inm1_y_end=ism1_inm2_y_end,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
//int  ism2acrsm1_inm1_x0=ism2_inm1_x0,      ism2acrsm1_inm1_x_end=ism1_inm2_x_end,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_y0=0,      ism2acrsm1_inm1_y_end=0,       ism2acrsm1_inm1_size_y=0;      //  ! ! ! the most important area ! ! !
int  ism2acrsm1_inm1_x0=0,      ism2acrsm1_inm1_x_end=0,       ism2acrsm1_inm1_size_x=0;      //  ! ! ! the most important area ! ! !

//int  ism1acrsm2_inm2_y0=;  ism1acrsm2_inm2_y_end=;       ism1acrsm2_inm2_size_y=;
//int  ism1acrsm2_inm2_x0=;  ism1acrsm2_inm2_x_end=;       ism1acrsm2_inm2_size_x=;

int  ism2acrsm1_puttom2_y0=0,  ism2acrsm1_puttom2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_puttom2_x0=0,  ism2acrsm1_puttom2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity

int  ism2acrsm1_fromm2_y0=0,  ism2acrsm1_fromm2_y_end=0;      // ism2acrsm1_inm2_size_y=;      //  there is not necessity
int  ism2acrsm1_fromm2_x0=0,  ism2acrsm1_fromm2_x_end=0;      // ism2acrsm1_inm2_size_x=;      //  there is not necessity



//  across section handling
//
//  case (1)  *  case (2)  *  case (3)  *  case (4)  *  case (5)  *   here
//            *        |   *            *            *            *
//            *        |   *            *            *            *
//            *            *        |   *            *            *  matr1->max_not_0_y_index
//      |     *      |     *      | |   *      |     *      |     *  |  
//      |     *      |     *      |     *      |     *      | |   *  |
//      |     *      |     *      |     *      |     *      | |   *  |                           ism2acrsm1_inm1_y_end  ***************** ism2acrsm1_inm1_y_end (cross area)
//      |     *      |     *      |     *      | |   *      |     *  |                           |                      ***************** ism2acrsm1_inm1_y0    (cross area)
//            *            *            *        |   *            *  matr1->min_not_0_y_index    |
//        |   *            *            *            *            *                              ism2acrsm1_inm1_y0 
//        |   *            *            *            *            *                           
//is_across_ar*is_across_ar*is_across_ar*is_across_ar*is_across_ar*
// =false     * =false     * =true      * =true      * =true      *

bool is_across_area=true;

  if ((ism2_inm1_size_y > 0) && (ism2_inm1_size_x > 0))        //  !!! ======== v3 M_div_2 line ======== !!!
  {
    if ((is_across_area == true) && (ism2_inm1_y0          <= m1_area_will_be_car_to_m2_y_end)) {ism2acrsm1_inm1_y0   =ism2_inm1_y0;    } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_y_end       >= m1_area_will_be_car_to_m2_y0   )) {ism2acrsm1_inm1_y_end=ism2_inm1_y_end; } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_x0          <= m1_area_will_be_car_to_m2_x_end)) {ism2acrsm1_inm1_x0   =ism2_inm1_x0;    } else {is_across_area=false;}
    if ((is_across_area == true) && (ism2_inm1_x_end       >= m1_area_will_be_car_to_m2_x0   )) {ism2acrsm1_inm1_x_end=ism2_inm1_x_end; } else {is_across_area=false;}

    if (is_across_area == true)      //  further limitation
    {
      if (ism2acrsm1_inm1_y0    < m1_area_will_be_car_to_m2_y0   ) {ism2acrsm1_inm1_y0   =m1_area_will_be_car_to_m2_y0;}
      if (ism2acrsm1_inm1_y_end > m1_area_will_be_car_to_m2_y_end) {ism2acrsm1_inm1_y_end=m1_area_will_be_car_to_m2_y_end;}  
    ism2acrsm1_inm1_size_y=ism2acrsm1_inm1_y_end-ism2acrsm1_inm1_y0+1;

      if (ism2acrsm1_inm1_x0    < m1_area_will_be_car_to_m2_x0   ) {ism2acrsm1_inm1_x0   =m1_area_will_be_car_to_m2_x0;}
      if (ism2acrsm1_inm1_x_end > m1_area_will_be_car_to_m2_x_end) {ism2acrsm1_inm1_x_end=m1_area_will_be_car_to_m2_x_end;}
    ism2acrsm1_inm1_x0=0;  ism2acrsm1_inm1_x_end=0;
    ism2acrsm1_inm1_size_x=1;

    ism2acrsm1_puttom2_y0=ism2acrsm1_inm1_y0+shf1_y;    ism2acrsm1_puttom2_y_end=ism2acrsm1_inm1_y_end+shf1_y;
    ism2acrsm1_puttom2_x0=ism2acrsm1_inm1_x0;           ism2acrsm1_puttom2_x_end=ism2acrsm1_inm1_x_end;

    //ism2acrsm1_fromm2_y0=matr2->min_not_0_y_index+ism2acrsm1_inm1_y0-ism2_inm1_y0;    ism2acrsm1_fromm2_y_end=ism2acrsm1_fromm2_y0+ism2acrsm1_inm1_size_y-1;
    //ism2acrsm1_fromm2_x0=matr2->min_not_0_x_index+ism2acrsm1_inm1_x0-ism2_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_fromm2_x0+ism2acrsm1_inm1_size_x-1;
    ism2acrsm1_fromm2_y0=ism2acrsm1_inm1_y0+shf2_y;    ism2acrsm1_fromm2_y_end=ism2acrsm1_inm1_y_end+shf2_y;
    ism2acrsm1_fromm2_x0=ism2acrsm1_inm1_x0;    ism2acrsm1_fromm2_x_end=ism2acrsm1_inm1_x_end;
    }
  }
  else
  {is_across_area=false;}


//  _(2)_   
  if ((ism2_inm1_size_y != 0) && (ism2_inm1_size_x != 0))
  {
  matr1_area_spec_add_to_matr2(matr1,  ism2_inm1_y0,          matr2,  matr2->min_not_0_y_index,  0,  ism2_inm1_size_y,  ism2_inm1_size_x,    ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  }

//  _(4)_
  //==DEBUG==//  cout<<"y0_matr1_to="<<ism2acrsm1_puttom2_y0<<", x0_matr1_to="<<ism2acrsm1_puttom2_x0<<endl;                  //==DEBUG==//
  //==DEBUG==//  cout<<"y0_matr2_from="<<ism2acrsm1_fromm2_y0<<", x0_matr2_from="<<ism2acrsm1_fromm2_x0<<endl;                //==DEBUG==//
  //==DEBUG==//  cout<<"y0_matr2_to="<<ism2acrsm1_inm1_y0<<", x0_matr2_to="<<ism2acrsm1_inm1_x0<<endl;                        //==DEBUG==//
  //==DEBUG==//  cout<<"subarea_len_y="<<ism2acrsm1_inm1_size_y<<"  "<<"subarea_len_x="<<ism2acrsm1_inm1_size_x<<"  "<<endl;  //==DEBUG==//
  if ((ism2acrsm1_inm1_size_y !=0) && (ism2acrsm1_inm1_size_x !=0))
  { 
  matr1_area_spec_mutual_add_swap(matr1, ism2acrsm1_puttom2_y0, ism2acrsm1_puttom2_x0,     matr2, ism2acrsm1_fromm2_y0, ism2acrsm1_fromm2_x0,                   ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,  ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  } 

//  _(5)_
  if (( ism1_inm2_size_y != 0) && (ism1_inm2_size_x != 0))
  {
  matr1_area_spec_add_to_matr2(matr2,  ism1_inm2_y0*matr2->size_x+ism1_inm2_x0,          matr1,  m1_area_will_be_car_to_m2_y0,  m1_area_will_be_car_to_m2_x0,  ism1_inm2_size_y,  ism1_inm2_size_x,    ism2acrsm1_inm1_y0, ism2acrsm1_inm1_x0,    ism2acrsm1_inm1_size_y,  ism2acrsm1_inm1_size_x);
  }

//  6:    update not-0-area indicators 

//==//matr2->min_not_0_y_index+=shf2_y;  matr2->max_not_0_y_index+=shf2_y;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)
//==//matr2->min_not_0_x_index+=shf2_x;  matr2->max_not_0_x_index+=shf2_x;      //  ! ! ! we don't do it because shift is realized inside of matr2->shift_elements_spec_fast_case(shf2_y,  shf2_x)


  if ((ism2_inm1_size_y != 0) && (ism2_inm1_size_x != 0))
  {
    if ((matr1->not_0_len_y != 0) && (matr1->not_0_len_x != 0))
    {
      if (ism2_inm1_y0                < matr1->min_not_0_y_index ) {matr1->min_not_0_y_index=ism2_inm1_y0;   }
      if (matr1->max_not_0_y_index    < ism2_inm1_y_end          ) {matr1->max_not_0_y_index=ism2_inm1_y_end;}
      if (ism2_inm1_x0                < matr1->min_not_0_x_index ) {matr1->min_not_0_x_index=ism2_inm1_x0;   }
      if (matr1->max_not_0_x_index    < ism2_inm1_x_end          ) {matr1->max_not_0_x_index=ism2_inm1_x_end;}
    }
    else {matr1->min_not_0_y_index=ism2_inm1_y0;  matr1->max_not_0_y_index=ism2_inm1_y_end;     matr1->min_not_0_x_index=ism2_inm1_x0;  matr1->max_not_0_x_index=ism2_inm1_x_end;}

  matr1->not_0_len_y=matr1->max_not_0_y_index-matr1->min_not_0_y_index+1;
  matr1->not_0_len_x=1;
  }


  if ((ism1_inm2_size_y != 0) && (ism1_inm2_size_x != 0))
  {
    if ((matr2->not_0_len_y != 0) && (matr2->not_0_len_x != 0))
    {
      if (ism1_inm2_y0                < matr2->min_not_0_y_index ) {matr2->min_not_0_y_index=ism1_inm2_y0;   }
      if (matr2->max_not_0_y_index    < ism1_inm2_y_end          ) {matr2->max_not_0_y_index=ism1_inm2_y_end;}
    matr2->min_not_0_x_index=0;   
    matr2->max_not_0_x_index=0;
    }
    else
    {matr2->min_not_0_y_index=ism1_inm2_y0;  matr2->max_not_0_y_index=ism1_inm2_y_end;    matr2->min_not_0_x_index=0;  matr2->max_not_0_x_index=0;}

  matr2->not_0_len_y=matr2->max_not_0_y_index-matr2->min_not_0_y_index+1;
  matr2->not_0_len_x=1;   
  }

matr1->min_not_0_x_index=0; matr1->max_not_0_x_index=0;  matr1->not_0_len_x=1;
matr2->min_not_0_x_index=0; matr2->max_not_0_x_index=0;  matr2->not_0_len_x=1;

//==DEBUG==//  cout<<"DEBUG  became:   matr1"<<endl;  matr1->display_on_screen_not_0_area_size_param();  //==DEBUG==//  
//==DEBUG==//  cout<<"DEBUG  became:   matr2"<<endl;  matr2->display_on_screen_not_0_area_size_param();  //==DEBUG==//  



return true;

}



//  --------------------------------------------------------------------------------------------------------------------------------------
//  ======================================================================================================================================
//  --------------------------------------------------------------------------------------------------------------------------------------




int optim_enumer_one_line_pbc(unsigned long long int & state, int len, bool & start_neo_count_and_do_multipl)
{

// int shM=bit1_count_row0+amount_of_bit1(row_1_state);


  if (state == 0 ) {state=1;  start_neo_count_and_do_multipl=true;  return len;}                                                // ......000001
  if (state == 1 ) {state=3;  start_neo_count_and_do_multipl=true;  if (len > 1) {return len-1;    } else {return 1;}  }        // ......000011
  if (state == 3 ) {state=5;  start_neo_count_and_do_multipl=true;  if (len > 2) {return len-2;    } else {return 1;}  }        // ......000101
  if (state == 5 ) {state=7;  start_neo_count_and_do_multipl=true;  if (len > 2) {return len-2;    } else {return 1;}  }        // ......000111
  if (state == 7 ) {state=9;  start_neo_count_and_do_multipl=true;  if (len > 3) {return len-3;    } else {return 1;}  }        // ......001001
  if (state == 9 ) {state=11; start_neo_count_and_do_multipl=true;  if (len > 3) {return 2*(len-3);} else {return 1;}  }        // ......001011
  if (state == 11) {state=15; start_neo_count_and_do_multipl=true;  if (len > 3) {return (len-3);  } else {return 1;}  }        // ......001111
  if (state == 15) {state=17; start_neo_count_and_do_multipl=true;  if (len > 4) {return (len-4);  } else {return 1;}  }        // ......010001

//  if (state >= 17) {...     already is

start_neo_count_and_do_multipl=false;



//  we nned to get info about state: (1) amount of bit1,  (2) pos of forw bit1,  (3) pos of bacw bit 1

int max_forw_range_len=20;
int range_len=3, forw_range=1,  backw_range=1,  mid_range=1;
int forw_bit_pos=0; 


unsigned long long int shifted_bit1=1;
shifted_bit1=(shifted_bit1 << (len-1));

//  define forw bit 1
  for (int pos=len-1; pos >= 0; pos--)
  {  if ((shifted_bit1 & state) != 0)  {forw_bit_pos=pos; break;}   shifted_bit1=(shifted_bit1 >> 1);}

//  ***
range_len=forw_bit_pos-1;
forw_range=range_len/2;    //  forw_range=(range_len-1)/2;
  if (forw_range > max_forw_range_len) {forw_range=max_forw_range_len;}
backw_range=forw_range;  mid_range=range_len-forw_range-backw_range; 

//  ***
unsigned long long int  forw_num=0, mid_num=0, backw_num=0;
unsigned long long int  forw_num_max=(((unsigned long long int ) 1) << (forw_range))-1,  mid_num_max=(((unsigned long long int ) 1) << (mid_range))-1,  backw_num_max=(((unsigned long long int ) 1) << (backw_range))-1;
unsigned long long int  backw_num_max_div_2=backw_num_max/2;


//  ***
write_bits_ullint(forw_num, forw_range-1,   state, forw_bit_pos-1,   forw_range);

//  ***
write_bits_ullint(mid_num,  mid_range-1,   state, forw_bit_pos-1-forw_range,  mid_range);

//  ***
write_bits_ullint(backw_num,  backw_range-1,   state, forw_bit_pos-1-forw_range-mid_range,   backw_range);



//  analizing

//  *** 
//  test mid_num

  if (mid_range > 0)
  {
  mid_num++;

    if (mid_num > mid_num_max) {mid_num=0;} 
    if (mid_num != 0) {write_bits_ullint(state, forw_bit_pos-1-forw_range,   mid_num, mid_range-1,   mid_range);  return 1;} 
  
  }



start_neo_count_and_do_multipl=true; 
//cout<<"forw_num="<<forw_num<<"  mid_num="<<mid_num<<"  backw_num="<<backw_num<<"  range_len="<<range_len<<endl;
//cout<<"forw_num_max="<<forw_num_max<<"  mid_num_max="<<mid_num_max<<"  backw_num_max="<<backw_num_max<<"  range_len="<<range_len<<endl;


//  *** 
//  test backw_num
int degen_factor=1; 
unsigned long long int  forw_invert_num=0;

  if (backw_range > 0)
  { 
  bool exit=false; 

    while (exit == false)
    {
    backw_num++;

    // (1) increment backw_num and forw_num
      if (backw_num > backw_num_max) 
      {
      backw_num=0;  forw_num++;  

        if (forw_num > forw_num_max) 
        {
        //  update   forw_range, mid_range, backw_range,  range_len,  forw_num, mid_num, backw_num
        forw_bit_pos++;
        state=1;  state=(state << forw_bit_pos);  state+=1;  exit=true;
        forw_num=0;  mid_num=0;  backw_num=0;  range_len++; 
        forw_range=range_len/2;  //  forw_range=(range_len-1)/2; 
        backw_range=forw_range;  mid_range=range_len-forw_range-backw_range;
          if (forw_range > max_forw_range_len) {forw_range=max_forw_range_len;}
        backw_range=forw_range;  mid_range=range_len-forw_range-backw_range;
          if (range_len+2 > len) {return 1;}      //  !!! end of enumeration  !!!  
        }

      }

    // (2) analising
    forw_invert_num=bit_inverse_ullint(forw_num, forw_range);

      if (forw_invert_num == backw_num) {exit=true;}
      else
      {
        if (forw_invert_num < backw_num) {degen_factor=2;  exit=true;}
      }

    }      //  while (exit == false)

  }  
  else
  {
  return (len-range_len-1);
  }


write_bits_ullint(state, forw_bit_pos-1-forw_range-mid_range,   backw_num, backw_range-1,   backw_range);
write_bits_ullint(state, forw_bit_pos-1-forw_range,   mid_num, mid_range-1,   mid_range);      //  already have done above
write_bits_ullint(state, forw_bit_pos-1,                        forw_num,  forw_range-1,    forw_range);


return degen_factor*(len-range_len-1);

}




void matr_symmetr_reflect(matr_t1 * matr_dest,  matr_t1 * matr_src,  bool is_refl_last_col)
{


#ifdef __MACR_MATR_TYPE_1__MPZ


int x_end=matr_src->max_not_0_x_index;

  if ((is_refl_last_col == false) && (matr_src->max_not_0_x_index+1 == matr_src->size_x))
  {x_end=matr_src->max_not_0_x_index-1;}
  

  for (int index_y=matr_src->min_not_0_y_index;  index_y <= matr_src->max_not_0_y_index;  index_y++)
  {

    for (int index_x=matr_src->min_not_0_x_index;  index_x <= x_end;  index_x++)
    {
    int matr_src_elem_index=index_y*matr_src->size_x+index_x;

      if (mpz_sgn(matr_src->element[matr_src_elem_index]) != 0)
      {
      int matr_dest_elem_index=index_y*matr_dest->size_x+index_x;
      mpz_set(matr_dest->element[matr_dest_elem_index],  matr_src->element[matr_src_elem_index]);
      mpz_set(matr_dest->element[matr_dest_elem_index-index_x+matr_dest->size_x-index_x-1],  matr_src->element[matr_src_elem_index]);
      }

    }

  } 


matr_dest->min_not_0_y_index=matr_src->min_not_0_y_index;  matr_dest->max_not_0_y_index=matr_dest->size_y-matr_dest->min_not_0_y_index-1;   
matr_dest->min_not_0_x_index=matr_src->min_not_0_x_index;  matr_dest->max_not_0_x_index=matr_dest->size_x-matr_dest->min_not_0_x_index-1;   
matr_dest->not_0_len_y=matr_dest->max_not_0_y_index-matr_dest->min_not_0_y_index+1;
matr_dest->not_0_len_x=matr_dest->max_not_0_x_index-matr_dest->min_not_0_x_index+1;

#endif 

}



void matr_yaxis_invers_symmetr_selfadd(matr_t1 * matr)
{


#ifdef __MACR_MATR_TYPE_1__MPZ

matr_t1 matr_temp; 
matr_temp.init(matr->size_y,  matr->size_x, 0);

//  matr  to   matr_temp

  for (int index_y=matr->min_not_0_y_index;  index_y <= matr->max_not_0_y_index;  index_y++)
  {
    for (int index_x=matr->min_not_0_x_index;  index_x <= matr->max_not_0_x_index;  index_x++)
    {
    int matr_index=index_y*matr->size_x+index_x;

      if (mpz_sgn(matr->element[matr_index]) != 0)
      {
      mpz_add(matr_temp.element[matr_index], matr_temp.element[matr_index], matr->element[matr_index]);
      mpz_add(matr_temp.element[matr_index+matr->size_x-index_x-index_x-1], matr_temp.element[matr_index+matr->size_x-index_x-index_x-1], matr->element[matr_index]);
      }

    }

  } 


int matr_distance_to_right_wall=matr->size_x-matr->max_not_0_x_index-1; 
int matr_distance_to_left_wall=matr->min_not_0_x_index;


matr_temp.min_not_0_y_index=matr->min_not_0_y_index;   matr_temp.max_not_0_y_index=matr->max_not_0_y_index;    matr_temp.not_0_len_y=matr->not_0_len_y;   
matr_temp.min_not_0_x_index=matr->min_not_0_x_index;   matr_temp.max_not_0_x_index=matr->max_not_0_x_index;    matr_temp.not_0_len_x=matr->not_0_len_x; 

  if (matr_distance_to_right_wall < matr_distance_to_left_wall) {matr_temp.min_not_0_x_index=matr_distance_to_right_wall;}
  if (matr_distance_to_right_wall > matr_distance_to_left_wall) {matr_temp.max_not_0_x_index=matr->size_x-matr_distance_to_left_wall-1;}

matr_temp.not_0_len_x=matr_temp.max_not_0_x_index-matr_temp.min_not_0_x_index+1;

//  matr_temp  to   matr

matr->set_val_for_all_el(0);
matr1_add(matr,  & matr_temp); 

//  matr->min_not_0_y_index=matr_temp.min_not_0_y_index;   matr->max_not_0_y_index=matr_temp.max_not_0_y_index;    matr->not_0_len_y= matr->max_not_0_y_index-matr->min_not_0_y_index+1; 
//  matr->min_not_0_x_index=matr_temp.min_not_0_x_index;   matr->max_not_0_x_index=matr_temp.max_not_0_x_index;    matr->not_0_len_x= matr->max_not_0_x_index-matr->min_not_0_x_index+1;

matr_temp.free();

#endif 

}



void matr_half_to_whole_fill(matr_t1 * matr_dest,  matr_t1 * matr_src)
{

#ifdef __MACR_MATR_TYPE_1__MPZ

bool is_even_len_x=false;
int len_x_div_2=matr_dest->size_x/2;

  if (matr_dest->size_x % 2 == 0) {is_even_len_x=true;}


//  matr  to   matr_temp


//  (1) copy matr_src
int index_x_end=matr_src->max_not_0_x_index;
  if (index_x_end > len_x_div_2-1) {index_x_end=len_x_div_2-1;  if (is_even_len_x == false) {index_x_end++;}  }

  for (int index_y=matr_src->min_not_0_y_index;  index_y <= matr_src->max_not_0_y_index;  index_y++)
  {
    for (int index_x=matr_src->min_not_0_x_index;  index_x <= index_x_end;  index_x++)
    {
    int matr_index_src =index_y*matr_src->size_x+index_x;
    int matr_index_dest=index_y*matr_dest->size_x+index_x;

      if (mpz_sgn(matr_src->element[matr_index_src]) != 0)
      {
      mpz_add(matr_dest->element[matr_index_dest], matr_dest->element[matr_index_dest], matr_src->element[matr_index_src]);     
      }

    }

  } 


//  (2) add-shift right in  x=1
//== here not ==//index_x_end=matr_src->max_not_0_x_index;
//== here not ==//  if (index_x_end > len_x_div_2-1) {index_x_end=len_x_div_2-1;  if (is_even_len_x == false) {index_x_end++;}  }
//== here not ==//  if (is_even_len_x == false) {if (index_x_end >= len_x_div_2) {index_x_end--;}} else {if (index_x_end >= len_x_div_2-1) {index_x_end--;}}

//== here not ==//  for (int index_y=matr_src->min_not_0_y_index;  index_y <= matr_src->max_not_0_y_index;  index_y++)
//== here not ==//  {
//== here not ==//    for (int index_x=matr_src->min_not_0_x_index;  index_x <= index_x_end;  index_x++)
//== here not ==//    {
//== here not ==//    int matr_index_src =index_y*matr_src->size_x+index_x;
//== here not ==//    int matr_index_dest=index_y*matr_dest->size_x+index_x+1;
//== here not ==//
//== here not ==//      if (mpz_sgn(matr_src->element[matr_index_src]) != 0)
//== here not ==//      {
//== here not ==//      mpz_add(matr_dest->element[matr_index_dest], matr_dest->element[matr_index_dest], matr_src->element[matr_index_src]);     
//== here not ==//      }
//== here not ==//
//== here not ==//    }
//== here not ==//
//== here not ==//  } 


//  (3) add-shift right in  x=1

//  if ((is_even_len_x  == false) && (matr_src->max_not_0_x_index+1 == matr_src->size_x))
//  {x_end=matr_src->max_not_0_x_index-1;}
  
  for (int index_y=matr_src->min_not_0_y_index;  index_y <= matr_src->max_not_0_y_index;  index_y++)
  {

    for (int index_x=matr_src->min_not_0_x_index;  index_x < len_x_div_2;  index_x++)
    {
    //int matr_src_elem_index=index_y*matr_dest->size_x+index_x;
    int matr_dest_elem_index=index_y*matr_dest->size_x+index_x; 

      if (mpz_sgn(matr_dest->element[matr_dest_elem_index]) != 0)
      {
      int simmetr_index=matr_dest_elem_index-index_x+matr_dest->size_x-index_x-1;
      mpz_add(matr_dest->element[simmetr_index],  matr_dest->element[simmetr_index],  matr_dest->element[matr_dest_elem_index]);
      }

    }

  } 


matr_dest->min_not_0_y_index=matr_src->min_not_0_y_index;  matr_dest->max_not_0_y_index=matr_src->max_not_0_y_index;   
matr_dest->min_not_0_x_index=matr_src->min_not_0_x_index;  matr_dest->max_not_0_x_index=matr_dest->size_x-matr_dest->min_not_0_x_index-1;   
matr_dest->not_0_len_y=matr_dest->max_not_0_y_index-matr_dest->min_not_0_y_index+1;
matr_dest->not_0_len_x=matr_dest->max_not_0_x_index-matr_dest->min_not_0_x_index+1;

#endif 


}








// ----------------------------------------------------
// in out procedures
// ----------------------------------------------------
//
//
//     [ O ]
//       \ \      p
//        \ \  \o/
//         \ \--'---_
//         /\ \   / ~~\_
//   ./---/__|=/_/------//~~~\
//  /___________________/O   O \
//  (===(\_________(===(Oo o o O)        
//   \~~~\____/     \---\Oo__o--
//     ~~~~~~~       ~~~~~~~~~~
//
//
//
//  format of out data
//
//  #  dim
//  2
//  #  lattice type
//  rect
//  #  bound cond
//  fbc
//  #  Ly
//  6
//  #  Lx
//  6
//  #  amount of threads 
//  4
//  #  cur thread number 
//  -1
//  #  last completed row main state
//  1636111
//  #  lattice struct
//  |    |    |    |
//  * -- * -- * -- * -- 
//  |    |    |    |
//  * -- * -- * -- * --
//  |    |    |    |
//  * -- * -- * -- * --
//  |    |    |    |
//  * -- * -- * -- * --
//  #  min energy
//  -125
//  #  J+ percent
//  20
//  #  J- percent
//  80
//  #  J table
//  #  s1_index_y	s1_index_x	s2_index_y	s2_index_x	bond_mean
//  5		        1		0		1		-1
//  5		        1		5		2		+1
//  5		        2		0		2		+1
//  #  J-matr
//    +   +   +   +   +   +   +   +
//  + * - * - * + * - * - * + * - * +
//    -   +   -   +   -   -   +   -
//  + * - * + * - * + * - * - * - * +
//    -   -   +   +   +   -   +   +
//  + * + * - * - * + * + * + * + * +
//    -   +   -   +   +   +   +   -
//  + * + * + * - * - * - * + * + * +
//    -   +   -   +   +   -   -   +
//  + * + * + * - * - * + * + * + * +
//    -   +   +   -   +   +   +   -
//  + * - * + * + * - * - * - * + * +
//    +   -   +   -   +   -   +   +
//  + * + * - * - * - * - * + * + * +
//    -   +   -   -   -   -   +   +
//  + * - * + * - * - * - * - * + * +
//    +   +   +   +   +   +   +   +
//  #  spin configuration for
//  ?????
//  #  spin configuration:
//  # s_y   s_x     spin
//  0       0       -1
//  0       1       -1
//  1       0       +1
//  1       1       -1
//  #  hist E M g 
//  #  ================================================================
//  -13   -19   971274912840123402394
//  -35   -62   43895723845244
//  -23   -65   19284078234982347823
//  -18   -90   87395782398752982347823
//  #  subhist for boundary state
//  0
//  #  subhist start 
//  #  ================================================================
//  -2   -9    092397856199010
//  -13  -73   3284
//  #  subhist end
//  #  subhist for boundary state
//  1
//  #  subhist start 
//  #  ================================================================
//  -3   -15    52311234
//  -24  -88    1122
//  # subhist end 
//  #  subhist for boundary state
//  2
//  #  subhist start 
//  #  ================================================================
//  -89   -75    121
//  -12   -67    8564543
//  # subhist end 
//
//  writing and extracting data
//  save J matr to file:       izing_EA_2D_rect_fbc_06x06__J                      +postfix    + .txt
//  save hist matr to file:    izing_EA_2D_rect_fbc_06x06___thrd_06_of_16__hist   +postfix    + .txt
//








int matr_t1_to_chstr_hist(matr_t1 * matr_hist,  char *&  chstr_text,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment,  int Ly,  int Lx)
{

  if (matr_hist == 0) {return 0;}

  if (chstr_text) {return 0;}


string string_var;


  if (print_comment) 
  {
  string_var.append((char *)  "#  hist E M g\n");
  //string_var.append((char *)  "#  Ly:\n");
  //add_int_to_string(& string_var, Ly);
  //string_var.append((char *)  "\n");
  //string_var.append((char *)  "#  Lx:\n");
  //add_int_to_string(& string_var, Lx);
  //string_var.append((char *)  "\n");
  string_var.append((char *)  "#  ================================================================\n");
  }


int index_y_start=matr_hist->min_not_0_y_index,  index_y_end=matr_hist->max_not_0_y_index,  len_y=index_y_end-index_y_start+1,  index_dy=1;
int index_x_start=matr_hist->min_not_0_x_index,  index_x_end=matr_hist->max_not_0_x_index,  len_x=index_x_end-index_x_start+1,  index_dx=1;

int E_y_start=E_down+dif_E_betw_y_cell*index_y_start,  E_y_end=E_y_start+(len_y-1)*dif_E_betw_y_cell,  E_dy=dif_E_betw_y_cell;
int M_x_start=M_left+dif_M_betw_x_cell*index_x_start,  M_x_end=M_x_start+(len_x-1)*dif_M_betw_x_cell,  M_dx=dif_M_betw_x_cell; 

  if (dif_E_betw_y_cell < 0) 
  {
  index_y_start=matr_hist->max_not_0_y_index;  index_y_end=matr_hist->min_not_0_y_index;  len_y=index_y_start-index_y_end+1;  index_dy=-1;
  E_y_start=E_down+dif_E_betw_y_cell*index_y_start;  E_y_end=E_y_start+(len_y-1)*dif_E_betw_y_cell;  E_dy=dif_E_betw_y_cell;
  }

  if (dif_M_betw_x_cell < 0) 
  {
  index_x_start=matr_hist->max_not_0_x_index;  index_x_end=matr_hist->min_not_0_x_index;  len_x=index_x_start-index_x_end+1;  index_dx=-1;
  M_x_start=M_left+dif_M_betw_x_cell*index_x_start;  M_x_end=M_x_start+(len_x-1)*dif_M_betw_x_cell;  M_dx=dif_M_betw_x_cell;
  }


char chstr_num[512];  chstr_num[0]='\0'; 

int E_y=E_y_start;
int M_x=M_x_start;

  for (int index_y=index_y_start;  index_y <= index_y_end;  index_y+=index_dy)
  {
  M_x=M_x_start;

    for (int index_x=index_x_start;  index_x <= index_x_end;  index_x+=index_dx)
    {
    int matr_src_elem_index=index_y*matr_hist->size_x+index_x;

      if (mpz_sgn(matr_hist->element[matr_src_elem_index]) != 0)
      { 
      add_int_to_string(& string_var, E_y);
      string_var.append((char *)  "  \t"); 
      add_int_to_string(& string_var, M_x);
      string_var.append((char *)  "  \t");
      matr_hist->get_element_as_chstr_stat(index_y,  index_x,  chstr_num, -1);
      string_var.append(chstr_num);
      chstr_num[0]='\0';
      string_var.append((char *)  "\n");
      }

    M_x+=M_dx;
    }

  E_y+=E_dy;
  } 



chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';

return strlen(chstr_text);

}




int matr_t1_to_file(matr_t1 * matr_hist,  char * chstr_filename,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment,  int Ly,  int Lx,  bool newfile)
{

char * chstr_text=0;
int text_size=matr_t1_to_chstr_hist(matr_hist,  chstr_text,  E_down,  M_left,  dif_E_betw_y_cell,  dif_M_betw_x_cell,  print_comment,  Ly,  Lx);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}




void matr_t1_to_screen(matr_t1 * matr_hist,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment,  int Ly,  int Lx)
{

char * chstr_text=0;
int text_size=matr_t1_to_chstr_hist(matr_hist,  chstr_text,  E_down,  M_left,  dif_E_betw_y_cell,  dif_M_betw_x_cell,  print_comment,  Ly,  Lx);
cout<<chstr_text<<endl;

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}




int matr_t1_to_short_display(matr_t1 * matr_hist,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment,  int Ly,  int Lx)
{

char *chstr_text=0;
  

  if (matr_hist == 0) {return 0;}

  if (chstr_text) {return 0;}


string string_var;



string_var.append((char *)  "#  hist E M g\n");
string_var.append((char *)  "#  ================================================================\n");


int index_y_start=matr_hist->min_not_0_y_index,  index_y_end=matr_hist->max_not_0_y_index,  len_y=index_y_end-index_y_start+1,  index_dy=1;
int index_x_start=matr_hist->min_not_0_x_index,  index_x_end=matr_hist->max_not_0_x_index,  len_x=index_x_end-index_x_start+1,  index_dx=1;

int E_y_start=E_down+dif_E_betw_y_cell*index_y_start,  E_y_end=E_y_start+(len_y-1)*dif_E_betw_y_cell,  E_dy=dif_E_betw_y_cell;
int M_x_start=M_left+dif_M_betw_x_cell*index_x_start,  M_x_end=M_x_start+(len_x-1)*dif_M_betw_x_cell,  M_dx=dif_M_betw_x_cell; 

  if (dif_E_betw_y_cell < 0) 
  {
  index_y_start=matr_hist->max_not_0_y_index;  index_y_end=matr_hist->min_not_0_y_index;  len_y=index_y_start-index_y_end+1;  index_dy=-1;
  E_y_start=E_down+dif_E_betw_y_cell*index_y_start;  E_y_end=E_y_start+(len_y-1)*dif_E_betw_y_cell;  E_dy=dif_E_betw_y_cell;
  }

  if (dif_M_betw_x_cell < 0) 
  {
  index_x_start=matr_hist->max_not_0_x_index;  index_x_end=matr_hist->min_not_0_x_index;  len_x=index_x_start-index_x_end+1;  index_dx=-1;
  M_x_start=M_left+dif_M_betw_x_cell*index_x_start;  M_x_end=M_x_start+(len_x-1)*dif_M_betw_x_cell;  M_dx=dif_M_betw_x_cell;
  }



int counter=0;
char chstr_num[512];  chstr_num[0]='\0';
int E_y=E_y_start;
int M_x=M_x_start;
string string_line;


  for (int index_y=index_y_start;  index_y <= index_y_end;  index_y+=index_dy)
  {
  M_x=M_x_start;

    for (int index_x=index_x_start;  index_x <= index_x_end;  index_x+=index_dx)
    {
    int matr_src_elem_index=index_y*matr_hist->size_x+index_x;

      if (mpz_sgn(matr_hist->element[matr_src_elem_index]) != 0) 
      { 
      string_line.clear();
      add_int_to_string(& string_line, E_y); 
      string_line.append((char *)  "\t");
      add_int_to_string(& string_line, M_x);
      string_line.append((char *)  "\t"); 
      matr_hist->get_element_as_chstr_stat(index_y,  index_x,  chstr_num);
      string_line.append(chstr_num);  
      chstr_num[0]='\0';      
      string_line.append((char *)  "\n"); 

        if (counter < 3) {string_var.append(string_line);}
  
      counter++;
      } 

    M_x+=M_dx;
    }

  E_y+=E_dy;
  } 


  if (counter >= 3) 
  {
  string_var.append((char *)  ".......................................\n");
  string_var.append(string_line);
  }
string_var.append((char *)  "#  ====================================\n");


chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';
cout<<chstr_text;
int len=strlen(chstr_text);

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}


return len;

}





bool extract_2D_rect_lattice_hist_from_file(char * filename, int Ly, int Lx, matr_t1 *matr, int E_min, int E_max, int M_min, int M_max, int dE_dcell, int dM_dcell)
{

  if (filename == 0) {return false;}
  if (matr == 0) {return false;}
  if ((E_max < E_min) || (M_max < M_min) || (dE_dcell < 1) || (dM_dcell < 1)) {return false;}


//  reading file
char * chstr_text=0;
int file_size=write_file_to_chstr(filename, chstr_text,  0);
  if (file_size < 1) {return 0;} 
bool ok=false;     



//  ----------------------------------------------------------------------------
//  defining:       extr_table_startpos,   extr_table_endpos, extr_table_len
//
int extr_table_startpos=0, extr_table_endpos=0, extr_table_len=0;

int  key_phr_line_number=0;                 //  n - next
int  key_phr_line_number_line_n=0;          //  n - next
int  key_phr_line_number_line_nn=0;         //  n - next

int line_end_pos_temp=0;


key_phr_line_number=sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, (char *) "#  hist E M g", 0); 


  if (key_phr_line_number >= 0) 
  {ok=sr_line_position_and_len_in_chstr(chstr_text,   key_phr_line_number, 0, & line_end_pos_temp, 0);}
      

  if (ok == true)
  {
  key_phr_line_number_line_n =sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, (char *) "#", line_end_pos_temp+1);              

    if (key_phr_line_number_line_n >= 0)            
    {ok=sr_line_position_and_len_in_chstr(chstr_text,  key_phr_line_number_line_n, 0, & line_end_pos_temp, 0);}      //  so we know that hist data starts with:    key_phr_line_number_line_n+1
    else {ok=false;}

    if (ok == true)
    {
    key_phr_line_number_line_nn=sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, (char *) "#", line_end_pos_temp+1);  

      if (key_phr_line_number_line_nn >= 0)
      {
      ok=sr_line_position_and_len_in_chstr(chstr_text,  key_phr_line_number_line_nn-1,  0,  & line_end_pos_temp,  0);  
      extr_table_endpos=line_end_pos_temp;
      }
      else
      {
      extr_table_endpos=strlen(chstr_text)-1; 
      }	

    }

  }
  

  if (ok == true) 
  {ok=sr_line_position_and_len_in_chstr(chstr_text,   key_phr_line_number_line_n+1, & extr_table_startpos, 0, 0);}

  if (ok == false) {  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}  return false;}
	
extr_table_len=extr_table_endpos-extr_table_startpos+1;
//  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}      //==// DEBUG

//  done ! ! !




//  ----------------------------------------------------------------------------
//  copying text chunk analizing and extracting phrases
//  format int extract_phrases_from_chstr(char *src_text, char ** & chstr_phrases, int & amount_of_extracted_phrases,  int & phrase_len,  int start_phrase_num,  int col_count=-1);

char remember_smb=chstr_text[extr_table_startpos+extr_table_len];      //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!
chstr_text[extr_table_startpos+extr_table_len]='\0';                   //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!

//  E phrases extracting
char ** chstr_phrases_E=0;
int  amount_of_extracted_phrases_E=0;
int  max_phrase_len_E=0;
extract_phrases_from_chstr(chstr_text+extr_table_startpos, chstr_phrases_E,  amount_of_extracted_phrases_E,  max_phrase_len_E,  0,  3);

//  M phrases extracting
char ** chstr_phrases_M=0;
int  amount_of_extracted_phrases_M=0;
int  max_phrase_len_M=0;
extract_phrases_from_chstr(chstr_text+extr_table_startpos, chstr_phrases_M,  amount_of_extracted_phrases_M,  max_phrase_len_M,  1,  3);

//  g phrases extracting
char ** chstr_phrases_g=0;      //  !  !  !
int  amount_of_extracted_phrases_g=0;
int  max_phrase_len_g=0;
extract_phrases_from_chstr(chstr_text+extr_table_startpos, chstr_phrases_g,  amount_of_extracted_phrases_g,  max_phrase_len_g,  2,  3);

chstr_text[extr_table_startpos+extr_table_len]=remember_smb;          //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}
  if (amount_of_extracted_phrases_E != amount_of_extracted_phrases_M) {ok=false;}
  if (amount_of_extracted_phrases_E != amount_of_extracted_phrases_g) {ok=false;}




//  ----------------------------------------------------------------------------
//  putting to matr

//  matr->free();
//  matr->init((E_max-E_min+1)/dE_dcell+1,  (M_max-M_min+1)/dM_dcell+1,  0);

  if (ok == true)
  { 
  int M_amount=(M_max-M_min+1)/dM_dcell;
    if (M_max-M_min+1-dM_dcell*M_amount > 0) {M_amount++;}

  //  cout<<"DEBUG:  M_amount="<<M_amount<<endl;
  int index_y=0,  index_x=0;
  int extracted_E=0, extracted_M=0;
  mpz_t temp_mpz_var;
  mpz_init2(temp_mpz_var, __MACR_MPZ_BIT_SIZE);

    for (int i=0;  i < amount_of_extracted_phrases_g;  i++)
    {
    extracted_E=atoi(chstr_phrases_E[i]);  extracted_M=atoi(chstr_phrases_M[i]);
    index_y=(extracted_E-E_min)/dE_dcell;  index_x=(extracted_M-M_min)/dM_dcell;
    mpz_set_str(temp_mpz_var, chstr_phrases_g[i], 10);
    mpz_add(matr->element[index_y*M_amount+index_x], matr->element[index_y*M_amount+index_x], temp_mpz_var);  
    }

  mpz_clear(temp_mpz_var);
  }



//  ----------------------------------------------------------------------------
//  cleaning

free_chstr_array(chstr_phrases_E, amount_of_extracted_phrases_E); 
free_chstr_array(chstr_phrases_M, amount_of_extracted_phrases_M); 
free_chstr_array(chstr_phrases_g, amount_of_extracted_phrases_g); 

  if (ok == false) {return false;}


return true;

}






bool read_J_bonds_from_file(char * chstr_filename_with_J_bonds,  bool * & J_ver,  bool * & J_hor,  int & Ly,  int & Lx, bool is_1pbc_0fbc)
{

  if (J_ver) {return false;}
  if (J_hor) {return false;}
  if (chstr_filename_with_J_bonds == 0) {return false;}


int posit_bond_perc=0,  negat_bond_perc=0;
int Ly_=0,  Lx_=0;        
bool is_extr_src_info=extract_2D_rect_fbc_lattice_serv_info_from_file_v15(chstr_filename_with_J_bonds,  & Ly_,  & Lx_, 0, 0, 0,  & posit_bond_perc,  & negat_bond_perc, 0);      //  !  !  !
Ly=Ly_;  Lx=Lx_;

//  J_ver=new bool[Lx*Ly+Lx];
//  J_hor=new bool[Lx*Ly+Ly];  


bool is_J_extracted=false;

  if ((is_extr_src_info) && (Ly > 1) && (Lx > 1))
  {
    if (is_1pbc_0fbc == true) {is_J_extracted=extract_2D_rect_lattice_J_table_from_file_dyn(chstr_filename_with_J_bonds, Ly, Lx, J_ver, J_hor, 'c', 'c', 'c', 'c');}
    else {is_J_extracted=extract_2D_rect_lattice_J_table_from_file_dyn(chstr_filename_with_J_bonds, Ly, Lx, J_ver, J_hor, 'n', 'n', 'n', 'n');}
  
    if (is_J_extracted == false)
    {  
      if (J_ver == 0) {delete[] J_ver;  J_ver=0;}  
      if (J_hor == 0) {delete[] J_hor;  J_hor=0;}  
    return false;
    } 

  }

return true;

}








void  st_spin_config_int_ar :: sort_by_M()
{

  if (current_g > 1)
  {
  bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
  bool ok=false;
  short int M1=0,  M2=0;

  //  bubble sorting
  //
  unsigned long long int permut_amount=0;
  unsigned short int state_buf=0;
  bool is_exit=false;

    for (unsigned long long int counter=0;  (counter < current_g-1) && (is_exit == false) ;  counter++)
    { 
    permut_amount=0;
  
      for (unsigned long long int i=0;  i < (current_g-1-counter);  i++)
      {
      // cout<<"i="<<i<<endl;
      convert_bit_row_state_to_bool_ar((unsigned long long int) i  ,  spin_array); 
      M1=calc_magnetiz_lat(spin_array,   __MACR_LATTICE_LY,   __MACR_LATTICE_LX);
      convert_bit_row_state_to_bool_ar((unsigned long long int) i+1,  spin_array); 
      M2=calc_magnetiz_lat(spin_array,   __MACR_LATTICE_LY,   __MACR_LATTICE_LX);
    
        if (M1 > M2) 
        {
          for (int row=0;  row < __MACR_LATTICE_LY;  row++)
          {state_buf=bit_row_state[i*__MACR_LATTICE_LY+row];  bit_row_state[i*__MACR_LATTICE_LY+row]=bit_row_state[(i+1)*__MACR_LATTICE_LY+row];  bit_row_state[(i+1)*__MACR_LATTICE_LY+row]=state_buf;}
        
        permut_amount++;   
        }      
  
      }
    
      if (permut_amount < 1) {is_exit=true;}
    }
  
  }

}










unsigned long long int  st_spin_config_int_ar :: replace_row_state_ar_with_nouv(const st_spin_config_int_ar & o_spin_config_int_ar_arg)
{

  if (o_spin_config_int_ar_arg.current_g < 1) {E=0;  current_g=0;  was_overfilled=false;  return 0;}

//  std::memcpy(bit_row_state, o_spin_config_int_ar_arg.bit_row_state, o_spin_config_int_ar_arg.current_g*sizeof(int)*__MACR_LATTICE_LY);      //  !  !  !
  
  for (unsigned long long int i=0;  i < o_spin_config_int_ar_arg.current_g*__MACR_LATTICE_LY;  i++)
  {bit_row_state[i]=o_spin_config_int_ar_arg.bit_row_state[i];}

E=o_spin_config_int_ar_arg.E;  current_g=o_spin_config_int_ar_arg.current_g;  was_overfilled=o_spin_config_int_ar_arg.was_overfilled;

return current_g; 

}





unsigned long long int st_spin_config_int_ar :: add_row_state_nouv(const st_spin_config_int_ar & o_spin_config_int_ar_arg)
{  //  current_cycle_out=0;  wished_cycle_out=0;

  if (o_spin_config_int_ar_arg.current_g < 1) {return 0;}


unsigned long long int  dif=treat_out_of_ullint_2_num_befor_add(current_g, o_spin_config_int_ar_arg.current_g, & was_overfilled);      //  insurance not to out of ULLINT type    //  unsigned long long int  dif=o_spin_config_int_ar_arg.current_g;

  if ((dif+current_g) > max_amount_of_spin_config) {dif=max_amount_of_spin_config-current_g;  was_overfilled=true;}

//std::memcpy(bit_row_state+current_g*sizeof(int)*__MACR_LATTICE_LY, o_spin_config_int_ar_arg.bit_row_state, dif*sizeof(int)*__MACR_LATTICE_LY);      //  !  !  !
//  current_cycle_out=0;  wished_cycle_out=0;

  for (unsigned long long int i=0;  i < dif*__MACR_LATTICE_LY;  i++)
  {bit_row_state[current_g*__MACR_LATTICE_LY+i]=o_spin_config_int_ar_arg.bit_row_state[i];}

current_g+=dif;  E=o_spin_config_int_ar_arg.E;  was_overfilled=(was_overfilled | o_spin_config_int_ar_arg.was_overfilled);

return current_g;

}




unsigned long long int st_spin_config_int_ar :: add_row_state_nouv_rand(const st_spin_config_int_ar & o_spin_config_int_ar_arg)
{

  if (o_spin_config_int_ar_arg.current_g < 1) {return 0;}

unsigned long long int  dif                =treat_out_of_ullint_2_num_befor_add(current_g, o_spin_config_int_ar_arg.current_g, & was_overfilled);      //  insurance not to out of ULLINT type         //  o_spin_config_int_ar_arg.current_g;
unsigned long long int  dif_chunk_until_end=dif;
unsigned long long int  dif_rest_to_carry=0;

  if ((dif+current_g) > max_amount_of_spin_config) 
  {
  dif_chunk_until_end=max_amount_of_spin_config-current_g;   dif_rest_to_carry=dif-dif_chunk_until_end;  
    if (dif_rest_to_carry > max_amount_of_spin_config) {dif_rest_to_carry=max_amount_of_spin_config;}  
  was_overfilled=true;
  }

//std::memcpy(bit_row_state+current_g*sizeof(int)*__MACR_LATTICE_LY, o_spin_config_int_ar_arg.bit_row_state, dif*sizeof(int)*__MACR_LATTICE_LY);      //  !  !  !
//  current_cycle_out=0;  wished_cycle_out=0;

  for (unsigned long long int i=0;  i < dif_chunk_until_end*__MACR_LATTICE_LY;  i++)
  {bit_row_state[current_g*__MACR_LATTICE_LY+i]=o_spin_config_int_ar_arg.bit_row_state[i];} 

current_g+=dif_chunk_until_end;  E=o_spin_config_int_ar_arg.E;  was_overfilled=(was_overfilled | o_spin_config_int_ar_arg.was_overfilled);


//  here we begin random adding
//
  if (dif_rest_to_carry > 0)
  {
  unsigned long long int ullint_rand=(unsigned long long int) rand(); 
  //  we start on point:  0...or...(wished_rand_interval-1)
  unsigned long long int wished_rand_interval=max_amount_of_spin_config-dif_rest_to_carry+1;      //  ! ! ! DANGEROUS! FOR MEMORY LEAK!
   
    if (wished_rand_interval > RAND_MAX) {wished_rand_interval=RAND_MAX;} 
    
    if (wished_rand_interval < 2) 
    {
      if (ullint_rand % 2 == 0)
      {  
        for (unsigned long long int i=0;  i < dif_rest_to_carry*__MACR_LATTICE_LY;  i++)
        {bit_row_state[i]=o_spin_config_int_ar_arg.bit_row_state[i];}
      }    
    }
    else
    {
    unsigned long long int rand_shift_from_0=ullint_rand % wished_rand_interval;    //  it's key random param  ! ! !
      if (ullint_rand % 2 == 0)
      {  
        for (unsigned long long int i=0;  i < dif_rest_to_carry*__MACR_LATTICE_LY;  i++)
        {bit_row_state[rand_shift_from_0*__MACR_LATTICE_LY+i]=o_spin_config_int_ar_arg.bit_row_state[i];}
      }        
    }
    

  }

return current_g;

}




bool st_spin_config_int_ar :: convert_bit_row_state_to_bool_ar(unsigned long long int bit_row_state_number,  bool * const spin_array, int *Ly_arg,  int *Lx_arg)
{


  if (spin_array == 0) {return false;}
  if (bit_row_state_number >= current_g) {return false;}

int Ly=__MACR_LATTICE_LY;
int Lx=__MACR_LATTICE_LX;


  if (Ly_arg) {Ly=(*Ly_arg);}
  if (Lx_arg) {Lx=(*Lx_arg);}  


  //  #define __MACR_LATTICE_LY 13
  //  #define __MACR_LATTICE_LX 13
  //  bit_row_state[bit_row_state_number][col]  

unsigned int shift_one=0;  

  for (int row=0;  row < Ly;  row++)
  { 
  shift_one=1;

    for (int col=0;  col < Lx;  col++)
    {
    //__MACR_LATTICE_LY 
      if ((shift_one & bit_row_state[bit_row_state_number*Ly+row]) > 0) {spin_array[row*Lx+Lx-col-1]=true;} else {spin_array[row* Lx+Lx-col-1]=false;}
    shift_one=(shift_one << 1); 
    }

  }

return true;

}




void st_spin_config_int_ar :: write_to_screen(unsigned long long int bit_row_state_number,  char ch_spin_up, char ch_spin_dn, bool add_info, bool is_pbc1_fbc0)
{

bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
bool ok=convert_bit_row_state_to_bool_ar(bit_row_state_number,  spin_array);

  if (ok == true)
  {
    if (add_info == true)
    {
    string string_var; 
    string_var.append((char *)  "#  spin config number=");  add_ullint_to_string(& string_var, bit_row_state_number);
      if (is_pbc1_fbc0 == true) {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX);}
      else {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX+__MACR_LATTICE_LY+__MACR_LATTICE_LX);}
    string_var.append((char *)  "\t\tM= ");                 add_int_to_string(& string_var, calc_magnetiz_lat(spin_array,   __MACR_LATTICE_LY,   __MACR_LATTICE_LX));
    string_var.append((char *)  "\n");
    cout<<string_var.c_str();
    string_var.clear();
    }
  
  print_to_screen_2D_spin_lattice(spin_array,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);      //  !  !  !
  }


}




int st_spin_config_int_ar :: write_to_file(char * chstr_filename, unsigned long long int bit_row_state_number,  bool is_new,  char ch_spin_up, char ch_spin_dn, bool add_info, bool is_pbc1_fbc0)
{

bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
bool ok=convert_bit_row_state_to_bool_ar(bit_row_state_number,  spin_array);

  if (ok == true)
  {
  char *text=0;

  string string_var;
 
    if (add_info == true)
    {
    string_var.append((char *)  "#  spin config number=");  add_ullint_to_string(& string_var, bit_row_state_number);
      if (is_pbc1_fbc0 == true) {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX);}
      else {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX+__MACR_LATTICE_LY+__MACR_LATTICE_LX);}
    string_var.append((char *)  "\t\tM= ");                 add_int_to_string(& string_var, calc_magnetiz_lat(spin_array,   __MACR_LATTICE_LY,   __MACR_LATTICE_LX));
    string_var.append((char *)  "\n");  
    }
  
  int file_size=print_to_chstr_2D_spin_lattice(text, spin_array,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);       //  !  !  !
    if (text) {string_var.append(text);  file_size=write_chstr_to_file(string_var.c_str(), chstr_filename, is_new);    delete[] text;  text=0;}

  string_var.clear();
  
  return file_size;
  }

return 0;

}


void st_spin_config_int_ar :: write_to_screen_all_config(char ch_spin_up, char ch_spin_dn, bool add_info, bool is_pbc1_fbc0)
{

bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];

  for (int i=0;  i < current_g;  i++)
  {
  bool ok=convert_bit_row_state_to_bool_ar((unsigned long long int) i,  spin_array);

    if (ok == true)
    {
      if (add_info == true)
      {
      string string_var; 
      string_var.append((char *)  "#  spin config number=");  add_int_to_string(& string_var, i);
        if (is_pbc1_fbc0 == true) {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX);}
        else {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX+__MACR_LATTICE_LY+__MACR_LATTICE_LX);}
      string_var.append((char *)  "\t\tM= ");                 add_int_to_string(& string_var, calc_magnetiz_lat(spin_array,   __MACR_LATTICE_LY,   __MACR_LATTICE_LX));
      string_var.append((char *)  "\n");
      cout<<string_var.c_str();
      string_var.clear();
      }
  
    print_to_screen_2D_spin_lattice(spin_array,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);      //  !  !  !
    }
  }

}



int st_spin_config_int_ar :: write_to_file_all_config(char * chstr_filename, bool is_new,  char ch_spin_up,  char ch_spin_dn, bool add_info, bool is_pbc1_fbc0)
{

bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
int file_size=0;
string string_var;

  for (int i=0;  i < current_g;  i++)
  {
  bool ok=convert_bit_row_state_to_bool_ar((unsigned long long int) i,  spin_array);

    if (ok == true)
    {
    char *text=0;
 
      if (add_info == true)
      {
      string_var.append((char *)  "#  spin config number=");  add_ullint_to_string(& string_var, i);
        if (is_pbc1_fbc0 == true) {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX);}
        else {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX+__MACR_LATTICE_LY+__MACR_LATTICE_LX);}
      string_var.append((char *)  "\t\tM= ");                 add_int_to_string(& string_var, calc_magnetiz_lat(spin_array,   __MACR_LATTICE_LY,   __MACR_LATTICE_LX));
      string_var.append((char *)  "\n");  
      } 
  
   file_size=print_to_chstr_2D_spin_lattice(text, spin_array,  __MACR_LATTICE_LY,   __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);       //  !  !  !
      if (text) {string_var.append(text);   delete[] text;  text=0;}
  
    return file_size;
    }
  }


file_size=write_chstr_to_file(string_var.c_str(), chstr_filename, is_new);
string_var.clear();

return file_size;

}










void  st_spin_config_int_ar :: write_to_screen_with_J_bonds_one_config(unsigned long long int config_num, bool * J_ver, bool * J_hor, char ch_spin_up, char ch_spin_dn, bool add_info, char ch_ver_frame, char ch_hor_frame, bool is_pbc1_fbc0, bool is_do_frame, char bond_bkg_symb, bool is_check_E, int *Ly_arg,  int *Lx_arg)
{



int  Ly=__MACR_LATTICE_LY,  Lx=__MACR_LATTICE_LX;
  if (Ly_arg) {Ly=(*Ly_arg);}
  if (Lx_arg) {Lx=(*Lx_arg);}  


bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
int E_=0;

  if (config_num < current_g)
  {
  
  bool ok=convert_bit_row_state_to_bool_ar((unsigned long long int) config_num,  spin_array, &Ly, &Lx);

    if (ok == true)
    {
    // (1)  draw E, M  
    //
      if (add_info == true)
      {
        if (is_check_E == true)
        {
          if (is_pbc1_fbc0 == true) {calc_pospair_and_upspin_EA_2D_lattice_PBC_v2(spin_array,  J_ver,  J_hor,  Ly,  Lx,  & E_,  0);}
          else {calc_pospair_and_upspin_EA_2D_lattice_FBC_v2(spin_array,  J_ver,  J_hor,  Ly,  Lx,  & E_,  0);}
        }
      
      
      string string_var; 
      string_var.append((char *)  "#  spin config number=");  add_int_to_string(& string_var, config_num);
        if (is_pbc1_fbc0 == true) {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E_-2*Ly*Lx);}
        else {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var,    E_-((2*Ly*Lx-Ly-Lx)-E_)     );}        
      string_var.append((char *)  "\t\tM= ");                 add_int_to_string(& string_var, calc_magnetiz_lat(spin_array,   Ly,  Lx));
      string_var.append((char *)  "\n");
      cout<<string_var.c_str();
      string_var.clear();
      }

    // (2)  draw lattice  
    //
      if (is_pbc1_fbc0 == true) {print_to_screen_2D_EA_with_spin_state_lat_struct(spin_array, Ly,  Lx, J_ver, J_hor, ch_spin_up, ch_spin_dn, ch_ver_frame, ch_hor_frame);}      //  !  !  !
      else 
      {
        if ((J_ver) && (J_hor))
        {
        // TEMP --- v
        st_char_str_lattice_console_draw oblat;
        //oblat.set(Ly,  Ly,  5,  9,  3,  3);
        //oblat.set(Ly,  Ly,  3,  5,  1,  3);        
        //oblat.set(Ly,  Ly,  1,  3,  1,  1);                
        oblat.set(Ly,  Lx,  1,  3,  1,  1);        
        //oblat.set(__MACR_LATTICE_LY,  __MACR_LATTICE_LX,  3,  5,  1,  3);
        //oblat.lat_draw_type_and_print_to_screen_fbc();
        oblat.draw_fbc_J_bond_and_spin_lattice(J_ver, J_hor, spin_array, is_do_frame, bond_bkg_symb, ch_spin_up, ch_spin_dn);
        //oblat.draw_fbc_spin_lattice(spin_array);
        oblat.free(); 
        // TEMP --- ^
        }
        else
        {
        // TEMP --- v
        st_char_str_lattice_console_draw oblat;
        //oblat.set(Ly,  Ly,  5,  9,  3,  3);
        //oblat.set(Ly,  Ly,  3,  5,  1,  3);
        //oblat.set(Ly,  Ly,  1,  3,  1,  1);
        oblat.set(Ly,  Lx,  1,  3,  1,  1);                
        //oblat.set(__MACR_LATTICE_LY,  __MACR_LATTICE_LX,  3,  5,  1,  3);
        //oblat.lat_draw_type_and_print_to_screen_fbc();
        oblat.draw_fbc_spin_lattice(spin_array, is_do_frame, ch_spin_up, ch_spin_dn);
        //oblat.draw_fbc_spin_lattice(spin_array);
        oblat.free(); 
        // TEMP --- ^
        }

      //print_to_screen_2D_EA_with_spin_state_lat_struct(spin_array, __MACR_LATTICE_LY, __MACR_LATTICE_LX, J_ver, J_hor, ch_spin_up, ch_spin_dn, ch_ver_frame, ch_hor_frame);
      
      }
      
    }
  }

}





void st_spin_config_int_ar :: write_to_screen_with_J_bonds_all_config(bool * J_ver, bool * J_hor, char ch_spin_up, char ch_spin_dn, bool add_info, char ch_ver_frame, char ch_hor_frame, bool is_pbc1_fbc0)
{

bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];

  for (int i=0;  i < current_g;  i++)
  {
  bool ok=convert_bit_row_state_to_bool_ar((unsigned long long int) i,  spin_array);

    if (ok == true)
    {
      if (add_info == true)
      {
      string string_var; 
      string_var.append((char *)  "#  spin config number=");  add_int_to_string(& string_var, i);
        if (is_pbc1_fbc0 == true) {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX);}
        else {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX+__MACR_LATTICE_LY+__MACR_LATTICE_LX);}
      string_var.append((char *)  "\t\tM= ");                 add_int_to_string(& string_var, calc_magnetiz_lat(spin_array,   __MACR_LATTICE_LY,   __MACR_LATTICE_LX));
      string_var.append((char *)  "\n");
      cout<<string_var.c_str();
      string_var.clear();
      }

    print_to_screen_2D_EA_with_spin_state_lat_struct(spin_array, __MACR_LATTICE_LY, __MACR_LATTICE_LX, J_ver, J_hor, ch_spin_up, ch_spin_dn, ch_ver_frame, ch_hor_frame);      //  !  !  !
    }
  }

}





int st_spin_config_int_ar :: write_to_file_with_J_bonds_all_config(char * chstr_filename, bool is_new, bool * J_ver, bool * J_hor, char ch_spin_up,  char ch_spin_dn, bool add_info, char ch_ver_frame, char ch_hor_frame, bool is_pbc1_fbc0)
{

bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
int file_size=0;
string string_var;  

  for (int i=0;  i < current_g;  i++)
  {
  bool ok=convert_bit_row_state_to_bool_ar((unsigned long long int) i,  spin_array); 

    if (ok == true)
    {
    char *text=0;
 
      if (add_info == true)
      {
      string_var.append((char *)  "#  spin config number=");  add_ullint_to_string(& string_var, i);
        if (is_pbc1_fbc0 == true) {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX);}
        else {string_var.append((char *)  "\t\tE= ");  add_int_to_string(& string_var, 2*E-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX+__MACR_LATTICE_LY+__MACR_LATTICE_LX);}
      string_var.append((char *)  "\t\tM= ");                 add_int_to_string(& string_var, calc_magnetiz_lat(spin_array,   __MACR_LATTICE_LY,   __MACR_LATTICE_LX));
      string_var.append((char *)  "\n");  
      }
  
    file_size=print_to_chstr_2D_EA_with_spin_state_lat_struct(text, spin_array, __MACR_LATTICE_LY, __MACR_LATTICE_LX, J_ver, J_hor, ch_spin_up, ch_spin_dn, ch_ver_frame, ch_hor_frame);       //  !  !  !
      if (text) {string_var.append(text);   delete[] text;  text=0;}
  
    }
  }


file_size=write_chstr_to_file(string_var.c_str(), chstr_filename, is_new);
string_var.clear();

return file_size;

} 



bool st_spin_config_int_ar :: write_to_screen_with_J_bonds_all_config_with_E_check(bool * J_ver, bool * J_hor, char ch_spin_up,  char ch_spin_dn, bool add_info, char ch_ver_frame, char ch_hor_frame, bool is_pbc1_fbc0)
{

  if (current_g < 1) {return false;}
  
  for (unsigned long long int gs_state=0;  gs_state < current_g;  gs_state++)
  {
  write_to_screen_with_J_bonds_one_config(gs_state, J_ver, J_hor, ch_spin_up, ch_spin_dn, add_info, ch_ver_frame, ch_hor_frame, is_pbc1_fbc0, true, ' ', true);
  }

return true;

}





bool st_spin_config_int_ar :: check_energy(bool * J_ver, bool * J_hor, short int E_must_be, bool is_pbc1_fbc0, unsigned long long int start_state, unsigned long long int amount_of_states)
{

  if (current_g < 1) {return false;}
  if (J_ver == 0) {return false;}
  if (J_hor == 0) {return false;}
  
  
  if (current_g <= start_state) {return false;}
  if (start_state+amount_of_states > current_g) {amount_of_states=current_g-start_state;}
  
bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
bool ok_is_extracted_spin_state_into_bool_ar=false;    
  
  for (unsigned long long int gs_state=start_state;  gs_state < start_state+amount_of_states;  gs_state++)
  {
  ok_is_extracted_spin_state_into_bool_ar=convert_bit_row_state_to_bool_ar(gs_state,  spin_array); 
  
    if (ok_is_extracted_spin_state_into_bool_ar == false) {return false;}  
  
  int E_=0;  
  
    if (is_pbc1_fbc0 == true) {calc_pospair_and_upspin_EA_2D_lattice_PBC_v2(spin_array,  J_ver,  J_hor,  __MACR_LATTICE_LY,  __MACR_LATTICE_LX,  & E_,  0);  E_=2*E_-2*__MACR_LATTICE_LY*__MACR_LATTICE_LX;}      //  !  !  !
    else {calc_pospair_and_upspin_EA_2D_lattice_FBC_v2(spin_array,  J_ver,  J_hor,  __MACR_LATTICE_LY,  __MACR_LATTICE_LX,  & E_,  0);  E_=E_-((2*__MACR_LATTICE_LY*__MACR_LATTICE_LX-__MACR_LATTICE_LY-__MACR_LATTICE_LX)-E_);}      //  !  !  !
    
    if (E_must_be != E_) {return false;}

  }

return true;

}





bool st_spin_config_int_ar :: check_energy_edge_v3(bool * J_ver, bool * J_hor, int  amount_of_null_edge_bonds,  bool * bool_edge_spin_exist_d_u_l_r_ar,  short int E_must_be, unsigned long long int start_state, unsigned long long int amount_of_states, int *Ly_arg,  int *Lx_arg)
{

  if (current_g < 1) {return false;}
  if (J_ver == 0) {return false;}
  if (J_hor == 0) {return false;}
  
  
  if (current_g <= start_state) {return false;}
  if (start_state+amount_of_states > current_g) {amount_of_states=current_g-start_state;}



  
int  Ly=__MACR_LATTICE_LY,  Lx=__MACR_LATTICE_LX;
  if (Ly_arg) {Ly=(*Ly_arg);}
  if (Lx_arg) {Lx=(*Lx_arg);}  
  
bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX];

bool ok_is_extracted_spin_state_into_bool_ar=false;    
int  pospair_amount=0, E_=0;
  
  for (unsigned long long int gs_state=start_state;  gs_state < start_state+amount_of_states;  gs_state++)
  {
  ok_is_extracted_spin_state_into_bool_ar=convert_bit_row_state_to_bool_ar(gs_state,  spin_array, & Ly, & Lx); 
  
    if (ok_is_extracted_spin_state_into_bool_ar == false) {return false;}  

  //  void calc_pospair_and_upspin_EA_2D_lattice_with_edges_without_down_edge_v3(bool * spin_array,  bool * bool_edge_spin_exist_d_u_l_r_ar,  bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  int * pospair_amount,  int * up_amount,  bool calc_edge_down)
  E_=0;  pospair_amount=0;
  calc_pospair_and_upspin_EA_2D_lattice_with_edges_without_down_edge_v3(spin_array,  bool_edge_spin_exist_d_u_l_r_ar,  J_ver,  J_hor,  Ly,  Lx,  & pospair_amount,  0,  true);      //  !  !  !  
  E_=2*pospair_amount-(2*Ly*Lx-amount_of_null_edge_bonds-Ly-Lx+2*Ly+2*Lx);
  
  //==DEBUG==//  cout<<endl<<"DEBUG:  "<<"  pospair_amount="<<pospair_amount<<"E_must_be="<<E_must_be<<"   E_="<<E_<<endl;      //==DEBUG==//    
      
      
    if (E_must_be != E_) {return false;}  
  }

return true;

}




void st_spin_config_int_ar :: write_to_screen_2_state_cmp(unsigned long long int bit_row_state_number,  st_spin_config_int_ar & o_spin_config_int_ar_inst_2,  unsigned long long int bit_row_state_number_inst_2)
{

bool spin_array_1[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
bool ok=convert_bit_row_state_to_bool_ar(bit_row_state_number,  spin_array_1);
 
  //if (ok == true)
  //{
  //print_to_screen_2D_spin_lattice(spin_array_1,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);      //  !  !  !
  //}


bool spin_array_2[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
ok=o_spin_config_int_ar_inst_2.convert_bit_row_state_to_bool_ar(bit_row_state_number_inst_2,  spin_array_2);

  //if (ok == true)
  //{
  //print_to_screen_2D_spin_lattice(spin_array,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);      //  !  !  !
  //}


print_to_screen_2D_spin_state_comparing_lattice(spin_array_1, spin_array_2,  __MACR_LATTICE_LY, __MACR_LATTICE_LX);  

//void print_to_screen_2D_spin_state_comparing_lattice(const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx);

//int print_to_file_2D_spin_state_comparing_lattice(char *chstr_filename, const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx,  bool newfile);

//int print_to_chstr_2D_spin_state_comparing_lattice(char * & chstr_text, const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx);

  

}





bool  st_spin_config_int_ar :: write_to_screen_2_state_cmp_2(char *chstr_filename_J_bond,  unsigned long long int bit_row_state_number,  st_spin_config_int_ar & o_spin_config_int_ar_inst_2,  unsigned long long int bit_row_state_number_inst_2,   int d_len, int u_len, int l_len, int r_len,   int start_y, int start_x,   bool is_flip_bond_1ver_0hor, int flip_bond_index_y, int flip_bond_index_x)
{

  if (chstr_filename_J_bond == 0) {return false;}   

//--DEBU--//  cout<<" d_len="<<d_len<<" u_len="<<u_len<<" l_len="<<l_len<<" r_len="<<r_len<<" start_y="<<start_y<<" start_x="<<start_x<<endl;      //--DEBU--//

      
bool * J_ver=0,  * J_hor=0;
int Ly_=0, Lx_=0;  
read_J_bonds_from_file(chstr_filename_J_bond, J_ver,  J_hor,  Ly_,  Lx_, false);
  if ((Ly_ < 2) || (Lx_ < 2) || (J_ver == 0) || (J_hor == 0)) {return false;} 




bool spin_array_1[Ly_*Lx_];
bool ok=convert_bit_row_state_to_bool_ar(bit_row_state_number,  spin_array_1);
 
  //if (ok == true)
  //{
  //print_to_screen_2D_spin_lattice(spin_array_1,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);      //  !  !  !
  //}


bool spin_array_2[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
ok=o_spin_config_int_ar_inst_2.convert_bit_row_state_to_bool_ar(bit_row_state_number_inst_2,  spin_array_2);

  //if (ok == true)
  //{
  //print_to_screen_2D_spin_lattice(spin_array,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);      //  !  !  !
  //}


    
st_char_str_lattice_console_draw oblat_1;
//oblat.set(Ly,  Ly,  5,  9,  3,  3);
oblat_1.set(Ly_,  Ly_,  3,  5,  1,  3);
//oblat_1.set(Ly_,  Ly_,  1,  3,  1,  1);
  
//oblat.lat_draw_type_and_print_to_screen_fbc();
//oblat.draw_fbc_J_bond_lattice(J_ver, J_hor);
//cout<<"g  "<<flip_index_y<<"  "<<flip_index_x<<"  ver=";  display_bool(flip_1ver_0hor);  cout<<endl;
//oblat.draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond(J_ver, J_hor, 0, flip_1ver_0hor, flip_index_y, flip_index_x, true, ' ', 'u', 'd');
//--DEBUG--//  cout<<"config 1 - initial"<<endl;
//  proc prototype
//  void draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame(bool *J_ver, bool *J_hor, bool *spin_array, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x, bool is_do_out_frame=true, char bond_bkg_symb=' ', char spin_up='u', char spin_down='d', char inner_frame_symb=':');
//

  if (is_flip_bond_1ver_0hor == true ) {J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]=(!J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]);}
  if (is_flip_bond_1ver_0hor == false) {J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]=(!J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]);}
  
oblat_1.draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame(J_ver, J_hor, spin_array_1, is_flip_bond_1ver_0hor, flip_bond_index_y, flip_bond_index_x, d_len, u_len, l_len, r_len, start_y, start_x, true,  ' ' , 'u', 'd');
  
cout<<endl;
oblat_1.free();

st_char_str_lattice_console_draw oblat_2;
//oblat.set(Ly,  Ly,  5,  9,  3,  3);
oblat_2.set(Ly_,  Ly_,  3,  5,  1,  3);
//oblat_2.set(Ly_,  Ly_,  1,  3,  1,  1);  
//--DEBUG--//  cout<<"config 2 - flipped J"<<endl;
//  proc prototype   
//  void draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame_cmp(bool *J_ver, bool *J_hor, bool *spin_array_to_cmp, bool *spin_array_cur, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x, bool is_do_out_frame, char bond_bkg_symb, char spin_up_same='u', char spin_down_same='d', char spin_up_diff='U', char spin_down_diff='D', char inner_frame_symb=':')
//

  if (is_flip_bond_1ver_0hor == true ) {J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]=(!J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]);}
  if (is_flip_bond_1ver_0hor == false) {J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]=(!J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]);}

oblat_2.draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame_cmp(J_ver, J_hor, spin_array_1, spin_array_2, is_flip_bond_1ver_0hor, flip_bond_index_y, flip_bond_index_x, d_len, u_len, l_len, r_len, start_y, start_x, true,  ' ' , 'u', 'd' , 'U', 'D');
cout<<endl;
oblat_2.free();  
    
  if (J_ver) {delete[] J_ver;  J_ver=0;}
  if (J_hor) {delete[] J_hor;  J_hor=0;}
    
return true;   
  
}






void st_spin_config_int_ar :: write_to_screen_with_J_bonds_config_edge_v3(bool * J_ver, bool * J_hor, int config_num, int amount_of_null_edge_bonds, bool *bool_edge_spin_exist_d_u_l_r_ar, bool check_E, char ch_spin_up, char ch_spin_dn, bool add_info, char ch_ver_frame, char ch_hor_frame, int *Ly_arg,  int *Lx_arg)
{

bool spin_array[__MACR_LATTICE_LY*__MACR_LATTICE_LX]; 

  //for (int i=0;  i < current_g;  i++)

int  Ly=__MACR_LATTICE_LY,  Lx=__MACR_LATTICE_LX;
  if (Ly_arg) {Ly=(*Ly_arg);}
  if (Lx_arg) {Lx=(*Lx_arg);}  
unsigned long long int i=config_num;  
  
  {
  bool ok=convert_bit_row_state_to_bool_ar((unsigned long long int) i,  spin_array, &Ly, &Lx);

    if (ok == true)
    {
    
      if (add_info == true)   
      {
      int pospair_amount_local=0;  
      string string_var; 
      string_var.append((char *)  "#  spin config number=");  add_int_to_string(& string_var, i);
      
        if (check_E == true)
        {
        calc_pospair_and_upspin_EA_2D_lattice_with_edges_without_down_edge_v3(spin_array,  bool_edge_spin_exist_d_u_l_r_ar, J_ver, J_hor, Ly, Lx, & pospair_amount_local, 0, true);
        string_var.append((char *)  "\t\tfact E= ");
        add_int_to_string(& string_var,   2*pospair_amount_local-(2*Ly*Lx-Ly-Lx+2*Ly+2*Lx-amount_of_null_edge_bonds));  
        }
        else
        {
        string_var.append((char *)  "\t\tsupposing E= ");       add_int_to_string(& string_var, 2*E_pospair-(2*Ly*Lx-Ly-Lx+2*Ly+2*Lx-amount_of_null_edge_bonds));
        }
        
      string_var.append((char *)  "\t\tM= ");                 add_int_to_string(& string_var, calc_magnetiz_lat(spin_array,   Ly,   Lx));
      
      string_var.append((char *)  "\n");
      cout<<string_var.c_str();
      string_var.clear();
      }

    st_char_str_lattice_console_draw oblat_1;
    oblat_1.set(Ly+2,  Lx+2,  3,  5,  1,  3);
    oblat_1.draw_fbc_J_bond_and_spin_lattice_frame_proj_v3(J_ver, J_hor, spin_array,  bool_edge_spin_exist_d_u_l_r_ar, true, ' ', 'u', 'd');
    oblat_1.free();
    }
  }




//  draw_fbc_J_bond_and_spin_lattice_frame_proj_v3(J_ver, J_hor, spin_array_, bool_edge_spin_exist_d_u_l_r_ar, true, ' ', 'u', 'd');
//int amount_of_positiv_E_par_revers_order_edge_v3(bool *J_line,  bool * bool_edge_spin_exist_d_u_l_r_ar, unsigned long long int spin_line_state, int len);
//void calc_pospair_and_upspin_EA_2D_lattice_with_edges_without_down_edge_v3(bool * spin_array,  bool * bool_edge_spin_exist_d_u_l_r_ar, bool * J_ver, bool * J_hor, int Ly, int Lx, int * pospair_amount=0, int * up_amount=0, bool calc_edge_down=false);


}






bool st_spin_config_int_ar :: write_to_screen_2_state_edge_cmp_2_v3( st_spin_config_int_ar * o_spin_config_int_ar_inst_to_cmp,  unsigned long long int number_of_gs_1,  unsigned long long int number_of_gs_2_to_cmp, bool * J_ver,  bool * J_hor,  bool * array_spin_exist_d_u_l_r,   bool is_flip_bond_1ver_0hor, int flip_bond_index_y, int flip_bond_index_x, bool is_draw_domain_wall)
{
   
  if (J_ver == 0) {return false;}
  if (J_hor == 0) {return false;}
  if (o_spin_config_int_ar_inst_to_cmp->current_g < 1) {return false;} 
  if (o_spin_config_int_ar_inst_to_cmp->current_g <= number_of_gs_2_to_cmp) {return false;} 
  if (this->current_g < 1) {return false;} 
  if (this->current_g <= number_of_gs_1) {return false;} 
    
  
//--DEBUG--//  cout<<" d_len="<<d_len<<" u_len="<<u_len<<" l_len="<<l_len<<" r_len="<<r_len<<" start_y="<<start_y<<" start_x="<<start_x<<endl;      //--DEBUG--//

      
int Ly_=__MACR_LATTICE_LY, Lx_=__MACR_LATTICE_LX;  
  if ((Ly_ < 2) || (Lx_ < 2) || (J_ver == 0) || (J_hor == 0)) {return false;} 




bool spin_array_1[Ly_*Lx_];
bool ok=convert_bit_row_state_to_bool_ar(number_of_gs_1,  spin_array_1);
 
  //if (ok == true)
  //{
  //print_to_screen_2D_spin_lattice(spin_array_1,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);      //  !  !  !
  //}


bool spin_array_2[__MACR_LATTICE_LY*__MACR_LATTICE_LX];
ok=o_spin_config_int_ar_inst_to_cmp->convert_bit_row_state_to_bool_ar(number_of_gs_2_to_cmp,  spin_array_2);

  //if (ok == true)
  //{
  //print_to_screen_2D_spin_lattice(spin_array,  __MACR_LATTICE_LY, __MACR_LATTICE_LX, ch_spin_up,  ch_spin_dn);      //  !  !  !
  //}


    
st_char_str_lattice_console_draw oblat_1;
//oblat.set(Ly,  Ly,  5,  9,  3,  3);
oblat_1.set(Ly_+2,  Lx_+2,  3,  5,  1,  3);
//oblat_1.set(Ly_,  Ly_,  1,  3,  1,  1);
  
//oblat.lat_draw_type_and_print_to_screen_fbc();
//oblat.draw_fbc_J_bond_lattice(J_ver, J_hor);
//cout<<"g  "<<flip_index_y<<"  "<<flip_index_x<<"  ver=";  display_bool(flip_1ver_0hor);  cout<<endl;
//oblat.draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond(J_ver, J_hor, 0, flip_1ver_0hor, flip_index_y, flip_index_x, true, ' ', 'u', 'd');
//--DEBUG--//  cout<<"config 1 - initial"<<endl;
//  proc prototype
//  void draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame(bool *J_ver, bool *J_hor, bool *spin_array, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x, bool is_do_out_frame=true, char bond_bkg_symb=' ', char spin_up='u', char spin_down='d', char inner_frame_symb=':');
//

  if (is_flip_bond_1ver_0hor == true ) {J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]=(!J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]);}
  if (is_flip_bond_1ver_0hor == false) {J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]=(!J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]);} 

oblat_1.draw_fbc_J_bond_and_spin_lattice_frame_proj_v3(J_ver, J_hor, spin_array_1, array_spin_exist_d_u_l_r, true, ' ', 'u', 'd');  

  
cout<<endl;
oblat_1.free();

st_char_str_lattice_console_draw oblat_2;
//oblat.set(Ly,  Ly,  5,  9,  3,  3);
oblat_2.set(Ly_+2,  Lx_+2,  3,  5,  1,  3);
//oblat_2.set(Ly_,  Ly_,  1,  3,  1,  1);  
//--DEBUG--//  cout<<"config 2 - flipped J"<<endl;
//  proc prototype   
//void st_char_str_lattice_console_draw :: draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame_cmp_frame_proj_v3(bool *J_ver, bool *J_hor, bool *spin_array_to_cmp_, bool *spin_array_cur_, bool * array_spin_exist_d_u_l_r, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int start_y,  int start_x, bool is_do_out_frame, char bond_bkg_symb, char spin_up_same='u', char spin_down_same='d', char spin_up_diff='U', char spin_down_diff='D', char inner_frame_symb=':',  char domain_wall_symb='.', bool is_draw_domain_wall=true)
//

  if (is_flip_bond_1ver_0hor == true ) {J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]=(!J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]);}
  if (is_flip_bond_1ver_0hor == false) {J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]=(!J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]);}

oblat_2.draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame_cmp_frame_proj_v3(J_ver, J_hor, spin_array_1, spin_array_2, array_spin_exist_d_u_l_r, is_flip_bond_1ver_0hor, flip_bond_index_y, flip_bond_index_x,  0, 0, true, ' ', 'u', 'd',  'U', 'D', ':',  '.', is_draw_domain_wall);
cout<<endl;
oblat_2.free();  

  if (is_flip_bond_1ver_0hor == true ) {J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]=(!J_ver[flip_bond_index_y*Lx_+flip_bond_index_x]);}
  if (is_flip_bond_1ver_0hor == false) {J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]=(!J_hor[flip_bond_index_y*(Lx_+1)+flip_bond_index_x]);}
    
    
return true;   

}





bool extract_J_bond_submatr_from_ver(bool * J_bond_ver_src, int Ly_src, int Lx_src, bool *& J_bond_ver_dest, int spin_start_pos_y_src,  int spin_start_pos_x_src,   int spin_dy,  int spin_dx)
{

  if (J_bond_ver_dest == 0) {J_bond_ver_dest=new bool[(spin_dy+1)*spin_dx];}
    
  
int index_y_src=0,  index_x_src=0;
  
  for (int y=0;  y < spin_dy+1;  y++)
  {
  index_y_src=spin_start_pos_y_src+y;  if (Ly_src < index_y_src) {index_y_src=index_y_src-Ly_src-1;}
  
    for (int x=0;  x < spin_dx;  x++)
    {
    index_x_src=spin_start_pos_x_src+x;  if (Lx_src <= index_x_src) {index_x_src=index_x_src-Lx_src;}
    J_bond_ver_dest[y*spin_dx+x]=J_bond_ver_src[index_y_src*Lx_src+index_x_src];
    }  
    
  }

return true;

}




bool extract_J_bond_submatr_from_hor(bool * J_bond_hor_src, int Ly_src, int Lx_src, bool *& J_bond_hor_dest, int spin_start_pos_y_src,  int spin_start_pos_x_src,   int spin_dy,  int spin_dx)
{

  if (J_bond_hor_dest) {J_bond_hor_dest=new bool[spin_dy*(spin_dx+1)];}
  
    
int index_y_src=0,  index_x_src=0;
  
  for (int y=0;  y < spin_dy;  y++)
  {
  index_y_src=spin_start_pos_y_src+y;  if (Ly_src <= index_y_src) {index_y_src=index_y_src-Ly_src;}
  
    for (int x=0;  x < spin_dx+1;  x++)
    {
    index_x_src=spin_start_pos_x_src+x;  if (Lx_src < index_x_src) {index_x_src=index_x_src-Lx_src-1;}
    J_bond_hor_dest[y*(spin_dx+1)+x]=J_bond_hor_src[index_y_src*(Lx_src+1)+index_x_src];
    }  
    
  }

return true;

}





//  spin_dy=3   int spin_dx=3
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

// (4)
//
bool extract_environment_4(bool *spin_array_src, int Ly, int Lx,  int down_y, int left_x, int dy, int dx,  bool *& bool_edge_spin_exist_d_u_l_r_ar, bool pbc_1_fbc_0)
{

  if (spin_array_src == 0) {return false;}
  if ((pbc_1_fbc_0 == false) && (left_x+dx > Lx)) {return false;}
  if ((pbc_1_fbc_0 == false) && (down_y+dy > Ly)) {return false;}
  if ((pbc_1_fbc_0 == false) && (left_x < 0)) {return false;}
  if ((pbc_1_fbc_0 == false) && (down_y < 0)) {return false;}


  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  bool_edge_spin_exist_d_u_l_r_ar=new bool[dy*4+dx*4];
  }
  
int y_src=0, x_src=0;
    
  if (pbc_1_fbc_0 == true)
  {
  //  pbc
  
    for (int x=0;  x < dx;  x++)
    {
    //  down-line 
    //    
    y_src=down_y-1;  cycl_correction(y_src, Ly);    
    x_src=left_x+x;  cycl_correction(x_src, Lx);
    bool_edge_spin_exist_d_u_l_r_ar[2*x]=spin_array_src[y_src*Lx+x_src];  bool_edge_spin_exist_d_u_l_r_ar[2*x+1]=true;

    //  up-line 
    //    
    y_src=down_y+dy; cycl_correction(y_src, Ly);    
    x_src=left_x+x;  cycl_correction(x_src, Lx);
    bool_edge_spin_exist_d_u_l_r_ar[2*dx+2*x]=spin_array_src[y_src*Lx+x_src];  bool_edge_spin_exist_d_u_l_r_ar[2*dx+2*x+1]=true;    
    }
      
    for (int y=0;  y < dy;  y++)
    {
    //  left-line 
    //        
    y_src=down_y+y;  cycl_correction(y_src, Ly);    
    x_src=left_x-1;  cycl_correction(x_src, Lx);
    bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y]=spin_array_src[y_src*Lx+x_src];  bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y]=true;

    //  right-line 
    //          
    y_src=down_y+y;     cycl_correction(y_src, Ly);    
    x_src=left_x+dx;    cycl_correction(x_src, Lx);
    bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*dy+2*y]=spin_array_src[y_src*Lx+x_src];  bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*dy+2*y+1]=true;    
    }
  
  }
  else
  {
  //  fbc
  
    for (int x=0;  x < dx;  x++)
    {
    //  down-line 
    //        
    y_src=down_y-1;   
    x_src=left_x+x;  
      if (down_y-1 > -1) {bool_edge_spin_exist_d_u_l_r_ar[2*x]=spin_array_src[y_src*Lx+x_src];  bool_edge_spin_exist_d_u_l_r_ar[2*x+1]=true;} 
      else {bool_edge_spin_exist_d_u_l_r_ar[2*x]=false;  bool_edge_spin_exist_d_u_l_r_ar[2*x+1]=false;}
    
    //  up-line 
    //        
    y_src=down_y+dy;     
    x_src=left_x+x;  
      if (down_y+dy < Ly) {bool_edge_spin_exist_d_u_l_r_ar[2*dx+2*x]=spin_array_src[y_src*Lx+x_src];  bool_edge_spin_exist_d_u_l_r_ar[2*dx+2*x+1]=true;} 
      else {bool_edge_spin_exist_d_u_l_r_ar[2*dx+2*x]=false;  bool_edge_spin_exist_d_u_l_r_ar[2*dx+2*x+1]=false;}    
    }
      
    for (int y=0;  y < dy;  y++)
    {
    //  left-line 
    //        
    y_src=down_y+y;  
    x_src=left_x-1;
      if (left_x-1 > -1) {bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y]=spin_array_src[y_src*Lx+x_src]; bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y+1]=true;} 
      else {bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y]=false; bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y+1]=false;}

    //  right-line 
    //        
    y_src=down_y+y;     
    x_src=left_x+dx;
      if (left_x+dx < Lx) {bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*dy+2*y]=spin_array_src[y_src*Lx+x_src];  bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*dy+2*y+1]=true;} 
      else {bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*dy+2*y]=false;  bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*dy+2*y+1]=false;}
        
    }
  
  }
  
return true;
  
}



// (5)
//
bool rand_generate_edges_spin_val(bool * bool_edge_spin_exist_d_u_l_r_ar, int dy, int dx)
{

  if (bool_edge_spin_exist_d_u_l_r_ar == 0) {return false;}
  
//  down-line and up-line 
//        

  for (int x=0;  x < 2*dx;  x++)
  {  if (bool_edge_spin_exist_d_u_l_r_ar[2*x+1] == true) {bool_edge_spin_exist_d_u_l_r_ar[2*x]=((bool) (rand() % 2));} else {bool_edge_spin_exist_d_u_l_r_ar[2*x]=false;}  }


//  left-line and right-line 
//        

  for (int y=0;  y < 2*dy;  y++)
  {  if (bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y+1] == true) {bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y]=((bool) (rand() % 2));} else {bool_edge_spin_exist_d_u_l_r_ar[4*dx+2*y]=false;}  }


return true;
  
}



// (6)
//
//
//   example:  dy=4,  dx=3          
//   down: spin_array_src[0]  spin_array_src[1]  spin_array_src[2]          up: spin_array_src[3]  spin_array_src[4]  spin_array_src[5]            left: spin_array_src[6]  spin_array_src[7]            right: spin_array_src[8]  spin_array_src[9]  
//
//
//   -------------------------
//   | *         *         * |
//   -------------------------
//   -----               -----
//   |   |               |   |
//   | * |               | * |
//   |   |               |   |
//   |   |               |   |
//   | * |               | * |
//   |   |               |   |
//   -----               ----- 
//   -------------------------
//   | *         *         * |
//   -------------------------

bool extract_frame_spin_value(bool *spin_array_src, int Ly, int Lx,  bool *& spin_frame_value_4, int down_y, int left_x, int dy, int dx, bool pbc_1_fbc_0)
{

  if (spin_frame_value_4) {return false;}
  if (spin_frame_value_4 == 0) {return false;}
  if ((pbc_1_fbc_0 == false) && (left_x+dx >= Lx)) {return false;}
  if ((pbc_1_fbc_0 == false) && (down_y+dy >= Ly)) {return false;}
  if ((pbc_1_fbc_0 == false) && (left_x <= 0)) {return false;}
  if ((pbc_1_fbc_0 == false) && (down_y <= 0)) {return false;}



spin_frame_value_4=new bool[4*dy+4*(dx-2)];
int y_src=0, x_src=0;


  for (int x=0;  x < dx;  x++)
  {
  //  down-line 
  //    
  y_src=down_y;    cycl_correction(y_src, Ly);    
  x_src=left_x+x;  cycl_correction(x_src, Lx);
  spin_frame_value_4[x]=spin_array_src[y_src*Lx+x_src]; 

  //  up-line 
  //    
  y_src=down_y+dy-1;  cycl_correction(y_src, Ly);    
  x_src=left_x+x;     cycl_correction(x_src, Lx);
  spin_frame_value_4[x+dx]=spin_array_src[y_src*Lx+x_src];     
  }
      
  for (int y=0;  y < dy;  y++)
  {
  //  left-line 
  //        
  y_src=down_y+y;  cycl_correction(y_src, Ly);    
  x_src=left_x;    cycl_correction(x_src, Lx);
  spin_frame_value_4[2*dx+y]=spin_array_src[y_src*Lx+x_src];

  //  right-line 
  //          
  y_src=down_y+y;       cycl_correction(y_src, Ly);    
  x_src=left_x+dx-1;    cycl_correction(x_src, Lx);
  spin_frame_value_4[2*dx+dy+y]=spin_array_src[y_src*Lx+x_src];    
  }
  

return true;

}




bool paste_spin_sublat_into_spin_lat(bool * spin_array_dest,  int Ly_dest,  int Lx_dest,  bool * spin_array_src, int dest_down_y, int dest_left_x, int dy, int dx, bool pbc_1_fbc_0)
{

  if (spin_array_dest == 0) {return false;}
  if (spin_array_src  == 0) {return false;}  
  
int y_dest=0, x_dest=0;

  for (int src_y=0;  src_y < dy;  src_y++)
  {
  y_dest=cycl_correction(dest_down_y+src_y, Ly_dest);

    for (int src_x=0;  src_x < dx;  src_x++)
    {
    x_dest=cycl_correction(dest_left_x+src_x, Lx_dest);
    spin_array_dest[y_dest*Lx_dest+x_dest]=spin_array_src[src_y*dx+src_x];
    }
    
  }  
  
return true;  

}






bool calc_amount_of_null_edge_bonds_edge_v3(int Ly,  int Lx,  bool * bool_edge_spin_exist_d_u_l_r_ar)
{

int amount_of_null_edge_bonds=0;      //  calc here

  {
    for (int y=0;  y < Ly;  y++)
    {  if (bool_edge_spin_exist_d_u_l_r_ar[Lx*4+2*y+1] == true)  {amount_of_null_edge_bonds++;}    if (bool_edge_spin_exist_d_u_l_r_ar[Lx*4+Ly*2+2*y+1] == true)  {amount_of_null_edge_bonds++;}  }

    for (int x=0;  x < Lx;  x++)
    {  if (bool_edge_spin_exist_d_u_l_r_ar[2*x+1] == true)  {amount_of_null_edge_bonds++;}         if (bool_edge_spin_exist_d_u_l_r_ar[Lx*2+2*x+1] == true)  {amount_of_null_edge_bonds++;}  }          
  }
  
return amount_of_null_edge_bonds;  

}





bool print_struct_with_J_edges_v3(bool *J_ver, bool *J_hor, int Ly,  int Lx,  bool * bool_edge_spin_exist_d_u_l_r_ar)
{

st_char_str_lattice_console_draw oblat_1;
oblat_1.set(Ly+2,  Lx+2,  3,  5,  1,  3);
oblat_1.draw_fbc_J_bond_and_spin_lattice_frame_proj_v3(J_ver, J_hor, 0,  bool_edge_spin_exist_d_u_l_r_ar, true, ' ', 'u', 'd');
oblat_1.free();

return true;

}





int match_test_for_0_radius_struct(int N)
{


int matc_cases_amount=0;
int E_1[4], E_2[4];


  for (int n=0;  n < N;  n++)
  {
  int min_E_1=100,  min_E_2=100;
  int min_E_counter_1=0, min_E_counter_2=0;
  int min_E_config_1[4], min_E_config_2[4];
  
  //  gen J
  //
  int J_ver[3], J_hor[4];

    for (int i=0;  i < 3;  i++)  {J_ver[i]=-1;}
    for (int i=0;  i < 4;  i++)  {J_hor[i]=-1;}  
    if ((rand() % 2) == 0) {J_ver[1]=+1;}
  
  
  //  bool_edge_spin_exist_d_u_l_r_ar
  //
  int edge_spins[6];

    for (int i=0;  i < 6;  i++)
    {edge_spins[i]=-1;}


    for (int k_ran=0;  k_ran < 3;  k_ran++)
    {
    int pos=0;
  
      for (int k=0;  k < 1000;  k++)
      {
      pos=rand() % 6;
        if (edge_spins[pos] == -1) {edge_spins[pos]=+1; break;}
      }
    }
    
  
  
  int spin_down=-1,  spin_up=-1;
  
    for (int config=0;  config < 4;  config++)
    {
    
      if (config == 0) {spin_down=-1;  spin_up=-1;}    
      if (config == 1) {spin_down=+1;  spin_up=-1;}
      if (config == 2) {spin_down=-1;  spin_up=+1;}
      if (config == 3) {spin_down=+1;  spin_up=+1;}      

    E_1[config]=spin_down*spin_up*J_ver[1]  +  spin_down*(edge_spins[0]*J_ver[0]+edge_spins[2]*J_hor[0]+edge_spins[3]*J_hor[1])  +  spin_up*(edge_spins[1]*J_ver[2]+edge_spins[4]*J_hor[2]+edge_spins[5]*J_hor[3]);    
    E_2[config]=E_1[config]-2*spin_down*spin_up*J_ver[1];
    
      if (E_1[config] <= min_E_1) {min_E_1=E_1[config];}
      if (E_2[config] <= min_E_2) {min_E_2=E_2[config];}      
      
    }
    
    for (int config=0;  config < 4;  config++)
    {
      if (E_1[config] == min_E_1) {min_E_config_1[min_E_counter_1]=config;  min_E_counter_1++;}
      if (E_2[config] == min_E_2) {min_E_config_2[min_E_counter_2]=config;  min_E_counter_2++;}      
    }
    

  int choiced_rand_config_1=rand() % min_E_counter_1;
  
    if (E_2[choiced_rand_config_1] == min_E_2) {matc_cases_amount++;}
    
    
  }      //  for (int n=0;  n < N;  n++)



return matc_cases_amount;

}




bool print_bool_edge_spin_exist_d_u_l_r_ar_for_debug(int Ly,  int Lx,  bool * bool_edge_spin_exist_d_u_l_r_ar)
{

  if (bool_edge_spin_exist_d_u_l_r_ar == 0) {return false;}


  
cout<<"down   bool_edge_spin_exist_d_u_l_r_ar:"<<endl;
cout<<":::"<<endl;

  for (int x=0;  x < Lx;  x++)
  {
    if (bool_edge_spin_exist_d_u_l_r_ar[2*x+1] == true) {cout<<"1";} else {cout<<"0";}    
  }
cout<<endl;
  for (int x=0;  x < Lx;  x++)
  {
    if (bool_edge_spin_exist_d_u_l_r_ar[2*x  ] == true) {cout<<"1";} else {cout<<"0";}    
  }
cout<<endl;

cout<<":::"<<endl;
cout<<"up   bool_edge_spin_exist_d_u_l_r_ar:"<<endl;
cout<<":::"<<endl;

  for (int x=0;  x < Lx;  x++)
  {
    if (bool_edge_spin_exist_d_u_l_r_ar[2*Lx+2*x+1] == true) {cout<<"1";} else {cout<<"0";}    
  }
cout<<endl;
  for (int x=0;  x < Lx;  x++)
  {
    if (bool_edge_spin_exist_d_u_l_r_ar[2*Lx+2*x  ] == true) {cout<<"1";} else {cout<<"0";}    
  }
cout<<endl;


cout<<":::"<<endl;  
cout<<"left   bool_edge_spin_exist_d_u_l_r_ar:"<<endl;
cout<<":::"<<endl;

  for (int y=0;  y < Ly;  y++)
  {
    if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*y+1] == true) {cout<<"1";} else {cout<<"0";}    
  }
cout<<endl;
  for (int y=0;  y < Ly;  y++)
  {
    if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*y  ] == true) {cout<<"1";} else {cout<<"0";}    
  }
cout<<endl;

cout<<":::"<<endl;
cout<<"right   bool_edge_spin_exist_d_u_l_r_ar:"<<endl;
cout<<":::"<<endl;

  for (int y=0;  y < Ly;  y++)
  {
    if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*Ly+2*y+1] == true) {cout<<"1";} else {cout<<"0";}    
  }
cout<<endl;
  for (int y=0;  y < Ly;  y++)
  {
    if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*Ly+2*y  ] == true) {cout<<"1";} else {cout<<"0";}    
  }
cout<<endl;


return true;

}




/*
void write_ullint_spin_sublat_state_to_big_lat(unsigned long long int spin_sublattice_state_code, bool * spin_array,  int Lx,  int Ly,  int sub_left_down_ind_x,  int sub_left_down_ind_y, int sub_len_x, int sub_len_y)
{

unsigned long long int staying_1_bit1=1;
unsigned long long int moving_1_bit1=1;
unsigned long long int one_spin_state=0;


int index_y=sub_left_down_ind_y,  index_x=sub_left_down_ind_x;


  for (int index_y=0;  index_y < sub_len_y;  index_y++)
  {
  index_cycl_shift(index_y, Ly);
  index_x=sub_left_down_ind_x;

    for (int index_x=0;  index_x < sub_len_x;  index_x++)
    {
    index_cycl_shift(index_x, Lx);
    one_spin_state=(spin_sublattice_state_code & moving_1_bit1);

      if (one_spin_state != 0) {spin_array[Lx*index_y+index_x]=true;} else {spin_array[Lx*index_y+index_x]=false;}

    moving_1_bit1=(moving_1_bit1 << 1);    
    index_x++;
    }

  index_y++;
  }


   *   *
   |   |
-- *-- *-- 
   |   |
-- *-- *-- 
   |   |
   *   *

64 random bit
((long long)rand() << 32) | rand()

}*/




/*




#include <iostream>
#include <random>
#include <cmath>

int main()
{
    std::random_device rd;

    std::mt19937_64 e2(rd());

    std::uniform_int_distribution<long long int> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));

    std::cout << std::llround(std::pow(2,61)) << std::endl; 
    std::cout << std::llround(std::pow(2,62)) << std::endl; 

    for (int n = 0; n < 10; ++n) {
            std::cout << dist(e2)<< ", " ;
    }
    std::cout << std::endl ;
}




*/




/*

   *   *   *
   |   |   |
---*---*---*---
   |   |   |
---*---*---*---
   |   |   |
---*---*---*---
   |   |   |
   *   *   *


   *   *   *
   |   |   |
---6---7---8---
   |   |   |
---3---4---5---
   |   |   |
---0---1---2---
   |   |   |
   *   *   *


J_ver

   *   *   *
   0   10  11
---*---*---*---
   6   7   8
---*---*---*---
   3   4   5
---*---*---*---
   0   1   2
   *   *   *

J_hor

   *   *   *
   |   |   |
-8-*-9-*-10*-11
   |   |   |
-4-*-5-*-6-*-7-
   |   |   |
-0-*-1-*-2-*-3-
   |   |   |
   *   *   *

*/









                                                                                                                                                                                     
                                                                                                                                                                                                        
                                                                                                                                                                                                        
                                                                                                                                                                                                        
//                                                                                                                                                                                                        
//                                                                                                                                                                                                        
//                                                                                                                                                                                                        
//        .,..                                                                                                                                                                                            
//       .*,&%%%&@@&%#/*,..                                                                                                                                                                               
//     ...  ##&@@&&&%&&&%%&&%%%%%%#%%/                                                                                                                                                                    
//      ./##(*.,*(%@@@(,(%%%(%%%%%%%%%%&%/                                                                                                                                                                
//        ..      .,*(&@@&&%%%%%%%###%%%%%%&%*                                                        *                                                                                                   
//                    .,*(&@&&&%&%%%%#%##%#%%%%&&/.            .@@%,.                               .&@%..                                                                                                
//                        .,/(&@@@&&%%%&&#%%%%####%%%%&#%/   #@@@@@%.                             .(#%%%@*                                                                                                
//                            .,*/((##%&@@@@@&&&%%%%%%%#&@@@@@@&&@@@@&..                          &%##@/#&@(        &@.                                                                                   
//                                   ...,,**///((#&&%%%@@@@&(&@@@@@@@&(.                          &%%**(%&&%/.   #%%@@**.                                                                                 
//                                                .,*%%%&@@@@@@@@&#&&&@@@@@&*.                    (@#/*##(/* #*/(&#@@@%*.                                                                                 
//                                                    ...,,*/(%@@@@@@@@@&&%%@@@@@@&.             .@&@@@&&&%@@#(#@@@&@%(*.                                                                                 
//                                                            ..,%%&#&@@@@@@@&&&%&@@@@@@&#. ,%((%&&&&(//(&@@@%%#&@@@%(*.                                                                                  
//                                                                ,**#@%#%&@@@@@@@@@&&%#%@@@@&%*/(#(%&&%&@@@&#/,./@((*.                                                                                   
//                                                                   /@&%@@@%@&&@@@@@@@@@@@&%%(%@&#&&%%###%&#/.   .,.                                                                                     
//                                                                    .,(@&@@@@@&@&%@@&%@@@@@@@&&&##(##&@@@@&#&&&/                                                                                        
//                                                                     (&#//&@@@@@#/*,,,,*#&@@@@@@@&&@@%%#((/&%###(%&&&%*                                                                                 
//                                                                   ,@%&&&%@&%#(/,.       .,*/((@&&%&%&&%##&%(/**//((###%&&&%(,                                                                          
//                                                                    (&@@%@&#/,               *@&%%%&@@@@%&@%%%#((((*//(((((%%&&&&&%/,                                                                   
//                                                                     .,*///,.              *@&&%#&@@%#(((#&#&&&%&%%###&%%#/*,((###%%#&&&(.                                                              
//                                                                                           @@&%(@@&#(*. .,**/#((,(@@&%%%%&%%%((((/ .(####&&&#                                                           
//                                                                                           .,(%&%#(*.         ,/#(((//(#%&@&%#%######(/**,/(##%%&&#/                                                    
//                                                                                              ..,,.             .*((#(,..*((*/#&@&%&&&&%%#(#(*,,/(%%%%&&&(*                                             
//                                                                                                                   ,/###/  .,*#(**(##%%&@@@@&&@&##((,..*((#%%%&#,. ##                                   
//                                                                                                                     .,/&#(*   .,*//((/*,,*//(#@@@%%%#%#(//. .*(((%@@%*                                 
//                                                                                                                        .,*/#%%/.     .,.,,,,/#%%%%%%##%%#/(%%%&&&&@&(@@@@@@&%#*                        
//                                                                                                                            ..,,,,/%@@@@&&&&&&&&&&&%&&&&&%&/&&&&&&&&@%&%,,*#%&@@@@@@@@@%(/,      .      
//                                                                                                                                  ..,*/(%(&&%%&%&&&&&&&&&&%(#%&&&&&&#,,*,@@&%%(,*/#%&@@@@@@@@@@@&@,.    
//                                                                                                                                        /(@*/*#@%(#%**/*//////**,,.,/*,#@&@@@@@@@&&#/,.,*((&&&&@@&(*.   
//                                                                                                                                         .*%@&%(**/#&&%.     (%#//%%###(*,,*/#@@@@@@@@@@&#%/**//@#/,    
//                                                                                                                                           .,*(%%&&&&%%##%@*#%##%#(/*.        .*/(%@@@@@@@@@&&@@#(,     
//                                                                                                                                               .,*((###%%##/#(*.(,                .,*/(&@@@@@@@#(,      
//                                                                                                                                                      .....,**,/,,,.                   .,*/(@&#/,       
//                                                                                                                                                             .....                        ./,/*.        
//                                                                                                                                                                                           ..          
//                                                                                                                                                                                                        
                                                                                                                                                                                                        

// ***************
