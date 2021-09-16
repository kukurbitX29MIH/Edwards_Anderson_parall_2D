// Last modified 09.02.2021 PadalkoMA
// how to compile:   g++ code.cpp  -o prg.exe -lgmp -lgmpxx
// how to launch:    ./prg.exe

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
#include <string.h>
#include <gmp.h>
#include <gmpxx.h>
//#include <mpfr.h>
#include "anyfunctions.cpp"
#include "rdk_io_formatting.cpp"
#include "rdk_io_formatting_ext.cpp"






//                    ___    ___________________________
//                   \    ________\_______   _______   /
//        ________    ______________________________/ /
//       \  ____________\_________   ____________/                                              
//            _______________   _______          
//     _____________________\___  ___________________\             
//                        
//    



using namespace std;
using namespace _NMSP_RDK_;







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

//int  matr_t1_to_chstr_hist(matr_t1 * matr_hist,  char *& chstr_text,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment=true,  int Ly=0,  int Lx=0);
//int  matr_t1_to_file(matr_t1 * matr_t1_var,  char * chstr_filename,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment=true,  int Ly=0,  int Lx=0,  bool newfile=true);
//void matr_t1_to_screen(matr_t1 * matr_hist,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment=true,  int Ly=0,  int Lx=0);
//int  matr_t1_to_short_display(matr_t1 * matr_hist,  int E_down,  int M_left,  int dif_E_betw_y_cell,  int dif_M_betw_x_cell,  bool print_comment=true,  int Ly=0,  int Lx=0);
//bool extract_2D_rect_lattice_hist_from_file(char * filename, int Ly, int Lx, matr_t1 *matr, int E_min, int E_max, int M_min, int M_max, int dE_dcell, int dM_dcell);


//       print common scheme of lattice with bonds.    1) ChInterBound__Row_Col_Srf=7=0b111 --> y and x PBC,  2) ChInterBound__Row_Col_Srf=0=0b000 --> y and x FBC
int  print_spin2D_lattice_with_bonds_to_chstr(char * & chstr_text, unsigned long long int N_Row, unsigned long long int N_Col, char elem='*', char ChInterBound__Row_Col_Srf=7);
void print_spin2D_lattice_with_bonds(unsigned long long int N_Row, unsigned long long int N_Col, char elem='*', char ChInterBound__Row_Col_Srf=7);


//  write Ising lattice with vacancies struct to var hstr_text.  vacan_array_  can be empty array so vacancies there are not 
int print_to_chstr_2D_Ising_with_vac_lat_struct(char * & chstr_text, bool * vacan_array_, int Ly, int Lx, char vac_=' ', char spin_='*', char ver_frame_=' ', char hor_frame_='=');

//  write Ising lattice with vacancies struct to screen.  vacan_array_  can be empty array so vacancies there are not 
void print_to_screen_2D_Ising_with_vac_lat_struct(bool * vacan_array_, int Ly,  int Lx,  char vac_=' ', char spin_='*', char ver_frame_=' ', char hor_frame_='=');

//  write Ising lattice with vacancies struct to file chstr_filename.  vacan_array_  can be empty array so vacancies there are not 
int print_to_file_2D_Ising_with_vac_lat_struct(char *chstr_filename, bool * vacan_array_,  int Ly,  int Lx, char vac_=' ', char spin_='*', char ver_frame_=' ', char hor_frame_='=',  bool newfile=false);


//  write Edwards-Anderson(EA) lattice with vacancies struct to var hstr_text
int print_to_chstr_2D_EA_with_vac_lat_struct(char * & chstr_text, bool * vacan_array_, int Ly, int Lx, bool * J_ver,  bool * J_hor,  char vac_=' ', char spin_='*', char ver_frame_=' ', char hor_frame_='=');

//  write Edwards-Anderson(EA) lattice with vacancies struct to screen
void print_to_screen_2D_EA_with_vac_lat_struct(bool * vacan_array_, int Ly,  int Lx,  bool * J_ver,  bool * J_hor,  char vac_=' ', char spin_='*', char ver_frame_=' ', char hor_frame_='=');

//  write Edwards-Anderson(EA) lattice with vacancies struct to file chstr_filename
int print_to_file_2D_EA_with_vac_lat_struct(char *chstr_filename, bool * vacan_array_,  int Ly,  int Lx, bool * J_ver,  bool * J_hor,  char vac_=' ', char spin_='*', char ver_frame_=' ', char hor_frame_='=',  bool newfile=false);


//  write Edwards-Anderson(EA) lattice with spin states and struct to var chstr_text
int print_to_chstr_2D_EA_with_spin_state_lat_struct(char * & chstr_text, bool * spin_array_, int Ly, int Lx, bool * J_ver, bool * J_hor, char spin_up_='u', char spin_dn_='d', char ver_frame_=' ', char hor_frame_='=');

//  write Edwards-Anderson(EA) lattice with spin states and struct to screen
void print_to_screen_2D_EA_with_spin_state_lat_struct(bool * spin_array_, int Ly, int Lx, bool * J_ver, bool * J_hor, char spin_up_='u', char spin_dn_='d', char ver_frame_=' ', char hor_frame_='=');

//  write Edwards-Anderson(EA) lattice with spin states and struct to file chstr_filename
int print_to_file_2D_EA_with_spin_state_lat_struct(char *chstr_filename, bool * spin_array_, int Ly, int Lx, bool * J_ver, bool * J_hor, char spin_up_='u', char spin_dn_='d', char ver_frame_=' ', char hor_frame_='=', bool newfile=false);


//  write spin configuration to var chstr_text:   j me
int print_to_chstr_2D_spin_lattice(char * & chstr_text, const bool * const spin_array,  int Ly,  int Lx, char up_spin_arg='1',  char down_spin_arg='0');

//  write spin configuration to screen
void print_to_screen_2D_spin_lattice(const bool * const spin_array,  int Ly,  int Lx, char up_spin_arg='1',  char down_spin_arg='0');

//  write spin configuration to file chstr_filename
int print_to_file_2D_spin_lattice(char *chstr_filename, const bool * const spin_array,  int Ly,  int Lx, char up_spin_arg='1',  char down_spin_arg='0', bool is_new=true);


int print_to_chstr_J_2D_lattice_PBC(char * & chstr_text, bool * J_ver,  bool * J_hor,  int Ly,  int Lx);

int print_to_chstr_J_2D_lattice_FBC(char * & chstr_text, bool * J_ver,  bool * J_hor,  int Ly,  int Lx);

void print_to_screen_J_2D_lattice_PBC(bool * J_ver,  bool * J_hor,  int Ly,  int Lx);

void print_to_screen_J_2D_lattice_FBC(bool * J_ver,  bool * J_hor,  int Ly,  int Lx);

int print_to_file_J_2D_lattice_PBC(char *chstr_filename, bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  bool newfile);

int print_to_file_J_2D_lattice_FBC(char *chstr_filename, bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  bool newfile);


void print_to_screen_2D_spin_state_comparing_lattice(const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx);

int print_to_file_2D_spin_state_comparing_lattice(char *chstr_filename, const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx,  bool newfile);

int print_to_chstr_2D_spin_state_comparing_lattice(char * & chstr_text, const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx);


int  write_2D_spin_config_array_to_chstr(char * & chstr_text, bool *spin_array, int Ly, int Lx, int subLy=-1, int subLx=-1, int start_y=0, int start_x=0,  int dy=0, int dx=0, bool print_comment=true);

int  write_2D_spin_config_array_to_screen(bool *spin_array, int Ly, int Lx, int subLy=-1, int subLx=-1, int start_y=0, int start_x=0,  int dy=0, int dx=0, bool print_comment=true);

int  write_spin_config_array_to_screen(char *chstr_filename, bool *spin_array, int Ly, int Lx, int subLy=-1, int subLx=-1, int start_y=0, int start_x=0,  int dy=0, int dx=0, bool print_comment=true);

int extract_2D_spin_config_array_from_file(char *chstr_filename, bool * const spin_array, int Ly, int Lx, int subLy=-1, int subLx=-1, int start_y=0, int start_x=0,  int dy=0, int dx=0);


// already implemented in rdk_io_formatting_ext
//bool write_mpf_array_to_chstr(mpf_t * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len,  int wished_len=-1);

// already implemented in rdk_io_formatting_ext
// bool extract_mpz_column_from_file(char * filename, mpz_t *& col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);
    
int write_2D_rect_lattice_info_to_chstr_v15(char * & text, int Ly, int Lx, unsigned long long int state=~0, int cur_thrd_num=-1, int thrd_amount=-1, int J_plus_perc=-1, int J_min_perc=-1, int min_energ=-999999999, bool pbc1fbc0=true);

int write_2D_rect_J_bond_table_to_chstr(char * &text, int array_size, int *spin1_index_y, int *spin1_index_x, int *spin2_index_y, int *spin2_index_x, bool *J_val, bool print_comment=true);

int write_2D_rect_J_bond_table_to_chstr(char * & text, int Ly, int Lx, bool *J_ver, bool *J_hor, char d_bonds_dcn='d', char u_bonds_ucn='u', char l_bonds_lcn='l', char r_bonds_rcn='r', int y_start=0, int y_end=-1, int x_start=0, int x_end=-1, int out_dx=0, int out_dy=0, bool print_comment=true);


int write_2D_rect_lattice_info_to_screen_v15(char * text, int Ly, int Lx, unsigned long long int state=~0, int cur_thrd_num=-1, int thrd_amount=1, int J_plus_perc=0, int J_min_perc=100, int min_energ=-999999999, bool pbc1fbc0=true);

int write_2D_rect_J_bond_table_to_screen(int array_size, int *spin1_index_y, int *spin1_index_x, int *spin2_index_y, int *spin2_index_x, bool *J_val, bool print_comment=true);

//int write_2D_rect_J_bond_table_to_screen(int Ly, int Lx, bool *J_ver, bool *J_hor, bool do_bottom_bond=true, bool do_top_bond=true, bool do_left_bond=true, bool do_right_bond=true, bool print_comment=true);
int write_2D_rect_J_bond_table_to_screen( int Ly, int Lx, bool *J_ver, bool *J_hor, char d_bonds_dcn='d', char u_bonds_ucn='u', char l_bonds_lcn='l', char r_bonds_rcn='r', int y_start=0, int y_end=-1, int x_start=0, int x_end=-1, int out_dx=0, int out_dy=0, bool print_comment=true);


int write_2D_rect_lattice_info_to_file_v15(char * filename, int Ly, int Lx, unsigned long long int state=~0, int cur_thrd_num=-1, int thrd_amount=1, int J_plus_perc=0, int J_min_perc=100, int min_energ=-999999999, bool pbc1fbc0=true, bool newfile=true);

int write_2D_rect_J_bond_table_to_file(char * filename, int array_size, int *spin1_index_y, int *spin1_index_x, int *spin2_index_y, int *spin2_index_x, bool *J_val, bool print_comment=true, bool newfile=true);

int write_2D_rect_J_bond_table_to_file(char * filename, int Ly, int Lx, bool *J_ver, bool *J_hor, char d_bonds_dcn='d', char u_bonds_ucn='u', char l_bonds_lcn='l', char r_bonds_rcn='r', int y_start=0, int y_end=-1, int x_start=0, int x_end=-1, int out_dx=0, int out_dy=0, bool print_comment=true, bool newfile=true);




bool extract_2D_rect_fbc_lattice_serv_info_from_file_v15(char * filename, int *Ly, int *Lx, unsigned long long int *state, int *cur_thrd_num, int *thrd_amount, int *J_plus_perc, int *J_min_perc, int *min_energ);

bool extract_2D_rect_lattice_J_table_from_file_dyn(char * filename, int Ly, int Lx, bool *& J_ver, bool *& J_hor, char d_bonds_dcn='d', char u_bonds_ucn='u', char l_bonds_lcn='l', char r_bonds_rcn='r', int y_start=0, int y_end=-1, int x_start=0, int x_end=-1, int index_dy=0, int index_dx=0);

bool extract_2D_rect_lattice_J_table_from_file(char *filename, int & array_size, int *& spin1_index_y, int *& spin1_index_x, int *& spin2_index_y, int *& spin2_index_x, bool *& J_val, bool *part_was_missed=0);  


//  clusters
//
int print_to_chstr_clusters_id_for_2d_lattice(char *& chstr_text, int Ly,  int Lx,  int *cluster_id);
void print_to_screen_clusters_id_for_2d_lattice(int Ly,  int Lx,  int *cluster_id);
int print_to_file_clusters_id_for_2d_lattice(char *chstr_filename,  int Ly,  int Lx,  int *cluster_id,  bool newfile=false);



//  ----------- v54
//
bool extract_2D_rect_gbc_lattice_ext_edge_and_up_gs_row_from_file_v54(char * filename, int Ly, int Lx, bool * J_ver,  bool * J_hor,  bool & is_up_gs_found,  bool *& bool_edge_spin_exist_d_u_l_r_ar,  bool * bool_up_row_gs_ar, int *min_energ=0);


bool print_to_file_2D_rect_gbc_lattice_info_and_ext_edge_and_up_gs_row_v54(int Ly, int Lx, bool * J_ver,  bool * J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, int *min_energ=0, bool *bool_up_row_gs_ar=0,  char * dest_filename_arg=0, int J_plus_perc=50, int J_min_perc=50,  bool is_write_lat_struct=false,  bool pbc1fbc0=false,  bool print_J_list=true, bool is_newfile=true);

bool extract_lattice_gbc_data_v54(char * filename, int &Ly, int &Lx, bool * & J_ver,  bool * &J_hor,  int  & J_plus_perc, int & J_min_perc, bool & is_edges_extracted, bool & is_up_gs_found, bool *& bool_edge_spin_exist_d_u_l_r_ar,  bool *& bool_up_row_gs_ar,  int *min_energ=0, bool pbc1fbc0=false,  bool is_print_comments=false,  bool is_to_extract_up_row_gs=false); 








//
//  lattice consists of:   cells, nodes, ver_lines, hor_lines, linefeeds 
//
 

struct st_char_str_lattice_console_draw
{

  st_char_str_lattice_console_draw() {chstr_lattice=0;    lattice_sz_y=0;  lattice_sz_x=0;  lattice_sz_incrrem_x=0;   size=0;}
  bool set(int cells_amount_y_,  int cells_amount_x_, int cell_sz_in_y_=7,  int cell_sz_in_x_=7,  int node_sz_y_=1,  int node_sz_x_=1);
  void free() {  if (chstr_lattice) {delete[] chstr_lattice;  chstr_lattice=0;}  cells_amount_y=0; cells_amount_x=0;   cell_sz_in_y=0; cell_sz_in_x=0;  node_sz_y=0; node_sz_x=0;   lattice_sz_y=0; lattice_sz_x=0;  lattice_sz_incrrem_x=0; size=0;}
  ~st_char_str_lattice_console_draw() {free();}
  
  char *chstr_lattice;
  char node_main_symb, ver_line_main_symb,  hor_line_main_symb;
  
  int cells_amount_y,  cells_amount_x;            //  !  !  !
  int cell_sz_in_y,  cell_sz_in_x;                //  !  !  !
  int node_sz_y,  node_sz_x;                      //  !  !  !

  int lattice_sz_y,  lattice_sz_x,  lattice_sz_incrrem_x;
  int size;

  //  proc
  //
  int get_cell_dl_y0(int cell_num_y) {int p=node_sz_y+ (cell_sz_in_y+node_sz_y)*cell_num_y;  if (p<lattice_sz_y) {return p;} return 0;}
  int get_cell_dl_x0(int cell_num_x) {int p=node_sz_x+ (cell_sz_in_x+node_sz_x)*cell_num_x;  if (p<lattice_sz_x) {return p;} return 0;}
  int get_cell_dl_pos(int cell_num_y, int cell_num_x    ) {int p=(node_sz_y+ (cell_sz_in_y+node_sz_y)*cell_num_y)*lattice_sz_incrrem_x  +  (node_sz_x+ (cell_sz_in_x+node_sz_x)*cell_num_x);  if (p<size) {return p;} return 0;}
  //int get_cell_middle_pos(int cell_num_y, int cell_num_x) {int p=(get_cell_dl_y0(cell_num_y)+cell_sz_in_y/2)*lattice_sz_incrrem_x  +  (get_cell_dl_x0(cell_num_x))+cell_sz_in_x/2;  if (p<size) {return p;} return 0;}  
  int get_cell_middle_pos(int cell_num_y, int cell_num_x) {int p=(node_sz_y+ (cell_sz_in_y+node_sz_y)*cell_num_y+cell_sz_in_y/2)*lattice_sz_incrrem_x  +  (node_sz_x+ (cell_sz_in_x+node_sz_x)*cell_num_x)+cell_sz_in_x/2;  if (p<size) {return p;} return 0;}
  int get_cell_dl_pos_with_shift(int cell_num_y, int cell_num_x, int shft_y, int shft_x) 
  {int p=(node_sz_y+ (cell_sz_in_y+node_sz_y)*cell_num_y/*+cell_sz_in_y/2*/)*lattice_sz_incrrem_x  +  (node_sz_x+ (cell_sz_in_x+node_sz_x)*cell_num_x)+/*cell_sz_in_x/2+ */shft_y*lattice_sz_incrrem_x+shft_x;  if (p<size) {return p;} return 0;} 
  
  int get_ver_line_l_y0(int cell_num_y) {return cell_num_y*(cell_sz_in_y+node_sz_y)+node_sz_y;}
  int get_ver_line_l_x0(int cell_num_x) {return cell_num_x*(cell_sz_in_x+node_sz_x);}
  int get_ver_line_r_y0(int cell_num_y) {return cell_num_y*(cell_sz_in_y+node_sz_y)+node_sz_y;}
  int get_ver_line_r_x0(int cell_num_x) {return cell_num_x*(cell_sz_in_x+node_sz_x)+node_sz_x+cell_sz_in_x;}

  int get_hor_line_d_y0(int cell_num_y)  {return cell_num_y*(cell_sz_in_y+node_sz_y);}
  int get_hor_line_d_x0(int cell_num_x)  {return cell_num_x*(cell_sz_in_x+node_sz_x)+node_sz_x;}
  int get_hor_line_u_y0(int cell_num_y) {return cell_num_y*(cell_sz_in_y+node_sz_y)+node_sz_y+cell_sz_in_y;}
  int get_hor_line_u_x0(int cell_num_x) {return cell_num_x*(cell_sz_in_x+node_sz_x)+node_sz_y;}


  int get_node_ld_y0(int node_num_y) {return node_num_y*(cell_sz_in_y+node_sz_y);}
  int get_node_ld_x0(int node_num_x) {return node_num_x*(cell_sz_in_x+node_sz_x);}
  int get_node_lu_y0(int node_num_y) {return node_num_y*(cell_sz_in_y+node_sz_y)+cell_sz_in_y+node_sz_y;}
  int get_node_lu_x0(int node_num_x) {return node_num_x*(cell_sz_in_x+node_sz_x);}

  int get_node_rd_y0(int node_num_y) {return node_num_y*(cell_sz_in_y+node_sz_y);}
  int get_node_rd_x0(int node_num_x) {return node_num_x*(cell_sz_in_x+node_sz_x)+cell_sz_in_x+node_sz_x;}
  int get_node_ru_y0(int node_num_y) {return node_num_y*(cell_sz_in_y+node_sz_y)+cell_sz_in_y+node_sz_y;}
  int get_node_ru_x0(int node_num_x) {return node_num_x*(cell_sz_in_x+node_sz_x)+cell_sz_in_x+node_sz_x;}


  int get_pos_ver_line_l_yx_with_shft(int cell_num_y, int cell_num_x, int shft_y, int shft_x) {int p=get_ver_line_l_y0(cell_num_y)*lattice_sz_incrrem_x+get_ver_line_l_x0(cell_num_x)+shft_y*lattice_sz_incrrem_x+shft_x; if (p<size) {return p;} return 0;}
  int get_pos_ver_line_r_yx_with_shft(int cell_num_y, int cell_num_x, int shft_y, int shft_x) {int p=get_ver_line_r_y0(cell_num_y)*lattice_sz_incrrem_x+get_ver_line_r_x0(cell_num_x)+shft_y*lattice_sz_incrrem_x+shft_x; if (p<size) {return p;} return 0;}
  int get_pos_hor_line_d_yx_with_shft(int cell_num_y, int cell_num_x, int shft_y, int shft_x) {int p=get_hor_line_d_y0(cell_num_y)*lattice_sz_incrrem_x+get_hor_line_d_x0(cell_num_x)+shft_y*lattice_sz_incrrem_x+shft_x; if (p<size) {return p;} return 0;}
  int get_pos_hor_line_u_yx_with_shft(int cell_num_y, int cell_num_x, int shft_y, int shft_x) {int p=get_hor_line_u_y0(cell_num_y)*lattice_sz_incrrem_x+get_hor_line_u_x0(cell_num_x)+shft_y*lattice_sz_incrrem_x+shft_x; if (p<size) {return p;} return 0;}

  int get_pos_node_ld_yx_with_shft(int node_num_y, int node_num_x, int shft_y, int shft_x) {int p=get_node_ld_y0(node_num_y)*lattice_sz_incrrem_x+get_node_ld_x0(node_num_x)+shft_y*lattice_sz_incrrem_x+shft_x; if (p<size) {return p;} return 0;}
  int get_pos_node_lu_yx_with_shft(int node_num_y, int node_num_x, int shft_y, int shft_x) {int p=get_node_lu_y0(node_num_y)*lattice_sz_incrrem_x+get_node_lu_x0(node_num_x)+shft_y*lattice_sz_incrrem_x+shft_x; if (p<size) {return p;} return 0;}
  int get_pos_node_rd_yx_with_shft(int node_num_y, int node_num_x, int shft_y, int shft_x) {int p=get_node_rd_y0(node_num_y)*lattice_sz_incrrem_x+get_node_rd_x0(node_num_x)+shft_y*lattice_sz_incrrem_x+shft_x; if (p<size) {return p;} return 0;}
  int get_pos_node_ru_yx_with_shft(int node_num_y, int node_num_x, int shft_y, int shft_x) {int p=get_node_ru_y0(node_num_y)*lattice_sz_incrrem_x+get_node_ru_x0(node_num_x)+shft_y*lattice_sz_incrrem_x+shft_x; if (p<size) {return p;} return 0;}
  
  void draw_node(int node_num_y, int node_num_x, char symb);    
  void draw_nodes();
  void draw_pbc__delete_some_ver_lines() {  if (chstr_lattice) {  for (int i=0;  i < cells_amount_y;  i++) {draw_ver_line(i, cells_amount_x, ' ');  draw_node(i, cells_amount_x, ' ');}  draw_node(cells_amount_y, cells_amount_x, ' ');}  }
  void draw_pbc__delete_some_hor_lines() {  if (chstr_lattice) {  for (int i=0;  i < cells_amount_x;  i++) {draw_hor_line(cells_amount_y, i, ' ');  draw_node(cells_amount_y, i, ' ');}  draw_node(cells_amount_y, cells_amount_x, ' ');}  }
  void draw_ver_line(int num_y, int num_x, char symb);
  void draw_hor_line(int num_y, int num_x, char symb);
  void draw_ver_lines();
  void draw_hor_lines();
  void draw_space() {  if (chstr_lattice) {  for (int i=0;  i < size;  i++) {chstr_lattice[i]=' ';}  }  }
  void draw_line_feeds() {  if (chstr_lattice) {  for (int i=0;  i < lattice_sz_y;  i++) {chstr_lattice[i*lattice_sz_incrrem_x+lattice_sz_x]='\n';}  }  }  

  void draw_frame(char hor_frame_symb='-',  char ver_frame_symb='|');  
  void draw_ver_lines_type1();
  void draw_hor_lines_type1();


  void draw_num_in_ver_lines_type1_fbc(int *ver_num_array);
  void draw_num_in_hor_lines_type1_fbc(int *hor_num_array);
  void draw_num_in_ver_lines_type1_pbc(int *ver_num_array);
  void draw_num_in_hor_lines_type1_pbc(int *hor_num_array);

  void draw_star_in_hor_lines_type1_fbc(int *transv_h_array);
  void draw_star_in_hor_lines_type1_pbc(int *transv_h_array);

  void draw_num_in_ver_lines_type2_fbc(int *ver_num_array);
  void draw_num_in_hor_lines_type2_fbc(int *hor_num_array);
  void draw_num_in_ver_lines_type2_pbc(int *ver_num_array);
  void draw_num_in_hor_lines_type2_pbc(int *hor_num_array);

  void draw_star_in_hor_lines_type2_fbc(int *transv_h_array);  
  void draw_star_in_hor_lines_type2_pbc(int *transv_h_array);

  
  void print_to_screen()  {    if (chstr_lattice) {  for (int y=lattice_sz_y-1; y >= 0; y--) {for (int x=0; x < lattice_sz_incrrem_x; x++) {cout<<chstr_lattice[y*lattice_sz_incrrem_x+x];}}  }      }
  
  void lat_draw_type_and_print_to_screen_fbc() {    if (chstr_lattice) {draw_space(); draw_line_feeds();  draw_ver_lines();  draw_hor_lines();   draw_nodes(); print_to_screen();}    }  
  void lat_draw_type_1_and_print_to_screen_fbc() {    if (chstr_lattice) {draw_space(); draw_line_feeds();  draw_ver_lines_type1();  draw_hor_lines_type1(); draw_nodes(); print_to_screen();}    }
  void lat_draw_with_nums_type_1_and_print_to_screen_fbc(int *ver_num_array, int *hor_num_array, int *transv_h_array=0);  
  void lat_draw_with_nums_type_2_and_print_to_screen_fbc(int *ver_num_array, int *hor_num_array, int *transv_h_array=0); 
  void lat_draw_with_nums_type_1_and_print_to_screen_pbc(int *ver_num_array, int *hor_num_array, int *transv_h_array=0);  
  void lat_draw_with_nums_type_2_and_print_to_screen_pbc(int *ver_num_array, int *hor_num_array, int *transv_h_array=0); 


  void draw_fbc_J_bond_lattice(bool *J_ver, bool *J_hor, bool is_do_out_frame=true, char bond_bkg_symb=' ', char mark_spin_by_symb='S');
  void draw_fbc_spin_lattice(bool *spin_array, bool is_do_out_frame, char spin_up='u', char spin_down='d');
  void draw_fbc_J_bond_and_spin_lattice(bool *J_ver, bool *J_hor, bool *spin_array, bool is_do_out_frame=true, char bond_bkg_symb=' ', char spin_up='u', char spin_down='d');
  //  spin_array can be 0
  void draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond(bool *J_ver, bool *J_hor, bool *spin_array, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, bool is_do_out_frame=true, char bond_bkg_symb=' ', char spin_up='u', char spin_down='d');
  //  spin_array can be 0        J_ver and J_hor can be 0        J_bond_y_index and  int J_bond_x_index  can be negative, it means there is no marking of flip bond   
  void draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame(bool *J_ver, bool *J_hor, bool *spin_array, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x, bool is_do_out_frame=true, char bond_bkg_symb=' ', char spin_up='u', char spin_down='d', char inner_frame_symb=':');
  void draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame_cmp(bool *J_ver, bool *J_hor, bool *spin_array_to_cmp, bool *spin_array_cur, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x, bool is_do_out_frame, char bond_bkg_symb, char spin_up_same='u', char spin_down_same='d', char spin_up_diff='U', char spin_down_diff='D', char inner_frame_symb=':');
  
  //
  void draw_fbc_J_bond_and_spin_lattice_frame_proj_v3(bool *J_ver, bool *J_hor, bool *spin_array_, bool * array_spin_exist_d_u_l_r, bool is_do_out_frame=true, char bond_bkg_symb=' ', char spin_up='u', char spin_down='d', bool draw_on_screen=true, bool draw_to_file=false, char *filename_dest=0, bool newfile=false);
  //
  void draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame_cmp_frame_proj_v3(bool *J_ver, bool *J_hor, bool *spin_array_to_cmp_, bool *spin_array_cur_, bool * array_spin_exist_d_u_l_r, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int start_y,  int start_x, bool is_do_out_frame, char bond_bkg_symb, char spin_up_same='u', char spin_down_same='d', char spin_up_diff='U', char spin_down_diff='D', char inner_frame_symb=':',  char domain_wall_symb='.', bool is_draw_domain_wall=true);
  
  
  
};




//
//
//
//
//
//
//
//















int print_spin2D_lattice_with_bonds_to_chstr(char * & chstr_text, unsigned long long int N_Row, unsigned long long int N_Col, char elem, char ChInterBound__Row_Col_Srf)
{

  if (chstr_text != 0) {return 0;}   


string string_var; 
char chstr_elem_loc[8];  
chstr_elem_loc[0]='*';  chstr_elem_loc[1]='\0';
chstr_elem_loc[0]=elem; chstr_elem_loc[1]='\0';


  for (unsigned long long int i=0;  i < N_Row-1; i++)
  {
    for (unsigned long long int j=0;  j < N_Col-1; j++)
    {string_var.append(chstr_elem_loc);  string_var.append((char *)  "--");}
  
  string_var.append(chstr_elem_loc);

    if ((ChInterBound__Row_Col_Srf == 2) || (ChInterBound__Row_Col_Srf == 3) || (ChInterBound__Row_Col_Srf == 6) || (ChInterBound__Row_Col_Srf == 7))
    {string_var.append((char *)  "--");}

  string_var.append((char *)  "\n");

    for (unsigned long long int j=0;  j < N_Col-1; j++)
    {string_var.append((char *)  "|  ");}
  string_var.append((char *)  "|\n");
  }


  for (unsigned long long int j=0;  j < N_Col-1; j++)
  {string_var.append(chstr_elem_loc);  string_var.append((char *)  "--");}
  
string_var.append(chstr_elem_loc);
  if ((ChInterBound__Row_Col_Srf == 2) || (ChInterBound__Row_Col_Srf == 3) || (ChInterBound__Row_Col_Srf == 6) || (ChInterBound__Row_Col_Srf == 7))
  {string_var.append((char *)  "--");}

string_var.append((char *)  "\n");

  if ((ChInterBound__Row_Col_Srf == 4) || (ChInterBound__Row_Col_Srf == 5) || (ChInterBound__Row_Col_Srf == 6) || (ChInterBound__Row_Col_Srf == 7))
  {
    for (unsigned long long int j=0;  j < N_Col-1; j++)
    {string_var.append((char *)  "|  ");}
  string_var.append((char *)  "|\n");
  }


//string_var.append((char *)  "\n");


chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';
int len=strlen(chstr_text);
string_var.clear();


return len;

}



void print_spin2D_lattice_with_bonds(unsigned long long int N_Row, unsigned long long int N_Col, char elem, char ChInterBound__Row_Col_Srf)
{

  for (unsigned long long int i=0;  i < N_Row-1; i++)
  {
    for (unsigned long long int j=0;  j < N_Col-1; j++)
    {cout<<elem<<"--";}
  
  cout<<elem;

    if ((ChInterBound__Row_Col_Srf == 2) || (ChInterBound__Row_Col_Srf == 3) || (ChInterBound__Row_Col_Srf == 6) || (ChInterBound__Row_Col_Srf == 7))
    {cout<<"--";}

  cout<<endl;

    for (unsigned long long int j=0;  j < N_Col-1; j++)
    {cout<<"|  ";}
  cout<<"|"<<endl;
  }


  for (unsigned long long int j=0;  j < N_Col-1; j++)
  {cout<<elem<<"--";}
  
cout<<elem;
  if ((ChInterBound__Row_Col_Srf == 2) || (ChInterBound__Row_Col_Srf == 3) || (ChInterBound__Row_Col_Srf == 6) || (ChInterBound__Row_Col_Srf == 7))
  {cout<<"--";}

cout<<endl;

  if ((ChInterBound__Row_Col_Srf == 4) || (ChInterBound__Row_Col_Srf == 5) || (ChInterBound__Row_Col_Srf == 6) || (ChInterBound__Row_Col_Srf == 7))
  {
    for (unsigned long long int j=0;  j < N_Col-1; j++)
    {cout<<"|  ";}
  cout<<"|"<<endl;
  }

cout<<endl;

}







//  write spin configuration to var chstr_text:   j me
int print_to_chstr_2D_Ising_with_vac_lat_struct(char * & chstr_text,bool * vacan_array_, int Ly,  int Lx,  char vac_, char spin_, char ver_frame_, char hor_frame_)
{

  if (chstr_text != 0) {return 0;}   
  if ((Ly < 1) || (Lx < 1)) {return 0;}

string string_var; 

bool * vacan_array=0;

  if (vacan_array_) {vacan_array=vacan_array_;}
  else
  {
  vacan_array=new bool[Ly*Lx];

    for (int i=0;  i < Ly*Lx;  i++)
    {vacan_array[i]=false;}

  }


char elem[]="*";                     elem[0]=spin_;         elem[1]='\0'; 
char elem_vac[]=" ";                 elem_vac[0]=vac_;      elem_vac[1]='\0'; 
char hor_bond[]="---";                
char not_hor_bond[]="   ";
char ver_bond[]="|";
char not_ver_bond[]=" ";
char interval_space[]="   ";
char frame_hor[]="e";               frame_hor[0]=hor_frame_;  frame_hor[1]='\0'; 
char frame_ver[]=" ";               frame_ver[0]=ver_frame_;  frame_ver[1]='\0';


int  main_spin_index=0,  neigh_right_spin_index=0,  neigh_up_spin_index=0;
bool main_spin=false, neigh_right_spin=false, neigh_up_spin=false;

//  main_spin_index       =index_y*Lx+index_x;     main_spin_index       =(Ly-1)*Lx+index_x;      main_spin_index       =index_y*Lx+Lx-1;          main_spin_index       =Ly*Lx-1;
//  neigh_up_spin_index   =index_y*Lx+index_x+Lx;  neigh_up_spin_index   =index_x;                neigh_up_spin_index   =index_y*Lx+Lx-1+Lx;       neigh_up_spin_index   =Lx-1;
//  neigh_right_spin_index=index_y*Lx+index_x+1;   neigh_right_spin_index=(Ly-1)*Lx+index_x+1;    neigh_right_spin_index=index_y*Lx;               neigh_right_spin_index=(Ly-1)*Lx; 
   

  for (int u=0;  u < Lx*4+2;  u++) {string_var.append(frame_hor);}  string_var.append((char *)  "\n");                       

string_var.append(frame_ver);
  for (int index_x=0;  index_x < Lx-1;  index_x++)
  {  if ((vacan_array[(Ly-1)*Lx+index_x] == false) && (vacan_array[index_x] == false)) {string_var.append(ver_bond);} else {string_var.append(not_ver_bond);}    string_var.append(interval_space);}
  if ((vacan_array[main_spin_index] == false) && (vacan_array[Lx-1] == false)) {string_var.append(ver_bond);} else {string_var.append(not_ver_bond);}    
string_var.append(interval_space);  string_var.append(frame_ver);  string_var.append((char *)  "\n");

string_var.append(frame_ver);
  if (vacan_array[Lx*(Ly-1)] == false) {string_var.append(elem);} else {string_var.append(elem_vac);} 
  for (int index_x=0;  index_x < Lx-1;  index_x++)
  {
    if ((vacan_array[(Ly-1)*Lx+index_x] == false) && (vacan_array[(Ly-1)*Lx+index_x+1] == false)) {string_var.append(hor_bond);} else {string_var.append(not_hor_bond);}
    if (vacan_array[Lx*(Ly-1)+index_x+1] == false) {string_var.append(elem);} else {string_var.append(elem_vac);}  
  }  
  if ((vacan_array[Ly*Lx-1] == false) && (vacan_array[(Ly-1)*Lx] == false)) {string_var.append(hor_bond);} else {string_var.append(not_hor_bond);} 
string_var.append(frame_ver);  string_var.append((char *)  "\n");

  for (int index_y=Ly-2;  index_y >= 0 ;  index_y--)
  {
  string_var.append(frame_ver);
    for (int index_x=0;  index_x < Lx-1;  index_x++)
    {  if ((vacan_array[index_y*Lx+index_x] == false) && (vacan_array[index_y*Lx+index_x+Lx] == false)) {string_var.append(ver_bond);} else {string_var.append(not_ver_bond);}    string_var.append(interval_space);}
    if ((vacan_array[index_y*Lx+Lx-1] == false) && (vacan_array[index_y*Lx+Lx-1+Lx] == false)) {string_var.append(ver_bond);} else {string_var.append(not_ver_bond);}    
  string_var.append(interval_space);  string_var.append(frame_ver);  string_var.append((char *)  "\n");

  string_var.append(frame_ver);
    if (vacan_array[Lx*index_y] == false) {string_var.append(elem);} else {string_var.append(elem_vac);} 
    for (int index_x=0;  index_x < Lx-1;  index_x++)
    {
      if ((vacan_array[index_y*Lx+index_x] == false) && (vacan_array[index_y*Lx+index_x+1] == false)) {string_var.append(hor_bond);} else {string_var.append(not_hor_bond);}
      if (vacan_array[Lx*index_y+index_x+1] == false) {string_var.append(elem);} else {string_var.append(elem_vac);}  
    }  
    if ((vacan_array[index_y*Lx+Lx-1] == false) && (vacan_array[index_y*Lx] == false)) {string_var.append(hor_bond);} else {string_var.append(not_hor_bond);} 
  string_var.append(frame_ver);  string_var.append((char *)  "\n");
  }

  for (int u=0;  u < Lx*4+2;  u++) {string_var.append(frame_hor);}  string_var.append((char *)  "\n");




chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';
int len=strlen(chstr_text);
string_var.clear();


  if (vacan_array_) {vacan_array=0;} else {  if (vacan_array) {delete[] vacan_array;  vacan_array=0;}  }

return len;

}



//  write spin configuration to screen
void print_to_screen_2D_Ising_with_vac_lat_struct(bool * vacan_array_, int Ly,  int Lx,  char vac_, char spin_, char ver_frame_, char hor_frame_)
{

char * chstr_text=0;
int text_size=print_to_chstr_2D_Ising_with_vac_lat_struct(chstr_text, vacan_array_, Ly,  Lx,  vac_, spin_, ver_frame_, hor_frame_);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}



//  write spin configuration to file chstr_filename
int print_to_file_2D_Ising_with_vac_lat_struct(char *chstr_filename, bool * vacan_array_, int Ly,  int Lx,  char vac_, char spin_, char ver_frame_, char hor_frame_,  bool newfile)
{

  if ((Ly < 1) || (Lx < 1)) {return 0;}

char * chstr_text=0;
int text_size=print_to_chstr_2D_Ising_with_vac_lat_struct(chstr_text, vacan_array_, Ly,  Lx,  vac_, spin_, ver_frame_, hor_frame_);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}



//  write Edwards-Anderson(EA) lattice with vacancies struct to var hstr_text
int print_to_chstr_2D_EA_with_vac_lat_struct(char * & chstr_text, bool * vacan_array_, int Ly, int Lx, bool * J_ver,  bool * J_hor,  char vac_, char spin_, char ver_frame_, char hor_frame_)
{


  if (chstr_text != 0) {return 0;}   
  if ((Ly < 1) || (Lx < 1)) {return 0;}
  if ((J_ver == 0) || (J_hor == 0)) {return 0;}

string string_var; 

bool * vacan_array=0;

  if (vacan_array_) {vacan_array=vacan_array_;}
  else
  {
  vacan_array=new bool[Ly*Lx];

    for (int i=0;  i < Ly*Lx;  i++)
    {vacan_array[i]=false;}

  }


char elem[]="*";                     elem[0]=spin_;         elem[1]='\0'; 
char elem_vac[]=" ";                 elem_vac[0]=vac_;      elem_vac[1]='\0'; 
char hor_bond[]="---";                
char not_hor_bond[]="   ";
char ver_bond[]="|";
char not_ver_bond[]=" ";
char interval_space[]="   ";
char frame_hor[]="e";               frame_hor[0]=hor_frame_;  frame_hor[1]='\0'; 
char frame_ver[]=" ";               frame_ver[0]=ver_frame_;  frame_ver[1]='\0';
char J_plus_ver[]="+";
char J_min_ver[]="-";  
char J_plus_hor[]=" + ";
char J_min_hor[]=" - ";  


int  main_spin_index=0,  neigh_right_spin_index=0,  neigh_up_spin_index=0;
bool main_spin=false, neigh_right_spin=false, neigh_up_spin=false;

//  main_spin_index       =index_y*Lx+index_x;     main_spin_index       =(Ly-1)*Lx+index_x;      main_spin_index       =index_y*Lx+Lx-1;          main_spin_index       =Ly*Lx-1;
//  neigh_up_spin_index   =index_y*Lx+index_x+Lx;  neigh_up_spin_index   =index_x;                neigh_up_spin_index   =index_y*Lx+Lx-1+Lx;       neigh_up_spin_index   =Lx-1;
//  neigh_right_spin_index=index_y*Lx+index_x+1;   neigh_right_spin_index=(Ly-1)*Lx+index_x+1;    neigh_right_spin_index=index_y*Lx;               neigh_right_spin_index=(Ly-1)*Lx; 
   

  for (int u=0;  u < Lx*4+2;  u++) {string_var.append(frame_hor);}  string_var.append((char *)  "\n");                       

string_var.append(frame_ver);
  for (int index_x=0;  index_x < Lx-1;  index_x++)
  {  
    if ((vacan_array[(Ly-1)*Lx+index_x] == false) && (vacan_array[index_x] == false)) 
    {if (J_ver[index_x] == true) {string_var.append(J_plus_ver);} else {string_var.append(J_min_ver);}} else {string_var.append(not_ver_bond);}    
  string_var.append(interval_space);
  }
  if ((vacan_array[main_spin_index] == false) && (vacan_array[Lx-1] == false)) 
  {if (J_ver[Lx-1] == true) {string_var.append(J_plus_ver);} else {string_var.append(J_min_ver);}} 
  else {string_var.append(not_ver_bond);}    
string_var.append(interval_space);  string_var.append(frame_ver);  string_var.append((char *)  "\n");

string_var.append(frame_ver);
  if (vacan_array[Lx*(Ly-1)] == false) {string_var.append(elem);} else {string_var.append(elem_vac);} 
  for (int index_x=0;  index_x < Lx-1;  index_x++)
  {
    if ((vacan_array[(Ly-1)*Lx+index_x] == false) && (vacan_array[(Ly-1)*Lx+index_x+1] == false)) 
    {if (J_hor[(Ly-1)*(Lx+1)+index_x+1] == true) {string_var.append(J_plus_hor);} else {string_var.append(J_min_hor);}} 
    else {string_var.append(not_hor_bond);}
    if (vacan_array[Lx*(Ly-1)+index_x+1] == false) {string_var.append(elem);} else {string_var.append(elem_vac);}  
  }  
  if ((vacan_array[Ly*Lx-1] == false) && (vacan_array[(Ly-1)*Lx] == false)) 
  {if (J_hor[(Ly-1)*(Lx+1)] == true) {string_var.append(J_plus_hor);} else {string_var.append(J_min_hor);}} 
  else {string_var.append(not_hor_bond);} 
string_var.append(frame_ver);  string_var.append((char *)  "\n");

  for (int index_y=Ly-2;  index_y >= 0 ;  index_y--)
  {
  string_var.append(frame_ver);
    for (int index_x=0;  index_x < Lx-1;  index_x++)
    {  
      if ((vacan_array[index_y*Lx+index_x] == false) && (vacan_array[index_y*Lx+index_x+Lx] == false))
      {if (J_ver[index_y*Lx+index_x+Lx] == true) {string_var.append(J_plus_ver);} else {string_var.append(J_min_ver);}} 
      else {string_var.append(not_ver_bond);}    
    string_var.append(interval_space);
    }
    if ((vacan_array[index_y*Lx+Lx-1] == false) && (vacan_array[index_y*Lx+Lx-1+Lx] == false)) 
    {if (J_ver[index_y*Lx+Lx-1+Lx] == true) {string_var.append(J_plus_ver);} else {string_var.append(J_min_ver);}} else {string_var.append(not_ver_bond);}    
  string_var.append(interval_space);  string_var.append(frame_ver);  string_var.append((char *)  "\n");

  string_var.append(frame_ver);
    if (vacan_array[Lx*index_y] == false) {string_var.append(elem);} else {string_var.append(elem_vac);} 
    for (int index_x=0;  index_x < Lx-1;  index_x++)
    {
      if ((vacan_array[index_y*Lx+index_x] == false) && (vacan_array[index_y*Lx+index_x+1] == false)) 
      {if (J_hor[index_y*(Lx+1)+index_x+1] == true) {string_var.append(J_plus_hor);} else {string_var.append(J_min_hor);}} 
      else {string_var.append(not_hor_bond);}
      if (vacan_array[Lx*index_y+index_x+1] == false) {string_var.append(elem);} else {string_var.append(elem_vac);}  
    }  
    if ((vacan_array[index_y*Lx+Lx-1] == false) && (vacan_array[index_y*Lx] == false)) {if (J_hor[index_y*(Lx+1)] == true) 
    {string_var.append(J_plus_hor);} else {string_var.append(J_min_hor);}} else {string_var.append(not_hor_bond);} 
  string_var.append(frame_ver);  string_var.append((char *)  "\n");
  }

  for (int u=0;  u < Lx*4+2;  u++) {string_var.append(frame_hor);}  string_var.append((char *)  "\n");




chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';
int len=strlen(chstr_text);
string_var.clear();


  if (vacan_array_) {vacan_array=0;} else {  if (vacan_array) {delete[] vacan_array;  vacan_array=0;}  }

return len;


}




//  write Edwards-Anderson(EA) lattice with spin states and struct to var chstr_text
int print_to_chstr_2D_EA_with_spin_state_lat_struct(char * & chstr_text, bool * spin_array_, int Ly, int Lx, bool * J_ver,  bool * J_hor,  char spin_up_, char spin_dn_, char ver_frame_, char hor_frame_)
{


  if (chstr_text != 0) {return 0;}   
  if ((Ly < 1) || (Lx < 1)) {return 0;}
  if ((J_ver == 0) || (J_hor == 0)) {return 0;}

string string_var; 

bool * spin_array=0;

  if (spin_array_) {spin_array=spin_array_;}
  else
  {
  spin_array=new bool[Ly*Lx];

    for (int i=0;  i < Ly*Lx;  i++)
    {spin_array[i]=false;}

  }


char spin_up_elem[]="u";           spin_up_elem[0]=spin_up_;  spin_up_elem[1]='\0'; 
char spin_dn_elem[]="d";           spin_dn_elem[0]=spin_dn_;  spin_dn_elem[1]='\0';
char hor_bond[]="---";                
char not_hor_bond[]="   ";
char ver_bond[]="|";
char not_ver_bond[]=" ";
char interval_space[]="   ";
char frame_hor[]="e";               frame_hor[0]=hor_frame_;  frame_hor[1]='\0'; 
char frame_ver[]=" ";               frame_ver[0]=ver_frame_;  frame_ver[1]='\0';
char J_plus_ver[]="+";
char J_min_ver[]="-";  
char J_plus_hor[]=" + ";
char J_min_hor[]=" - ";  


int  main_spin_index=0,  neigh_right_spin_index=0,  neigh_up_spin_index=0;
bool main_spin=false, neigh_right_spin=false, neigh_up_spin=false;

//  main_spin_index       =index_y*Lx+index_x;     main_spin_index       =(Ly-1)*Lx+index_x;      main_spin_index       =index_y*Lx+Lx-1;          main_spin_index       =Ly*Lx-1;
//  neigh_up_spin_index   =index_y*Lx+index_x+Lx;  neigh_up_spin_index   =index_x;                neigh_up_spin_index   =index_y*Lx+Lx-1+Lx;       neigh_up_spin_index   =Lx-1;
//  neigh_right_spin_index=index_y*Lx+index_x+1;   neigh_right_spin_index=(Ly-1)*Lx+index_x+1;    neigh_right_spin_index=index_y*Lx;               neigh_right_spin_index=(Ly-1)*Lx; 
   

  for (int u=0;  u < Lx*4+2;  u++) {string_var.append(frame_hor);}  string_var.append((char *)  "\n");                       

string_var.append(frame_ver);
  for (int index_x=0;  index_x < Lx-1;  index_x++)
  {  
    //--//if ((spin_array[(Ly-1)*Lx+index_x] == false) && (spin_array[index_x] == false)) 
    {if (J_ver[index_x] == true) {string_var.append(J_plus_ver);} else {string_var.append(J_min_ver);}} //--//else {string_var.append(not_ver_bond);}    
  string_var.append(interval_space);
  }
  //--//if ((spin_array[main_spin_index] == false) && (spin_array[Lx-1] == false)) 
  {if (J_ver[Lx-1] == true) {string_var.append(J_plus_ver);} else {string_var.append(J_min_ver);}} 
  //--//else {string_var.append(not_ver_bond);}    
string_var.append(interval_space);  string_var.append(frame_ver);  string_var.append((char *)  "\n");

string_var.append(frame_ver);
  if (spin_array[Lx*(Ly-1)] == false) {string_var.append(spin_dn_elem);} else {string_var.append(spin_up_elem);}  //
  for (int index_x=0;  index_x < Lx-1;  index_x++)
  {
    //--//if ((spin_array[(Ly-1)*Lx+index_x] == false) && (spin_array[(Ly-1)*Lx+index_x+1] == false))  
    {if (J_hor[(Ly-1)*(Lx+1)+index_x+1] == true) {string_var.append(J_plus_hor);} else {string_var.append(J_min_hor);}} 
    //--//else {string_var.append(not_hor_bond);}
    if (spin_array[Lx*(Ly-1)+index_x+1] == false) {string_var.append(spin_dn_elem);} else {string_var.append(spin_up_elem);}   //
  }  
  //--//if ((spin_array[Ly*Lx-1] == false) && (spin_array[(Ly-1)*Lx] == false)) 
  {if (J_hor[(Ly-1)*(Lx+1)] == true) {string_var.append(J_plus_hor);} else {string_var.append(J_min_hor);}} 
  //--//else {string_var.append(not_hor_bond);} 
string_var.append(frame_ver);  string_var.append((char *)  "\n");

  for (int index_y=Ly-2;  index_y >= 0 ;  index_y--)
  {
  string_var.append(frame_ver);
    for (int index_x=0;  index_x < Lx-1;  index_x++) 
    {  
      //--//if ((spin_array[index_y*Lx+index_x] == false) && (spin_array[index_y*Lx+index_x+Lx] == false))
      {if (J_ver[index_y*Lx+index_x+Lx] == true) {string_var.append(J_plus_ver);} else {string_var.append(J_min_ver);}} 
      //--//else {string_var.append(not_ver_bond);}    
    string_var.append(interval_space);
    }
    //if ((spin_array[index_y*Lx+Lx-1] == false) && (spin_array[index_y*Lx+Lx-1+Lx] == false)) 
    {if (J_ver[index_y*Lx+Lx-1+Lx] == true) {string_var.append(J_plus_ver);} else {string_var.append(J_min_ver);}} //--//else {string_var.append(not_ver_bond);}    
  string_var.append(interval_space);  string_var.append(frame_ver);  string_var.append((char *)  "\n");

  string_var.append(frame_ver);
    if (spin_array[Lx*index_y] == false) {string_var.append(spin_dn_elem);} else {string_var.append(spin_up_elem);}  //
    for (int index_x=0;  index_x < Lx-1;  index_x++)
    {
      //--//if ((spin_array[index_y*Lx+index_x] == false) && (spin_array[index_y*Lx+index_x+1] == false)) 
      {if (J_hor[index_y*(Lx+1)+index_x+1] == true) {string_var.append(J_plus_hor);} else {string_var.append(J_min_hor);}} 
      //--//else {string_var.append(not_hor_bond);}
      if (spin_array[Lx*index_y+index_x+1] == false) {string_var.append(spin_dn_elem);} else {string_var.append(spin_up_elem);}    //
    }  
    //--//if ((spin_array[index_y*Lx+Lx-1] == false) && (spin_array[index_y*Lx] == false)) 
    {if (J_hor[index_y*(Lx+1)] == true) {string_var.append(J_plus_hor);} else {string_var.append(J_min_hor);}} 
    //--//else {string_var.append(not_hor_bond);} 
  string_var.append(frame_ver);  string_var.append((char *)  "\n");
  }

  for (int u=0;  u < Lx*4+2;  u++) {string_var.append(frame_hor);}  string_var.append((char *)  "\n");




chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';
int len=strlen(chstr_text);
string_var.clear();


  if (spin_array_) {spin_array=0;} else {  if (spin_array) {delete[] spin_array;  spin_array=0;}  }

return len;


}




//  write Edwards-Anderson(EA) lattice with spin states and struct to screen
void print_to_screen_2D_EA_with_spin_state_lat_struct(bool * spin_array_, int Ly, int Lx, bool * J_ver, bool * J_hor, char spin_up_, char spin_dn_, char ver_frame_, char hor_frame_)
{

char * chstr_text=0;
int text_size=print_to_chstr_2D_EA_with_spin_state_lat_struct(chstr_text, spin_array_, Ly, Lx, J_ver, J_hor, spin_up_, spin_dn_, ver_frame_, hor_frame_);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}


//  write Edwards-Anderson(EA) lattice with spin states and struct to file chstr_filename
int print_to_file_2D_EA_with_spin_state_lat_struct(char *chstr_filename, bool * spin_array_, int Ly, int Lx, bool * J_ver, bool * J_hor, char spin_up_, char spin_dn_, char ver_frame_, char hor_frame_, bool newfile)
{

  if ((Ly < 1) || (Lx < 1)) {return 0;}

char * chstr_text=0;
int text_size=print_to_chstr_2D_EA_with_vac_lat_struct(chstr_text, spin_array_, Ly, Lx, J_ver, J_hor, spin_up_, spin_dn_, ver_frame_, hor_frame_);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}





//  write Edwards-Anderson(EA) lattice with vacancies struct to screen
void print_to_screen_2D_EA_with_vac_lat_struct(bool * vacan_array_, int Ly,  int Lx,  bool * J_ver,  bool * J_hor,  char vac_, char spin_, char ver_frame_, char hor_frame_)
{

char * chstr_text=0;
int text_size=print_to_chstr_2D_EA_with_vac_lat_struct(chstr_text, vacan_array_, Ly,  Lx, J_ver, J_hor, vac_, spin_, ver_frame_, hor_frame_);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}



//  write Edwards-Anderson(EA) lattice with vacancies struct to file chstr_filename
int print_to_file_2D_EA_with_vac_lat_struct(char *chstr_filename, bool * vacan_array_,  int Ly,  int Lx, bool * J_ver,  bool * J_hor,  char vac_, char spin_, char ver_frame_, char hor_frame_,  bool newfile)
{

  if ((Ly < 1) || (Lx < 1)) {return 0;}

char * chstr_text=0;
int text_size=print_to_chstr_2D_EA_with_vac_lat_struct(chstr_text, vacan_array_, Ly,  Lx, J_ver, J_hor, vac_, spin_, ver_frame_, hor_frame_);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}





int print_to_chstr_2D_spin_lattice(char * & chstr_text, const bool * const spin_array,  int Ly,  int Lx, char up_spin_arg,  char down_spin_arg)
{


  if (chstr_text == 0)
  {
  char up_spin[]="1";
  char down_spin[]="0";
    if (up_spin_arg   != '1') {up_spin[0]=up_spin_arg;}
    if (down_spin_arg != '0') {down_spin[0]=down_spin_arg;}
  char ver_bond[]="|";
  char hor_bond[]="--";
  char space[]="  ";


  int chstr_text_size=8*strlen(up_spin)*Lx*Ly*4+1024;
  chstr_text=new char[chstr_text_size];
  chstr_text[0]='\0';


    for (int index_y_=0;  index_y_ < Ly; index_y_++)
    {
    int index_y=Ly-index_y_-1;

      for (int index_x=0;  index_x < Lx; index_x++)
      {
        if (spin_array[index_y*Lx+index_x] == false) {strcat(chstr_text, down_spin);} else {strcat(chstr_text, up_spin);}
        if (index_x < Lx-1) {strcat(chstr_text, hor_bond);}
      }

    strcat(chstr_text, (char *) "\n\0");


      if (index_y_ < Ly-1)
      {
    
        for (int index_x=0;  index_x < Lx; index_x++)
        {
        strcat(chstr_text, ver_bond);  strcat(chstr_text, space);
        }

      strcat(chstr_text, (char *) "\n\0");
      }
  
    }
  }

return strlen(chstr_text);

}



void print_to_screen_2D_spin_lattice(const bool * const spin_array,  int Ly,  int Lx, char up_spin_arg,  char down_spin_arg)
{

char * chstr_text=0;
int text_size=print_to_chstr_2D_spin_lattice(chstr_text, spin_array,  Ly,  Lx,  up_spin_arg,  down_spin_arg);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}



int print_to_file_2D_spin_lattice(char *chstr_filename, const bool * const spin_array,  int Ly,  int Lx, char up_spin_arg,  char down_spin_arg, bool is_new)
{

char * chstr_text=0;
int text_size=print_to_chstr_2D_spin_lattice(chstr_text, spin_array,  Ly,  Lx,  up_spin_arg,  down_spin_arg);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, is_new);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}





int print_to_chstr_J_2D_lattice_PBC(char * & chstr_text, bool * J_ver,  bool * J_hor,  int Ly,  int Lx)
{

  if (chstr_text) {return 0;}
  if (J_ver == 0) {return 0;}
  if (J_hor == 0) {return 0;} 

char hor_sp[]=" ";
char spin[]="*";
char posit_bond[]="+";
char negat_bond[]="-";


string string_var;

//string_var.append((char *)  "# ground state search with using Nefedev algorithm v1 for EA 2D lattice\n");
//  string_var.append((char *)  "# total lattice dim:  Lx: \0");
//  add_int_to_string(& string_var, Lx);
//  string_var.append((char *)  "\n# total lattice dim:  Ly: \0");
//  add_int_to_string(& string_var, Ly);
//  string_var.append((char *)  "\n\0");
string_var.append((char *)  "#  J-matr\n");



// ver up bonds

string_var.append(hor_sp);
string_var.append(hor_sp);

  for (int index_x=0;  index_x < Lx; index_x++)
  { 

    if (J_ver[Ly*Lx+index_x] == false) {string_var.append(negat_bond);} else {string_var.append(posit_bond);}
    if (index_x == Lx-1) {string_var.append((char *)  "\n");} else {string_var.append(hor_sp);  string_var.append(hor_sp);  string_var.append(hor_sp);}
   
  }


//  hor bonds
//  ver bonds
//  ...
//  hor bonds
//  ver bonds

  for (int index_y_=0;  index_y_ < Ly; index_y_++)
  {
  int index_y=Ly-index_y_-1;

    // hor left bond

    if (J_hor[index_y*(Lx+1)] == false) {string_var.append(negat_bond);} else {string_var.append(posit_bond);}

  string_var.append(hor_sp);


  // hor bonds
    for (int index_x=1;  index_x <= Lx; index_x++)
    { 
    string_var.append(spin);
    string_var.append(hor_sp);

      if (J_hor[index_y*(Lx+1)+index_x] == false) {string_var.append(negat_bond);} else {string_var.append(posit_bond);}
      if (index_x == Lx) {string_var.append((char *)  "\n");} else {string_var.append(hor_sp);}

    }

  // ver bonds
  string_var.append(hor_sp);
  string_var.append(hor_sp);

    for (int index_x=0;  index_x < Lx; index_x++)
    { 

      if (J_ver[index_y*Lx+index_x] == false) {string_var.append(negat_bond);} else {string_var.append(posit_bond);}

      if (index_x == Lx-1) {string_var.append((char *)  "\n");} else {string_var.append(hor_sp);  string_var.append(hor_sp);  string_var.append(hor_sp);}
   
    }

  }



chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';

return strlen(chstr_text);

}







int print_to_chstr_J_2D_lattice_FBC(char * & chstr_text, bool * J_ver,  bool * J_hor,  int Ly,  int Lx)
{

  if (chstr_text) {return 0;}
  if (J_ver == 0) {return 0;}
  if (J_hor == 0) {return 0;} 

char hor_sp[]=" ";
char spin[]="*";
char posit_bond[]="+";
char negat_bond[]="-";


string string_var;

//string_var.append((char *)  "# ground state search with using Nefedev algorithm v1 for EA 2D lattice\n");
//  string_var.append((char *)  "# total lattice dim:  Lx: \0");
//  add_int_to_string(& string_var, Lx);
//  string_var.append((char *)  "\n# total lattice dim:  Ly: \0");
//  add_int_to_string(& string_var, Ly);
//  string_var.append((char *)  "\n\0");
string_var.append((char *)  "#  J-matr\n");



// ver up bonds

string_var.append(hor_sp);
string_var.append(hor_sp);

  for (int index_x=0;  index_x < Lx; index_x++)
  { 

    if (J_ver[Ly*Lx+index_x] == false) {string_var.append(negat_bond);} else {string_var.append(posit_bond);}

    if (index_x == Lx-1) {string_var.append((char *)  "\n");} else {string_var.append(hor_sp);  string_var.append(hor_sp);  string_var.append(hor_sp);}
   
  }


//  hor bonds
//  ver bonds
//  ...
//  hor bonds
//  ver bonds

  for (int index_y_=0;  index_y_ < Ly; index_y_++)
  {
  int index_y=Ly-index_y_-1;

    // hor left bond

    if (J_hor[index_y*(Lx+1)] == false) {string_var.append(negat_bond);} else {string_var.append(posit_bond);}

  string_var.append(hor_sp);


  // hor bonds
    for (int index_x=1;  index_x <= Lx; index_x++)
    { 
    string_var.append(spin);
    string_var.append(hor_sp);

      if (J_hor[index_y*(Lx+1)+index_x] == false) {string_var.append(negat_bond);} else {string_var.append(posit_bond);}
      if (index_x == Lx) {string_var.append((char *)  "\n");} else {string_var.append(hor_sp);}

    }

  // ver bonds
  string_var.append(hor_sp);
  string_var.append(hor_sp);

    for (int index_x=0;  index_x < Lx; index_x++)
    { 

      if (J_ver[index_y*Lx+index_x] == false) {string_var.append(negat_bond);} else {string_var.append(posit_bond);}

      if (index_x == Lx-1) {string_var.append((char *)  "\n");} else {string_var.append(hor_sp);  string_var.append(hor_sp);  string_var.append(hor_sp);}
   
    }

  }



chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';


return strlen(chstr_text);

}







void print_to_screen_J_2D_lattice_PBC(bool * J_ver,  bool * J_hor,  int Ly,  int Lx)
{

char * chstr_text=0;
int text_size=print_to_chstr_J_2D_lattice_PBC(chstr_text, J_ver,  J_hor,  Ly,  Lx);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}






void print_to_screen_J_2D_lattice_FBC(bool * J_ver,  bool * J_hor,  int Ly,  int Lx)
{

char * chstr_text=0;
int text_size=print_to_chstr_J_2D_lattice_FBC(chstr_text, J_ver,  J_hor,  Ly,  Lx);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}





int print_to_file_J_2D_lattice_PBC(char *chstr_filename, bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  bool newfile)
{

char * chstr_text=0;
int text_size=print_to_chstr_J_2D_lattice_PBC(chstr_text, J_ver,  J_hor,  Ly,  Lx);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}





int print_to_file_J_2D_lattice_FBC(char *chstr_filename, bool * J_ver,  bool * J_hor,  int Ly,  int Lx,  bool newfile)
{

char * chstr_text=0;
int text_size=print_to_chstr_J_2D_lattice_FBC(chstr_text, J_ver,  J_hor,  Ly,  Lx);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}




int print_to_chstr_2D_spin_state_comparing_lattice(char * & chstr_text, const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx)
{

char up_spin_arg='u';  
char down_spin_arg='d';


  if (chstr_text == 0)
  {
  char up_spin[]="1";
  char down_spin[]="0";
    if (up_spin_arg   != '1') {up_spin[0]=up_spin_arg;}
    if (down_spin_arg != '0') {down_spin[0]=down_spin_arg;}
  char ver_bond[]="|";
  char hor_bond[]="--";
  char space[]="  ";


  int chstr_text_size=2*(8*strlen(up_spin)*Lx*Ly*4+1024)+128;
  chstr_text=new char[chstr_text_size];
  chstr_text[0]='\0';


  // ***********************************************
  // 1
  strcat(chstr_text, (char *) "comparing spin state 1:\n\0");


    for (int index_y_=0;  index_y_ < Ly; index_y_++)
    {
    int index_y=Ly-index_y_-1;

      for (int index_x=0;  index_x < Lx; index_x++)
      {
        if (spin_array_1[index_y*Lx+index_x] == false) {strcat(chstr_text, down_spin);} else {strcat(chstr_text, up_spin);}
        if (index_x < Lx-1) {strcat(chstr_text, hor_bond);}
      }

    strcat(chstr_text, (char *) "\n\0");


      if (index_y_ < Ly-1)
      {
    
        for (int index_x=0;  index_x < Lx; index_x++)
        {
        strcat(chstr_text, ver_bond);  strcat(chstr_text, space);
        }

      strcat(chstr_text, (char *) "\n\0");
      }
  
    }


  // ***********************************************
  // 2
  char up_spin_arg_2='U';  
  char down_spin_arg_2='D';

  char up_spin_2[]="1";
  char down_spin_2[]="0";
    if (up_spin_arg_2   != '1') {up_spin_2[0]=up_spin_arg_2;}
    if (down_spin_arg_2 != '0') {down_spin_2[0]=down_spin_arg_2;}
  
  strcat(chstr_text, (char *) "\ncomparing spin state 2:\n\0");


    for (int index_y_=0;  index_y_ < Ly;  index_y_++)
    {
    int index_y=Ly-index_y_-1;

      for (int index_x=0;  index_x < Lx; index_x++)
      {  
        if (spin_array_1[index_y*Lx+index_x] == spin_array_2[index_y*Lx+index_x])
        {
          if (spin_array_2[index_y*Lx+index_x] == false) {strcat(chstr_text, down_spin);} else {strcat(chstr_text, up_spin);}
        }
        else
        {
          if (spin_array_2[index_y*Lx+index_x] == false) {strcat(chstr_text, down_spin_2);} else {strcat(chstr_text, up_spin_2);}
        }
                
        if (index_x < Lx-1) {strcat(chstr_text, hor_bond);}
      }

    strcat(chstr_text, (char *) "\n\0");


      if (index_y_ < Ly-1)
      {
    
        for (int index_x=0;  index_x < Lx; index_x++)
        {
        strcat(chstr_text, ver_bond);  strcat(chstr_text, space);
        }

      strcat(chstr_text, (char *) "\n\0");
      }
  
    } 
    
    
  }

return strlen(chstr_text);

}




void print_to_screen_2D_spin_state_comparing_lattice(const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx)
{

char * chstr_text=0;
int text_size=print_to_chstr_2D_spin_state_comparing_lattice(chstr_text, spin_array_1, spin_array_2,  Ly,  Lx);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}





int print_to_file_2D_spin_state_comparing_lattice(char *chstr_filename, const bool * const spin_array_1, const bool * const spin_array_2,  int Ly,  int Lx,  bool newfile)
{

char * chstr_text=0;
int text_size=print_to_chstr_2D_spin_state_comparing_lattice(chstr_text, spin_array_1, spin_array_2,  Ly,  Lx);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}






int  write_2D_spin_config_array_to_chstr(char * & chstr_text, bool *spin_array, int Ly, int Lx, int subLy, int subLx, int start_y, int start_x, int dy, int dx, bool print_comment)
{

  if (chstr_text) {return 0;}
  if (spin_array) {return 0;}

  if ((Ly < 1) || (Lx < 1)) {return 0;}
  if (subLy < 0) {subLy=Ly;}
  if (subLx < 0) {subLx=Lx;}
  if (subLy > Ly) {return 0;}
  if (subLx > Lx) {return 0;}

  if ((start_y < 0) || (start_x < 0)) {return 0;}

int end_y=start_y+subLy-1, end_x=start_x+subLx-1;

  if (end_y > Ly-1) {end_y=Ly-1;}
  if (end_x > Lx-1) {end_x=Lx-1;}

subLy=end_y-start_y+1;  subLx=end_x-start_x+1;

  if (dy < 0) {dy=0;}
  if (dx < 0) {dx=0;}



//  format data for writing
//
//  #  spin configuration for
//  ?????
//  #  spin configuration:
//  # s_y   s_x     spin
//  0       0       -1
//  0       1       -1
//  1       0       +1
//  1       1       -1

string string_var;

  if (print_comment == true) 
  {
  string_var.append((char *)  "#  spin configuration for\n");
  string_var.append((char *)  "?\n");
  string_var.append((char *)  "#  spin configuration:\n");
  string_var.append((char *)  "# s_y\ts_x\tspin\n");
  }


  {
  int index_y_loc=0, index_x_loc=0;
 
    for (int index_y=start_y;  index_y <= end_y;  index_y++)
    {
      for (int index_x=start_x;  index_x <= end_x;  index_x++)
      {
      index_y_loc=index_y;  index_y_loc=index_y % Ly;
      index_x_loc=index_x;  index_x_loc=index_x % Lx;

      add_int_to_string(& string_var, index_y_loc+dy);
      string_var.append((char *)  "\t");
      add_int_to_string(& string_var, index_x_loc+dx);
      string_var.append((char *)  "\t");
      int spin_value=-1; 
        if (spin_array[index_y_loc*Ly+index_x_loc] > 0) {spin_value=+1;}
        if (spin_value == 0) {string_var.append((char *)  " ");}  else {  if (spin_value > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, spin_value);
      string_var.append((char *)  "\n");
      }
    }
  }


chstr_text=new char[string_var.length()+1];
chstr_text[0]='\0';
strcpy(chstr_text, string_var.c_str());

//  for (int i=0;  i < string_var.length(); i++)
//  {chstr_text[i]=string_var.c_str()[i];}
//chstr_text[string_var.length()]='\0';

string_var.clear();


return strlen(chstr_text);

}



int  write_2D_spin_config_array_to_screen(bool *spin_array, int Ly, int Lx, int subLy, int subLx, int start_y, int start_x,  int dy, int dx, bool print_comment)
{


char * chstr_text=0;
int text_size=write_2D_spin_config_array_to_chstr(chstr_text, spin_array, Ly, Lx, subLy, subLx, start_y, start_x,  dy, dx, print_comment);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return 0;

}



int  write_spin_config_array_to_file(char *filename, bool *spin_array, int Ly, int Lx, int subLy=-1, int subLx=-1, int start_y=0, int start_x=0,  int dy=0, int dx=0, bool print_comment=true)
{

char * chstr_text=0;
int chstr_text_size=write_2D_spin_config_array_to_chstr(chstr_text, spin_array, Ly, Lx, subLy, subLx, start_y, start_x,  dy, dx, print_comment);
int file_size=0;

  if ((chstr_text_size > 0) && (chstr_text))  {file_size=write_chstr_to_file(chstr_text, filename, true);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;


}



int extract_2D_spin_config_array_from_file(char *filename, bool * const spin_array, int Ly, int Lx, int subLy, int subLx, int start_y, int start_x,  int dy, int dx)
{

//  format data for writing
//
//  #  spin configuration for
//  ?????
//  #  spin configuration:
//  # s_y   s_x     spin
//  0       0       -1
//  0       1       -1
//  1       0       +1
//  1       1       -1


  if (filename == 0) {return 0;}
  if ((Ly < 1) || (Lx < 1)) {return 0;}
  if (spin_array == 0) {return 0;}
  if (subLy < 0) {subLy=Ly;}
  if (subLx < 0) {subLx=Lx;}

  if (subLy > Ly) {return 0;}
  if (subLx > Lx) {return 0;}

  if ((start_y < 0) || (start_x < 0)) {return 0;}

int end_y=start_y+subLy-1, end_x=start_x+subLx-1;

  if (end_y > Ly-1) {end_y=Ly-1;}
  if (end_x > Lx-1) {end_x=Lx-1;}

subLy=end_y-start_y+1;  subLx=end_x-start_x+1;

  if (dy < 0) {dy=0;}
  if (dx < 0) {dx=0;}

bool was_error=false;
int number_of_extracted_spins=0;



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


key_phr_line_number=sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, (char *) "#  spin configuration for", 0); 

  if (key_phr_line_number >= 0) 
  {ok=sr_line_position_and_len_in_chstr(chstr_text,   key_phr_line_number, 0, & line_end_pos_temp, 0);}
      

  if (ok == true)
  {
  key_phr_line_number_line_n =sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, (char *) "#", line_end_pos_temp+3);              

    if (key_phr_line_number_line_n >= 0)            
    {ok=sr_line_position_and_len_in_chstr(chstr_text,  key_phr_line_number_line_n, 0, & line_end_pos_temp, 0);}      //  so we know that hist data starts with:    key_phr_line_number_line_n+3
    else {ok=false;}

    if (ok == true)
    {
    key_phr_line_number_line_nn=sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, (char *) "#", line_end_pos_temp+3);  

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

//==DEBUG==//cout<<"DEBUG: "<<"  extr_table_startpos="<<extr_table_startpos<<"  extr_table_endpos="<<extr_table_endpos<<"  extr_table_len="<<extr_table_len<<endl;
char remember_smb=chstr_text[extr_table_startpos+extr_table_len];      //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!
chstr_text[extr_table_startpos+extr_table_len]='\0';                   //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!

char ** chstr_phrases=0;
int  amount_of_extracted_phrases=0;
int  max_phrase_len=0;
extract_phrases_from_chstr(chstr_text+extr_table_startpos, chstr_phrases,  amount_of_extracted_phrases,  max_phrase_len,  0,  1);

chstr_text[extr_table_startpos+extr_table_len]=remember_smb;          //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!

  if (chstr_text) 
  {delete[] chstr_text;  chstr_text=0;}
  if (amount_of_extracted_phrases < 1) 
  {ok=false; was_error=true;}
  
  if (amount_of_extracted_phrases % 3 != 0) 
  {ok=false; was_error=true;}



//  ----------------------------------------------------------------------------
//  putting to matr

  //if (part_was_missed) {*part_was_missed=true;}

  if (ok == true)
  {
  int amount_of_lines=amount_of_extracted_phrases/3;

  int extracted_index_y=0,  extracted_index_x=0,  spin_value=0;
  
  //array_size=amount_of_lines; 

  //spin1_index_y=new int[array_size];
  //spin1_index_x=new int[array_size];
  //  spin_array  already implemented

    //for (int i=0;  i < array_size;  i++)
    //{
    //spin1_index_y[i]=0;  spin1_index_x[i]=0;
    //  spin_array[i]=false;
    //}

    for (int i=0;  i < (amount_of_lines) && (ok == true);  i++)
    {
    bool put_data_ok=true;

      if (is_chstr_valid_num(chstr_phrases[i*3+0]) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_phrases[i*3+1]) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_phrases[i*3+2]) == false) {put_data_ok=false;}

      if (put_data_ok == true)
      {
      extracted_index_y=atoi(chstr_phrases[i*3+0]);  extracted_index_x=atoi(chstr_phrases[i*3+1]);
      spin_value=atoi(chstr_phrases[i*3+2]);

      //  treatment of extracted indeces
        if (extracted_index_y <    0) {extracted_index_y=Ly-(-extracted_index_y % Ly);}
        if (extracted_index_y > Ly-1) {extracted_index_y=extracted_index_y % Ly;}
        if (extracted_index_x <    0) {extracted_index_x=Lx-(-extracted_index_x % Lx);}
        if (extracted_index_x > Lx-1) {extracted_index_x=extracted_index_x % Lx;}

        if ((extracted_index_y < start_y) || (end_y < extracted_index_y)) {put_data_ok=false;  was_error=true;}
        if ((extracted_index_x < start_x) || (end_x < extracted_index_x)) {put_data_ok=false;  was_error=true;}

        if (put_data_ok == true)
        {
        int index=(extracted_index_y+dy)*Lx+extracted_index_x+dx;
          if ((index >= 0) && (index <= Ly*Lx-1))
          {  if (spin_value >= 0) {spin_array[index]=true;} else {spin_array[index]=false;}   number_of_extracted_spins++;}
          else
          {put_data_ok=false;  was_error=true;}
        }

        if (extracted_index_y < start_y) {extracted_index_y=start_y;}
      }      //  if (put_data_ok == true)

    }      //  for (int i=0;  i < (amount_of_extracted_phrases) && (ok == true);  i++) 

  }      //  if (ok == true)



//  ----------------------------------------------------------------------------
//  cleaning

  if (chstr_phrases) {free_chstr_array(chstr_phrases, amount_of_extracted_phrases);}
 

return number_of_extracted_spins;

}






// already implemented in rdk_io_formatting_ext
//bool write_mpf_array_to_chstr(mpf_t * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len,  int wished_len=-1);





// already implemented in rdk_io_formatting_ext
// bool extract_mpz_column_from_file(char * filename, mpz_t *& col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);




int write_2D_rect_lattice_info_to_chstr_v15(char * & text, int Ly, int Lx, unsigned long long int state, int cur_thrd_num, int thrd_amount, int J_plus_perc, int J_min_perc, int min_energ, bool pbc1fbc0)
{

  if (text) {return 0;}
  if (Ly < 1) {return 0;}
  if (Lx < 1) {return 0;} 


string string_var;

string_var.append((char *)  "#  dim\n");
add_int_to_string(& string_var, 2);
string_var.append((char *)  "\n");

string_var.append((char *)  "#  lattice type\n");
string_var.append((char *)  "rect\n");

string_var.append((char *)  "#  bond cond\n");
  if (pbc1fbc0 == true) {string_var.append((char *)  "pbc\n");} else {string_var.append((char *)  "fbc\n");}

string_var.append((char *)  "#  Ly\n");
add_int_to_string(& string_var, Ly);
string_var.append((char *)  "\n");

string_var.append((char *)  "#  Lx\n");
add_int_to_string(& string_var, Lx);
string_var.append((char *)  "\n");

  if (cur_thrd_num >= 0)
  {
  string_var.append((char *)  "#  cur_thrd_num\n");
  add_int_to_string(& string_var, cur_thrd_num);
  string_var.append((char *)  "\n");
  }

  if (thrd_amount >= 0)
  {
  string_var.append((char *)  "#  amount of threads\n");
  add_int_to_string(& string_var, thrd_amount);
  string_var.append((char *)  "\n");
  }

  if (state != ~0)
  {
  string_var.append((char *)  "#  last completed row main state\n");
  add_ullint_to_string(& string_var, state);  //  else {string_var.append((char *)  "?\n");}
  string_var.append((char *)  "\n");
  }

  {
  char * chstr_temp=0;
  char elem='*';
  int size=0;
  string_var.append((char *)  "#  lattice and bond structure\n");
    if (pbc1fbc0 == true) {size=print_spin2D_lattice_with_bonds_to_chstr(chstr_temp, Ly, Lx, elem, 7);} else {size=print_spin2D_lattice_with_bonds_to_chstr(chstr_temp, Ly, Lx, elem, 0);}
    if (size > 0) {string_var.append(chstr_temp);}
    if (chstr_temp) {delete[] chstr_temp;  chstr_temp=0;}
  }

  if (min_energ != -999999999)
  {
  string_var.append((char *)  "#  min energy\n");
  add_int_to_string(& string_var, min_energ);  // else {string_var.append((char *)  "?\n");}
  string_var.append((char *)  "\n");
  }

  if (J_plus_perc >= 0)
  {
  string_var.append((char *)  "#  J_plus_perc\n");
  add_int_to_string(& string_var, J_plus_perc);
  string_var.append((char *)  "\n");
  }

  if (J_min_perc >= 0)
  {
  string_var.append((char *)  "#  J_min_perc\n");
  add_int_to_string(& string_var, J_min_perc);
  string_var.append((char *)  "\n");
  }


text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {text[i]=string_var.c_str()[i];}

text[string_var.length()]='\0';


return strlen(text);

}



int write_2D_rect_J_bond_table_to_chstr(char * & text, int array_size, int *spin1_index_y, int *spin1_index_x, int *spin2_index_y, int *spin2_index_x, bool *J_val, bool print_comment)
{

  if (array_size < 1) {return 0;}
  if (text) {return 0;}
  
  if (spin1_index_y == 0) {return 0;}
  if (spin1_index_x == 0) {return 0;}
  if (spin2_index_y == 0) {return 0;}
  if (spin2_index_x == 0) {return 0;}
  if (J_val         == 0) {return 0;}


string string_var;

  if (print_comment == true)
  {
  string_var.append("#  J table\n");
  string_var.append("#  s1_index_x\ts1_index_y\ts2_index_x\ts2_index_y\tbond_mean\n");
  }


int int_J_val=0;

  for (int i=0;  i < array_size;  i++)
  {
  add_int_to_string(& string_var, spin1_index_y[i]);
  string_var.append((char *)  "\t\t");
  add_int_to_string(& string_var, spin1_index_x[i]);
  string_var.append((char *)  "\t\t");
  add_int_to_string(& string_var, spin2_index_y[i]);
  string_var.append((char *)  "\t\t");
  add_int_to_string(& string_var, spin2_index_x[i]);
  string_var.append((char *)  "\t\t");
    if (J_val[i] == true) {int_J_val=1;} else {int_J_val=-1;}
    if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
  add_int_to_string(& string_var, (int) int_J_val);
  string_var.append((char *)  "\n");  
  }


text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {text[i]=string_var.c_str()[i];}

text[string_var.length()]='\0';


return strlen(text);

}






int write_2D_rect_J_bond_table_to_chstr(char * & text, int Ly, int Lx, bool *J_ver, bool *J_hor, char d_bonds_dcn, char u_bonds_ucn, char l_bonds_lcn, char r_bonds_rcn, int y_start, int y_end, int x_start, int x_end, int out_dx, int out_dy, bool print_comment)
{

  if (text) {return 0;}
  if ((Ly < 1) && (Lx < 1)) {return 0;}
  if (J_ver == 0) {return 0;}
  if (J_hor == 0) {return 0;}
  if (y_end < 0) {y_end=Ly-1;}
  if (x_end < 0) {x_end=Lx-1;}
  if (y_start < 0) {y_start=0;}
  if (x_start < 0) {x_start=0;}

  if (y_start > y_end) {return 0;}
  if (x_start > x_end) {return 0;}


string string_var;
int start_down_index_y=0, start_down_index_x=0; 

  if (print_comment == true)
  {
  string_var.append("#  J table\n");
  string_var.append("#  s1_index_y\ts1_index_x\ts2_index_y\ts2_index_x\tbond_mean\n");
  } 




int int_J_val=0;

// ver bonds
//
  {
  int J_ver_index_y=0, J_ver_index_x=0;

    for (int index_y=y_start;  index_y < y_end;  index_y++)
    {

      for (int index_x=x_start;  index_x <= x_end;  index_x++)
      {
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+1+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");

      J_ver_index_y=index_y+1;  J_ver_index_x=index_x;
        if (J_ver[Lx*J_ver_index_y+J_ver_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }
  }


// hor bonds
//
  {
  int J_hor_index_y=0, J_hor_index_x=0;

    for (int index_y=y_start;  index_y <= y_end;  index_y++)
    {

      for (int index_x=x_start;  index_x < x_end;  index_x++)
      {
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+1+out_dx);
      string_var.append((char *)  "\t\t");

      J_hor_index_y=index_y;  J_hor_index_x=index_x+1;
        if (J_hor[(Lx+1)*J_hor_index_y+J_hor_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

  }




//  - - - - - - - - - - - - - - - -


//  down ver bonds    |  |  |  |  |
//
  if ((d_bonds_dcn == 'd') || (d_bonds_dcn == 'c'))
  {
  int J_ver_index_y=0, J_ver_index_x=0;
  int index_y=0;

    if (d_bonds_dcn == 'd')
    {
    index_y=y_start;
      for (int index_x=x_start;  index_x <= x_end;  index_x++)
      {
      add_int_to_string(& string_var, index_y-1+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");

      J_ver_index_y=index_y;  J_ver_index_x=index_x;
        if (J_ver[Lx*J_ver_index_y+J_ver_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

    if ((d_bonds_dcn == 'c') && (Ly > 2))
    {
    index_y=0;
      for (int index_x=x_start;  index_x <= x_end;  index_x++)
      {
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, Ly-1);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");

      J_ver_index_y=Ly;  J_ver_index_x=index_x; 
        if (J_ver[Lx*J_ver_index_y+J_ver_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

  }      //  if ((d_bonds_dun == 'd') || (u_bonds_dun == 'u'))



//  up ver bonds    |  |  |  |  |
//
  if ((u_bonds_ucn == 'u') || (u_bonds_ucn == 'c'))
  {
  int J_ver_index_y=0, J_ver_index_x=0;
  int index_y=y_end;

    if (u_bonds_ucn == 'u')
    {
    index_y=y_end;
      for (int index_x=x_start;  index_x <= x_end;  index_x++)
      {
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+1+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");

      J_ver_index_y=index_y+1;  J_ver_index_x=index_x;
        if (J_ver[Lx*J_ver_index_y+J_ver_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

    if ((u_bonds_ucn == 'c') && (Ly > 2))
    {
    index_y=y_end;
      for (int index_x=x_start;  index_x <= x_end;  index_x++)
      {
      add_int_to_string(& string_var, 0);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");

      J_ver_index_y=0;  J_ver_index_x=index_x;        
        if (J_ver[Lx*J_ver_index_y+J_ver_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

  }      //  if ((d_bonds_dun == 'd') || (u_bonds_dun == 'u'))



//  left hor bonds    --  --  --  --  --
//
  if ((l_bonds_lcn == 'l') || (l_bonds_lcn == 'c'))
  {
  int J_hor_index_y=0, J_hor_index_x=0;
  int index_x=0;

    if (l_bonds_lcn == 'l')
    {
    index_x=x_start;
      for (int index_y=y_start;  index_y <= y_end;  index_y++)
      {
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x-1+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");

      J_hor_index_y=index_y;  J_hor_index_x=index_x;
        if (J_hor[(Lx+1)*J_hor_index_y+J_hor_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

    if ((l_bonds_lcn == 'c') && (Lx > 2))
    {
    index_x=x_start;
      for (int index_y=y_start;  index_y <= y_end;  index_y++)
      {
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, Lx-1);
      string_var.append((char *)  "\t\t");

      J_hor_index_y=index_y;  J_hor_index_x=Lx;
        if (J_hor[(Lx+1)*J_hor_index_y+J_hor_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

  }      //  if ((d_bonds_dun == 'd') || (u_bonds_dun == 'u'))




//  right hor bonds    --  --  --  --  --
//
  if ((r_bonds_rcn == 'r') || (r_bonds_rcn == 'c'))
  {
  int J_hor_index_y=0, J_hor_index_x=0;
  int index_x=0;

    if (r_bonds_rcn == 'r')
    {
    index_x=x_end;
      for (int index_y=y_start;  index_y <= y_end;  index_y++)
      {
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+1+out_dx);
      string_var.append((char *)  "\t\t");

      J_hor_index_y=index_y;  J_hor_index_x=index_x+1;
        if (J_hor[(Lx+1)*J_hor_index_y+J_hor_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

    if ((l_bonds_lcn == 'c') && (Lx > 2))
    {
    index_x=x_end;
      for (int index_y=y_start;  index_y <= y_end;  index_y++)
      {
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, 0);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_y+out_dy);
      string_var.append((char *)  "\t\t");
      add_int_to_string(& string_var, index_x+out_dx);
      string_var.append((char *)  "\t\t");

      J_hor_index_y=index_y;  J_hor_index_x=Lx;
        if (J_hor[(Lx+1)*J_hor_index_y+J_hor_index_x] == true) {int_J_val=1;} else {int_J_val=-1;}
        if (int_J_val == 0) {string_var.append((char *)  " ");}  else {  if (int_J_val > 0) {string_var.append((char *)  "+");}  }
      add_int_to_string(& string_var, (int) int_J_val);
      string_var.append((char *)  "\n");    
      }
    }

  }      //  if ((d_bonds_dun == 'd') || (u_bonds_dun == 'u'))


text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {text[i]=string_var.c_str()[i];}

text[string_var.length()]='\0';


return strlen(text);

}




int write_2D_rect_lattice_info_to_screen_v15(int Ly, int Lx, unsigned long long int state, int cur_thrd_num, int thrd_amount, int J_plus_perc, int J_min_perc, int min_energ, bool pbc1fbc0)
{

char * chstr_text=0;
int text_size=write_2D_rect_lattice_info_to_chstr_v15(chstr_text, Ly, Lx, state, cur_thrd_num, thrd_amount, J_plus_perc, J_min_perc, min_energ, pbc1fbc0);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return 0;

}




int write_2D_rect_J_bond_table_to_screen(int array_size, int *spin1_index_y, int *spin1_index_x, int *spin2_index_y, int *spin2_index_x, bool *J_val, bool print_comment)
{

char * chstr_text=0;
int text_size=write_2D_rect_J_bond_table_to_chstr(chstr_text, array_size, spin1_index_y, spin1_index_x, spin2_index_y, spin2_index_x, J_val, print_comment);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return 0;

}



int write_2D_rect_J_bond_table_to_screen(int Ly, int Lx, bool *J_ver, bool *J_hor, char d_bonds_dcn, char u_bonds_ucn, char l_bonds_lcn, char r_bonds_rcn, int y_start, int y_end, int x_start, int x_end, int out_dx, int out_dy, bool print_comment)
{

char * chstr_text=0;
int text_size=write_2D_rect_J_bond_table_to_chstr(chstr_text, Ly, Lx, J_ver, J_hor, d_bonds_dcn, u_bonds_ucn, l_bonds_lcn, r_bonds_rcn, y_start, y_end, x_start, x_end, out_dx, out_dy, print_comment);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return 0;

}



int write_2D_rect_lattice_info_to_file_v15(char * filename, int Ly, int Lx, unsigned long long int state, int cur_thrd_num, int thrd_amount, int J_plus_perc, int J_min_perc, int min_energ, bool pbc1fbc0, bool newfile)
{

char * chstr_text=0;
int chstr_text_size=write_2D_rect_lattice_info_to_chstr_v15(chstr_text, Ly, Lx, state, cur_thrd_num, thrd_amount, J_plus_perc, J_min_perc, min_energ, pbc1fbc0);
int file_size=0;

  if ((chstr_text_size > 0) && (chstr_text))  {file_size=write_chstr_to_file(chstr_text, filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}




int write_2D_rect_J_bond_table_to_file(char * filename, int array_size, int *spin1_index_y, int *spin1_index_x, int *spin2_index_y, int *spin2_index_x, bool *J_val, bool print_comment, bool newfile)
{

char * chstr_text=0;
int chstr_text_size=write_2D_rect_J_bond_table_to_chstr(chstr_text, array_size, spin1_index_y, spin1_index_x, spin2_index_y, spin2_index_x, J_val, print_comment);
int file_size=0;

  if ((chstr_text_size > 0) && (chstr_text))  {file_size=write_chstr_to_file(chstr_text, filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}



int write_2D_rect_J_bond_table_to_file(char * filename, int Ly, int Lx, bool *J_ver, bool *J_hor, char d_bonds_dcn, char u_bonds_ucn, char l_bonds_lcn, char r_bonds_rcn, int y_start, int y_end, int x_start, int x_end, int out_dx, int out_dy, bool print_comment, bool newfile)
{

char * chstr_text=0;
int chstr_text_size=write_2D_rect_J_bond_table_to_chstr(chstr_text, Ly, Lx, J_ver, J_hor, d_bonds_dcn, u_bonds_ucn, l_bonds_lcn, r_bonds_rcn, y_start, y_end, x_start, x_end, out_dx, out_dy, print_comment);
int file_size=0;

  if ((chstr_text_size > 0) && (chstr_text))  {file_size=write_chstr_to_file(chstr_text, filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}



bool extract_2D_rect_fbc_lattice_serv_info_from_file_v15(char * filename, int *Ly, int *Lx, unsigned long long int *state, int *cur_thrd_num, int *thrd_amount, int *J_plus_perc, int *J_min_perc, int *min_energ)
{

  if (filename == 0) {return false;}


//  reading file
char * chstr_text=0;
//  func format: int write_file_to_chstr(char *filename, char *& text,  char *chstr_forbid_1st_line_symb_list,  int start_line_number=0,  bool not_take_empty_lines=true);
write_file_to_chstr(filename, chstr_text,  0);
char * chstr_extracting_phrase=0;
bool ok=false;
int returning_len=0;


//  -----------------
//  #  Ly
//  6

  if (Ly)
  {
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  Ly",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;}
    if (ok == true)
    {       
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *Ly=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }


//  -----------------
//  #  Lx
//  6

  if (Lx)
  {
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  Lx",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;}
    if (ok == true)
    {       
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *Lx=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }


//  -----------------
//  #  last completed row main state"
//  1636111

  if (state)
  {
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  last completed row main state",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;}
    if (ok == true)
    {       
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *state=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }


//  -----------------
//  #  cur thread number
//  -1

  if (cur_thrd_num)
  {
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  cur thread number",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;}
    if (ok == true)
    {       
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *cur_thrd_num=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }


//  -----------------
//  #  amount of threads
//  4

  if (thrd_amount)
  {
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  amount of threadsr",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;}
    if (ok == true)
    {       
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *thrd_amount=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }


//  -----------------
//  #  J+ percent
//  20

  if (J_plus_perc)
  {
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  J_plus_perc",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;}
    if (ok == true)
    {       
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *J_plus_perc=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }


//  -----------------
//  #  J- percent
//  80

  if (J_min_perc)
  { 
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  J_min_perc",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;}
    if (ok == true)
    { 
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *J_min_perc=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }


//  -----------------
//  #  min energy
//  -125

  if (min_energ)
  {
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  min energy",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;} 
    if (ok == true)
    {       
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *min_energ=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }



  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}


return true;

}







bool extract_2D_rect_lattice_J_table_from_file_dyn(char * filename, int Ly, int Lx, bool *& J_ver, bool *& J_hor, char d_bonds_dcn, char u_bonds_ucn, char l_bonds_lcn, char r_bonds_rcn, int y_start, int y_end, int x_start, int x_end, int index_dy, int index_dx)
{    

  if (filename == 0) {return false;}            
  if ((Ly < 1) || (Lx < 1)) {return false;}       
  if ((J_ver) || (J_hor)) {return false;}          

  if (y_end < 0) {y_end=Ly-1;}        
  if (x_end < 0) {x_end=Lx-1;}                
  if (y_start < 0) {y_start=0;}            
  if (x_start < 0) {x_start=0;}      

  if (index_dy < 0) {index_dy=0;}           
  if (index_dx < 0) {index_dx=0;}        

bool was_error=false; 



char * chstr_text=0;
int  file_size=write_file_to_chstr(filename, chstr_text, 0); 
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


key_phr_line_number=sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, (char *) "#  J table", 0); 

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
//
char *chstr_phrases=0;
int *len_array=0;
int *start_pos_array=0;
int max_extracted_phrase_len=0;

int  amount_of_extracted_phrases=extract_phrases_from_chstr(chstr_text, chstr_phrases, len_array, start_pos_array, & max_extracted_phrase_len,  0,  1,  extr_table_startpos, extr_table_endpos); 

  if (amount_of_extracted_phrases < 1) {ok=false;}      
  if (amount_of_extracted_phrases % 5 != 0) {ok=false;}     
  if (max_extracted_phrase_len < 1) {ok=false;}          
  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}   



//  ----------------------------------------------------------------------------
//  putting to matr

  if (ok == true)
  {
  int line_amount=amount_of_extracted_phrases/5;
  int Ly_loc=y_end-y_start+1, Lx_loc=x_end-x_start+1;
  int max_index_J_ver=(Ly_loc+1)*Lx_loc-1;
  int max_index_J_hor=Ly_loc*(Lx_loc+1)-1;
  J_ver=new bool[(Ly_loc+1)*Lx_loc];
  J_hor=new bool[Ly_loc*(Lx_loc+1)];

    for (int i=0;  i < (Ly+1)*Lx;  i++)
    {J_ver[i]=false;}

    for (int i=0;  i < Ly*(Ly+1);  i++)
    {J_hor[i]=false;}


  int spin1_index_y=0,      spin1_index_x=0,      spin2_index_y=0,      spin2_index_x=0;
  int loc_spin1_index_y=0,  loc_spin1_index_x=0,  loc_spin2_index_y=0,  loc_spin2_index_x=0;
  bool spin_value=0;
  int amount_of_lines=amount_of_extracted_phrases/5;
  char  *chstr_spin1_index_y=0,  *chstr_spin1_index_x=0,  *chstr_spin2_index_y=0, *chstr_spin2_index_x=0;
  char  *chstr_J_value=0;
  chstr_spin1_index_y=new char[max_extracted_phrase_len+8];  chstr_spin1_index_y[0]='\0';     chstr_spin1_index_x=new char[max_extracted_phrase_len+8];  chstr_spin1_index_x[0]='\0';  
  chstr_spin2_index_y=new char[max_extracted_phrase_len+8];  chstr_spin2_index_y[0]='\0';     chstr_spin2_index_x=new char[max_extracted_phrase_len+8];  chstr_spin2_index_x[0]='\0';
  chstr_J_value=new char[max_extracted_phrase_len+8];        chstr_J_value[0]='\0';

    for (int i=0;  i < (line_amount) && (ok == true);  i++)
    {
    bool put_data_ok=true;
    int index=0;
 
    chstr_spin1_index_y[0]='\0';  chstr_spin1_index_x[0]='\0';  chstr_spin2_index_y[0]='\0';  chstr_spin2_index_x[0]='\0';        chstr_J_value[0]='\0';

   //=DEBUG=//  cout<<"DEBUG: tot_len="<<strlen(chstr_phrases)<<"  start_pos_array[i*5+0]="<<start_pos_array[i*5+0]<<"  "<<"len_array[i*5+0]="<<len_array[i*5+0]<<endl;  //=DEBUG=// 
   //=DEBUG=//  cout<<"DEBUG: tot_len="<<strlen(chstr_phrases)<<"  start_pos_array[i*5+1]="<<start_pos_array[i*5+1]<<"  "<<"len_array[i*5+1]="<<len_array[i*5+1]<<endl;  //=DEBUG=// 
   //=DEBUG=//  cout<<"DEBUG: tot_len="<<strlen(chstr_phrases)<<"  start_pos_array[i*5+2]="<<start_pos_array[i*5+2]<<"  "<<"len_array[i*5+2]="<<len_array[i*5+2]<<endl;  //=DEBUG=// 
   //=DEBUG=//  cout<<"DEBUG: tot_len="<<strlen(chstr_phrases)<<"  start_pos_array[i*5+3]="<<start_pos_array[i*5+3]<<"  "<<"len_array[i*5+3]="<<len_array[i*5+3]<<endl;  //=DEBUG=// 
   //=DEBUG=//  cout<<"DEBUG: tot_len="<<strlen(chstr_phrases)<<"  start_pos_array[i*5+4]="<<start_pos_array[i*5+4]<<"  "<<"len_array[i*5+4]="<<len_array[i*5+4]<<endl;  //=DEBUG=// 

      if (copy_replace(chstr_phrases,  start_pos_array[i*5+0],  len_array[i*5+0], chstr_spin1_index_y,  0, true) == false)  {put_data_ok=false;}    //  cout<<"DEBUG: chstr_spin1_index_y="<<chstr_spin1_index_y<<"  ";
      if (copy_replace(chstr_phrases,  start_pos_array[i*5+1],  len_array[i*5+1], chstr_spin1_index_x,  0, true) == false)  {put_data_ok=false;}    //  cout<<"DEBUG: chstr_spin1_index_x="<<chstr_spin1_index_x<<"  ";
      if (copy_replace(chstr_phrases,  start_pos_array[i*5+2],  len_array[i*5+2], chstr_spin2_index_y,  0, true) == false)  {put_data_ok=false;}    //  cout<<"DEBUG: chstr_spin2_index_y="<<chstr_spin2_index_y<<"  ";
      if (copy_replace(chstr_phrases,  start_pos_array[i*5+3],  len_array[i*5+3], chstr_spin2_index_x,  0, true) == false)  {put_data_ok=false;}    //  cout<<"DEBUG: chstr_spin2_index_x="<<chstr_spin2_index_x<<"  ";
      if (copy_replace(chstr_phrases,  start_pos_array[i*5+4],  len_array[i*5+4], chstr_J_value,        0, true) == false)  {put_data_ok=false;}    //  cout<<"DEBUG: chstr_J_value="<<chstr_J_value<<"  "<<endl;

      if (is_chstr_valid_num(chstr_spin1_index_y) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_spin1_index_x) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_spin2_index_y) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_spin2_index_y) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_J_value)       == false) {put_data_ok=false;}

    spin1_index_y=atoi(chstr_spin1_index_y);  spin1_index_x=atoi(chstr_spin1_index_x);
    spin2_index_y=atoi(chstr_spin2_index_y);  spin2_index_x=atoi(chstr_spin2_index_x); 
      if (atoi(chstr_J_value) > 0) {spin_value=true;} else {spin_value=false;}

    //  all values are extracted,    now we need sort it and put to J_ver and J_hor

      if (spin2_index_y < spin1_index_y) {swap(spin2_index_y, spin1_index_y);}
      if (spin2_index_x < spin1_index_x) {swap(spin2_index_x, spin1_index_x);}

    loc_spin1_index_y=spin1_index_y+index_dy;  loc_spin1_index_x=spin1_index_x+index_dx;  loc_spin2_index_y=spin2_index_y+index_dy;  loc_spin2_index_x=spin2_index_x+index_dx;
    
      //if ((loc_spin1_index_y < 0) || (loc_spin1_index_x < 0) || (loc_spin2_index_y < 0) || (loc_spin2_index_x < 0)) {put_data_ok=false;}
      //if ((loc_spin1_index_y > Ly_loc) || (loc_spin1_index_x > Lx_loc) || (loc_spin2_index_y > Ly_loc) || (loc_spin2_index_x > Lx_loc)) {put_data_ok=false;}
 
      if ((loc_spin1_index_y < y_start-1) || (loc_spin2_index_y < y_start-1) || (loc_spin1_index_y > y_end+1) || (loc_spin2_index_y > y_end+1)) {put_data_ok=false;}
      if ((loc_spin1_index_x < x_start-1) || (loc_spin2_index_x < x_start-1) || (loc_spin1_index_x > x_end+1) || (loc_spin2_index_x > x_end+1)) {put_data_ok=false;}
  
  
 
      if (put_data_ok == true) 
      { 

        if (((spin2_index_y-spin1_index_y == 1) || (spin2_index_y-spin1_index_y == Ly_loc-1)) && (spin2_index_x-spin1_index_x == 0))
        {  
        // ver bond
          if ((y_start <= spin1_index_y) && (spin2_index_y <= y_end) && (x_start <= spin1_index_x) && (spin2_index_x <= x_end))
          {
          //  inner spins        
            if (spin2_index_y-spin1_index_y == 1) 
            {index=(loc_spin1_index_y+1)*Lx_loc+loc_spin1_index_x;  if ((index >= 0) && (index <= max_index_J_ver)) {J_ver[index]=spin_value;} else {was_error=true;}  }  //  adjacent spins      
            //**//else {
            if ((spin2_index_y-spin1_index_y != 1) || ((spin2_index_y-spin1_index_y == 1) && (Ly_loc <= 2)))
            {
            // ver pbc spins
              if (u_bonds_ucn == 'c') {index=(loc_spin2_index_y+1)*Lx_loc+loc_spin1_index_x;  if ((index >= 0) && (index <= max_index_J_ver)) {J_ver[index]=spin_value;} else {was_error=true;}  }
              if (d_bonds_dcn == 'c') {index=loc_spin1_index_x;                     if ((index >= 0) && (index <= max_index_J_ver)) {J_ver[index]=spin_value;} else {was_error=true;}  }
            }
          }
          else
          {
          //  edge spins
            if ((y_end == spin1_index_y) && (spin2_index_y == y_end+1) && (x_start <= spin1_index_x) && (spin2_index_x <= x_end) && (u_bonds_ucn == 'u'))      //  up edge
            {index=(loc_spin2_index_y)*Lx_loc+loc_spin1_index_x;    if ((index >= 0) && (index <= max_index_J_ver)) {J_ver[index]=spin_value;} else {was_error=true;}  }
          //  edge spins
            if ((y_start-1 == spin1_index_y) && (spin2_index_y == y_start) && (x_start <= spin1_index_x) && (spin2_index_x <= x_end) && (d_bonds_dcn == 'd'))      //  down edge
            {index=(spin2_index_y)*Lx_loc+loc_spin1_index_x;     if ((index >= 0) && (index <= max_index_J_ver)) {J_ver[index]=spin_value;} else {was_error=true;}  }
          }
        }

        if (((spin2_index_x-spin1_index_x == 1) || (spin2_index_x-spin1_index_x == Lx_loc-1)) && (spin2_index_y-spin1_index_y == 0))
        {
        //  hor bond
          if ((y_start <= spin1_index_y) && (spin2_index_y <= y_end) && (x_start <= spin1_index_x) && (spin2_index_x <= x_end))
          {
          //  inner spins        
            if (spin2_index_x-spin1_index_x == 1)  
            {index=loc_spin1_index_y*(Lx_loc+1)+loc_spin1_index_x+1;  if ((index >= 0) && (index <= max_index_J_hor)) {J_hor[index]=spin_value;} else {was_error=true;}  }    //  adjacent spins      
            //**//else {
            if ((spin2_index_x-spin1_index_x != 1) || ((spin2_index_x-spin1_index_x == 1) && (Lx_loc <= 2)))
            {
            // ver pbc spins
              if (l_bonds_lcn == 'c') {index=loc_spin1_index_y*(Lx_loc+1)+loc_spin1_index_x  ;    if ((index >= 0) && (index <= max_index_J_hor)) {J_hor[index]=spin_value;} else {was_error=true;}  }
              if (r_bonds_rcn == 'c') {index=loc_spin1_index_y*(Lx_loc+1)+loc_spin2_index_x+1;    if ((index >= 0) && (index <= max_index_J_hor)) {J_hor[index]=spin_value;} else {was_error=true;}  }
            }
          }
          else
          {
          //  edge spins
            if ((x_end == spin1_index_x) && (spin2_index_x == x_end+1) && (y_start <= spin1_index_y) && (spin2_index_y <= y_end) && (r_bonds_rcn == 'r'))      //  right edge
            {index=loc_spin1_index_y*(Lx_loc+1)+loc_spin2_index_x;    if ((index >= 0) && (index <= max_index_J_hor)) {J_hor[index]=spin_value;} else {was_error=true;}  }
          //  edge spins
            if ((x_start-1 == spin1_index_x) && (spin2_index_x == x_start) && (y_start <= spin1_index_y) && (spin2_index_y <= y_end) && (l_bonds_lcn == 'l'))      //  left edge
            {index=loc_spin1_index_y*(Lx_loc+1)+loc_spin2_index_x;    if ((index >= 0) && (index <= max_index_J_hor)) {J_hor[index]=spin_value;} else {was_error=true;}  }
          } 
        }
      
      }      //  if (put_data_ok == true)

    }      //  for (int i=0;  i < (amount_of_extracted_phrases) && (ok == true);  i++) 

    if (chstr_spin1_index_y) {delete[] chstr_spin1_index_y;  chstr_spin1_index_y=0;}        if (chstr_spin1_index_x) {delete[] chstr_spin1_index_x;  chstr_spin1_index_x=0;}
    if (chstr_spin2_index_y) {delete[] chstr_spin2_index_y;  chstr_spin2_index_y=0;}        if (chstr_spin2_index_x) {delete[] chstr_spin2_index_x;  chstr_spin2_index_x=0;}
    if (chstr_J_value      ) {delete[] chstr_J_value;        chstr_J_value=0;}

  }      //  if (ok == true)


//  ----------------------------------------------------------------------------
//  cleaning

  if (chstr_phrases  ) {delete[] chstr_phrases;    chstr_phrases=0;}
  if (len_array      ) {delete[] len_array;        len_array=0;}
  if (start_pos_array) {delete[] start_pos_array;  start_pos_array=0;}


  if (ok == false) 
  {
    if (J_ver) {delete[] J_ver;  J_ver=0;}
    if (J_hor) {delete[] J_hor;  J_hor=0;}

  return false;
  }


return true;

}

  





bool extract_2D_rect_lattice_J_table_from_file(char *filename, int & array_size, int *& spin1_index_y, int *& spin1_index_x, int *& spin2_index_y, int *& spin2_index_x, bool *& J_val, bool *part_was_missed)
{

  if (filename == 0) {return false;}
  if ((spin1_index_y) || (spin2_index_y) || (spin1_index_x) || (spin2_index_x) || (J_val)) {return false;}



char * chstr_text=0;
int file_size=write_file_to_chstr(filename, chstr_text,  0);
  if (file_size < 1) {return 0;} 
bool ok=true;    
array_size=0;  





//  ----------------------------------------------------------------------------
//  defining:       extr_table_startpos,   extr_table_endpos, extr_table_len
//
int extr_table_startpos=0, extr_table_endpos=0, extr_table_len=0;

int  key_phr_line_number=0;                 //  n - next
int  key_phr_line_number_line_n=0;          //  n - next
int  key_phr_line_number_line_nn=0;         //  n - next

int line_end_pos_temp=0;


key_phr_line_number=sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, (char *) "#  J table", 0); 

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
//==DEBUG==//  for (int i=extr_table_startpos;  i <= extr_table_endpos;  i++) {cout<<chstr_text[i];} cout<<endl;  //==//DEBUG
//  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}      //==// DEBUG

//  done ! ! !




//  ----------------------------------------------------------------------------
//  copying text chunk analizing and extracting phrases
//  format int extract_phrases_from_chstr(char *src_text, char ** & chstr_phrases, int & amount_of_extracted_phrases,  int & phrase_len,  int start_phrase_num,  int col_count=-1);

//==DEBUG==//cout<<"DEBUG: "<<"  extr_table_startpos="<<extr_table_startpos<<"  extr_table_endpos="<<extr_table_endpos<<"  extr_table_len="<<extr_table_len<<endl;
char remember_smb=chstr_text[extr_table_startpos+extr_table_len];      //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!
chstr_text[extr_table_startpos+extr_table_len]='\0';                   //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!

char ** chstr_phrases=0;
int  amount_of_extracted_phrases=0;
int  max_phrase_len=0;
extract_phrases_from_chstr(chstr_text+extr_table_startpos, chstr_phrases,  amount_of_extracted_phrases,  max_phrase_len,  0,  1);

chstr_text[extr_table_startpos+extr_table_len]=remember_smb;          //  ! ! ! ! DANGER. CLOSE TO MEM LEAK! BE CAREFULLY!

  if (chstr_text) 
  {delete[] chstr_text;  chstr_text=0;}
  if (amount_of_extracted_phrases < 1) 
  {ok=false;}
  
  if (amount_of_extracted_phrases % 5 != 0) 
  {ok=false;}



//  ----------------------------------------------------------------------------
//  putting to matr

  if (part_was_missed) {*part_was_missed=true;}

  if (ok == true)
  {

  int amount_of_lines=amount_of_extracted_phrases/5;
  array_size=amount_of_lines;

  spin1_index_y=new int[array_size];
  spin1_index_x=new int[array_size];
  spin2_index_y=new int[array_size];
  spin2_index_x=new int[array_size];
  J_val=new bool[array_size];

    for (int i=0;  i < array_size;  i++)
    {
    spin1_index_y[i]=0;  spin1_index_x[i]=0;
    spin2_index_y[i]=0;  spin2_index_x[i]=0;
    J_val[i]=false;
    }

    for (int i=0;  i < (amount_of_lines) && (ok == true);  i++)
    {
    bool put_data_ok=true;

      if (is_chstr_valid_num(chstr_phrases[i*5+0]) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_phrases[i*5+1]) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_phrases[i*5+2]) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_phrases[i*5+3]) == false) {put_data_ok=false;}
      if (is_chstr_valid_num(chstr_phrases[i*5+4]) == false) {put_data_ok=false;}

      if (put_data_ok == true)
      {
      spin1_index_y[i]=atoi(chstr_phrases[i*5+0]);  spin1_index_x[i]=atoi(chstr_phrases[i*5+1]);
      spin2_index_y[i]=atoi(chstr_phrases[i*5+2]);  spin2_index_x[i]=atoi(chstr_phrases[i*5+3]);
        if (atoi(chstr_phrases[i*5+4]) > 0) {J_val[i]=true;} else {J_val[i]=false;}
      }      //  if (put_data_ok == true)

      if (put_data_ok == false)  {  if (part_was_missed) {*part_was_missed=false;}  }

    }      //  for (int i=0;  i < (amount_of_extracted_phrases) && (ok == true);  i++) 

  }      //  if (ok == true)



//  ----------------------------------------------------------------------------
//  cleaning

  if (chstr_phrases) {free_chstr_array(chstr_phrases, amount_of_extracted_phrases);}
 
  if (ok == false) 
  {
  array_size=0;

    if (spin1_index_y) {delete[] spin1_index_y;  spin1_index_y=0;}
    if (spin1_index_x) {delete[] spin1_index_x;  spin1_index_x=0;}
    if (spin2_index_y) {delete[] spin2_index_y;  spin2_index_y=0;}
    if (spin2_index_x) {delete[] spin2_index_x;  spin2_index_x=0;}
    if (J_val        ) {delete[] J_val;          J_val=0;        }

  return false;
  }


return true;

}  




int print_to_chstr_clusters_id_for_2d_lattice(char *& chstr_text, int Ly,  int Lx,  int *cluster_id)
{

  if (cluster_id == 0) {return 0;}
  if (chstr_text) {return 0;}

string string_var;

string_var.append((char *)  "\n\n");    
  
  
char chstr_num[64];  chstr_num[0]='\0';
sprintf(chstr_num, "%d", Ly*Lx);
int cluster_id_len=strlen(chstr_num);
  if (cluster_id_len % 2 == 0) {cluster_id_len++;}
  
  for (int y=Ly-1;  y >= 0;  y--)
  {
  
    for (int x=0;  x < Lx;  x++)  
    {
    chstr_num[0]='\0';
    sprintf(chstr_num, "% *d", cluster_id_len, cluster_id[y*Lx+x]);
    string_var.append(chstr_num);
      if (x < Lx-1) {string_var.append((char *)  "-");} else {string_var.append((char *)  "\n");}
    }
    
    if (y > 0)
    {
      {
      for (int x=0;  x < Lx;  x++)  
      {
        for (int i=0;  i < (cluster_id_len-1)/2; i++) {string_var.append((char *)  " ");}      string_var.append((char *)  "|");    for (int i=0;  i < (cluster_id_len-1)/2; i++) {string_var.append((char *)  " ");}  
          if (x < Lx-1) {string_var.append((char *)  " ");} else {string_var.append((char *)  "\n");}
        }
      }
    }
    
  }  
      
      

chstr_text=new char[string_var.length()+1];

  for (int i=0;  i < string_var.length(); i++)
  {chstr_text[i]=string_var.c_str()[i];}

chstr_text[string_var.length()]='\0';
string_var.clear();

return strlen(chstr_text);
      
}





void print_to_screen_clusters_id_for_2d_lattice(int Ly,  int Lx,  int *cluster_id)
{

char * chstr_text=0;
int text_size=print_to_chstr_clusters_id_for_2d_lattice(chstr_text, Ly,  Lx,  cluster_id);

  if (text_size > 0) {cout<<chstr_text<<endl;}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

}




//  write spin configuration to file chstr_filename
int print_to_file_clusters_id_for_2d_lattice(char *chstr_filename,  int Ly,  int Lx,  int *cluster_id,  bool newfile)
{

  if ((Ly < 1) || (Lx < 1)) {return 0;}

char * chstr_text=0;
int text_size=print_to_chstr_clusters_id_for_2d_lattice(chstr_text, Ly,  Lx,  cluster_id);
int file_size=0;

  if (text_size > 0) {file_size=write_chstr_to_file(chstr_text, chstr_filename, newfile);}

  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

return file_size;

}










//  ============================     ============================    ============================     ============================
//
//  =============     ============================     ============================    ============================     ============================ 
//
//  =================     ============================    ============================     ============================






bool  st_char_str_lattice_console_draw :: set(int cells_amount_y_, int cells_amount_x_,  int cell_sz_in_y_,  int cell_sz_in_x_,  int node_sz_y_,  int node_sz_x_)
{

  if (chstr_lattice) {return false;}
  if (cells_amount_y_ < 0) {return false;}
  if (cells_amount_x_ < 0) {return false;}

cells_amount_y=cells_amount_y_;  cells_amount_x=cells_amount_x_;
cell_sz_in_y=cell_sz_in_y_;  cell_sz_in_x=cell_sz_in_x_;  
node_sz_y=node_sz_y_;   node_sz_x=node_sz_x_;

lattice_sz_y=cells_amount_y*cell_sz_in_y+(cells_amount_y+1)*node_sz_y;    lattice_sz_x=cells_amount_x*cell_sz_in_x+(cells_amount_x+1)*node_sz_x;    lattice_sz_incrrem_x=lattice_sz_x+1;
size=lattice_sz_y*lattice_sz_x+lattice_sz_y;      //  !  !  !  !  !  ! very important part.    wrong estimation can lead to memory leaks.        here   lattice_sz_y   takes into account  feedline symbol (or lattice_sz_incrrem_x already)    // !  !  !

node_main_symb='o'; ver_line_main_symb='|'; hor_line_main_symb='-';

  if (size <= 0) {return 0;}

  if (chstr_lattice == 0) {chstr_lattice=new char[size+1];} else {return 0;}

  for (int i=0;  i < size;  i++)
  {chstr_lattice[i]=' ';}
  
chstr_lattice[size]='\0';

return true;

}




void   st_char_str_lattice_console_draw :: draw_node(int node_num_y, int node_num_x, char symb)
{
  for (int y=0;  y < node_sz_y;  y++)
  {
    for (int x=0;  x < node_sz_x;  x++)
    {
    chstr_lattice[get_pos_node_ld_yx_with_shft(node_num_y, node_num_x, y, x)]=symb;
    }
  }
}

    
void  st_char_str_lattice_console_draw :: draw_nodes()
{

  if (chstr_lattice)
  {
    for (int y=0;  y < cells_amount_y+1;  y++)
    {
      for (int x=0;  x < cells_amount_x+1;  x++)
      {
      draw_node(y, x, node_main_symb);
      }
    }
  }

}


void  st_char_str_lattice_console_draw :: draw_ver_line(int num_y, int num_x, char symb)
{
  for (int y=0;  y < cell_sz_in_y;  y++)
  {
    for (int x=0;  x < node_sz_x;  x++)
    {
    chstr_lattice[get_pos_ver_line_l_yx_with_shft(num_y, num_x, y, x)]=symb;
    }
  } 
}


void  st_char_str_lattice_console_draw :: draw_hor_line(int num_y, int num_x, char symb)
{
  for (int y=0;  y < node_sz_y;  y++)
  {
    for (int x=0;  x < cell_sz_in_x;  x++)
    {
    chstr_lattice[get_pos_hor_line_d_yx_with_shft(num_y, num_x, y, x)]=symb;
    }
  } 
}


void  st_char_str_lattice_console_draw :: draw_ver_lines()
{

  if (chstr_lattice)
  {
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x+1;  x++)
      {
      draw_ver_line(y, x, ver_line_main_symb);
      } 
    }
  }

}
 

void  st_char_str_lattice_console_draw :: draw_hor_lines()
{

  if (chstr_lattice)
  {
    for (int y=0;  y < cells_amount_y+1;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      draw_hor_line(y, x, hor_line_main_symb);
      } 
    }
  }

}






void  st_char_str_lattice_console_draw :: draw_frame(char hor_frame_symb,  char ver_frame_symb)
{

  if (chstr_lattice)
  {
    for (int x=0;  x < lattice_sz_x;  x++)
    {chstr_lattice[x]= hor_frame_symb;  chstr_lattice[(lattice_sz_y-1)*(lattice_sz_incrrem_x)+x]=hor_frame_symb;}

    for (int y=1;  y < lattice_sz_y-1;  y++)
    {chstr_lattice[(lattice_sz_incrrem_x)*y]=ver_frame_symb;  chstr_lattice[(lattice_sz_incrrem_x)*y+lattice_sz_x-1]=ver_frame_symb;}
  }

}






void  st_char_str_lattice_console_draw :: draw_ver_lines_type1()
{

char symb='*';

  if (chstr_lattice)
  {
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x+1;  x++)
      {
      symb='*';
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
          chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb;          
          }
          //if (symb == '*') {symb=' ';} else {symb='*';}  
        }
                     
      } 
    }
  }

}



void  st_char_str_lattice_console_draw :: draw_hor_lines_type1()
{

char symb='*';

  if (chstr_lattice)
  {
    for (int y=0;  y < cells_amount_y+1;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {

        for (int y_=0;  y_ < node_sz_y;  y_++)
        {    
        symb=' ';
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;
            if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   
              
      } 
    }
  }

}


void  st_char_str_lattice_console_draw :: draw_num_in_ver_lines_type1_fbc(int *ver_num_array)
{

char chstr_num[64];  chstr_num[0]='\0';

  if ((chstr_lattice) && (ver_num_array))
  {
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x+1;  x++)
      {
      chstr_num[0]='\0';
      sprintf(chstr_num, "%d", ver_num_array[y*(cells_amount_x+1)+x]);
      chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2, 0)]=chstr_num[0];
        if ((ver_num_array[y*(cells_amount_x)+x] > 99) && (node_sz_x > 2)) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2, 2)]=chstr_num[2];}
        if ((ver_num_array[y*(cells_amount_x)+x] > 9 ) && (node_sz_x > 1)) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2, 1)]=chstr_num[1];}
      } 
    }
  }

chstr_num[0]='\0';

}


void  st_char_str_lattice_console_draw :: draw_num_in_hor_lines_type1_fbc(int *hor_num_array)
{

char chstr_num[64];  chstr_num[0]='\0';

  if ((chstr_lattice) && (hor_num_array))
  {
    for (int y=0;  y < cells_amount_y+1;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      chstr_num[0]='\0';
      sprintf(chstr_num, "%d", hor_num_array[y*cells_amount_x+x]);
      //chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2)]=chstr_num[0];
        if (hor_num_array[y*cells_amount_x+x] < 10) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2)]=chstr_num[0];}
        if ((hor_num_array[y*cells_amount_x+x] > 9) && (cell_sz_in_x > 1))  {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2-1)]=chstr_num[0];}      
        if ((hor_num_array[y*cells_amount_x+x] > 99) && (cell_sz_in_x > 2)) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2)+1]=chstr_num[2];}
        if ((hor_num_array[y*cells_amount_x+x] > 9 ) && (cell_sz_in_x > 1)) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2+0)]=chstr_num[1];}
      } 
    }
  }

chstr_num[0]='\0';

}



void  st_char_str_lattice_console_draw :: draw_num_in_ver_lines_type1_pbc(int *ver_num_array)
{

char chstr_num[64];  chstr_num[0]='\0';

  if ((chstr_lattice) && (ver_num_array))
  {
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      chstr_num[0]='\0';
      sprintf(chstr_num, "%d", ver_num_array[y*(cells_amount_x)+x]);
      chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2, 0)]=chstr_num[0];
        if ((ver_num_array[y*(cells_amount_x)+x] > 99) && (node_sz_x > 2)) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2, 2)]=chstr_num[2];}
        if ((ver_num_array[y*(cells_amount_x)+x] > 9 ) && (node_sz_x > 1)) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2, 1)]=chstr_num[1];}
      } 
    }
  }

chstr_num[0]='\0';

}


void  st_char_str_lattice_console_draw :: draw_num_in_hor_lines_type1_pbc(int *hor_num_array)
{

char chstr_num[64];  chstr_num[0]='\0';

  if ((chstr_lattice) && (hor_num_array))
  {
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      chstr_num[0]='\0';
      sprintf(chstr_num, "%d", hor_num_array[y*cells_amount_x+x]);
        if (hor_num_array[y*cells_amount_x+x] < 10) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2)]=chstr_num[0];}
        if ((hor_num_array[y*cells_amount_x+x] > 9) && (cell_sz_in_x > 1))  {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2-1)]=chstr_num[0];}      
        if ((hor_num_array[y*cells_amount_x+x] > 99) && (cell_sz_in_x > 2)) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2)+1]=chstr_num[2];}
        if ((hor_num_array[y*cells_amount_x+x] > 9 ) && (cell_sz_in_x > 1)) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, 0, cell_sz_in_x/2+0)]=chstr_num[1];}
      } 
    }
  }

chstr_num[0]='\0';

}




void st_char_str_lattice_console_draw :: draw_star_in_hor_lines_type1_fbc(int *transv_h_array)
{

char chstr_num[64];  chstr_num[0]='\0';

  for (int y=0;  y < cells_amount_y+1;  y++)
  {
    for (int x=0;  x < cells_amount_x;  x++)
    {
    //  for (int y_node=0;  y_node < node_sz_y;  y_node++)
    //  {
    //    for (int x_node=0;  x_node < node_sz_x;  x_node++)
    //    {
    //    chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, y_node, x_node)]=' ';
    //    }
    //  }  
    
    
    chstr_num[0]='\0';
    sprintf(chstr_num, "%d", transv_h_array[y*cells_amount_x+x]);    
      
      if ((transv_h_array[y*cells_amount_x+x] > 0) && ((transv_h_array[y*cells_amount_x+x] < 10)))
      {  
        if (node_sz_x > 0) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[0];}          
      }      

      if ((transv_h_array[y*cells_amount_x+x] > 9) && ((transv_h_array[y*cells_amount_x+x] < 100)))
      {  
        if (node_sz_x > 1) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 1)]=chstr_num[0]; chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[1];}  
        if (node_sz_x == 1) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[0];}          
      }      
            
      if ((transv_h_array[y*cells_amount_x+x] > 99))
      {  
        if (node_sz_x > 2) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 2)]=chstr_num[0]; chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 1)]=chstr_num[1];  chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[2];}  
        if (node_sz_x == 2) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 1)]=chstr_num[0]; chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[1];}  
        if (node_sz_x == 1) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[0];}          
      }
            
    }
  }    

chstr_num[0]='\0';  

}






void st_char_str_lattice_console_draw :: draw_star_in_hor_lines_type1_pbc(int *transv_h_array)
{

char chstr_num[64];  chstr_num[0]='\0';

  for (int y=0;  y < cells_amount_y;  y++)
  {
    for (int x=0;  x < cells_amount_x;  x++)
    {
    //  for (int y_node=0;  y_node < node_sz_y;  y_node++)
    //  {
    //    for (int x_node=0;  x_node < node_sz_x;  x_node++)
    //    {
    //    chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, y_node, x_node)]=' ';
    //    }
    //  }  
    
    
    chstr_num[0]='\0';
    sprintf(chstr_num, "%d", transv_h_array[y*cells_amount_x+x]);    
      
      if ((transv_h_array[y*cells_amount_x+x] > 0) && ((transv_h_array[y*cells_amount_x+x] < 10)))
      {  
        if (node_sz_x > 0) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[0];}          
      }      

      if ((transv_h_array[y*cells_amount_x+x] > 9) && ((transv_h_array[y*cells_amount_x+x] < 100)))
      {  
        if (node_sz_x > 1) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 1)]=chstr_num[0]; chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[1];}  
        if (node_sz_x == 1) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[0];}          
      }      
            
      if ((transv_h_array[y*cells_amount_x+x] > 99))
      {  
        if (node_sz_x > 2) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 2)]=chstr_num[0]; chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 1)]=chstr_num[1];  chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[2];}  
        if (node_sz_x == 2) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 1)]=chstr_num[0]; chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[1];}  
        if (node_sz_x == 1) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]=chstr_num[0];}          
      }
            
    }
  }    

chstr_num[0]='\0';  

}




void  st_char_str_lattice_console_draw :: draw_num_in_ver_lines_type2_fbc(int *ver_num_array)
{

char symb;

  if ((chstr_lattice) && (ver_num_array))
  {
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x+1;  x++)
      {

        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {    
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
            if (x_ >= ver_num_array[y*(cells_amount_x+1)+x]) {symb=' ';} else {symb='|';}
          chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb; 
          }
        }
        
      } 
    }
  }

}


void  st_char_str_lattice_console_draw :: draw_num_in_hor_lines_type2_fbc(int *hor_num_array)
{

char symb;

  if ((chstr_lattice) && (hor_num_array))
  {
    for (int y=0;  y < cells_amount_y+1;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      char symb='*';

        for (int y_=0;  y_ < node_sz_y;  y_++)
        {    
          if (y_*2 >= hor_num_array[y*cells_amount_x+x]) {symb=' ';} else {  if (y_*2+1 == hor_num_array[y*cells_amount_x+x]) {symb='-';} else {symb='=';}  }
           
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;
          }
        } 
          
      } 
    }
  }

}
  




void  st_char_str_lattice_console_draw :: draw_num_in_ver_lines_type2_pbc(int *ver_num_array)
{

char symb;

  if ((chstr_lattice) && (ver_num_array))
  {
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {

        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {    
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
            if (x_ >= ver_num_array[y*(cells_amount_x)+x]) {symb=' ';} else {symb='|';}
          chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb; 
          }
        }
        
      } 
    }
  }

}


void  st_char_str_lattice_console_draw :: draw_num_in_hor_lines_type2_pbc(int *hor_num_array)
{

char symb;

  if ((chstr_lattice) && (hor_num_array))
  {
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      char symb='*';

        for (int y_=0;  y_ < node_sz_y;  y_++)
        {    
          if (y_*2 >= hor_num_array[y*cells_amount_x+x]) {symb=' ';} else {  if (y_*2+1 == hor_num_array[y*cells_amount_x+x]) {symb='-';} else {symb='=';}  }
           
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;
          }
        } 
          
      } 
    }
  }

}




void st_char_str_lattice_console_draw :: draw_star_in_hor_lines_type2_fbc(int *transv_h_array)
{

char chstr_num[64];  chstr_num[0]='\0';

  for (int y=0;  y < cells_amount_y;  y++)
  {
    for (int x=0;  x < cells_amount_x+1;  x++)
    {
      for (int y_node=0;  y_node < node_sz_y;  y_node++)
      {
        for (int x_node=0;  x_node < node_sz_x;  x_node++)
        {
        chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, y_node, x_node)]=' ';
        }
      }  
        
      if (transv_h_array[y*cells_amount_x+x] > 0) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]='X';}      //  !  !  !
      
    }
  }

chstr_num[0]='\0';      

}





void st_char_str_lattice_console_draw :: draw_star_in_hor_lines_type2_pbc(int *transv_h_array)
{

char chstr_num[64];  chstr_num[0]='\0';

  for (int y=0;  y < cells_amount_y;  y++)
  {
    for (int x=0;  x < cells_amount_x;  x++)
    {
      for (int y_node=0;  y_node < node_sz_y;  y_node++)
      {
        for (int x_node=0;  x_node < node_sz_x;  x_node++)
        {
        chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, y_node, x_node)]=' ';
        }
      }  
        
      if (transv_h_array[y*cells_amount_x+x] > 0) {chstr_lattice[get_pos_node_ld_yx_with_shft(y, x, 0, 0)]='X';}      //  !  !  !
      
    }
  }

chstr_num[0]='\0';      

}




void st_char_str_lattice_console_draw :: lat_draw_with_nums_type_1_and_print_to_screen_fbc(int *ver_num_array, int *hor_num_array, int *transv_h_array) 
{ 

  if ((chstr_lattice) && (ver_num_array) && (hor_num_array))
  {
  draw_space(); draw_line_feeds();  draw_ver_lines_type1();  draw_hor_lines_type1(); draw_nodes();   draw_num_in_ver_lines_type1_fbc(ver_num_array); draw_num_in_hor_lines_type1_fbc(hor_num_array);
  //  if (transv_h_array != 0) {draw_star_in_hor_lines_type1_fbc(transv_h_array);}    
  print_to_screen();
  }

ver_num_array=0;  hor_num_array=0;      

}



void st_char_str_lattice_console_draw :: lat_draw_with_nums_type_2_and_print_to_screen_fbc(int *ver_num_array, int *hor_num_array, int *transv_h_array) 
{
    
  if ((chstr_lattice) && (ver_num_array) && (hor_num_array))
  {
  draw_space(); draw_line_feeds();  draw_ver_lines_type1();  draw_hor_lines_type1(); draw_nodes();   draw_num_in_ver_lines_type2_fbc(ver_num_array); draw_num_in_hor_lines_type2_fbc(hor_num_array);    
  //  if (transv_h_array != 0) {draw_star_in_hor_lines_type2_fbc(transv_h_array);}        
  print_to_screen();
  }    

ver_num_array=0;  hor_num_array=0;      

}



void st_char_str_lattice_console_draw :: lat_draw_with_nums_type_1_and_print_to_screen_pbc(int *ver_num_array, int *hor_num_array, int *transv_h_array) 
{ 

  if ((chstr_lattice) && (ver_num_array) && (hor_num_array))
  {
  draw_space(); draw_line_feeds();  draw_ver_lines_type1();  draw_hor_lines_type1(); draw_nodes();   draw_num_in_ver_lines_type1_pbc(ver_num_array); draw_num_in_hor_lines_type1_pbc(hor_num_array);  
  draw_pbc__delete_some_ver_lines();  draw_pbc__delete_some_hor_lines();
    if (transv_h_array != 0) {draw_star_in_hor_lines_type1_pbc(transv_h_array);}    
  print_to_screen();
  }

ver_num_array=0;  hor_num_array=0;      

}



void st_char_str_lattice_console_draw :: lat_draw_with_nums_type_2_and_print_to_screen_pbc(int *ver_num_array, int *hor_num_array, int *transv_h_array) 
{
    
  if ((chstr_lattice) && (ver_num_array) && (hor_num_array))
  {
  draw_space(); draw_line_feeds();  draw_ver_lines_type1();  draw_hor_lines_type1(); draw_nodes();   draw_num_in_ver_lines_type2_pbc(ver_num_array); draw_num_in_hor_lines_type2_pbc(hor_num_array);  
  draw_pbc__delete_some_ver_lines();  draw_pbc__delete_some_hor_lines();
    if (transv_h_array != 0) {draw_star_in_hor_lines_type2_pbc(transv_h_array);}        
  print_to_screen();
  }    

ver_num_array=0;  hor_num_array=0;      

}






void st_char_str_lattice_console_draw :: draw_fbc_J_bond_lattice(bool *J_ver, bool *J_hor, bool is_do_out_frame, char bond_bkg_symb, char mark_spin_by_symb)
{
//  good confugs for working    Ly,  Lx,  5,  9,  3,  3
//  good confugs for working    Ly,  Lx,  3,  5,  1,  3

  if ((chstr_lattice) && (J_ver) && (J_hor))
  {
  draw_space();
  
  
  // ver lines (=hor bonds)
  //
  char symb='*';

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=1;  x < cells_amount_x;  x++)
      {
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
          chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb;          
          }
 
        }
                     
        if (J_hor[y*(cell_sz_in_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
        {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}                               
      } 
    }


  // hor lines (=ver bonds)
  //    
    for (int y=1;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {

        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   
     
        if (J_ver[y*(cell_sz_in_x)+cell_sz_in_x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
        {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       
      }
         
    }     

    
  //  cells (=spins)
  //
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      //cout<<" y= "<<get_cell_dl_y0(y)<<" x= "<<get_cell_dl_x0(x)<<"  ";
      chstr_lattice[get_cell_middle_pos(y, x)]=mark_spin_by_symb;  
      }
    cout<<endl;  
    }  

      
      
  draw_line_feeds();   
    if (is_do_out_frame == true) {draw_frame();}   
  print_to_screen();      // draw_nodes();    
      
  }      //  if ((chstr_lattice) && (J_ver) && (J_hor))

}





void st_char_str_lattice_console_draw :: draw_fbc_spin_lattice(bool *spin_array, bool is_do_out_frame, char spin_up, char spin_down)
{
//  good confugs for working    Ly,  Lx,  5,  9,  3,  3
//  good confugs for working    Ly,  Lx,  3,  5,  1,  3

  if ((chstr_lattice) && (spin_array))
  {
  draw_space();
    
  //  cells (=spins)
  //
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      //cout<<" y= "<<get_cell_dl_y0(y)<<" x= "<<get_cell_dl_x0(x)<<"  ";
        if (spin_array[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down;}        
      }
    }  
     
  draw_line_feeds();   
    if (is_do_out_frame == true) {draw_frame();}   
  print_to_screen();      // draw_nodes();    
      
  }      //  if ((chstr_lattice) && (J_ver) && (J_hor))

}



void st_char_str_lattice_console_draw :: draw_fbc_J_bond_and_spin_lattice(bool *J_ver, bool *J_hor, bool *spin_array, bool is_do_out_frame, char bond_bkg_symb, char spin_up, char spin_down)
{

//  good confugs for working    Ly,  Lx,  5,  9,  3,  3
//  good confugs for working    Ly,  Lx,  3,  5,  1,  3

  if ((chstr_lattice) && (J_ver) && (J_hor) && (spin_array))
  {
  draw_space();
  
  
  // ver lines (=hor bonds)
  //
  char symb='*';

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=1;  x < cells_amount_x;  x++)
      {
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
          chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb;          
          }
 
        }
                     
//        if (J_hor[y*(cell_sz_in_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
//        {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}
        if (J_hor[y*(cells_amount_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
        {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}
                                       
      } 
    }




  // hor lines (=ver bonds)
  //    
    for (int y=1;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {

        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   
     
//        if (J_ver[y*(cell_sz_in_x)+cell_sz_in_x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
//        {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       
        if (J_ver[y*(cells_amount_x)+x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
        {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       

      }
         
    }     

    
  //  cells (=spins)
  //
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      //cout<<" y= "<<get_cell_dl_y0(y)<<" x= "<<get_cell_dl_x0(x)<<"  ";
        if (spin_array[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down;}        
      }
    }  

      
      
  draw_line_feeds();   
    if (is_do_out_frame == true) {draw_frame();}   
  print_to_screen();      // draw_nodes();    
      
  }      //  if ((chstr_lattice) && (J_ver) && (J_hor))

}








void  st_char_str_lattice_console_draw :: draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond(bool *J_ver, bool *J_hor, bool *spin_array, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, bool is_do_out_frame, char bond_bkg_symb, char spin_up, char spin_down)
{

//  good confugs for working    Ly,  Lx,  5,  9,  3,  3
//  good confugs for working    Ly,  Lx,  3,  5,  1,  3

char bond_separ_symb='^';

  if ((chstr_lattice) && (J_ver) && (J_hor))
  {
  draw_space();
  
  
  // ver lines (=hor bonds)
  //
  char symb='*';

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=1;  x < cells_amount_x;  x++)
      {
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
          chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb;          
          }
 
        }
                     
        if (J_hor[y*(cell_sz_in_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
        {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}                               
      } 
    }


  // hor lines (=ver bonds)
  //    
    for (int y=1;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {

        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   
     
        if (J_ver[y*(cell_sz_in_x)+cell_sz_in_x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
        {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       
      }
         
    }     

    
  //  cells (=spins)
  //
    

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      //cout<<" y= "<<get_cell_dl_y0(y)<<" x= "<<get_cell_dl_x0(x)<<"  ";
        if (spin_array)
        {  if (spin_array[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down;}  }
        else 
        {chstr_lattice[get_cell_middle_pos(y, x)]='S';}        
      }
    }  

     
    
  // marking of flip bond  
  //    
    if (is_1ver_0hor == true)
    {
    //  ver lines (=hor bonds)
    //    
    //  J_bond_y_index, int J_bond_x_index
    //
    int y=J_bond_y_index;
    int x=J_bond_x_index;
    
    //for (int y=0;  y < cells_amount_y;  y++)
      {
        //for (int x=1;  x < cells_amount_x;  x++)
        {
        symb=bond_bkg_symb;
      
          for (int y_=0;  y_ < cell_sz_in_y;  y_++)
          {      
            for (int x_=0;  x_ < node_sz_x;  x_++)
            {
            chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=bond_separ_symb;          
            }
          }
                     
          if (J_hor[y*(cell_sz_in_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}
                                         
        } 
      }
    }      
     

    if (is_1ver_0hor == false)
    {
    // hor lines (=ver bonds)
    //
    //  J_bond_y_index, int J_bond_x_index
    //           
    int y=J_bond_y_index;
    int x=J_bond_x_index;
    
      //for (int y=1;  y < cells_amount_y;  y++)
      {
        //for (int x=0;  x < cells_amount_x;  x++)
        {

          for (int y_=0;  y_ < node_sz_y;  y_++)  
          {    
          symb=bond_bkg_symb;
          
            for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
            {
            chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=bond_separ_symb;  
            //  if (symb == '*') {symb=' ';} else {symb='*';}
            }
          }   
     
          if (J_ver[y*(cell_sz_in_x)+cell_sz_in_x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
          {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}     
                                            
        }
         
      }
      
    }     
      
      
      
  draw_line_feeds();   
    if (is_do_out_frame == true) {draw_frame();}   
  print_to_screen();      // draw_nodes();    
      
  }      //  if ((chstr_lattice) && (J_ver) && (J_hor))


}







  //   J_bond_y_index  and   int J_bond_x_index  can be negative
void st_char_str_lattice_console_draw :: draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame(bool *J_ver, bool *J_hor, bool *spin_array, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x, bool is_do_out_frame, char bond_bkg_symb, char spin_up, char spin_down, char inner_frame_symb)
{

//is_1ver_0hor=true; J_bond_y_index=8; J_bond_x_index=0;      //--DEBUG--//  
//d_len=6;  u_len=0;  l_len=0;  r_len=6;  start_y=3; start_x=0;      //--DEBUG--//

//cout<<"is_1ver_0hor="<<((int) (is_1ver_0hor))<<"  J_bond_y_index="<<J_bond_y_index<<"  J_bond_x_index="<<J_bond_x_index<<endl;      //--DEBUG--//
//cout<<"d_len="<<d_len<<"  u_len="<<u_len<<"  l_len="<<l_len<<"  r_len="<<r_len<<"  start_y="<<start_y<<" start_x="<<start_x<<endl;      //--DEBUG--//


//  good confugs for working    Ly,  Lx,  5,  9,  3,  3
//  good confugs for working    Ly,  Lx,  3,  5,  1,  3

char bond_separ_symb='^';

  if ((chstr_lattice) && (J_ver) && (J_hor))
  {
  draw_space();
  
    
  //  cells (=spins)
  //
    

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      //cout<<" y= "<<get_cell_dl_y0(y)<<" x= "<<get_cell_dl_x0(x)<<"  ";
        if (spin_array)
        {  if (spin_array[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down;}  }
        else 
        {chstr_lattice[get_cell_middle_pos(y, x)]='S';}        
      }
    }  

            
 
 
    
    
     
    if (is_do_out_frame == true) {draw_frame();}  
      
 
 
 
      
  //  drawing inner frame
  //
  //  int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x,
  //
  
  
  
  //  bottom
  //
    if (d_len > 0)
    {
    // down part
    draw_node(start_y-1, start_x, inner_frame_symb);     
    
      for (int cell_num_x=start_x;  cell_num_x < start_x+d_len;  cell_num_x++)
      {
      draw_node(start_y-1, cell_num_x+1, inner_frame_symb);
      draw_hor_line(start_y-1, cell_num_x, inner_frame_symb);
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)
        


    // up part 
    draw_node(start_y+1-1, start_x, inner_frame_symb);     
    
      for (int cell_num_x=start_x;  cell_num_x < start_x+d_len;  cell_num_x++)
      {
      draw_node(start_y+1-1, cell_num_x+1, inner_frame_symb);
        if (!(((l_len > 0) && (cell_num_x == start_x)) || ((r_len > 0) && (cell_num_x == start_x+d_len-1)))) {draw_hor_line(start_y+1-1, cell_num_x, inner_frame_symb);}
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     
       
       
       if (l_len > 0) {draw_ver_line(start_y-1, start_x,                inner_frame_symb);} 
       if (r_len > 0) {draw_ver_line(start_y-1, start_x+d_len,          inner_frame_symb);}
        
    }   



  //  top
  //
    if (u_len > 0)
    {
    int y_row_not_null=start_y+l_len;
      if (r_len > 0) {y_row_not_null=start_y+r_len;}    
    int row_len_not_null=d_len;
      if (u_len > 0) {row_len_not_null=u_len;}    
      
    
    // up part
    draw_node(y_row_not_null+1, start_x, inner_frame_symb);     
    
      for (int cell_num_x=start_x;  cell_num_x < start_x+u_len;  cell_num_x++)
      {
      draw_node(y_row_not_null+1, cell_num_x+1, inner_frame_symb);
      draw_hor_line(y_row_not_null+1, cell_num_x, inner_frame_symb);
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     


    //  down part 
    draw_node(y_row_not_null, start_x, inner_frame_symb);     
    
      for (int cell_num_x=start_x;  cell_num_x < start_x+u_len;  cell_num_x++)
      {
      draw_node(y_row_not_null, cell_num_x+1, inner_frame_symb);
        if (!(((l_len > 0) && (cell_num_x == start_x)) || ((r_len > 0) && (cell_num_x == start_x+u_len-1)))) {draw_hor_line(y_row_not_null, cell_num_x, inner_frame_symb);}
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     
       
       
       if (l_len > 0) {draw_ver_line(y_row_not_null, start_x,                        inner_frame_symb);} 
       if (r_len > 0) {draw_ver_line(y_row_not_null, start_x+row_len_not_null,       inner_frame_symb);}
        
    }

 
  
  //  left
  //
    if (l_len > 0)
    {
    // left part
    draw_node(start_y, start_x, inner_frame_symb);     
    
      for (int cell_num_y=start_y;  cell_num_y < start_y+l_len;  cell_num_y++)
      {
      draw_node(cell_num_y, start_x, inner_frame_symb);
      draw_ver_line(cell_num_y, start_x, inner_frame_symb);
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     


    // right part 
    draw_node(start_y, start_x+1, inner_frame_symb);     
    
      for (int cell_num_y=start_y;  cell_num_y < start_y+l_len;  cell_num_y++)
      {
      draw_node(cell_num_y, start_x+1, inner_frame_symb);
      draw_ver_line(cell_num_y, start_x+1, inner_frame_symb);
      } 
             
       
      //if (d_len > 0) {draw_hor_line(start_y,         start_x,          inner_frame_symb);} 
      //if (u_len > 0) {draw_hor_line(start_y+l_len, start_x,          inner_frame_symb);}
        
    }

 
  
  //  right
  //
    if (r_len > 0)
    { 
    int row_len_not_null=d_len;
      if (u_len > 0) {row_len_not_null=u_len;}        
    int x_col=start_x+d_len-1;
      if (u_len > 0) {x_col=start_x+u_len-1;}

    
    // right part
    draw_node(start_y, x_col+1, inner_frame_symb);     
    
      for (int cell_num_y=start_y;  cell_num_y < start_y+r_len;  cell_num_y++)
      {
      draw_node(cell_num_y, x_col+1, inner_frame_symb);
      draw_ver_line(cell_num_y, x_col+1, inner_frame_symb);
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     


    // left part 
    draw_node(start_y, x_col, inner_frame_symb);     
    
      for (int cell_num_y=start_y;  cell_num_y < start_y+r_len;  cell_num_y++)
      {
      draw_node(cell_num_y, x_col, inner_frame_symb);
      draw_ver_line(cell_num_y, x_col, inner_frame_symb);
      }          
       
      //if (d_len > 0) {draw_hor_line(start_y,         start_x+r_len+1,          inner_frame_symb);} 
      //if (u_len > 0) {draw_hor_line(start_y+l_len, start_x+r_len+1,          inner_frame_symb);}
        
    }





  
  // ver lines (=hor bonds)
  //
  char symb='*';

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=1;  x < cells_amount_x;  x++)
      {
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
          //--//chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb;          
          }
 
        }
            
        if (J_hor)
        {             
          if (J_hor[y*(cells_amount_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}
        }
                                         
      } 
    }


  // hor lines (=ver bonds)
  //    
    for (int y=1;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {

        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          //--//chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   

        if (J_ver)
        {             
          if (J_ver[y*(cells_amount_x)+x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
          {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       
        }

      }
         
    }     




     
  // marking of flip bond  
  //  
    if (is_1ver_0hor == false)
    {
    //  ver lines (=hor bonds)
    //    
    //  J_bond_y_index, int J_bond_x_index
    //
    int y=J_bond_y_index;
    int x=J_bond_x_index;
    
    //for (int y=0;  y < cells_amount_y;  y++)
      {
        //for (int x=1;  x < cells_amount_x;  x++)
        {
        symb=bond_bkg_symb;
      
          for (int y_=0;  y_ < cell_sz_in_y;  y_++)
          {      
            for (int x_=0;  x_ < node_sz_x;  x_++)
            {
            chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=bond_separ_symb;          
            }
          }
                     
          if (J_hor[y*(cells_amount_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}
                                         
        } 
      }
    }      
     
  // marking of flip bond  
  //
    if (is_1ver_0hor == true)
    {
    // hor lines (=ver bonds)
    //
    //  J_bond_y_index, int J_bond_x_index
    //           
    int y=J_bond_y_index;
    int x=J_bond_x_index;
    
      //for (int y=1;  y < cells_amount_y;  y++)
      {
        //for (int x=0;  x < cells_amount_x;  x++)
        {

          for (int y_=0;  y_ < node_sz_y;  y_++)  
          {    
          symb=bond_bkg_symb;
          
            for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
            {
            chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=bond_separ_symb;  
            //  if (symb == '*') {symb=' ';} else {symb='*';}
            }
          }   
     
          if (J_ver[y*(cells_amount_x)+x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
          {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}     
                                            
        }
         
      }
      
    }     




 
       
      
  draw_line_feeds();   
  print_to_screen();      // draw_nodes();    
      
  }      //  if ((chstr_lattice) && (J_ver) && (J_hor))


}






  //   J_bond_y_index  and   int J_bond_x_index  can be negative
void st_char_str_lattice_console_draw :: draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame_cmp(bool *J_ver, bool *J_hor, bool *spin_array_to_cmp, bool *spin_array_cur, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x, bool is_do_out_frame, char bond_bkg_symb, char spin_up_same, char spin_down_same, char spin_up_diff, char spin_down_diff, char inner_frame_symb)
{

//is_1ver_0hor=true; J_bond_y_index=8; J_bond_x_index=0;      //--DEBUG--//  
//d_len=6;  u_len=0;  l_len=0;  r_len=6;  start_y=3; start_x=0;      //--DEBUG--//

//cout<<"is_1ver_0hor="<<((int) (is_1ver_0hor))<<"  J_bond_y_index="<<J_bond_y_index<<"  J_bond_x_index="<<J_bond_x_index<<endl;      //--DEBUG--//
//cout<<"d_len="<<d_len<<"  u_len="<<u_len<<"  l_len="<<l_len<<"  r_len="<<r_len<<"  start_y="<<start_y<<" start_x="<<start_x<<endl;      //--DEBUG--//


//  good confugs for working    Ly,  Lx,  5,  9,  3,  3
//  good confugs for working    Ly,  Lx,  3,  5,  1,  3

char bond_separ_symb='^';

  if ((chstr_lattice) && (J_ver) && (J_hor) && (spin_array_to_cmp) && (spin_array_cur))
  {
  draw_space();
  
    
  //  cells (=spins)
  //
    

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      //cout<<" y= "<<get_cell_dl_y0(y)<<" x= "<<get_cell_dl_x0(x)<<"  ";
      
        if (spin_array_to_cmp[y*cells_amount_x+x] == spin_array_cur[y*cells_amount_x+x])
        {  
          if (spin_array_cur[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up_same;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down_same;} 
        }
        else
        {
          if (spin_array_cur[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up_diff;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down_diff;} 
        }
        
      }
    }  

            
 
 
    
    
     
    if (is_do_out_frame == true) {draw_frame();}  
      
 
 
 
      
  //  drawing inner frame
  //
  //  int d_len, int u_len, int l_len, int r_len,   int start_y,  int start_x,
  //
  
  
  
  //  bottom
  //
    if (d_len > 0)
    {
    // down part
    draw_node(start_y-1, start_x, inner_frame_symb);     
    
      for (int cell_num_x=start_x;  cell_num_x < start_x+d_len;  cell_num_x++)
      {
      draw_node(start_y-1, cell_num_x+1, inner_frame_symb);
      draw_hor_line(start_y-1, cell_num_x, inner_frame_symb);
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)
        


    // up part 
    draw_node(start_y+1-1, start_x, inner_frame_symb);     
    
      for (int cell_num_x=start_x;  cell_num_x < start_x+d_len;  cell_num_x++)
      {
      draw_node(start_y+1-1, cell_num_x+1, inner_frame_symb);
        if (!(((l_len > 0) && (cell_num_x == start_x)) || ((r_len > 0) && (cell_num_x == start_x+d_len-1)))) {draw_hor_line(start_y+1-1, cell_num_x, inner_frame_symb);}
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     
       
       
       if (l_len > 0) {draw_ver_line(start_y-1, start_x,                inner_frame_symb);} 
       if (r_len > 0) {draw_ver_line(start_y-1, start_x+d_len,          inner_frame_symb);}
        
    }   



  //  top
  //
    if (u_len > 0)
    {
    int y_row_not_null=start_y+l_len;
      if (r_len > 0) {y_row_not_null=start_y+r_len;}    
    int row_len_not_null=d_len;
      if (u_len > 0) {row_len_not_null=u_len;}    
      
    
    // up part
    draw_node(y_row_not_null+1, start_x, inner_frame_symb);     
    
      for (int cell_num_x=start_x;  cell_num_x < start_x+u_len;  cell_num_x++)
      {
      draw_node(y_row_not_null+1, cell_num_x+1, inner_frame_symb);
      draw_hor_line(y_row_not_null+1, cell_num_x, inner_frame_symb);
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     


    //  down part 
    draw_node(y_row_not_null, start_x, inner_frame_symb);     
    
      for (int cell_num_x=start_x;  cell_num_x < start_x+u_len;  cell_num_x++)
      {
      draw_node(y_row_not_null, cell_num_x+1, inner_frame_symb);
        if (!(((l_len > 0) && (cell_num_x == start_x)) || ((r_len > 0) && (cell_num_x == start_x+u_len-1)))) {draw_hor_line(y_row_not_null, cell_num_x, inner_frame_symb);}
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     
       
       
       if (l_len > 0) {draw_ver_line(y_row_not_null, start_x,                        inner_frame_symb);} 
       if (r_len > 0) {draw_ver_line(y_row_not_null, start_x+row_len_not_null,       inner_frame_symb);}
        
    }

 
  
  //  left
  //
    if (l_len > 0)
    {
    // left part
    draw_node(start_y, start_x, inner_frame_symb);     
    
      for (int cell_num_y=start_y;  cell_num_y < start_y+l_len;  cell_num_y++)
      {
      draw_node(cell_num_y, start_x, inner_frame_symb);
      draw_ver_line(cell_num_y, start_x, inner_frame_symb);
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     


    // right part 
    draw_node(start_y, start_x+1, inner_frame_symb);     
    
      for (int cell_num_y=start_y;  cell_num_y < start_y+l_len;  cell_num_y++)
      {
      draw_node(cell_num_y, start_x+1, inner_frame_symb);
      draw_ver_line(cell_num_y, start_x+1, inner_frame_symb);
      } 
             
       
      //if (d_len > 0) {draw_hor_line(start_y,         start_x,          inner_frame_symb);} 
      //if (u_len > 0) {draw_hor_line(start_y+l_len, start_x,          inner_frame_symb);}
        
    }

 
  
  //  right
  //
    if (r_len > 0)
    { 
    int row_len_not_null=d_len;
      if (u_len > 0) {row_len_not_null=u_len;}        
    int x_col=start_x+d_len-1;
      if (u_len > 0) {x_col=start_x+u_len-1;}

    
    // right part
    draw_node(start_y, x_col+1, inner_frame_symb);     
    
      for (int cell_num_y=start_y;  cell_num_y < start_y+r_len;  cell_num_y++)
      {
      draw_node(cell_num_y, x_col+1, inner_frame_symb);
      draw_ver_line(cell_num_y, x_col+1, inner_frame_symb);
      }      //  for (int cell_num_x=0;  cell_num_x < cells_amount_x;  cell_num_x++)     


    // left part 
    draw_node(start_y, x_col, inner_frame_symb);     
    
      for (int cell_num_y=start_y;  cell_num_y < start_y+r_len;  cell_num_y++)
      {
      draw_node(cell_num_y, x_col, inner_frame_symb);
      draw_ver_line(cell_num_y, x_col, inner_frame_symb);
      }          
       
      //if (d_len > 0) {draw_hor_line(start_y,         start_x+r_len+1,          inner_frame_symb);} 
      //if (u_len > 0) {draw_hor_line(start_y+l_len, start_x+r_len+1,          inner_frame_symb);}
        
    }




  
  // ver lines (=hor bonds)
  //
  char symb='*';

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=1;  x < cells_amount_x;  x++)
      {
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
          //--//chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb;          
          }
 
        }
            
        if (J_hor)
        {             
          if (J_hor[y*(cells_amount_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}
        }
                                         
      } 
    }


  // hor lines (=ver bonds)
  //    
    for (int y=1;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {

        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          //--//chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   

        if (J_ver)
        {             
          if (J_ver[y*(cells_amount_x)+x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
          {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       
        }

      }
         
    }     




     
  // marking of flip bond  
  //  
    if (is_1ver_0hor == false)
    {
    //  ver lines (=hor bonds)
    //    
    //  J_bond_y_index, int J_bond_x_index
    //
    int y=J_bond_y_index;  
    int x=J_bond_x_index; 
    
    //for (int y=0;  y < cells_amount_y;  y++)
      {
        //for (int x=1;  x < cells_amount_x;  x++)
        {
        symb=bond_bkg_symb;
      
          for (int y_=0;  y_ < cell_sz_in_y;  y_++)
          {      
            for (int x_=0;  x_ < node_sz_x;  x_++)
            {
            chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=bond_separ_symb;          
            }
          }
                     
          if (J_hor[y*(cells_amount_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}
                                         
        } 
      }
    }      
     
  // marking of flip bond  
  //
    if (is_1ver_0hor == true)
    {
    // hor lines (=ver bonds)
    //
    //  J_bond_y_index, int J_bond_x_index
    //           
    int y=J_bond_y_index;
    int x=J_bond_x_index;
    
      //for (int y=1;  y < cells_amount_y;  y++)
      {
        //for (int x=0;  x < cells_amount_x;  x++)
        {

          for (int y_=0;  y_ < node_sz_y;  y_++)  
          {    
          symb=bond_bkg_symb;
          
            for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
            {
            chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=bond_separ_symb;  
            //  if (symb == '*') {symb=' ';} else {symb='*';}
            }
          }   
     
          if (J_ver[y*(cells_amount_x)+x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
          {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}     
                                            
        }
         
      }
      
    }     




 
       
      
  draw_line_feeds();   
  print_to_screen();      // draw_nodes();    
      
  }      //  if ((chstr_lattice) && (J_ver) && (J_hor))


}






void st_char_str_lattice_console_draw :: draw_fbc_J_bond_and_spin_lattice_frame_proj_v3(bool *J_ver, bool *J_hor, bool *spin_array_, bool * array_spin_exist_d_u_l_r, bool is_do_out_frame, char bond_bkg_symb, char spin_up, char spin_down, bool draw_on_screen, bool draw_to_file, char *filename_dest, bool newfile)
{

//  good confugs for working    Ly,  Lx,  5,  9,  3,  3
//  good confugs for working    Ly,  Lx,  3,  5,  1,  3

char bond_separ_symb='^';

  if ((chstr_lattice) && (J_ver) && (J_hor))
  {
  draw_space();

      
  //  creating spin_array
  //  
  bool *spin_array=0;
  spin_array=new bool[cells_amount_y*cells_amount_x];
  
    for (int i=0;  i < cells_amount_y*cells_amount_x;  i++)
    {spin_array[i]=false;}

    if (spin_array_)
    {
      for (int y=0;  y < cells_amount_y-2;  y++)
      {
        for (int x=0;  x < cells_amount_x-2;  x++)
        {spin_array[cells_amount_x+y*cells_amount_x+ 1+x]=spin_array_[y*(cells_amount_x-2)+x];}
      }
    } 
    
    for (int x=0;  x < cells_amount_x-2;  x++) 
    {
      if (array_spin_exist_d_u_l_r[2*x                     +1] == true) {spin_array[                                  1+x]=array_spin_exist_d_u_l_r[                     2*x];  /*cout<<"FFFF:  "<<array_spin_exist_d_u_l_r[                     2*x]<<" ";*/}  
      if (array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*x+1] == true) {spin_array[(cells_amount_y-1)*cells_amount_x+1+x]=array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*x];  /*cout<<"FFFF:  "<<array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*x]<<" ";*/}
    } 
    for (int y=0;  y < cells_amount_y-2;  y++) 
    {
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)                     +2*y+1] == true) {spin_array[(1+y)*cells_amount_x  ]=array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+                     2*y];}
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(cells_amount_y-2)+2*y+1] == true) {spin_array[(2+y)*cells_amount_x-1]=array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(cells_amount_y-2)+2*y];}
    }  


  
  // ver lines (=hor bonds)
  //
  char symb='*';

    for (int y=1;  y < cells_amount_y-1;  y++)
    {
      for (int x=1;  x < cells_amount_x;  x++)
      {
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
          chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb;          
          }
 
        }
                     
      //int node_shift_y=0, node_shift_y=0  
                     
        if (J_hor[(y-1)*(cells_amount_x-1)+ x-1] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
        {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}                               
      } 
    }   


  // hor lines (=ver bonds)
  //    
    for (int y=1;  y < cells_amount_y;  y++)
    {
      for (int x=1;  x < cells_amount_x-1;  x++)
      {

        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   
     
        if (J_ver[(y-1)*(cells_amount_x-2)+x-1] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
        {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       
      }
         
    }   

    
    
    
  //  cells (=spins)
  //

    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      //cout<<" y= "<<get_cell_dl_y0(y)<<" x= "<<get_cell_dl_x0(x)<<"  ";
        if ((spin_array_) || (y == 0) || (y == cells_amount_y-1) || (x == 0) || (x == cells_amount_x-1))
        {  if (spin_array[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down;}  }
        else 
        {chstr_lattice[get_cell_middle_pos(y, x)]='S';}        
      }
    }  






  //  edge empties (down, up)
  //
    for (int x=1;  x < cells_amount_x-1;  x++)
    {
    
      if (array_spin_exist_d_u_l_r[2*(x-1)+1] == false )
      {
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {chstr_lattice[get_cell_dl_pos_with_shift(0, x, y_, x_)]=' ';}
        }      
        
        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(1, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   
      }
      
      if (array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*(x-1)+1] == false )
      {
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {chstr_lattice[get_cell_dl_pos_with_shift(cells_amount_y-1, x, y_, x_)]=' ';}
        }
        
              
        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(1+cells_amount_y-2, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }           
      }      
      
    }


  //  edge empties (left, right)
  //
    for (int y=1;  y < cells_amount_y-1;  y++)
    {
    
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(y-1)+1] == false)
      {
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {chstr_lattice[get_cell_dl_pos_with_shift(y, 0, y_, x_)]=' ';}
        }
        
       symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, 1, y_, x_)]=symb;}
        }       
      }
      
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(cells_amount_y-2)+2*(y-1)+1] == false )
      {
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {chstr_lattice[get_cell_dl_pos_with_shift(y, cells_amount_x-1, y_, x_)]=' ';}
        }
        
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, 1+cells_amount_x-2, y_, x_)]=symb;}
        }       
                    
      }      
      
    }





     
      
  //  draw empties in corners   
  //
  
    for (int y_=0;  y_ < cell_sz_in_y;  y_++)
    {      
      for (int x_=0;  x_ < cell_sz_in_x;  x_++)
      { 
      chstr_lattice[get_cell_dl_pos_with_shift(0,                0, y_, x_)]=' ';  chstr_lattice[get_cell_dl_pos_with_shift(0,                cells_amount_x-1, y_, x_)]=' '; 
      chstr_lattice[get_cell_dl_pos_with_shift(cells_amount_y-1, 0, y_, x_)]=' ';  chstr_lattice[get_cell_dl_pos_with_shift(cells_amount_y-1, cells_amount_x-1, y_, x_)]=' ';                
      }
    }
// get_cell_dl_pos_with_shift(int cell_num_y, int cell_num_x, int shft_y, int shft_x)


    if (is_do_out_frame == true) {draw_frame();} 
  draw_line_feeds();
      
      
      
  //  result out 
  //  args:        (... ,  bool draw_on_screen,  bool draw_to_file,  char *filename_dest,  bool newfile)    
  //    
    if (draw_on_screen == true)
    {  
    print_to_screen();      // draw_nodes();      //    char *chstr_lattice;
    }    
    
    if ((draw_to_file == true) && (filename_dest) && (strlen(filename_dest) > 0))
    {
    write_chstr_to_file(chstr_lattice, filename_dest, newfile);
    }
    
  
  
    if (spin_array) {delete[] spin_array;  spin_array=0;}  
      
  }      //  if ((chstr_lattice) && (J_ver) && (J_hor))


}










  //   J_bond_y_index  and   int J_bond_x_index  can be negative
void st_char_str_lattice_console_draw :: draw_fbc_J_bond_and_spin_lattice_with_marking_flip_bond_and_inner_frame_cmp_frame_proj_v3(bool *J_ver, bool *J_hor, bool *spin_array_to_cmp_, bool *spin_array_cur_, bool * array_spin_exist_d_u_l_r, bool is_1ver_0hor, int J_bond_y_index, int J_bond_x_index, int start_y,  int start_x, bool is_do_out_frame, char bond_bkg_symb, char spin_up_same, char spin_down_same, char spin_up_diff, char spin_down_diff, char inner_frame_symb,  char domain_wall_symb, bool is_draw_domain_wall)
{

//is_1ver_0hor=true; J_bond_y_index=8; J_bond_x_index=0;      //--DEBUG--//  
//d_len=6;  u_len=0;  l_len=0;  r_len=6;  start_y=3; start_x=0;      //--DEBUG--//

//cout<<"is_1ver_0hor="<<((int) (is_1ver_0hor))<<"  J_bond_y_index="<<J_bond_y_index<<"  J_bond_x_index="<<J_bond_x_index<<endl;      //--DEBUG--//
//cout<<"d_len="<<d_len<<"  u_len="<<u_len<<"  l_len="<<l_len<<"  r_len="<<r_len<<"  start_y="<<start_y<<" start_x="<<start_x<<endl;      //--DEBUG--//




//  good confugs for working    Ly,  Lx,  5,  9,  3,  3
//  good confugs for working    Ly,  Lx,  3,  5,  1,  3

char bond_separ_symb='^';
char symb=' ';

  if ((chstr_lattice) && (J_ver) && (J_hor) && (spin_array_to_cmp_) && (spin_array_cur_))
  {
  draw_space();
  


  //  creating  spin_array_cur  and  spin_array_to_cmp
  //  
  bool *spin_array_cur=0;                                          bool *spin_array_to_cmp=0;
  spin_array_cur=new bool[cells_amount_y*cells_amount_x];          spin_array_to_cmp=new bool[cells_amount_y*cells_amount_x];
  
    for (int i=0;  i < cells_amount_y*cells_amount_x;  i++)
    {spin_array_cur[i]=false;  spin_array_to_cmp[i]=false;}

    for (int y=0;  y < cells_amount_y-2;  y++)
    {
      for (int x=0;  x < cells_amount_x-2;  x++)
      {spin_array_cur[cells_amount_x+y*cells_amount_x+ 1+x]=spin_array_cur_[y*(cells_amount_x-2)+x];  spin_array_to_cmp[cells_amount_x+y*cells_amount_x+ 1+x]=spin_array_to_cmp_[y*(cells_amount_x-2)+x];}
    }
    
    
    for (int x=0;  x < cells_amount_x-2;  x++) 
    {
      if (array_spin_exist_d_u_l_r[2*x                     +1] == true) {spin_array_cur[                                  1+x]=array_spin_exist_d_u_l_r[                     2*x];}  
      if (array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*x+1] == true) {spin_array_cur[(cells_amount_y-1)*cells_amount_x+1+x]=array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*x];}
      
      if (array_spin_exist_d_u_l_r[2*x                     +1] == true) {spin_array_to_cmp[                                  1+x]=array_spin_exist_d_u_l_r[                     2*x];}  
      if (array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*x+1] == true) {spin_array_to_cmp[(cells_amount_y-1)*cells_amount_x+1+x]=array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*x];}      
    } 
    for (int y=0;  y < cells_amount_y-2;  y++) 
    {
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)                     +2*y+1] == true) {spin_array_cur[(1+y)*cells_amount_x  ]=array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+                     2*y];}
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(cells_amount_y-2)+2*y+1] == true) {spin_array_cur[(2+y)*cells_amount_x-1]=array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(cells_amount_y-2)+2*y];}
      
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)                     +2*y+1] == true) {spin_array_to_cmp[(1+y)*cells_amount_x  ]=array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+                     2*y];}
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(cells_amount_y-2)+2*y+1] == true) {spin_array_to_cmp[(2+y)*cells_amount_x-1]=array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(cells_amount_y-2)+2*y];}
    } 
    






  //  domains
  //
    if (is_draw_domain_wall == true)
    {
    
      for (int y=0;  y < cells_amount_y;  y++)
      {
        for (int x=0;  x < cells_amount_x;  x++)
        {

          if (spin_array_to_cmp[y*cells_amount_x+x] != spin_array_cur[y*cells_amount_x+x])
          {
            for (int y_=0;  y_ < cell_sz_in_y;  y_++)
            {    
            char symb=domain_wall_symb;
              for (int x_=0;  x_ < cell_sz_in_x;  x_++)
              {
              chstr_lattice[get_cell_dl_pos_with_shift(y, x, y_, x_)]=symb;     
              }
            }   
          }
              
        }
      }
      
    }
  


  // ver lines (=hor bonds)
  //
  char symb='*';

    for (int y=1;  y < cells_amount_y-1;  y++)
    {
      for (int x=1;  x < cells_amount_x;  x++)
      {
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {
          chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=symb;          
          }
 
        }
                     
      //int node_shift_y=0, node_shift_y=0  
                     
        if (J_hor[(y-1)*(cells_amount_x-1)+ x-1] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
        {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}                               
      } 
    }   


  // hor lines (=ver bonds)
  //    
    for (int y=1;  y < cells_amount_y;  y++)
    {
      for (int x=1;  x < cells_amount_x-1;  x++)
      {

        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   
     
        if (J_ver[(y-1)*(cells_amount_x-2)+x-1] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
        {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       
      }
         
    }   



  //  cells (=spins)
  //    
    for (int y=0;  y < cells_amount_y;  y++)
    {
      for (int x=0;  x < cells_amount_x;  x++)
      {
      //cout<<" y= "<<get_cell_dl_y0(y)<<" x= "<<get_cell_dl_x0(x)<<"  ";
        if ((spin_array_cur) || (y == 0) || (y == cells_amount_y-1) || (x == 0) || (x == cells_amount_x-1))
        {  
        
          if (spin_array_to_cmp[y*cells_amount_x+x] == spin_array_cur[y*cells_amount_x+x])
          {  
            if (spin_array_cur[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up_same;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down_same;} 
          }
          else
          {
            if (spin_array_cur[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up_diff;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down_diff;} 
          }
                
          //  if (spin_array[y*cells_amount_x+x] == true) {chstr_lattice[get_cell_middle_pos(y, x)]=spin_up;} else {chstr_lattice[get_cell_middle_pos(y, x)]=spin_down;}  
        }
        else 
        {chstr_lattice[get_cell_middle_pos(y, x)]='S';}        
      }
    }  







  //  edge empties (down, up)
  //
    for (int x=1;  x < cells_amount_x-1;  x++)
    {
    
      if (array_spin_exist_d_u_l_r[2*(x-1)+1] == false )
      {
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {chstr_lattice[get_cell_dl_pos_with_shift(0, x, y_, x_)]=' ';}
        }      
        
        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(1, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }   
      }
      
      if (array_spin_exist_d_u_l_r[2*(cells_amount_x-2)+2*(x-1)+1] == false )
      {
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {chstr_lattice[get_cell_dl_pos_with_shift(cells_amount_y-1, x, y_, x_)]=' ';}
        }
        
              
        for (int y_=0;  y_ < node_sz_y;  y_++)  
        {    
        symb=bond_bkg_symb;
          
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
          {
          chstr_lattice[get_pos_hor_line_d_yx_with_shft(1+cells_amount_y-2, x, y_, x_)]=symb;  
          //  if (symb == '*') {symb=' ';} else {symb='*';}
          }
        }           
      }      
      
    }


  //  edge empties (left, right)
  //
    for (int y=1;  y < cells_amount_y-1;  y++)
    {
    
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(y-1)+1] == false)
      {
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {chstr_lattice[get_cell_dl_pos_with_shift(y, 0, y_, x_)]=' ';}
        }
        
       symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, 1, y_, x_)]=symb;}
        }       
      }
      
      if (array_spin_exist_d_u_l_r[4*(cells_amount_x-2)+2*(cells_amount_y-2)+2*(y-1)+1] == false )
      {
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < cell_sz_in_x;  x_++)
          {chstr_lattice[get_cell_dl_pos_with_shift(y, cells_amount_x-1, y_, x_)]=' ';}
        }
        
      symb=bond_bkg_symb;
      
        for (int y_=0;  y_ < cell_sz_in_y;  y_++)
        {      
          for (int x_=0;  x_ < node_sz_x;  x_++)
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, 1+cells_amount_x-2, y_, x_)]=symb;}
        }       
                    
      }      
      
    }



     
  // marking of flip bond  
  //    
    if (is_1ver_0hor == false)
    {
    //  ver lines (=hor bonds)
    //    
    //  J_bond_y_index, int J_bond_x_index
    //
    int y=J_bond_y_index+1; 
    int x=J_bond_x_index+1;     cout<<"================================================================================aaaaaaaaaaaaaaaaaaaaaaa"<<endl<<"YYY="<<y<<"  "<<"XXX="<<x<<"  "<<endl;
    
    //for (int y=0;  y < cells_amount_y;  y++)
      {
        //for (int x=1;  x < cells_amount_x;  x++)
        {
        symb=bond_bkg_symb;
      
          for (int y_=0;  y_ < cell_sz_in_y;  y_++)
          {      
            for (int x_=0;  x_ < node_sz_x;  x_++)
            {
            chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, y_, x_)]=bond_separ_symb;          
            }
          }
                     
          //if (J_hor_[y*(cell_sz_in_x+1)+ x] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
          //{chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}

          if (J_hor[J_bond_y_index*(cells_amount_x-1)+ J_bond_x_index] == true) {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='+';} else 
          {chstr_lattice[get_pos_ver_line_l_yx_with_shft(y, x, cell_sz_in_y/2,  node_sz_x/2)]='-';}                               
                                         
        } 
      }
    }      


    if (is_1ver_0hor == true)
    {
    // hor lines (=ver bonds)
    //
    //  J_bond_y_index, int J_bond_x_index
    //           
    int y=J_bond_y_index+1;
    int x=J_bond_x_index+1;
    
      //for (int y=1;  y < cells_amount_y;  y++)
      {
        //for (int x=0;  x < cells_amount_x;  x++)
        {

          for (int y_=0;  y_ < node_sz_y;  y_++)  
          {    
          symb=bond_bkg_symb;
          
            for (int x_=0;  x_ < cell_sz_in_x;  x_++)  
            {
            chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, y_, x_)]=bond_separ_symb;  
            //  if (symb == '*') {symb=' ';} else {symb='*';}
            }
          }   
     
          //if (J_ver_[y*(cell_sz_in_x)+cell_sz_in_x] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
          //{chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}
               
          if (J_ver[J_bond_y_index*(cells_amount_x-2)+J_bond_x_index] == true) {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='+';} else 
          {chstr_lattice[get_pos_hor_line_d_yx_with_shft(y, x, node_sz_y/2,  cell_sz_in_x/2)]='-';}                                       
                                            
        }
         
      }
      
    }     




      
  //  draw empties in corners   
  //
  
    for (int y_=0;  y_ < cell_sz_in_y;  y_++)
    {      
      for (int x_=0;  x_ < cell_sz_in_x;  x_++)
      { 
      chstr_lattice[get_cell_dl_pos_with_shift(0,                0, y_, x_)]=' ';  chstr_lattice[get_cell_dl_pos_with_shift(0,                cells_amount_x-1, y_, x_)]=' '; 
      chstr_lattice[get_cell_dl_pos_with_shift(cells_amount_y-1, 0, y_, x_)]=' ';  chstr_lattice[get_cell_dl_pos_with_shift(cells_amount_y-1, cells_amount_x-1, y_, x_)]=' ';                
      }
    }
// get_cell_dl_pos_with_shift(int cell_num_y, int cell_num_x, int shft_y, int shft_x)


    if (is_do_out_frame == true) {draw_frame();} 
      
 
       
      
  draw_line_feeds();   
  print_to_screen();      // draw_nodes();    

    if (spin_array_cur   ) {delete[] spin_array_cur;     spin_array_cur=0;   }  
    if (spin_array_to_cmp) {delete[] spin_array_to_cmp;  spin_array_to_cmp=0;}  

      
  }      //    if ((chstr_lattice) && (J_ver) && (J_hor) && (spin_array_to_cmp_) && (spin_array_cur_))


}










bool extract_2D_rect_gbc_lattice_ext_edge_and_up_gs_row_from_file_v54(char * filename, int Ly, int Lx, bool * J_ver,  bool * J_hor,  bool & is_up_gs_found,  bool *& bool_edge_spin_exist_d_u_l_r_ar,  bool * bool_up_row_gs_ar, int *min_energ)
{ 

  if (filename == 0) {return false;}
  if (bool_edge_spin_exist_d_u_l_r_ar == 0) {return false;}


//  reading file
char * chstr_text=0;
//  func format: int write_file_to_chstr(char *filename, char *& text,  char *chstr_forbid_1st_line_symb_list,  int start_line_number=0,  bool not_take_empty_lines=true);
write_file_to_chstr(filename, chstr_text,  0);
char * chstr_extracting_phrase=0;
bool ok=false;
int returning_len=0;



//  -----------------
//  #  min energy
//  -125

  if (min_energ)
  {
  returning_len=extract_phrase_after_key_phr_in_x_lines_from_chstr(chstr_text,  (char *) "#  min energy",  chstr_extracting_phrase,  0,  1,  0);
    if (returning_len > 0) {ok=true;} else {ok=false;} 
    if (ok == true)
    {       
    ok=is_chstr_valid_num(chstr_extracting_phrase);
      if (ok == true)
      {
      *min_energ=atoi(chstr_extracting_phrase);
      }     
    }

    if (chstr_extracting_phrase) {delete[] chstr_extracting_phrase;  chstr_extracting_phrase=0;}
  }



//  -----------------
//  #  up row gs
//
  /*if (bool_up_row_gs_ar)
  {
  int line_number=0;  
  int amount_of_extracted_elements=0;
   
  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  up row gs", bool_up_row_gs_ar, Lx);
  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Lx)) {is_up_gs_found=true;} else {is_up_gs_found=false;}
   }*/






//  -----------------
//  #  bool *& bool_edge_spin_exist_d_u_l_r_ar
//
//
  {
  bool *bool_array_temp=0;
  bool_array_temp=new bool[Lx+Ly+1];
  
    for (int i=0;  i < Lx+Ly+1;  i++)
    {bool_array_temp[i]=false;}
    
    if (bool_edge_spin_exist_d_u_l_r_ar == 0) {bool_edge_spin_exist_d_u_l_r_ar=new bool[Ly*4+Lx*4];} 
      

  int amount_of_extracted_elements=0;


     
  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  down edge spin state", bool_array_temp, Lx);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Lx)) {      for (int i=0;  i < Lx;  i++) {bool_edge_spin_exist_d_u_l_r_ar[2*i]=bool_array_temp[i];}  }        

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  up edge spin state", bool_array_temp, Lx);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Lx)) {      for (int i=0;  i < Lx;  i++) {bool_edge_spin_exist_d_u_l_r_ar[2*Lx+2*i]=bool_array_temp[i];}      }

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  left edge spin state", bool_array_temp, Ly);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Ly)) {      for (int i=0;  i < Lx;  i++) {bool_edge_spin_exist_d_u_l_r_ar[2*Ly+2*Lx+2*i]=bool_array_temp[i];}      }

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  right edge spin state", bool_array_temp, Ly);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Ly)) {      for (int i=0;  i < Lx;  i++) {bool_edge_spin_exist_d_u_l_r_ar[2*Ly+4*Lx+2*i]=bool_array_temp[i];}      }
    
    

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  down edge spin existing", bool_array_temp, Lx);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Lx)) {      for (int i=0;  i < Lx;  i++) {bool_edge_spin_exist_d_u_l_r_ar[2*i+1]=bool_array_temp[i];}      }

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  up edge spin existing", bool_array_temp, Lx);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Lx)) {      for (int i=0;  i < Lx;  i++) {bool_edge_spin_exist_d_u_l_r_ar[2*Lx+2*i+1]=bool_array_temp[i];}      }

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  left edge spin existing", bool_array_temp, Ly);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Ly)) {      for (int i=0;  i < Lx;  i++) {bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*i+1]=bool_array_temp[i];}      }

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  right edge spin existing", bool_array_temp, Ly);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Ly)) {      for (int i=0;  i < Lx;  i++) {bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*Ly+2*i+1]=bool_array_temp[i];}      }
  

  
  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  down edge J_ver", bool_array_temp, Lx);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Lx)) {      for (int i=0;  i < Lx;  i++) {J_ver[i]=bool_array_temp[i];}      }

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  up edge J_ver", bool_array_temp, Lx);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Lx)) {      for (int i=0;  i < Lx;  i++) {J_ver[i+Lx*Ly]=bool_array_temp[i];}      }

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  left edge J_hor", bool_array_temp, Ly);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Ly)) {      for (int i=0;  i < Ly;  i++) {J_hor[i*(Lx+1)]=bool_array_temp[i];}      }

  amount_of_extracted_elements=extract_bool_array_from_chstr_text_in_line_after_key_phrase(chstr_text, (char *) "#  right edge J_hor", bool_array_temp, Ly);  
    if ((amount_of_extracted_elements > 0) && (amount_of_extracted_elements == Ly)) {      for (int i=0;  i < Ly;  i++) {J_hor[i*(Lx+1)+Lx]=bool_array_temp[i];}      }
  
    
    
    if (bool_array_temp) {delete[] bool_array_temp;  bool_array_temp=0;}
    
  }




  if (chstr_text) {delete[] chstr_text;  chstr_text=0;}

//  int extract_phrases_from_chstr(char *src_text, char ** & chstr_phrases, int & amount_of_extracted_phrases,  int & phrase_len,  int start_phrase_num,  int col_count=-1);
return true;

}






bool print_to_file_2D_rect_gbc_lattice_info_and_ext_edge_and_up_gs_row_v54(int Ly, int Lx, bool * J_ver,  bool * J_hor,  bool * bool_edge_spin_exist_d_u_l_r_ar, int *min_energ, bool *bool_up_row_gs_ar,  char * dest_filename_arg, int J_plus_perc, int J_min_perc,  bool is_write_lat_struct,  bool pbc1fbc0,  bool print_J_list, bool is_newfile)
{



char dest_filename[128];  dest_filename[0]='\0';  


//  0)  edit filename
//
  if (dest_filename_arg == 0)
  { 
  string result_data_filename;  result_data_filename.append((char *) "Ising_EA_2D_gbc_");
  result_data_filename.append((char *) "bond_pl");
  add_int_to_string(& result_data_filename, J_plus_perc);  
  result_data_filename.append((char *) "_min");
  add_int_to_string(& result_data_filename, J_min_perc);
  result_data_filename.append((char *) "_");  
 
    if (Ly < 10  ) {result_data_filename.append((char *) "0");}
    if (Ly < 100 ) {result_data_filename.append((char *) "0");}
    if (Ly < 1000) {result_data_filename.append((char *) "0");}             
  add_int_to_string(& result_data_filename, Ly);
    
  result_data_filename.append((char *) "x");
    
    if (Lx < 10  ) {result_data_filename.append((char *) "0");}
    if (Lx < 100 ) {result_data_filename.append((char *) "0");}
    if (Lx < 1000) {result_data_filename.append((char *) "0");}          
  add_int_to_string(& result_data_filename, Lx);
   
  result_data_filename.append((char *) ".txt");
  
  strcpy(dest_filename, result_data_filename.c_str());   
  result_data_filename.clear();
  }
  else
  {strcpy(dest_filename, dest_filename_arg);}



string string_var;


//  1)  write Ly
//  ------------
//
string_var.append((char *)  "#  Ly\n");
add_int_to_string(& string_var, Ly);
string_var.append((char *)  "\n");


//  2)  write Lx
//  ------------
//
string_var.append((char *)  "#  Lx\n");
add_int_to_string(& string_var, Lx);
string_var.append((char *)  "\n");


//  3)  write Ly
//  ------------
//
string_var.append((char *)  "#  J_plus_perc\n");
add_int_to_string(& string_var, J_plus_perc);
string_var.append((char *)  "\n");


//  4)  write Lx
//  ------------
//
string_var.append((char *)  "#  J_min_perc\n");
add_int_to_string(& string_var, J_min_perc);
string_var.append((char *)  "\n");

      
//  5)  lat and bond struct
//  ------------
//
  if (is_write_lat_struct == true)
  {
  char * chstr_temp=0;
  char elem='*';
  int size=0;
  string_var.append((char *)  "#  lattice and bond structure\n");
    if (pbc1fbc0 == true) {size=print_spin2D_lattice_with_bonds_to_chstr(chstr_temp, Ly, Lx, elem, 7);} else {size=print_spin2D_lattice_with_bonds_to_chstr(chstr_temp, Ly, Lx, elem, 0);}
    if (size > 0) {string_var.append(chstr_temp);}
    if (chstr_temp) {delete[] chstr_temp;  chstr_temp=0;}
  }


//  6)  lat and bond struct
//  ------------
//
  if (min_energ)
  {
  string_var.append((char *)  "#  min energy\n");
  add_int_to_string(& string_var, *min_energ);  // else {string_var.append((char *)  "?\n");}
  string_var.append((char *)  "\n");
  } 


//  7)  up row gs
//  ------------
//
  if (bool_up_row_gs_ar)
  {
  string_var.append((char *)  "#  up row gs\n");
    for (int index_x=0;  index_x < Lx;  index_x++)
    {    if (bool_up_row_gs_ar[index_x] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index_x < Lx-1) {string_var.append((char *)  " ");}    }
  string_var.append((char *)  "\n");
  }




// (char *) "#  down edge spin state"      2*index
// (char *) "#  up edge spin state"        2*Lx+2*index
// (char *) "#  left edge spin state"      4*Lx+2*index
// (char *) "#  right edge spin state"     4*Lx+2*Ly+2*index
 
// (char *) "#  down edge spin existing"   Lx+2*index+1
// (char *) "#  up edge spin existing"     2*Lx+2*index+1
// (char *) "#  left edge spin existing"   4*Lx+2*index+1
// (char *) "#  right edge spin existing"  4*Lx+2*Ly+2*index+1    
 
// (char *) "#  down edge J_ver"    index
// (char *) "#  up edge J_ver"      index+Ly*Lx
// (char *) "#  left edge J_hor"    index*(Lx+1)   
// (char *) "#  right edge J_hor"   index*(Lx+1)+Lx



//  8)  down edge spin state
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  down edge spin state\n");
    for (int index=0;  index < Lx;  index++)
    {    if (bool_edge_spin_exist_d_u_l_r_ar[2*index] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Lx-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  9)  up edge spin state
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  up edge spin state\n");
    for (int index=0;  index < Lx;  index++)
    {    if (bool_edge_spin_exist_d_u_l_r_ar[2*Lx+2*index] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Lx-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  10)  left edge spin state
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  left edge spin state\n");
    for (int index=0;  index < Ly;  index++)
    {    if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*index] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Ly-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  11)  right edge spin state
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  right edge spin state\n");
    for (int index=0;  index < Ly;  index++)
    {    if (bool_edge_spin_exist_d_u_l_r_ar[4*Lx+2*Ly+2*index] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Ly-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }





//  12)  down edge spin existing
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  down edge spin existing\n");
    for (int index=0;  index < Lx;  index++)
    {    if (bool_edge_spin_exist_d_u_l_r_ar[2*index+1] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Lx-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  13)  up edge spin existing
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  up edge spin existing\n");
    for (int index=0;  index < Lx;  index++)
    {    if (bool_edge_spin_exist_d_u_l_r_ar[2*Lx+2*index+1] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Lx-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  14)  left edge spin existing
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  left edge spin existing\n");
    for (int index=0;  index < Ly;  index++)
    {    if (bool_edge_spin_exist_d_u_l_r_ar[4*Ly+2*index+1] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Ly-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  15)  right edge spin existing
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  right edge spin existing\n");
    for (int index=0;  index < Ly;  index++)
    {    if (bool_edge_spin_exist_d_u_l_r_ar[4*Ly+2*Lx+2*index+1] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Ly-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }



//  16)  down edge J_ver
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  down edge J_ver\n");
    for (int index=0;  index < Lx;  index++)
    {    if (J_ver[index] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Lx-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  17)  up edge spin existing
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  up edge J_ver\n");
    for (int index=0;  index < Lx;  index++)
    {    if (J_ver[index+Ly*Lx] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Lx-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  14)  left edge spin existing
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  left edge J_hor\n");
    for (int index=0;  index < Ly;  index++)
    {    if (J_hor[index*(Lx+1)] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Ly-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }


//  15)  right edge spin existing
//  ------------
//
  if (bool_edge_spin_exist_d_u_l_r_ar)
  {
  string_var.append((char *)  "#  right edge J_hor\n");
    for (int index=0;  index < Ly;  index++)
    {    if (J_hor[index*(Lx+1)+Lx] == true) {string_var.append((char *)  "+1");} else {string_var.append((char *)  "-1");}    if (index < Ly-1) {string_var.append((char *)  " ");}  }
  string_var.append((char *)  "\n");
  }




write_chstr_to_file(string_var.c_str(), dest_filename, is_newfile);      //  !  !  !


//  16)  J-bond list
//  ------------
//
  if ((print_J_list == true) && (J_ver) && (J_hor))
  {
    if (pbc1fbc0  == true)
    {
    write_2D_rect_J_bond_table_to_file(dest_filename, Ly, Lx, J_ver, J_hor, 'c', 'c', 'c', 'c', 0, -1, 0, -1, 0, 0, true, false);  //  ! ! !  recording J array to file
    }
    else
    {
    write_2D_rect_J_bond_table_to_file(dest_filename, Ly, Lx, J_ver, J_hor, 'n', 'u', 'n', 'n', 0, -1, 0, -1, 0, 0, true, false);  //  ! ! !  recording J array to file    
    }
  }
  
  
  
//  cleaning
//
string_var.clear();



return true;

}






bool extract_lattice_gbc_data_v54(char * filename, int &Ly, int &Lx, bool * & J_ver,  bool * &J_hor,  int  & J_plus_perc, int & J_min_perc, bool & is_edges_extracted, bool & is_up_gs_found, bool *& bool_edge_spin_exist_d_u_l_r_ar,  bool *& bool_up_row_gs_ar,  int *min_energ, bool pbc1fbc0,  bool is_print_comments,  bool is_to_extract_up_row_gs)
{

  if (is_print_comments == true) {cout<<endl<<"start extracting lattice parameters from file ..."<<endl;}
  if (filename == 0) {return false;}
  
  
bool is_J_ver_initially_mem_allocated=false, is_J_hor_initially_mem_allocated=false, is_edge_inited=false,  is_up_row_gs_initially_mem_allocated=false;

  if (J_ver) {is_J_ver_initially_mem_allocated=true;}
  if (J_hor) {is_J_hor_initially_mem_allocated=true;}
  if (bool_edge_spin_exist_d_u_l_r_ar) {is_edge_inited=true;} 
  if (bool_up_row_gs_ar) {is_up_row_gs_initially_mem_allocated=true;} 
  




bool is_len_extracted=false;
//bool is_bond_ratio_extracted=false;
bool is_J_extracted=false;
//bool is_edges_extracted=false; 


//  extract Ly, Lx
//
  {
  int J_plus_perc_loc=50, J_min_perc_loc=50;
  is_len_extracted=extract_2D_rect_fbc_lattice_serv_info_from_file_v15(filename, & Ly, & Lx, 0, 0, 0, & J_plus_perc_loc, & J_min_perc_loc, 0);      //  !  !  !
  
    if (is_print_comments == true) {  if (is_len_extracted == true) {cout<<"is_len_extracted is extracted successfully: yes"<<endl;} else {cout<<"is_len_extracted is extracted successfully: no"<<endl;}  }
  J_plus_perc=J_plus_perc_loc;    J_min_perc=J_min_perc_loc;
  }


//  extract J_ver, J_hor
//
  if (is_len_extracted == true)
  {
  bool * J_ver_temp=0,  * J_hor_temp=0;         
  //J_ver_temp=new bool[(Ly+1)*Lx];  J_hor_temp=new bool[(Lx+1)*Ly];
  
  //  for (int k=0; k < (Ly+1)*Lx; k++) {J_ver_temp[k]=false;}
  //  for (int k=0; k < (Lx+1)*Ly; k++) {J_hor_temp[k]=false;}  
  
    if (J_ver == 0) {J_ver=new bool[(Ly+1)*Lx];}
    if (J_hor == 0) {J_hor=new bool[(Lx+1)*Ly];}

    for (int k=0; k < (Ly+1)*Lx; k++) {J_ver[k]=false;}
    for (int k=0; k < (Lx+1)*Ly; k++) {J_hor[k]=false;}  
 
    if (pbc1fbc0 == false)  {is_J_extracted=extract_2D_rect_lattice_J_table_from_file_dyn(filename, Ly, Lx, J_ver_temp, J_hor_temp, 'n', 'n', 'n', 'n');}      //  !  !  !
    if (pbc1fbc0 == true )  {is_J_extracted=extract_2D_rect_lattice_J_table_from_file_dyn(filename, Ly, Lx, J_ver_temp, J_hor_temp, 'c', 'c', 'c', 'c');}      //  !  !  !    
  
      if (is_J_extracted == true)
      {            
        for (int k=0; k < (Ly+1)*Lx; k++) {J_ver[k]=J_ver_temp[k];}
        for (int k=0; k < (Lx+1)*Ly; k++) {J_hor[k]=J_hor_temp[k];}
      }
      else
      {
        if ((is_J_ver_initially_mem_allocated == false) && (J_ver)) {delete[] J_ver;  J_ver=0;}
        if ((is_J_hor_initially_mem_allocated == false) && (J_hor)) {delete[] J_hor;  J_hor=0;}   
      }
    
    
      if (J_ver_temp) {delete[] J_ver_temp;  J_ver_temp=0;}
      if (J_hor_temp) {delete[] J_hor_temp;  J_hor_temp=0;}
          
  
    //if ((is_J_extracted == false) && (is_J_ver_initially_mem_allocated == false) && (J_ver)) {delete[] J_ver;  J_ver=0;}
    //if ((is_J_extracted == false) && (is_J_hor_initially_mem_allocated == false) && (J_hor)) {delete[] J_hor;  J_hor=0;}    

    if (is_print_comments == true) {  if (is_J_extracted == true) {cout<<"is_J_extracted successfully: yes"<<endl;} else {cout<<"is_J_extracted successfully: no"<<endl;}  }
  
  }


//  extract edge
//  
  if (is_J_extracted == true)
  {  
  
    if (bool_edge_spin_exist_d_u_l_r_ar == 0) {bool_edge_spin_exist_d_u_l_r_ar=new bool[Ly*4+Lx*4];}
    
    for (int i=0;  i < (Ly*4+Lx*4);  i++) 
    {bool_edge_spin_exist_d_u_l_r_ar[i]=false;}


    if (is_to_extract_up_row_gs == true)
    {  if (bool_up_row_gs_ar == 0) {bool_up_row_gs_ar=new bool[Lx];}  }
    
    if (bool_up_row_gs_ar)
    { 
      for (int i=0;  i < Lx;  i++) 
      {bool_up_row_gs_ar[i]=false;}
    }
    
  bool is_up_gs_found_loc=false;
  
  is_edges_extracted=extract_2D_rect_gbc_lattice_ext_edge_and_up_gs_row_from_file_v54(filename, Ly, Lx,  J_ver,  J_hor,   is_up_gs_found_loc, bool_edge_spin_exist_d_u_l_r_ar,  bool_up_row_gs_ar,  min_energ);

  
      if (is_print_comments == true) {  if (is_edges_extracted == true) {cout<<"is_edges_extracted: yes"<<endl;} else {cout<<"is_edges_extracted: no"<<endl;}  }
  
  is_up_gs_found=is_up_gs_found_loc;

  
    if (is_edges_extracted == false) 
    { 
      if (bool_edge_spin_exist_d_u_l_r_ar) {delete[] bool_edge_spin_exist_d_u_l_r_ar;  bool_edge_spin_exist_d_u_l_r_ar=0;}   
      if (is_up_row_gs_initially_mem_allocated == false) {  if (bool_up_row_gs_ar) {delete[] bool_up_row_gs_ar;  bool_up_row_gs_ar=0;}  }
    }
    
  }

  
  
  if (is_print_comments == true) {cout<<"process of extracting lattice parameters from file completed!"<<endl<<endl;}
  
return is_J_extracted;

}





