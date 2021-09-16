#ifndef __RDK_IO_FORMATTING_CPP
#define __RDK_IO_FORMATTING_CPP

//  Last modified 03.11.2020 PadalkoMA

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

#include "anyfunctions.cpp"


using namespace std;




//
//                        _________________________|                 
//  ---------------                    
//                |                      
//                 \____________|        
//                     |                               |
//     _______ _______________          ---------------
//                     




// ==========================================================================
// ==========================================================================

namespace _NMSP_RDK_
{

void free_chstr_array(char ** & chstr_var, int ext_array_size); 
  //  is ch in chstr,   if (len == -1) {strlen will be completed in proc}
inline bool is_matching(char ch,  const char *const chstr, int len=-1);
inline bool copy_replace(char * chstr_src,  int src_start_pos,  int len,     char * chstr_dest,  int dest_start_pos=0,    bool add_end_str_symb=false);
inline int calc_max_len_of_phrase_in_chstr(char ** text, int amount);
inline int chstr_fill(char *chstr, char symb, int start_pos, int len, bool put_end_line_symb=true);      //  proc returns cur_pos  
inline int string_append_x_times(string *str, char* add, int x);      //  proc returns length
int delete_symbs_in_chstr(char *chstr, char *excess_symbs_list);


inline bool is_ch_natural_num(char ch);
inline bool is_ch_integer_num(char ch);
inline bool is_ch_real_num_elem(char ch, char ch_mantis_separ1='.',  char ch_mantis_separ2='.');
bool is_chstr_valid_num(char *chstrnum, char ch_mantis_separ1='.',  char ch_mantis_separ2='.', int max_len=32);
int put_mantis_exp_to_chstrnum_stat(char * chstr_dest, char *chstr_mantis, char *chstr_exp);

int sr_altern_1sp_2phrase__pos2phrase(const char *const chstr,  int start_pos, bool * is_found=0);
int sr_altern_1phrase_2sp__pos2sp(const char *const chstr,  int start_pos, bool * is_found=0);
int sr_altern_1phrase_2sp__lastpos1phrase(const char *const chstr,  int start_pos, bool * is_found=0);

int sr_altern_1sp_2num__pos2num(const char *const chstr,  int start_pos, bool * is_found=0,  char *symbs_consid_as_num=(char *) "+-.,E\0");
int sr_altern_1num_2sp__pos2sp(const char *const chstr, int start_pos, bool * is_found=0,  char *symbs_consid_as_num=(char *) "+-.,E\0");
int sr_altern_1num_2sp__lastpos1num(const char * const chstr, int start_pos, bool * is_found=0,  char *symbs_consid_as_num=(char *) "+-.,E\0");    

bool sr_given_phrase_in_chstr(const char * const chstr, char * chstr_phrase, int start_pos=0, int *returning_start_pos=0, int *returning_end_pos=0);

int calc_amount_of_lines_in_chstr(char *text);
int calc_amount_of_lines_in_chstr(char *text,  char *chstr_forbid_1st_line_symb_list);
int calc_amount_of_lines_in_chstr(char *text,  char *chstr_forbid_1st_line_symb_list, int start_line_number);
int calc_amount_of_not_empt_lines_in_chstr(char *text);
int calc_amount_of_not_empt_lines_in_chstr(char *text, int start_line_number);
int calc_max_line_size_in_chstr(char *text);
int calc_amount_of_lines_in_file(char *filename);
int calc_amount_of_not_empt_lines_in_file(char *filename);
int calc_amount_of_not_empt_lines_in_file(char *filename, int start_line_number);
int calc_max_line_size_in_file(char *filename); 
int calc_num(char *filename);


  //  line_size is calculated and returned by proc  
bool sr_line_position_and_len_in_chstr(char *chstr,  int line_number, int * line_start_pos=0, int * line_end_pos=0, int * line_len=0);
int  sr_line_number_for_line_with_given_phrase_from_chstr(char *chstr, char *chstr_key_phrase, int start_pos=0);      // returns -1 if not found
int  extract_line_with_given_number_from_chstr(char *chstr, int line_number, char * & chstr_extracting_line);
int  extract_phrase_in_given_area_from_chstr(char *chstr,  int given_start_pos,  int given_end_pos,  char * & chstr_extracting_phrase, int *extr_phr_start_pos=0, int *extr_phr_end_pos=0);
int  extract_phrase_after_key_phr_in_x_lines_from_chstr(char *chstr,  char *chstr_key_phrase,  char * & chstr_extracting_phrase,  int start_pos=0,  int d_x_lines=0,  int phrase_number_in_line=0);
//int  extract_text_from_chstr_after_key_phrase_untill_next_key_phr(char *chstr,  char *chstr_key_phrase,  char * & chstr_extracting_phrase,  int start_pos=0);      //  not done
int extract_line_after_line_with_given_phrase_from_chstr(char *chstr, char *chstr_key_phrase, char *& extracting_phrase, int start_pos=0);      //  still not checked!
int extract_lines_from_file(char *filename, char ** & chstr_lines, int & amount_of_lines,  int & max_line_size,  int start_line_number, int line_wished_amount_be_extracted);
int extract_lines_as_text_from_file(char *filename, char * & text, int & amount_of_lines,  int start_line_number, int line_wished_amount_be_extracted);
int extract_lines_from_file(char *filename, char * & chstr_lines, int & amount_of_lines,  int & line_size,  int start_line_number=0);
int extract_not_empty_lines_from_file(char *filename, char *chstr_lines, int & amount_of_not_empt_lines,  int & line_size);
int extract_lines_from_file(char *filename, char ** & chstr_lines, int & amount_of_lines,  int & line_size,  int start_line_number=0);
int extract_lines_from_file(char *filename, char ** & chstr_lines, int & amount_of_lines,  int & line_size,  char *chstr_forbid_1st_line_symb_list,  int start_line_number=0);

int calc_phrases_after_key_phr_in_x_lines_from_chstr(char *chstr,  char *chstr_key_phrase,  int start_pos=0,  int d_x_lines=0); 
int calc_amount_of_phrases_in_chstr(char * text);
int calc_max_len_of_phrase_in_chstr(char * text);
int calc_max_len_of_phrases_in_phrase_ar(char ** chstr_phrases, int phrases_count);
int calc_max_len_of_phrases_in_chstr_by_col(char ** chstr_phrases, int phrases_count, int col_count, int * & max_phrase_len);      // proc returns sum of max len of phrases in row
bool sr_if_there_signs_bef_phrases_in_chstr_by_col(char ** chstr_phrases, int phrases_count, int col_count, bool * & sign_is);      // proc returns true if there is at least 1 + or -, it's needed for sign align
int extract_phrases_from_file(char *filename, char ** & chstr_phrase, int & amount_of_phrases,  int & phrase_len,  int start_line_number=0);

int calc_amount_of_anytypenums_in_chstr(const char * const text, char * symbs_consid_as_num=(char *) "+-.,E\0");
int calc_max_len_of_anytypenums_in_chstr(char * text, char *symbs_consid_as_num=(char *) "+-.,E\0");
int extract_anytypenums_from_file(char *filename, char ** & chstr_nums, int & amount_of_nums,  int & num_len,  int start_line_number=0, char *symbs_consid_as_num=(char *) "+-.,E\0");


void display_chstr_lines(char * chstr_var,  int amount_of_lines,  int & line_size);
void display_chstr_lines(char ** chstr_var,  int amount_of_lines);
void display_chstr_phrases_by_col(char ** chstr_var,  int phrases_count, int col_count, char *chstr_separ=(char *) "\t\t");
int write_chstr_phrases_to_file_by_col(char * filename, char ** chstr_var,  int phrases_count, int col_count, char *chstr_separ=(char *) "\t\t", bool B1New0Add=true);

int write_chstr_phrases_to_chstr_by_align_col_1(char ** chstr_phrases, char * & chstr_dest,  int phrases_count, int col_count, int hor_interval_len, int ver_interval_len);                  // dif cols
int write_chstr_phrases_to_file_by_align_col_1(char * filename,  char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, bool B1New0Add=true);  // dif cols
int write_chstr_phrases_to_screen_by_align_col_1(char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len);                                       // dif cols
int write_chstr_phrases_to_chstr_by_align_col_2(char ** chstr_phrases, char * & chstr_dest,  int phrases_count, int col_count, int hor_interval_len, int ver_interval_len);                  // same cols
int write_chstr_phrases_to_file_by_align_col_2(char * filename,  char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, bool B1New0Add=true);  // same cols
int write_chstr_phrases_to_screen_by_align_col_2(char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len);                                       // same cols
int write_chstr_phrases_to_chstr_by_align_col_3(char ** chstr_phrases, char * & chstr_dest,  int phrases_count, int col_count, int hor_interval_len, int ver_interval_len);                  // same cols, sign align 
int write_chstr_phrases_to_file_by_align_col_3(char * filename,  char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, bool B1New0Add=true);  // same cols, sign align
int write_chstr_phrases_to_screen_by_align_col_3(char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len);                                       // same cols, sign align
int write_chstr_phrases_to_chstr_by_align_col_4(char ** chstr_phrases, char * & chstr_dest,  int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, char type=1);     // with frame
int write_chstr_phrases_to_file_by_align_col_4(char * filename,  char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, bool B1New0Add=true, char type=1); // with frame 
int write_chstr_phrases_to_screen_by_align_col_4(char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, char type=1);                          // with frame


int calc_phrases_amount_in_line(char *chstr_text, int line_number);
int extract_phrases_from_chstr(char *src_text, char ** & chstr_phrases, int & amount_of_extracted_phrases,  int & phrase_len,  int start_phrase_num,  int col_count=-1); 
//  super SAFE and CHECKED proc for phrase extracting.  3 main arrays for filling: phrases, lens, startposes.   max_phrase_len - if you want to know max phrase len   all phrases we consider like table with cols  
int extract_phrases_from_chstr(char *src_text, char * &chstr_phrases, int *&len_array, int *&start_pos_array, int * max_phrase_len=0,  int start_phrase_num=0,  int col_count=1, int src_start_pos=0, int src_end_pos=-1);

//  the most needed
//  -------------______________****** ------------>  
bool extract_column_from_file(char * filename, char **&col_phrases, int & amount_in_col,  int & max_phrase_len, int col_start_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_number=0);
bool write_columns_to_screen_3(int amount_in_col, int dy, int dx, bool B1New0Add, char ** col1_phrases=0, char ** col2_phrases=0, char ** col3_phrases=0, char ** col4_phrases=0, char ** col5_phrases=0);
bool write_columns_to_screen_4(int amount_in_col, int dy, int dx, bool B1New0Add, char ** col1_phrases=0, char ** col2_phrases=0, char ** col3_phrases=0, char ** col4_phrases=0, char ** col5_phrases=0);
int write_columns_to_file_3(char * filename, int amount_in_col, int dy, int dx, bool B1New0Add, char ** col1_phrases=0, char ** col2_phrases=0, char ** col3_phrases=0, char ** col4_phrases=0, char ** col5_phrases=0);
int write_columns_to_file_4(char * filename, int amount_in_col, int dy, int dx, bool B1New0Add, char ** col1_phrases=0, char ** col2_phrases=0, char ** col3_phrases=0, char ** col4_phrases=0, char ** col5_phrases=0);


template<typename T> bool extract_num_array_from_chstr(char * chstr_text,  T *& num_array, int & amount_of_num,  bool type_dbl1_int0=false);
bool extract_bool_array_from_chstr(char * chstr_text,  bool *& num_array, int & amount_of_num);


template< typename T > bool extract_dbl_column_from_file(char * filename, T *&col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);
template< typename T > bool extract_int_column_from_file(char * filename, T *&col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);


template<typename T> bool write_dbl_array_to_chstr(T * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len, int prec=3, int wished_len=8);
template<typename T> bool write_int_array_to_chstr(T * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len, int wished_len=-1);


void write_chstr_phrases_to_screen(char ** col_phrases, int amount_of_elem, int start_elem_num=0, bool hor1ver0=true, bool do_enum=false);
template<typename T> void write_dbl_array_to_screen(T * col_phrases, int amount_of_elem, int start_elem_num=0, bool hor1ver0=true, bool do_enum=false);
template<typename T> void write_int_array_to_screen(T * col_phrases, int amount_of_elem, int start_elem_num=0, bool hor1ver0=true, bool do_enum=false);


//  adding coulomn at the position col_num near other columns,  very important procs!!!
//   col_num-=1 means to put at the end
// 
int add_chstr_array_column_to_file(char ** col_phrases, char * filename, int amount_in_col,  int start_line_num=0, char *forbid_1st_symbs=(char *) "#/", int col_num=-1, char type=0);
template<typename T> int add_dbl_array_column_to_file(T *col_phrases, char *filename, int amount_in_col, int prec=3, int wished_len=8, int start_line_num=0, char *forbid=(char *) "#/", int col_num=-1, char type=0);
template<typename T> int add_int_array_column_to_file(T *col_phrases, char * filename, int amount_in_col, int wished_len=-1, int start_line_num=0, char *forbid_1st_symbs=(char *) "#/", int col_num=-1, char type=0);





int extract_bool_array_from_chstr_text_in_line_after_key_phrase(char * chstr_text, char * key_phrase, bool * const bool_array, int given_amount_of_array_elements=-1);







//  ========        ========          ========           =======        ======
//  --------------------------------------------------------------------------
//
//
//
//
//
//                                   ===========================================
//   ------------------------------------------------------------------
//                   \_____________
//                                 /
//                                /
// ------------------------------/
//
//
//
//              ========        ========          ========           =======        ======













void free_chstr_array(char ** & chstr_var,  int ext_array_size)
{

  if (chstr_var)
  {
    for (int i=0;  i < ext_array_size;  i++)
    {
      if (chstr_var[i]) {delete[] chstr_var[i];  chstr_var[i]=0;}
    }
 
  delete[] chstr_var;  chstr_var=0;
  }

}




// ==========================================================================
// ==========================================================================



inline bool is_matching(char ch,  const char * const chstr, int len)
{

  if (chstr == 0) {return false;}

  if (len < 0) {len=strlen(chstr);}

  for (int i=0;  i < len;  i++)
  {  if (chstr[i]  == ch) {return true;}  }

return false;

}



// ==========================================================================
// ==========================================================================



inline bool copy_replace(char * chstr_src,  int src_start_pos,  int len,     char * chstr_dest,  int dest_start_pos,    bool add_end_str_symb)
{

  if ((chstr_src == 0) || (chstr_dest == 0)) {return false;} 

  for (int pos=0;  pos < len;  pos++)
  {chstr_dest[dest_start_pos+pos]=chstr_src[src_start_pos+pos];}

  if (add_end_str_symb == true) {chstr_dest[dest_start_pos+len]='\0';}

return true;

}



// ==========================================================================
// ==========================================================================



inline int calc_max_len_of_phrase_in_chstr(char ** text, int amount)
{

  if ((text == 0) || (amount < 1)) {return 0;}

int len=0, maxlen=0;

  for (int i=0;  i < amount;  i++)
  {
  len=strlen(text[i]);
    if (maxlen < len) {maxlen=len;}
  }

return maxlen;

}



// ==========================================================================
// ==========================================================================



inline int chstr_fill(char *chstr, char symb, int start_pos, int len, bool put_end_line_symb)
{

  if ((chstr == 0) || (start_pos < 0) || (len <= 0)) {return start_pos;}

  for (int i=0;  i < len;  i++)
  {chstr[start_pos+i]=symb;}

  if (put_end_line_symb == true) {chstr[start_pos+len]='\0';}

return start_pos+len;

}



// ==========================================================================
// ==========================================================================



inline int string_append_x_times(string *str, char* add, int x)
{

  if (str == 0) {return 0;}
  if (add == 0) {return str->length();}

  for (int i=0;  i < x;  i++)
  {str->append(add);}

return str->length();

}



// ==========================================================================
// ==========================================================================



int delete_symbs_in_chstr(char *chstr, char *excess_symbs_list)
{

  if ((chstr == 0) || (excess_symbs_list == 0)) {return 0;}

int chstr_len=strlen(chstr);
  if (chstr_len < 1) {return 0;}
int excess_symbs_list_len=strlen(excess_symbs_list);
  if (excess_symbs_list_len < 1) {return strlen(chstr);}


int shift=0;

  for (int pos=0;  pos < chstr_len;  pos++)
  {
    if (is_matching(chstr[pos], excess_symbs_list, excess_symbs_list_len) == true)  {shift++;}
    else {chstr[pos-shift]=chstr[pos];}
  }

chstr_len=chstr_len-shift;
chstr[chstr_len]='\0';

return chstr_len;

}


// ==========================================================================
// ==========================================================================



inline bool is_ch_natural_num(char ch)
{

  if ((ch > 47) && (ch < 58)) {return true;}

return false;

}



// ==========================================================================
// ==========================================================================



inline bool is_ch_integer_num(char ch)
{

  if (((ch > 47) && (ch < 58)) || (ch == '+') || (ch == '-')) {return true;}

return false;

}



// ==========================================================================
// ==========================================================================



inline bool is_ch_real_num_elem(char ch, char ch_mantis_separ1,  char ch_mantis_separ2)
{

  if (((ch > 47) && (ch < 58)) || (ch == '.') || (ch == ch_mantis_separ1) || (ch == ch_mantis_separ2) || (ch == '+') || (ch == 'e') || (ch == 'E')) {return true;}

return false;

}


// ==========================================================================
// ==========================================================================


bool is_chstr_valid_num(char *chstrnum, char ch_mantis_separ1,  char ch_mantis_separ2, int max_len)
{

  if (chstrnum ==0) {return false;}

int len=strlen(chstrnum);
int sepsymbcount=0, esymbcount=0, sign1symbcount=0, sign2symbcount=0;
int sepsymbpos=0, esymbpos=0, sign1symbpos=0, sign2symbpos=0;
int cyfsymbcount=0;
bool isbadsymbol=false;

  if (len > max_len) {return false;}
  if (len < 1 ) {return false;}

  for (int i=0; i < len; i++)
  {
    if ((chstrnum[i] <= 47) || (chstrnum[i] >= 58))
    {
      if ((chstrnum[i] == '-') || (chstrnum[i] == '+'))
      {
        if (sign1symbcount < 1) {sign1symbpos=i;  sign1symbcount++;}
        else {sign2symbpos=i;  sign2symbcount++;}
      }
      else
      {
        if ((chstrnum[i] == ch_mantis_separ1) || (chstrnum[i] == ch_mantis_separ2))
        {
        sepsymbpos=i;  sepsymbcount++;
        }
        else
        {
          if ((chstrnum[i] == 'e') || (chstrnum[i] == 'E'))
          {
          esymbpos=i;  esymbcount++;
          }
          else
          {
          //isbadsymbol=true;
          return false;
          }
        }
      }
    }
    else
    {cyfsymbcount++;}
  } 


  if (cyfsymbcount   < 1) {return false;}
  if (sign1symbcount > 1) {return false;}
  if (sign2symbcount > 1) {return false;}
  if (sepsymbcount   > 1) {return false;}
  if (esymbcount     > 1) {return false;}

  if ((len > 1) && (sepsymbpos == len-1))   {return false;}
  if ((esymbpos != 0) && (sepsymbpos > esymbpos)) {return false;}
  if ((sepsymbpos > esymbpos) && (esymbcount > 0)) {return false;}
  if ((sign1symbcount > 0) && (sign1symbpos != 0)) {return false;}
  if ((sign2symbcount > 0) && (esymbcount > 0) && (sign2symbpos != esymbcount+1)) {return false;}
  if ((esymbcount     > 0) && (esymbpos == 0)) {return false;}
  if ((esymbcount     > 0) && (esymbpos == len-1)) {return false;}

return true;

}



// ==========================================================================
// ==========================================================================



int sr_altern_1sp_2phrase__pos2phrase(const char * const chstr,  int start_pos, bool * is_found)
{

bool is_found_loc=false;
int pos=start_pos;
bool is_prev_symb_sp=false;
//int chstr_len=strlen(chstr);

  if (start_pos == 0) {is_prev_symb_sp=true;}

  while (chstr[pos] != '\0')  //&& (chstr_len > pos))
  {
    if ((chstr[pos] == ' ') || (chstr[pos] == '\n') || (chstr[pos] == '\t')) 
    {is_prev_symb_sp=true;}
    else
    {  if (is_prev_symb_sp == true) {is_found_loc=true;  break;}   is_prev_symb_sp=false;}

  pos++;
  }

  if (is_found) {*is_found=is_found_loc;}

return pos; 

}



// ==========================================================================
// ==========================================================================



bool sr_given_phrase_in_chstr(const char * const chstr, char * given_phrase, int start_pos, int *returning_start_pos, int *returning_end_pos)
{

  if (chstr == 0) {return false;}

  if ((returning_start_pos == 0) && (returning_end_pos == 0)) {return false;}

  if (given_phrase == 0) {return false;}

  if (chstr == 0) {  if (returning_start_pos) {*returning_start_pos=-1;}   if (returning_end_pos) {*returning_end_pos=-2;}  return false;}

  if (start_pos < 0) {start_pos=0;}


int pos=start_pos;
int match_symb_counter=0;
int phrase_len=0;
phrase_len=strlen(given_phrase);
int proposed_start_pos_of_phrase=start_pos;

  if (returning_start_pos) {*returning_start_pos=-1;}   
  if (returning_end_pos) {*returning_end_pos=-2;} 

  if (phrase_len < 1) {return false;}


  while (chstr[pos] != '\0')
  {
    if (chstr[pos] != given_phrase[match_symb_counter]) {match_symb_counter=0;}
    else 
    {
      if (match_symb_counter == 0) {proposed_start_pos_of_phrase=pos;}
    match_symb_counter++;
      if (match_symb_counter >= phrase_len) 
      {
        if (returning_start_pos != 0) {*returning_start_pos=proposed_start_pos_of_phrase;}
        if (returning_end_pos   != 0) {*returning_end_pos  =pos;                         }
      return true;
      }      
    }
  pos++;  
  }


return false;

}



// ==========================================================================
// ==========================================================================



int put_mantis_exp_to_chstrnum_stat(char * chstr_dest, char *chstr_mantis, char *chstr_exp)
{

  if (chstr_dest == 0) {return 0;}
  if ((chstr_mantis == 0) && (chstr_exp == 0)) {return 0;}

chstr_dest[0]='\0';
strcpy(chstr_dest,  chstr_mantis);

  if (chstr_exp)
  {   
  int exp_len=strlen(chstr_exp);
  bool exp_is=false;
    if ((exp_len > 0) && (chstr_exp[0] == 'e')) {exp_is=true;}
    if ((exp_len > 0) && (chstr_exp[0] == 'E')) {exp_is=true;}
    if ((exp_len > 1) && (chstr_exp[1] == 'e')) {exp_is=true;}
    if ((exp_len > 1) && (chstr_exp[1] == 'E')) {exp_is=true;}

    if (exp_is == false) 
    {
    strcat(chstr_dest,  (char *) "E");
 
      if ((chstr_exp[0] != '+') && (chstr_exp[0] != '-')) {strcat(chstr_dest,  (char *) "+");}

    }

  }

strcat(chstr_dest,  chstr_exp);

return strlen(chstr_dest);

}



// ==========================================================================
// ==========================================================================



int sr_altern_1phrase_2sp__pos2sp(const char * const chstr,  int start_pos, bool * is_found)
{

bool is_found_loc=false;
int pos=start_pos;
bool is_prev_symb_phrase=false;
//int chstr_len=strlen(chstr);

  while (chstr[pos] != '\0')  //&& (chstr_len > pos))
  {
    if ((chstr[pos] != ' ') && (chstr[pos] != '\n') && (chstr[pos] != '\t')) 
    {is_prev_symb_phrase=true;}
    else
    {  if (is_prev_symb_phrase == true) {is_found_loc=true;  break;}   is_prev_symb_phrase=false;}

  pos++;
  }

  if (is_found) {*is_found=is_found_loc;}

return pos; 

}



// ==========================================================================
// ==========================================================================



int sr_altern_1phrase_2sp__lastpos1phrase(const char * const chstr,  int start_pos, bool * is_found)
{

bool is_found_loc=false;
int pos=start_pos;
int chstr_len=strlen(chstr);

  while (chstr[pos] != '\0')  //&& (chstr_len > pos))
  {
    if ((chstr[pos] != ' ') && (chstr[pos] != '\n') && (chstr[pos] != '\t')) 
    {
      if ((chstr[pos+1] == '\0') || ((chstr[pos+1] == ' ') || (chstr[pos+1] == '\n') || (chstr[pos+1] == '\t'))) {is_found_loc=true;  break;}
    }
  pos++;
  }

  if (is_found) {*is_found=is_found_loc;}

return pos; 

}




// ==========================================================================
// ==========================================================================



int sr_altern_1sp_2num__pos2num(const char *const chstr,  int start_pos, bool * is_found, char *symbs_consid_as_num)
{

bool is_found_loc=false;
int pos=start_pos;
bool is_prev_symb_sp=false;
int symbs_consid_as_num_len=0;
//int chstr_len=strlen(chstr);

  if (symbs_consid_as_num) {symbs_consid_as_num_len=strlen(symbs_consid_as_num);}
  if (start_pos == 0) {is_prev_symb_sp=true;}

  while (chstr[pos] != '\0')  //&& (chstr_len > pos))
  {
    if ((is_ch_natural_num(chstr[pos]) == false) && (is_matching(chstr[pos], symbs_consid_as_num, symbs_consid_as_num_len) == false))
    {is_prev_symb_sp=true;}
    else
    {  if (is_prev_symb_sp == true) {is_found_loc=true;  break;}   is_prev_symb_sp=false;}

  pos++;
  }

  if (is_found) {*is_found=is_found_loc;}

return pos; 

}



// ==========================================================================
// ==========================================================================



int sr_altern_1num_2sp__pos2sp(const char *const chstr, int start_pos, bool * is_found, char *symbs_consid_as_num)
{

bool is_found_loc=false;
int pos=start_pos;
bool is_prev_symb_num=false;
int symbs_consid_as_num_len=0;
//int chstr_len=strlen(chstr);

  if (symbs_consid_as_num) {symbs_consid_as_num_len=strlen(symbs_consid_as_num);}
  if (start_pos == 0) {is_prev_symb_num=true;}


  while (chstr[pos] != '\0')  //&& (chstr_len > pos))
  {
    if ((is_ch_natural_num(chstr[pos]) == true) || (is_matching(chstr[pos], symbs_consid_as_num, symbs_consid_as_num_len) == true))
    {is_prev_symb_num=true;}
    else
    {  if (is_prev_symb_num == true) {is_found_loc=true;  break;}   is_prev_symb_num=false;}

  pos++;
  }

  if (is_found) {*is_found=is_found_loc;}

return pos; 

}



// ==========================================================================
// ==========================================================================



int sr_altern_1num_2sp__lastpos1num(const char *const chstr, int start_pos, bool * is_found,  char *symbs_consid_as_num)
{

bool is_found_loc=false;
int pos=start_pos;
int symbs_consid_as_num_len=0;
//int chstr_len=strlen(chstr);

  if (symbs_consid_as_num) {symbs_consid_as_num_len=strlen(symbs_consid_as_num);}

  while (chstr[pos] != '\0')  //&& (chstr_len > pos))
  {
    if ((is_ch_natural_num(chstr[pos]) == true) || (is_matching(chstr[pos], symbs_consid_as_num, symbs_consid_as_num_len) == true)) 
    {
      if ((chstr[pos+1] == '\0') || ((is_ch_natural_num(chstr[pos+1]) == false) && (is_matching(chstr[pos+1], symbs_consid_as_num, symbs_consid_as_num_len) == false))) 
      {is_found_loc=true;   break;}
    }
  pos++;
  } 

  if (is_found) {*is_found=is_found_loc;}

return pos; 

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_lines_in_chstr(char *text)
{

  if (text == 0) {return 0;} 

int text_size=strlen(text); 
int line_counter=0, lastpos=0;

  for (int pos=0;  pos < text_size;  pos++)
  {
    if (text[pos] == '\n') {lastpos=pos;  line_counter++;} 
  }

  if ((text_size > 0) && (text[text_size-1] != '\n')) {line_counter++;}

return line_counter;

}



// ==========================================================================
// ==========================================================================




int calc_amount_of_lines_in_chstr(char *text,  char *chstr_forbid_1st_line_symb_list)
{

  if (text == 0) {return 0;} 

int text_size=strlen(text); 
int line_counter=0, lastpos=0;
int amount_forbid_symb=0,  amount_forbid_lines=0;

  if (chstr_forbid_1st_line_symb_list) {amount_forbid_symb=strlen(chstr_forbid_1st_line_symb_list);}

  if (text_size > 0)
  {        for (int i=0;  i < amount_forbid_symb;  i++)  {  if (text[0] == chstr_forbid_1st_line_symb_list[i]) {amount_forbid_lines++;  break;}  }        }

  for (int pos=0;  pos < text_size;  pos++)
  {
    if (text[pos] == '\n') 
    {
    lastpos=pos;  line_counter++;
      if (pos+1 <= text_size-1)
      {        for (int i=0;  i < amount_forbid_symb;  i++)  {  if (text[pos+1] == chstr_forbid_1st_line_symb_list[i]) {amount_forbid_lines++;  break;}  }        }
    }
  }

  if ((text_size > 0) && (text[text_size-1] != '\n')) {line_counter++;}

return line_counter-amount_forbid_lines;

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_lines_in_chstr(char *text,  char *chstr_forbid_1st_line_symb_list, int start_line_number)
{

  if (text == 0) {return 0;} 

int text_size=strlen(text); 
int line_counter=0, lastpos=0;
int amount_forbid_symb=0,  amount_forbid_lines=0;

  if (chstr_forbid_1st_line_symb_list) {amount_forbid_symb=strlen(chstr_forbid_1st_line_symb_list);}

  if ((text_size > 0) && (line_counter >= start_line_number))
  {        for (int i=0;  i < amount_forbid_symb;  i++)  {  if (text[0] == chstr_forbid_1st_line_symb_list[i]) {amount_forbid_lines++;  break;}  }        }

  for (int pos=0;  pos < text_size;  pos++)
  {
    if (text[pos] == '\n') 
    {
    lastpos=pos; 
      if (line_counter >= start_line_number)
      { 
      line_counter++;
        if (pos+1 <= text_size-1)
        {        for (int i=0;  i < amount_forbid_symb;  i++)  {  if (text[pos+1] == chstr_forbid_1st_line_symb_list[i]) {amount_forbid_lines++;  break;}  }        }
      }
    }
  }

  if ((text_size > 0) && (text[text_size-1] != '\n')) {line_counter++;}

return line_counter-amount_forbid_lines;

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_not_empt_lines_in_chstr(char *text)
{

  if (text == 0) {return 0;} 

int text_size=strlen(text); 
int line_counter=0, lastpos=0;

  for (int pos=0;  pos < text_size;  pos++)
  {
    if ((text[pos] == '\n') && (pos-1 >= 0) && (text[pos-1] != '\n')) {lastpos=pos;  line_counter++;} 
  }

  if ((text_size > 0) && (text[text_size-1] != '\n')) {line_counter++;}
 
return line_counter;

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_not_empt_lines_in_chstr(char *text, int start_line_number)
{

  if (text == 0) {return 0;} 

int text_size=strlen(text); 
int notempt_line_counter=0, line_counter=0, lastpos=0; 

  for (int pos=0;  pos < text_size;  pos++)
  {
    if ((text[pos] == '\n') || (text[pos] == '\0'))
    {
    line_counter++;
      if ((pos-1 >= 0) && (text[pos-1] != '\n') && (line_counter >= start_line_number)) {lastpos=pos;  notempt_line_counter++;}
    } 
  }

  if ((text[text_size-1] != '\n') && (text[text_size-1] != '\0')) {line_counter++;}

return notempt_line_counter;

}



// ==========================================================================
// ==========================================================================



int calc_max_line_size_in_chstr(char *text)
{

  if (text == 0) {return 0;} 

int text_size=strlen(text); 
int line_start_pos=0, line_len=0, line_len_max=0;

  for (int pos=0;  pos < text_size;  pos++)
  {
    if ((text[pos] == '\n') || ((text[pos] == '\0'))) 
    {
    line_len=pos-line_start_pos;    
      if (line_len_max < line_len) {line_len_max=line_len;}    
    line_start_pos=pos+1;
    } 
  }

  if ((text_size > 0) && (text[text_size-1] != '\n') && (text[text_size-1] != '\0')) 
  {
  line_len=text_size-line_start_pos;    
    if (line_len_max < line_len) {line_len_max=line_len;}
  }

return line_len_max;

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_lines_in_file(char *filename)
{

  if (filename == 0) {return 0;}

char * text=0;
int text_size=write_file_to_chstr(filename, text);  //  cout<<"DEBUG: "<<" text_size="<<text_size<<endl;
int line_counter=0, lastpos=0;

  for (int pos=0;  pos < text_size;  pos++)
  {
    if (text[pos] == '\n') {lastpos=pos;  line_counter++;} 
  }

  if (text[text_size-1] != '\n') {line_counter++;}
  if (text) {delete[] text;  text=0;}

return line_counter;

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_not_empt_lines_in_file(char *filename)
{

  if (filename == 0) {return 0;}

char * text=0;
int text_size=write_file_to_chstr(filename, text);
int line_counter=0, lastpos=0;

  for (int pos=0;  pos < text_size;  pos++)
  {
    if ((text[pos] == '\n') && (pos-1 >= 0) && (text[pos-1] != '\n')) {lastpos=pos;  line_counter++;} 
  }

  if (text[text_size-1] != '\n') {line_counter++;}
  if (text) {delete[] text;  text=0;} 

return line_counter;

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_not_empt_lines_in_file(char *filename, int start_line_number)
{

  if (filename == 0) {return 0;}

char * text=0;
int text_size=write_file_to_chstr(filename, text);
int notempt_line_counter=calc_amount_of_not_empt_lines_in_chstr(text, start_line_number);

  if (text) {delete[] text;  text=0;} 

return notempt_line_counter;

}



// ==========================================================================
// ==========================================================================



int calc_max_line_size_in_file(char *filename) 
{

  if (filename == 0) {return 0;}

char * text=0;
int text_size=write_file_to_chstr(filename, text);  //  cout<<"DEBUG: "<<" text_size="<<text_size<<endl;
int line_start_pos=0, line_len=0, line_len_max=0;

  for (int pos=0;  pos < text_size;  pos++)
  {
    if (text[pos] == '\n') 
    {
    line_len=pos-line_start_pos;    
      if (line_len_max < line_len) {line_len_max=line_len;}    
    line_start_pos=pos+1;
    } 
  }

  if ((text_size > 0) && (text[text_size-1] != '\n')) 
  {
  line_len=text_size-line_start_pos;    
    if (line_len_max < line_len) {line_len_max=line_len;}
  }

  if (text) {delete[] text;  text=0;}

return line_len_max;

}



// ==========================================================================
// ==========================================================================



bool sr_line_position_and_len_in_chstr(char *chstr,  int line_number, int * line_start_pos, int * line_end_pos, int * line_len)
{

  if ((chstr == 0) || (line_number < 0)) {return false;}
  if (line_len) {*line_len=0;}


int line_counter=0;
int pos=0, line_start_pos_loc=0, line_end_pos_loc=0,  line_len_loc=0;


  while (chstr[pos] != '\0')
  {

    if ((chstr[pos] == '\n') || ((chstr[pos+1] == '\0') && (chstr[pos] != '\n')))
    {    

      if (line_counter == line_number)
      {
      line_end_pos_loc=pos-1;
        if (chstr[pos] != '\n') {line_end_pos_loc=pos;}
      line_len_loc=line_end_pos_loc-line_start_pos_loc+1;
      
      
        if (line_start_pos) 
        {*line_start_pos=line_start_pos_loc;}
        if (line_end_pos ) {*line_end_pos=line_end_pos_loc;}
        if (line_len     ) {*line_len=line_len_loc;}

      return true;
      }

    line_start_pos_loc=pos+1;  line_counter++;
    } 

  pos++;
  }

  if (line_start_pos) {*line_start_pos=-1;}
  if (line_end_pos ) {*line_end_pos=-2;}
  if (line_len     ) {*line_len=0;}

return false;

}



// ==========================================================================
// ==========================================================================



int sr_line_number_for_line_with_given_phrase_from_chstr(char *chstr, char *chstr_key_phrase, int start_pos)      // returns -1 if not found
{

  if (chstr == 0) {return -1;}
  if (chstr_key_phrase == 0) {return -1;}
  if (strlen(chstr_key_phrase) < 1) {return -1;}

int returning_start_pos=0;
bool ok=sr_given_phrase_in_chstr(chstr, chstr_key_phrase, start_pos, & returning_start_pos, 0);

  if (ok == true) 
  {
  int line_feed_smb_counter=0;

    for (int pos=0;  pos < returning_start_pos;  pos++)
    {
      if (chstr[pos] == '\n') {line_feed_smb_counter++;}
    }    

  return line_feed_smb_counter;
  }

return -1;

}



// ==========================================================================
// ==========================================================================



int  extract_line_with_given_number_from_chstr(char *chstr, int line_number, char * & chstr_extracting_line)
{

  if (chstr == 0) {return false;}
  if (line_number < 0) {return false;}
  if (chstr_extracting_line) {return false;}

int line_start_pos=0, line_end_pos=0, line_len=0;

bool ok=sr_line_position_and_len_in_chstr(chstr,  line_number, & line_start_pos, & line_end_pos, & line_len);

  if (ok == true) 
  {
  chstr_extracting_line=new char[line_len+1];  chstr_extracting_line[0]='\0';
  ok=copy_replace(chstr,  line_start_pos,  line_len,     chstr_extracting_line,  0,    true);
    if (ok == true) {return line_len;}
  }

return false;

}



// ==========================================================================
// ==========================================================================



int  extract_phrase_in_given_area_from_chstr(char *chstr,  int given_start_pos,  int given_end_pos,  char * & chstr_extracting_phrase, int *extr_phr_start_pos, int *extr_phr_end_pos)
{

  if (chstr == 0) {return -1;}
  if (chstr_extracting_phrase) {return -1;}
  if (given_start_pos < 0) {return -1;}


bool is_found=false;

int extr_phr_start_pos_loc=0, extr_phr_end_pos_loc=0;


extr_phr_start_pos_loc=sr_altern_1sp_2phrase__pos2phrase(chstr,  given_start_pos,  & is_found);

  if (is_found == true)
  {
  extr_phr_end_pos_loc=sr_altern_1sp_2phrase__pos2phrase(chstr,  extr_phr_start_pos_loc,  & is_found);
    if (is_found == true)
    {
    int len=extr_phr_end_pos_loc-extr_phr_start_pos_loc+1;
    chstr_extracting_phrase=new char[len+1];
    chstr_extracting_phrase[0]='\0';
    bool ok=copy_replace(chstr,  extr_phr_start_pos_loc,  len,    chstr_extracting_phrase,  0,    true);    // ! ! !
      if (ok == true) 
      {

        if (extr_phr_start_pos) {*extr_phr_start_pos=extr_phr_start_pos_loc;}
        if (extr_phr_end_pos  ) {*extr_phr_end_pos  =extr_phr_end_pos_loc;  }

      return len;
      }
    }
  }


return -1;

}



// ==========================================================================
// ==========================================================================



int  extract_phrase_after_key_phr_in_x_lines_from_chstr(char *chstr,  char *chstr_key_phrase,  char * & chstr_extracting_phrase,  int start_pos,  int d_x_lines,  int phrase_number_in_line)
{

  if (chstr == 0) {return -1;}
  if (chstr_key_phrase == 0) {return -1;}
  if (chstr_extracting_phrase) {return -1;}

  if (start_pos < 0) {start_pos=0;}
  if (d_x_lines < 0) {d_x_lines=0;}
  if (phrase_number_in_line < 0) {phrase_number_in_line=0;}


int key_phrase_start_pos=0, key_phrase_end_pos=0;
int line_start_pos=0, line_end_pos=0, line_len=0;                  
int extract_phrase_start_pos=0, extract_phrase_end_pos=0;
int line_number=0;
bool ok=false;
int returning_len=0;


line_number=sr_line_number_for_line_with_given_phrase_from_chstr(chstr, chstr_key_phrase, start_pos);
ok=sr_line_position_and_len_in_chstr(chstr,  line_number+d_x_lines,  & line_start_pos,  & line_end_pos, & line_len);

  if (ok == true)
  {
  int phrase_counter=0;
  bool is_last_smb_sp=true;    
  bool is_cur_smb_sp=true;    
  bool is_next_smb_sp=true;    

    for (int pos=line_start_pos;  pos <= line_end_pos;  pos++)
    {

      if ((chstr[pos] != ' ') && (chstr[pos] != '\t') && (chstr[pos] != '\n')) {is_cur_smb_sp=false;} else {is_cur_smb_sp=true;}
      if ((pos == line_end_pos) || (chstr[pos+1] == ' ') || (chstr[pos+1] == '\t') || (chstr[pos+1] == '\n')) {is_next_smb_sp=true;} else {is_next_smb_sp=false;}

      if ((is_last_smb_sp == true) && (is_cur_smb_sp == false)) {extract_phrase_start_pos=pos;}
      if ((is_cur_smb_sp == false) && (is_next_smb_sp == true)) {extract_phrase_end_pos=pos;}

      if ((is_cur_smb_sp == false) && (is_next_smb_sp == true))
      {
        if (phrase_counter == phrase_number_in_line)
        {
        int phr_len=extract_phrase_end_pos-extract_phrase_start_pos+1;
        chstr_extracting_phrase=new char[phr_len+1];
        chstr_extracting_phrase[0]='\0';
        ok=copy_replace(chstr,  extract_phrase_start_pos,  phr_len,    chstr_extracting_phrase,  0,    true);    // ! ! !
          if (ok == true) {return true;} else {return false;}
        }
      phrase_counter++;
      }

    is_last_smb_sp=is_cur_smb_sp;
    }      //  for (int pos=line_start_pos;  pos <= line_end_pos;  pos++)

  }      //  if (ok == true)


return false;

}



// ==========================================================================
// ==========================================================================



int extract_line_after_line_with_given_phrase_from_chstr(char *chstr, char *chstr_key_phrase, char *& extracting_phrase, int start_pos)
{

  if (chstr == 0) {return -1;}
  if (chstr_key_phrase == 0) {return -1;}
  if (strlen(chstr_key_phrase) < 1) {return -1;}
  if (extracting_phrase != 0) {return -1;}

  if (start_pos < 0) {start_pos=0;}


int key_phrase_start_pos=0, key_phrase_end_pos=0;
int extracting_phrase_start_pos=0, extracting_phrase_end_pos=0;


bool ok=sr_given_phrase_in_chstr(chstr,  chstr_key_phrase,  start_pos,  & key_phrase_start_pos,  & key_phrase_end_pos);    // sr key phrase

  if (ok == true)
  {
  bool is_found=false;
  extracting_phrase_start_pos=sr_altern_1sp_2phrase__pos2phrase(chstr,  key_phrase_end_pos+1,  & is_found);    // sr extracting value(phrase)
    if (is_found == true)
    {
    is_found=false;
    extracting_phrase_end_pos=sr_altern_1phrase_2sp__lastpos1phrase(chstr,  extracting_phrase_start_pos, & is_found);
      if (is_found == true)    //  value is found!  now copy it!
      {
      int extr_phrase_len=extracting_phrase_end_pos-extracting_phrase_start_pos+1;
        if (extr_phrase_len > 0)
        {
        extracting_phrase=new char[extr_phrase_len+1];
        extracting_phrase[0]='\0';
        ok=copy_replace(chstr,  extracting_phrase_start_pos, extr_phrase_len,   extracting_phrase,  0,   true);
        return true;
        }
      }     
    }
  }

return false;

}



// ==========================================================================
// ==========================================================================



int extract_lines_as_text_from_file(char *filename, char * & text, int & amount_of_lines,  int start_line_number, int line_wished_amount_be_extracted)
{

amount_of_lines=0;

  if (filename == 0) {return 0;}
  if (text != 0) {return 0;}
  if (line_wished_amount_be_extracted < 1) {return 0;}


int line_amount_will_be_extracted=0;

char *src_text=0;
int src_size=write_file_to_chstr(filename, src_text, start_line_number);      //==DEBUG==//  cout<<"src_size="<<src_size<<endl;      //==DEBUG==//
amount_of_lines=calc_amount_of_lines_in_chstr(src_text);
  if (amount_of_lines <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
  if (line_amount_will_be_extracted > amount_of_lines)  {line_amount_will_be_extracted=amount_of_lines;}
  if (line_amount_will_be_extracted <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}

int start_pos_front_line=0,  end_pos_back_line=0;
bool ok2=sr_line_position_and_len_in_chstr(src_text,  line_amount_will_be_extracted-1, 0, & end_pos_back_line, 0);
  if (ok2 == false) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
  if ((end_pos_back_line+1 < src_size) && (src_text[end_pos_back_line+1] == '\n')) {end_pos_back_line++;}

int dest_size=end_pos_back_line+1;
text=new char[dest_size+1];
copy_replace(src_text, start_pos_front_line, dest_size, text,  0, true);


  if (src_text) {delete[] src_text;  src_text=0;}

return dest_size;

}



// ==========================================================================
// ==========================================================================



int extract_lines_from_file(char *filename, char ** & chstr_lines, int & amount_of_lines,  int & max_line_size,  int start_line_number, int line_wished_amount_be_extracted)
{

amount_of_lines=0;  max_line_size=0;

  if (filename == 0) {return 0;}
  if (chstr_lines != 0) {return 0;}
  if (line_wished_amount_be_extracted < 1) {return 0;}


int line_amount_will_be_extracted=0;

char *src_text=0;
int src_size=write_file_to_chstr(filename, src_text, start_line_number);      //==DEBUG==//  cout<<"src_size="<<src_size<<endl;      //==DEBUG==//
amount_of_lines=calc_amount_of_lines_in_chstr(src_text);
  if (amount_of_lines <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
  if (line_amount_will_be_extracted > amount_of_lines)  {line_amount_will_be_extracted=amount_of_lines;}
  if (line_amount_will_be_extracted <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
max_line_size=calc_max_line_size_in_chstr(src_text);    //==DEBUG==//  cout<<"line_size="<<line_size<<"  amount_of_lines="<<amount_of_lines<<endl;  //==DEBUG==//

chstr_lines=new char *[line_amount_will_be_extracted];

  for (int line_num=0;  line_num < line_amount_will_be_extracted;  line_num++)
  {chstr_lines[line_num]=0;  chstr_lines[line_num]=new char[max_line_size+1];  chstr_lines[line_num][0]='\0';}


int line_start_pos=0,  line_end_pos=0, cur_line_len=0;
bool ok=false;

  for (int line_num=0;  line_num < line_amount_will_be_extracted;  line_num++)
  {
  ok=sr_line_position_and_len_in_chstr(src_text,  line_num, &line_start_pos,  &line_end_pos, &cur_line_len);
    if (ok == true) {copy_replace(src_text, line_start_pos, cur_line_len, chstr_lines[line_num], 0, true);}
  }

  if (src_text) {delete[] src_text;  src_text=0;}

return line_amount_will_be_extracted;

}


// ==========================================================================
// ==========================================================================



int extract_lines_from_file(char *filename, char * & chstr_lines, int & amount_of_lines,  int & line_size,  int start_line_number)
{

amount_of_lines=0;  line_size=0;

  if (filename == 0) {return 0;}
  if (chstr_lines != 0) {return 0;}


char *src_text=0;
int src_size=write_file_to_chstr(filename, src_text);      //==DEBUG==//  cout<<"src_size="<<src_size<<endl;      //==DEBUG==//
amount_of_lines=calc_amount_of_lines_in_chstr(src_text)-start_line_number;
  if (amount_of_lines <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
line_size=calc_max_line_size_in_chstr(src_text);    //==DEBUG==//  cout<<"line_size="<<line_size<<"  amount_of_lines="<<amount_of_lines<<endl;  //==DEBUG==//

int dest_size=amount_of_lines*(line_size+1);    //==DEBUG==//      cout<<"dest_size="<<dest_size<<endl;  
chstr_lines=new char[dest_size];

  for (int pos=0;  pos < dest_size;  pos++)
  {chstr_lines[pos]='\0';}

  for (int src_pos=0, dest_pos=0, line_number=0;  src_pos < src_size;  src_pos++)
  {
    if ((src_text[src_pos] == '\n') || (src_text[src_pos] == '\0')) 
    {  if ((src_text[src_pos] == '\n') || ((src_text[src_pos] == '\0') && (src_pos > 0) && (src_text[src_pos-1] != '\n'))) {line_number++;  dest_pos=0;}  }
    else
    {  
      if (start_line_number <= line_number) {chstr_lines[(line_number-start_line_number)*(line_size+1)+dest_pos]=src_text[src_pos];  dest_pos++;}
    }
  }

  if (src_text) {delete[] src_text;  src_text=0;}

return amount_of_lines;

}



// ==========================================================================
// ==========================================================================



int extract_not_empty_lines_from_file(char *filename,  char * chstr_lines,  int & amount_of_not_empt_lines,  int & line_size)
{

  if (filename == 0) {return 0;}
  if (chstr_lines != 0) {return 0;}


char *src_text=0;

int src_size=write_file_to_chstr(filename, src_text);
amount_of_not_empt_lines=calc_amount_of_not_empt_lines_in_chstr(src_text);
line_size=calc_max_line_size_in_chstr(src_text); 

int dest_size=amount_of_not_empt_lines*(line_size+1);  
chstr_lines=new char[dest_size];

  for (int pos=0;  pos < dest_size;  pos++)
  {chstr_lines[pos]='\0';}

  for (int src_pos=0, dest_pos=0, line_number=0;  src_pos < src_size;  src_pos++)
  {
    if ((src_text[src_pos] == '\n') || (src_text[src_pos] == '\0')) 
    {  if ((src_text[src_pos] == '\n') || ((src_text[src_pos] == '\0') && (src_pos > 0) && (src_text[src_pos-1] != '\n'))) {line_number++;  dest_pos=0;}  }
    else
    {chstr_lines[line_number*(line_size+1)+dest_pos]=src_text[src_pos];  dest_pos++;}
  }


  if (src_text) {delete[] src_text;  src_text=0;}

return amount_of_not_empt_lines;

}



// ==========================================================================
// ==========================================================================



int extract_lines_from_file(char *filename, char ** & chstr_lines, int & amount_of_lines,  int & line_size,  int start_line_number)
{

  if (filename == 0) {return 0;}
  if (chstr_lines != 0) {return 0;}

line_size=0; 

char *src_text=0;
int src_size=write_file_to_chstr(filename, src_text);
amount_of_lines=calc_amount_of_lines_in_chstr(src_text)-start_line_number;
  if (amount_of_lines <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
line_size=calc_max_line_size_in_chstr(src_text);

chstr_lines=new char *[amount_of_lines];
  for (int line_num=0;  line_num < amount_of_lines;  line_num++)
  {chstr_lines[line_num]=new char[line_size+1];  chstr_lines[line_num][0]='\0';}

  for (int src_pos=0, dest_pos=0, line_number=0;  src_pos < src_size;  src_pos++)
  {
    if ((src_text[src_pos] == '\n') || (src_text[src_pos] == '\0')) 
    {  if ((src_text[src_pos] == '\n') || ((src_text[src_pos] == '\0') && (src_pos > 0) && (src_text[src_pos-1] != '\n'))) {line_number++;  dest_pos=0;}  }
    else
    {  
      if (start_line_number <= line_number) 
      {chstr_lines[line_number-start_line_number][dest_pos]=src_text[src_pos];  chstr_lines[line_number-start_line_number][dest_pos+1]='\0';  dest_pos++;}  
    }
  }

  if (src_text) {delete[] src_text;  src_text=0;}

return amount_of_lines;

}



// ==========================================================================
// ==========================================================================




int extract_lines_from_file(char *filename, char ** & chstr_lines, int & amount_of_lines,  int & line_size,  char *chstr_forbid_1st_line_symb_list,  int start_line_number)
{

  if ((filename == 0) || (chstr_forbid_1st_line_symb_list == 0)) {return 0;}
  if (chstr_lines != 0) {return 0;}

line_size=0; 
int record_line_number=0,  amount_forbid_symb=strlen(chstr_forbid_1st_line_symb_list);
bool line_allowed=true;


char *src_text=0;
int src_size=write_file_to_chstr(filename, src_text, start_line_number);
amount_of_lines=calc_amount_of_lines_in_chstr(src_text,  chstr_forbid_1st_line_symb_list);
  if (amount_of_lines <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
line_size=calc_max_line_size_in_chstr(src_text);

chstr_lines=new char *[amount_of_lines];
  for (int line_num=0;  line_num < amount_of_lines;  line_num++)
  {chstr_lines[line_num]=new char[line_size+1];  chstr_lines[line_num][0]='\0';}


  for (int i=0;  i < amount_forbid_symb;  i++)  {  if (src_text[0] == chstr_forbid_1st_line_symb_list[i]) {line_allowed=false;  break;}  }

  for (int src_pos=0, dest_pos=0, line_number=0;  src_pos < src_size;  src_pos++)
  {
    if ((src_text[src_pos] == '\n') || (src_text[src_pos] == '\0')) 
    {
      if ((src_text[src_pos] == '\n') || ((src_text[src_pos] == '\0') && (src_pos > 0) && (src_text[src_pos-1] != '\n'))) 
      {
        if (line_allowed == true) {record_line_number++;}
        if (src_pos+1 < src_size) {line_allowed=true;      for (int i=0;  i < amount_forbid_symb;  i++)  {  if (src_text[src_pos+1] == chstr_forbid_1st_line_symb_list[i]) {line_allowed=false;  break;}  }      }
      line_number++;  dest_pos=0;       
      }
    }
    else
    {  
      if (line_allowed == true) 
      {chstr_lines[record_line_number][dest_pos]=src_text[src_pos];  chstr_lines[record_line_number][dest_pos+1]='\0';  dest_pos++;}  
    }
  }

  if (src_text) {delete[] src_text;  src_text=0;}

return amount_of_lines;

}




// ==========================================================================
// ==========================================================================



int  calc_phrases_after_key_phr_in_x_lines_from_chstr(char *chstr,  char *chstr_key_phrase,  int start_pos,  int d_x_lines)
{

  if (chstr == 0) {return 0;}
  if (chstr_key_phrase == 0) {return 0;}

  if (start_pos < 0) {start_pos=0;}
  if (d_x_lines < 0) {d_x_lines=0;}


int key_phrase_start_pos=0, key_phrase_end_pos=0;
int line_start_pos=0, line_end_pos=0, line_len=0;                  
int phrase_start_pos=0, phrase_end_pos=0;
int line_number=0;
bool ok=false;
int returning_len=0;
int phrase_counter=0;


line_number=sr_line_number_for_line_with_given_phrase_from_chstr(chstr, chstr_key_phrase, start_pos);
ok=sr_line_position_and_len_in_chstr(chstr,  line_number+d_x_lines,  & line_start_pos,  & line_end_pos, & line_len);

  if (ok == true)
  {
  bool is_last_smb_sp=true;    

    for (int pos=line_start_pos;  pos <= line_end_pos;  pos++)
    {

      if ((chstr[pos] != ' ') && (chstr[pos] != '\t'))
      {
        if (is_last_smb_sp == true) {phrase_counter++;  phrase_start_pos=pos;  is_last_smb_sp=false;}  
      }
      else
      {
        if ((is_last_smb_sp == false) || (pos == line_end_pos)) 
        {phrase_end_pos=pos;  is_last_smb_sp=true;}
      }

    }      //  for (int pos=line_start_pos;  pos <= line_end_pos;  pos++)

  }      //  if (ok == true)


return phrase_counter;

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_phrases_in_chstr(char * text)
{ 

  if (text == 0) {return 0;}

bool is_found=true;
int phrase_counter=0;
int pos=0;
int whole_text_len=strlen(text);

  while ((is_found == true) && (pos < whole_text_len))
  {
  pos=sr_altern_1phrase_2sp__lastpos1phrase(text,  pos, & is_found);
    if (is_found == true) {phrase_counter++;}
  pos++;
  }

return phrase_counter;

}




// ==========================================================================
// ==========================================================================



int calc_max_len_of_phrase_in_chstr(char * text)
{

  if (text == 0) {return 0;}

bool is_found=true;
int pos_start=0, pos_end=0;
int whole_text_len=strlen(text);
int maxlen=0; 

  while ((is_found == true) && (pos_start < whole_text_len))
  {
  pos_start=sr_altern_1sp_2phrase__pos2phrase(text,  pos_start, & is_found);
    if (is_found == true) 
    {
    pos_end=sr_altern_1phrase_2sp__lastpos1phrase(text,  pos_start, & is_found);
      if ((is_found == true) && (pos_end-pos_start+1 > maxlen)) {maxlen=pos_end-pos_start+1;} 
    }
  pos_start=pos_end+1;
  }       

return maxlen;

}



// ==========================================================================
// ==========================================================================



int calc_max_len_of_phrases_in_phrase_ar(char ** chstr_phrases, int phrases_count)
{

  if (chstr_phrases == 0) {return 0;}

int max_phrase_len=0,  len=0; 

  for (int i=0;  i < phrases_count;  i++)
  {
  len=strlen(chstr_phrases[i]);
 
   if (len > max_phrase_len) {max_phrase_len=len;} 

  }

return max_phrase_len;

}



// ==========================================================================
// ==========================================================================



int calc_max_len_of_phrases_in_chstr_by_col(char ** chstr_phrases, int phrases_count, int col_count, int * & max_phrase_len)
{

  if ((chstr_phrases == 0) || (col_count <= 0)|| (max_phrase_len)) {return 0;}

max_phrase_len=new int[col_count];  

  for (int i=0;  i < col_count;  i++)
  {
  max_phrase_len[i]=0;
  int len=0;

    for (int counter=i;  counter < phrases_count;  counter+=col_count)
    {
    len=strlen(chstr_phrases[counter]);
      if (len > max_phrase_len[i]) {max_phrase_len[i]=len;} 
    }

  }


int sum_max_row_sz=0;

  for (int i=0;  i < col_count;  i++)
  {sum_max_row_sz+=max_phrase_len[i];}  

return sum_max_row_sz;

}



// ==========================================================================
// ==========================================================================



bool sr_if_there_signs_bef_phrases_in_chstr_by_col(char ** chstr_phrases, int phrases_count, int col_count, bool * & sign_is)
{

  if ((chstr_phrases == 0) || (col_count <= 0) || (sign_is)) {return 0;}

sign_is=new bool[col_count];  
bool any_sign_is=false;

  for (int i=0;  i < col_count;  i++)
  {
  sign_is[i]=false;

    for (int counter=i;  counter < phrases_count;  counter+=col_count)
    {
      if ((chstr_phrases[counter][0] == '+') || (chstr_phrases[counter][0] == '-')) {sign_is[i]=true;  any_sign_is=true;  break;}
    }

  }

return any_sign_is;

}



// ==========================================================================
// ==========================================================================



int extract_phrases_from_file(char *filename, char ** & chstr_phrase, int & amount_of_phrases,  int & phrase_len,  int start_line_number)
{

  if (filename == 0) {return 0;}
  if (chstr_phrase != 0) {return 0;}

phrase_len=0; 

char *src_text=0;
int src_size=write_file_to_chstr(filename, src_text, start_line_number);
amount_of_phrases=calc_amount_of_phrases_in_chstr(src_text);
  if (amount_of_phrases <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
phrase_len=calc_max_len_of_phrase_in_chstr(src_text);

chstr_phrase=new char *[amount_of_phrases];
  for (int phrase_num=0;  phrase_num < amount_of_phrases;  phrase_num++)
  {chstr_phrase[phrase_num]=new char[phrase_len+1];  chstr_phrase[phrase_num][0]='\0';}


int pos_start=0, pos_end=0, len=0;
bool is_found=true; 
int phrase_counter=0;
  
  while (is_found == true)
  {
  pos_start=sr_altern_1sp_2phrase__pos2phrase(src_text,  pos_end, & is_found);
    if (is_found == true) 
    {
    pos_end=sr_altern_1phrase_2sp__lastpos1phrase(src_text,  pos_start, & is_found);
      if (is_found == true) {len=pos_end-pos_start+1;  copy_replace(src_text,  pos_start,  len,  chstr_phrase[phrase_counter],  0,  true);  phrase_counter++;} 
    }
  }

  if (src_text) {delete[] src_text;  src_text=0;}

return phrase_counter;

}



// ==========================================================================
// ==========================================================================



int calc_amount_of_anytypenums_in_chstr(const char *const text, char *symbs_consid_as_num)
{

  if (text == 0) {return 0;}
int text_len=strlen(text);

bool is_found=true;
int num_counter=0;
int pos=0; 

  while ((text_len > pos) && (is_found == true))
  {
  pos=sr_altern_1num_2sp__lastpos1num(text,  pos, & is_found, symbs_consid_as_num);
    if (is_found == true) {num_counter++;  pos++;  }
  }

return num_counter;

}



// ==========================================================================
// ==========================================================================



int calc_max_len_of_anytypenums_in_chstr(char *text, char *symbs_consid_as_num)
{

  if (text == 0) {return 0;}

bool is_found=true;
int pos_start=0, pos_end=0;
int whole_text_len=strlen(text);
int maxlen=0;

  while ((is_found == true) && (pos_end < whole_text_len))
  {
  pos_start=sr_altern_1sp_2num__pos2num(text, pos_end, & is_found, symbs_consid_as_num);
    if (is_found == true) 
    {
    pos_end=sr_altern_1num_2sp__lastpos1num(text, pos_start, & is_found, symbs_consid_as_num);
      if ((is_found == true) && (pos_end-pos_start+1 > maxlen)) {maxlen=pos_end-pos_start+1;  pos_end++;} 
    }
  }

return maxlen;

}


// ==========================================================================
// ==========================================================================



int extract_anytypenums_from_file(char *filename, char ** & chstr_nums, int & amount_of_nums,  int & num_len,  int start_line_number, char *symbs_consid_as_num)
{

  if (filename == 0) {return 0;}
  if (chstr_nums != 0) {return 0;}


num_len=0; 

char *src_text=0;
int src_size=write_file_to_chstr(filename, src_text, start_line_number);
amount_of_nums=calc_amount_of_anytypenums_in_chstr(src_text);  
  if (amount_of_nums <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;}
num_len=calc_max_len_of_anytypenums_in_chstr(src_text)+64;

chstr_nums=new char *[amount_of_nums];
  for (int number_of_num=0;  number_of_num < amount_of_nums;  number_of_num++)
  {chstr_nums[number_of_num]=new char[num_len+1];  chstr_nums[number_of_num][0]='\0';}


int pos_start=0, pos_end=0, len=0;
bool is_found=true;
int num_counter=0;
  
  while ((is_found == true) && (pos_start < src_size))
  {
  pos_start=sr_altern_1sp_2num__pos2num(src_text,  pos_end, & is_found, symbs_consid_as_num); 
    if (is_found == true) 
    {
    pos_end=sr_altern_1num_2sp__lastpos1num(src_text,  pos_start, & is_found, symbs_consid_as_num);
      if (is_found == true) {len=pos_end-pos_start+1; copy_replace(src_text,  pos_start,  len,  chstr_nums[num_counter],  0,  true);  pos_end++; num_counter++;} 
    }
  pos_start=pos_end+1;
  }

  if (src_text) {delete[] src_text;  src_text=0;}

return num_counter;


}



// ==========================================================================
// ==========================================================================



void display_chstr_lines(char * chstr_lines,  int amount_of_lines,  int & line_size)
{

  if (chstr_lines == 0) {cout<<"=====  @@ empty chstr var @@ ===="<<endl;}
  else 
  {
  //-------------------------------------------------------------
    for (int line_num=0;  line_num < amount_of_lines;  line_num++)
    {

      for (int pos=0;  pos < line_size;  pos++) 
      { 
        if (chstr_lines[(line_size+1)*line_num+pos] == '\0') {break;}  
      cout<<chstr_lines[(line_size+1)*line_num+pos];
      }

    cout<<endl;
    }  
  //-------------------------------------------------------------
  }

}



// ==========================================================================
// ==========================================================================



void display_chstr_lines(char ** chstr_lines,  int amount_of_lines)
{

  if (chstr_lines == 0) {cout<<"=====  @@ empty chstr var @@ ===="<<endl;}
  else 
  {
  //-------------------------------------------------------------
    for (int line_num=0;  line_num < amount_of_lines;  line_num++)
    {
    int pos=0;

      while (chstr_lines[line_num][pos] != '\0')
      {cout<<chstr_lines[line_num][pos];  pos++;}

    cout<<endl;
    }  
  //-------------------------------------------------------------
  }

}



// ==========================================================================
// ==========================================================================



void display_chstr_phrases_by_col(char ** chstr_var,  int phrases_count, int col_count, char *chstr_separ)
{

  if ((chstr_var) && (phrases_count > 0) && (col_count > 0))
  {
  int counter=0;

    for (int counter=0;  counter < phrases_count;  counter++)
    {
    cout<<chstr_var[counter];

      if ((counter % col_count) < col_count-1) {  if (chstr_separ) {cout<<chstr_separ;}  } else {cout<<endl;}
 
    }

  }

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_file_by_col(char * filename, char ** chstr_var,  int phrases_count, int col_count, char *chstr_separ, bool B1New0Add)
{

  if ((filename == 0) || (chstr_var == 0) || (phrases_count < 0) || (col_count < 1)) {return 0;}


string string_var;

int counter=0;

  for (int counter=0;  counter < phrases_count;  counter++)
  {
  string_var.append(chstr_var[counter]);

    if (counter % col_count < col_count-1) {  if (chstr_separ) {string_var.append(chstr_separ);}  } else {string_var.append((char *) "\n\0");} 
  }


int file_sz=write_chstr_to_file(string_var.c_str(), filename, B1New0Add);

string_var.clear();

return file_sz;

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_chstr_by_align_col_1(char ** chstr_phrases, char * & chstr_dest,  int phrases_count, int col_count, int hor_interval_len, int ver_interval_len)
{

  if ((chstr_phrases == 0) || (phrases_count < 0) || (col_count < 1) || (chstr_dest != 0) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


int max_phrase_len=calc_max_len_of_phrase_in_chstr(chstr_phrases, phrases_count);
int chstr_sz=phrases_count*(max_phrase_len+hor_interval_len+ver_interval_len+1)+1;
  if (chstr_sz <= 0) {return 0;}
char *chstr=0;
chstr=new char[chstr_sz+1];
chstr[0]='\0';
char interval_symb=' ';


int counter=0, len=0,  rest_fill_sp_len=0, cur_pos=0;

  for (int counter=0;  counter < phrases_count;  counter++)
  {
  len=strlen(chstr_phrases[counter]);   
  rest_fill_sp_len=max_phrase_len-len;
  strcat(chstr, chstr_phrases[counter]);
  cur_pos=strlen(chstr);
  cur_pos=chstr_fill(chstr, interval_symb, cur_pos, rest_fill_sp_len);

    if (counter % col_count < col_count-1) 
    {cur_pos=chstr_fill(chstr, interval_symb, cur_pos, col_count);} 
    else 
    {cur_pos=chstr_fill(chstr, (char) '\n', cur_pos, ver_interval_len+1);}
 
  }


int result_sz=strlen(chstr);
chstr_dest=new char[result_sz+1];
chstr_dest[0]='\0';
strcpy(chstr_dest, chstr);

  if (chstr) {delete[] chstr;  chstr=0;}

return result_sz+1;

}



// ==========================================================================
// ==========================================================================


int write_chstr_phrases_to_file_by_align_col_1(char * filename,  char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, bool B1New0Add)
{

  if ((filename == 0) || (chstr_phrases == 0) || (phrases_count < 0) || (col_count < 1) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


char *chstr_dest=0;
int file_sz=0;
file_sz=write_chstr_phrases_to_chstr_by_align_col_1(chstr_phrases, chstr_dest, phrases_count, col_count, hor_interval_len, ver_interval_len);
file_sz=write_chstr_to_file(chstr_dest, filename, B1New0Add);

  if (chstr_dest) {delete[] chstr_dest;  chstr_dest=0;}

return file_sz;

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_screen_by_align_col_1(char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len)
{

  if ((chstr_phrases == 0) || (phrases_count <= 0) || (col_count < 1) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


char *chstr_dest=0;
int chstr_dest_sz=0;
chstr_dest_sz=write_chstr_phrases_to_chstr_by_align_col_1(chstr_phrases, chstr_dest, phrases_count, col_count, hor_interval_len, ver_interval_len);
cout<<chstr_dest<<endl;

  if (chstr_dest) {delete[] chstr_dest;  chstr_dest=0;}

return chstr_dest_sz;

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_chstr_by_align_col_2(char ** chstr_phrases, char * & chstr_dest,  int phrases_count, int col_count, int hor_interval_len, int ver_interval_len)
{

  if ((chstr_phrases == 0) || (phrases_count < 0) || (col_count < 1) || (chstr_dest != 0) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}

string string_var;
int *max_phrase_len=0; 
calc_max_len_of_phrases_in_chstr_by_col(chstr_phrases, phrases_count, col_count, max_phrase_len);   
  if (max_phrase_len == 0) {  if (max_phrase_len) {delete[] max_phrase_len;  max_phrase_len=0;} return 0;} 

char chstrinterval_str[]=" ";
int counter=0, len=0,  rest_fill_sp_len=0, col_number=0;

  for (int counter=0;  counter < phrases_count;  counter++)
  {
  col_number=counter % col_count;  len=strlen(chstr_phrases[counter]); 
  rest_fill_sp_len=max_phrase_len[col_number]-len;
  string_var.append(chstr_phrases[counter]); 
  string_append_x_times(& string_var, chstrinterval_str, rest_fill_sp_len);

    if (col_number < col_count-1) {string_append_x_times(& string_var, chstrinterval_str, hor_interval_len);} else {string_append_x_times(& string_var, (char *) "\n", ver_interval_len+1);}
  }


int result_len=string_var.length();
chstr_dest=new char[result_len+1];   chstr_dest[0]='\0';
strcpy(chstr_dest, string_var.c_str());
string_var.clear();


  if (max_phrase_len) {delete[] max_phrase_len;  max_phrase_len=0;}

return result_len+1;

}



// ==========================================================================
// ==========================================================================


int write_chstr_phrases_to_file_by_align_col_2(char * filename,  char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, bool B1New0Add)
{

  if ((filename == 0) || (chstr_phrases == 0) || (phrases_count < 0) || (col_count < 1) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


char *chstr_dest=0;
int file_sz=0;
file_sz=write_chstr_phrases_to_chstr_by_align_col_2(chstr_phrases, chstr_dest, phrases_count, col_count, hor_interval_len, ver_interval_len);
file_sz=write_chstr_to_file(chstr_dest, filename, B1New0Add);

  if (chstr_dest) {delete[] chstr_dest;  chstr_dest=0;}

return file_sz;

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_screen_by_align_col_2(char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len)
{

  if ((chstr_phrases == 0) || (phrases_count <= 0) || (col_count < 1) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


char *chstr_dest=0;
int chstr_dest_sz=0;
chstr_dest_sz=write_chstr_phrases_to_chstr_by_align_col_2(chstr_phrases, chstr_dest, phrases_count, col_count, hor_interval_len, ver_interval_len);
cout<<chstr_dest<<endl;

  if (chstr_dest) {delete[] chstr_dest;  chstr_dest=0;}

return chstr_dest_sz;

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_chstr_by_align_col_3(char ** chstr_phrases, char * & chstr_dest,  int phrases_count, int col_count, int hor_interval_len, int ver_interval_len)
{ 

  if ((chstr_phrases == 0) || (phrases_count < 0) || (col_count < 1) || (chstr_dest != 0) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


string string_var;
int *max_phrase_len=0; 
calc_max_len_of_phrases_in_chstr_by_col(chstr_phrases, phrases_count, col_count, max_phrase_len);   
  if (max_phrase_len == 0) {  if (max_phrase_len) {delete[] max_phrase_len;  max_phrase_len=0;} return 0;} 
bool *sign_is=0;
sr_if_there_signs_bef_phrases_in_chstr_by_col(chstr_phrases, phrases_count, col_count, sign_is); 
  if (sign_is == 0) {  if (sign_is) {delete[] sign_is;  sign_is=0;} return 0;}
  for (int col=0;  col < col_count;  col++) {  if (sign_is[col] == true) {max_phrase_len[col]++;}  }


char chstrinterval_str[]=" ";
int counter=0, len=0,  rest_fill_sp_len=0, col_number=0;

  for (int counter=0;  counter < phrases_count;  counter++)
  {
  col_number=counter % col_count;  len=strlen(chstr_phrases[counter]); 
    if ((sign_is[col_number] == true) && (chstr_phrases[counter][0] != '+') && (chstr_phrases[counter][0] != '-')) {string_var.append((char *) " ");  len++;}      // sign in col alignment
  rest_fill_sp_len=max_phrase_len[col_number]-len;
  string_var.append(chstr_phrases[counter]); 
  string_append_x_times(& string_var, chstrinterval_str, rest_fill_sp_len);

    if (col_number < col_count-1) {string_append_x_times(& string_var, chstrinterval_str, hor_interval_len);} else {string_append_x_times(& string_var, (char *) "\n", ver_interval_len+1);}
  }


int result_len=string_var.length();
chstr_dest=new char[result_len+1];   chstr_dest[0]='\0';
strcpy(chstr_dest, string_var.c_str());
string_var.clear();


  if (max_phrase_len) {delete[] max_phrase_len;  max_phrase_len=0;}
  if (sign_is) {delete[] sign_is;  sign_is=0;}

return result_len+1;
  
}



// ==========================================================================
// ==========================================================================


int write_chstr_phrases_to_file_by_align_col_3(char * filename,  char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, bool B1New0Add)
{ 

  if ((filename == 0) || (chstr_phrases == 0) || (phrases_count < 0) || (col_count < 1) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


char *chstr_dest=0;
int chstr_dest_sz=0;
chstr_dest_sz=write_chstr_phrases_to_chstr_by_align_col_3(chstr_phrases, chstr_dest, phrases_count, col_count, hor_interval_len, ver_interval_len);
chstr_dest_sz=write_chstr_to_file(chstr_dest, filename, B1New0Add);

  if (chstr_dest) {delete[] chstr_dest;  chstr_dest=0;}

return chstr_dest_sz;

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_screen_by_align_col_3(char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len)
{

  if ((chstr_phrases == 0) || (phrases_count <= 0) || (col_count < 1) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


char *chstr_dest=0;
int chstr_dest_sz=0;
chstr_dest_sz=write_chstr_phrases_to_chstr_by_align_col_3(chstr_phrases, chstr_dest, phrases_count, col_count, hor_interval_len, ver_interval_len);
cout<<chstr_dest<<endl;

  if (chstr_dest) {delete[] chstr_dest;  chstr_dest=0;}

return chstr_dest_sz;

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_chstr_by_align_col_4(char ** chstr_phrases, char * & chstr_dest,  int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, char type)
{

  if ((chstr_phrases == 0) || (phrases_count < 0) || (col_count < 1) || (chstr_dest != 0) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}
  if (type < 0) {type=0;}     if (type > 1) {type=1;}

string string_var;
int *max_phrase_len=0; 
int max_row_len=calc_max_len_of_phrases_in_chstr_by_col(chstr_phrases, phrases_count, col_count, max_phrase_len);   
  if (max_phrase_len == 0) {  if (max_phrase_len) {delete[] max_phrase_len;  max_phrase_len=0;} return 0;} 
bool *sign_is=0;
sr_if_there_signs_bef_phrases_in_chstr_by_col(chstr_phrases, phrases_count, col_count, sign_is); 
  if (sign_is == 0) {  if (sign_is) {delete[] sign_is;  sign_is=0;} return 0;}
  for (int col=0;  col < col_count;  col++) {  if (sign_is[col] == true) {max_phrase_len[col]++;  max_row_len++;}  }
char chstrinterval_str[]=" ",  frame_hor_str[]="*",  frame_ver_str[]="*";
  if (type == 1) {frame_hor_str[0]='-';   frame_ver_str[0]='|';}
int max_row_len_full=max_row_len+hor_interval_len*col_count*2+(col_count+1)*strlen(frame_ver_str);


int counter=0, len=0,  rest_fill_sp_len=0, col_number=0;


string_append_x_times(& string_var, frame_hor_str, max_row_len_full);  string_var.append((char *) "\n");     // hor line                  caption


bool mean_line_new=true;


  for (int counter=0;  counter < phrases_count;  counter++)      //  main
  {

    if (mean_line_new == true)
    {

      for (int i=0;  i < ver_interval_len;  i++)      //  cell up or down part
      {
      string_var.append(frame_ver_str);  
  
        for (int j=0;  j < col_count;  j++)
        {
        string_append_x_times(& string_var, chstrinterval_str, hor_interval_len);  string_append_x_times(& string_var, chstrinterval_str, max_phrase_len[j]);
        string_append_x_times(& string_var, chstrinterval_str, hor_interval_len);  string_var.append(frame_ver_str);  
        }
      string_var.append((char *) "\n");
      } 

    string_var.append(frame_ver_str);
    }  
    
  string_append_x_times(& string_var, chstrinterval_str, hor_interval_len);
  col_number=counter % col_count;  len=strlen(chstr_phrases[counter]); 
    if ((sign_is[col_number] == true) && (chstr_phrases[counter][0] != '+') && (chstr_phrases[counter][0] != '-')) {string_var.append((char *) " ");  len++;}      // sign in col alignment
  rest_fill_sp_len=max_phrase_len[col_number]-len;
  string_var.append(chstr_phrases[counter]); 
  string_append_x_times(& string_var, chstrinterval_str, rest_fill_sp_len);
  string_append_x_times(& string_var, chstrinterval_str, hor_interval_len);
  string_var.append(frame_ver_str);

    if (col_number >= col_count-1) 
    {
    string_append_x_times(& string_var, (char *) "\n", ver_interval_len);

      for (int i=0;  i < ver_interval_len;  i++)      //  cell up or down part
      {
      string_var.append(frame_ver_str);  
  
        for (int j=0;  j < col_count;  j++)
        {
        string_append_x_times(& string_var, chstrinterval_str, hor_interval_len);  string_append_x_times(& string_var, chstrinterval_str, max_phrase_len[j]);
        string_append_x_times(& string_var, chstrinterval_str, hor_interval_len);  string_var.append(frame_ver_str);  
        }
      string_var.append((char *) "\n");
      }  

    string_append_x_times(& string_var, frame_hor_str, max_row_len_full);  string_var.append((char *) "\n");     // hor line

    mean_line_new=true;
    }
    else
    {mean_line_new=false;}
  }


int result_len=string_var.length();
chstr_dest=new char[result_len+1];   chstr_dest[0]='\0';
strcpy(chstr_dest, string_var.c_str());
string_var.clear();


  if (max_phrase_len) {delete[] max_phrase_len;  max_phrase_len=0;}
  if (sign_is) {delete[] sign_is;  sign_is=0;}

return result_len+1;
  
}



// ==========================================================================
// ==========================================================================


int write_chstr_phrases_to_file_by_align_col_4(char * filename,  char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, bool B1New0Add, char type)
{ 

  if ((filename == 0) || (chstr_phrases == 0) || (phrases_count < 0) || (col_count < 1) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}

char *chstr_dest=0;
int chstr_dest_sz=0;
chstr_dest_sz=write_chstr_phrases_to_chstr_by_align_col_4(chstr_phrases, chstr_dest, phrases_count, col_count, hor_interval_len, ver_interval_len, type);
chstr_dest_sz=write_chstr_to_file(chstr_dest, filename, B1New0Add);

  if (chstr_dest) {delete[] chstr_dest;  chstr_dest=0;}

return chstr_dest_sz;

}



// ==========================================================================
// ==========================================================================



int write_chstr_phrases_to_screen_by_align_col_4(char ** chstr_phrases, int phrases_count, int col_count, int hor_interval_len, int ver_interval_len, char type)
{

  if ((chstr_phrases == 0) || (phrases_count <= 0) || (col_count < 1) || (hor_interval_len < 0) || (ver_interval_len < 0)) {return 0;}


char *chstr_dest=0;
int chstr_dest_sz=0;
chstr_dest_sz=write_chstr_phrases_to_chstr_by_align_col_4(chstr_phrases, chstr_dest, phrases_count, col_count, hor_interval_len, ver_interval_len, type);
cout<<chstr_dest<<endl;

  if (chstr_dest) {delete[] chstr_dest;  chstr_dest=0;}

return chstr_dest_sz;

}



// ==========================================================================
// ==========================================================================



int calc_phrases_amount_in_line(char *chstr_text, int line_number)
{

  if ((chstr_text == 0) || (line_number < 0)) {return 0;} 

int text_sz=strlen(chstr_text); 

  if (text_sz < 1) {return 0;}


bool is_match_found=true;
int counter=0;
int line_start_pos=0,  line_end_pos=0,  line_len=0;  
bool ok=sr_line_position_and_len_in_chstr(chstr_text,  line_number, & line_start_pos, & line_end_pos, & line_len);
int pos=line_start_pos; 
  if (pos > 0) {pos--;} 

  if (ok == true)
  { 
    while ((pos <= line_end_pos) && (is_match_found == true))
    {  
    pos=sr_altern_1sp_2phrase__pos2phrase(chstr_text,  pos,  & is_match_found); 
      if ((is_match_found == true) && (pos <= line_end_pos)) {counter++;} else { return counter;}
    pos++;
    }
  }        

return counter;
	
}



// ==========================================================================
// ==========================================================================



int extract_phrases_from_chstr(char *src_text, char ** & chstr_phrases, int & amount_of_extracted_phrases,  int & phrase_len,  int start_phrase_num,  int col_count)
{     

  if (src_text == 0) {return 0;}
  if (chstr_phrases != 0) {return 0;}
  if (col_count == 0) {return 0;}   

phrase_len=0;

int src_size=strlen(src_text); 
int amount_of_phrases=calc_amount_of_phrases_in_chstr(src_text);     
  if (amount_of_phrases <= 0) {  if (src_text) {delete[] src_text;  src_text=0;}  return 0;} 

  if (col_count < 0)     //  if (col_count < 0)   means we don't know num of columns, so we must find it  
  {
  col_count=calc_phrases_amount_in_line(src_text, 0);
    if (col_count < 1) {  if (src_text) {delete[] src_text;  src_text=0;}   return 0;}
  }    

phrase_len=calc_max_len_of_phrase_in_chstr(src_text);

int amount_of_must_be_extracted_phrases=(amount_of_phrases-start_phrase_num)/col_count;  
  if ((amount_of_phrases-start_phrase_num)-(amount_of_must_be_extracted_phrases)*col_count > 0) {amount_of_must_be_extracted_phrases++;} 


chstr_phrases=new char *[amount_of_must_be_extracted_phrases];
  for (int phrase_num=0;  phrase_num < amount_of_must_be_extracted_phrases;  phrase_num++)
  {chstr_phrases[phrase_num]=new char[phrase_len+1];  chstr_phrases[phrase_num][0]='\0';}


int pos_start=0, pos_end=0, len=0;
bool is_found=true; 
int phrase_counter=0;
  
  while ((is_found == true) && (amount_of_extracted_phrases < amount_of_must_be_extracted_phrases) && (pos_start <src_size))
  { 
  pos_start=sr_altern_1sp_2phrase__pos2phrase(src_text,  pos_start, & is_found);
    if (is_found == true) 
    { 
    pos_end=sr_altern_1phrase_2sp__lastpos1phrase(src_text,  pos_start, & is_found);
      if (is_found == true) 
      { 
        if ((phrase_counter == start_phrase_num) || ((phrase_counter-start_phrase_num) % col_count == 0)) 
        {len=pos_end-pos_start+1;    copy_replace(src_text,  pos_start,  len,  chstr_phrases[amount_of_extracted_phrases],  0,  true);    amount_of_extracted_phrases++;}
      phrase_counter++; 
      } 
    }  
  pos_start=pos_end+1;
  } 


//  if (amount_of_extracted_phrases < amount_of_must_be_extracted_phrases) {cout<<"problem in extract_phrases_from_chstr(...)!"<<endl;}

return amount_of_must_be_extracted_phrases;

}


// ==========================================================================
// ==========================================================================



int extract_phrases_from_chstr(char *src_text, char * &  chstr_phrases, int * & len_array, int * & start_pos_array, int * max_phrase_len,  int start_phrase_num,  int col_count, int src_start_pos, int src_end_pos)
{
//  proc returns array size

  if (src_text == 0) {return 0;}
  if ((chstr_phrases) || (start_pos_array) || (len_array)) {return 0;}
  if ((src_end_pos >= 0) && (src_end_pos >= strlen(src_text)))  {return 0;}
  if (src_end_pos < 0) {src_end_pos=strlen(src_text)-1;}
  if (src_end_pos < src_start_pos) {return 0;}
  if (col_count < 1) {return 0;}
  if (start_phrase_num+1 > col_count) {return 0;}



//  * * * * * * * * * * * * *
//  define amount of phrases
//
int amount_of_found_phrases=0;
int amount_of_found_phrases_for_giv_column=0;

  {
  bool is_last_sp=true;
    for (int pos=src_start_pos;  pos <= src_end_pos;  pos++)
    {
      if ((src_text[pos] == ' ') || (src_text[pos] == '\t') || (src_text[pos] == '\n')) {is_last_sp=true;} else {  if (is_last_sp == true) {amount_of_found_phrases++;}  is_last_sp=false;}
    }

    if (amount_of_found_phrases < 1) {return 0;}

  amount_of_found_phrases_for_giv_column=amount_of_found_phrases/col_count+(amount_of_found_phrases-(amount_of_found_phrases/col_count)*col_count)/(start_phrase_num+1);
  }


//  * * * * * * * * * * * * *
//  define max len of phrases
//
  if (max_phrase_len)
  {
  bool is_last_sp=true;
  int cur_len=0;
  int phr_counter=0;

    for (int pos=src_start_pos;  pos <= src_end_pos;  pos++)
    {
      if ((src_text[pos] == ' ') || (src_text[pos] == '\t') || (src_text[pos] == '\n')) 
      {  
        if ((is_last_sp == false) && ((*max_phrase_len) < cur_len)) {    if ((phr_counter % col_count) == start_phrase_num) {*max_phrase_len=cur_len;}    phr_counter++;}
      cur_len=0; is_last_sp=true;
      } 
      else {cur_len++;  is_last_sp=false;}
    }  

    if ((is_last_sp == false) && ((*max_phrase_len) < cur_len)) {    if ((phr_counter % col_count) == start_phrase_num) {*max_phrase_len=cur_len;}    phr_counter++;} 
  }


//  * * * * * * * * * * * * *
//  define total sum of phrase sizes 
//
int total_phr_size_sum=0;
  {
  bool is_last_sp=true;
  int cur_len=0;
  int phr_counter=0;

    for (int pos=src_start_pos;  pos <= src_end_pos;  pos++)
    {
      if ((src_text[pos] == ' ') || (src_text[pos] == '\t') || (src_text[pos] == '\n')) 
      {  
        if (is_last_sp == false) {    if ((phr_counter % col_count) == start_phrase_num) {total_phr_size_sum+=cur_len;}    phr_counter++;}
      cur_len=0; is_last_sp=true; 
      } 
      else {cur_len++;  is_last_sp=false;}
    }  

    if (is_last_sp == false) {    if ((phr_counter % col_count) == start_phrase_num) {total_phr_size_sum+=cur_len;}    phr_counter++;} 

    if (total_phr_size_sum < 1) {return 0;}
  }


//  * * * * * * * * * * * * *
//  extracting phrases  MAIN part 
//
  { 
  int phr_counter=0;
  chstr_phrases=new char[total_phr_size_sum+1];  chstr_phrases[0]='\0';  //chstr_phrases[total_phr_size_sum]='\0';
  len_array=new int[amount_of_found_phrases_for_giv_column];                            for (int i=0;  i < amount_of_found_phrases_for_giv_column;  i++) {len_array[i]=0;}
  start_pos_array=new int[amount_of_found_phrases_for_giv_column];                      for (int i=0;  i < amount_of_found_phrases_for_giv_column;  i++) {start_pos_array[i]=0;}
  
  bool is_last_sp=true, is_cur_sp=true, is_next_sp=true;
  int phr_start_pos=0,  phr_end_pos=0,  phr_len=0;
  int amount_of_written_phrases=0;      //  will be compared to amount_of_found_phrases_for_giv_column
  int cur_dest_pos=0;
  
    for (int pos=src_start_pos;  pos <= src_end_pos;  pos++)
    {
      if ((src_text[pos] == ' ') || (src_text[pos] == '\t') || (src_text[pos] == '\n')) {is_cur_sp=true;} else {is_cur_sp=false;}
      if (pos+1 <= src_end_pos) {  if ((src_text[pos+1] == ' ') || (src_text[pos+1] == '\t') || (src_text[pos+1] == '\n')) {is_next_sp=true;} else {is_next_sp=false;}  }
      if (pos == src_end_pos) {is_next_sp=true;}
      
      //  analize
      if ((is_last_sp == true) && (is_cur_sp == false)) {phr_start_pos=pos;}
      if ((is_cur_sp == false) && (is_next_sp == true)) {phr_end_pos=pos;  phr_len=phr_end_pos-phr_start_pos+1;}

     // writing
      if ((is_cur_sp == false) && (is_next_sp == true) && (amount_of_written_phrases < amount_of_found_phrases_for_giv_column))
      {  
        if ((phr_counter % col_count) == start_phrase_num)
        {
        //  phrase is localised ! ! !
          if (cur_dest_pos < total_phr_size_sum) {start_pos_array[amount_of_written_phrases]=cur_dest_pos;  len_array[amount_of_written_phrases]=phr_len;}      //  put start posand len into array

          for (int src_rec_pos=phr_start_pos;  src_rec_pos <= phr_end_pos;  src_rec_pos++)
          {  
            if (cur_dest_pos < total_phr_size_sum)      //    ! ! ! ensurence for mem leak and main part ! ! ! 
            {chstr_phrases[cur_dest_pos]=src_text[src_rec_pos];  cur_dest_pos++;} 
          }   
          
        amount_of_written_phrases++;
        }   
      phr_counter++;
      }

    is_last_sp=is_cur_sp;
    }      //  for (int pos=src_start_pos;  pos <= src_end_pos;  pos++)  


    for (int i=cur_dest_pos;  i <= total_phr_size_sum;  i++)      //  for safety     line end smb putting
    {chstr_phrases[i]='\0';}                                      //  for safety     

  }


return amount_of_found_phrases_for_giv_column;

}


// ==========================================================================
// ==========================================================================



bool extract_column_from_file(char * filename, char ** & col_phrases, int & amount_in_col,  int & max_phrase_len, int col_start_number, int col_count, char *forbid_1st_symbs,  int start_line_number)
{
 
  if ((filename == 0) || (col_phrases != 0)) {return false;}


max_phrase_len=0;  amount_in_col=0;

char *extracted_file_text=0;
int extracted_file_text_sz=0;

extracted_file_text_sz=write_file_to_chstr(filename, extracted_file_text,  forbid_1st_symbs,  start_line_number,  true);    //  file text extracting ! ! !
  if (extracted_file_text_sz <= 0) {  if (extracted_file_text) {delete[] extracted_file_text;  extracted_file_text=0;}  return 0;}

int amount_in_col_loc=0, max_phrase_len_loc=0;   
extract_phrases_from_chstr(extracted_file_text, col_phrases, amount_in_col_loc, max_phrase_len_loc,  col_start_number,  col_count);  
amount_in_col=amount_in_col_loc;  max_phrase_len=max_phrase_len_loc;


  if (extracted_file_text) {delete[] extracted_file_text;  extracted_file_text=0;} 
  if (amount_in_col < 1) {return false;} 

return true;

}



// ==========================================================================
// ==========================================================================




bool write_columns_to_screen_3(int amount_in_col, int dy, int dx, bool B1New0Add, char ** col1_phrases, char ** col2_phrases, char ** col3_phrases, char ** col4_phrases, char ** col5_phrases)
{

  if ((dy < 0) || (dx < 0)) {return false;}
  if ((col1_phrases == 0) && (col2_phrases == 0) && (col3_phrases == 0) && (col4_phrases == 0) && (col5_phrases == 0)) {return false;}


int len=0, max_len=0;

  for (int line_num=0;  line_num < amount_in_col;  line_num++)
  {
    if ((col1_phrases) && (col1_phrases[line_num])) {len=strlen(col1_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col2_phrases) && (col2_phrases[line_num])) {len=strlen(col2_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col3_phrases) && (col3_phrases[line_num])) {len=strlen(col3_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col4_phrases) && (col4_phrases[line_num])) {len=strlen(col4_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col5_phrases) && (col5_phrases[line_num])) {len=strlen(col5_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
  }

int amount_of_not_null_cols=0;

  if (col1_phrases) {amount_of_not_null_cols++;}
  if (col2_phrases) {amount_of_not_null_cols++;}
  if (col3_phrases) {amount_of_not_null_cols++;}
  if (col4_phrases) {amount_of_not_null_cols++;}
  if (col5_phrases) {amount_of_not_null_cols++;}


int amount_of_elements=amount_in_col*amount_of_not_null_cols;  

char **chstr_intermed=0;
chstr_intermed=new char *[amount_of_elements];

  for (int i=0;  i < amount_of_elements;  i++)
  {
  chstr_intermed[i]=new char[max_len+1];
  chstr_intermed[i][0]='\0';
  }


int counter=0;

  for (int line_num=0;  line_num < amount_in_col;  line_num++)   
  {
    if ((col1_phrases) && (col1_phrases[line_num])) {chstr_intermed[counter][0]='\0';  strcpy(chstr_intermed[counter], col1_phrases[line_num]);  counter++;}
    if ((col2_phrases) && (col2_phrases[line_num])) {chstr_intermed[counter][0]='\0';  strcpy(chstr_intermed[counter], col2_phrases[line_num]);  counter++;}
    if ((col3_phrases) && (col3_phrases[line_num])) {chstr_intermed[counter][0]='\0';  strcpy(chstr_intermed[counter], col3_phrases[line_num]);  counter++;}
    if ((col4_phrases) && (col4_phrases[line_num])) {chstr_intermed[counter][0]='\0';  strcpy(chstr_intermed[counter], col4_phrases[line_num]);  counter++;}
    if ((col5_phrases) && (col5_phrases[line_num])) {chstr_intermed[counter][0]='\0';  strcpy(chstr_intermed[counter], col5_phrases[line_num]);  counter++;}
  }

write_chstr_phrases_to_screen_by_align_col_3(chstr_intermed, amount_of_elements, amount_of_not_null_cols, dx, dy); 

free_chstr_array(chstr_intermed,  amount_of_elements);

return true;

}



// ==========================================================================
// ==========================================================================



bool write_columns_to_screen_4(int amount_in_col, int dy, int dx, bool B1New0Add, char ** col1_phrases, char ** col2_phrases, char ** col3_phrases, char ** col4_phrases, char ** col5_phrases)
{


  if ((dy < 0) || (dx < 0)) {return false;}
  if ((col1_phrases == 0) && (col2_phrases == 0) && (col3_phrases == 0) && (col4_phrases == 0) && (col5_phrases == 0)) {return false;}


int len=0, max_len=0;

  for (int line_num=0;  line_num < amount_in_col;  line_num++)
  {
    if ((col1_phrases) && (col1_phrases[line_num])) {len=strlen(col1_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col2_phrases) && (col2_phrases[line_num])) {len=strlen(col2_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col3_phrases) && (col3_phrases[line_num])) {len=strlen(col3_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col4_phrases) && (col4_phrases[line_num])) {len=strlen(col4_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col5_phrases) && (col5_phrases[line_num])) {len=strlen(col5_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
  }

int amount_of_not_null_cols=0;

  if (col1_phrases) {amount_of_not_null_cols++;}
  if (col2_phrases) {amount_of_not_null_cols++;}
  if (col3_phrases) {amount_of_not_null_cols++;}
  if (col4_phrases) {amount_of_not_null_cols++;}
  if (col5_phrases) {amount_of_not_null_cols++;}


int amount_of_elements=amount_in_col*amount_of_not_null_cols;

char **chstr_intermed=0;
chstr_intermed=new char *[amount_of_elements];

  for (int i=0;  i < amount_of_elements;  i++)
  {
  chstr_intermed[i]=new char[max_len+1];
  chstr_intermed[i][0]='\0';
  }


int counter=0;

  for (int line_num=0;  line_num < amount_in_col;  line_num++)   
  {
    if ((col1_phrases) && (col1_phrases[line_num])) {strcpy(chstr_intermed[counter], col1_phrases[line_num]);  counter++;}
    if ((col2_phrases) && (col2_phrases[line_num])) {strcpy(chstr_intermed[counter], col2_phrases[line_num]);  counter++;}
    if ((col3_phrases) && (col3_phrases[line_num])) {strcpy(chstr_intermed[counter], col3_phrases[line_num]);  counter++;}
    if ((col4_phrases) && (col4_phrases[line_num])) {strcpy(chstr_intermed[counter], col4_phrases[line_num]);  counter++;}
    if ((col5_phrases) && (col5_phrases[line_num])) {strcpy(chstr_intermed[counter], col5_phrases[line_num]);  counter++;}
  }

write_chstr_phrases_to_screen_by_align_col_4(chstr_intermed, amount_of_elements, amount_of_not_null_cols, dx, dy); 

free_chstr_array(chstr_intermed,  amount_of_elements);

return true;

}



// ==========================================================================
// ==========================================================================



int write_columns_to_file_3(char * filename, int amount_in_col, int dy, int dx, bool B1New0Add, char ** col1_phrases, char ** col2_phrases, char ** col3_phrases, char ** col4_phrases, char ** col5_phrases)
{

  if ((filename == 0) || (dy < 0) || (dx < 0)) {return 0;}
  if ((col1_phrases == 0) && (col2_phrases == 0) && (col3_phrases == 0) && (col4_phrases == 0) && (col5_phrases == 0)) {return false;}


int len=0, max_len=0;

  for (int line_num=0;  line_num < amount_in_col;  line_num++)
  {
    if ((col1_phrases) && (col1_phrases[line_num])) {len=strlen(col1_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col2_phrases) && (col2_phrases[line_num])) {len=strlen(col2_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col3_phrases) && (col3_phrases[line_num])) {len=strlen(col3_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col4_phrases) && (col4_phrases[line_num])) {len=strlen(col4_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col5_phrases) && (col5_phrases[line_num])) {len=strlen(col5_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
  }

int amount_of_not_null_cols=0;

  if (col1_phrases) {amount_of_not_null_cols++;}
  if (col2_phrases) {amount_of_not_null_cols++;}
  if (col3_phrases) {amount_of_not_null_cols++;}
  if (col4_phrases) {amount_of_not_null_cols++;}
  if (col5_phrases) {amount_of_not_null_cols++;}


int amount_of_elements=amount_in_col*amount_of_not_null_cols;

char **chstr_intermed=0;
chstr_intermed=new char *[amount_of_elements];

  for (int i=0;  i < amount_of_elements;  i++)
  {
  chstr_intermed=new char *[max_len+1];
  chstr_intermed[i][0]='\0';
  }


int counter=0;

  for (int line_num=0;  line_num < amount_in_col;  line_num++)   
  {
    if ((col1_phrases) && (col1_phrases[line_num])) {strcpy(chstr_intermed[counter], col1_phrases[line_num]);  counter++;}
    if ((col2_phrases) && (col2_phrases[line_num])) {strcpy(chstr_intermed[counter], col2_phrases[line_num]);  counter++;}
    if ((col3_phrases) && (col3_phrases[line_num])) {strcpy(chstr_intermed[counter], col3_phrases[line_num]);  counter++;}
    if ((col4_phrases) && (col4_phrases[line_num])) {strcpy(chstr_intermed[counter], col4_phrases[line_num]);  counter++;}
    if ((col5_phrases) && (col5_phrases[line_num])) {strcpy(chstr_intermed[counter], col5_phrases[line_num]);  counter++;}
  }

int file_size=write_chstr_phrases_to_file_by_align_col_3(filename,  chstr_intermed, amount_of_elements, amount_of_not_null_cols, dx, dy, B1New0Add); 


free_chstr_array(chstr_intermed,  amount_of_elements);

return file_size;

}



// ==========================================================================
// ==========================================================================



int write_columns_to_file_4(char * filename, int amount_in_col, int dy, int dx, bool B1New0Add, char ** col1_phrases, char ** col2_phrases, char ** col3_phrases, char ** col4_phrases, char ** col5_phrases)
{

  if ((filename == 0) || (dy < 0) || (dx < 0)) {return 0;}
  if ((col1_phrases == 0) && (col2_phrases == 0) && (col3_phrases == 0) && (col4_phrases == 0) && (col5_phrases == 0)) {return false;}


int len=0, max_len=0;

  for (int line_num=0;  line_num < amount_in_col;  line_num++)
  {
    if ((col1_phrases) && (col1_phrases[line_num])) {len=strlen(col1_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col2_phrases) && (col2_phrases[line_num])) {len=strlen(col2_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col3_phrases) && (col3_phrases[line_num])) {len=strlen(col3_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col4_phrases) && (col4_phrases[line_num])) {len=strlen(col4_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
    if ((col5_phrases) && (col5_phrases[line_num])) {len=strlen(col5_phrases[line_num]);  if (len > max_len) {max_len=len;}  }
  }

int amount_of_not_null_cols=0;

  if (col1_phrases) {amount_of_not_null_cols++;}
  if (col2_phrases) {amount_of_not_null_cols++;}
  if (col3_phrases) {amount_of_not_null_cols++;}
  if (col4_phrases) {amount_of_not_null_cols++;}
  if (col5_phrases) {amount_of_not_null_cols++;}


int amount_of_elements=amount_in_col*amount_of_not_null_cols;

char **chstr_intermed=0;
chstr_intermed=new char *[amount_of_elements];

  for (int i=0;  i < amount_of_elements;  i++)
  {
  chstr_intermed=new char *[max_len+1];
  chstr_intermed[i][0]='\0';
  }


int counter=0;

  for (int line_num=0;  line_num < amount_in_col;  line_num++)
  {
    if ((col1_phrases) && (col1_phrases[line_num])) {strcpy(chstr_intermed[counter], col1_phrases[line_num]);  counter++;}
    if ((col2_phrases) && (col2_phrases[line_num])) {strcpy(chstr_intermed[counter], col2_phrases[line_num]);  counter++;}
    if ((col3_phrases) && (col3_phrases[line_num])) {strcpy(chstr_intermed[counter], col3_phrases[line_num]);  counter++;}
    if ((col4_phrases) && (col4_phrases[line_num])) {strcpy(chstr_intermed[counter], col4_phrases[line_num]);  counter++;}
    if ((col5_phrases) && (col5_phrases[line_num])) {strcpy(chstr_intermed[counter], col5_phrases[line_num]);  counter++;}
  }

int file_size=write_chstr_phrases_to_file_by_align_col_4(filename,  chstr_intermed, amount_of_elements, amount_of_not_null_cols, dx, dy, B1New0Add, 0); 


free_chstr_array(chstr_intermed,  amount_of_elements);

return file_size;

}



// ==========================================================================
// ==========================================================================



template<typename T> bool extract_num_array_from_chstr(char * chstr_text,  T *& num_array, int & amount_of_num,  bool type_dbl1_int0)
{

  if (num_array) {return false;}


char * chstr_phrases=0;
int  * len_array=0; 
int  * start_pos_array=0;

int max_phrase_len=0;


amount_of_num=extract_phrases_from_chstr(chstr_text, chstr_phrases, len_array, start_pos_array,  & max_phrase_len,  0,  1,  0,  -1);      //  !  !  !


  if (amount_of_num > 0) 
  {
  char *chstr_temp=0;
  chstr_temp=new char[max_phrase_len+32];
  chstr_temp[0]='\0';
  
  
  num_array=new T[amount_of_num]; 


    if (type_dbl1_int0 == false)
    {       
      for (int i=0;  i < amount_of_num; i++)
      {
        for (int k=0;  k < len_array[i];  k++)  {chstr_temp[k]=chstr_phrases[start_pos_array[i]+k];}    chstr_temp[len_array[i]]='\0';  
      num_array[i]=atoi(chstr_temp);         
      }
    }
    else
    {       
      for (int i=0;  i < amount_of_num; i++)
      {
        for (int k=0;  k < len_array[i];  k++)  {chstr_temp[k]=chstr_phrases[start_pos_array[i]+k];}     chstr_temp[len_array[i]]='\0';  
      num_array[i]=atof(chstr_temp);     
      }
    }
  
    if (chstr_temp) {delete[]  chstr_temp;  chstr_temp=0;}
  
  }


  if (start_pos_array) {delete[] start_pos_array;  start_pos_array=0;}
  if (len_array      ) {delete[] len_array;        len_array=0;}
  if (chstr_phrases  ) {delete[] chstr_phrases;    chstr_phrases=0;}
  
  
  if (amount_of_num < 1) {return false;}
  
return true;    
  

}



// ==========================================================================
// ==========================================================================



bool extract_bool_array_from_chstr(char * chstr_text,  bool *& num_array, int & amount_of_num)
{

  if (num_array) {return false;}


char * chstr_phrases=0;
int  * len_array=0; 
int  * start_pos_array=0;

int max_phrase_len=0;


amount_of_num=extract_phrases_from_chstr(chstr_text, chstr_phrases, len_array, start_pos_array,  & max_phrase_len,  0,  1,  0,  -1);      //  !  !  !


  if (amount_of_num > 0) 
  {
  char *chstr_temp=0;
  chstr_temp=new char[max_phrase_len+32];
  chstr_temp[0]='\0';
  
  num_array=new bool[amount_of_num]; 
           
    for (int i=0;  i < amount_of_num; i++)
    {
      for (int k=0;  k < len_array[i];  k++)  {chstr_temp[k]=chstr_phrases[start_pos_array[i]+k];}    chstr_temp[len_array[i]]='\0';  
      if (atoi(chstr_temp) == 1) {num_array[i]=true;} else {num_array[i]=false;}         
    }
    
  
    if (chstr_temp) {delete[]  chstr_temp;  chstr_temp=0;}
  
  }


  if (start_pos_array) {delete[] start_pos_array;  start_pos_array=0;}
  if (len_array      ) {delete[] len_array;        len_array=0;}
  if (chstr_phrases  ) {delete[] chstr_phrases;    chstr_phrases=0;}
  
  
  if (amount_of_num < 1) {return false;}
  
return true;    
  

}



// ==========================================================================
// ==========================================================================




template<typename T> bool extract_dbl_column_from_file(char * filename,  T *& col_phrases,  int & amount_in_col,  int col_number,  int col_count,  char *forbid_1st_symbs,  int start_line_num)
{

  if ((filename == 0) || (col_phrases) || (col_count == 0) || ((col_count >= 0) && (col_number+1 > col_count))) {return false;}
 

char **chstr_col_phrases=0;
int max_phrase_len=0;
bool ok=extract_column_from_file(filename,  chstr_col_phrases,  amount_in_col,  max_phrase_len,  col_number,  col_count,  forbid_1st_symbs,  start_line_num);

  if ((ok == false) || (amount_in_col < 1)) {free_chstr_array(chstr_col_phrases, amount_in_col);  return false;}

col_phrases=new T[amount_in_col]; 

  for (int i=0;  i < amount_in_col; i++)
  {col_phrases[i]=atof(chstr_col_phrases[i]);}


free_chstr_array(chstr_col_phrases, amount_in_col);

return true;

}


// ==========================================================================
// ==========================================================================




template<typename T> bool extract_int_column_from_file(char * filename,  T *& col_phrases,  int & amount_in_col,  int col_number,  int col_count,  char *forbid_1st_symbs,  int start_line_num)
{

  if ((filename == 0) || (col_phrases) || (col_count == 0) || ((col_count >= 0) && (col_number+1 > col_count))) {return false;}
 

char **chstr_col_phrases=0;
int max_phrase_len=0;
bool ok=extract_column_from_file(filename,  chstr_col_phrases,  amount_in_col,  max_phrase_len,  col_number,  col_count,  forbid_1st_symbs,  start_line_num);

  if ((ok == false) || (amount_in_col < 1)) {free_chstr_array(chstr_col_phrases, amount_in_col);  return false;}

col_phrases=new T[amount_in_col]; 

  for (int i=0;  i < amount_in_col; i++)
  {col_phrases[i]=atoi(chstr_col_phrases[i]);}


free_chstr_array(chstr_col_phrases, amount_in_col);

return true;

}



// ==========================================================================
// ==========================================================================



template< typename T > bool write_dbl_array_to_chstr(T * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len, int prec, int wished_len)
{

  if ((col_phrases == 0) || (chstr_phrases) || (amount_of_elem < 1)) {return false;}

char chstr_dbl[64];  chstr_dbl[0]='\0';


//  define max len size
//
max_elem_len=0;
int len=0;

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_dbl[0]='\0';
    if (wished_len < 1) {sprintf(chstr_dbl, "%.*E", prec, col_phrases[i]);} else {sprintf(chstr_dbl, "% *.*E", wished_len, prec, col_phrases[i]);}
  //  if (wished_len < 1) {sprintf(chstr_dbl, "%.*E", prec, col_phrases[i]);} else {sprintf(chstr_dbl, "% .*E", wished_len, prec, col_phrases[i]);}
  len=strlen(chstr_dbl);
    if (max_elem_len < len) {max_elem_len=len;}
  }


//  memory allocation
//
chstr_phrases=new char*[amount_of_elem];

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_phrases[i]=new char[max_elem_len+1];
  chstr_phrases[i][0]='\0';
  }


//  arary writing to chstr
//
  for (int i=0;  i < amount_of_elem;  i++)
  {
    if (wished_len < 1) {sprintf(chstr_phrases[i], "%.*E", prec, col_phrases[i]);} else {sprintf(chstr_phrases[i], "% *.*E", wished_len, prec, col_phrases[i]);}
  //  if (wished_len < 1) {sprintf(chstr_phrases[i], "%.*E", prec, col_phrases[i]);} else {sprintf(chstr_phrases[i], "% .*E", wished_len, prec, col_phrases[i]);}
  }


return true;

}



// ==========================================================================
// ==========================================================================




template< typename T > bool write_int_array_to_chstr(T * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len, int wished_len)
{


  if ((col_phrases == 0) || (chstr_phrases) || (amount_of_elem < 1)) {return false;}

char chstr_int[128];  chstr_int[0]='\0';


//  define max len size
//
max_elem_len=0;
int len=0;

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_int[0]='\0';
    if (wished_len < 1) 
    {sprintf(chstr_int, "%lli", (long long int) col_phrases[i]);} 
    else  
    {sprintf(chstr_int, "% *lli", wished_len, (long long int) col_phrases[i]);}
  len=strlen(chstr_int);
    if (max_elem_len < len) {max_elem_len=len;}
  }


//  memory allocation
//
chstr_phrases=new char*[amount_of_elem];

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_phrases[i]=new char[max_elem_len+1];
  chstr_phrases[i][0]='\0';
  }


//  arary writing to chstr
//
  for (int i=0;  i < amount_of_elem;  i++)
  {
    if (wished_len < 1) 
    {sprintf(chstr_phrases[i], "%lli", (long long int) col_phrases[i]);} 
    else 
    {sprintf(chstr_phrases[i], "% *lli", wished_len, (long long int) col_phrases[i]);}
  }


return true;

}



// ==========================================================================
// ==========================================================================



void write_chstr_phrases_to_screen(char ** col_phrases, int amount_of_elem, int start_elem_num, bool hor1ver0, bool do_enum)
{

  if ((col_phrases) && (amount_of_elem > 0))
  {

    for (int i=0;  i < amount_of_elem;  i++)
    {
      if (hor1ver0 == true)
      {
        if (do_enum == true) {cout<<(start_elem_num+i)<<") ";}
      cout<<col_phrases[start_elem_num+i];
        if (i < amount_of_elem-1) {cout<<"   ";} else {cout<<endl;} 
      }
      else
      {
        if (do_enum == true) {cout<<(start_elem_num+i)<<")\t";}
      cout<<col_phrases[start_elem_num+i];
        if (i < amount_of_elem-1) {cout<<"\n";} else {cout<<endl;}
      }  
    }

  }
  else
  {
  cout<<"==== ! ! PROBLEM in write_chstr_phrases_to_screen:   empty array ! !   ----------------------------------"<<endl;
  }

}



// ==========================================================================
// ==========================================================================



template< typename T > void write_dbl_array_to_screen(T * col_phrases, int amount_of_elem, int start_elem_num, bool hor1ver0, bool do_enum)
{

  if ((col_phrases) && (amount_of_elem > 0))
  {

    for (int i=0;  i < amount_of_elem;  i++)
    {
      if (hor1ver0 == true)
      {
        if (do_enum == true) {cout<<(start_elem_num+i)<<") ";}
      cout<<col_phrases[start_elem_num+i];
        if (i < amount_of_elem-1) {cout<<"   ";} else {cout<<endl;} 
      }
      else
      {
        if (do_enum == true) {cout<<(start_elem_num+i)<<")\t";}
      cout<<col_phrases[start_elem_num+i];
        if (i < amount_of_elem-1) {cout<<"\n";} else {cout<<endl;}
      }  
    }

  }
  else
  {
  cout<<"==== ! ! PROBLEM in write_dbl_array_to_screen:   empty array ! !   ----------------------------------"<<endl;
  }

}




// ==========================================================================
// ==========================================================================




template< typename T > void write_int_array_to_screen(T * col_phrases, int amount_of_elem, int start_elem_num, bool hor1ver0, bool do_enum)
{

  if ((col_phrases) && (amount_of_elem > 0))
  {

    for (int i=0;  i < amount_of_elem;  i++)
    {
      if (hor1ver0 == true)
      {
        if (do_enum == true) {cout<<(start_elem_num+i)<<") ";}
      cout<<col_phrases[start_elem_num+i];
        if (i < amount_of_elem-1) {cout<<"   ";} else {cout<<endl;} 
      }
      else
      {
        if (do_enum == true) {cout<<(start_elem_num+i)<<")\t";}
      cout<<col_phrases[start_elem_num+i];
        if (i < amount_of_elem-1) {cout<<"\n";} else {cout<<endl;}
      }  
    }

  }
  else
  {
  cout<<"==== ! ! PROBLEM in write_int_array_to_screen:   empty array ! !   ----------------------------------"<<endl;
  }

}


// ==========================================================================
// ==========================================================================



int add_chstr_array_column_to_file(char ** col_phrases, char * filename, int amount_in_col,  int start_line_num, char *forbid_1st_symbs, int col_num, char type)
{

  if ((col_phrases == 0) || (filename == 0) || (amount_in_col < 1) || (start_line_num < 0)) {return 0;}
  if (type < 0) {type=0;}
  if (type > 1) {type=1;}


char *src_pre_text_caption=0;
char *src_text=0;
bool not_take_empty_lines=true;


//  1) extract src file text
//
int src_pre_text_caption_sz=0;
int temp1_not_using=0;
  if (start_line_num > 0) {extract_lines_as_text_from_file(filename, src_pre_text_caption, temp1_not_using, 0, start_line_num);}
int extracted_src_file_sz=write_file_to_chstr(filename,  src_text,  forbid_1st_symbs,  start_line_num,  not_take_empty_lines);

  if (src_text == 0)  // || (extracted_src_file_sz < 1)) 
  {
    if (src_text) {delete[] src_text;  src_text=0;}  
    if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  return 0;
  }  


//  2) extracted src phrases
//
int col_count=calc_phrases_amount_in_line(src_text, start_line_num);
char **chstr_extracted_phrases=0;
int amount_of_extracted_phrases=0;
int max_extracted_phrase_len=0;                                                                                       //    start phrase num v     v col_count   
int src_phrase_amount=0;
  if (extracted_src_file_sz > 1) {src_phrase_amount=extract_phrases_from_chstr(src_text, chstr_extracted_phrases, amount_of_extracted_phrases,  max_extracted_phrase_len,  0,    1);} 
int dest_phrase_amount=dest_phrase_amount=amount_of_extracted_phrases+amount_in_col;

  if ((col_count > 0) && (amount_of_extracted_phrases/col_count != amount_in_col)) 
  {
  free_chstr_array(chstr_extracted_phrases, amount_of_extracted_phrases);
    if (src_text) {delete[] src_text;  src_text=0;}  
    if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  return 0;
  }  


//  3) define new max len
//
int added_col_phrase_len_max=calc_max_len_of_phrases_in_phrase_ar(col_phrases, amount_in_col);      // proc returns sum of max len of phrases in row
int phrase_len_max=0;

  if (added_col_phrase_len_max < max_extracted_phrase_len) {phrase_len_max=max_extracted_phrase_len;} else {phrase_len_max=added_col_phrase_len_max;}


//  4) memory allocation for new phrase array
//
char **chstr_result_phrases=0;
chstr_result_phrases=new char *[dest_phrase_amount];

  for (int i=0;  i < dest_phrase_amount;  i++)
  {
  chstr_result_phrases[i]=new char[phrase_len_max+1];
  chstr_result_phrases[i][0]='\0';
  }


//  5) writing-distribution all (аfrom extracted and col_phrases) phrases to chstr_result_phrases 
//
  if (col_num < 0) {col_num=col_count;}
  if (col_num > col_count) {col_num=col_count;}
int counter_in_col_phrases=0;
int total_phrase_counter=0;

  for (int i=0;  i < dest_phrase_amount;  i++)
  {

    if ((i % (col_count+1))  == col_num)
    {
    chstr_result_phrases[i][0]='\0';

      if ((col_phrases[counter_in_col_phrases]) && (strlen(col_phrases[counter_in_col_phrases]) > 0)) {strcpy(chstr_result_phrases[i], col_phrases[counter_in_col_phrases]);}  
    
    counter_in_col_phrases++;  total_phrase_counter++;
    }  
    else
    {
    chstr_result_phrases[i][0]='\0'; 

      if ((chstr_extracted_phrases[i-counter_in_col_phrases]) && (strlen(chstr_extracted_phrases[i-counter_in_col_phrases]) > 0)) 
      {strcpy(chstr_result_phrases[i], chstr_extracted_phrases[i-counter_in_col_phrases]);}  

    total_phrase_counter++;
    }

  }


//  6) writing chstr_result_phrases to file 
//
int result_file_sz=0;

  if ((src_pre_text_caption_sz > 0) && (src_pre_text_caption )) {result_file_sz=write_chstr_to_file(src_pre_text_caption, filename, true);}
  if (total_phrase_counter > 0) 
  {
  int hor_interval_len=2, ver_interval_len=0;
  bool is_new_file=false;
    if (result_file_sz < 1) {is_new_file=true;}

    if (type == 0)
    {result_file_sz=write_chstr_phrases_to_file_by_align_col_3(filename,  chstr_result_phrases, total_phrase_counter, col_count+1, hor_interval_len, ver_interval_len, is_new_file);}
    if (type == 1)
    {result_file_sz=write_chstr_phrases_to_file_by_align_col_4(filename,  chstr_result_phrases, total_phrase_counter, col_count+1, hor_interval_len, ver_interval_len, is_new_file, 0);}

  }


// cleaning
free_chstr_array(chstr_result_phrases, dest_phrase_amount);
  if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  if (src_text) {delete[] src_text;  src_text=0;}


return result_file_sz;

}



// ==========================================================================
// ==========================================================================



template<typename T> int add_dbl_array_column_to_file(T *dbl_col_phrases, char *filename, int amount_in_col, int prec, int wished_len, int start_line_num, char *forbid_1st_symbs, int col_num, char type)
{

  if ((dbl_col_phrases == 0) || (filename == 0) || (amount_in_col < 1) || (start_line_num < 0)) {return 0;}
  if (type < 0) {type=0;}
  if (type > 1) {type=1;}


char **col_phrases=0;
int temp_0=0;                                                        // max_len     prec,  wished_len
write_dbl_array_to_chstr(dbl_col_phrases, col_phrases, amount_in_col,   temp_0,     3,     wished_len);



char *src_pre_text_caption=0;
char *src_text=0;
bool not_take_empty_lines=true;


//  1) extract src file text
//
int src_pre_text_caption_sz=0;
int temp1_not_using=0;
  if (start_line_num > 0) {extract_lines_as_text_from_file(filename, src_pre_text_caption, temp1_not_using, 0, start_line_num);}
int extracted_src_file_sz=write_file_to_chstr(filename,  src_text,  forbid_1st_symbs,  start_line_num,  not_take_empty_lines);

  if (src_text == 0)  // || (extracted_src_file_sz < 1)) 
  {
    if (src_text) {delete[] src_text;  src_text=0;}  
    if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  return 0;
  }  


//  2) extracted src phrases
//
int col_count=calc_phrases_amount_in_line(src_text, start_line_num);
char **chstr_extracted_phrases=0;
int amount_of_extracted_phrases=0;
int max_extracted_phrase_len=0;                                                                                       //    start phrase num v     v col_count   
int src_phrase_amount=0;
  if (extracted_src_file_sz > 1) {src_phrase_amount=extract_phrases_from_chstr(src_text, chstr_extracted_phrases, amount_of_extracted_phrases,  max_extracted_phrase_len,  0,    1);} 
int dest_phrase_amount=dest_phrase_amount=amount_of_extracted_phrases+amount_in_col;

  if ((col_count > 0) && (amount_of_extracted_phrases/col_count != amount_in_col)) 
  {
  free_chstr_array(chstr_extracted_phrases, amount_of_extracted_phrases);
    if (src_text) {delete[] src_text;  src_text=0;}  
    if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  return 0;
  }  


//  3) define new max len
//
int added_col_phrase_len_max=calc_max_len_of_phrases_in_phrase_ar(col_phrases, amount_in_col);      // proc returns sum of max len of phrases in row
int phrase_len_max=0;

  if (added_col_phrase_len_max < max_extracted_phrase_len) {phrase_len_max=max_extracted_phrase_len;} else {phrase_len_max=added_col_phrase_len_max;}


//  4) memory allocation for new phrase array
//
char **chstr_result_phrases=0;
chstr_result_phrases=new char *[dest_phrase_amount];

  for (int i=0;  i < dest_phrase_amount;  i++)
  {
  chstr_result_phrases[i]=new char[phrase_len_max+1];
  chstr_result_phrases[i][0]='\0';
  }


//  5) writing-distribution all (аfrom extracted and col_phrases) phrases to chstr_result_phrases 
//
  if (col_num < 0) {col_num=col_count;}
  if (col_num > col_count) {col_num=col_count;}
int counter_in_col_phrases=0;
int total_phrase_counter=0;

  for (int i=0;  i < dest_phrase_amount;  i++)
  {

    if ((i % (col_count+1))  == col_num)
    {
    chstr_result_phrases[i][0]='\0';

      if ((col_phrases[counter_in_col_phrases]) && (strlen(col_phrases[counter_in_col_phrases]) > 0)) {strcpy(chstr_result_phrases[i], col_phrases[counter_in_col_phrases]);}  
    
    counter_in_col_phrases++;  total_phrase_counter++;
    }  
    else
    {
    chstr_result_phrases[i][0]='\0'; 

      if ((chstr_extracted_phrases[i-counter_in_col_phrases]) && (strlen(chstr_extracted_phrases[i-counter_in_col_phrases]) > 0)) 
      {strcpy(chstr_result_phrases[i], chstr_extracted_phrases[i-counter_in_col_phrases]);}  

    total_phrase_counter++;
    }

  }


//  6) writing chstr_result_phrases to file 
//
int result_file_sz=0;

  if ((src_pre_text_caption_sz > 0) && (src_pre_text_caption )) {result_file_sz=write_chstr_to_file(src_pre_text_caption, filename, true);}
  if (total_phrase_counter > 0) 
  {
  int hor_interval_len=2, ver_interval_len=0;
  bool is_new_file=false;
    if (result_file_sz < 1) {is_new_file=true;}

    if (type == 0)
    {result_file_sz=write_chstr_phrases_to_file_by_align_col_3(filename,  chstr_result_phrases, total_phrase_counter, col_count+1, hor_interval_len, ver_interval_len, is_new_file);}
    if (type == 1)
    {result_file_sz=write_chstr_phrases_to_file_by_align_col_4(filename,  chstr_result_phrases, total_phrase_counter, col_count+1, hor_interval_len, ver_interval_len, is_new_file, 0);}

  }


// cleaning
free_chstr_array(chstr_result_phrases, dest_phrase_amount);
  if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  if (src_text) {delete[] src_text;  src_text=0;}


return result_file_sz;

}



// ==========================================================================
// ==========================================================================



template<typename T> int add_int_array_column_to_file(T * int_col_phrases, char * filename, int amount_in_col, int wished_len, int start_line_num, char *forbid_1st_symbs, int col_num, char type)
{

  if ((int_col_phrases == 0) || (filename == 0) || (amount_in_col < 1) || (start_line_num < 0)) {return 0;}
  if (type < 0) {type=0;}
  if (type > 1) {type=1;}


char **col_phrases=0;
int temp_0=0;                                                        // max_len      wished_len
write_int_array_to_chstr(int_col_phrases, col_phrases, amount_in_col,   temp_0,      wished_len);


char *src_pre_text_caption=0;
char *src_text=0;
bool not_take_empty_lines=true;


//  1) extract src file text
//
int src_pre_text_caption_sz=0;
int temp1_not_using=0;
  if (start_line_num > 0) {extract_lines_as_text_from_file(filename, src_pre_text_caption, temp1_not_using, 0, start_line_num);}
int extracted_src_file_sz=write_file_to_chstr(filename,  src_text,  forbid_1st_symbs,  start_line_num,  not_take_empty_lines);

  if (src_text == 0)  // || (extracted_src_file_sz < 1)) 
  {
    if (src_text) {delete[] src_text;  src_text=0;}  
    if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  return 0;
  }  


//  2) extracted src phrases
//
int col_count=calc_phrases_amount_in_line(src_text, start_line_num);
char **chstr_extracted_phrases=0;
int amount_of_extracted_phrases=0;
int max_extracted_phrase_len=0;                                                                                       //    start phrase num v     v col_count   
int src_phrase_amount=0;
  if (extracted_src_file_sz > 1) {src_phrase_amount=extract_phrases_from_chstr(src_text, chstr_extracted_phrases, amount_of_extracted_phrases,  max_extracted_phrase_len,  0,    1);} 
int dest_phrase_amount=dest_phrase_amount=amount_of_extracted_phrases+amount_in_col;

  if ((col_count > 0) && (amount_of_extracted_phrases/col_count != amount_in_col)) 
  {
  free_chstr_array(chstr_extracted_phrases, amount_of_extracted_phrases);
    if (src_text) {delete[] src_text;  src_text=0;}  
    if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  return 0;
  }  


//  3) define new max len
//
int added_col_phrase_len_max=calc_max_len_of_phrases_in_phrase_ar(col_phrases, amount_in_col);      // proc returns sum of max len of phrases in row
int phrase_len_max=0;

  if (added_col_phrase_len_max < max_extracted_phrase_len) {phrase_len_max=max_extracted_phrase_len;} else {phrase_len_max=added_col_phrase_len_max;}


//  4) memory allocation for new phrase array
//
char **chstr_result_phrases=0;
chstr_result_phrases=new char *[dest_phrase_amount];

  for (int i=0;  i < dest_phrase_amount;  i++)
  {
  chstr_result_phrases[i]=new char[phrase_len_max+1];
  chstr_result_phrases[i][0]='\0';
  }


//  5) writing-distribution all (аfrom extracted and col_phrases) phrases to chstr_result_phrases 
//
  if (col_num < 0) {col_num=col_count;}
  if (col_num > col_count) {col_num=col_count;}
int counter_in_col_phrases=0;
int total_phrase_counter=0;

  for (int i=0;  i < dest_phrase_amount;  i++)
  {

    if ((i % (col_count+1))  == col_num)
    {
    chstr_result_phrases[i][0]='\0';

      if ((col_phrases[counter_in_col_phrases]) && (strlen(col_phrases[counter_in_col_phrases]) > 0)) {strcpy(chstr_result_phrases[i], col_phrases[counter_in_col_phrases]);}  
    
    counter_in_col_phrases++;  total_phrase_counter++;
    }  
    else
    {
    chstr_result_phrases[i][0]='\0'; 

      if ((chstr_extracted_phrases[i-counter_in_col_phrases]) && (strlen(chstr_extracted_phrases[i-counter_in_col_phrases]) > 0)) 
      {strcpy(chstr_result_phrases[i], chstr_extracted_phrases[i-counter_in_col_phrases]);}  

    total_phrase_counter++;
    }

  }


//  6) writing chstr_result_phrases to file 
//
int result_file_sz=0;

  if ((src_pre_text_caption_sz > 0) && (src_pre_text_caption )) {result_file_sz=write_chstr_to_file(src_pre_text_caption, filename, true);}
  if (total_phrase_counter > 0) 
  {
  int hor_interval_len=2, ver_interval_len=0;
  bool is_new_file=false;
    if (result_file_sz < 1) {is_new_file=true;}

    if (type == 0)
    {result_file_sz=write_chstr_phrases_to_file_by_align_col_3(filename,  chstr_result_phrases, total_phrase_counter, col_count+1, hor_interval_len, ver_interval_len, is_new_file);}
    if (type == 1)
    {result_file_sz=write_chstr_phrases_to_file_by_align_col_4(filename,  chstr_result_phrases, total_phrase_counter, col_count+1, hor_interval_len, ver_interval_len, is_new_file, 0);}

  }


// cleaning
free_chstr_array(chstr_result_phrases, dest_phrase_amount);
  if (src_pre_text_caption) {delete[] src_pre_text_caption;  src_pre_text_caption=0;}
  if (src_text) {delete[] src_text;  src_text=0;}


return result_file_sz;

}



// ==========================================================================
// ==========================================================================



int extract_bool_array_from_chstr_text_in_line_after_key_phrase(char * chstr_text, char * key_phrase, bool * const bool_array, int given_amount_of_array_elements)
{ 

  if (bool_array == 0) {return 0;}
  if (chstr_text == 0) {return 0;}
  if (key_phrase == 0) {return 0;}  
  if (given_amount_of_array_elements == 0) {return 0;}
  

bool *bool_array_temp=0;  
char *chstr_extracted_text=0;

int   returning_len=0;
int   line_number=0;
int  amount_of_extracted_num=0;


line_number=sr_line_number_for_line_with_given_phrase_from_chstr(chstr_text, key_phrase, 0);
                
  if (line_number >= 0) {returning_len=extract_line_with_given_number_from_chstr(chstr_text, line_number+1, chstr_extracted_text);}  
  
  
  
  
  if (returning_len > 0)
  {
  int amount_of_extracted_num_loc=0;
  bool ok_extracted=extract_bool_array_from_chstr(chstr_extracted_text,  bool_array_temp, amount_of_extracted_num_loc);     //  !  !  !        false means here type int instead of dbl  
  
  amount_of_extracted_num=amount_of_extracted_num_loc;   
  
  
  
    if (given_amount_of_array_elements > -1)       //  protecting
    {
      if (given_amount_of_array_elements != amount_of_extracted_num) {amount_of_extracted_num=0; }    
    }  
    
    if (ok_extracted == false) {amount_of_extracted_num=0;}       //  protecting
    

    if (amount_of_extracted_num > 0)
    {
      for (int i=0;  i < amount_of_extracted_num;  i++)
      {bool_array[i]=bool_array_temp[i];}
    }  
    
  }


  
  if (chstr_extracted_text) {delete[] chstr_extracted_text;  chstr_extracted_text=0;}
  if (bool_array_temp) {delete[] bool_array_temp;  bool_array_temp=0;}
    
    
return amount_of_extracted_num;

}    
    
  


// ==========================================================================
// ==========================================================================







// ==========================================================================
// ==========================================================================



}



//
//                        _________________________|                 
//  ---------------                    
//                |                      
//                 \____________|        
//                     |
//     _______ _______________          ---------------
//                     
















//int main(int argn,  char **argc)
//{ 

//cout<<endl<<"=============================================================== start =============================================================== "<<endl;
//cout<<endl; 

//clock_t start=0, end=0;  
//start = clock();
//srand( (unsigned) time(NULL) ); 

/*
// * * * * *
// example 1

char chstr_filename_1[]="L9x9_E_M_Histogram__fer_cycly_.txt\0";
char *chstr_lines_1=0;
int amount_of_lines_1=0;
int line_size_1=0;
amount_of_lines_1=extract_lines_from_file(chstr_filename_1,  chstr_lines_1,  amount_of_lines_1, line_size_1,  0);
cout<<"amount_of_lines="<<amount_of_lines_1<<"  "<<"line_size="<<line_size_1<<endl;
display_chstr_lines(chstr_lines_1,  amount_of_lines_1, line_size_1);
  if (chstr_lines_1) {delete[] chstr_lines_1;  chstr_lines_1=0;}


// * * * * *
// example 2

char chstr_filename_2[]="L9x9_E_M_Histogram__fer_cycly_.txt\0";
char **chstr_lines_2=0;
int amount_of_lines_2=0;
int line_size_2=0;
char chstr_forbid_1st_line_symb_list_2[]="3\0";
amount_of_lines_2=extract_lines_from_file(chstr_filename_2,  chstr_lines_2,  amount_of_lines_2, line_size_2,  chstr_forbid_1st_line_symb_list_2,  0);
cout<<"amount_of_lines="<<amount_of_lines_2<<"  "<<"line_size="<<line_size_2<<endl;
display_chstr_lines(chstr_lines_2,  amount_of_lines_2);
free_chstr_array(chstr_lines_2, amount_of_lines_2);
*/


/*
// * * * * *
// example 3

char chstr_filename_3[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
char **chstr_nums_3=0;
int amount_of_nums_3=0;
int max_num_size_3=0;
amount_of_nums_3=extract_anytypenums_from_file(chstr_filename_3,  chstr_nums_3,  amount_of_nums_3, max_num_size_3,  0);
cout<<"amount_of_nums="<<amount_of_nums_3<<"  "<<"max_num_size="<<max_num_size_3<<endl;
display_chstr_phrases_by_col(chstr_nums_3,  amount_of_nums_3, 3, (char *) "\t\t");
cout<<"writiting in file:  ";
write_chstr_phrases_to_file_by_col((char *) "result_3.txt\0",  chstr_nums_3,  amount_of_nums_3, 3, (char *) "\t\t", true);
free_chstr_array(chstr_nums_3, amount_of_nums_3);
*/


/*
// * * * * *
// example 4

char chstr_filename_4[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
char **chstr_nums_4=0;
int amount_of_nums_4=0;
int max_num_size_4=0;
amount_of_nums_4=extract_anytypenums_from_file(chstr_filename_4,  chstr_nums_4,  amount_of_nums_4, max_num_size_4,  0);
cout<<"amount_of_nums="<<amount_of_nums_4<<"  "<<"num_size="<<max_num_size_4<<endl;
write_chstr_phrases_to_screen_by_align_col_1(chstr_nums_4, amount_of_nums_4, 3, 3, 0);
cout<<"writing in file:  "<<write_chstr_phrases_to_file_by_align_col_1((char *) "result_4.txt\0",  chstr_nums_4,  amount_of_nums_4, 3, 3, 0, true)<<endl;
free_chstr_array(chstr_nums_4, amount_of_nums_4);
*/

/*
// * * * * *
// example 5

char chstr_filename_5[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
char **chstr_nums_5=0;
int amount_of_nums_5=0;
int max_num_size_5=0;
amount_of_nums_5=extract_anytypenums_from_file(chstr_filename_5,  chstr_nums_5,  amount_of_nums_5, max_num_size_5,  0);
cout<<"amount_of_nums="<<amount_of_nums_5<<"  "<<"num_size="<<max_num_size_5<<endl;
write_chstr_phrases_to_screen_by_align_col_3(chstr_nums_5, amount_of_nums_5, 3, 3, 0);
cout<<"writing in file:  "<<write_chstr_phrases_to_file_by_align_col_3((char *) "result_5.txt\0",  chstr_nums_5,  amount_of_nums_5, 3, 1, 0, true)<<endl;

free_chstr_array(chstr_nums_5, amount_of_nums_5);
*/


/*
// * * * * *
// example 6

char chstr_filename_6[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
char **chstr_nums_6=0;
int amount_of_nums_6=0;
int max_num_size_6=0;
amount_of_nums_6=extract_anytypenums_from_file(chstr_filename_6,  chstr_nums_6,  amount_of_nums_6, max_num_size_6,  0);
cout<<"amount_of_nums="<<amount_of_nums_6<<"  "<<"num_size="<<max_num_size_6<<endl;
write_chstr_phrases_to_screen_by_align_col_2(chstr_nums_6, amount_of_nums_6, 3, 3, 0);
cout<<"writing in file:  "<<write_chstr_phrases_to_file_by_align_col_2((char *) "result_6.txt\0",  chstr_nums_6,  amount_of_nums_6, 3, 1, 0, true)<<endl;
free_chstr_array(chstr_nums_6, amount_of_nums_6);
*/


/*
// * * * * *
// example 7

char chstr_filename_7[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
char *chstr_var_7=0;
write_file_to_chstr(chstr_filename_7, chstr_var_7);   int h=strlen(chstr_var_7);
cout<<chstr_var_7<<endl;
  if (chstr_var_7) {delete[] chstr_var_7;   chstr_var_7=0;}
*/


/*
// * * * * *
// example 8

char chstr_filename_8[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
char **chstr_nums_8=0;
int amount_of_nums_8=0;
int max_num_size_8=0;
amount_of_nums_8=extract_anytypenums_from_file(chstr_filename_8,  chstr_nums_8,  amount_of_nums_8, max_num_size_8,  0);
cout<<"amount_of_nums="<<amount_of_nums_8<<"  "<<"num_size="<<max_num_size_8<<endl;
write_chstr_phrases_to_screen_by_align_col_4(chstr_nums_8, amount_of_nums_8, 3, 3, 1, 1);
cout<<"writing in file:  "<<write_chstr_phrases_to_file_by_align_col_4((char *) "result_8.txt\0",  chstr_nums_8,  amount_of_nums_8, 3, 3, 1, true, 1)<<endl;
free_chstr_array(chstr_nums_8, amount_of_nums_8);
*/


// * * * * *
// example 9
/*
char chstr_filename_9[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
char *chstr_text_9=0;
cout<<"writing in file:  "<<write_file_to_chstr(chstr_filename_9, chstr_text_9,  (char *) "012345678",  0,  true)<<"  bytes"<<endl;
cout<<chstr_text_9<<endl;
  if (chstr_text_9) {delete[] chstr_text_9;  chstr_text_9=0;}
*/


/*
// * * * * *
// example 10

char chstr_filename_10[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
char **chstr_phrases_10=0;
char **chstr_phrases_10_2=0;
int amount_of_phrases_10=0;
int max_phrase_len_10=0;
int col_start_number_10=2;
extract_column_from_file(chstr_filename_10, chstr_phrases_10  , amount_of_phrases_10,  max_phrase_len_10, 0, -1, (char *) "#/",  0);
extract_column_from_file(chstr_filename_10, chstr_phrases_10_2, amount_of_phrases_10,  max_phrase_len_10, 2, -1, (char *) "#/",  0);
cout<<" amount_of_phrases_10 "<< amount_of_phrases_10<<endl;
write_columns_to_screen_3(amount_of_phrases_10, 0, 2, false, chstr_phrases_10, chstr_phrases_10_2);
free_chstr_array(chstr_phrases_10, amount_of_phrases_10);
free_chstr_array(chstr_phrases_10_2, amount_of_phrases_10);
*/


/*
// * * * * *
// example 11

//  test: template< typename T > bool extract_dbl_column_from_file(char * filename, T *&col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);

char chstr_filename_11[]="L9x9_E_M_Histogram__fer_cycly.txt\0";
double * dbl_phrases_11=0;
int amount_of_phrases_11=0;
int max_phrase_len_11=0;
int start_line_num_11=0;
extract_dbl_column_from_file(chstr_filename_11, dbl_phrases_11,  amount_of_phrases_11,  0, 3, (char *) "#/",  0);      //  ! ! !
cout<<" amount_of_phrases_11 "<< amount_of_phrases_11<<endl;
write_dbl_array_to_screen(dbl_phrases_11, amount_of_phrases_11, 0, true, true);
  if (dbl_phrases_11) {delete[] dbl_phrases_11;  dbl_phrases_11=0;}
*/


/*
// * * * * *
// example 12

// bool extract_column_from_file(char * filename, char **&col_phrases, int & amount_in_col,  int & max_phrase_len, int col_start_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/", int start_line_number=0);
// test: int add_chstr_array_column_to_file(char ** col_phrases, char * filename, int amount_in_col,  int start_line_num=0, char *forbid_1st_symbs=(char *) "#/",  char type=0, int col_num=-1);

char chstr_filename_12[]="L9x9_E_M_Histogram__fer_cycly_.txt\0";
char ** chstr_extracted_phrases_12=0;
int amount_of_phrases_12=0;
int max_phrase_len_12=0;
int start_line_num_12=0;  
extract_column_from_file(chstr_filename_12,  chstr_extracted_phrases_12, amount_of_phrases_12,  max_phrase_len_12, 2, -1,  (char *) "#/",  0);
int file_size_12=add_chstr_array_column_to_file(chstr_extracted_phrases_12, chstr_filename_12, amount_of_phrases_12,  0, (char *) "#/", 0, 0);
free_chstr_array(chstr_extracted_phrases_12, amount_of_phrases_12);
*/


/*
// * * * * *
// example 13

char chstr_filename_13[]="L9x9_E_M_Histogram__fer_cycly_.txt\0";
char ** chstr_extracted_phrases_13=0;
int amount_of_phrases_13=0;
int amount_in_col_13=3;
int max_phrase_len_13=0;
int start_line_num_13=0;  

mpf_t *x_mpf=0;
x_mpf=new mpf_t[amount_in_col_13];
  for (int i=0;  i < amount_in_col_13;  i++)
  {mpf_init(x_mpf[i]);  mpf_set_d(x_mpf[i], (double) i);}

//int add_mpf_array_column_to_file(T *col_phrases, char *filename, int amount_in_col, int prec=3, int wished_len=8, int start_line_num=0, char *forbid=(char *) "#/", int col_num=-1, char type=0);

int file_size_13=add_mpf_array_column_to_file(x_mpf, chstr_filename_13, amount_in_col_13,  3, 0, (char *) "#/", 0, 0);
cout<<"written file size = "<<file_size_13<<endl;

  if (x_mpf)
  {

    for (int i=0;  i < amount_in_col_13;  i++)
    {mpf_clear(x_mpf[i]);}

  delete[] x_mpf;  x_mpf=0;
  }

*/



/*
// * * * * *
// example 15

//  bool extract_mpf_column_from_file(char * filename, mpf_t *& col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);
//  int add_mpz_array_column_to_file(mpz_t * mpz_col_phrases, char * filename, int amount_in_col,  int wished_len=-1,  int start_line_num=0, char *forbid_1st_symbs=(char *) "#/", int col_num=-1, char type=0);

char chstr_filename_15[]="L9x9_E_M_Histogram__fer_cycly_.txt\0";
char ** chstr_extracted_phrases_15=0;
int amount_of_phrases_15=0;
int amount_in_col_15=3;
int max_phrase_len_15=0;
int start_line_num_15=0;  
mpz_t *col_phrases_15;

extract_mpz_column_from_file(chstr_filename_15, col_phrases_15, amount_in_col_15, 1, -1, (char *) "#/",  0);
  for (int i=0;  i < amount_in_col_15;  i++)
  {cout<<"  "<<col_phrases_15[i];}
cout<<endl;
int file_size_15=add_mpz_array_column_to_file(col_phrases_15, chstr_filename_15, amount_in_col_15,  -1, 0, (char *) "#/", 0, 0);
cout<<"written file size = "<<file_size_15<<endl;

  if (col_phrases_15)
  {

    for (int i=0;  i < amount_in_col_15;  i++)
    {mpz_clear(col_phrases_15[i]);}

  delete[] col_phrases_15;  col_phrases_15=0;
  }
*/



//end = clock();
//printf("work took %d milliseconds\n", (int)((end-start)*1E3/CLOCKS_PER_SEC));
//cout<<endl;



//cout<<endl<<"=============================================================== end =============================================================== "<<endl<<endl;

//return 0;

//}








  



#endif

//
