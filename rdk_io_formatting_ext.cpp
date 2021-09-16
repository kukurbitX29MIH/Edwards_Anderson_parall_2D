#ifndef __RDK_IO_FORMATTING_EXT_CPP
#define __RDK_IO_FORMATTING_EXT_CPP

// Last modified 03.11.2020 PadalkoMA
// Last modified 03.11.2020 PadalkoMA

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
#include <mpfr.h>

#include "anyfunctions.cpp"
#include "rdk_io_formatting.cpp"

#define _MACR_GMP_CHSTR_MAX_SZ_ 1024

using namespace std;



namespace _NMSP_RDK_
{


#ifdef __GMP_H__
bool extract_mpf_column_from_file(char * filename, mpf_t *& col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);
bool extract_mpz_column_from_file(char * filename, mpz_t *& col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);
#endif
#ifdef __MPFR_H
bool extract_mpfr_column_from_file(char *filename, mpfr_t *& col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);
#endif


#ifdef __GMP_H__
bool write_mpf_array_to_chstr(mpf_t * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len,  int wished_len=-1);
bool write_mpz_array_to_chstr(mpz_t * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len,  int wished_len=-1);
#endif
#ifdef __MPFR_H
int mpfr_to_chstr_1_num(mpfr_t  *mpfr_num,  char * chstr_num, int wished_len=-1);
bool write_mpfr_array_to_chstr(mpfr_t * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len,  int wished_len=-1);
#endif


#ifdef __GMP_H__
void write_mpf_array_to_screen(mpf_t * col_phrases, int amount_of_elem, int start_elem_num=0, bool hor1ver0=true, bool do_enum=true);
void write_mpz_array_to_screen(mpz_t * col_phrases, int amount_of_elem, int start_elem_num=0, bool hor1ver0=true, bool do_enum=true);
#endif
#ifdef __MPFR_H
void write_mpfr_1_num_to_screen(mpfr_t  *mpfr_num,  int wished_len=-1);
void write_mpfr_array_to_screen(mpfr_t * col_phrases, int amount_of_elem, int start_elem_num=0, bool hor1ver0=true, bool do_enum=true, int wished_len=3);
#endif


#ifdef __GMP_H__
int add_mpf_array_column_to_file(mpf_t * mpf_col_phrases, char * filename, int amount_in_col,  int wished_len=-1,  int start_line_num=0, char *forbid_1st_symbs=(char *) "#/", int col_num=-1, char type=0);
int add_mpz_array_column_to_file(mpz_t * mpz_col_phrases, char * filename, int amount_in_col,  int wished_len=-1,  int start_line_num=0, char *forbid_1st_symbs=(char *) "#/", int col_num=-1, char type=0);
#endif
#ifdef __MPFR_H
int add_mpfr_array_column_to_file(mpfr_t * mpf_col_phrases, char * filename, int amount_in_col,  int wished_len=-1,  int start_line_num=0, char *forbid_1st_symbs=(char *) "#/", int col_num=-1, char type=0);
#endif






//  u   u   u   u   u   u   u   u   u   u   u   u   u   u   u   u
//
//  u   u   u   u   u   u   u   u   u   u   u   u   u   u   u   u






// ==========================================================================
// ==========================================================================



#ifdef __GMP_H__

bool extract_mpf_column_from_file(char * filename,  mpf_t *& col_phrases,  int & amount_in_col,  int col_number,  int col_count,  char *forbid_1st_symbs,  int start_line_num)
{

  if ((filename == 0) || (col_phrases) || (col_count == 0) || ((col_count >= 0) && (col_number+1 > col_count))) {return false;}
 

amount_in_col=0;
char **chstr_col_phrases=0;
int max_phrase_len=0;
bool ok=extract_column_from_file(filename,  chstr_col_phrases,  amount_in_col,  max_phrase_len,  col_number,  col_count,  forbid_1st_symbs,  start_line_num);

  if ((ok == false) || (amount_in_col < 1)) {free_chstr_array(chstr_col_phrases, amount_in_col);  return false;}


double bit_amount_per_1_el__double=log((double) max_phrase_len)/log(2.0)+1.0;
int    bit_amount_per_1_el__int=floor(bit_amount_per_1_el__int)+10+1;      // mantissa bit len + exp bit len + sign

col_phrases=new mpf_t[amount_in_col]; 

  for (int i=0;  i < amount_in_col; i++) 
  {mpf_init2(col_phrases[i], bit_amount_per_1_el__int);}


  for (int i=0;  i < amount_in_col; i++)
  {mpf_set_str(col_phrases[i],  chstr_col_phrases[i], 10);}


free_chstr_array(chstr_col_phrases, amount_in_col);

return true;

}

#endif



// ==========================================================================
// ==========================================================================



#ifdef __GMP_H__

bool extract_mpz_column_from_file(char * filename, mpz_t *& col_phrases,  int & amount_in_col,  int col_number, int col_count, char *forbid_1st_symbs,  int start_line_num)
{

  if ((filename == 0) || (col_phrases) || (col_count == 0) || ((col_count >= 0) && (col_number+1 > col_count))) {return false;}


amount_in_col=0;
char **chstr_col_phrases=0;
int max_phrase_len=0;
bool ok=extract_column_from_file(filename,  chstr_col_phrases,  amount_in_col,  max_phrase_len,  col_number,  col_count,  forbid_1st_symbs,  start_line_num);

  if ((ok == false) || (amount_in_col < 1)) {free_chstr_array(chstr_col_phrases, amount_in_col);  return false;}


double bit_amount_per_1_el__double=log((double) max_phrase_len)/log(2.0)+1.0;
int    bit_amount_per_1_el__int=floor(bit_amount_per_1_el__int)+3;      // mantissa bit len + exp bit len + sign

col_phrases=new mpz_t[amount_in_col]; 

  for (int i=0;  i < amount_in_col; i++)
  {mpz_init2(col_phrases[i], bit_amount_per_1_el__int);}


  for (int i=0;  i < amount_in_col; i++)
  {mpz_set_str(col_phrases[i],  chstr_col_phrases[i], 10);}


free_chstr_array(chstr_col_phrases, amount_in_col);

return true;

}

#endif



// ==========================================================================
// ==========================================================================



#ifdef __MPFR_H

bool extract_mpfr_column_from_file(char *filename, mpfr_t *& col_phrases, int & amount_in_col, int col_number, int col_count, char *forbid_1st_symbs,  int start_line_num)
{

  if ((filename == 0) || (col_phrases) || (col_count == 0) || ((col_count >= 0) && (col_number+1 > col_count))) {return false;}
 

amount_in_col=0;
char **chstr_col_phrases=0;
int max_phrase_len=0;
bool ok=extract_column_from_file(filename,  chstr_col_phrases,  amount_in_col,  max_phrase_len,  col_number,  col_count,  forbid_1st_symbs,  start_line_num);

  if ((ok == false) || (amount_in_col < 1)) {free_chstr_array(chstr_col_phrases, amount_in_col);  return false;}


double bit_amount_per_1_el__double=log((double) max_phrase_len)/log(2.0)+1.0;
int    bit_amount_per_1_el__int=floor(bit_amount_per_1_el__int)+10+1;      // mantissa bit len + exp bit len + sign

col_phrases=new mpfr_t[amount_in_col]; 

  for (int i=0;  i < amount_in_col; i++) 
  {mpfr_init2(col_phrases[i], bit_amount_per_1_el__int);}


  for (int i=0;  i < amount_in_col; i++)
  {mpfr_set_str(col_phrases[i],  chstr_col_phrases[i], 10, MPFR_RNDZ);}


free_chstr_array(chstr_col_phrases, amount_in_col);

return true;

}

#endif



// ==========================================================================
// ==========================================================================



#ifdef __GMP_H__

bool write_mpf_array_to_chstr(mpf_t * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len,  int wished_len)
{


  if ((col_phrases == 0) || (chstr_phrases) || (amount_of_elem < 1)) {return false;}


char chstr_mantis1[_MACR_GMP_CHSTR_MAX_SZ_], chstr_mantis2[_MACR_GMP_CHSTR_MAX_SZ_], chstr_exp[128], chstr_result_num[_MACR_GMP_CHSTR_MAX_SZ_*2];
chstr_mantis1[0]='\0';  chstr_mantis2[0]='\0';  chstr_exp[0]='\0';  chstr_result_num[0]='\0';      //  mantis1 - only digits: 8253818,      mantis2 - digits with 0.:   0.8253818
  if (wished_len < 1) {wished_len=1;}


//  define max len size
//
int  elem_len=0;  
max_elem_len=0;
int  mantis_elem_len=0,  exp_elem_len=0;
mp_exp_t  mp_exp_t_var;

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_mantis1[0]='\0';  chstr_mantis2[0]='\0';  chstr_exp[0]='\0';  chstr_result_num[0]='\0';
  mpf_get_str(chstr_mantis1,  & mp_exp_t_var,  10,  wished_len,  col_phrases[i]);
  sprintf(chstr_exp, "%li", mp_exp_t_var);

  strcpy(chstr_mantis2, (char *) "0.");
  strcat(chstr_mantis2, chstr_mantis1);
  mantis_elem_len=strlen(chstr_mantis2);
  exp_elem_len=strlen(chstr_exp);

  strcat(chstr_result_num, chstr_mantis2);
  strcat(chstr_result_num, (char *) "E");

    if (((exp_elem_len > 0) && (chstr_exp[0] != '+') && (chstr_exp[0] != '-')) || (exp_elem_len < 1)) {strcat(chstr_result_num, (char *) "+");}

  strcat(chstr_result_num, chstr_exp);
  elem_len=strlen(chstr_result_num);
 
    if (max_elem_len < elem_len) {max_elem_len=elem_len;}

  }


//  memory allocation
//
chstr_phrases=new char*[amount_of_elem];

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_phrases[i]=new char[max_elem_len+1];
  chstr_phrases[i][0]='\0';
  }



//  array writing to chstr
//
  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_mantis1[0]='\0';  chstr_mantis2[0]='\0';  chstr_exp[0]='\0';  chstr_result_num[0]='\0';
  mpf_get_str(chstr_mantis1,  & mp_exp_t_var,  10,  wished_len,  col_phrases[i]);
  sprintf(chstr_exp, "%li", mp_exp_t_var);

  strcpy(chstr_mantis2, (char *) "0.");
  strcat(chstr_mantis2, chstr_mantis1);

  strcat(chstr_result_num, chstr_mantis2);
  strcat(chstr_result_num, (char *) "E");

    if (((exp_elem_len > 0) && (chstr_exp[0] != '+') && (chstr_exp[0] != '-')) || (exp_elem_len < 1)) {strcat(chstr_result_num, (char *) "+");}

  strcat(chstr_result_num, chstr_exp);
  elem_len=strlen(chstr_result_num);
  strcat(chstr_result_num, chstr_exp);

  chstr_phrases[i][0]='\0';  
  strcpy(chstr_phrases[i], chstr_result_num);
  }


return true;

}

#endif



// ==========================================================================
// ==========================================================================



#ifdef __GMP_H__

bool write_mpz_array_to_chstr(mpz_t * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len,  int wished_len)
{

  if ((col_phrases == 0) || (chstr_phrases) || (amount_of_elem < 1)) {return false;}


char  chstr_result_num[2048];
chstr_result_num[0]='\0';
  if (wished_len < 1) {wished_len=1;}


//  define max len size
//
int  elem_len=0;
max_elem_len=0;

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_result_num[0]='\0';
  mpz_get_str(chstr_result_num, 10, col_phrases[i]);
  elem_len=strlen(chstr_result_num);
 
    if (max_elem_len < elem_len) {max_elem_len=elem_len;}

  }


//  memory allocation
//
chstr_phrases=new char*[amount_of_elem];

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_phrases[i]=new char[max_elem_len+1];
  chstr_phrases[i][0]='\0';
  }


//  array writing to chstr
//
  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_result_num[0]='\0';
  mpz_get_str(chstr_result_num, 10, col_phrases[i]);
  elem_len=strlen(chstr_result_num);
 
    if (max_elem_len < elem_len) {max_elem_len=elem_len;}

  strcpy(chstr_phrases[i], chstr_result_num);
  }



return true;

}

#endif



// ==========================================================================
// ==========================================================================



#ifdef __MPFR_H

int mpfr_to_chstr_1_num(mpfr_t  *mpfr_num,  char * chstr_num,  int wished_len)
{

  if ((mpfr_num == 0) || (chstr_num == 0)) {return 0;}

char chstr_mantis1[_MACR_GMP_CHSTR_MAX_SZ_], chstr_mantis2[_MACR_GMP_CHSTR_MAX_SZ_], chstr_exp[128];
chstr_mantis1[0]='\0';  chstr_mantis2[0]='\0';  chstr_exp[0]='\0';  chstr_num[0]='\0';      //  mantis1 - only digits: 8253818,      mantis2 - digits with 0.:   0.8253818
  if (wished_len < 1) {wished_len=1;}

int  elem_len=0;  
int  mantis_elem_len=0,  exp_elem_len=0;
mp_exp_t  mp_exp_t_var;

mpfr_get_str(chstr_mantis1,  & mp_exp_t_var,  10,  wished_len,  *mpfr_num, MPFR_RNDD);
sprintf(chstr_exp, "%li", mp_exp_t_var);

int shift=0;
  if ((strlen(chstr_mantis1) > 0) && ((chstr_mantis1[0] == '+') || (chstr_mantis1[0] == '-'))) {chstr_mantis2[0]=chstr_mantis1[0];  chstr_mantis2[1]='\0';  shift=1;}
strcat(chstr_mantis2, (char *) "0.");
int mantissa1_len=strlen(chstr_mantis1);
int startpos_mant2=strlen(chstr_mantis2);
  for (int i=0;  i < (mantissa1_len-shift); i++)
  {chstr_mantis2[startpos_mant2+i]=chstr_mantis1[i+shift];}
chstr_mantis2[mantissa1_len+2]='\0';
//strcat(chstr_mantis2, chstr_mantis1);
mantis_elem_len=strlen(chstr_mantis2);
exp_elem_len=strlen(chstr_exp);

strcat(chstr_num, chstr_mantis2);
strcat(chstr_num, (char *) "E");

  if (((exp_elem_len > 0) && (chstr_exp[0] != '+') && (chstr_exp[0] != '-')) || (exp_elem_len < 1)) {strcat(chstr_num, (char *) "+");}

strcat(chstr_num, chstr_exp);
elem_len=strlen(chstr_num);

return elem_len;

}



bool write_mpfr_array_to_chstr(mpfr_t * col_phrases, char ** & chstr_phrases, int amount_of_elem,  int & max_elem_len,  int wished_len)
{

  if ((col_phrases == 0) || (chstr_phrases) || (amount_of_elem < 1)) {return false;}


char chstr_mantis1[_MACR_GMP_CHSTR_MAX_SZ_], chstr_mantis2[_MACR_GMP_CHSTR_MAX_SZ_], chstr_exp[128], chstr_result_num[_MACR_GMP_CHSTR_MAX_SZ_*1];
chstr_mantis1[0]='\0';  chstr_mantis2[0]='\0';  chstr_exp[0]='\0';  chstr_result_num[0]='\0';      //  mantis1 - only digits: 8253818,      mantis2 - digits with 0.:   0.8253818
  if (wished_len < 1) {wished_len=1;}


//  define max len size
//
int  elem_len=0;  
max_elem_len=0;
int  mantis_elem_len=0,  exp_elem_len=0;
mp_exp_t  mp_exp_t_var;

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_mantis1[0]='\0';  chstr_mantis2[0]='\0';  chstr_exp[0]='\0';  chstr_result_num[0]='\0';
  mpfr_get_str(chstr_mantis1,  & mp_exp_t_var,  10,  wished_len,  col_phrases[i], MPFR_RNDD);
  sprintf(chstr_exp, "%li", mp_exp_t_var);

  strcpy(chstr_mantis2, (char *) "0.");
  strcat(chstr_mantis2, chstr_mantis1);
  mantis_elem_len=strlen(chstr_mantis2);
  exp_elem_len=strlen(chstr_exp);

  strcat(chstr_result_num, chstr_mantis2);
  strcat(chstr_result_num, (char *) "E");

    if (((exp_elem_len > 0) && (chstr_exp[0] != '+') && (chstr_exp[0] != '-')) || (exp_elem_len < 1)) {strcat(chstr_result_num, (char *) "+");}

  strcat(chstr_result_num, chstr_exp);
  elem_len=strlen(chstr_result_num);
 
    if (max_elem_len < elem_len) {max_elem_len=elem_len;}

  }


//  memory allocation
//
chstr_phrases=new char*[amount_of_elem];

  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_phrases[i]=new char[max_elem_len+1];
  chstr_phrases[i][0]='\0';
  }



//  array writing to chstr
//
  for (int i=0;  i < amount_of_elem;  i++)
  {
  chstr_mantis1[0]='\0';  chstr_mantis2[0]='\0';  chstr_exp[0]='\0';  chstr_result_num[0]='\0';
  mpfr_get_str(chstr_mantis1,  & mp_exp_t_var,  10,  wished_len,  col_phrases[i], MPFR_RNDD);
  sprintf(chstr_exp, "%li", mp_exp_t_var);

  strcpy(chstr_mantis2, (char *) "0.");
  strcat(chstr_mantis2, chstr_mantis1);

  strcat(chstr_result_num, chstr_mantis2);
  strcat(chstr_result_num, (char *) "E");

    if (((exp_elem_len > 0) && (chstr_exp[0] != '+') && (chstr_exp[0] != '-')) || (exp_elem_len < 1)) {strcat(chstr_result_num, (char *) "+");}

  strcat(chstr_result_num, chstr_exp);
  elem_len=strlen(chstr_result_num);
  strcat(chstr_result_num, chstr_exp);

  chstr_phrases[i][0]='\0';  
  strcpy(chstr_phrases[i], chstr_result_num);
  }


return true;

}

#endif



// ==========================================================================
// ==========================================================================



#ifdef __GMP_H__

void write_mpf_array_to_screen(mpf_t * col_phrases, int amount_of_elem, int start_elem_num, bool hor1ver0, bool do_enum)
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
  cout<<"==== ! ! PROBLEM in write_mpf_array_to_screen:   empty array ! !   ----------------------------------"<<endl;
  }

}

#endif



// ==========================================================================
// ==========================================================================



#ifdef __GMP_H__

void write_mpz_array_to_screen(mpz_t * col_phrases, int amount_of_elem, int start_elem_num, bool hor1ver0, bool do_enum)
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
  cout<<"==== ! ! PROBLEM in write_mpz_array_to_screen:   empty array ! !   ----------------------------------"<<endl;
  }

}

#endif



// ==========================================================================
// ==========================================================================



void write_mpfr_1_num_to_screen(mpfr_t  *mpfr_num,  int wished_len)
{

  if (mpfr_num)
  {
  char chstr_num[_MACR_GMP_CHSTR_MAX_SZ_];  chstr_num[0]='\0';
  mpfr_to_chstr_1_num(mpfr_num,  chstr_num, wished_len);
  cout<<chstr_num;
  }

}




// ==========================================================================
// ==========================================================================



#ifdef __MPFR_H

void write_mpfr_array_to_screen(mpfr_t * col_phrases, int amount_of_elem, int start_elem_num, bool hor1ver0, bool do_enum, int wished_len)
{

  if ((col_phrases) && (amount_of_elem > 0))
  {
  char chstr_num[_MACR_GMP_CHSTR_MAX_SZ_];  chstr_num[0]='\0';

    for (int i=0;  i < amount_of_elem;  i++)
    {
    chstr_num[0]='\0';

      if (hor1ver0 == true)
      {
        if (do_enum == true) {cout<<(start_elem_num+i)<<") ";}
      mpfr_to_chstr_1_num(& col_phrases[start_elem_num+i],  chstr_num, wished_len);
      cout<<chstr_num;
        if (i < amount_of_elem-1) {cout<<"   ";} else {cout<<endl;} 
      }
      else
      {
        if (do_enum == true) {cout<<(start_elem_num+i)<<")\t";}
      mpfr_to_chstr_1_num(& col_phrases[start_elem_num+i],  chstr_num, wished_len);
      cout<<chstr_num;
        if (i < amount_of_elem-1) {cout<<"\n";} else {cout<<endl;}
      }  
    }

  }
  else
  {
  cout<<"==== ! ! PROBLEM in write_mpf_array_to_screen:   empty array ! !   ----------------------------------"<<endl;
  }

}

#endif



// ==========================================================================
// ==========================================================================



#ifdef __GMP_H__

int add_mpf_array_column_to_file(mpf_t * mpf_col_phrases, char * filename, int amount_in_col,  int wished_len,  int start_line_num, char *forbid_1st_symbs, int col_num, char type)
{

  if ((mpf_col_phrases == 0) || (filename == 0) || (amount_in_col < 1) || (start_line_num < 0)) {return 0;}
  if (type < 0) {type=0;}
  if (type > 1) {type=1;}


char **col_phrases=0;
int temp_0=0;                                                         // max_len      wished_len
write_mpf_array_to_chstr(mpf_col_phrases, col_phrases, amount_in_col,   temp_0,     wished_len); 


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

#endif



// ==========================================================================
// ==========================================================================



#ifdef __GMP_H__

int add_mpz_array_column_to_file(mpz_t * mpz_col_phrases, char * filename, int amount_in_col,  int wished_len,  int start_line_num, char *forbid_1st_symbs, int col_num, char type)
{


  if ((mpz_col_phrases == 0) || (filename == 0) || (amount_in_col < 1) || (start_line_num < 0)) {return 0;}
  if (type < 0) {type=0;}
  if (type > 1) {type=1;}


char **col_phrases=0;
int temp_0=0;                                                         // max_len        wished_len
write_mpz_array_to_chstr(mpz_col_phrases, col_phrases, amount_in_col,   temp_0,       wished_len);


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

#endif



// ==========================================================================
// ==========================================================================



#ifdef __MPFR_H

int add_mpfr_array_column_to_file(mpfr_t * mpfr_col_phrases, char * filename, int amount_in_col,  int wished_len,  int start_line_num, char *forbid_1st_symbs, int col_num, char type)
{

  if ((mpfr_col_phrases == 0) || (filename == 0) || (amount_in_col < 1) || (start_line_num < 0)) {return 0;}
  if (type < 0) {type=0;}
  if (type > 1) {type=1;}


char **col_phrases=0;
int temp_0=0;                                                         // max_len      wished_len
write_mpfr_array_to_chstr(mpfr_col_phrases, col_phrases, amount_in_col,   temp_0,     wished_len); 


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

#endif




// ==========================================================================
// ==========================================================================

}












//int main(int argn,  char **argc)
//{ 

//cout<<endl<<"=============================================================== start =============================================================== "<<endl;
//cout<<endl; 

//clock_t start=0, end=0;  
//start = clock();
//srand( (unsigned) time(NULL) ); 

//using namespace _NMSP_RDK_;

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
// example 14

//  int add_int_array_column_to_file(T *col_phrases, char * filename, int amount_in_col, int wished_len=-1, int start_line_num=0, char *forbid_1st_symbs=(char *) "#/", int col_num=-1, char type=0);

char chstr_filename_14[]="L9x9_E_M_Histogram__fer_cycly_.txt\0";
char ** chstr_extracted_phrases_14=0;
int amount_of_phrases_14=0;
int amount_in_col_14=3;
int max_phrase_len_14=0;
int start_line_num_14=0;  

int num_14_ar[3];
num_14_ar[0]=6;  num_14_ar[1]=7;  num_14_ar[2]=8;

int file_size_14=add_int_array_column_to_file(num_14_ar, chstr_filename_14, amount_in_col_14,  3, 0, (char *) "#/", 0, 0);
cout<<"written file size = "<<file_size_14<<endl;

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
mpz_t *col_phrases_15=0;

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




/*
// * * * * *
// example 16

//  bool extract_mpf_column_from_file(char * filename, mpf_t *& col_phrases, int & amount_in_col, int col_number, int col_count=-1, char *forbid_1st_symbs=(char *) "#/",  int start_line_num=0);
//  int add_mpz_array_column_to_file(mpz_t * mpz_col_phrases, char * filename, int amount_in_col,  int wished_len=-1,  int start_line_num=0, char *forbid_1st_symbs=(char *) "#/", int col_num=-1, char type=0);

char chstr_filename_16[]="L9x9_E_M_Histogram__fer_cycly_.txt\0";
char ** chstr_extracted_phrases_16=0;
int amount_of_phrases_16=0;
int amount_in_col_16=3;
int max_phrase_len_16=0;
int start_line_num_16=0;  
mpfr_t *col_phrases_16=0;

extract_mpfr_column_from_file(chstr_filename_16, col_phrases_16, amount_in_col_16, 1, -1, (char *) "#/",  0);
//  for (int i=0;  i < amount_in_col_16;  i++)
//  {cout<<"  "<<col_phrases_16[i];}
cout<<endl;
int file_size_16=add_mpfr_array_column_to_file(col_phrases_16, chstr_filename_16, amount_in_col_16,  -1, 0, (char *) "#/", 0, 0);
cout<<"written file size = "<<file_size_16<<endl;

  if (col_phrases_16)
  {

    for (int i=0;  i < amount_in_col_16;  i++)
    {mpfr_clear(col_phrases_16[i]);}

  delete[] col_phrases_16;  col_phrases_16=0;
  }
*/



//end = clock();
//printf("work took %d milliseconds\n", (int)((end-start)*1E3/CLOCKS_PER_SEC));
//cout<<endl;



//cout<<endl<<"=============================================================== end =============================================================== "<<endl<<endl;

//return 0;

//}


//
//  LD_LIBRARY_PATH=/usr/local/lib
//  export LD_LIBRARY_PATH
//  echo ${LD_LIBRARY_PATH}
//


#endif

//
