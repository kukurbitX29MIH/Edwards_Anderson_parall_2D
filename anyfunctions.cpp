#ifndef __ANYFUNCTIONS_CPP
#define __ANYFUNCTIONS_CPP

//  Last modified 25.12.2020 PadalkoMA

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <limits.h>
#include <fstream>

#define __MACR_ULLINT_MAX_VALUE 18446744073709551615



using namespace :: std;

inline void display_bool(bool x, bool is_endline=false) {  if (x == true) {cout<<"true";} else {cout<<"false";}  if (is_endline == true) {cout<<endl;}  }

//                              // 1                  // 2                    // 3                    // 4                    // 5                     // 6
// N={1, 11}
int concat11_int(int N,  char * c,  char * c1, int i2,  char * c3=0, int i4=0,  char * c5=0, int i6=0,  char * c7=0, int i8=0,  char * c9=0, int i10=0,  char * c11=0);
int concat11_llint(int N,  char * c,  char * c1, long long int i2,  char * c3=0, long long int i4=0,  char * c5=0, long long int i6=0,  char * c7=0, long long int i8=0,  char * c9=0, long long int i10=0,  char * c11=0);
int concat11_dbl(int N,  char * c,  char * c1, double d2,  char * c3=0, double d4=0,  char * c5=0, double d6=0,  char * c7=0, double d8=0,  char * c9=0, double d10=0,  char * c11=0);

int write_chstr_to_file(const char * const text, const char * const filename, bool B1New0Add=true);
int write_int_to_file(int num, char *filename, bool B1New0Add=true);
int write_llint_to_file(long long int num, char *filename, bool B1New0Add=true);
int write_dbl_to_file(double num, char *filename, bool B1New0Add=true);
//int write_file_to_chstr(char *filename, char * & text);                             // linux bug because of file end problem,  return read char size without text end symbol
//int write_file_to_chstr(char *filename, char *& text,  int start_line_number=0);      // recomeneded instead of the previous
int write_file_to_chstr(char *filename, char * & chstr_text); 
int write_file_to_chstr(char *filename, char *& chstr_text,  int start_line_number); 
int write_file_to_chstr_with_last_emp_lines_filter(char *filename, char *& chstr_text,  int start_line_number);
int write_file_to_chstr(char *filename, char *& text,  char *chstr_forbid_1st_line_symb_list,  int start_line_number=0,  bool not_take_empty_lines=true);
int write_dbl1couln_to_file(char *filename, int N_Point, double *data);  
int write_dbl2couln_to_file(char *filename, int N_Point, double *data_x, double *data_y);
bool to_empty_file(char *filename);
  
int calc_amount_of_lines(char *filename,  int & max_size_of_line,  bool calc_only_not_emp_lines=true);    //  analizing of file
int extract_all_lines_from_file(char *filename, char ** & chstr_line, bool calc_only_not_emp_lines=true, int *line_size=0);  //  extract lines from file. Return size of file. Each line = each separate ptr chstr_line
int extract_all_lines_from_file(char *filename, char * & chstr_line, int & line_size, bool calc_only_not_emp_lines=true);
unsigned long long int extract_line_from_file(char *filename, char *& chstr_line, int desired_line_number);      //  extract line from file. Returns size of file.
int write_file_to_int_array(char * filename,  int *& int_array,  int linenumber_we_start_with=0);      //  extract all int nums and puts it to int_array

unsigned char ChULLIntToBinChStr(unsigned long long int num, char * const chstr, bool fillzero=true, bool littleendian=true);
unsigned char ChUIntToBinChStr(unsigned int num, char * const chstr, bool fillzero=true, bool littleendian=true);
unsigned char ChUSIntToBinChStr(unsigned short int num, char * const chstr, bool fillzero=true, bool littleendian=true);
unsigned char ChUCharToBinChStr(unsigned char num, char * const chstr, bool fillzero=true, bool littleendian=true);

void VDisplayULLIntToBin(unsigned long long int num, bool fillzero=true, bool littleendian=true);
void VDisplayUIntToBin(unsigned int num, bool fillzero=true, bool littleendian=true);
void VDisplayUSIntToBin(unsigned short int num, bool fillzero=true, bool littleendian=true);
void VDisplayUCharToBin(unsigned char num, bool fillzero=true, bool littleendian=true);
void VDisplayFixLenULLIntToBin(unsigned long long int num, int len, bool fillzero, bool littleendian=true);


bool is_ch_ref_to_dbl_num(char ch);
bool is_int_or_llint_num(char *chstrnum);
bool is_fl_or_dbl_num(char *chstrnum);

void add_chstr_to_string(string * string_ptr_dest,  char *chstr_added_str,  int given_len);     //  concatenation
void add_string_to_string(string * string_ptr_dest,  string * string_ptr_src,  int given_len);  //  concatenation
void add_int_to_string(string * string_ptr, int num);  //  concatenation
void add_uint_to_string(string * string_ptr, unsigned int num);  //  concatenation
void add_llint_to_string(string * string_ptr, long long int num);  //  concatenation
void add_ullint_to_string(string * string_ptr, unsigned long long int num);  //  concatenation
void add_double_to_string(string * string_ptr, double num);  //  concatenation
void add_int_to_chstr(char * chstr, int num);
void add_ullint_to_chstr(char * chstr, unsigned long long int num);

void display_mem_size(long long int mem_size_bit, bool is_line_feed=true);
void print_time(int time_ms) {cout<<(time_ms/(3600000*24))<<" days  "<<((time_ms % (24*60*60*1000))/(60*60*1000))<<" h  "<<((time_ms % (60*60*1000))/(60*1000))<<" m  "<<((time_ms % (60*1000))/(1000))<<" s  "<<((time_ms % (1000))/1)<<" ms";}

void  display_to_screen_ullint_fix_len(unsigned long long int x,  unsigned int fix_len);
void  display_to_screen_llint_fix_len(long long int x,  unsigned int fix_len);


inline unsigned long long int treat_out_of_ullint_2_num_befor_add(unsigned long long int x1,  unsigned long long int x2,  bool *is_out_of_ullint=0);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *




//inline void display_bool(bool x, bool is_endline)
//{  if (x == true) {cout<<"true";} else {cout<<"false";}  if (is_endline == true) {cout<<endl;}  }



int concat11_int(int N,  char * c, char * c1, int i2,  char * c3, int i4,  char * c5, int i6,  char * c7, int i8,  char * c9, int i10,  char * c11)
{

  if (c1 == 0) {return 0;}
  if (N < 1) {return 0;}

char chstrmum[32];  chstrmum[0]='\0';

  if (c1 != 0) {strcat(c, c1);}
  if (N == 1 ) {return (int) strlen(c);}

sprintf(chstrmum, "%d", i2);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 2 ) {return (int) strlen(c);}

  if (c3 != 0) {strcat(c, c3);}
  if (N == 3 ) {return strlen(c);}

sprintf(chstrmum, "%d", i4);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 4 ) {return strlen(c);}

  if (c5 != 0) {strcat(c, c5);}
  if (N == 5 ) {return strlen(c);}

sprintf(chstrmum, "%d", i6);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 6 ) {return strlen(c);}

  if (c7 != 0) {strcat(c, c7);}
  if (N == 7 ) {return strlen(c);}

sprintf(chstrmum, "%d", i8);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 8 ) {return strlen(c);}

  if (c9 != 0) {strcat(c, c9);}
  if (N == 9 ) {return strlen(c);}

sprintf(chstrmum, "%d", i10);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 10) {return strlen(c);}

  if (c11 != 0) {strcat(c, c11);}
return strlen(c);

}



int concat11_llint(int N,  char * c, char * c1, long long int i2,  char * c3, long long int i4,  char * c5, long long int i6,  char * c7, long long int i8,  char * c9, long long int i10,  char * c11)
{

  if (c1 == 0) {return 0;}
  if (N < 1) {return 0;}

char chstrmum[32];  chstrmum[0]='\0';

  if (c1 != 0) {strcat(c, c1);}
  if (N == 1 ) {return strlen(c);}

sprintf(chstrmum, "%lld", i2);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 2 ) {return strlen(c);}

  if (c3 != 0) {strcat(c, c3);}
  if (N == 3 ) {return strlen(c);}

sprintf(chstrmum, "%lld", i4);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 4 ) {return strlen(c);}

  if (c5 != 0) {strcat(c, c5);}
  if (N == 5 ) {return strlen(c);}

sprintf(chstrmum, "%lld", i6);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 6 ) {return strlen(c);}

  if (c7 != 0) {strcat(c, c7);}
  if (N == 7 ) {return strlen(c);}

sprintf(chstrmum, "%lld", i8);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 8 ) {return strlen(c);}

  if (c9 != 0) {strcat(c, c9);}
  if (N == 9 ) {return strlen(c);}

sprintf(chstrmum, "%lld", i10);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 10) {return strlen(c);}

  if (c11 != 0) {strcat(c, c11);}
return strlen(c);

}



int concat11_int(int N,  char * c, char * c1, double d2,  char * c3, double d4,  char * c5, double d6,  char * c7, double d8,  char * c9, double d10,  char * c11)
{

  if (c1 == 0) {return 0;}
  if (N < 1) {return 0;}

char chstrmum[32];  chstrmum[0]='\0';

  if (c1 != 0) {strcat(c, c1);}
  if (N == 1 ) {return strlen(c);}

sprintf(chstrmum, "%E", d2);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 2 ) {return strlen(c);}

  if (c3 != 0) {strcat(c, c3);}
  if (N == 3 ) {return strlen(c);}

sprintf(chstrmum, "%E", d4);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 4 ) {return strlen(c);}

  if (c5 != 0) {strcat(c, c5);}
  if (N == 5 ) {return strlen(c);}

sprintf(chstrmum, "%E", d6);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 6 ) {return strlen(c);}

  if (c7 != 0) {strcat(c, c7);}
  if (N == 7 ) {return strlen(c);}

sprintf(chstrmum, "%E", d8);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 8 ) {return strlen(c);}

  if (c9 != 0) {strcat(c, c9);}
  if (N == 9 ) {return strlen(c);}

sprintf(chstrmum, "%E", d10);   strcat(c, chstrmum);  chstrmum[0]='\0';
  if (N == 10) {return strlen(c);}

  if (c11 != 0) {strcat(c, c11);}
return strlen(c);

}



int write_chstr_to_file(const char * const text, const char * const filename, bool B1New0Add)
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

fwrite(text, textlen, 1, file);
fclose(file);

return textlen+1;

}




int write_int_to_file(int num, char *filename, bool B1New0Add)
{

char text[32];  text[0]='\0';
sprintf(text, "%d", num); 
return write_chstr_to_file(text, filename, B1New0Add);

}



int write_llint_to_file(long long int num, char *filename, bool B1New0Add)
{

char text[32];  text[0]='\0';
sprintf(text, "%lld", num); 
return write_chstr_to_file(text, filename, B1New0Add);

}



int write_dbl_to_file(double num, char *filename, bool B1New0Add)
{

char text[32];  text[0]='\0';
sprintf(text, "%E", num); 
return write_chstr_to_file(text, filename, B1New0Add);

}


//  special case, it has linux bug with file end: memory leak happens, so it's strongly recommended the next variant
//
//spec//  int write_file_to_chstr(char *filename, char *&text)
//spec//  {
//spec//  
//spec//    if (filename == 0) {return 0;}
//spec//    if (text != 0) {return 0;}
//spec//  
//spec//  FILE *file=0;
//spec//  file=fopen(filename, "rt");
//spec//  
//spec//    if (file == 0) {return 0;}
//spec//  
//spec//  */char ch=0;
//spec//  int filesize=0;
//spec//  
//spec//    while (feof(file) == false)
//spec//    {fread(&ch, 1, 1, file);  filesize++;}
//spec//  
//spec//    if (filesize == 0) {return 0;}
//spec//  
//spec//  fseek(file, 0, SEEK_SET);
//spec//  text=new char[filesize+1];
//spec//  text[0]='\0';
//spec//  fread(text, filesize, 1, file);
//spec//  text[filesize]='\0';
//spec//  
//spec//  fclose(file);
//spec//  
//spec//  return filesize;
//spec//  
//spec//    }




//spec//  int write_file_to_chstr(char *filename, char *& text,  int start_line_number)
//spec//  {
//spec//  
//spec//    if (filename == 0) {return 0;}
//spec//    if (text != 0) {return 0;}
//spec//  
//spec//  FILE *file=0;
//spec//  file=fopen(filename, "rt");
//spec//  
//spec//    if (file == 0) {return 0;}
//spec//  
//spec//  char ch=0;
//spec//  int filesize=0;
//spec//  int start_line_pos=0, line_counter=0,  pos=0;
//spec//  
//spec//  
//spec//  // file analizing
//spec//  //  calc 4 quantities:   1) line_counter,  2) start_line_pos,  3) end_line_feed_counter,  4) filesize    
//spec//  // * *
//spec//  int end_line_feed_counter=0;      //  for Linux BUG problem solving!!! about end of file
//spec//  
//spec//    while (feof(file) == false)
//spec//    {
//spec//    filesize++;
//spec//    fread(&ch, 1, 1, file);   
//spec//      if (ch == '\n') 
//spec//      {
//spec//      line_counter++;    
//spec//        if (line_counter == start_line_number) {start_line_pos=pos+1;}   
//spec//      end_line_feed_counter++;
//spec//      }
//spec//      else {end_line_feed_counter=0;}
//spec//    pos++;
//spec//    }
//spec//  
//spec//  
//spec//    if (filesize == 0) {fclose(file);  return 0;}
//spec//  
//spec//  //  last line handling
//spec//    if ((filesize > 0) && (ch != '\n')) {line_counter++;}  
//spec//  
//spec//  //  correcting  start_line_number  with taking into account end_line_feed_counter   
//spec//  int not_empt_last_line_counter=0;
//spec//    if (end_line_feed_counter > 0) {not_empt_last_line_counter=line_counter-end_line_feed_counter+1;} else {not_empt_last_line_counter=line_counter;}
//spec//    if (start_line_number+1 > not_empt_last_line_counter) {fclose(file);  return 0;}
//spec//    if (start_line_pos+1 > filesize) {fclose(file);  return 0;}
//spec//  
//spec//  //  correcting  filesize
//spec//    if (end_line_feed_counter > 0) {filesize=filesize-end_line_feed_counter+1;}
//spec//  
//spec//  
//spec//  // recording
//spec//  // * *
//spec//  fseek(file, start_line_pos, SEEK_SET);
//spec//  int recording_size=filesize-start_line_pos;  
//spec//  //==DEBUG==//  cout<<"filesize="<<filesize<<"  recording_size="<<recording_size<<"  start_line_pos="<<start_line_pos<<" end_line_feed_counter="<<end_line_feed_counter<<endl;  //==DEBUG==//
//spec//  text=new char[recording_size+1];
//spec//  text[0]='\0';
//spec//  fread(text, recording_size, 1, file); 
//spec//  text[recording_size]='\0';
//spec//  
//spec//  
//spec//  fclose(file);
//spec//  
//spec//  return recording_size;
//spec//  
//spec//  }



int write_file_to_chstr(char *filename, char *& chstr_text)
{

  if (filename == 0) {return 0;}
  if (chstr_text != 0) {return 0;}

int file_size=0;

//streampos size;

ifstream file(filename, ios::in|ios::ate);

  if (file.is_open())
  {
    file_size=file.tellg(); 
    chstr_text=new char[file_size+1];
    file.seekg(0, ios::beg);
    file.read(chstr_text, file_size);
    chstr_text[file_size]='\0';
    file.close();

  return file_size;
  }

return 0;

}



int write_file_to_chstr(char *filename, char *& chstr_text,  int start_line_number)
{ 

  if (filename == 0) {return 0;}
  if (chstr_text != 0) {return 0;}
  if (start_line_number < 0) {start_line_number=0;}

int file_size=0;
string line;
string string_extracted_text;

//streampos size;

ifstream file(filename);  

  if (file.is_open())
  {       
  int line_counter=0;
  bool is_just_started=true;

    while (getline(file, line))
    {
      if (start_line_number <= line_counter) 
      {
        if (is_just_started == true) {is_just_started=false;} else {string_extracted_text.append((char *) "\n");}
      string_extracted_text.append(line);
      }  
    
    line_counter++;
    }

  file.close();

  int text_size=string_extracted_text.length();
  chstr_text=new char[text_size+1];  chstr_text[0]='\0';
  strcat(chstr_text, string_extracted_text.c_str());
  string_extracted_text.clear();
  
  return text_size;
  }


line.clear();

return 0;

}

/*int write_file_to_chstr(char *filename, char *& chstr_text,  int start_line_number)
{  

  if (filename == 0) {return 0;}
  if (chstr_text != 0) {return 0;}
  if (start_line_number < 0) {start_line_number=0;}

int file_size=0;
string line;
string string_extracted_text;

//streampos size;

ifstream file(filename);       

  if (file.is_open())
  {
  int line_counter=0;    
  bool is_just_started=true;

    while (getline(file, line))
    {
      if (start_line_number <= line_counter) 
      {
        if (is_just_started == true) {is_just_started=false;} else {string_extracted_text.append((char *) "\n");}  
      }  
    
    line_counter++;
    }

  file.close();

  int text_size=string_extracted_text.length();
  chstr_text=new char[text_size+1];  chstr_text[0]='\0';
  strcat(chstr_text, string_extracted_text.c_str());
  string_extracted_text.clear();
  
  return text_size;
  }


line.clear();

return 0;

}*/




int write_file_to_chstr_with_last_emp_lines_filter(char *filename, char *& chstr_text,  int start_line_number)
{

  if (filename == 0) {return 0;}
  if (chstr_text != 0) {return 0;}
  if (start_line_number < 0) {start_line_number=0;}

int file_size=0;
string line;
string string_extracted_text;

//streampos size;

ifstream file(filename);

  if (file.is_open())
  {
  int line_counter=0;
  bool is_just_started=true;

    while (getline(file, line))
    {
      if (start_line_number <= line_counter) 
      {
        if (is_just_started == true) {is_just_started=false;} else {string_extracted_text.append((char *) "\n");}
      string_extracted_text.append(line);
      }  
    
    line_counter++;
    }

  file.close();

  int text_size=string_extracted_text.length();
  chstr_text=new char[text_size+1];  chstr_text[0]='\0';
  strcat(chstr_text, string_extracted_text.c_str());
  string_extracted_text.clear();
  
  return text_size;
  }


line.clear();

return 0;

}



//  this proc uses prev     int write_file_to_chstr(char *filename, char *& text,  int start_line_number)

int write_file_to_chstr(char *filename, char *& text,  char *chstr_forbid_1st_line_symb_list,  int start_line_number,  bool not_take_empty_lines)
{

  if ((filename == 0) || (text != 0)) {return 0;} 

int forbid_len=0;
  if (chstr_forbid_1st_line_symb_list) {forbid_len=strlen(chstr_forbid_1st_line_symb_list);}
char *text_temp=0;

int read_file_sz=write_file_to_chstr(filename, text_temp, start_line_number);      //    ! ! !
  if (read_file_sz <= 0) {  if (text_temp) {delete[] text_temp;  text_temp=0;}  }

string string_var;
int line_start_pos=0,  line_end_pos=0, line_len=0,  line_counter=0;
bool is_empty=false, is_forbid=false;  


  for (int pos=0;  pos < read_file_sz;  pos++)
  {
    if ((text_temp[pos] == '\n') || ((pos+1 == read_file_sz) && (pos > 0) && (text_temp[pos] != '\n')))
    { 
      if (text_temp[pos] == '\n') {line_end_pos=pos-1;} else {line_end_pos=pos;}

    line_counter++;    line_len=line_end_pos-line_start_pos+1;  is_forbid=false;
      if (line_len <= 0) {is_empty=true;} else {is_empty=false;}
    
      if (is_empty == false)
      {
        for (int i=0;  i < forbid_len;  i++)
        {  if (text_temp[line_start_pos] == chstr_forbid_1st_line_symb_list[i]) {is_forbid=true;  break;}  }
      }

      //if (start_line_number < line_counter)     //  already all lines are valid under condit  (start_line_number < line_counter)   see: int read_file_sz=write_file_to_chstr(filename, text_temp, start_line_number); 
      {
        if (is_empty == true)
        {  if (not_take_empty_lines == false) {string_var.append((char *)  "\n");}  }
        else
        { 
          if (is_forbid == false)
          {
            for (int k=0;  k < line_len;  k++)                      //  eventually copy line  !!!
            {string_var.push_back(text_temp[line_start_pos+k]);}
          string_var.append((char *)  "\n");
          }
        }
      }

      if (pos+1 < read_file_sz) {line_start_pos=pos+1;}
    }
  }


int text_sz=string_var.length();
text=new char[text_sz+1];    text[0]='\0';
strcpy(text,  string_var.c_str());

  if (text_temp) {delete[] text_temp;  text_temp=0;}
string_var.clear();

return text_sz;

}



int write_dbl1couln_to_file(char *filename, int N_Point, double *data)  
{

  if (filename == 0) {return 0;}
  if (data == 0) {return 0;}
  if (N_Point == 0) {return 0;}

FILE *file=0;
file=fopen(filename, "wr");

  if (file == 0) {return 0;}

char chstrnum[128];    chstrnum[0]='\0';
char *text=0;
text=new char[N_Point*(32+1)];

  for (int i=0; i < N_Point; i++)
  {
  sprintf(chstrnum, "%E", data[i]);
  strcat(text, chstrnum);
  strcat(text, (char *) "\n");
  }

int filesize=strlen(text);
filesize=fwrite(text, filesize, 1, file);
fclose(file);

  if (text) {delete[] text;  text=0;}

return filesize;

}



int write_dbl2couln_to_file(char *filename, int N_Point, double *data_x, double *data_y)  
{

  if (filename == 0) {return 0;}
  if (data_x == 0) {return 0;}
  if (data_y == 0) {return 0;}
  if (N_Point == 0) {return 0;}

FILE *file=0;
file=fopen(filename, "wt");

  if (file == 0) {return 0;}

char chstrnum[128];    chstrnum[0]='\0';
char *text=0;
text=new char[N_Point*(32*2+2)];
text[0]='\0';

  for (int i=0; i < N_Point; i++)
  {
  sprintf(chstrnum, "%E", data_x[i]);
  strcat(text, chstrnum);
  strcat(text, (char *) "\t");

  sprintf(chstrnum, "%E", data_y[i]);
  strcat(text, chstrnum);
  strcat(text, (char *) "\n");
  }

int filesize=strlen(text);
filesize=fwrite(text, filesize, 1, file);
fclose(file);

  if (text) {delete[] text;  text=0;}

return filesize;

}





bool to_empty_file(char *filename)
{

  if (filename == 0) {return 0;}

FILE *file=0;
file=fopen(filename, "wt");

  if (file == 0) {return false;}

fclose(file);

  if (file == 0) {return 0;}

return true;

}




//  analizing of file
//
int calc_amount_of_lines(char *filename,  int & max_size_of_line, bool calc_only_not_emp_lines)
{  

char * text=0;
int text_size=write_file_to_chstr(filename, text);  //  cout<<"DEBUG: "<<" text_size="<<text_size<<endl;


  if (text == 0) {return 0;}

int line_counter=0;
int size_of_line=0;      //  max_size_of_line
// int pos_line_start=0, pos_line_end=0;
bool is_last_step_line_finished=true;

  if (calc_only_not_emp_lines == false)
  {

    for (int pos=0;  pos < text_size;  pos++)
    {

      if (text[pos] == '\n')
      {
      is_last_step_line_finished=true;
        if (size_of_line > max_size_of_line) {max_size_of_line=size_of_line;}
      line_counter++;
      size_of_line=0;
      }
      else
      {

        if (is_last_step_line_finished == true) 
        {
        is_last_step_line_finished=false;  
        //  pos_line_start=pos;
        }

      size_of_line++;
      }     

    }

    if ((text_size > 0) && (is_last_step_line_finished == false))
    {
    is_last_step_line_finished=true;
      if (size_of_line > max_size_of_line) {max_size_of_line=size_of_line;}
    line_counter++;
    size_of_line=0;
    }

  }
  else
  {
  bool is_not_emp_symb_in_line=false;

    for (int pos=0;  pos < text_size;  pos++)
    {

      if (text[pos] == '\n')
      {
      is_last_step_line_finished=true;
        if (size_of_line > max_size_of_line) {max_size_of_line=size_of_line;}
        if (is_not_emp_symb_in_line == true) {line_counter++;}
      size_of_line=0;  is_not_emp_symb_in_line=false;
      }
      else
      {

        if (is_last_step_line_finished == true) 
        {
        is_last_step_line_finished=false;  
        //  pos_line_start=pos;
        }

        if ((text[pos] != '\n') && (text[pos] != '\t') && (text[pos] != ' ')) {is_not_emp_symb_in_line=true;}

      size_of_line++;
      }     

    }

    if ((text_size > 0) && (is_last_step_line_finished == false))
    {
    is_last_step_line_finished=true;
      if (size_of_line > max_size_of_line) {max_size_of_line=size_of_line;}
    line_counter++;
    size_of_line=0;
    }

  }


  if (text) {delete[] text;  text=0;}

return line_counter;

}




//  extract line from file. Return size of file.
//
int extract_all_lines_from_file(char *filename, char ** & chstr_line, bool calc_only_not_emp_lines, int *line_size)
{

  if (chstr_line != 0) {return 0;}

int max_size_of_line=0;
int amount_of_lines=calc_amount_of_lines(filename,  max_size_of_line,  calc_only_not_emp_lines); 

  if (amount_of_lines <= 0) {return 0;}
  if (max_size_of_line <= 0) {return 0;}


// * * * * * *
// memory allocation
//
chstr_line=new char *[amount_of_lines];

  for (int i=0;  i < amount_of_lines;  i++)
  {
  chstr_line[i]=0;
  chstr_line[i]=new char[max_size_of_line+1];      //  +1 symbol for line-end-symbol
  chstr_line[i][0]='\0';
  }
  

// * * * * * *
// start extracting
//
char * text=0;
int text_size=write_file_to_chstr(filename, text);


  if (text == 0) {return 0;}

int line_counter=0;
int size_of_line=0;      //  max_size_of_line
int pos_line_start=0;
bool is_last_step_line_finished=true;

  if (calc_only_not_emp_lines == false)
  {

    for (int pos=0;  pos < text_size;  pos++)
    {

      if (text[pos] == '\n')
      {
      is_last_step_line_finished=true;

        //  ===
        //  writing V
        if (size_of_line <= 0) {chstr_line[line_counter][0]='\0';}
        else
        {
          if (size_of_line > max_size_of_line) {size_of_line=max_size_of_line;}

          if (line_counter <= amount_of_lines)
          { 
            for (int pos_copy=0;  pos_copy < size_of_line;  pos_copy++)
            {
            chstr_line[line_counter][pos_copy]=text[pos_line_start+pos_copy];
            }
          chstr_line[line_counter][size_of_line]='\0';
          }
        }
        //  writing ^
        //  ===

      line_counter++;
      size_of_line=0;
      }
      else
      {

        if (is_last_step_line_finished == true) 
        {
        is_last_step_line_finished=false;  
        pos_line_start=pos;
        }

      size_of_line++;
      }     

    }

    if ((text_size > 0) && (is_last_step_line_finished == false))
    {
    is_last_step_line_finished=true;

      //  ===
      //  writing V
      if (size_of_line <= 0) {chstr_line[line_counter][0]='\0';}
      else
      {
        if (size_of_line > max_size_of_line) {size_of_line=max_size_of_line;}

        if (line_counter <= amount_of_lines)
        { 
          for (int pos_copy=0;  pos_copy < size_of_line;  pos_copy++)
          {
          chstr_line[line_counter][pos_copy]=text[pos_line_start+pos_copy];
          }
        chstr_line[line_counter][size_of_line]='\0';
        }
      }
      //  writing ^
      //  ===

    line_counter++;
    size_of_line=0;
    }

  }
  else
  {
  bool is_not_emp_symb_in_line=false;

    for (int pos=0;  pos < text_size;  pos++)
    {

      if (text[pos] == '\n')
      {
      is_last_step_line_finished=true;

        if (is_not_emp_symb_in_line == true) 
        {

          //  ===
          //  writing V
          if (line_counter <= amount_of_lines)
          { 
            if (size_of_line <= 0) {chstr_line[line_counter][0]='\0';}
            else
            {
              if (size_of_line > max_size_of_line) {size_of_line=max_size_of_line;}

              for (int pos_copy=0;  pos_copy < size_of_line;  pos_copy++)
              {chstr_line[line_counter][pos_copy]=text[pos_line_start+pos_copy];}
            }
          chstr_line[line_counter][size_of_line]='\0';
          }
          //  writing ^
          //  ===

        line_counter++;
        }
      size_of_line=0;  is_not_emp_symb_in_line=false;
      }
      else
      {

        if (is_last_step_line_finished == true) 
        {
        is_last_step_line_finished=false;  
        pos_line_start=pos;
        }

        if ((text[pos] != '\n') && (text[pos] != '\t') && (text[pos] != ' ')) {is_not_emp_symb_in_line=true;}

      size_of_line++;
      }     

    }

    if ((text_size > 0) && (is_last_step_line_finished == false))
    {
    is_last_step_line_finished=true;

      if (is_not_emp_symb_in_line == true) 
      {

        //  ===
        //  writing V
        if (line_counter <= amount_of_lines)
        { 
          if (size_of_line <= 0) {chstr_line[line_counter][0]='\0';}
          else
          {
            if (size_of_line > max_size_of_line) {size_of_line=max_size_of_line;}

            for (int pos_copy=0;  pos_copy < size_of_line;  pos_copy++)
            {chstr_line[line_counter][pos_copy]=text[pos_line_start+pos_copy];}
          }
        chstr_line[line_counter][size_of_line]='\0';
        }
        //  writing ^
        //  ===

      line_counter++;
      }
    size_of_line=0;  is_not_emp_symb_in_line=false;
    }

  }

  

  if (line_size != 0) {*line_size=max_size_of_line;}

  if (text) {delete[] text;  text=0;}


return line_counter;

}





//  extract line from file. Return size of file.
//
int extract_all_lines_from_file(char *filename, char * & chstr_line, int & line_size, bool calc_only_not_emp_lines)
{

  if (chstr_line != 0) {return 0;}

int max_size_of_line=0;
int amount_of_lines=calc_amount_of_lines(filename,  max_size_of_line,  calc_only_not_emp_lines); 

  if (amount_of_lines <= 0) {return 0;}
  if (max_size_of_line <= 0) {return 0;}


// * * * * * *
// memory allocation
//
int max_size_of_line_with_lineend=max_size_of_line+1;
int alloc_size=amount_of_lines*max_size_of_line_with_lineend;
chstr_line=new char[alloc_size];

  for (int i=0;  i < alloc_size; i++)
  {chstr_line[i]='\0';}


// * * * * * *
// start extracting
//
char * text=0;
int text_size=write_file_to_chstr(filename, text);


  if (text == 0) {return 0;}

int line_counter=0;
int size_of_line=0;      //  max_size_of_line
int pos_line_start=0;
bool is_last_step_line_finished=true;

  if (calc_only_not_emp_lines == false)
  {

    for (int pos=0;  pos < text_size;  pos++)
    {

      if (text[pos] == '\n')
      {
      is_last_step_line_finished=true;

        //  ===
        //  writing V
        if (size_of_line > 0)
        {
          if (size_of_line > max_size_of_line) {size_of_line=max_size_of_line;}

          if (line_counter <= amount_of_lines)
          { 
            for (int pos_copy=0;  pos_copy < size_of_line;  pos_copy++)
            {
            chstr_line[line_counter*max_size_of_line_with_lineend+pos_copy]=text[pos_line_start+pos_copy];
            //  chstr_line[line_counter*max_size_of_line_with_lineend+size_of_line]='\0';      //  already is
            }
          }
        }
        //  writing ^
        //  ===

      line_counter++;
      size_of_line=0;
      }
      else
      {

        if (is_last_step_line_finished == true) 
        {
        is_last_step_line_finished=false;  
        pos_line_start=pos;
        }

      size_of_line++;
      }     

    }      //  for (int pos=0;  pos < text_size;  pos++)

    if ((text_size > 0) && (is_last_step_line_finished == false))
    {
    is_last_step_line_finished=true;

      //  ===
      //  writing V
      if (size_of_line > 0)
      {
        if (size_of_line > max_size_of_line) {size_of_line=max_size_of_line;}

        if (line_counter <= amount_of_lines)
        { 
          for (int pos_copy=0;  pos_copy < size_of_line;  pos_copy++)
          {
          chstr_line[line_counter*max_size_of_line_with_lineend+pos_copy]=text[pos_line_start+pos_copy];
          }
        //  chstr_line[line_counter*max_size_of_line_with_lineend+size_of_line]='\0';      //  already is
        }
      }
        //  writing ^
        //  ===

    line_counter++;
    size_of_line=0;
    }

  }
  else
  {
  bool is_not_emp_symb_in_line=false;

    for (int pos=0;  pos < text_size;  pos++)
    {

      if (text[pos] == '\n')
      {
      is_last_step_line_finished=true;

        if (is_not_emp_symb_in_line == true) 
        {

          //  ===
          //  writing V
          if (line_counter <= amount_of_lines)
          { 
            if (size_of_line > 0) 
            {
              if (size_of_line > max_size_of_line) {size_of_line=max_size_of_line;}

              for (int pos_copy=0;  pos_copy < size_of_line;  pos_copy++)
              {chstr_line[line_counter*max_size_of_line_with_lineend+pos_copy]=text[pos_line_start+pos_copy];}
            }
          }
          //  writing ^
          //  ===

        line_counter++;
        }

      size_of_line=0;  is_not_emp_symb_in_line=false;
      }
      else
      {

        if (is_last_step_line_finished == true) 
        {
        is_last_step_line_finished=false;  
        pos_line_start=pos;
        }

        if ((text[pos] != '\n') && (text[pos] != '\t') && (text[pos] != ' ')) {is_not_emp_symb_in_line=true;}

      size_of_line++;
      }     

    }      //  for (int pos=0;  pos < text_size;  pos++)

    if ((text_size > 0) && (is_last_step_line_finished == false))
    {
    is_last_step_line_finished=true;

      if (is_not_emp_symb_in_line == true) 
      {

        //  ===
        //  writing V
        if (line_counter <= amount_of_lines)
        { 
          if (size_of_line > 0) 
          {
            if (size_of_line > max_size_of_line) {size_of_line=max_size_of_line;}

            for (int pos_copy=0;  pos_copy < size_of_line;  pos_copy++)
            {chstr_line[line_counter*max_size_of_line_with_lineend+pos_copy]=text[pos_line_start+pos_copy];}
          }
        //  chstr_line[line_counter*max_size_of_line_with_lineend+size_of_line]='\0';      //  already is
        }
        //  writing ^
        //  ===

      line_counter++;
      }

    size_of_line=0;  is_not_emp_symb_in_line=false;
    }

  }

  

line_size=max_size_of_line_with_lineend;

  if (text) {delete[] text;  text=0;}


return line_counter;

}




unsigned long long int extract_line_from_file(char *filename, char *& chstr_line, int desired_line_number)
{

  if (filename == 0) {return 0;}
  if (chstr_line != 0) {return 0;}

FILE *file=0;
file=fopen(filename, "rt");

  if (file == 0) {return 0;}

char ch=0;
int line_size=0;
int line_number_counter=0;
int start_pos_of_line=0;
bool is_line_found=false;
bool line_just_finished=false;

  while ((feof(file) == false) && (is_line_found == false))
  {

    if (line_just_finished == true) {start_pos_of_line+=line_size+1;  line_size=0;  line_just_finished=false;}

  fread(&ch, 1, 1, file);  
  
    if (ch == '\n') 
    {
    line_number_counter++;

      if (line_number_counter > desired_line_number) {is_line_found=true;}

    line_just_finished=true;
    } 
    else 
    {
    line_size++;
    } 
    

  }

  if ((is_line_found == false) && (line_number_counter-1 == desired_line_number)) {is_line_found=true;}
 
  if (line_size == 0) {return 0;}


fseek(file, start_pos_of_line, SEEK_SET);
chstr_line=new char[line_size+1];
chstr_line[0]='\0';
fread(chstr_line, line_size, 1, file);
chstr_line[line_size]='\0';

fclose(file);

return line_size;

}



int write_file_to_int_array(char * filename,  int *& int_array,  int linenumber_we_start_with)
{

  if (!filename) {return 0;}
  if (int_array) {return 0;}
  if (linenumber_we_start_with < 0) {linenumber_we_start_with=0;}

int text_size=0;
char * text=0; 
text_size=write_file_to_chstr(filename, text);

  if (text_size < 1) {return 0;}
  if (text == 0) {return 0;}

text_size=strlen(text);


//  count the amount of nums

int num_count=0;
int line_number=0;
int start_pos_num=-1;
char sign='n'; 

  for (int pos=0;  pos < text_size; pos++)
  {

    if (line_number >= linenumber_we_start_with)
    {

      if (((text[pos] >= 48) && (text[pos] <= 57)) || (text[pos] == '-') || (text[pos] == '+')) 
      {
        if ((text[pos] == '-') || (text[pos] == '+')) {sign=text[pos];  start_pos_num=-1;} else {start_pos_num=pos;}
      }
      else
      {

        if (start_pos_num != -1) {num_count++;}

      start_pos_num=-1;  sign='n';
      }

    }
    else
    {
      if (text[pos] == '\n') {line_number++;}
    }
    
  } 

  if (start_pos_num != -1) {num_count++;} 

  if (num_count < 1) {return 0;}


//  allocate memory

int_array=new int[num_count];


//  read  nums

num_count=0;
line_number=0;
start_pos_num=-1;
sign='n';
int last_pos_num=-1;
int num_counter=0;
char chstr_num[128];  chstr_num[0]='\0';

  for (int pos=0;  pos < text_size;  pos++)
  {
  last_pos_num=-1;

    if (line_number >= linenumber_we_start_with)
    {

      if (((text[pos] >= 48) && (text[pos] <= 57)) || (text[pos] == '-') || (text[pos] == '+')) 
      {

        if ((text[pos] == '-') || (text[pos] == '+')) {sign=text[pos];  start_pos_num=-1;} 
        else 
        {
          if (start_pos_num == -1) {start_pos_num=pos;}  
          if (pos+1 == text_size) {last_pos_num=pos;}          
        }

      }
      else
      {

        if (start_pos_num != -1) {last_pos_num=pos;} else {start_pos_num=-1;  sign='n';}

      }

    }
    else
    {   
      if (text[pos] == '\n') {line_number++;}
    }

    
    if ((last_pos_num > -1) && ((last_pos_num - start_pos_num) < 11))
    {
      
      for (int i=start_pos_num;  i <= last_pos_num;  i++)
      {chstr_num[i-start_pos_num]=text[i];}

    chstr_num[last_pos_num-start_pos_num+1]='\0';
    int_array[num_counter]=atoi(chstr_num);  
      if (sign == '-') {int_array[num_counter]=-int_array[num_counter];}
    num_counter++;
    chstr_num[0]='\0';
    start_pos_num=-1;  sign='n';
    }

  }
    
  if (text) {delete[] text;  text=0;}

return num_counter;

}





unsigned char ChULLIntToBinChStr(unsigned long long int num, char * const chstr, bool fillzero, bool littleendian)
{

  if (chstr == 0) {return 0;}


char curpos=0;

unsigned long long int comparedivisor=(ULLONG_MAX-1)/2+1;


  if (fillzero == true)
  {
    for (int i=64-1; i > 0; i--)
    {
      if (num >= comparedivisor) {chstr[curpos]='1'; num-=comparedivisor;} else {chstr[curpos]='0';}
    comparedivisor=(comparedivisor >> 1);  curpos++;
    }  

    if (num == 1) {chstr[curpos]='1'; curpos++;} else {chstr[curpos]='0'; curpos++;}
  }
  else
  {
  bool wasone=false;

    for (int i=64-1; i > 0; i--)
    {
      if (num >= comparedivisor) {chstr[curpos]='1';  wasone=true;  num-=comparedivisor; curpos++;} 
      else 
      {  if (wasone == true) {chstr[curpos]='0';  curpos++;}  }
    comparedivisor=(comparedivisor >> 1); 
    }  

    if (num == 1) {chstr[curpos]='1';  curpos++;} else {chstr[curpos]='0';  curpos++;}
  
  }


  if (littleendian == false)
  {
  char ch=0;
  char lendiv2=curpos/2;
    for (int i=0; i < lendiv2; i++)
    {ch=chstr[i]; chstr[i]=chstr[curpos-i-1]; chstr[curpos-i-1]=ch;}
  }


chstr[curpos]='\0';

return curpos;

}



unsigned char ChUIntToBinChStr(unsigned int num, char * const chstr, bool fillzero, bool littleendian)
{

  if (chstr == 0) {return 0;}


char curpos=0;

unsigned int comparedivisor=(UINT_MAX-1)/2+1;


  if (fillzero == true)
  {
    for (int i=32-1; i > 0; i--)
    {
      if (num >= comparedivisor) {chstr[curpos]='1'; num-=comparedivisor;} else {chstr[curpos]='0';}
    comparedivisor=(comparedivisor >> 1);  curpos++;
    }  

    if (num == 1) {chstr[curpos]='1'; curpos++;} else {chstr[curpos]='0'; curpos++;}
  }
  else
  {
  bool wasone=false;

    for (int i=32-1; i > 0; i--)
    {
      if (num >= comparedivisor) {chstr[curpos]='1';  wasone=true;  num-=comparedivisor; curpos++;} 
      else 
      {  if (wasone == true) {chstr[curpos]='0';  curpos++;}  }
    comparedivisor=(comparedivisor >> 1); 
    }  

    if (num == 1) {chstr[curpos]='1';  curpos++;} else {chstr[curpos]='0';  curpos++;}
  
  }



  if (littleendian == false)
  {
  char ch=0;
  char lendiv2=curpos/2;
    for (int i=0; i < lendiv2; i++)
    {ch=chstr[i]; chstr[i]=chstr[curpos-i-1]; chstr[curpos-i-1]=ch;}
  }


chstr[curpos]='\0';

return curpos;

}



unsigned char ChUSIntToBinChStr(unsigned short int num, char * const chstr, bool fillzero, bool littleendian)
{

  if (chstr == 0) {return 0;}


char curpos=0;

unsigned int comparedivisor=(USHRT_MAX-1)/2+1;


  if (fillzero == true)
  {
    for (int i=16-1; i > 0; i--)
    {
      if (num >= comparedivisor) {chstr[curpos]='1'; num-=comparedivisor;} else {chstr[curpos]='0';}
    comparedivisor=(comparedivisor >> 1);  curpos++;
    }  

    if (num == 1) {chstr[curpos]='1'; curpos++;} else {chstr[curpos]='0'; curpos++;}
  }
  else
  {
  bool wasone=false;

    for (int i=16-1; i > 0; i--)
    {
      if (num >= comparedivisor) {chstr[curpos]='1';  wasone=true;  num-=comparedivisor; curpos++;} 
      else 
      {  if (wasone == true) {chstr[curpos]='0';  curpos++;}  }
    comparedivisor=(comparedivisor >> 1); 
    }  

    if (num == 1) {chstr[curpos]='1';  curpos++;} else {chstr[curpos]='0';  curpos++;}
  
  }


  if (littleendian == false)
  {
  char ch=0;
  char lendiv2=curpos/2;
    for (int i=0; i < lendiv2; i++)
    {ch=chstr[i]; chstr[i]=chstr[curpos-i-1]; chstr[curpos-i-1]=ch;}
  }


chstr[curpos]='\0';

return curpos;

}



unsigned char ChUCharToBinChStr(unsigned char num, char * const chstr, bool fillzero, bool littleendian)
{

  if (chstr == 0) {return 0;}


char curpos=0;

unsigned int comparedivisor=(UCHAR_MAX-1)/2+1;


  if (fillzero == true)
  {
    for (int i=8-1; i > 0; i--)
    {
      if (num >= comparedivisor) {chstr[curpos]='1'; num-=comparedivisor;} else {chstr[curpos]='0';}
    comparedivisor=(comparedivisor >> 1);  curpos++;
    }  

    if (num == 1) {chstr[curpos]='1'; curpos++;} else {chstr[curpos]='0'; curpos++;}
  }
  else
  {
  bool wasone=false;

    for (int i=8-1; i > 0; i--)
    {
      if (num >= comparedivisor) {chstr[curpos]='1';  wasone=true;  num-=comparedivisor; curpos++;} 
      else 
      {  if (wasone == true) {chstr[curpos]='0';  curpos++;}  }
    comparedivisor=(comparedivisor >> 1); 
    }  

    if (num == 1) {chstr[curpos]='1';  curpos++;} else {chstr[curpos]='0';  curpos++;}
  
  }


  if (littleendian == false)
  {
  char ch=0;
  char lendiv2=curpos/2;
    for (int i=0; i < lendiv2; i++)
    {ch=chstr[i]; chstr[i]=chstr[curpos-i-1]; chstr[curpos-i-1]=ch;}
  }


chstr[curpos]='\0';

return curpos;

}



void VDisplayULLIntToBin(unsigned long long int num, bool fillzero, bool littleendian)
{

char chstr[128];   chstr[0]='\0';
ChULLIntToBinChStr(num, chstr, fillzero, littleendian);
cout<<chstr;

}



void VDisplayUIntToBin(unsigned int num, bool fillzero, bool littleendian)
{

char chstr[128];   chstr[0]='\0';
ChUIntToBinChStr(num, chstr, fillzero, littleendian);
cout<<chstr;

}



void VDisplayUSIntToBin(unsigned short int num, bool fillzero, bool littleendian)
{

char chstr[128];   chstr[0]='\0';
ChUSIntToBinChStr(num, chstr, fillzero, littleendian);
cout<<chstr;

}



void VDisplayUCharToBin(unsigned char num, bool fillzero, bool littleendian)
{

char chstr[128];   chstr[0]='\0';
ChUCharToBinChStr(num, chstr, fillzero, littleendian);
cout<<chstr;

}



//  for debug: representation in bin format to screen
void VDisplayFixLenULLIntToBin(unsigned long long int num, int len, bool fillzero, bool littleendian)
{

  if (len > 64) {len=64;}

  if (len > 0)
  {
  unsigned long long int one=1;
  int shifted_one=0;

    for (int i=0;  i < len;  i++) 
    {
    int shifted_one=(one << (len-i-1));  

      if (littleendian == false) {shifted_one=(one << i);}
    
      if ((shifted_one & num) != 0)  {cout<<"1 ";} else  {cout<<"0 ";} 

    }

  }

}



bool is_ch_ref_to_dbl_num(char ch)
{

  if (((ch > 47) && (ch < 58)) || (ch == '.') || (ch == ',') || (ch == '-') || (ch == '+') || (ch == 'e') || (ch == 'E'))
  {return true;}

return false;

}



bool is_int_or_llint_num(char *chstrnum)
{

  if (chstrnum == 0) {return false;}

int len=strlen(chstrnum);
int signsymbcount=0;
int signsymbpos=0;
int cyfsymbcount=0;
bool isbadsymbol=false;

  if (len > 32) {return false;}
  if (len < 1 ) {return false;}

  for (int i=0; i < len; i++)
  {
    if ((chstrnum[i] <= 47) || (chstrnum[i] >= 58))
    {
      if ((chstrnum[i] == '-') || (chstrnum[i] == '+'))
      {
      signsymbcount++;  signsymbpos=i;

        if (signsymbcount > 1) {return false;}

      }
      else
      {
      return false;
      }
    }
    else
    {cyfsymbcount++;}
  } 


  if (signsymbcount > 1) {return false;}

  if (len > 1)    {return false;}
  if ((signsymbcount > 0) && (signsymbpos != 0)) {return false;}

return true;

}



bool is_fl_or_dbl_num(char *chstrnum)
{

  if (chstrnum == 0) {return false;}

int len=strlen(chstrnum);
int sepsymbcount=0, esymbcount=0, sign1symbcount=0, sign2symbcount=0;
int sepsymbpos=0, esymbpos=0, sign1symbpos=0, sign2symbpos=0;
int cyfsymbcount=0;
bool isbadsymbol=false;

  if (len > 32) {return false;}
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
        if ((chstrnum[i] == '.') || (chstrnum[i] == ','))
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



void add_chstr_to_string(string * string_ptr_dest,  char *chstr_added_str,  int given_len)
{

int start_pos=string_ptr_dest->length();
int added_len=0;
  if (chstr_added_str) {added_len=strlen(chstr_added_str);}

char chstr[2];  chstr[0]=' ';  chstr[1]='\0';

  if (string_ptr_dest != 0)
  {
  string_ptr_dest->append(chstr_added_str);

    for (int i=0; i < given_len-added_len;  i++)
    {string_ptr_dest->append(chstr);} 
  }

}




void add_string_to_string(string * string_ptr_dest,  string * string_ptr_src,  int given_len)
{

int start_pos=string_ptr_dest->length();
int added_len=string_ptr_src->length();
char chstr[2];  chstr[0]=' ';  chstr[1]='\0';

  if (string_ptr_dest != 0)
  {
  string_ptr_dest->append(string_ptr_src->c_str());

    for (int i=0; i < given_len-added_len;  i++)
    {string_ptr_dest->append(chstr);}
  }

}



void add_int_to_string(string * string_ptr, int num)
{

char chstr_num[64];  chstr_num[0]='\0';
sprintf(chstr_num, "%d", num);
string_ptr->append(chstr_num);

}



void add_uint_to_string(string * string_ptr, unsigned int num)
{

char chstr_num[64];  chstr_num[0]='\0';
sprintf(chstr_num, "%u", num);
string_ptr->append(chstr_num);

}



void add_llint_to_string(string * string_ptr, long long int num)
{

char chstr_num[64];  chstr_num[0]='\0';
sprintf(chstr_num, "%lld", num);
string_ptr->append(chstr_num);

}



void add_ullint_to_string(string * string_ptr, unsigned long long int num)
{

char chstr_num[64];  chstr_num[0]='\0';
sprintf(chstr_num, "%llu", num);
string_ptr->append(chstr_num);

}



void add_double_to_string(string * string_ptr, double num)
{

char chstr_num[64];  chstr_num[0]='\0';
sprintf(chstr_num, "%E", num);
string_ptr->append(chstr_num);

}



void add_int_to_chstr(char * chstr, int num)
{

char chstr_num[64];  chstr_num[0]='\0';
sprintf(chstr_num, "%d", num);
strcat(chstr, chstr_num);

}



void add_ullint_to_chstr(char * chstr, unsigned long long int num)
{

char chstr_num[64];  chstr_num[0]='\0';
sprintf(chstr_num, "%llu", num);
strcat(chstr, chstr_num);

}




void display_mem_size(long long int mem_size_bit, bool is_line_feed)
{

long long int mem_size_Gb=mem_size_bit/8/1024/1024/1024;
mem_size_bit-=mem_size_Gb*1024*1024*1024*8;
long long int mem_size_Mb=mem_size_bit/8/1024/1024;
mem_size_bit-=mem_size_Mb*1024*1024*8;
long long int mem_size_kb=mem_size_bit/8/1024;
mem_size_bit-=mem_size_kb*1024*8;
long long int mem_size_b =mem_size_bit/8;
cout<<mem_size_Mb<<" Mb + "<<mem_size_kb<<" kb + "<<mem_size_b<<" b";

  if (is_line_feed == true) {cout<<endl;}

}




void  display_to_screen_ullint_fix_len(unsigned long long int x,  unsigned int fix_len)
{

char chstr_x[64];
chstr_x[0]='\0';
sprintf(chstr_x,  "% *lli", fix_len, x);
cout<<chstr_x;

}




void  display_to_screen_llint_fix_len(long long int x,  unsigned int fix_len)
{

char chstr_x[64];
chstr_x[0]='\0';
sprintf(chstr_x,  "% *lli", fix_len, x);
cout<<chstr_x;

}







inline unsigned long long int treat_out_of_ullint_2_num_befor_add(unsigned long long int x1,  unsigned long long int x2,  bool *is_out_of_ullint)
{

unsigned long long int rest=((unsigned long long int) __MACR_ULLINT_MAX_VALUE)-x1;

  if (rest >= x2) 
  {  
  //  if (is_out_of_ullint) {*is_out_of_ullint=false;}  
  return x2;
  }

  if (is_out_of_ullint) {*is_out_of_ullint=true;}
  
return rest;  
    
}






#endif

