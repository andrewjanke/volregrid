/* lex.l - scanner for volregrid config files */

DIGIT     [0-9]
LIST_ID   [A-Z]+(_[A-Z]+)*
REAL      "-"?{DIGIT}+(\.{DIGIT}*)?([eE][-+]?{DIGIT}+)?
INT       "-"?{DIGIT}+
COORD     {REAL}\,{REAL}\,{REAL}
RCOORD    "@["{COORD}"]"

%{

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "arb_path_io.h"

/* function prototypes */
int      get_coord_from_str(char *str);
int      get_reps_from_str(char *str);
void     yyerror(void);

/* global parsing vars */
int      line_num = 1;
int      c_list_id = 0;
int      id;
char     list_name[256];
double   x, y, z, reps;

int      mat_pos = 0;
double   MAT[12] = { 1, 0, 0, 0,
   0, 1, 0, 0,
   0, 0, 1, 0
   };
 
%}

/* States */
%x LIST
%x MATRIX

%%

<*>#.*                 /* gobble up comments  */
<*>[ \t]+              /* eat up whitespace   */
<*>\n                  { line_num++; }  /* count lines */

 /* co-ordinates */
{COORD}                { if(get_coord_from_str(yytext) && 
                            add_coord_to_list(0, x, y, z)){
                            return 1;                                
                            }                                        
                         else{                                       
                            yyerror();                               
                            }                                        
                         }                                           
{RCOORD}               { yytext+=2;                                         
                         if(get_coord_from_str(yytext) &&       
                            add_rcoord_to_list(0, x, y, z, 1)){
                            return 1;                                       
                            }                                               
                         else{                                              
                            yyerror();                                      
                            }                                               
                         }                                                  
{RCOORD}\*{REAL}       { yytext+=2;                                                          
                         if(get_coord_from_str(yytext) &&                        
                            add_rcoord_to_list(0, x, y, z, get_reps_from_str(yytext))){
                            return 1;                                                        
                            }                                                                
                         else{                                                               
                            yyerror();                                                       
                            }                                                                
                         }                                                                   

 /* Matrix */
"M = ["                { mat_pos = 0;
                         BEGIN(MATRIX);
                         }
<MATRIX>{
   {REAL}              { MAT[mat_pos++] = atof(yytext); }
   ]                   { if(mat_pos != 12){
                            fprintf(stdout, "Less than 12 elements found in Matrix!\n");
                            yyerror();
                            }
                         set_curr_matrix(MAT);
                         BEGIN(INITIAL);
                         }
   }

 /* list definition */
LIST\ {LIST_ID}        { yytext += 5;   
                         c_list_id = new_coord_list(0, yytext);
                         if(c_list_id == 0){
                            fprintf(stderr, "List %s cannot be redefined\n", yytext);
                            yyerror();
                            }
                         BEGIN(LIST);
                         }
<LIST>{
   \{                  
   {COORD}             { if(!get_coord_from_str(yytext) ||           
                            !add_coord_to_list(c_list_id, x, y, z)){
                            yyerror();                                          
                            }                                                   
                         }                                                      
   {RCOORD}            { yytext+=2;                                              
                         if(!get_coord_from_str(yytext) ||            
                            !add_rcoord_to_list(c_list_id, x, y, z, 1)){
                            yyerror();                                           
                            }                                                    
                         }                                                       
   {RCOORD}\*{REAL}    { yytext+=2;
                         if(!get_coord_from_str(yytext) ||
                            !add_rcoord_to_list(c_list_id, x, y, z, get_reps_from_str(yytext))){
                            yyerror();                                                                   
                            }                                                                            
                         }                                                                               
   \}                  { BEGIN(INITIAL);
                         }
   }

 /* list call */
CALL\ {LIST_ID}        { yytext += 5;
                         if(call_list(yytext)){
                            return 1;
                            }
                         else{
                            yyerror();
                            }
                         }

 /* default case */
<*>.                   { yyerror(); }

%%

void yyerror(void){
   char tmp[100]; /* Haaaaack! :) */
   int c;
   
   /* slurp in rest of the offending line */
   c = 0;
   while((tmp[c++] = input()) != '\n' && c < 99){
      }
   tmp[c++] = '\0';
   
   /* pretty print out the error */
   fprintf(stderr, "\nError at line %d:\n", line_num);
   fprintf(stderr, "%s%s", yytext, tmp);
   for(c=0; c<strlen(yytext); c++){
      fputc(' ', stderr);
      }
   fprintf(stderr, "^\nparse error\n");
   
   exit(EXIT_FAILURE);
   }
 
int yywrap(){  
   return 1;
   }  

int get_coord_from_str(char *str){
   char *ptr;
   
   /* get the first double */
   x = strtod(str, &ptr);
   
   /* skip a "," */
   if(ptr[0] != ','){
      return 0;
      }
   ptr++;
   
   y = strtod(ptr, &ptr);
   
   if(ptr[0] != ','){
      return 0;
      }
   ptr++;
   
   z = strtod(ptr, &ptr);
   
   return 1;
   }

int get_reps_from_str(char *str){
   
   /* get to the number in question */
   while(str[0] != '*'){
      str++;
      }
   str++;
   
   return atoi(str);
   }

int yyflex_init(char *fname){
   yyin = fopen(fname, "r");
   
   return (yyin == NULL) ? 0 : 1;
   }

int yyflex_end(void){
   fclose(yyin);
   return 1;
   }
