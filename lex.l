/* lex.l - scanner for volregrid config files */

DIGIT     [0-9]
LIST_ID   [A-Z]+(_[A-Z]+)*
REAL      "-"?{DIGIT}+(\.{DIGIT}*)?([eE][-+]?{DIGIT}+)?
INT       "-"?{DIGIT}+
COORD     {REAL}\,{REAL}\,{REAL}
RCOORD    "@["{COORD}"]"

%{ 
   #include <math.h>
   #include <string.h>
   #include "arb_path_io.h"
   
   /* function prototypes */
   int get_coord_from_str(char *str, double *x, double *y, double *z);
   int get_reps_from_str(char *str);
   void yyerror(void);
   
   /* global parsing vars */
   int line_num = 1;
   int c_list_id = 0;
   int id;
   char list_name[256];
   double x,y,z,reps;
   
   int    mat_pos = 0;
   double MAT[12] = {1,0,0,0,
                     0,1,0,0,
                     0,0,1,0};
   
%}

/* States */
%x LIST
%x MATRIX
%x PARAM
%x CALL
%x INT_KWD
%x REAL_KWD
%x STR_KWD

%%

<*>#.*                 /* gobble up comments  */
<*>[ \t]+              /* eat up whitespace   */
<*>\n                  { line_num++; }  /* count lines */

 /* co-ordinates */
{COORD}                { if(get_coord_from_str(yytext, &x, &y, &z) && 
                            add_coord_to_list(0, x, y, z)){
                            return 1;                                
                            }                                        
                         else{                                       
                            yyerror();                               
                            }                                        
                         }                                           
{RCOORD}               { yytext+=2;                                         
                         if(get_coord_from_str(yytext, &x, &y, &z) &&       
                            add_rcoord_to_list(0, x, y, z, 1)){
                            return 1;                                       
                            }                                               
                         else{                                              
                            yyerror();                                      
                            }                                               
                         }                                                  
{RCOORD}\*{REAL}       { yytext+=2;                                                          
                         if(get_coord_from_str(yytext, &x, &y, &z) &&                        
                            add_rcoord_to_list(0, x, y, z, get_reps_from_str(yytext))){
                            return 1;                                                        
                            }                                                                
                         else{                                                               
                            yyerror();                                                       
                            }                                                                
                         }                                                                   

 /* parameters */
DIMENSIONS\:           |
SLICES\:               printf("***| INT KEYWORD        %s\n", yytext);
SLICES_START\:         |
SLICES_STEP\:          printf("***| REAL KEYWORD       %s\n", yytext);
SLICES_DIM\:           printf("***| STRING KEYWORD     %s\n", yytext);

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
                         printf("***| LIST:  %s\n", yytext);    
                         c_list_id = new_coord_list(0, yytext);
                         if(c_list_id == 0){
                            fprintf(stderr, "List %s cannot be redefined\n", yytext);
                            yyerror();
                            }
                        
                         BEGIN(LIST);
                         }
<LIST>{
   \{                  
   {COORD}             { printf("***| LCoord   %d", c_list_id);
                         if(get_coord_from_str(yytext, &x, &y, &z) &&           
                            add_coord_to_list(c_list_id, x, y, z)){             
                            printf("  -- %g %g %g\n", x, y, z);                 
                            }                                                   
                          else{                                                 
                            yyerror();                                          
                            }                                                   
                         }                                                      
   {RCOORD}            { printf("***| LRCoord  %d", c_list_id);
                         yytext+=2;                                              
                         if(get_coord_from_str(yytext, &x, &y, &z) &&            
                            add_rcoord_to_list(c_list_id, x, y, z, 1)){          
                            printf("  -- %g %g %g\n", x, y, z);                  
                            }                                                    
                          else{                                                  
                            yyerror();                                           
                            }                                                    
                         }                                                       
   {RCOORD}\*{REAL}    { printf("***| LRRCoord %d", c_list_id);
                         yytext+=2;                                                                      
                         if(get_coord_from_str(yytext, &x, &y, &z) &&                                    
                            add_rcoord_to_list(c_list_id, x, y, z, get_reps_from_str(yytext))){          
                            printf("  -- %g %g %g\n", x, y, z);                                          
                            }                                                                            
                          else{                                                                          
                            yyerror();                                                                   
                            }                                                                            
                         }                                                                               
   \}                  { printf("***|   END LIST\n");
                         BEGIN(INITIAL);
                         }
   }

 /* list call */
CALL\ {LIST_ID}        { yytext += 5;
                         printf("***| Calling list %s\n", yytext);
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
   fprintf(stderr, "Error at line %d: %s\n", line_num, yytext);
   exit(EXIT_FAILURE);
   }
 
int yywrap(){  
   return 1;
   }  

int get_coord_from_str(char *str, double *x, double *y, double *z){
   char *ptr;
   
   /* get the first double */
   *x = strtod(str, &ptr);
   
   /* skip a "," */
   if(ptr[0] != ','){
      return 0;
      }
   ptr++;
   
   *y = strtod(ptr, &ptr);
   
   if(ptr[0] != ','){
      return 0;
      }
   ptr++;
   
   *z = strtod(ptr, &ptr);
   
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
