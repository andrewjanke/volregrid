/* lex.l - scanner for volregrid config files */

DIGIT     [0-9]
LIST_ID   [A-Z]+(_[A-Z]+)*
REAL      "-"?{DIGIT}+(\.{DIGIT}*)?([eE][-+]?{DIGIT}+)?
INT       "-"?{DIGIT}+

%{ 
   #include <math.h>
   #include <string.h>
   #include "arb_path_io.h"
   
   /* function prototypes */
   int get_coord_from_str(char *str, double *x, double *y, double *z);
   int get_num_reps_from_str(char *str);
   
   /* global parsing vars */
   int yylineno;
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

<*>#.*\n               /* gobble up comments  */
<*>[ \t\n]+            /* eat up whitespace   */

 /* co-ordinates */
{REAL},{REAL},{REAL}         { printf("***| Coord  ");
                               if(get_coord_from_str(yytext, &x, &y, &z) && 
                                  add_coord_to_list(0, x, y, z)){
                                  printf("  -- %g %g %g\n", x, y, z);
                                  return 1;
                                  }
                               }
@{REAL},{REAL},{REAL}        { printf("***| RCoor  ");
                               yytext++;
                               if(get_coord_from_str(yytext, &x, &y, &z) && 
                                  add_rcoord_to_list(0, x, y, z)){
                                  printf("  -- %g %g %g\n", x, y, z);
                                  return 1;
                                  }
                               }
@{REAL},{REAL},{REAL}:{REAL} { printf("***| RRCoo   %s\n", yytext);
                               yytext++;
                               if(get_coord_from_str(yytext, &x, &y, &z)){
                                  
                                  reps = get_num_reps_from_str(yytext);
                                  add_rrcoord_to_list(0, x, y, z, reps);
                                  printf("  -- %g %g %g\n", x, y, z);
                                  return 1;
                                  }
                               }

 /* parameters */
DIMENSIONS\:      |
SLICES\:          printf("***| INT KEYWORD        %s\n", yytext);
SLICES_START\:    |
SLICES_STEP\:     printf("***| REAL KEYWORD       %s\n", yytext);
SLICES_DIM\:      printf("***| STRING KEYWORD     %s\n", yytext);

 /* Matrix */
"M = ["         { mat_pos = 0;
                  BEGIN(MATRIX);
                  }
<MATRIX>{
   {REAL}       { MAT[mat_pos++] = atof(yytext); }
   ]            { if(mat_pos != 12){
                     fprintf(stdout, "Less than 12 elements found in Matrix!\n");
                     return 0;
                     }
                  set_curr_matrix(MAT);
                  BEGIN(INITIAL);
                  }
   }

 /* list definition */
LIST                   { BEGIN(LIST);                        }
<LIST>{
   {LIST_ID}           { printf("\n***| FOUND LIST:  %s\n", yytext);    
                         c_list_id = new_coord_list(0, yytext);
                         printf("\n***| GOT ID:  %d\n", c_list_id); 
                         }
   \{                  
   {REAL},{REAL},{REAL}          { printf("***| LCoord  %d %s\n", c_list_id, yytext);
                                   if(!get_coord_from_str(yytext, &x, &y, &z) || 
                                      !add_coord_to_list(c_list_id, x, y, z)){
                                      printf("  -- %g %g %g\n", x, y, z);
                                      return 0;
                                      }
                                   }
   @{REAL},{REAL},{REAL}         { printf("***| LRCoor   %s\n", yytext);
                                   yytext++;
                                   if(!get_coord_from_str(yytext, &x, &y, &z) || 
                                      !add_rcoord_to_list(c_list_id, x, y, z)){
                                      printf("  -- %g %g %g\n", x, y, z);
                                      return 0;
                                      }
                                   }
   @{REAL},{REAL},{REAL}:{REAL}  { printf("***| LRRCoo   %s\n", yytext);
                                   yytext++;
                                   if(get_coord_from_str(yytext, &x, &y, &z)){

                                      reps = get_num_reps_from_str(yytext);
                                      add_rrcoord_to_list(c_list_id, x, y, z, reps);
                                      printf("  -- %g %g %g\n", x, y, z);
                                      return 1;
                                      }
                                   }
   \}                  { printf("***|   END LIST\n\n");
                         BEGIN(INITIAL);
                         }
   }

 /* list call */
CALL                   { BEGIN(CALL);                           }
<CALL>{LIST_ID}        { printf("***| Calling list %s\n", yytext); 
                         id = get_list_id_from_name(yytext);
                         
                         printf("Got list ID: %d (%s)\n", id, yytext); 
                         if(id != 0 && call_list(id)){
                            BEGIN(INITIAL);
                            return 1;
                            }
                         else{
                            printf("Couldn't find list %s\n", yytext);
                            }
                         
                         BEGIN(INITIAL);
                         }

 /* default case */
<*>.                   { printf( "***| Unrecognized character: %s\n", yytext );
                         }

%%

void yyerror(char *s){
   fprintf(stderr, "line %d: %s\n", yylineno, s);
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

int get_num_reps_from_str(char *str){
   
   /* get to the number in question */
   while(str[0] != ':'){
      str++;
      }
   str++;
   
   return atof(str);
   }

int yyflex_init(char *fname){
   yyin = fopen(fname, "r");
   
   return (yyin == NULL) ? 0 : 1;
   }

int yyflex_end(void){
   fclose(yyin);
   return 1;
   }
