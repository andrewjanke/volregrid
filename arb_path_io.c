/* arb_path_io.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arb_path_io.h"

/* global vars for arb_path_regrid */
/* list 0 is a special list, (the working buffer one) */
int      num_lists;
Coord_list *coord_lists;

/* current affine matrix */
int      aff_M_identity = 1;
double   aff_M[3][4] = { {1, 0, 0, 0},
{0, 1, 0, 0},
{0, 0, 1, 0}
};

/* function prototypes */

/* use realloc to set the size of the point data */
int      set_coord_list_size(Coord_list list, size_t size);

/* lex bits */
int      yyflex_init(char *fname);
int      yyflex_end(void);
int      yylex(void);

/* init the whole schebang */
int init_arb_path(char *fname)
{
   int id;

   /* initialise the parser with the config file */
   if(!yyflex_init(fname)){
      fprintf(stderr, "Failed to init the parser\n");
      return 0;
      }

   /* init the data store - list[0] */
   num_lists = 0;
   id = new_coord_list(0, "default");
   fprintf(stdout, "Returned list id of %d\n", id);

   return 1;
   }

/* finalize the parser */
void end_arb_path(void){
   
   /* tidy up */
   yyflex_end();
   }

/* get some more co-ordinates, returns an empty list when done */
Coord_list get_some_arb_path_coords(int max_buf_size)
{
   int i, j, c;
   double tmp[3];
   Coord_list list = coord_lists[0];

   fprintf(stdout, "Getting %d coords\n", max_buf_size);

   /* get some co-ordinates */
   c = 0;
   while(list->n_pts < max_buf_size && c < max_buf_size){
      if(!yylex()){
         fprintf(stdout, "Error in parsing, this isn't good!\n");
         exit(EXIT_FAILURE);
         }
      c++;
      }
   
   /* multiply by the current affine matrix */
   if(!aff_M_identity){
      fprintf(stdout, "Multiplying by aff matrix\n");
      
      for(c=0; c<list->n_pts; c++){
         for(i = 0; i < 3; i++){
            tmp[i] = 0;
            for(j = 0; j < 3; j++){
               tmp[i] += aff_M[i][j]* list->pts[c].coord[j];
               }
            }
         
         list->pts[c].coord[0] = tmp[0] + aff_M[0][3];
         list->pts[c].coord[1] = tmp[1] + aff_M[1][3];
         list->pts[c].coord[2] = tmp[2] + aff_M[2][3];
         }
      
      }

   return (list);
   }

/* return a new Coord_list's id */
int new_coord_list(size_t size, char *name)
{
   int new_id;
   
   new_id = num_lists;
   
   /* grow the list of list pointers */
   coord_lists = (Coord_list *) realloc(coord_lists, num_lists * sizeof(*coord_lists));

   /* allocate space for the struct */
   coord_lists[new_id] = (Coord_list) malloc(sizeof(*coord_lists[new_id]));

   /* init the structure */
   coord_lists[new_id]->id = num_lists;
   strcpy(coord_lists[new_id]->name, name);
   coord_lists[new_id]->n_pts = 0;
   coord_lists[new_id]->pts = NULL;

   set_coord_list_size(coord_lists[new_id], size);

   num_lists++;
   return (new_id);
   }

int set_coord_list_size(Coord_list list, size_t size)
{

   list->n_pts = size;
   list->pts = (Coord *) realloc(list->pts, list->n_pts * sizeof(*list->pts));

   return 1;
   }

/* update the current affine matrix */
void set_curr_matrix(double mat[12])
{
   int      i, j, c;

   aff_M_identity = 1;
   c = 0;
   for(i = 0; i < 3; i++){
      for(j = 0; j < 4; j++){
         aff_M[i][j] = mat[c];

         /* check for "identity-ness" */
         if(i == j){
            if(mat[c] != 1.0){
               aff_M_identity = 0;
               }
            }
         else{
            if(mat[c] != 0.0){
               aff_M_identity = 0;
               }
            }

         c++;
         }
      }

   print_curr_matrix();
   }

int add_coord_to_list(int list_id, double x, double y, double z){
   Coord_list list = coord_lists[list_id];
   int c_pt = list->n_pts;
   
   set_coord_list_size(list, c_pt+1);
   
   list->pts[c_pt].coord[0] = x;
   list->pts[c_pt].coord[1] = y;
   list->pts[c_pt].coord[2] = z;
   
   return 1;
   }

int add_rcoord_to_list(int list_id, double x, double y, double z){
   Coord_list list = coord_lists[list_id];
   int c_pt = list->n_pts;
   
   set_coord_list_size(list, c_pt+1);
   
   list->pts[c_pt].coord[0] = list->pts[c_pt-1].coord[0] + x;
   list->pts[c_pt].coord[1] = list->pts[c_pt-1].coord[1] + y;
   list->pts[c_pt].coord[2] = list->pts[c_pt-1].coord[2] + z;
   
   return 1;
   }

int add_rrcoord_to_list(int list_id, double x, double y, double z, double rep){
   Coord_list list = coord_lists[list_id];
   int c_pt = list->n_pts;
   int c;
   
   set_coord_list_size(list, c_pt+rep);
   
   for(c=0; c<rep; c++){
      list->pts[c_pt+c].coord[0] = list->pts[c_pt+c-1].coord[0] + x;
      list->pts[c_pt+c].coord[1] = list->pts[c_pt+c-1].coord[1] + y;
      list->pts[c_pt+c].coord[2] = list->pts[c_pt+c-1].coord[2] + z;
      }
   
   return 1;
   }

int get_list_id_from_name(char *name){
   int c;
   
   for(c=0; c<num_lists; c++){
      if(strcmp(coord_lists[c]->name, name) == 0){
         return coord_lists[c]->id;
         }
      }
   
   return 0;
   }


int call_list(int list_id){
   int c;
   Coord_list list = coord_lists[list_id];
   
   /* copy the co-ords into the main (0) list */
   for(c=0; c<list->n_pts; c++){
      add_coord_to_list(0, list->pts[c].coord[0],
                           list->pts[c].coord[1],
                           list->pts[c].coord[2]);
      }
   
   return 1;
   }

/* debugging function to print a coord list */
void print_coord_list(Coord_list list)
{
   int      i;

   fprintf(stdout, "List[%d]: %s (%d pts)\n", list->id, list->name, list->n_pts);
   for(i = 0; i < list->n_pts; i++){
      fprintf(stdout, "| [%d] %g %g %g\n", i, list->pts[i].coord[0],
              list->pts[i].coord[1], list->pts[i].coord[2]);
      }
   }

void print_curr_matrix(void)
{
   int      i;

   fprintf(stdout, "Curr Matrix - %sidentity\n",
           (aff_M_identity) ? "" : "non-");
   for(i = 0; i < 3; i++){
      fprintf(stdout, "  %6.2f %6.2f %6.2f %6.2f\n",
              aff_M[i][0], aff_M[i][1], aff_M[i][2], aff_M[i][3]);
      }
   }
