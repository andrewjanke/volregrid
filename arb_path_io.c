/* arb_path_io.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arb_path_io.h"

/* global vars for arb_path_regrid */
/* list 0 is a special list, (the working buffer one) */
int      num_lists = 0;
Coord_list *coord_lists;

/* current affine matrix */
int      aff_M_identity = 1;
double   aff_M[3][4] = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0} };

/* current position in space */
double   c_pos[3];

/* file pointer for data file */
FILE    *data_fp = NULL;

/* function prototypes */
int      set_coord_list_size(Coord_list list, size_t size);
int      yyflex_init(char *fname);
int      yyflex_end(void);
int      yylex(void);


/* init the whole schebang */
int init_arb_path(char *coord_fname, char *data_fname)
{
   /* initialise the parser with the coordinates config file */
   if(!yyflex_init(coord_fname)){
      fprintf(stderr, "Failed to init the parser\n");
      return 0;
      }

   /* init the primary coord list - list[0] */
   if(new_coord_list(0, "DEFAULT") != 0){
      fprintf(stdout, "Something is amiss, the initial list is not 0!\n");
      exit(EXIT_FAILURE);
      }
   
   /* open the data file */
   data_fp = fopen(data_fname, "r");
   if(data_fp == NULL){
      fprintf(stderr, "Failed to open the data file %s\n", data_fname);
      return 0;
      }
    
   return 1;
   }

/* finalize the parser */
void end_arb_path(void)
{

   /* tidy up */
   yyflex_end();
   }

/* get some more co-ordinates, returns an empty list when done */
Coord_list get_some_arb_path_coords(size_t max_buf_size)
{
   Coord_list list = coord_lists[0];

   /* first empty out the list */
   if(!set_coord_list_size(list, 0)){
      fprintf(stderr, "Error clearing out list[%d] %s\n", list->id, list->name);
      }

   /* then fill it up again */
   while (list->n_pts < max_buf_size){
      if(!yylex()){

         /* exit early if we are done */
         fprintf(stdout, "yylex() only got %d coords\n", list->n_pts);
         return (list);
         }
      }

   return (list);
   }

/* get some arb_path data, returns NULL on fail                */
int get_some_arb_path_data(double *data_buf, size_t n_pts, size_t vector_size)
{
   

   return 1;
   }

/* return a new Coord_list's id */
int new_coord_list(size_t size, char *name)
{
   int      new_id;
   int      c;

   new_id = num_lists;

   /* first check if a list with that name exists already */
   for(c = 0; c < num_lists; c++){
      if(strcmp(coord_lists[c]->name, name) == 0){
         return 0;
         }
      }

   /* grow the list of list pointers */
   num_lists++;
   coord_lists = (Coord_list *) realloc(coord_lists, num_lists * sizeof(*coord_lists));

   /* allocate space for the struct */
   coord_lists[new_id] = (Coord_list) malloc(sizeof(*coord_lists[new_id]));

   /* init the structure */
   coord_lists[new_id]->id = new_id;
   strcpy(coord_lists[new_id]->name, name);
   coord_lists[new_id]->pts = NULL;

   coord_lists[new_id]->n_pts = 0;
   coord_lists[new_id]->alloc_size = 0;
   set_coord_list_size(coord_lists[new_id], size);

   return (new_id);
   }

/* use realloc to set the size of the point data */
int set_coord_list_size(Coord_list list, size_t size)
{
   /* only grow a list, never shrink */
   if(size > list->alloc_size){
      list->alloc_size = size;
      list->pts = (Coord *) realloc(list->pts, list->alloc_size * sizeof(*list->pts));
      }

   list->n_pts = size;

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
   }

int add_coord_to_list(int list_id, double x, double y, double z)
{
   Coord_list list = coord_lists[list_id];
   int      c_pt = list->n_pts;

   set_coord_list_size(list, list->n_pts + 1);

   /* only apply the matrix if we have to             */
   /* all lists are referenced to the identity matrix */
   if(list_id > 0 || aff_M_identity){
      list->pts[c_pt].coord[0] = x;
      list->pts[c_pt].coord[1] = y;
      list->pts[c_pt].coord[2] = z;
      }
   else{
      int      i;

      fprintf(stdout, "Multiplying");
      for(i = 0; i < 3; i++){
         list->pts[c_pt].coord[i] =
            (aff_M[i][0] * x) + (aff_M[i][1] * y) + (aff_M[i][2] * z) + aff_M[i][3];
         }
      }

   return 1;
   }

int add_rcoord_to_list(int list_id, double x, double y, double z, double rep)
{
   Coord_list list = coord_lists[list_id];
   int      c_pt = list->n_pts;
   int      c;

   set_coord_list_size(list, c_pt + rep);

   /* multiply by the current affine matrix */
   if(list_id > 0 || aff_M_identity){
      for(c = 0; c < rep; c++){
         list->pts[c_pt + c].coord[0] = list->pts[c_pt + c - 1].coord[0] + x;
         list->pts[c_pt + c].coord[1] = list->pts[c_pt + c - 1].coord[1] + y;
         list->pts[c_pt + c].coord[2] = list->pts[c_pt + c - 1].coord[2] + z;
         }
      }
   else{
      int      i;

      fprintf(stdout, "Multiplying");
      for(c = 0; c < rep; c++){

         for(i = 0; i < 3; i++){
            list->pts[c_pt + c].coord[i] = list->pts[c_pt + c - 1].coord[i] +
               (aff_M[i][0] * x) + (aff_M[i][1] * y) + (aff_M[i][2] * z) + aff_M[i][3];
            }
         }
      }

   return 1;
   }


int call_list(char *name)
{
   int      c = 0;
   Coord_list list = NULL;

   /* get the list by name */
   while (c < num_lists && list == NULL){
      if(strcmp(coord_lists[c]->name, name) == 0){
         list = coord_lists[c];
         }
      c++;
      }

   if(list == NULL){
      fprintf(stderr, "Couldn't find list called '%s'\n", name);
      return 0;
      }

   /* copy the co-ords into the main (0) list */
   for(c = 0; c < list->n_pts; c++){
      add_coord_to_list(0, list->pts[c].coord[0],
                        list->pts[c].coord[1], list->pts[c].coord[2]);
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

   fprintf(stdout, "Curr Matrix - %sidentity\n", (aff_M_identity) ? "" : "non-");
   for(i = 0; i < 3; i++){
      fprintf(stdout, "  %6.2f %6.2f %6.2f %6.2f\n",
              aff_M[i][0], aff_M[i][1], aff_M[i][2], aff_M[i][3]);
      }
   }
