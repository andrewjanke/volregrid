/* arb_path_io.h */

#ifndef ARB_PATH_IO_H
#define ARB_PATH_IO_H

#include <netcdf.h>

/* Structures for regrid information */
typedef struct {
   double   coord[3];
} Coord;

struct coord_list_struct {
   int      id;
   char     name[256];
   size_t   n_pts;
   size_t   alloc_size;
   Coord   *pts;
};
typedef struct coord_list_struct *Coord_list;

/* initialise the parser with a filename */
int      init_arb_path(char *coord_fname, char *data_fname);

/* finalize the parser */
void     end_arb_path(void);

/* get some co-ordinates, returns an empty list when done      */
/* a max buffer size can be input (but might not be obeyed..)  */
Coord_list get_some_arb_path_coords(size_t max_buf_size);

/* get some arb_path data, places data in data_buf             */
int      get_some_arb_path_data(double *data_buf, nc_type dtype, int is_signed,
                                size_t n_pts, size_t vector_size);

/* functions called by lex */
/* update the current affine matrix */
void     set_curr_matrix(double mat[12]);

/* return a new Coord_list's id */
int      new_coord_list(size_t size, char *name);

int      add_coord_to_list(int list_id, double x, double y, double z);
int      add_rcoord_to_list(int list_id, double x, double y, double z, double rep);

/* call a list by name */
int      call_list(char *name);



/* debugging functions */
void     print_coord_list(Coord_list list);
void     print_curr_matrix(void);

#endif
