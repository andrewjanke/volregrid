/* volregrid.c                                                               */
/*                                                                           */
/* Performs regridding on a MINC volume or raw input data                    */
/*                                                                           */
/* Andrew Janke - rotor@cmr.uq.edu.au                                        */
/* Mark Griffin - mark@cmr.uq.edu.au                                         */
/* Center for Magnetic Resonance                                             */
/* University of Queensland                                                  */
/*                                                                           */
/* Copyright Andrew Janke, The University of Queensland.                     */
/* Permission to use, copy, modify, and distribute this software and its     */
/* documentation for any purpose and without fee is hereby granted,          */
/* provided that the above copyright notice appear in all copies.  The       */
/* author and the University of Queensland make no representations about the */
/* suitability of this software for any purpose.  It is provided "as is"     */
/* without express or implied warranty.                                      */
/*                                                                           */
/* Wed Dec 11 14:57:44 EST 2002 - inital version for sodium k-space data     */
/* Wed Dec 18 09:20:02 EST 2002 - added lex parser for arbitrary path data   */


#include <float.h>
#include <sys/stat.h>
#include <ctype.h>
#include <volume_io.h>
#include <ParseArgv.h>
#include <time_stamp.h>
#include "arb_path_io.h"

#define NO_VALUE DBL_MAX               /* Constant to flag fact that value not set */
#define DEF_BOOL -1
#define X_IDX 0
#define Y_IDX 1
#define Z_IDX 2
#define V_IDX 3

/* function prototypes */
int      read_config_file(char *filename, char *args[]);

int      verbose = FALSE;
int      debug = FALSE;
int      clobber = FALSE;
int      regrid_dim = 3;
char    *arb_path_cfg_fn = NULL;
int      arb_path_buff_size = 100;
nc_type  in_dtype = NC_FLOAT;
int      in_is_signed = FALSE;

/* output file parameters */
char    *out_config_fn = NULL;
nc_type  out_dtype = NC_UNSPECIFIED;
int      out_is_signed = DEF_BOOL;
double   valid_range[2] = { -DBL_MAX, DBL_MAX };
double   real_range[2] = { 0, 1 };
double   out_start[MAX_VAR_DIMS] = { -50, -50, -50 };
double   out_step[MAX_VAR_DIMS] = { 1, 1, 1 };
int      out_length[MAX_VAR_DIMS] = { 100, 100, 100, 1 };
double   out_dircos[3][3] = {
   1, 0, 0,
   0, 1, 0,
   0, 0, 1
};

static ArgvInfo argTable[] = {
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "Print out extra information."},
   {"-debug", ARGV_CONSTANT, (char *)TRUE, (char *)&debug,
    "Spew copious amounts of debugging info."},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "Clobber existing files."},

   {NULL, ARGV_HELP, NULL, NULL, "\nRaw Infile Options"},
   {"-byte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&in_dtype,
    "Input data is byte data."},
   {"-short", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&in_dtype,
    "Input data is short integer data."},
   {"-int", ARGV_CONSTANT, (char *)NC_INT, (char *)&in_dtype,
    "Input data is 32-bit integer"},
   {"-float", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&in_dtype,
    "Input data is single-precision data. (Default)"},
   {"-double", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&in_dtype,
    "Input data is double-precision data."},
   {"-signed", ARGV_CONSTANT, (char *)TRUE, (char *)&in_is_signed,
    "Input data is signed integer data."},
   {"-ounsigned", ARGV_CONSTANT, (char *)FALSE, (char *)&in_is_signed,
    "Input data is unsigned integer data."},
   {"-range", ARGV_FLOAT, (char *)2, (char *)valid_range,
    "Valid range of input values (default = full range)."},
   {"-real_range", ARGV_FLOAT, (char *)2, (char *)real_range,
    "Real range of input values (ignored for floating-point types)."},


   {NULL, ARGV_HELP, NULL, NULL, "\nOutfile Options"},
   {"-outconfig", ARGV_STRING, (char *)1, (char *)&out_config_fn,
    "Get the output geometry from the input filename (overrides args below)"},
   {"-obyte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&out_dtype,
    "Write out byte data."},
   {"-oshort", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&out_dtype,
    "Write out short integer data."},
   {"-int", ARGV_CONSTANT, (char *)NC_INT, (char *)&in_dtype,
    "Input data is 32-bit integer"},
   {"-ofloat", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&out_dtype,
    "Write out single-precision data. (Default)"},
   {"-odouble", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&out_dtype,
    "Write out double-precision data."},
   {"-osigned", ARGV_CONSTANT, (char *)TRUE, (char *)&out_is_signed,
    "Write signed integer data."},
   {"-ounsigned", ARGV_CONSTANT, (char *)FALSE, (char *)&out_is_signed,
    "Write unsigned integer data."},
   {"-xstart", ARGV_FLOAT, (char *)1, (char *)&out_start[X_IDX],
    "Starting coordinate for x dimension."},
   {"-ystart", ARGV_FLOAT, (char *)1, (char *)&out_start[Y_IDX],
    "Starting coordinate for y dimension."},
   {"-zstart", ARGV_FLOAT, (char *)1, (char *)&out_start[Z_IDX],
    "Starting coordinate for z dimension."},
   {"-xstep", ARGV_FLOAT, (char *)1, (char *)&out_step[X_IDX],
    "Step size for x dimension."},
   {"-ystep", ARGV_FLOAT, (char *)1, (char *)&out_step[Y_IDX],
    "Step size for y dimension."},
   {"-zstep", ARGV_FLOAT, (char *)1, (char *)&out_step[Z_IDX],
    "Step size for z dimension."},
   {"-xlength", ARGV_INT, (char *)1, (char *)&out_length[X_IDX],
    "Number of samples in x dimension."},
   {"-ylength", ARGV_INT, (char *)1, (char *)&out_length[Y_IDX],
    "Number of samples in y dimension."},
   {"-zlength", ARGV_INT, (char *)1, (char *)&out_length[Z_IDX],
    "Number of samples in z dimension."},
   {"-vector_length", ARGV_INT, (char *)1, (char *)&out_length[V_IDX],
    "Number of samples in the vector dimension."},
   {"-xdircos", ARGV_FLOAT, (char *)3, (char *)out_dircos[X_IDX],
    "Direction cosines for x dimension."},
   {"-ydircos", ARGV_FLOAT, (char *)3, (char *)out_dircos[Y_IDX],
    "Direction cosines for y dimension."},
   {"-zdircos", ARGV_FLOAT, (char *)3, (char *)out_dircos[Z_IDX],
    "Direction cosines for z dimension."},


   {NULL, ARGV_HELP, NULL, NULL, "\nRegridding options"},
   {"-2D", ARGV_CONSTANT, (char *)2, (char *)&regrid_dim,
    "Regrid slice by slice (Default 3D)."},

   {"-arb_path", ARGV_STRING, (char *)1, (char *)&arb_path_cfg_fn,
    "Regrid data using an arbitrary path from the input filename"},
   {"-arb_path_coord_buffer", ARGV_INT, (char *)1, (char *)&arb_path_buff_size,
    "Size of arbitrary path co-ordinate buffer"},


   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

char    *std_dimorder[] = { MIzspace, MIyspace, MIxspace };
char    *std_dimorder_v[] = { MIzspace, MIyspace, MIxspace, MIvector_dimension };

main(int argc, char *argv[])
{
   char    *in_fn, *out_fn;
   char    *history;
   Status   status;
   Volume   out_vol;

   /* get the history string */
   history = time_stamp(argc, argv);

   /* get args */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 3)){
      fprintf(stderr, "\nUsage: %s [<options>] <infile.raw> <outfile.mnc>\n", argv[0]);
      fprintf(stderr, "       %s [-help]\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }
   in_fn = argv[1];
   out_fn = argv[2];

   /* check for infile and outfile */
   if(!file_exists(in_fn)){
      fprintf(stderr, "%s: Couldn't find input file %s.\n\n", argv[0], in_fn);
      exit(EXIT_FAILURE);
      }
   if(!clobber && file_exists(out_fn)){
      fprintf(stderr, "%s: File %s exists, use -clobber to overwrite.\n\n", argv[0],
              out_fn);
      exit(EXIT_FAILURE);
      }

   /* set up parameters for reconstruction */
   if(out_dtype == NC_UNSPECIFIED){
      out_dtype = in_dtype;
      }
   if(out_is_signed == DEF_BOOL){
      out_is_signed = in_is_signed;
      }

   /* check vector dimension size */
   if(out_length[V_IDX] < 1){
      fprintf(stderr, "%s: -vector_length must be 1 or greater.\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }


   /* read in the output file config from a file is specified */
   if(out_config_fn != NULL){
      int      ext_args_c;
      char    *ext_args[32];           /* max possible is 32 arguments */

      ext_args_c = read_config_file(out_config_fn, ext_args);
      if(ParseArgv(&ext_args_c, ext_args, argTable,
                    ARGV_DONT_SKIP_FIRST_ARG || ARGV_NO_LEFTOVERS || ARGV_NO_DEFAULTS)){
         fprintf(stderr, "\nError in parameters in %s\n", out_config_fn);
         exit(EXIT_FAILURE);
         }
      }

   /* create the output volume */
   out_vol = create_volume((out_length[V_IDX] > 1) ? 4 : 3,
                           (out_length[V_IDX] > 1) ? std_dimorder_v : std_dimorder,
                           out_dtype, out_is_signed, 0.0, 0.0);
   set_volume_sizes(out_vol, out_length);
   set_volume_starts(out_vol, out_start);
   set_volume_separations(out_vol, out_step);
   set_volume_real_range(out_vol, valid_range[0], valid_range[1]);
   alloc_volume_data(out_vol);

   /* print some pretty output */
   if(verbose){

      fprintf(stdout, " | Input data:      %s\n", in_fn);
      fprintf(stdout, " | Arb path:        %s\n", arb_path_cfg_fn);
      fprintf(stdout, " | Valid Range:    [%8.3f:%8.3f]\n", valid_range[0],
              valid_range[1]);
      fprintf(stdout, " | Real Range:     [%8.3f:%8.3f]\n", real_range[0], real_range[1]);
      fprintf(stdout, " | Output file:     %s\n", out_fn);
      }


   /* arbitrary path */
   if(arb_path_cfg_fn != NULL){
      Coord_list coord_buf;
      double  *data_buf = NULL;
      size_t   data_buf_alloc_size = 0;

      /* check for the config file */
      if(!file_exists(arb_path_cfg_fn)){
         fprintf(stderr, "%s: Couldn't find config file %s.\n\n", argv[0],
                 arb_path_cfg_fn);
         exit(EXIT_FAILURE);
         }

      /* initialise the parser with the config file */
      if(!init_arb_path(arb_path_cfg_fn, in_fn)){
         fprintf(stderr, "%s: Failed to init arb_path, this isn't good\n", argv[0]);
         exit(EXIT_FAILURE);
         }

      /* get some co-ordinates */
      coord_buf = get_some_arb_path_coords(arb_path_buff_size);
      while(coord_buf->n_pts != 0){
         
      //   print_coord_list(coord_buf);
         fprintf(stdout, "Got %d co-ords\n", coord_buf->n_pts);

         /* grow data_buf if we have to */
         if(coord_buf->n_pts > data_buf_alloc_size){
            data_buf_alloc_size = coord_buf->n_pts;
            data_buf = (double *)realloc(data_buf, data_buf_alloc_size * out_length[V_IDX]);
            fprintf(stdout, "Allocated %d x %d = %d bytes to data_buf\n", data_buf_alloc_size,
                 out_length[V_IDX], data_buf_alloc_size * out_length[V_IDX]);
            }
         
         /* get the data */
         if(!get_some_arb_path_data(data_buf, coord_buf->n_pts, out_length[V_IDX])){
            fprintf(stderr, "failed getting data\n");
            }
      

         /* regrid (do the nasty) */
//      if(regrid_arbitrary_path(out_vol, in_fn, path_data) != OK){
//         print_error("%s: Problems regridding data\n", argv[0]);
//         exit(EXIT_FAILURE);
//         }
         
         /* get the next lot of co-ordinates */
         coord_buf = get_some_arb_path_coords(arb_path_buff_size);
         }

      end_arb_path();
      }

   /* output the result */
   if(verbose){
      fprintf(stdout, "Outputting %s...\n", out_fn);
      }
   status = output_volume(out_fn, out_dtype, out_is_signed, 0, 0, out_vol, history, NULL);
   if(status != OK){
      print_error("Problems outputing: %s", out_fn);
      }

   delete_volume(out_vol);
   return (status);
   }


/* regrid a volume with an arbitrary path */



/* convert input file to equiv C/L, return number of args read */
int read_config_file(char *filename, char *args[])
{
   FILE    *fp;
   int      ch;
   char     tmp[256];
   int      nargs = 0;
   int      ichar = 0;
   int      in_comment = FALSE;

   /* Open the file */
   if((fp = fopen(filename, "r")) == NULL){
      fprintf(stderr, "Unable to open options file %s\n\n", filename);
      exit(EXIT_FAILURE);
      }

   /* Read in the arguments, skipping comments and poking them in an array */
   while ((ch = getc(fp)) != EOF){
      switch (ch){
      case '#':
         in_comment = TRUE;
         break;

      case '\n':
         in_comment = FALSE;
      default:
         if(!in_comment){
            if(isspace(ch)){
               /* we have a complete arg, copy it over */
               if(ichar != 0){
                  tmp[ichar++] = '\0';
                  args[nargs] = (char *)malloc(ichar * sizeof(*args[nargs]));
                  strncpy(args[nargs], tmp, ichar);

                  ichar = 0;
                  nargs++;
                  }
               }
            else{
               tmp[ichar++] = ch;
               }
            }
         }
      }
   fclose(fp);
   return nargs;
   }
