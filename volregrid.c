/* volregrid.c                                                               */
/*                                                                           */
/* Performs regridding on a series of input MINC volumes or raw input data   */
/*                                                                           */
/* Andrew Janke - rotor@cmr.uq.edu.au                                        */
/* Mark Griffin - mark@cmr.uq.edu.au                                         */
/* Center for Magnetic Resonance                                             */
/* University of Queensland                                                  */
/*                                                                           */
/* Copyright Andrew Janke & Mark Griffin The University of Queensland.       */
/* Permission to use, copy, modify, and distribute this software and its     */
/* documentation for any purpose and without fee is hereby granted,          */
/* provided that the above copyright notice appear in all copies.  The       */
/* author and the University of Queensland make no representations about the */
/* suitability of this software for any purpose.  It is provided "as is"     */
/* without express or implied warranty.                                      */
/*                                                                           */
/* Wed Dec 11 14:57:44 EST 2002 - inital version for sodium k-space data     */
/* Wed Dec 18 09:20:02 EST 2002 - added lex parser for arbitrary path data   */
/* Fri Jan 10 15:21:30 EST 2003 - First working version                      */

#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include <ctype.h>
#include <volume_io.h>
#include <ParseArgv.h>
#include <time_stamp.h>
#include <voxel_loop.h>
#include <gsl/gsl_sf_bessel.h>
#include "arb_path_io.h"
#include "minc_support.h"

#define SQR2(x) ((x)*(x))
#define SQR3(x) ((x)*(x)*(x))

#define DEF_BOOL -1
#define X_IDX 2
#define Y_IDX 1
#define Z_IDX 0
#define V_IDX 3

/* permutation array for IDX's */
static int perm[3] = { X_IDX, Y_IDX, Z_IDX };

static char *std_dimorder[] = { MIzspace, MIyspace, MIxspace };
static char *std_dimorder_v[] = { MIzspace, MIyspace, MIxspace, MIvector_dimension };

/* enum for regridding types */
typedef enum {
   UNSPECIFIED_FUNC = 0,
   KAISERBESSEL_FUNC,
   GAUSSIAN_FUNC,
   NEAREST_FUNC
   } Regrid_op;

/* function prototypes */
int      read_config_file(char *filename, char *args[]);
int      get_model_file_info(char *dst, char *key, char *nextArg);
void     scale_volume(Volume * vol, double o_min, double o_max, double min, double max);
void     regrid_point(Volume * totals, Volume * counts,
                      double x, double y, double z, int v_size, double *data_buf);
void     regrid_arb_path(char *coord_fn, char *data_fn, int buff_size,
                         Volume * totals, Volume * counts, int v_size);
void     regrid_loop(void *caller_data, long num_voxels,
                     int input_num_buffers, int input_vector_length,
                     double *input_data[],
                     int output_num_buffers, int output_vector_length,
                     double *output_data[], Loop_Info * loop_info);
void     regrid_minc(char *in_fn, Volume * totals, Volume * counts, int v_size);

/* argument variables and table */
static int verbose = FALSE;
static int clobber = FALSE;
static nc_type in_dtype = NC_FLOAT;
static int in_is_signed = FALSE;
static int max_buffer_size_in_kb = 4 * 1024;
static int vect_size = 1;

/* arb path variables */
char    *ap_coord_fn = NULL;

/* regridding options */
double   regrid_radius = 2.0;
Regrid_op regrid_type = GAUSSIAN_FUNC;
double   regrid_sigma = 1.0;

/* output file parameters */
char    *out_config_fn = NULL;
nc_type  out_dtype = NC_UNSPECIFIED;
int      out_is_signed = DEF_BOOL;
double   out_range[2] = { -DBL_MAX, DBL_MAX };

Volume_Definition out_inf = {
   3,

   {MIxspace, MIyspace, MIzspace},     /* dimnames */

   {Z_IDX, Y_IDX, X_IDX},              /* space to dim */
   {X_IDX, Y_IDX, Z_IDX},              /* dim to space */

   {100, 100, 100},                    /* nelem */
   {-50.0, -50.0, -50.0},              /* start */
   {1.0, 1.0, 1.0},                    /* step */
   {{1.0, 0.0, 0.0},                   /* direction cosines */
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}},

   {{1.0, 0.0, 0.0, -50.0},            /* voxel to world */
    {0.0, 1.0, 0.0, -50.0},
    {0.0, 0.0, 1.0, -50.0}}
   };

static ArgvInfo argTable[] = {
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "Print out extra information."},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "Overwrite existing files."},
   {"-max_buffer_size_in_kb", ARGV_INT, (char *)1, (char *)&max_buffer_size_in_kb,
    "maximum size of internal buffers."},

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
   {"-unsigned", ARGV_CONSTANT, (char *)FALSE, (char *)&in_is_signed,
    "Input data is unsigned integer data. (Default)"},
   {"-vector", ARGV_INT, (char *)1, (char *)&vect_size,
    "Size of vector dimension of Input data."},

   {NULL, ARGV_HELP, NULL, NULL, "\nOutfile Options"},
   {"-outconfig", ARGV_STRING, (char *)1, (char *)&out_config_fn,
    "Get the output geometry from the input filename (overrides args below)"},
   {"-like", ARGV_FUNC, (char *)get_model_file_info, (char *)&out_inf,
    "Specifies a model file for the output geometry."},
   {"-obyte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&out_dtype,
    "Write out byte data."},
   {"-oshort", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&out_dtype,
    "Write out short integer data."},
   {"-oint", ARGV_CONSTANT, (char *)NC_INT, (char *)&out_dtype,
    "Write out 32-bit integer data"},
   {"-ofloat", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&out_dtype,
    "Write out single-precision data. (Default)"},
   {"-odouble", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&out_dtype,
    "Write out double-precision data."},
   {"-osigned", ARGV_CONSTANT, (char *)TRUE, (char *)&out_is_signed,
    "Write signed integer data."},
   {"-ounsigned", ARGV_CONSTANT, (char *)FALSE, (char *)&out_is_signed,
    "Write unsigned integer data."},
   {"-range", ARGV_FLOAT, (char *)2, (char *)out_range,
    "Range to scale output values between (default = range of input data)."},
   {"-xnelements", ARGV_INT, (char *)1, (char *)&out_inf.nelem[0],
    "Number of samples in x dimension."},
   {"-ynelements", ARGV_INT, (char *)1, (char *)&out_inf.nelem[1],
    "Number of samples in y dimension."},
   {"-znelements", ARGV_INT, (char *)1, (char *)&out_inf.nelem[2],
    "Number of samples in z dimension."},
   {"-xstart", ARGV_FLOAT, (char *)1, (char *)&out_inf.start[0],
    "Starting coordinate for x dimension."},
   {"-ystart", ARGV_FLOAT, (char *)1, (char *)&out_inf.start[1],
    "Starting coordinate for y dimension."},
   {"-zstart", ARGV_FLOAT, (char *)1, (char *)&out_inf.start[2],
    "Starting coordinate for z dimension."},
   {"-xstep", ARGV_FLOAT, (char *)1, (char *)&out_inf.step[0],
    "Step size for x dimension."},
   {"-ystep", ARGV_FLOAT, (char *)1, (char *)&out_inf.step[1],
    "Step size for y dimension."},
   {"-zstep", ARGV_FLOAT, (char *)1, (char *)&out_inf.step[2],
    "Step size for z dimension."},
   {"-xdircos", ARGV_FLOAT, (char *)3, (char *)out_inf.dircos[0],
    "Direction cosines along the x dimension"},
   {"-ydircos", ARGV_FLOAT, (char *)3, (char *)out_inf.dircos[1],
    "Direction cosines along the y dimension"},
   {"-zdircos", ARGV_FLOAT, (char *)3, (char *)out_inf.dircos[2],
    "Direction cosines along the z dimension"},

   {NULL, ARGV_HELP, NULL, NULL, "\nRegridding options"},
   {"-regrid_radius", ARGV_FLOAT, (char *)1, (char *)&regrid_radius,
    "Defines the Window radius for regridding (in mm)."},
   {"-kaiser_bessel", ARGV_CONSTANT, (char *)KAISERBESSEL_FUNC, (char *)&regrid_type,
    "Use a Kaiser-Bessel convolution kernel for the reconstruction."},
   {"-gaussian", ARGV_CONSTANT, (char *)GAUSSIAN_FUNC, (char *)&regrid_type,
    "Use a Gaussian convolution kernel for the reconstruction."},
   {"-nearest", ARGV_CONSTANT, (char *)NEAREST_FUNC, (char *)&regrid_type,
    "Use nearest neighbour reconstruction."},
   {"-sigma", ARGV_FLOAT, (char *)1, (char *)&regrid_sigma,
    "Value of sigma for -gaussian and -kaiser_bessel function"},

   {NULL, ARGV_HELP, NULL, NULL, "\nArbitrary path Regridding options"},
   {"-arb_path", ARGV_STRING, (char *)1, (char *)&ap_coord_fn,
    "<file> Regrid data using an arbitrary path from the input file"},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
   };

int main(int argc, char *argv[])
{
   char   **infiles;
   int      n_infiles;
   char    *out_fn;
   char    *history;
   progress_struct progress;
   Status   status;
   Volume   totals, counts;
   int      i, j, k, v;
   double   min, max;
   long     num_missed;
   double   weight, value;
   
   int sizes[MAX_VAR_DIMS];
   double starts[MAX_VAR_DIMS];
   double steps[MAX_VAR_DIMS];
   Real tmp_dircos[WORLD_NDIMS];

   /* get the history string */
   history = time_stamp(argc, argv);

   /* get args */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 3)){
      fprintf(stderr, "\nUsage: %s [options] <infile.raw> <out.mnc>\n", argv[0]);
      fprintf(stderr, "       %s [options] <in1.mnc> [<in2.mnc> [...]] <out.mnc>\n",
              argv[0]);
      fprintf(stderr, "       %s [-help]\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   /* get file names */
   n_infiles = argc - 2;
   infiles = (char **)malloc(sizeof(char *) * n_infiles);
   for(i = 0; i < n_infiles; i++){
      infiles[i] = argv[i + 1];
      }
   out_fn = argv[argc - 1];

   /* check for infiles and outfile */
   for(i = 0; i < n_infiles; i++){
      if(!file_exists(infiles[i])){
         fprintf(stderr, "%s: Couldn't find input file %s.\n\n", argv[0], infiles[i]);
         exit(EXIT_FAILURE);
         }
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
   if(vect_size < 1){
      fprintf(stderr, "%s: -vector (%d) must be 1 or greater.\n\n", argv[0], vect_size);
      exit(EXIT_FAILURE);
      }

   /* check sigma */
   if(regrid_sigma <= 0){
      fprintf(stderr, "%s: -sigma must be greater than 0.\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   /* read in the output file config from a file is specified */
   if(out_config_fn != NULL){
      int      ext_args_c;
      char    *ext_args[32];           /* max possible is 32 arguments */

      ext_args_c = read_config_file(out_config_fn, ext_args);
      if(ParseArgv(&ext_args_c, ext_args, argTable,
                   ARGV_DONT_SKIP_FIRST_ARG | ARGV_NO_LEFTOVERS | ARGV_NO_DEFAULTS)){
         fprintf(stderr, "\nError in parameters in %s\n", out_config_fn);
         exit(EXIT_FAILURE);
         }
      }

   if(verbose){
      fprintf_vol_def(stdout, &out_inf);
      }
   
   /* transpose the geometry arrays */
   for(i=0; i<WORLD_NDIMS; i++){
      sizes[i] = out_inf.nelem[perm[i]];
      starts[i] = out_inf.start[perm[i]];
      steps[i] = out_inf.step[perm[i]];
      }
   sizes[WORLD_NDIMS] = vect_size; 
   
   /* create the totals volume */
   totals = create_volume((vect_size > 1) ? 4 : 3,
                          (vect_size > 1) ? std_dimorder_v : std_dimorder,
                          out_dtype, out_is_signed, 0.0, 0.0);
   set_volume_sizes(totals, sizes);
   set_volume_starts(totals, starts);
   set_volume_separations(totals, steps);
   for(i=0; i<WORLD_NDIMS; i++){
      for(j=0; j<WORLD_NDIMS; j++){
         tmp_dircos[j] = (Real)out_inf.dircos[perm[i]][j];
         fprintf(stdout, "[%d] Dircos[%d] %g\n", i, j, tmp_dircos[j]); 
         }
      set_volume_direction_cosine(totals, i, tmp_dircos);
      }
   alloc_volume_data(totals);
   
   /* create the "counts" volume */
   counts = create_volume(3, std_dimorder, out_dtype, out_is_signed, 0.0, 0.0);
   set_volume_sizes(counts, sizes);
   set_volume_starts(counts, starts);
   set_volume_separations(counts, steps);
   alloc_volume_data(counts);

   /* initialize both of them */
   for(k = sizes[Z_IDX]; k--;){
      for(j = sizes[Y_IDX]; j--;){
         for(i = sizes[X_IDX]; i--;){
            set_volume_real_value(counts, k, j, i, 0, 0, 0.0);
            for(v = vect_size; v--;){
               set_volume_real_value(totals, k, j, i, v, 0, 0.0);
               }
            }
         }
      }

   /* if regridding via an arbitrary path */
   if(ap_coord_fn != NULL){

      if(n_infiles > 1){
         fprintf(stderr, "%s: arb_path only works for one input file (so far).\n\n",
                 argv[0]);
         exit(EXIT_FAILURE);
         }

      /* print some pretty output */
      if(verbose){
         fprintf(stdout, " | Input data:      %s\n", infiles[0]);
         fprintf(stdout, " | Arb path:        %s\n", ap_coord_fn);
         fprintf(stdout, " | Output range:    [%g:%g]\n", out_range[0], out_range[1]);
         fprintf(stdout, " | Output file:     %s\n", out_fn);
         }

      regrid_arb_path(ap_coord_fn, infiles[0], max_buffer_size_in_kb, &totals, &counts,
                      vect_size);
      }

   /* else if regridding via a series of input minc file(s) */
   else {
      int      c;

      for(c = 0; c < n_infiles; c++){

         /* print some pretty output */
         if(verbose){
            fprintf(stdout, " | Input file:      %s\n", infiles[c]);
            }

         regrid_minc(infiles[c], &totals, &counts, vect_size);
         }
      }

   /* divide totals/counts volume for result */
   min = DBL_MAX;
   max = -DBL_MAX;
   num_missed = 0;
   initialize_progress_report(&progress, FALSE, out_inf.nelem[Z_IDX], "Dividing through");
   for(k = sizes[Z_IDX]; k--;){
      for(j = sizes[Y_IDX]; j--;){
         for(i = sizes[X_IDX]; i--;){
            weight = get_volume_real_value(counts, k, j, i, 0, 0);
            if(weight != 0){
               for(v = vect_size; v--;){
                  value = get_volume_real_value(totals, k, j, i, v, 0) / weight;
                  if(value < min){
                     min = value;
                     }
                  if(value > max){
                     max = value;
                     }

                  set_volume_real_value(totals, k, j, i, v, 0, value);
                  }
               }
            else {
               num_missed++;
               }
            }
         }
      update_progress_report(&progress, k + 1);
      }
   terminate_progress_report(&progress);

   /* set the volumes (initial) range */
   if(verbose){
      fprintf(stdout, " + data range: [%g:%g]\n", min, max);
      }
   set_volume_real_range(totals, min, max);

   if(num_missed > 0 && verbose){
      int      nvox;

      nvox = out_inf.nelem[X_IDX] * out_inf.nelem[Y_IDX] * out_inf.nelem[Z_IDX];
      fprintf(stdout,
              "\n-regrid_radius possibly too small, no data in %ld/%d[%2.2f%%] voxels\n\n",
              num_missed, nvox, ((float)num_missed / nvox * 100));
      }

   /* rescale data if required */
   if(out_range[0] != -DBL_MAX && out_range[1] != DBL_MAX){
      double   o_min, o_max;

      /* get the existing range */
      get_volume_real_range(totals, &o_min, &o_max);

      /* rescale it */
      scale_volume(&totals, o_min, o_max, out_range[0], out_range[1]);
      }

   /* output the result */
   if(verbose){
      fprintf(stdout, "Outputting %s...\n", out_fn);
      }
   status =
      output_volume(out_fn, out_dtype, out_is_signed, 0.0, 0.0, totals, history, NULL);
   if(status != OK){
      print_error("Problems outputing: %s", out_fn);
      }

   delete_volume(totals);
   delete_volume(counts);

   return (EXIT_SUCCESS);
   }

/* re-scale a volume between a new range                                     */
/* rescaling is done based upon the formula:                                 */
/*                                                                           */
/*                         (max - min)                                       */
/*    x' = (x - o_min) * --------------- + min                               */
/*                       (o_max - o_min)                                     */
/*                                                                           */
/* or more succinctly (and speedily)                                         */
/*                                                                           */
/*    x' = ax + b                                                            */
/*                                                                           */
/* where                                                                     */
/*          (max - min)                                                      */
/*    a = ---------------                                                    */
/*        (o_max - o_min)                                                    */
/*                                                                           */
/*    b = min - (o_min * a)                                                  */
/*                                                                           */
void scale_volume(Volume * vol, double o_min, double o_max, double min, double max)
{
   double   value, a, b;
   int      sizes[MAX_VAR_DIMS];
   int      i, j, k, v;
   progress_struct progress;

   get_volume_sizes(*vol, sizes);

   /* rescale the volume */
   a = (max - min) / (o_max - o_min);
   b = min - (o_min * a);
   initialize_progress_report(&progress, FALSE, sizes[Z_IDX], "Rescaling Volume");
   for(k = sizes[Z_IDX]; k--;){
      for(j = sizes[Y_IDX]; j--;){
         for(i = sizes[X_IDX]; i--;){
            for(v = sizes[V_IDX]; v--;){
               value = (get_volume_real_value(*vol, k, j, i, v, 0) * a) + b;

               set_volume_real_value(*vol, k, j, i, v, 0, value);
               }
            }
         }
      update_progress_report(&progress, k + 1);
      }
   terminate_progress_report(&progress);

   if(verbose){
      fprintf(stdout, " + rescaled data range: [%g:%g]\n", min, max);
      }

   set_volume_real_range(*vol, min, max);
   }

/* struct for regriding using a minc volume */
typedef struct {

   int      file_ndims;
   int      space_to_dim[WORLD_NDIMS];
   int      dim_to_space[MAX_VAR_DIMS];
   double   voxel_to_world[WORLD_NDIMS][WORLD_NDIMS + 1];

   double  *data_buf;

   Volume  *totals;
   Volume  *counts;

   } Loop_Data;

/* voxel loop function for regrid_minc */
void regrid_loop(void *caller_data, long num_voxels,
                 int input_num_buffers, int input_vector_length,
                 double *input_data[],
                 int output_num_buffers, int output_vector_length,
                 double *output_data[], Loop_Info * loop_info)
{
   long     ivox;
   long     idx[MAX_VAR_DIMS];
   int      v, idim;
   double   voxel_coord[WORLD_NDIMS];
   double   world_coord[WORLD_NDIMS];

   /* get pointer to loop data */
   Loop_Data *ld = (Loop_Data *) caller_data;

   /* shut the compiler up - yes I _know_ I don't use these */
   (void)input_num_buffers;
   (void)output_num_buffers;
   (void)output_vector_length;
   (void)output_data;

   /* for each (vector) voxel */
   for(ivox = 0; ivox < num_voxels * (long)input_vector_length;
       ivox += (long)input_vector_length){

      /* figure out where we are in space */
      get_info_voxel_index(loop_info, ivox, ld->file_ndims, idx);

      /* convert voxel index to world co-ordinate */
      for(idim = 0; idim < WORLD_NDIMS; idim++){
         voxel_coord[idim] = idx[ld->space_to_dim[idim]];
         }
      transform_coord(&world_coord[0], ld->voxel_to_world, &voxel_coord[0]);

      /* get the data */
      for(v = 0; v < input_vector_length; v++){
         ld->data_buf[v] = input_data[0][(ivox * (long)input_vector_length) + (long)v];
         }

      /* then regrid the point */
      regrid_point(ld->totals, ld->counts,
                   world_coord[0], world_coord[1], world_coord[2],
                   input_vector_length, ld->data_buf);
      }

   return;
   }

/* regrid using a minc volume */
void regrid_minc(char *in_fn, Volume * totals, Volume * counts, int v_size)
{
   Loop_Data ld;
   Loop_Options *loop_opt;
   int      mincid;

   /* Open the file to get some information */
   mincid = miopen(in_fn, NC_NOWRITE);
   ld.file_ndims = get_minc_ndims(mincid);
   get_minc_spatial_dims(mincid, ld.space_to_dim, ld.dim_to_space);
   get_minc_voxel_to_world(mincid, ld.voxel_to_world);

   /* alloc space for data_buf */
   ld.data_buf = (double *)malloc(sizeof(double) * v_size);

   ld.totals = totals;
   ld.counts = counts;

   /* set up and do voxel_loop */
   loop_opt = create_loop_options();
   set_loop_first_input_mincid(loop_opt, mincid);
   set_loop_verbose(loop_opt, verbose);
   set_loop_buffer_size(loop_opt, (long)1024 * max_buffer_size_in_kb);
   voxel_loop(1, &in_fn, 0, NULL, NULL, loop_opt, regrid_loop, (void *)&ld);
   free_loop_options(loop_opt);

   /* tidy up */
   free(ld.data_buf);
   }

/* regrid a point in a file using the input co-ordinate and data */
void regrid_point(Volume * totals, Volume * counts,
                  double x, double y, double z, int v_size, double *data_buf)
{

   int      sizes[MAX_VAR_DIMS];
   Real     steps[MAX_VAR_DIMS];
   Real     starts[MAX_VAR_DIMS];
   int      start_idx[3];
   int      stop_idx[3];
   double   value, weight;
   double   euc_dist;
   double   euc[3];
   double   c_pos[3];
   int      i, j, k, v;
   double   coord[3];

   coord[0] = x;
   coord[1] = y;
   coord[2] = z;

   get_volume_sizes(*totals, sizes);
   get_volume_separations(*totals, steps);
   get_volume_starts(*totals, starts);

   /* figure out the neighbouring voxels start and stop (in voxel co-ordinates) */
   for(i = 0; i < 3; i++){
      start_idx[i] =
         (int)rint((coord[i] - starts[perm[i]] - regrid_radius) / steps[perm[i]]);
      stop_idx[i] = start_idx[i] + rint((regrid_radius * 2) / steps[perm[i]]);

      /* check that we aren't off the edge */
      if(start_idx[i] < 0){
         start_idx[i] = 0;
         }
      if(stop_idx[i] > sizes[perm[i]]){
         stop_idx[i] = sizes[perm[i]];
         }
      }

   /* loop over the neighbours, getting euclidian distance */
   c_pos[0] = starts[perm[0]] + (start_idx[0] * steps[perm[0]]);
   for(i = start_idx[0]; i < stop_idx[0]; i++){
      euc[0] = fabs(c_pos[0] - x);

      c_pos[1] = starts[perm[1]] + (start_idx[1] * steps[perm[1]]);
      for(j = start_idx[1]; j < stop_idx[1]; j++){
         euc[1] = fabs(c_pos[1] - y);

         c_pos[2] = starts[perm[2]] + (start_idx[2] * steps[perm[2]]);
         for(k = start_idx[2]; k < stop_idx[2]; k++){
            euc[2] = fabs(c_pos[2] - z);

            euc_dist = sqrt(SQR2(euc[0]) + SQR2(euc[1]) + SQR2(euc[2]));
            if(euc_dist <= regrid_radius){

               /* calculate the weighting factor */
               switch (regrid_type){
               default:
                  fprintf(stderr, "Erk! unknown regrid_type. File: %s Line: %d\n",
                          __FILE__, __LINE__);
                  exit(EXIT_FAILURE);
                  break;

               case NEAREST_FUNC:
                  fprintf(stderr, "Erk! -nearest is not implemented yet\n");
                  exit(EXIT_FAILURE);
                  break;

               case KAISERBESSEL_FUNC:
                  weight =
                     gsl_sf_bessel_I0(regrid_sigma *
                                      sqrt(1 - SQR2(euc[0] / regrid_radius))) *
                     gsl_sf_bessel_I0(regrid_sigma *
                                      sqrt(1 - SQR2(euc[1] / regrid_radius))) *
                     gsl_sf_bessel_I0(regrid_sigma *
                                      sqrt(1 - SQR2(euc[2] / regrid_radius))) /
                     SQR3(regrid_radius);

                  break;

               case GAUSSIAN_FUNC:
                  weight = exp(-SQR2(euc_dist) / SQR2(regrid_sigma));
                  break;
                  }

               /* set data values */
               for(v = 0; v < v_size; v++){
                  value = get_volume_real_value(*totals, k, j, i, v, 0);
                  set_volume_real_value(*totals, k, j, i, v, 0,
                                        value + (data_buf[0 + v] * weight));
                  }

               /* increment count value */
               value = get_volume_real_value(*counts, k, j, i, 0, 0);
               set_volume_real_value(*counts, k, j, i, 0, 0, value + weight);
               }

            c_pos[2] += steps[perm[2]];
            }

         c_pos[1] += steps[perm[1]];
         }

      c_pos[0] += steps[perm[0]];
      }

   }

/* regrid a volume with an arbitrary path     */
/* return resulting totals and counts volumes */
void regrid_arb_path(char *coord_fn, char *data_fn, int buff_size,
                     Volume * totals, Volume * counts, int v_size)
{
   Coord_list coord_buf;
   double  *data_buf = NULL;
   size_t   data_buf_alloc_size = 0;
   int      sizes[MAX_VAR_DIMS];
   int      c;
   int      total_pts;

   /* get volume info */
   get_volume_sizes(*totals, sizes);

   /* check for the config file */
   if(!file_exists(coord_fn)){
      fprintf(stderr, "Couldn't find config file %s.\n\n", coord_fn);
      exit(EXIT_FAILURE);
      }

   /* initialise the parser with the config file */
   if(!init_arb_path(coord_fn, data_fn)){
      fprintf(stderr, "Failed to init arb_path, this isn't good\n");
      exit(EXIT_FAILURE);
      }

   /* get some co-ordinates */
   if(verbose){
      fprintf(stdout, " + Doing arbitrary path (vector: %d)\n", v_size);
      }
   total_pts = 0;
   coord_buf = get_some_arb_path_coords(buff_size);
   while(coord_buf->n_pts != 0){

      /* grow data_buf if we have to */
      if(coord_buf->n_pts * v_size > data_buf_alloc_size){
         data_buf_alloc_size = coord_buf->n_pts * v_size;
         data_buf = realloc(data_buf, data_buf_alloc_size * sizeof(*data_buf));
         }

      /* get the data */
      if(!get_some_arb_path_data
         (data_buf, in_dtype, in_is_signed, coord_buf->n_pts, v_size)){
         fprintf(stderr, "failed getting data\n");
         exit(EXIT_FAILURE);
         }

      total_pts += coord_buf->n_pts;
      if(verbose){
         fprintf(stdout, " | got %d co-ords   total: %d\n", coord_buf->n_pts, total_pts);
         }

      /* regrid (do the nasty) */
      for(c = 0; c < coord_buf->n_pts; c++){

         regrid_point(totals, counts,
                      coord_buf->pts[c].coord[0],
                      coord_buf->pts[c].coord[1],
                      coord_buf->pts[c].coord[2], v_size, &data_buf[c * v_size]);
         }

      /* get the next lot of co-ordinates */
      coord_buf = get_some_arb_path_coords(buff_size);
      }

   /* finish up */
   end_arb_path();
   }

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
   while((ch = getc(fp)) != EOF){
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
                  strncpy(args[nargs], tmp, (size_t) ichar);

                  ichar = 0;
                  nargs++;
                  }
               }
            else {
               tmp[ichar++] = ch;
               }
            }
         }
      }
   fclose(fp);
   return nargs;
   }

/* get info from a model file (steps, starts, etc */
int get_model_file_info(char *dst, char *key, char *nextArg)
{
   int      mincid;
   char    *fname;

   /* Get pointer to volume definition structure */
   Volume_Definition *vd = (Volume_Definition *) dst;

   /* Check for following argument */
   if(nextArg == NULL){
      (void)fprintf(stderr, "\"%s\" option requires an additional argument\n", key);
      exit(EXIT_FAILURE);
      }
   fname = nextArg;

   /* check that the file exists */
   if(!file_exists(fname)){
      fprintf(stderr, "Couldn't find -like file %s.\n\n", fname);
      exit(EXIT_FAILURE);
      }

   /* get volume definition from input filename */
   mincid = miopen(fname, NC_NOWRITE);
   get_vol_def(mincid, vd);

   /* Close the file */
   miclose(mincid);

   return TRUE;
   }
