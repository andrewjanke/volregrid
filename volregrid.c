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

int      verbose = FALSE;
int      debug   = FALSE;
int      clobber = FALSE;
int      regrid_dim = 3;
char    *arb_path_cfg_fn = NULL;
nc_type  dtype = NC_FLOAT;
int      is_signed = FALSE;

static ArgvInfo argTable[] = {
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "Print out extra information."},
   {"-debug", ARGV_CONSTANT, (char *)TRUE, (char *)&debug,
    "Spew copious amounts of debugging info."},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "Clobber existing files."},

   {NULL, ARGV_HELP, NULL, NULL, "\nOutfile Options"},
   {"-byte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&dtype,
    "Write out byte data."},
   {"-short", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&dtype,
    "Write out short integer data."},
   {"-long", ARGV_CONSTANT, (char *)NC_LONG, (char *)&dtype,
    "Write out long integer data."},
   {"-float", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&dtype,
    "Write out single-precision data. (Default)"},
   {"-double", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&dtype,
    "Write out double-precision data."},
   {"-signed", ARGV_CONSTANT, (char *)TRUE, (char *)&is_signed,
    "Write signed integer data."},
   {"-unsigned", ARGV_CONSTANT, (char *)FALSE, (char *)&is_signed,
    "Write unsigned integer data."},

   {NULL, ARGV_HELP, NULL, NULL, "\nRegridding options"},
   {"-2D", ARGV_CONSTANT, (char *)2, (char *)&regrid_dim,
    "Regrid slice by slice (Default 3D)."},

   {"-arb_path_cfg", ARGV_STRING, (char *)1, (char *)&arb_path_cfg_fn,
    "Regrid data using an arbitrary path from the input filename"},


   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

char    *output_dimorder[] = { MIzspace, MIyspace, MIxspace };
char    *output_dimorder_v[] = { MIzspace, MIyspace, MIxspace, MIvector_dimension };

main(int argc, char *argv[])
{
   char    *in_fn, *out_fn;
   char    *history;
   Status   status;
   Volume   out_vol;
   int      c;
   
   Coord_list tmp;

//   Real     min, max;

//   minc_input_options in_ops;

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
   
   /* arbitrary path */
   if(arb_path_cfg_fn != NULL){
      
      /* check for the config file */
      if(!file_exists(arb_path_cfg_fn)){
         fprintf(stderr, "%s: Couldn't find config file %s.\n\n", argv[0], arb_path_cfg_fn);
         exit(EXIT_FAILURE);
         }
      
      /* initialise the parser with the config file */
      if(!init_arb_path(arb_path_cfg_fn)){
         fprintf(stderr, "%s: Failed to init arb_path, this isn't good\n", argv[0]);
         exit(EXIT_FAILURE);
         }
   
      /* get some data */
      tmp = get_some_arb_path_coords(100);
      printf("CO-ORDS\n");
      print_coord_list(tmp);
      
   //   tmp = get_some_arb_path_coords(10);
      
   //   print_coord_list(tmp);
      
      end_arb_path();
      
      }
   exit(0);

//   if(status != OK){
//      fprintf(stderr, "Problems reading: %s\n", in_fn);
//      exit(EXIT_FAILURE);
//      }

//   if(verbose){
//      Real     min_value, max_value;
//
//      get_volume_real_range(data, &min_value, &max_value);
//
//      fprintf(stdout, " | Input file:     %s\n", in_fn);
//      fprintf(stdout, " | Input ndims:    %d\n", in_ndims);
//      fprintf(stdout, " | min/max:        [%8.3f:%8.3f]\n", min_value, max_value);
//      fprintf(stdout, " | Output files:\n");
//      for(c = 0; c < MAX_OUTFILES; c++){
//         if(outfiles[c] != NULL){
//            fprintf(stdout, " |   [%d]:         %s => %s\n", c, out_names[c],
//                    outfiles[c]);
//            }
//         }
//      fprintf(stdout, " | FFT order:      %d\n", fft_dim);
//      }


   /* create the output volume */
   {
      Real     min, max;

      int      sizes[4];
      Real     starts[4];
      Real     separations[4];

      for(c = 0; c < 3; c++){
         sizes[c] = 10;
         starts[c] = -5;
         separations[c] = 1;
         }

      /* setup frequency dimension */
      sizes[3] = 2;
      starts[3] = 0;
      separations[3] = 1;

      min = 0;
      max = 10;

      /* define new out_vol volume  */
      out_vol = create_volume(4, output_dimorder_v, NC_FLOAT, TRUE, 0.0, 0.0);
      set_volume_sizes(out_vol, sizes);
      set_volume_starts(out_vol, starts);
      set_volume_separations(out_vol, separations);
      set_volume_real_range(out_vol, min, max);

      /* allocate space for out_vol */
      alloc_volume_data(out_vol);
      }

   /* regrid (do the nasty) */
//   if(regrid_arbitrary_path(out_vol, in_fn, path_data) != OK){
//      print_error("%s: Problems regridding data\n", argv[0]);
//      exit(EXIT_FAILURE);
//      }

   /* output the result */
   if(verbose){
      fprintf(stdout, "Outputting %s...\n", out_fn);
      }
   status = output_volume(out_fn, dtype, is_signed, 0, 0, out_vol, history, NULL);
   if(status != OK){
      print_error("Problems outputing: %s", out_fn);
      }

   delete_volume(out_vol);
   return (status);
   }


/* regrid a volume with an arbitrary path */
