/* minc_support.c  -  these are functions that should be in MINC itself*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minc_support.h"

/* Get the total number of image dimensions in a minc file */
int get_minc_ndims(int mincid)
{
   int      imgid;
   int      ndims;

   imgid = ncvarid(mincid, MIimage);
   (void)ncvarinq(mincid, imgid, NULL, NULL, &ndims, NULL, NULL);

   return ndims;
   }

/* Get a double attribute from a minc file */
void get_minc_attribute(int mincid, char *varname, char *attname,
                        int maxvals, double vals[])
{
   int      varid;
   int      old_ncopts;
   int      att_length;

   if(!mivar_exists(mincid, varname)){
      return;
      }
   varid = ncvarid(mincid, varname);
   old_ncopts = ncopts;
   ncopts = 0;
   (void)miattget(mincid, varid, attname, NC_DOUBLE, maxvals, vals, &att_length);
   ncopts = old_ncopts;
   }

/* Get the permutation from spatial (x, y, z) - to file dimensions
   and vice-versa. */
void get_minc_spatial_dims(int mincid, int space_to_dim[], int dim_to_space[])
{
   int      imgid, dim[MAX_VAR_DIMS];
   int      idim, ndims, world_index;
   char     dimname[MAX_NC_NAME];

   /* Set default values */
   for(idim = 0; idim < 3; idim++){
      space_to_dim[idim] = -1;
      }
   for(idim = 0; idim < MAX_VAR_DIMS; idim++){
      dim_to_space[idim] = -1;
      }

   /* Get the dimension ids for the image variable */
   imgid = ncvarid(mincid, MIimage);
   (void)ncvarinq(mincid, imgid, NULL, NULL, &ndims, dim, NULL);

   /* Loop over them to find the spatial ones */
   for(idim = 0; idim < ndims; idim++){

      /* Get the name and check that this is a spatial dimension */
      (void)ncdiminq(mincid, dim[idim], dimname, NULL);
      if((dimname[0] == '\0') || (strcmp(&dimname[1], "space") != 0)){
         continue;
         }

      /* Look for the spatial dimensions */
      switch (dimname[0]){
      case 'x':
         world_index = 0;
         break;
      case 'y':
         world_index = 1;
         break;
      case 'z':
         world_index = 2;
         break;
      default:
         world_index = 0;
         break;
         }
      space_to_dim[world_index] = idim;
      dim_to_space[idim] = world_index;
      }
   }

/* Get the voxel to world transform (for column vectors) */
void get_minc_voxel_to_world(int mincid,
                             double voxel_to_world[WORLD_NDIMS][WORLD_NDIMS + 1])
{
   int      idim, jdim;
   double   dircos[WORLD_NDIMS];
   double   step, start;
   char    *dimensions[WORLD_NDIMS];

   dimensions[0] = MIxspace;
   dimensions[1] = MIyspace;
   dimensions[2] = MIzspace;

   /* Zero the matrix */
   for(idim = 0; idim < WORLD_NDIMS; idim++){
      for(jdim = 0; jdim < WORLD_NDIMS + 1; jdim++)
         voxel_to_world[idim][jdim] = 0.0;
      }

   for(jdim = 0; jdim < WORLD_NDIMS; jdim++){

      /* Set default values */
      step = 1.0;
      start = 0.0;
      for(idim = 0; idim < WORLD_NDIMS; idim++)
         dircos[idim] = 0.0;
      dircos[jdim] = 1.0;

      /* Get the attributes */
      get_minc_attribute(mincid, dimensions[jdim], MIstart, 1, &start);
      get_minc_attribute(mincid, dimensions[jdim], MIstep, 1, &step);
      get_minc_attribute(mincid, dimensions[jdim], MIdirection_cosines,
                         WORLD_NDIMS, dircos);
      normalize_vector(dircos);

      /* Put them in the matrix */
      for(idim = 0; idim < WORLD_NDIMS; idim++){
         voxel_to_world[idim][jdim] = step * dircos[idim];
         voxel_to_world[idim][WORLD_NDIMS] += start * dircos[idim];
         }
      }
   }

void normalize_vector(double vector[])
{
   int      idim;
   double   magnitude;

   magnitude = 0.0;
   for(idim = 0; idim < WORLD_NDIMS; idim++){
      magnitude += (vector[idim] * vector[idim]);
      }
   magnitude = sqrt(magnitude);
   if(magnitude > 0.0){
      for(idim = 0; idim < WORLD_NDIMS; idim++){
         vector[idim] /= magnitude;
         }
      }
   }

/* Transforms a coordinate through a linear transform */
void transform_coord(double out_coord[],
                     double transform[WORLD_NDIMS][WORLD_NDIMS + 1], double in_coord[])
{
   int      idim, jdim;
   double   homogeneous_coord[WORLD_NDIMS + 1];

   for(idim = 0; idim < WORLD_NDIMS; idim++){
      homogeneous_coord[idim] = in_coord[idim];
      }
   homogeneous_coord[WORLD_NDIMS] = 1.0;

   for(idim = 0; idim < WORLD_NDIMS; idim++){
      out_coord[idim] = 0.0;
      for(jdim = 0; jdim < WORLD_NDIMS + 1; jdim++){
         out_coord[idim] += transform[idim][jdim] * homogeneous_coord[jdim];
         }
      }
   }

/* get volume definition from a mincid */
void get_vol_def(int mincid, Volume_Definition * vd)
{
   int      i, j;
   long     tmp_long;

   /* number of dimensions */
   vd->ndims = get_minc_ndims(mincid);

   /* fill in permutation arrays */
   get_minc_spatial_dims(mincid, vd->space_to_dim, vd->dim_to_space);

   /* get nelems, start, step and dircos */
   for(j = 0; j < WORLD_NDIMS; j++){

      /* Set default values */
      vd->step[j] = 1.0;
      vd->start[j] = 0.0;
      for(i = 0; i < WORLD_NDIMS; i++){
         vd->dircos[j][i] = 0.0;
         }
      vd->dircos[j][j] = 1.0;

      /* dimension size */
      (void)ncdiminq(mincid, vd->dim_to_space[j], NULL, &tmp_long);
      vd->nelem[j] = (int)tmp_long;

      /* start, step, dircos */
      get_minc_attribute(mincid, vd->dimname[j], MIstart, 1, &(vd->start[j]));
      get_minc_attribute(mincid, vd->dimname[j], MIstep, 1, &(vd->step[j]));
      get_minc_attribute(mincid, vd->dimname[j], MIdirection_cosines, WORLD_NDIMS,
                         vd->dircos[j]);
      normalize_vector(vd->dircos[j]);
      }

   /* voxel to world transform */
   for(i = 0; i < WORLD_NDIMS; i++){
      for(j = 0; j < WORLD_NDIMS + 1; j++)
         vd->voxel_to_world[i][j] = 0.0;
      }

   /* Put them in the matrix */
   for(j = 0; j < WORLD_NDIMS; j++){
      for(i = 0; i < WORLD_NDIMS; i++){
         vd->voxel_to_world[i][j] = vd->step[j] * vd->dircos[j][i];
         vd->voxel_to_world[i][WORLD_NDIMS] += vd->start[j] * vd->dircos[j][i];
         }
      }
   }

void fprintf_vol_def(FILE * fp, Volume_Definition * vd)
{
   int      i, idx;

   fprintf(fp, "+ ndims: %d\n", vd->ndims);
   
   fprintf(fp, "| The following are in file dimension order\n");
   
   fprintf(fp, "|  dimname s2d d2s  nelem  start      step        dircos\n");
   for(i = 0; i < vd->ndims; i++){
      idx = vd->dim_to_space[i];
      
      if(idx != -1){
         fprintf(fp, "   %s   %d   %d    %d  %6f  %6f  %6f %6f %6f\n", 
                 vd->dimname[idx],
                 vd->space_to_dim[idx], vd->dim_to_space[idx],
                 vd->nelem[idx], vd->start[idx], vd->step[idx], 
                 vd->dircos[idx][0], vd->dircos[idx][1], vd->dircos[idx][2]);
         }
      else{
         fprintf(fp, "  non-spatial dimension --- -- - - - - - - - --\n");
         }

      }
   }
