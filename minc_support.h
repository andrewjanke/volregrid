/* minc_support.h  - these are functions that should be in MINC itself */

#ifndef MINC_SUPPORT_H
#define MINC_SUPPORT_H

#include <netcdf.h>
#include <minc.h>

#define WORLD_NDIMS 3
#define NO_AXIS -1

/* structure for volume info -- essentially the good bits of volume_io */
typedef struct {
   int ndims;
   
   char *dimname[WORLD_NDIMS];
   
   /* permutation arrays from dimension order */
   /* to spatial order and vice-versa         */
   int space_to_dim[WORLD_NDIMS];
   int dim_to_space[MAX_VAR_DIMS];
   
   /* These are subscripted by X, Y and Z */
   int      nelem[WORLD_NDIMS];
   double   start[WORLD_NDIMS];
   double   step[WORLD_NDIMS];
   double   dircos[WORLD_NDIMS][WORLD_NDIMS];
   
   /* tranformation from voxel to world co-ordinates */
   double voxel_to_world[WORLD_NDIMS][WORLD_NDIMS + 1];
   } Volume_Definition;

/* get various dimension info -- nicked from mincstats */
void     get_minc_attribute(int mincid, char *varname, char *attname,
                            int maxvals, double vals[]);
int      get_minc_ndims(int mincid);
void     get_minc_spatial_dims(int mincid, int space_to_dim[], int dim_to_space[]);
void     get_minc_voxel_to_world(int mincid,
                                 double voxel_to_world[WORLD_NDIMS][WORLD_NDIMS + 1]);
void     normalize_vector(double vector[]);
void     transform_coord(double out_coord[],
                         double transform[WORLD_NDIMS][WORLD_NDIMS + 1],
                         double in_coord[]);
void     get_vol_def(int mincid, Volume_Definition * vd);

void     fprintf_vol_def(FILE *fp, Volume_Definition *vd);

#endif
