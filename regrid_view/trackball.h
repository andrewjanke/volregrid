/* (c) Copyright 1993, 1994, Silicon Graphics, Inc.                          */
/* ALL RIGHTS RESERVED                                                       */
/* Permission to use, copy, modify, and distribute this software for         */
/* any purpose and without fee is hereby granted, provided that the above    */
/* copyright notice appear in all copies and that both the copyright notice  */
/* and this permission notice appear in supporting documentation, and that   */
/* the name of Silicon Graphics, Inc. not be used in advertising             */
/* or publicity pertaining to distribution of the software without specific, */
/* written prior permission.                                                 */
/*                                                                           */
/* THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"         */
/* AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,          */
/* INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR          */
/* FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON              */
/* GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,           */
/* SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY             */
/* KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,            */
/* LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF         */
/* THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN            */
/* ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON            */
/* ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE         */
/* POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.                          */
/*                                                                           */
/* US Government Users Restricted Rights                                     */
/* Use, duplication, or disclosure by the Government is subject to           */
/* restrictions set forth in FAR 52.227.19(c)(2) or subparagraph             */
/* (c)(1)(ii) of the Rights in Technical Data and Computer Software          */
/* clause at DFARS 252.227-7013 and/or in similar or successor               */
/* clauses in the FAR or the DOD or NASA FAR Supplement.                     */
/* Unpublished-- rights reserved under the copyright laws of the             */
/* United States.  Contractor/manufacturer is Silicon Graphics,              */
/* Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.             */
/*                                                                           */
/* OpenGL(TM) is a trademark of Silicon Graphics, Inc.                       */

/* trackball.h                                                               */
/* A virtual trackball implementation                                        */
/* Written by Gavin Bell for Silicon Graphics, November 1988.                */
/* Modified by Andrew Janke - 2001 - quat_to_axis, vects_to_quat             */


/* Pass the x and y coordinates of the last and current positions of         */
/* the mouse, scaled so they are from (-1.0 ... 1.0).                        */
/*                                                                           */
/* The resulting rotation is returned as a quaternion rotation in the        */
/* first paramater.                                                          */
void trackball(double q[4], double p1x, double p1y, double p2x, double p2y);

/* Given two quaternions, add them together to get a third quaternion.       */
/* Adding quaternions to get a compound rotation is analagous to adding      */
/* translations to get a compound translation.  When incrementally           */
/* adding rotations, the first argument here should be the new               */
/* rotation, the second and third the total rotation (which will be          */
/* over-written with the resulting new total rotation).                      */
void add_quats(double *q1, double *q2, double *dest);

/* Builds a rotation matrix based on a given quaternion.                     */
void build_rotmatrix(double m[4][4], double quat[4]);

/* This function computes a quaternion based on an axis (defined by          */
/* the given vector) and an angle about which to rotate.  The angle is       */
/* expressed in radians.  The result is put into the third argument.         */
void axis_to_quat(double vec[3], double  phi, double quat[4]);

/* extra functions  -- Andrew Janke                                          */
void quat_to_axis(double vec[3], double *phi, double quat[4]);
void quat_to_xvec(double vec[3], double quat[4]);
void quat_to_yvec(double vec[3], double quat[4]);
void quat_to_zvec(double vec[3], double quat[4]);
/* creates a quat from 2 up vectors                                          */
void vects_to_quat(double v1[3], double v2[3], double q[4]);


#define SQR(x) ((x)*(x))

/* inits a vector                                                   */
void vset(double *v, double x, double y, double z);

/* copy a vector                                                    */
void vcopy(double *copy, double *v);

/* compute the addition of two vectors                              */
void vadd(double *add, double *v1, double *v2);

/* compute the subtraction of two vectors (v1 - v2)                 */
void vsub(double *sub, double *v1, double *v2);

/* multiply a vector by a constant                                  */
void vscale(double *v, double scale);

/* returns the length of a vector                                   */
double vlength(double *v);

/* returns the euc distance between two points                      */
double veuc(double *p1, double *p2);

/* normalise a vector                                               */
void vnormal(double *v);

/* computes the vector cross product of two vectors                 */
void vcross(double *cross, double *v1, double *v2);
 
/* returns the vector dot product of two vectors                    */
double vdot(double *v1, double *v2);

