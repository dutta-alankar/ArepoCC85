/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/mesh/criterion_refinement.c
 * \date        05/2018
 * \brief       Criteria for the refinement of a cell.
 * \details     Routines which are checking whether a cell should be refined.
 *              contains functions:
 *                int should_this_cell_be_split(int i)
 *                static int can_this_cell_be_split(int i)
 *                static int refine_criterion_default(int i)
 *                static int refine_criterion_jeans_ref(int i)
 *                static int jeans_refinement_criteria(int i)
 *                static int refine_criterion_volume(int i)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#if defined(REFINEMENT_SPLIT_CELLS) && !defined(ONEDIMS)
static int can_this_cell_be_split(int i);
static int refine_criterion_default(int i);
static int refine_criterion_jeans_ref(int i);
static int jeans_refinement_criteria(int i);

#ifdef REFINEMENT_VOLUME_LIMIT
static int refine_criterion_volume(int i);
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */

#ifdef REFINEMENT_MERGE_CELLS
char *FlagDoNotRefine;
#endif /* #ifdef REFINEMENT_MERGE_CELLS */

/*! \brief Should this cell be refined?
 *
 *  This function signals whether a cell needs further refinement. This needs
 *  to be adjusted according to the needs of the simulation in question.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be split.
 */
int should_this_cell_be_split(int i)
{
#ifdef REFINEMENT_MERGE_CELLS
  if(FlagDoNotRefine[i])
    return 0;
#endif /* #ifdef REFINEMENT_MERGE_CELLS */

  if(P[i].Mass == 0 && P[i].ID == 0) /* skip cells that have been swallowed or dissolved */
    return 0;

#ifdef AGNWIND_FLAG
  if(SphP[i].AGNFlag == 1 || SphP[i].AGNFlag == 2)
    return 0;
#endif /* #ifdef AGNWIND_FLAG */

#ifdef STICKYFLAGS
  if(SphP[i].StickyFlag > 0)
    return 0;
#endif /* #ifdef STICKYFLAGS */

#if defined(REFINEMENT_VOLUME_LIMIT)
  double maxvolume = All.MaxVolume;
  double minvolume = All.MinVolume;

  if(SphP[i].Volume > 2. * maxvolume)
    if(can_this_cell_be_split(i))
      return 1;

  if(SphP[i].Volume < 2. * minvolume)
    return 0;

  if(refine_criterion_volume(i))
    if(can_this_cell_be_split(i))
      return 1;
#endif /* #if defined(REFINEMENT_VOLUME_LIMIT) */

  switch(All.RefinementCriterion) /* select the function that evaluates the refinement criterion */
    {
      case 0:
        return 0;
        break;

      case 1:
        return refine_criterion_default(i);
        break;

      case 2:
        return refine_criterion_jeans_ref(i);
        break;

      default:
        terminate("invalid refinement criterion specified");
        break;
    }

  return 0;
}

/*
 * static functions; i.e. functions that are only called within this file
 */

/*! \brief Is cell round enough to be refined?
 *
 *  This function signals whether a cell is allowed refinement. A cell that
 *  is supposed to be refined needs to match certain roundness criteria, which
 *  are specified in this function.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell is allowed to be refined.
 */
static int can_this_cell_be_split(int i)
{
#ifdef REGULARIZE_MESH_FACE_ANGLE
  if(SphP[i].MaxFaceAngle < 1.5 * All.CellMaxAngleFactor)
    return 1;

#else  /* #ifdef REGULARIZE_MESH_FACE_ANGLE */
  double dx      = nearest_x(P[i].Pos[0] - SphP[i].Center[0]);
  double dy      = nearest_y(P[i].Pos[1] - SphP[i].Center[1]);
  double dz      = nearest_z(P[i].Pos[2] - SphP[i].Center[2]);
  double d       = sqrt(dx * dx + dy * dy + dz * dz);
  double cellrad = get_cell_radius(i);

  if(d < 2.0 * All.CellShapingFactor * cellrad) /* only refine cells which are reasonably 'round' */
    return 1;
#endif /* #ifdef REGULARIZE_MESH_FACE_ANGLE #else */

  return 0;
}

/*! \brief Default refinement criterion.
 *
 *  Checks if cell is within a factor of 2 of the target gas mass.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be refined.
 */
static int refine_criterion_default(int i)
{
#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement != 0) {
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */
#ifdef AGNWIND_FLAG
  locate_AGN_sphere();
  double xpos, ypos, zpos;
  double r, theta, phi;

  xpos = P[i].Pos[0] - SpherePosX;
  ypos = P[i].Pos[1] - SpherePosY;
  zpos = P[i].Pos[2] - SpherePosZ;

  r = sqrt(xpos*xpos+ypos*ypos+zpos*zpos);
  theta = acos(zpos/r);
  phi = atan2(ypos,xpos);
  phi = (phi<0)?(phi+2*M_PI):phi;

  /* Assumes cloud is along positive y-axis */
  int Nprp, Nprl;
  double delTheta, delPhi;

  Nprp = 10;
  Nprl = 300;

  /* For now check ic-gen-relax1-cloud.py */
  double Rcl = 10.;
  double Rini = 300.;

  double Rmax = Nprl*Rcl;

  double enhance = All.EnhanceResFactor;
#ifdef CLOUD_PRESENT
  double enhanceCloud = All.EnhanceResCloudNeigh;
  double neighColor = SphP[i].AGNWindFraction;
  if ( (neighColor>1.e-3) || ( (ypos>0 && ypos<20*Rcl) && (fabs(xpos)<=4*Rcl) && (fabs(zpos)<=4*Rcl)) ){
    if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * enhanceCloud *All.TargetGasMass)
          return 1;
    else return 0;
  }
#endif /* #ifdef CLOUD_PRESENT */

  delTheta = Nprp*(Rcl/Rini);
  delPhi = delTheta;

  // if(r>1.5*Rmax) return 0;
  if (r<=Rmax){
    if ( (theta>=(M_PI/2.-delTheta)) && (theta<=(M_PI/2.+delTheta)) ){
      if ( (phi>=(M_PI/2.-delPhi)) && (phi<=(M_PI/2.+delPhi)) ){
#if defined(REFINE_AGNWIND)
        if(can_this_cell_be_split(i) && SphP[i].AGNFlag == 3) {
          if(P[i].Mass > 2.0 * enhance * All.TargetGasMass/All.ResolutionBoostAGNWind)
            return 1;
        }
        else return 0;
#else /* #if defined(REFINE_AGNWIND) */
        if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * enhance *All.TargetGasMass)
          return 1;
        else return 0;
#endif /* #if defined(REFINE_AGNWIND) #else */
      }
    }
  }
  else {
#if defined(REFINE_AGNWIND)
    if(can_this_cell_be_split(i) && SphP[i].AGNFlag == 3) {
       if(P[i].Mass > 2.0 * All.TargetGasMass/All.ResolutionBoostAGNWind)
         return 1;
    }
  }
#endif /* #if defined(REFINE_AGNWIND) */
#endif /* #ifdef AGNWIND_FLAG */
    if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * All.TargetGasMass)
      return 1;
#ifdef REFINEMENT_HIGH_RES_GAS
  }
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */
  return 0; /* default is not to refine */
}

/*! \brief Jeans refinement criterion additional target mass criterion
 *
 *  Resolving the Jeans length is an additional criterion, apart from obeying
 *  the usual factor of 2 within a target mass criterion.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be refined.
 */
static int refine_criterion_jeans_ref(int i)
{
#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement != 0)
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */
    if(can_this_cell_be_split(i))
      {
        if(P[i].Mass > 2.0 * All.TargetGasMass)
          return 1;

#ifdef JEANS_REFINEMENT
        return jeans_refinement_criteria(i);
#else  /* #ifdef JEANS_REFINEMENT */
      return 0;
#endif /* #ifdef JEANS_REFINEMENT #else */
      }

  return 0;
}

/*! \brief Refinement criterion according to Jeans stability of a cell.
 *
 *  The cell will be refined if the Jeans length is not resolved by
 *  JEANS_REFINEMENT cells.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be refined.
 */
static int jeans_refinement_criteria(int i)
{
#ifdef JEANS_REFINEMENT
  if(can_this_cell_be_split(i))
    {
      double jeans_number, jeans_length, sound_speed, dx;

      sound_speed  = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
      jeans_length = sqrt(M_PI / All.G / SphP[i].Density) * sound_speed;
      dx           = 2.0 * get_cell_radius(i);
      jeans_number = jeans_length / dx;

      if(jeans_number < JEANS_REFINEMENT)
        {
          return 1;
        }
    }
#endif /* #ifdef JEANS_REFINEMENT */

  return 0;
}

#ifdef REFINEMENT_VOLUME_LIMIT
/*! \brief Refinement criterion for based on the minimum volume of a
 *  neighboring cell.
 *
 *  This criterion is supposed to avoid sudden jumps in resolution which lead
 *  to an inaccurate result. Each cell that has a volume larger than a
 *  specified factor times the minimum volume of all neighboring cells will be
 *  refined. This also includes a global absolute minimum and maximum volume.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be refined.
 */
static int refine_criterion_volume(int i)
{
#ifdef AGNWIND_FLAG
  locate_AGN_sphere();
  double xpos, ypos, zpos;
  double r, theta, phi;

  xpos = P[i].Pos[0] - SpherePosX;
  ypos = P[i].Pos[1] - SpherePosY;
  zpos = P[i].Pos[2] - SpherePosZ;

  r = sqrt(xpos*xpos+ypos*ypos+zpos*zpos);
  theta = acos(zpos/r);
  phi = atan2(ypos,xpos);
  phi = (phi<0)?(phi+2*M_PI):phi;

  double rHres = 1.5*(All.AGNWindSphereRad*PARSEC/All.UnitLength_in_cm);
  double enhance = 1.; //(r<=rHres)?(2./3):1.0;
  if(All.MaxVolumeDiff > 0 && SphP[i].Volume > enhance * All.MaxVolumeDiff * SphP[i].MinNgbVolume)
#else /* #ifdef AGNWIND_FLAG */
  if(All.MaxVolumeDiff > 0 && SphP[i].Volume > All.MaxVolumeDiff * SphP[i].MinNgbVolume)
#endif /* #ifdef AGNWIND_FLAG #else */
    {
#ifdef REGULARIZE_MESH_FACE_ANGLE
      if(SphP[i].MaxFaceAngle < 1.5 * All.CellMaxAngleFactor)
        return 1;
#else  /* #ifdef REGULARIZE_MESH_FACE_ANGLE */

      double dx      = nearest_x(P[i].Pos[0] - SphP[i].Center[0]);
      double dy      = nearest_y(P[i].Pos[1] - SphP[i].Center[1]);
      double dz      = nearest_z(P[i].Pos[2] - SphP[i].Center[2]);
      double d       = sqrt(dx * dx + dy * dy + dz * dz);
      double cellrad = get_cell_radius(i);

      if(d < 2.0 * All.CellShapingFactor * cellrad) /* only refine cells which are reasonably 'round' */
        return 1;
#endif /* #ifdef REGULARIZE_MESH_FACE_ANGLE #else */
    }

  return 0;
}
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */

#endif /* #if defined(REFINEMENT_SPLIT_CELLS) && !defined(ONEDIMS) */
