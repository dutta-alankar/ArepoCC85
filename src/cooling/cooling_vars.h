/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              Alankar Dutta (alankardutta@IISC.AC.IN).
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
 * \file        src/cooling/cooling_vars.h
 * \date        05/2021
 * \brief       Variables for cooling.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 27.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#define MAX_TABLESIZE 20000 /* Max # of lines in TREECOOL */
#define MAX_SUBCYCLE 50 /* Max number of sub-cycling steps */

/* data for gas state */
typedef struct
{
  double XH, YHe, Zmet;
  double ethmin; /* minimum internal energy for neutral gas */
  double mu, mue, mui;
} GasState;

/* tabulated rates */
typedef struct
{
  double Temperature[MAX_TABLESIZE], LAMBDA[MAX_TABLESIZE], townY[MAX_TABLESIZE];
  double townY_inv[MAX_TABLESIZE], Temperature_invY[MAX_TABLESIZE];
} RateTable;

/* cooling data */
typedef struct
{
  double u_old_input, rho_input, dt_input;
} DoCoolData;
