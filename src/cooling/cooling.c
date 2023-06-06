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
 * \file        src/cooling/cooling.c
 * \date        05/2018
 * \brief       Module for gas radiative cooling
 * \details     contains functions:
 *                double DoCooling(double u_old, double rho, double dt)
 *                double convert_u_to_temp(double u, double rho)
 *                double CoolingRateFromU(double u, double rho)
 *                double CoolingRate(double T, double rho)
 *                void MakeRateTable(char *fname)
 *                void InitCool(void)
 *                void cooling_only(void)
 *                void cool_cell(int i)
 *                double GetCoolingTime (double u, double rho)
 *                double lambda_interp(double temperature)
 *                double Y_interp(double temperature)
 *                double Y_interp(double temperature)
 *                double invY_interp(double townY)
 *                double interp1D(double* x_data, double* y_data, double x_request, char* msg)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef COOLING

static double Tmin = 10.0;     /*!< min temperature */
//static double Tmax = 1.e9;     /*!< max temperature */
static GasState gs;           /*!< gas state */
static RateTable *RateT;      /*!< tabulated rates */
static int NCOOLTAB = MAX_TABLESIZE;

/*! \brief Computes the new internal energy per unit mass.
 *
 *  The function solves for the new internal energy per unit mass of the gas
 *  by integrating the equation for the internal energy with an implicit
 *  Euler scheme. The root of resulting non linear equation,
 *  which gives tnew internal energy, is found with the bisection method.
 *  Arguments are passed in internal units.
 *
 *  \param[in] u_old the initial (before cooling is applied) internal energy
 *             per unit mass of the gas cell.
 *  \param[in] rho   the proper density of the gas cell.
 *  \param[in] dt    the duration of the time step.
 *  \param[in] ne_guess electron number density relative to hydrogen number
 *             density (for molecular weight computation).
 *
 *  \return The new internal energy per unit mass of the gas cell.
 */
double DoCooling(double u_old, double rho, double dt)
{
  double u, dt_cool;

  if(!gsl_finite(u_old))
    terminate("invalid input: u_old=%g\n", u_old);

  if(u_old < 0 || rho < 0)
    terminate("invalid input: task=%d u_old=%g  rho=%g  dt=%g  All.MinEgySpec=%g\n", ThisTask, u_old, rho, dt, All.MinEgySpec);

  dt_cool = GetCoolingTime(u_old, rho)*All.UnitTime_in_s / All.HubbleParam; /* convert to physical cgs units */
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;

  /* -----------------------------------------------------
   *     Subcycling starts here
   * ----------------------------------------------------- */
  u = u_old;
  double dt_sub = 1.0/(1.0/(dt_cool/MAX_SUBCYCLE) + 2./dt); /* Subcycle limit */
  int sub_steps = (int)ceil(dt/dt_sub);
  dt_sub = dt/sub_steps;
  int n_sub = 0;
  double prs, tcool;
  double T0, T1, Tref = 1.e9;
  for (n_sub = 0; n_sub<sub_steps; n_sub++ )
  {
    prs = rho*u*GAMMA_MINUS1;
    T0  = convert_u_to_temp (u, rho);

    tcool = GetCoolingTime( u/(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs),
	     rho/(All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam) )*All.UnitTime_in_s / All.HubbleParam;

    T1 = invY_interp( Y_interp(T0) + (T0/Tref)* (CoolingRate(Tref, rho)/CoolingRate(T0, rho))* (dt_sub/tcool)); /* Townsend cooling update */
    prs = T1*rho*BOLTZMANN/(PROTONMASS*gs.mu);

    /* ------------------------------------------------------
     *       If gas cooled below lower threshold, redefine
     *               pressure so that T = MinGasTemp.
     * ------------------------------------------------------ */
    if (T1 < All.MinGasTemp && T0 > All.MinGasTemp)  prs = All.MinGasTemp*rho*BOLTZMANN/(PROTONMASS*gs.mu);

    u = prs/rho/GAMMA_MINUS1;
  }   /* -- End Subcycle Loop   -- */

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs; /* to internal units */

  return u;
}

/*! \brief Compute gas temperature from internal energy per unit mass.
 *
 *   This function calculates the gas temperature.
 *   Inputs are all in physical cgs units.
 *
 *  \param[in] u   internal energy per unit mass.
 *  \param[in] rho gas density.
 *
 *  \return The gas temperature in Kelvin.
 */
double convert_u_to_temp(double u, double rho)
{
  double temp = GAMMA_MINUS1 * gs.mu * u * (PROTONMASS / BOLTZMANN);
  if (temp<0) terminate ("negative temperature (u,rho): (%12.6e  %12.6e)\n",u,rho);

  return temp;
}

/*! \brief Get cooling rate from gas internal energy.
 *
 *  This function first calculates
 *  (cooling rate-heating rate)/(n_e n_i) in cgs units.
 *  All inputs in physical cgs units.
 *
 *  \param[in] u Gas internal energy per unit mass.
 *  \param[in] rho Gas density.
 *
 *  \return Cooling rate in physical cgs units.
 */
double CoolingRateFromU(double u, double rho)
{
  double temp = convert_u_to_temp(u, rho);
  return CoolingRate(temp, rho);
}

/*! \brief  Calculate (cooling rate-heating rate)/(n_e n_i) in cgs units.
 *
 *  All inputs in physical cgs units.
 *  \param[in] T Gas temperature.
 *  \param[in] rho Gas density.
 *
 *  \return (cooling rate-heating rate)/(n_e n_i).
 */
double CoolingRate(double T, double rho)
{
  double Lambda, Heat;
  Heat = 0.;
  Lambda = lambda_interp (T);
  return (Lambda - Heat);
}

/*! \brief Make cooling rates interpolation table.
 *
 *  Set up interpolation tables in T for cooling rates given in
 *  Cloudy generated cooling table from PLUTO example and corresponding Townsend Y
 *  for integrating cooling function.
 *
 *  \param[in] fname Filename of the Table.
 *
 *  \return void
 */
void MakeRateTable(char *fname)
{
  FILE *fcool;

  gs.XH = 0.7154, gs.YHe = 0.2703, gs.Zmet = 0.0143;
  gs.mu = 1./(2*gs.XH+0.75*gs.YHe+(9./16)*gs.Zmet);
  gs.mue = 2./(1+gs.XH);
  gs.mui = 1./(1/gs.mu-1/gs.mue);

  if(All.MinGasTemp > 0.0)
    Tmin = All.MinGasTemp;
  else
    Tmin = 10.0;
  gs.ethmin =  Tmin * BOLTZMANN / (PROTONMASS * gs.mu * GAMMA_MINUS1);
  /* minimum internal energy for neutral gas */

  mpi_printf (" > Reading table %s from disk...\n",fname);
  fcool = fopen(fname,"r");
  if (fcool == NULL)
      terminate ("cooling table file %s could not be found.\n",fname);

  NCOOLTAB = 0;
  while (fscanf(fcool, "%lf %lf %lf %lf %lf\n", &RateT->Temperature[ NCOOLTAB ],
                               &RateT->LAMBDA[ NCOOLTAB ], &RateT->townY[ NCOOLTAB ],
							   &RateT->townY_inv[ NCOOLTAB ], &RateT->Temperature_invY[ NCOOLTAB ])!= EOF)
	NCOOLTAB++;

}

/*! \brief Initialize the cooling module.
 *
 *  This function initializes the cooling module. In particular,
 *  it allocates the memory for the cooling rate and ionization tables
 *  and initializes them.
 *
 *  \return void
 */
void InitCool(void)
{
  /* allocate and construct rate table */
  mpi_printf("Loading cooling lookup table to memory.\n");
  RateT = (RateTable *)mymalloc("RateT", sizeof(RateTable));
  mpi_printf("Memory reserved %d bytes.\n",sizeof(RateTable));
  MakeRateTable(All.TreecoolFile);

  mpi_printf("GFM_COOLING: time, time begin = %le\t%le\n", All.Time, All.TimeBegin);
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
}

/*! \brief Apply the isochoric cooling to all the active gas cells.
 *
 *  \return void
 */
void cooling_only(void) /* normal cooling routine when star formation is disabled */
{
  int idx, i;

  CPU_Step[CPU_MISC] += measure_time();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

          cool_cell(i);  /* rewrites gas pressure based on cooling rate */
        }
    }
  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

/*! \brief Apply the isochoric cooling to a given gas cell.
 *
 *  This function applies the normal isochoric cooling to a single gas cell.
 *  Once the cooling has been applied according to one of the cooling models
 *  implemented, the internal energy per unit mass, the total energy and the
 *  pressure of the cell are updated.
 *
 *  \param[in] i Index of the gas cell to which cooling is applied.
 *
 *  \return void
 */
void cool_cell(int i)
{
  double dt, dtime;
  double unew, dens, dtcool;

#ifdef AGNWIND_FLAG
  if(SphP[i].AGNFlag == 1)
    return;
#endif /* #ifdef AGNWIND_FLAG */

  dens = SphP[i].Density;

  dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

  dtime = All.cf_atime * dt / All.cf_time_hubble_a;

  dtcool = dtime;

  unew       = DoCooling(dmax(All.MinEgySpec, SphP[i].Utherm), dens * All.cf_a3inv, dtcool);
  SphP[i].Ne = 1;

  if(unew < 0)
    terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

  double du = unew - SphP[i].Utherm;

  if(unew < All.MinEgySpec)
    du = All.MinEgySpec - SphP[i].Utherm;

  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

#ifdef OUTPUT_COOLHEAT
  if(dtime > 0)
    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif /* #ifdef OUTPUT_COOLHEAT */

  set_pressure_of_cell(i);
}

#endif /* #ifdef COOLING */



double lambda_interp(double temperature)
{
  return interp1D(&RateT->Temperature[0], &RateT->LAMBDA[0], temperature, "lambda_interp");
}

double Y_interp(double temperature)
{
  return interp1D(&RateT->Temperature[0], &RateT->townY[0], temperature, "Y_interp");
}

double invY_interp(double townY)
{
  return interp1D(&RateT->townY_inv[0], &RateT->Temperature_invY[0], townY, "invY_interp");
}

/*! \brief A generic one dimensional linear interpolating function.
 *
 *  This function applies the Table lookup by binary search
 *  used for interpolating y as a function of x.
 *  x data in table is assumed to be arranged in ascending order
 *
 *  \param[in] x_data Independent variable of interpolation.
 *  \param[in] y_data Dependent variable of interpolation.
 *  \param[in] x_request Independent variable value where interpolation is needed.
 *  \param[in] msg Some random string message to identify from where function is called (useful for debugging).
 *
 *  \return Interpolated value.
 *
*/
double interp1D(double* x_data, double* y_data, double x_request, char* msg)
{
  int ntab = NCOOLTAB; /*(int)(sizeof(x_data)/sizeof(x_data[0]));  number of entries */
  int klo = 0;
  int khi = ntab - 1;
  int kmid;
  double xmid, dx, scrh;

  if (x_request > x_data[khi] || x_request < x_data[klo])
  {
    mpi_printf ("called from %s\n",msg);
    terminate ("requested value out of range: %12.6e\n",x_request);
  }
  while (klo != (khi - 1))
  {
    kmid = (klo + khi)/2;
    xmid = x_data[kmid];
    if (x_request <= xmid)
		khi = kmid;
    else if (x_request > xmid)
		klo = kmid;
  }
  /*  Compute r.h.s */
  dx       = x_data[khi] - x_data[klo];
  scrh     = y_data[klo]*(x_data[khi] - x_request)/dx + y_data[khi]*(x_request - x_data[klo])/dx;

  return scrh;
}

/*! \brief Computes the cooling time of gas.
 *
 *  The function calculates the cooling time of the gas in internal units
 *  from the provided cooling table. Arguments passed in internal units.
 *
 *  \param[in] u the internal energy per unit mass of the gas cell.
 *  \param[in] rho   the proper density of the gas cell.
 *
 *  \return The cooling time of the gas cell in internal units.
 */
double GetCoolingTime (double u, double rho)
{
  double necgs, nicgs, ratefact, temperature, tcool;

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* convert to physical cgs units */
  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs; /* density in cgs units */

  necgs =  rho / (gs.mue*PROTONMASS); /* electron number dens in cgs units */
  nicgs =  rho / (gs.mui*PROTONMASS); /* ion number dens in cgs units */
  ratefact = necgs * nicgs / rho;

  temperature = convert_u_to_temp(u, rho);
  tcool = u/(ratefact*CoolingRate(temperature, rho));
  tcool *= All.HubbleParam / All.UnitTime_in_s; /* cooling time in internal units */
  return tcool;
}
