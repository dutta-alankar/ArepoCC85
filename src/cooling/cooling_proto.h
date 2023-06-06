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
 * \file        src/cooling/cooling_proto.h
 * \date        05/2018
 * \brief       Header for cooling functions.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 27.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef INLINE_FUNC
#define INLINE_FUNC
#endif /* #ifndef INLINE_FUNC */

double convert_u_to_temp(double u, double rho);
double CoolingRate(double T, double rho);
double CoolingRateFromU(double u, double rho);
double DoCooling(double u_old, double rho, double dt);
double GetCoolingTime(double u, double rho);

double interp1D(double* x_data, double* y_data, double x_request, char* msg);
double lambda_interp(double temperature);
double Y_interp(double temperature);
double invY_interp(double townY);

void InitCool(void);
