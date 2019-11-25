/*
*  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/


#ifndef GalacticBinaryData_h
#define GalacticBinaryData_h

#include <stdio.h>

void GalacticBinaryReadData(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags);
void GalacticBinarySimulateData(struct Data *data);
void GalacticBinaryInjectVerificationSource(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags);
void GalacticBinaryInjectSimulatedSource(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags);
void GalacticBinaryCatalogSNR(struct Data *data, struct Orbit *orbit, struct Flags *flags);

#endif /* GalacticBinaryData_h */
