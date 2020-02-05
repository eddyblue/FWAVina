/*
        FWAVina version 1.0                     Date: 26/11/2018

        This file is revised from mutate.h in AutoDock Vina.

        Authors: Marcus C. K. Ng  <marcus.ckng@gmail.com>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_fwa_MUTATE_H
#define VINA_fwa_MUTATE_H

#include "fwa.h"
#include "model.h"
#include "quasi_newton.h"

// does not set model
void fwa_mutate_conf(conf& c, const model& m, fl amplitude, rng& generator);
void fwa_mutate_conf_position(conf& c, const model& m, fl amplitude, rng& generator);
void fwa_mutate_conf_oritation(conf& c, const model& m, fl amplitude, rng& generator);
void fwa_mutate_conf_torsion(conf& c, const model& m, fl amplitude, rng& generator);




#endif
