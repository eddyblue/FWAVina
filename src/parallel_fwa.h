/*
        FWAVina                       Date: 26/11/2018

        This file is revised from parallel_mc.h in AutoDock Vina.

        Authors: Marcus, C. K. Ng  <marcus.ckng@gmail.com>

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

#ifndef VINA_PARALLEL_fwa_H
#define VINA_PARALLEL_fwa_H

#include "fwa_search.h"

struct parallel_fwa {
	fwa_search fwa;
	sz num_tasks;
	sz num_threads;
	bool display_progress;
	parallel_fwa() : num_tasks(8), num_threads(1), display_progress(true) {}
	void operator()(const model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generatorint, int num_fwa, double w, double c1, double c2) const;
};

#endif
