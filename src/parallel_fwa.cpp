/*
        fwaVina                       Date: 26/11/2018

        This file is revised from parallel_fwa.cpp in AutoDock Vina.

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

#include "parallel.h"
#include "parallel_fwa.h"
#include "coords.h"
#include "parallel_progress.h"

struct parallel_fwa_task {
	model m;
	output_container out;
	rng generator;
	parallel_fwa_task(const model& m_, int seed) : m(m_), generator(static_cast<rng::result_type>(seed)) {}
};

typedef boost::ptr_vector<parallel_fwa_task> parallel_fwa_task_container;

struct parallel_fwa_aux {
	const fwa_search* fwa;
	const precalculate* p;
	const igrid* ig;
	const precalculate* p_widened;
	const igrid* ig_widened;
	const vec* corner1;
	const vec* corner2;
	const int num_fwa;
	const double w,c1,c2;
	parallel_progress* pg;
	parallel_fwa_aux(const fwa_search* fwa_, const precalculate* p_, const igrid* ig_, const precalculate* p_widened_, const igrid* ig_widened_, const vec* corner1_, const vec* corner2_, parallel_progress* pg_, int num_fwa_, double w_, double c1_, double c2_)
		: fwa(fwa_), p(p_), ig(ig_), p_widened(p_widened_), ig_widened(ig_widened_), corner1(corner1_), corner2(corner2_), pg(pg_), num_fwa(num_fwa_), w(w_), c1(c1_), c2(c2_) {}
	void operator()(parallel_fwa_task& t) const {
		(*fwa)(t.m, t.out, *p, *ig, *p_widened, *ig_widened, *corner1, *corner2, pg, t.generator, num_fwa, w, c1, c2);
	}
};

void merge_output_containers(const output_container& in, output_container& out, fl min_rmsd, sz max_size) {
	VINA_FOR_IN(i, in)
		add_to_output_container(out, in[i], min_rmsd, max_size);
}

void merge_output_containers(const parallel_fwa_task_container& many, output_container& out, fl min_rmsd, sz max_size) {
	min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
	VINA_FOR_IN(i, many)
		merge_output_containers(many[i].out, out, min_rmsd, max_size);
	out.sort();
}

void parallel_fwa::operator()(const model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator, int num_fwa, double w, double c1, double c2) const {
	parallel_progress pp;
	parallel_fwa_aux parallel_fwa_aux_instance(&fwa, &p, &ig, &p_widened, &ig_widened, &corner1, &corner2, (display_progress ? (&pp) : NULL), num_fwa, w, c1, c2);
	parallel_fwa_task_container task_container;
	VINA_FOR(i, num_tasks)
		task_container.push_back(new parallel_fwa_task(m, random_int(0, 1000000, generator)));
	if(display_progress) 
		pp.init(num_tasks * fwa.num_steps);
	parallel_iter<parallel_fwa_aux, parallel_fwa_task_container, parallel_fwa_task, true> parallel_iter_instance(&parallel_fwa_aux_instance, num_threads);
	parallel_iter_instance.run(task_container);
	merge_output_containers(task_container, out, fwa.min_rmsd, fwa.num_saved_mins);
}
