/*
        fwaVina version 1.0                     Date: 26/11/2014

        This file is revised from monte_carlo.cpp in AutoDock Vina.

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


#include "fwa_search.h"
#include "conf.h"
#include "coords.h"
#include "quasi_newton.h"
#include "fwa_mutate.h"
#include <algorithm>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<time.h>
#include <vector>
#include<string>
using namespace std;



output_type fwa_search::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,int num_fwa, double w, double c1, double c2) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator, num_fwa, w, c1, c2); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}


 void swap(int *a,int *b)
 {
    int temp;
    temp=*a;
    *a=*b;
    *b=temp;
 }

 int rand_X(int x)
{
    return rand()%x;
}

 void GetRandomArray(int arry[],int number)
{
    int i,k;
    srand((int)time(0));
    for(i=0;i<number-1;i++)
    {
        k=rand()%(number-i)+i;
        swap(&arry[i],&arry[k]); 
    }
}

// out is sorted
void fwa_search::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,int num_fwa,double w,double c1,double c2) const {
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	
//《改》定义烟花数组和火花数组；
std::vector<output_type> Fwa;
std::vector<sort_type> Esort;
sort_type sort_a;
const int LENGTH=54;
int index[LENGTH];
int b;
	conf_size s = m.get_size();
	change g(s);
	output_type tmp(s, 0);

//《改》初始化n个烟花
int y;
for(y=0;y<num_fwa;y++)
	{
     tmp.c.randomize(corner1, corner2, generator); 
	 Fwa.push_back(tmp);
	}

	fl best_e = max_fl;
	quasi_newton quasi_newton_par; 
		
	double energy=0;
	int count=0;


VINA_U_FOR(step, num_steps/80) 
{

		if(increment_me)
			++(*increment_me);
		
//《改》每个烟花先做一次拟牛顿	
	for(y=0;y<num_fwa;y++)
		{
		tmp=Fwa[y];
	    output_type	candidate=tmp;
		quasi_newton_par.max_steps =0;
		quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
		if(step == 0 || metropolis_accept(tmp.e, candidate.e, temperature, generator))
		Fwa[y] = candidate;
        }

 //《改》将每个烟花按e从小到大排序

	sort(Fwa.begin(),Fwa.end());
		//for each particle loop

//烟花开始爆炸
	std::vector<output_type> Spark;
	int j;
	 for(y=0;y<num_fwa;y++)
	{
		    if(y<=2)
			for(j=0;j<num_fwa-y;j++)
			{
			 output_type candidate=Fwa[y];
			 fwa_mutate_conf_torsion(candidate.c, m, mutation_amplitude, generator); 
    		 quasi_newton_par.max_steps=ssd_par.evals;
		//	 quasi_newton_par.max_steps=0;
			 quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
		//	 m.set(candidate.c);
			 Spark.push_back(candidate);
			}
			if(y>2&&y<6)
			for(j=0;j<num_fwa-y;j++)
			{
			 output_type candidate=Fwa[y];
			 if(j<(num_fwa-y)/2)
			 fwa_mutate_conf_oritation(candidate.c, m, mutation_amplitude, generator); 
			 else
			 fwa_mutate_conf_torsion(candidate.c, m, mutation_amplitude, generator); 
    		 quasi_newton_par.max_steps=ssd_par.evals;
		//	 quasi_newton_par.max_steps=0;
			 quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
		//	 m.set(candidate.c);
			 Spark.push_back(candidate);
			}
			if(y>=6)
			for(j=0;j<num_fwa-y;j++)
			{
			 output_type candidate=Fwa[y];
			 fwa_mutate_conf_oritation(candidate.c, m, mutation_amplitude, generator); 
    		 quasi_newton_par.max_steps=ssd_par.evals;
		//	 quasi_newton_par.max_steps=0;
			 quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
		//	 m.set(candidate.c);
			 Spark.push_back(candidate);
			}
		
	}

	for(y=0;y<num_fwa;y++)
	{
	 Spark.push_back(Fwa[y]);
	}

//将火花入栈
	for(b=0;b<LENGTH;b++)
	{
	index[b]=b;
	sort_a.index=b;
	sort_a.e=Spark[b].e;
	Esort.push_back(sort_a);
	}

	GetRandomArray(index,LENGTH);

	for(b=0;b<num_fwa;b++)//变异烟花，随机选9个出来变异
	{ 
	if(index[b]!=0)
	{
	output_type candidate=Spark[index[b]];
	fwa_mutate_conf(candidate.c, m, mutation_amplitude, generator); 
	quasi_newton_par.max_steps=0;
    quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
 //   m.set(candidate.c);
    Spark[index[b]]=candidate;
	}
	}



//将所有的火花进行排序，选出优质烟花；
//	sort(Spark.begin(),Spark.end());
	sort(Esort.begin(),Esort.end());


//	for(y=num_fwa-1;y>=0;y--)
	for(y=0;y<num_fwa;y++)
	{
//	Fwa[y]=Spark[y];
	Fwa[y]=Spark[Esort[y].index];
	
	tmp=Fwa[y];
		//tmp = candidate;
//		m.set(tmp.c); // FIXME? useless?

			// FIXME only for very promising ones
		if(tmp.e < best_e || out.size() < num_saved_mins) 
		{  
			quasi_newton_par.max_steps = ssd_par.evals;
//			quasi_newton_par.max_steps =0;
			quasi_newton_par(m, p, ig, tmp, g, authentic_v);
			m.set(tmp.c); // FIXME? useless?
			tmp.coords = m.get_heavy_atom_movable_coords();
			add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
			if(tmp.e < best_e)
			{
					best_e = tmp.e;
					Fwa[y]=tmp;     //如果优质烟花做拟牛顿后值更小，则替换优质烟花。
			}
		}

	

		/***Criteria defined by fwaVina***/

		if(std::abs(best_e - energy) < 0.0001)
		{

			count +=1;
			if(count > 450)   //测试了100，结果又不稳定了《改》
			{

//输出迭代次数	
FILE *fplj;
int iter=step; 
fplj=fopen("F:\\research\\NO3\\docking\\FWAvinaOut\\interation.txt","a");
fprintf(fplj,"real %d\n",iter);
fclose(fplj);
//结束
				//printf("Terminated: %d \n",step);
				step = num_steps; //break the loop
				count = 0;


			}
			
		}else{
			energy = best_e;
			count =0;
		}

	}
	
	
}


	
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order

}
