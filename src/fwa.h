/*
        fwaVina version 1.0                     Date: 26/11/2014

        Authors: Marcus C. K. Ng  <marcus.ckng@gmail.com>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
#ifndef fwa_H_
#define fwa_H_

#include "common.h"
#include "conf.h"
#include <vector>

class fwa
{	
public:
	struct bird{
		     vec velocity,vO;  //velocity of position, velocity of orientation
		     
		     vec pbest_pos;
		     vec current_position;
 
		     qt pbest_orientation; 
		     qt current_orientation;
		     
		     fl* pbest_torsion;
		     fl* current_torsion;
		     fl* vT;    //velocity of torsion
		     
		     double pbest_fit,tmp_fit;  

	};
	
	int torsionSize;
	double w,c1,c2;						// weight, learning coefficient(1&2)
	rng g;

	int number;							//number of fwa
	vec corner1,corner2;				//corners of search box
	
	static vec gbest_position;			//global best of degree of freedom vector
	static qt gbest_orientation;		//global best of degree of freedom vector
	static fl* gbest_torsion;			//global best of degree of freedom vector
	
	static double gbest_fit;			//global best value

	std::vector<bird> particle;
	
	
	  double R1Max_;
	  double R1Min_;
	  double R2Max_;
	  double R2Min_;
	
	
	fwa(int,double,double,double,const vec,const vec,rng&,conf&);
	void init(rng&,conf&);
    
    /*Update Velocity for vector*/
    void updateVelocity(rng&,int);
    void updateVelocityO(rng&,int);
    void updateVelocityT(rng&,int,sz);
    
    
    /*Compute the new pose*/
    void computeNewPositions(int);
	void computeNewOrientation(int);
	void computeNewTorsion(int,rng&,sz);
    
    
    /*Update personal best value*/
	void updatePersonalBest(int,double);
    
    /*Update global best value*/
	void updateGlobalBest(int);
    
    /*Get personal best value*/
	double getPersonalBest(int);	//return the personal best value
    
    /*Set personal best vector*/
    void updateBestPosition(int,vec);
    void updateBestOrientation(int,qt);
    void updateBestTorsion(int,fl,sz);


    /*Get current vector*/
	vec getCurrentPosition(int);
    qt getCurrentOrientation(int);
    fl getCurrentTorsion(int,sz);
    
};


#endif /*fwa_H_*/
