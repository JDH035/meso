
# include "OPFZoneTempFD.hpp"
# include <iostream>
# include <fstream>
# include <cmath>

// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

OPFZoneTempFD::OPFZoneTempFD(const CommonParams& pin,
							 const GetPot& input_params) : p(pin), c(p), tempL(p)
{

    // ---------------------------------------
    // set needed parameters:
    // ---------------------------------------

    nxyz = p.nx*p.ny*p.nz;
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    deli = (nz+2)*(ny+2);
	delj = (nz+2);
	delk = 1;
    co = input_params("PFApp/C0",0.5);
	initNoise = input_params("PFApp/initNoise",0.1);
    M = input_params("PFApp/Mob",1.0);
    noiseStr = input_params("PFApp/noiseStr",0.1);
    wzone = input_params("PFApp/wzone",10.0);
	vzone = input_params("PFApp/vzone",1.0);
	Tmin = input_params("PFApp/Tmin",373);
	Tmax = input_params("PFApp/Tmax",473);
	N = input_params("PFApp/Poly",288);
	alpha = input_params("PFApp/alpha",3.9);
	beta = input_params("PFApp/beta",0.028);
	templating = input_params("PFApp/templating",0);
	tempSpace1 = input_params("PFApp/tempSpace1",4);
	tempSpace2 = input_params("PFApp/tempSpace2",2);
	tempSpace3 = input_params("PFApp/tempSpace3",18);
	tempSpace4 = input_params("PFApp/tempSpace4",24);
	tempSpace5 = input_params("PFApp/tempSpace5",3);
	bX = input_params("PFApp/bX",0);
	bY = input_params("PFApp/bY",1);
	bZ = input_params("PFApp/bZ",0);
	topWetting = input_params("PFApp/topWetting",0);
	bcp = input_params("PFApp/bcp",1);
	chiNs =  70123.549*pow(co,8)
			-280494.19*pow(co,7)
			+485599.11*pow(co,6)
			-475067.64*pow(co,5)
			+287860.20*pow(co,4)
			-111184.23*pow(co,3)
			+27041.564*pow(co,2)
			-3878.3585*co
			+271.87882;
}

// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

OPFZoneTempFD::~OPFZoneTempFD()
{

}

// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void OPFZoneTempFD::initPhaseField()
{

    //	---------------------------------------
    // initialize the concentration field:
    //	---------------------------------------

    srand(time(NULL)*(p.rank+1));   // set the random seed
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double r = (double)rand()/RAND_MAX;
                double val = co + initNoise*(r-0.5);
                c.setValue(ndx,val);
            }//i
        }//j
    }//k
	
	//ensure all processors reach this point before proceeding
	MPI::COMM_WORLD.Barrier();
	
	if (templating != 0)
	{
		//determine z position for (1)2d or (0)3d simulation
		int k = (nz==1);

		//apply templating type
		switch(templating)
		{
			case 1: // vertical striped templating 
					// int 		diameter of templating = tempSpace1
					// double 	x-separation = tempSpace3
			{
				double L = tempSpace3;
				double w = tempSpace1*0.5;
				double b = 2*M_PI/L;
				double shift = cos(b*w);
				double height = (1-shift);
				for (int i=1; i<nx+1; i++) {
					double zDot = cos(b*(i+p.xOff-1));
					zDot -= shift;
					zDot /= height;
					if (zDot > 0.0) {
						for (int j=1; j<ny+1; j++) {					
							int ndx = i*deli + j*delj + k*delk;
							tempL.setValue(ndx,zDot);
						}//j
					}//if
				}//i				
			}//case 1
			break;
			
			case 2: /// leading edge horizontal line array
					// int 		line diameter = tempSpace1
					// int		x-length = tempSpace2
					// double 	y-separation = tempSpace3
			{
				double L = tempSpace3;
				int endAt = min(tempSpace2-p.xOff,nx);
				if (endAt > 0) {
					double w = tempSpace1*0.5;
					double b = 2*M_PI/L;
					double shift = cos(b*w);
					double height = (1-shift);
					for (int j=1; j<ny+1; j++) {
						double zCyl = cos(b*(j-1));
						zCyl -= shift;
						zCyl /= height;
						if (zCyl > 0.0) {
							for (int i=1; i<endAt+1; i++) {
								int ndx = i*deli + j*delj + k*delk;
								tempL.setValue(ndx,zCyl);
							}//i
						}//if
					}//j					
				}//if
			}//case 2
			break;
			
			case 3: // rectangle dot array
					// int 		dot diameter = tempSpace1
					// int 		number of columns = tempSpace2
					// double 	x-separation distance = tempSpace3
					// double	y-separation distance = tempSpace4
					// int 		additional columns at regular spacing = tempSpace5
			{
				int changeAt = min(floor(tempSpace3*(tempSpace2-0.5))-p.xOff,nx*1.0);
				double w = tempSpace1*0.5;
				double periodX = 2*M_PI/tempSpace3;
				double periodY = 2*M_PI/tempSpace4;
				double shiftX = cos(periodX*w);
				double shiftY = cos(periodY*w);
				double heightX = (1-shiftX);
				double heightY = (1-shiftY);
				
				for (int j=1; j<ny+1; j++) {
					double yDot = (cos(periodY *(j-1)) - shiftY) / heightY;
					for (int i=1; i<nx+1; i++) {
						double zDot = 0;
						if (i<=changeAt) {
							zDot = (cos(periodX *(i+p.xOff-1)) - shiftX) / heightX + yDot;
							if (zDot > 0.0) {
								int ndx = i*deli + j*delj + k*delk;
								tempL.setValue(ndx,0.75);
							}//if
						} else {
							if (fmod((i-(tempSpace3*(tempSpace2-0.5))+p.xOff),(tempSpace5*tempSpace3))>=((tempSpace5-1)*tempSpace3)) {
								zDot = (cos(periodX *(i+p.xOff-1)) - shiftX) / heightX + yDot;
								if (zDot > 0.0) {
									int ndx = i*deli + j*delj + k*delk;
									tempL.setValue(ndx,0.75);
								}//if
							}//if
						}//else
							
							
							
						
					}//i
				}//j					
				
				
			}//case 3
			break;
			
			case 4: // square dot array 45
					// int 		dot diameter = tempSpace1
					// int 		number of columns = tempSpace2
					// double 	separation distance = tempSpace3
			{
				double L = tempSpace3;
				int endAt = min(floor(L*0.5*(0.25+(tempSpace2-1)))-p.xOff,nx*1.0);
				double w = tempSpace1*0.5;
				double b = 2*M_PI/L;
				double shift = cos(b*w)*cos(b);
				double height = (1-shift);
				if (endAt > 0) {
					for (int j=1; j<ny+1; j++) {
						for (int i=1; i<endAt+1; i++) {
							double zDot = cos(b*(i+p.xOff-1))*cos(b*(j-1));
							zDot -= shift;
							zDot /= height;
							if (zDot > 0.0) {
								int ndx = i*deli + j*delj + k*delk;
								tempL.setValue(ndx,zDot);
							}//if
						}//i
					}//j					
				}//if
			}//case 4
			break;
			
			case 5: // square dot array plus single row
					// int 		dot diameter = tempSpace1
					// int 		number of columns = tempSpace2
					// double 	separation distance = tempSpace3
			{
				double L = tempSpace3;
				int endAt = min(floor(L*(0.5+(tempSpace2-1)))-p.xOff,nx*1.0);
				double w = tempSpace1*0.5;
				double b = 2*M_PI/L;
				double shift = cos(b*w)+cos(b);
				double height = (2-shift);
				//loop through to create full column template
				if (endAt > 0) {
					for (int j=1; j<ny+1; j++) {
						for (int i=1; i<endAt+1; i++) {
							double zDot = cos(b*(i+p.xOff-1))+cos(b*(j-1));
							zDot -= shift;
							zDot /= height;
							if (zDot > 0.0) {
								int ndx = i*deli + j*delj + k*delk;
								tempL.setValue(ndx,zDot);
							}//if
						}//i
					}//j					
				}//if
				//loop through single row template
				for (int j=floor(0.5*tempSpace3); j<ceil(1.5*tempSpace3)+1; j++) {
					for (int i=1; i<nx+1; i++) {
						double zDot = cos(b*(i+p.xOff-1))+cos(b*(j-1));
						zDot -= shift;
						zDot /= height;
						if (zDot > 0.0) {
							int ndx = i*deli + j*delj + k*delk;
							tempL.setValue(ndx,zDot);
						}//if
					}//i
				}//j
			}//case 5
			break;
			
			case 6: // square dot array 45 plus single row
					// int 		dot diameter = tempSpace1
					// int 		number of columns = tempSpace2
					// double 	separation distance = tempSpace3
			{
				double L = tempSpace3;
				int endAt = min(floor(L*0.5*(0.25+(tempSpace2-1)))-p.xOff,nx*1.0);
				double w = tempSpace1*0.5;
				double b = 2*M_PI/L;
				double shift = cos(b*w)*cos(b);
				double height = (1-shift);
				if (endAt > 0) {
					for (int j=1; j<ny+1; j++) {
						for (int i=1; i<endAt+1; i++) {
							double zDot = cos(b*(i+p.xOff-1))*cos(b*(j-1));
							zDot -= shift;
							zDot /= height;
							if (zDot > 0.0) {
								int ndx = i*deli + j*delj + k*delk;
								tempL.setValue(ndx,zDot);
							}//if
						}//i
					}//j					
				}//if
				for (int j=floor(tempSpace3*0.75); j<ceil(tempSpace3*1.75)+1; j++) {
					for (int i=1; i<nx+1; i++) {
						double zDot = cos(b*(i+p.xOff-1))*cos(b*(j-1));
						zDot -= shift;
						zDot /= height;
						if (zDot > 0.0) {
							int ndx = i*deli + j*delj + k*delk;
							tempL.setValue(ndx,zDot);
						}//if
					}//i
				}//j						
			}//case 6
			break;
			
		}//switch templating
	}//if templating
}

// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------
double calculateC2(double x, double x2, double g2, double g4)
{
	double c_2 = -(  5.920 
					-14.31	  		*g2
					-398.5	  		*g4		//+398.5
					+2.025	  *x
					-4.285	  *x	*g2
					+39.47	  *x	*g4
					+0.005522 *x2
					-0.1915	  *x2	*g2
					-1.003	  *x2	*g4 );
	return c_2;
}
	
double calculateC3(double x, double x2, double g, double g3, double g5)
{
	double c_3 = -(  9.741		*g
					-39.46		*g3
					-999.0		*g5
					+9.224	*x	*g
					+14.69	*x	*g3
					-510.3	*x	*g5
					+0.06433*x2	*g
					-1.281	*x2	*g3
					+15.70	*x2	*g5 );
	return c_3;
}
double calculateC4(double x, double x2, double g2, double g4)
{
	double c_4 =  (  9.686 
					+53.00		*g2
					-1775.0		*g4
					+3.600	*x
				  //+0.00	*x	*g2
				  //+0.00	*x	*g4
					+0.02068*x2
					-0.2385	*x2	*g2
					-0.4559	*x2	*g4 );
	return c_4;
}
double calculateC5(double x, double g2, double g4)
{
	double c_5 =(( (0.7885	+0.2119	*x)				//0.7853 //0.2356
				  +(-5.654	-1.170	*x)	*g2
				  +(-16.22	+3.659	*x)	*g4 ) 
				/(1	+(	 0.139		*x				//0.1185
						-0.7423		*x	*g2
						+5.481		*x	*g4 )) );
	return c_5;
}

double calculateC6(double x, double x2, double g2, double g4, double c_2)
{
	double c_6 = -c_2/(  0.5 
					  //+0.0			*g2
					  //+0.0			*g4
						+0.2230		*x			//0.2357
						+-1.956		*x	*g2
						+7.147		*x	*g4
						+0.0006666	*x2
						+-0.02858	*x2	*g2
						+0.05316	*x2	*g4 );					
	return c_6;
}



void OPFZoneTempFD::updatePhaseField()
{
	
	//allocate storage space
    SfieldFD mu(p);
    SfieldFD mob(p);
	SfieldFD C_5(p);
	SfieldFD C_6(p);
	
    // ---------------------------------------
    // calculate chemical potential & mobility
    // ---------------------------------------

    c.updateBoundaries(bX,bY,bZ);

	//calculate local coefficient values
	double Mc = M;						//assume constant mobility
	double Temp = Tmax;					//assume constant temperature
	double chiN  = (beta + alpha/Temp)*N;
	std::cout << chiN <<" "<< beta <<" "<< alpha  <<" "<< Temp <<" "<< N << std::endl;
	double x = chiN - chiNs;
	double g = abs(0.5 - co);
	double g2 = pow(g,2);
	double g3 = pow(g,3);
	double g4 = pow(g,4);
	double g5 = pow(g,5);
	double x2 = pow(x,2);
	
	double c_2 = calculateC2(x,x2,g2,g4);
	double c_3 = calculateC3(x,x2,g,g3,g5);
	double c_4 = calculateC4(x,x2,g2,g4);
	double c_5 = calculateC5(x,g2,g4);
	double c_6 = calculateC6(x,x2,g2,g4,c_2);
		
	
    for (int i=1; i<nx+1; i++) {

		//if zone annealing is active
		if (vzone != 0) { 
			//calculate time 
			double time = current_step*p.dt;
		
			//calculate offset
			int xf = int(time*vzone - 0.5*wzone - p.xOff);
		
			//calculate local mobility
			Mc = 0.5*(1 - tanh(6.0*double(i-xf)/wzone));
			
			//calculate local temperature
			Temp = Tmin + Mc*(Tmax - Tmin);

			//update local coefficient values
			chiN  = (beta + alpha/Temp)*N;
			x = chiN - chiNs;
			g = abs(0.5 - co);
			g2 = pow(g,2);
			g3 = pow(g,3);
			g4 = pow(g,4);
			g5 = pow(g,5);
			x2 = pow(x,2);
		
			c_2 = calculateC2(x,x2,g2,g4);
			c_3 = calculateC3(x,x2,g,g3,g5);
			c_4 = calculateC4(x,x2,g2,g4);
			c_5 = calculateC5(x,g2,g4);
			c_6 = calculateC6(x,x2,g2,g4,c_2);
		}
		

		for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                
				//get current position
				int ndx = i*deli + j*delj + k*delk;
				
				//get current concentration 
				double cc = c.getValue(ndx);
				
				//calculate chemical potential
                double dFdc = 2*c_2*(cc-0.5) 
							+ 3*c_3*(cc-0.5)*(cc-0.5) 
							+ 4*c_4*(cc-0.5)*(cc-0.5)*(cc-0.5); 
				
				//write mobility values
				mob.setValue(ndx, Mc);
				
				//write potential values
				mu.setValue(ndx,dFdc);
				
				//write coefficients needed outside current scope
				C_5.setValue(ndx,c_5);
				C_6.setValue(ndx,c_6);
            }//k
        }//j
    }//i
	
	//ensure all processors reach this point before proceeding
	MPI::COMM_WORLD.Barrier();
	
	//update boundary conditions for mu and mobility
	mu.updateBoundaries(bX,bY,bZ);
    mob.updateBoundaries(bX,bY,bZ);
	
	// ---------------------------------------
    // Apply surface templating:
    // ---------------------------------------

	if (templating)
	{
		//determine z position for (1)2d or (0)3d simulation
		int k = (nz==1);
		// apply templating
		for (int j=1; j<ny+1; j++) {
			for (int i=1; i<nx+1; i++) {
				int ndx = i*deli + j*delj + k*delk;
				double zDot = tempL.getValue(ndx);
				double cc = c.getValue(ndx);
				c.setValue(ndx,max(cc,zDot));
			}//i
		}//j
	}//if templating
	
	if (topWetting)
	{
		//determine z position for (1)2d or (nz+1)3d simulation
		int k = nz*(nz>1)+1;
		// apply top surface
		for (int j=1; j<ny+1; j++) {
			for (int i=1; i<nx+1; i++) {
				int ndx = i*deli + j*delj + k*delk;
				c.setValue(ndx,0.0);
			}//i
		}//j
	}//if topWetting
	
		
    // ---------------------------------------
    // update CH equation:
    // ---------------------------------------

	//calculate first laplacian of CH right hand side
	SfieldFD inside = mu - C_5*c.Laplacian();
	inside.updateBoundaries(bX,bY,bZ);
	
	//calculate total right hand side of CH equation
	SfieldFD RHS = inside.Laplacian(mob); 	
	// Apply energy barrier for Di-block system:
	RHS -= C_6*(c-co);
	
	//add thermal noise
	for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
				
				// Add random fluctuations:
				double r = (double)rand()/RAND_MAX;
                double val = noiseStr*(r-0.5);
                RHS.addValue(ndx,val);
            }
        }
    }
	
	RHS.updateBoundaries(bX,bY,bZ);
	c += p.dt*RHS;
	std::cout << c_2 << " " << c_3 << " " << c_4 << " " << c_5 << " " << c_6 <<  std::endl;
}//updatePhaseField

// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void OPFZoneTempFD::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
}
