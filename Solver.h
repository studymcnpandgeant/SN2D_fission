#ifndef SOLVER_H_
#define SOLVER_H_


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <memory.h>


#include "Quadrature.h"


class Solver {

private:
    /* Geometry */
	int _Nx;
	int _Ny;
	int _Lx;
	int _Ly;
	int _deltax;
	int _deltay;
	int _idx1;
	int _idx2;
	int _idy1;
	int _idy2;
	double _S0;
	double pi;

	/* Material */	
	double _sigma_t;
	double _sigma_s0;
	double _sigma_s1;
	double _sigma_s2;
	double _externalsource;
	double _nu_sigma_f;
	
	/* Quadrature */
	int _n; /* Number of directions */
	int _M4;
	int _M1;
	int _M2;
	int _M3;
	Quadrature* _quadrature;

	/* Scalar flux */
	double*** _flux_x;
	double*** _flux_y;
	double*** _flux_c;//three dimension arry
	double** _flux;
	double** _flux_pre;
	
	/* Angular flux */
	//double* _angular_flux;
	
	/* Boundary flux */
	//double* _boundary_flux;

	/* Source */
	double*** _tot_source;
	double*** _fiss_source;
	//double*** _ex_source;
	//double* _fiss_source;

	double _keff;
	
	/* Iteration */
	int _max_iter;
	double _converge_thresh;
	//     double _res_keff;
	double _res_flux;
	double _res_keff;

public:

	Solver(int Nx=100,int Ny=100, int n=4);
	virtual ~Solver();

	void setQuadrature(Quadrature* quad);

	void startIteration();//key project

	void InitializeFlux();
	void InitializeSource();
	void InputAnVoidFlux();

	//void computeTotalSource();
	//void computeFissionSource();

	void Quadrat_1();
	void Quadrat_2();
	void Quadrat_3();
	void Quadrat_4();

	void storeFlux();
	//void zeroFlux();

	void sourceIteration();
	void externalsource();//pull-in a fixed source

	void computeResidual();
	void computeTotalSource();
	void computeKeff();
	void computeFissionSource();

	double getResidual() { return _res_flux; };

};

Solver::Solver(int Nx,int Ny, int n) {

	_Nx = Nx;
	_Ny = Ny;
	_n = n;
	_M4 = n*(n+2)/2;
	_M1 = n*(n+2)/8;
	_M2 = n*(n+2)/4;
	_M3 = 3*n*(n+2)/8;
	pi = 3.14159;

	_Lx = 100;
	_Ly = 100;
	_deltax = _Lx/_Nx;
	_deltay = _Ly/_Ny;
	_sigma_t = 0.25;
	_sigma_s0 = 0.15;
	_sigma_s1 = 0.01;
	_sigma_s2 = 0.0025;
	_nu_sigma_f = 0.0225;
	_keff =1.0;
	_max_iter = 1000;
	_S0=0.0;//no ex source

	_converge_thresh = 0.00001;
	_res_keff = 1.0;
	//_res_flux = 1.0;
	_quadrature = NULL;
	_idx1=1;
	_idx2=25/_deltax;
	_idy1=25/_deltay+1;
	_idy2 = _idy1 + 25 / _deltay - 1;
	
}

Solver::~Solver() {
	/*if (_quadrature != NULL)
		delete _quadrature;

	for (int i=0;i<_Ny+1;i++){
		for (int j=0;j<_Nx;j++){
			if (_flux_x[i][j] != NULL)
			{
				delete [] _flux_x[i][j];
			}
		}
	}

	for (int i=0;i<_Ny+1;i++){
		if (_flux_x[i] != NULL)
			{
				delete [] _flux_x[i];
			}
	}

	if (_flux_x != NULL)
			{
				delete [] _flux_x;
			}


	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx+1;j++){
			if (_flux_y[i][j] != NULL)
			{
				delete [] _flux_y[i][j];
			}
		}
	}

	for (int i=0;i<_Ny;i++){
		if (_flux_y[i] != NULL)
			{
				delete [] _flux_y[i];
			}
	}

	if (_flux_y != NULL)
			{
				delete [] _flux_y;
			}

    for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			if (_flux_c[i][j] != NULL)
			{
				delete [] _flux_c[i][j];
			}
		}
	}

	for (int i=0;i<_Ny;i++){
		if (_flux_c[i] != NULL)
			{
				delete [] _flux_c[i];
			}
	}

	if (_flux_c != NULL)
			{
				delete [] _flux_c;
			}

	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			if (_tot_source[i][j] != NULL)
			{
				delete [] _tot_source[i][j];
			}
		}
	}

	for (int i=0;i<_Ny;i++){
		if (_tot_source != NULL)
			{
				delete [] _tot_source;
			}
	}

	if (_tot_source != NULL)
			{
				delete [] _tot_source;
			}
    

	for (int i=0;i<_Ny;i++){
		if (_flux[i] != NULL)
			{
				delete [] _flux[i];
			}
	}

	if (_flux != NULL)
			{
				delete [] _flux;
			}

    for (int i=0;i<_Ny;i++){
		if (_flux_pre[i] != NULL)
			{
				delete [] _flux_pre[i];
			}
	}

	if (_flux_pre != NULL)
			{
				delete [] _flux_pre;
			}*/
}//should change delete in three dimension 

void Solver::setQuadrature(Quadrature* quad) {
	_quadrature = quad;
}

void Solver::startIteration() {
	int num = 0;
	int i = 0;
	_keff = 1.0;
	InitializeFlux();
    std::cout << "Flux array created!!!!!--!..." << std::endl;
    InputAnVoidFlux();    
	InitializeSource();
        std::cout << "Flux array created!!!..." << std::endl;
	externalsource();
    std::cout << "Flux array created!!!!!!..." << std::endl;    
	//storeFlux();//if there should be a storeFlux?????
	std::cout << "Flux array created!!!!!22222!..." << std::endl; 
	//0();
	//computeTotalSource();
	computeFissionSource();




	for (int i=0; i<_max_iter; i++) {
		computeTotalSource();//actually this function only add fission source

		Quadrat_1();
		
		Quadrat_2();

		Quadrat_3();
		std::cout << "Flux array created!!!!!!11111..." << std::endl;
		Quadrat_4();
		std::cout << "Flux array created!!!!!!11111..." << std::endl;
		sourceIteration();
		
		externalsource();
		storeFlux();
		
		computeKeff();//important
		computeResidual();
		




		//normalize();
		//computeTotalSource();
		//for (int j=0; j<_n/2; j++)
			//sweepMuMinus(j);
		//for (int j=_n/2-1; j>=0; j--)
			//sweepMuPlus(j);	
		//computeKeff();
		//computeResidual();
		//storeFlux();
		num++;
		if (i>1 && _res_flux < _converge_thresh)
			break;
		std::cout << "# " << num << '\t' << "resflux = " << _res_flux << std:: endl;
		std::cout << "# " << num << '\t' << "Keff = " << _keff << std:: endl;
	}

	std::cout << "Scalar Flux: " << std::endl;
	for (int i=0; i<_Nx; i++) {
		for (int j=0; j<_Ny; j++){
			std::cout << "Mesh #" <<i<<j<< " = " << _flux[i][j] << '\t';
		if ((i+1)%4 ==0)
			std::cout << '\n';
		}
	}
}

void Solver::InputAnVoidFlux(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){	
		  _flux[i][j]=1.0;	
		}
	}
}

void Solver::computeFissionSource(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_fiss_source[i][j][m]=(_nu_sigma_f/_keff)*_flux[i][j]/(4*pi);
			}
		}
	}
}

void Solver::computeKeff(){
	double old_keff = _keff;
	double old_fiss_source = 0.0;
	double new_fiss_source = 0.0;
	for (int i=0; i<_Ny; i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				old_fiss_source =old_fiss_source + _fiss_source[i][j][m];
			}
		}
	}
	computeFissionSource();
	for (int i=0; i<_Ny; i++){
		for (int j=0;j<_Nx; j++){
			for (int m=0;m<_M4;m++){
				new_fiss_source =new_fiss_source + _fiss_source[i][j][m];
			}
		}
	}
	_keff = _keff*new_fiss_source/old_fiss_source;
	_res_keff = fabs((_keff-old_keff)/_keff);
	std::cout << "Show old fission source" <<'\t'<<old_fiss_source<< std::endl;
	std::cout << "Show new fission source" <<'\t'<<new_fiss_source<< std::endl;
}


void Solver::computeTotalSource(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m]=_tot_source[i][j][m]+_fiss_source[i][j][m];
			}
		}
	}
}

void Solver::InitializeFlux(){
	std::cout << "Flux array created..." << std::endl;
	_flux_x = (double***) new double**[_Ny + 1];
	for (int i=0;i<_Ny+1;i++){
		_flux_x[i] = (double**) new double* [_Nx];
	}
	for (int i=0;i<_Ny+1;i++){
		for (int j=0;j<_Nx;j++){
			_flux_x[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny+1;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_flux_x[i][j][m]=0.0;
			}
		}
	}

	_flux_y = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_y[i] = (double**) new double* [_Nx+1];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx+1;j++){
			_flux_y[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx+1;j++){
			for (int m=0;m<_M4;m++){
				_flux_y[i][j][m]=0.0;
			}
		}
	}


	_flux_c = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_c[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_flux_c[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_flux_c[i][j][m]=0.0;
			}
		}
	}


	_flux = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux[i] = new double [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
		
			_flux[i][j]=0.0;
		
		}
	}


	_flux_pre = (double**) new double*[_Ny];
	for (int i=0;i<_Ny;i++){
		_flux_pre[i] = new double [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
		
			_flux_pre[i][j]=0.0;
		
		}
	}


}

void Solver::InitializeSource(){
	_tot_source = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_tot_source[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_tot_source[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m]=0.0;
			}
		}
	}
	
	_fiss_source = (double***) new double**[_Ny];
	for (int i=0;i<_Ny;i++){
		_fiss_source[i] = (double**) new double* [_Nx];
	};
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_fiss_source[i][j] = new double[_M4];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_fiss_source[i][j][m]=0.0;
			}
		}
	}
}

void Solver::externalsource(){
	for (int i=_idy1-1;i<_idy2;i++){
		for (int j=_idx1-1;j<_idx2;j++){
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m]=_tot_source[i][j][m]+_S0/(2*pi);
			}

		}
	}
}

void Solver::Quadrat_1(){
	for (int j=0;j<_Nx;j++){
		for (int m=0;m<_M1;m++){
			_flux_x[0][j][m]=0.0;//there is 0
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int m=0;m<_M1;m++){
			_flux_y[i][0][m]=0.0;
		}
	}
	for (int m=0;m<_M1;m++){
		double E = fabs(_quadrature->getMux(m))*(2 * _deltay);	
		double F = fabs(_quadrature->getMuy(m))*(2 * _deltay);
		double C = E + F + _sigma_t*_deltax*_deltay;
		for (int i=0;i<_Ny;i++){
			for (int j=0;j<_Nx;j++){
				_flux_c[i][j][m]=(E*_flux_y[i][j][m]+F*_flux_x[i][j][m]+_tot_source[i][j][m]*_deltax*_deltay)/C;
				_flux_x[i+1][j][m]=2*_flux_c[i][j][m]-_flux_x[i][j][m];
				_flux_y[i][j+1][m]=2*_flux_c[i][j][m]-_flux_y[i][j][m];
				if (_flux_x[i+1][j][m]<0){
					_flux_x[i+1][j][m]=0.0;
				}
				if (_flux_y[i][j+1][m]<0)
				{
					_flux_y[i][j+1][m]=0.0;
				}
			}
		}
	}
}

void Solver::Quadrat_2(){
	for (int j=0;j<_Nx;j++){
		for (int m=_M1;m<_M2;m++){
			_flux_x[0][j][m]=0.0;
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int m=_M1;m<_M2;m++){
			_flux_y[i][_Nx][m]=0.0;
		}
	}
	for (int m=_M1;m<_M2;m++){
		double E=fabs(_quadrature->getMux(m))*(2*_deltay);
		double F=fabs(_quadrature->getMuy(m))*(2*_deltax);
		double C=E+F+_sigma_t*_deltax*_deltay;
		for (int i=0;i<_Ny;i++){
			for (int j=_Nx-1;j>-1;j--){
				_flux_c[i][j][m]=(F*_flux_x[i][j][m]+E*_flux_y[i][j+1][m]+_tot_source[i][j][m]*_deltax*_deltay)/C;
				_flux_x[i+1][j][m]=2*_flux_c[i][j][m]-_flux_x[i][j][m];
				_flux_y[i][j][m]=2*_flux_c[i][j][m]-_flux_y[i][j+1][m];
				if (_flux_x[i+1][j][m]<0)
				{
					_flux_x[i+1][j][m]=0.0;
				}
				if (_flux_y[i][j][m]<0)
				{
					_flux_y[i][j][m] = 0.0;
				}
			}
		}
	}
}

void Solver::Quadrat_3(){
	for (int j=0;j<_Nx;j++){
		for (int m=_M2;m<_M3;m++){
			_flux_x[_Ny][j][m]=0.0;
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int m=_M2;m<_M3;m++){
			_flux_y[i][_Nx][m]=0.0;
		}
	}

	for (int m=_M2;m<_M3;m++){
		double E=fabs(_quadrature->getMux(m))*(2*_deltay);
		double F=fabs(_quadrature->getMuy(m))*(2*_deltax);
		double C=E+F+_sigma_t*_deltax*_deltay;

		for (int i=_Ny-1;i>-1;i--){
			for (int j=_Nx-1;j>-1;j--){
				_flux_c[i][j][m]=(F*_flux_x[i+1][j][m]+E*_flux_y[i][j+1][m]+_tot_source[i][j][m]*_deltax*_deltay)/C;
				_flux_x[i][j][m]=2*_flux_c[i][j][m]-_flux_x[i+1][j][m];
				_flux_y[i][j][m]=2*_flux_c[i][j][m]-_flux_y[i][j+1][m];
				if (_flux_x[i][j][m]<0)
				{
					_flux_x[i][j][m]=0.0;
				}
				if (_flux_y[i][j][m]<0)
				{
					_flux_y[i][j][m] = 0.0;
				}
			}
		}
	}
}

void Solver::Quadrat_4(){
	for (int j=0;j<_Nx;j++){
		for (int m=_M3;m<_M4;m++){
			_flux_x[_Ny][j][m]=0.0;
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int m=_M3;m<_M4;m++){
			_flux_y[i][0][m]=0.0;
		}
	}
	for (int m=_M3;m<_M4;m++){
		double E=fabs(_quadrature->getMux(m))*(2*_deltay);
		double F=fabs(_quadrature->getMuy(m))*(2*_deltax);
		double C=E+F+_sigma_t*_deltax*_deltay;
		for (int i=_Ny-1;i>-1;i--){
			for (int j=0;j<_Nx;j++){
				_flux_c[i][j][m]=(F*_flux_x[i+1][j][m]+E*_flux_y[i][j][m]+_tot_source[i][j][m]*_deltax*_deltay)/C;
				_flux_x[i][j][m]=2*_flux_c[i][j][m]-_flux_x[i+1][j][m];
				_flux_y[i][j+1][m]=2*_flux_c[i][j][m]-_flux_y[i][j][m];
				if (_flux_x[i][j][m]<0)
				{
					_flux_x[i][j][m]=0.0;
				}
				if (_flux_y[i][j+1][m]<0)
				{
					_flux_y[i][j+1][m] = 0.0;
				}
			}
		}
	}
}

void Solver::sourceIteration(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m] = 0.0;
			}
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			double C00 = 0.0;
			double C1_1 = 0.0;
			double C11 = 0.0;
			double C2_2 = 0.0;
			double C20 = 0.0;
			double C22 = 0.0;
			for (int m=0;m<_M4;m++){
				C00=C00+_quadrature->getY00(m)*_flux_c[i][j][m]/_M4;
				C1_1=C1_1+_quadrature->getY1_1(m)*_flux_c[i][j][m]/_M4;
				C11=C11+_quadrature->getY11(m)*_flux_c[i][j][m]/_M4;
				C2_2=C2_2+_quadrature->getY2_2(m)*_flux_c[i][j][m]/_M4;
				C20=C20+_quadrature->getY20(m)*_flux_c[i][j][m]/_M4;
				C22=C22+_quadrature->getY22(m)*_flux_c[i][j][m]/_M4;
			}
			for (int m=0;m<_M4;m++){
				_tot_source[i][j][m]=_sigma_s0*_quadrature->getY00(m)*C00+3*_sigma_s1*(_quadrature->getY11(m)*C11+_quadrature->getY1_1(m)*C1_1)+5*_sigma_s2*(_quadrature->getY2_2(m)*C2_2/12+_quadrature->getY20(m)*C20/12+_quadrature->getY22(m)*C22/12);
				_tot_source[i][j][m]=_tot_source[i][j][m]/(2*pi);
			}
		}
	}
}

void Solver::storeFlux(){
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_flux_pre[i][j]=_flux[i][j];
		}
	}
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			_flux[i][j]=0.0;
		}
	}
	double XX=0.0;
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			for (int m=0;m<_M4;m++){
				XX=XX+4*pi*_flux_c[i][j][m]/_M4;
			}
			_flux[i][j]=XX;
			XX=0.0;
		}
	}
}

void Solver::computeResidual(){
	double max_scalarflux = 0.0;
	double mm;
	for (int i=0;i<_Ny;i++){
		for (int j=0;j<_Nx;j++){
			mm = fabs(_flux[i][j]-_flux_pre[i][j]);
			if (mm > max_scalarflux)
					{
						max_scalarflux = mm;
					}		
		 }
	}
	_res_flux = max_scalarflux;
}


#endif
