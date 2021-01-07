#include<iostream>
#include<string>
#include<fstream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<complex>
#include"vsrc.h"
#include"isrc.h"
#include"Branch.h"
using namespace std;
using namespace Eigen;
using namespace std::complex_literals;


void simplification(Vsrc a, complex<double> b) 
{
	complex<double> x = (a.getcomplex()) / b;
	double p = x.real() * x.real() +x.imag() * x.imag();;
	double mag = sqrt(p);
	double phase = tanf(x.imag() / x.real());
	Isrc c (a.getname(), a.getnode1(), a.getnode2(), phase, mag);
	Branch q((a.getnode1()), a.getnode2());
	q.setcurrent(c.getcomplex());
	c.zerosetters();
}


enum component
{
	w,
	vsrc,
	isrc,
	vcvs,
	ccvs,
	vccs,
	cccs,
	res,
	cap,
	ind
};

int main()
{
	double c,a = 2;
	double b = 2;
	complex <double> c2(a,b);
	c = c2.real()+1;
	c2.imag(b);
	complex<double> c1;
	c1.imag(1);
	c1.real(0);
	complex<double> c3;
	c3 = c2 / c1;
	Eigen::MatrixXcd m;
	Eigen::Matrix3cd f;
	f <<  c1,c2,c3, 7.0 +5i, 8, 9,10,11,12;
	cout << f<<endl;

	

	ifstream inputfile;
	ofstream outputfile;
	outputfile.open("circuit.txt");
	inputfile.open("netlist.txt");
	if (!(inputfile.is_open()))
	{
		cout << "file is not opened , try again" << endl;
		return 0;
	}
	string name;
	while (!(inputfile.eof()))
	{
		inputfile >> name;
		if (name== "w")
		{
			double w;
			inputfile >> w;
			outputfile << name << " " << w << endl;
			
		}
		else if (name == "vsrc")
		{
			string Name;
			int n1, n2;
			double phase,magintude;
			inputfile >> Name >> n1 >> n2 >> magintude >> phase;
			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << magintude << " " << phase << endl;
		}
		else if (name == "isrc")
		{

			string Name;
			int n1, n2;
			double phase, magintude;
			inputfile >> Name >> n1 >> n2 >> magintude >> phase;

			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << magintude << " " << phase << endl;
		}
		else if (name == "vcvs")
		{
			string Name;
			int n1, n2;
			int n3, n4;
			string Name2;
			double val;
			inputfile >> Name >> n1 >> n2 >> n3 >> n4 >> Name2 >> val;
			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << n3 << " " << n4 <<" "<< Name2<<" "<<val<< endl;
		}
		else if (name == "ccvs")
		{
			string Name;
			int n1, n2;
			int n3, n4;
			string Name2;
			double val;
			inputfile >> Name >> n1 >> n2 >> n3 >> n4 >> Name2 >> val;
			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << n3 << " " << n4 << " " << Name2 << " " << val << endl;

		}
		else if (name == "vccs")
		{
			string Name;
			int n1, n2;
			int n3, n4;
			string Name2;
			double val;
			inputfile >> Name >> n1 >> n2 >> n3 >> n4 >> Name2 >> val;
			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << n3 << " " << n4 << " " << Name2 << " " << val << endl;

		}
		else if (name == "cccs")
		{
			string Name;
			int n1, n2;
			int n3, n4;
			string Name2;
			double val;
			inputfile >> Name >> n1 >> n2 >> n3 >> n4 >> Name2 >> val;
			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << n3 << " " << n4 << " " << Name2 << " " << val << endl;

		}
		else if (name == "res")
		{
			string Name;
			int n1, n2;
			double val;
			inputfile >> Name >> n1 >> n2 >> val;
			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << val << endl;

		}
		else if (name == "cap")
		{
			string Name;
			int n1, n2;
			double val;
			inputfile >> Name >> n1 >> n2 >> val;
			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << val << endl;
		}
		else if (name == "ind")
		{
			string Name;
			int n1, n2;
			double val;
			inputfile >> Name >> n1 >> n2 >> val;
			outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << val << endl;
		}
		else
		{

		}
	}

	int n = 2;
	Branch** nbranch = new Branch * [n];
	complex<double> zero;
	zero.real(0);
	zero.imag(0);
	for (int i = 0; i < n; i++)
	{
		if (nbranch[i]->getcurrent().real() == zero.real() && nbranch[i]->getcurrent().imag() == zero.imag() && nbranch[i]->getz().real() > zero.real() && nbranch[i]->getz().imag() > zero.imag() && nbranch[i]->getvolt().getcomplex().real() > zero.real() && nbranch[i]->getvolt().getcomplex().imag() > zero.imag()) {
			Branch q = *nbranch[i];
			simplification(nbranch[i]->getvolt(), nbranch[i]->getz());
			nbranch[i]->getvolt().zerosetters();
		}
	}

	outputfile.close();
	inputfile.close();
	return 0;
}
	/*
			switch (component)
			{
			case w:
				break;
			case vsrc:
				break;
			case isrc:
				break;
			case vcvs:
				break;
			case ccvs:
				break;
			case vccs:
				break;
			case cccs:
				break;
			case res:
				break;
			case cap:
				break;
			case ind:
				break;
			default:
				break;
			}
		}*/
	
	