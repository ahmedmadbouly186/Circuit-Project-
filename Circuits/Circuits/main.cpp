#include<iostream>
#include<string>
#include<fstream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<complex>
#include"vsrc.h"
#include"Capacaitor.h"
#include"Resistance.h"
#include"Inductor.h"
#include"isrc.h"
#include"Branch.h"
using namespace std;
using namespace Eigen;
using namespace std::complex_literals;
#define n_max 10
Component* complist[n_max] = { NULL };
Component* searcher(string x);
int counter = 0;


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

MatrixXcd CalculateVoltNode(MatrixXcd G, MatrixXcd I);
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
	Eigen::MatrixXcd s;
	Eigen::Matrix3cd f;
	f <<  c1,c2,c3, 7.0 +5i, 8, 9,10,11,12;
	cout << f<<endl;

	/*
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

	outputfile.close();
	inputfile.close();
	*/
	int n = 2;
	Branch** nbranch = new Branch * [n];
	complex<double> zero;
	zero.real(0);
	zero.imag(0);
	for (int i = 0; i < n; i++)
	{
		if (nbranch[i]->getcurrent().real() == zero.real() && nbranch[i]->getcurrent().imag() == zero.imag() && nbranch[i]->getz().real() > zero.real() 
			&& nbranch[i]->getz().imag() > zero.imag() && nbranch[i]->getvolt().getcomplex().real() > zero.real() && nbranch[i]->getvolt().getcomplex().imag() > zero.imag()) {
			Branch q = *nbranch[i];
			simplification(nbranch[i]->getvolt(), nbranch[i]->getz());
			nbranch[i]->getvolt().zerosetters();
		}
	}

	///////////////////hany
	int num_non_sim_nodes = 3;
	int num = num_non_sim_nodes - 1;
	Node* arr[2] = { 1,2 };
	//Eigen::MatrixXd G(num, num);
	Eigen::MatrixXcd m(num, num);
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < num; j++)
		{
			m(i, j).real(0);
			m(i, j).imag(0);
		}
	}
	Eigen::MatrixXcd I(1,num);
	/*m(0, 0).real(3);
	m(0, 0).imag(5);
	m(1, 0).real(3);
	m(1, 0).imag(5);
	m(2, 0).real(3);
	m(2, 0).imag(5);
	m(2, 2).real(3);
	m(2, 2).imag(5);
	cout << m;*/
	for (int j = 0; j < num; j++)
	{
		for (int k = 0; k < num; k++)
		{
			complex <double> Admittance_s(0, 0);
			complex <double> Admittance_d(0, 0);
			for (int i = 0; i < counter; i++)
			{
				Component* comp = complist[i];
				Resistance* resistance = dynamic_cast<Resistance*>(comp);
				Inductor* inductor = dynamic_cast<Inductor*>(comp);
				Capacaitor* capacaitor = dynamic_cast<Capacaitor*>(comp);
				if (resistance != NULL)
				{
					if ((resistance->get_node1() == arr[k]) || (resistance->get_node2() == arr[k]))
					{
						complex <double> r;
						r = 1.0 / resistance->get_Impedance();
						Admittance_s += r;
					}
					if ((resistance->get_node1() == arr[j]) && (resistance->get_node2() == arr[k]))
					{
						complex <double> r;
						r = 1.0 / resistance->get_Impedance();
						Admittance_d += r;
					}
				}
				else if (inductor != NULL)
				{
					if ((inductor->get_node1() == arr[k]) || (inductor->get_node2() == arr[k]))
					{
						complex <double> in;
						in = inductor->get_Admittance();
						Admittance_s += in;
					}
					if ((inductor->get_node1() == arr[j]) && (inductor->get_node2() == arr[k]))
					{
						complex <double> r;
						r = 1.0 / inductor->get_Impedance();
						Admittance_d += r;
					}
				}
				else if (capacaitor != NULL)
				{
					if ((capacaitor->get_node1() == arr[k]) || (capacaitor->get_node2() == arr[k]))
					{
						complex <double> ca;
						ca = capacaitor->get_Admittance();
						Admittance_s += ca;
					}
					if ((capacaitor->get_node1() == arr[j]) && (capacaitor->get_node2() == arr[k]))
					{
						complex <double> r;
						r = 1.0 / capacaitor->get_Impedance();
						Admittance_d += r;
					}

				}
			}
			if (j == k)
			{
				m(k, k) = Admittance_s;
			}
			else
			{
				m(j, k) = -Admittance_d;
			}
		}
	}
	MatrixXcd V = CalculateVoltNode(m , I); // function calculate volt of each node
	Node* nodes = new Node[num];  // should be created above 
	for (int i = 0; i < num; i++)
	{
		nodes[i].setvoltage(V(1, i));
	}

	Branch* branches = new Branch[num]; // should be created above  
	int number_of_branchess; // should be created above 
	for (int i = 0; i < number_of_branchess; i++)
	{
		branches->Calculatecurrent();
	}
	// now we have voltage od each node and current in each branch
	
	return 0;
}

MatrixXcd CalculateVoltNode(MatrixXcd G, MatrixXcd I)
{
	MatrixXcd V;
	V = G.inverse() * I;
	return V;
}

Component* searcher(string x)
{
	Component* comp = NULL;
	for (int i = 0; i < counter; i++)
	{
		string name = complist[i]->get_name();
		if ((name.compare(x) == 0))
		{
			comp = complist[i];
		}
	}
	return comp;
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
	
	