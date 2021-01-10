#include<iostream>
#include <cmath>
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
#include "Node.h"
#include"Dependent_sources.h"
#include"Current_C_Current.h"
using namespace std;
using namespace Eigen;
using namespace std::complex_literals;
#include"Component.h"
double W;
int searcher(string x);
#define n_max 10
Component* complist[n_max] = { NULL };
int comcount = 0;
/////////////nodes&branches/////////////
int nodecount = 0;
int nonsimple_nodecount = 0;
int simple_nodecount = 0;
int branchcount = 0;
int nodes_rannk[10] = { -1 ,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
Node** Nodelist;
Node** nonsimple_Nodelist;
Node** simple_Nodelist;
Branch** fullBranchlist;
Branch** Branchlist;
/////////////nodes&branches/////////////
Isrc* simplification(Vsrc* a, complex<double> b, Node* N1, Node* N2)
{
	complex<double> x = (a->getcomplex()) / b;
	double p = x.real() * x.real() + x.imag() * x.imag();
	double mag = sqrt(p);
	double phase = atan(x.imag() / x.real()) * 180 / 3.14;
	Isrc* c = NULL;
	if(N2->getrank() == a->get_node2() || N1->getrank() == a->get_node1())
	{
		c = new Isrc(a->get_name(), N1->getrank(), N2->getrank(), phase, mag);
	}
	else
	{
		c = new Isrc(a->get_name(), N2->getrank(), N1->getrank(), phase, mag);

	}
	//Branch q(Nodelist[a->get_node1()], Nodelist[a->get_node2()]);
	Branchlist[branchcount] = new Branch(N1, N2);
	Branchlist[branchcount]->setcurentsource(c);
	//q.setcurrent(c->getcomplex());
	//	c->zerosetters();
	branchcount++;
	return c;
}
MatrixXcd CalculateVoltNode(MatrixXcd G, MatrixXcd I);

int main()
{
	ifstream inputfile;
	inputfile.open("netlist3.txt");
	if (!(inputfile.is_open()))
	{
		cout << "file is not opened , try again" << endl;
		return 0;
	}
	string name;

	while (!(inputfile.eof()))
	{
		int n1 = -1, n2 = -1;
		inputfile >> name;
		if (name == "w"|| name =="W")
		{
			double w;
			inputfile >> w;
			W = w;
		}
		else if (name == "vsrc" || name == "Vsrc" || name == "VSrc")
		{
			string Name;
			double phase, magintude;
			inputfile >> Name >> n1 >> n2 >> magintude >> phase;
			Vsrc* v1 = new Vsrc(Name, n1, n2, phase, magintude);
			complist[comcount++] = v1;
		}
		else if (name == "isrc" || name == "Isrc" || name == "ISrc")
		{
			string Name;
			double phase, magintude;
			inputfile >> Name >> n1 >> n2 >> magintude >> phase;
			Isrc* v1 = new Isrc(Name, n1, n2, phase, magintude);
			complist[comcount++] = v1;
		}
		else if (name == "vcvs" || name == "Vcvs" || name == "VCvs")
		{
			/*string Name;
			int n3, n4;
			string Name2;
			double val;
			inputfile >> Name >> n1 >> n2 >> n3 >> n4 >> Name2 >> val;*/
			//Vcvs* v1 = new Vcvs(Name, n1, n2, phase, magintude);
			//complist[comcount++] = v1;
			//outputfile << name << " " << Name << " " << n1 << " " << n2 << " " << n3 << " " << n4 <<" "<< Name2<<" "<<val<< endl;
		}
		else if (name == "cccs" || name == "Cccs" || name == "CCcs")
		{
			string Name;
			int n3, n4;
			string Name2;
			double val;
			inputfile >> Name >> n1 >> n2 >> n3 >> n4 >> Name2 >> val;
			Current_C_Current* cccs = new Current_C_Current(n1,n2,Name,n3,n4,Name2,val);
			complist[comcount++] = cccs;
		}
		else if (name == "res")
		{
			string Name;
			double val;
			inputfile >> Name >> n1 >> n2 >> val;
			Resistance* v1 = new Resistance(val, Name, n1, n2);
			complist[comcount++] = v1;
		}
		else if (name == "cap")
		{
			string Name;
			double val;
			inputfile >> Name >> n1 >> n2 >> val;
			Capacaitor* v1 = new Capacaitor(val, W, Name, n1, n2);
			complist[comcount++] = v1;
		}
		else if (name == "ind")
		{
			string Name;
			double val;
			inputfile >> Name >> n1 >> n2 >> val;
			Inductor* v1 = new Inductor(val, W, Name, n1, n2);
			complist[comcount++] = v1;
		}
		bool findN1 = false;
		bool findN2 = false;
		for (int i = 0; i < 10; i++)
		{
			if (n1 == nodes_rannk[i])
			{
				findN1 = true;
			}
			if (n2 == nodes_rannk[i])
			{
				findN2 = true;
			}
		}
		if (findN1 == false)
		{
			nodes_rannk[nodecount++] = n1;
		}
		if (findN2 == false)
		{
			nodes_rannk[nodecount++] = n2;
		}

	}
	inputfile.close();
	/////////////////////////////////// Sabry ////////////////////////////////////////////
	Current_C_Current* cccs = NULL;
	for (int i = 0; i < comcount; i++)
	{
		cccs = dynamic_cast<Current_C_Current*>(complist[i]);
		if (cccs != NULL)
		{
			break;
		}
	}
	if (cccs != NULL)
	{
		for (int j = 0; j < comcount; j++)
		{
			if (   cccs->get_Noded1() == complist[j]->get_node1() && cccs->get_Noded2() == complist[j]->get_node2()
				|| cccs->get_Noded2() == complist[j]->get_node1() && cccs->get_Noded1() == complist[j]->get_node2())
			{
				Resistance* r1 = dynamic_cast <Resistance*>(complist[j]);
				Capacaitor* c1 = dynamic_cast <Capacaitor*>(complist[j]);
				Inductor* in1 = dynamic_cast <Inductor*>(complist[j]);
				if (r1 != NULL)
				{
					cccs->set_Impedance(r1->get_Impedance());
				}
				else if(c1 !=NULL)
				{
					cccs->set_Impedance(c1->get_Impedance());
				}
				else if(in1 !=NULL)
				{
					cccs->set_Impedance(in1->get_Impedance());
				}

			}
		}
	}
	////////////////////madbouly el gamed statr /////////////////////////
	Nodelist = new Node * [nodecount];
	for (int i = 0; i < nodecount; i++)
	{
		Nodelist[i] = new Node(i);
	}

	for (int i = 0; i < comcount; i++)
	{
		int N1 = complist[i]->get_node1();
		int N2 = complist[i]->get_node2();
		int indexnode1, indexnode2;
		for (int j = 0; j < nodecount; j++)
		{
			if (Nodelist[j]->getrank() == N1)
			{
				indexnode1 = j;
			}
			else if (Nodelist[j]->getrank() == N2)
			{
				indexnode2 = j;
			}
		}
		Nodelist[indexnode1]->increment_connections(Nodelist[indexnode2]);
		Nodelist[indexnode2]->increment_connections(Nodelist[indexnode1]);
	}
	for (int j = 0; j < nodecount; j++)
	{
		if (!(Nodelist[j]->ifsimple()))
		{
			nonsimple_nodecount++;
		}
	}
	simple_nodecount = nodecount - nonsimple_nodecount;
	nonsimple_Nodelist = new Node * [nonsimple_nodecount];
	simple_Nodelist = new Node * [simple_nodecount];

	{
		int nonsimple = 0; int simple = 0;
		for (int i = 0; i < nodecount; i++)
		{
			if (Nodelist[i]->ifsimple())
			{
				simple_Nodelist[simple++] = Nodelist[i];
			}
			else
			{
				nonsimple_Nodelist[nonsimple++] = Nodelist[i];
			}
		}
	}

	////////////////////////////// creat branches /////////////////////////////////////
	//fullBranchlist = new Branch * [comcount];
	branchcount = comcount - simple_nodecount;
	Branchlist = new Branch * [comcount];
	for (int i = 0; i < comcount; i++)
	{
		Branchlist[i] = NULL;
	}
	{
		int counter = 0;
		for (int i = 0; i < comcount; i++)
		{
			int N1 = complist[i]->get_node1();
			int N2 = complist[i]->get_node2();
			//fullBranchlist[i] = new Branch(Nodelist[N1], Nodelist[N2]);
			bool find = false;
			for (int j = 0; j < simple_nodecount; j++)
			{
				if (N1 == simple_Nodelist[j]->getrank() || N2 == simple_Nodelist[j]->getrank())
				{
					find = true;
				}
			}
			if (find == false)
			{
				Vsrc* v1 = dynamic_cast <Vsrc*>(complist[i]);
				Isrc* i1 = dynamic_cast <Isrc*>(complist[i]);
				Current_C_Current* D_I = dynamic_cast<Current_C_Current*>(complist[i]);
				Resistance* r1 = dynamic_cast <Resistance*>(complist[i]);
				Capacaitor* c1 = dynamic_cast <Capacaitor*>(complist[i]);
				Inductor* in1 = dynamic_cast <Inductor*>(complist[i]);

				Branchlist[counter] = new Branch(Nodelist[N1], Nodelist[N2]);
				if (v1 != NULL)
				{
					Branchlist[counter]->setvolt(v1);
				}
				else if (i1 != NULL)
				{
					Branchlist[counter]->setcurentsource(i1);
				}
				else if (r1 != NULL)
				{
					Branchlist[counter]->setz(r1->get_Impedance());
				}
				else if (c1 != NULL)
				{
					Branchlist[counter]->setz(c1->get_Impedance());
				}
				else if (in1 != NULL)
				{
					Branchlist[counter]->setz(in1->get_Impedance());
				}
				else if(D_I != NULL)
				{
					Branchlist[counter]->setDependentCsource(D_I);
				}
				counter++;
			}
		}

		for (int i = 0; i < simple_nodecount; i++)
		{
			Node** connection = simple_Nodelist[i]->getNodeconnections();
			int N1 = connection[0]->getrank();        // the node that before the simplenode that we work on
			int N = simple_Nodelist[i]->getrank();   // the simplenode that we work on
			int N2 = connection[1]->getrank();       // the node that after the simplenode that we work on
			Branchlist[counter] = new Branch(Nodelist[N1], Nodelist[N2]);
			bool find = false;
			for (int j = 0; j < comcount; j++)
			{
				int N3 = complist[j]->get_node1();
				int N4 = complist[j]->get_node2();
				if ((N3 == N && N4 == N1) || (N3 == N1 && N4 == N) || (N3 == N && N4 == N2) || (N3 == N2 && N4 == N))
				{
					find = true;
				}


				if (find == true)
				{
					Vsrc* v1 = dynamic_cast <Vsrc*>(complist[j]);
					Isrc* i1 = dynamic_cast <Isrc*>(complist[j]);

					Resistance* r1 = dynamic_cast <Resistance*>(complist[j]);
					Capacaitor* c1 = dynamic_cast <Capacaitor*>(complist[j]);
					Inductor* in1 = dynamic_cast <Inductor*>(complist[j]);
					if (v1 != NULL)
					{
						Branchlist[counter]->setvolt(v1);
					}
					else if (i1 != NULL)
					{
						Branchlist[counter]->setcurentsource(i1);
					}
					else if (r1 != NULL)
					{
						Branchlist[counter]->setz(r1->get_Impedance());
					}
					else if (c1 != NULL)
					{
						Branchlist[counter]->setz(c1->get_Impedance());
					}
					else if (in1 != NULL)
					{
						Branchlist[counter]->setz(in1->get_Impedance());
					}

				}
				find = false;
			}

			// this code work when N2 is also a simple Node  
			// then we should add the component between N2 and the next Node to the same btanch
			// and Repeat this step if the Node that next to N2 is alse asimple
			// and so , if all nodes is asimple Node we must make exactly one branch
			// this branch will take all the components
			for (int k = i; k < simple_nodecount; k++)
			{
				N2;// ÚÇíÒ ÇáÝ ÚáíåÇ ÇÊÇßÏ åíÇ simple æáÇ áÇ
				bool find2 = false;
				for (int k1 = i; k1 < simple_nodecount; k1++)
				{
					if (N2 == simple_Nodelist[k1]->getrank())
					{
						find2 = true;
					}
				}
				if (find2 == true)
				{
					i++;
					N = N2;
					connection = simple_Nodelist[i]->getNodeconnections();
					N2 = connection[1]->getrank();
					bool f = false;
					for (int k2 = 0; k2 < comcount; k2++)
					{
						int N3 = complist[k2]->get_node1();
						int N4 = complist[k2]->get_node2();
						if ((N3 == N && N4 == N2) || (N3 == N2 && N4 == N))
						{
							f = true;
						}
						if (f == true)
						{
							Vsrc* v1 = dynamic_cast <Vsrc*>(complist[k2]);
							Isrc* i1 = dynamic_cast <Isrc*>(complist[k2]);

							Resistance* r1 = dynamic_cast <Resistance*>(complist[k2]);
							Capacaitor* c1 = dynamic_cast <Capacaitor*>(complist[k2]);
							Inductor* in1 = dynamic_cast <Inductor*>(complist[k2]);
							if (v1 != NULL)
							{
								Branchlist[counter]->setvolt(v1);
							}
							else if (i1 != NULL)
							{
								Branchlist[counter]->setcurentsource(i1);
							}
							else if (r1 != NULL)
							{
								Branchlist[counter]->setz(r1->get_Impedance());
							}
							else if (c1 != NULL)
							{
								Branchlist[counter]->setz(c1->get_Impedance());
							}
							else if (in1 != NULL)
							{
								Branchlist[counter]->setz(in1->get_Impedance());
							}
						}
					}
				}
				else
				{

				}

			}
			counter++;
		}
	}

	/////////////////////////madbouly el gamed end //////////////////////////

	complex<double> zero;
	zero.real(0);
	zero.imag(0);
	for (int i = 0; i < branchcount; i++)
	{
		cout << Branchlist[i]->getnode1()->getrank();
		cout << Branchlist[i]->getnode2()->getrank();

	}
	Isrc* c;
	for (int i = 0; i < branchcount; i++)
	{
		Branch* b1 = Branchlist[i];
		if (Branchlist[i] != NULL)
		{
			if ((Branchlist[i]->getvolt() != NULL))
			{

				/*if (Branchlist[i]->getcsource()->getcomplex().real() == zero.real() && Branchlist[i]->getcsource()->getcomplex().imag() == zero.imag() &&
					(Branchlist[i]->getz().real() > zero.real() || Branchlist[i]->getz().imag() > zero.imag())
					&& Branchlist[i]->getvolt()->getcomplex().real() > zero.real() && Branchlist[i]->getvolt()->getcomplex().imag() > zero.imag())*/
				if (Branchlist[i]->getcsource() == NULL && Branchlist[i]->getz() != zero && Branchlist[i]->getvolt()->getcomplex() != zero)
				{
					//Branch q = *Branchlist[i];
					int n = searcher(Branchlist[i]->getvolt()->get_name());
					c = simplification(Branchlist[i]->getvolt(), Branchlist[i]->getz(), Branchlist[i]->getnode1(), Branchlist[i]->getnode2());
					complist[n] = c;
					Branchlist[i]->getvolt()->zerosetters();

				}
			}

		}
	}
	/////////////////////////////////////Ahmed hany
	int num_non_sim_nodes = nonsimple_nodecount;
	int num = num_non_sim_nodes - 1;
	/*Node n1, n2;
	Node* arr[2] = { &n1,&n2 };*/

	int* arr;
	arr = new int[nonsimple_nodecount - 1];
	for (int i = 1; i < nonsimple_nodecount; i++)
	{
		arr[i - 1] = nonsimple_Nodelist[i]->getrank();
	}
	cout << arr[0] << " " << arr[1] << " " << arr[2] << endl;
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
	Eigen::MatrixXcd I(num, 1);
	for (int i = 0; i < num; i++)
	{
		I(i, 0) = 0;
	}
	/*for (int j = 1; j < num; j++)
	{
		for (int k = 0; k < num; k++)
		{
			complex <double> Admittance_s(0, 0);
			complex <double> Admittance_d(0, 0);
			for (int i = 0; i < branchcount; i++)
			{
				complex <double> admit = Branchlist[i]->getz();
				if ((Branchlist[i]->getnode1()->getrank() == arr[k]) || (Branchlist[i]->getnode2()->getrank() == arr[k]))
				{
					complex <double> r;
					r = Branchlist[i]->get_Admittance();
					Admittance_s += r;
				}
			}
			for (int l = 0; l < branchcount; l++)
			{
				if (((Branchlist[l]->getnode1()->getrank() == arr[j - 1]) || (Branchlist[l]->getnode2()->getrank() == arr[j - 1])) &&
					((Branchlist[l]->getnode1()->getrank() == arr[j])) || (Branchlist[l]->getnode2()->getrank() == arr[j]))
				{
					Admittance_d += Branchlist[l]->get_Admittance();
				}
			}
			m(k, k) = Admittance_s;
			cout << m(k, k) << endl;
			if (j!=k) {
				m(j, k) = -Admittance_d;
				m(k, j) = -Admittance_d;
			}
		}
	}*/
	for (int k = 0; k < num; k++)
	{
		complex <double> Admittance_s(0, 0);
		for (int i = 0; i < branchcount; i++)
		{
			complex <double> admit = Branchlist[i]->getz();
			if ((Branchlist[i]->getnode1()->getrank() == arr[k]) || (Branchlist[i]->getnode2()->getrank() == arr[k]))
			{
				complex <double> r;
				r = Branchlist[i]->get_Admittance();
				Admittance_s += r;
			}
		}
		m(k, k) = Admittance_s;
	}
	complex <double>* Glist = new complex <double>[num];
	for (int v = 0; v < num; v++)
	{
		for (int k = v+1; k < num; k++)
		{
			complex <double> Admittance_d(0, 0);
			for (int l = 0; l < branchcount; l++)
			{
				if (((Branchlist[l]->getnode1()->getrank() == arr[v]) || (Branchlist[l]->getnode2()->getrank() == arr[v])) 
					&&((Branchlist[l]->getnode1()->getrank() == arr[k]) || (Branchlist[l]->getnode2()->getrank() == arr[k])))
				{
					Admittance_d += Branchlist[l]->get_Admittance();
				}
			}
			m(v, k) = -Admittance_d;
			m(k, v) = -Admittance_d;	
		}
	}
	for (int j = 0; j < num; j++)
	{
		for (int i = 0; i < comcount; i++)
		{
			Component* comp = complist[i];
			Isrc* current = dynamic_cast<Isrc*>(comp);
			if (current != NULL)
			{
				if (current->get_node1() == arr[j])
				{
					complex <double> curr = current->getcomplex();
					I(j, 0) += curr;

				}
				if (current->get_node2() == arr[j])
				{
					complex <double> curr = -current->getcomplex();
					I(j, 0) += curr;
				}
			}
		}
	}
	
	//for (int i=0;i<)
	cout << "................" << endl;
	cout << m<< endl;
	cout << "................" << endl;
	cout << I << endl;
	cout << "................" << endl;
	///////////////////////////sabry el gamed///////////////////////

	// Dependent Source CCCs
	int cccs_node1;
	double sign = 1;
	int index;
	if (cccs != NULL)
	{
		for (int i = 1; i < nonsimple_nodecount; i++)
		{
			if (cccs->get_node1() == nonsimple_Nodelist[i]->getrank())
			{
				index = i;
				sign *= 1;
			}
			else if (cccs->get_node2() == nonsimple_Nodelist[i]->getrank())
			{
				index = i;
				sign *= -1;
			}

		}
		for (int i = 1; i < nonsimple_nodecount; i++)
		{
			if(nonsimple_Nodelist[i]->getrank() == cccs->get_Noded1())
			{
				m(index - 1, i - 1) += sign * cccs->get_Admittance() * cccs->get_cofficient();
			}
			else if (nonsimple_Nodelist[i]->getrank() == cccs->get_Noded2())
			{
				m(index - 1, i - 1) -= sign * cccs->get_Admittance() * cccs->get_cofficient();
			}
		}
	}
	MatrixXcd V = CalculateVoltNode(m, I); // function calculate volt of each node
	for (int i = 1; i < nonsimple_nodecount; i++)
	{
		nonsimple_Nodelist[i]->setvoltage(V(i-1, 0));
	}  
	for (int i = 0; i < branchcount; i++)
	{
		Branchlist[i]->Calculatecurrent();
	}
	cout << V;
	return 0;
}

///////////////////////////sabry el gamed///////////////////////
MatrixXcd CalculateVoltNode(MatrixXcd G, MatrixXcd I)
{
	MatrixXcd V;
	V = G.inverse() * I;
	return V;
}
//////////////////////////////Ahmed hany
int searcher(string x)
{
	int c;
	Component* comp = NULL;
	for (int i = 0; i < comcount; i++)
	{
		string name = complist[i]->get_name();
		if ((name.compare(x) == 0))
		{
			comp = complist[i];
			c = i;
		}
	}
	return c;
}