#include <iostream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <vector>


#include "plif.h"

int main()
{
	using namespace std;
	using namespace plif;
	typedef piecewise_linear_function PLF;
	typedef piecewise_linear_interval_function PLIF;
	
	interval A(-5, 5);
	interval B(6, 10);
	interval C(0, 7);
	interval D(3, 7);
	interval E(8, 11);
	vector<interval> list;
	list.push_back(A);
	list.push_back(B);
	list.push_back(C);
	list.push_back(D);
	list.push_back(E);
	vector<interval> partitioned;
	partitioned = partition(list);
	for(int i=0;i<partitioned.size();i++)
	{
		cout << "[" << partitioned[i].lower() << ", " << partitioned[i].upper() << "]" << endl;
	}
	
		cout << endl << endl;
	
	PLF f, g;
	f.set_domain(interval(-5, 5));
	g.set_domain(interval(-3, 3));
	f.set_left_bounded();
	g.set_left_bounded();
	f.set_right_bounded();
	g.set_right_bounded();
	f.set_left_closed(true);
	f.set_right_closed(true);
	g.set_left_closed(true);
	g.set_right_closed(true);
	f.insert_point(-5, 1, 1, 1);
	f.insert_point(0, 2, 2, 2);
	f.insert_point(5, 1, 1, 1);
	g.insert_point(-3, 2, 2, 2);
	g.insert_point(0, 1, 1, 1);
	g.insert_point(3, 2, 2, 2);
	PLF h;

	h= f+g;
	h.display();
	
	cout << endl << endl; 
	
	PLF l, u;
	l.set_domain(interval(1, 4));
	u.set_domain(interval(1, 4));
	l.set_left_closed(true);
	u.set_left_closed(true);
	l.set_right_closed(true);
	u.set_right_closed(true);
	l.set_left_bounded();
	l.set_right_bounded();
	u.set_left_bounded();
	u.set_right_bounded();
	l.insert_point(1, 2, 2, 2);
	l.insert_point(2, 1, 1, 1);
	l.insert_point(3, 1, 1, 1);
	l.insert_point(4, 2, 2, 2);
	u.insert_point(1, 3, 3, 3);
	u.insert_point(2.5, 1.2, 1.3, 1.4);
	u.insert_point(4, 3, 3, 3);	
	shortest_path(l, u, breakpoint(1, 2, 2, 2), breakpoint(4, 2, 2, 2)).display();
	cout << endl << dist(l, u) << endl;
	
	PLIF p(l, u);
	p.display();
	cout << endl;
	std::vector<precision_type> points;
	points.push_back(1.5);
	points.push_back(2.5);
	points.push_back(3.5);


	std::vector<PLIF> dissected = dissect(p, points);
	int i;
	for(i=0;i<dissected.size();i++)
	{
		dissected[i].display();
		cout << endl;
	}
	cout << endl;
//	std::vector<PLIF> list_of_plif = split(p);
//	cout << list_of_plif.size() << endl;

	std::vector<interval> subdomains = above_threshold_subdomains(l, 1.5);
//	cout << subdomains.size() << endl;
	
	for(i=0;i<subdomains.size();i++)
	{
		cout << "[" << subdomains[i].lower() << "," << subdomains[i].upper() << "]" << endl;
	}
}
