/*
 * sample_plotter.h
 *
 *  Created on: Aug 9, 2011
 *      Author: ray
 */

#ifndef SAMPLE_PLOTTER_H_
#define SAMPLE_PLOTTER_H_

#include <fstream>

template<typename scalar_type>
class sample_plotter {
public:
	struct sample_point {
		scalar_type x;
		scalar_type f_x;
	};

	sample_plotter(){};
	~sample_plotter(){};

	/**
	 * Prints the sample point to the console
	 */
	void print_samples() const {
		std::cout << "Sampling points in the sequence of request:" << std::endl;

		for(typename std::list<sample_point>::const_iterator it = my_samples_list.begin(); it!=my_samples_list.end(); it++){
			std::cout << "[" <<(*it).x << ", " << (*it).f_x << "]" << std::endl;
		}

	};
	void add_sample(scalar_type x, scalar_type f_x){
		sample_point p;
		p.x = x;
		p.f_x = f_x;
		my_samples_set.insert(x);
		my_samples_list.push_back(p);
	};

	/**
	 * Returns the number of distinct samples requested
	 *
	 * @return size of the samples set.
	 */
	size_t get_size() const {
		return my_samples_set.size();
	}
	/**
	 * Plots the samples as blue points
	 */
	void plot_graph() const {
		std::ofstream temp_file;
		temp_file.open("/tmp/lb_samples.txt");
		for(typename std::list<sample_point>::const_iterator it = my_samples_list.begin(); it!=my_samples_list.end();it++){
			temp_file << (*it).x << " " << (*it).f_x << std::endl;
		}
		temp_file.close();
		int res = std::system("graph -TX -C -S 3 -m -3 /tmp/lb_samples.txt");
	}
private:
	std::list<sample_point> my_samples_list; // List preserve the order of sample request
	std::set<scalar_type> my_samples_set; // Keeps track of distinct sample point without function values.

};

#endif /* SAMPLE_PLOTTER_H_ */
