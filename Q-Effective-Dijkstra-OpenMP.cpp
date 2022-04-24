// Pavel Manatin, South Ural State University 
// 2022
#include <vector>
#include <omp.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
std::vector<std::vector<unsigned short>> input(std::string filename, int n)
{
	std::vector<std::vector<unsigned short>> a;
	a.resize(n);
	std::string line;
	std::ifstream in(filename);
	if (in.is_open())
	{
		unsigned short x;
		while (!in.eof())
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					in >> x;
					a[i].push_back(x);
				}
			break;
		}
	}
	in.close();
	return a;

}
void output(std::vector<unsigned short> r, double start, double finish)
{
#pragma omp parallel
#pragma omp single
	std::cout << "Number of threads:" << omp_get_num_threads() << std::endl;
	std::cout << "Shortest paths from first node:" << std::endl;
	for (int v = 0; v < r.size() - 1; v++)
		std::cout << r[v] << ",";
	std::cout << r.back() << std::endl;
	std::cout << "Time: " << finish - start << "." << std::endl;
}

std::vector<std::vector<unsigned short> > RandomMatrix(int size)
{
	std::vector<std::vector<unsigned short> > matrix(size, std::vector<unsigned short>(size, 0));
	srand(int(time(NULL)));
#pragma omp parallel  for 
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i != j)
				matrix[i][j] = rand() % 100;
	return matrix;
}
void dijkstra(std::vector<std::vector<unsigned short> > a)
{
	double e_time, s_time, f_time;
	std::vector<unsigned short> r(a.size(), 0);
	std::vector<bool> p(a.size(), 0);
	s_time = omp_get_wtime();
#pragma omp parallel for 
	for (int i = 1; i < a.size(); i++)
	{
		r[i] = 101;
		p[i] = 0;
	}
#pragma omp parallel for 
	for (int i = 0; i < a.size(); i++)
		if (a[0][i] != 0)
			r[i] = a[0][i];
	for (int t = 1; t < a.size(); t++)
	{
		int temp = 0;
		unsigned short min = 101;
		for (int i = 1; i < a.size(); i++)
			if (r[i] <= min && p[i] == 0)
			{
				min = r[i];
				temp = i;
			}
		p[temp] = 1;
#pragma omp parallel for
		for (int i = 1; i < a.size(); i++)
			if (p[temp] == 1 && i != temp)
				if (a[temp][i] != 0)
					if ((r[temp] + a[temp][i]) < r[i])
						r[i] = r[temp] + a[temp][i];
	}

	f_time = omp_get_wtime();
	e_time = f_time - s_time;
	output(r, s_time, f_time);
}
int main(int argc, char* argv[])
{
	//unsigned short tries = 10;
	//std::vector<double> times(tries, 0);
	omp_set_dynamic(0);
	omp_set_num_threads(std::stoi(argv[1]));
	std::vector<std::vector<unsigned short> > a = input((std::string)argv[2], std::stoi(argv[3]));
	dijkstra(a);
	return 0;
}
