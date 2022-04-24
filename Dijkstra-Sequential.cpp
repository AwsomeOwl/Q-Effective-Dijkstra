// Manatin Pavel, South Ural State University 
// 2022
#include <vector>
#include <omp.h>
#include <time.h>
#include <iostream>
using namespace std;
vector<vector<unsigned short> > RandomMatrix(int size)
{
	srand(time(NULL));
	vector<vector<unsigned short> > matrix(size, vector<unsigned short>(size, 0));
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i != j)
				matrix[i][j] = rand() % 100;
	return matrix;
}
int main()
{
	vector < double> times(0,0);
	for (int k = 0; k < 10; k++)
	{
		for (int t = 1000; t < 2000; t += 1000)
		{
			vector<vector<unsigned short> > a = RandomMatrix(t);
			double e_time, s_time, f_time;
			int i;
			std::vector<unsigned short> r (1,0);
			std::vector<bool> p (1,0);
			std::vector<bool> p_flag (1,0);
			omp_set_dynamic(0);
			omp_set_num_threads(1);
			s_time = omp_get_wtime();
			for (i = 1; i < a.size(); i++)
			{
				r.push_back(USHRT_MAX);
				p.push_back(0);
				p_flag.push_back(1);
			}
			i = 0;
			while (i < a.size())
			{
				if (a[0][i] != 0)
					r[i] = a[0][i];
				i += 1;
			}

			while (p != p_flag)
			{
				int temp = 0;
				unsigned short min = USHRT_MAX;
				i = 1;
				while (i < a.size())
				{
					if (r[i] <= min and p[i] == 0)
					{
						min = r[i];
						temp = i;
					}
					i += 1;
				}
				p[temp] = 1;
				i = 1;
				while (i < a.size() and temp < a.size())
				{
					if (p[temp] == 1 and i != temp)
						if (a[temp][i] != 0)
							if ((r[temp] + a[temp][i]) < r[i])
								r[i] = r[temp] + a[temp][i];
					i += 1;
				}
			}
			f_time = omp_get_wtime();
			e_time = f_time - s_time;
			times.push_back(e_time);
		}
		for (int v = 0; v < times.size(); v++)
			std::cout << times[v] << ",";
		std::cout << endl;
		std::cout << endl;
		times.clear();
	}
}
