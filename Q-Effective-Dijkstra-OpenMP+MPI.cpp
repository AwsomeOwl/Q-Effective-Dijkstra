// Manatin Pavel, South Ural State University 
// 2022
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "mpi.h"
#include "omp.h"
std::vector<short> input(std::string filename)
{
	std::vector<short> a;
	std::string line;
	std::ifstream in(filename);
	if (in.is_open())
	{
		short x;
		while (!in.eof())
		{
			in >> x;
			a.push_back(x);
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
std::vector<short> RandomMatrix(int size)
{
	std::vector<short>  matrix(size*size, 0);
	srand(int(time(NULL)));
#pragma omp parallel  for 
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i != j)
				matrix[i * size + j] = (rand()+omp_get_thread_num()) % 100;
	return matrix;
}
/*vector<vector<double>> e = {{0.0, 10.0, 0.0,0.0,0.0,0.0,3.0,6.0,12.0},
	{ 10.0, 0.0, 18.0, 0.0, 0.0, 0.0, 2.0, 0.0, 13.0 },
	{ 0.0, 18.0, 0.0, 25.0, 0.0, 20.0, 0.0, 0.0, 7.0 },
	{ 0.0, 0.0, 25.0, 0.0, 5.0, 16.0, 4.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0, 5.0, 0.0, 10.0, 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 20.0, 0.0, 10.0, 0.0, 14.0, 15.0, 9.0 },
	{ 0.0, 2.0, 0.0, 4.0, 0.0, 14.0, 0.0, 0.0, 24.0 },
	{ 6.0, 0.0, 0.0, 0.0, 23.0, 15.0, 0.0, 0.0, 5.0 },
	{ 12.0, 13.0, 0.0, 0.0, 0.0, 9.0, 24.0, 5.0, 0.0 } };
/* { {0, 1, 0, 1}, {0,0,1,1}, {0,1,0,0}, {1,0,1,0}};*/
/*vector <double> eTransformed(pow(e.size(), 2), 0);*/
void dijkstra(int world_size, int world_rank, std::vector<short> &eTransformed)
{
	int size;
	size = pow(eTransformed.size(), 0.5);
	//size = size - size % world_size;
	/* { {0, 1, 0, 1}, {0,0,1,1}, {0,1,0,0}, {1,0,1,0}};*/
#pragma omp parallel for
	for (int i = 0; i < size; i++)
		eTransformed[i*size+i] = -1;
#pragma omp parallel for
	for (int i = 0; i < size*size; i++)
		if (eTransformed[i]==0)
			eTransformed[i] = 999;
	MPI_Barrier(MPI_COMM_WORLD);
	double start_time = MPI_Wtime();
	double finish_time;
	//MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	std::vector<unsigned short> distances_prev(size, 0);
	std::vector<unsigned short> distances_curr(size, 0);
	bool FinishInitialized = false;
#pragma omp parallel  for 
	for (int i = 1; i < size; i++) 
		distances_curr[i] = eTransformed[i]; //заполним начальную матрицу расстояний
	//vector<double>distances_prev_splitted(e.size()*world_size, 0);
	//vector<double>distances_curr_splitted(e.size() * world_size, INFINITY);
	std::vector<short>distancesBuffer(size*size / world_size, 0);
	std::vector<unsigned short>distancesReducedBuffer(size, 999);
	while (true)
	{
		/* {
			cout << "current " << world_rank << endl;
			for (int i = 0; i < e[0].size(); i++)
				cout << distances_curr[i] << " ";
			cout << endl;
			cout << "prev " << world_rank << endl;
			for (int i = 0; i < e[0].size(); i++)
				cout << distances_prev[i] << " ";
			cout << endl;
		}*/
		MPI_Barrier(MPI_COMM_WORLD);
		if (world_rank == 0)
		{
			if (distances_curr == distances_prev)
			{
				finish_time = MPI_Wtime();
				FinishInitialized = true;
			}
			else
				distances_prev = distances_curr;
		}
		MPI_Bcast(&FinishInitialized, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (FinishInitialized)
			break;
		MPI_Bcast(distances_curr.data(), size, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
			//MPI_Bcast(distances_prev.data(), distances_prev.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter(eTransformed.data(), size*size/world_size, MPI_SHORT, distancesBuffer.data(), size*size / world_size, MPI_SHORT, 0, MPI_COMM_WORLD);
		for (int k = 0; k < size / world_size; k++) //проходим по строкам
		{
			short rootDistance = 0;
			int rootNumber = -1;
#pragma omp parallel for
			for (int i = 0; i < size; i++)//найти корневую вершину
				if (distancesBuffer[k * size + i] == -1)
				{
					distancesBuffer[k * size + i] = distances_curr[i];
					rootDistance = distances_curr[i];
					rootNumber = i;
				}
#pragma omp parallel for
				for (int i = 0; i < size; i++)
					distancesBuffer[k * size + i] += rootDistance; // преобразуем буфер
#pragma omp parallel for
			for (int i = 0; i < size; i++)//проходим по строке и редуцируем в буфер размером, равным кол-вом вершин
				distancesReducedBuffer[i] = ((distancesBuffer[k * size + i] < distancesReducedBuffer[i]))
				? (distancesBuffer[k * size + i])
				: distancesReducedBuffer[i];
		}
		MPI_Reduce(distancesReducedBuffer.data(), distances_curr.data(), size, MPI_UNSIGNED_SHORT, MPI_MIN,0, MPI_COMM_WORLD);
		//cout << "reduced " << world_rank << endl;
		//for (int i = 0; i < e[0].size(); i++)
			//cout << distances_curr[i] << " ";
		//cout << endl;
		//distancesBuffer.clear();
		//distancesBuffer.resize(eTransformed.size()/world_size);
	}
	if (world_rank == 0)
	{
		output(distances_curr, start_time, finish_time);
	}
}
int main(int argc, char** argv)
{
	int tries = 10;
	int world_size = 0;
	int world_rank = 0;
	int size;
	std::vector<short> eTransformed = input(argv[2]);
	omp_set_dynamic(0);
	omp_set_num_threads(std::stoi(argv[1]));
	size = pow(eTransformed.size(), 0.5);
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	dijkstra(world_size, world_rank, eTransformed);
	MPI_Finalize();
		
}
