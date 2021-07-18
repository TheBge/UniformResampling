/***********************************************************************************************
author: Bge

Description:
		A resample program for point cloud file.
		The point file(.txt) is from a slice program.
input:
		X coordinate      Y coordinate
ouput:
		X coordinate      Y coordinate      Z coordinate(layer height)

************************************************************************************************
History				Name				Action
18-Jul-2021			Bge					Initial
************************************************************************************************/

#include <string>
#include <array>
#include <vector>

#include <cmath>
#include <math.h>

#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/types.hpp>

using namespace std;
using namespace cv;

template<typename T, typename V>
void CoordinateSplit(const vector<T>& pl, 
	vector<V>& contourx, 
	vector<V>& contoury,
	vector<V>& contourz
)
{
	int n = 0;
	contourx.resize( pl.size() / 3 );
	contoury.resize( pl.size() / 3 );
	contourz.resize( pl.size() / 3 );

	for (int j = 0; j < pl.size();)
	{
		contourx[n] = (T)(pl[j]);
		contoury[n] = (T)(pl[j+1]);
		contourz[n] = (T)(pl[j + 2]);

		j = j + 3;
		n++;
	}
}

template<typename T, typename V>
void PolyLineSplit(const vector<Point_<T> >& pl, vector<V>& contourx, vector<V>& contoury)
{
	contourx.resize(pl.size());
	contoury.resize(pl.size());
	for (int j = 0; j < pl.size(); j++)
	{
		contourx[j] = (V)(pl[j].x);
		contoury[j] = (V)(pl[j].y);
	}
}

template<typename T, typename V>
void PolyLineMerge(vector<Point_<T> >& pl, const vector<V>& contourx, const vector<V>& contoury) 
{
	assert(contourx.size() == contoury.size());
	pl.resize(contourx.size());

	for (int j = 0; j < contourx.size(); j++) 
	{
		pl[j].x = (V)(contourx[j]);
		pl[j].y = (V)(contoury[j]);
	}
}

template<typename T>
double pl_arcLength(vector <Point_<T>>xycoordinates, bool isOpen) 
{
	double length = 0;
	for (int i = 0; i < xycoordinates.size() - 1; i++) 
	{
		length = length + norm(xycoordinates[i] - xycoordinates[i + 1]);
	}

	if (isOpen) 
	{
		length += norm(xycoordinates.front() - xycoordinates.back());
	}
	return length;
}

// perform the linear interpolation for inputed open contours.
void ResampleCurve(const vector<double>& curvex, 
	const vector<double>& curvey,
	const vector<double>& curvez,
	vector<double>& resampleX, 
	vector<double>& resampleY,
	int N,
	bool isOpen,
	vector<double>& sumXCoordinates,
	vector<double>& sumYCoordinates,
	vector<double>& sumZCoordinates
) 
{
	assert(curvex.size() > 0 && curvey.size() > 0 && curvex.size() == curvey.size());

	vector<Point2d> resamplepl(N);
	resamplepl[0].x = curvex[0];
	resamplepl[0].y = curvey[0];
	vector<Point2d> pl;
	PolyLineMerge(pl, curvex, curvey);
	double pl_length = pl_arcLength(pl, isOpen);
	double resample_size = pl_length / (double)N;
	int curr = 0;
	double dist = 0.0;

	for (int i = 1; i < N; ) 
	{
		assert(curr < pl.size() - 1);
		double last_dist = norm(pl[curr] - pl[curr + 1]);
		dist += last_dist;
		// cout << curr << " and " << curr+1 << "\t\t" << last_dist << " ("<<dist<<")"<<endl;
		if (dist >= resample_size) 
		{
			//put a point on line
			double _d = last_dist - (dist - resample_size);
			Point2d cp(pl[curr].x, pl[curr].y), cp1(pl[curr + 1].x, pl[curr + 1].y);
			Point2d dirv = cp1 - cp; dirv = dirv * (1.0 / norm(dirv));
			// cout << "point " << i << " between " << curr << " and " << curr+1 << " remaining " << dist << endl;
			assert(i < resamplepl.size());
			resamplepl[i] = cp + dirv * _d;
			i++;
			dist = last_dist - _d; //remaining dist
			//if remaining dist to next point needs more sampling... (within some epsilon)
			while (dist - resample_size > 1e-10) 
			{
				// cout << "point " << i << " between " << curr << " and " << curr+1 << " remaining " << dist << endl;
				assert(i < resamplepl.size());
				resamplepl[i] = resamplepl[i - 1] + dirv * resample_size;
				dist -= resample_size;
				i++;
			}
		}
		curr++;
	}

	PolyLineSplit(resamplepl, resampleX, resampleY);

	for (auto sampleXPoint : resampleX)
	{
		sumXCoordinates.push_back(sampleXPoint);
	}

	for (auto sampleYPoint : resampleY)
	{
		sumYCoordinates.push_back(sampleYPoint);
	}

	double currentLayer = curvez[0];
	for (int i = 0; i < resampleX.size(); ++i)
	{
		sumZCoordinates.push_back(currentLayer);
	}
}

int main()
{
	Eigen::Vector2d row_vec;
	ifstream in("D:/workdir/TheBge/UniformResampling/result/Path.txt");
	string line;
	vector<double> coordinates;
	vector<double> x_coordinates;
	vector<double> y_coordinates;
	vector<double> z_coordinates;

	vector<double> resamplex;
	vector<double> resampley;

	// input
	for (string s; getline(in, s); )
	{
		istringstream sin(s);     
		for (double ia; sin >> ia; )      
		{
			coordinates.push_back(ia);      
		}
	}

	CoordinateSplit(coordinates, x_coordinates, y_coordinates, z_coordinates);

	vector<double> xCoordinates;
	vector<double> yCoordinates;
	vector<double> zCoordinates;
	vector<double>::iterator iter = z_coordinates.begin();

	int i = 0;
	vector<double> sumXCoordinates{};
	vector<double> sumYCoordinates{};
	vector<double> sumZCoordinates{};

	for (vector<double>::iterator iter= z_coordinates.begin(); iter != z_coordinates.end(); ++iter)
	{
		if ((iter+1) != z_coordinates.end() && z_coordinates[i] == z_coordinates[i+1])
		{
			xCoordinates.push_back(x_coordinates[i]);
			yCoordinates.push_back(y_coordinates[i]);
			zCoordinates.push_back(z_coordinates[i]);
			i++;
			continue;
		}
		else
		{
			xCoordinates.push_back(x_coordinates[i]);
			yCoordinates.push_back(y_coordinates[i]);
			zCoordinates.push_back(z_coordinates[i]);
			i++;
			ResampleCurve(xCoordinates, yCoordinates, zCoordinates, 
				resamplex, resampley, 500, false, sumXCoordinates, sumYCoordinates, sumZCoordinates);
			xCoordinates.clear();
			yCoordinates.clear();
			zCoordinates.clear();
		}
	}

	fstream output;
	output.open("D:/workdir/TheBge/UniformResampling/result/ResamplePath.txt", fstream::out);
	for (int i = 0; i < sumXCoordinates.size(); i++)
	{
		output << sumXCoordinates[i] << "     " << sumYCoordinates[i] << "     " << sumZCoordinates[i] << endl;
	}
	output.close();

	return 0;
}
