#include <vector>
#include <tuple>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

//-----------------------------------------------------------

void save3D(std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> data, std::string filename)
{

    std::vector<double> xX = std::get<0>(data);
    std::vector<double> yY = std::get<1>(data);
    std::vector<double> zZ = std::get<2>(data);

    std::ofstream ofs(filename, std::ofstream::out);

    double fx, fy, fz;
    for (int i = 0; i <= xX.size() - 1; i++)
    {
        fx = xX.at(i);
        fy = yY.at(i);
        fz = zZ.at(i);

        ofs << std::fixed << fx << " ";
        ofs << std::fixed << fy << " ";
        ofs << std::fixed << fz;
        ofs << std::endl;
    }

    ofs.flush();
    ofs.close();
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> read3D(std::string filename)
{
    // The input file is expected to have three floating point values per line seperated by white space
    std::vector<double> pathX;
    std::vector<double> pathY;
    std::vector<double> pathZ;

    std::fstream in(filename);
    for (std::string line; std::getline(in, line);) // read stream line by line
    {
        std::stringstream ss(line); // make a stream for the line
        double x, y, z;
        ss >> x >> y >> z; // read the whitespace-separated floats
        pathX.push_back(x);
        pathY.push_back(y);
        pathZ.push_back(z);
    }
    return std::make_tuple(pathX, pathY, pathZ);
}
//-----------------------------------------------------------

double computeBinominal(int n, int k)
{
    if (k > n)
        return 0;
    double value = 1.0;
    for (int i = 1; i <= k; i++)
    {
        value *= n--;
        value /= i;
    }
    return value;
}

//-----------------------------------------------------------

std::vector<double> bezierCurve(std::vector<double> lineData, unsigned int outputSize)
{
    std::vector<double> bCurve;
    int n = lineData.size() - 1;
    for (int t_index = 0; t_index <= outputSize-1; t_index++)
    {
        double t = double(t_index) / double(outputSize-1);
        double bCurvet{0};
        for (int i = 0; i <= n; ++i)
            bCurvet += computeBinominal(n, i) * std::pow((1 - t), (n - i)) * std::pow(t, i) * lineData[i];
        bCurve.push_back(bCurvet);
    }
    return bCurve;
    // Note the curve shape can be modified by adding rational weights to the numerator and denominator as per
    // https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Rational_B%C3%A9zier_curves
}

//-----------------------------------------------------------

std::vector<double> bezierCurve(std::vector<double> lineData, unsigned int outputSize, std::vector<double> localWeight, double globalWeight = 1.0)
{
    /*
    localWeight: 
                    This vector contains the rational weighting values for the 
                    bezier curve. It is the same length as lineData and will 
                    emphasize the importance of point i in proportion to the 
                    value of its weight compared to the sum of all the weights.
    globalWeight:
                    When globalWeight = 1, this function produces a rationally 
                    weighted bezier curve.
                    When globalWeight > 1, the curve will be "pushed" toward
                    the original data points creating a tighter but less smooth
                    fit.
                    When globalWeight < 1, the curve behaviour is undocumented.
    A useful description of Rational Weighted beziers is here:
    https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Rational_B%C3%A9zier_curves

    */                  
    std::vector<double> bCurve;
    int n = lineData.size() - 1;
    for (int t_index = 0; t_index <= outputSize-1; t_index++)
    {
        double t = double(t_index) / double(outputSize-1);
        double bCurvetNum{0};
        double bCurvetDen{0};
        for (int i = 0; i <= n; ++i)
        {
            double binomial = computeBinominal(n, i);
            double globalWeightFactor = std::pow(globalWeight, -(std::abs(t*n - double(i))) );
            bCurvetNum += binomial * std::pow((1 - t), (n - i)) * std::pow(t, i) * lineData[i] * localWeight[i] * globalWeightFactor;
            bCurvetDen += binomial * std::pow((1 - t), (n - i)) * std::pow(t, i) * localWeight[i] * globalWeightFactor;
        }
        bCurve.push_back(bCurvetNum/bCurvetDen);
    }
    return bCurve;
}

//-----------------------------------------------------------

int main()
{
    // Load the point data from a file
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> inData = read3D("UEPath.dat");

    // The bezier operation happens on one axis at a time and does not depend on the value of the
    // other axes. For this reason it is useful to store the x, y, and z values each in their own array.
    std::vector<double> xX = std::get<0>(inData); // Put all the X values into their own array
    std::vector<double> yY = std::get<1>(inData); // Put all the Y values into their own array
    std::vector<double> zZ = std::get<2>(inData); // Put all the Z values into their own array

    std::vector<double> weight;
    for (int i = 0; i < xX.size(); ++i)
        {
            weight.push_back(1.0);
        }

    // This is the important stuff
    const unsigned int outputSize = xX.size()*2;           // The number of points in the output curve. This could be a useful parameter to alter at the blueprint level or maybe compute a reasonable value based on how long the original path is.
    std::vector<double> xXb = bezierCurve(xX, outputSize); // Compute all the bezier X values
    std::vector<double> yYb = bezierCurve(yY, outputSize); // Compute all the bezier Y values
    std::vector<double> zZb = bezierCurve(zZ, outputSize); // Compute all the bezier Z values
    // The important stuff is done

    // Store the results to a file
    save3D(std::make_tuple(xXb, yYb, zZb), "ue_bz.dat");

    double gw = 1.5;
    std::vector<double> xXbw = bezierCurve(xX, outputSize, weight, gw); // Compute all the bezier X values
    std::vector<double> yYbw = bezierCurve(yY, outputSize, weight, gw); // Compute all the bezier Y values
    std::vector<double> zZbw = bezierCurve(zZ, outputSize, weight, gw); // Compute all the bezier Z values
                                                           // The important stuff is done

    // Store the results to a file
    save3D(std::make_tuple(xXbw, yYbw, zZbw), "ue_wbz.dat");
}