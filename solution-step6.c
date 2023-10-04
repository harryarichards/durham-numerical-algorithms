// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2018.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <omp.h>

double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0001;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
    NumberOfBodies = (argc-2) / 7;

    x    = new double*[NumberOfBodies];
    v    = new double*[NumberOfBodies];
    mass = new double [NumberOfBodies];

    int readArgument = 1;

    tPlotDelta  = std::stof(argv[readArgument]); readArgument++;
    tFinal      = std::stof(argv[readArgument]); readArgument++;

    for (int i=0; i<NumberOfBodies; i++) {
        x[i] = new double[3];
        v[i] = new double[3];

        x[i][0] = std::stof(argv[readArgument]); readArgument++;
        x[i][1] = std::stof(argv[readArgument]); readArgument++;
        x[i][2] = std::stof(argv[readArgument]); readArgument++;

        v[i][0] = std::stof(argv[readArgument]); readArgument++;
        v[i][1] = std::stof(argv[readArgument]); readArgument++;
        v[i][2] = std::stof(argv[readArgument]); readArgument++;

        mass[i] = std::stof(argv[readArgument]); readArgument++;

        if (mass[i]<=0.0 ) {
            std::cerr << "invalid mass for body " << i << std::endl;
            exit(-2);
        }
    }

    std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

    if (tPlotDelta<=0.0) {
        std::cout << "plotting switched off" << std::endl;
        tPlot = tFinal + 1.0;
    }
    else {
        std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
        tPlot = 0.0;
    }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
    videoFile.open( "result.pvd" );
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
              << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
    static int counter = -1;
    counter++;
    std::stringstream filename;
    filename << "result-" << counter <<  ".vtp";
    std::ofstream out( filename.str().c_str() );
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    for (int i=0; i<NumberOfBodies; i++) {
        out << x[i][0]
            << " "
            << x[i][1]
            << " "
            << x[i][2]
            << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>"  << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
    //We need to parallelise this.
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();

    // Declares and initialises an array that stores the forces acting on each particle.
    double force[NumberOfBodies][3];
    #pragma omp parallel
    {
        #pragma omp for
        for (int current = 0; current < NumberOfBodies; current++) {
            force[current][0] = 0.0;
            force[current][1] = 0.0;
            force[current][2] = 0.0;
        }

    }

    /* For each particle we calculate the force between it and each of the other particles.
     * We do this using a nested for loop and both for loops are parallelised.
     */
    #pragma omp parallel
    {
        #pragma omp for
        for (int current = 0; current < NumberOfBodies; current++) {
            #pragma omp parallel
            {
                #pragma omp for
                for (int other = current + 1; other < NumberOfBodies; other++) {
                    // Calculate the distance along each axis between the current particle and the other particle.
                    const double zero_distance = x[other][0] - x[current][0];
                    const double one_distance = x[other][1] - x[current][1];
                    const double two_distance = x[other][2] - x[current][2];

                    // Calculates the distance between the two particles.
                    const double distance = sqrt(
                            (zero_distance * zero_distance) +
                            (one_distance * one_distance) +
                            (two_distance * two_distance)
                    );

                    /* Calculates the product of the mass of the two particles divided by the distance between them squared.
                     * This is a scalar value used repeatedly in the force calculations below.
                     */
                    double mass_product_over_distance_cubed = (
                            mass[current] * mass[other] / (distance * distance * distance)
                    );

                    // Calculates the magnitude of the force between the two particles in each direction.
                    const double symmetric_force_zero = zero_distance * mass_product_over_distance_cubed;
                    const double symmetric_force_one = one_distance * mass_product_over_distance_cubed;
                    const double symmetric_force_two = two_distance * mass_product_over_distance_cubed;


                    // Use critical to avoid race conditions.
                    #pragma omp critical
                    {
                        // Sums the magnitude of the force for the current particle (in each direction).
                        force[current][0] += symmetric_force_zero;
                        force[current][1] += symmetric_force_one;
                        force[current][2] += symmetric_force_two;
                        force[other][0] -= symmetric_force_zero;
                        force[other][1] -= symmetric_force_one;
                        force[other][2] -= symmetric_force_two;
                    }

                    minDx = std::min(minDx, distance);

                }
            }
        }
    }

    /* Here we update the position and velocity of each particle using the forces calculated above.
     * We do this in parallel.
     * If we were to update the position of the particles in parallel in the above loop it would lead to inaccuracies in the calculation of velocity.
     */
    #pragma omp parallel
    {
        #pragma omp for
        for (int current = 0; current < NumberOfBodies; current++) {
            // Calculates the new coordinates of the current particle.
            x[current][0] += timeStepSize * v[current][0];
            x[current][1] += timeStepSize * v[current][1];
            x[current][2] += timeStepSize * v[current][2];
            const double time_step_over_mass = timeStepSize / mass[current];
            v[current][0] += time_step_over_mass * force[current][0];
            v[current][1] += time_step_over_mass * force[current][1];
            v[current][2] += time_step_over_mass * force[current][2];
            maxV = std::max(maxV, std::sqrt(
                    v[current][0] * v[current][0] + v[current][1] * v[current][1] + v[current][2] * v[current][2]));
        }
    }

    t += timeStepSize;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
    if (argc==1) {
        std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time objects" << std::endl
                  << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
                  << "  final-time      simulated time (greater 0)" << std::endl
                  << std::endl
                  << "Examples:" << std::endl
                  << "0.01  100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
                  << "0.01  100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
                  << "0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
                  << std::endl
                  << "In this naive code, only the first body moves" << std::endl;

        return -1;
    }
    else if ( (argc-3)%7!=0 ) {
        std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
        return -2;
    }

    setUp(argc,argv);

    //openParaviewVideoFile();

    int snapshotCounter = 0;
    if (t > tPlot) {
        //printParaviewSnapshot();
        std::cout << "plotted initial setup" << std::endl;
        tPlot = tPlotDelta;
    }

    int timeStepCounter = 0;
    while (t<=tFinal) {
        updateBody();
        timeStepCounter++;
        if (t >= tPlot) {
            //printParaviewSnapshot();
            std::cout << "plot next snapshot"
                      << ",\t time step=" << timeStepCounter
                      << ",\t t="         << t
                      << ",\t dt="        << timeStepSize
                      << ",\t v_max="     << maxV
                      << ",\t dx_min="    << minDx
                      << std::endl;

            tPlot += tPlotDelta;
        }
    }

    //closeParaviewVideoFile();

    return 0;
}
