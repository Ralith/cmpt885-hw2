/** Implementation of an ant colony optimization for the travelling salesman problem.  
    Drawn HEAVILY from Jones' "Artificial Intelligence:  A Systems Approach", pp. 423-430.

    @author Adrian Kwok
*/
#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include "omp.h"

using namespace std;

int iterations = 1000;
int numAnts = 1000;


//ant colony parameters:
double alpha = 0.7; //favour pheromone level over distance (pp 430)
double beta = 1.0;  //favour distance over pheromone level (pp 430)
double p = 0.8;     //intensification/evaporation value.  means 10% evaporation per tour
double Q = 100;     //constant affecting how much base pheromone to put into an edge at each iteration
double BASE_PHEROMONE = 1;
double MAX_PHEROMONE = 10000.0;

/** Overall algorithm (from "Jones: Artificial Intelligence:  A Systems Approach")
    for each iteration:  
        a.)  distribute ants evenly amongst the graph
        b.)  performing a tour (based on pheromone strength at each ant step) with every ant:
             i.)  at the end of an ant's tour, calculate the distance travelled
             ii.)  distribute pheromone in the tour (amount of pheromone deposited at each edge is inversely proportional to distance of tour)
        c.)  evaporating some of the pheromone on all edges
*/

struct ant
{  vector<bool> visitedCities;
   vector<int> tour;

   double tourLength;
   int tourIndex;
   ant(){}
   ant(int numCities, int initCity)
   {  visitedCities.resize(numCities);
      tour.resize(numCities);
      tourLength = 0;
      tourIndex = 0;

      //init vectors
      for (int i = 0; i < numCities; i++)
        visitedCities[i] = false;
      for (int i = 0; i < numCities; i++)
        tour[i] = -1;  

      visitedCities[initCity] = true;
      tour[0] = initCity;
   }

    void evaluateDistance(double** distMatrix)
    {   //calculate distance of the path 
        tourLength = 0;
        for (int i = 0; i < tour.size(); i++)
          tourLength+= distMatrix[tour[i]][tour[(i+1)%tour.size()]]; //mod for wraparound
    }
};
   
struct city {
  float x, y;
  string name;
  int index;
};

double** genDistMatrix(const vector<city> cities);
double** genPheromMatrix(int numCities);
void chooseNextCity(ant &curAnt, int numCities, double **pheromMatrix, double **distMatrix);
void intensifyPheromoneTrails(ant &curAnt, int numCities, double **pheromMatrix);
void antTSP(vector<city> cities, unsigned timeout);
double genRandom(); //generate random double between 0 and 1.

ostream& operator<<(ostream& os, const city& c) {
  os << c.name << " (" << c.x << ", " << c.y << ")";
  return os;
}

ostream& operator<<(ostream& os, const vector<int> &p) {
  os << "(";
  for (unsigned i = 0; i < p.size(); i++)
  {  if (i != p.size()-1)
       os << p[i] << ",";
     else
       os << p[i] << ")";
  }
  return os;
}
ostream& operator<<(ostream& os, const vector<bool> &p) {
  os << "(";
  for (unsigned i = 0; i < p.size(); i++)
  {  if (i != p.size()-1)
       os << p[i] << ",";
     else
       os << p[i] << ")";
  }
  return os;
}

istream& operator>>(istream& is, city &c) {
  is >> c.x;
  is >> c.y;
  is.get();
  getline(is, c.name);
  return is;
}

int main(int argc, char **argv) {
  string path;
  unsigned cores = 0, timeout;
  if(argc == 4) {
    path = argv[3];
    timeout = atoi(argv[2]);
    cores = atoi(argv[1]);
  } else if(argc == 3) {
    timeout = atoi(argv[1]);
    path = argv[2];
  } else {
    cerr << "Usage: " << argv[0] << " [threads] <timeout> <datafile>" << endl;
    return 1;
  }
  cores = cores ? cores : omp_get_num_procs();
  omp_set_num_threads(cores);
  cout << "Number of threads: " << cores << endl;
  cout << "Data file: " << path << endl;

  srand ( time(NULL) );

  ifstream datafile;
  datafile.open(path.c_str(), ios::in);

  if(!datafile.is_open()) {
    cerr << "Failed to open datafile: " << strerror(errno) << endl;
    return 2;
  }

  size_t length;
  if(path.substr(path.size()-3, path.size()) == "tsp") {
    string line;
    getline(datafile, line);
    getline(datafile, line);
    getline(datafile, line);
    getline(datafile, line);
    length = atoi(line.substr(line.find(':')+1).c_str());
    vector<city> cities(length);
    getline(datafile, line);
    if(line != "EDGE_WEIGHT_TYPE : EUC_2D") {
      cerr << "Noneuclidian inputs not supported!" << endl;
      return 6;
    }
    getline(datafile, line);
    for(size_t i = 0; i < length; ++i) {
      datafile >> cities[i].index;
      datafile >> cities[i].x;
      datafile >> cities[i].y;
    }
    antTSP(cities, timeout);
    for(vector<city>::iterator i = cities.begin(); i != cities.end(); ++i) {
      cout << "City " << (*i).index << ":" << *i << endl;
    }

  } 
  else 
  {  datafile >> length;
     vector<city> cities(length);
     for(size_t i = 0; i < length; ++i)
     {  datafile >> cities[i];
	    cities[i].index = i;    //encode each city as a numerical index for easy retrieval and path representation
   
	    if(datafile.eof()) 
        {  cerr << "Unexpected EOF in datafile!" << endl;
	       return 3;
	    } 
        else if(datafile.bad()) 
        { cerr << "Error reading from datafile: " << strerror(errno) << endl;
	      return 4;
	    } 
        else if(datafile.fail()) 
        { cerr << "Malformed datafile!" << endl;
	      return 5;
	    }
    }
     antTSP(cities, timeout);
    for(vector<city>::iterator i = cities.begin(); i != cities.end(); ++i)
	  cout << "City " << (*i).index << ":" << *i << endl;
      
  }

  return 0;
}

//return ordered path corresponding to best path found
//note that this is easily parallelizable since each ant works independently.  
void antTSP(vector<city> cities, unsigned timeout)
{   cout << "Got " << cities.size() << " cities.  Setting one ant per city." << endl;
    
    
    int numCities = cities.size();
    numAnts = numCities;
    vector<ant> ants(numAnts);
    BASE_PHEROMONE = (double)1.0 /numCities;

    //generate distance matrix for the complete tsp graph
    double** distMatrix = genDistMatrix(cities);
    //initialise pheromone matrix to BASE_PHEROM
    double** pheromMatrix = genPheromMatrix(numCities);
    
    //from here on:  refer to paths as encoded indices of cities only (optimization) -- don't need x, y, or names anymore
    
    vector<int> bestPathSoFar(numCities);
    double bestPathDistanceSoFar=HUGE_VAL;

    //for POSIX:
    //struct timespec zero;
    //clock_gettime(CLOCK_MONOTONIC, &zero);

    clock_t beginTicks = clock();
    clock_t endTicks = clock() + timeout*CLK_TCK;
    
    for (int iteration = 1; true; ++iteration)
      {
        //cout << "a" << endl;
        //a.)  distribute ants evenly amongst the graph.  assign one ant to each city.
        #pragma omp parallel for
        for (int i = 0; i < numAnts; i++)
        {   ants[i] = ant(numCities, i);
        }

        //cout << "b" << endl;
        //b.) get all ants to perform a tour.  this corresponds to telling them to go to a nextCity() NUMCITIES-1 times (since a tour == NUMCITIES, and we already chose the first city for them)
        for (int i = 0; i < numCities-1; i++)
        {   //parallelize the ants (they create their own tours independently)
            #pragma omp parallel for
            for (int antNum = 0; antNum < numAnts; antNum++)
            {   chooseNextCity(ants[antNum], numCities, pheromMatrix, distMatrix);
            }
            #pragma omp barrier
        }

        //cout << "c" << endl;
        //c.)  distribute pheromone across the tour (amount deposited = inversely proportional to distance of tour)
        /***** THIS HAS A CRITICAL SECTION */
        #pragma omp parallel for
        for (int i = 0; i < numAnts; i++)
        {   intensifyPheromoneTrails(ants[i], numCities, pheromMatrix);
        }

        //cout << "d" << endl;
        //d.)  evaporate pheromone trails for all edges in the graph since it's the end of an iteration
        #pragma omp parallel for
        for (int i = 0; i < numCities; i++)
        {   for (int j = 0; j < numCities; j++)
            {   pheromMatrix[i][j] = pheromMatrix[i][j] * (1-p);
                if (pheromMatrix[i][j] < 0)
                    pheromMatrix[i][j] = BASE_PHEROMONE;
            }
        }

        //cout << "e" << endl;
        //e.)  finally, find the best path from this iteration by going through the pheromone trail.  start at first city
        vector<int> bestPath(numCities);
        vector<bool> visited(numCities);
        for (int i = 0; i < numCities; i++)
        {  bestPath[i] = -1;
           visited[i] = false;
        }
        bestPath[0] = 0;
        visited[0] = true;
        for (int pathIndex = 0; pathIndex < numCities-1; pathIndex++)
        {   int bestCity = -1;
            double bestPherom = -1;
            int from = bestPath[pathIndex];

            for (int to = 0; to < numCities; to++)
            {  if (visited[to] == false)
               {  if (bestPherom < pheromMatrix[from][to])
                  {  bestPherom = pheromMatrix[from][to];
                     bestCity = to;
                  }
               }
            }
            bestPath[pathIndex+1] = bestCity;
            visited[bestCity] = true;
        }
        double dist = 0.0;
        for (int i = 0; i < bestPath.size(); i++)
        {   dist+= distMatrix[bestPath[i]][bestPath[(i+1)%bestPath.size()]]; //mod for wraparound
        }

        if (dist < bestPathDistanceSoFar)
        {  bestPathDistanceSoFar = dist;
           bestPathSoFar = bestPath;
        }

        #ifndef _WIN32
	    struct timespec now;
	    clock_gettime(CLOCK_MONOTONIC, &now);
        float dt = (now.tv_sec - zero.tv_sec) + 1e-9*(now.tv_nsec - zero.tv_nsec);
        cout << iteration << "," << dt << "," << bestPathDistanceSoFar << endl;
        if(dt > timeout) {  exit(0); }
        #endif

        #ifdef _WIN32
        clock_t curTicks = clock();
        
        double dt = (curTicks-beginTicks)/(double)CLOCKS_PER_SEC;
        cout << "Iteration:" << iteration << " Time: " << dt << " Best Distance Found:" << bestPathDistanceSoFar << endl;

        if (curTicks >= endTicks)
        {   cout << "Iterations completed:" << iteration << endl;
            exit(0);
        }
        #endif
    }
}

void chooseNextCity(ant &curAnt, int numCities, double **pheromMatrix, double **distMatrix)
{  int from = curAnt.tour[curAnt.tourIndex];
   double denominator = 0.0;

   if (curAnt.tourIndex == numCities-1)
   {  cout << " this shouldn't happen " << endl;
   }

   //calculate denominator in probabilistic equation.  this corresponds to the sum of all non visited cities' edges' pheromones * the edge's visibility (1/dist)
   for (int to = 0; to < numCities; to++)
   {  if (curAnt.visitedCities[to] == false)
      {  denominator += pow(pheromMatrix[from][to], alpha) * pow(1.0 / distMatrix[from][to], beta);
      }
   }

   //go through all cities.  calculate probability for taking path to next city, and roll 'dice' to see if we should take it.  we try this 10 times, and if we still don't have a result, just choose the path with the highest probability (to prevent deadlocks)
   double bestProb = -1;
   int bestCity = -1;

   int nextCity = -1;
   bool foundCity = false;
   double prob = 0.0;
   for (int tries = 0; tries < 5; tries++)
   {  for (int to = 0; to < numCities; to++)
      {  if (curAnt.visitedCities[to] == false) //only consider cities for which we haven't visited yet...
         {  
            if (denominator == 0 || distMatrix[from][to] == 0) //resolve divide by zero problem
               prob = 1.0;
            else
               prob = min(1.0,(pow(pheromMatrix[from][to],alpha) * pow(1.0 / distMatrix[from][to], beta)) / denominator);
                
            if (bestProb <= prob) //to prevent deadlocking situation, keep track of the best cities found 
            {  bestProb = prob;
               bestCity = to;
            }
            if (genRandom() <= prob)
            {  foundCity = true;
               bestCity = to;
               break;
            }
         }
      }
      if (foundCity)
         break;
   }
   
   nextCity = bestCity;

   if (nextCity < 0)
   {  cout << "Incorrect nextCity" << endl;
      cout << "Ant tourIndex:" << curAnt.tourIndex << endl;
      bool hasUnvisited = false;
      for (int i = 0; i < curAnt.visitedCities.size(); i++)
      { if (curAnt.visitedCities[i] == false)
        {  hasUnvisited = true;
           cout << "Unvisited city:" << i << endl;
           break;
        }
      }
      cout << "Ant has unvisited cities:" << hasUnvisited << endl;
      cout << "Calculated probability:"  << prob << endl;

   }

   curAnt.tourIndex+= 1;
   curAnt.tour[curAnt.tourIndex] = nextCity;
   curAnt.tourLength+= distMatrix[from][nextCity];
   curAnt.visitedCities[nextCity] = true;

   if (curAnt.tourIndex == numCities-1) //i.e. we're at the end of the tour, so update the last edge anyways
   {   curAnt.tourLength+= distMatrix[curAnt.tour[numCities-1]][curAnt.tour[0]];
   }
}

void intensifyPheromoneTrails(ant &curAnt, int numCities, double **pheromMatrix)
{   for (int i = 0; i < numCities; i++)
    { //go through each edge in curAnt's path
        int from = curAnt.tour[i];
        int to = curAnt.tour[(i+1)%numCities];

        //increase pheromone at edge in path
        /** THIS HAS TO BE IN A CRITICAL SECTION!!!!!!!! **/
        #pragma omp critical
        {  
           #pragma omp flush(pheromMatrix) //flush so that all threads' view of pheromMatrix is consistent
           pheromMatrix[from][to] += (Q/curAnt.tourLength)*p;
           if (pheromMatrix[from][to] > MAX_PHEROMONE)
              pheromMatrix[from][to] = MAX_PHEROMONE;
           pheromMatrix[to][from] = pheromMatrix[from][to];
        }
    }
}

//helper func: each entry in the returned distance matrix corresponds to the distance between city(i,j)
double** genDistMatrix(vector <city> cities)
{   int numCities = cities.size();   
    double** retDistMatrix = new double*[numCities];
    
    //can parallelize this if we want to using openMP (independent for)
    #pragma omp parallel for
    for (int i = 0; i < numCities; i++)
    {   retDistMatrix[i] = new double[numCities];
        //calculate distance from city i to all other cities:
        for (int j = 0; j < numCities; j++)
        {   retDistMatrix[i][j] = sqrt(pow(fabs(cities[i].x - cities[j].x),2) + pow(fabs(cities[i].y - cities[j].y),2));
        }
    }
    return retDistMatrix;
}
//helper func: initialise pherome matrix
double** genPheromMatrix(int numCities)
{   double **pheromMatrix = new double*[numCities];
    #pragma omp parallel for
    for (int i = 0; i < numCities; i++)
    {  pheromMatrix[i] = new double[numCities];
       for (int j = 0; j < numCities; j++)
          pheromMatrix[i][j] = BASE_PHEROMONE;
    }
    return pheromMatrix;
}
double genRandom()
{   return (double)rand()/(double)RAND_MAX;
}
        
