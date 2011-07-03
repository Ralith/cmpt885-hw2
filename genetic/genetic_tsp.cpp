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

//define these as user input later if necessary
int population_size = 20000; //population should always be an even number
int generations = 300;
double mutation_likelihood = 0.1;

//tournament selection
int tournamentSize = 5;
double tournamentProb = 0.6;

struct city {
  float x, y;
  string name;
  int index;
};

struct calculatedpath {
    vector<int> path; //vector of cities, correspnding a path (note: each city is referred to by its int index)
    double distance;

    calculatedpath(){}
    calculatedpath(vector<int> &p)
    {   path = p;
        distance = 0;
    }
    
    void evaluateDistance(double** distMatrix)
    {   //calculate distance of the path 
        distance = 0;
        for (int i = 0; i < path.size(); i++)
        {   distance+= distMatrix[path[i]][path[(i+1)%path.size()]]; //mod for wraparound
        }
    }
    
    bool operator < (const calculatedpath &other) const
    {    return distance < other.distance;
    }
    calculatedpath& operator = (const calculatedpath &other) 
    {    path = other.path;
         distance = other.distance;
         return *this;
    }
};


calculatedpath geneticTSP(vector<city> &cities, unsigned timeout);
double** genDistMatrix(const vector<city> &cities);
void genInitialPopulation(vector<calculatedpath> &retPop, unsigned numCities);
void crossover(calculatedpath &child, const calculatedpath &parent1, const calculatedpath &parent2, double **distMatrix);

ostream& operator<<(ostream& os, const city& c) {
  os << c.name << " (" << c.x << ", " << c.y << ")";
  return os;
}

ostream& operator<<(ostream& os, const calculatedpath& p) {
  os << "Path dist:" << p.distance; // ", Path: (";
  /*for (unsigned i = 0; i < p.path.size(); i++)
  {  if (i != p.path.size()-1)
       os << p.path[i] << ",";
     else
       os << p.path[i] << ")";
  }*/
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
    geneticTSP(cities, timeout);
    for(vector<city>::iterator i = cities.begin(); i != cities.end(); ++i) {
      cout << "City " << (*i).index << ":" << *i << endl;
    }

  } else {
    datafile >> length;
    vector<city> cities(length);
    for(size_t i = 0; i < length; ++i)
      {  datafile >> cities[i];
	cities[i].index = i;    //encode each city as a numerical index for easy retrieval and path representation
 
	if(datafile.eof()) {
	  cerr << "Unexpected EOF in datafile!" << endl;
	  return 3;
	} else if(datafile.bad()) {
	  cerr << "Error reading from datafile: " << strerror(errno) << endl;
	  return 4;
	} else if(datafile.fail()) {
	  cerr << "Malformed datafile!" << endl;
	  return 5;
	}
      }
    geneticTSP(cities, timeout);
      for(vector<city>::iterator i = cities.begin(); i != cities.end(); ++i) {
	cout << "City " << (*i).index << ":" << *i << endl;
      }
  }
  
  return 0;
}

//return ordered path corresponding to best path found
//note that this is easily parallelizable since during each generation, evaluating distances for each candidate path is independent of other paths.  
//similarly, crossing over to generate new children is also independent from child to child.
calculatedpath geneticTSP(vector<city> &cities, unsigned timeout)
{
  cout << "Got " << cities.size() << " cities." << endl;
    //generate distance matrix for the complete tsp graph
    double** distMatrix = genDistMatrix(cities);
    
    //from here on:  refer to paths as encoded indices of cities only (optimization) -- don't need x, y, or names anymore
    //generate an initial population (i.e. paths) to seed the genetic algorithm
    vector<calculatedpath> population(population_size);
    genInitialPopulation(population, cities.size());

    //calculate all of the distances for the initial population first.
    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
    {   population[i].evaluateDistance(distMatrix);
    }
    
    struct timespec zero;
    clock_gettime(CLOCK_MONOTONIC, &zero);
    //"evolve" the seeded population for a specified number of generations.
    for (unsigned generation = 1; true; generation++)
    {
 /***** 1.)  SELECT only the best half of the population for 'mating' */

       //get minimum path for this generation and make sure it's in the bestpop.  we don't want to lose any good paths to the selection process!
       float min = HUGE_VAL;
       int minIndex = -1;
       for (int i = 0; i < population.size(); i++)
       {  if (population[i].distance < min)
          {  min = population[i].distance;
             minIndex = i;
          }
       }
       struct timespec now;
       clock_gettime(CLOCK_MONOTONIC, &now);
       float dt = (now.tv_sec - zero.tv_sec) + 1e-9*(now.tv_nsec - zero.tv_nsec);
       cout << generation << "," << dt
	    << "," << population[minIndex].distance << endl;

       if(dt > timeout) {
	 exit(0);
       }

       //tournament elitist selection
       //http://en.wikipedia.org/wiki/Tournament_selection
       vector <calculatedpath> bestpop(population.size()/2);

       int numSelected = 0;
       while (numSelected != population.size()/2)
       {  vector<calculatedpath> curTournament(tournamentSize);
          for (unsigned i = 0; i < tournamentSize; i++)
          {  curTournament[i] = population[rand()%population.size()];
          }
          sort(curTournament.begin(), curTournament.end());
          for (int i = 0; i < tournamentSize; i++)
          {   if (rand()%10000 < (double)(tournamentProb*pow((double)(1-tournamentProb),i)*10000))
              {  bestpop[numSelected] = curTournament[i];
                 numSelected++;
                 if (numSelected == population.size()/2)
                    break;
              }
          }                  
       }
     
 /***** 2.)  CROSSOVER the best half of the population to reproduce children */
       vector<calculatedpath> children(bestpop.size()); //every two parent pairs creates one child, i.e. #children == #bestpop, and #children + #bestpop == population
       #pragma omp parallel for
       for (int i = 0; i < children.size(); i++)
       {  crossover(children[i], bestpop[i], bestpop[(i+1)%bestpop.size()], distMatrix); //mod for wraparound
       }
       

 /***** 3.)  Randomly MUTATE some of the children.  This corresponds to just flipping two nodes in the path, and keeping it only if the path is an improvement. */
       #pragma omp parallel for
       for (int i = 0; i < children.size(); i++)
       {   //mutate
           if (rand()%100 < (mutation_likelihood*100)) //random number from 0-99.  accurate to 2 decimal places
           {  for (int j = 0; j < 5; j++) //swap a few cities
              {  swap(children[i].path[rand()%children[i].path.size()], children[i].path[rand()%children[i].path.size()]);
              }
           }
       }

 /***** 4.)  EVALUATE the distance for each newly created child
             Can easily parallelize this */
       #pragma omp parallel for
       for (int i = 0; i < children.size(); i++)
       {  //only compute distances if we haven't done so yet -- this is an optimization step since the top half of the population's distances will have been computed already
          if (children[i].distance == 0)
            children[i].evaluateDistance(distMatrix);
          else
            cout << " WTF " << endl;
       }
       
 /***** 5.)  POPULATE the new population by having the top half as the best selected parents, and the bottom half as the children of those parents */
      copy(bestpop.begin(), bestpop.end(), population.begin());
      copy(children.begin(), children.end(), population.begin() + (int)(bestpop.size()));
      


      //onto the next generation...
    }

    //finished algorithm.  sort to find path with lowest distance (lazy)
    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
       population[i].evaluateDistance(distMatrix);
    sort(population.begin(), population.end());
    
    return population[0];
}

/* Crossover function to generate a new child from two parents.
   Use modified Grefenstette Greedy Crossover: 
      Child path initially only has parent1's first city.  Then:
         for each city x in Child's path:
            select next cities a,b from x in both parents' paths.
                a.)  choose closest a,b, and if city DNE in Child's path then extend path with it
                b.)  o/w, if one city already exists in Child's path, choose other one that DNE and extend with it
                c.)  o/w, if both cities exist in Child's path, choose random unchosen city and extend with it
*/
void crossover(calculatedpath &child, const calculatedpath &parent1, const calculatedpath &parent2, double **distMatrix)
{    vector<int> path(parent1.path.size());
     
     //hash map to quickly check if a city already exists in child's path
	 vector<bool> hasUsedCity(path.size());
     for (unsigned i=0; i < hasUsedCity.size(); i++)
     {   hasUsedCity[i] = false;
     }
     
     //initially, child takes first city of parent
     path[0] = parent1.path[0];
     hasUsedCity[parent1.path[0]] = true;
     
     for (unsigned i = 0; i < path.size()-1; i++)
     {         
         //find candidate cities connected to current city in child         
         int parent1Index = distance(parent1.path.begin(), find(parent1.path.begin(), parent1.path.end(), path[i]));
         int nextCity1 = parent1.path[(parent1Index+1)%parent1.path.size()]; //wraparound in-case
         int parent2Index = distance(parent2.path.begin(), find(parent2.path.begin(), parent2.path.end(), path[i]));
         int nextCity2 = parent2.path[(parent2Index+1)%parent2.path.size()]; //wraparound in-case
         
         //see which candidate city has a shorter distance, simple lookup
         int closerCity = (distMatrix[path[i]][nextCity1] <= distMatrix[path[i]][nextCity2])?nextCity1:nextCity2;
         int fartherCity = (distMatrix[path[i]][nextCity1] <= distMatrix[path[i]][nextCity2])?nextCity2:nextCity1;
         
         //now we try setting the next city in the child to the closer city, if possible
         if (!hasUsedCity[closerCity]) //haven't used it yet, so set it as next
         {  path[i+1] = closerCity;
            hasUsedCity[closerCity] = true;
         }
         else if (!hasUsedCity[fartherCity]) //closerCity has been used, use other
         {  path[i+1] = fartherCity;
            hasUsedCity[fartherCity] = true;
         } //both cities have been used, randomly choose one that hasn't then
         else
         {   vector<int> availCities;
             for (unsigned x = 0; x < path.size(); x++)
             {   if (!hasUsedCity[x])
                    availCities.push_back(x);
             }
             int randCity = availCities[rand()%availCities.size()];
             path[i+1] = randCity;
             hasUsedCity[randCity] = true;
         }
     }     
     child = calculatedpath(path); //encapsulate newly created path into our data structure for easier dist calc/comparison later.
}
     
//helper func: each entry in the returned distance matrix corresponds to the distance between city(i,j)
double** genDistMatrix(const vector<city> &cities)
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
//helper func:  generate an initial population (i.e. paths) to seed the genetic algorithm
void genInitialPopulation(vector<calculatedpath> &retPop, unsigned numCities)
{   vector<int> unshuffledPath(numCities);
    for (unsigned i = 0; i < numCities; i++)
        unshuffledPath[i] = i;
    
    #pragma omp parallel for
    for (int i = 0; i < retPop.size(); i++)
    {   vector<int> randPath = unshuffledPath;
        random_shuffle(randPath.begin(), randPath.end());
        retPop[i] = calculatedpath(randPath);
    }
}
        
