#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "time.h"
#include "omp.h"

using namespace std;

//define these as user input later if necessary
int population_size = 20000; //population should always be an even number
int generations = 100;
double mutation_likelihood = 0.03;


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
        for (int i = 0; i < path.size(); i++)
        {   distance+= distMatrix[path[i]][path[(i+1)%path.size()]]; //mod for wraparound
        }
    }
    
    bool operator < (const calculatedpath &other) const
    {    return distance < other.distance;
    }
};


calculatedpath geneticTSP(vector<city> &cities);
double** genDistMatrix(const vector<city> &cities);
void genInitialPopulation(vector<calculatedpath> &retPop, unsigned numCities);
void crossover(calculatedpath &child, const calculatedpath &parent1, const calculatedpath &parent2, double **distMatrix);

ostream& operator<<(ostream& os, const city& c) {
  os << c.name << " (" << c.x << ", " << c.y << ")";
  return os;
}

ostream& operator<<(ostream& os, const calculatedpath& p) {
  os << "Path dist:" << p.distance << ", Path: (";
  for (unsigned i = 0; i < p.path.size(); i++)
  {  if (i != p.path.size()-1)
       os << p.path[i] << ",";
     else
       os << p.path[i] << ")";
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
  if(argc != 2) {
    cerr << "Usage: " << argv[0] << " <dataset>" << endl;
    return 1;
  }
  
  ifstream datafile;
  datafile.open(argv[1], ios::in);

  if(!datafile.is_open()) {
    cerr << "Failed to open datafile: " << strerror(errno) << endl;
    return 2;
  }

  size_t length;
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

  srand ( time(NULL) );

  geneticTSP(cities);
  
  for(vector<city>::iterator i = cities.begin(); i != cities.end(); ++i) {
    cout << "City " << (*i).index << ":" << *i << endl;
  }
  
  return 0;
}

//return ordered path corresponding to best path found
//note that this is easily parallelizable since during each generation, evaluating distances for each candidate path is independent of other paths.  
//similarly, crossing over to generate new children is also independent from child to child.
calculatedpath geneticTSP(vector<city> &cities)
{   /**initialise**/

    int cores = omp_get_num_procs();
    cout << "Number of cores: " << cores << endl;
    omp_set_num_threads(6);
    
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
    
    //"evolve" the seeded population for a specified number of generations.
    for (int generation = 1; generation <= generations ; generation++)
    {
     /* 1.)  Sort the population and SELECT only the best (top) half of the current
             population for 'mating' */
       sort(population.begin(), population.end());
       vector<calculatedpath> bestpop(population.size()/2);
       copy(population.begin(), (population.begin() + (int)(population.size()/2)), bestpop.begin());
       cout << "Best path in generation " << generation << ": " << bestpop[0] << endl;
       
     /* 2.)  CROSSOVER the best half of the population to reproduce children */
       vector<calculatedpath> children(bestpop.size()); //every two parent pairs creates one children, i.e. # of children == # of bestpop, and #children + #bestpop == population
       #pragma omp parallel for
       for (int i = 0; i < children.size(); i++)
       {  crossover(children[i], bestpop[i], bestpop[(i+1)%bestpop.size()], distMatrix); //mod for wraparound
       }
       
     /* 3.)  Randomly MUTATE some of the children.  This corresponds to just flipping two nodes in the path, and keeping it only if the path is an improvement. */  
       #pragma omp parallel for
       for (int i = 0; i < children.size(); i++)
       {   //mutate
           if (rand()%100 < (mutation_likelihood*100)) //random number from 0-99.  accurate to 2 decimal places
           {  calculatedpath old = children[i];
              for (int j = 0; j < 2; j++) //swap a few cities
              {  swap(children[i].path[rand()%children[i].path.size()], children[i].path[rand()%children[i].path.size()]);
              }
              children[i].evaluateDistance(distMatrix);
              if (old.distance < children[i].distance)
              {   children[i] = old;
              }
           }
       }

     /* 4.)  EVALUATE the distance for each newly created child
             Can easily parallelize this */
       #pragma omp parallel for
       for (int i = 0; i < children.size(); i++)
       {  //only compute distances if we haven't done so yet -- this is an optimization step since the top half of the population's distances will have been computed already
          if (children[i].distance == 0)
            children[i].evaluateDistance(distMatrix);
          else
            cout << " WTF " << endl;
       }
       
     /* 5.)  POPULATE the new population by replacing poor half of the current population with newly created children */
      copy(children.begin(), children.end(), population.begin() + (int)(population.size()-children.size()));

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
        
