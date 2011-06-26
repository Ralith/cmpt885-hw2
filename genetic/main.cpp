#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

struct city {
  float x, y;
  string name;
  int index;
};

struct calculatedpath {
    vector<int> path; //vector of cities, correspnding a path (note: each city is referred to by its int index)
    double distance;
    calculatedpath(vector<int> p, double** distMatrix)
    {   distance = 0;
        path = p;
        //calculate distance of the path
        for (int i = 0; i < path.size(); i++)
        {   if (i != path.size()-1)
                distance+= distMatrix[path[i]][path[i+1]];
            else //wraparound since a tsp path is a cycle
                distance+= distMatrix[path[i]][path[0]];
        }
        
    }
    calculatedpath(){}
};

vector<city> geneticTSP(vector<city> cities);
double** genDistMatrix(vector<city> cities);
vector<calculatedpath> genInitialPopulation(int numCities, double** distMatrix);

//define these as user input later if necessary
int population_size = 500;
int generations = 500;

ostream& operator<<(ostream& os, const city& c) {
  os << c.name << " (" << c.x << ", " << c.y << ")";
  return os;
}

ostream& operator<<(ostream& os, const calculatedpath& p) {
  os << "Path dist:" << p.distance << ", Path: (";
  for (int i = 0; i < p.path.size(); i++)
  {  os << p.path[i] << ",";
  }
  os << ")";
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
  }

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
  
  srand ( time(NULL) );
  geneticTSP(cities);
  
 /* for(vector<city>::iterator i = cities.begin(); i != cities.end(); ++i) {
    cout << *i << endl;
  }*/
  
  return 0;
}

//return ordered path corresponding to best path found
vector<city> geneticTSP(vector<city> cities)
{   /**initialise**/
    
    //generate distance matrix for the complete tsp graph
    double** distMatrix = genDistMatrix(cities);
    
    //from here on:  refer to paths as encoded indices of cities only (optimization) -- don't need x, y, or names anymore
    //generate an initial population (i.e. paths) to seed the genetic algorithm
    vector<calculatedpath> population = genInitialPopulation(cities.size(), distMatrix);
    
    for (int i = 0; i < population.size(); i++)
    {   cout << population[i] << endl;
    }
    //evalute   //evaluate the distance for each path in the population
    //select    //select the best half of the current population for mating
    //crossover //mate each pair in the selected best population
    //mutate    //with mutation_likelihood, flip two cities in each trip
    //populate  //generate new population by taking parents and newly formed children
    
    return cities;
}

//helper func: each entry in the returned distance matrix corresponds to the distance between city(i,j)
double** genDistMatrix(vector<city> cities)
{   int numCities = cities.size();   
    double** retDistMatrix = new double*[numCities];
    
    //can parallelize this if we want to using openMP (independent for)
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
vector<calculatedpath> genInitialPopulation(int numCities, double** distMatrix)
{   vector<calculatedpath> initialPop(population_size);

    vector<int> unshuffledPath(numCities);
    for (int i = 0; i < numCities; i++)
        unshuffledPath[i] = i;
    
    for (int i = 0; i < population_size; i++)
    {   vector<int> randPath= unshuffledPath;
        random_shuffle(randPath.begin(), randPath.end());
        initialPop[i] = calculatedpath(randPath, distMatrix);
    }
    return initialPop;
}
        
