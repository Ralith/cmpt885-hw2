#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

//define these as user input later if necessary
int population_size = 50000; //population should always be an even number
int generations = 100;

struct city {
  float x, y;
  string name;
  int index;
};

struct calculatedpath {
    vector<int> path; //vector of cities, correspnding a path (note: each city is referred to by its int index)
    double distance;

    calculatedpath(){}
    calculatedpath(vector <int> p)
    {   path = p;
    }
    
    void evaluateDistance(double** distMatrix)
    {   distance = 0;
        //calculate distance of the path
        for (int i = 0; i < path.size(); i++)
        {   if (i != path.size()-1)
                distance+= distMatrix[path[i]][path[i+1]];
            else //wraparound since a tsp path is a cycle
                distance+= distMatrix[path[i]][path[0]];
        }
    }
    
    bool operator < (const calculatedpath &other) const
    {    return distance < other.distance;
    }
};


vector<city> geneticTSP(vector<city> cities);
double** genDistMatrix(vector<city> cities);
void genInitialPopulation(vector<calculatedpath> &retPop, int numCities, double** distMatrix);
void crossover(calculatedpath &child, calculatedpath parent1, calculatedpath parent2, double **distMatrix);

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
  
  for(vector<city>::iterator i = cities.begin(); i != cities.end(); ++i) {
    cout << "City " << (*i).index << ":" << *i << endl;
  }
  
  return 0;
}

//return ordered path corresponding to best path found
vector<city> geneticTSP(vector<city> cities)
{   /**initialise**/
    
    //generate distance matrix for the complete tsp graph
    double** distMatrix = genDistMatrix(cities);
    
    //from here on:  refer to paths as encoded indices of cities only (optimization) -- don't need x, y, or names anymore
    //generate an initial population (i.e. paths) to seed the genetic algorithm
    vector<calculatedpath> population(population_size);
    genInitialPopulation(population, cities.size(), distMatrix);
    
    //"evolve" the seeded population for a specified number of times.
    for (int generation = 1; generation <= generations ; generation++)
    {/* 1.)  EVALUATE the distance for each path in the population
             Can easily parallelize this */
       for (int i = 0; i < population.size(); i++)
       {   population[i].evaluateDistance(distMatrix);
       }
               
     /* 2.)  Sort the population and SELECT only the best half of the current
             population for 'mating' */
       
       sort(population.begin(), population.end());
       vector<calculatedpath> bestpop(population.size()/2);
       copy(population.begin(), (population.end() - (int)(population.size()/2)), bestpop.begin());
       cout << "Best population in generation " << generation << ": " << bestpop[0] << endl;
       
       
     /* 3.)  CROSSOVER the best half to reproduce children 
             Use modified Grefenstette Greedy Crossover 
             a.)  Child[i] takes Parent[i]'s first city.  Then:
                    for each city x in child[i]:
                        select cities a,b from x in either parent[i]|[i+1] path.
                               -  choose closest a,b, and if DNE in child[i] then
                                  extend
                               -  o/w, if a or b exists, choose other that DNE and extend
                               -  o/w, if both exist in child[i], choose random unchosen city
             Loop across all children.  */
             
       vector<calculatedpath> children(population.size()/2); 
       for (int i = 0; i < children.size(); i++)
       {   
           if (i != children.size()-1)
              crossover(children[i], bestpop[i], bestpop[i+1], distMatrix);
           else
              crossover(children[i], bestpop[i], bestpop[0], distMatrix);
       }
       
     /* 4.)  Randomly MUTATE some of the children */
     
     /* 5.)  POPULATE the new population by replacing poor population with children*/
      
      copy(children.begin(), children.end(), population.begin() + (int)population.size()/2);
        
       /*
       vector<calculatedpath> worstpop(population.size()/2);
       copy(population.begin()+(int)(population.size()/2), population.end(), worstpop.begin());
*/
       
   /*            
      for (int i = 0; i < bestpop.size(); i++)
       {   cout << bestpop[i] << endl;
       }
       cout << "children: " << endl;
       for (int i = 0; i < children.size(); i++)
       {   cout << children[i] << endl;
       }
       */
               
    }          
       
       
               
    
    

    //mutate    //with mutation_likelihood, flip two cities in each trip
    //populate  //generate new population by taking parents and newly formed children
    
    

        
    return cities;
}
//Crossover algorithm based on Grefenstette's modified Greedy Crossover.  See above for pseudocode
void crossover(calculatedpath &child, calculatedpath parent1, calculatedpath parent2, double **distMatrix)
{    vector<int> path(parent1.path.size());
     
     bool hasUsedCity[path.size()];
     for (int i=0; i < path.size(); i++)
     {   hasUsedCity[i] = false;
     }
     
     //initially, child takes first city of parent
     path[0] = parent1.path[0];
     hasUsedCity[parent1.path[0]] = true;
     
     for (int i = 0; i < path.size()-1; i++)
     {   //find canditate cities connected to current city in child         
         int parent1Index = distance(parent1.path.begin(), find(parent1.path.begin(), parent1.path.end(), path[i]));
         int nextCity1 = (parent1Index != parent1.path.size()-1)?parent1.path[parent1Index+1]:parent1.path[0]; //wraparound in-case
         int parent2Index = distance(parent2.path.begin(), find(parent2.path.begin(), parent2.path.end(), path[i]));
         int nextCity2 = (parent2Index != parent2.path.size()-1)?parent2.path[parent2Index+1]:parent2.path[0]; //wraparound in-case
         
         //see which has a shorter distance, simple lookup
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
             for (int x = 0; x < path.size(); x++)
             {   if (!hasUsedCity[x])
                    availCities.push_back(x);
             }
             int randCity = availCities[rand()%availCities.size()];
             path[i+1] = randCity;
             hasUsedCity[randCity] = true;
         }
     }
     
     child = calculatedpath(path);
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
void genInitialPopulation(vector<calculatedpath> &retPop, int numCities, double** distMatrix)
{   vector<int> unshuffledPath(numCities);
    for (int i = 0; i < numCities; i++)
        unshuffledPath[i] = i;
            
    for (int i = 0; i < retPop.size(); i++)
    {   vector<int> randPath = unshuffledPath;
        random_shuffle(randPath.begin(), randPath.end());
        retPop[i] = calculatedpath(randPath);
    }
    
    delete &unshuffledPath;
}
        
