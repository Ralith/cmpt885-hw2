#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>

using namespace std;

struct city {
  float x, y;
  string name;
};

ostream& operator<<(ostream& os, const city& dt) {
  os << dt.name << " (" << dt.x << ", " << dt.y << ")";
  return os;
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
  char namebuf[256];
  for(size_t i = 0; i < length; ++i) {
    datafile >> cities[i].x;
    datafile >> cities[i].y;
    datafile.get();
    datafile.get(namebuf, 256);
    cities[i].name = namebuf;

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

  for(vector<city>::iterator i = cities.begin(); i != cities.end(); ++i) {
    cout << *i << endl;
  }

  return 0;
}
