/***********************************************************
 * Copyright (C) 2013 PYXAID group
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#include <boost/python.hpp>
#include <iostream>
#include <iomanip>
#include "wfc_export.h"
#include "namd_export.h"
using namespace std;
using namespace boost::python;

class info {
public:
  // Constructor
  info() {}

  // Class members
  void version() {
    std::cout
        << "================================================================================\n";
    std::cout << "PYXAID2: PYthon eXtension for Ab Inition Dynamics version 2.0\n";
    std::cout << "/***********************************************************\n";
    std::cout << " * Copyright (C) 2017 PYXAID2 group\n";
    std::cout << " * This program is free software distributed under the terms of the\n";
    std::cout << " * GNU General Public License as published by the\n";
    std::cout << " * Free Software Foundation; either version 3 of the\n";
    std::cout << " * License, or (at your option) any later version.\n";
    std::cout << " * http://www.gnu.org/copyleft/gpl.txt\n";
    std::cout << "***********************************************************/\n";
  }

  void developers() {
    std::cout << "===== Name ========================== Contact info ==============\n";
    std::cout << "Alexey V. Akimov              e-mail: alexvakimov@gmail.com      \n";
    std::cout << "Wei Li                        e-mail: liwei0099@gmail.com\n";
    std::cout << "...\n";
  }

  void documentation() { std::cout << "Coming soon...\n"; }
};

BOOST_PYTHON_MODULE(pyxaid_core) {
  class_<info>("info", init<>())
      .def("version", &info::version)
      .def("developers", &info::developers)
      .def("documentation", &info::documentation);

  export_wfc();
  export_namd();
}
