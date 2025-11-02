/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#ifndef pyxaid_state_h
#define pyxaid_state_h

#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/python.hpp>
using namespace boost::python;
using namespace std;

#include "liblibra_core.h"
using namespace liblibra;
using namespace liblibra::libutil;

int ext2int(int, vector<int>&);
//int delta(vector<int>& A,vector<int>& B,int& a,int& b);

// Multi-electron state
// Basically it is a Slater product
class me_state {
public:
  std::string name;  // label of the determinant
  // Data
  std::vector<int> active_space;  // All occupied and virtual orbitals involved in dynamics
  std::vector<int> actual_state;  // This is an actual(current) state

  double Eshift;  // this correction goes to Exc and is read directly from input - this is
                  // to simplify parameter development
  double Exc;     // correlation correction, to better describe the energy
                  // of the state with 2 electrons on the same orbital

  // NAC scaling constants for different pairs of the states
  // the sizes of the below vectors should be equal and elements be in correspondense
  std::vector<int> nac_scl_indx;
  std::vector<double> nac_scl;

  me_state() = default;
  me_state(std::vector<int>& as_, std::vector<int>& cs_)
      : active_space{as_}, actual_state{cs_}, Exc{0.} {}

  void set_me_state(std::vector<int>& as_, std::vector<int>& cs_) {
    active_space = as_;
    actual_state = cs_;
    Exc = 0.0;
  }

  ~me_state() = default;

  int calculate_Exc(std::vector<int>&,
                    std::vector<int>&,
                    std::vector<double>&,
                    std::vector<int>&,
                    std::vector<double>&);
  void show_state();
};

void input_iconds(boost::python::dict params,
                  int me_numstates,
                  std::vector<std::vector<int> >& icond);
void input_states(boost::python::dict params, std::vector<me_state>& states);

#endif  // state_h
