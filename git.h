// git.h --- Gauss Integrals Tuned (GIT) library code
//           Original C-code by Peter Roegen, Copyright (C) 2011
//           Please cite: 
//           P. R{\o}gen, Evaluating protein structure descriptors and tuning Gauss integral
//           based descriptors: Journal of Physics Condensed Matter, vol: 17, pages: 1523-1538, 2005.
//           
//           C++ version by Tim Harder and Wouter Boomsma, Copyright (C) 2011
//
// This file is part of Git
//
// Git is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Git.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef GIT_H
#define GIT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <cstdlib>
#include <cassert>

namespace git {

//! Vector class for 3D coordinates
class Vec3 {

     //! Internal data
     double data[3];
public:

     //! Constructor (default)
     Vec3() {}

     //! Constructor (x,y,z coordinates)
     Vec3(double x, double y, double z) {
          data[0] = x;
          data[1] = y;
          data[2] = z;
     }

     //! Constructor (from array)
     Vec3(double array[3]) {
          data[0] = array[0];
          data[1] = array[1];
          data[2] = array[2];
     }

     //! Overload [] indexing operator (const)
     double operator[](const int index) const {
          return data[index];
     }

     //! Overload indexing operator (non-const)
     double& operator[](const int index) {
          return data[index];
     }

     //! Calculate L2 norm (length) of vector
     //! \return length of vector
     double norm() const {
          return sqrt((*this)*(*this));
     }

     //! Return normalized vector
     //! \return new vector
     Vec3 normalize() const {
          return (*this)/norm();
     }

     //! Overload * operator (dot product)
     friend double operator*(const Vec3 &v1, const Vec3 &v2) {
          return (v1.data[0] * v2.data[0] + 
                  v1.data[1] * v2.data[1] + 
                  v1.data[2] * v2.data[2]);     
     }

     //! Overload * operator
     friend Vec3 operator*(const Vec3 &v, const double value) {
          return Vec3(v.data[0] * value,
                      v.data[1] * value,
                      v.data[2] * value);
     }

     //! Overload * operator
     friend Vec3 operator*(const double value, const Vec3 &v) {
          return operator*(v,value);
     }
     
     //! Overload / operator
     friend Vec3 operator/(const Vec3 &v, const double value) {
          return Vec3(v.data[0] / value,
                      v.data[1] / value,
                      v.data[2] / value);
     }

     //! Overload + operator
     friend Vec3 operator+(const Vec3 &v, const double value) {
          return Vec3(v.data[0] + value,
                      v.data[1] + value,
                      v.data[2] + value);
     }

     //! Overload + operator
     friend Vec3 operator+(const double value, const Vec3 &v) {
          return operator+(v,value);
     }

     //! Overload + operator
     friend Vec3 operator+(const Vec3 &v1, const Vec3 &v2) {
          return Vec3(v1.data[0] + v2.data[0],
                      v1.data[1] + v2.data[1],
                      v1.data[2] + v2.data[2]);
     }     

     //! Overload - operator
     friend Vec3 operator-(const Vec3 &v1, const Vec3 &v2) {
          return Vec3(v1.data[0] - v2.data[0],
                      v1.data[1] - v2.data[1],
                      v1.data[2] - v2.data[2]);
     }     

     //! Cross product
     friend Vec3 cross_product(const Vec3 &v1, const Vec3 &v2) {
          return Vec3((v1.data[1]*v2.data[2])-(v1.data[2]*v2.data[1]),
                      (v1.data[2]*v2.data[0])-(v1.data[0]*v2.data[2]),
                      (v1.data[0]*v2.data[1])-(v1.data[1]*v2.data[0]));
     }

     //! Overload output operator
     friend std::ostream& operator<<(std::ostream& out, const Vec3 &v) {
          out << "[" << v.data[0] << ", " << v.data[1]<< ", " << v.data[2] << "]";
          return out;
     }
     
};

//! Gauss Integral Tuned class 
class Git {
private:

     //! Whether to run in smoothen_backbone mode
     const bool smoothen_backbone_mode;

     //! Pointer to reference (average) Gauss integral data
     std::vector<std::vector<double> > *reference_gauss_integrals;

     //! Maintain map of reference gauss integral filenames - and their corresponding
     //! parameters. This ensures that Git objects can be constructed quickly also
     //! when a non-default parameter file is used.
     static std::map<std::string, std::vector<std::vector<double> > > reference_gauss_integrals_cache;

     //! Parsed version of default reference gauss integrals parameters - smoothen_backbone case.
     //! Note that this could be part of reference_gauss_integrals_cache, but is included in a separate
     //! variable to avoid the map lookip
     static std::vector<std::vector<double> > git_default_reference_gauss_integrals_smoothen;

     //! Parsed version of default reference gauss integrals parameters - non-smoothen_backbone case.
     //! Note that this could be part of reference_gauss_integrals_cache, but is included in a separate
     //! variable to avoid the map lookip
     static std::vector<std::vector<double> > git_default_reference_gauss_integrals;


     // Initialize reference Gauss Integral data
     void init_reference_gauss_integrals(std::string parameter_filename="");

     //! Draw a uniformly distributed random value
     //! Note: this function can be overwritten in a derived class
     //! if you wish to use a different random number generator.
     virtual double sample_uniform_01() {
          return rand()/(0.5*RAND_MAX)-1.0;
     }

     //! Make pertubations to the backbone coordinates
     void perturb_backbone(std::vector<Vec3> &polyl, double ampitude);

     //! Smoothens the backbone if the deformation of the backbone 
     //! can be done without selfintersections
     void smoothen_backbone(std::vector<Vec3> &polyl);

public:

     //! Constructor
     //!
     //! \param parameter_filename Filename containing average gauss
     //! integrals parameters. If not specified, a meaningful default
     //! is used.
     //! \param smoothen_backbone_mode Whether to smoothen the backbone
     Git(const std::string parameter_filename="", const bool smoothen_backbone_mode=1);

     //! Generate Gauss Integrals
     //!
     //! \param polyl Coordinates of protein chain
     //! \return 30 dimensional GIT vector
     std::vector<double> generate_gauss_integrals(std::vector<Vec3> &polyl);

     //! Get the maximum length (number of residues) supported by this version of GIT
     int get_max_length() const;

     //! Retrieve the current smoothen_backbone_mode setting
     bool get_smoothen_backbone_mode() const;

};

}


#endif
