/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
#include "multiv_gauss.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles

  std::default_random_engine gen; // Random engine for sampling
  // Creating Normal (Gaussian) distributions
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i) {
    Particle particle;

    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1;

    particles.push_back(particle);

  }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  for(int i = 0; i < num_particles; ++i){
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    // Assuming bicycle model
    if(yaw_rate == 0){
      x +=velocity*delta_t*cos(theta); 
      y +=velocity*delta_t*sin(theta); 
    }
    else{
      x += (velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
      y += (velocity/yaw_rate)*(-cos(particles[i].theta+yaw_rate*delta_t)+cos(particles[i].theta));
      theta += yaw_rate*delta_t;
    }
    std::default_random_engine gen; // Random engine for sampling
    // Creating Normal (Gaussian) distributions
    normal_distribution<double> dist_x(x, std_pos[0]);
    normal_distribution<double> dist_y(y, std_pos[1]);
    normal_distribution<double> dist_theta(theta, std_pos[2]);
    
    // Set the new particle predicted location
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  // Given all landmarks within the sensor range of a particle, find the nearest ones to
  // the sensor measurements (using nearest-neighbors)
  // NOTE: these observation associations are stored in the Particle `particle[i].association` vector
  // use setAssociations method for this
  // NOTE: observations do not contain an ID, since these are only measurements to unknown landmarks

  for(int i=0; i<predicted.size(); ++i){
    double min_dist = 0; 
    for(int j=0; j<observations.size(); ++j){
      double current_dist = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
      // Update associations based on nearest neighbour
      if(current_dist <= min_dist){
        min_dist = current_dist;
        observations[j].id = predicted[i].id;
      }
    }
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  // `observations` are the same for all particles since these are the actual sensor measurements
  // create a vector called `predicted` to store all landmarks which are within sensor range for each particle
  // use dataAssociation to perform a nearest-neighbors search
  for(int i=0; i<num_particles; ++i){
    float x_max, x_min;
    float y_max, y_min;

    vector<LandmarkObs> predicted;
    // Checking for landmarks within sensor range of each particle
    x_max = particles[i].x + sensor_range;
    y_max = particles[i].y + sensor_range;
    x_min = particles[i].x - sensor_range;
    y_min = particles[i].y - sensor_range;

    for(int j=0; j<map_landmarks.landmark_list.size(); ++j){
      if((x_min <= map_landmarks.landmark_list[j].x_f <= x_max) && 
          (y_min <= map_landmarks.landmark_list[j].y_f <= y_max)){
            LandmarkObs landmark;
            landmark.id = map_landmarks.landmark_list[j].id_i;
            landmark.x = map_landmarks.landmark_list[j].x_f;
            landmark.y = map_landmarks.landmark_list[j].y_f;
            predicted.push_back(landmark);
          }
    }
    // Copy observations vector for mutability in the dataAssociation method
    // for each partile
    vector<LandmarkObs> particleObs = observations;
    for(int j=0; j<particleObs.size(); ++j){
      // transform observations to map coordinates
      double x_map, y_map;
      particleObs[j].x = particles[i].x + (cos(particles[i].theta) * particleObs[j].x) - (sin(particles[i].theta) * particleObs[j].y);
      particleObs[j].y = particles[i].y + (sin(particles[i].theta) * particleObs[j].x) + (cos(particles[i].theta) * particleObs[j].y);

    }
    // Find observation associations
    dataAssociation(predicted, particleObs);

    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;
    for(int j=0; j<particleObs.size(); ++j){
      associations.push_back(particleObs[j].id);
      sense_x.push_back(particleObs[j].x);
      sense_y.push_back(particleObs[j].y);
    }
    // Update particle landmark associations
    SetAssociations(particles[i], associations, sense_x, sense_y);
    
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}