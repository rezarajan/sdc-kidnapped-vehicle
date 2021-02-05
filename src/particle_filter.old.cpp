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

  std::default_random_engine gen;
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
  is_initialized = true;

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
  std::default_random_engine gen; // Random engine for sampling
  for(int i = 0; i < num_particles; ++i){
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    // std::cout << "Actual X: " << x << "\n";
    // std::cout << "Actual Y: " << y << "\n";
    // std::cout << "Actual Theta: " << theta << "\n";
    // std::cout << "Velocity: " << velocity << "\n";
    // std::cout << "Yaw Rate: " << yaw_rate << "\n";

  // No division by zero
    if (yaw_rate > 0.0001) {
        x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
        y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
        theta += yaw_rate * delta_t;
    } else {
        x += velocity * cos(particles[i].theta) * delta_t;
        y += velocity * sin(particles[i].theta) * delta_t;
        theta += yaw_rate * delta_t;
    }


    // Creating Normal (Gaussian) distributions
    normal_distribution<double> dist_x(x, std_pos[0]);
    normal_distribution<double> dist_y(y, std_pos[1]);
    normal_distribution<double> dist_theta(theta, std_pos[2]);
    
    // Set the new particle predicted location
    // particles[i].x = dist_x(gen);
    // particles[i].y = dist_y(gen);
    // particles[i].theta = std::fmod(dist_theta(gen), 2*M_PI);
    particles[i].x = x;
    particles[i].y = y;
    particles[i].theta = std::fmod(theta, 2*M_PI);;
    // std::cout << "X: " << particles[i].x << "\n";
    // std::cout << "Y: " << particles[i].y << "\n";
    // std::cout << "Theta: " << particles[i].theta << "\n";
    // std::cout << "Weight: " << particles[i].weight << "\n";
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

  for(auto & o: observations){
    bool first_check = true;
    double min_dist = 0;
    for(auto const & p: predicted){
      // distance from landmark (predicted) to sensor reading (observations)
      double current_dist = dist(p.x, p.y, o.x, o.y);
      if(!first_check){
        if(current_dist <= min_dist){
          min_dist = current_dist;
          o.id = p.id;
        }
      }
      else{
        min_dist = current_dist;
        o.id = p.id;
        first_check = false;
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
  
  // Find landmark associations for each particle 
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
      double theta = particles[i].theta;
      if(particles[i].theta < 0){
        theta = particles[i].theta + 2*M_PI; // -90 degrees
      }
      // double theta = 0; // -90 degrees
      particleObs[j].x = particles[i].x + (cos(theta) * particleObs[j].x) - (sin(theta) * particleObs[j].y);
      particleObs[j].y = particles[i].y + (sin(theta) * particleObs[j].x) + (cos(theta) * particleObs[j].y);
 
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

    // Update particle weights
  for(int i=0; i<num_particles; ++i){
    particles[i].weight = 1;
    for(int j=0; j<particles[i].associations.size(); ++j){
      int id = particles[i].associations[j]-1; // Map landmark id is the index + 1
      double landmark_x, landmark_y;
      landmark_x = map_landmarks.landmark_list[id].x_f;
      landmark_y = map_landmarks.landmark_list[id].y_f;
      // std::cout << "Landmark Coordinates: (" << landmark_x << "," << landmark_y << ")\n";
      // std::cout << "Observation Coordinates: (" << particles[i].sense_x[j] << "," << particles[i].sense_y[j] << ")\n";

      // particles[i].weight *= multiv_prob(std_landmark[0], std_landmark[1], particles[i].sense_x[j], particles[i].sense_y[j],
      //                                   landmark_x, landmark_y);
      // std::cout << "New Weight for Particle[" << i << "]: " << multiv_prob(std_landmark[0], std_landmark[1], particles[i].sense_x[j], particles[i].sense_y[j],
      //                                   landmark_x, landmark_y) << "\n";
    }
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  // Create a sampling distribution from the particle weights
  vector<double> weights;
  for (auto const & p: particles){
    weights.push_back(p.weight);
  }

  std::default_random_engine gen; // Random engine for sampling
  std::discrete_distribution<> dist_weights(weights.begin(), weights.end());
  
  // Resample
  std::vector<Particle> new_particles;
  for(int i=0; i<num_particles; ++i){
      new_particles.push_back(particles[dist_weights(gen)]);
  }
  // Update new particle set
  particles = new_particles;

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