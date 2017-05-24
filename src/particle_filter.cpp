/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	std::default_random_engine gen_init;
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);
	particles.resize(num_particles);
	weights.resize(num_particles);
	for (int i=0; i<num_particles; i++){
    particles[i].id = i;
    particles[i].x = dist_x(gen_init);
    particles[i].y = dist_y(gen_init);
    particles[i].theta = dist_theta(gen_init);
    particles[i].weight = 1.0;
    weights[i] = 1.0;
	}
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  std::default_random_engine gen_pred;
  std::normal_distribution<double> dist_x(0,std_pos[0]);
  std::normal_distribution<double> dist_y(0,std_pos[1]);
  std::normal_distribution<double> dist_theta(0,std_pos[2]);
  for (int i=0; i<num_particles; i++){
    //Updating theta
    double init_theta = particles[i].theta;
    particles[i].theta += yaw_rate*delta_t;
    // updating x and y on the basis of yaw rate
    if (yaw_rate > 1e-3){
      particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta) - sin(init_theta));
      particles[i].y += (velocity/yaw_rate)*(-cos(particles[i].theta) + cos(init_theta));
    }
    else{
      particles[i].x += velocity*cos(init_theta)*delta_t;
      particles[i].y += velocity*sin(init_theta)*delta_t;
    }
    // Adding Noise
    particles[i].x += dist_x(gen_pred);
    particles[i].y += dist_y(gen_pred);
    particles[i].theta += dist_theta(gen_pred);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
  //float dist=INFINITY;
  std::vector<double> euDist;
  if(!observations.empty() && !predicted.empty()){
    for (size_t i=0; i<observations.size(); i++){
      for (size_t j=0; j<predicted.size(); j++){
        euDist.push_back(dist(observations[i].x,observations[i].y, predicted[j].x, predicted[j].y));
      }
      std::vector<double>::iterator min_val_ptr = std::min_element(std::begin(euDist), std::end(euDist));
      int ind_lm = std::distance(std::begin(euDist),min_val_ptr);
      // For DEBUGGING ONLY ->
      //for (auto elem : euDist){
      //  std::cout<< elem<<" ";
      //}
      //std::cout<<"\n";
      //std::cout<<"The minimum index is "<< ind_lm<< std::endl;
      observations[i].id = predicted[ind_lm].id;
      euDist.erase(euDist.begin(),euDist.end());
    }
  }
  else{
    std::cout<<"Data Association cannot be performed if one of the vectors is empty!\n";
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	//   std::cout <<"Observation vector size="<< observations.size()<<std::endl;

	std::vector <LandmarkObs> transObs;
	std::vector<LandmarkObs> predLandmarks;
	// Browse through every particle
	for (int i=0; i<num_particles; i++){
    double theta = particles[i].theta;
    double x = particles[i].x;
    double y = particles[i].y;
    // Create a list of Landmarks within the range of the particle
    for (size_t k=0; k<map_landmarks.landmark_list.size(); k++){
      LandmarkObs lm;
      Map::single_landmark_s lm_map = map_landmarks.landmark_list[k];
      lm.id = lm_map.id_i;
      lm.x = (double) lm_map.x_f;
      lm.y = (double) lm_map.y_f;
      if (dist(x,y,lm.x,lm.y) <= sensor_range){
        predLandmarks.push_back(lm);
      }
    }
    // Transform the observation coordinates to Map coordinates
    for (size_t j=0; j<observations.size(); j++){
      LandmarkObs t_obs;
      double obs_x = observations[j].x;
      double obs_y = observations[j].y;
      double cos_theta = cos(theta);
      double sin_theta = sin(theta);
      t_obs.id = observations[j].id;
      t_obs.x = x + cos_theta*obs_x - sin_theta*obs_y;
      t_obs.y = y + cos_theta*obs_y + sin_theta*obs_x;
      transObs.push_back(t_obs);
    }
    // Associate Actual Landmarks with Transformed Observations.
    dataAssociation(predLandmarks,transObs);
    //for(auto& obs:transObs){
    //  obs.id = predLandmarks[0].id;
    //}
    // Update weights based on Landmark and observation distances.
    double wt = 1;
    while (!transObs.empty()){
      LandmarkObs lastObs = transObs.back();
      transObs.pop_back();
      for (size_t m=0; m<predLandmarks.size(); m++){
        if (lastObs.id == predLandmarks[m].id){
          wt *= gaussian_prob(lastObs, predLandmarks[m], std_landmark);
          break;
        }
      }
    }
    particles[i].weight = wt;
    weights[i] = wt;
    // Empty the predLandmark vector
    predLandmarks.clear();
    // Check if vectors are empty
    if (!transObs.empty() || !predLandmarks.empty()){
      std::cout<< "The vectors are not empty!!"<<std::endl;
    }
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	//Setup random bits
	//std::random_device rd;
	//std::mt19937 gen(rd());
  std::default_random_engine gen;
  //std::default_random_engine gen;
  std::discrete_distribution<> discDist(weights.begin(),weights.end());
  std::vector<Particle> particlesResampled;
  particlesResampled.resize(num_particles);
  for (int i=0; i<num_particles; i++){
    int index = discDist(gen);
    particlesResampled[i] = particles[index];
  }
  particles = particlesResampled;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
