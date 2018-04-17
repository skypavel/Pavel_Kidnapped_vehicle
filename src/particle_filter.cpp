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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    default_random_engine gen;

    //set number of particles
    num_particles = 50;

    //Resize vector of particles
    particles.resize(num_particles);

    //Resize vector of weights
    weights.resize(num_particles);

    //Create Normal distribution for x, y and theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    std::cout<<" gps_x = "<<x<<" gps_y = "<<y<<endl;
    //
    for (int i = 0; i < num_particles; ++i) {

        double sample_x = dist_x(gen);
        double sample_y = dist_y(gen);
        double sample_theta = dist_theta(gen);

        particles[i].id = i;
        particles[i].x = sample_x;
        particles[i].y = sample_y;
        particles[i].theta = sample_theta;
        particles[i].weight = 1.0;
        }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;

    //Create Normal distribution for noise
    normal_distribution<double> dist_x(0.0, std_pos[0]);
    normal_distribution<double> dist_y(0.0, std_pos[1]);
    normal_distribution<double> dist_theta(0.0, std_pos[2]);

    for (int i = 0; i < num_particles; i++){
        if (fabs(yaw_rate) < 0.0001) {
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
        }
        else {
            particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t))-sin(particles[i].theta));
            particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
            particles[i].theta += yaw_rate * delta_t;
        }

        //add noise to the particles
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);

      //  std::cout<<" x_pred = "<<particles[i].x<<" y_pred = "<<particles[i].y<<" v = "<< velocity<<endl;

    }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html





    double std_x = std_landmark[0];
    double std_y = std_landmark[1];
    double gauss_norm = (1/(2 * std_x * std_y));
    double exponent = 0.0;

    double weights_sum = 0.0;

    //Consider all particles
    for (int i = 0; i < num_particles; ++i) {

        double cur_w = 1.0;

        //Consider each observation
        for (int j = 0; j < observations.size(); ++j) {

            LandmarkObs current_obs = observations[j];
            LandmarkObs transformed_obs;

            //Transform the observation point from vehicle coordinates to mapp coordinates
            transformed_obs.x = observations[j].x * cos(particles[i].theta) -
                    observations[j].y * sin(particles[i].theta)+particles[i].x;
            transformed_obs.y = observations[j].x * sin(particles[i].theta) +
                    observations[j].y * cos(particles[i].theta)+particles[i].y;

            //Find all landmark which are within range of sensor
            vector<double> landmark_observation_dist (map_landmarks.landmark_list.size());

            double min_distance = 1000000.0;
            int min_lm = 0;

            //find the nearest landmark
            for (int k = 0; k < map_landmarks.landmark_list.size(); ++k) {
                //calculate distance between current particle and landmarks
                double landmark_particle_dist = sqrt(pow(particles[i].x - map_landmarks.landmark_list[k].x_f, 2) +
                                                     pow(particles[i].y - map_landmarks.landmark_list[k].y_f, 2));

                // take in account landmarks, which are located within sensor range
                if (landmark_particle_dist < sensor_range) {
                    landmark_observation_dist[k] = sqrt(pow(transformed_obs.x - map_landmarks.landmark_list[k].x_f, 2) +
                                                     pow(transformed_obs.y - map_landmarks.landmark_list[k].y_f, 2));
                    // find the nearest landmark
                    if (landmark_observation_dist[k] < min_distance){
                        min_distance = landmark_observation_dist[k];
                        min_lm = k;
                    }
                    } else {
                    landmark_observation_dist[k] = 1000000.0;
                }


            }


            double x_diff = transformed_obs.x - map_landmarks.landmark_list[min_lm].x_f;
            double y_diff = transformed_obs.y - map_landmarks.landmark_list[min_lm].y_f;
            exponent = 0.5 * ((x_diff*x_diff)/(std_x*std_x) +
                              (y_diff*y_diff)/(std_y*std_y));
            cur_w = exp(-exponent)*gauss_norm;
            particles[i].weight *= cur_w;

            }


        weights_sum += particles[i].weight;


        }
        for (int i = 0; i < num_particles; i++) {
            particles[i].weight /= weights_sum;
            weights[i] = particles[i].weight;
          }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> new_particles(num_particles);

    random_device rd;
    default_random_engine gen(rd());

    for (int i = 0; i < num_particles; ++i) {
        discrete_distribution<int> index(weights.begin(), weights.end());
        new_particles[i] = particles[index(gen)];
    }

    particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
