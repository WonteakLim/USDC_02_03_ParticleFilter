/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */
#include <cstdlib>
#include <ctime>
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

const int NUMBER_OF_PARTICLES = 50; //300
const double INITIAL_WEIGHT = 1.0;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    this->num_particles = NUMBER_OF_PARTICLES;

    default_random_engine gen;
    normal_distribution<double> particle_x(x, std[0]);
    normal_distribution<double> particle_y(y, std[1]);
    normal_distribution<double> particle_theta(theta, std[2]);

        for (int i = 0; i < NUMBER_OF_PARTICLES; ++i) {

            Particle p;
            p.x = particle_x(gen);
            p.y = particle_y(gen);
            p.theta = particle_theta(gen);
            p.weight = 1;
            weights.push_back(p.weight);
            this->particles.push_back(p);

        }
    this->is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

      default_random_engine gen;
      normal_distribution<double> noise_x(0, std_pos[0]);
      normal_distribution<double> noise_y(0, std_pos[1]);
      normal_distribution<double> noise_theta(0, std_pos[2]);

      for (int i = 0;  i < NUMBER_OF_PARTICLES; ++i) {

      this->particles[i].x += velocity * delta_t * cos(this->particles[i].theta)
                                                            + noise_x(gen);
      this->particles[i].y += velocity * delta_t * sin(this->particles[i].theta)
                                                            + noise_y(gen);
      this->particles[i].theta += yaw_rate * delta_t + noise_theta(gen);
        }
}

std::vector<LandmarkObs> ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
        // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
        //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.


        std::vector<LandmarkObs> association;
    LandmarkObs NearLandmark;

    for (int i = 0; i < observations.size(); ++i){
        LandmarkObs obs = observations[i];
        double Compare = 10000.0;

        for (int j = 0; j < predicted.size(); ++j){
            LandmarkObs pred = predicted[j];
            double distance = dist(obs.x, obs.y, pred.x, pred.y);

            if (distance < Compare){
                Compare = distance;
                NearLandmark = pred;
            }
        }
        // printf("Distance = %f\n",Compare);
        association.push_back(NearLandmark);
        // printf("association %d = %f\n",i,association[i].x);
        // printf("transform %d = %f\n",i,observations[i].x);
    }


    return association;
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

    double sum = 0;
    std::vector<double> likelihood;
    for (int i = 0; i < NUMBER_OF_PARTICLES; ++i){
        Particle p = this->particles[i];
        std::vector<LandmarkObs> predicted_landmark;
        std::vector<LandmarkObs> transformed_observations;
        for (int j = 0; j < observations.size(); ++j){
            LandmarkObs observation = observations[j];
            LandmarkObs transformed_observation;
            transformed_observation.x = p.x + observation.x * cos(p.theta) - observation.y * sin(p.theta);
            transformed_observation.y = p.y + observation.x * sin(p.theta) + observation.y * cos(p.theta);
            transformed_observations.push_back(transformed_observation);
            // printf("%f\n",transformed_observation.x);
        }

        // printf("x_p = %f\n",p.x);
        // printf("y_p = %f\n",p.y);

        for (auto landmark: map_landmarks.landmark_list){

            double distance = dist(p.x,p.y,landmark.x_f,landmark.y_f);
            if (distance < sensor_range) {
                LandmarkObs Add_landmark;
                Add_landmark.id = landmark.id_i;
                Add_landmark.x = landmark.x_f;
                Add_landmark.y = landmark.y_f;
                predicted_landmark.push_back(Add_landmark);
            }
        }
        // printf("old = %f\n",predicted_landmark[0].x);
        std::vector<LandmarkObs> associated_landmark;
        associated_landmark = dataAssociation(predicted_landmark, transformed_observations);

        // printf("new = %f\n",associated_landmark[0].x);
        double probability = 1;
        double x;
        for (int j=0; j < associated_landmark.size(); ++j){
            // printf("predicted = %f\n",associated_landmark[1].x);
            // printf("transform = %f\n",transformed_observations[1].x);
            // printf("transformed = %f\n",transformed_observations[j].x);
            double dx = associated_landmark[j].x - transformed_observations[j].x;
            // printf("dx = %f\n",dx);
            double dy = associated_landmark[j].y - transformed_observations[j].y;
            probability *= 1.0/(2*M_PI*std_landmark[0]*std_landmark[1]) * exp(-dx*dx / (2*std_landmark[0]*std_landmark[0])/1000)* exp(-dy*dy / (2*std_landmark[1]*std_landmark[1])/1000);
            // probability *= 1/(dx*dy);
            // x = -dx*dx / (2*std_landmark[0]*std_landmark[0]);
            // printf("%f\n",x);
            // printf("%f\n",probability);
            // double inv = 1/probability;
            // printf("value = %f\n",inv);
        }
        weights.push_back(probability);
        this->particles[i].weight = probability;
				sum += this->particles[i].weight;

    }

    // for (int i = 0; i < weights.size(); i++)
    // {
    //     sum += weights.at(i);
    // }
  //
    // for (int i = 0; i < weights.size(); i++){
    //     weights.at(i) /= sum;
    //     // this->particles[i].weight /= sum;
    //     // printf("%f\n",this->particles[i].weight);
    //     // x += this->particles[i].weight;
    //     // printf("sum = %f\n",x);
    // }
  //

    // printf("sum = %f\n",sum);
    double x = 0;
    for (int i = 0; i < NUMBER_OF_PARTICLES; ++i){
        // weights.at(i) /= sum;
        this->particles[i].weight /= sum;
        // printf("%f\n",this->particles[i].weight);
        // x += this->particles[i].weight;
        // printf("sum = %f\n",x);
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    // double highest_weight = -1.0;
    // Particle best_particle;

    // for (int i = 0; i < num_particles; ++i) {
    //     if (particles[i].weight > highest_weight) {
    //         highest_weight = particles[i].weight;
    //         best_particle = particles[i];
    //     }
    // }
  //
    // for (int i = 0; i < num_particles; ++i) {
    //     particles[i] = best_particle;
    // }

    default_random_engine gen;
    std::vector<Particle> weighted_sample(num_particles);
    std::uniform_real_distribution<float> distribution(0.0, 1.0);

    for(int i = 0; i < num_particles; ++i){
        float myrand = distribution(gen);
        double weightsum = 0;
        for(int j = 0; j < num_particles; ++j){
            weightsum += this->particles[j].weight;

            if(weightsum > myrand){
                weighted_sample.at(i) = this->particles[j];
                break;
            }
        }
    }
    this->particles = weighted_sample;



// std::discrete_distribution<int> d(weights.begin(), weights.end());
// std::vector<Particle> weighted_sample(num_particles);
// default_random_engine gen;
// for(int i = 0; i < num_particles; i++){
//     int j = d(gen);
//     weighted_sample.at(i) = particles.at(j);
// }
//
// particles = weighted_sample;

return;

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
