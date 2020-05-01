#include <accelmath.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm> // std::copy
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

//#include <curand.h>

using namespace std;

#define threshold 0.1
// #pragma acc routine(sqrt) seq

double** kmeans(double **x, double **initial_centroids, int num_samples, int num_features, int k);
void writeCentroidsToFile(double **final_centroid, int k, int num_features);
void writeLabelsToFile(double **x, double *labels, int num_samples, int num_features);

enum Color { red,
              green,
              blue,
              yellow,
              purple,
              orange,
              pink,
              brown 
            };

int main()
{
  clock_t tStart = clock();
  srand(2); // seed 1
  setlocale(LC_ALL, "en_US.UTF-8");
  // Synthetic data
  int k = 8;              // Number of clusters
  int num_features = 10;   // x1, x2
  int num_samples = 6500; // total number of data points
  double **x = new double *[num_samples];

  for (int i = 0; i < num_samples; i++)
  {
    x[i] = new double[num_features];
  }

  std::string line;
  std::ifstream myfile("synthetic_dataset.txt");
  std::string delimiter = " ";
  int linenum = 0;
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      std::vector<std::string> result;
      std::istringstream iss(line);
      for (std::string s; iss >> s;)
        result.push_back(s);
      for (int i = 0; i < result.size(); i++)
      {
        stringstream num(result[i]);
        double x_temp = 0;
        num >> x_temp;
        x[linenum][i] = x_temp;
      }
      linenum++;
    }
  }
  myfile.close();

  for (int i = 0; i < num_samples; i++)
    for (int j = 2; j < num_features; j++)
    {
      x[i][j] = 1.0;
    }
  // Do some preprocessing
  // Initialize clusters randomly, but only within the min-max range
  double mins[num_features];
  double maxes[num_features];

  for (int j = 0; j < num_features; j++)
  {
    mins[0] = 0.0;
    maxes[0] = 0.0;
  }

  // Find min/max of each feature
  for (int i = 0; i < num_samples; i++)
  {
    for (int j = 0; j < num_features; j++)
    {
      if (x[i][j] < mins[j])
        mins[j] = x[i][j];
      if (x[i][j] > maxes[j])
        maxes[j] = x[i][j];
    }
  }
  
  double **initial_centroids = new double *[k];
  for (int i = 0; i < k; i++)
  {
    initial_centroids[i] = new double[num_features];
    for (int j = 0; j < num_features; j++)
    {
      if(j < 2){
        initial_centroids[i][j] = (mins[j] +
                                        rand()) % (int)maxes[j];
      }
      else{
        initial_centroids[i][j] = 0;
      }
       
      }
  }
  
  printf("Initial Centroids\n");
  for (int i = 0; i < k; i++)
  {
    printf("Centroid %d: (", i);
    for (int j = 0; j < num_features; j++)
    {
      printf("%f ", initial_centroids[i][j]);
    }
    printf(")\n");
  }

  int random[num_samples][num_features];
  for (int i = 0; i < num_samples; i++)
  {
    for(int j = 0; j <num_features; j++)
      random[i][j] = mins[j] +
                rand() % (int)maxes[j];
  }

  // double **final_centroid = kmeans(x, initial_centroids, num_samples, num_features, k);
  double **distances = new double *[num_samples];
  for (int i = 0; i < num_samples; i++)
  {
    distances[i] = new double[k];
    for (int j = 0; j < k; j++)
      distances[i][j] = INT_MAX;
  }
  double *clusters = new double[num_samples];
  for(int i = 0; i < num_samples; i++){
    clusters[i] = 0;
  }
  double **centroids = new double*[k];
  for(int i = 0; i < k ; i++){
    centroids[i] = new double[num_features];
    for(int j = 0; j < num_features;j++){
      if(j >= 2){
        centroids[i][j] =0;
      }
      else{
        centroids[i][j] = initial_centroids[i][j];
      }
    }
  }
  int rand_counter = 0;
  for(int i = 0; i < k; i++){
    for(int j = 0; j < 2; j++){
      cout<<centroids[i][j]<<"  ";
    }
    cout<<endl;
  }

  //#pragma acc data create(clusters[0:num_samples]) copy(centroids [0:k] [0:num_features]) copyin(x[0:num_samples][0:num_features])
  double thresh = 1;
  //#pragma acc routine seq
  while(thresh >= 0.1)
  {
    thresh = 0;
    int* clusters_old = new int[num_samples];
    for(int i = 0; i < num_samples; i++){
      clusters_old[i] = clusters[i];
    }
    #pragma acc loop independent
    for (int i = 0; i < num_samples; i++)
    {
      for (int j = 0; j < k; j++){
        distances[i][j] = 0.0;
      }
    }
    double closest_dist;
    //Loop through each sample
    //Loop through cluster for the sample and find closest centroid
    for (int i = 0; i < num_samples; i++)
    {
      closest_dist = INT_MAX;
      #pragma acc loop independent
      for (int c = 0; c < k; c++)
      {
        //Calculate l2 distance from each cluster
        //This is a data independet loop so I should be able to do a parallelization
        for (int j = 0; j < 2; j++)
        {
          distances[i][c] += (x[i][j] - centroids[c][j]) * (x[i][j] - centroids[c][j]);
        }
        // distances[i][c] = sqrt(distances[i][c]);
      
      closest_dist = distances[i][0];
      if (distances[i][c] < closest_dist)
      {
        //cout<<"Hello there"<<endl;
        closest_dist = distances[i][c];
        clusters[i] = c;
      }
      
      //cout<<"went through a loop fine"<<endl;
    }
    }
    double temp =0;
    for(int i =0; i < num_samples; i++){
      if(clusters_old[i] - clusters[i] != 0){
        temp+=1;
      }
    }
    //cout<<"This is temp "<< temp << endl;
    thresh = 1.0*temp/num_samples;

    // counts holds the number of data points currently in the cluster
    int *counts = new int[k];
    for (int c = 0; c < k; c++)
    {
      // new_centroids[c] = new double[num_features];
      counts[c] = 0.0;

      #pragma acc loop independent
      for (int i = 0; i < num_samples; i++)
      {
        if (clusters[i] == c)
        {
          counts[c] += 1;
          // #pragma acc parallel loop
          for (int j = 0; j < num_features; j++)
          {
            centroids[c][j] += x[i][j];
          }
        }
      }
    }
    #pragma acc loop independent
    for (int c = 0; c < k; c++)
    {
      // Divide by number of data points in cluster
      // This is the new centroid (average)
      for (int j = 0; j < num_features; j++)
      {
        if (counts[c] == 0)
          centroids[c][j] = random[rand_counter][j] % 500000; // If no data points in group, then reinitialize
        else
          centroids[c][j] = centroids[c][j] / counts[c];
      }
      rand_counter = (rand_counter * 7 + 31) % num_samples;
    }
    writeLabelsToFile(x, clusters, num_samples, num_features);
    cout << thresh << endl;
  }
  
  writeCentroidsToFile(centroids, k, num_features);
  printf("Time taken for clustering parallel : %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  writeLabelsToFile(x, clusters, num_samples, num_features);
  system("gnuplot -p plotCluster");
}

std::string getColor(int val){
  std::string color;
  switch (val)
  {
  case Color::red:
    color = "red";
    break;
  case Color::green:
    color = "green";
    break;
  case Color::blue:
    color = "blue";
    break;
  case Color::yellow:
    color = "yellow";
    break;
  case Color::purple:
    color = "purple";
    break;
  case Color::orange:
    color = "orange";
    break;
  case Color::pink:
    color = "pink";
    break;
  case Color::brown:
    color = "brown";
    break;
  }
  return color;
}

void writeCentroidsToFile(double **final_centroid, int k, int num_features)
{
  std::ofstream outfile;
  outfile.open("centroids.txt");
  printf("\n");
  printf("Final Centroids\n");
  std::string color;
  for (int i = 0; i < k; i++)
  {
    color = getColor(i);
    printf("Centroid %d: (", i);
    for (int j = 0; j < num_features; j++)
    {
      printf("%f ", final_centroid[i][j]);
      outfile << final_centroid[i][j] << " ";
    }
    printf(")\n");
    outfile << i << "\n";
  }
  outfile.close();
}

void writeLabelsToFile(double **x, double *labels, int num_samples, int num_features)
{
  std::string color;
  std::ofstream outfile;
  outfile.open("data_colors.txt");
  for (int i = 0; i < num_samples; i++)
  {
    color = getColor(labels[i]);
    for (int j = 0; j < num_features; j++)
    {
      outfile << x[i][j] << " ";
    }
    outfile << labels[i] << "\n";
  }
  outfile.close();

}
