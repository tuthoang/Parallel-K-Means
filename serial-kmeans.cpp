#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm> // std::copy
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <iostream>
using namespace std;

double **kmeans(double **x, double **initial_centroids, int num_samples, int num_features, int k, double *mins, double *maxes);
int *getLabels(double **x, double **centroids, int num_samples, int num_features, int k);
double **getCentroids(double **x, double **centroids, int *clusters, int num_samples, int num_features, int k, double *mins, double *maxes);
void writeCentroidsToFile(string centroids_file, double **final_centroid, int k, int num_features);
void writeLabelsToFile(string filename, double **x, int *labels, int num_samples, int num_features);

enum Color
{
  red,
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
  srand(157); // seed 1
  // srand(time(NULL));
  setlocale(LC_ALL, "en_US.UTF-8");
  // Synthetic data
  int k = 20;              // Number of clusters
  int num_features = 2;   // x1, x2
  int num_samples = 3000; // total number of data points
  double **x = new double *[num_samples];

  for (int i = 0; i < num_samples; i++)
  {
    x[i] = new double[num_features];
  }

  string filename = "a1.txt";
  string ground_truth_filename = "a1-ga-cb.txt";

  std::string line;
  std::ifstream myfile(filename);
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
        x[linenum][i] = std::stod(result[i]);
      }
      linenum++;
    }
  }
  myfile.close();

  double **centroids = new double*[k];
  double **min_centroids = new double*[k];
  double **initial_centroids = new double*[k];
  int *clusters = new int[k];
  double **distances = new double *[num_samples];
  double **old_centroids = new double *[k];
  double *mins = new double[num_features];
  double *maxes = new double[num_features];
  int *counts = new int[k];
  // Initialize clusters randomly, but only within the min-max range
  mins = new double[num_features];
  mins[0] = x[0][0];
  mins[1] = x[0][1];
  maxes = new double[num_features];
  maxes[0] = x[0][0];
  maxes[1] = x[0][1];

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
  cout << maxes[0] << "\t" << maxes[1] << endl;

  for (int loop = 0; loop < 10000; loop++)
  {
    // Do some preprocessing


    // initial_centroids = new double *[k];
    // centroids = new double *[k];
    // min_centroids = new double *[k];
    
    for (int i = 0; i < k; i++)
    {
      initial_centroids[i] = new double[num_features];
      centroids[i] = new double[num_features];
      min_centroids[i] = new double[num_features];
      for (int j = 0; j < num_features; j++)
      {
        initial_centroids[i][j] = mins[j] +
                                  rand() % (int)maxes[j];
        centroids[i][j] = initial_centroids[i][j];
        min_centroids[i][j] = initial_centroids[i][j];
      }
    }
    // centroids = initial_centroids;
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

      distances = new double *[num_samples];
      for (int i = 0; i < num_samples; i++)
      {
        distances[i] = new double[k];
        for (int j = 0; j < k; j++)
          distances[i][j] = 0.0;
      }

      clusters = new int[num_samples];
      for (int i = 0; i < num_samples; i++)
      {
        clusters[i] = -1;
      }

      old_centroids = new double *[k];

      for (int i = 0; i < k; i++)
      {
        old_centroids[i] = new double[num_features];
        for (int j = 0; j < num_features; j++)
        {
          old_centroids[i][j] = centroids[i][j];
        }
      }
      // double **final_centroid = kmeans(x, initial_centroids, num_samples, num_features, k, mins, maxes);

      double thresh = 0.05;

      int thresh_met_counter = 0;
      bool thresh_met_max = false;

      int min_WCSS = INT_MAX;

      int count = 0;
      
      while (thresh_met_counter < (k * num_features) && count < 10000)
      // while(count<1000)
      {
          thresh_met_counter = 0;

          //Loop through each sample
          //Loop through cluster for the sample and find closest centroid
          double closest_dist;
          for (int i = 0; i < num_samples; i++)
          {
            closest_dist = INT_MAX;
            //#pragma acc loop independent
            for (int c = 0; c < k; c++)
            {
              //Calculate l2 distance from each cluster
              //This is a data independet loop so I should be able to do a parallelization
              for (int j = 0; j < num_features; j++)
              {
                distances[i][c] += ((x[i][j] - centroids[c][j]) * (x[i][j] - centroids[c][j]));
              }
              distances[i][c] = sqrt(distances[i][c]);
              // cout << distances[i][c] <<endl;
              // closest_dist = distances[i][0];
              if (distances[i][c] < closest_dist)
              {
                // cout<< c << endl;
                closest_dist = distances[i][c];
                clusters[i] = c;
              }
            }
          }


          // counts holds the number of data points currently in the cluster
          counts = new int[k];
          for (int c = 0; c < k; c++)
          {
            counts[c] = 0;
            for (int i = 0; i < num_samples; i++)
            {
              if (clusters[i] == c)
              {
                counts[c] += 1;
                for (int j = 0; j < num_features; j++)
                {
                  centroids[c][j] += x[i][j];
                }
              }
            }
          }

          for (int c = 0; c < k; c++)
          {
            // Divide by number of data points in cluster
            // This is the new centroid (average)
            for (int j = 0; j < num_features; j++)
            {
              if (counts[c] == 0)
              {
                // cout<< c << "  Has no data points.  " << mins[j] << " \t" << maxes[j]<<endl;
                centroids[c][j] = mins[j] + rand() % (int)maxes[j]; // If no data points in group, then reinitialize
              }
              else
                centroids[c][j] = centroids[c][j] / counts[c];
            }
          }

          for (int i = 0; i < k; i++)
          {
            for (int j = 0; j < num_features; j++)
            {
              // cout << abs(centroids[i][j] - old_centroids[i][j]) <<endl;
              if (abs(centroids[i][j] - old_centroids[i][j]) / abs(centroids[i][j]) < thresh)
              {
                thresh_met_counter++;
              }
            }
          }

          for (int i = 0; i < k; i++)
          {
            for (int j = 0; j < num_features; j++)
            {
              old_centroids[i][j] = centroids[i][j];
            }
          }


          count++;
          //writeLabelsToFile(x, clusters, num_samples, num_features);
          cout << "This is the iter " << count << " this is the number of points that changed centroids " << thresh_met_counter << endl;
      }

      double t1 = (double)(clock() - tStart) / CLOCKS_PER_SEC;
      printf("Time taken for clustering serially : %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

      
      //WITH IN SUM OF SQUARES
      double wcss = 0;
      double wcss_cluster = 0;
      for(int c = 0; c < k; c++){
        wcss_cluster = 0;
        for(int i = 0; i < num_samples; i++){
          if(clusters[i] == c){
            for(int j = 0; j < num_features; j++){
              wcss_cluster = wcss_cluster + (x[i][j] - centroids[c][j]) * (x[i][j] - centroids[c][j]);
            }
          }
        }
        wcss += sqrt(wcss_cluster);
      }
      // cout <<wcss<<endl;
      wcss = wcss/num_samples;
      // cout << wcss << endl;

      if(wcss < min_WCSS){
        min_WCSS = wcss;
        for(int c = 0; c < k; c++){
            for(int j = 0; j < num_features; j++){
              min_centroids[c][j] = centroids[c][j];
            }
        }
      }
      cout << "end" << endl;
  }

  printf("InitialCentroids\n");
  for (int i = 0; i < k; i++)
  {
    printf("Centroid %d: (", i);
    for (int j = 0; j < num_features; j++)
    {
      printf("%f ", initial_centroids[i][j]);
    }
    printf(")\n");
  }

  string out = "data_colors.txt";
  string centroids_file = "centroids.txt";
  writeCentroidsToFile(centroids_file, min_centroids, k, num_features);
  int *labels = getLabels(x, min_centroids, num_samples, num_features, k);
  writeLabelsToFile(out, x, labels, num_samples, num_features);
  std::stringstream fmt;
  fmt << "gnuplot -p -e \"k=" << k << "; data_labels= '" << out
      << "'; centroids_file= '" << centroids_file
      << "'; Title='K-Means'\" plotCluster.gp";
  system(fmt.str().c_str());


  double **ground_truth = new double*[k];
  line;
  std::ifstream myfile1(ground_truth_filename);
  delimiter = " ";
  linenum = 0;
  if (myfile1.is_open())
  {
    while (getline(myfile1, line))
    {
      std::vector<std::string> result;
      std::istringstream iss(line);
      for (std::string s; iss >> s;)
        result.push_back(s);
      ground_truth[linenum] = new double[num_features];
      for (int i = 0; i < result.size(); i++)
      {
        ground_truth[linenum][i] = std::stod(result[i]);
      }
      linenum++;
    }
  }
  myfile1.close();
  printf("--------------------\n");

  for(int i =0 ; i<k; i++){
    for(int j = 0; j < num_features; j++){
      cout<< ground_truth[i][j] << "  ";
    }
    cout<<endl;
  }
  printf("asdasdsaasd\n");
  string ground_truth_labels_file = "ground_truth_centroids.txt";
  labels = getLabels(x, ground_truth, num_samples, num_features, k);
  writeLabelsToFile(ground_truth_labels_file, x, labels, num_samples, num_features);
  printf("2132132121asdasdsaasd\n");

  std::stringstream fmt1;
  fmt1 << "gnuplot -p -e \"k=" << k << "; data_labels= '" << ground_truth_labels_file 
        << "'; centroids_file= '" << ground_truth_filename 
        << "'; Title='Ground Truth' \" plotCluster.gp";
  system(fmt1.str().c_str());
}

double **kmeans(double **x, double **centroids,
                int num_samples, int num_features, int k,
                double *mins, double *maxes)
{

  double **old_centroids = new double *[k];

  for (int i = 0; i < k; i++){
    old_centroids[i] = new double[num_features];
    for (int j = 0; j < num_features; j++){
      old_centroids[i][j] = centroids[i][j];
    }
  }

  int *clusters;

  int count = 0;
  double thresh = 0.05;


  int thresh_met_counter = 0;
  bool thresh_met_max = false;


  while (thresh_met_counter < (k * num_features) && count < 10000)
  // while(count < 10000)
  {
    // std::copy(&centroids[0][0], // Used for convergence checking
    //           &centroids[0][0] + num_samples * num_features,
    //           old_centroids[0][0]);
    // printf("\n");
    // for (int i = 0; i < k; i++)
    // {
    //   printf("Centroid %d: (", i);
    //   for (int j = 0; j < num_features; j++)
    //   {
    //     printf("%f ", centroids[i][j]);
    //   }
    //   printf(")\n");
    // }

    thresh_met_counter = 0;
    clusters = getLabels(x, centroids, num_samples, num_features, k);
    centroids = getCentroids(x, centroids, clusters, num_samples, num_features, k, mins, maxes);

    for (int i = 0; i < k; i++)
    {
      for (int j = 0; j < num_features; j++)
      {
        // cout << abs(centroids[i][j] - old_centroids[i][j]) <<endl;
        if (abs(centroids[i][j] - old_centroids[i][j]) / abs(centroids[i][j]) < thresh)
        {
          thresh_met_counter++;
        }
      }
    }

    for (int i = 0; i < k; i++)
    {
      for (int j = 0; j < num_features; j++)
      {
        old_centroids[i][j] = centroids[i][j];
      }
    }

    count++;
    cout << count << "    " << thresh_met_counter << endl;
  }

  return centroids;
}

/**
 * Assigns clusters on the given data by
 * calculating the closest distance to current centroids
 */
int *getLabels(double **x, double **centroids,
               int num_samples, int num_features, int k)
{
  int *clusters = new int[num_samples];
  double l2_dist;
  double closest_dist;
  // Loop through each sample
  // Loop through cluster for the sample and find closest centroid
  for (int i = 0; i < num_samples; i++)
  {
    closest_dist = INT_MAX;
    for (int c = 0; c < k; c++)
    {

      // Calculate l2 distance from each cluster
      l2_dist = 0.0;
      for (int j = 0; j < num_features; j++)
      {
        l2_dist += pow(x[i][j] - centroids[c][j], 2);
      }
      l2_dist = sqrt(l2_dist);

      // Assign closest centroid to data point
      if (l2_dist < closest_dist)
      {
        closest_dist = l2_dist;
        clusters[i] = c;
      }
    }
  }
  return clusters;
}

/**
 * Updates the centroids by calculating
 * the mean of the data points belonging to that cluster
 */
double **getCentroids(double **x, double **centroids, int *clusters,
                      int num_samples, int num_features, int k, 
                      double *mins, double *maxes)
{
  //Initializing everything to 0
  double **new_centroids = new double *[k];
  for (int i = 0; i < k; i++)
  {
    new_centroids[i] = new double[num_features];
    for (int j = 0; j < num_features; j++)
    {
      new_centroids[i][j] = 0;
      //cout<<centroids[i][j]<<endl;
    }
  }

  // counts holds the number of data points currently in the cluster
  int *counts = new int[k];

  for (int c = 0; c < k; c++)
  {
    // new_centroids[c] = new double[num_features];
    counts[c] = 0.0;
    for (int i = 0; i < num_samples; i++)
    {
      if (clusters[i] == c)
      {
        counts[c] += 1;
        for (int j = 0; j < num_features; j++)
        {
          new_centroids[c][j] += x[i][j];
        }
      }
    }

    // Divide by number of data points in cluster
    // This is the new centroid (average)
    for (int j = 0; j < num_features; j++)
    {
      if (counts[c] == 0)
        new_centroids[c][j] = mins[j] +
                              rand() % (int)maxes[j]; //600000; // If no data points in group, then reinitialize
      else
        new_centroids[c][j] = new_centroids[c][j] / counts[c];
    }
  }
  return new_centroids;
}

std::string getColor(int val)
{
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

void writeCentroidsToFile(string centroids_file, double **final_centroid, int k, int num_features)
{
  std::ofstream outfile;
  outfile.open(centroids_file);
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

void writeLabelsToFile(string filename, double **x, int *labels, int num_samples, int num_features)
{
  std::string color;
  std::ofstream outfile;
  outfile.open(filename);
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