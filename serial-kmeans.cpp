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

double **kmeans(double **x, double **initial_centroids, int num_samples, int num_features, int k);
int *getLabels(double **x, double **centroids, int num_samples, int num_features, int k);
void getCentroids(double **x, double **centroids, int *clusters, int num_samples, int num_features, int k);
void writeCentroidsToFile(double **final_centroid, int k, int num_features);
void writeLabelsToFile(double **x, int *labels, int num_samples, int num_features);

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
  srand(50); // seed 1
  setlocale(LC_ALL, "en_US.UTF-8");
  // Synthetic data
  int k = 8;              // Number of clusters
  int num_features = 2;   // x1, x2
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
        x[linenum][i] = std::stod(result[i]);
      }
      linenum++;
    }
  }
  myfile.close();

  // Do some preprocessing
  // Initialize clusters randomly, but only within the min-max range
  double *mins = new double[num_features];
  mins[0] = x[0][0];
  mins[1] = x[0][1];
  double *maxes = new double[num_features];
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

  double **initial_centroids = new double *[k];
  for (int i = 0; i < k; i++)
  {
    initial_centroids[i] = new double[num_features];
    for (int j = 0; j < num_features; j++)
    {
      initial_centroids[i][j] = mins[j] +
                                rand() % (int)maxes[j];
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
  double **final_centroid = kmeans(x, initial_centroids, num_samples, num_features, k);
  double t1 = (double)(clock() - tStart)/CLOCKS_PER_SEC;
  printf("Time taken for clustering serially : %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);  


  writeCentroidsToFile(final_centroid, k, num_features);
  int* labels = getLabels(x, final_centroid, num_samples, num_features, k);
  writeLabelsToFile(x, labels, num_samples, num_features);
  system("gnuplot -p plotCluster");
}

double **kmeans(double **x, double **centroids,
                int num_samples, int num_features, int k)
{

  double **old_centroids = centroids;
  int *clusters;

  int count = 0;
  while (count < 10000)
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
    clusters = getLabels(x, centroids, num_samples, num_features, k);
    getCentroids(x, centroids, clusters, num_samples, num_features, k);
    count++;
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
void getCentroids(double **x, double **centroids, int *clusters,
                    int num_samples, int num_features, int k)
{
  // double **new_centroids = new double*[k];

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
        counts[c]++;
        for (int j = 0; j < num_features; j++)
        {
          centroids[c][j] += x[i][j];
        }
      }
    }

    // Divide by number of data points in cluster
    // This is the new centroid (average)
    for (int j = 0; j < num_features; j++)
    {
      if (counts[c] == 0)
        centroids[c][j] = rand() % 500000; // If no data points in group, then reinitialize
      else
        centroids[c][j] = centroids[c][j] / counts[c];
    }
  }
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

void writeLabelsToFile(double **x, int *labels, int num_samples, int num_features)
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