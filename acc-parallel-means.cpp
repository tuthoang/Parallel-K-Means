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

#include <curand.h>

using namespace std;

#define threshold 0.1
// #pragma acc routine(sqrt) seq

double** kmeans(double **x, double **initial_centroids, int num_samples, int num_features, int k);
int* getLabels(double **x, double **centroids, int num_samples, int num_features, int k);
double** getCentroids(double **x, double **centroids, int *clusters, int num_samples, int num_features, int k);
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
      distances[i][j] = 0.0;
  }
  double *clusters = new double[num_samples];
  for(int i = 0; i < num_samples; i++){
    clusters[i] = 0;
  }
  double **centroids = initial_centroids;
  int count = 0;
  int rand_counter = 0;


  #pragma acc data create(clusters[0:num_samples]) copy(centroids [0:k] [0:num_features]) copyin(x[0:num_samples][0:num_features])
  for (count =0; count < 10000; count++)
  {
    #pragma acc loop independent
    for (int i = 0; i < num_samples; i++)
    {
      for (int j = 0; j < k; j++)
        distances[i][j] = 0.0;
    }
    double closest_dist;
    // Loop through each sample
    // Loop through cluster for the sample and find closest centroid
    for (int i = 0; i < num_samples; i++)
    {
      closest_dist = INT_MAX;
      // #pragma acc parallel loop reduction(+:l2_dist)
      #pragma acc loop independent
      // #pragma acc data create(distances [0:num_samples] [0:k]) 
      for (int c = 0; c < k; c++)
      {
        // Calculate l2 distance from each cluster
        //This is a data independet loop so I should be able to do a parallelization
        // #pragma acc parallel loop reduction(+:l2_dist)
        for (int j = 0; j < num_features; j++)
        {
          distances[i][c] += (x[i][j] - centroids[c][j]) * (x[i][j] - centroids[c][j]);
        }
        // distances[i][c] = sqrt(distances[i][c]);
      
      //cout<<"This is one of the clusters calculation " << distances[i][0]<<endl;
      closest_dist = distances[i][0];
      // #pragma acc data copy(distances [0:num_samples] [0:k])
      // #pragma acc parallel loop
      
      if (distances[i][c] < closest_dist)
      {
        closest_dist = distances[i][c];
        clusters[i] = c;
      }
      
      //cout<<"went through a loop fine"<<endl;
    }
  }
  



// ====

    //Initializing everything to 0
    // double **new_centroids = new double *[k];

    // for (int i = 0; i < k; i++)
    // {
    //   new_centroids[i] = new double[num_features];
    //   for (int j = 0; j < num_features; j++)
    //   {
    //     new_centroids[i][j] = 0;
    //     //cout<<centroids[i][j]<<endl;
    //   }
    // }


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

    // centroids = new_centroids;
    // for(int p = 0; p<k;p++){
    //   for(int q =0;q<num_features;q++)
    //     cout<<centroids[p][q] << " ";
    //   cout<<endl;
    // }
    // count++;
    cout << count << endl;
  }
  
  writeCentroidsToFile(centroids, k, num_features);

  printf("Time taken for clustering parallel : %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  int *labels = getLabels(x, centroids, num_samples, num_features, k);
  writeLabelsToFile(x, labels, num_samples, num_features);
  system("gnuplot -p plotCluster");
}

// double **kmeans(double **x, double **centroids,
//                 int num_samples, int num_features, int k)
// {

//   double** new_centroids = new double*[k];
//   int *clusters;

//   int count = 0;
//   while (count < 10000)
//   {
//     // std::copy(&centroids[0][0], // Used for convergence checking
//     //           &centroids[0][0] + num_samples * num_features,
//     //           old_centroids[0][0]);
//     // printf("\n");
//     // for (int i = 0; i < k; i++)
//     // {
//     //   printf("Centroid %d: (", i);
//     //   for (int j = 0; j < num_features; j++)
//     //   {
//     //     printf("%f ", centroids[i][j]);
//     //   }
//     //   printf(")\n");
//     // }
//     clusters = getLabels(x, centroids, num_samples, num_features, k);
//     new_centroids = getCentroids(x, centroids, clusters, num_samples, num_features, k);
//     for(int i = 0; i < k; i++){
//         for(int j = 0; j < num_features; j++){
//             centroids[i][j] = new_centroids[i][j];
//         }
//     }
//     count++;
//     cout<<count<<endl;
//   }

//   return centroids;
// }

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
  double** distances = new double*[num_samples];

  // #pragma acc parallel
  for (int i = 0; i < num_samples; i++)
  {
    distances[i] = new double[k];
    for (int j = 0; j < k; j++)
      distances[i][j] = 0.0;
  }
  // #pragma acc data copyin(centroids[0:k][0:num_features]), create(distances[0:num_samples][0:k])
  for (int i = 0; i < num_samples; i++)
  {
    closest_dist = INT_MAX;
    for (int c = 0; c < k; c++)
    {
      // Calculate l2 distance from each cluster
      // l2_dist = 0.0;
      //This is a data independet loop so I should be able to do a parallelization
      // #pragma acc parallel loop reduction(+:l2_dist)
      for (int j = 0; j < num_features; j++)
      {
        l2_dist += (x[i][j] - centroids[c][j]) * (x[i][j] - centroids[c][j]);
      }
      distances[i][c] = sqrt(l2_dist);
      // Assign closest centroid to data point
      // distances[i][c] = l2_dist;
    }
    //cout<<"This is one of the clusters calculation " << distances[i][0]<<endl;
    closest_dist = distances[i][0];
    // #pragma acc kernels
    for (int c = 1; c < k; c++)
    {
      if (distances[i][c] < closest_dist)
      {
        closest_dist = distances[i][c];
        clusters[i] = c;
      }
    }
    //cout<<"went through a loop fine"<<endl;
  }
    //cout<< "got here "<<endl;
    return clusters;
}

  /**
 * Updates the centroids by calculating
 * the mean of the data points belonging to that cluster
 * 
 * Idea now is to somehow stop updating some clusters if their mean hasn't shifted too much, therefore we need a way to store and compare the old centroids 
 * with new ones and then somehow update the dataset. 
 */
  double **getCentroids(double **x, double **centroids, int *clusters,
                        int num_samples, int num_features, int k)
  {
    //Initializing everything to 0
    double **new_centroids = new double *[k];

    for (int i = 0; i < k; i++)
    {
      new_centroids[i] = new double[num_features];
      for(int j = 0; j < num_features; j++){
          new_centroids[i][j] = 0;
          //cout<<centroids[i][j]<<endl;
      }
  }
  // int random[num_samples];
  // for (int i = 0; i < num_samples; i++){
  //   random[i] = rand();
  // }

  // counts holds the number of data points currently in the cluster
  int *counts = new int[k];
  // #pragma acc data copy(random[0:num_samples])
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
        new_centroids[c][j] = rand() % 500000; // If no data points in group, then reinitialize
      else
        new_centroids[c][j] = new_centroids[c][j] / counts[c];
    }
  }

  return new_centroids;
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