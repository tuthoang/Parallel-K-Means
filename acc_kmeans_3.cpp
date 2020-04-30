// #include <accelmath.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
// #include <algorithm> // std::copy
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <iostream>

// using namespace std;

// #define k 8
// #define num_samples 6500
// #define num_features 2
// #define iterations 10000

int *getLabels(double **x, double **centroids, int num_samples, int num_features, int k);
double **getCentroids(double **x, double **centroids, int *clusters, int num_samples, int num_features, int k, double *mins, double *maxes);
void writeCentroidsToFile(std::string centroids_file, double **final_centroid, int k, int num_features);
void writeLabelsToFile(std::string filename, double **x, int *labels, int num_samples, int num_features);

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
  #pragma acc init

  clock_t tStart = clock();
  srand(157); // seed 1
  // srand(time(NULL));
  setlocale(LC_ALL, "en_US.UTF-8");
  // Synthetic data
  int const k = 8;
  int const num_samples = 6500;
  int const num_features = 2;
  int const iterations = 10000;
  double **x = new double *[num_samples];

  for (int i = 0; i < num_samples; i++)
  {
    x[i] = new double[num_features];
    for(int j = 0; j < num_features; j++)
      x[i][j] = 0;
  }

  std::string filename = "synthetic_dataset.txt";
  std::string ground_truth_filename = "synthetic_dataset_gt.txt";

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
  // **this is my issue here ^^^
  // for (int i = 0; i < num_samples; i++)
  // {
  //   for (int j = 0; j < num_features; j++)
  //   {
  //     printf("%f ", x[i][j]);
  //   }
  //   printf("\n");
  // }

  std::cout<<linenum << "   "<<num_samples<<std::endl;
  //All the neccessarry array declarations
  double **centroids = new double *[k];
  double **min_centroids = new double *[k];
  double **initial_centroids = new double *[k];
  double *clusters = new double[num_samples];
  double **distances = new double *[num_samples];
  int *mins = new int[num_features];
  int *maxes = new int[num_features];
  double **old_centroids = new double *[k];
  int *counts = new int[num_features];
  int *wcss_array = new int[iterations];
  double **final_centroids = new double *[k];

  // Initialize clusters randomly, but only within the min-max range
  mins = new int[num_features];
  mins[0] = x[0][0];
  mins[1] = x[0][1];
  maxes = new int[num_features];
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
  printf("this is max %d and %d \n",maxes[0],maxes[1]);
  printf("this is min %d and %d \n",mins[0],mins[1]);
  int **random = new int*[num_samples];
  for (int i = 0; i < num_samples; i++)
  {
    random[i] = new int[num_features];
    for(int j = 0; j < num_features; j++){
        random[i][j] = (mins[j] + rand())%maxes[j];
    }
    
  }
   for (int i = 0; i < 100; i++)
  {
    
    for(int j = 0; j < num_features; j++){
        printf("%d  ",random[i][j]);
    }
    printf("\n");
    
  }
  
  for(int z = 0; z < iterations; z++){
    wcss_array[z] = 0;
  }
  for(int i = 0; i < k; i++){
    final_centroids[i] = new double[num_features];
    centroids[i] = new double[num_features];
    initial_centroids[i] = new double[num_features];
    old_centroids[i] = new double[num_features];
    min_centroids[i] = new double[num_features];
    for (int j = 0; j < num_features; j++)
    {
      final_centroids[i][j] = 0;
      centroids[i][j] = 0;
      initial_centroids[i][j] = 0;
      old_centroids[i][j] = 0;  
      min_centroids[i][j] = 0;
      }
  }
  
  printf("init random 2\n");

  for(int i = 0; i <num_samples; i++){
    distances[i] = new double[k];
    for(int j = 0; j < k ; j++){
      distances[i][k] =0;
    }
  }


  int counter = 0;

  printf("begin big loop\n");

  int thresh_met_counter = 0;
  double min_WCSS = INT_MAX;
  int minIndex = 0;
  // int min_WCSS = INT_MAX;

  // printf("%d\n", min_WCSS);
  //This big for loop is sort of independent we can just copy all of the clusters into an array and see which ones have the min wcss and then select the one with the min


    #pragma acc data copyin(x [0:num_samples] [0:num_features], \
                          random [0:num_samples][0:num_features], \
                          maxes [0:num_features],                  \
                          mins [0:num_features])      \
        // copy(centroids [0:k] [0:num_features], clusters[0:num_samples], distances[0:num_samples][0:k], counts[0:k])
    for (int loop1 = 0; loop1 < iterations; loop1++)
    {
      int counter1 = 0;
      clock_t tStart1 = clock();

      printf("Counter: %d\n", counter);
      // fflush(stdout);

      // printf("1\n");
      // fflush(stdout);
      for (int i = 0; i < num_samples; i++)
      {
        clusters[i] = 0;
      }
      // #pragma acc parallel loop vector_length(k)
      // #pragma acc parallel loop
      // #pragma acc parallel pcopy(counts [0:k])
      for (int i = 0; i < k; i++)
      {
        counts[i] = 0;
      }

      // Do some preprocessing
      counter = 0;

      // #pragma acc parallel loop present(random[:num_samples][:num_features]) collapse(2)
        for (int i = 0; i < k; i++)
        {
          for (int j = 0; j < num_features; j++)
          {
            //  printf("%d \n",counter);
              counter = counter%num_samples;
              centroids[i][j] = (random[counter][j] + mins[j])%maxes[j];
              // printf("%f :         %d                  %f \n", random[counter][j], counter, random[counter][j]);
              // fflush(stdout);

              counter = (counter + loop1*17)%num_samples;

              // counter++;  
          }
          // printf("\n");
          // fflush(stdout);
        }

        // #pragma acc loop collapse(2)
        for (int i = 0; i < k; i++)
        {
          for (int j = 0; j < num_features; j++)
          {
            old_centroids[i][j] = centroids[i][j];
          }
        }
        
        // printf("111\n");
        // fflush(stdout);

        // double **final_centroid = kmeans(x, initial_centroids, num_samples, num_features, k, mins, maxes);

        double thresh = 0.05;
        thresh_met_counter = 0;
       clock_t tStartb = clock();

        int count = 0;

        // #pragma acc loop private(thresh_met_counter, min_WCSS,minIndex)

        // #pragma acc data copy(centroids[0:k][0:num_features]) copy(distances[0:num_samples][0:k])
        // #pragma acc data copyin(x [0:num_samples] [0:num_features], \
        //                         centroids[0:k][0:2]) \
        //                         create(distances[0:num_samples][0:k]) 
                          // create(counts[0:k])
        while((thresh_met_counter < (k * num_features) && count < 10000) )
        // while(count<1000)
        {
          // printf("count: %d\n", count);
            thresh_met_counter = 0;

            //Loop through each sample
            //Loop through cluster for the sample and find closest centroid
            double closest_dist;
            // #pragma acc loop independent

            // #pragma acc loop gang
            double min_dist = (double)INT_MAX; 
            // #pragma acc loop
            for (int i = 0; i < num_samples; i++)
            {
              min_dist = (double)INT_MAX;
              // #pragma acc parallel loop
              for (int c = 0; c < k; c++)
              {
                //Calculate l2 distance from each cluster
                //This is a data independet loop so I should be able to do a parallelization

                distances[i][c] = sqrt(pow(x[i][0] - centroids[c][0], 2) + pow(x[i][1] - centroids[c][1], 2));
                // min_dist = fmin(min_dist, distances[i][c]);
                if (distances[i][c] <= min_dist)
                {
                  min_dist = distances[i][c];
                  // cout<< c << endl;
                  clusters[i] = c;
                }
              }
            }

            // for (int i = 0; i< num_samples; i++){
            //   closest_dist = INT_MAX;
            //   for(int c = 0; c < k; c++)
            //     if(distances[i][c] < closest_dist)
            //     {
            //       // cout<< c << endl;
            //       closest_dist = distances[i][c];
            //       clusters[i] = c;
            //     }
            // }
            // #pragma acc parallel loop present(centroids[0:k][0:num_features], x[0:num_samples][0:num_features], counts[0:k])
            // for (int c = 0; c < k; c++)
            // {
            //     counts[c] = 0;
            //     // #pragma acc loop independent
            //     for (int i = 0; i < num_samples; i++)
            //     {
            //       if (clusters[i] == c)
            //       {
            //           counts[c] += 1;
            //           // #pragma acc loop independent
            //           for (int j = 0; j < num_features; j++)
            //           {
            //             centroids[c][j] += x[i][j];
            //           }
            //       }
            //     }
            // }
            // counts holds the number of data points currently in the cluster
            // counts = new int[k];
            // #pragma acc loop
            #pragma acc parallel loop vector
            for (int c = 0; c < k; c++){
              // #pragma acc atomic write
              counts[c] = 0;
              // #pragma acc atomic write
              centroids[c][0] = 0.0;
              // #pragma acc atomic write
              centroids[c][1] = 0.0;
            }
            // #pragma acc update device(centroids[0:k][0:num_features])

            // Get number of points in each cluster.
            #pragma acc parallel loop copy(counts[0:k])
            for (int i = 0; i < num_samples; i++)
            {
              #pragma acc atomic update
              counts[(int)clusters[i]]++;
            }

            // If no data points in any cluster, then reinitialize the centroid
            for(int c = 0; c < k; c++){
              if (counts[c] == 0)
              {
                counts[c] = 1;
                for (int j = 0; j < num_features; j++)
                {
                  centroids[c][j] = (random[counter1][j] + mins[j]) % maxes[j]; 
                }
                counter1 = (abs(counter1 + 39 + c * 82)) % num_samples;
              }
            }
            // #pragma acc update device(centroids[0:k][0:num_features])


            #pragma acc parallel loop copy(centroids[0:k][0:num_features]) 
            for (int i = 0; i < num_samples; i++)
            {
              int c = clusters[i]; // label of data point i

              #pragma acc atomic update
              centroids[c][0] += (x[i][0]);

              #pragma acc atomic update
              centroids[c][1] += (x[i][1]);
            }
            // // #pragma acc update self(centroids[0:k][0:num_features])
            // // #pragma acc update device(centroids[0:k][0:num_features])

            // // #pragma acc parallel loop present(counts[0:k], random[0:num_samples][0:num_features], centroids[0:k][0:num_features]) 

            #pragma acc parallel loop copy(centroids[0:k][0:2]) copyin(counts[0:k])
            for (int c = 0; c < k; c++)
            {
                // Divide by number of data points in cluster
                // This is the new centroid (average)
                centroids[c][0] = centroids[c][0] / counts[c];

                centroids[c][1] = centroids[c][1] / counts[c];

            }
            // #pragma acc update self(centroids[0:k][0:num_features])
            // #pragma acc update device(centroids[0:k][0:num_features])

            // fflush(stdout);
            #pragma acc parallel loop reduction(+:thresh_met_counter)
            for (int i = 0; i < k; i++)
            {
              #pragma acc loop reduction(+: thresh_met_counter)
              for (int j = 0; j < num_features; j++)
              {
                // cout << abs(centroids[i][j] - old_centroids[i][j]) <<endl;
                double temp = fabs(centroids[i][j] - old_centroids[i][j]) / fabs(centroids[i][j]);
                // temp = temp*temp;
                // temp = pow(temp,1/2);
                if (temp < thresh)
                {
                  thresh_met_counter++;
                }
              }
            }
          
            #pragma acc parallel loop collapse(2)
            for (int i = 0; i < k; i++)
            {
              for (int j = 0; j < num_features; j++)
              {
                old_centroids[i][j] = centroids[i][j];
                }
            }
            
            count++;
          
        }

        // double t1 = (double)(clock() - tStart) / CLOCKS_PER_SEC;
        // printf("Time taken for clustering serially : %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        
        //WITH IN SUM OF SQUARES
        double wcss = 0;
        double wcss_cluster = 0;

        // for(int c = 0; c < k; c++){
        //     wcss_cluster = 0;

            // #pragma acc loop reduction(+: wcss_cluster)
            for(int i = 0; i < num_samples; i++){
              // if(clusters[i] == c){
                  // #pragma acc loop reduction(+: wcss_cluster)
                  for(int j = 0; j < num_features; j++){
                    // printf("clusters: %d                x: %d             centroids: %d\n", c, x[i][j], centroids[c][j]);
                    wcss += (x[i][j] - centroids[(int)clusters[i]][j]) * (x[i][j] - centroids[(int)clusters[i]][j]);
                  }
              }
            // }
            wcss =  sqrt(wcss);
        // }
       
        wcss =  wcss/num_samples;

        // printf("WCSS : %f\n", wcss);

        if(wcss < min_WCSS){

          min_WCSS = wcss;
          // printf("MIN WCSS : %f\n", min_WCSS);
          // fflush(stdout);
          // printf("This is the loop that resulted in a change %d %.9f %d \n", loop1, wcss,count);
          // minIndex = loop1;
          // fflush(stdout);

          for (int c = 0; c < k; c++)
          {
            for (int j = 0; j < num_features; j++)
            {
              min_centroids[c][j] = centroids[c][j];
            }
          }

        }

        // final_centroids[loop] = centroids;
        printf("COUNT: %d \t\t Time taken for clustering acc : %.2fs\n", count, (double)(clock() - tStartb) / CLOCKS_PER_SEC);
        // printf("This is the current wcss %.9f \n", min_WCSS);
         printf("This is the loop %d \n", loop1);
        fflush(stdout);


    }
      // std::cout << "end" << endl;
  // #pragma acc update self(centroids[0:k][0:num_features])
  printf("This is the min index and this is the minWCSS %d %.9f \n", minIndex,min_WCSS);  
  printf("Time taken for clustering acc : %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

  // printf("This is the loop %d \n", loop);
  // #pragma acc exit data copyout(min_centroids[0:k][0:num_features], clusters[0:num_samples])

  printf("min centroids\n");
  for(int i = 0; i < k ; i++){
    for(int j = 0; j < num_features; j++){
      // printf("000\n");

      printf("%f ", min_centroids[i][j]);
    }
    printf("\n");
  }

  std::string out = "acc_data_colors.txt";
  std::string centroids_file = "acc_centroids.txt";
  writeCentroidsToFile(centroids_file, min_centroids, k, num_features);
  int *labels = getLabels(x, min_centroids, num_samples, num_features, k);
  writeLabelsToFile(out, x, labels, num_samples, num_features);
  std::stringstream fmt;
  fmt << "gnuplot -p -e \"k=" << k << "; data_labels= '" << out
      << "'; centroids_file= '" << centroids_file
      << "'; Title='K-Means'\" plotCluster.gp";

  // PLOT WITH GNUPLOT
  // system(fmt.str().c_str());

  double **ground_truth = new double *[k];
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
      std::cout << ground_truth[i][j] << "  ";
    }
    std::cout << std::endl;
  }
  // printf("asdasdsaasd\n");
  std::string ground_truth_labels_file = "ground_truth_centroids.txt";
  labels = getLabels(x, ground_truth, num_samples, num_features, k);
  writeLabelsToFile(ground_truth_labels_file, x, labels, num_samples, num_features);
  // printf("2132132121asdasdsaasd\n");

  std::stringstream fmt1;
  fmt1 << "gnuplot -p -e \"k=" << k << "; data_labels= '" << ground_truth_labels_file 
        << "'; centroids_file= '" << ground_truth_filename 
        << "'; Title='Ground Truth' \" plotCluster.gp";

  // PLOT GROUND TRUTH
  // system(fmt1.str().c_str());
}

/**
 * Assigns clusters on the given data by
 * calculating the closest distance to current centroids
 */
int *getLabels(double **x, double **centroids, int num_samples, int num_features, int k)
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

void writeCentroidsToFile(std::string centroids_file, double **final_centroid, int l, int features)
{
  std::ofstream outfile;
  outfile.open(centroids_file);
  printf("\n");
  printf("Final Centroids\n");
  std::string color;
  for (int i = 0; i < l; i++)
  {
    color = getColor(i);
    printf("Centroid %d: (", i);
    for (int j = 0; j < features; j++)
    {
      printf("%f ", final_centroid[i][j]);
      outfile << final_centroid[i][j] << " ";
    }
    printf(")\n");
    outfile << i << "\n";
  }
  outfile.close();
}

void writeLabelsToFile(std::string filename, double **x, int *labels, int samples, int features)
{
  std::string color;
  std::ofstream outfile;
  outfile.open(filename);
  for (int i = 0; i < samples; i++)
  {
    color = getColor(labels[i]);
    for (int j = 0; j < features; j++)
    {
      outfile << x[i][j] << " ";
    }
    outfile << labels[i] << "\n";
  }
  outfile.close();
}
