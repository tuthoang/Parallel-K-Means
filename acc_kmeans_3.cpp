#include <accelmath.h>
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

#define k 20
#define num_samples 3000
#define num_features 2
#define iterations 10000

int *getLabels(double **x, double **centroids);
void writeCentroidsToFile(std::string centroids_file, double **final_centroid, int l, int features);
void writeLabelsToFile(std::string filename, double **x, int *labels, int samples, int features);

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
  
  double **x = new double *[num_samples];

  for (int i = 0; i < num_samples; i++)
  {
    x[i] = new double[num_features];
    for(int j = 0; j < num_features; j++)
      x[i][j] = 0.0;
  }

  std::string filename = "a1.txt";
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
        // std::cout<<x[linenum][i] << std::endl;
      }
      linenum++;
    }
  }
  myfile.close();

  // **this is my issue here ^^^
  for (int i = 0; i < num_samples; i++)
  {
    for (int j = 0; j < num_features; j++)
    {
      printf("%d ", x[i][j]);
    }
    printf("\n");
  }

  std::cout<<num_samples<<std::endl;
  //All the neccessarry array declarations 
  double **centroids = new double*[k];
  double **min_centroids = new double*[k];
  double **initial_centroids = new double*[k];
  int **clusters = new int*[iterations];
  double **distances = new double *[num_samples];
  double *mins = new double[num_features];
  double *maxes = new double[num_features];
  double **old_centroids = new double*[k];
  int **counts = new int*[iterations];
  double* wcss_array = new double[iterations];
  double*** final_centroids = new double**[iterations];
  
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
  printf("this is max %d and %d \n",maxes[0],maxes[1]);
  printf("this is min %d and %d \n",mins[0],mins[1]);
  int ***random = new int**[iterations];
  for (int i = 0; i < iterations; i++)
  {
    random[i] = new int*[num_samples];
    for(int j = 0; j <num_samples; j++){
        random[i][j] = new int[num_features];
        for(int l = 0; l < num_features; l++){
            random[i][j][l] = mins[l] +
                rand()%(int)maxes[l];
        }
    }
  }
  printf("init random\n");

  // for(int i = 0; i < num_samples; i++){
  //   for(int j = 0; j < num_features; j++){
  //     printf("%d ", random[0][i][j]);
  //   }
  //   printf("\n");
  // }
  for(int z = 0; z < iterations; z++){
    final_centroids[z] = new double*[k];
    counts[z] = new int[k];
    wcss_array[z] = 0;
    for(int i = 0; i < k; i++){
        final_centroids[z][i] = new double[num_features];
        counts[z][i] = 0;
        centroids[i] = new double[num_features];
        initial_centroids[i] = new double[num_features];
        old_centroids[i] = new double[num_features];
        for(int j = 0; j < num_features; j++){
            final_centroids[z][i][j] = 0;
            centroids[i][j] = 0;
            initial_centroids[i][j] = 0;
            old_centroids[i][j] = 0;  
        }
    }
  }
  printf("init random 2\n");

  for(int i = 0; i <num_samples; i++){
    distances[i] = new double[k];
    for(int j = 0; j < k ; j++){
      distances[i][k] =0;
    }
  }

  for(int i = 0; i < iterations; i++){
    clusters[i] = new int[num_samples];
    for(int j = 0; j < num_samples; j++){
      clusters[i][j] = 0;
    }
  }

  int counter = 0;
  printf("begin big loop\n");

  //This big for loop is sort of independent we can just copy all of the clusters into an array and see which ones have the min wcss and then select the one with the min
    #pragma data copyin(x [0:num_samples] [0:num_features]),                                                       \
    copyin(initial_centroids [0:k] [0:num_features]),                                                                 \
    copyin(final_centroids [0:iterations] [0:k] [0:num_features]),                                             \
    copyin(wcss_array [0:iterations]),                                                                               \
    copyin(random [0:iterations] [0:num_samples] [0:num_features]),                                                  \
    copyin(counts [0:iterations] [0:k], centroids [0:k] [0:num_features]),                                           \
    copyin((oldcentroids [0:k] [0:num_features]),                                                                    \
    copyin(clusters [0:iterations] [0:num_samples]),                                                         \
    copyin(distances [0:num_samples] [0:k]))
  {
    // #pragma acc loop vector reduction(min:)
    for (int loop = 0; loop < iterations; loop++)
    {
      printf("Counter: %d\n", counter);

      // Do some preprocessing
        for (int i = 0; i < k; i++)
        {
          #pragma acc loop vector
          for (int j = 0; j < num_features; j++)
          {
              centroids[i][j] = mins[j] + random[loop][counter][j] % (int)maxes[j];
              counter = (random[loop][counter][j] * random[(loop + loop)%iterations][counter][j])%num_samples;
              // counter++;  
          }
        }
        // centroids = initial_centroids;
        printf("Initial Centroids\n");
        for (int i = 0; i < k; i++)
        {
          printf("Centroid %d: (", i);
          for (int j = 0; j < num_features; j++)
          {
            printf("%f ", centroids[i][j]);
          }
          printf(")\n");
        }

        #pragma acc loop vector
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < num_features; j++)
            {
            old_centroids[i][j] = centroids[i][j];
            }
        }
        // double **final_centroid = kmeans(x, initial_centroids, num_samples, num_features, k, mins, maxes);

        double thresh = 0.05;
        int thresh_met_counter = 0;
        // int min_WCSS = INT_MAX;

        int count = 0;
        
        #pragma acc loop vector
       while (thresh_met_counter < (k * num_features) && count < 10000)
        // while(count<1000)
        {
            thresh_met_counter = 0;

            //Loop through each sample
            //Loop through cluster for the sample and find closest centroid
            double closest_dist;
            #pragma acc parallel loop vector_length(num_samples)
            for (int i = 0; i < num_samples; i++)
            {
                closest_dist = INT_MAX;
                #pragma acc loop independent
                for (int c = 0; c < k; c++)
                {
                  //Calculate l2 distance from each cluster
                  //This is a data independet loop so I should be able to do a parallelization
                  double sum_dist = 0;
                  #pragma acc loop vector reduction(+:sum_dist)
                  for (int j = 0; j < num_features; j++)
                  {
                      sum_dist += ((x[i][j] - centroids[c][j]) * (x[i][j] - centroids[c][j]));
                  }
                  distances[i][c] = sqrt(sum_dist);
                  // cout << distances[i][c] <<endl;
                  // closest_dist = distances[i][0];
                  if (distances[i][c] < closest_dist)
                  {
                      // cout<< c << endl;
                      closest_dist = distances[i][c];
                      clusters[loop][i] = c;
                  }
                }
            }


            // counts holds the number of data points currently in the cluster
            #pragma acc parallel loop vector_length(k)
            for (int c = 0; c < k; c++)
            {
                counts[loop][c] = 0;
                #pragma acc loop independent
                for (int i = 0; i < num_samples; i++)
                {
                  if (clusters[loop][i] == c)
                  {
                      counts[loop][c] += 1;
                      #pragma acc loop independent
                      for (int j = 0; j < num_features; j++)
                      {
                        centroids[c][j] += x[i][j];
                      }
                  }
                }
            }


            #pragma acc parallel loop vector_length(k)
            for (int c = 0; c < k; c++)
            {
                // Divide by number of data points in cluster
                // This is the new centroid (average)
                int counter1 =0 ;


                #pragma acc loop independent
                for (int j = 0; j < num_features; j++)
                {
                  if (counts[loop][c] == 0)
                  {
                      // cout<< c << "  Has no data points.  " << mins[j] << " \t" << maxes[j]<<endl;
                      centroids[c][j] = mins[j] + random[loop][counter1][j] % (int)maxes[j]; // If no data points in group, then reinitialize
                  }
                  else{
                      centroids[c][j] = centroids[c][j] / counts[loop][c];
                  }
                  // counter = random[loop][counter][j]*42%num_samples;
                  counter1++;
                }
                
            }

            #pragma acc parallel loop reduction(+:thresh_met_counter)
            for (int i = 0; i < k; i++)
            {
              #pragma acc loop reduction(+: thresh_met_counter)
              for (int j = 0; j < num_features; j++)
              {
                // cout << abs(centroids[i][j] - old_centroids[i][j]) <<endl;
                double temp = abs(centroids[i][j] - old_centroids[i][j]) / abs(centroids[i][j]);
                // temp = temp*temp;
                // temp = pow(temp,1/2);
                if (temp < thresh)
                {
                  thresh_met_counter++;
                }
              }
            }

            #pragma acc parallel loop vector_length(k)
            for (int i = 0; i < k; i++)
            {
              #pragma acc loop independent
              for (int j = 0; j < num_features; j++)
              {
                old_centroids[i][j] = centroids[i][j];
                }
            }


            count++;
            //writeLabelsToFile(x, clusters, num_samples, num_features);
            // cout << "This is the iter " << count << " this is the number of points that changed centroids " << thresh_met_counter << endl;
        }

        // double t1 = (double)(clock() - tStart) / CLOCKS_PER_SEC;
        // printf("Time taken for clustering serially : %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        // std::cout<<"This is the end of "<<count<<" iter"<<endl;  

        
        //WITH IN SUM OF SQUARES
        double wcss = 0;
        double wcss_cluster = 0;

        #pragma acc parallel loop reduction(+: wcss)
        for(int c = 0; c < k; c++){
            wcss_cluster = 0;

            #pragma acc loop reduction(+: wcss_cluster)
            for(int i = 0; i < num_samples; i++){
              if(clusters[loop][i] == c){
                  for(int j = 0; j < num_features; j++){
                    wcss_cluster += (x[i][j] - centroids[c][j]) * (x[i][j] - centroids[c][j]);
                  }
              }
            }
            wcss += sqrt(wcss_cluster);
        }
        // cout <<wcss<<endl;
        wcss = wcss/num_samples;
        // cout << wcss << endl;

        // if(wcss < min_WCSS){
        //     min_WCSS = wcss;
        //     for(int c = 0; c < k; c++){
        //         for(int j = 0; j < num_features; j++){
        //         min_centroids[c][j] = centroids[c][j];
        //         }
        //     }
        // }
        wcss_array[loop] = wcss;
        final_centroids[loop] = centroids;
         printf("This is the loop %d \n", loop);
        }
      // std::cout << "end" << endl;
     
  } 
  #pragma acc exit data copyout(final_centroids[0:iterations][0:k][0:num_features], wcss_array[0:iterations], clusters[0:iterations][0:num_samples])

  double min = wcss_array[0];
  int index = 0;
  for(int i = 0; i < 10000; i++){
    if(min > wcss_array[i]){
        index= i;
        min = wcss_array[i];
    }
  }
  min_centroids = final_centroids[index];
  for(int i = 0; i < k ; i++){
    for(int j = 0; j < num_features; i++){
      std::cout << min_centroids[i][j] << " ";
    }
    std::cout << std::endl;
  }
  int* final_cluster = new int[num_samples];
  final_cluster = clusters[index];
  // printf("InitialCentroids\n");
  // for (int i = 0; i < k; i++)
  // {
  //   printf("Centroid %d: (", i);
  //   for (int j = 0; j < num_features; j++)
  //   {
  //     printf("%f ", initial_centroids[i][j]);
  //   }
  //   printf(")\n");
  // }

  std::string out = "acc_data_colors.txt";
  std::string centroids_file = "acc_centroids.txt";
  writeCentroidsToFile(centroids_file, min_centroids, k, num_features);
  int *labels = getLabels(x, min_centroids);//, num_samples, num_features, k);
  writeLabelsToFile(out, x, labels, num_samples, num_features);
  std::stringstream fmt;
  fmt << "gnuplot -p -e \"k=" << k << "; data_labels= '" << out
      << "'; centroids_file= '" << centroids_file
      << "'; Title='K-Means'\" plotCluster.gp";

  // PLOT WITH GNUPLOT
  // system(fmt.str().c_str());


  double **ground_truth = new double*[k];
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
  labels = getLabels(x, ground_truth);//, num_samples, num_features, k);
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
int *getLabels(double **x, double **centroids)
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
