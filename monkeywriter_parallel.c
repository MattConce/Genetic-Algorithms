#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <ctype.h>
#include <unistd.h>

#define alphalen 58

// Dna is equivalent to an individual
struct dna {
   int len;
   double score;
   char *gene;
};

typedef struct dna *DNA;

void calcFitness(DNA *population, int size, char *target, double record[]) {
   int best, i = 0; 
   double bestScore = 0.0;
   #pragma omp parallel for shared(size, best) private(i) \
   reduction(max : bestScore)
   for (i = 0; i < size; i++) {
      for (int j = 0; j < population[i]->len; j++) {
         if (population[i]->gene[j] == target[j])
            population[i]->score += 1.0;
      }
      population[i]->score = (population[i]->score / (double) population[i]->len);
      population[i]->score = pow(population[i]->score, 2.0);

      if (bestScore < population[i]->score) {
         bestScore = population[i]->score;
         best = i; 
      }
   }
   record[0] = (double) best, record[1] = bestScore;
}

DNA* initPopulation(int size, int len, char *target, unsigned int *seed) {
   DNA *population = malloc(size * sizeof(DNA));

   int i = 0;

   #pragma omp parallel for private(i) 
   for (i = 0; i < size; i++) {
      population[i] = malloc(sizeof(struct dna));
      population[i]->len = len;
      population[i]->score = 0.0;
      population[i]->gene = malloc(len * sizeof(char));
      for (int j = 0; j < len; j++) {
         population[i]->gene[j] = (char) 65 + rand_r(seed) % alphalen;;
      }
   }
   return population;
}

void destroyPopulation(DNA *population, int size) {
   for (int i = 0; i < size; i++) {
      free(population[i]->gene);
      free(population[i]);
   }
   free(population);
}

void showDNA(DNA *population, int size, int generations, double record[]) {
   printf(" \nGeneration %d\n", generations);
   printf(" \nPopulation size %d\n", size);
   int best = (int) record[0];
   for (int i = 0; i < size; i++) {
      for (int j = 0; j < population[i]->len; j++)
         printf("%c ", population[i]->gene[j]);
      printf(" score: %f \n", population[i]->score);
      printf("\n");
   }
   printf("best: ");
   for (int i = 0; i < population[best]->len; i++)
      printf("%c ", population[best]->gene[i]);
   printf("\n");
}

int binarySearch(double arr[], int l, int r, double x) {
   if (r >= l) { 
      int mid = l + (r - l) / 2; 

      if (arr[mid] == x) 
         return mid; 

      if (arr[mid] > x) 
         return binarySearch(arr, l, mid - 1, x); 

      return binarySearch(arr, mid + 1, r, x); 
   } 

   return l; 
} 

DNA findParent(DNA *population, int size, unsigned int *seed) {
   double *arr = calloc(size, sizeof(double));
      
   for (int i = 0; i < size; i++)
      arr[i] = population[i]->score;

   for (int i = 1; i < size; i++)
      arr[i] += arr[i-1];

   double x = rand_r(seed) / (RAND_MAX + 1.0);
   double k = x * arr[size-1];
   int index = binarySearch(arr, 0, size-1, k);
   free(arr);
   return population[index];
}

DNA crossover(DNA a, DNA b, int len) {
   DNA offspring = malloc(sizeof(struct dna));
   offspring->len = len;
   offspring->score = 0.0;
   offspring->gene = malloc(len * sizeof(char));
   int mid = floor(len/2);
   for (int i = 0; i < len; i++) {
      if (i > mid)
         offspring->gene[i] = a->gene[i];
      else
         offspring->gene[i] = b->gene[i];
   }
   return offspring;
}

void mutation(DNA c, double r, unsigned int *seed) {
   for (int i = 0; i < c->len; i++) {
      double rr = rand_r(seed) / (RAND_MAX + 1.0); 
      if (rr < r) {
         c->gene[i] = (char) 65 + rand_r(seed) % alphalen;;
      }
   }
}


DNA* generate(DNA *population, int size, int len, char *target, 
      double r, unsigned int *seed) {

   DNA *newPopulation = initPopulation(size, len, target, seed);
   for (int i = 0; i < size; i++) {
      DNA parentA = findParent(population, size, seed);
      DNA parentB = findParent(population, size, seed);
      DNA offspring = crossover(parentA, parentB, len);
      mutation(offspring, r, seed);
      newPopulation[i] = offspring;
   }
   destroyPopulation(population, size);
   return newPopulation;
}

int evaluate(DNA eval, int size) {
   for (int i = 0; i < size; i++)
      if (eval->score >= 1.0)
         return 1;
   return 0;
}

char* readInput(int *length) {
   int len;
   scanf("%d",&len);
   getchar();
   char *input = malloc((len+1) * sizeof(char));
   int i;
   for (i = 0; i < len; i++) {
      input[i] = getchar();
   }
   input[i] = '\0';
   *length = len;
   return input;
}


int main(int argc, char *argv[]) {

  if (argc != 4) {
    printf("Type ./monkeywriter_parallel2 population_size=int debug={0,1} nthreads=int\n");
    printf("%d",argc);
    return 0;
  }

  int size  = atoi(argv[1]);
  int debug = atoi(argv[2]);

  int done = 0, generations = 0;

  // Number of Threads
  omp_set_num_threads(atoi(argv[3]));

  /*
  int alphaLen = 56;
  char *alphabet = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
                    'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
                    'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g',
                    'h', 'i', 'j', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's',
                    't', 'u', 'v', 'w', 'x', 'y', 'z', ' ', ',', '.', '!',
                    '?'};
  */
   
   int len;
   char *target = readInput(&len);

   unsigned int seed = 0;

   DNA* population = malloc(size*sizeof(DNA));

   // Set the seed for the random numbers
   seed = time(NULL);

   population = initPopulation(size, len, target, &seed);
   // Standard genetic algorithm
   while(!done && generations < 20000) {
      double max_fitness = -INFINITY;
      double record[] = {0.0, 0.0};
      // Calculate the fitness
      calcFitness(population, size, target, record);
      if (debug)
        showDNA(population, size, generations, record);
      int best = (int) record[0];
      // Evaluate the population
      done = evaluate(population[best], size);
      if (record[1] > max_fitness) {
        max_fitness  = record[1];
      }
      // Generate the next population based on the parents
      population = generate(population, size, len, target, 0.01, &seed);
      generations++;
   }
   free(target);
   printf("Generations: %d\n", generations);
   return 0;
}
