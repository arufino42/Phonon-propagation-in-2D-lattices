
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dt 1e-4

typedef struct 
{
  double** position;
  double** velocity;
  double mass;
  double l;
  double f0;
  double k;
  int m;
  int n;
} state;

void import(char* filename, double*** matriz, int m, int n)
{
  FILE *f1;
  f1=fopen(filename,"rt");
  int i,j;
  for(i=0;i<m;++i)
  {
    for(j=0;j<n;++j)
    {
      fscanf(f1,"%lf\n",&((*matriz)[i][j]));
      printf("%lf  ",(*matriz)[i][j]);
    }
    printf("\n");
  }
  fclose(f1);
}

double** create_matrix(int m,int n)
{
  double **matrix;
  matrix = (double**) malloc (m*sizeof(double*));
  int i,j;
  for(i=0;i<m;++i)
  {
    matrix[i]=(double*) calloc(n,sizeof(double));
  }
  return matrix;
}

void run_sim_square_z(state *system,int no_iterations)
{
  double **temp_pos,**temp_vel;
  temp_pos=create_matrix(system->m,system->n);
  temp_vel=create_matrix(system->m,system->n);
  int i,j,k; 
  double accel;
  double omega=(system->f0)/((system->m)*(system->l));
  for(k=0;k<no_iterations;k++)
  {
    for(i=0;i<(system->m);++i)
    {
      for(j=0;j<(system->n);++j)
      {
        accel=0;
        accel-=4*system->position[i][j];
        if(i>0)
          accel+=system->position[i-1][j];
        if(i<(system->m)-1)
          accel+=system->position[i+1][j];
        if(j>0)
          accel+=system->position[i][j-1];
        if(j<(system->n)-1)
          accel+=system->position[i][j+1];
        accel *= omega;
        temp_vel[i][j] = system->velocity[i][j] + accel*dt;
        temp_pos[i][j] = system->position[i][j] + temp_vel[i][j]*dt;
      }
    }
  }
  
  system->velocity=temp_vel;
  system->position=temp_pos;
  for(i=0;i<system->m;++i)
  {
    free(temp_vel[i]);
    free(temp_pos[i]);
  }
}

void run_sim_square_xy(state *system,int no_iterations)
{
  int i,j,k; 
  double accel;
  double omega=(system->f0)/((system->m)*(system->l));
  for(k=0;k<no_iterations;k++)
  {
    for(i=0;i<(system->m);++i)
    {
      for(j=0;j<(system->n);++j)
      {
        accel=0;
        accel-=4*system->position[i][j];
        if(i>0)
          accel+=system->position[i-1][j];
        if(i<(system->m)-1)
          accel+=system->position[i+1][j];
        if(j>0)
          accel+=system->position[i][j-1];
        if(j<(system->n)-1)
          accel+=system->position[i][j+1];
        accel *= omega;
        system->velocity[i][j] = system->velocity[i][j] + accel*dt;
        system->position[i][j] = system->position[i][j] + system->velocity[i][j]*dt;
      }
    }
  }
}

int main()
{
  FILE *f1,*f2,*f3,*f4;
  f1 = fopen("Ficheiros/file00.txt","wt");
  f2 = fopen("Ficheiros/file44.txt","wt");
  f3 = fopen("Ficheiros/file23.txt","wt");
  f4 = fopen("Ficheiros/file32.txt","wt");
  state grid;
  grid.position = create_matrix(5,5);
  grid.velocity = create_matrix(5,5);
  grid.velocity[0][0]=1;
  grid.l=1;
  grid.f0=1;
  grid.k=1;
  grid.mass=1;
  grid.m=5;
  grid.n=5;
  int i,j,k;
  for(k=0;k<2000;k++)
  {
  run_sim_square_z(&grid,1);
  /*for(i=0;i<10;++i)
  {
    for(j=0;j<10;++j)
    {
      printf("%g  ",grid.position[i][j]);
    }
    printf("\n");  
  }*/
  fprintf(f1,"%g  \n",grid.position[0][0]);
  fprintf(f2,"%g  \n",grid.position[4][4]);
  fprintf(f3,"%g  \n",grid.position[2][3]);
  fprintf(f4,"%g  \n",grid.position[0][4]);
  }
  return 0;
  fclose(f1);
  fclose(f2);
}
