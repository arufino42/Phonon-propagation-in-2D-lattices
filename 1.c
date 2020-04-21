#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dt 1e-9

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

void run_sim_square_z(state *system)
{
  double **temp;
  temp_pos=create_matrix(system->m,system->n);
  tem_vel=create_matrix(system->m,system->n);
  int i,j; 
  double accel;
  double omega=1;
  printf("---");
  for(i=0;system->m;++i)
  {
    for(j=0;system->n;++j)
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
  system->velocity=temp_vel;
  system->position=temp_pos;
  for(i=0;i<system->m;++i)
  {
    free(temp_vel[i]);
    free(temp_pos[i]);
  }
}

int main()
{
  state grid;
  grid.position = create_matrix(3,3);
  grid.velocity = create_matrix(3,3);
  grid.l=1;
  grid.f0=1;
  grid.k=1;
  grid.mass=1;
  grid.m=3;
  grid.n=3;
  int i,j;
  run_sim_square_z(&grid);
  for(i=0;i<3;++i)
  {
    for(j=0;j<3;++j)
    {
      printf("%g  ",grid.position[i][j]);
    }
    printf("\n");  
  }
  return 0;
}