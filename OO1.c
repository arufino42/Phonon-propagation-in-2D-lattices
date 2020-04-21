
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dt 1e-4
#define M 50
#define N 50

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

double **temp_pos,**temp_vel;

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
    for(i=0;i<(system->m);++i)
    {
      for(j=0;j<(system->n);++j)
      {
        system->velocity[i][j]=temp_vel[i][j];
        system->position[i][j]=temp_pos[i][j];
      }
    }
  }
}

void run_sim_square_xy(state *system,int no_iterations)
{
  int i,j,k; 
  double accel;
  double omega1=(system->f0)/((system->m)*(system->l));
  double omega0=(system->k)/(system->m);
  for(k=0;k<no_iterations;k++)
  {
    for(i=0;i<(system->m);++i)
    {
      for(j=0;j<(system->n);++j)
      {
        accel=0;
        accel=-2*omega0*system->position[i][j]-2*omega1*system->position[i][j];
        if(i>0)
          accel+=omega1*system->position[i-1][j];
        if(i<(system->m)-1)
          accel+=omega1*system->position[i+1][j];
        if(j>0)
          accel+=omega0*system->position[i][j-1];
        if(j<(system->n)-1)
          accel+=omega0*system->position[i][j+1];
        temp_vel[i][j] = system->velocity[i][j] + accel*dt;
        temp_pos[i][j] = system->position[i][j] + temp_vel[i][j]*dt;
      }
    }  
    for(i=0;i<(system->m);++i)
    {
      for(j=0;j<(system->n);++j)
      {
        system->velocity[i][j]=temp_vel[i][j];
        system->position[i][j]=temp_pos[i][j];
      }
    }
  }
}

int main()
{
  temp_pos=create_matrix(M,N);
  temp_vel=create_matrix(M,N);
  FILE *f1,*f2,*f3,*f4;
  f1 = fopen("Ficheiros/file30_10.txt","wt");
  f2 = fopen("Ficheiros/file30_30.txt","wt");
  //f3 = fopen("Ficheiros/file23.txt","wt");
  //f4 = fopen("Ficheiros/file32.txt","wt");
  state grid;
  grid.position = create_matrix(M,N);
  grid.velocity = create_matrix(M,N);
  int i,j,k;
  for(i=0;i<M;i++)
  {
    for(j=0;j<N;++j)
    {
      grid.position[i][j]=1;
    }
  }
  grid.l=1;
  grid.f0=1;
  grid.k=1;
  grid.mass=1;
  grid.m=M;
  grid.n=N;
  for(k=0;k<2000;k++)
  {
  run_sim_square_xy(&grid,1000);
  /*for(i=0;i<10;++i)
  {
    for(j=0;j<10;++j)
    {
      printf("%g  ",grid.position[i][j]);
    }
    printf("\n");  
  }*/
  fprintf(f1,"%g  \n",grid.position[30][10]);
  fprintf(f2,"%g  \n",grid.position[30][30]);
  //fprintf(f3,"%g  \n",grid.position[2][3]);
  //fprintf(f4,"%g  \n",grid.position[0][4]);
  }
  fclose(f1);
  fclose(f2);
  //fclose(f3);
  //fclose(f4);
  return 0;
}
