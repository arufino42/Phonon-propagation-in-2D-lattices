
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dt 1e-4
#define M 41
#define N 41
#define Mh 40
#define Nh 40

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

double **temp_pos,**temp_vel,**temp_posh,**temp_velh;



//Criação ficheiros e matrizes

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




// Simulações Rede Quadrada

//z

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

//x y

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



//Simulações Rede Hexagonal


//z

void run_sim_square_zh(state *system,int no_iterations)
{
  
  int i,j,k,r1,r2; 
  double accel;
  double omega=(system->f0)/((system->m)*(system->l));
  for(k=0;k<no_iterations;k++)
  {
    for(i=1;i<=(system->m);++i)
    {
      for(j=1;j<=(system->n);++j)
      {
        accel=0;

	r1= i%2;
	r2= j%2;
	if (r1!=r2)
	  {
	    if  ((((j==1)||(j==system->n)) && (r1!=0)) || (i==1) || (i==system->m))
	      {
		temp_vel[i][j] =0;
		temp_pos[i][j] =0;
	      }
	    else 
	      {
		accel-=3*system->position[i][j];
		accel+=(system->position[i+1][j]+ system->position[i-1][j] + system->position[i][j+1]);
	      }
	  }
      	if (r1==r2)
	  {
	    if ((((j==1)||(j==system->n)) && (r1!=0)) || (i==1) || (i==system->m))
	      {
		temp_vel[i][j] =0;
		temp_pos[i][j] =0;
	      }
	    else 
	      {
		accel-=3*system->position[i][j];
		accel+=(system->position[i+1][j]+ system->position[i-1][j] + system->position[i][j-1]);
	      }
	  }
	 
        accel *= omega;
        temp_velh[i][j] = system->velocity[i][j] + accel*dt;
        temp_posh[i][j] = system->position[i][j] + temp_vel[i][j]*dt;
      }
    }  
    for(i=1;i<=(system->m);++i)
    {
      for(j=1;j<=(system->n);++j)
      {
        system->velocity[i][j]=temp_velh[i][j];
        system->position[i][j]=temp_posh[i][j];
      }
    }
  }
}

//xy

void run_sim_square_xyh(state *system,int no_iterations,double c)
{
  
  int i,j,k,r1,r2; 
  double accel;
  double omega=c/((system->k)*(system->m));
  for(k=0;k<no_iterations;k++)
  {
    for(i=1;i<=(system->m);++i)
    {
      for(j=1;j<=(system->n);++j)
      {
        accel=0;

	r1= i%2;
	r2= j%2;
	if (r1!=r2)
	  {
	    if  ((((j==1)||(j==system->n)) && (r1!=0)) || (i==1) || (i==system->m))
	      {
		temp_vel[i][j] =0;
		temp_pos[i][j] =0;
	      }
	    else 
	      {
		accel-=3*system->position[i][j];
		accel+=(system->position[i+1][j]+ system->position[i-1][j] + system->position[i][j+1]);
	      }
	  }
      	if (r1==r2)
	  {
	    if ((((j==1)||(j==system->n)) && (r1!=0)) || (i==1) || (i==system->m))
	      {
		temp_vel[i][j] =0;
		temp_pos[i][j] =0;
	      }
	    else 
	      {
		accel-=3*system->position[i][j];
		accel+=(system->position[i+1][j]+ system->position[i-1][j] + system->position[i][j-1]);
	      }
	  }
	 
        accel *= omega;
        temp_velh[i][j] = system->velocity[i][j] + accel*dt;
        temp_posh[i][j] = system->position[i][j] + temp_vel[i][j]*dt;
      }
    }  
    for(i=1;i<=(system->m);++i)
    {
      for(j=1;j<=(system->n);++j)
      {
        system->velocity[i][j]=temp_velh[i][j];
        system->position[i][j]=temp_posh[i][j];
      }
    }
  }
}







int main()
{
  temp_pos=create_matrix(M,N);
  temp_vel=create_matrix(M,N);

  temp_posh=create_matrix(Mh+1,Nh+1);
  temp_velh=create_matrix(Mh+1,Nh+1);


  FILE *f1,*f2,*f3,*f4, *f5, *f6, *f7, *f8;

  //Ficheiros Quadrados

  //z
  
  f1 = fopen("Ficheiros/file00.txt","wt");
  f2 = fopen("Ficheiros/file44.txt","wt");
  f3 = fopen("Ficheiros/file23.txt","wt");
  f4 = fopen("Ficheiros/file32.txt","wt");

  //xy

  f8=fopen("Ficheiros/filexy.txt","wt");

  
  //Ficheiros hexagonos
  
  f5 = fopen("Ficheiros/filehz.txt","wt");
  f6 = fopen("Ficheiros/filehx.txt","wt");
  f7 = fopen("Ficheiros/filehy.txt","wt");

  
  state grid, gridxy;
  state gridh,gridhx,gridhy;

  //Hexagonos

  //z

  gridh.position = create_matrix(Mh+1,Nh+1);
  gridh.velocity = create_matrix(Mh+1,Nh+1);
  gridh.velocity[0][0]=1;
  gridh.l=1;
  gridh.f0=1;
  gridh.k=1;
  gridh.mass=1;
  gridh.m=Mh;
  gridh.n=Nh;

  //x

  gridhx.position = create_matrix(Mh+1,Nh+1);
  gridhx.velocity = create_matrix(Mh+1,Nh+1);
  gridhx.velocity[2][2]=1;
  gridhx.l=1;
  gridhx.f0=1;
  gridhx.k=1;
  gridhx.mass=1;
  gridhx.m=Mh;
  gridhx.n=Nh;

  //y

  gridhy.position = create_matrix(Mh+1,Nh+1);
  gridhy.velocity = create_matrix(Mh+1,Nh+1);
  gridhy.velocity[2][2]=1;
  gridhy.l=1;
  gridhy.f0=1;
  gridhy.k=1;
  gridhy.mass=1;
  gridhy.m=Mh;
  gridhy.n=Nh;

  //Quadrada

  //z
  
  grid.position = create_matrix(M,N);
  grid.velocity = create_matrix(M,N);
  grid.velocity[2][2]=1;
  grid.l=1;
  grid.f0=1;
  grid.k=1;
  grid.mass=1;
  grid.m=M;
  grid.n=N;

  //xy

  gridxy.position = create_matrix(M,N);
  gridxy.velocity = create_matrix(M,N);
  gridxy.velocity[2][2]=1;
  gridxy.l=1;
  gridxy.f0=1;
  gridxy.k=1;
  gridxy.mass=1;
  gridxy.m=M;
  gridxy.n=N;
  
  int i,j,k;

  
  for(k=0;k<2000;k++)
  {
  run_sim_square_z(&grid,1000);
  run_sim_square_xy(&gridxy,1000);
  run_sim_square_zh(&gridh,1000);
  run_sim_square_xyh(&gridhx,1000,1/4);
  run_sim_square_xyh(&gridhy,1000,3/4);

  fprintf(f5,"%g  \n",gridh.position[2][2]);
  fprintf(f6,"%g  \n",gridhx.position[2][2]);
  fprintf(f7,"%g  \n",gridhy.position[2][2]);


  fprintf(f1,"%g  \n",grid.position[0][0]);
  fprintf(f2,"%g  \n",grid.position[4][4]);
  fprintf(f3,"%g  \n",grid.position[2][3]);
  fprintf(f4,"%g  \n",grid.position[0][4]);

  fprintf(f8,"%g  \n",gridxy.position[2][2]);

  }

  fclose(f1);
  fclose(f2);

  return 0;
}
