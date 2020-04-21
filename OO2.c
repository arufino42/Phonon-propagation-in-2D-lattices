
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define dt 1e-2
#define M 100
#define N 100
#define Mh 100
#define Nh 100

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

double delta=0.1, p=2;

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
        if  ((((j==1)||(j==system->n)) && (r1!=0)) || (i==1) || (i==system->m))
        {
    	  temp_velh[i][j] =0;
		  temp_posh[i][j] =0;
	    }
	    else
	    {  
	      if (r1!=r2)
	      {
		    accel-=3*system->position[i][j];
		    accel+=(system->position[i+1][j]+ system->position[i-1][j] + system->position[i][j+1]);
	      }
	      if (r1==r2)
	      {
		    accel-=3*system->position[i][j];
		    accel+=(system->position[i+1][j]+ system->position[i-1][j] + system->position[i][j-1]);
	      }
        }
        accel *= omega;
        temp_velh[i][j] = system->velocity[i][j] + accel*dt;
        temp_posh[i][j] = system->position[i][j] + temp_velh[i][j]*dt;
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

void run_sim_square_xyh(state *system,int no_iterations,double c,int sel)
{
  
  int i,j,k,r1,r2; 
  double accel;
  double omega=1/(double)(system->m);
  double o;
  if (sel == 0)
  	{
  	o=((system->f0)/(system->l));
  	}
  if (sel==1)
  {
  	o=system->k;
  	}
  for(k=0;k<no_iterations;k++)
  {
    for(i=1;i<=(system->m);++i)
    {
      for(j=1;j<=(system->n);++j)
      {
        accel=0;
		r1= i%2;
		r2= j%2;
		if ((((j==1)||(j==system->n)) && (r1!=0)) || (i==1) || (i==system->m))
	    {
		  temp_velh[i][j] =0;
		  temp_posh[i][j] =0;
		}
	    else
		{
		  if (r1!=r2)
	      {
		  	accel-=(2*c*(system->k)+o)*system->position[i][j];
		  	accel+=c*system->position[i+1][j]+ c*system->position[i-1][j] + o*system->position[i][j+1];
	      }
      	  if (r1==r2)
      	  {
			accel-=(2*c*(system->k)+o)*system->position[i][j];
		  	accel+=c*system->k*system->position[i+1][j]+ c*system->k*system->position[i-1][j] + o*system->position[i][j-1];
	      }
	  	}
        accel *= omega;
        temp_velh[i][j] = system->velocity[i][j] + accel*dt;
        temp_posh[i][j] = system->position[i][j] + temp_velh[i][j]*dt;
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

double grid_max(state *system)
{
  double max=0;
  int i,j;
  for(i=1;i<(system->m);++i)
  {
    for(j=1;j<(system->n);++j)
    {
      if(system->position[i][j]>max)
      {
        max=system->position[i][j];
      }
    }
  }
  return max;
}

double correct_x(int i, int j)
{
  double x;
  if(i%2==j%2)
  {
    x=j-1/(double)6;
  }
  else
  {
    x=j+1/(double)6;
  }
  return x*sqrt(3);
}

int main()
{
  temp_pos=create_matrix(M,N);
  temp_vel=create_matrix(M,N);

  temp_posh=create_matrix(Mh+1,Nh+1);
  temp_velh=create_matrix(Mh+1,Nh+1);


  FILE *f1,*f2,*f3,*f4, *f5, *f6, *f7, *f8;

  //Ficheiros
  
  state grid, gridxy;
  state gridh,gridhx,gridhy;

  //Hexagonos

  //z

  gridh.position = create_matrix(Mh+1,Nh+1);
  gridh.velocity = create_matrix(Mh+1,Nh+1);
  gridh.position[50][50]=p;
  gridh.l=1;
  gridh.f0=0.5;
  gridh.k=1;
  gridh.mass=1;
  gridh.m=Mh;
  gridh.n=Nh;

  //x

  gridhx.position = create_matrix(Mh+1,Nh+1);
  gridhx.velocity = create_matrix(Mh+1,Nh+1);
  gridhx.position[50][50]=p;
  gridhx.l=1;
  gridhx.f0=0.5;
  gridhx.k=1;
  gridhx.mass=1;
  gridhx.m=Mh;
  gridhx.n=Nh;

  //y

  gridhy.position = create_matrix(Mh+1,Nh+1);
  gridhy.velocity = create_matrix(Mh+1,Nh+1);
  gridhy.position[50][50]=p;
  gridhy.l=1;
  gridhy.f0=0.5;
  gridhy.k=1;
  gridhy.mass=1;
  gridhy.m=Mh;
  gridhy.n=Nh;

  //Quadrada

  //z
  
  grid.position = create_matrix(M,N);
  grid.velocity = create_matrix(M,N);
  grid.position[50][50]=p;
  grid.l=1;
  grid.f0=0.5;
  grid.k=1;
  grid.mass=1;
  grid.m=M;
  grid.n=N;

  //xy

  gridxy.position = create_matrix(M,N);
  gridxy.velocity = create_matrix(M,N);
  gridxy.position[50][50]=p;
  gridxy.l=1;
  gridxy.f0=0.5;
  gridxy.k=1;
  gridxy.mass=1;
  gridxy.m=M;
  gridxy.n=N;
  
  //f1=fopen("grafico_s_z.txt","wt");
  f2=fopen("grafico_s_y.txt","wt");
  //f3=fopen("grafico_h_z.txt","wt");
  //f4=fopen("grafico_h_x.txt","wt");
  //f5=fopen("grafico_h_y.txt","wt");
  int t;
  for(t=0;t<100;++t)
  {
    //run_sim_square_z(&grid,10e2);
    run_sim_square_xy(&gridxy,10e2);
    //run_sim_square_zh(&gridh,10e2);
    //run_sim_square_xyh(&gridhx,10e2,0.25,1);
    //run_sim_square_xyh(&gridhy,10e2,0.75,0);
    //fprintf(f1,"%lf\n",grid.position[25][50]);
    fprintf(f2,"%lf\n",gridxy.position[25][50]);
    //fprintf(f3,"%lf\n",gridh.position[25][50]);
    /*fprintf(f4,"%lf\n",gridhx.position[25][50]);
    fprintf(f5,"%lf\n",gridhy.position[25][50]);*/
    printf("%d\n",t);
  }
  
  //fclose(f1);
  fclose(f2);
  //fclose(f3);
  //fclose(f4);
  //fclose(f5);
  return 0;
}
