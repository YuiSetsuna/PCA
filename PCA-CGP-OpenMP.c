#include "stdio.h"
#include "time.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "omp.h"

#define C 100
#define maxgen 2000 
#define pop 160
#define pcross 1 
#define pmutation 1 
#define MAXFIT 204.081633

int rank=0;
int size=0;
MPI_Status status;

int **allocmatrix(int n1, int n2);
void initial(int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C]);
void bcastmvalue(int mvalue[pop][49+C], int **localmvalue);
void ffunction(double fitness[], int **mvalue, int xvalue[49+C], int yvalue[49+C]);
void choice(double fitness[], int **mvalue);
void cross(int **mvalue);
void mutation(int **mvalue);
void localsearch(double fitness[], int **mvalue, int xvalue[49+C], int yvalue[49+C]);
int maxfitness(double fitness[], double *maxfitnessvalue);
double linefitness(int linenum, int **mvalue, int xvalue[49+C], int yvalue[49+C]);
void commmax(int **mvalue, double fitness[], double *maxfitnessvalue);
void freematrix(int **m);

int **allocmatrix(int n1, int n2)
{
    int **m;
  	int i=0;
    if ((m=(int **)calloc(n1,sizeof(int *))) == NULL) 
		return NULL;
    if ((m[0]=(int *)calloc(n1*n2,sizeof(int))) == NULL) 
		return NULL;
    for (i=1; i<n1; i++) 
		m[i]=m[i-1]+n2;
    return m;
}

void freematrix(int **m)
{
    free(m[0]);
    free(m);
}

void initial(int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C])
{
  if(rank==0)
  {
	printf("Initial Start\n");
	int i=0;
	int j=0;
	srand((unsigned)time(NULL));
	#pragma omp parallel for 
	for(i=0;i<=pop-1;i++)
	{
		for(j=0;j<=49+C-1;j++)
		{
		mvalue[i][j]=rand()%2;
		//mvalue[i][j]=1;
		}
	}
	#pragma omp parallel for 
	for(i=0;i<49;i++)
	{
		xvalue[i]=(i%7)*41+20;
		yvalue[i]=(i/7)*41+20;
	}
	#pragma omp parallel for 
	for(i=49;i<=49+C-1;i++)
	{
		xvalue[i]=rand()%286;
		yvalue[i]=rand()%286;
	}
  printf("Initial Done\n");
  }
  MPI_Bcast(xvalue,49+C,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(yvalue,49+C,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
}

void bcastmvalue(int mvalue[pop][49+C], int **localmvalue)
{
  MPI_Bcast(mvalue,pop*(49+C),MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  
  int i=0;
  int j=0;
  int k=0;
  
  int quotient=pop/size;
  int remainder=pop%size;

  if(rank<remainder)
  { 
    for(i=0,j=rank*(quotient+1);(i<=quotient+1-1)&&(j<=(rank+1)*(quotient+1)-1);i++,j++)
    {
		for(k=0;k<=49+C-1;k++)
		{
			localmvalue[i][k]=mvalue[j][k];
		}	
    }
  }
  else
  {
	for(i=0,j=rank*quotient+remainder;(i<=quotient-1)&&(j<=(rank+1)*quotient+remainder-1);i++,j++)
    {
      	for(k=0;k<=49+C-1;k++)
		{
			localmvalue[i][k]=mvalue[j][k];
		}
    }
  }
}

void ffunction(double fitness[], int **mvalue, int xvalue[49+C], int yvalue[49+C])
{
	//printf("Fitness Start\n");
	int i=0;
	int j=0;
	int k=0;
	int temppop=0;
	int count=0;
	int countm=0;
	double countcover=0;
	double fxvalue=0;
	int tempcover[287][287];

	int quotient=pop/size;
	int remainder=pop%size;

  if(rank<remainder)
  {
	#pragma omp parallel for 
    for(temppop=0;temppop<=quotient+1-1;temppop++)
    {
      fitness[temppop]=linefitness(temppop, mvalue, xvalue, yvalue);
    }
  }
  else
  {
	#pragma omp parallel for 
    for(temppop=0;temppop<=quotient-1;temppop++)
    {
      fitness[temppop]=linefitness(temppop, mvalue, xvalue, yvalue);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  /*if(rank==0)
  {
  	for(temppop=0;temppop<=pop-1;temppop++)
	{
		fitness[temppop]=linefitness(temppop, mvalue, xvalue, yvalue);
	}
  }*/
	//for(i=0;i<=pop-1;i++)
		//printf("fitness[%d]=%f\n",i,fitness[i]);
	//printf("Fitness Done\n");
}

double linefitness(int linenum, int **mvalue, int xvalue[49+C], int yvalue[49+C])
{
	int i=0;
	int j=0;
	int k=0;
	int tempcover[287][287];
	int count=0;
	int countm=0;
	double countcover=0;
	double fxvalue=0;
	for(i=0;i<=286;i++)
		{
			for(j=0;j<=286;j++)
			{
				tempcover[i][j]=0;
			}
		}
		for(i=0;i<=49+C-1;i++)
		{
			if(mvalue[linenum][i]==1)
			{
			countm++;
			//printf("temppop=%d,i=%d\n",temppop,i);
			if(xvalue[i]<20&&yvalue[i]<20)
			{
			//printf("xvalue=%d<20,yvalue=%d<20\n",xvalue[i],yvalue[i]);
			//sleep(0.2);
			for(j=0;j<=xvalue[i]+20;j++)
			{
				for(k=0;k<=yvalue[i]+20;k++)
				{
					tempcover[j][k]=1;
				}
			}
			}
			else if(xvalue[i]<20&&yvalue[i]>266)
			{
			//printf("xvalue=%d<20,yvalue=%d>266\n",xvalue[i],yvalue[i]);
			//sleep(0.2);
			for(j=0;j<=xvalue[i]+20;j++)
			{
				for(k=yvalue[i]-20;k<=286;k++)
				{
					tempcover[j][k]=1;
				}
			}
			}
			else if(xvalue[i]<20&&yvalue[i]>=20&&yvalue[i]<=266)
			{
			//printf("xvalue=%d<20,yvalue=%d in [20,266]\n",xvalue[i],yvalue[i]);
			//sleep(0.2);
			for(j=0;j<=xvalue[i]+20;j++)
			{
				for(k=yvalue[i]-20;k<=yvalue[i]+20;k++)
				{
					tempcover[j][k]=1;
				}
			}
			}
			else if(xvalue[i]>266&&yvalue[i]<20)
			{
			//printf("xvalue=%d>266,yvalue=%d<20\n",xvalue[i],yvalue[i]);
			//sleep(0.2);
			for(j=xvalue[i]-20;j<=286;j++)
			{
				for(k=0;k<=yvalue[i]+20;k++)
				{
					tempcover[j][k]=1;
				}
			}
			}
			else if(xvalue[i]>266&&yvalue[i]>266)
			{
			//printf("xvalue=%d>266,yvalue=%d>266\n",xvalue[i],yvalue[i]);
			//sleep(0.2);
			for(j=xvalue[i]-20;j<=286;j++)
			{
			for(k=yvalue[i]-20;k<=286;k++)
			{
				tempcover[j][k]=1;
			}
			}
			}
			else if(xvalue[i]>266&&yvalue[i]>=20&&yvalue[i]<=266)
			{
			//printf("xvalue=%d>266,yvalue=%d\n",xvalue[i],yvalue[i]);
			//sleep(0.2);
			for(j=xvalue[i]-20;j<=286;j++)
			{	
				for(k=yvalue[i]-20;k<=yvalue[i]+20;k++)
				{
					tempcover[j][k]=1;
				}
			}
			}
			else if(yvalue[i]>266&&xvalue[i]>=20&&xvalue[i]<=266)
			{
			//printf("xvalue=%d in [20,266],yvalue=%d>266\n",xvalue[i],yvalue[i]);
			//sleep(0.2);/
			for(j=xvalue[i]-20;j<=xvalue[i]+20;j++)
			{
				for(k=yvalue[i]-20;k<=286;k++)
				{
					tempcover[j][k]=1;
				}
			}
			}
			else if(yvalue[i]<20&&xvalue[i]>=20&&xvalue[i]<=266)
			{
			//printf("xvalue=%d in [20,266],yvalue=%d<20\n",xvalue[i],yvalue[i]);
			//sleep(0.2);
			for(j=xvalue[i]-20;j<=xvalue[i]+20;j++)
			{
				for(k=0;k<=yvalue[i]+20;k++)
				{
					tempcover[j][k]=1;
				}
			}
			}
			else
			{
			//printf("xvalue=%d in [20,266],yvalue=%d in [20,266]\n",xvalue[i],yvalue[i]);
			//sleep(0.2);
			for(j=xvalue[i]-20;j<=xvalue[i]+20;j++)
			{
				for(k=yvalue[i]-20;k<=yvalue[i]+20;k++)
				{
					tempcover[j][k]=1;

				}
			}
			}
			}		
		}
		for(i=0;i<=286;i++)
		{
			for(j=0;j<=286;j++)
			{
				if(tempcover[i][j]==1)
				{
					count++;
				}	
			}
		}
		//printf("count=%d\n",count);
		//printf("countm=%d\n",countm);
		countcover=(double)count/82369.0;
		countcover*=100.0;
		fxvalue=countcover*countcover/(double)countm;
	return fxvalue;
}

void choice(double fitness[], int **mvalue)
{
	//printf("Choice Start\n");
	int i=0;
	int j=0;
	int k=0;
	double sum=0;
	double chooseprob=0;

  int **choice;
  double *fitprob;

	int quotient=pop/size;
	int remainder=pop%size;
	
	if(rank<remainder)
	{
		choice=allocmatrix(quotient+1,49+C); 
		fitprob=(double*)calloc(quotient+1,sizeof(double));
		for(i=0;i<=quotient+1-1;i++)
		{
			sum+=fitness[i];
    	}
		for(i=0;i<=quotient+1-1;i++)
		{
			fitprob[i]=fitness[i]/sum; //probability array
			//printf("fitprob[%d]=%f\n",i,fitprob[i]);
			//sleep(1);
		}
		for(i=0;i<=quotient+1-1;i++)
		{
			chooseprob=((double)rand())/RAND_MAX; // random 0 to 1  
			for(j=0;j<=quotient+1-1;j++)
			{
            	chooseprob=chooseprob-fitprob[j];
            	if(chooseprob<=0)
            	{
                	for(k=0;k<=49+C-1;k++)
                    choice[i][k]=mvalue[j][k]; //choose
                	break;    
            	}
        	}

		}
		for(i=0;i<=quotient+1-1;i++)
		{
			for(j=0;j<=49+C-1;j++)
				mvalue[i][j]=choice[i][j];  //switch
		}
	}
	else
	{
		choice=allocmatrix(quotient,49+C); 
		fitprob=(double*)calloc(quotient,sizeof(double));
		for(i=0;i<=quotient-1;i++)
		{
			sum+=fitness[i];
    	}
		for(i=0;i<=quotient-1;i++)
		{
			fitprob[i]=fitness[i]/sum; //probability array
			//printf("fitprob[%d]=%f\n",i,fitprob[i]);
			//sleep(1);
		}
		for(i=0;i<=quotient-1;i++)
		{
			chooseprob=((double)rand())/RAND_MAX; // random 0 to 1  
			for(j=0;j<=quotient-1;j++)
			{
            	chooseprob=chooseprob-fitprob[j];
            	if(chooseprob<=0)
            	{
                	for(k=0;k<=49+C-1;k++)
                    choice[i][k]=mvalue[j][k]; //choose
                	break;    
            	}
        	}

		}
		for(i=0;i<=quotient-1;i++)
		{
			for(j=0;j<=49+C-1;j++)
				mvalue[i][j]=choice[i][j];  //switch
		}
	}
  freematrix(choice);
  free(fitprob);
  //printf("Choice Done\n");
}

void cross(int **mvalue)
{
	//printf("Cross Start\n");
	int i=0;
	int j=0;
	int parent0=0;
	int parent1=0;
	int temp=0;
	int forb=0;
	int crosspos=0;	//cross position
	int pos=0;	//population position
	int localpop=0;
	double crossprob=0;

	srand((unsigned)time(NULL));

	int quotient=pop/size;
	int remainder=pop%size;
	
	if(rank<remainder)
	{
		localpop=quotient-1;
	}
	else
	{
		localpop=quotient;
	}
	while(pos<=localpop-2)
	{
        parent0 = pos; 
        parent1 = pos+1; 
	
		crosspos=rand()%(49+C);
		forb=rand()%2;
		//printf("crosspos=%d,forb=%d\n",crosspos,forb);
		if(forb==0)
		{
			//printf("Switch Front\n");
			for(i=0;i<=crosspos;i++)
        		{
				temp=mvalue[parent0][i];
				mvalue[parent0][i]=mvalue[parent1][i];
				mvalue[parent1][i]=temp; 
        		}
        	//printf("Switch Done\n");
		}
		else
		{
			//printf("Switch Behind\n");
			for(i=crosspos;i<=49+C-1;i++)
        		{
				temp=mvalue[parent0][i];
				mvalue[parent0][i]=mvalue[parent1][i];
				mvalue[parent1][i]=temp; 
        		}
        	//printf("Switch Done\n");
		}
        pos+=2;
        //printf("pos=%d\n",pos);
	}
	//printf("Cross Done\n");
}

void mutation(int **mvalue)
{
	//printf("Mutation Start\n");
	int i=0;	
	int mutationpos=0;
	int temp=0;
  int localpop=0;
	double mutationprob;
	
	int quotient=pop/size;
	int remainder=pop%size;
	
	if(rank<remainder)
	{
		localpop=quotient-1;
	}
	else
	{
		localpop=quotient;
	}
	for(i=0;i<=localpop-1;i++)
	{	
		mutationpos=rand()%(49+C); 
        mvalue[i][mutationpos]=(mvalue[i][mutationpos]+1)%2;
	}
	//printf("Mutation Done\n");
}

void localsearch(double fitness[], int **mvalue, int xvalue[49+C], int yvalue[49+C])
{
	//printf("Local Search Start\n");
	int i=0;
	int j=0;
	int k=0;
	int times=0;
	int temppop=0;
	int temptimes=6; 
	int moved=0;
  int localpop=0;
	double tempfit[pop];
	double localfit[pop];
	
	int quotient=pop/size;
	int remainder=pop%size;
	
	if(rank<remainder)
	{
		localpop=quotient-1;
	}
	else
	{
		localpop=quotient;
	}
	
	for(i=0;i<=localpop-1;i++)
		tempfit[i]=fitness[i];
	for(temppop=0;temppop<=temptimes-1;temppop++)
	{
		i=rand()%localpop;
		//j=rand()%(49+C);
		for(j=0;j<=49+C-1;j++)
		{
		if(mvalue[i][j]==1)
		{
			moved=0;
			times=0;
			for(k=j+1;k<=49+C-1;k++)
			{
				if(mvalue[i][k]==0)
				{
					if(times>=3)
						break;
					times++;
					mvalue[i][k]=1;
					mvalue[i][j]=0;
					localfit[i]=linefitness(i, mvalue, xvalue, yvalue);
					if(localfit[i]>tempfit[i])
					{
						tempfit[i]=localfit[i];
						fitness[i]=localfit[i]; 
						moved++;
					}
					else
					{
						mvalue[i][k]=0;
						mvalue[i][j]=1;
						moved=0;
					}
				}
			}
			if(!moved)
			{
				mvalue[i][j]=0;
				localfit[i]=linefitness(i, mvalue, xvalue, yvalue);
				if(localfit[i]>tempfit[i])
				{
					tempfit[i]=localfit[i];
					fitness[i]=localfit[i];
				}
				else
					mvalue[i][j]=1;
			}
		}
		}
	}
	//printf("Local Search Done\n");
}

int maxfitness(double fitness[], double *maxfitnessvalue)
{
	//printf("Max Start\n");
	int i=0;
	int index=0;
  int localpop=0;
	double max=0;
	
	int quotient=pop/size;
	int remainder=pop%size;
	
	if(rank<remainder)
	{
		localpop=quotient-1;
	}
	else
	{
		localpop=quotient;
	}
	
	for(i=0;i<=localpop-1;i++)
	{
		//printf("fitness[%d]=%f\n",i,fitness[i]);
		if(fitness[i]>max)
		{
			max=fitness[i];
			index=i;
		}
	}
	*maxfitnessvalue=max;
	
	//sleep(1);
	//printf("Max Done\n");
	return index;
}

int minfitness(double fitness[])
{
	//printf("Min Start\n");
	int i=0;
	int index=0;
  int localpop=0;
	double min=0;
	
	int quotient=pop/size;
	int remainder=pop%size;
	
	if(rank<remainder)
	{
		localpop=quotient-1;
	}
	else
	{
		localpop=quotient;
	}
	
	min=fitness[0];
	index=0;
	
	for(i=1;i<=localpop-1;i++)
	{
		//printf("fitness[%d]=%f\n",i,fitness[i]);
		if(fitness[i]<min)
		{
			min=fitness[i];
			index=i;
		}
	}
	//printf("Min Done\n");
	return index;
}

void commmax(int **mvalue, double fitness[], double *maxfitnessvalue)
{
	int minindex=0;
	int maxindex=0;
  int localpop=0;
  
	int quotient=pop/size;
	int remainder=pop%size;
	
	if(rank<remainder)
	{
		localpop=quotient+1;
	}
	else
	{
		localpop=quotient;
	}
	
	minindex=minfitness(fitness);
	maxindex=maxfitness(fitness, maxfitnessvalue);

	MPI_Sendrecv(&mvalue[maxindex][0], 49+C, MPI_INT, (rank+1)%size, 0, &mvalue[minindex][0], 49+C, MPI_INT, (rank+size-1)%size, 0, MPI_COMM_WORLD, &status);
  MPI_Barrier(MPI_COMM_WORLD);
	MPI_Sendrecv(&fitness[maxindex], 1, MPI_DOUBLE, (rank+1)%size, 0, &fitness[minindex], 1, MPI_DOUBLE, (rank+size-1)%size, 0, MPI_COMM_WORLD, &status);
	MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
	MPI_Comm_size(MPI_COMM_WORLD,&size);
  
	int xvalue[49+C];
	int yvalue[49+C];
	int mvalue[pop][49+C];
	double maxfitnessvalue=0;

	int i=0;
	int j=0;
	int k=0;
	int m=0;
	int l=0;
	int index=0;
	int bestfitgen=0;
	int oldbestfitgen=0;
	int result[49+C];
	int oldresult[49+C];
	int oldindex;
	int exitrank=0;
	int exit=0;
  int **localmvalue;
  double *fitness;
	double oldmaxvalue=0;
	double time;

	time_t start,finish;
	start = clock(); 

	int quotient=pop/size;
	int remainder=pop%size;

	if(rank<remainder)
	{
		localmvalue=allocmatrix(quotient+1,49+C); 
		fitness=(double*)calloc(quotient+1,sizeof(double));
	}
	else
	{
		localmvalue=allocmatrix(quotient,49+C); 
		fitness=(double*)calloc(quotient,sizeof(double));
	}
	
	initial(mvalue, xvalue, yvalue);
	bcastmvalue(mvalue,localmvalue);
	ffunction(fitness, localmvalue, xvalue, yvalue);

	oldindex=maxfitness(fitness, &maxfitnessvalue);
	oldmaxvalue=maxfitnessvalue;
	for(i=0;i<=49+C-1;i++)
		oldresult[i] = localmvalue[oldindex][i]; 
	for(i=0;i<=maxgen-1;i++) 
	{
		if(rank==0)
			printf("Generation:%d\n",i);
		choice(fitness, localmvalue); 
		cross(localmvalue); 
		mutation(localmvalue); 
		localsearch(fitness, localmvalue, xvalue, yvalue);
		ffunction(fitness, localmvalue, xvalue, yvalue);
		if((i+1)%20==0)
			commmax(localmvalue, fitness, &maxfitnessvalue);
		index=maxfitness(fitness, &maxfitnessvalue);
		if(maxfitnessvalue>=oldmaxvalue)
		{
			//printf("%f>%f\n",maxfitnessvalue,oldmaxvalue);
			for(j=0;j<=49+C-1;j++)
			{
				result[j]=localmvalue[index][j]; 
				oldresult[j]=result[j];
			}
			bestfitgen=i+1; 
			oldindex=index;
			oldmaxvalue=maxfitnessvalue;
			oldbestfitgen=bestfitgen;
		}
		else
		{
			//printf("%f<%f\n",maxfitnessvalue,oldmaxvalue);
			for(j=0;j<=49+C-1;j++)
			{
				result[j] = oldresult[j]; 
				localmvalue[index][j]=result[j];
			}
			bestfitgen=oldbestfitgen;
			maxfitnessvalue=oldmaxvalue;
		}
		if(rank==0)
			printf("oldmaxvalue=%f,maxfitnessvalue=%f,index=%d\n",oldmaxvalue,maxfitnessvalue,index);
		if(fabs(maxfitnessvalue-MAXFIT)<1e-6)
		{
			exitrank=rank;
			exit=1;
		}
		MPI_Bcast(&exitrank,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&exit,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if(exit==1)
			break;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==exitrank)
		printf("Evolution Done\n");
		finish=clock();
		time=((double)(finish-start))/CLOCKS_PER_SEC;
	if(rank==exitrank)
	{
		for(i=0;i<=49+C-1;i++)
			printf("result[%d]=%d\n",i,result[i]);
		printf("time=%f,best fitness generation=%d,max fitness=%f\n",time,bestfitgen,maxfitnessvalue);
	}
  freematrix(localmvalue);
  free(fitness);
	MPI_Finalize();
}
