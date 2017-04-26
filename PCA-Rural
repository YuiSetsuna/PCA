#include "stdio.h"
#include "time.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "omp.h"

#define C 100
#define maxgen 2000 
#define pop 160
#define pcross 1 
#define pmutation 1 
#define MAXFIT 204.081633

int initial(int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C]);
void ffunction(double fitness[pop], int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C]);
void choice(double fitness[pop], int mvalue[pop][49+C]);
void cross(int mvalue[pop][49+C]);
void mutation(int mvalue[pop][49+C]);
void localsearch(double fitness[pop], int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C]);
int maxfitness(double fitness[pop], double *maxfitnessvalue);
double linefitness(int linenum, int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C]);
int userx[20];
int usery[20];

int initial(int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C])
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
	for(i=0;i<=49+C-1;i++)
	{
		xvalue[i]=rand()%286;
		yvalue[i]=rand()%286;
	}
	for(i=0;i<=19;i++)
	{
		userx[i]=rand()%286;
		usery[i]=rand()%286;
	}
	//for(i=0;i<=pop-1;i++)
	//{
		//for(j=0;j<=49+C-1;j++)
		//{
		//printf("mvalue[%d][%d]=%d\n",i,j,mvalue[i][j]);
		//}
	//}
	//for(i=0;i<49;i++)
	//{
		//printf("xvalue[%d]=%d,yvalue[%d]=%d\n",i,xvalue[i],i,yvalue[i]);
	//}
	printf("Initial Done\n");
}

void ffunction(double fitness[pop], int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C])
{
	printf("Fitness Start\n");
	int i=0;
	int j=0;
	int k=0;
	int temppop=0;
	int count=0;
	int countm=0;
	double countcover=0;
	double fxvalue=0;
	int tempcover[287][287];
	#pragma omp parallel for 
	for(temppop=0;temppop<=pop-1;temppop++)
	{
		fitness[temppop]=linefitness(temppop, mvalue, xvalue, yvalue);
	}
	//for(i=0;i<=pop-1;i++)
		//printf("fitness[%d]=%f\n",i,fitness[i]);
	printf("Fitness Done\n");
}

double linefitness(int linenum, int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C])
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
		for(i=0;i<=19;i++)
		{
				if(tempcover[userx[i]][usery[i]]==1)
				{
					count++;
				}	
		}
		//printf("count=%d\n",count);
		//printf("countm=%d\n",countm);
		countcover=(double)count;
		fxvalue=countcover*countcover/(double)countm;
	return fxvalue;
}

void choice(double fitness[pop], int mvalue[pop][49+C])
{
	printf("Choice Start\n");
	int i=0;
	int j=0;
	int k=0;
	double choice[pop][49+C];
	double fitprob[pop];
	double sum=0;
	double chooseprob=0;

	for(i=0;i<=pop-1;i++)
 	{
		sum+=fitness[i];
    	}
	//printf("sum=%f\n",sum);
	//sleep(1);
	for(i=0;i<=pop-1;i++)
	{
		fitprob[i]=fitness[i]/sum; //probability array
		//printf("fitprob[%d]=%f\n",i,fitprob[i]);
		//sleep(1);
	}
	for(i=0;i<=pop-1;i++)
	{
		chooseprob=((double)rand())/RAND_MAX; // random 0 to 1  
		for(j=0;j<=pop-1;j++)
        	{
					//printf("chooseprob1=%f\n",chooseprob);
            		chooseprob=chooseprob-fitprob[j];
            		//printf("chooseprob2=%f\n",chooseprob);
            		if(chooseprob<=0)
            		{
                		for(k=0;k<=49+C-1;k++)
                    			choice[i][k]=mvalue[j][k]; //choose
                    	//printf("Choose Done...\n");
                		break;    
            		}
        	}

    	}
	for(i=0;i<=pop-1;i++)
	{
		for(j=0;j<=49+C-1;j++)
            		mvalue[i][j]=choice[i][j];  //switch
    	}
    printf("Choice Done\n");
}

void cross(int mvalue[pop][49+C])
{
	printf("Cross Start\n");
	int i=0;
	int j=0;
	int parent0=0;
	int parent1=0;
	int temp=0;
	int forb=0;
	int crosspos=0;	//cross position
	int pos=0;	//population position
	double crossprob=0;

	srand((unsigned)time(NULL));

	while(pos<=pop-2)
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
	printf("Cross Done\n");
}

void mutation(int mvalue[pop][49+C])
{
	printf("Mutation Start\n");
	int i=0;	
	int mutationpos=0;
	int temp=0;
	double mutationprob;
	for(i=0;i<=pop-1;i++)
	{	
		mutationpos=rand()%(49+C); 
        mvalue[i][mutationpos]=(mvalue[i][mutationpos]+1)%2;
	}
	printf("Mutation Done\n");
}

void localsearch(double fitness[pop], int mvalue[pop][49+C], int xvalue[49+C], int yvalue[49+C])
{
	printf("Local Search Start\n");
	int i=0;
	int j=0;
	int k=0;
	int times=0;
	int temppop=0;
	int temptimes=6; 
	int moved=0;
	double tempfit[pop];
	double localfit[pop];
	for(i=0;i<=pop-1;i++)
		tempfit[i]=fitness[i];
	for(temppop=0;temppop<=temptimes-1;temppop++)
	{
		i=rand()%(49+C);
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
	printf("Local Search Done\n");
}

int maxfitness(double fitness[pop], double *maxfitnessvalue)
{
	printf("Max Start\n");
	int i=0;
	int index=0;
	double max=0;
	for(i=0;i<=pop-1;i++)
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
	printf("Max Done\n");
	return index;
}

int main()
{
	int xvalue[49+C];
	int yvalue[49+C];
	int mvalue[pop][49+C];
	double fitness[pop];
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
	double oldmaxvalue=0;
	double time;

	time_t start,finish;
	start = clock(); 

	initial(mvalue, xvalue, yvalue);
	ffunction(fitness, mvalue, xvalue, yvalue);

	oldindex=maxfitness(fitness, &maxfitnessvalue);
	oldmaxvalue=maxfitnessvalue;
	for(i=0;i<=49+C-1;i++)
		oldresult[i] = mvalue[oldindex][i]; 
	for(i=0;i<=maxgen-1;i++) 
	{
		printf("Generation:%d\n",i);
		//printf("oldmaxvalue=%f\n",oldmaxvalue);
		/*for(m=0;m<=pop-1;m++)
		{
		for(l=0;l<=49+C-1;l++)
		{
		printf("mvalue[%d][%d]=%d\n",m,l,mvalue[m][l]);
		}
		}*/
		//for(k=0;k<=pop-1;k++)
			//printf("fitness[%d]=%f\n",k,fitness[k]);
		//sleep(10);
		choice(fitness, mvalue); 
		cross(mvalue); 
		mutation(mvalue); 
		localsearch(fitness, mvalue, xvalue, yvalue);
		ffunction(fitness, mvalue, xvalue, yvalue);
        index=maxfitness(fitness, &maxfitnessvalue); 

        if(maxfitnessvalue>=oldmaxvalue)
		{
			//printf("%f>%f\n",maxfitnessvalue,oldmaxvalue);
			for(j=0;j<=49+C-1;j++)
			{
                		result[j]=mvalue[index][j]; 
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
				mvalue[index][j]=result[j];
			}
			bestfitgen=oldbestfitgen;
			maxfitnessvalue=oldmaxvalue;
		}
		printf("oldmaxvalue=%f,maxfitnessvalue=%f,index=%d\n",oldmaxvalue,maxfitnessvalue,index);
		if(fabs(maxfitnessvalue-MAXFIT)<1e-6)
			break;
	}
	finish=clock(); 
	time=((double)(finish-start))/CLOCKS_PER_SEC;
	for(i=0;i<=19;i++)
		printf("userx[%d]=%d,usery[%d]=%d\n",i,userx[i],i,usery[i]);
	for(i=0;i<=49+C-1;i++)
		{
			if(result[i]==1)
			{
				printf("xvalue[%d]=%d,yvalue[%d]=%d\n",i,xvalue[i],i,yvalue[i]);
			}
		}
	printf("time=%f,best fitness generation=%d,max fitness=%f\n",time,bestfitgen,maxfitnessvalue);
}
