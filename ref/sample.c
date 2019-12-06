#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h>
#include <mpi.h>

void MPI_Get_stealing_profile_data(int **, int **);

static void thread_sample_data(){
	int i;
	const int cores = 36;
	long long timestep = 0;
	int *priv_shm_idle_process;
	int *priv_shm_avail_tasks;
	int total_idle_process;
	int total_avail_tasks;

	MPI_Get_stealing_profile_data(&priv_shm_idle_process, &priv_shm_avail_tasks);

	FILE *fp = fopen("timestep.log","w");
	do{
		total_avail_tasks = total_idle_process = 0;
		for(i=0;i<cores;++i){
			total_idle_process += priv_shm_idle_process[i];
			total_avail_tasks += priv_shm_avail_tasks[i];
		}
		fprintf(fp, "%lld %d %d\n", timestep, total_idle_process, total_avail_tasks);
		usleep(100);
		timestep += 100;
	} while(1);
}

void sample_process_and_task_(){
	int err;
	pthread_attr_t attr;
	pthread_t thread_id;
	pthread_attr_init(&attr);
	err = pthread_create(&thread_id, &attr, &thread_sample_data, NULL);
	if(err != 0){
		printf("pthread create error\n");
		exit(1);
	}
	pthread_attr_destroy(&attr);
	return;
}
