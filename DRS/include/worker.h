#ifndef ALP_INCLUDE_WORKER
#define ALP_INCLUDE_WORKER

#include "index.h"
#include "DRS.h"
#include <pthread.h>
#include <atomic>

#include <iostream>

extern std::atomic<int> iter;
extern pthread_barrier_t barrier;

void async_worker(int thread_id, DRS drs, Parameter& param, Result& res)
{
	int max_iter = param.max_iter;
	int id;
	int epoch;
	pthread_barrier_wait(&barrier);
	if (!thread_id)
		res.time = get_wall_time();
	epoch = 0;

	for (int i = 0; ; i++, iter++)
	{
		if (param.stop || i > max_iter)
		{
			param.stop = true;
			break;
		}

		id = drs.get_id();
		drs.update(id);
		if (iter > (epoch+1) * param.check_step)
		{
			epoch++;
			pthread_barrier_wait(&barrier);
			if (!thread_id)
				res.check_time -= get_wall_time();
			drs.parallel_getmin(thread_id);
			pthread_barrier_wait(&barrier);
			drs.parallel_stop(thread_id);
			pthread_barrier_wait(&barrier);
			if (!thread_id)
			{
				drs.check_stop();
			}
			if (!thread_id)
				res.check_time += get_wall_time();
			pthread_barrier_wait(&barrier);
		}
	}
}


#endif // !ALP_INCLUDE_WORKER
