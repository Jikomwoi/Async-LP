#ifndef ALP_INCLUDE_WORKER
#define ALP_INCLUDE_WORKER

#include "index.h"
#include "primdual.h"
#include <pthread.h>
#include <atomic>
#include "util.h"
#include <iostream>

extern std::atomic<int> iter;
extern pthread_barrier_t barrier;

void async_worker(int thread_id, PD pd, Parameter& param, Result& res)
{
	int max_iter = param.max_iter;
	int id;
	int epoch = 0;

	pthread_barrier_wait(&barrier);
	if (!thread_id)
	{
		res.time = get_wall_time();
		res.time_without = 0;
	}

	for (int i = 0; ; i++, iter++)
	{
		if (param.stop || i > max_iter)
		{
			param.stop = true;
			break;
		}

		id = pd.get_id();
		pd.update(id, epoch+1);
		// stop check
		if (iter > (epoch+1)*param.check_step)
		{
			double temp;
			epoch ++;
			pthread_barrier_wait(&barrier);
			if (!thread_id)
			{
				temp = get_wall_time();
			}
			pd.parallel_stop(thread_id);
			pthread_barrier_wait(&barrier);
			if (!thread_id)
			{
				pd.check_stop();
				res.time_without += get_wall_time() - temp;
			}
			pthread_barrier_wait(&barrier);
		}
	}
}


#endif // !ALP_INCLUDE_WORKER
