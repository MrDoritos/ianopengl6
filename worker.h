#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <deque>
#include <algorithm>

#define TASK_COMPLETED 1
#define TASK_NOT_STARTED 0
#define TASK_FAILED 2
#define TASK_RUNNING 3

typedef void *voidarg;
typedef void (*voidfunc)(voidarg);

struct task {
    task() {
        function = 0;
        argument = 0;
        status = TASK_NOT_STARTED;
    }
    task(voidfunc func, voidarg arg)
    :task() {
        this->function = func;
        this->argument = arg;
    }
    voidfunc function;
    voidarg argument;
    int status;
};

#define WORKER_NOT_STARTED 0
#define WORKER_RUNNING 1
#define WORKER_STOPPING 2

void worker_function();
void worker_sleep(int milli);

struct worker {
    worker() {
        threads = 1;
        join_timeout = 500;
        state = WORKER_NOT_STARTED;
        INSTANCE = this;
    }
    worker(int thread_count) 
        :worker() {
            this->threads = thread_count;
    }

    void stop() {
        puts("Stopping worker");
        {
            std::unique_lock<std::mutex> lock(queue_lock);
            puts("Obtained lock");

            state = WORKER_STOPPING;
            task_queue.empty();
            puts("Emptied task queue and set state to stopping... waiting on threads to finish");
        }
        puts("Unlock queue");
        printf("Joining %i threads after %i ms\n", threads, join_timeout);
        //std::this_thread::sleep_for(std::chrono::milliseconds(join_timeout));
        worker_sleep(join_timeout);
        for (std::thread *thread : workers) {
            if (!thread->joinable())
                puts("Thread not joinable");
        }
        puts("Joining all threads");
        for (std::thread *thread : workers) {
            thread->join();
            puts("Deleting");
            delete thread;
        }
    }

    void start() {
        if (threads < 1)
            puts("Warning: threads not greater than 0");
        puts("Set state to running");
        state = WORKER_RUNNING;
        for (int i = 0; i < threads; i++) {
            printf("Creating thread %i\n", i);
            std::thread *new_worker = new std::thread(worker_function);
            workers.push_back(new_worker);
        }
    }

    void add_task(task* newtask) {
        {
            std::unique_lock<std::mutex> lock(queue_lock);
            //if (std::find(task_queue.begin(), task_queue.end(), newtask) == task_queue.end()) {
                task_queue.push_back(newtask);
            //} else {
            //    printf("This task %p exists\n", newtask);
            //}
        }
    }

    int state;
    int threads;
    int join_timeout;
    std::mutex queue_lock;
    std::deque<task*> task_queue;
    std::vector<std::thread*> workers;
    static worker *INSTANCE;
};

worker *worker::INSTANCE;

void worker_sleep(int milli) {
    std::this_thread::sleep_for(std::chrono::milliseconds(milli));
}

void worker_function() {
    const int wait_timeout = 10;
    while (worker::INSTANCE->state == WORKER_RUNNING) {
        task *new_task = 0;
        if (worker::INSTANCE->task_queue.size() > 0) {
            {
                std::unique_lock<std::mutex> lock(worker::INSTANCE->queue_lock);
                if (worker::INSTANCE->task_queue.size() > 0) {
                    new_task = worker::INSTANCE->task_queue.front();
                    //printf("New Task: %p\n", new_task);
                    //worker::INSTANCE->task_queue.pop_front();
                    worker::INSTANCE->task_queue.erase(worker::INSTANCE->task_queue.begin());
                }
            }
        }

        if (new_task && new_task->function) {
            new_task->status = TASK_RUNNING;
            new_task->function(new_task->argument);
            new_task->status = TASK_COMPLETED;
            //printf("Finished Task: %p\n", new_task);
        } else {
            worker_sleep(wait_timeout);
        }
    }

    puts("Worker stopped");
}