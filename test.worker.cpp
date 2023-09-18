#include "worker.h"

struct do_some_math_args {
    int start;
    int skip;
    int end;
};

void do_some_math(voidarg void_math_args) {
    struct do_some_math_args *math_args = (do_some_math_args*)void_math_args;
    for (int i = math_args->start; i < math_args->end; i += math_args->skip) {
        printf("i = %i\n", i);
    }
}


int main() {
    worker new_worker(2);
    do_some_math_args a, b;
    a.start = 10; a.skip = 10; a.end = 500;
    b.start = -5; b.skip = 3; b.end = 50;
    new_worker.start();
    
    task _a, _b;
    _a.argument = (voidarg)&a;
    _b.argument = (voidarg)&b;
    _a.function = (voidfunc)&do_some_math;
    _b.function = (voidfunc)&do_some_math;

    new_worker.add_task(&_a);
    new_worker.add_task(&_b);

    worker_sleep(1000);

    new_worker.stop();

    return 0;
}