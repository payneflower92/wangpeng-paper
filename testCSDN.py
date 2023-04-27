import math
import datetime
import multiprocessing as mp


def train_on_parameter(name, param,haha,hello):
    result = 0
    for num in param:
        result += math.sqrt(num * math.tanh(num) / math.log2(num) / math.log10(num))
    print(haha)
    print(hello)
    return {name: result}


if __name__ == '__main__':

    start_t = datetime.datetime.now()
    haha="haha"
    hello = "hello"

    num_cores = int(mp.cpu_count())
    print("本地计算机有: " + str(num_cores) + " 核心")
    pool = mp.Pool(num_cores)
    param_dict = {'task1': list(range(10, 3000)),
                  'task2': list(range(3000, 600000)),
                  'task3': list(range(600000, 900000)),
                  'task4': list(range(900000, 1200000)),
                  'task5': list(range(1200000, 1500000)),
                  'task6': list(range(1500000, 1800000)),
                  'task7': list(range(1800000, 2100000)),
                  'task8': list(range(2100000, 2400000))}
    results = [pool.apply_async(train_on_parameter, args=(name, param,haha,hello)) for name, param in param_dict.items()]
    results = [p.get() for p in results]

    end_t = datetime.datetime.now()
    elapsed_sec = (end_t - start_t).total_seconds()
    print("多进程计算 共消耗: " + "{:.2f}".format(elapsed_sec) + " 秒")