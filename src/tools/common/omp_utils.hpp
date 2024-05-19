//
// Created by anton on 7/27/20.
//
#pragma once
#include "logging.hpp"
#include <functional>
#include <parallel/algorithm>
#include <omp.h>
#include <utility>
#include <numeric>
#include <wait.h>
#include "unistd.h"

//用于在多线程环境下进行并行计数
template<typename T>
class UniversalParallelCounter {
    std::vector<T> cnt;
public:
    explicit UniversalParallelCounter(size_t thread_num) : cnt(thread_num){
    }

    void operator++() {
        cnt[omp_get_thread_num()] += 1;
    }

    void operator+=(const T val) {
        cnt[omp_get_thread_num()] += val;
    }

    size_t get() const {
        return std::accumulate(cnt.begin(), cnt.end(), size_t(0));
    }
};

typedef UniversalParallelCounter<size_t> ParallelCounter;

template<class T>
class ParallelRecordCollector {
    std::vector<std::vector<T>> recs;
public:
    friend class Iterator;
    class Iterator : public std::iterator<std::forward_iterator_tag, T, size_t,  T*, T&>{
    private:
        ParallelRecordCollector<T> &data;
        size_t row;
        size_t col;
    public:
        explicit Iterator(ParallelRecordCollector<T> &_data, size_t _row = 0, size_t _col = 0) : data(_data), row(_row), col(_col) {
            while(row < data.recs.size() && col == data.recs[row].size()) {
                row += 1;
                col = 0;
            }
        }

        void operator++() {
            col += 1;
            while(row < data.recs.size() && col == data.recs[row].size()) {
                row += 1;
                col = 0;
            }
        }

        T &operator *() {
            return data.recs[row][col];
        }

        bool operator==(const Iterator &other) {
            return row == other.row && col == other.col;
        }
        bool operator!=(const Iterator &other) {
            return row != other.row || col != other.col;
        }

    };
    explicit ParallelRecordCollector(size_t thread_num) : recs(thread_num){
    }

    void add(const T &rec) {
        recs[omp_get_thread_num()].emplace_back(rec);
    }

    template<class I>
    void addAll(I begin, I end) {
        recs[omp_get_thread_num()].insert(recs[omp_get_thread_num()].end(), begin, end);
    }

    template< class... Args >
    void emplace_back( Args&&... args ) {
        recs[omp_get_thread_num()].emplace_back(args...);
    }

    Iterator begin() {
        return Iterator(*this, 0, 0);
    }

    Iterator end() {
        return Iterator(*this, recs.size(), 0);
    }

    size_t size() const {
        size_t res = 0;
        for (const std::vector<T> & row : recs) {
            res += row.size();
        }
        return res;
    }

    bool empty() const {
        return size() == 0;
    }

    std::vector<T> collect() {
        std::vector<T> res;
        for(std::vector<T> &row : recs) {
            for(T &val : row)
                res.emplace_back(std::move(val));
            row.clear();
        }
        return std::move(res);
    }

    void clear() {
        for(std::vector<T> &row : recs) {
            row.clear();
        }
    }

    std::vector<T> collectUnique() {
        std::vector<T> res = collect();
        __gnu_parallel::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        return std::move(res);
    }
};

template<class T>
std::ostream& operator<<(std::ostream& out, const ParallelRecordCollector<T>& tree) {
    if(tree.size() == 0) {
        return out << "[]" << std::endl;
    }
    out << "[";
    for(const T & item: tree) {
        out << item << ", ";
    }
    return out << "]";
}


template<class V>
class ParallelProcessor {
public:
    std::function<void(size_t, V &)> task = [] (size_t, V &) {};
    std::function<void ()> doBefore = [] () {};
    std::function<void ()> doAfter = [] () {};
    std::function<void ()> doInParallel = [] () {};
    std::function<void (V&)> doInOneThread = [] (V &) {};
    std::function<void ()> doInTheEnd = [] () {};
    logging::Logger &logger;
    size_t threads;

    /**
     * @brief 并行处理器
     *
     * 使用指定的任务函数、日志记录器和线程数创建一个并行处理器对象。
     *
     * @param _task 任务函数，用于处理每个元素
     * @param _logger 日志记录器对象，用于记录日志信息
     * @param _threads 线程数，指定并行处理时使用的线程数量
     */
    ParallelProcessor(std::function<void(size_t, V &)> _task, logging::Logger & _logger, size_t _threads) :
                    task(_task), logger(_logger), threads(_threads) {
    }

//This method expects iterator to be a generator, i.e. it returns temporary objects. Thus we have to store them in a buffer and
//keep track of total size of stored objects.
    template<class I>
    /**
     * @brief 处理记录
     *
     * 将指定范围内的记录进行并行处理。
     * 作者采用的并行模型为openmp 并行编程框架,把底层抽象了,那么上层只需要调用函数
     * @param begin 记录范围的起始迭代器
     * @param end 记录范围的结束迭代器
     * @param bucket_length 每个桶的长度，默认为 1MB
     */
    void processRecords(I begin, I end, size_t bucket_length = 1024 * 1024) {
        //设置openmp的线程个数
        omp_set_num_threads(threads);
        //default(none):强制程序员明确指定所有变量的数据共享属性,
        //default(shared): 默认情况下所有变量都是共享的，这意味着它们对所有线程都是可见的，并且线程可以读取和修改它们。但是，这种使用需要格外小心，因为共享变量的修改可能导致数据竞争和其他并发问题。
#pragma omp parallel default(none)
        {//下面single相当于锁住了,只允许1个线程执行
#pragma omp single
            {
                logger.trace() << "Starting parallel calculation using " << omp_get_num_threads() << " threads" << std::endl;
            }
        };
//        size_t bucket_length = 1024 * 1024;
        size_t buffer_size = 1024 * 1024;
        size_t max_length = 1024 * 1024 * 1024;
        ParallelProcessor<V> &self = *this;
        size_t total = 0;
        size_t total_len = 0;
        while(begin != end) {
            size_t clen = 0;
            std::vector<V> items;
            items.reserve(buffer_size);//预留buffer_size的容器内存
            doBefore();
#pragma omp parallel default(none) shared(begin, end, items, buffer_size, clen, max_length, bucket_length, self, std::cout, total)
            {//除了begin, end, items, buffer_size, clen,
            // max_length, bucket_length, self, std::cout, total全部共享
#pragma omp single
                {
#pragma omp task default(none) shared(self, std::cout)
                    {
                        self.doInParallel();//1个线程负责执行任务
                    }
                    //剩下的1个线程负责执行下面代码
                    size_t cur_pos = 0;
                    size_t left = cur_pos;
                    size_t right = cur_pos;
                    size_t cur_length = 0;
                    while (begin != end && items.size() < buffer_size && clen < max_length) {
                        items.emplace_back(*begin);
                        ++begin;
                        clen += items.back().size();
                        cur_length += items.back().size();
                        right += 1;
                        if(cur_length >= bucket_length || begin == end || items.size() >= buffer_size || clen >= max_length) {
#pragma omp task default(none) shared(items, self, std::cout) firstprivate(total, left, right)
                            {
                                for(size_t i = left; i < right; i++)
                                    self.task(total + i, items[i]);
                            }
                            left = right;
                            cur_length = 0;
                        }
                    }
                }
            }
            doAfter();
            logger.trace() << items.size() << " items of total length "<< clen << " processed " << std::endl;
            total += items.size();
            items.clear();
            total_len += clen;
        }
        doInTheEnd();
        logger.trace() << "Finished parallel processing. Processed " << total <<
               " items with total length " << total_len << std::endl;
    }


    //This method expects that iterators return references to objects instead of temporary objects.
    template<class I>
    void processObjects(I begin, I end, size_t bucket_size = 1024) {
        logger.trace() << "Starting parallel calculation" << std::endl;
        omp_set_num_threads(threads);
        ParallelProcessor<V> &self = *this;
        size_t buffer_size = 1024 * 1024;
        size_t total = 0;
        while(begin != end) {
            std::vector<V*> items;
            items.reserve(buffer_size);
            doBefore();
#pragma omp parallel default(none) shared(begin, end, items, buffer_size, bucket_size, self, total)
            {
#pragma omp single
                {
#pragma omp task default(none) shared(self)
                    {
                        self.doInParallel();
                    }
                    while (begin != end && items.size() < buffer_size) {
                        size_t left = items.size();
                        size_t right = items.size();
                        while (begin != end && items.size() < buffer_size && right - left < bucket_size) {
                            items.push_back(&(*begin));
                            self.doInOneThread(*items.back());
                            ++begin;
                            right += 1;
                        }
#pragma omp task default(none) shared(items, self) firstprivate(total, left, right)
                        {
                            for(size_t i = left; i < right; i++)
                                self.task(total + i, *items[i]);
                        }
                    }
                }
            }
            doAfter();
            logger.trace() << "Processed " << items.size() << " items" << std::endl;
            total += items.size();
            items.clear();
        }
        doInTheEnd();
        logger.trace() << "Finished parallel processing. Processed " << total << " items " << std::endl;
    }

};

//This method expects that iterators return references to objects instead of temporary objects.
template<class I>
void processObjects(I begin, I end, logging::Logger &logger, size_t threads, std::function<void(size_t, typename I::value_type &)> task,
                    size_t bucket_size = 1024) {
    typedef typename I::value_type V;
    ParallelProcessor<V>(task, logger, threads).processObjects(begin, end, bucket_size);
}

//This method expects iterator to be a generator, i.e. it returns temporary objects. Thus we have to store them in a buffer and
//keep track of total size of stored objects.
template<class I>
/**
 * @brief 处理记录
 *
 * 使用给定的参数在并行处理器中处理记录。
 *
 * @param begin 记录范围的起始迭代器
 * @param end 记录范围的结束迭代器
 * @param logger 日志记录器
 * @param threads 线程数
 * @param task 处理任务的函数
 * @param bucket_length 存储桶的长度，默认为 1MB（1024 * 1024）
 */
void processRecords(I begin, I end, logging::Logger &logger, size_t threads, std::function<void(size_t, typename I::value_type &)> task,
                    size_t bucket_length = 1024 * 1024) {
    typedef typename I::value_type V;
    ParallelProcessor<V>(task, logger, threads).processRecords(begin, end, bucket_length);
}

/**
 * @brief 在子进程中运行函数
 *
 * 在子进程中运行给定的函数，并在子进程结束后等待其退出。
 *
 * @param f 要在子进程中运行的函数
 */
inline void runInFork(const std::function<void()>& f) {
    pid_t p = fork();
    if (p < 0) {
        std::cout << "Fork failed" << std::endl;
        exit(1);
    }
    if(p == 0) {
        f();
        exit(0);
    } else {
        int status = 0;
        waitpid(p, &status, 0);
        if (WEXITSTATUS(status) || WIFSIGNALED(status)) {
            std::cout << "Child process crashed" << std::endl;
            exit(1);
        }
    }
}
