/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef THREADPOOL_HPP_
#define THREADPOOL_HPP_

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

/**
 * A simple, reusable fixed-size thread pool built on the C++11/17 standard library.
 *
 * The pool launches a fixed number of persistent worker threads at construction. Tasks are
 * submitted with Enqueue(), which returns a std::future so that the caller can wait for the
 * task and retrieve any exception it threw (exceptions raised inside a task are captured in
 * the returned future rather than terminating the program). The worker threads are shut down
 * and joined by the destructor.
 *
 * Typical use is a "parallel for": enqueue one task per worker, each of which pulls work items
 * from a shared std::atomic index counter, then wait on the returned futures. Results should be
 * written into per-item slots (distinct indices) so no locking is required for the results.
 *
 * The pool does not itself provide any thread-safety for the work the tasks do: tasks must not
 * touch shared mutable state or non-thread-safe singletons without their own synchronisation.
 */
class ThreadPool
{
public:
    /**
     * Constructor. Launches @p numThreads persistent worker threads.
     *
     * @param numThreads  the number of worker threads (must be >= 1).
     */
    explicit ThreadPool(unsigned numThreads)
            : mStop(false)
    {
        if (numThreads < 1u)
        {
            throw std::invalid_argument("ThreadPool requires at least one thread.");
        }
        mWorkers.reserve(numThreads);
        for (unsigned i = 0; i < numThreads; ++i)
        {
            mWorkers.emplace_back([this] { this->WorkerLoop(); });
        }
    }

    /**
     * Destructor. Signals the workers to stop once the queue is drained, and joins them.
     */
    ~ThreadPool()
    {
        {
            std::unique_lock<std::mutex> lock(mQueueMutex);
            mStop = true;
        }
        mCondition.notify_all();
        for (std::thread& r_worker : mWorkers)
        {
            if (r_worker.joinable())
            {
                r_worker.join();
            }
        }
    }

    // Non-copyable and non-movable (owns threads).
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /**
     * Submit a task to be run by one of the worker threads.
     *
     * @param task  a nullary callable to run. Any exception it throws is stored in the returned
     *              future and re-thrown when future.get() is called.
     * @return a std::future<void> that becomes ready when the task has finished.
     */
    std::future<void> Enqueue(std::function<void()> task)
    {
        auto p_packaged = std::make_shared<std::packaged_task<void()> >(std::move(task));
        std::future<void> result = p_packaged->get_future();
        {
            std::unique_lock<std::mutex> lock(mQueueMutex);
            if (mStop)
            {
                throw std::runtime_error("Enqueue called on a stopped ThreadPool.");
            }
            mTasks.emplace([p_packaged] { (*p_packaged)(); });
        }
        mCondition.notify_one();
        return result;
    }

    /**
     * @return the number of worker threads in the pool.
     */
    unsigned GetNumThreads() const
    {
        return (unsigned)mWorkers.size();
    }

private:
    /**
     * The loop each worker thread runs: wait for a task, run it, repeat until stopped and drained.
     */
    void WorkerLoop()
    {
        while (true)
        {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(mQueueMutex);
                mCondition.wait(lock, [this] { return mStop || !mTasks.empty(); });
                if (mStop && mTasks.empty())
                {
                    return;
                }
                task = std::move(mTasks.front());
                mTasks.pop();
            }
            task();
        }
    }

    /** The worker threads. */
    std::vector<std::thread> mWorkers;

    /** The queue of pending tasks. */
    std::queue<std::function<void()> > mTasks;

    /** Protects #mTasks and #mStop. */
    std::mutex mQueueMutex;

    /** Signals workers that a task is available or that the pool is stopping. */
    std::condition_variable mCondition;

    /** Set true by the destructor to tell workers to finish. */
    bool mStop;
};

#endif // THREADPOOL_HPP_
