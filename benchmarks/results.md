# Results and Performance Analysis

## Introduction
This document presents the performance results obtained from the parallelization of the LU Decomposition algorithm using POSIX threads. Our primary focus was to evaluate the computational efficiency and numerical stability enhancements achieved through parallel execution on multicore systems.

## Experimental Setup
Experiments were conducted on two hardware platforms to assess the scalability and efficiency of the parallel algorithm:
- **Linux System:** Equipped with an Intel Xeon E312xx (Sandy Bridge) processor featuring 32 CPUs.
- **MacBook Pro:** Housing an M3 Pro chip with 11 cores.

Matrices ranging from 3x3 to 10000x10000, encompassing both sparse and dense configurations, were used to test the algorithm under varied conditions.

## Results

### Performance Evaluation on Different Systems
The parallel LU Decomposition algorithm demonstrated significant speed improvements across all tested systems. Below are the highlights:

#### Linux Platform (32 CPUs)
- **5000x5000 Matrix:** Achieved a maximum speedup of 7.934753 with 27 threads.
- **10000x10000 Matrix:** Speedup peaked at 10.781217 with 31 threads, showcasing efficient utilization of computational resources.

#### MacBook Pro (11 Cores)
- Despite hardware limitations, notable speedups were observed, with a **5000x5000 matrix** reaching a speedup of 4.940405 with 10 threads.

### Speedup and Efficiency Analysis
Our analysis indicates that the speedup and efficiency of the parallel LU Decomposition are directly influenced by the number of threads utilized. Optimal performance was observed when the number of threads closely matched the physical core count of the execution environment.

### Impact of Partial Pivoting on Performance
The incorporation of partial pivoting slightly extended computation times but was critical for enhancing numerical stability, particularly for ill-conditioned matrices. The trade-off between a negligible increase in computation time and significant gains in accuracy justifies the use of partial pivoting.

## Discussion
The experimental results underscore the potential of parallel computing to significantly enhance the efficiency of numerical algorithms like LU Decomposition. The successful application of the BSP model and POSIX threads in our project illustrates a viable pathway to tackling large-scale computational problems more effectively.

## Future Work
Further optimizations and exploration of sophisticated parallelization strategies could yield even greater performance improvements. Additionally, applying the algorithm to real-world problems and extending its applicability to different matrix types remain promising avenues for future research.

## Conclusion
The parallel LU Decomposition algorithm represents a significant step forward in our ability to solve systems of linear equations efficiently. This project not only demonstrated substantial speedups but also provided insights into the critical balance between computational speed and numerical accuracy.

