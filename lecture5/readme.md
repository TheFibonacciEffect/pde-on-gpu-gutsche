# Lecture 5
## Exercise 2 task 2
![l5e2t2](figs/memcopy.png)

I have measured the peak memory bandwith of my [Intel® Core™ i5-10210U Processor](https://ark.intel.com/content/www/us/en/ark/products/195436/intel-core-i510210u-processor-6m-cache-up-to-4-20-ghz.html) by copying an $n \times n$ array for different $n$. The max memory bandwidth by the vendor is 45.8 GB/s. The peak memory throughput found using the benchmark tools is 15.33 GB/s for array programming and 5.89 GB/s using kernel programming. The array programming approach seems to be better suited for most array sizes.