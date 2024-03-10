#!/usr/bin/perl

use strict;
use warnings;

# Output file path
#my $output_file = 'benchmarks/1000x1000_mac.txt';
#my $output_file = 'benchmarks/2000x2000_mac.txt';
#my $output_file = 'benchmarks/5000x5000_mac.txt';
#my $output_file = 'benchmarks/5000x5000_linux_16.txt';
my $output_file = 'benchmarks/10000x10000_linux_32.txt';

# Open the file for writing
open(my $fh, '>', $output_file) or die "Could not open file '$output_file' $!";

# Print the header to the file
printf $fh "%-20s %-15s %-15s\n", "Number of Threads", "Time (Seconds)", "Speedup";

# Variable to store time for single-thread execution
my $single_thread_time;

# Run the program with different numbers of threads (1 to 16)
for (my $num_threads = 1; $num_threads <= 32; $num_threads++) {

    my $output = `./bin/parallel_solver_linux matrices/py_generated/10000x10000.txt $num_threads`; 

    # Extract the time from the output.
    my ($time_in_seconds) = $output =~ /Time taken: ([\d\.]+) seconds/;

    if (defined $time_in_seconds) {
        # Calculate speedup. Avoid division by zero by checking if $num_threads is 1.
        my $speedup = $num_threads == 1 ? 1 : ($single_thread_time / $time_in_seconds);
        $single_thread_time = $time_in_seconds if $num_threads == 1;

        # Write the results to the output file
        printf $fh "%-20d %-15f %-15f\n", $num_threads, $time_in_seconds, $speedup;
    } else {
        warn "Could not extract execution time for $num_threads threads";
    }
}

# Close the file
close $fh;
print "Results saved in $output_file\n";
