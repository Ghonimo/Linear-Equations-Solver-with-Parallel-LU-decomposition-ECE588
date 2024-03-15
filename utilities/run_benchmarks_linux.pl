#!/usr/bin/perl

use strict;
use warnings;

# Step 1: Parse command-line arguments
my ($matrix_file_path, $max_threads) = @ARGV;
die "Usage: $0 <matrix_file_path> <max_threads>\n" unless defined $matrix_file_path && defined $max_threads;

# Validate max_threads
die "Maximum number of threads must be a positive integer\n" unless $max_threads =~ /^\d+$/ && $max_threads > 0;

# Step 2: Extract matrix size from the given path
my ($matrix_size) = $matrix_file_path =~ m|/(\d+x\d+)\.txt$|;
die "Could not extract matrix size from file path\n" unless defined $matrix_size;

# Construct the output file name dynamically
my $output_file = "benchmarks/linux_data/${matrix_size}_mac_${max_threads}.txt";

# Step 3: Open the output file for writing
open(my $fh, '>', $output_file) or die "Could not open file '$output_file' $!";

# Print the header to the file
printf $fh "%-20s %-15s %-15s\n", "Number of Threads", "Time (Seconds)", "Speedup";

# Variable to store time for single-thread execution
my $single_thread_time;

# Step 4: Run the program with different numbers of threads (1 to max_threads)
for (my $num_threads = 1; $num_threads <= $max_threads; $num_threads++) {

    my $output = `./bin/linux/parallel_solver $matrix_file_path $num_threads`;

    # Extract the time from the output
    my ($time_in_seconds) = $output =~ /Time taken: ([\d\.]+) seconds/;

    if (defined $time_in_seconds) {
        # Calculate speedup
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
