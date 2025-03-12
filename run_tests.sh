#!/bin/bash

# Define the thread counts to test
thread_counts=(1 2 4 8)

# Loop through each thread count
for threads in "${thread_counts[@]}"
do
    echo "Running tests with OMP_NUM_THREADS=$threads"
    export OMP_NUM_THREADS=$threads

    # Loop through each test file
    for i in {1..5}
    do
        exec_summary_file="test/exec_time/test0$i_summary.err"
        # Define the input, output, expected, and error file names
        input_file="test/in/test0$i.in"
        output_file="test/out/test0$i_${threads}threads.out"
        expected_file="test/expected/test0$i.out"
        exec_time_file="test/exec_time/test0$i_${threads}threads.err"

        
        touch "$exec_time_file"


        # Read input values from the input file
        # Assuming the input file has the values in the order: seed side ncside n_part time_steps
        read -r seed side ncside n_part time_steps < "$input_file"

        # Run the program with the input values as command-line arguments
        # Redirect stdout to the output file and stderr to the error file
        ./build/bin/parsim "$seed" "$side" "$ncside" "$n_part" "$time_steps" > "$output_file" 2> "$exec_time_file"

        # Compare the output to the expected output
        if diff -q "$output_file" "$expected_file" > /dev/null; then
            echo "Test $i with $threads threads: PASSED"
        else
            echo "Test $i with $threads threads: FAILED"
            echo "Differences:"
            diff "$output_file" "$expected_file"
        fi

        # Print the execution time (stderr content)
        echo -n "Test $i with $threads threads -> " >> "$exec_summary_file"
        cat "$exec_time_file" >> "$exec_summary_file"
    done
done