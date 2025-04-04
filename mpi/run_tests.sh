#!/bin/bash

lines="83b2a6e80b806169ba35b56bf557eac9d12ebcd6" # Commit previous to 'implement blocks'
blocks="29b80082ee62a06fa670129aa4a6c565c0122326" # Commit previous to 'implement load balancing'
balancing="80deed8fd845680bb4e32b025ea099ae74c59a14" # Commit 'implement load balancing'
balancingv2="3cf326e99d304f7aba0cad0ff09f01fa54e52ebd" # Commit 'implement load balancing v2'
# openmp

# Commits to test
hashes=("$lines" "$blocks" "$balancing" "$balancingv2")
commits=("lines" "blocks" "balancing" "balancingv2")

# Stashing any changes before changing commits
git stash

# Loop through each thread count
for c in "${!hashes[@]}"
do
    echo "Running tests for commit ${commits[$c]}"
    git checkout "${hashes[$c]}" include src
    make clean && make

    exec_summary_file="test/exec_time/test${commits[$c]}.err"
    # Loop through each test file
    for i in $(seq -f "%02g" 1 12)
    do

        input_file="test/in/test$i.in"
        output_file="test/out/test${i}_${commits[$c]}.out"
        expected_file="test/expected/test$i.out"
        exec_time_file="test/exec_time/test${i}_${commits[$c]}.err"

        # Read input values from the input file
        # Assuming the input file has the values in the order: seed side ncside n_part time_steps
        read -r seed side ncside n_part time_steps < "$input_file"

        # Run the program with the input values as command-line arguments
        # Redirect stdout to the output file and stderr to the error file
        mpirun -n 8 parsim "$seed" "$side" "$ncside" "$n_part" "$time_steps" > "$output_file" 2> "$exec_time_file"

        # Compare the output to the expected output
        if diff -q "$output_file" "$expected_file" > /dev/null; then
            echo "Test $i for commit ${commits[$c]}: PASSED"
        else
            echo "Test $i for commit ${commits[$c]}: FAILED"
            echo "Differences:"
            diff "$output_file" "$expected_file"
        fi

        # Print the execution time (stderr content) in the desired format
        echo -n "TEST $i" >> "$exec_summary_file"
        echo -n "s" >> "$exec_summary_file"
        echo -n " -> " >> "$exec_summary_file"
        cat "$exec_time_file" >> "$exec_summary_file"
        rm "$exec_time_file"
    done
    git restore --staged .
    git restore .
done

# Restoring changes
git stash pop
