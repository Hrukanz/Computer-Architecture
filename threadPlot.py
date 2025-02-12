import subprocess
import matplotlib.pyplot as plt
import re

n = 51 # Grid size
iterations = 10000
max_threads = 10  # Maximum number of threads 

thread_counts = range(1, max_threads + 1)

threads_list = []
times_list = []

compile_command = "gcc -o poisson poisson.c -lpthread"
compilation = subprocess.run(compile_command, shell=True)
if compilation.returncode != 0:
    print("Compilation failed.")

for threads in thread_counts:
    command = [
        "./poisson",
        f"-n {n}",
        f"-i {iterations}",
        f"-t {threads}"
    ]
    print(f"Running with {threads} thread(s)...")
    result = subprocess.run(' '.join(command), shell=True, capture_output=True, text=True)

    match = re.search(r"Execution time: ([\d\.]+) seconds", result.stdout)
    if match:
        execution_time = float(match.group(1))
        threads_list.append(threads)
        times_list.append(execution_time)
        print(f"Execution time: {execution_time} seconds")

plt.figure(figsize=(10, 6))
plt.plot(threads_list, times_list, marker='o')
plt.title('Number of Threads vs. Runtime')
plt.xlabel('Number of Threads')
plt.ylabel('Runtime (seconds)')
plt.xticks(thread_counts)
plt.grid(True)
plt.show()