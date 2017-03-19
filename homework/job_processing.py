#!/usr/bin/python3.5

def diff_cost(job):
  return job[0] - job[1]

def ratio_cost(job):
  return job[0] / job[1]

def job_processing(cost_fn):
  # input_file_name = 'job_processing_test_results_are_68615_67247.txt'
  input_file_name = 'job_processing_input.txt'

  jobs = []
  time = 0
  weighted_completion_time = 0

  with open(input_file_name) as f:
    for line in f:
      job = [int(x) for x in line.rstrip().split(' ')]

      if len(job) == 2:
        jobs.append((cost_fn(job), job[0], job[1]))

  while len(jobs):
    max_job = jobs[0]
    max_index = 0

    for index, job in enumerate(jobs):
      if job[0] > max_job[0] or (job[0] == max_job[0] and job[1] >= max_job[1]):
        max_job = job
        max_index = index

    jobs.pop(max_index)

    time += max_job[2]
    weighted_completion_time += max_job[1] * time

  print('Done, ' + str(weighted_completion_time))

if __name__ == '__main__':
  job_processing(diff_cost)
  job_processing(ratio_cost)
