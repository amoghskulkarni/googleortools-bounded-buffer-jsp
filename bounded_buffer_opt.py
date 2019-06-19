from __future__ import print_function
import json
# from math import inf

from datetime import datetime
import time as time_

def millis():
    return int(round(time_.time() * 1000))

# Import Python wrapper for or-tools constraint solver.
from ortools.constraint_solver import pywrapcp

def main(buffer_capacities):
  with open('jobfile.json', 'r') as f:
    jobfile = json.load(f)

  # Upper limit for cycle time 
  horizon = jobfile['requirements']['max_cycle_time']

  # Get the job types
  job_types = []
  for job_type in jobfile['job_types']:
    task_cycle_times = []
    for task in job_type['tasks']:
      task_cycle_times.append(int(task['cycle_time']['mean']))
    job_types.append(task_cycle_times)
  
  # Get machines count 
  machines_count = jobfile['system']['machines']

  # Get buffers count
  buffers_count = jobfile['system']['buffers']

  # Get jobs count
  job_counts = []
  total_jobs_count = 0
  for job_req in jobfile['requirements']['jobs']:
    job_counts.append(job_req['quantity'])
  total_jobs_count = sum(job_counts)
  
  # Get buffer capacities 
  buffer_capacities = []
  for buf_info in jobfile['requirements']['buffer_capacities']:
    buffer_capacities.append(buf_info['capacity'])

  # Create the solver.
  solver = pywrapcp.Solver('jobshop')

  all_machines = range(machines_count)
  all_jobs = range(total_jobs_count)
  
  # Data
  machines = []
  for job_n in range(total_jobs_count):
    machines.append(range(machines_count))

  processing_times = []
  # for job_n in range(total_jobs_count):
  #   processing_times.append(job_types[job_n % len(job_types)])
  for job_n in range(len(job_types)):
    # Assumes that len(job_types) == len(job_counts) and the requirements are given in the same order
    for _ in range(job_counts[job_n]):
      processing_times.append(job_types[job_n])

  # Create jobs
  all_tasks = {}
  for i in all_jobs:
    for j in range(0, len(machines[i])):
      all_tasks[(i, j)] = solver.FixedDurationIntervalVar(0,
                                                          horizon,
                                                          processing_times[i][j],
                                                          False,
                                                          'Job_%i_%i' % (i, j))

  # Creates sequence variables and add disjunctive constraints.
  all_sequences = []
  for i in all_machines:
    machines_jobs = []
    for j in all_jobs:
      for k in range(0, len(machines[j])):
        if machines[j][k] == i:
          machines_jobs.append(all_tasks[(j, k)])
    disj = solver.DisjunctiveConstraint(machines_jobs, 'machine %i' % i)
    all_sequences.append(disj.SequenceVar())
    solver.Add(disj)

  # Add conjunctive contraints.
  for i in all_jobs:
    for j in range(0, len(machines[i]) - 1):
      solver.Add(all_tasks[(i, j + 1)].StartsAfterEnd(all_tasks[(i, j)]))
  
  ###
  job_present_at_enqueue = {}
  job_present_at_dequeue = {}
  # Create the Integer variables for every task's completion and starting times
  for i in range(len(all_jobs)):
    for j in range(len(all_machines) - 1):
      for k in range(len(all_jobs)):
        for l in range(len(all_machines) - 1):
          if i != k and j == l:
            # print("B_%i_%i__E_%i_%i" % (i, j, k, l))
            enq = solver.IntVar(0, 1, "B_%i_%i__E_%i_%i" % (i, j, k, l))
            job_present_at_enqueue[(i, j, k, l)] = enq
            
            # all_tasks[(p, q)].EndExpr() ---- E_pq
            # all_tasks[(p, q)].StartExpr() -- D_pq

            # enq = 1 if (E_ij <= E_kl) and (E_kl < D_ij+1)
            # AND
            # enq = 0 otherwise

            enq_value = solver.ConditionalExpression(all_tasks[(i, j)].EndExpr() <= all_tasks[(k, l)].EndExpr() and all_tasks[(k, l)].EndExpr() < all_tasks[(i, j+1)].StartExpr(), solver.IntConst(1), 0)
            solver.Add(enq == enq_value)

  for i in range(len(all_jobs)):
    for j in range(len(all_machines) - 1):
      for k in range(len(all_jobs)):
        for l in range(1, len(all_machines)):
          if i != k and j == l:
            # print("B_%i_%i__D_%i_%i" % (i, j, k, l))
            deq = solver.IntVar(0, 1, "B_%i_%i__D_%i_%i" % (i, j, k, l))
            job_present_at_dequeue[(i, j, k, l)] = deq
            
            # all_tasks[(p, q)].EndExpr() ---- E_pq
            # all_tasks[(p, q)].StartExpr() -- D_pq

            # deq = 1 if (E_ij <= D_kl) and (D_kl < D_ij+1)
            # AND
            # deq = 0 otherwise

            deq_value = solver.ConditionalExpression(all_tasks[(i, j)].EndExpr() <= all_tasks[(k, l)].StartExpr() and all_tasks[(k, l)].StartExpr() < all_tasks[(i, j+1)].StartExpr(), solver.IntConst(1), 0)
            solver.Add(deq == deq_value)
  
  for l in range(len(all_machines) - 1):
    for j in range(len(all_machines) - 1):
      for k in range(len(all_jobs)):
        sum_constraint_over = []
        for i in range(len(all_jobs)):
          if i != k and j == l:
            sum_constraint_over.append(job_present_at_enqueue[(i, j, k, l)])
        solver.Add(solver.Sum(sum_constraint_over) < buffer_capacities[j])

  # for l in range(1, len(all_machines)):
  #   for j in range(1, len(all_machines) - 1):
  #     for k in range(len(all_jobs)):
  #       sum_constraint_over = []
  #       for i in range(len(all_jobs)):
  #         if i != k and j == l:
  #           sum_constraint_over.append(job_present_at_dequeue[(i, j, k, l)])
  #       solver.Add(solver.Sum(sum_constraint_over) > 0)
  ##

  # Set the objective.
  obj_var = solver.Max([all_tasks[(i, len(machines[i])-1)].EndExpr()
                        for i in all_jobs])
  objective_monitor = solver.Minimize(obj_var, 1)
  # Create search phases.
  sequence_phase = solver.Phase([all_sequences[i] for i in all_machines],
                                solver.SEQUENCE_DEFAULT)
  vars_phase = solver.Phase([obj_var],
                            solver.CHOOSE_FIRST_UNBOUND,
                            solver.ASSIGN_MIN_VALUE)
  main_phase = solver.Compose([sequence_phase, vars_phase])
  # Create the solution collector.
  collector = solver.LastSolutionCollector()

  # Add the interesting variables to the SolutionCollector.
  collector.Add(all_sequences)
  collector.AddObjective(obj_var)

  for i in all_machines:
    sequence = all_sequences[i]
    sequence_count = sequence.Size()
    for j in range(0, sequence_count):
      t = sequence.Interval(j)
      collector.Add(t.StartExpr().Var())
      collector.Add(t.EndExpr().Var())
  # Solve the problem.
  disp_col_width = 10
  before = millis()
  if solver.Solve(main_phase, [objective_monitor, collector]):
    after = millis()
    print("\nTime taken: ", str(after - before))
    print("\nOptimal Schedule Length:", collector.ObjectiveValue(0), "\n")
    sol_line = ""
    sol_line_tasks = ""
    print("Optimal Schedule", "\n")
    
    schedulefile = open("schedulefile.json", "w")
    schedule_data = {
      "machine_count": machines_count,
      "job_count": total_jobs_count,
      "optimal_schedule_length": collector.ObjectiveValue(0),
      "machine_schedules": []
    }

    for i in all_machines:
      seq = all_sequences[i]
      sol_line += "Machine " + str(i) + ": "
      sol_line_tasks += "Machine " + str(i) + ": "
      sequence = collector.ForwardSequence(0, seq)
      seq_size = len(sequence)

      for j in range(0, seq_size):
        t = seq.Interval(sequence[j])
        # Add spaces to output to align columns.
        sol_line_tasks +=  t.Name() + " " * (disp_col_width - len(t.Name()))

      for j in range(0, seq_size):
        t = seq.Interval(sequence[j])
        sol_tmp = "[" + str(collector.Value(0, t.StartExpr().Var())) + ","
        sol_tmp += str(collector.Value(0, t.EndExpr().Var())) + "] "
        # Add spaces to output to align columns.
        sol_line += sol_tmp + " " * (disp_col_width - len(sol_tmp))
      
      ### Populate schedule data
      task_types_sequence = []
      for j in sequence:
        j_copy = j
        # determine which kind of job this task is of 
        for k, job_count in enumerate(job_counts):
          if j_copy >= job_count:
            j_copy -= job_count
          else:
            task_types_sequence.append(jobfile['job_types'][k]['tasks'][i])
            break

      schedule_data["machine_schedules"].append({
        "machine": i,
        "tasks": sequence,
        "task_types": task_types_sequence
      })
      ###

      sol_line += "\n"
      sol_line_tasks += "\n"
    
    json.dump(schedule_data, schedulefile, indent=4)
    schedulefile.close()

    print(sol_line_tasks)
    print("Time Intervals for Tasks\n")
    print(sol_line)
    return collector.ObjectiveValue(0)
  else:
    print('No valid solution')
    return horizon

if __name__ == '__main__':
  sched = main([4, 4])

  # Plot buffer capacities vs optimal sched length
  # import plotly.plotly as py
  # import plotly.graph_objs as go

  # X_data = []
  # Y_data = []
  # Z_data = []
  # for buffer_1_capacity in range(3, 7):
  #   for buffer_2_capacity in range(3, 7):
  #     optimal_sched_length = main([buffer_1_capacity, buffer_2_capacity])
  #     X_data.append(buffer_1_capacity)
  #     Y_data.append(buffer_2_capacity)
  #     Z_data.append(optimal_sched_length)
  
  # print(X_data)
  # print(Y_data)
  # print(Z_data)

  # data = [
  #   go.Surface(
  #       x=X_data,
  #       y=Y_data,
  #       z=Z_data
  #   )
  # ]
  # layout = go.Layout(
  #     title='Optimal schedules vs. Buffer capacities',
  #     autosize=True,
  #     margin=dict(
  #         l=65,
  #         r=50,
  #         b=65,
  #         t=90
  #     )
  # )
  # fig = go.Figure(data=data, layout=layout)
  # py.iplot(fig, filename='Optimal schedules')

