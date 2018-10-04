#!/usr/bin/python3.4

import sys
import math
import re
import subprocess

import numpy as np
from scipy.optimize import minimize

# physical params:
Ez = 0

#cache init:
cache_idx = 0
subprocess.call("[ -d cache ] && rm -r cache", shell=True)
subprocess.call("mkdir cache", shell=True)


def extract_classical_energy(log_file_path, target, target_nk):
  r = re.compile(r"\[DATA   \] \[NUMERICAL\] ground state classical energy     : (.*)")
  with open(log_file_path) as f:
   for line in f:
     line = line.strip()
     result = re.match(r, line)
     if result:
       print("[DEBUG] [f(x)] match line: ", line)
       energy = float(result.group(1).strip())
       return energy 
  return None


def extract_numerical_energy(log_file_path, target , target_nk):
  target_to_pattern = {
    'ground' : r"^\[DATA   \] \[NUMERICAL\] state nk, n_stars, n_stars_A, n_stars_B, energy: (.*),(.*),(.*),(.*),(.*)",
    'excited' : r"^\[DATA   \] \[NUMERICAL\] state nk, n_stars, n_stars_A, n_stars_B, energy, excitation_energy: (.*),(.*),(.*),(.*),(.*),(.*)"}
  r = re.compile(target_to_pattern[target])
  with open(log_file_path) as f:
   for line in f:
     line = line.strip()
     result = re.match(r, line)
     if result:
       #print("[DEBUG] [f(x)] pre-match line: ", line)
       energy = float(result.group(5).strip())
       nk = opt_int_str_to_int(result.group(1).strip())       
       if target == 'ground' or (target == 'excited' and nk == target_nk):
         print("[DEBUG] [f(x)] match line: ", line)
         print("[DEBUG] [f(x)] energy: ", energy)
         return energy 
  return None


def opt_int_str_to_int(s):
  if (s == "--"):
    return float("nan")
  else:
    return int(s)

def f(x, target : '("ground"|"excited")' = "ground", target_nk : '(int|None)' = None ):
  global cache_idx
  print("[DEBUG] [f(x)] [Enter]")
  print("[DEBUG] [f(x)] target    : ", target)
  print("[DEBUG] [f(x)] target_nk : ", target_nk)
  # init data:
  theta_0 = x[0]
  delta_theta = 0.0 if len(x) < 2 else x[1]
  print("[DEBUG] [f(x)] theta_0:     ", theta_0)
  print("[DEBUG] [f(x)] delta_theta: ", delta_theta)
  # run starsring:
  log_file_name = "log_{cache_idx:04d}_{theta_0:+6.9f}_{delta_theta:+6.9f}.log".format(cache_idx = cache_idx, theta_0 = theta_0, delta_theta = delta_theta)
  log_file_path = "cache/" + log_file_name
  target_to_switch = {
    'ground' : '--omit_one_star_space_calculations',
    'excited' : '--omit_zero_star_space_calculations'}
  command = './stars_ring -H jabcdzx --phi0 {theta_0} --delta_phi {delta_theta} -z{Ez} --n_extra_states_one_star_space 4'.format(theta_0 = theta_0, delta_theta = delta_theta, Ez = Ez) + ' ' + target_to_switch[target]
  command =  command + " > " + log_file_path + " 2>&1"
  print("[DEBUG] [f(x)] command: " + command)
  subprocess.check_call(command, shell=True)
  # create link, this may be useful when looking into the cache directory:
  link_file_name = "link_{cache_idx:04d}.log".format(cache_idx = cache_idx)
  link_file_path = "cache/" + link_file_name
  command2 = "ln -s \"{from_path}\" \"{to_path}\"".format(from_path = log_file_name, to_path = link_file_path)
  print("[DEBUG] [f(x)] command2: " + command2)
  subprocess.check_call(command2, shell=True)
  # now advance cach idx:
  cache_idx += 1
  # parse lig file:
  #energy = extract_classical_energy(log_file_path, target, target_nk)
  energy = extract_numerical_energy(log_file_path, target , target_nk)
  if (energy == None):
    print("[ERROR] [f(x)] fail to grep the output energy")
    raise RuntimeError("[f(x)]: Fail to grep the output energy");
  return energy

init_guesses = []
#init_guesses.append([0, 0])
#init_guesses.append([math.pi, 0])
#init_guesses.append([0, math.pi / 2.0])
#init_guesses.append([math.pi/2, math.pi/2])
x=0.5
y=0.00
#init_guesses.append([math.pi/2-x - y, 2*x])
init_guesses.append([math.pi/2-x - y, 2*x*0.90])

#init_guesses.append([math.pi/2, 0.1])
#init_guesses.append([1.0, 0.1])

#ctx = {"target" : "ground"}
ctx = {"target" : "excited", "target_nk" : 2 } 

results = [ minimize(lambda x : f(x, **ctx), init_guess, method='nelder-mead') for init_guess in init_guesses] #options={'xtol': 1e-8, 'disp': True}
the_best_result = None
for result in results :
  print("result:")
  print(result)
  print()
 
successful_results = [result for result in results if result.success]
the_best_result = min(successful_results, key = lambda result: result.fun)

print("the best result:")
print(the_best_result)
print()

print("refactor the result:")
if (the_best_result.x[0] < math.pi/2):
  print("opt x: {phi_0:+8.6f} {delta_phi:+8.6f}".format(phi_0 = the_best_result.x[0], delta_phi = the_best_result.x[-1]))
else:
  print("mirror!")
  print("opt x: {phi_0:+8.6f} {delta_phi:+8.6f}".format(phi_0 =math.pi-the_best_result.x[0], delta_phi = -the_best_result.x[-1]))
