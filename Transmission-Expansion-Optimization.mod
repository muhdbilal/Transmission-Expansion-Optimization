reset;

set BUS;    # set of buses
set BRANCH; # set of branches
set NBRANCH; # set of new branches
set GEND;   # Gen Data


#@@@@@@@@@@@@@@@
#### PARAMETERS:
# Bus Data
param bus_num		{BUS}; # Bus Number
param bus_Pd		{BUS}; # Real Power Demand
param bus_RPd		{BUS}; # Reactive Power Demand

# GEN Data
param genD_bus		{GEND}; # GEN location
param genD_Pmax		{GEND}; # Max gen production
param genC_Cost		{GEND}; # Linear Cost Term

# Branch Data
param branch_fbus	{BRANCH}; # from bus for line
param branch_tbus	{BRANCH}; # to bus for line
param branch_length	{BRANCH}; # length of line
param branch_b		{BRANCH}; # line susceptance
param branch_rateA	{BRANCH}; # long term thermal rating

# New Branch Data
param nbranch_fbus		{NBRANCH}; # from bus for line
param nbranch_tbus		{NBRANCH}; # to bus for line
param nbranch_length	{NBRANCH}; # length of line
param nbranch_b			{NBRANCH}; # line susceptance
param nbranch_rateA		{NBRANCH}; # long term thermal rating
param nbranch_cost		{NBRANCH}; # installation cost of the branch

param MBase; let MBase:=100; # the MVA Base
param M; let M:=1000; # the big-M parameter

#@@@@@@@@@@@@@@@


#@@@@@@@@@@@@@@@
#### VARIABLES:
var gen_supply {g in GEND};     # Variable for GEN Supply
var bus_theta {b in BUS};		#bus phase
var branch_load {k in BRANCH};  # Load on each branch
var nbranch_load {k in NBRANCH};  # Load on each new branch
var W1{k in NBRANCH} binary; # wheter a new branch is to be installed
#@@@@@@@@@@@@@@@


#@@@@@@@@@@@@@@@
#### OBJECTIVE:
# Objective is to Minimize Cost
minimize COST: sum{g in GEND}(gen_supply[g]*genC_Cost[g]) + sum{k in NBRANCH}(W1[k]*nbranch_cost[k]*nbranch_length[k]); 
#@@@@@@@@@@@@@@@


#@@@@@@@@@@@@@@@
###CONSTRAINTS:

subject to PGen{g in GEND}:
	0 <= gen_supply[g] <= genD_Pmax[g];

subject to Branch_limit{k in BRANCH}:
	-branch_rateA[k] <= branch_load[k] <= branch_rateA[k];

subject to b_thetha{k in BRANCH}:
	branch_load[k] = MBase*branch_b[k]*(bus_theta[branch_tbus[k]] - bus_theta[branch_fbus[k]]);


subject to new_Branch_limit_A{k in NBRANCH}:
	nbranch_load[k] >= -1*(nbranch_rateA[k])*(W1[k]);

subject to new_Branch_limit_B{k in NBRANCH}:
	nbranch_load[k] <= nbranch_rateA[k]*(W1[k]);

subject to new_b_thetha_1{k in NBRANCH}:
	-(1 - (W1[k]))*M <= nbranch_load[k] - MBase*nbranch_b[k]*(bus_theta[nbranch_tbus[k]] - bus_theta[nbranch_fbus[k]]);

subject to new_b_thetha_2{k in NBRANCH}:
	nbranch_load[k] - MBase*nbranch_b[k]*(bus_theta[nbranch_tbus[k]] - bus_theta[nbranch_fbus[k]]) <= (1 - (W1[k]))*M;


subject to supplydemand{n in BUS}:
	sum{x in NBRANCH: nbranch_fbus[x]==n}(nbranch_load[x]) + sum{s in BRANCH: branch_fbus[s]==n}(branch_load[s]) - sum{r in BRANCH: branch_tbus[r]==n}(branch_load[r]) - sum{y in NBRANCH: nbranch_tbus[y]==n}(nbranch_load[y]) = sum{g in GEND: genD_bus[g]==n}(gen_supply[g]) - bus_Pd[n];


#@@@@@@@@@@@@@@@

#### Load data:
data;

param: BUS: bus_num bus_Pd bus_RPd:= include bus_data.dat;
param: GEND: genD_bus genD_Pmax genC_Cost:= include gen_data.dat;
param: BRANCH: branch_fbus branch_tbus branch_length branch_b branch_rateA:= include branch_data.dat;
param: NBRANCH: nbranch_fbus nbranch_tbus nbranch_length nbranch_b nbranch_rateA nbranch_cost:= include nbranch_data.dat;

#@@@@@@@@@@@@@@@

for {n in BUS}{
	let bus_Pd[n] := 3*bus_Pd[n];
	let bus_RPd[n] := 3*bus_RPd[n];
};

for {k in BRANCH}{
	let branch_b[k] := 1/branch_b[k]
};

for {k in NBRANCH}{
	let nbranch_b[k] := 1/nbranch_b[k]
};

#@@@@@@@@@@@@@@@



#@@@@@@@@@@@@@@@

option solver gurobi;
option gurobi_options('mipgap=0.001 timelim=1000');
solve;
display _total_solve_elapsed_time; 
option show_stats 1;


for{n in NBRANCH}{
printf "%8.2f\n", W1[n] > W.Result1;
printf "%8.2f\n", nbranch_load[n] > new_branch_loads.Result2;
};

for{n in BRANCH}{
printf "%8.2f\n", branch_load[n] > original_branch_loads.Result3;
};

for{g in GEND}{
printf "%8.2f\n", gen_supply[g] > gen_output.Result4;
};



