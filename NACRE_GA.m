function nga = NACRE_GA( initial_aircraft_sequence, ...
    deterministic_permutations, ground_database, runway_name )
%NACRE_GA Performs GA based search for an optimal sequencing solution
%   The GA works by nearest-neighbor reordering.  No aircraft can be out of
%   position by more than one place.  This creates a number of Abelian 
%   one-flip configurations that grows as N(L) = Fibonacci(L), or about
%   10^8 for L = 40 aircraft and 10^10 for L = 50 aircraft. This is still 
%   too large of a search space to search exhaustively and quickly.
%
%   The GA takes as input an initial schedule that does not incorporate
%   wake vortex constraints.  It uses that schedule to generate sequences
%   that do satisfy wake-vortex constraints and evolves their population.

%  Bruce Sawhill, PhD
%  NextGen AeroSciences, LLC
%  2014

%  Convert PGM input to canonical form (bks data format, 
%  schedled runway use ordering)

iac = PGM_converter(initial_aircraft_sequence, ground_database, runway_name);
igs = PGM_Ground_Converter(initial_aircraft_sequence,...
    ground_database, runway_name);
num_ac = size(iac,1);

%  Establish baseline of aircraft in initial scheduled order of runway
%  usage, modified for wake vortex to be physically flyable. Order by
%  plane number for comparison with new solutions.
 
init_runway_order = zeros(num_ac,2);

%  Keep track of the aircraft indices and their associated runway times for
%  later delay calculation

unsequenced_runway_times = zeros(num_ac,2);
unsequenced_runway_times = sortrows(iac(:,1:2),1);

%  Run initial sequence through sequencer to get vortex separated times
%  init_runway_order orders by index for comparison with initial conditions

init_runway_order(:,1) = iac(:,1);
report = sequencer(iac);
init_runway_order(:,2) = report.start_times;
init_runway_order = sortrows(init_runway_order, 1);

% Introduce small problem size exhaustive search here as an alternative to GA
if num_ac < 11

    exs = exhaustive_searchf(iac, deterministic_permutations);
    nga.delays = exs.delays;
    nga.fitness = exs.fitness;
    nga.runway_times = exs.runway_times;
    nga.gate_release_times = exs.gate_release_times;
    nga.spot_begin_end_times = exs.spot_begin_end_times;
    
    for j = 1: num_ac
        if iac(j,6) ~= -1
        end
    end
    
else
    
%  Global constants and arrays
rng(1);  %  Control random number generator
POPULATION = 2 * num_ac;
SELECTION_STRENGTH = 0.8;% Lower is stronger
NUM_GEN = 40; %  This may or may not be used depending on stopping criteria
%  Generate objects that will be combined genetically
%  These are single flip mutations.
%  Generate all nearest neighbors to the initial sequence
%  The genome is indexed by the the position of the first member of a flip
%  This is sort of an inside-out GA, that starts with mutations

of = oneflips(iac);

initial_scores = zeros(num_ac - 1, 3); % only N-1 single flips are possible
for i = 1: num_ac - 1
    
    sequence_stats = sequencer(of(:,:,i));
    initial_scores(i,:) = [i, sum(sequence_stats.delays), ...
        sequence_stats.start_times(end)-sequence_stats.start_times(1)];
    
end
    
ranked_scores = sortrows(initial_scores, 3);

% Generate new population of proper size to initialize the evolutionary loop
previous_generation = of;
cutoff = floor(SELECTION_STRENGTH * length(previous_generation));
ranked_scores(1: cutoff, :);

%  Create some random index pairs
indices = randi([1,num_ac-1], 3 * POPULATION,2);

sorted_indices = indices;

% Normal order the tuples
for j = 1: 2 * POPULATION
  sorted_indices(j,:) = sort(indices(j,:));  
end

% Find legitimate pairs (unique, no nearest-neighbor flips)
non_redundant = unique(sorted_indices,'rows');
lnr = length(non_redundant);

%  Filter out the illegal combinations of flips
non_redundant_filtered = [];

for g = 1: lnr  
    if abs(non_redundant(g,2)-non_redundant(g,1)) > 1
        
        non_redundant_filtered =...
            [non_redundant_filtered; non_redundant(g,:)];
        
    end
    
end

%  Generate a random subset of non_redundant_filtered to select next gen
lnrf = length(non_redundant_filtered);
new_indices = non_redundant_filtered(randperm(lnrf),:);
new_pop_indices = new_indices(1:POPULATION - cutoff,:);

genome = cell(POPULATION, 1);

%  Fill in the first slots with the winners of the 1-flip competition
for k = 1: cutoff
    genome{k} = ranked_scores(k,1);
end

%  Add some scores with double flips to make up a full population for
%  evolving
for k = cutoff + 1: POPULATION
    genome{k} = new_pop_indices(k - cutoff,:);
end

%  ================
%  ================


%  The start of the big evolutionary loop

for gen = 1: NUM_GEN
%  Make and evaluate new generation
new_generation = zeros(num_ac, 6, POPULATION);

for i = 1: POPULATION
    
    % Initialize start sequence to baseline ordering by runway usage time
    start_sequence = iac;
    
    for j = 1: length(genome{i})
        
        new_generation(:,:,i) = flipper(start_sequence, genome{i}(j));
        start_sequence = new_generation(:,:, i);
        
    end
    
end

%  Sequence and rank
evol_scores = zeros(POPULATION, 3);
ordering = zeros(num_ac, 2, POPULATION);
delays = zeros(num_ac,1, POPULATION);

for i = 1: POPULATION
    
    %  Rectify and score each member of the new generation
  
    seq = sequencer(new_generation(:,:,i));
    
      ordering(:,1,i) = new_generation(:,1,i);
      ordering(:,2,i) = seq.start_times;
    
    evol_scores(i,:) = [i, sum(seq.delays.^2), ...
        seq.start_times(end)-seq.start_times(1)];
    
end

%  Compute delays

delays = zeros(num_ac,1, POPULATION);
sorted = zeros(num_ac,2, POPULATION);
gate_release_times = -1 * ones(num_ac, 1, POPULATION);
spot_begin_end_times = -1 * ones(num_ac, 2, POPULATION);

for i = 1: POPULATION
    
  sorted(:,:,i) = sortrows(ordering(:,:,i),1); % put in original sequence
  delays(:,1,i) = sorted(:,1,i);
  delays(:,2,i) = sorted(:,2,i) - unsequenced_runway_times(:,2);
  
end
% sorti = sortrows(ordering(:,:,3),1)

%  Rank by sum of squares of delays (inverse of throughput)
rankings = sortrows(evol_scores, 2);
previous_generation = new_generation;
cutoff = floor(SELECTION_STRENGTH * length(previous_generation));
% best_makespan = evol_scores(1,3)
genome_source_list = cell(cutoff, 2);

for i = 1: cutoff
    genome_source_list{i} = genome{rankings(i,1)};
end

%  Generational reportage
% nga.genes = genome;
nga.fitness = evol_scores(:,2);
nga.delays = delays;
nga.runway_times = sorted;
nga.gate_release_times = gate_release_times;
nga.spot_begin_end_times = spot_begin_end_times;

%  Generate random pairs of genome indices for parents of next generation
indices = randi([1,cutoff], 2 * POPULATION,2);
sorted_indices = indices;
%  Sort them
for j = 1: 2 * POPULATION
  sorted_indices(j,:) = sort(indices(j,:));  
end

%  Find legitimate pairs (unique, no nearest-neighbor flips)
%  First, unique pairs of indices for parents
non_redundant = unique(sorted_indices,'rows');
lnr = length(non_redundant);

%  Filter out the illegal combos of flips
%  This requires mapping indices to genes
non_redundant_filtered = [];
for g = 1: lnr
    
    %  Get indices from list
    index1 = non_redundant(g,1);
    index2 = non_redundant(g,2);
    %  Get genes from parents
    aa = genome_source_list{index1};
    bb = genome_source_list{index2};
    
    %  Combine them and remove redundancies
    testcase = unique([aa,bb]);
    stc = sort(testcase);
    
    %  Check for illegality(flips too close)
    parity = 1;
    
    for c = 1: length(stc)-1
        if abs(stc(c+1)-stc(c)) >1
            parity = parity * 1;
        else
            parity = 0;
        end
           
    end
    
    %  Put good candidates on the list
    if parity == 1
            non_redundant_filtered =...
            [non_redundant_filtered; [index1 index2]];
            
    end
    
end

non_redundant_filtered = unique(non_redundant_filtered,'rows');
lnrf = length(non_redundant_filtered);
% Randomize indices to facilitate random selection from sorted list
new_indices = non_redundant_filtered(randperm(lnrf),:);
new_pop_indices = new_indices(1:POPULATION - cutoff,:);

%  Build the genome for the next generation
for k = 1: cutoff
    genome{k} = genome_source_list{k};
end

for k = cutoff + 1: POPULATION
    
    npi1 = new_pop_indices(k - cutoff, 1);
    npi2 = new_pop_indices(k - cutoff, 2);
    gsl1 = genome_source_list{npi1};
    gsl2 = genome_source_list{npi2};
    combined = unique([gsl1,gsl2]);
    
    genome{k} = combined;
end

end %  The end of the big evolutionary loop

end

unimpeded_offset_times = zeros(num_ac, 1);
for a = 1: num_ac
    unimpeded_offset_times(a,1) = igs(a).total_time;
end

for b = 1: length(nga.fitness) % loop over population
    
    for a = 1: num_ac  %  loop over aircraft in each solution
        
        if strcmp(igs(a).opType,'DEP') == 1
            % compute unimpeded gate release times to match up with
            % optimized runway times by offsetting with unimpeded
            % taxi times from SOSS-derived database
            nga.gate_release_times(a,1,b) = nga.runway_times(a,2,b) - ...
                unimpeded_offset_times(a,1);
           
        else
            nga.gate_release_times(a,1,b) = -1;
            
        end
    end
end

end


%------------------------------------------------------------------------
%------------------------------------------------------------------------

%  All the functions that are used by the master GA loop


function of = oneflips( initial_aircraft_sequence )
%oneflips Generates the starter set for the GA of one flip neighbors

%  Initialize array
size = length(initial_aircraft_sequence);
of = zeros(size, 6, size - 1);

    % Generate the one flip neighbors
    for i = 1: size - 1

        of(:,:,i) = flipper(initial_aircraft_sequence, i);

    end

end
% ==================================================================
% ==================================================================
function fl = flipper( aircraft_sequence_data, forward_mover_index )
%flipper performs a pairwise nearest-neighbor flip of two aircraft
%   The forward_mover_index identifies wich aircaft moves forward one
%   position.  

first = aircraft_sequence_data(forward_mover_index, :);
second = aircraft_sequence_data(forward_mover_index + 1, :);
fl = aircraft_sequence_data;
fl(forward_mover_index, :) = second;
fl(forward_mover_index+1, :) = first;

end
% ==================================================================
% ==================================================================
function exsearch = exhaustive_searchf( converted_ground_schedule,...
    permutation_array )
%exhaustive_search exhaustively searches up-to-one-position-out-of-place solutions
%   When the solution space is *too small* to allow for a sufficiently 
%   large population to allow a genetic algorithm to work, just enumerate
%   and evaluate all possible aircraft orderings where no aircraft is more
%   than one position out of place.

seq_length = size(converted_ground_schedule, 1);
numperms = length(permutation_array{seq_length});
init_numbering = converted_ground_schedule(:,1);

%  Allocate space for results
stats = zeros(seq_length, 2, numperms);
fitness = zeros(numperms,2);

for j = 1: numperms
   
    order = permutation_array{seq_length}(j,:);
    trial = converted_ground_schedule(order,:);
    seqinfo = sequencer(trial);
    stats(:, 1, j) = seqinfo.start_times;
    stats(:, 2, j) = seqinfo.delays;
    
    %  Sum of squares of delays
    fitness(j, 1) = sum(seqinfo.delays.^2);
    
    % Makespan
    fitness(j, 2) = seqinfo.start_times(end) - ...
        seqinfo.start_times(1);
        
end

if numperms <= 20
    population = numperms;
else
    population = 20;
end

exsearch.fitness = zeros(population, 1);
exsearch.delays = zeros(seq_length, 2, population);
exsearch.runway_times = zeros(seq_length, 2, population);
gate_release_times = -1 * ones(seq_length, 1, population);
spot_begin_end_times = -1 * ones(seq_length, 2, population);

[fitness_ranking,permutation_index] =...
    sort(fitness(:,1));
%  Sort by fitness ranking
%  Note: Change index above from '1' to other to change fitness criteria

exsearch.fitness = fitness_ranking(1: population,1);
exsearch.delays = stats(:,2,permutation_index(1: population));
exsearch.runway_times = stats(:,1,permutation_index(1: population));
exsearch.gate_release_times = gate_release_times;
exsearch.spot_begin_end_times = spot_begin_end_times;

pind = permutation_index(1: population,1);
new_ordering = init_numbering(permutation_array{seq_length}...
    (permutation_index(1: population), :));

for p = 1: population
    
    delays(:,:,p) = sortrows([ new_ordering(p,:)',...
    stats(:, 2, pind(p))],1);
    runway_times(:,:,p) = sortrows([ new_ordering(p,:)',...
    stats(:, 1, pind(p))],1);

end

exsearch.delays = delays;
exsearch.runway_times = runway_times;

end
% ==================================================================
% ==================================================================
function gf = PGM_Ground_Converter( PGM_Parser_Output, ground_database, runway_name )
% PGM_Ground_Converter creates 4DT ground paths from node sequences 
% This function performs what script buildnodetracks.m does, but more
% generally.

num_ac = length(PGM_Parser_Output.aircraftData);
runway_ops_counter = zeros(num_ac,1);

for i = 1: num_ac  %  Filter for ops that use a particular runway
    runway_ops_counter(i) =...
        strcmp(PGM_Parser_Output.aircraftData(i).runway, runway_name);
end
indices = find(runway_ops_counter == 1); % pick out relevant flight nos.
num_ops = length(indices); % total ops of interest

for j = 1: num_ops
    
    gf(j).original_flight_index = indices(j);
    gf(j).opType = PGM_Parser_Output.aircraftData(indices(j)).opType;
    gf(j).taxi_route = PGM_Parser_Output.aircraftData(indices(j)).fullTaxiRoute;
    soss_times = soss_unimpeded_times(PGM_Parser_Output.aircraftData(indices(j)).numId, ...
        ground_database);
    gf(j).sossTimes = soss_times;
    gf(j).total_time = soss_times(end, 1);
    
end

end
% ==================================================================
% ==================================================================
function sust = soss_unimpeded_times( aircraft_ID_num,...
    soss_database )
%soss_unimpeded_times mines database of unimpeded surface traffic trajectories
%  The traffic scenarios used in this project are a combination of
%  aircraft trajectories drawn from a larger set in various combinations.

%  pick out unimpeded (nomSchTimes) times and normalize so start is 0.
%  Note that flight numbers start with zero and require an offset

aci = aircraft_ID_num + 1;
sust = zeros(length(soss_database(aci).taxiRoute), 2);
sust(:, 1) = soss_database(aci).nomSchTimes - ...
    soss_database(aci).nomSchTimes(1);

%  Load up reverse cumulative times into array so end is 0, time reversed
sust(:, 2) = soss_database(aci).nomSchTimes(end) - ...
    soss_database(aci).nomSchTimes;

end
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
function cf = PGM_converter( PGM_Parser_Output, ground_database, runway_name )
%  PGM_converter Converts PGM input into canonical NACRE_GA format
%  Also sorts by order of initial scheduled usage of runway

sizedata = size( PGM_Parser_Output.aircraftData, 2);
runway_ops_counter = zeros(sizedata,1);
ground_array = PGM_Ground_Converter( PGM_Parser_Output, ground_database, runway_name );

for i = 1: sizedata  %  Filter for ops that use a particular runway
    runway_ops_counter(i) =...
        strcmp(PGM_Parser_Output.aircraftData(i).runway, runway_name);
end
num_ops = sum(runway_ops_counter); %  Find total ops on relevant runway

cf = zeros(num_ops,6); % Preallocate for output
j = 0;  % Initialize ops index

for i = 1: sizedata

    %  Only take entries that use the particular runway
    if strcmp(PGM_Parser_Output.aircraftData(i).runway, runway_name) == 1
    j = j+ 1;   
        cf(j,1) = i; %  Original flight index
       
        if PGM_Parser_Output.aircraftData(i).earliestRunwayArrivalTime - ...
               PGM_Parser_Output.aircraftData(i).earliestGateArrivalTime < ...
               ground_array(j).total_time && ...
               strcmp(PGM_Parser_Output.aircraftData(i).opType,'DEP') == 1
            cf(j,2) = ground_array(j).total_time + ...
                   PGM_Parser_Output.aircraftData(i).earliestGateArrivalTime;
              
        else
            cf(j,2) = PGM_Parser_Output.aircraftData(i).earliestRunwayArrivalTime;
        end
        
        % cf(j,2) = PGM_Parser_Output.aircraftData(i).earliestRunwayArrivalTime;
        cf(j,4) = strcmp(PGM_Parser_Output.aircraftData(i).opType,'DEP');
        cf(j,5) = PGM_Parser_Output.aircraftData(i).fix;
        %  Assign weight classes
        if PGM_Parser_Output.aircraftData(i).acWtClass == 3
            cf(j,3) = 1; %  Small
        elseif PGM_Parser_Output.aircraftData(i).acWtClass == 1
            cf(j,3) = 2; %  Large
        else 
            cf(j,3) = 3; %  Heavy
        end
        
        %  For departures, load in earliest release time
        if strcmp(PGM_Parser_Output.aircraftData(i).startTaxiRouteCode , 'GATE')
            cf(j,6) =  PGM_Parser_Output.aircraftData(i).earliestGateArrivalTime;
        else
            cf(j,6) = -1;
        end
        
    end
    %  Move up one space if previous flight used 
    %  chosen runway.
   
             
end
cf = sortrows(cf, 2); %  Sort by initially planned runway usage time

end
% ==================================================================
% ==================================================================
function sequence_report = sequencer( aircraft_ground_schedule_data)
%sequencer creates min makespan schedule incorporating wake vortex

%  sequencer takes a set of aircraft sith scheduled departure times  and 
%  creates the minimum makespan sequence for the aircraft in the specified 
%  order. 

%  Generate constants

num_ac = size(aircraft_ground_schedule_data, 1);

%  Create arrays to store results of sequencer
sequence_report.delays = zeros(num_ac,1);
sequence_report.start_times = zeros(num_ac,1);

% Create matrix for wake vortex separation
% Order of rows and columns: {Small, Large, Heavy, B757}
% Column is leader and row is follower: (later, earlier) or (2nd, 1st))

wake_vortex_matrix_dd =...
    [60	60	120	120;... 
     60	60	120	120;...
     60	60	90	90;...
     60	60	90	90];

wake_vortex_matrix_ad = ... % Arrival after a departure
    [50	60	70	60;...
     50	60	70	60;...
     50	60	70	60;...
     50	60	70	60];

wake_vortex_matrix_da = ... %  Departure after an arrival
    [40	40	40	40;...
     28	28	28	28;...
     24	24	24	24;...
     28	28	28	28];

 wake_vortex_matrix_aa = ... %  Consecutive arrivals, in nautical miles of separation
    [3	5	6;...            %  Does not include 'Super = A380" category.
     3	3	5;...            %  These are ICAO standards.
     3  3   4];
 
%  Initialize with first aircraft time
sequence_report.start_times(1) = aircraft_ground_schedule_data(1,2);

%  Pack aircraft as close as possible, allowing for wake vortex
for i = 1: num_ac - 1
 
    weightclass_first = aircraft_ground_schedule_data(i,3);
    weightclass_second = aircraft_ground_schedule_data(i+1,3);
    optype_first = aircraft_ground_schedule_data(i,4);
    optype_second = aircraft_ground_schedule_data(i+1,4);
    fix_first = aircraft_ground_schedule_data(i,5);
    fix_second = aircraft_ground_schedule_data(i+1,5);
      
    if optype_first == 1 && optype_second == 1 &&...
            weightclass_first < 3 && fix_first == fix_second  
        
        min_sep = 80;
        %  Some weight combos of A/C going to the same departure fix need
        %  more time than specified by wake vortex rules because of 
        %  separation rules if they're headed in the same direction
        
    elseif optype_first == 1 && optype_second == 1 &&...
            weightclass_first >= 3 && fix_first == fix_second
        
        min_sep = wake_vortex_matrix_dd(weightclass_second,weightclass_first);
        %  If the lead aircraft is bigger and they're both going to the
        %  same fix, wake vortex rules dominate
    
    elseif optype_first == 1 && optype_second == 1 &&...
            fix_first ~= fix_second
        
        min_sep = wake_vortex_matrix_dd(weightclass_second,weightclass_first);
        %  If they're not going to the same fix, vortex rules apply
        
    elseif optype_first == 1 && optype_second == 0
       
        min_sep = wake_vortex_matrix_da(weightclass_second,weightclass_first);
        
    elseif optype_first == 0 && optype_second == 1
        
        min_sep = wake_vortex_matrix_ad(weightclass_second,weightclass_first);
        
    else  %  In the rare case of consecutive arrivals
        
        min_sep = 80;
        
        % wake_vortex_matrix_aa(weightclass_second,weightclass_first) * 30;
        %  SOSS only has a requirement when coming from same fix, but it
        %  seems this needs to be extended to any fix.  ICAO distance-based
        %  separation gives a minimum separation from 72-90 sec.
        %  Assume flight speed of 120 knots to convert distance separation
        %  standards into time separations.  150 knots nominal speed would 
        %  have a multiplier of 24 seconds/nm.
       
    end
    
    %  Check to see if next flight is constrained by schedule or by
    %  separation, increment start_times appropriately
    if sequence_report.start_times(i) + min_sep >...
            aircraft_ground_schedule_data(i+1,2)
        
        sequence_report.start_times(i+1) =...
            sequence_report.start_times(i) + min_sep;
        
    else
        
        sequence_report.start_times(i + 1) = ...
            aircraft_ground_schedule_data(i + 1, 2);
        
    end
    
end

for t = 1: num_ac
    
    sequence_report.delays(t) = sequence_report.start_times(t) - ...
        aircraft_ground_schedule_data(t, 2);
    
end

end



