import pandas as pd
import numpy as np
import re
import subprocess
import argparse
import os
import random
import json
import time
from typing import List, Dict, Tuple
from exportable_custom_functions import setup_experimental_data

class GeneticOptimizer:
    def __init__(
        self,
        population_size=50,
        generations=30,
        mutation_rate=0.2,
        crossover_rate=0.8,
        tournament_size=3,
        n_processes=os.cpu_count()
    ):
        self.population_size = population_size
        self.generations = generations
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.tournament_size = tournament_size
        self.n_processes = n_processes
        
        # Constants
        self.initial_O2 = 5.7E22
        self.flux = 2.2E14
        
        # Setup experimental data
        self.exp_data = setup_experimental_data()
        
        # Calculate weights that increase linearly from first to last point
        self.calculate_weights()
        
        # Load parameter ranges
        self.model_inp_ranges = pd.read_csv("./model_inp_values/model_ranges.dat")
        self.photo_processes_ranges = pd.read_csv("./photo_processes_values/photo_processes_ranges.dat")
        self.variables = pd.concat([self.model_inp_ranges, self.photo_processes_ranges], ignore_index=True)
        
        # Ensure work directories exist
        self.setup_directories()
    
    def calculate_weights(self):
        """Calculate increasing weights for experimental data points"""
        num_points = len(self.exp_data)
        # Linear weight increase from 1.0 to 5.0
        self.weights = np.linspace(1.0, 5.0, num_points)
        # Normalize weights so they sum to num_points (to keep RMSD comparable)
        self.weights = self.weights * (num_points / np.sum(self.weights))
        
        # Add weights to experimental data for reference
        self.exp_data['weight'] = self.weights
        print(f"Calculated weights for {num_points} experimental points:")
        for i, w in enumerate(self.weights):
            print(f"  Point {i+1}: {w:.2f}")
        
    def setup_directories(self):
        """Create necessary directories for parallel processing"""
        # Create a main work directory
        os.makedirs("ga_workdir", exist_ok=True)
        
        # Verify monaco exists and has execute permission
        monaco_path = os.path.abspath("./monaco")
        if not os.path.exists(monaco_path):
            print(f"ERROR: Monaco executable not found at {monaco_path}")
            print(f"Current directory: {os.getcwd()}")
            print(f"Directory contents: {os.listdir('.')}")
            raise FileNotFoundError(f"Monaco executable not found: {monaco_path}")
        
        if not os.access(monaco_path, os.X_OK):
            print(f"Making monaco executable: {monaco_path}")
            os.chmod(monaco_path, 0o755)  # Add execute permission
        
        # List all files in the current directory
        print("Files in current directory:")
        all_files = os.listdir('.')
        print(all_files)
        
        # List of files that monaco might need
        essential_files = [
            "model.inp", 
            "photo_processes.dat",
            "network.dat",
            "rd_eff.txt",
            "init_gas_ab.inp",   # Added as requested
            "init_bulk_ab.inp",  # Added as requested
            "init_surf_ab.inp"   # Added as requested
        ]
        
        # Look for any additional data files that might be needed
        data_file_extensions = ['.dat', '.txt', '.inp', '.ini', '.cfg']
        for file in all_files:
            _, ext = os.path.splitext(file)
            if ext.lower() in data_file_extensions and file not in essential_files:
                essential_files.append(file)
                
        print(f"Essential files to copy: {essential_files}")
        
        # Create individual process directories
        for i in range(self.n_processes):
            process_dir = f"ga_workdir/process_{i}"
            os.makedirs(process_dir, exist_ok=True)
            
            # Create csv directory within process directory
            os.makedirs(f"{process_dir}/csv", exist_ok=True)
            
            # Copy all essential files
            for file in essential_files:
                if os.path.exists(file):
                    try:
                        import shutil
                        dest_path = os.path.join(process_dir, file)
                        shutil.copy2(file, dest_path)
                        print(f"Copied {file} to {process_dir}")
                    except Exception as e:
                        print(f"Error copying {file}: {e}")
                else:
                    print(f"Warning: Essential file {file} not found")
            
            # Copy additional required directories
            required_dirs = []  # Add directory names monaco might need
            for dir_name in required_dirs:
                if os.path.exists(dir_name):
                    try:
                        import shutil
                        dest_dir = os.path.join(process_dir, dir_name)
                        if os.path.exists(dest_dir):
                            shutil.rmtree(dest_dir)
                        shutil.copytree(dir_name, dest_dir)
                        print(f"Copied directory {dir_name} to {process_dir}")
                    except Exception as e:
                        print(f"Error copying directory {dir_name}: {e}")
            
            # Copy the monaco executable directly (safer than symlinks)
            try:
                import shutil
                monaco_dest = os.path.join(process_dir, "monaco")
                shutil.copy2(monaco_path, monaco_dest)  # preserves permissions
                print(f"Copied monaco to {process_dir}")
            except Exception as e:
                print(f"Error copying monaco: {e}")
                
        print("Directory setup complete")
        
    def generate_individual(self) -> List[float]:
        """Generate a random individual (normalized values 0-1)"""
        return [random.random() for _ in range(len(self.variables))]
    
    def generate_population(self, size: int) -> List[List[float]]:
        """Generate initial population"""
        return [self.generate_individual() for _ in range(size)]
    
    def scale_parameters(self, individual: List[float]) -> Dict[str, float]:
        """Scale normalized parameters (0-1) to their actual ranges"""
        scaled_params = {}
        for i, (_, row) in enumerate(self.variables.iterrows()):
            min_val = np.log10(row['minValue'])
            max_val = np.log10(row['maxValue'])
            scaled_val = 10 ** (min_val + individual[i] * (max_val - min_val))
            scaled_params[row['toModify']] = scaled_val
        return scaled_params

    def update_model_files(self, params: Dict[str, float], process_dir: str):
        """Update model.inp and photo_processes.dat with new parameters"""
        # Update model.inp
        filepath = f"{process_dir}/model.inp"
        with open(filepath, "r") as file:
            lines = file.readlines()
        with open(filepath, "w") as file:
            for line in lines:
                for param_name, new_value in params.items():
                    if re.search(rf"\b{param_name}\b", line):
                        formatted_value = f"{new_value:.4E}"
                        line = re.sub(r"=\s+[-+]?\d+\.\d+E[-+]?\d+", 
                                    f"= {formatted_value}", line)
                file.write(line)
        
        # Update photo_processes.dat
        filepath = f"{process_dir}/photo_processes.dat"
        with open(filepath, "r") as file:
            lines = file.readlines()
        with open(filepath, "w") as file:
            for line in lines:
                for param_name, new_value in params.items():
                    if re.search(rf"\b{param_name}\b", line):
                        formatted_value = f"{new_value:.2E}"
                        start_col = 107
                        end_col = 118
                        line = (line[:start_col - 1] + formatted_value + 
                               line[end_col:] + "\n")
                file.write(line)

    def evaluate_individual(self, individual: List[float], process_id: int) -> float:
        """Evaluate fitness of an individual using weighted RMSD"""
        process_dir = f"ga_workdir/process_{process_id}"
        
        # Scale parameters to their actual ranges
        params = self.scale_parameters(individual)
        
        # Update model files with new parameters
        self.update_model_files(params, process_dir)
        
        # Run monaco in the process directory
        current_dir = os.getcwd()
        try:
            # Full diagnostics of the environment
            print(f"\nProcess {process_id} diagnostics:")
            print(f"  Current directory before chdir: {current_dir}")
            print(f"  Process directory: {process_dir}")
            print(f"  Process directory exists: {os.path.exists(process_dir)}")
            print(f"  Process directory contents: {os.listdir(process_dir)}")
            
            os.chdir(process_dir)
            print(f"  Current directory after chdir: {os.getcwd()}")
            
            # Verify monaco exists and is executable
            monaco_exe = "./monaco"
            if not os.path.exists(monaco_exe):
                print(f"  ERROR: Monaco not found in {os.getcwd()}")
                raise FileNotFoundError(f"Monaco executable not found in {os.getcwd()}")
                
            if not os.access(monaco_exe, os.X_OK):
                print(f"  Making monaco executable in {os.getcwd()}")
                os.chmod(monaco_exe, 0o755)
            
            # Run monaco with timeout (30 seconds per evaluation)
            print(f"  Executing: {monaco_exe} in {os.getcwd()}")
            try:
                result = subprocess.run(monaco_exe, capture_output=True, text=True, timeout=30)
                
                # Log the execution results
                print(f"  Monaco return code: {result.returncode}")
                print(f"  Monaco stdout (truncated): {result.stdout[:500]}")
                print(f"  Monaco stderr (truncated): {result.stderr[:500]}")
            except subprocess.TimeoutExpired:
                print(f"  Monaco execution timed out after 30 seconds")
                # Don't raise an exception, we'll check for output files anyway
            
            # Check directory structure after execution  
            print(f"  Directory contents after execution: {os.listdir('.')}")
            print(f"  CSV directory exists: {os.path.exists('csv')}")
            if os.path.exists('csv'):
                print(f"  CSV directory contents: {os.listdir('csv')}")
            
            # Ensure the csv directory exists
            os.makedirs("csv", exist_ok=True)
            
            # Check if the bO3.csv file exists
            if not os.path.exists("csv/bO3.csv"):
                raise FileNotFoundError("csv/bO3.csv not found")
                
            # Read and process the CSV file safely
            try:
                # First check file content
                with open("csv/bO3.csv", 'r') as f:
                    content = f.read()
                    lines = content.strip().split('\n')
                    print(f"  bO3.csv has {len(lines)} lines")
                    print(f"  First few lines: {lines[:3]}")
                
                # If file has only header or is empty, raise an exception
                if len(lines) <= 2:  # Need at least header + data row
                    raise ValueError(f"bO3.csv has insufficient data: {len(lines)} lines")
                
                # Read CSV with pandas - flexible header handling
                if len(lines) > 2:
                    calc_data = pd.read_csv("csv/bO3.csv", header=1, 
                                           names=["Fluence", "Abundance"])
                else:
                    # Try without header
                    calc_data = pd.read_csv("csv/bO3.csv", header=None, 
                                           names=["Fluence", "Abundance"])
                
                # Verify we have data
                if len(calc_data) == 0:
                    raise ValueError("Empty dataframe after reading bO3.csv")
                    
                print(f"  Successfully read {len(calc_data)} rows from bO3.csv")
                
                calc_data["Abundance"] = (calc_data["Abundance"] / self.initial_O2) * 100.0
                
                # Calculate weighted RMSD
                weighted_deviations = []
                for idx in self.exp_data.index:
                    scalar = self.exp_data.loc[idx, "expX"]
                    closest_index = (calc_data["Fluence"] - scalar).abs().idxmin()
                    closest_value = calc_data.loc[closest_index, "Abundance"]
                    
                    # Get the weight for this experimental point
                    weight = self.weights[idx]
                    
                    # Calculate weighted squared deviation
                    weighted_deviation = weight * ((closest_value - self.exp_data.loc[idx, "expY"]) ** 2)
                    weighted_deviations.append(weighted_deviation)
                    
                    # Debug print for this comparison
                    print(f"  Point {idx}: Model={closest_value:.2f}, Exp={self.exp_data.loc[idx, 'expY']:.2f}, "
                          f"Weight={weight:.2f}, Weighted deviation={weighted_deviation:.2f}")
                
                # Calculate weighted RMSD
                weighted_rmsd = (sum(weighted_deviations) / sum(self.weights)) ** 0.5
                print(f"  Process {process_id} evaluation successful. Weighted RMSD: {weighted_rmsd:.6f}")
                
            except Exception as e:
                print(f"  Error processing CSV: {e}")
                import traceback
                traceback.print_exc()
                weighted_rmsd = float('inf')
            
        except Exception as e:
            print(f"Error during evaluation in process {process_id}: {e}")
            import traceback
            traceback.print_exc()
            weighted_rmsd = float('inf')  # Penalty value for failed evaluations
        finally:
            os.chdir(current_dir)
            
        return weighted_rmsd
        
    def evaluate_population_batch(self, batch_id: int, individuals: List[List[float]]):
        """Evaluate a batch of individuals in a separate process"""
        results = []
        for i, ind in enumerate(individuals):
            fitness = self.evaluate_individual(ind, batch_id)
            results.append((i, fitness))
            
        # Write results to a temporary file
        with open(f"ga_workdir/results_{batch_id}.json", "w") as f:
            json.dump(results, f)
    
    def batch_individuals(self, population):
        """Split population into batches for parallel processing"""
        batch_size = len(population) // self.n_processes
        batches = []
        
        for i in range(self.n_processes):
            start_idx = i * batch_size
            end_idx = start_idx + batch_size if i < self.n_processes - 1 else len(population)
            batches.append(population[start_idx:end_idx])
            
        return batches
    
    def tournament_selection(self, population, fitnesses, tournament_size):
        """Select an individual using tournament selection"""
        indices = random.sample(range(len(population)), tournament_size)
        tournament = [(population[i], fitnesses[i]) for i in indices]
        winner = min(tournament, key=lambda x: x[1])
        return winner[0]
    
    def crossover(self, parent1, parent2):
        """Perform two-point crossover"""
        if random.random() > self.crossover_rate:
            return parent1, parent2
            
        size = len(parent1)
        point1 = random.randint(1, size-2)
        point2 = random.randint(point1, size-1)
        
        child1 = parent1[:point1] + parent2[point1:point2] + parent1[point2:]
        child2 = parent2[:point1] + parent1[point1:point2] + parent2[point2:]
        
        return child1, child2
    
    def mutate(self, individual):
        """Perform Gaussian mutation"""
        mutated = individual.copy()
        for i in range(len(mutated)):
            if random.random() < self.mutation_rate:
                # Add Gaussian noise
                mutated[i] += random.gauss(0, 0.2)
                # Ensure values stay within 0-1 range
                mutated[i] = max(0, min(1, mutated[i]))
        return mutated
    
    def parallel_evaluate(self, population):
        """Evaluate population in parallel using subprocess"""
        # Split population into batches
        batches = self.batch_individuals(population)
        
        # Launch evaluation processes
        processes = []
        for i, batch in enumerate(batches):
            # Write batch to file
            with open(f"ga_workdir/batch_{i}.json", "w") as f:
                json.dump(batch, f)
            
            # Launch the process using the same Python interpreter
            cmd = [
                "python3", "-c",
                f"""
import sys, json
from {os.path.splitext(os.path.basename(__file__))[0]} import GeneticOptimizer
# Load the batch
with open("ga_workdir/batch_{i}.json", "r") as f:
    batch = json.load(f)
# Create optimizer with same parameters
optimizer = GeneticOptimizer(
    population_size={self.population_size},
    generations={self.generations},
    mutation_rate={self.mutation_rate},
    crossover_rate={self.crossover_rate},
    tournament_size={self.tournament_size},
    n_processes={self.n_processes}
)
# Evaluate the batch
optimizer.evaluate_population_batch({i}, batch)
                """
            ]
            
            process = subprocess.Popen(cmd)
            processes.append(process)
        
        # Wait for all processes to complete
        for p in processes:
            p.wait()
        
        # Collect results
        all_results = []
        for i in range(len(batches)):
            with open(f"ga_workdir/results_{i}.json", "r") as f:
                results = json.load(f)
                
            # Calculate original indices
            batch_size = len(population) // self.n_processes
            start_idx = i * batch_size
            
            # Adjust indices
            adjusted_results = [(start_idx + idx, fitness) for idx, fitness in results]
            all_results.extend(adjusted_results)
        
        # Sort by original index and extract just the fitness values
        all_results.sort(key=lambda x: x[0])
        fitnesses = [result[1] for result in all_results]
        
        return fitnesses
    
    def run_optimization(self):
        """Run the genetic algorithm optimization"""
        print("Initializing population...")
        population = self.generate_population(self.population_size)
        
        # Track the best solution
        best_individual = None
        best_fitness = float('inf')
        
        # Main evolution loop
        for gen in range(self.generations):
            print(f"\nGeneration {gen+1}/{self.generations}")
            start_time = time.time()
            
            # Evaluate population
            print(f"Evaluating {len(population)} individuals...")
            fitnesses = self.parallel_evaluate(population)
            
            # Update best solution
            min_fitness_idx = fitnesses.index(min(fitnesses))
            if fitnesses[min_fitness_idx] < best_fitness:
                best_fitness = fitnesses[min_fitness_idx]
                best_individual = population[min_fitness_idx]
                
            print(f"Best fitness: {best_fitness:.6f}")
            
            # Create next generation
            if gen < self.generations - 1:  # Skip for last generation
                new_population = []
                
                # Elitism: keep the best individual
                elite_idx = fitnesses.index(min(fitnesses))
                new_population.append(population[elite_idx])
                
                # Fill the rest of the population
                while len(new_population) < self.population_size:
                    # Selection
                    parent1 = self.tournament_selection(population, fitnesses, self.tournament_size)
                    parent2 = self.tournament_selection(population, fitnesses, self.tournament_size)
                    
                    # Crossover
                    child1, child2 = self.crossover(parent1, parent2)
                    
                    # Mutation
                    child1 = self.mutate(child1)
                    child2 = self.mutate(child2)
                    
                    # Add to new population
                    new_population.append(child1)
                    if len(new_population) < self.population_size:
                        new_population.append(child2)
                
                population = new_population
            
            gen_time = time.time() - start_time
            print(f"Generation completed in {gen_time:.2f} seconds")
            
        # Final results
        best_params = self.scale_parameters(best_individual)
        
        # Save results
        results_df = pd.DataFrame([best_params])
        results_df['weighted_RMSD'] = best_fitness
        results_df.to_csv("genetic_algorithm_results.csv")
        
        return best_params, best_fitness

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Run genetic algorithm optimization for model parameters')
    
    parser.add_argument('--pop-size', type=int, default=50,
                        help='Population size (default: 50)')
    parser.add_argument('--generations', type=int, default=30,
                        help='Number of generations (default: 30)')
    parser.add_argument('--mutation-rate', type=float, default=0.2,
                        help='Mutation rate (default: 0.2)')
    parser.add_argument('--crossover-rate', type=float, default=0.8,
                        help='Crossover rate (default: 0.8)')
    parser.add_argument('--tournament-size', type=int, default=3,
                        help='Tournament size for selection (default: 3)')
    parser.add_argument('--processes', type=int, default=os.cpu_count(),
                        help=f'Number of parallel processes (default: {os.cpu_count()})')
    parser.add_argument('--output', type=str, default='genetic_algorithm_results.csv',
                        help='Output file name (default: genetic_algorithm_results.csv)')
    
    return parser.parse_args()

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_arguments()
    
    print("Starting optimization with parameters:")
    print(f"Population size: {args.pop_size}")
    print(f"Generations: {args.generations}")
    print(f"Mutation rate: {args.mutation_rate}")
    print(f"Crossover rate: {args.crossover_rate}")
    print(f"Tournament size: {args.tournament_size}")
    print(f"Number of processes: {args.processes}")
    print(f"Output file: {args.output}")
    
    # Initialize optimizer with command line arguments
    optimizer = GeneticOptimizer(
        population_size=args.pop_size,
        generations=args.generations,
        mutation_rate=args.mutation_rate,
        crossover_rate=args.crossover_rate,
        tournament_size=args.tournament_size,
        n_processes=args.processes
    )
    
    # Run optimization
    print("\nStarting optimization...")
    best_params, best_rmsd = optimizer.run_optimization()
    
    print(f"\nOptimization complete!")
    print(f"Best weighted RMSD: {best_rmsd}")
    print("\nBest parameters:")
    for param, value in best_params.items():
        print(f"{param}: {value:.4E}")
