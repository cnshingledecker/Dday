import pandas as pd
import numpy as np
import subprocess
import os
import time
from exportable_custom_functions import setup_experimental_data

def evaluate_current_model():
    """
    Run monaco executable in the current directory and calculate 
    weighted RMSD compared to experimental data.
    """
    print("Evaluating current model parameters...")
    
    # Constants
    initial_O2 = 5.7E22
    
    # Setup experimental data
    exp_data = setup_experimental_data()
    
    # Calculate weights that increase linearly from first to last point
    num_points = len(exp_data)
    # Linear weight increase from 1.0 to 5.0
    weights = np.linspace(1.0, 5.0, num_points)
    # Normalize weights so they sum to num_points (to keep RMSD comparable)
    weights = weights * (num_points / np.sum(weights))
    
    # Add weights to experimental data for reference
    exp_data['weight'] = weights
    print(f"Weights for {num_points} experimental points:")
    for i, w in enumerate(weights):
        print(f"  Point {i+1}: {w:.2f}")
    
    # Verify monaco exists and is executable
    monaco_exe = "./monaco"
    if not os.path.exists(monaco_exe):
        raise FileNotFoundError(f"Monaco executable not found in {os.getcwd()}")
    
    if not os.access(monaco_exe, os.X_OK):
        print(f"Making monaco executable")
        os.chmod(monaco_exe, 0o755)
    
    # Ensure csv directory exists
    os.makedirs("csv", exist_ok=True)
    
    print(f"Executing monaco in {os.getcwd()}")
    start_time = time.time()
    
    try:
        # Run monaco with timeout (60 seconds)
        result = subprocess.run(monaco_exe, capture_output=True, text=True, timeout=60)
        
        # Log the execution results
        print(f"Monaco return code: {result.returncode}")
        print(f"Monaco execution time: {time.time() - start_time:.2f} seconds")
        
        if result.returncode != 0:
            print("Warning: Monaco execution returned non-zero exit code")
            print(f"Monaco stderr: {result.stderr[:500]}")
    except subprocess.TimeoutExpired:
        print(f"Monaco execution timed out after 60 seconds")
        # Continue to check for output files anyway
    
    # Check if the bO3.csv file exists
    if not os.path.exists("csv/bO3.csv"):
        raise FileNotFoundError("csv/bO3.csv not found after running monaco")
    
    # Read and process the CSV file
    try:
        # First check file content
        with open("csv/bO3.csv", 'r') as f:
            content = f.read()
            lines = content.strip().split('\n')
            print(f"bO3.csv has {len(lines)} lines")
        
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
            
        print(f"Successfully read {len(calc_data)} rows from bO3.csv")
        
        # Convert abundance to percentage
        calc_data["Abundance"] = (calc_data["Abundance"] / initial_O2) * 100.0
        
        # Calculate weighted RMSD
        weighted_deviations = []
        unweighted_deviations = []
        results_data = []
        
        for idx in exp_data.index:
            scalar = exp_data.loc[idx, "expX"]
            closest_index = (calc_data["Fluence"] - scalar).abs().idxmin()
            closest_value = calc_data.loc[closest_index, "Abundance"]
            
            # Get the weight for this experimental point
            weight = weights[idx]
            
            # Calculate deviations
            deviation = (closest_value - exp_data.loc[idx, "expY"]) ** 2
            weighted_deviation = weight * deviation
            
            weighted_deviations.append(weighted_deviation)
            unweighted_deviations.append(deviation)
            
            # Store results for detailed output
            results_data.append({
                "Point": idx + 1,
                "Fluence": scalar,
                "Model_Value": closest_value,
                "Exp_Value": exp_data.loc[idx, "expY"],
                "Difference": closest_value - exp_data.loc[idx, "expY"],
                "Weight": weight,
                "Weighted_Deviation": weighted_deviation
            })
            
        # Calculate RMSDs
        weighted_rmsd = (sum(weighted_deviations) / sum(weights)) ** 0.5
        unweighted_rmsd = (sum(unweighted_deviations) / len(exp_data)) ** 0.5
        
        # Create and display results table
        results_df = pd.DataFrame(results_data)
        
        print("\nDetailed comparison results:")
        print(results_df.to_string(index=False))
        
        print("\nSummary:")
        print(f"Unweighted RMSD: {unweighted_rmsd:.6f}")
        print(f"Weighted RMSD: {weighted_rmsd:.6f}")
        
        # Save results to CSV
        results_df.to_csv("current_model_evaluation.csv", index=False)
        
        # Create a summary dataframe with both RMSD values
        summary_df = pd.DataFrame({
            'Metric': ['Unweighted_RMSD', 'Weighted_RMSD'],
            'Value': [unweighted_rmsd, weighted_rmsd]
        })
        summary_df.to_csv("current_model_summary.csv", index=False)
        
        return weighted_rmsd, unweighted_rmsd, results_df
        
    except Exception as e:
        print(f"Error processing results: {e}")
        import traceback
        traceback.print_exc()
        return float('inf'), float('inf'), None

if __name__ == "__main__":
    try:
        weighted_rmsd, unweighted_rmsd, results = evaluate_current_model()
        print("\nEvaluation complete!")
        
        # Plot results if matplotlib is available
        try:
            import matplotlib.pyplot as plt
            
            if results is not None:
                plt.figure(figsize=(12, 8))
                
                # Plot model vs experimental values
                plt.subplot(2, 1, 1)
                plt.plot(results['Fluence'], results['Model_Value'], 'b-o', label='Model')
                plt.plot(results['Fluence'], results['Exp_Value'], 'r-x', label='Experimental')
                plt.xlabel('Fluence')
                plt.ylabel('Abundance (%)')
                plt.title('Model vs Experimental Values')
                plt.legend()
                plt.grid(True)
                
                # Plot weights and deviations
                plt.subplot(2, 1, 2)
                plt.bar(range(1, len(results) + 1), results['Weight'], alpha=0.3, label='Weight')
                plt.bar(range(1, len(results) + 1), results['Weighted_Deviation'], alpha=0.7, label='Weighted Deviation')
                plt.xlabel('Data Point')
                plt.title(f'Weights and Deviations (Weighted RMSD: {weighted_rmsd:.4f})')
                plt.legend()
                plt.grid(True)
                
                plt.tight_layout()
                plt.savefig('model_evaluation.png')
                print("Created visualization: model_evaluation.png")
        except ImportError:
            print("Matplotlib not available - skipping visualization")
        
    except Exception as e:
        print(f"Evaluation failed: {e}")
        import traceback
        traceback.print_exc()
