import os
import pandas as pd
import sys
import time

# Set the directory containing the CSV files
csv_dir = sys.argv[1]

# Get a list of all CSV files in the directory
csv_files = [f for f in os.listdir(csv_dir) if f.endswith("_counts.csv")]

# Initialize an empty dictionary to hold the total counts
total_counts = {}

# Loop through each CSV file
for csv_file in csv_files:
    # Read in the CSV file as a DataFrame
    df = pd.read_csv(os.path.join(csv_dir, csv_file), sep="\t")
    
    # Loop through each row in the DataFrame
    for index, row in df.iterrows():
        # Get the keyword and count
        keyword = row["keyword"]
        count = row["count"]
        
        # Add the count to the total for this keyword
        if keyword in total_counts:
            total_counts[keyword] += count
        else:
            total_counts[keyword] = count

# Print the total counts for each keyword
for keyword, count in total_counts.items():
    print(f"{keyword}: {count}")

timestr = time.strftime("%Y%m%d-%H%M%S")
print(total_counts)

counts_df = pd.DataFrame.from_dict(total_counts, orient='index', columns=['count'])
counts_df.to_csv(f'keyword_counts_{timestr}.csv')

print(f'Written to file: keyword_counts_{timestr}.csv')

