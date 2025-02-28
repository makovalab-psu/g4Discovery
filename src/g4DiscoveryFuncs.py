from Bio import SeqIO
import numpy as np
import os
import pandas as pd
import re
import subprocess

def run_docker(input_file_path, output_file_path, output_file_name, pqs_min_score, container_name='kxk302/pqsfinder:1.0.0'):
    '''
    Run the pqsfinder Docker container with the specified input and output files.
    '''
    try:
        command = [
            "docker", "run", 
            "-v", f"{os.path.abspath(input_file_path)}:/input",  # Mount the input file
            "-v", f"{os.path.abspath(output_file_path)}:/output",  # Mount the output file
            container_name,  # Docker container name
            "/input",  # Path to input file inside the container
            f"/output/{output_file_name}",  # Path to output file inside the container
            f"{pqs_min_score}",
            "1" #overlapping = True
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=True)

    except subprocess.CalledProcessError as e:
        print(e.stderr)  # Standard error from the container

# regex expressions derived from Guedin A. et al. (2010), Mukundan V.T. et al. (2013), and Maity A. et al. (2020)
pos_pattern = r'((G{3,}[ATCG]{1,12}){3,}G{3,})|((G([ATC]{0,1})G([ATC]{0,1})G([ATCG]{1,3})){3,}G([ATC]{0,1})G([ATC]{0,1})G)|((G{1,2}[ATC]{1,2}){7,}G{1,2})'
neg_pattern = r'((C{3,}[ATCG]{1,12}){3,}C{3,})|((C([ATG]{0,1})C([ATG]{0,1})C([ATCG]{1,3})){3,}C([ATG]{0,1})C([ATG]{0,1})C)|((C{1,2}[ATG]{1,2}){7,}C{1,2})'

def filterG4s(fasta_file, chr, min_tetrad, min_score, min_g4hunterscore, pos_pattern=pos_pattern, neg_pattern=neg_pattern):
    '''
    Filter G4s from a fasta file based on the following criteria:
    1. Minimum number of tetrads
    2. Minimum G4Hunter score
    3. Minimum PQS score
    4. Sequence must be a valid G4 motif
    '''
    choose = []
    for record in SeqIO.parse(fasta_file, "fasta"):

        params = record.id.split(";")
        sequence = (str(record.seq)).upper()

        start = int(params[2].split("=")[-1]) - 1 #convert to 0-based
        end = int(params[3].split("=")[-1]) #convert to 0-based exclusive [start, end)
        strand = params[4].split("=")[-1]
        score = float(params[5].split("=")[-1])
        nt = int(params[6].split("=")[-1]) #number of tetrads
        length = end-start

        if strand == "+":
            regex = re.match(pos_pattern, sequence)
        elif strand == "-":
            regex = re.match(neg_pattern, sequence)

        g4HunterScore = np.round(float(CalScore(BaseScore(sequence)[1],length)[0]),2) #G4 Hunter score calculation

        if nt >= min_tetrad and score >= min_score and regex != None and regex.span() == (0, len(sequence)) and abs(g4HunterScore) >= min_g4hunterscore:
            choose.append([f"chr{chr}", start, end, score, length, strand, g4HunterScore])

    choose = pd.DataFrame(choose, columns=["chr","start", "end", "score", "length", "strand", "g4HunterScore"])
    return choose

def filterNonOverlappingG4sCmplx(df):
    '''
    Filter non-overlapping G4s in a region by 
    choosing the one with the highest pqsscore 
    and shortest length.
    '''
    grouped = df.groupby("start") #group by start position

    filteredG4s = [] 
    for _, group in grouped: #iterate thorught the start positions
        group = group.sort_values(by=["score","length"],ascending=[False,True]) #choose the shortest and best scoring G4
        row = group.iloc[0] #the first one
        filteredG4s.append(row.name) #indicate where in the dataframe the row is
    df = df.loc[filteredG4s] #filter the dataframe
    df = df.sort_values(by=["start"],ascending=[True]) #sort by start position
    
    if df.empty:
        return df
    else:
        df_workon = df.copy() #copy the dataframe to work on
        choose = [] 
        finalData = []
        rowendgroup = 0 #initialize the end of the region
        while df_workon.empty == False and df_workon.iloc[-1]["start"] > rowendgroup: # while the last start position is greater than the end of the region
            for idx, row in df_workon.iterrows(): #iterate through the dataframe
                rowendgroup = df_workon.iloc[0]["end"] #set the end of the region to the end of the first row
                if row["start"] >= rowendgroup: #if the start position is greater than the end of the region
                    subset = df_workon.loc[choose] #choose the rows
                    subset = subset.sort_values(by=["score","length"],ascending=[False,True]) #sort by score and length
                    choosedRow = subset.iloc[0] #choose the first row
                    df_workon = df[df["start"] >= choosedRow["end"]] #filter the dataframe
                    finalData.append(choosedRow.name) #append the index
                    choose = [] #reset the choose list
                    break #break the for loop
                else:
                    choose.append(idx) #append the index until if condition is met

        #For the last subset
        if not df_workon.empty:
            subset = df_workon.loc[choose] 
            subset = subset.sort_values(by=["score","length"],ascending=[False,True])
            choosedRow = subset.iloc[0]
        finalData.append(choosedRow.name)

        df = df.loc[finalData]
        df.drop_duplicates(inplace=True) #to account for duplicate if last subset result is same
        
        return df
