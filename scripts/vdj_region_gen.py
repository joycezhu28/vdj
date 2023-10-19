import pandas as pd
import numpy as np
from difflib import SequenceMatcher
import argparse
pd.options.mode.chained_assignment = None  # default='warn'

#How to use in command line:
#python3 vdj_region_gen.py (vj_file) -v (v_file) -j (j_file) -o (output name)

#Function creates ID column matching that of the vj file
def combine_columns(df, col1, col2):
    '''
    Inputs:
        df = v or j file
        col1 = accession number ('AN')
        col2 = 'V.Name' or 'J.Name'
    Output:
        df updated with VID/JID column
    '''
    col_one = df[col1].astype(str)
    col_two = df[col2].astype(str)
    new_column = col_one + '|' + col_two + '|Homo'
    df['ID'] = new_column
    return df

#Function matches VID and JID to extract corresponding AA sequences
def id_match(vj, v, j):
    '''
    Inputs:
        vj = df containing both V and J IDs in 'VID' and 'JID' columns
        v = df containing VIDs and corresponding AA sequences
        j = df containing JIDs and corresponding AA sequences
    Output:
        new df updated with AA sequences
    '''
    vj['VID'] = vj['VID'].astype(str)
    vj['JID'] = vj['JID'].astype(str)
    v_seq_list = []
    j_seq_list = []
    for i, row in vj.iterrows():
        v_id = row['VID']
        j_id = row['JID']
        try:
            filtered_v = v[v['ID'] == v_id].reset_index()
            v_seq = filtered_v.loc[0]['V.AA.String']
        except:
            v_seq = 'NA'
        try:
            filtered_j = j[j['ID'] == j_id].reset_index()
            j_seq = filtered_j.loc[0]['J.AA.String']
        except:
            j_seq = 'NA'
        v_seq_list.append(v_seq)
        j_seq_list.append(j_seq)
    new_df = vj.copy()
    new_df['V_AA_Seq'] = v_seq_list
    new_df['J_AA_Seq'] = j_seq_list
    return new_df

#Finds longest common sequence between V & CDR3 and keeps only the preceding V sequence
#Requires longest common sequence to start with a 'C'
def v_trimming(v_seq, cdr3_seq):
    '''
    Inputs:
        v_seq = V sequence
        cdr3_seq = CDR3 sequence
    Output:
        v_seq_trimmed where the longest common sequence between V & CDR3 starting with C is removed along with the sequence following it
    '''
    match = SequenceMatcher(None, v_seq, cdr3_seq).find_longest_match()
    match_seq = v_seq[match.a:match.a + match.size]
    if match_seq[0] == 'C':
        v_seq_parts = v_seq.split(match_seq)[:-1]
        v_seq_trimmed = ''.join(v_seq_parts)
    else:
        all_v_strings = [v_seq[i:k+1] for i in range(len(v_seq)) for k in range(len(v_seq))]
        common_substrings = [x for x in all_v_strings if x in v_seq and x in cdr3_seq and x.startswith('C')]
        ordered_cs = sorted(common_substrings, key=len, reverse=True)
        new_match_seq = ordered_cs[0]
        v_seq_parts = v_seq.split(new_match_seq)[:-1]
        v_seq_trimmed = ''.join(v_seq_parts)
    return v_seq_trimmed

#Finds motif in J sequence to keep the following sequence
def j_trimming(j_seq):
    '''
    Input: j_seq = J sequence
    Output: j_seq_trimmed which is the sequence following F/W motif ((F/W)G(N)G)
    '''
    j_dict = {'F':2, 'W':2, 'G':1}
    value_list = []
    for n in j_seq:
        try:
            value = j_dict[n]
            value_list.append(value)
        except:
            value = 0
            value_list.append(value)
    value_string = ''.join(str(num) for num in value_list)
    match = SequenceMatcher(None, value_string, '2101').find_longest_match()
    string_start = match.a + 1
    j_seq_trimmed = j_seq[string_start:]
    return j_seq_trimmed

#Applies v_trimming and j_trimming functions to create full region sequence
def full_region_gen(vj):
    '''
    Input: vj file containing V and J AA sequences (output from id_match function)
    Output: df updated with V-CDR3-J full AA sequence column
    '''
    v_trim_list = []
    j_trim_list = []
    for i, row in vj.iterrows():
        v = row['V_AA_Seq']
        c = row['CDR3']
        j = row['J_AA_Seq']
        v_trimmed = v_trimming(v,c)
        j_trimmed = j_trimming(j)
        v_trim_list.append(v_trimmed)
        j_trim_list.append(j_trimmed)
    vj['V_AA_Trimmed'] = v_trim_list
    vj['J_AA_Trimmed'] = j_trim_list
    vj['V_CDR3_J_Sequence'] = vj['V_AA_Trimmed'] + vj['CDR3'] + vj['J_AA_Trimmed']
    return vj

### MAIN FUNCTION
def vdj_region_gen(vj_file, v_file, j_file, output):
    '''
    Input: vj file, v data file, j data file, output csv file name
    Output: vj file with full V-CDR3-J AA sequences generated
    '''
    #Read in files
    vj_df = pd.read_csv(vj_file)
    v_df = pd.read_csv(v_file)
    j_df = pd.read_csv(j_file)
    #Convert accession numbers into strings
    v_df['V.Accession.Number'] = v_df['V.Accession.Number'].astype(str)
    j_df['J.Accession.Number'] = j_df['J.Accession.Number'].astype(str)
    #Remove '>' from accession number for v
    v_an = []
    for x in v_df['V.Accession.Number']:
        accession_number = x.split('>')[1]
        v_an.append(accession_number)
    v_df['AN'] = v_an
    #Remove '>' from accession number for j
    j_an = []
    for x in j_df['J.Accession.Number']:
        accession_number = x.split('>')[1]
        j_an.append(accession_number)
    j_df['AN'] = j_an
    #Create VID/JID columns
    combine_columns(v_df, 'AN', 'V.Name')
    combine_columns(j_df, 'AN', 'J.Name')
    #Add corresponding V/J AA sequences to the main df
    vj_with_seq = id_match(vj_df, v_df, j_df)
    #Remove rows where AA sequence for VID could not be found
    vj_with_seq2 = vj_with_seq[(vj_with_seq['V_AA_Seq'] != 'NA') & (vj_with_seq['J_AA_Seq'] != 'NA')]
    #Generate column containing full V-CDR3-J AA sequences
    final_df = full_region_gen(vj_with_seq2)
    #Save output to csv
    final_df.to_csv(output, header=True, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates full V-CDR3-J AA sequences')
    parser.add_argument('vj_file', type=str, help='CSV file with VID, JID, CDR3 sequence, etc.')
    parser.add_argument('-v', '--v_file', type=str, help='V_data file containing corresponding V AA sequences')
    parser.add_argument('-j', '--j_file', type=str, help='J_data file containing corresponding J AA sequences')
    parser.add_argument('-o', '--output', default='vdj_output.csv', help='Assign name to output CSV file')
    args = parser.parse_args()
    vdj_region_gen(vj_file=args.vj_file, v_file=args.v_file, j_file=args.j_file, output=args.output)
