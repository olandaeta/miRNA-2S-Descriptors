# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 10:10:51 2018
version 02
@author: Omar Gabriel Landaeta Herrera
To run this script I use the python Shell prompt where I call the descriptors
function like this:
    descriptors()
    
    This will print the prompt on the screen like this:
            
    Path to file and name: miRNA_Dataset.xlsx
    File type i.e.(excel or csv): excel
    number or name ofthe sheet: 0
    Columns names that will be loaded: Y1_3,Y4_6,Expression level,mRNA SS
    minimun descriptor size: 2
    maximun descriptor size: 3
    Path to file and name: ./results/test_2-3
    Type of the ouput file i.e.(excel or csv): excel

"""

import pandas as pd
from itertools import product


def read_file(l,sn,columns_names,f_type):
    """Function that take as input the file location(l), sheet name(sn), columns 
    names(columns_names) and the file type(f_type) and return a dataframe with
    the columns requested"""
    if f_type == 'excel':
        file_df_original = pd.read_excel(r''+l+'',sheet_name=sn)
    elif f_type == 'csv':
        file_df_original = pd.read_csv(r''+l+'')
    
    miRNA_df = file_df_original[columns_names].copy()
    return miRNA_df

def descriptor_maker(n=2,m=0):
    """Function that take two inputs, the minimun window size and a maximun window
    siz and return a descriptors list. The default window size is 2.    """     
    lista_array= []
    if (n == 2 and m == 0) or (n == m):        
        lista_array = map(''.join,product('.()', repeat=n))
    elif n > 2 and m == 0:
        for i in range(2,n+1):
            lista_array += map(''.join,product('.()', repeat=i))
    else:
        for i in range(n,m+1):
            lista_array += map(''.join,product('.()', repeat=i))
    return list(lista_array)

def horspool_P(text, pattern):
    """Function that implement the Boyer Moore Horspool algorithm for exact string
    matching. It has two inputs, the characters chain and the pattern """
    m = len(pattern)
    n = len(text)    
    ocurr = 0
    shift = {}
    for i in range(m):
        shift[pattern[i]] = m
    for i in range(m-1):
        shift[pattern[i]] = m-1-i
    k = m-1    
    while k < n:
        j = m-1
        i = k        
        while j >= 0 and text[i] == pattern[j]:
            j -= 1
            i -= 1
        if j == -1:
            ocurr += 1
        k += shift.get(text[k], m)
    return ocurr

def descriptors():
    """The main function that will request user inputs to process a file and 
    return a process file"""
    location = input('Path to file and name: ')
    file_type = input("File type i.e.(excel or csv): ")
    if file_type == 'excel':
        sheet_name = int(input('number or name ofthe sheet: '))
    columns_names = input('Columns names that will be loaded: ').split(',')
        
    min_dex = int(input('minimun descriptor size: '))
    if min_dex < 2:
        min_dex = 2
    max_dex = int(input('maximun descriptor size: '))
    if max_dex == None:
        max_dex = 0
        
    output_file_name = input('Path to file and name: ')
    output_file_type = input('Type of the ouput file i.e.(excel or csv): ')
    
    ss_df = read_file(location,sheet_name,columns_names,file_type)    
    dex = descriptor_maker(min_dex,max_dex)
    dex_length = len(dex)
    rows_size = ss_df.shape[0]
    header = [i for i in ss_df.columns[0:4]]
    
    for i in range(dex_length):
        ss_df['Fr'+str(len(dex[i]))+'_'+str(i+1)+'_'+dex[i]] = 0
        for j in range(rows_size):
            ss_df.iat[j,i+4]=horspool_P(ss_df.iat[j,3], dex[i])
    
    for v in range(dex_length):
        if ss_df.iloc[:,v+4].sum() > 0:
            header.append(ss_df.columns[v+4])
             
    ss_df_output = ss_df[header]
    
    if output_file_type == 'excel':
        ss_df_output.to_excel(output_file_name+'.xlsx', index = False)
    elif output_file_type == 'csv':
        ss_df_output.to_csv(output_file_name+'.csv', index = False)
    
    
