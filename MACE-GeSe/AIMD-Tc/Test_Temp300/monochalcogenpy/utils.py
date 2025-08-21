import numpy as np

def fill_abc(string):
    assert len(string) ==2
    missing_char = [char for char in 'abc' if char not in string][0]
    return string + missing_char

def sort_axis_indices(orientation, ref_orientation):
    #goal - turn indices currently along abc into current orientation
    #understand this as (the first part gives mapping from 'abc' to ref_orientation.)
    #then the second maps orientation to bc 
    #so in effect we are mapping from orientation to abc obeying order from ref_orientation
    return sort_string_indices('abc',ref_orientation)[sort_string_indices(orientation, 'abc')]

def sort_string_indices(string, ref_string):
    """
    Sorts the indices of string1 based on the sorted order of characters in string2.
    
    Args:
    string: The string whose indices need to be sorted.
    ref_string: The string based on whose sorted order the indices will be sorted.
    
    Returns:
    A list of integers representing the sorted indices of string.
    """
    indices = {c: i for i, c in enumerate(list(string))}
    sort_index = np.array([indices[char] for char in list(ref_string)])
    return sort_index
